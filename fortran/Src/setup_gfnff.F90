  subroutine setup_gfnff
!
!  Initialisation routine for GFNFF in GULP
!
!   2/21 GFNFF reference configuration added
!   2/21 Option to increment pi system added
!   2/21 Huckel solve separated into a subroutine
!   3/21 Use of reference geometry added
!   3/21 Correction to extra torsion for sp3-sp3 case
!   3/21 Change to unit conversion of torsions to avoid problem where ntors is zero
!   3/21 pibond initialised to zero to avoid problem for N in pi systems
!   3/21 Reference charges added
!   3/21 Corrections to torsions added from version 6.4.0 - fij mods moved outside of k/l loop
!   3/21 Overflow due to exponential of qfac trapped
!   3/21 cuts set to a small value to avoid problems with electrostatics
!   3/21 Trap for shell models added
!   3/21 Topological distances corrected for PBC and 1+12 vs 1d+12 error
!   4/21 Parallel modifications added
!   4/21 K point setup added
!   5/21 pibond correction to bonding now uses a finite tolerance to avoid rounding errors
!   5/21 Setting of ngfnff_current_cfg added
!   5/21 Wolf sum in topo added as an option
!   5/21 Exclusion of bonds added as an option
!   6/21 Use of xclat/yclat/zclat changed to xalat/yalat/zalat for benefit of MD
!   7/21 In non-MD case xclat/yclat/zclat are now used to allow for symmetry
!   7/21 Atomic number check added
!   7/21 lgfnff_fragment_bond added
!   7/21 Setting of gfnff_ref_q added
!   7/21 Coordinates changed to always use xclat/yclat/zclat as MD should matter for setup
!  10/21 gfnff_ref arrays moved to configurations
!  10/21 Renamed from gfnff_setup to setup_gfnff
!   1/22 Cell list updated if reference cell used
!  11/23 Trap for zero atoms added
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2023
!
!  Julian Gale, Curtin University, November 2023
!
  use datatypes
  use configurations,  only : lgfnff_ref_rv, lgfnff_ref_crd
  use configurations,  only : gfnff_ref_rv, gfnff_ref_crd, gfnff_ref_q
  use control,         only : keyword
  use current,         only : numat, nat, qf, ndim, ncf, maxat
  use current,         only : nftype, nasym, nrela2f
  use current,         only : xalat, yalat, zalat
  use current,         only : xclat, yclat, zclat
  use current,         only : r1x, r1y, r1z, rv
  use current,         only : r2x, r2y, r2z, kv
  use current,         only : r3x, r3y, r3z
  use element,         only : atsym
  use gulp_gfnff
  use m_gfnff_c6
  use m_gfnff_nbr
  use molecule,        only : nnobo, nobond, nobotyp
  use parallel,        only : ioproc
  use shells,          only : cuts, nshell
  use times,           only : tgfnffsetup
  implicit none
!
!  Local variables
!
  character(len=2), dimension(:),  allocatable, save :: atsym_ele
  integer(i4)                                        :: i
  integer(i4)                                        :: ii
  integer(i4)                                        :: maxelein
  integer(i4)                                        :: status
  logical                                            :: lbondsok
  logical                                            :: lverbose
  real(dp),    dimension(:),       allocatable, save :: qsave
  real(dp)                                           :: g_cpu_time
  real(dp)                                           :: rvsave(3,3)
  real(dp)                                           :: t1
  real(dp)                                           :: t2
  real(dp),    dimension(:,:),     allocatable, save :: xyzsave
!
  t1 = g_cpu_time()
!
  lverbose = ((index(keyword,'verb').ne.0.or.index(keyword,'gver').ne.0).and.ioproc)
!
!  Store the configuration whose parameters are being setup
!
  ngfnff_current_cfg = ncf
!
!  Prevent runs with shells for now
!
  if (nshell.gt.0) then
    call outerror('GFNFF is current incompatible with the shell model',0_i4)
    call stopnow('setup_gfnff')
  endif
!
!  Check that there are atoms in the structure otherwise the setup will fail!
!
  if (numat.eq.0) then
    call outerror('GFNFF parameter setup requires a number of atoms > 0',0_i4)
    call stopnow('setup_gfnff')
  endif
!
!  Allocate local memory
!
  allocate(atsym_ele(numat),stat=status)
  if (status/=0) call outofmemory('setup_gfnff','atsym_ele')
!
!  Setup up atomic symbols
!
  do i = 1,numat
    atsym_ele(i) = atsym(nat(i))
  enddo
!
!  If reference geometry hasn't been specified then save the current one as the reference
!  If it has been specified then switch in here
!
  if (ndim.gt.0) then
    rvsave(1:3,1:3) = rv(1:3,1:3)
    if (lgfnff_ref_rv(ncf)) then
      rv(1:3,1:3) = gfnff_ref_rv(1:3,1:3,ncf)
      r1x = rv(1,1)
      r1y = rv(2,1)
      r1z = rv(3,1)
      r2x = rv(1,2)
      r2y = rv(2,2)
      r2z = rv(3,2)
      r3x = rv(1,3)
      r3y = rv(2,3)
      r3z = rv(3,3)
!
!  If cell is changed then update vector list
!
      call rlist2
    else
      gfnff_ref_rv(1:3,1:3,ncf) = rv(1:3,1:3)
    endif
  endif
  if (lgfnff_ref_crd(ncf)) then
    allocate(xyzsave(3,numat),stat=status)
    if (status/=0) call outofmemory('setup_gfnff','xyzsave')
    allocate(qsave(numat),stat=status)
    if (status/=0) call outofmemory('setup_gfnff','qsave')
!
    do i = 1,numat
      xyzsave(1,i) = xclat(i)
      xyzsave(2,i) = yclat(i)
      xyzsave(3,i) = zclat(i)
    enddo
    qsave(1:numat) = qf(1:numat)
  endif
  if (lgfnff_ref_crd(ncf)) then
    do i = 1,numat
      xclat(i) = gfnff_ref_crd(1,i,ncf)
      yclat(i) = gfnff_ref_crd(2,i,ncf)
      zclat(i) = gfnff_ref_crd(3,i,ncf)
    enddo
    do i = 1,nasym
      ii = nrela2f(i)
      xalat(i) = xclat(ii)
      yalat(i) = yclat(ii)
      zalat(i) = zclat(ii)
    enddo
    qf(1:numat) = gfnff_ref_q(1:numat,ncf)
  else
    do i = 1,numat
      gfnff_ref_crd(1,i,ncf) = xclat(i)
      gfnff_ref_crd(2,i,ncf) = yclat(i)
      gfnff_ref_crd(3,i,ncf) = zclat(i)
    enddo
    gfnff_ref_q(1:numat,ncf) = qf(1:numat)
  endif
!################################################################################
!  Pre-loop setup                                                               #
!################################################################################
!
!  Build atom in cell coordinates
!
  call getinbox
!
!  Set cuts to a small value otherwise electrostatics can go wrong for distances below this value
!
  cuts = 0.01_dp
!################################################################################
!  Build the neighbour lists                                                    #
!################################################################################
#ifndef NOGFNFF
!
!  Find the maximum cutoff for the neighbour list
!
  call pgfnff_get_max_cutoff(numat,nat,gfnff_cut2_nbr)
!
!  Build the neighbour list
!
  call gfnff_getnbr(gfnff_cut2_nbr,.true.)
!
!  Check that all atoms have neighbours where they should have
!
  call pgfnff_check_atom_with_nobonds(numat,nat,nnbr,lbondsok)
!
!  If bonding is not OK then try a larger cutoff radius
!
  if (.not.lbondsok) then
    gfnff_cut2_nbr = 2.0_dp*gfnff_cut2_nbr
    call gfnff_getnbr(gfnff_cut2_nbr,.true.)
  endif
#endif
!**********************************************
!  Call routine to generate GFNFF parameters  *
!**********************************************
#ifndef NOGFNFF
  call pgfnff_pargen(ndim,kv,numat,nat,nftype,xclat,yclat,zclat,qf,nfrag,nfraglist,qfrag, &
                     nnobo,nobond,nobotyp,maxnbr,nnbr,nbrno,ncnbr,rnbr,xnbr,ynbr,znbr,nnbr_bond, &
                     nbrno_bond,ncnbr_bond,rbnbr,xbnbr,ybnbr,zbnbr,lverbose)
#else
  call outerror('GFNFF cannot be setup without linking the pGFNFF library',0_i4)
  call stopnow('setup_gfnff')
#endif
!
!  Ensure that arrays are large enough for numat + nfrag for use later in EEM
!
  if (numat+nfrag.gt.maxat) then
    maxat = numat + nfrag
    call changemaxat
  endif
!
!  Ensure arrays that depend on the number of atoms have the right size
!
  maxat_gfnff = numat
  call changemaxatgfnff
!
!  Reset geometry if a reference one was used
!
  if (ndim.gt.0) then
    if (lgfnff_ref_rv(ncf)) then
      rv(1:3,1:3) = rvsave(1:3,1:3)
      r1x = rv(1,1)
      r1y = rv(2,1)
      r1z = rv(3,1)
      r2x = rv(1,2)
      r2y = rv(2,2)
      r2z = rv(3,2)
      r3x = rv(1,3)
      r3y = rv(2,3)
      r3z = rv(3,3)
!
!  If cell is changed then update vector list
!
      call rlist2
    endif
  endif
  if (lgfnff_ref_crd(ncf)) then
    do i = 1,numat
      xclat(i) = xyzsave(1,i)
      yclat(i) = xyzsave(2,i)
      zclat(i) = xyzsave(3,i)
    enddo
    qf(1:numat) = qsave(1:numat)
!
    deallocate(qsave,stat=status)
    if (status/=0) call deallocate_error('setup_gfnff','qsave')
    deallocate(xyzsave,stat=status)
    if (status/=0) call deallocate_error('setup_gfnff','xyzsave')
  endif
!***************************************
!  Set up nb3atm and associated lists  *
!***************************************
  call gfnff_setb3atm(ndim,numat,nnbr_bond,maxnbr,nbrno_bond,xbnbr,ybnbr,zbnbr)
#ifndef NOGFNFF
!*********************************************************
!  Transfer GFNFF force field parameters to GULP arrays  *
!*********************************************************
!
!  Find the maximum element number supported by the library
!
  call pgfnff_get_maxele(max_gfnff_ele)
!
  call changemaxgfnffele
!
!  Find the maximum atomic number in the actual system
!
  maxelein = 0
  do i = 1,numat
    maxelein = max(maxelein,nat(i))
  enddo
  maxc6ele = maxelein
!
!  General parameters
!
  call pgfnff_get_maxref(maxc6ref)
  call pgfnff_get_cnpar(gfnff_kn_cn,gfnff_kn_hb,gfnff_cnmax)
!
!  Coordination number thresholds
!
  call pgfnff_get_coordination_thresholds(gfnff_cnthr,gfnff_cnhbthr,gfnff_cnerfcut_hb,gfnff_cnerfcut_cn)
!
!  Check array sizes for dispersion
!
  call changemaxd4c6
!
!  Dispersion parameters
!
  call pgfnff_get_dispersion_threshold(gfnff_dispthr)
  t_dispthr = gfnff_dispthr*gfnff_taper
  call pgfnff_get_dispersion(numat,nat,maxelein,maxc6ref,d4_nref,d4_cn,d4_c6,d4_c9,d4_zeta_c6,d4_r0,d4_sqrtZr4r2)
!
!  Covalent radii
!
  call pgfnff_get_covalent_radii(maxelein,gfnff_rcov)
!
!  General radii
!
  call pgfnff_get_general_radii(maxelein,gfnff_rad)
!
!  Coordination number parameters for radii
!
  call pgfnff_get_cn_radii(maxelein,gfnff_rad_cn)
!
!  EEQ parameters
!
  call pgfnff_get_eeq(numat,gfnff_eeq_alp,gfnff_eeq_chi,gfnff_eeq_gam,gfnff_eeq_cnf)
!
!  Repulsion parameters
!
  call pgfnff_get_repulsion_threshold(gfnff_repthr)
  t_repthr = gfnff_repthr*gfnff_taper
  call pgfnff_get_repulsion_scale(gfnff_repscale_b,gfnff_repscale_n,gfnff_repscale_13,gfnff_repscale_14)
  call pgfnff_get_repulsion(numat,maxelein,gfnff_repulsion_a,gfnff_repulsion_z,gfnff_repulsion_p)
!
!  Bond parameters
!
  call pgfnff_get_bond_parameters(numat,maxnbr,nnbr_bond,par_gfnff_bond)
  call pgfnff_get_bond_scale(gfnff_bondscale)
!
!  Angle bends
!
  call pgfnff_get_number_of_angles(n_gfnff_angles)
  max_gfnff_angles = n_gfnff_angles
  call changemaxgfnffangles
  call pgfnff_get_angles(max_gfnff_angles,n_gfnff_angleptr,par_gfnff_angle,gfnff_angle_damp)
!
!  Torsions
!
  call pgfnff_get_number_of_torsions(n_gfnff_torsions)
  max_gfnff_torsions = n_gfnff_torsions
  call changemaxgfnfftorsions
  call pgfnff_get_torsions(max_gfnff_torsions,n_gfnff_torsionptr,par_gfnff_torsion,gfnff_torsion_damp)
!
!  Hydrogen bond thresholds
!
  call pgfnff_get_hydrogen_bond_thresholds(gfnff_hbthr1,gfnff_hbthr2)
!
!  Hydrogen bond cutoffs
!
  call pgfnff_get_hydrogen_bond_cutoffs(gfnff_hb_a_cut,gfnff_hb_long_cut,gfnff_hb_nb_cut, &
                                        gfnff_hb_short_cut,gfnff_hb_alp)
!
!  Halogen bond cutoffs
!
  call pgfnff_get_halogen_bond_cutoffs(gfnff_xb_a_cut,gfnff_xb_short_cut,gfnff_xb_long_cut)
!
!  Hydrogen bond potential AB atoms
!
  call pgfnff_get_number_hydrogen_bond_AB(n_gfnff_hb_AB)
  call pgfnff_get_hydrogen_bond_AB(numat,n_gfnff_hb_AB,n_gfnff_hb_ABptr,gfnff_hb_ABq, &
                                   gfnff_hb_acid,gfnff_hb_base)
!
!  Hydrogen bond potential H atoms
!
  call pgfnff_get_number_hydrogen_bond_H(n_gfnff_hb_H)
  if (n_gfnff_hb_H.gt.0) then
    if (n_gfnff_hb_H.gt.maxathbH) then
      maxathbH = n_gfnff_hb_H
      call changemaxathbH
    endif
    call pgfnff_get_hydrogen_bond_H(numat,n_gfnff_hb_H,n_gfnff_hb_Hptr,n_gfnff_hb_Hrptr)
  endif
!
!  Hydrogen bond scale factors
!
  call pgfnff_get_hydrogen_bond_scale(gfnff_hb_scale_gen,gfnff_hb_scale_coh)
!
!  Hydrogen bond parameters
!
  call pgfnff_get_hydrogen_bond_parameters(gfnff_bend_hb,gfnff_tors_hb,gfnff_hbabmix)
!
!  Halogen bond potential AB atoms
!
  call pgfnff_get_number_halogen_bond_AB(n_gfnff_xb_AB)
  if (n_gfnff_xb_AB.gt.0) then
    if (n_gfnff_xb_AB.gt.maxatxbAB) then
      maxatxbAB = n_gfnff_xb_AB
      call changemaxatxbAB
    endif
    call pgfnff_get_halogen_bond_AB(numat,n_gfnff_xb_AB,n_gfnff_xb_ABptr,gfnff_xb_ABq)
  endif
!
!  Halogen bond scale factors
!
  call pgfnff_get_halogen_bond_scale(max_gfnff_ele,gfnff_xb_scale)
#endif
!
!  Deallocate local memory
!
  deallocate(atsym_ele,stat=status)
  if (status/=0) call deallocate_error('setup_gfnff','atsym_ele')
!
  t2 = g_cpu_time()
  tgfnffsetup = tgfnffsetup + t2 - t1
!
  end subroutine setup_gfnff
