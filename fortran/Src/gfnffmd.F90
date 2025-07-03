  subroutine gfnffmd(egfnff,lgrad1)
!
!  Calculates the energy and derivatives for the GFNFF force field.
!  Energy and first derivatives only
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
!
!  On exit :
!
!  egfnff          = the value of the energy contribution
!
!   9/20 Created
!  10/20 First derivatives added
!  10/20 Handling of surfaces and regions added
!  11/20 Energy/force only version created
!   2/21 Cut no longer passed to coordination number routines
!   4/21 Parallelisation added
!   4/21 Energy format for output changed
!   6/21 Damped ATM dispersion for GFNFF added
!   7/21 Hydrogen bond setup call added
!   7/21 Neighbour list now regenerated 
!  10/21 Generic damping routine introduced
!   2/22 Changes to parallel handling of the energy output
!   2/22 Arguments added to bond_hb_set_AHB
!   3/22 egfnff_tot now set for a single processor
!   5/22 Call to gfnff_setc6 moved
!   5/22 C6 derivatives rearranged for improved performance
!   7/23 nregionnocfg added
!   7/23 lupdateinbox added
!   9/23 Energy sum logic corrected
!   9/23 gfnff_cn called with lgrad1 true if ldcharge is true
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
!  Julian Gale, CIC, Curtin University, September 2023
!
  use datatypes
  use configurations, only : nregions, lsliceatom, lregionrigid
  use control,        only : keyword, lseok, ldcharge
  use current
  use derivatives
  use EDIPdata
  use energies,       only : eattach, esregion12, esregion2
  use energies,       only : eregion2region, siteenergy
  use g_constants,    only : pi
  use gulp_gfnff
  use iochannels
  use m_gfnff_c6
  use m_gfnff_nbr
  use m_gfnff_nbr3
  use m_gfnff_pairs
  use m_strain,       only : real1strterm
  use m_strain,       only : realstrterms
  use neighbours
  use numbers,        only : third
  use parallel,       only : ioproc, nprocs, procid
  use spatial,        only : lupdateinbox
  use spatialbo
  use symmetry,       only : lstr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(out)                         :: egfnff
  logical,     intent(in)                          :: lgrad1
!
!  Local variables
!
  integer(i4)                                      :: hbA
  integer(i4)                                      :: hbH
  integer(i4)                                      :: i
  integer(i4)                                      :: ia
  integer(i4)                                      :: imin
  integer(i4)                                      :: ind
  integer(i4)                                      :: indn
  integer(i4)                                      :: indij
  integer(i4)                                      :: indil
  integer(i4)                                      :: indjl
  integer(i4)                                      :: j
  integer(i4)                                      :: jmax
  integer(i4)                                      :: k
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: l
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: natl
  integer(i4)                                      :: ni
  integer(i4)                                      :: nij
  integer(i4)                                      :: nik
  integer(i4)                                      :: nil
  integer(i4)                                      :: nj
  integer(i4)                                      :: njk
  integer(i4)                                      :: nk
  integer(i4)                                      :: nkl
  integer(i4)                                      :: nl
  integer(i4)                                      :: np
  integer(i4)                                      :: nregioni           ! Region number for i
  integer(i4)                                      :: nregionj           ! Region number for j
  integer(i4)                                      :: nregionk           ! Region number for k
  integer(i4)                                      :: nregionl           ! Region number for l
  integer(i4)                                      :: nrot
  integer(i4)                                      :: nshell
  integer(i4)                                      :: status
  logical                                          :: lattach
  logical                                          :: lbonded
  logical                                          :: ldohbcn
  logical                                          :: linear
  logical                                          :: lreg12
  logical                                          :: lreg2one
  logical                                          :: lreg2pair
  logical                                          :: lreg2qtet
  logical                                          :: lreg2trio
  logical                                          :: lslicei
  logical                                          :: lslicej
  logical                                          :: lslicek
  logical                                          :: lslicel
  logical                                          :: lverbose
  real(dp)                                         :: c6
  real(dp)                                         :: dc6dcni
  real(dp)                                         :: dc6dcnj
  real(dp)                                         :: d2c6dcni2
  real(dp)                                         :: d2c6dcnidcnj
  real(dp)                                         :: d2c6dcnj2
  real(dp)                                         :: c9
  real(dp),    dimension(:),     allocatable, save :: cn                ! log of coordination number
  real(dp),    dimension(:),     allocatable, save :: cn_hb             ! coordination number for HB only
  real(dp)                                         :: cos1
  real(dp)                                         :: dcos1drij
  real(dp)                                         :: dcos1dril
  real(dp)                                         :: dcos1drjl
  real(dp)                                         :: cos2
  real(dp)                                         :: dcos2drij
  real(dp)                                         :: dcos2dril
  real(dp)                                         :: dcos2drjl
  real(dp)                                         :: cos3
  real(dp)                                         :: dcos3drij
  real(dp)                                         :: dcos3dril
  real(dp)                                         :: dcos3drjl
  real(dp)                                         :: cosang
  real(dp)                                         :: cosfct
  real(dp)                                         :: cosphi
  real(dp)                                         :: cosphi0
  real(dp)                                         :: costheta0
  real(dp)                                         :: cutdr
  real(dp)                                         :: damp
  real(dp)                                         :: dampij
  real(dp)                                         :: dampik
  real(dp)                                         :: dampjk
  real(dp)                                         :: dampjl
  real(dp)                                         :: ddampij
  real(dp)                                         :: ddampik
  real(dp)                                         :: ddampjk
  real(dp)                                         :: ddampjl
  real(dp)                                         :: d2dampij
  real(dp)                                         :: d2dampik
  real(dp)                                         :: d2dampjk
  real(dp)                                         :: d2dampjl
  real(dp)                                         :: dcosdrij
  real(dp)                                         :: dcosdrik
  real(dp)                                         :: dcosdrjk
  real(dp)                                         :: disp
  real(dp)                                         :: ddispdr
  real(dp),    dimension(:),     allocatable, save :: dcn               ! First derivatives of log of coordination number w.r.t. coordination number
  real(dp),    dimension(:),     allocatable, save :: d2cn              ! Second derivatives of log of coordination number w.r.t. coordination number
  real(dp)                                         :: dr0dcni
  real(dp)                                         :: dr0dcnj
  real(dp)                                         :: dtheta
  real(dp)                                         :: dr2ds(6)
  real(dp)                                         :: d2r2dx2(6)
  real(dp)                                         :: d2r2ds2(6,6)
  real(dp)                                         :: d2r2dsdx(6,3)
  real(dp)                                         :: dr2ds6(6,6)
  real(dp)                                         :: d2r2dx26(6,6)
  real(dp)                                         :: d2r2ds26(6,6,6)
  real(dp)                                         :: d2r2dsdx6(6,3,6)
  real(dp)                                         :: e1d(6)
  real(dp)                                         :: e2d(21)
  real(dp)                                         :: eangle            ! Angle bending energy
  real(dp)                                         :: ebatm             ! Non-bond triples energy
  real(dp)                                         :: ebond             ! Bonding energy
  real(dp)                                         :: ecoulomb          ! Coulomb energy
  real(dp)                                         :: edispersion       ! Dispersion energy
  real(dp)                                         :: ehb               ! Hydrogen bonding energy
  real(dp)                                         :: erepulsion        ! Repulsive energy
  real(dp)                                         :: etorsion          ! Torsional energy
  real(dp)                                         :: exb               ! XB energy
  real(dp)                                         :: eterm6th
  real(dp)                                         :: egfnff_tot
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: dr
  real(dp)                                         :: deijdcn           ! Derivative of eij w.r.t. coordination number
  real(dp)                                         :: deijdcn_sum       ! Derivative of eij w.r.t. coordination number - sum of terms
  real(dp)                                         :: d2eijdcn2         ! Second derivative of eij w.r.t. coordination number (dummy)
  real(dp)                                         :: deijddr           ! Derivative of eij w.r.t. dr
  real(dp)                                         :: deijdr            ! (1/r)(dEij/dr)
  real(dp)                                         :: eij
  real(dp)                                         :: eijk
  real(dp)                                         :: eijknodamp
  real(dp)                                         :: deijkdcos
  real(dp)                                         :: deijkdr
  real(dp)                                         :: deijkdrij
  real(dp)                                         :: deijkdril
  real(dp)                                         :: deijkdrjl
  real(dp)                                         :: eijkl
  real(dp)                                         :: eijklnod
  real(dp)                                         :: deijkldphi
  real(dp)                                         :: deijkldrij
  real(dp)                                         :: deijkldrjk
  real(dp)                                         :: deijkldrjl
  real(dp)                                         :: fct
  real(dp)                                         :: ff23
  real(dp)                                         :: fil
  real(dp),    dimension(:,:),   allocatable, save :: gw
  real(dp),    dimension(:,:),   allocatable, save :: dgwdcn
  real(dp),    dimension(:,:),   allocatable, save :: d2gwdcn2
  real(dp)                                         :: phi
  real(dp)                                         :: dphi
  real(dp)                                         :: fphi
  real(dp)                                         :: phi0
  real(dp)                                         :: phi1dx(3,3)
  real(dp)                                         :: phi1ds(6)
  real(dp)                                         :: phi2ds(6,6)
  real(dp)                                         :: phi2dsdx(6,3,3)
  real(dp)                                         :: phi2dx(3,3,6)
  real(dp)                                         :: r
  real(dp)                                         :: r0
  real(dp)                                         :: r03
  real(dp)                                         :: r04
  real(dp)                                         :: r15
  real(dp)                                         :: r2
  real(dp)                                         :: r2ij
  real(dp)                                         :: r2ik
  real(dp)                                         :: r2jk
  real(dp)                                         :: r4r2ij
  real(dp)                                         :: rcut2ij
  real(dp)                                         :: rcut2ik
  real(dp)                                         :: rcut2jk
  real(dp)                                         :: rcut2jl
  real(dp)                                         :: repab
  real(dp)                                         :: rep_rho
  real(dp)                                         :: rep_rholoc
  real(dp)                                         :: rij0
  real(dp)                                         :: ril0
  real(dp)                                         :: rjl0
  real(dp)                                         :: rij03
  real(dp)                                         :: ril03
  real(dp)                                         :: rjl03
  real(dp)                                         :: rij
  real(dp)                                         :: rij2
  real(dp)                                         :: rrijk9
  real(dp)                                         :: rik
  real(dp)                                         :: rik2
  real(dp)                                         :: ril
  real(dp)                                         :: ril2
  real(dp)                                         :: rjk
  real(dp)                                         :: rjk2
  real(dp)                                         :: rjl
  real(dp)                                         :: rjl2
  real(dp)                                         :: rkl
  real(dp)                                         :: rkang
  real(dp)                                         :: rn
  real(dp)                                         :: rrij
  real(dp)                                         :: rril
  real(dp)                                         :: rrik
  real(dp)                                         :: rrjl
  real(dp)                                         :: rrij3
  real(dp)                                         :: rril3
  real(dp)                                         :: rrjl3
  real(dp)                                         :: sinphi
  real(dp)                                         :: sum1(9)
  real(dp)                                         :: sum0(9)
  real(dp)                                         :: t            ! Taper function
  real(dp)                                         :: dtdr         ! Taper function first derivative w.r.t. r divided by r
  real(dp)                                         :: dtdr2        ! Taper function first derivative w.r.t. r2
  real(dp)                                         :: d2tdr22      ! Taper function second derivative w.r.t. r2
  real(dp)                                         :: d3tdr23      ! Taper function third derivative w.r.t. r2
  real(dp)                                         :: p2
  real(dp)                                         :: p3
  real(dp)                                         :: trm6
  real(dp)                                         :: trm8
  real(dp)                                         :: tmp
  real(dp)                                         :: tmprr
  real(dp)                                         :: theta
  real(dp)                                         :: theta0
  real(dp)                                         :: time1
  real(dp)                                         :: time2
  real(dp)                                         :: time3
  real(dp)                                         :: time4
  real(dp)                                         :: vij(3)
  real(dp)                                         :: vik(3)
  real(dp)                                         :: vil(3)
  real(dp)                                         :: vjk(3)
  real(dp)                                         :: vjl(3)
  real(dp)                                         :: vkl(3)
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xcom(6)
  real(dp)                                         :: ycom(6)
  real(dp)                                         :: zcom(6)
  real(dp)                                         :: xv4(6)
  real(dp)                                         :: yv4(6)
  real(dp)                                         :: zv4(6)
  real(dp)                                         :: zetaprod
#ifdef TRACE
  call trace_in('egfnff_energy')
#endif
!
  time1 = g_cpu_time()
!
!  Initialise energies
!
  egfnff = 0.0_dp
  eangle = 0.0_dp
  ebatm = 0.0_dp
  ebond = 0.0_dp
  ecoulomb = 0.0_dp
  edispersion = 0.0_dp
  ehb = 0.0_dp
  erepulsion = 0.0_dp
  etorsion = 0.0_dp
  exb = 0.0_dp
!
  d2eijdcn2 = 0.0_dp
!
  lverbose = ((index(keyword,'verb').ne.0.or.index(keyword,'gver').ne.0).and.ioproc)
!
  if (lgrad1) then
!
!  Initialise centre of mass arrays to zero
!
    xcom(1:6) = 0.0_dp
    ycom(1:6) = 0.0_dp
    zcom(1:6) = 0.0_dp
  endif
!
!  Set flag as to whether hydrogen bond coordination number is needed
!
  ldohbcn = (n_gfnff_hb_AB.gt.0.and.n_gfnff_hb_H.gt.0)
!
!  Allocate local memory 
!
  allocate(cn(numat),stat=status)
  if (status/=0) call outofmemory('gfnff_energy','cn')
  allocate(dcn(numat),stat=status)
  if (status/=0) call outofmemory('gfnff_energy','dcn')
  allocate(d2cn(1),stat=status)
  if (status/=0) call outofmemory('gfnff_energy','d2cn')
  allocate(gw(maxc6ref,numat),stat=status)
  if (status/=0) call outofmemory('gfnff_energy','gw')
!
  if (lgrad1) then
    allocate(dgwdcn(maxc6ref,numat),stat=status)
    if (status/=0) call outofmemory('gfnff_energy','dgwdcn')
  else
    allocate(dgwdcn(maxc6ref,1),stat=status)
    if (status/=0) call outofmemory('gfnff_energy','dgwdcn')
  endif
  allocate(d2gwdcn2(maxc6ref,1),stat=status)
  if (status/=0) call outofmemory('gfnff_energy','d2gwdcn2')
!
  if (ldohbcn) then
    allocate(cn_hb(numat),stat=status)
    if (status/=0) call outofmemory('gfnff_energy','cn_hb')
  endif
!
!  Set coordinates in the box
!
  if (lupdateinbox) call getinbox
!********************************************************
!  Update the neighbour lists for the current geometry  *
!********************************************************
  if (lverbose) time3 = g_cpu_time()
!
!  Re-build the neighbour list
!
  call gfnff_getnbr(gfnff_cut2_nbr,.false.)
  call gfnff_update_nbr_bond
  if (lverbose) then 
    time4 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for neighbour list = '',f12.6)') time4 - time3
  endif
!****************************************
!  Compute the log coordination number  *
!****************************************
!
!  Need the first derivatives if ldcharge, even if lgrad1 is not true
!
  call gfnff_cn(numat,gfnff_cnthr,cn,dcn,d2cn,lgrad1.or.ldcharge,.false.)
!*********************************************************
!  Compute the coordination number for hydrogen bonding  *
!*********************************************************
  if (ldohbcn) then
!
!  Set up AHB trio list for cn_hb
!
    call bond_hb_set_AHB(numat,nat,nnbr_bond,maxnbr,nbrno_bond,xbnbr,ybnbr,zbnbr)
!
    call gfnff_cn_hb(cn_hb)
  endif
!
  if (lverbose) then 
    time3 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for coordination number = '',f12.6)') time3 - time4
  endif
!***********************************************************
!  Compute the coordination-dependent electrostatic terms  *
!***********************************************************
  call eem_gfnff(ecoulomb,cn,dcn,d2cn,lgrad1,.false.,.false.)
!
  if (lverbose) then 
    time4 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for charges = '',f12.6)') time4 - time3
  endif
!********************************************************
!  Compute the coordination-dependent dispersion terms  *
!********************************************************
  call gfnff_cnc6(numat,nat,cn,gw,dgwdcn,d2gwdcn2,lgrad1,.false.)
!
  if (lverbose) then 
    time3 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for c6 terms from cn = '',f12.6)') time3 - time4
  endif
!****************
!  Bond energy  *
!****************
  do i = procid+1,numat,nprocs
    nati = nat(i)
    nregioni = nregionno(nrelf2a(i))
    lslicei = lsliceatom(nsft+nrelf2a(i))
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
!
!  Only compute half of bonds to avoid duplication between i->j and j->i
!
      if (j.lt.i) cycle
!
!  Set up quantities for this pair
!
      natj = nat(j)
      nregionj = nregionno(nrelf2a(j))
      lslicej = lsliceatom(nsft+nrelf2a(j))
!
      lreg2one  = .false.
      lreg2pair = .false.
      if (lseok.and.nregions(ncf).gt.1) then
        lreg2pair = (nregioni.eq.nregionj.and.lregionrigid(nregioni,ncf))
        if (.not.lreg2pair) lreg2one = (lregionrigid(nregioni,ncf).neqv.lregionrigid(nregionj,ncf))
      endif
      lattach = (lslicei.and..not.lslicej.or.lslicej.and..not.lslicei)
!
!  Correction for i = j
!
      if (i.eq.j) then
        fct = 0.5_dp
      else
        fct = 1.0_dp
      endif
      r  = rbnbr(ni,i)
      r0 = par_gfnff_bond(1,ni,i)
      p2 = par_gfnff_bond(2,ni,i)
      p3 = par_gfnff_bond(3,ni,i)*fct
!
      if (lgrad1) then
        xji = xbnbr(ni,i)
        yji = ybnbr(ni,i)
        zji = zbnbr(ni,i)
      endif
!
!  Estimate corrected r0 based on coordination numbers and charges
!
      call gfnff_radnoq(nati,natj,cn(i),cn(j),r0,dr0dcni,dr0dcnj,lgrad1)
!
      dr = r - r0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Bond energy for bonded pairs  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (nbrno_hb(ni,i).ge.1) then
!------------------------------------
!  Special case for hydrogen bonds  !
!------------------------------------
        if (nati.eq.1) then
          hbH = i
          hbA = j
        elseif (natj.eq.1) then
          hbH = j
          hbA = i
        else
          call outerror('no H-atom found in bond where one was expected',0_i4)
          call stopnow('gfnff_energy')
        endif
!
!  Energy
!
        tmp = (1.0_dp - (1.0_dp - gfnff_bondscale)*cn_hb(hbH))*p2
        eij = p3*exp(-tmp*dr**2)
!
!  First derivatives
!
        if (lgrad1) then
          deijddr = - 2.0_dp*tmp*dr*eij   ! (dE/dr)
!
!  Coordination number derivatives
!
          deijdcn = - deijddr*dr0dcni
          call gfnff_drv_dcn(i,deijdcn,cn,dcn)
          deijdcn = - deijddr*dr0dcnj
          call gfnff_drv_dcn(j,deijdcn,cn,dcn)
!
          deijddr = deijddr/r   ! (1/r(dE/dr))
!
          xdrv(i) = xdrv(i) - deijddr*xji
          ydrv(i) = ydrv(i) - deijddr*yji
          zdrv(i) = zdrv(i) - deijddr*zji
!
          xdrv(j) = xdrv(j) + deijddr*xji
          ydrv(j) = ydrv(j) + deijddr*yji
          zdrv(j) = zdrv(j) + deijddr*zji
!
          if (lstr) then
            call real1strterm(ndim,xji,yji,zji,xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijddr*dr2ds(ks)
            enddo
          endif
!
!  Get HB coordination number derivatives for hbH
!
          deijdcn = dr*dr*eij*p2*(1.0_dp - gfnff_bondscale)
          call gfnff_drv_dcn_hb(hbH,deijdcn)
        endif
      else
!-----------------
!  General case  !
!-----------------
!
!  Energy
!
        eij = p3*exp(-p2*dr**2)
!
!  First derivatives
!
        if (lgrad1) then
          deijddr = - 2.0_dp*p2*dr*eij   ! (dE/dr)
!
!  Coordination number derivatives
!
          deijdcn = - deijddr*dr0dcni
          call gfnff_drv_dcn(i,deijdcn,cn,dcn)
          deijdcn = - deijddr*dr0dcnj
          call gfnff_drv_dcn(j,deijdcn,cn,dcn)
!
          deijddr = deijddr/r   ! (1/r(dE/dr))
!
          xdrv(i) = xdrv(i) - deijddr*xji
          ydrv(i) = ydrv(i) - deijddr*yji
          zdrv(i) = zdrv(i) - deijddr*zji
!
          xdrv(j) = xdrv(j) + deijddr*xji
          ydrv(j) = ydrv(j) + deijddr*yji
          zdrv(j) = zdrv(j) + deijddr*zji
!
          if (lstr) then
            call real1strterm(ndim,xji,yji,zji,xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijddr*dr2ds(ks)
            enddo
          endif
        endif
      endif
!
!  Handle regions
!
      if (lreg2one) then
        esregion12 = esregion12 + eij
      elseif (lreg2pair) then
        esregion2 = esregion2 + eij
      else
        ebond = ebond + eij
      endif
!
!  Region - region energy
!
      eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eij
!
!  Attachment energy
!
      if (lattach) eattach = eattach + eij
!
!  Site energies
!
      siteenergy(i) = siteenergy(i) + 0.5_dp*eij
      siteenergy(j) = siteenergy(j) + 0.5_dp*eij
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Repulsive energy for bonded pairs  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rep_rho = sqrt(gfnff_repulsion_a(nati)*gfnff_repulsion_a(natj))
      repab = gfnff_repulsion_z(nati)*gfnff_repulsion_z(natj)*gfnff_repscale_b
!
!  Energy
!
      r15 = r**1.5_dp
      tmp = exp(-rep_rho*r15)*repab*fct
      eij = tmp/r
!
!  Handle regions
!
      if (lreg2one) then
        esregion12 = esregion12 + eij
      elseif (lreg2pair) then
        esregion2 = esregion2 + eij
      else
        erepulsion = erepulsion + eij
      endif
!
!  Region - region energy
!
      eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eij
!
!  Attachment energy
!
      if (lattach) eattach = eattach + eij
!
!  Site energies
!
      siteenergy(i) = siteenergy(i) + 0.5_dp*eij
      siteenergy(j) = siteenergy(j) + 0.5_dp*eij
!
!  First derivatives
!
      if (lgrad1) then
        deijdr = - eij/r**2 - 1.5_dp*rep_rho*eij/r**0.5_dp     !  (1/r)(dE/dr)
!
        xdrv(i) = xdrv(i) - xji*deijdr
        ydrv(i) = ydrv(i) - yji*deijdr
        zdrv(i) = zdrv(i) - zji*deijdr
!
        xdrv(j) = xdrv(j) + xji*deijdr
        ydrv(j) = ydrv(j) + yji*deijdr
        zdrv(j) = zdrv(j) + zji*deijdr
!
        if (lstr) then
          call real1strterm(ndim,xji,yji,zji,xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijdr*dr2ds(ks)
          enddo
        endif
      endif
    enddo
  enddo
!
  if (lverbose) then 
    time4 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for bond energy = '',f12.6)') time4 - time3
  endif
!**************************************************
!  Long-range terms for repulsion and dispersion  *
!**************************************************
!
!  Initialise lattice vectors for this cutoff
!
  cutdr = max(gfnff_dispthr,gfnff_repthr)
  call gfnff_setpaircell(cutdr,lr_paircell)
!
  ind = 0
  if (ndim.eq.0) then
    imin = 2
  else
    imin = 1
  endif
  do i = procid+imin,numat,nprocs
    nati = nat(i)
    nregioni = nregionno(nrelf2a(i))
    lslicei = lsliceatom(nsft+nrelf2a(i))
    if (ndim.eq.0) then
      jmax = i - 1
    else
      jmax = i
    endif
!
!  In order to handle H-H pairs in repulsion where they are connected by 2/3 bonds set extended neighbour list
!
    if (nati.eq.1) then
      call gfnff_get_n3atoms(numat,i,3_i4)
    else
      call gfnff_get_n3atoms(numat,i,1_i4)
    endif
!
    do j = 1,jmax
      ind = i*(i-1)/2 + j
      natj = nat(j)
      nregionj = nregionno(nrelf2a(j))
      lslicej = lsliceatom(nsft+nrelf2a(j))
      if (nati.ge.natj) then
        indn = nati*(nati-1)/2 + natj
      else
        indn = natj*(natj-1)/2 + nati
      endif
!
!  Correction for i = j
!
      if (i.eq.j) then
        fct = 0.5_dp
      else
        fct = 1.0_dp
      endif
!
!  Set up region related quantities
!
      lreg2one  = .false.
      lreg2pair = .false.
      if (lseok.and.nregions(ncf).gt.1) then
        lreg2pair = (nregioni.eq.nregionj.and.lregionrigid(nregioni,ncf))
        if (.not.lreg2pair) lreg2one = (lregionrigid(nregioni,ncf).neqv.lregionrigid(nregionj,ncf))
      endif
      lattach = (lslicei.and..not.lslicej.or.lslicej.and..not.lslicei)
!
!  Find interactions for this pair
!
      call gfnff_getpairs(i,j,cutdr,lr_paircell,lr_pairs)
!
!  Set distance-independent constants for dispersion
!
      r4r2ij = 3.0_dp*d4_sqrtZr4r2(nati)*d4_sqrtZr4r2(natj)
      r0 = d4_r0(indn)
!
!  Set distance-independent constants for repulsion
!
      rep_rho = gfnff_repulsion_p(ind)
      repab = gfnff_repulsion_z(nati)*gfnff_repulsion_z(natj)*gfnff_repscale_n
!
!  Set c6 for this pair
!
      zetaprod = d4_zeta_c6(i)*d4_zeta_c6(j)
      call gfnff_setc6(numat,nati,natj,i,j,gw,dgwdcn,d2gwdcn2,c6,dc6dcni,dc6dcnj, &
                       d2c6dcni2,d2c6dcnj2,d2c6dcnidcnj,lgrad1,.false.)
!
      if (lgrad1) then
        deijdcn_sum = 0.0_dp
      endif
!
!  Loop over valid interactions (self term is excluded by getpairs)
!
      do np = 1,lr_pairs%npair
        r2 = lr_pairs%r2pair(np)
        xji = lr_pairs%xpair(np)
        yji = lr_pairs%ypair(np)
        zji = lr_pairs%zpair(np)
!!!!!!!!!!!!!!!!!!!!!!!
!  Dispersion energy  !
!!!!!!!!!!!!!!!!!!!!!!!
        if (r2.lt.gfnff_dispthr) then
!
!  Taper
!
          if (r2.gt.t_dispthr) then
            call mdftaper(r2,t_dispthr,gfnff_dispthr,t,dtdr2,d2tdr22,d3tdr23,lgrad1,.false.,.false.)
            if (lgrad1) dtdr = 2.0_dp*dtdr2
          else
            t = 1.0_dp
            dtdr = 0.0_dp
          endif
!
!  Compute dispersion terms
!
          r03 = r0**3
          r04 = r0**4
          trm6 = 1.0_dp/(r2**3 + r03)
          trm8 = 1.0_dp/(r2**4 + r04)
!
          disp = (trm6 + 2.0_dp*r4r2ij*trm8)*zetaprod
          eij = - fct*c6*disp*t
!
!  Handle regions
!
          if (lreg2one) then
            esregion12 = esregion12 + eij
          elseif (lreg2pair) then
            esregion2 = esregion2 + eij
          else
            edispersion = edispersion + eij
          endif
!
!  Region - region energy
!
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eij
!
!  Attachment energy
!
          if (lattach) eattach = eattach + eij
!
!  Site energies
!
          siteenergy(i) = siteenergy(i) + 0.5_dp*eij
          siteenergy(j) = siteenergy(j) + 0.5_dp*eij
!
          if (lgrad1) then
!
!  Coordination number derivatives of c6
!
            deijdcn_sum = deijdcn_sum - disp*t
!
!  Distance derivatives
!
            ddispdr = - r2*r2*(6.0_dp*trm6*trm6 + 16.0_dp*r4r2ij*r2*trm8*trm8)*zetaprod
            deijdr = - fct*c6*(ddispdr*t + disp*dtdr)   ! (1/r(dE/dr))
!
            xdrv(i) = xdrv(i) - deijdr*xji
            ydrv(i) = ydrv(i) - deijdr*yji
            zdrv(i) = zdrv(i) - deijdr*zji
!
            xdrv(j) = xdrv(j) + deijdr*xji
            ydrv(j) = ydrv(j) + deijdr*yji
            zdrv(j) = zdrv(j) + deijdr*zji
!
            if (lstr) then
              call real1strterm(ndim,xji,yji,zji,xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
              do kl = 1,nstrains
                ks = nstrptr(kl)
                rstrd(kl) = rstrd(kl) + deijdr*dr2ds(ks)
              enddo
            endif
          endif
        endif
!!!!!!!!!!!!!!!!!!!!!!
!  Repulsive energy  !
!!!!!!!!!!!!!!!!!!!!!!
        if (r2.lt.gfnff_repthr) then
!
!  Taper
!
          if (r2.gt.t_repthr) then
            call mdftaper(r2,t_repthr,gfnff_repthr,t,dtdr2,d2tdr22,d3tdr23,lgrad1,.false.,.false.)
            if (lgrad1) dtdr = 2.0_dp*dtdr2
          else
            t = 1.0_dp
            dtdr = 0.0_dp
          endif
!
!  Excluded bonded pairs which are handled separately
!
          call gfnff_bond_check(j,xji,yji,zji,1_i4,lbonded)
          if (lbonded) cycle
!
          ff23 = 1.0_dp
          if (natj.eq.1) then
!
!  Check whether j is 2 or 3 bonds away from i
!
            call gfnff_bond_shell(j,xji,yji,zji,nshell)
            if (nshell.eq.2) then
              ff23 = gfnff_repscale_13
            elseif (nshell.eq.3) then
              ff23 = gfnff_repscale_14
            endif
          endif
          rep_rholoc = rep_rho*ff23
!
!  Energy
!
          r = sqrt(r2)
          r15 = r**1.5_dp
          tmp = exp(-rep_rholoc*r15)*repab
          tmprr = tmp/r
          eij = fct*tmprr*t
!
!  Handle regions
!
          if (lreg2one) then
            esregion12 = esregion12 + eij
          elseif (lreg2pair) then
            esregion2 = esregion2 + eij
          else
            erepulsion = erepulsion + eij
          endif
!
!  Region - region energy
!
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eij
!
!  Attachment energy
!
          if (lattach) eattach = eattach + eij
!
!  Site energies
!
          siteenergy(i) = siteenergy(i) + 0.5_dp*eij
          siteenergy(j) = siteenergy(j) + 0.5_dp*eij
!
!  First derivatives
!
          if (lgrad1) then
            deijdr = - tmprr*t/r**2 - 1.5_dp*rep_rholoc*tmprr*t/r**0.5_dp + tmprr*dtdr  !  (1/r)(dE/dr)
            deijdr = fct*deijdr
!
            xdrv(i) = xdrv(i) - xji*deijdr
            ydrv(i) = ydrv(i) - yji*deijdr
            zdrv(i) = zdrv(i) - zji*deijdr
!
            xdrv(j) = xdrv(j) + xji*deijdr
            ydrv(j) = ydrv(j) + yji*deijdr
            zdrv(j) = zdrv(j) + zji*deijdr
!
            if (lstr) then
              call real1strterm(ndim,xji,yji,zji,xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
              do kl = 1,nstrains
                ks = nstrptr(kl)
                rstrd(kl) = rstrd(kl) + deijdr*dr2ds(ks)
              enddo
            endif
          endif
        endif
      enddo
!
!  Collective derivatives summed over all pairs
!
      if (lgrad1) then
!
!  Coordination number derivatives of c6
!
        deijdcn = fct*dc6dcni*deijdcn_sum
        call gfnff_drv_dcn(i,deijdcn,cn,dcn)
        deijdcn = fct*dc6dcnj*deijdcn_sum
        call gfnff_drv_dcn(j,deijdcn,cn,dcn)
      endif
    enddo
  enddo
!
  if (lverbose) then 
    time3 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for long range energy = '',f12.6)') time3 - time4
  endif
!******************
!  Angle bending  *
!******************
  do ia = procid+1,n_gfnff_angles,nprocs
    i  = n_gfnff_angleptr(1,ia)
    nj = n_gfnff_angleptr(2,ia)
    nk = n_gfnff_angleptr(3,ia)
    j = nbrno_bond(nj,i)
    k = nbrno_bond(nk,i)
!
    nregioni = nregionno(nrelf2a(i))
    nregionj = nregionno(nrelf2a(j))
    nregionk = nregionno(nrelf2a(k))
    lslicei = lsliceatom(nsft+nrelf2a(i))
    lslicej = lsliceatom(nsft+nrelf2a(j))
    lslicek = lsliceatom(nsft+nrelf2a(k))
!
!  Setup region quantities
!
    lreg12    = .false.
    lreg2trio = .false.
    if (lseok.and.nregions(ncf).gt.1) then
      lreg2trio = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1)
      if (.not.lreg2trio) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1)
    endif
    lattach = .true.
    if (lslicei.and.lslicej.and.lslicek) lattach = .false.
    if (.not.lslicei.and..not.lslicej.and..not.lslicek) lattach = .false.
!
!  Set parameters
!
    theta0 = par_gfnff_angle(1,ia)
    rkang  = par_gfnff_angle(2,ia)
!
    vij(1) = xbnbr(nj,i)
    vij(2) = ybnbr(nj,i)
    vij(3) = zbnbr(nj,i)
!
    vik(1) = xbnbr(nk,i)
    vik(2) = ybnbr(nk,i)
    vik(3) = zbnbr(nk,i)
!
    vjk(1) = vik(1) - vij(1)
    vjk(2) = vik(2) - vij(2)
    vjk(3) = vik(3) - vij(3)
!
    rij = rbnbr(nj,i)
    rik = rbnbr(nk,i)
!
    r2ij = rij**2
    r2ik = rik**2
    r2jk = vjk(1)**2 + vjk(2)**2 + vjk(3)**2
!
!  Compute the cosine of the angle and the angle
!
    cosang = 0.5_dp*(r2ij + r2ik - r2jk)/(rij*rik)
    cosang  = dble(min(1.0_dp,max(-1.0_dp,cosang)))
    theta = dacos(cosang)
!
!  Compute cutoff factors for damping
!
    rcut2ij = (gfnff_angle_damp*(gfnff_rcov(nat(i))+gfnff_rcov(nat(j)))**2)**2
    rcut2ik = (gfnff_angle_damp*(gfnff_rcov(nat(i))+gfnff_rcov(nat(k)))**2)**2
!
!  Damping factors
!
    call gfnff_damp(r2ij,rcut2ij,dampij,ddampij,d2dampij,lgrad1,.false.)
    call gfnff_damp(r2ik,rcut2ik,dampik,ddampik,d2dampik,lgrad1,.false.)
    damp = dampij*dampik
!
    linear = (pi-theta0.lt.1.d-6)
    if (linear) then 
!
!  Linear
!
      dtheta = theta - theta0
      eijknodamp = rkang*dtheta**2
      eijk = damp*eijknodamp
! DEBUG - need to improve handling of angles approaching pi
      if ((pi-theta).gt.1.0d-6) then
        deijkdcos = - 2.0_dp*damp*rkang*dtheta/sin(theta)
      else
        deijkdcos = 0.0_dp
      endif
    else
      costheta0 = cos(theta0)
      eijknodamp = rkang*(cosang - costheta0)**2
      eijk = damp*eijknodamp
      deijkdcos = 2.0_dp*damp*rkang*(cosang - costheta0)
    endif
!
!  Handle regions
!
    if (lreg2trio) then
      esregion2 = esregion2 + eijk
    elseif (lreg12) then
      esregion12 = esregion12 + eijk
    else
      eangle = eangle + eijk
    endif
!
!  Region - region energy
!
    eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*eijk
    eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + third*eijk
    eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + third*eijk
!
!  Attachment energy
!
    if (lattach) eattach = eattach + eijk
!
!  Site energies
!
    siteenergy(i) = siteenergy(i) + third*eijk
    siteenergy(j) = siteenergy(j) + third*eijk
    siteenergy(k) = siteenergy(k) + third*eijk
!
!  First derivatives
!
    if (lgrad1) then
      rrij = 1.0_dp/rij
      rrik = 1.0_dp/rik
!
      dcosdrij = rrij*rrik - cosang*rrij**2 
      dcosdrik = rrij*rrik - cosang*rrik**2
      dcosdrjk = - rrij*rrik
!
!  i-j
!
      deijkdr = deijkdcos*dcosdrij      ! (1/r(dEdr))
      deijkdr = deijkdr + eijknodamp*dampik*ddampij   ! Add on damping derivatives for i-j
      xdrv(i) = xdrv(i) - deijkdr*vij(1)
      ydrv(i) = ydrv(i) - deijkdr*vij(2)
      zdrv(i) = zdrv(i) - deijkdr*vij(3)
      xdrv(j) = xdrv(j) + deijkdr*vij(1)
      ydrv(j) = ydrv(j) + deijkdr*vij(2)
      zdrv(j) = zdrv(j) + deijkdr*vij(3)
!
      if (lstr) then
        call real1strterm(ndim,vij(1),vij(2),vij(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkdr*dr2ds(ks)
        enddo
      endif
!
!  i-k
!
      deijkdr = deijkdcos*dcosdrik      ! (1/r(dEdr))
      deijkdr = deijkdr + eijknodamp*dampij*ddampik   ! Add on damping derivatives for i-k
      xdrv(i) = xdrv(i) - deijkdr*vik(1)
      ydrv(i) = ydrv(i) - deijkdr*vik(2)
      zdrv(i) = zdrv(i) - deijkdr*vik(3)
      xdrv(k) = xdrv(k) + deijkdr*vik(1)
      ydrv(k) = ydrv(k) + deijkdr*vik(2)
      zdrv(k) = zdrv(k) + deijkdr*vik(3)
!
      if (lstr) then
        call real1strterm(ndim,vik(1),vik(2),vik(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkdr*dr2ds(ks)
        enddo
      endif
!
!  j-k
!
      deijkdr = deijkdcos*dcosdrjk      ! (1/r(dEdr))
      xdrv(j) = xdrv(j) - deijkdr*vjk(1)
      ydrv(j) = ydrv(j) - deijkdr*vjk(2)
      zdrv(j) = zdrv(j) - deijkdr*vjk(3)
      xdrv(k) = xdrv(k) + deijkdr*vjk(1)
      ydrv(k) = ydrv(k) + deijkdr*vjk(2)
      zdrv(k) = zdrv(k) + deijkdr*vjk(3)
!
      if (lstr) then
        call real1strterm(ndim,vjk(1),vjk(2),vjk(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkdr*dr2ds(ks)
        enddo
      endif
    endif
  enddo
!
  if (lverbose) then 
    time4 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for angle energy = '',f12.6)') time4 - time3
  endif
!*************
!  Torsions  *
!*************
  do ia = procid+1,n_gfnff_torsions,nprocs
    nrot = n_gfnff_torsionptr(5,ia)
!
!  Only do proper torsions in this loop
!
    if (nrot.le.0) cycle
!
    rn = dble(nrot)
!
!  Set parameters
!
    phi0 = par_gfnff_torsion(1,ia)
!
    j  = n_gfnff_torsionptr(1,ia)
    njk = n_gfnff_torsionptr(2,ia)
    nij = n_gfnff_torsionptr(3,ia)
    nkl = n_gfnff_torsionptr(4,ia)
!
    i = nbrno_bond(nij,j)
    k = nbrno_bond(njk,j)
    l = nbrno_bond(nkl,k)
!
    nregioni = nregionno(nrelf2a(i))
    nregionj = nregionno(nrelf2a(j))
    nregionk = nregionno(nrelf2a(k))
    nregionl = nregionno(nrelf2a(l))
    lslicei = lsliceatom(nsft+nrelf2a(i))
    lslicej = lsliceatom(nsft+nrelf2a(j))
    lslicek = lsliceatom(nsft+nrelf2a(k))
    lslicel = lsliceatom(nsft+nrelf2a(l))
!
!  Set region 2 quartet flag
!
    lreg12    = .false.
    lreg2qtet = .false.
    if (lseok.and.nregions(ncf).gt.1) then
      lreg2qtet = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1.and.nregionl.gt.1)
      if (.not.lreg2qtet) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1.or.nregionl.gt.1)
    endif
    lattach = .true.
    if (lslicei.and.lslicej.and.lslicek.and.lslicel) lattach = .false.
    if (.not.lslicei.and..not.lslicej.and..not.lslicek.and..not.lslicel) lattach = .false.
!
    vjl(1) = xbnbr(njk,j) + xbnbr(nkl,k)
    vjl(2) = ybnbr(njk,j) + ybnbr(nkl,k)
    vjl(3) = zbnbr(njk,j) + zbnbr(nkl,k)
!
    vik(1) = - xbnbr(nij,j) + xbnbr(njk,j)
    vik(2) = - ybnbr(nij,j) + ybnbr(njk,j)
    vik(3) = - zbnbr(nij,j) + zbnbr(njk,j)
!
    rjl2 = vjl(1)**2 + vjl(2)**2 + vjl(3)**2
    rjk = rbnbr(njk,j)
    rjk2 = rjk**2
    rik2 = vik(1)**2 + vik(2)**2 + vik(3)**2
!
!  Compute cutoff factors for damping
!
    rcut2ik = (gfnff_torsion_damp*(gfnff_rcov(nat(i))+gfnff_rcov(nat(k)))**2)**2
    rcut2jk = (gfnff_torsion_damp*(gfnff_rcov(nat(j))+gfnff_rcov(nat(k)))**2)**2
    rcut2jl = (gfnff_torsion_damp*(gfnff_rcov(nat(j))+gfnff_rcov(nat(l)))**2)**2
!
!  Damping factors
!
    call gfnff_damp(rik2,rcut2ik,dampik,ddampik,d2dampik,lgrad1,.false.)
    call gfnff_damp(rjk2,rcut2jk,dampjk,ddampjk,d2dampjk,lgrad1,.false.)
    call gfnff_damp(rjl2,rcut2jl,dampjl,ddampjl,d2dampjl,lgrad1,.false.)
    damp = dampik*dampjk*dampjl
!
    vij(1) = - xbnbr(nij,j)
    vij(2) = - ybnbr(nij,j)
    vij(3) = - zbnbr(nij,j)
!
    vil(1) = vij(1) + vjl(1)
    vil(2) = vij(2) + vjl(2)
    vil(3) = vij(3) + vjl(3)
    ril2 = vil(1)**2 + vil(2)**2 + vil(3)**2
!
    rij = rbnbr(nij,j)
    rik = sqrt(rik2)
    ril = sqrt(ril2)
    rjl = sqrt(rjl2)
    rkl = rbnbr(nkl,k)
!
    call gfnff_torsion(1_i4,par_gfnff_torsion(2,ia),rn,phi0,rij,rik,ril,rjk,rjl,rkl,eijklnod,e1d,e2d,lgrad1,.false.)
    eijkl = eijklnod*damp
!
!  Handle regions
!
    if (lreg2qtet) then
      esregion2 = esregion2 + eijkl
    elseif (lreg12) then
      esregion12 = esregion12 + eijkl
    else
      etorsion = etorsion + eijkl
    endif
!
!  Region - region energy
!
    eterm6th = eijkl/6.0_dp
    eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eterm6th
    eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + eterm6th
    eregion2region(nregionl,nregioni) = eregion2region(nregionl,nregioni) + eterm6th
    eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + eterm6th
    eregion2region(nregionl,nregionj) = eregion2region(nregionl,nregionj) + eterm6th
    eregion2region(nregionl,nregionk) = eregion2region(nregionl,nregionk) + eterm6th
!
!  Attachment energy
!
    if (lattach) eattach = eattach + eijkl
!
!  Site energies
!
    siteenergy(i) = siteenergy(i) + 0.25_dp*eijkl
    siteenergy(j) = siteenergy(j) + 0.25_dp*eijkl
    siteenergy(k) = siteenergy(k) + 0.25_dp*eijkl
    siteenergy(l) = siteenergy(l) + 0.25_dp*eijkl
!
!  First derivatives
!
    if (lgrad1) then
!
!  Finish interatomic vectors
!
      vjk(1) = xbnbr(njk,j)
      vjk(2) = ybnbr(njk,j)
      vjk(3) = zbnbr(njk,j)
!
      vkl(1) = xbnbr(nkl,k)
      vkl(2) = ybnbr(nkl,k)
      vkl(3) = zbnbr(nkl,k)
!
!  Scale current derivatives by damping
!
      e1d(1:6) = e1d(1:6)*damp
!
!  Add damping derivatives
!
      e1d(2) = e1d(2) + eijklnod*dampjl*dampjk*ddampik
      e1d(4) = e1d(4) + eijklnod*dampjl*ddampjk*dampik
      e1d(5) = e1d(5) + eijklnod*ddampjl*dampjk*dampik
!
      xdrv(i) = xdrv(i) - vij(1)*e1d(1) - vik(1)*e1d(2) - vil(1)*e1d(3)
      ydrv(i) = ydrv(i) - vij(2)*e1d(1) - vik(2)*e1d(2) - vil(2)*e1d(3)
      zdrv(i) = zdrv(i) - vij(3)*e1d(1) - vik(3)*e1d(2) - vil(3)*e1d(3)
!
      xdrv(j) = xdrv(j) - vjk(1)*e1d(4) + vij(1)*e1d(1) - vjl(1)*e1d(5)
      ydrv(j) = ydrv(j) - vjk(2)*e1d(4) + vij(2)*e1d(1) - vjl(2)*e1d(5)
      zdrv(j) = zdrv(j) - vjk(3)*e1d(4) + vij(3)*e1d(1) - vjl(3)*e1d(5)
!
      xdrv(k) = xdrv(k) + vjk(1)*e1d(4) - vkl(1)*e1d(6) + vik(1)*e1d(2)
      ydrv(k) = ydrv(k) + vjk(2)*e1d(4) - vkl(2)*e1d(6) + vik(2)*e1d(2)
      zdrv(k) = zdrv(k) + vjk(3)*e1d(4) - vkl(3)*e1d(6) + vik(3)*e1d(2)
!
      xdrv(l) = xdrv(l) + vkl(1)*e1d(6) + vjl(1)*e1d(5) + vil(1)*e1d(3)
      ydrv(l) = ydrv(l) + vkl(2)*e1d(6) + vjl(2)*e1d(5) + vil(2)*e1d(3)
      zdrv(l) = zdrv(l) + vkl(3)*e1d(6) + vjl(3)*e1d(5) + vil(3)*e1d(3)
!
      if (lstr) then
!
!  Set up strain products
!
        xv4(1) = vij(1)
        yv4(1) = vij(2)
        zv4(1) = vij(3)
        xv4(2) = vik(1)
        yv4(2) = vik(2)
        zv4(2) = vik(3)
        xv4(3) = vil(1)
        yv4(3) = vil(2)
        zv4(3) = vil(3)
        xv4(4) = vjk(1)
        yv4(4) = vjk(2)
        zv4(4) = vjk(3)
        xv4(5) = vjl(1)
        yv4(5) = vjl(2)
        zv4(5) = vjl(3)
        xv4(6) = vkl(1)
        yv4(6) = vkl(2)
        zv4(6) = vkl(3)
!
        call realstrterms(ndim,6_i4,6_i4,xv4,yv4,zv4,xcom,ycom,zcom,dr2ds6,d2r2dx26,d2r2dsdx6,d2r2ds26,.false.)
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + e1d(1)*dr2ds6(ks,1)
          rstrd(kl) = rstrd(kl) + e1d(2)*dr2ds6(ks,2)
          rstrd(kl) = rstrd(kl) + e1d(3)*dr2ds6(ks,3)
          rstrd(kl) = rstrd(kl) + e1d(4)*dr2ds6(ks,4)
          rstrd(kl) = rstrd(kl) + e1d(5)*dr2ds6(ks,5)
          rstrd(kl) = rstrd(kl) + e1d(6)*dr2ds6(ks,6)
        enddo
      endif
    endif
  enddo
!**********************
!  Improper torsions  *
!**********************
  do ia = procid+1,n_gfnff_torsions,nprocs
    nrot = n_gfnff_torsionptr(5,ia)
!
!  Only do improper torsions and out of planes
!
    if (nrot.gt.0) cycle
!
    rn = dble(nrot)
!
!  Set parameters
!
    phi0 = par_gfnff_torsion(1,ia)
!
    i  = n_gfnff_torsionptr(1,ia)
    nij = n_gfnff_torsionptr(2,ia)
    nik = n_gfnff_torsionptr(3,ia)
    nil = n_gfnff_torsionptr(4,ia)
!
    j = nbrno_bond(nij,i)
    k = nbrno_bond(nik,i)
    l = nbrno_bond(nil,i)
!
    nregioni = nregionno(nrelf2a(i))
    nregionj = nregionno(nrelf2a(j))
    nregionk = nregionno(nrelf2a(k))
    nregionl = nregionno(nrelf2a(l))
    lslicei = lsliceatom(nsft+nrelf2a(i))
    lslicej = lsliceatom(nsft+nrelf2a(j))
    lslicek = lsliceatom(nsft+nrelf2a(k))
    lslicel = lsliceatom(nsft+nrelf2a(l))
!
!  Set region 2 quartet flag
!
    lreg12    = .false.   
    lreg2qtet = .false.
    if (lseok.and.nregions(ncf).gt.1) then
      lreg2qtet = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1.and.nregionl.gt.1)
      if (.not.lreg2qtet) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1.or.nregionl.gt.1)
    endif
    lattach = .true.
    if (lslicei.and.lslicej.and.lslicek.and.lslicel) lattach = .false.
    if (.not.lslicei.and..not.lslicej.and..not.lslicek.and..not.lslicel) lattach = .false.
!
    vjk(1) = xbnbr(nik,i) - xbnbr(nij,i)
    vjk(2) = ybnbr(nik,i) - ybnbr(nij,i)
    vjk(3) = zbnbr(nik,i) - zbnbr(nij,i)
!
    vjl(1) = xbnbr(nil,i) - xbnbr(nij,i)
    vjl(2) = ybnbr(nil,i) - ybnbr(nij,i)
    vjl(3) = zbnbr(nil,i) - zbnbr(nij,i)
!
    rij = rbnbr(nij,i)
    rij2 = rij**2
    rjk2 = vjk(1)**2 + vjk(2)**2 + vjk(3)**2
    rjl2 = vjl(1)**2 + vjl(2)**2 + vjl(3)**2
!
!  Compute cutoff factors for damping
!
    rcut2ij = (gfnff_torsion_damp*(gfnff_rcov(nat(i))+gfnff_rcov(nat(j)))**2)**2
    rcut2jk = (gfnff_torsion_damp*(gfnff_rcov(nat(j))+gfnff_rcov(nat(k)))**2)**2
    rcut2jl = (gfnff_torsion_damp*(gfnff_rcov(nat(j))+gfnff_rcov(nat(l)))**2)**2
!
    call gfnff_damp(rij2,rcut2ij,dampij,ddampij,d2dampij,lgrad1,.false.)
    call gfnff_damp(rjk2,rcut2jk,dampjk,ddampjk,d2dampjk,lgrad1,.false.)
    call gfnff_damp(rjl2,rcut2jl,dampjl,ddampjl,d2dampjl,lgrad1,.false.)
    damp = dampij*dampjk*dampjl
!
    call gfnff_inversion(ndim,numat,nnbr_bond,maxnbr,rbnbr,xbnbr,ybnbr,zbnbr,i,nij,nik,nil,phi,phi1dx, &
                         phi1ds,phi2dx,phi2ds,phi2dsdx,lgrad1,.false.)
!
    if (nrot.eq.0) then
      dphi = phi - phi0
      fphi = dphi + pi
      cosphi = cos(fphi)
      eijklnod = (1.0_dp + cosphi)*par_gfnff_torsion(2,ia)
      if (lgrad1) then
        sinphi = sin(fphi)
        deijkldphi = - par_gfnff_torsion(2,ia)*damp*sinphi
      endif
    else
      cosphi  = cos(phi)
      cosphi0 = cos(phi0)
      eijklnod = par_gfnff_torsion(2,ia)*(cosphi - cosphi0)**2
      if (lgrad1) then
        sinphi  = sin(phi)
        deijkldphi = - 2.0_dp*par_gfnff_torsion(2,ia)*damp*sinphi*(cosphi - cosphi0)
      endif
    endif
!
    eijkl = eijklnod*damp
!
!  Handle regions
!
    if (lreg2qtet) then
      esregion2 = esregion2 + eijkl 
    elseif (lreg12) then
      esregion12 = esregion12 + eijkl
    else
      etorsion = etorsion + eijkl
    endif 
!
!  Region - region energy
!
    eterm6th = eijkl/6.0_dp
    eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eterm6th
    eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + eterm6th
    eregion2region(nregionl,nregioni) = eregion2region(nregionl,nregioni) + eterm6th
    eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + eterm6th
    eregion2region(nregionl,nregionj) = eregion2region(nregionl,nregionj) + eterm6th
    eregion2region(nregionl,nregionk) = eregion2region(nregionl,nregionk) + eterm6th
!
!  Attachment energy
!
    if (lattach) eattach = eattach + eijkl
!
!  Site energies
!
    siteenergy(i) = siteenergy(i) + 0.25_dp*eijkl
    siteenergy(j) = siteenergy(j) + 0.25_dp*eijkl
    siteenergy(k) = siteenergy(k) + 0.25_dp*eijkl
    siteenergy(l) = siteenergy(l) + 0.25_dp*eijkl
!
!  First derivatives
!
    if (lgrad1) then
!
!  Finish interatomic vectors
!
      vij(1) = xbnbr(nij,i)
      vij(2) = ybnbr(nij,i)
      vij(3) = zbnbr(nij,i)
!
!  Derivatives of phi
!
      xdrv(i) = xdrv(i) - deijkldphi*phi1dx(1,1)
      ydrv(i) = ydrv(i) - deijkldphi*phi1dx(2,1)
      zdrv(i) = zdrv(i) - deijkldphi*phi1dx(3,1)
!
      xdrv(j) = xdrv(j) + deijkldphi*phi1dx(1,1)
      ydrv(j) = ydrv(j) + deijkldphi*phi1dx(2,1)
      zdrv(j) = zdrv(j) + deijkldphi*phi1dx(3,1)
!
      xdrv(j) = xdrv(j) - deijkldphi*phi1dx(1,2)
      ydrv(j) = ydrv(j) - deijkldphi*phi1dx(2,2)
      zdrv(j) = zdrv(j) - deijkldphi*phi1dx(3,2)
!
      xdrv(k) = xdrv(k) + deijkldphi*phi1dx(1,2)
      ydrv(k) = ydrv(k) + deijkldphi*phi1dx(2,2)
      zdrv(k) = zdrv(k) + deijkldphi*phi1dx(3,2)
!
      xdrv(i) = xdrv(i) - deijkldphi*phi1dx(1,3)
      ydrv(i) = ydrv(i) - deijkldphi*phi1dx(2,3)
      zdrv(i) = zdrv(i) - deijkldphi*phi1dx(3,3)
!
      xdrv(l) = xdrv(l) + deijkldphi*phi1dx(1,3)
      ydrv(l) = ydrv(l) + deijkldphi*phi1dx(2,3)
      zdrv(l) = zdrv(l) + deijkldphi*phi1dx(3,3)
!
      if (lstr) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkldphi*phi1ds(ks)
        enddo
      endif
!
!  Add damping derivatives
!
!  i-j
!
      deijkldrij = eijklnod*ddampij*dampjk*dampjl
      xdrv(i) = xdrv(i) - deijkldrij*vij(1)
      ydrv(i) = ydrv(i) - deijkldrij*vij(2)
      zdrv(i) = zdrv(i) - deijkldrij*vij(3)
!
      xdrv(j) = xdrv(j) + deijkldrij*vij(1)
      ydrv(j) = ydrv(j) + deijkldrij*vij(2)
      zdrv(j) = zdrv(j) + deijkldrij*vij(3)
!
      if (lstr) then
        call real1strterm(ndim,vij(1),vij(2),vij(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkldrij*dr2ds(ks)
        enddo
      endif
!
!  j-k
!
      deijkldrjk = eijklnod*dampij*ddampjk*dampjl
      xdrv(j) = xdrv(j) - deijkldrjk*vjk(1)
      ydrv(j) = ydrv(j) - deijkldrjk*vjk(2)
      zdrv(j) = zdrv(j) - deijkldrjk*vjk(3)
!
      xdrv(k) = xdrv(k) + deijkldrjk*vjk(1)
      ydrv(k) = ydrv(k) + deijkldrjk*vjk(2)
      zdrv(k) = zdrv(k) + deijkldrjk*vjk(3)
!
      if (lstr) then
        call real1strterm(ndim,vjk(1),vjk(2),vjk(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkldrjk*dr2ds(ks)
        enddo
      endif
!
!  j-l
!
      deijkldrjl = eijklnod*dampij*dampjk*ddampjl
      xdrv(j) = xdrv(j) - deijkldrjl*vjl(1)
      ydrv(j) = ydrv(j) - deijkldrjl*vjl(2)
      zdrv(j) = zdrv(j) - deijkldrjl*vjl(3)
!
      xdrv(l) = xdrv(l) + deijkldrjl*vjl(1)
      ydrv(l) = ydrv(l) + deijkldrjl*vjl(2)
      zdrv(l) = zdrv(l) + deijkldrjl*vjl(3)
!
      if (lstr) then
        call real1strterm(ndim,vjl(1),vjl(2),vjl(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkldrjl*dr2ds(ks)
        enddo
      endif
    endif
  enddo
!
  if (lverbose) then 
    time3 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for torsion energy = '',f12.6)') time3 - time4
  endif
!***********************************
!  Bonded atoms 3 body dispersion  *
!***********************************
!
!  Loop over 1-4 (3 bond) atom pairs
!
  do ia = procid+1,nb3atm,nprocs
    i = nb3list(1,ia)
    ni = nb3list(2,ia)
    nj = nb3list(3,ia)
    nk = nb3list(4,ia)
!
    j = nbrno_bond(ni,i)
    k = nbrno_bond(nj,j)
    l = nbrno_bond(nk,k)
!
    nati = nat(i)
    natl = nat(l)
!
    nregioni = nregionno(nrelf2a(i))
    nregionl = nregionno(nrelf2a(l))
    lslicei = lsliceatom(nsft+nrelf2a(i))
    lslicel = lsliceatom(nsft+nrelf2a(l))
!
!  Construct vector and distance from i to l
!
    vil(1) = xbnbr(ni,i) + xbnbr(nj,j) + xbnbr(nk,k)
    vil(2) = ybnbr(ni,i) + ybnbr(nj,j) + ybnbr(nk,k)
    vil(3) = zbnbr(ni,i) + zbnbr(nj,j) + zbnbr(nk,k)
    ril2 = vil(1)**2 + vil(2)**2 + vil(3)**2
    ril = sqrt(ril2)
    rril = 1.0_dp/ril
!
    fil = d4_c9(i)*d4_c9(l)
!
!  Damping terms
!
    if (nati.ge.natl) then
      indil = nati*(nati-1)/2 + natl
    else
      indil = natl*(natl-1)/2 + nati
    endif
    ril0 = atm_alpha1*d4_r0(indil)
    ril03 = ril0**1.5_dp
!------------------------------
!  Loop over neighbours of i  !
!------------------------------
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
      natj = nat(j)
!
!  Compute distances for the triad
!
      rij = rbnbr(ni,i)
      rij2 = rij**2
      vjl(1) = vil(1) - xbnbr(ni,i)
      vjl(2) = vil(2) - ybnbr(ni,i)
      vjl(3) = vil(3) - zbnbr(ni,i)
      rjl2 = vjl(1)**2 + vjl(2)**2 + vjl(3)**2
!
!  Check that j and l are different atoms
!
      if (rjl2.lt.1.0d-2) cycle
!
!  Setup region quantities
!
      nregionj = nregionno(nrelf2a(j))
      lslicej = lsliceatom(nsft+nrelf2a(j))
      lreg12    = .false.
      lreg2trio = .false.
      if (lseok.and.nregions(ncf).gt.1) then
        lreg2trio = (nregioni.gt.1.and.nregionj.gt.1.and.nregionl.gt.1)
        if (.not.lreg2trio) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionl.gt.1)
      endif
      lattach = .true.
      if (lslicei.and.lslicej.and.lslicel) lattach = .false.
      if (.not.lslicei.and..not.lslicej.and..not.lslicel) lattach = .false.
!
!  Form C9 coefficient
!
      c9 = fil*d4_c9(j)
!
!  Damping terms
!
      if (nati.ge.natj) then
        indij = nati*(nati-1)/2 + natj
      else
        indij = natj*(natj-1)/2 + nati
      endif
      if (natj.ge.natl) then
        indjl = natj*(natj-1)/2 + natl
      else
        indjl = natl*(natl-1)/2 + natj
      endif
      rij0 = atm_alpha1*d4_r0(indij)
      rjl0 = atm_alpha1*d4_r0(indjl)
      rij03 = rij0**1.5_dp
      rjl03 = rjl0**1.5_dp
!
      rjl = sqrt(rjl2)
!
      rrij = 1.0_dp/rij
      rrjl = 1.0_dp/rjl
!
      cos1 = 0.5_dp*(rjl2 + rij2 - ril2)*rrij*rrjl
      cos2 = 0.5_dp*(ril2 + rij2 - rjl2)*rrij*rril
      cos3 = 0.5_dp*(ril2 + rjl2 - rij2)*rril*rrjl
!
      rrij3 = 1.0_dp/(rij2*rij + rij03)
      rril3 = 1.0_dp/(ril2*ril + ril03)
      rrjl3 = 1.0_dp/(rjl2*rjl + rjl03)
      rrijk9 = rrij3*rril3*rrjl3
!
      cosfct = (3.0_dp*cos1*cos2*cos3 + 1.0_dp)*rrijk9
      eijk = c9*cosfct ! energy
!
!  Handle regions
!
      if (lreg2trio) then
        esregion2 = esregion2 + eijk
      elseif (lreg12) then
        esregion12 = esregion12 + eijk
      else
        ebatm = ebatm + eijk
      endif
!
!  Region - region energy
!
      eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*eijk
      eregion2region(nregionl,nregioni) = eregion2region(nregionl,nregioni) + third*eijk
      eregion2region(nregionl,nregionj) = eregion2region(nregionl,nregionj) + third*eijk
!
!  Attachment energy
!
      if (lattach) eattach = eattach + eijk
!
!  Site energies
!
      siteenergy(i) = siteenergy(i) + third*eijk
      siteenergy(j) = siteenergy(j) + third*eijk
      siteenergy(l) = siteenergy(l) + third*eijk
!
!  First derivatives
!
      if (lgrad1) then
        dcos1drij = rrij*rrjl - cos1*rrij**2
        dcos1drjl = rrij*rrjl - cos1*rrjl**2
        dcos1dril = - rrij*rrjl
!
        dcos2drij = rrij*rril - cos2*rrij**2
        dcos2dril = rrij*rril - cos2*rril**2
        dcos2drjl = - rrij*rril
!
        dcos3dril = rril*rrjl - cos3*rril**2
        dcos3drjl = rril*rrjl - cos3*rrjl**2
        dcos3drij = - rril*rrjl
!
        deijkdrij = - 3.0_dp*eijk*rij*rrij3 + (3.0_dp*c9*rrijk9)*(dcos1drij*cos2*cos3 + &
                                                                  dcos2drij*cos1*cos3 + &
                                                                  dcos3drij*cos1*cos2)
        deijkdril = - 3.0_dp*eijk*ril*rril3 + (3.0_dp*c9*rrijk9)*(dcos1dril*cos2*cos3 + &
                                                                  dcos2dril*cos1*cos3 + &
                                                                  dcos3dril*cos1*cos2)
        deijkdrjl = - 3.0_dp*eijk*rjl*rrjl3 + (3.0_dp*c9*rrijk9)*(dcos1drjl*cos2*cos3 + &
                                                                  dcos2drjl*cos1*cos3 + &
                                                                  dcos3drjl*cos1*cos2)
!
!  i-j
!
        xdrv(i) = xdrv(i) - deijkdrij*xbnbr(ni,i)
        ydrv(i) = ydrv(i) - deijkdrij*ybnbr(ni,i)
        zdrv(i) = zdrv(i) - deijkdrij*zbnbr(ni,i)
!
        xdrv(j) = xdrv(j) + deijkdrij*xbnbr(ni,i)
        ydrv(j) = ydrv(j) + deijkdrij*ybnbr(ni,i)
        zdrv(j) = zdrv(j) + deijkdrij*zbnbr(ni,i)
!
        if (lstr) then
          call real1strterm(ndim,xbnbr(ni,i),ybnbr(ni,i),zbnbr(ni,i),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds(ks)
          enddo
        endif
!
!  i-l
!
        xdrv(i) = xdrv(i) - deijkdril*vil(1)
        ydrv(i) = ydrv(i) - deijkdril*vil(2)
        zdrv(i) = zdrv(i) - deijkdril*vil(3)
!
        xdrv(l) = xdrv(l) + deijkdril*vil(1)
        ydrv(l) = ydrv(l) + deijkdril*vil(2)
        zdrv(l) = zdrv(l) + deijkdril*vil(3)
!     
        if (lstr) then
          call real1strterm(ndim,vil(1),vil(2),vil(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdril*dr2ds(ks)
          enddo
        endif
!
!  j-l
!
        xdrv(j) = xdrv(j) - deijkdrjl*vjl(1)
        ydrv(j) = ydrv(j) - deijkdrjl*vjl(2)
        zdrv(j) = zdrv(j) - deijkdrjl*vjl(3)
!
        xdrv(l) = xdrv(l) + deijkdrjl*vjl(1)
        ydrv(l) = ydrv(l) + deijkdrjl*vjl(2)
        zdrv(l) = zdrv(l) + deijkdrjl*vjl(3)
!     
        if (lstr) then
          call real1strterm(ndim,vjl(1),vjl(2),vjl(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdrjl*dr2ds(ks)
          enddo
        endif
      endif
    enddo
!------------------------------
!  Loop over neighbours of l  !
!------------------------------
    do nl = 1,nnbr_bond(l)
      j = nbrno_bond(nl,l)
      natj = nat(j)
!
!  Compute distances for the triad
!
      rjl = rbnbr(nl,l)
      rjl2 = rjl**2
      vij(1) = vil(1) + xbnbr(nl,l)
      vij(2) = vil(2) + ybnbr(nl,l)
      vij(3) = vil(3) + zbnbr(nl,l)
      rij2 = vij(1)**2 + vij(2)**2 + vij(3)**2
!
!  Check that i and j are different atoms
!
      if (rij2.lt.1.0d-2) cycle
!
!  Setup region quantities
!
      nregionj = nregionno(nrelf2a(j))
      lslicej = lsliceatom(nsft+nrelf2a(j))
      lreg12    = .false.
      lreg2trio = .false.
      if (lseok.and.nregions(ncf).gt.1) then
        lreg2trio = (nregioni.gt.1.and.nregionj.gt.1.and.nregionl.gt.1)
        if (.not.lreg2trio) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionl.gt.1)
      endif
      lattach = .true.
      if (lslicei.and.lslicej.and.lslicel) lattach = .false.
      if (.not.lslicei.and..not.lslicej.and..not.lslicel) lattach = .false.
!
!  Form C9 coefficient
!
      c9 = fil*d4_c9(j)
!
!  Damping terms
!
      if (nati.ge.natj) then
        indij = nati*(nati-1)/2 + natj
      else
        indij = natj*(natj-1)/2 + nati
      endif
      if (natj.ge.natl) then
        indjl = natj*(natj-1)/2 + natl
      else
        indjl = natl*(natl-1)/2 + natj
      endif
      rij0 = atm_alpha1*d4_r0(indij)
      rjl0 = atm_alpha1*d4_r0(indjl)
      rij03 = rij0**1.5_dp
      rjl03 = rjl0**1.5_dp
!
      rij = sqrt(rij2)
!
      rrij = 1.0_dp/rij
      rrjl = 1.0_dp/rjl
!
      cos1 = 0.5_dp*(rjl2 + rij2 - ril2)*rrij*rrjl
      cos2 = 0.5_dp*(ril2 + rij2 - rjl2)*rrij*rril
      cos3 = 0.5_dp*(ril2 + rjl2 - rij2)*rril*rrjl
!
      rrij3 = 1.0_dp/(rij2*rij + rij03)
      rril3 = 1.0_dp/(ril2*ril + ril03)
      rrjl3 = 1.0_dp/(rjl2*rjl + rjl03)
      rrijk9 = rrij3*rril3*rrjl3
!
      cosfct = (3.0_dp*cos1*cos2*cos3 + 1.0_dp)*rrijk9
      eijk = c9*cosfct ! energy
!
!  Handle regions
!
      if (lreg2trio) then
        esregion2 = esregion2 + eijk
      elseif (lreg12) then
        esregion12 = esregion12 + eijk
      else
        ebatm = ebatm + eijk
      endif
!
!  Region - region energy
!
      eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*eijk
      eregion2region(nregionl,nregioni) = eregion2region(nregionl,nregioni) + third*eijk
      eregion2region(nregionl,nregionj) = eregion2region(nregionl,nregionj) + third*eijk
!
!  Attachment energy
!
      if (lattach) eattach = eattach + eijk
!
!  Site energies
!
      siteenergy(i) = siteenergy(i) + third*eijk
      siteenergy(j) = siteenergy(j) + third*eijk
      siteenergy(l) = siteenergy(l) + third*eijk
!
!  First derivatives
!
      if (lgrad1) then
        dcos1drij = rrij*rrjl - cos1*rrij**2
        dcos1drjl = rrij*rrjl - cos1*rrjl**2
        dcos1dril = - rrij*rrjl
!
        dcos2drij = rrij*rril - cos2*rrij**2
        dcos2dril = rrij*rril - cos2*rril**2
        dcos2drjl = - rrij*rril
!
        dcos3dril = rril*rrjl - cos3*rril**2
        dcos3drjl = rril*rrjl - cos3*rrjl**2
        dcos3drij = - rril*rrjl
!
        deijkdrij = - 3.0_dp*eijk*rij*rrij3 + (3.0_dp*c9*rrijk9)*(dcos1drij*cos2*cos3 + &
                                                                  dcos2drij*cos1*cos3 + &
                                                                  dcos3drij*cos1*cos2)
        deijkdril = - 3.0_dp*eijk*ril*rril3 + (3.0_dp*c9*rrijk9)*(dcos1dril*cos2*cos3 + &
                                                                  dcos2dril*cos1*cos3 + &
                                                                  dcos3dril*cos1*cos2)
        deijkdrjl = - 3.0_dp*eijk*rjl*rrjl3 + (3.0_dp*c9*rrijk9)*(dcos1drjl*cos2*cos3 + &
                                                                  dcos2drjl*cos1*cos3 + &
                                                                  dcos3drjl*cos1*cos2)
!
!  i-j
!
        xdrv(i) = xdrv(i) - deijkdrij*vij(1)
        ydrv(i) = ydrv(i) - deijkdrij*vij(2)
        zdrv(i) = zdrv(i) - deijkdrij*vij(3)
!
        xdrv(j) = xdrv(j) + deijkdrij*vij(1)
        ydrv(j) = ydrv(j) + deijkdrij*vij(2)
        zdrv(j) = zdrv(j) + deijkdrij*vij(3)
!
        if (lstr) then
          call real1strterm(ndim,vij(1),vij(2),vij(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds(ks)
          enddo
        endif
!
!  i-l
!
        xdrv(i) = xdrv(i) - deijkdril*vil(1)
        ydrv(i) = ydrv(i) - deijkdril*vil(2)
        zdrv(i) = zdrv(i) - deijkdril*vil(3)
!
        xdrv(l) = xdrv(l) + deijkdril*vil(1)
        ydrv(l) = ydrv(l) + deijkdril*vil(2)
        zdrv(l) = zdrv(l) + deijkdril*vil(3)
!
        if (lstr) then
          call real1strterm(ndim,vil(1),vil(2),vil(3),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdril*dr2ds(ks)
          enddo
        endif
!
!  j-l
!
        xdrv(l) = xdrv(l) - deijkdrjl*xbnbr(nl,l)
        ydrv(l) = ydrv(l) - deijkdrjl*ybnbr(nl,l)
        zdrv(l) = zdrv(l) - deijkdrjl*zbnbr(nl,l)
!
        xdrv(j) = xdrv(j) + deijkdrjl*xbnbr(nl,l)
        ydrv(j) = ydrv(j) + deijkdrjl*ybnbr(nl,l)
        zdrv(j) = zdrv(j) + deijkdrjl*zbnbr(nl,l)
!
        if (lstr) then
          call real1strterm(ndim,xbnbr(nl,l),ybnbr(nl,l),zbnbr(nl,l),xcom(1),ycom(1),zcom(1),dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdrjl*dr2ds(ks)
          enddo
        endif
      endif
    enddo
  enddo
!
  if (lverbose) then 
    time4 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for bond 3 dispersion energy = '',f12.6)') time4 - time3
  endif
!*********************************
!  Hydrogen and halogen bonding  *
!*********************************
  call gfnff_ehb(ehb,exb,nnbr_bond,maxnbr,nbrno_bond,rbnbr,xbnbr,ybnbr,zbnbr,lgrad1,.false.)
!
  if (lverbose) then 
    time3 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for hydrogen bond energy = '',f12.6)') time3 - time4
  endif
!
!  Add up energy components
!
  egfnff = ebond + eangle + etorsion + erepulsion + ecoulomb + edispersion + ehb + exb + ebatm
!
!  Add components
!
  if (nprocs.gt.1.and.index(keyword,'verb').ne.0) then
    sum0(1) = ebond
    sum0(2) = eangle
    sum0(3) = etorsion
    sum0(4) = erepulsion
    sum0(5) = ecoulomb
    sum0(6) = edispersion
    sum0(7) = ehb
    sum0(8) = exb
    sum0(9) = ebatm
!
    call sumall(sum0,sum1,9_i4,"gfnffmd","energies")
!
    ebond = sum1(1)
    eangle = sum1(2)
    etorsion = sum1(3)
    erepulsion = sum1(4)
    ecoulomb = sum1(5)
    edispersion = sum1(6)
    ehb = sum1(7)
    exb = sum1(8)
    ebatm = sum1(9)
    egfnff_tot = ebond + eangle + etorsion + erepulsion + ecoulomb + edispersion + ehb + exb + ebatm
  else
    egfnff_tot = ebond + eangle + etorsion + erepulsion + ecoulomb + edispersion + ehb + exb + ebatm
  endif
!
  if (lverbose) then
    write(ioout,'(/,''  GFNFF energy components: '',/)')
    write(ioout,'(''  Bonding        '',1x,f14.6,'' eV '',1x,f16.8,'' au'')') ebond,ebond/gfnff_autoev
    write(ioout,'(''  Angle bending  '',1x,f14.6,'' eV '',1x,f16.8,'' au'')') eangle,eangle/gfnff_autoev
    write(ioout,'(''  Torsions       '',1x,f14.6,'' eV '',1x,f16.8,'' au'')') etorsion,etorsion/gfnff_autoev
    write(ioout,'(''  Repulsion      '',1x,f14.6,'' eV '',1x,f16.8,'' au'')') erepulsion,erepulsion/gfnff_autoev
    write(ioout,'(''  Electrostatics '',1x,f14.6,'' eV '',1x,f16.8,'' au'')') ecoulomb,ecoulomb/gfnff_autoev
    write(ioout,'(''  Dispersion     '',1x,f14.6,'' eV '',1x,f16.8,'' au'')') edispersion,edispersion/gfnff_autoev
    write(ioout,'(''  Hydrogen bonds '',1x,f14.6,'' eV '',1x,f16.8,'' au'')') ehb,ehb/gfnff_autoev
    write(ioout,'(''  XB energy      '',1x,f14.6,'' eV '',1x,f16.8,'' au'')') exb,exb/gfnff_autoev
    write(ioout,'(''  Bonded atom    '',1x,f14.6,'' eV '',1x,f16.8,'' au'')') ebatm,ebatm/gfnff_autoev
    write(ioout,'(''  TOTAL          '',1x,f14.6,'' eV '',1x,f16.8,'' au'',/)') egfnff_tot,egfnff_tot/gfnff_autoev
  endif
!
!  Free local memory
!
  if (ldohbcn) then
    deallocate(cn_hb,stat=status)
    if (status/=0) call deallocate_error('gfnff_energy','cn_hb')
  endif
!
  deallocate(d2gwdcn2,stat=status)
  if (status/=0) call deallocate_error('gfnff_energy','d2gwdcn2')
  deallocate(dgwdcn,stat=status)
  if (status/=0) call deallocate_error('gfnff_energy','dgwdcn')
  deallocate(gw,stat=status)
  if (status/=0) call deallocate_error('gfnff_energy','gw')
  deallocate(d2cn,stat=status)
  if (status/=0) call deallocate_error('gfnff_energy','d2cn')
  deallocate(dcn,stat=status)
  if (status/=0) call deallocate_error('gfnff_energy','dcn')
  deallocate(cn,stat=status)
  if (status/=0) call deallocate_error('gfnff_energy','cn')
!
  time2 = g_cpu_time()
  tgfnff = tgfnff + time2 - time1
#ifdef TRACE
  call trace_out('egfnff_energy')
#endif
!
  return
  end
