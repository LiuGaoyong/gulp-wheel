  subroutine gfnff_energyd_s(egfnff,lgrad1,lgrad2)
!
!  Calculates the energy and derivatives for the GFNFF force field.
!  Distributed memory version. Algorithm improves scaling for dispersion term.
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!
!  On exit :
!
!  egfnff          = the value of the energy contribution
!
!  NB: For second derivatives, add to d2(j,i) j > i and symmetrise
!      at the end.
!
!   5/22 Created from gfnff_energyd
!   7/23 nregionnocfg added
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
!  Julian Gale, CIC, Curtin University, September  2023
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
  use parallel
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
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: hbA
  integer(i4)                                      :: hbH
  integer(i4)                                      :: i
  integer(i4)                                      :: ia
  integer(i4)                                      :: ii
  integer(i4)                                      :: iloc
  integer(i4)                                      :: imin
  integer(i4)                                      :: ind
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: indn
  integer(i4)                                      :: indij
  integer(i4)                                      :: indil
  integer(i4)                                      :: indjl
  integer(i4)                                      :: indx(2,6)
  integer(i4)                                      :: indy(2,6)
  integer(i4)                                      :: indz(2,6)
  integer(i4)                                      :: indxf(2,6)
  integer(i4)                                      :: indyf(2,6)
  integer(i4)                                      :: indzf(2,6)
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: ixf
  integer(i4)                                      :: iyf
  integer(i4)                                      :: izf
  integer(i4)                                      :: j
  integer(i4)                                      :: jj
  integer(i4)                                      :: jloc
  integer(i4)                                      :: jmax
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: jxf
  integer(i4)                                      :: jyf
  integer(i4)                                      :: jzf
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: l
  integer(i4)                                      :: lloc
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: lxf
  integer(i4)                                      :: lyf
  integer(i4)                                      :: lzf
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
  logical                                          :: lilocal
  logical                                          :: ljlocal
  logical                                          :: lklocal
  logical                                          :: lllocal
  logical                                          :: locals(2,6)
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
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
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
  real(dp)                                         :: d2cos1(6)
  real(dp)                                         :: d2cos2(6)
  real(dp)                                         :: d2cos3(6)
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
  real(dp)                                         :: d2cosdrij2
  real(dp)                                         :: d2cosdrik2
  real(dp)                                         :: d2cosdrijdrjk
  real(dp)                                         :: d2cosdrijdrik
  real(dp)                                         :: d2cosdrikdrjk
  real(dp)                                         :: disp
  real(dp)                                         :: dispt
  real(dp)                                         :: ddispdr
  real(dp)                                         :: d2dispdr2
  real(dp),    dimension(:),     allocatable, save :: dcn               ! First derivatives of log of coordination number w.r.t. coordination number
  real(dp),    dimension(:),     allocatable, save :: d2cn              ! Second derivatives of log of coordination number w.r.t. coordination number
  real(dp)                                         :: dr0dcni
  real(dp)                                         :: dr0dcnj
  real(dp)                                         :: dtheta
  real(dp)                                         :: dtmps(6)
  real(dp)                                         :: dtmps2(6)
  real(dp)                                         :: dtrm
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
  real(dp)                                         :: deijdcnh          ! Derivative of eij w.r.t. HB coordination number
  real(dp)                                         :: deijdcni          ! Derivative of eij w.r.t. coordination number of i
  real(dp)                                         :: d2eijdcni2        ! Second derivative of eij w.r.t. coordination number of i
  real(dp)                                         :: d2eijdcnh2        ! Second derivative of eij w.r.t. HB coordination number
  real(dp)                                         :: d2eijdcnhdr       ! Second derivative of eij w.r.t. HB coordination number and r
  real(dp)                                         :: d2eijdcnhdcni     ! Second derivative of eij w.r.t. HB coordination number and CN of i
  real(dp)                                         :: d2eijdcnhdcnj     ! Second derivative of eij w.r.t. HB coordination number and CN of j
  real(dp)                                         :: d2eijdcnidcnj     ! Second derivative of eij w.r.t. coordination numbers of i and j
  real(dp)                                         :: d2eijdcnidrij     ! Second derivative of eij w.r.t. coordination number of i and rij
  real(dp)                                         :: d2eijdcnjdrij     ! Second derivative of eij w.r.t. coordination number of j and rij
  real(dp)                                         :: deijdcnj          ! Derivative of eij w.r.t. coordination number of j
  real(dp)                                         :: d2eijdcnj2        ! Second derivative of eij w.r.t. coordination number of j
  real(dp)                                         :: deijddr           ! Derivative of eij w.r.t. dr
  real(dp)                                         :: d2eijddr2         ! Second derivative of eij w.r.t. dr
  real(dp)                                         :: deijdr            ! (1/r)(dEij/dr)
  real(dp)                                         :: d2eijdr2          ! (d2Eij/dr2)
  real(dp)                                         :: eij
  real(dp)                                         :: eijk
  real(dp)                                         :: eijknodamp
  real(dp)                                         :: deijkdcos
  real(dp)                                         :: d2eijkdcos2
  real(dp)                                         :: deijkdrij
  real(dp)                                         :: deijkdrik
  real(dp)                                         :: deijkdril
  real(dp)                                         :: deijkdrjk
  real(dp)                                         :: deijkdrjl
  real(dp)                                         :: eijkl
  real(dp)                                         :: eijklnod
  real(dp)                                         :: deijkldphi
  real(dp)                                         :: d2eijkldphi2
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
  real(dp)                                         :: ofctij
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
  real(dp)                                         :: repab_b
  real(dp)                                         :: repab_n
  real(dp)                                         :: rep_rho_b
  real(dp)                                         :: rep_rho_n
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
  real(dp)                                         :: sintheta
  real(dp)                                         :: sum0(9)
  real(dp)                                         :: sum1(9)
  real(dp)                                         :: t            ! Taper function
  real(dp)                                         :: dtdr         ! Taper function first derivative w.r.t. r divided by r
  real(dp)                                         :: dtdr2        ! Taper function first derivative w.r.t. r2
  real(dp)                                         :: d2tdr2       ! Taper function second derivative w.r.t. r
  real(dp)                                         :: d2tdr22      ! Taper function second derivative w.r.t. r2
  real(dp)                                         :: d3tdr23      ! Taper function third derivative w.r.t. r2
  real(dp)                                         :: tmp
  real(dp)                                         :: tmprr
  real(dp)                                         :: theta
  real(dp)                                         :: theta0
  real(dp)                                         :: time1
  real(dp)                                         :: time2
  real(dp)                                         :: time3
  real(dp)                                         :: time4
  real(dp)                                         :: trm6
  real(dp)                                         :: trm8
  real(dp)                                         :: p2
  real(dp)                                         :: p3
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
  call trace_in('gfnff_energyd_s')
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
  if (status/=0) call outofmemory('gfnff_energyd_s','cn')
  allocate(dcn(numat),stat=status)
  if (status/=0) call outofmemory('gfnff_energyd_s','dcn')
  allocate(d2cn(numat),stat=status)
  if (status/=0) call outofmemory('gfnff_energyd_s','d2cn')
  allocate(gw(maxc6ref,numat),stat=status)
  if (status/=0) call outofmemory('gfnff_energyd_s','gw')
!
  if (lgrad1) then
    allocate(dgwdcn(maxc6ref,numat),stat=status)
    if (status/=0) call outofmemory('gfnff_energyd_s','dgwdcn')
  else
    allocate(dgwdcn(maxc6ref,1),stat=status)
    if (status/=0) call outofmemory('gfnff_energyd_s','dgwdcn')
  endif
!
  if (lgrad2) then
    allocate(d2gwdcn2(maxc6ref,numat),stat=status)
    if (status/=0) call outofmemory('gfnff_energyd_s','d2gwdcn2')
  else
    allocate(d2gwdcn2(maxc6ref,1),stat=status)
    if (status/=0) call outofmemory('gfnff_energyd_s','d2gwdcn2')
  endif
!
  if (ldohbcn) then
    allocate(cn_hb(numat),stat=status)
    if (status/=0) call outofmemory('gfnff_energyd_s','cn_hb')
  endif
!
!  Set coordinates in the box
!
  call getinbox
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
  call gfnff_cn(numat,gfnff_cnthr,cn,dcn,d2cn,lgrad1.or.ldcharge,lgrad2)
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
  call eem_gfnff(ecoulomb,cn,dcn,d2cn,lgrad1,lgrad2,.false.)
!
  if (lverbose) then 
    time4 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for charges = '',f12.6)') time4 - time3
  endif
!********************************************************
!  Compute the coordination-dependent dispersion terms  *
!********************************************************
  call gfnff_cnc6(numat,nat,cn,gw,dgwdcn,d2gwdcn2,lgrad1,lgrad2)
!
  if (lverbose) then 
    time3 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for c6 terms from cn = '',f12.6)') time3 - time4
  endif
!****************
!  Bond energy  *
!****************
  ixf = - 2
  iyf = - 1
  izf =   0
  do i = 1,numat
    iloc = atom2local(i)
    lilocal = (iloc.ne.0)
!
!  Only do energy and first derivatives if i is local
!
    if (.not.lilocal.and..not.lgrad2) cycle
!
    nati = nat(i)
    nregioni = nregionno(nrelf2a(i))
    lslicei = lsliceatom(nsft+nrelf2a(i))
!
    if (lgrad2) then
      ixf = ixf + 3
      iyf = iyf + 3
      izf = izf + 3
      if (lilocal) then
        indi = 3*(iloc - 1)
        ix = indi + 1
        iy = indi + 2
        iz = indi + 3
      endif
    endif
!
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
!
!  Only compute half of bonds to avoid duplication between i->j and j->i
!
      if (j.lt.i) cycle
!
!  Set up quantities for this pair
!
      jloc = atom2local(j)
      ljlocal = (jloc.ne.0)
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
!  Derivative setup tasks
!
      if (lgrad1) then
        xji = xbnbr(ni,i)
        yji = ybnbr(ni,i)
        zji = zbnbr(ni,i)
        if (lgrad2) then
          indj = 3*(j-1)
          jxf = indj + 1
          jyf = indj + 2
          jzf = indj + 3
          if (ljlocal) then
            indj = 3*(jloc - 1)
            jx = indj + 1
            jy = indj + 2
            jz = indj + 3
          endif
        endif
      endif
!
!  Estimate corrected r0 based on coordination numbers and charges
!  NB: Second derivatives w.r.t. cni and cnj are zero since terms are linear
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
        if (lgrad1) then
          deijddr = - 2.0_dp*tmp*dr*eij   ! (dE/dr)
!
!  Coordination number derivatives
!
          if (lgrad2) then
            dtrm = tmp*eij*(2.0_dp - 4.0_dp*tmp*dr*dr)
!
            deijdcni = - deijddr*dr0dcni
            deijdcnj = - deijddr*dr0dcnj
            d2eijdcni2 = - dtrm*dr0dcni*dr0dcni
            d2eijdcnj2 = - dtrm*dr0dcnj*dr0dcnj
            d2eijdcnidcnj = - dtrm*dr0dcni*dr0dcnj
            d2eijdcnidrij = dtrm*dr0dcni/r
            d2eijdcnjdrij = dtrm*dr0dcnj/r
!
            call gfnff_drv2_dcnd(lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,j,jx,jy,jz,jxf,jyf,jzf,xji,yji,zji, &
                                 deijdcni,deijdcnj,d2eijdcni2,d2eijdcnj2,d2eijdcnidcnj, &
                                 d2eijdcnidrij,d2eijdcnjdrij,dcn,d2cn,lgrad2)
          else
            deijdcni = - deijddr*dr0dcni
            deijdcnj = - deijddr*dr0dcnj
            call gfnff_drv_dcn(i,deijdcni,cn,dcn)
            call gfnff_drv_dcn(j,deijdcnj,cn,dcn)
          endif
!
          deijddr = deijddr/r   ! (1/r(dE/dr))
!
          if (lilocal) then
            xdrv(i) = xdrv(i) - deijddr*xji
            ydrv(i) = ydrv(i) - deijddr*yji
            zdrv(i) = zdrv(i) - deijddr*zji
!
            xdrv(j) = xdrv(j) + deijddr*xji
            ydrv(j) = ydrv(j) + deijddr*yji
            zdrv(j) = zdrv(j) + deijddr*zji
          endif
!
          if (lstr.or.lgrad2) then
            call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
          endif
!
          if (lstr) then
            if (lilocal) then
              do kl = 1,nstrains
                ks = nstrptr(kl)
                rstrd(kl) = rstrd(kl) + deijddr*dr2ds(ks)
              enddo
            endif
          endif
!
!  Second derivatives
!
          if (lgrad2) then
            d2eijddr2 = (4.0_dp*tmp*dr*dr + 2.0_dp*dr/r - 2.0_dp)*tmp*eij/r**2   ! (1/r.d(1/r.dE/dr)/dr)
!
            if (lilocal) then
              derv2(jxf,ix) = derv2(jxf,ix) - d2eijddr2*d2r2dx2(1)
              derv2(jyf,ix) = derv2(jyf,ix) - d2eijddr2*d2r2dx2(6)
              derv2(jzf,ix) = derv2(jzf,ix) - d2eijddr2*d2r2dx2(5)
              derv2(jxf,iy) = derv2(jxf,iy) - d2eijddr2*d2r2dx2(6)
              derv2(jyf,iy) = derv2(jyf,iy) - d2eijddr2*d2r2dx2(2)
              derv2(jzf,iy) = derv2(jzf,iy) - d2eijddr2*d2r2dx2(4)
              derv2(jxf,iz) = derv2(jxf,iz) - d2eijddr2*d2r2dx2(5)
              derv2(jyf,iz) = derv2(jyf,iz) - d2eijddr2*d2r2dx2(4)
              derv2(jzf,iz) = derv2(jzf,iz) - d2eijddr2*d2r2dx2(3)
              derv2(jxf,ix) = derv2(jxf,ix) - deijddr
              derv2(jyf,iy) = derv2(jyf,iy) - deijddr
              derv2(jzf,iz) = derv2(jzf,iz) - deijddr
            endif
!
            if (ljlocal) then
              derv2(ixf,jx) = derv2(ixf,jx) - d2eijddr2*d2r2dx2(1)
              derv2(iyf,jx) = derv2(iyf,jx) - d2eijddr2*d2r2dx2(6)
              derv2(izf,jx) = derv2(izf,jx) - d2eijddr2*d2r2dx2(5)
              derv2(ixf,jy) = derv2(ixf,jy) - d2eijddr2*d2r2dx2(6)
              derv2(iyf,jy) = derv2(iyf,jy) - d2eijddr2*d2r2dx2(2)
              derv2(izf,jy) = derv2(izf,jy) - d2eijddr2*d2r2dx2(4)
              derv2(ixf,jz) = derv2(ixf,jz) - d2eijddr2*d2r2dx2(5)
              derv2(iyf,jz) = derv2(iyf,jz) - d2eijddr2*d2r2dx2(4)
              derv2(izf,jz) = derv2(izf,jz) - d2eijddr2*d2r2dx2(3)
              derv2(ixf,jx) = derv2(ixf,jx) - deijddr
              derv2(iyf,jy) = derv2(iyf,jy) - deijddr
              derv2(izf,jz) = derv2(izf,jz) - deijddr
            endif
!
            if (lstr) then
!
!  Mixed derivatives
!
              if (lilocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(ix,kl) = derv3(ix,kl) - deijddr*d2r2dsdx(ks,1) - xji*d2eijddr2*dr2ds(ks)
                  derv3(iy,kl) = derv3(iy,kl) - deijddr*d2r2dsdx(ks,2) - yji*d2eijddr2*dr2ds(ks)
                  derv3(iz,kl) = derv3(iz,kl) - deijddr*d2r2dsdx(ks,3) - zji*d2eijddr2*dr2ds(ks)
                enddo
              endif
              if (ljlocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(jx,kl) = derv3(jx,kl) + deijddr*d2r2dsdx(ks,1) + xji*d2eijddr2*dr2ds(ks)
                  derv3(jy,kl) = derv3(jy,kl) + deijddr*d2r2dsdx(ks,2) + yji*d2eijddr2*dr2ds(ks)
                  derv3(jz,kl) = derv3(jz,kl) + deijddr*d2r2dsdx(ks,3) + zji*d2eijddr2*dr2ds(ks)
                enddo
              endif
!
!  Strain-strain
!
              if (lilocal) then
                do kk = 1,nstrains
                  ks = nstrptr(kk)
                  do kl = 1,nstrains
                    kt = nstrptr(kl)
                    sderv2(kl,kk) = sderv2(kl,kk) + d2eijddr2*dr2ds(kt)*dr2ds(ks) + deijddr*d2r2ds2(kt,ks)
                  enddo
                enddo
              endif
            endif
          endif
!
!  Get HB coordination number derivatives for hbH
!
          deijdcnh = dr*dr*p2*eij*(1.0_dp - gfnff_bondscale)
          if (lgrad2) then
            d2eijdcnh2 = dr*dr*p2*deijdcnh*(1.0_dp - gfnff_bondscale)
            dtrm = 2.0_dp*eij*p2*dr*(1.0_dp - dr*dr*tmp)*(1.0_dp - gfnff_bondscale)
            d2eijdcnhdcni = - dtrm*dr0dcni
            d2eijdcnhdcnj = - dtrm*dr0dcnj
            d2eijdcnhdr = dtrm/r
            call gfnff_drv2_dcnd_hb(hbH,i,j,xji,yji,zji,deijdcnh,d2eijdcnh2,d2eijdcnhdr,d2eijdcnhdcni, &
                                    d2eijdcnhdcnj,dcn,lgrad2)
          else
            call gfnff_drv_dcn_hb(hbH,deijdcnh)
          endif
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
          if (lgrad2) then
            dtrm = p2*eij*(2.0_dp - 4.0_dp*p2*dr*dr)
!
            deijdcni = - deijddr*dr0dcni
            deijdcnj = - deijddr*dr0dcnj
            d2eijdcni2 = - dtrm*dr0dcni*dr0dcni
            d2eijdcnj2 = - dtrm*dr0dcnj*dr0dcnj
            d2eijdcnidcnj = - dtrm*dr0dcni*dr0dcnj
            d2eijdcnidrij = dtrm*dr0dcni/r
            d2eijdcnjdrij = dtrm*dr0dcnj/r
!
            call gfnff_drv2_dcnd(lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,j,jx,jy,jz,jxf,jyf,jzf,xji,yji,zji, &
                                 deijdcni,deijdcnj,d2eijdcni2,d2eijdcnj2,d2eijdcnidcnj, &
                                 d2eijdcnidrij,d2eijdcnjdrij,dcn,d2cn,lgrad2)
          else
            deijdcni = - deijddr*dr0dcni
            deijdcnj = - deijddr*dr0dcnj
            call gfnff_drv_dcn(i,deijdcni,cn,dcn)
            call gfnff_drv_dcn(j,deijdcnj,cn,dcn)
          endif
!
          deijddr = deijddr/r   ! (1/r(dE/dr))
!
          if (lilocal) then
            xdrv(i) = xdrv(i) - deijddr*xji
            ydrv(i) = ydrv(i) - deijddr*yji
            zdrv(i) = zdrv(i) - deijddr*zji
!
            xdrv(j) = xdrv(j) + deijddr*xji
            ydrv(j) = ydrv(j) + deijddr*yji
            zdrv(j) = zdrv(j) + deijddr*zji
          endif
!
          if (lstr.or.lgrad2) then
            call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
          endif
          if (lstr) then
            if (lilocal) then
              do kl = 1,nstrains
                ks = nstrptr(kl)
                rstrd(kl) = rstrd(kl) + deijddr*dr2ds(ks)
              enddo
            endif
          endif
!
!  Second derivatives
!
          if (lgrad2) then
            d2eijddr2 = (4.0_dp*p2*dr*dr + 2.0_dp*dr/r - 2.0_dp)*p2*eij/r**2   ! (1/r.d(1/r.dE/dr)/dr)
!
            if (lilocal) then
              derv2(jxf,ix) = derv2(jxf,ix) - d2eijddr2*d2r2dx2(1)
              derv2(jyf,ix) = derv2(jyf,ix) - d2eijddr2*d2r2dx2(6)
              derv2(jzf,ix) = derv2(jzf,ix) - d2eijddr2*d2r2dx2(5)
              derv2(jxf,iy) = derv2(jxf,iy) - d2eijddr2*d2r2dx2(6)
              derv2(jyf,iy) = derv2(jyf,iy) - d2eijddr2*d2r2dx2(2)
              derv2(jzf,iy) = derv2(jzf,iy) - d2eijddr2*d2r2dx2(4)
              derv2(jxf,iz) = derv2(jxf,iz) - d2eijddr2*d2r2dx2(5)
              derv2(jyf,iz) = derv2(jyf,iz) - d2eijddr2*d2r2dx2(4)
              derv2(jzf,iz) = derv2(jzf,iz) - d2eijddr2*d2r2dx2(3)
              derv2(jxf,ix) = derv2(jxf,ix) - deijddr
              derv2(jyf,iy) = derv2(jyf,iy) - deijddr
              derv2(jzf,iz) = derv2(jzf,iz) - deijddr
            endif
!
            if (ljlocal) then
              derv2(ixf,jx) = derv2(ixf,jx) - d2eijddr2*d2r2dx2(1)
              derv2(iyf,jx) = derv2(iyf,jx) - d2eijddr2*d2r2dx2(6)
              derv2(izf,jx) = derv2(izf,jx) - d2eijddr2*d2r2dx2(5)
              derv2(ixf,jy) = derv2(ixf,jy) - d2eijddr2*d2r2dx2(6)
              derv2(iyf,jy) = derv2(iyf,jy) - d2eijddr2*d2r2dx2(2)
              derv2(izf,jy) = derv2(izf,jy) - d2eijddr2*d2r2dx2(4)
              derv2(ixf,jz) = derv2(ixf,jz) - d2eijddr2*d2r2dx2(5)
              derv2(iyf,jz) = derv2(iyf,jz) - d2eijddr2*d2r2dx2(4)
              derv2(izf,jz) = derv2(izf,jz) - d2eijddr2*d2r2dx2(3)
              derv2(ixf,jx) = derv2(ixf,jx) - deijddr
              derv2(iyf,jy) = derv2(iyf,jy) - deijddr
              derv2(izf,jz) = derv2(izf,jz) - deijddr
            endif
!
            if (lstr) then
!
!  Mixed derivatives
!
              if (lilocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(ix,kl) = derv3(ix,kl) - deijddr*d2r2dsdx(ks,1) - xji*d2eijddr2*dr2ds(ks)
                  derv3(iy,kl) = derv3(iy,kl) - deijddr*d2r2dsdx(ks,2) - yji*d2eijddr2*dr2ds(ks)
                  derv3(iz,kl) = derv3(iz,kl) - deijddr*d2r2dsdx(ks,3) - zji*d2eijddr2*dr2ds(ks)
                enddo
              endif
              if (ljlocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(jx,kl) = derv3(jx,kl) + deijddr*d2r2dsdx(ks,1) + xji*d2eijddr2*dr2ds(ks)
                  derv3(jy,kl) = derv3(jy,kl) + deijddr*d2r2dsdx(ks,2) + yji*d2eijddr2*dr2ds(ks)
                  derv3(jz,kl) = derv3(jz,kl) + deijddr*d2r2dsdx(ks,3) + zji*d2eijddr2*dr2ds(ks)
                enddo
              endif
!
!  Strain-strain
!
              if (lilocal) then
                do kk = 1,nstrains
                  ks = nstrptr(kk)
                  do kl = 1,nstrains
                    kt = nstrptr(kl)
                    sderv2(kl,kk) = sderv2(kl,kk) + d2eijddr2*dr2ds(kt)*dr2ds(ks) + deijddr*d2r2ds2(kt,ks)
                  enddo
                enddo
              endif
            endif
          endif
        endif
      endif
      if (lilocal) then
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
      endif
    enddo
  enddo
!
  if (lverbose) then 
    time4 = g_cpu_time()
    write(ioout,'(''  GFNFF: Time for bond energy = '',f12.6)') time4 - time3
  endif
!************************************
!  Long-range terms for dispersion  *
!************************************
!
!  Initialise lattice vectors for this cutoff
!
  cutdr = gfnff_dispthr
  call gfnff_setpaircell(cutdr,lr_paircell)
!
!  Dispersion part 1
!
  if (lgrad1) then
    if (ndim.eq.0) then
      imin = 2
    else
      imin = 1
    endif
    do i = imin,numat
      nati = nat(i)
      iloc = atom2local(i)
      lilocal = (iloc.ne.0)
!
!  Only do energy and first derivatives if i is local
!
      if (.not.lilocal.and..not.lgrad2) cycle
!
      if (ndim.eq.0) then
        jmax = i - 1
      else
        jmax = i
      endif
!
      if (lgrad2) then
        indi = 3*(i-1)
        ixf = indi + 1
        iyf = indi + 2
        izf = indi + 3
        if (lilocal) then
          indi = 3*(iloc-1)
          ix = indi + 1
          iy = indi + 2
          iz = indi + 3
        endif
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
        jloc = atom2local(j)
        ljlocal = (jloc.ne.0)
        natj = nat(j)
        if (nati.ge.natj) then
          indn = nati*(nati-1)/2 + natj
        else
          indn = natj*(natj-1)/2 + nati
        endif
!
        if (lgrad2) then
          indj = 3*(j-1)
          jxf = indj + 1
          jyf = indj + 2
          jzf = indj + 3
          if (ljlocal) then
            indj = 3*(jloc-1)
            jx = indj + 1
            jy = indj + 2
            jz = indj + 3
          endif
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
!  Find interactions for this pair
!
        call gfnff_getpairs(i,j,cutdr,lr_paircell,lr_pairs)
!
!  Set distance-independent constants for dispersion
!
        r4r2ij = 3.0_dp*d4_sqrtZr4r2(nati)*d4_sqrtZr4r2(natj)
        r0 = d4_r0(indn)
        r03 = r0**3
        r04 = r0**4
!
!  Set c6 for this pair
!
        zetaprod = d4_zeta_c6(i)*d4_zeta_c6(j)
        call gfnff_setc6(numat,nati,natj,i,j,gw,dgwdcn,d2gwdcn2,c6,dc6dcni,dc6dcnj, &
                         d2c6dcni2,d2c6dcnj2,d2c6dcnidcnj,.true.,lgrad2)
        dispt = 0.0_dp
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
!
!  Taper
!
          if (r2.gt.t_dispthr) then
            call mdftaper(r2,t_dispthr,gfnff_dispthr,t,dtdr2,d2tdr22,d3tdr23,lgrad1,lgrad2,.false.)
            dtdr = 2.0_dp*dtdr2
            if (lgrad2) d2tdr2 = 4.0_dp*d2tdr22
          else
            t = 1.0_dp
            dtdr = 0.0_dp
            d2tdr2 = 0.0_dp
          endif
!
!  Compute dispersion terms
!
          trm6 = 1.0_dp/(r2**3 + r03)
          trm8 = 1.0_dp/(r2**4 + r04)
!
          disp = (trm6 + 2.0_dp*r4r2ij*trm8)*zetaprod
          dispt = dispt + disp*t
!
          if (lgrad2) then
            ddispdr = - r2*r2*(6.0_dp*trm6*trm6 + 16.0_dp*r4r2ij*r2*trm8*trm8)*zetaprod
            d2eijdcnidrij = - fct*(ddispdr*t + disp*dtdr)*dc6dcni
            d2eijdcnjdrij = - fct*(ddispdr*t + disp*dtdr)*dc6dcnj
            if (abs(d2eijdcnidrij).gt.gfnff_cnc6tol.or.abs(d2eijdcnjdrij).gt.gfnff_cnc6tol) then
              call gfnff_drv2_dcnd_x2(i,ixf,iyf,izf,j,jxf,jyf,jzf,xji,yji,zji, &
                                      d2eijdcnidrij,d2eijdcnjdrij,dcn)
            endif
          endif
        enddo
!
!  Collective derivative sums
!
        if (lgrad2) then
          deijdcni = - fct*dc6dcni*dispt
          deijdcnj = - fct*dc6dcnj*dispt
          d2eijdcni2 = - fct*d2c6dcni2*dispt
          d2eijdcnj2 = - fct*d2c6dcnj2*dispt
          d2eijdcnidcnj = - fct*d2c6dcnidcnj*dispt
!
!  Tolerance check
!
          if (abs(d2eijdcni2).gt.gfnff_cnc6tol.or.abs(d2eijdcnj2).gt.gfnff_cnc6tol) then
            call gfnff_drv2_dcnd_nox(lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,j,jx,jy,jz,jxf,jyf,jzf, &
                                     deijdcni,deijdcnj,d2eijdcni2,d2eijdcnj2,d2eijdcnidcnj, &
                                     dcn,d2cn,lgrad2)
          else
            call gfnff_drv_dcn(i,deijdcni,cn,dcn)
            call gfnff_drv_dcn(j,deijdcnj,cn,dcn)
          endif
        else
          deijdcni = - fct*dc6dcni*dispt
          call gfnff_drv_dcn(i,deijdcni,cn,dcn)
          deijdcnj = - fct*dc6dcnj*dispt
          call gfnff_drv_dcn(j,deijdcnj,cn,dcn)
        endif
      enddo
    enddo
  endif
!
!  Dispersion part 2
!
  ix = - 2
  iy = - 1
  iz =   0
  do iloc = 1,natomsonnode
    i = node2atom(iloc)
    nati = nat(i)
    nregioni = nregionno(nrelf2a(i))
    lslicei = lsliceatom(nsft+nrelf2a(i))
!
    if (lgrad2) then
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
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
    jx = - 2
    jy = - 1
    jz =   0
    do j = 1,numat
      natj = nat(j)
      nregionj = nregionno(nrelf2a(j))
      lslicej = lsliceatom(nsft+nrelf2a(j))
      if (nati.ge.natj) then
        indn = nati*(nati-1)/2 + natj
      else
        indn = natj*(natj-1)/2 + nati
      endif
!
      if (lgrad2) then
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
      endif
!
!  Correction for i = j
!
      if (i.eq.j) then
        fct = 0.5_dp
        ofctij = 1.0_dp
      else
        fct = 1.0_dp
        ofctij = 0.5_dp
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
      r03 = r0**3
      r04 = r0**4
!
!  Set c6 for this pair
!
      zetaprod = d4_zeta_c6(i)*d4_zeta_c6(j)
      call gfnff_setc6(numat,nati,natj,i,j,gw,dgwdcn,d2gwdcn2,c6,dc6dcni,dc6dcnj, &
                       d2c6dcni2,d2c6dcnj2,d2c6dcnidcnj,lgrad1,lgrad2)
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
!
!  Taper
!
        if (r2.gt.t_dispthr) then
          call mdftaper(r2,t_dispthr,gfnff_dispthr,t,dtdr2,d2tdr22,d3tdr23,lgrad1,lgrad2,.false.)
          if (lgrad1) dtdr = 2.0_dp*dtdr2
          if (lgrad2) d2tdr2 = 4.0_dp*d2tdr22
        else
          t = 1.0_dp
          dtdr = 0.0_dp
          d2tdr2 = 0.0_dp
        endif
!
!  Compute dispersion terms
!
        trm6 = 1.0_dp/(r2**3 + r03)
        trm8 = 1.0_dp/(r2**4 + r04)
!
        disp = (trm6 + 2.0_dp*r4r2ij*trm8)*zetaprod
        eij = - fct*c6*disp*t
!
!  Handle regions
!
        if (lreg2one) then
          esregion12 = esregion12 + eij*ofctij
        elseif (lreg2pair) then
          esregion2 = esregion2 + eij*ofctij
        else
          edispersion = edispersion + eij*ofctij
        endif
!
!  Region - region energy
!
        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eij*ofctij
!
!  Attachment energy
!
        if (lattach) eattach = eattach + eij*ofctij
!
!  Site energies
!
        siteenergy(i) = siteenergy(i) + 0.5_dp*eij*ofctij
        siteenergy(j) = siteenergy(j) + 0.5_dp*eij*ofctij
!
        if (lgrad1) then
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
          if (lstr.or.lgrad2) then
            call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
          endif
          if (lstr) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + ofctij*deijdr*dr2ds(ks)
            enddo
          endif
!
!  Second derivatives
!
          if (lgrad2) then
            d2dispdr2 = r2*(24.0_dp*trm6*trm6*(3.0_dp*r2*r2*r2*trm6 - 1.0_dp) + &
                            16.0_dp*r4r2ij*r2*trm8*trm8*(16.0_dp*trm8*r2*r2*r2*r2 - 6.0_dp))*zetaprod
            d2eijdr2 = - fct*c6*(d2dispdr2*t + 2.0_dp*ddispdr*dtdr + disp*d2tdr2)   ! (1/r.d(1/r.dE/dr)/dr)
!
            derv2(jx,ix) = derv2(jx,ix) - d2eijdr2*d2r2dx2(1)
            derv2(jy,ix) = derv2(jy,ix) - d2eijdr2*d2r2dx2(6)
            derv2(jz,ix) = derv2(jz,ix) - d2eijdr2*d2r2dx2(5)
            derv2(jx,iy) = derv2(jx,iy) - d2eijdr2*d2r2dx2(6)
            derv2(jy,iy) = derv2(jy,iy) - d2eijdr2*d2r2dx2(2)
            derv2(jz,iy) = derv2(jz,iy) - d2eijdr2*d2r2dx2(4)
            derv2(jx,iz) = derv2(jx,iz) - d2eijdr2*d2r2dx2(5)
            derv2(jy,iz) = derv2(jy,iz) - d2eijdr2*d2r2dx2(4)
            derv2(jz,iz) = derv2(jz,iz) - d2eijdr2*d2r2dx2(3)
            derv2(jx,ix) = derv2(jx,ix) - deijdr
            derv2(jy,iy) = derv2(jy,iy) - deijdr
            derv2(jz,iz) = derv2(jz,iz) - deijdr
!
            if (lstr) then
!
!  Mixed derivatives
!
              if (i.ne.j) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(ix,kl) = derv3(ix,kl) - deijdr*d2r2dsdx(ks,1) - xji*d2eijdr2*dr2ds(ks)
                  derv3(iy,kl) = derv3(iy,kl) - deijdr*d2r2dsdx(ks,2) - yji*d2eijdr2*dr2ds(ks)
                  derv3(iz,kl) = derv3(iz,kl) - deijdr*d2r2dsdx(ks,3) - zji*d2eijdr2*dr2ds(ks)
                enddo
              endif
!
!  Strain-strain
!
              do kk = 1,nstrains
                ks = nstrptr(kk)
                do kl = 1,nstrains
                  kt = nstrptr(kl)
                  sderv2(kl,kk) = sderv2(kl,kk) + ofctij*(d2eijdr2*dr2ds(kt)*dr2ds(ks) + deijdr*d2r2ds2(kt,ks))
                enddo
              enddo
            endif
          endif
!
!  Coordination number derivatives of c6
!
          if (lgrad2) then
            d2eijdcnidrij = - fct*(ddispdr*t + disp*dtdr)*dc6dcni
            d2eijdcnjdrij = - fct*(ddispdr*t + disp*dtdr)*dc6dcnj
!
!  Tolerance check
!
            if (abs(d2eijdcnidrij).gt.gfnff_cnc6tol.or.abs(d2eijdcnjdrij).gt.gfnff_cnc6tol) then
              call gfnff_drv2_dcnd_x1(i,ix,iy,iz,j,jx,jy,jz,xji,yji,zji, &
                                      d2eijdcnidrij,d2eijdcnjdrij,dcn,ofctij)
            endif
          endif
        endif
      enddo
    enddo
  enddo
!***********************************
!  Long-range terms for repulsion  *
!***********************************
!
!  Initialise lattice vectors for this cutoff
!
  cutdr = gfnff_repthr
  call gfnff_setpaircell(cutdr,lr_paircell)
!
  do ii = 1,natomsonnode
    i = node2atom(ii)
    nati = nat(i)
    nregioni = nregionno(nrelf2a(i))
    lslicei = lsliceatom(nsft+nrelf2a(i))
!
    if (lgrad2) then
      indi = 3*(ii-1)
      ix = indi + 1
      iy = indi + 2
      iz = indi + 3
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
    do j = 1,numat
      if (i.ge.j) then
        ind = i*(i-1)/2 + j
      else
        ind = j*(j-1)/2 + i
      endif
      natj = nat(j)
      nregionj = nregionno(nrelf2a(j))
      lslicej = lsliceatom(nsft+nrelf2a(j))
!
      if (lgrad2) then
        indj = 3*(j-1)
        jx = indj + 1
        jy = indj + 2
        jz = indj + 3
      endif
!
!  Correction for i = j
!
      if (i.eq.j) then
        fct = 0.5_dp
        ofctij = 1.0_dp
      else
        fct = 1.0_dp
        ofctij = 0.5_dp
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
!  Set distance-independent constants for repulsion
!
      rep_rho_b = sqrt(gfnff_repulsion_a(nati)*gfnff_repulsion_a(natj))
      repab_b = gfnff_repulsion_z(nati)*gfnff_repulsion_z(natj)*gfnff_repscale_b
      rep_rho_n = gfnff_repulsion_p(ind)
      repab_n = gfnff_repulsion_z(nati)*gfnff_repulsion_z(natj)*gfnff_repscale_n
!
!  Loop over valid interactions (self term is excluded by getpairs)
!
      do np = 1,lr_pairs%npair
        r2 = lr_pairs%r2pair(np)
        xji = lr_pairs%xpair(np)
        yji = lr_pairs%ypair(np)
        zji = lr_pairs%zpair(np)
!!!!!!!!!!!!!!!!!!!!!!
!  Repulsive energy  !
!!!!!!!!!!!!!!!!!!!!!!
!
!  Taper
!
        if (r2.gt.t_repthr) then
          call mdftaper(r2,t_repthr,gfnff_repthr,t,dtdr2,d2tdr22,d3tdr23,lgrad1,lgrad2,.false.)
          if (lgrad1) dtdr = 2.0_dp*dtdr2
          if (lgrad2) d2tdr2 = 4.0_dp*d2tdr22
        else
          t = 1.0_dp
          dtdr = 0.0_dp
          d2tdr2 = 0.0_dp
        endif
        ff23 = 1.0_dp
!
!  Check whether this is a bonded pair
!
        call gfnff_bond_check(j,xji,yji,zji,1_i4,lbonded)
        if (lbonded) then
          repab = repab_b
          rep_rholoc = rep_rho_b
        else
          repab = repab_n
!
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
          rep_rholoc = rep_rho_n*ff23
        endif
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
          esregion12 = esregion12 + eij*ofctij
        elseif (lreg2pair) then
          esregion2 = esregion2 + eij*ofctij
        else
          erepulsion = erepulsion + eij*ofctij
        endif
!
!  Region - region energy
!
        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + eij*ofctij
!
!  Attachment energy
!
        if (lattach) eattach = eattach + eij*ofctij
!
!  Site energies
!
        siteenergy(i) = siteenergy(i) + 0.5_dp*eij*ofctij
        siteenergy(j) = siteenergy(j) + 0.5_dp*eij*ofctij
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
          if (lstr.or.lgrad2) then
            call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
          endif
          if (lstr) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijdr*dr2ds(ks)*ofctij
            enddo
          endif
!
!  Second derivatives
!
          if (lgrad2) then
            d2eijdr2 = tmprr*t*(3.0_dp/r**4 + (15.0_dp/4.0_dp)*rep_rholoc/r**2.5_dp + (9.0_dp/4.0_dp)*rep_rholoc*rep_rholoc/r) &
                       - tmprr*dtdr*(2.0_dp/r**2 + 3.0_dp*rep_rholoc/r**0.5_dp) + tmprr*d2tdr2 ! (1/r.d(1/r.dE/dr)/dr)
            d2eijdr2 = fct*d2eijdr2
!
            derv2(jx,ix) = derv2(jx,ix) - d2eijdr2*d2r2dx2(1)
            derv2(jy,ix) = derv2(jy,ix) - d2eijdr2*d2r2dx2(6)
            derv2(jz,ix) = derv2(jz,ix) - d2eijdr2*d2r2dx2(5)
            derv2(jx,iy) = derv2(jx,iy) - d2eijdr2*d2r2dx2(6)
            derv2(jy,iy) = derv2(jy,iy) - d2eijdr2*d2r2dx2(2)
            derv2(jz,iy) = derv2(jz,iy) - d2eijdr2*d2r2dx2(4)
            derv2(jx,iz) = derv2(jx,iz) - d2eijdr2*d2r2dx2(5)
            derv2(jy,iz) = derv2(jy,iz) - d2eijdr2*d2r2dx2(4)
            derv2(jz,iz) = derv2(jz,iz) - d2eijdr2*d2r2dx2(3)
            derv2(jx,ix) = derv2(jx,ix) - deijdr
            derv2(jy,iy) = derv2(jy,iy) - deijdr
            derv2(jz,iz) = derv2(jz,iz) - deijdr
!
            if (lstr) then
!
!  Mixed derivatives
!
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(ix,kl) = derv3(ix,kl) - deijdr*d2r2dsdx(ks,1) - xji*d2eijdr2*dr2ds(ks)
                derv3(iy,kl) = derv3(iy,kl) - deijdr*d2r2dsdx(ks,2) - yji*d2eijdr2*dr2ds(ks)
                derv3(iz,kl) = derv3(iz,kl) - deijdr*d2r2dsdx(ks,3) - zji*d2eijdr2*dr2ds(ks)
              enddo
!
!  Strain-strain
!
              do kk = 1,nstrains
                ks = nstrptr(kk)
                do kl = 1,nstrains
                  kt = nstrptr(kl)
                  sderv2(kl,kk) = sderv2(kl,kk) + ofctij*(d2eijdr2*dr2ds(kt)*dr2ds(ks) + deijdr*d2r2ds2(kt,ks))
                enddo
              enddo
            endif
          endif
        endif
      enddo
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
  do ia = 1,n_gfnff_angles
    i  = n_gfnff_angleptr(1,ia)
    iloc = atom2local(i)
    lilocal = (iloc.ne.0)
!
!  Only do energy and first derivatives if i is local
!
    if (.not.lilocal.and..not.lgrad2) cycle
!
    nj = n_gfnff_angleptr(2,ia)
    nk = n_gfnff_angleptr(3,ia)
    j = nbrno_bond(nj,i)
    k = nbrno_bond(nk,i)
    jloc = atom2local(j)
    kloc = atom2local(k)
    ljlocal = (jloc.ne.0)
    lklocal = (kloc.ne.0)
!
!  Are any of the atoms local to this node?
!
    if (.not.lilocal.and..not.ljlocal.and..not.lklocal) cycle
!
    if (lgrad2) then
      indi = 3*(i-1)
      ixf = indi + 1
      iyf = indi + 2
      izf = indi + 3
      if (lilocal) then
        indi = 3*(iloc-1)
        ix = indi + 1
        iy = indi + 2
        iz = indi + 3
      endif
!
      indj = 3*(j-1)
      jxf = indj + 1
      jyf = indj + 2
      jzf = indj + 3
      if (ljlocal) then
        indj = 3*(jloc-1)
        jx = indj + 1
        jy = indj + 2
        jz = indj + 3
      endif
!
      indk = 3*(k-1)
      kxf = indk + 1
      kyf = indk + 2
      kzf = indk + 3
      if (lklocal) then
        indk = 3*(kloc-1)
        kx = indk + 1
        ky = indk + 2
        kz = indk + 3
      endif
    endif
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
    call gfnff_damp(r2ij,rcut2ij,dampij,ddampij,d2dampij,lgrad1,lgrad2)
    call gfnff_damp(r2ik,rcut2ik,dampik,ddampik,d2dampik,lgrad1,lgrad2)
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
        sintheta = sin(theta)
        deijkdcos = - 2.0_dp*rkang*dtheta/sintheta
        d2eijkdcos2 = 2.0_dp*rkang*(1.0_dp/sintheta**2 - dtheta/cos(theta))
      else
        deijkdcos = 0.0_dp
        d2eijkdcos2 = 0.0_dp
      endif
    else
      costheta0 = cos(theta0)
      eijknodamp = rkang*(cosang - costheta0)**2
      eijk = damp*eijknodamp
      deijkdcos = 2.0_dp*rkang*(cosang - costheta0)
      d2eijkdcos2 = 2.0_dp*rkang
    endif
    if (lilocal) then
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
    endif
!----------------------
!  First derivatives  |
!----------------------
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
      deijkdrij = damp*deijkdcos*dcosdrij      ! (1/r(dEdr))
      deijkdrij = deijkdrij + eijknodamp*dampik*ddampij   ! Add on damping derivatives for i-j
!
      if (lilocal) then
        xdrv(i) = xdrv(i) - deijkdrij*vij(1)
        ydrv(i) = ydrv(i) - deijkdrij*vij(2)
        zdrv(i) = zdrv(i) - deijkdrij*vij(3)
        xdrv(j) = xdrv(j) + deijkdrij*vij(1)
        ydrv(j) = ydrv(j) + deijkdrij*vij(2)
        zdrv(j) = zdrv(j) + deijkdrij*vij(3)
      endif
!
      if (lstr.or.lgrad2) then
        call real1strterm(ndim,vij(1),vij(2),vij(3),0.0_dp,0.0_dp,0.0_dp,dr2ds6(1,1),d2r2dx26(1,1), &
                          d2r2dsdx6(1,1,1),d2r2ds26(1,1,1),lgrad2)
      endif
      if (lilocal) then
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds6(ks,1)
          enddo
        endif
      endif
!
!  i-k
!
      deijkdrik = damp*deijkdcos*dcosdrik      ! (1/r(dEdr))
      deijkdrik = deijkdrik + eijknodamp*dampij*ddampik   ! Add on damping derivatives for i-k
!
      if (lilocal) then
        xdrv(i) = xdrv(i) - deijkdrik*vik(1)
        ydrv(i) = ydrv(i) - deijkdrik*vik(2)
        zdrv(i) = zdrv(i) - deijkdrik*vik(3)
        xdrv(k) = xdrv(k) + deijkdrik*vik(1)
        ydrv(k) = ydrv(k) + deijkdrik*vik(2)
        zdrv(k) = zdrv(k) + deijkdrik*vik(3)
      endif
!
      if (lstr.or.lgrad2) then
        call real1strterm(ndim,vik(1),vik(2),vik(3),0.0_dp,0.0_dp,0.0_dp,dr2ds6(1,2),d2r2dx26(1,2), &
                          d2r2dsdx6(1,1,2),d2r2ds26(1,1,2),lgrad2)
      endif
!
      if (lilocal) then
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdrik*dr2ds6(ks,2)
          enddo
        endif
      endif
!
!  j-k
!
      deijkdrjk = damp*deijkdcos*dcosdrjk      ! (1/r(dEdr))
      if (lilocal) then
        xdrv(j) = xdrv(j) - deijkdrjk*vjk(1)
        ydrv(j) = ydrv(j) - deijkdrjk*vjk(2)
        zdrv(j) = zdrv(j) - deijkdrjk*vjk(3)
        xdrv(k) = xdrv(k) + deijkdrjk*vjk(1)
        ydrv(k) = ydrv(k) + deijkdrjk*vjk(2)
        zdrv(k) = zdrv(k) + deijkdrjk*vjk(3)
      endif
!
      if (lstr.or.lgrad2) then
        call real1strterm(ndim,vjk(1),vjk(2),vjk(3),0.0_dp,0.0_dp,0.0_dp,dr2ds6(1,3),d2r2dx26(1,3), &
                          d2r2dsdx6(1,1,3),d2r2ds26(1,1,3),lgrad2)
      endif
!
      if (lilocal) then
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdrjk*dr2ds6(ks,3)
          enddo
        endif
      endif
!
      if (lgrad2) then
!-----------------------
!  Second derivatives  |
!-----------------------
        d2cosdrij2 = - rrij*rrij*(dcosdrij + rrij*rrik - 2.0_dp*cosang*rrij**2)
        d2cosdrijdrik = - rrij*rrik**3 - rrij*rrij*dcosdrik
        d2cosdrijdrjk = rrik*rrij**3
        d2cosdrik2 = - rrik*rrik*(dcosdrik + rrij*rrik - 2.0_dp*cosang*rrik**2)
        d2cosdrikdrjk = rrij*rrik**3
!
!  i-j / i-j
!
        d2trm = dampik*(dampij*(d2eijkdcos2*dcosdrij**2 + deijkdcos*d2cosdrij2) + 2.0_dp*ddampij*deijkdcos*dcosdrij + &
                d2dampij*eijknodamp)
        call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,vij(1),vij(2),vij(3), &
                          deijkdrij,d2trm,dr2ds6(1,1),d2r2dx26(1,1),d2r2ds26(1,1,1),d2r2dsdx6(1,1,1),lilocal)
!
!  i-j / i-k
!
        d2trm = dampij*dampik*(d2eijkdcos2*dcosdrij*dcosdrik + deijkdcos*d2cosdrijdrik) + &
                dampij*ddampik*deijkdcos*dcosdrij + dampik*ddampij*deijkdcos*dcosdrik + &
                eijknodamp*ddampij*ddampik
        call add_drv2_2dm(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,vij(1),vij(2),vij(3),vik(1),vik(2),vik(3), &
                          d2trm,dr2ds6(1,1),dr2ds6(1,2),lilocal)
!
!  i-j / j-k
!
        d2trm = dampij*dampik*(d2eijkdcos2*dcosdrij*dcosdrjk + deijkdcos*d2cosdrijdrjk) + &
                dampik*ddampij*deijkdcos*dcosdrjk
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,vij(1),vij(2),vij(3),vjk(1),vjk(2),vjk(3), &
                          d2trm,dr2ds6(1,1),dr2ds6(1,3),lilocal)
!
!  i-k / i-k
!
        d2trm = dampij*(dampik*(d2eijkdcos2*dcosdrik**2 + deijkdcos*d2cosdrik2) + 2.0_dp*ddampik*deijkdcos*dcosdrik + &
                d2dampik*eijknodamp)
        call add_drv2_1dm(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,vik(1),vik(2),vik(3), &
                          deijkdrik,d2trm,dr2ds6(1,2),d2r2dx26(1,2),d2r2ds26(1,1,2),d2r2dsdx6(1,1,2),lilocal)
!
!  i-k / j-k
!
        d2trm = dampij*dampik*(d2eijkdcos2*dcosdrik*dcosdrjk + deijkdcos*d2cosdrikdrjk) + &
                dampij*ddampik*deijkdcos*dcosdrjk
        call add_drv2_2dm(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,vik(1),vik(2),vik(3),vjk(1),vjk(2),vjk(3), &
                          d2trm,dr2ds6(1,2),dr2ds6(1,3),lilocal)
!
!  j-k / j-k
!
        d2trm = dampij*dampik*d2eijkdcos2*dcosdrjk**2
        call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,vjk(1),vjk(2),vjk(3), &
                          deijkdrjk,d2trm,dr2ds6(1,3),d2r2dx26(1,3),d2r2ds26(1,1,3),d2r2dsdx6(1,1,3),lilocal)
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
  do ia = 1,n_gfnff_torsions
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
    jloc = atom2local(j)
    ljlocal = (jloc.gt.0)
!
!  Only do energy and first derivatives if j is local
!
    if (.not.ljlocal.and..not.lgrad2) cycle
!
    njk = n_gfnff_torsionptr(2,ia)
    nij = n_gfnff_torsionptr(3,ia)
    nkl = n_gfnff_torsionptr(4,ia)
!
    i = nbrno_bond(nij,j)
    k = nbrno_bond(njk,j)
    l = nbrno_bond(nkl,k)
!
    iloc = atom2local(i)
    kloc = atom2local(k)
    lloc = atom2local(l)
    lilocal = (iloc.gt.0)
    lklocal = (kloc.gt.0)
    lllocal = (lloc.gt.0)
!
!  Only do if one of the atoms is local
!
    if (.not.lilocal.and..not.ljlocal.and..not.lklocal.and..not.lllocal) cycle
!
    if (lgrad2) then
      indi = 3*(i-1)
      ixf = indi + 1 
      iyf = indi + 2 
      izf = indi + 3 
      if (lilocal) then
        indi = 3*(iloc-1)
        ix = indi + 1 
        iy = indi + 2 
        iz = indi + 3 
      endif
!     
      indj = 3*(j-1)
      jxf = indj + 1
      jyf = indj + 2
      jzf = indj + 3
      if (ljlocal) then
        indj = 3*(jloc-1)
        jx = indj + 1 
        jy = indj + 2 
        jz = indj + 3 
      endif
!   
      indk = 3*(k-1)
      kxf = indk + 1
      kyf = indk + 2
      kzf = indk + 3
      if (lklocal) then
        indk = 3*(kloc-1)
        kx = indk + 1 
        ky = indk + 2 
        kz = indk + 3 
      endif
!   
      indl = 3*(l-1)
      lxf = indl + 1
      lyf = indl + 2
      lzf = indl + 3
      if (lllocal) then
        indl = 3*(lloc-1)
        lx = indl + 1 
        ly = indl + 2 
        lz = indl + 3 
      endif
!
!  Place into arrays to make looping over vectors possible
!
      locals(1,1) = lilocal
      locals(2,1) = ljlocal
!
      indx(1,1) = ix
      indy(1,1) = iy
      indz(1,1) = iz
      indx(2,1) = jx
      indy(2,1) = jy
      indz(2,1) = jz
!
      indxf(1,1) = ixf
      indyf(1,1) = iyf
      indzf(1,1) = izf
      indxf(2,1) = jxf
      indyf(2,1) = jyf
      indzf(2,1) = jzf
!
      locals(1,2) = lilocal
      locals(2,2) = lklocal
!
      indx(1,2) = ix
      indy(1,2) = iy
      indz(1,2) = iz
      indx(2,2) = kx
      indy(2,2) = ky
      indz(2,2) = kz
!
      indxf(1,2) = ixf
      indyf(1,2) = iyf
      indzf(1,2) = izf
      indxf(2,2) = kxf
      indyf(2,2) = kyf
      indzf(2,2) = kzf
!
      locals(1,3) = lilocal
      locals(2,3) = lllocal
!
      indx(1,3) = ix
      indy(1,3) = iy
      indz(1,3) = iz
      indx(2,3) = lx
      indy(2,3) = ly
      indz(2,3) = lz
!
      indxf(1,3) = ixf
      indyf(1,3) = iyf
      indzf(1,3) = izf
      indxf(2,3) = lxf
      indyf(2,3) = lyf
      indzf(2,3) = lzf
!
      locals(1,4) = ljlocal
      locals(2,4) = lklocal
!
      indx(1,4) = jx
      indy(1,4) = jy
      indz(1,4) = jz
      indx(2,4) = kx
      indy(2,4) = ky
      indz(2,4) = kz
!
      indxf(1,4) = jxf
      indyf(1,4) = jyf
      indzf(1,4) = jzf
      indxf(2,4) = kxf
      indyf(2,4) = kyf
      indzf(2,4) = kzf
!
      locals(1,5) = ljlocal
      locals(2,5) = lllocal
!
      indx(1,5) = jx
      indy(1,5) = jy
      indz(1,5) = jz
      indx(2,5) = lx
      indy(2,5) = ly
      indz(2,5) = lz
!
      indxf(1,5) = jxf
      indyf(1,5) = jyf
      indzf(1,5) = jzf
      indxf(2,5) = lxf
      indyf(2,5) = lyf
      indzf(2,5) = lzf
!
      locals(1,6) = lklocal
      locals(2,6) = lllocal
!
      indx(1,6) = kx
      indy(1,6) = ky
      indz(1,6) = kz
      indx(2,6) = lx
      indy(2,6) = ly
      indz(2,6) = lz
!
      indxf(1,6) = kxf
      indyf(1,6) = kyf
      indzf(1,6) = kzf
      indxf(2,6) = lxf
      indyf(2,6) = lyf
      indzf(2,6) = lzf
    endif
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
    call gfnff_damp(rik2,rcut2ik,dampik,ddampik,d2dampik,lgrad1,lgrad2)
    call gfnff_damp(rjk2,rcut2jk,dampjk,ddampjk,d2dampjk,lgrad1,lgrad2)
    call gfnff_damp(rjl2,rcut2jl,dampjl,ddampjl,d2dampjl,lgrad1,lgrad2)
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
    call gfnff_torsion(1_i4,par_gfnff_torsion(2,ia),rn,phi0,rij,rik,ril,rjk,rjl,rkl,eijklnod,e1d,e2d,lgrad1,lgrad2)
    eijkl = eijklnod*damp
    if (ljlocal) then
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
    endif
!----------------------
!  First derivatives  |
!----------------------
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
!----------------------------
!  Add damping derivatives  |
!----------------------------
!
!  Do second derivatives first so that e1d is not yet modified
!
      if (lgrad2) then
!
!  Scale current derivatives by damping
!
        e2d(1:21) = e2d(1:21)*damp
!
!  Indexing in e2d
!
!  1 = 11  7 = 22 12 = 33 16 = 44 19 = 55 21 = 66
!  2 = 21  8 = 32 13 = 43 17 = 54 20 = 65
!  3 = 31  9 = 42 14 = 53 18 = 64
!  4 = 41 10 = 52 15 = 63
!  5 = 51 11 = 62
!  6 = 61
!
        e2d(2)  = e2d(2)  + e1d(1)*dampjl*dampjk*ddampik
        e2d(7)  = e2d(7)  + e1d(2)*dampjl*dampjk*ddampik*2.0_dp
        e2d(8)  = e2d(8)  + e1d(3)*dampjl*dampjk*ddampik
        e2d(9)  = e2d(9)  + e1d(4)*dampjl*dampjk*ddampik
        e2d(10) = e2d(10) + e1d(5)*dampjl*dampjk*ddampik
        e2d(11) = e2d(11) + e1d(6)*dampjl*dampjk*ddampik
!
        e2d(4)  = e2d(4)  + e1d(1)*dampjl*ddampjk*dampik
        e2d(9)  = e2d(9)  + e1d(2)*dampjl*ddampjk*dampik
        e2d(13) = e2d(13) + e1d(3)*dampjl*ddampjk*dampik
        e2d(16) = e2d(16) + e1d(4)*dampjl*ddampjk*dampik*2.0_dp
        e2d(17) = e2d(17) + e1d(5)*dampjl*ddampjk*dampik
        e2d(18) = e2d(18) + e1d(6)*dampjl*ddampjk*dampik
!
        e2d(5)  = e2d(5)  + e1d(1)*ddampjl*dampjk*dampik
        e2d(10) = e2d(10) + e1d(2)*ddampjl*dampjk*dampik
        e2d(14) = e2d(14) + e1d(3)*ddampjl*dampjk*dampik
        e2d(17) = e2d(17) + e1d(4)*ddampjl*dampjk*dampik
        e2d(19) = e2d(19) + e1d(5)*ddampjl*dampjk*dampik*2.0_dp
        e2d(20) = e2d(20) + e1d(6)*ddampjl*dampjk*dampik
!
        e2d(7)  = e2d(7)  + eijklnod*dampjl*dampjk*d2dampik
        e2d(16) = e2d(16) + eijklnod*dampjl*d2dampjk*dampik
        e2d(19) = e2d(19) + eijklnod*d2dampjl*dampjk*dampik
!
        e2d(9)  = e2d(9)  + eijklnod*dampjl*ddampjk*ddampik
        e2d(10) = e2d(10) + eijklnod*ddampjl*dampjk*ddampik
        e2d(17) = e2d(17) + eijklnod*ddampjl*ddampjk*dampik
      endif
!
!  Scale current derivatives by damping
!
      e1d(1:6) = e1d(1:6)*damp
!
!  First derivatives of damping
!
      e1d(2) = e1d(2) + eijklnod*dampjl*dampjk*ddampik
      e1d(4) = e1d(4) + eijklnod*dampjl*ddampjk*dampik
      e1d(5) = e1d(5) + eijklnod*ddampjl*dampjk*dampik
!
      if (ljlocal) then
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
      endif
!
      if (lstr.or.lgrad2) then
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
        call realstrterms(ndim,6_i4,6_i4,xv4,yv4,zv4,xcom,ycom,zcom,dr2ds6,d2r2dx26,d2r2dsdx6,d2r2ds26,lgrad2)
      endif
!
      if (ljlocal.and.lstr) then
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
!-----------------------
!  Second derivatives  |
!-----------------------
      if (lgrad2) then
        ind = 0
        do ii = 1,6
          do jj = ii,6
            ind = ind + 1
            if (ii.eq.jj) then
              call add_drv2_1dm(locals(1,ii),locals(2,ii),indx(1,ii),indy(1,ii),indz(1,ii),indxf(1,ii),indyf(1,ii),indzf(1,ii), &
                                indx(2,ii),indy(2,ii),indz(2,ii),indxf(2,ii),indyf(2,ii),indzf(2,ii),xv4(ii),yv4(ii),zv4(ii), &
                                e1d(ii),e2d(ind),dr2ds6(1,ii),d2r2dx26(1,ii),d2r2ds26(1,1,ii),d2r2dsdx6(1,1,ii),ljlocal)
            else
              call add_drv2_2dm(locals(1,ii),locals(2,ii),locals(1,jj),locals(2,jj), &
                                indx(1,ii),indy(1,ii),indz(1,ii),indxf(1,ii),indyf(1,ii),indzf(1,ii), &
                                indx(2,ii),indy(2,ii),indz(2,ii),indxf(2,ii),indyf(2,ii),indzf(2,ii), &
                                indx(1,jj),indy(1,jj),indz(1,jj),indxf(1,jj),indyf(1,jj),indzf(1,jj), &
                                indx(2,jj),indy(2,jj),indz(2,jj),indxf(2,jj),indyf(2,jj),indzf(2,jj), &
                                xv4(ii),yv4(ii),zv4(ii),xv4(jj),yv4(jj),zv4(jj),e2d(ind), &
                                dr2ds6(1,ii),dr2ds6(1,jj),ljlocal)
            endif
          enddo
        enddo
      endif
    endif
  enddo
!**********************
!  Improper torsions  *
!**********************
  do ia = 1,n_gfnff_torsions
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
    iloc = atom2local(i)
    lilocal = (iloc.ne.0)
!
!  Only do energy and first derivatives if i is local
!
    if (.not.lilocal.and..not.lgrad2) cycle
!
    nij = n_gfnff_torsionptr(2,ia)
    nik = n_gfnff_torsionptr(3,ia)
    nil = n_gfnff_torsionptr(4,ia)
!
    j = nbrno_bond(nij,i)
    k = nbrno_bond(nik,i)
    l = nbrno_bond(nil,i)
!
    jloc = atom2local(j)
    kloc = atom2local(k)
    lloc = atom2local(l)
!
    ljlocal = (jloc.ne.0)
    lklocal = (kloc.ne.0)
    lllocal = (lloc.ne.0)
!
!  Are any of the atoms local to this node?
!
    if (.not.lilocal.and..not.ljlocal.and..not.lklocal.and..not.lllocal) cycle
!
    if (lgrad2) then
      indi = 3*(i-1)
      ixf = indi + 1
      iyf = indi + 2
      izf = indi + 3
      if (lilocal) then
        indi = 3*(iloc-1)
        ix = indi + 1
        iy = indi + 2
        iz = indi + 3
      endif
!
      indj = 3*(j-1)
      jxf = indj + 1
      jyf = indj + 2
      jzf = indj + 3
      if (ljlocal) then
        indj = 3*(jloc-1)
        jx = indj + 1
        jy = indj + 2
        jz = indj + 3
      endif
!
      indk = 3*(k-1)
      kxf = indk + 1
      kyf = indk + 2
      kzf = indk + 3
      if (lklocal) then
        indk = 3*(kloc-1)
        kx = indk + 1
        ky = indk + 2
        kz = indk + 3
      endif
!
      indl = 3*(l-1)
      lxf = indl + 1
      lyf = indl + 2
      lzf = indl + 3
      if (lllocal) then
        indl = 3*(lloc-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
      endif
    endif
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
    call gfnff_damp(rij2,rcut2ij,dampij,ddampij,d2dampij,lgrad1,lgrad2)
    call gfnff_damp(rjk2,rcut2jk,dampjk,ddampjk,d2dampjk,lgrad1,lgrad2)
    call gfnff_damp(rjl2,rcut2jl,dampjl,ddampjl,d2dampjl,lgrad1,lgrad2)
    damp = dampij*dampjk*dampjl
!
    call gfnff_inversion(ndim,numat,nnbr_bond,maxnbr,rbnbr,xbnbr,ybnbr,zbnbr,i,nij,nik,nil,phi,phi1dx, &
                         phi1ds,phi2dx,phi2ds,phi2dsdx,lgrad1,lgrad2)
!
    if (nrot.eq.0) then
      dphi = phi - phi0
      fphi = dphi + pi
      cosphi = cos(fphi)
      eijklnod = (1.0_dp + cosphi)*par_gfnff_torsion(2,ia)
      if (lgrad1) then
        sinphi = sin(fphi)
        deijkldphi = - par_gfnff_torsion(2,ia)*sinphi
        if (lgrad2) then
          d2eijkldphi2 = - par_gfnff_torsion(2,ia)*cosphi
        endif
      endif
    else
      cosphi  = cos(phi)
      cosphi0 = cos(phi0)
      eijklnod = par_gfnff_torsion(2,ia)*(cosphi - cosphi0)**2
      if (lgrad1) then
        sinphi  = sin(phi)
        deijkldphi = - 2.0_dp*par_gfnff_torsion(2,ia)*sinphi*(cosphi - cosphi0)
        if (lgrad2) then
          d2eijkldphi2 = - 2.0_dp*par_gfnff_torsion(2,ia)*(cosphi*(cosphi - cosphi0) - sinphi*sinphi)
        endif
      endif
    endif
!
    eijkl = eijklnod*damp
    if (lilocal) then
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
    endif
!**********************
!  First derivatives  *
!**********************
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
      d1trm = damp*deijkldphi
      if (lilocal) then
        xdrv(i) = xdrv(i) - d1trm*phi1dx(1,1)
        ydrv(i) = ydrv(i) - d1trm*phi1dx(2,1)
        zdrv(i) = zdrv(i) - d1trm*phi1dx(3,1)
!
        xdrv(j) = xdrv(j) + d1trm*phi1dx(1,1)
        ydrv(j) = ydrv(j) + d1trm*phi1dx(2,1)
        zdrv(j) = zdrv(j) + d1trm*phi1dx(3,1)
!
        xdrv(j) = xdrv(j) - d1trm*phi1dx(1,2)
        ydrv(j) = ydrv(j) - d1trm*phi1dx(2,2)
        zdrv(j) = zdrv(j) - d1trm*phi1dx(3,2)
!
        xdrv(k) = xdrv(k) + d1trm*phi1dx(1,2)
        ydrv(k) = ydrv(k) + d1trm*phi1dx(2,2)
        zdrv(k) = zdrv(k) + d1trm*phi1dx(3,2)
!
        xdrv(i) = xdrv(i) - d1trm*phi1dx(1,3)
        ydrv(i) = ydrv(i) - d1trm*phi1dx(2,3)
        zdrv(i) = zdrv(i) - d1trm*phi1dx(3,3)
!
        xdrv(l) = xdrv(l) + d1trm*phi1dx(1,3)
        ydrv(l) = ydrv(l) + d1trm*phi1dx(2,3)
        zdrv(l) = zdrv(l) + d1trm*phi1dx(3,3)
!
        if (lstr) then
          do ks = 1,nstrains
            rstrd(ks) = rstrd(ks) + d1trm*phi1ds(ks)
          enddo
        endif
      endif
!
!  Add damping derivatives
!
!  i-j
!
      deijkldrij = eijklnod*ddampij*dampjk*dampjl
      if (lilocal) then
        xdrv(i) = xdrv(i) - deijkldrij*vij(1)
        ydrv(i) = ydrv(i) - deijkldrij*vij(2)
        zdrv(i) = zdrv(i) - deijkldrij*vij(3)
!
        xdrv(j) = xdrv(j) + deijkldrij*vij(1)
        ydrv(j) = ydrv(j) + deijkldrij*vij(2)
        zdrv(j) = zdrv(j) + deijkldrij*vij(3)
      endif
!
      if (lstr.or.lgrad2) then
        call real1strterm(ndim,vij(1),vij(2),vij(3),0.0_dp,0.0_dp,0.0_dp,dr2ds6(1,1),d2r2dx26(1,1), &
                          d2r2dsdx6(1,1,1),d2r2ds26(1,1,1),lgrad2)
      endif
      if (lilocal) then
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkldrij*dr2ds6(ks,1)
          enddo
        endif
      endif
!
!  j-k
!
      deijkldrjk = eijklnod*dampij*ddampjk*dampjl
      if (lilocal) then
        xdrv(j) = xdrv(j) - deijkldrjk*vjk(1)
        ydrv(j) = ydrv(j) - deijkldrjk*vjk(2)
        zdrv(j) = zdrv(j) - deijkldrjk*vjk(3)
!
        xdrv(k) = xdrv(k) + deijkldrjk*vjk(1)
        ydrv(k) = ydrv(k) + deijkldrjk*vjk(2)
        zdrv(k) = zdrv(k) + deijkldrjk*vjk(3)
      endif
!
      if (lstr.or.lgrad2) then
        call real1strterm(ndim,vjk(1),vjk(2),vjk(3),0.0_dp,0.0_dp,0.0_dp,dr2ds6(1,2),d2r2dx26(1,2), &
                          d2r2dsdx6(1,1,2),d2r2ds26(1,1,2),lgrad2)
      endif
      if (lilocal) then
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkldrjk*dr2ds6(ks,2)
          enddo
        endif
      endif
!
!  j-l
!
      deijkldrjl = eijklnod*dampij*dampjk*ddampjl
      if (lilocal) then
        xdrv(j) = xdrv(j) - deijkldrjl*vjl(1)
        ydrv(j) = ydrv(j) - deijkldrjl*vjl(2)
        zdrv(j) = zdrv(j) - deijkldrjl*vjl(3)
!
        xdrv(l) = xdrv(l) + deijkldrjl*vjl(1)
        ydrv(l) = ydrv(l) + deijkldrjl*vjl(2)
        zdrv(l) = zdrv(l) + deijkldrjl*vjl(3)
      endif
!
      if (lstr.or.lgrad2) then
        call real1strterm(ndim,vjl(1),vjl(2),vjl(3),0.0_dp,0.0_dp,0.0_dp,dr2ds6(1,3),d2r2dx26(1,3), &
                          d2r2dsdx6(1,1,3),d2r2ds26(1,1,3),lgrad2)
      endif
      if (lilocal) then
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkldrjl*dr2ds6(ks,3)
          enddo
        endif
      endif
      if (lgrad2) then
!***********************
!  Second derivatives  *
!***********************
!----------------------------
!  Derivatives of phi only  |
!----------------------------
        d2trm = damp*d2eijkldphi2
!
!  i-j / i-j
!
        call add_drv2_1bdm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,d1trm,d2trm, &
                           phi1dx(1,1),phi1ds,phi2dx(1,1,1),phi2dsdx(1,1,1))
!
!  i-j / j-k
!
        call add_drv2_2bdm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,jx,jy,jz,jxf,jyf,jzf, &
                           kx,ky,kz,kxf,kyf,kzf,d1trm,d2trm,phi1dx(1,1),phi1dx(1,2),phi2dx(1,1,2))
!
!  i-j / i-l
!
        call add_drv2_2bdm(lilocal,ljlocal,lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf, &
                           lx,ly,lz,lxf,lyf,lzf,d1trm,d2trm,phi1dx(1,1),phi1dx(1,3),phi2dx(1,1,3))
!
!  j-k / j-k
!
        call add_drv2_1bdm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,d1trm,d2trm, &
                           phi1dx(1,2),phi1ds,phi2dx(1,1,4),phi2dsdx(1,1,2))
!
!  j-k / i-l
!
        call add_drv2_2bdm(ljlocal,lklocal,lilocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,ix,iy,iz,ixf,iyf,izf, &
                           lx,ly,lz,lxf,lyf,lzf,d1trm,d2trm,phi1dx(1,2),phi1dx(1,3),phi2dx(1,1,5))
!
!  i-l / i-l
!
        call add_drv2_1bdm(lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,d1trm,d2trm, &
                           phi1dx(1,3),phi1ds,phi2dx(1,1,6),phi2dsdx(1,1,3))
!
        if (lilocal) then
          if (lstr) then
!
!  Strain-strain
!
            do ks = 1,nstrains
              do kt = 1,nstrains
                sderv2(kt,ks) = sderv2(kt,ks) + d1trm*phi2ds(kt,ks) + d2trm*phi1ds(kt)*phi1ds(ks)
              enddo
            enddo
          endif
        endif
!-------------------------------------------------
!  Second derivatives of damping functions only  |
!-------------------------------------------------
!
!  i-j / i-j
!
        d2trm = eijklnod*d2dampij*dampjk*dampjl
        call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,vij(1),vij(2),vij(3),deijkldrij, &
                          d2trm,dr2ds6(1,1),d2r2dx26(1,1),d2r2ds26(1,1,1),d2r2dsdx6(1,1,1),lilocal)
!
!  i-j / j-k
!
        d2trm = eijklnod*ddampij*ddampjk*dampjl
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,vij(1),vij(2),vij(3),vjk(1),vjk(2),vjk(3),d2trm, &
                          dr2ds6(1,1),dr2ds6(1,2),lilocal)
!
!  i-j / j-l
!
        d2trm = eijklnod*ddampij*dampjk*ddampjl
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,vij(1),vij(2),vij(3),vjl(1),vjl(2),vjl(3),d2trm, &
                          dr2ds6(1,1),dr2ds6(1,3),lilocal)
!
!  j-k / j-k
!
        d2trm = eijklnod*dampij*d2dampjk*dampjl
        call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,vjk(1),vjk(2),vjk(3),deijkldrjk, &
                          d2trm,dr2ds6(1,2),d2r2dx26(1,2),d2r2ds26(1,1,2),d2r2dsdx6(1,1,2),lilocal)
!
!  j-k / j-l
!
        d2trm = eijklnod*dampij*ddampjk*ddampjl
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,vjk(1),vjk(2),vjk(3),vjl(1),vjl(2),vjl(3),d2trm, &
                          dr2ds6(1,2),dr2ds6(1,3),lilocal)
!
!  j-l / j-l
!
        d2trm = eijklnod*dampij*dampjk*d2dampjl
        call add_drv2_1dm(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,vjl(1),vjl(2),vjl(3),deijkldrjl, &
                          d2trm,dr2ds6(1,3),d2r2dx26(1,3),d2r2ds26(1,1,3),d2r2dsdx6(1,1,3),lilocal)
!---------------------------------------------------------------
!  Second derivatives with a mix of phi and damping functions  |
!---------------------------------------------------------------
!
!  Create copy of strain derivatives for phi divided by 3 to correct for being added in 3 different cases
!
        dtmps2 = phi1ds/3.0_dp
!
!  phi / damp
!
!  i-j / i-j
!
        d2trm = deijkldphi*ddampij*dampjk*dampjl
        do ks = 1,nstrains
          kl = nstrptr(ks)
          dtmps(ks) = dr2ds6(kl,1)
        enddo
        call add_drv2_1bxdm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,d2trm, &
                            phi1dx(1,1),vij,dtmps2,dtmps,lilocal)
!
!  j-k / i-j
!
        call add_drv2_2bxdm(ljlocal,lklocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,d2trm,phi1dx(1,2),vij,dtmps2,dtmps,lilocal)
!
!  i-l / i-j
!
        call add_drv2_2bxdm(lilocal,lllocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,d2trm,phi1dx(1,3),vij,dtmps2,dtmps,lilocal)
!
!  i-j / j-k
!
        d2trm = deijkldphi*dampij*ddampjk*dampjl
        do ks = 1,nstrains
          kl = nstrptr(ks)
          dtmps(ks) = dr2ds6(kl,2)
        enddo
        call add_drv2_2bxdm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,d2trm,phi1dx(1,1),vjk,dtmps2,dtmps,lilocal)
!
!  j-k / j-k
!
        call add_drv2_1bxdm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,d2trm, &
                            phi1dx(1,2),vjk,dtmps2,dtmps,lilocal)
!
!  i-l / j-k
!
        call add_drv2_2bxdm(lilocal,lllocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,d2trm,phi1dx(1,3),vjk,dtmps2,dtmps,lilocal)
!
!  i-j / j-l
!
        d2trm = deijkldphi*dampij*dampjk*ddampjl
        do ks = 1,nstrains
          kl = nstrptr(ks)
          dtmps(ks) = dr2ds6(kl,3)
        enddo
        call add_drv2_2bxdm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d2trm,phi1dx(1,1),vjl,dtmps2,dtmps,lilocal)
!
!  j-k / j-l
!
        call add_drv2_2bxdm(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d2trm,phi1dx(1,2),vjl,dtmps2,dtmps,lilocal)
!
!  i-l / j-l
!
        call add_drv2_2bxdm(lilocal,lllocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d2trm,phi1dx(1,3),vjl,dtmps2,dtmps,lilocal)
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
  do ia = 1,nb3atm
    i = nb3list(1,ia)
    iloc = atom2local(i)
    lilocal = (iloc.gt.0)
!
!  Only do energy and first derivatives if i is local
!
    if (.not.lilocal.and..not.lgrad2) cycle
!
    ni = nb3list(2,ia)
    nj = nb3list(3,ia)
    nk = nb3list(4,ia)
!
    j = nbrno_bond(ni,i)
    k = nbrno_bond(nj,j)
    l = nbrno_bond(nk,k)
!
    lloc = atom2local(l)
    lllocal = (lloc.gt.0)
!
    nati = nat(i)
    natl = nat(l)
!
    if (lgrad2) then
      indi = 3*(i-1)
      ixf = indi + 1
      iyf = indi + 2
      izf = indi + 3
      if (lilocal) then
        indi = 3*(iloc-1)
        ix = indi + 1
        iy = indi + 2
        iz = indi + 3
      endif
!
      indl = 3*(l-1)
      lxf = indl + 1
      lyf = indl + 2
      lzf = indl + 3
      if (lllocal) then
        indl = 3*(lloc-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
      endif
    endif
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
    ril03 = ril0**1.5
!------------------------------
!  Loop over neighbours of i  !
!------------------------------
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
      jloc = atom2local(j)
      ljlocal = (jloc.gt.0)
!
!  Only do if one of the atoms for derivatives is local
!
      if (.not.lilocal.and..not.ljlocal.and..not.lllocal) cycle
!
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
      if (lgrad2) then
        indj = 3*(j-1)
        jxf = indj + 1
        jyf = indj + 2
        jzf = indj + 3
        if (ljlocal) then
          indj = 3*(jloc-1)
          jx = indj + 1
          jy = indj + 2
          jz = indj + 3
        endif
      endif
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
      if (lilocal) then
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
      endif
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
        if (lilocal) then
          xdrv(i) = xdrv(i) - deijkdrij*xbnbr(ni,i)
          ydrv(i) = ydrv(i) - deijkdrij*ybnbr(ni,i)
          zdrv(i) = zdrv(i) - deijkdrij*zbnbr(ni,i)
!
          xdrv(j) = xdrv(j) + deijkdrij*xbnbr(ni,i)
          ydrv(j) = ydrv(j) + deijkdrij*ybnbr(ni,i)
          zdrv(j) = zdrv(j) + deijkdrij*zbnbr(ni,i)
        endif
!
        if (lstr.or.lgrad2) then
          call real1strterm(ndim,xbnbr(ni,i),ybnbr(ni,i),zbnbr(ni,i),xcom(1),ycom(1),zcom(1),dr2ds6(1,1),d2r2dx26(1,1), &
                            d2r2dsdx6(1,1,1),d2r2ds26(1,1,1),lgrad2)
        endif
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds6(ks,1)
            enddo
          endif
        endif
!
!  i-l
!
        if (lilocal) then
          xdrv(i) = xdrv(i) - deijkdril*vil(1)
          ydrv(i) = ydrv(i) - deijkdril*vil(2)
          zdrv(i) = zdrv(i) - deijkdril*vil(3)
!
          xdrv(l) = xdrv(l) + deijkdril*vil(1)
          ydrv(l) = ydrv(l) + deijkdril*vil(2)
          zdrv(l) = zdrv(l) + deijkdril*vil(3)
        endif
!     
        if (lstr.or.lgrad2) then
          call real1strterm(ndim,vil(1),vil(2),vil(3),xcom(1),ycom(1),zcom(1),dr2ds6(1,2),d2r2dx26(1,2), &
                            d2r2dsdx6(1,1,2),d2r2ds26(1,1,2),lgrad2)
        endif
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijkdril*dr2ds6(ks,2)
            enddo
          endif
        endif
!
!  j-l
!
        if (lilocal) then
          xdrv(j) = xdrv(j) - deijkdrjl*vjl(1)
          ydrv(j) = ydrv(j) - deijkdrjl*vjl(2)
          zdrv(j) = zdrv(j) - deijkdrjl*vjl(3)
!
          xdrv(l) = xdrv(l) + deijkdrjl*vjl(1)
          ydrv(l) = ydrv(l) + deijkdrjl*vjl(2)
          zdrv(l) = zdrv(l) + deijkdrjl*vjl(3)
        endif
!     
        if (lstr.or.lgrad2) then
          call real1strterm(ndim,vjl(1),vjl(2),vjl(3),xcom(1),ycom(1),zcom(1),dr2ds6(1,3),d2r2dx26(1,3), &
                            d2r2dsdx6(1,1,3),d2r2ds26(1,1,3),lgrad2)
        endif
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijkdrjl*dr2ds6(ks,3)
            enddo
          endif
        endif
        if (lgrad2) then
!-----------------------
!  Second derivatives  |
!-----------------------
!
!  Order of distances:
!  d1 : 1 = ij / 2 = il / 3 = jl
!  d2 : 1 = 11 / 2 = 12 / 3 = 13 / 4 = 22 / 5 = 23 / 6 = 33
!
          d2cos1(1) = - rrij*rrij*(dcos1drij + rrij*rrjl - 2.0_dp*cos1*rrij**2)
          d2cos1(2) = rrjl*rrij**3
          d2cos1(3) = - rrij*rrjl**3 - rrij*rrij*dcos1drjl
          d2cos1(4) = 0.0_dp
          d2cos1(5) = rrij*rrjl**3
          d2cos1(6) = - rrjl*rrjl*(dcos1drjl + rrij*rrjl - 2.0_dp*cos1*rrjl**2)
!
          d2cos2(1) = - rrij*rrij*(dcos2drij + rrij*rril - 2.0_dp*cos2*rrij**2)
          d2cos2(2) = - rrij*rril**3 - rrij*rrij*dcos2dril
          d2cos2(3) = rril*rrij**3
          d2cos2(4) = - rril*rril*(dcos2dril + rrij*rril - 2.0_dp*cos2*rril**2)
          d2cos2(5) = rrij*rril**3
          d2cos2(6) = 0.0_dp
!
          d2cos3(1) = 0.0_dp
          d2cos3(2) = rrjl*rril**3
          d2cos3(3) = rril*rrjl**3
          d2cos3(4) = - rril*rril*(dcos3dril + rril*rrjl - 2.0_dp*cos3*rril**2)
          d2cos3(5) = - rril*rrjl**3 - rril*rril*dcos3drjl
          d2cos3(6) = - rrjl*rrjl*(dcos3drjl + rril*rrjl - 2.0_dp*cos3*rrjl**2)
!
!  i-j / i-j
!
          if (lilocal.or.ljlocal) then
            d2trm = eijk*rrij3*(18.0_dp*rrij3*rij2 - 3.0_dp*rrij) &
                    - (18.0_dp*c9*rrijk9*rrij3*rij)*(dcos1drij*cos2*cos3 + cos1*dcos2drij*cos3 + dcos3drij*cos1*cos2) &
                    + (3.0_dp*c9*rrijk9)*(d2cos1(1)*cos2*cos3 + cos1*d2cos2(1)*cos3 + cos1*cos2*d2cos3(1) &
                                         + dcos1drij*(dcos2drij*cos3 + cos2*dcos3drij) &
                                         + dcos2drij*(dcos1drij*cos3 + cos1*dcos3drij) &
                                         + dcos3drij*(dcos1drij*cos2 + cos1*dcos2drij))
!
            call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xbnbr(ni,i),ybnbr(ni,i),zbnbr(ni,i), &
                              deijkdrij,d2trm,dr2ds6(1,1),d2r2dx26(1,1),d2r2ds26(1,1,1),d2r2dsdx6(1,1,1),lilocal)
          endif
!
!  i-j / i-l
!
          d2trm = 9.0_dp*eijk*rij*ril*rrij3*rril3 &
                  - (9.0_dp*c9*rrijk9*rrij3*rij)*(dcos1dril*cos2*cos3 + cos1*dcos2dril*cos3 + dcos3dril*cos1*cos2) &
                  - (9.0_dp*c9*rrijk9*rril3*ril)*(dcos1drij*cos2*cos3 + cos1*dcos2drij*cos3 + dcos3drij*cos1*cos2) &
                  + (3.0_dp*c9*rrijk9)*(d2cos1(2)*cos2*cos3 + cos1*d2cos2(2)*cos3 + cos1*cos2*d2cos3(2) &
                                       + dcos1drij*(dcos2dril*cos3 + cos2*dcos3dril) &
                                       + dcos2drij*(dcos1dril*cos3 + cos1*dcos3dril) &
                                       + dcos3drij*(dcos1dril*cos2 + cos1*dcos2dril))
!
          call add_drv2_2dm(lilocal,ljlocal,lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf, &
                            lx,ly,lz,lxf,lyf,lzf,xbnbr(ni,i),ybnbr(ni,i),zbnbr(ni,i),vil(1),vil(2),vil(3),d2trm, &
                            dr2ds6(1,1),dr2ds6(1,2),lilocal)
!
!  i-j / j-l
!
          d2trm = 9.0_dp*eijk*rij*rjl*rrij3*rrjl3 &
                  - (9.0_dp*c9*rrijk9*rrij3*rij)*(dcos1drjl*cos2*cos3 + cos1*dcos2drjl*cos3 + dcos3drjl*cos1*cos2) &
                  - (9.0_dp*c9*rrijk9*rrjl3*rjl)*(dcos1drij*cos2*cos3 + cos1*dcos2drij*cos3 + dcos3drij*cos1*cos2) &
                  + (3.0_dp*c9*rrijk9)*(d2cos1(3)*cos2*cos3 + cos1*d2cos2(3)*cos3 + cos1*cos2*d2cos3(3) &
                                       + dcos1drij*(dcos2drjl*cos3 + cos2*dcos3drjl) &
                                       + dcos2drij*(dcos1drjl*cos3 + cos1*dcos3drjl) &
                                       + dcos3drij*(dcos1drjl*cos2 + cos1*dcos2drjl))
!
          call add_drv2_2dm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,jx,jy,jz,jxf,jyf,jzf, &
                            lx,ly,lz,lxf,lyf,lzf,xbnbr(ni,i),ybnbr(ni,i),zbnbr(ni,i),vjl(1),vjl(2),vjl(3),d2trm, &
                            dr2ds6(1,1),dr2ds6(1,3),lilocal)
!
!  i-l / i-l
!
          if (lilocal.or.lllocal) then
            d2trm = eijk*rril3*(18.0_dp*rril3*ril2 - 3.0_dp*rril) &
                    - (18.0_dp*c9*rrijk9*rril3*ril)*(dcos1dril*cos2*cos3 + cos1*dcos2dril*cos3 + dcos3dril*cos1*cos2) &
                    + (3.0_dp*c9*rrijk9)*(d2cos1(4)*cos2*cos3 + cos1*d2cos2(4)*cos3 + cos1*cos2*d2cos3(4) &
                                         + dcos1dril*(dcos2dril*cos3 + cos2*dcos3dril) &
                                         + dcos2dril*(dcos1dril*cos3 + cos1*dcos3dril) &
                                         + dcos3dril*(dcos1dril*cos2 + cos1*dcos2dril))
!
            call add_drv2_1dm(lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,vil(1),vil(2),vil(3), &
                              deijkdril,d2trm,dr2ds6(1,2),d2r2dx26(1,2),d2r2ds26(1,1,2),d2r2dsdx6(1,1,2),lilocal)
          endif
!
!  i-l / j-l
!
          d2trm = 9.0_dp*eijk*ril*rjl*rril3*rrjl3 &
                  - (9.0_dp*c9*rrijk9*rril3*ril)*(dcos1drjl*cos2*cos3 + cos1*dcos2drjl*cos3 + dcos3drjl*cos1*cos2) &
                  - (9.0_dp*c9*rrijk9*rrjl3*rjl)*(dcos1dril*cos2*cos3 + cos1*dcos2dril*cos3 + dcos3dril*cos1*cos2) &
                  + (3.0_dp*c9*rrijk9)*(d2cos1(5)*cos2*cos3 + cos1*d2cos2(5)*cos3 + cos1*cos2*d2cos3(5) &
                                       + dcos1dril*(dcos2drjl*cos3 + cos2*dcos3drjl) &
                                       + dcos2dril*(dcos1drjl*cos3 + cos1*dcos3drjl) &
                                       + dcos3dril*(dcos1drjl*cos2 + cos1*dcos2drjl))
!
          call add_drv2_2dm(lilocal,lllocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,jx,jy,jz,jxf,jyf,jzf, &
                            lx,ly,lz,lxf,lyf,lzf,vil(1),vil(2),vil(3),vjl(1),vjl(2),vjl(3),d2trm, &
                            dr2ds6(1,2),dr2ds6(1,3),lilocal)
!
!  j-l / j-l
!
          if (ljlocal.or.lllocal) then
            d2trm = eijk*rrjl3*(18.0_dp*rrjl3*rjl2 - 3.0_dp*rrjl) &
                    - (18.0_dp*c9*rrijk9*rrjl3*rjl)*(dcos1drjl*cos2*cos3 + cos1*dcos2drjl*cos3 + dcos3drjl*cos1*cos2) &
                    + (3.0_dp*c9*rrijk9)*(d2cos1(6)*cos2*cos3 + cos1*d2cos2(6)*cos3 + cos1*cos2*d2cos3(6) &
                                         + dcos1drjl*(dcos2drjl*cos3 + cos2*dcos3drjl) &
                                         + dcos2drjl*(dcos1drjl*cos3 + cos1*dcos3drjl) &
                                         + dcos3drjl*(dcos1drjl*cos2 + cos1*dcos2drjl))
!
            call add_drv2_1dm(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,vjl(1),vjl(2),vjl(3), &
                              deijkdrjl,d2trm,dr2ds6(1,3),d2r2dx26(1,3),d2r2ds26(1,1,3),d2r2dsdx6(1,1,3),lilocal)
          endif
        endif
      endif
    enddo
!------------------------------
!  Loop over neighbours of l  !
!------------------------------
    do nl = 1,nnbr_bond(l)
      j = nbrno_bond(nl,l)
      jloc = atom2local(j)
      ljlocal = (jloc.gt.0)
!
!  Only do if one of the atoms for derivatives is local
!
      if (.not.lilocal.and..not.ljlocal.and..not.lllocal) cycle
!
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
      if (lgrad2) then
        indj = 3*(j-1)
        jxf = indj + 1
        jyf = indj + 2
        jzf = indj + 3
        if (ljlocal) then
          indj = 3*(jloc-1)
          jx = indj + 1
          jy = indj + 2
          jz = indj + 3
        endif
      endif
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
      if (lilocal) then
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
      endif
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
        if (lilocal) then
          xdrv(i) = xdrv(i) - deijkdrij*vij(1)
          ydrv(i) = ydrv(i) - deijkdrij*vij(2)
          zdrv(i) = zdrv(i) - deijkdrij*vij(3)
!
          xdrv(j) = xdrv(j) + deijkdrij*vij(1)
          ydrv(j) = ydrv(j) + deijkdrij*vij(2)
          zdrv(j) = zdrv(j) + deijkdrij*vij(3)
        endif
!
        if (lstr.or.lgrad2) then
          call real1strterm(ndim,vij(1),vij(2),vij(3),xcom(1),ycom(1),zcom(1),dr2ds6(1,1),d2r2dx26(1,1), &
                            d2r2dsdx6(1,1,1),d2r2ds26(1,1,1),lgrad2)
        endif
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds6(ks,1)
            enddo
          endif
        endif
!
!  i-l
!
        if (lilocal) then
          xdrv(i) = xdrv(i) - deijkdril*vil(1)
          ydrv(i) = ydrv(i) - deijkdril*vil(2)
          zdrv(i) = zdrv(i) - deijkdril*vil(3)
!
          xdrv(l) = xdrv(l) + deijkdril*vil(1)
          ydrv(l) = ydrv(l) + deijkdril*vil(2)
          zdrv(l) = zdrv(l) + deijkdril*vil(3)
        endif
!
        if (lstr.or.lgrad2) then
          call real1strterm(ndim,vil(1),vil(2),vil(3),xcom(1),ycom(1),zcom(1),dr2ds6(1,2),d2r2dx26(1,2), &
                            d2r2dsdx6(1,1,2),d2r2ds26(1,1,2),lgrad2)
        endif
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijkdril*dr2ds6(ks,2)
            enddo
          endif
        endif
!
!  j-l
!
        if (lilocal) then
          xdrv(l) = xdrv(l) - deijkdrjl*xbnbr(nl,l)
          ydrv(l) = ydrv(l) - deijkdrjl*ybnbr(nl,l)
          zdrv(l) = zdrv(l) - deijkdrjl*zbnbr(nl,l)
!
          xdrv(j) = xdrv(j) + deijkdrjl*xbnbr(nl,l)
          ydrv(j) = ydrv(j) + deijkdrjl*ybnbr(nl,l)
          zdrv(j) = zdrv(j) + deijkdrjl*zbnbr(nl,l)
        endif
!
        if (lstr.or.lgrad2) then
          call real1strterm(ndim,xbnbr(nl,l),ybnbr(nl,l),zbnbr(nl,l),xcom(1),ycom(1),zcom(1),dr2ds6(1,3),d2r2dx26(1,3), &
                            d2r2dsdx6(1,1,3),d2r2ds26(1,1,3),lgrad2)
        endif
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijkdrjl*dr2ds6(ks,3)
            enddo
          endif
        endif
        if (lgrad2) then
!-----------------------
!  Second derivatives  |
!-----------------------
!
!  Order of distances:
!  d1 : 1 = ij / 2 = il / 3 = jl
!  d2 : 1 = 11 / 2 = 12 / 3 = 13 / 4 = 22 / 5 = 23 / 6 = 33
!
          d2cos1(1) = - rrij*rrij*(dcos1drij + rrij*rrjl - 2.0_dp*cos1*rrij**2)
          d2cos1(2) = rrjl*rrij**3
          d2cos1(3) = - rrij*rrjl**3 - rrij*rrij*dcos1drjl
          d2cos1(4) = 0.0_dp
          d2cos1(5) = rrij*rrjl**3
          d2cos1(6) = - rrjl*rrjl*(dcos1drjl + rrij*rrjl - 2.0_dp*cos1*rrjl**2)
!
          d2cos2(1) = - rrij*rrij*(dcos2drij + rrij*rril - 2.0_dp*cos2*rrij**2)
          d2cos2(2) = - rrij*rril**3 - rrij*rrij*dcos2dril
          d2cos2(3) = rril*rrij**3
          d2cos2(4) = - rril*rril*(dcos2dril + rrij*rril - 2.0_dp*cos2*rril**2)
          d2cos2(5) = rrij*rril**3
          d2cos2(6) = 0.0_dp
!
          d2cos3(1) = 0.0_dp
          d2cos3(2) = rrjl*rril**3
          d2cos3(3) = rril*rrjl**3
          d2cos3(4) = - rril*rril*(dcos3dril + rril*rrjl - 2.0_dp*cos3*rril**2)
          d2cos3(5) = - rril*rrjl**3 - rril*rril*dcos3drjl
          d2cos3(6) = - rrjl*rrjl*(dcos3drjl + rril*rrjl - 2.0_dp*cos3*rrjl**2)
!
!  i-j / i-j
!
          if (lilocal.or.ljlocal) then
            d2trm = eijk*rrij3*(18.0_dp*rrij3*rij2 - 3.0_dp*rrij) &
                    - (18.0_dp*c9*rrijk9*rrij3*rij)*(dcos1drij*cos2*cos3 + cos1*dcos2drij*cos3 + dcos3drij*cos1*cos2) &
                    + (3.0_dp*c9*rrijk9)*(d2cos1(1)*cos2*cos3 + cos1*d2cos2(1)*cos3 + cos1*cos2*d2cos3(1) &
                                         + dcos1drij*(dcos2drij*cos3 + cos2*dcos3drij) &
                                         + dcos2drij*(dcos1drij*cos3 + cos1*dcos3drij) &
                                         + dcos3drij*(dcos1drij*cos2 + cos1*dcos2drij))
!
            call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,vij(1),vij(2),vij(3), &
                              deijkdrij,d2trm,dr2ds6(1,1),d2r2dx26(1,1),d2r2ds26(1,1,1),d2r2dsdx6(1,1,1),lilocal)
          endif
!
!  i-j / i-l
!
          d2trm = 9.0_dp*eijk*rij*ril*rrij3*rril3 &
                  - (9.0_dp*c9*rrijk9*rrij3*rij)*(dcos1dril*cos2*cos3 + cos1*dcos2dril*cos3 + dcos3dril*cos1*cos2) &
                  - (9.0_dp*c9*rrijk9*rril3*ril)*(dcos1drij*cos2*cos3 + cos1*dcos2drij*cos3 + dcos3drij*cos1*cos2) &
                  + (3.0_dp*c9*rrijk9)*(d2cos1(2)*cos2*cos3 + cos1*d2cos2(2)*cos3 + cos1*cos2*d2cos3(2) &
                                       + dcos1drij*(dcos2dril*cos3 + cos2*dcos3dril) &
                                       + dcos2drij*(dcos1dril*cos3 + cos1*dcos3dril) &
                                       + dcos3drij*(dcos1dril*cos2 + cos1*dcos2dril))
!
          call add_drv2_2dm(lilocal,ljlocal,lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf, &
                            lx,ly,lz,lxf,lyf,lzf,vij(1),vij(2),vij(3),vil(1),vil(2),vil(3),d2trm, &
                            dr2ds6(1,1),dr2ds6(1,2),lilocal)
!
!  i-j / l-j
!
          d2trm = 9.0_dp*eijk*rij*rjl*rrij3*rrjl3 &
                  - (9.0_dp*c9*rrijk9*rrij3*rij)*(dcos1drjl*cos2*cos3 + cos1*dcos2drjl*cos3 + dcos3drjl*cos1*cos2) &
                  - (9.0_dp*c9*rrijk9*rrjl3*rjl)*(dcos1drij*cos2*cos3 + cos1*dcos2drij*cos3 + dcos3drij*cos1*cos2) &
                  + (3.0_dp*c9*rrijk9)*(d2cos1(3)*cos2*cos3 + cos1*d2cos2(3)*cos3 + cos1*cos2*d2cos3(3) &
                                       + dcos1drij*(dcos2drjl*cos3 + cos2*dcos3drjl) &
                                       + dcos2drij*(dcos1drjl*cos3 + cos1*dcos3drjl) &
                                       + dcos3drij*(dcos1drjl*cos2 + cos1*dcos2drjl))
!
          call add_drv2_2dm(lilocal,ljlocal,lllocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,vij(1),vij(2),vij(3),xbnbr(nl,l),ybnbr(nl,l),zbnbr(nl,l),d2trm, &
                            dr2ds6(1,1),dr2ds6(1,3),lilocal)
!
!  i-l / i-l
!
          if (lilocal.or.lllocal) then
            d2trm = eijk*rril3*(18.0_dp*rril3*ril2 - 3.0_dp*rril) &
                    - (18.0_dp*c9*rrijk9*rril3*ril)*(dcos1dril*cos2*cos3 + cos1*dcos2dril*cos3 + dcos3dril*cos1*cos2) &
                    + (3.0_dp*c9*rrijk9)*(d2cos1(4)*cos2*cos3 + cos1*d2cos2(4)*cos3 + cos1*cos2*d2cos3(4) &
                                         + dcos1dril*(dcos2dril*cos3 + cos2*dcos3dril) &
                                         + dcos2dril*(dcos1dril*cos3 + cos1*dcos3dril) &
                                         + dcos3dril*(dcos1dril*cos2 + cos1*dcos2dril))
!
            call add_drv2_1dm(lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,vil(1),vil(2),vil(3), &
                              deijkdril,d2trm,dr2ds6(1,2),d2r2dx26(1,2), d2r2ds26(1,1,2),d2r2dsdx6(1,1,2),lilocal)
          endif
!
!  i-l / l-j
!
          d2trm = 9.0_dp*eijk*ril*rjl*rril3*rrjl3 &
                  - (9.0_dp*c9*rrijk9*rril3*ril)*(dcos1drjl*cos2*cos3 + cos1*dcos2drjl*cos3 + dcos3drjl*cos1*cos2) &
                  - (9.0_dp*c9*rrijk9*rrjl3*rjl)*(dcos1dril*cos2*cos3 + cos1*dcos2dril*cos3 + dcos3dril*cos1*cos2) &
                  + (3.0_dp*c9*rrijk9)*(d2cos1(5)*cos2*cos3 + cos1*d2cos2(5)*cos3 + cos1*cos2*d2cos3(5) &
                                       + dcos1dril*(dcos2drjl*cos3 + cos2*dcos3drjl) &
                                       + dcos2dril*(dcos1drjl*cos3 + cos1*dcos3drjl) &
                                       + dcos3dril*(dcos1drjl*cos2 + cos1*dcos2drjl))
!
          call add_drv2_2dm(lilocal,lllocal,lllocal,ljlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,vil(1),vil(2),vil(3),xbnbr(nl,l),ybnbr(nl,l),zbnbr(nl,l),d2trm, &
                            dr2ds6(1,2),dr2ds6(1,3),lilocal)
!
!  j-l / j-l
!
          if (ljlocal.or.lllocal) then
            d2trm = eijk*rrjl3*(18.0_dp*rrjl3*rjl2 - 3.0_dp*rrjl) &
                    - (18.0_dp*c9*rrijk9*rrjl3*rjl)*(dcos1drjl*cos2*cos3 + cos1*dcos2drjl*cos3 + dcos3drjl*cos1*cos2) &
                    + (3.0_dp*c9*rrijk9)*(d2cos1(6)*cos2*cos3 + cos1*d2cos2(6)*cos3 + cos1*cos2*d2cos3(6) &
                                         + dcos1drjl*(dcos2drjl*cos3 + cos2*dcos3drjl) &
                                         + dcos2drjl*(dcos1drjl*cos3 + cos1*dcos3drjl) &
                                         + dcos3drjl*(dcos1drjl*cos2 + cos1*dcos2drjl))
!
            call add_drv2_1dm(lllocal,ljlocal,lx,ly,lz,lxf,lyf,lzf,jx,jy,jz,jxf,jyf,jzf,xbnbr(nl,l),ybnbr(nl,l),zbnbr(nl,l), &
                              deijkdrjl,d2trm,dr2ds6(1,3),d2r2dx26(1,3),d2r2ds26(1,1,3),d2r2dsdx6(1,1,3),lilocal)
          endif
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
  call gfnff_ehbd(ehb,exb,nnbr_bond,maxnbr,nbrno_bond,rbnbr,xbnbr,ybnbr,zbnbr,lgrad1,lgrad2)
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
    call sumall(sum0,sum1,9_i4,"gfnff_energyd_s","energies")
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
  endif
!
  if (lverbose) then
    if (ioproc) then
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
  endif
!
  if (lgrad2) then
!
!  Sum diagonals
!
    ix = -2
    iy = -1
    iz =  0
    do ii = 1,natomsonnode
      i = node2atom(ii)
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      indi = 3*(i-1)
      ixf = indi + 1
      iyf = indi + 2
      izf = indi + 3
      derv2(ixf,ix) = derv2d(ix)
      derv2(ixf,iy) = 0.0_dp
      derv2(ixf,iz) = 0.0_dp
      derv2(iyf,ix) = 0.0_dp
      derv2(iyf,iy) = derv2d(iy)
      derv2(iyf,iz) = 0.0_dp
      derv2(izf,ix) = 0.0_dp
      derv2(izf,iy) = 0.0_dp
      derv2(izf,iz) = derv2d(iz)
      do j = 1,numat
        if (i.ne.j) then
          indj = 3*(j-1)
          jx = indj + 1
          jy = indj + 2
          jz = indj + 3
          derv2(ixf,ix) = derv2(ixf,ix) - derv2(jx,ix)
          derv2(ixf,iy) = derv2(ixf,iy) - derv2(jy,ix)
          derv2(ixf,iz) = derv2(ixf,iz) - derv2(jz,ix)
          derv2(iyf,ix) = derv2(iyf,ix) - derv2(jx,iy)
          derv2(iyf,iy) = derv2(iyf,iy) - derv2(jy,iy)
          derv2(iyf,iz) = derv2(iyf,iz) - derv2(jz,iy)
          derv2(izf,ix) = derv2(izf,ix) - derv2(jx,iz)
          derv2(izf,iy) = derv2(izf,iy) - derv2(jy,iz)
          derv2(izf,iz) = derv2(izf,iz) - derv2(jz,iz)
        endif
      enddo
    enddo
  endif
!
!  Free local memory
!
  if (ldohbcn) then
    deallocate(cn_hb,stat=status)
    if (status/=0) call deallocate_error('gfnff_energyd_s','cn_hb')
  endif
!
  deallocate(d2gwdcn2,stat=status)
  if (status/=0) call deallocate_error('gfnff_energyd_s','d2gwdcn2')
  deallocate(dgwdcn,stat=status)
  if (status/=0) call deallocate_error('gfnff_energyd_s','dgwdcn')
  deallocate(gw,stat=status)
  if (status/=0) call deallocate_error('gfnff_energyd_s','gw')
  deallocate(d2cn,stat=status)
  if (status/=0) call deallocate_error('gfnff_energyd_s','d2cn')
  deallocate(dcn,stat=status)
  if (status/=0) call deallocate_error('gfnff_energyd_s','dcn')
  deallocate(cn,stat=status)
  if (status/=0) call deallocate_error('gfnff_energyd_s','cn')
!
  time2 = g_cpu_time()
  tgfnff = tgfnff + time2 - time1
#ifdef TRACE
  call trace_out('gfnff_energyd_s')
#endif
!
  return
  end
