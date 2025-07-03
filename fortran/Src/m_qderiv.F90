module m_qderiv
!
!  Module containing routines for charge derivatives
!
!   4/23 Subroutine added specifically for d2chargeself term
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
!  Julian Gale, Curtin University, April 2023
!
  public d2charge, d2charged, d2charge2D, d2charge2Dd, d2charge3
  public d2chargep, d2chargepd, d2charge3p, d2chargefc, d2charge3fc
  public d2chargep2D, d2chargep2Dd, d2chargeself, d2chargeselfp
  public d2charge_finish, d2charged_finish, d2chargepd_finish
  public d2chargeselfpd, d2chargeselfd

CONTAINS

  subroutine d2charge(i,j,nor,ix,iy,iz,jx,jy,jz,lopi,lopj,dei,dej,d1xi,d1yi,d1zi, &
                      d1xj,d1yj,d1zj,ds1i,ds1j,d2i2,d2ij,d2j2,d2self,dei0,dej0, &
                      lreal,lDoSelf,nbosptr)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from variable charge models
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!
!   2/01 Structure of rpd changed to suit 2-D case
!  10/02 ReaxFF modifications added
!  11/02 Second strain derivatives corrected
!  11/02 Error in derv3 calculation for QEq(H) case corrected
!   9/04 Order of terms in dqds switched
!   9/04 Modified to handle general charge derivatives and not just EEM
!   9/04 Contribution from d2q/dalpha.dbeta added
!   9/04 dei0/dej0 arguments added
!   9/04 xtmp/ytmp/ztmp removed from arguments and replace by d1ix, d1jx
!        etc being precomputed prior to entry to routine and passed in
!        as d1xi, d1xj etc
!   9/04 rpd/mdis removed & ds1i/ds1j made into sum as per above
!   9/04 lDoSelf flag added to control addition of self term for EEM
!  10/04 Position of adding to d2qds2 to sderv2 corrected
!   7/05 Streitz and Mintmire modifications added
!   7/05 lReverse set to true for k.le.(i.or.j)
!  11/07 Unused variables removed
!   3/08 ReaxFF charge scheme removed as this is handled elsewhere
!   9/12 Pacha added
!   2/18 Trace added
!   5/18 Multiple qranges added
!   1/19 Correction to transfer of ds1i/ds1j 
!   1/19 Correction to volume related term at finite strain
!   1/19 Error in lReverse flag setting corrected
!   2/19 Corrections to alpha/beta order for derv2
!   7/19 lReverse and alpha/beta changes reversed
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   2/21 Modifications added for GFNFF
!   6/21 Modified so that nor sums are outside loops
!   6/21 lexact_d2q added
!   6/21 Exact d2q terms rearranged to speed up
!   7/21 Timing added
!   7/21 lexact_d2q replaced by loldd2q
!   7/21 New algorithm to avoid 4th power scaling added
!   7/21 QEq terms corrected
!  10/21 lgfnff moved to control
!  12/22 Modified to handle reaxFF
!  12/22 Optional argument, nbosptr added
!   9/23 Trap for shells added
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
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use optimisation
  use m_strain,       only : strainddetds, straindet
  use reaxFFdata,     only : reaxFFmu
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)           :: i
  integer(i4), intent(in)           :: ix
  integer(i4), intent(in)           :: iy
  integer(i4), intent(in)           :: iz
  integer(i4), intent(in)           :: j
  integer(i4), intent(in)           :: jx
  integer(i4), intent(in)           :: jy
  integer(i4), intent(in)           :: jz
  integer(i4), intent(in)           :: nor
  integer(i4), intent(in), optional :: nbosptr(*)       ! Only required for reaxFF calls - points to reaxFF species from atoms
  logical,     intent(in)           :: lopi
  logical,     intent(in)           :: lopj
  logical,     intent(in)           :: lreal
  logical,     intent(in)           :: lDoSelf
  real(dp),    intent(in)           :: d1xi          ! d2E/dqi.dx
  real(dp),    intent(in)           :: d1yi          ! d2E/dqi.dy
  real(dp),    intent(in)           :: d1zi          ! d2E/dqi.dz
  real(dp),    intent(in)           :: d1xj          ! d2E/dqj.dx
  real(dp),    intent(in)           :: d1yj          ! d2E/dqj.dy
  real(dp),    intent(in)           :: d1zj          ! d2E/dqj.dz
  real(dp),    intent(in)           :: d2i2(*)       ! d2E/dqi2
  real(dp),    intent(in)           :: d2ij(*)       ! d2E/dqi.dqj
  real(dp),    intent(in)           :: d2j2(*)       ! d2E/dqj2
  real(dp),    intent(in)           :: d2self        ! Self term
  real(dp),    intent(in)           :: dei(*)        ! dE/dqi
  real(dp),    intent(in)           :: dej(*)        ! dE/dqj
  real(dp),    intent(in)           :: dei0          ! dE(selfterm)/dqi
  real(dp),    intent(in)           :: dej0          ! dE(selfterm)/dqj
  real(dp),    intent(in)           :: ds1i(*)       ! d2E/dqi.d(epsilon)
  real(dp),    intent(in)           :: ds1j(*)       ! d2E/dqj.d(epsilon)
!
!  Local variables
!
  integer(i4)                       :: iix
  integer(i4)                       :: iiy
  integer(i4)                       :: iiz
  integer(i4)                       :: ind
  integer(i4)                       :: indi
  integer(i4)                       :: indj
  integer(i4)                       :: indk
  integer(i4)                       :: indl
  integer(i4)                       :: iv
  integer(i4)                       :: ixl
  integer(i4)                       :: iyl
  integer(i4)                       :: izl
  integer(i4)                       :: jjx
  integer(i4)                       :: jjy
  integer(i4)                       :: jjz
  integer(i4)                       :: jxl
  integer(i4)                       :: jyl
  integer(i4)                       :: jzl
  integer(i4)                       :: k
  integer(i4)                       :: kk
  integer(i4)                       :: kl
  integer(i4)                       :: kx
  integer(i4)                       :: ky
  integer(i4)                       :: kz
  integer(i4)                       :: kxl
  integer(i4)                       :: kyl
  integer(i4)                       :: kzl
  integer(i4)                       :: l
  integer(i4)                       :: lx
  integer(i4)                       :: ly
  integer(i4)                       :: lz
  integer(i4)                       :: lxl
  integer(i4)                       :: lyl
  integer(i4)                       :: lzl
  integer(i4)                       :: m
  integer(i4)                       :: mnxx
  integer(i4)                       :: mnxy
  integer(i4)                       :: mnxz
  integer(i4)                       :: mnyx
  integer(i4)                       :: mnyy
  integer(i4)                       :: mnyz
  integer(i4)                       :: mnzx
  integer(i4)                       :: mnzy
  integer(i4)                       :: mnzz
  integer(i4)                       :: mx
  integer(i4)                       :: my
  integer(i4)                       :: mz
  integer(i4)                       :: n
  integer(i4)                       :: nqr
  integer(i4)                       :: nx
  integer(i4)                       :: nk
  logical                           :: lopk
  logical                           :: lopl
  logical                           :: lNonIJQDeriv
  logical                           :: lReverse
  real(dp)                          :: d1ix
  real(dp)                          :: d1iy
  real(dp)                          :: d1iz
  real(dp)                          :: d1jx
  real(dp)                          :: d1jy
  real(dp)                          :: d1jz
  real(dp)                          :: d2i2s
  real(dp)                          :: d2ijs
  real(dp)                          :: d2j2s
  real(dp)                          :: d2qk
  real(dp)                          :: d1si
  real(dp)                          :: d2si
  real(dp)                          :: d3si
  real(dp)                          :: d4si
  real(dp)                          :: d5si
  real(dp)                          :: d6si
  real(dp)                          :: d1sj
  real(dp)                          :: d2sj
  real(dp)                          :: d3sj
  real(dp)                          :: d4sj
  real(dp)                          :: d5sj
  real(dp)                          :: d6sj
  real(dp)                          :: deisum
  real(dp)                          :: dejsum
  real(dp)                          :: dikx
  real(dp)                          :: diky
  real(dp)                          :: dikz
  real(dp)                          :: d2ikx
  real(dp)                          :: d2iky
  real(dp)                          :: d2ikz
  real(dp)                          :: dilx
  real(dp)                          :: dily
  real(dp)                          :: dilz
  real(dp)                          :: dis1
  real(dp)                          :: dis2
  real(dp)                          :: dis3
  real(dp)                          :: dis4
  real(dp)                          :: dis5
  real(dp)                          :: dis6
  real(dp)                          :: djkx
  real(dp)                          :: djky
  real(dp)                          :: djkz
  real(dp)                          :: d2jkx
  real(dp)                          :: d2jky
  real(dp)                          :: d2jkz
  real(dp)                          :: djlx
  real(dp)                          :: djly
  real(dp)                          :: djlz
  real(dp)                          :: djs1
  real(dp)                          :: djs2
  real(dp)                          :: djs3
  real(dp)                          :: djs4
  real(dp)                          :: djs5
  real(dp)                          :: djs6
  real(dp)                          :: dkix
  real(dp)                          :: dkiy
  real(dp)                          :: dkiz
  real(dp)                          :: dkjx
  real(dp)                          :: dkjy
  real(dp)                          :: dkjz
  real(dp)                          :: dk1
  real(dp)                          :: dk2
  real(dp)                          :: dk3
  real(dp)                          :: dk4
  real(dp)                          :: dk5
  real(dp)                          :: dk6
  real(dp)                          :: dsi(6)
  real(dp)                          :: dsi2(6)
  real(dp)                          :: dsj(6)
  real(dp)                          :: dsj2(6)
  real(dp)                          :: g_cpu_time
  real(dp)                          :: time1
  real(dp)                          :: time2
  real(dp)                          :: ock
  real(dp)                          :: qlk
  real(dp)                          :: zetah0
  real(dp), parameter               :: tsqrt2pi = 0.797884560802866_dp
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charge')
#endif
!
  time1 = g_cpu_time()
!
!  For ReaxFF check that nbosptr is present
!
  if (lreaxFF) then
    if (.not.present(nbosptr)) then
      call outerror('d2charge called for reaxFF without passing nbosptr',0_i4)
      call stopnow('d2charge')
    endif
  endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  if (lopi.and.lopj) then
    iix = ix
    iiy = iy
    iiz = iz
    jjx = jx
    jjy = jy
    jjz = jz
  elseif (lopi) then
    iix = ix
    iiy = iy
    iiz = iz
    jjx = ix
    jjy = iy
    jjz = iz
  elseif (lopj) then
    iix = jx
    iiy = jy
    iiz = jz
    jjx = jx
    jjy = jy
    jjz = jz
  endif
  indi = 3*(i - 1)
  indj = 3*(j - 1)
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d2ijs = d2ijs + d2ij(iv)
    d2i2s = d2i2s + d2i2(iv)
    d2j2s = d2j2s + d2j2(iv)
  enddo
  d1ix = d1xi
  d1iy = d1yi
  d1iz = d1zi
  d1jx = d1xj
  d1jy = d1yj
  d1jz = d1zj
  if (leem) then
    if (nat(i).gt.maxele) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    endif
  endif
  if (.not.loldd2q) then
!******************
!  New algorithm  *
!******************
    d2edqdq(j,i) = d2edqdq(j,i) + d2ijs
    d2edq2(i) = d2edq2(i) + d2i2s
    if (i.ne.j) then
      d2edqdq(i,j) = d2edqdq(i,j) + d2ijs
      d2edq2(j) = d2edq2(j) + d2j2s
    endif
  endif
!
  if (lstr) then
    if (ndim.eq.3) then
      dis1 = dqds(1,i)
      dis2 = dqds(2,i)
      dis3 = dqds(3,i)
      dis4 = dqds(4,i)
      dis5 = dqds(5,i)
      dis6 = dqds(6,i)
      djs1 = dqds(1,j)
      djs2 = dqds(2,j)
      djs3 = dqds(3,j)
      djs4 = dqds(4,j)
      djs5 = dqds(5,j)
      djs6 = dqds(6,j)
    elseif (ndim.eq.2) then
      dis1 = dqds(1,i)
      dis2 = dqds(2,i)
      dis3 = dqds(3,i)
      djs1 = dqds(1,j)
      djs2 = dqds(2,j)
      djs3 = dqds(3,j)
    elseif (ndim.eq.1) then
      dis1 = dqds(1,i)
      djs1 = dqds(1,j)
    endif
    if (lreal) then
!
!  Real space
!
      do kl = 1,nstrains
        dsi(kl) = ds1i(kl)
        dsj(kl) = ds1j(kl)
        dsi2(kl) = dsi(kl)
        dsj2(kl) = dsj(kl)
      enddo
    else
!
!  Reciprocal space
!
      do kl = 1,nstrains
        dsi(kl) = ds1i(kl)
        dsj(kl) = ds1j(kl)
        dsi2(kl) = dsi(kl)
        dsj2(kl) = dsj(kl)
      enddo
      deisum = 0.0_dp
      dejsum = 0.0_dp
      do iv = 1,nor
        deisum = deisum + dei(iv)
        dejsum = dejsum + dej(iv)
      enddo
!
      if (ndim.eq.3) then
        if (lfinitestrain) then
          do kl = 1,3
            dsi(kl) = dsi(kl) - deisum*strainddetds(kl)*straindet
            dsj(kl) = dsj(kl) - dejsum*strainddetds(kl)*straindet
          enddo
        else
          do kl = 1,3
            dsi(kl) = dsi(kl) - deisum
            dsj(kl) = dsj(kl) - dejsum
          enddo
        endif
      elseif (ndim.eq.2) then
        if (lfinitestrain) then
          do kl = 1,3
            dsi(kl) = dsi(kl) - deisum*strainddetds(kl)*straindet
            dsj(kl) = dsj(kl) - dejsum*strainddetds(kl)*straindet
          enddo
        else
          do kl = 1,2
            dsi(kl) = dsi(kl) - deisum
            dsj(kl) = dsj(kl) - dejsum
          enddo
        endif
      elseif (ndim.eq.1) then
        dsi(1) = dsi(1) - deisum
        dsj(1) = dsj(1) - dejsum
      endif
      if (leem.or.lgfnff) then
        do kl = 1,nstrains
          dsi2(kl) = dsi(kl)
          dsj2(kl) = dsj(kl)
        enddo
      endif
    endif
  endif
!
!  The following flag causes the second derivatives of the charge distribution 
!  with respect to two atoms neither of which are the i/j of the energy term
!  to be calculated. At the moment this contribution is not used in any method
!  and so the flag is set to false for efficiency.
!
  lNonIJQDeriv = .false.
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff.and..not.lreaxFF) then
    deisum = dei0
    dejsum = dej0
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!   
!  Strain - strain contribution
!   
    if (lstr.and.ndim.gt.0) then
      ind = 0
      do kk = 1,nstrains
        do kl = 1,kk-1
          ind = ind + 1
          sderv2(kl,kk) = sderv2(kl,kk) + deisum*d2qds2(ind,i) + dejsum*d2qds2(ind,j)
          sderv2(kk,kl) = sderv2(kk,kl) + deisum*d2qds2(ind,i) + dejsum*d2qds2(ind,j)
        enddo
        ind = ind + 1
        sderv2(kk,kk) = sderv2(kk,kk) + deisum*d2qds2(ind,i) + dejsum*d2qds2(ind,j)
      enddo
    endif
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(i)
      k = nqatomptr(m,i)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (i.gt.k) then
        derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,i)
        derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,i)
        derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,i)
        derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,i)
        derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,i)
        derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,i)
        derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,i)
        derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,i)
        derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,i)
      elseif (i.lt.k) then
        derv2(ix,kx) = derv2(ix,kx) + deisum*d2qdxyz2(mnxx,i)
        derv2(iy,kx) = derv2(iy,kx) + deisum*d2qdxyz2(mnxy,i)
        derv2(iz,kx) = derv2(iz,kx) + deisum*d2qdxyz2(mnxz,i)
        derv2(ix,ky) = derv2(ix,ky) + deisum*d2qdxyz2(mnxy,i)
        derv2(iy,ky) = derv2(iy,ky) + deisum*d2qdxyz2(mnyy,i)
        derv2(iz,ky) = derv2(iz,ky) + deisum*d2qdxyz2(mnyz,i)
        derv2(ix,kz) = derv2(ix,kz) + deisum*d2qdxyz2(mnxz,i)
        derv2(iy,kz) = derv2(iy,kz) + deisum*d2qdxyz2(mnyz,i)
        derv2(iz,kz) = derv2(iz,kz) + deisum*d2qdxyz2(mnzz,i)
      endif
!
!  Mixed strain contribution
!
      if (lstr.and.ndim.gt.0) then
        do kk = 1,nstrains
          derv3(kx,kk) = derv3(kx,kk) + deisum*d2qdxyzs(kk,mx,i)
          derv3(ky,kk) = derv3(ky,kk) + deisum*d2qdxyzs(kk,my,i)
          derv3(kz,kk) = derv3(kz,kk) + deisum*d2qdxyzs(kk,mz,i)
          derv3(ix,kk) = derv3(ix,kk) - deisum*d2qdxyzs(kk,mx,i)
          derv3(iy,kk) = derv3(iy,kk) - deisum*d2qdxyzs(kk,my,i)
          derv3(iz,kk) = derv3(iz,kk) - deisum*d2qdxyzs(kk,mz,i)
        enddo
      endif
    enddo
!
    do m = 1,nqatoms(j)
      k = nqatomptr(m,j)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (j.gt.k) then
        derv2(kx,jx) = derv2(kx,jx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(ky,jx) = derv2(ky,jx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(kz,jx) = derv2(kz,jx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(kx,jy) = derv2(kx,jy) + dejsum*d2qdxyz2(mnxy,j)
        derv2(ky,jy) = derv2(ky,jy) + dejsum*d2qdxyz2(mnyy,j)
        derv2(kz,jy) = derv2(kz,jy) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kx,jz) = derv2(kx,jz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(ky,jz) = derv2(ky,jz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kz,jz) = derv2(kz,jz) + dejsum*d2qdxyz2(mnzz,j)
      elseif (j.lt.k) then
        derv2(jx,kx) = derv2(jx,kx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(jy,kx) = derv2(jy,kx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jz,kx) = derv2(jz,kx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jx,ky) = derv2(jx,ky) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jy,ky) = derv2(jy,ky) + dejsum*d2qdxyz2(mnyy,j)
        derv2(jz,ky) = derv2(jz,ky) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jx,kz) = derv2(jx,kz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jy,kz) = derv2(jy,kz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jz,kz) = derv2(jz,kz) + dejsum*d2qdxyz2(mnzz,j)
      endif
!
!  Mixed strain contribution
!
      if (lstr.and.ndim.gt.0) then
        do kk = 1,nstrains
          derv3(kx,kk) = derv3(kx,kk) + dejsum*d2qdxyzs(kk,mx,j)
          derv3(ky,kk) = derv3(ky,kk) + dejsum*d2qdxyzs(kk,my,j)
          derv3(kz,kk) = derv3(kz,kk) + dejsum*d2qdxyzs(kk,mz,j)
          derv3(jx,kk) = derv3(jx,kk) - dejsum*d2qdxyzs(kk,mx,j)
          derv3(jy,kk) = derv3(jy,kk) - dejsum*d2qdxyzs(kk,my,j)
          derv3(jz,kk) = derv3(jz,kk) - dejsum*d2qdxyzs(kk,mz,j)
        enddo
      endif
    enddo
!
!  Terms involving Q derivatives for non-i/j case
!
    if (lNonIJQDeriv) then
      do m = 2,nqatoms(i)
        k = nqatomptr(m,i)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,i)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,i)
          derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,i)
          derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,i)
          derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,i)
          derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,i)
          derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,i)
          derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,i)
          derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,i)
          derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,i)
        enddo
      enddo
!
      do m = 2,nqatoms(j)
        k = nqatomptr(m,j)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,j)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,j)
          derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,j)
          derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,j)
          derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,j)
          derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,j)
          derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,j)
          derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,j)
          derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,j)
          derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,j)
        enddo
      enddo
    endif
  endif
!
!  Loop over atoms to add ij-k contributions
!
!  Note indices must be set to cope with the fact that only
!  half of derv2 is constructed in upper triangle and the
!  effects of lopi/lopj/lopk if freezing is being used.
!
  kx = - 2
  ky = - 1
  kz =   0
  do k = 1,numat
    lopk = (.not.lfreeze.or.lopf(k))
    indk = 3*(k-1)
    if (lopk) then
      kx = kx + 3
      ky = ky + 3
      kz = kz + 3
    endif
!
    dikx = dqdxyz(indk+1,i)
    diky = dqdxyz(indk+2,i)
    dikz = dqdxyz(indk+3,i)
    djkx = dqdxyz(indk+1,j)
    djky = dqdxyz(indk+2,j)
    djkz = dqdxyz(indk+3,j)
!
    if (i.ne.k.and.(lopi.or.lopk)) then
      if (lopi.and.lopk) then
        if (k.le.i) then
          ixl = ix
          iyl = iy
          izl = iz
          kxl = kx
          kyl = ky
          kzl = kz
          lReverse = .true.
        else
          ixl = kx
          iyl = ky
          izl = kz
          kxl = ix
          kyl = iy
          kzl = iz
          lReverse = .true.
        endif
      elseif (lopi) then
        ixl = ix
        iyl = iy
        izl = iz
        kxl = ix
        kyl = iy
        kzl = iz
        lReverse = .false.
      elseif (lopk) then
        ixl = kx
        iyl = ky
        izl = kz
        kxl = kx
        kyl = ky
        kzl = kz
        lReverse = .true.
      endif
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-k
!
      if (lReverse) then
        derv2(kxl,ixl) = derv2(kxl,ixl) - d1ix*dikx
        derv2(kyl,ixl) = derv2(kyl,ixl) - d1iy*dikx
        derv2(kzl,ixl) = derv2(kzl,ixl) - d1iz*dikx
        derv2(kxl,iyl) = derv2(kxl,iyl) - d1ix*diky
        derv2(kyl,iyl) = derv2(kyl,iyl) - d1iy*diky
        derv2(kzl,iyl) = derv2(kzl,iyl) - d1iz*diky
        derv2(kxl,izl) = derv2(kxl,izl) - d1ix*dikz
        derv2(kyl,izl) = derv2(kyl,izl) - d1iy*dikz
        derv2(kzl,izl) = derv2(kzl,izl) - d1iz*dikz
!
        derv2(kxl,ixl) = derv2(kxl,ixl) - d1jx*djkx
        derv2(kyl,ixl) = derv2(kyl,ixl) - d1jy*djkx
        derv2(kzl,ixl) = derv2(kzl,ixl) - d1jz*djkx
        derv2(kxl,iyl) = derv2(kxl,iyl) - d1jx*djky
        derv2(kyl,iyl) = derv2(kyl,iyl) - d1jy*djky
        derv2(kzl,iyl) = derv2(kzl,iyl) - d1jz*djky
        derv2(kxl,izl) = derv2(kxl,izl) - d1jx*djkz
        derv2(kyl,izl) = derv2(kyl,izl) - d1jy*djkz
        derv2(kzl,izl) = derv2(kzl,izl) - d1jz*djkz
      else
        derv2(kxl,ixl) = derv2(kxl,ixl) - d1ix*dikx
        derv2(kyl,ixl) = derv2(kyl,ixl) - d1ix*diky
        derv2(kzl,ixl) = derv2(kzl,ixl) - d1ix*dikz
        derv2(kxl,iyl) = derv2(kxl,iyl) - d1iy*dikx
        derv2(kyl,iyl) = derv2(kyl,iyl) - d1iy*diky
        derv2(kzl,iyl) = derv2(kzl,iyl) - d1iy*dikz
        derv2(kxl,izl) = derv2(kxl,izl) - d1iz*dikx
        derv2(kyl,izl) = derv2(kyl,izl) - d1iz*diky
        derv2(kzl,izl) = derv2(kzl,izl) - d1iz*dikz
!
        derv2(kxl,ixl) = derv2(kxl,ixl) - d1jx*djkx
        derv2(kyl,ixl) = derv2(kyl,ixl) - d1jx*djky
        derv2(kzl,ixl) = derv2(kzl,ixl) - d1jx*djkz
        derv2(kxl,iyl) = derv2(kxl,iyl) - d1jy*djkx
        derv2(kyl,iyl) = derv2(kyl,iyl) - d1jy*djky
        derv2(kzl,iyl) = derv2(kzl,iyl) - d1jy*djkz
        derv2(kxl,izl) = derv2(kxl,izl) - d1jz*djkx
        derv2(kyl,izl) = derv2(kyl,izl) - d1jz*djky
        derv2(kzl,izl) = derv2(kzl,izl) - d1jz*djkz
      endif
    endif
    if (j.ne.k.and.(lopj.or.lopk)) then
      if (lopj.and.lopk) then
        if (k.le.j) then
          jxl = jx
          jyl = jy
          jzl = jz
          kxl = kx
          kyl = ky
          kzl = kz
          lReverse = .true.
        else
          jxl = kx
          jyl = ky
          jzl = kz
          kxl = jx
          kyl = jy
          kzl = jz
          lReverse = .true.
        endif
      elseif (lopj) then
        jxl = jx
        jyl = jy
        jzl = jz
        kxl = jx
        kyl = jy
        kzl = jz
        lReverse = .false.
      elseif (lopk) then
        jxl = kx
        jyl = ky
        jzl = kz
        kxl = kx
        kyl = ky
        kzl = kz
        lReverse = .true.
      endif
!
!  d2E/(d(alpha).dq) x dq/d(beta) : j-k
!
      if (lReverse) then
        derv2(kxl,jxl) = derv2(kxl,jxl) + d1jx*djkx
        derv2(kyl,jxl) = derv2(kyl,jxl) + d1jy*djkx
        derv2(kzl,jxl) = derv2(kzl,jxl) + d1jz*djkx
        derv2(kxl,jyl) = derv2(kxl,jyl) + d1jx*djky
        derv2(kyl,jyl) = derv2(kyl,jyl) + d1jy*djky
        derv2(kzl,jyl) = derv2(kzl,jyl) + d1jz*djky
        derv2(kxl,jzl) = derv2(kxl,jzl) + d1jx*djkz
        derv2(kyl,jzl) = derv2(kyl,jzl) + d1jy*djkz
        derv2(kzl,jzl) = derv2(kzl,jzl) + d1jz*djkz
!
        derv2(kxl,jxl) = derv2(kxl,jxl) + d1ix*dikx
        derv2(kyl,jxl) = derv2(kyl,jxl) + d1iy*dikx
        derv2(kzl,jxl) = derv2(kzl,jxl) + d1iz*dikx
        derv2(kxl,jyl) = derv2(kxl,jyl) + d1ix*diky
        derv2(kyl,jyl) = derv2(kyl,jyl) + d1iy*diky
        derv2(kzl,jyl) = derv2(kzl,jyl) + d1iz*diky
        derv2(kxl,jzl) = derv2(kxl,jzl) + d1ix*dikz
        derv2(kyl,jzl) = derv2(kyl,jzl) + d1iy*dikz
        derv2(kzl,jzl) = derv2(kzl,jzl) + d1iz*dikz
      else
        derv2(kxl,jxl) = derv2(kxl,jxl) + d1jx*djkx
        derv2(kyl,jxl) = derv2(kyl,jxl) + d1jx*djky
        derv2(kzl,jxl) = derv2(kzl,jxl) + d1jx*djkz
        derv2(kxl,jyl) = derv2(kxl,jyl) + d1jy*djkx
        derv2(kyl,jyl) = derv2(kyl,jyl) + d1jy*djky
        derv2(kzl,jyl) = derv2(kzl,jyl) + d1jy*djkz
        derv2(kxl,jzl) = derv2(kxl,jzl) + d1jz*djkx
        derv2(kyl,jzl) = derv2(kyl,jzl) + d1jz*djky
        derv2(kzl,jzl) = derv2(kzl,jzl) + d1jz*djkz
!
        derv2(kxl,jxl) = derv2(kxl,jxl) + d1ix*dikx
        derv2(kyl,jxl) = derv2(kyl,jxl) + d1ix*diky
        derv2(kzl,jxl) = derv2(kzl,jxl) + d1ix*dikz
        derv2(kxl,jyl) = derv2(kxl,jyl) + d1iy*dikx
        derv2(kyl,jyl) = derv2(kyl,jyl) + d1iy*diky
        derv2(kzl,jyl) = derv2(kzl,jyl) + d1iy*dikz
        derv2(kxl,jzl) = derv2(kxl,jzl) + d1iz*dikx
        derv2(kyl,jzl) = derv2(kyl,jzl) + d1iz*diky
        derv2(kzl,jzl) = derv2(kzl,jzl) + d1iz*dikz
      endif
    endif
    if ((lopi.or.lopj).and.lreal.and.lDoSelf) then
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
      nk = nat(k)
      ock = occuf(k)
      if (leem.and.nk.le.maxele) then
        if (lmultiqrange.and.neemrptr(k).ne.0) then
          nqr = nqrnow(neemrptr(k))
        else
          nqr = 1
        endif
        if (lqeq) then
          if (nk.ne.1) then
            d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
          else
            qlk = qf(k)
            zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
            d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp + 2.0_dp*qlk/zetah0)
          endif
        else
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        endif
      elseif (lgfnff) then
        d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
      elseif (lreaxFF) then
        d2qk = 2.0_dp*ock*reaxFFmu(nbosptr(k))
      else
        d2qk = 0.0_dp
      endif
      dkix = dqdxyz(indi+1,k)
      dkiy = dqdxyz(indi+2,k)
      dkiz = dqdxyz(indi+3,k)
      dkjx = dqdxyz(indj+1,k)
      dkjy = dqdxyz(indj+2,k)
      dkjz = dqdxyz(indj+3,k)
      derv2(jjx,iix) = derv2(jjx,iix) + d2qk*dkjx*dkix
      derv2(jjy,iix) = derv2(jjy,iix) + d2qk*dkjx*dkiy
      derv2(jjz,iix) = derv2(jjz,iix) + d2qk*dkjx*dkiz
      derv2(jjx,iiy) = derv2(jjx,iiy) + d2qk*dkjy*dkix
      derv2(jjy,iiy) = derv2(jjy,iiy) + d2qk*dkjy*dkiy
      derv2(jjz,iiy) = derv2(jjz,iiy) + d2qk*dkjy*dkiz
      derv2(jjx,iiz) = derv2(jjx,iiz) + d2qk*dkjz*dkix
      derv2(jjy,iiz) = derv2(jjy,iiz) + d2qk*dkjz*dkiy
      derv2(jjz,iiz) = derv2(jjz,iiz) + d2qk*dkjz*dkiz
      if (lstr.and.i.eq.j) then
!
!  Mixed strain terms from self energy - only need to do once for each atom - hence check on i = j
!
        if (lopi) then
          do kl = 1,nstrains
            derv3(iix,kl) = derv3(iix,kl) + d2qk*dqds(kl,k)*dkix
            derv3(iiy,kl) = derv3(iiy,kl) + d2qk*dqds(kl,k)*dkiy
            derv3(iiz,kl) = derv3(iiz,kl) + d2qk*dqds(kl,k)*dkiz
          enddo
        endif
!
!  Strain-strain derivatives - only need to do if i=j=k
!
        if (i.eq.k) then
          if (ndim.eq.3) then
            dk1 = dqds(1,k)
            dk2 = dqds(2,k)
            dk3 = dqds(3,k)
            dk4 = dqds(4,k)
            dk5 = dqds(5,k)
            dk6 = dqds(6,k)
            sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
            sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
            sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
            sderv2(4,1) = sderv2(4,1) + d2qk*dk4*dk1
            sderv2(5,1) = sderv2(5,1) + d2qk*dk5*dk1
            sderv2(6,1) = sderv2(6,1) + d2qk*dk6*dk1
            sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
            sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
            sderv2(4,2) = sderv2(4,2) + d2qk*dk4*dk2
            sderv2(5,2) = sderv2(5,2) + d2qk*dk5*dk2
            sderv2(6,2) = sderv2(6,2) + d2qk*dk6*dk2
            sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
            sderv2(4,3) = sderv2(4,3) + d2qk*dk4*dk3
            sderv2(5,3) = sderv2(5,3) + d2qk*dk5*dk3
            sderv2(6,3) = sderv2(6,3) + d2qk*dk6*dk3
            sderv2(4,4) = sderv2(4,4) + d2qk*dk4*dk4
            sderv2(5,4) = sderv2(5,4) + d2qk*dk5*dk4
            sderv2(6,4) = sderv2(6,4) + d2qk*dk6*dk4
            sderv2(5,5) = sderv2(5,5) + d2qk*dk5*dk5
            sderv2(6,5) = sderv2(6,5) + d2qk*dk6*dk5
            sderv2(6,6) = sderv2(6,6) + d2qk*dk6*dk6
          elseif (ndim.eq.2) then
            dk1 = dqds(1,k)
            dk2 = dqds(2,k)
            dk3 = dqds(3,k)
            sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
            sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
            sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
            sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
            sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
            sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
          elseif (ndim.eq.1) then
            dk1 = dqds(1,k)
            sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
          endif
        endif
      endif
    endif
    if (lstr.and.lopk) then
!
!  Mix strain-internal contribution
!
!  d2E/d(strain).dq x dq/d(alpha)
!
      do kl = 1,nstrains
        derv3(kx,kl) = derv3(kx,kl) + dsi2(kl)*dikx + dsj2(kl)*djkx
        derv3(ky,kl) = derv3(ky,kl) + dsi2(kl)*diky + dsj2(kl)*djky
        derv3(kz,kl) = derv3(kz,kl) + dsi2(kl)*dikz + dsj2(kl)*djkz
      enddo
!
!  d2E/dq.dq x dq/d(alpha) x dq/d(strain)
!
      do kl = 1,nstrains
        derv3(kx,kl) = derv3(kx,kl) + d2ijs*(dqds(kl,i)*djkx + dqds(kl,j)*dikx)
        derv3(ky,kl) = derv3(ky,kl) + d2ijs*(dqds(kl,i)*djky + dqds(kl,j)*diky)
        derv3(kz,kl) = derv3(kz,kl) + d2ijs*(dqds(kl,i)*djkz + dqds(kl,j)*dikz)
      enddo
      if (abs(d2i2s).gt.1.0d-8) then
        do kl = 1,nstrains
          derv3(kx,kl) = derv3(kx,kl) + d2i2s*dqds(kl,i)*dikx
          derv3(ky,kl) = derv3(ky,kl) + d2i2s*dqds(kl,i)*diky
          derv3(kz,kl) = derv3(kz,kl) + d2i2s*dqds(kl,i)*dikz
          derv3(ix,kl) = derv3(ix,kl) - d2i2s*dqds(kl,i)*dikx
          derv3(iy,kl) = derv3(iy,kl) - d2i2s*dqds(kl,i)*diky
          derv3(iz,kl) = derv3(iz,kl) - d2i2s*dqds(kl,i)*dikz
        enddo
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        do kl = 1,nstrains
          derv3(kx,kl) = derv3(kx,kl) + d2j2s*dqds(kl,j)*djkx
          derv3(ky,kl) = derv3(ky,kl) + d2j2s*dqds(kl,j)*djky
          derv3(kz,kl) = derv3(kz,kl) + d2j2s*dqds(kl,j)*djkz
          derv3(jx,kl) = derv3(jx,kl) - d2j2s*dqds(kl,j)*djkx
          derv3(jy,kl) = derv3(jy,kl) - d2j2s*dqds(kl,j)*djky
          derv3(jz,kl) = derv3(jz,kl) - d2j2s*dqds(kl,j)*djkz
        enddo
      endif
    endif
    if (loldd2q) then
!******************
!  Old algorithm  *
!******************
!
!  Loop over atoms to add ij-kl correction
!
      lx = - 2
      ly = - 1
      lz =   0
!
!  Scale charge derivatives for k by d2ijs
!
      d2ikx = d2ijs*dikx     
      d2iky = d2ijs*diky     
      d2ikz = d2ijs*dikz     
      d2jkx = d2ijs*djkx     
      d2jky = d2ijs*djky     
      d2jkz = d2ijs*djkz     
!
      do l = 1,k-1
        lopl = (.not.lfreeze.or.lopf(l))
        indl = 3*(l-1)
        if (lopl) then
          lx = lx + 3
          ly = ly + 3
          lz = lz + 3
        endif
        if (lopk.or.lopl) then
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
          dilx = dqdxyz(indl+1,i)
          dily = dqdxyz(indl+2,i)
          dilz = dqdxyz(indl+3,i)
          djlx = dqdxyz(indl+1,j)
          djly = dqdxyz(indl+2,j)
          djlz = dqdxyz(indl+3,j)
!
          if (lopk.and.lopl) then
            kxl = kx
            kyl = ky
            kzl = kz
            lxl = lx
            lyl = ly
            lzl = lz
          elseif (lopk) then
            kxl = kx
            kyl = ky
            kzl = kz
            lxl = kx
            lyl = ky
            lzl = kz
          elseif (lopl) then
            kxl = lx
            kyl = ly
            kzl = lz
            lxl = lx
            lyl = ly
            lzl = lz
          endif
!
          derv2(lxl,kxl) = derv2(lxl,kxl) + dilx*d2jkx + djlx*d2ikx
          derv2(lyl,kxl) = derv2(lyl,kxl) + dilx*d2jky + djlx*d2iky
          derv2(lzl,kxl) = derv2(lzl,kxl) + dilx*d2jkz + djlx*d2ikz
          derv2(lxl,kyl) = derv2(lxl,kyl) + dily*d2jkx + djly*d2ikx
          derv2(lyl,kyl) = derv2(lyl,kyl) + dily*d2jky + djly*d2iky
          derv2(lzl,kyl) = derv2(lzl,kyl) + dily*d2jkz + djly*d2ikz
          derv2(lxl,kzl) = derv2(lxl,kzl) + dilz*d2jkx + djlz*d2ikx
          derv2(lyl,kzl) = derv2(lyl,kzl) + dilz*d2jky + djlz*d2iky
          derv2(lzl,kzl) = derv2(lzl,kzl) + dilz*d2jkz + djlz*d2ikz
!
!  d2Edi2/d2Edj2 
!
          if (abs(d2i2s).gt.1.0d-8) then
            derv2(lxl,kxl) = derv2(lxl,kxl) + d2i2s*dilx*dikx
            derv2(lyl,kxl) = derv2(lyl,kxl) + d2i2s*dilx*diky
            derv2(lzl,kxl) = derv2(lzl,kxl) + d2i2s*dilx*dikz
            derv2(lxl,kyl) = derv2(lxl,kyl) + d2i2s*dily*dikx
            derv2(lyl,kyl) = derv2(lyl,kyl) + d2i2s*dily*diky
            derv2(lzl,kyl) = derv2(lzl,kyl) + d2i2s*dily*dikz
            derv2(lxl,kzl) = derv2(lxl,kzl) + d2i2s*dilz*dikx
            derv2(lyl,kzl) = derv2(lyl,kzl) + d2i2s*dilz*diky
            derv2(lzl,kzl) = derv2(lzl,kzl) + d2i2s*dilz*dikz
          endif
          if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
            derv2(lxl,kxl) = derv2(lxl,kxl) + d2j2s*djlx*djkx
            derv2(lyl,kxl) = derv2(lyl,kxl) + d2j2s*djlx*djky
            derv2(lzl,kxl) = derv2(lzl,kxl) + d2j2s*djlx*djkz
            derv2(lxl,kyl) = derv2(lxl,kyl) + d2j2s*djly*djkx
            derv2(lyl,kyl) = derv2(lyl,kyl) + d2j2s*djly*djky
            derv2(lzl,kyl) = derv2(lzl,kyl) + d2j2s*djly*djkz
            derv2(lxl,kzl) = derv2(lxl,kzl) + d2j2s*djlz*djkx
            derv2(lyl,kzl) = derv2(lyl,kzl) + d2j2s*djlz*djky
            derv2(lzl,kzl) = derv2(lzl,kzl) + d2j2s*djlz*djkz
          endif
        endif
!
!  End of loop over l
!
      enddo
    endif
!
!  End of loop over k
!
  enddo
!
!  Strain - charge derivatives terms
!
  if (lstr) then
!
!  Mixed strain-internal term
!
!  d2E/d(alpha).dq x dq/d(strain)
!
    if (lopi) then
      do kl = 1,nstrains
        derv3(ix,kl) = derv3(ix,kl) - d1ix*dqds(kl,i)
        derv3(iy,kl) = derv3(iy,kl) - d1iy*dqds(kl,i)
        derv3(iz,kl) = derv3(iz,kl) - d1iz*dqds(kl,i)
        derv3(ix,kl) = derv3(ix,kl) - d1jx*dqds(kl,j)
        derv3(iy,kl) = derv3(iy,kl) - d1jy*dqds(kl,j)
        derv3(iz,kl) = derv3(iz,kl) - d1jz*dqds(kl,j)
      enddo
    endif
    if (lopj) then
      do kl = 1,nstrains
        derv3(jx,kl) = derv3(jx,kl) + d1jx*dqds(kl,j)
        derv3(jy,kl) = derv3(jy,kl) + d1jy*dqds(kl,j)
        derv3(jz,kl) = derv3(jz,kl) + d1jz*dqds(kl,j)
        derv3(jx,kl) = derv3(jx,kl) + d1ix*dqds(kl,i)
        derv3(jy,kl) = derv3(jy,kl) + d1iy*dqds(kl,i)
        derv3(jz,kl) = derv3(jz,kl) + d1iz*dqds(kl,i)
      enddo
    endif
!
!  Strain-strain terms for charge derivatives
!
    if (ndim.eq.3) then
      d1si = dsi(1)
      d2si = dsi(2)
      d3si = dsi(3)
      d4si = dsi(4)
      d5si = dsi(5)
      d6si = dsi(6)
      d1sj = dsj(1)
      d2sj = dsj(2)
      d3sj = dsj(3)
      d4sj = dsj(4)
      d5sj = dsj(5)
      d6sj = dsj(6)
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
      sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2
      sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1
      sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3
      sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1
      sderv2(4,1) = sderv2(4,1) + d1si*dis4 + d1sj*djs4
      sderv2(4,1) = sderv2(4,1) + d4si*dis1 + d4sj*djs1
      sderv2(5,1) = sderv2(5,1) + d1si*dis5 + d1sj*djs5
      sderv2(5,1) = sderv2(5,1) + d5si*dis1 + d5sj*djs1
      sderv2(6,1) = sderv2(6,1) + d1si*dis6 + d1sj*djs6
      sderv2(6,1) = sderv2(6,1) + d6si*dis1 + d6sj*djs1
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2+ d2sj*djs2)
      sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3
      sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2
      sderv2(4,2) = sderv2(4,2) + d2si*dis4 + d2sj*djs4
      sderv2(4,2) = sderv2(4,2) + d4si*dis2 + d4sj*djs2
      sderv2(5,2) = sderv2(5,2) + d2si*dis5 + d2sj*djs5
      sderv2(5,2) = sderv2(5,2) + d5si*dis2 + d5sj*djs2
      sderv2(6,2) = sderv2(6,2) + d2si*dis6 + d2sj*djs6
      sderv2(6,2) = sderv2(6,2) + d6si*dis2 + d6sj*djs2
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3)
      sderv2(4,3) = sderv2(4,3) + d3si*dis4 + d3sj*djs4
      sderv2(4,3) = sderv2(4,3) + d4si*dis3 + d4sj*djs3
      sderv2(5,3) = sderv2(5,3) + d3si*dis5 + d3sj*djs5
      sderv2(5,3) = sderv2(5,3) + d5si*dis3 + d5sj*djs3
      sderv2(6,3) = sderv2(6,3) + d3si*dis6 + d3sj*djs6
      sderv2(6,3) = sderv2(6,3) + d6si*dis3 + d6sj*djs3
      sderv2(4,4) = sderv2(4,4) + 2.0_dp*(d4si*dis4 + d4sj*djs4)
      sderv2(5,4) = sderv2(5,4) + d4si*dis5 + d4sj*djs5
      sderv2(5,4) = sderv2(5,4) + d5si*dis4 + d5sj*djs4
      sderv2(6,4) = sderv2(6,4) + d4si*dis6 + d4sj*djs6
      sderv2(6,4) = sderv2(6,4) + d6si*dis4 + d6sj*djs4
      sderv2(5,5) = sderv2(5,5) + 2.0_dp*(d5si*dis5 + d5sj*djs5)
      sderv2(6,5) = sderv2(6,5) + d5si*dis6 + d5sj*djs6
      sderv2(6,5) = sderv2(6,5) + d6si*dis5 + d6sj*djs5
      sderv2(6,6) = sderv2(6,6) + 2.0_dp*(d6si*dis6 + d6sj*djs6)
!
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
      sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
      sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
      sderv2(4,1) = sderv2(4,1) + d2ijs*(dis1*djs4 + dis4*djs1)
      sderv2(5,1) = sderv2(5,1) + d2ijs*(dis1*djs5 + dis5*djs1)
      sderv2(6,1) = sderv2(6,1) + d2ijs*(dis1*djs6 + dis6*djs1)
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
      sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
      sderv2(4,2) = sderv2(4,2) + d2ijs*(dis2*djs4 + dis4*djs2)
      sderv2(5,2) = sderv2(5,2) + d2ijs*(dis2*djs5 + dis5*djs2)
      sderv2(6,2) = sderv2(6,2) + d2ijs*(dis2*djs6 + dis6*djs2)
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
      sderv2(4,3) = sderv2(4,3) + d2ijs*(dis3*djs4 + dis4*djs3)
      sderv2(5,3) = sderv2(5,3) + d2ijs*(dis3*djs5 + dis5*djs3)
      sderv2(6,3) = sderv2(6,3) + d2ijs*(dis3*djs6 + dis6*djs3)
      sderv2(4,4) = sderv2(4,4) + 2.0_dp*d2ijs*dis4*djs4
      sderv2(5,4) = sderv2(5,4) + d2ijs*(dis4*djs5 + dis5*djs4)
      sderv2(6,4) = sderv2(6,4) + d2ijs*(dis4*djs6 + dis6*djs4)
      sderv2(5,5) = sderv2(5,5) + 2.0_dp*d2ijs*dis5*djs5
      sderv2(6,5) = sderv2(6,5) + d2ijs*(dis5*djs6 + dis6*djs5)
      sderv2(6,6) = sderv2(6,6) + 2.0_dp*d2ijs*dis6*djs6
      if (abs(d2i2s).gt.1.0d-8) then
        sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
        sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
        sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
        sderv2(4,1) = sderv2(4,1) + d2i2s*dis1*dis4
        sderv2(5,1) = sderv2(5,1) + d2i2s*dis1*dis5
        sderv2(6,1) = sderv2(6,1) + d2i2s*dis1*dis6
        sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
        sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
        sderv2(4,2) = sderv2(4,2) + d2i2s*dis2*dis4
        sderv2(5,2) = sderv2(5,2) + d2i2s*dis2*dis5
        sderv2(6,2) = sderv2(6,2) + d2i2s*dis2*dis6
        sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
        sderv2(4,3) = sderv2(4,3) + d2i2s*dis3*dis4
        sderv2(5,3) = sderv2(5,3) + d2i2s*dis3*dis5
        sderv2(6,3) = sderv2(6,3) + d2i2s*dis3*dis6
        sderv2(4,4) = sderv2(4,4) + d2i2s*dis4*dis4
        sderv2(5,4) = sderv2(5,4) + d2i2s*dis4*dis5
        sderv2(6,4) = sderv2(6,4) + d2i2s*dis4*dis6
        sderv2(5,5) = sderv2(5,5) + d2i2s*dis5*dis5
        sderv2(6,5) = sderv2(6,5) + d2i2s*dis5*dis6
        sderv2(6,6) = sderv2(6,6) + d2i2s*dis6*dis6
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
        sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
        sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
        sderv2(4,1) = sderv2(4,1) + d2j2s*djs1*djs4
        sderv2(5,1) = sderv2(5,1) + d2j2s*djs1*djs5
        sderv2(6,1) = sderv2(6,1) + d2j2s*djs1*djs6
        sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
        sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
        sderv2(4,2) = sderv2(4,2) + d2j2s*djs2*djs4
        sderv2(5,2) = sderv2(5,2) + d2j2s*djs2*djs5
        sderv2(6,2) = sderv2(6,2) + d2j2s*djs2*djs6
        sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
        sderv2(4,3) = sderv2(4,3) + d2j2s*djs3*djs4
        sderv2(5,3) = sderv2(5,3) + d2j2s*djs3*djs5
        sderv2(6,3) = sderv2(6,3) + d2j2s*djs3*djs6
        sderv2(4,4) = sderv2(4,4) + d2j2s*djs4*djs4
        sderv2(5,4) = sderv2(5,4) + d2j2s*djs4*djs5
        sderv2(6,4) = sderv2(6,4) + d2j2s*djs4*djs6
        sderv2(5,5) = sderv2(5,5) + d2j2s*djs5*djs5
        sderv2(6,5) = sderv2(6,5) + d2j2s*djs5*djs6
        sderv2(6,6) = sderv2(6,6) + d2j2s*djs6*djs6
      endif
    elseif (ndim.eq.2) then
      d1si = dsi(1)
      d2si = dsi(2)
      d3si = dsi(3)
      d1sj = dsj(1)
      d2sj = dsj(2)
      d3sj = dsj(3)
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
      sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2
      sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1
      sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3
      sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2 + d2sj*djs2)
      sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3
      sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3)
!
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
      sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
      sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
      sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
      if (abs(d2i2s).gt.1.0d-8) then
        sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
        sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
        sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
        sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
        sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
        sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
        sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
        sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
        sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
        sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
        sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
      endif
    elseif (ndim.eq.1) then
      d1si = dsi(1)
      d1sj = dsj(1)
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
      if (abs(d2i2s).gt.1.0d-8) then
        sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
      endif
    endif
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charge')
#endif
!
  return
  end

  subroutine d2charged(iloc,i,j,nor,ix,iy,iz,jx,jy,jz,dei,dej,d1xi,d1yi,d1zi, &
                       d1xj,d1yj,d1zj,ds1i,ds1j,d2i2,d2ij,d2j2,d2self,dei0,dej0, &
                       lreal,lDoSelf,nbosptr)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from variable charge models.
!  This version is for calls from distributed memory routines.
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!  NB: loldd2q = .true. not supported
!
!   1/22 Created from d2charge
!   1/22 lNonIJQDeriv removed as not currently used
!   1/22 lopi/lopj removed from arguments as freezing is not allowed for this case
!   1/22 dedqc and d2edqc added for variable charge second derivatives
!  12/22 Optional argument, nbosptr added
!   9/23 Trap for shells added
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
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use m_strain,       only : strainddetds, straindet
  use optimisation
  use parallel
  use reaxFFdata,     only : reaxFFmu
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)           :: i       ! Atom i
  integer(i4), intent(in)           :: iloc    ! Atom i number on local node
  integer(i4), intent(in)           :: ix
  integer(i4), intent(in)           :: iy
  integer(i4), intent(in)           :: iz
  integer(i4), intent(in)           :: j
  integer(i4), intent(in)           :: jx
  integer(i4), intent(in)           :: jy
  integer(i4), intent(in)           :: jz
  integer(i4), intent(in)           :: nor
  integer(i4), intent(in), optional :: nbosptr(*)       ! Only required for reaxFF calls - points to reaxFF species from atoms
  logical,     intent(in)           :: lreal
  logical,     intent(in)           :: lDoSelf
  real(dp),    intent(in)           :: d1xi          ! d2E/dqi.dx
  real(dp),    intent(in)           :: d1yi          ! d2E/dqi.dy
  real(dp),    intent(in)           :: d1zi          ! d2E/dqi.dz
  real(dp),    intent(in)           :: d1xj          ! d2E/dqj.dx
  real(dp),    intent(in)           :: d1yj          ! d2E/dqj.dy
  real(dp),    intent(in)           :: d1zj          ! d2E/dqj.dz
  real(dp),    intent(in)           :: d2i2(*)       ! d2E/dqi2
  real(dp),    intent(in)           :: d2ij(*)       ! d2E/dqi.dqj
  real(dp),    intent(in)           :: d2j2(*)       ! d2E/dqj2
  real(dp),    intent(in)           :: d2self        ! Self term
  real(dp),    intent(in)           :: dei(*)        ! dE/dqi
  real(dp),    intent(in)           :: dej(*)        ! dE/dqj
  real(dp),    intent(in)           :: dei0          ! dE(selfterm)/dqi
  real(dp),    intent(in)           :: dej0          ! dE(selfterm)/dqj
  real(dp),    intent(in)           :: ds1i(*)       ! d2E/dqi.d(epsilon)
  real(dp),    intent(in)           :: ds1j(*)       ! d2E/dqj.d(epsilon)
!
!  Local variables
!
  integer(i4)                       :: ind
  integer(i4)                       :: indi
  integer(i4)                       :: indj
  integer(i4)                       :: ixf
  integer(i4)                       :: iyf
  integer(i4)                       :: izf
  integer(i4)                       :: iv
  integer(i4)                       :: k
  integer(i4)                       :: kk
  integer(i4)                       :: kl
  integer(i4)                       :: kx
  integer(i4)                       :: ky
  integer(i4)                       :: kz
  integer(i4)                       :: m
  integer(i4)                       :: mnxx
  integer(i4)                       :: mnxy
  integer(i4)                       :: mnxz
  integer(i4)                       :: mnyy
  integer(i4)                       :: mnyz
  integer(i4)                       :: mnzz
  integer(i4)                       :: mx
  integer(i4)                       :: my
  integer(i4)                       :: mz
  integer(i4)                       :: nqr
  integer(i4)                       :: nk
  real(dp)                          :: d1ix
  real(dp)                          :: d1iy
  real(dp)                          :: d1iz
  real(dp)                          :: d1jx
  real(dp)                          :: d1jy
  real(dp)                          :: d1jz
  real(dp)                          :: d2i2s
  real(dp)                          :: d2ijs
  real(dp)                          :: d2j2s
  real(dp)                          :: d2qk
  real(dp)                          :: d1si
  real(dp)                          :: d2si
  real(dp)                          :: d3si
  real(dp)                          :: d4si
  real(dp)                          :: d5si
  real(dp)                          :: d6si
  real(dp)                          :: d1sj
  real(dp)                          :: d2sj
  real(dp)                          :: d3sj
  real(dp)                          :: d4sj
  real(dp)                          :: d5sj
  real(dp)                          :: d6sj
  real(dp)                          :: deisum
  real(dp)                          :: dejsum
  real(dp)                          :: dikx
  real(dp)                          :: diky
  real(dp)                          :: dikz
  real(dp)                          :: dis1
  real(dp)                          :: dis2
  real(dp)                          :: dis3
  real(dp)                          :: dis4
  real(dp)                          :: dis5
  real(dp)                          :: dis6
  real(dp)                          :: djkx
  real(dp)                          :: djky
  real(dp)                          :: djkz
  real(dp)                          :: djs1
  real(dp)                          :: djs2
  real(dp)                          :: djs3
  real(dp)                          :: djs4
  real(dp)                          :: djs5
  real(dp)                          :: djs6
  real(dp)                          :: dkix
  real(dp)                          :: dkiy
  real(dp)                          :: dkiz
  real(dp)                          :: dkjx
  real(dp)                          :: dkjy
  real(dp)                          :: dkjz
  real(dp)                          :: dk1
  real(dp)                          :: dk2
  real(dp)                          :: dk3
  real(dp)                          :: dk4
  real(dp)                          :: dk5
  real(dp)                          :: dk6
  real(dp)                          :: dsi(6)
  real(dp)                          :: dsi2(6)
  real(dp)                          :: dsj(6)
  real(dp)                          :: dsj2(6)
  real(dp)                          :: g_cpu_time
  real(dp)                          :: time1
  real(dp)                          :: time2
  real(dp)                          :: ock
  real(dp)                          :: qlk
  real(dp)                          :: zetah0
  real(dp), parameter               :: tsqrt2pi = 0.797884560802866_dp
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charged')
#endif
!
  time1 = g_cpu_time()
! 
!  For ReaxFF check that nbosptr is present
!   
  if (lreaxFF) then
    if (.not.present(nbosptr)) then
      call outerror('d2charged called for reaxFF without passing nbosptr',0_i4)
      call stopnow('d2charged')
    endif
  endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  indi = 3*(i - 1)
  indj = 3*(j - 1)
  ixf = indi + 1
  iyf = indi + 2
  izf = indi + 3
!
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d2ijs = d2ijs + d2ij(iv)
    d2i2s = d2i2s + d2i2(iv)
    d2j2s = d2j2s + d2j2(iv)
  enddo
  d1ix = d1xi
  d1iy = d1yi
  d1iz = d1zi
  d1jx = d1xj
  d1jy = d1yj
  d1jz = d1zj
  if (leem) then
    if (nat(i).gt.maxele) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    endif
  endif
!******************
!  New algorithm  *
!******************
  d2edqdq(j,iloc) = d2edqdq(j,iloc) + d2ijs
  d2edq2(i) = d2edq2(i) + d2i2s
!
!  Add d1ix / d1iy / d1iz to sums for globalisation later
!
  dedqc(1,i) = dedqc(1,i) + d1ix
  dedqc(2,i) = dedqc(2,i) + d1iy
  dedqc(3,i) = dedqc(3,i) + d1iz
!
!  Add d1jx / d1iy / d1jz to sums for globalisation later
!
  d2edqc(ixf,j) = d2edqc(ixf,j) + d1jx
  d2edqc(iyf,j) = d2edqc(iyf,j) + d1jy
  d2edqc(izf,j) = d2edqc(izf,j) + d1jz
!
  if (lstr) then
    if (ndim.eq.3) then
      dis1 = dqds(1,i)
      dis2 = dqds(2,i)
      dis3 = dqds(3,i)
      dis4 = dqds(4,i)
      dis5 = dqds(5,i)
      dis6 = dqds(6,i)
      djs1 = dqds(1,j)
      djs2 = dqds(2,j)
      djs3 = dqds(3,j)
      djs4 = dqds(4,j)
      djs5 = dqds(5,j)
      djs6 = dqds(6,j)
    elseif (ndim.eq.2) then
      dis1 = dqds(1,i)
      dis2 = dqds(2,i)
      dis3 = dqds(3,i)
      djs1 = dqds(1,j)
      djs2 = dqds(2,j)
      djs3 = dqds(3,j)
    elseif (ndim.eq.1) then
      dis1 = dqds(1,i)
      djs1 = dqds(1,j)
    endif
    if (lreal) then
!
!  Real space
!
      do kl = 1,nstrains
        dsi(kl) = ds1i(kl)
        dsj(kl) = ds1j(kl)
        dsi2(kl) = dsi(kl)
        dsj2(kl) = dsj(kl)
      enddo
    else
!
!  Reciprocal space
!
      do kl = 1,nstrains
        dsi(kl) = ds1i(kl)
        dsj(kl) = ds1j(kl)
        dsi2(kl) = dsi(kl)
        dsj2(kl) = dsj(kl)
      enddo
      deisum = 0.0_dp
      dejsum = 0.0_dp
      do iv = 1,nor
        deisum = deisum + dei(iv)
        dejsum = dejsum + dej(iv)
      enddo
!
      if (ndim.eq.3) then
        if (lfinitestrain) then
          do kl = 1,3
            dsi(kl) = dsi(kl) - deisum*strainddetds(kl)*straindet
            dsj(kl) = dsj(kl) - dejsum*strainddetds(kl)*straindet
          enddo
        else
          do kl = 1,3
            dsi(kl) = dsi(kl) - deisum
            dsj(kl) = dsj(kl) - dejsum
          enddo
        endif
      elseif (ndim.eq.2) then
        if (lfinitestrain) then
          do kl = 1,3
            dsi(kl) = dsi(kl) - deisum*strainddetds(kl)*straindet
            dsj(kl) = dsj(kl) - dejsum*strainddetds(kl)*straindet
          enddo
        else
          do kl = 1,2
            dsi(kl) = dsi(kl) - deisum
            dsj(kl) = dsj(kl) - dejsum
          enddo
        endif
      elseif (ndim.eq.1) then
        dsi(1) = dsi(1) - deisum
        dsj(1) = dsj(1) - dejsum
      endif
      if (leem.or.lgfnff) then
        do kl = 1,nstrains
          dsi2(kl) = dsi(kl)
          dsj2(kl) = dsj(kl)
        enddo
      endif
    endif
!
!  Add contributions to arrays for globalisation and use later
!
    if (i.eq.j) then
      do kl = 1,nstrains
        ds2g(kl,i)   = ds2g(kl,i)   + 2.0_dp*dsi2(kl)
        ds2gs(kl,i)  = ds2gs(kl,i)  + 2.0_dp*d2ijs*dqds(kl,j)
        d2s2gs(kl,i) = d2s2gs(kl,i) + d2i2s*dqds(kl,i)
      enddo
    else
      do kl = 1,nstrains
        ds2g(kl,i)   = ds2g(kl,i)   + dsi2(kl)
        ds2gs(kl,i)  = ds2gs(kl,i)  + d2ijs*dqds(kl,j)
        d2s2gs(kl,i) = d2s2gs(kl,i) + d2i2s*dqds(kl,i)
      enddo
    endif
  endif
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff.and..not.lreaxFF) then
    deisum = dei0
    dejsum = dej0
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!   
!  Strain - strain contribution
!   
    if (lstr.and.ndim.gt.0) then
      ind = 0
      do kk = 1,nstrains
        do kl = 1,kk-1
          ind = ind + 1
          sderv2(kl,kk) = sderv2(kl,kk) + deisum*d2qds2(ind,iloc)
          sderv2(kk,kl) = sderv2(kk,kl) + deisum*d2qds2(ind,iloc)
        enddo
        ind = ind + 1
        sderv2(kk,kk) = sderv2(kk,kk) + deisum*d2qds2(ind,iloc)
      enddo
    endif
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(iloc)
      k = nqatomptr(m,iloc)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,iloc)
      derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,iloc)
      derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,iloc)
      derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,iloc)
      derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,iloc)
      derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,iloc)
      derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,iloc)
      derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,iloc)
      derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,iloc)
!
!  Mixed strain contribution
!
      if (lstr.and.ndim.gt.0) then
        do kk = 1,nstrains
! DEBUG - k terms not computed and so may not be correct
! DEBUG - bond order charges not currently enabled and so needs checking in future
          derv3(ix,kk) = derv3(ix,kk) - deisum*d2qdxyzs(kk,mx,iloc)
          derv3(iy,kk) = derv3(iy,kk) - deisum*d2qdxyzs(kk,my,iloc)
          derv3(iz,kk) = derv3(iz,kk) - deisum*d2qdxyzs(kk,mz,iloc)
        enddo
      endif
    enddo
  endif
!
!  Loop over atoms to add ij-k contributions
!
  kx = - 2
  ky = - 1
  kz =   0
  do k = 1,numat
    kx = kx + 3
    ky = ky + 3
    kz = kz + 3
!
    dikx = dqdxyz(kx,i)
    diky = dqdxyz(ky,i)
    dikz = dqdxyz(kz,i)
    djkx = dqdxyz(kx,j)
    djky = dqdxyz(ky,j)
    djkz = dqdxyz(kz,j)
    if (lreal.and.lDoSelf) then
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
      nk = nat(k)
      ock = occuf(k)
      if (leem.and.nk.le.maxele) then
        if (lmultiqrange.and.neemrptr(k).ne.0) then
          nqr = nqrnow(neemrptr(k))
        else
          nqr = 1
        endif
        if (lqeq) then
          if (nk.ne.1) then
            d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
          else
            qlk = qf(k)
            zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
            d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp + 2.0_dp*qlk/zetah0)
          endif
        else
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        endif
      elseif (lgfnff) then
        d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
      elseif (lreaxFF) then
        d2qk = 2.0_dp*ock*reaxFFmu(nbosptr(k))
      else
        d2qk = 0.0_dp
      endif
      dkix = dqdxyz(indi+1,k)
      dkiy = dqdxyz(indi+2,k)
      dkiz = dqdxyz(indi+3,k)
      dkjx = dqdxyz(indj+1,k)
      dkjy = dqdxyz(indj+2,k)
      dkjz = dqdxyz(indj+3,k)
      derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
      derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
      derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
      derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
      derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
      derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
      derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
      derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
      derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
      if (lstr.and.i.eq.j) then
!
!  Mixed strain terms from self energy - only need to do once for each atom - hence check on i = j
!
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) + d2qk*dqds(kl,k)*dkix
          derv3(iy,kl) = derv3(iy,kl) + d2qk*dqds(kl,k)*dkiy
          derv3(iz,kl) = derv3(iz,kl) + d2qk*dqds(kl,k)*dkiz
        enddo
!
!  Strain-strain derivatives - only need to do if i=j=k
!
        if (i.eq.k) then
          if (ndim.eq.3) then
            dk1 = dqds(1,k)
            dk2 = dqds(2,k)
            dk3 = dqds(3,k)
            dk4 = dqds(4,k)
            dk5 = dqds(5,k)
            dk6 = dqds(6,k)
            sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
            sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
            sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
            sderv2(4,1) = sderv2(4,1) + d2qk*dk4*dk1
            sderv2(5,1) = sderv2(5,1) + d2qk*dk5*dk1
            sderv2(6,1) = sderv2(6,1) + d2qk*dk6*dk1
            sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
            sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
            sderv2(4,2) = sderv2(4,2) + d2qk*dk4*dk2
            sderv2(5,2) = sderv2(5,2) + d2qk*dk5*dk2
            sderv2(6,2) = sderv2(6,2) + d2qk*dk6*dk2
            sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
            sderv2(4,3) = sderv2(4,3) + d2qk*dk4*dk3
            sderv2(5,3) = sderv2(5,3) + d2qk*dk5*dk3
            sderv2(6,3) = sderv2(6,3) + d2qk*dk6*dk3
            sderv2(4,4) = sderv2(4,4) + d2qk*dk4*dk4
            sderv2(5,4) = sderv2(5,4) + d2qk*dk5*dk4
            sderv2(6,4) = sderv2(6,4) + d2qk*dk6*dk4
            sderv2(5,5) = sderv2(5,5) + d2qk*dk5*dk5
            sderv2(6,5) = sderv2(6,5) + d2qk*dk6*dk5
            sderv2(6,6) = sderv2(6,6) + d2qk*dk6*dk6
          elseif (ndim.eq.2) then
            dk1 = dqds(1,k)
            dk2 = dqds(2,k)
            dk3 = dqds(3,k)
            sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
            sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
            sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
            sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
            sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
            sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
          elseif (ndim.eq.1) then
            dk1 = dqds(1,k)
            sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
          endif
        endif
      endif
    endif
    if (lstr) then
!
!  Mix strain-internal contribution
!
      if (abs(d2i2s).gt.1.0d-8) then
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) - d2i2s*dqds(kl,i)*dikx
          derv3(iy,kl) = derv3(iy,kl) - d2i2s*dqds(kl,i)*diky
          derv3(iz,kl) = derv3(iz,kl) - d2i2s*dqds(kl,i)*dikz
        enddo
      endif
    endif
!
!  End of loop over k
!
  enddo
!
!  Strain - charge derivatives terms
!
  if (lstr) then
!
!  Mixed strain-internal term
!
!  d2E/d(alpha).dq x dq/d(strain)
!
!  NB: Exclude i = j below as it doesn't cancel as it should
!
    if (i.ne.j) then
      do kl = 1,nstrains
        derv3(ix,kl) = derv3(ix,kl) - d1ix*dqds(kl,i)
        derv3(iy,kl) = derv3(iy,kl) - d1iy*dqds(kl,i)
        derv3(iz,kl) = derv3(iz,kl) - d1iz*dqds(kl,i)
        derv3(ix,kl) = derv3(ix,kl) - d1jx*dqds(kl,j)
        derv3(iy,kl) = derv3(iy,kl) - d1jy*dqds(kl,j)
        derv3(iz,kl) = derv3(iz,kl) - d1jz*dqds(kl,j)
      enddo
    endif
!
!  Strain-strain terms for charge derivatives
!
!  NB: Only add if j is less than or equal to i to avoid double counting
!
    if (j.le.i) then
      if (ndim.eq.3) then
        d1si = dsi(1)
        d2si = dsi(2)
        d3si = dsi(3)
        d4si = dsi(4)
        d5si = dsi(5)
        d6si = dsi(6)
        d1sj = dsj(1)
        d2sj = dsj(2)
        d3sj = dsj(3)
        d4sj = dsj(4)
        d5sj = dsj(5)
        d6sj = dsj(6)
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
        sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2
        sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1
        sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3
        sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1
        sderv2(4,1) = sderv2(4,1) + d1si*dis4 + d1sj*djs4
        sderv2(4,1) = sderv2(4,1) + d4si*dis1 + d4sj*djs1
        sderv2(5,1) = sderv2(5,1) + d1si*dis5 + d1sj*djs5
        sderv2(5,1) = sderv2(5,1) + d5si*dis1 + d5sj*djs1
        sderv2(6,1) = sderv2(6,1) + d1si*dis6 + d1sj*djs6
        sderv2(6,1) = sderv2(6,1) + d6si*dis1 + d6sj*djs1
        sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2+ d2sj*djs2)
        sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3
        sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2
        sderv2(4,2) = sderv2(4,2) + d2si*dis4 + d2sj*djs4
        sderv2(4,2) = sderv2(4,2) + d4si*dis2 + d4sj*djs2
        sderv2(5,2) = sderv2(5,2) + d2si*dis5 + d2sj*djs5
        sderv2(5,2) = sderv2(5,2) + d5si*dis2 + d5sj*djs2
        sderv2(6,2) = sderv2(6,2) + d2si*dis6 + d2sj*djs6
        sderv2(6,2) = sderv2(6,2) + d6si*dis2 + d6sj*djs2
        sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3)
        sderv2(4,3) = sderv2(4,3) + d3si*dis4 + d3sj*djs4
        sderv2(4,3) = sderv2(4,3) + d4si*dis3 + d4sj*djs3
        sderv2(5,3) = sderv2(5,3) + d3si*dis5 + d3sj*djs5
        sderv2(5,3) = sderv2(5,3) + d5si*dis3 + d5sj*djs3
        sderv2(6,3) = sderv2(6,3) + d3si*dis6 + d3sj*djs6
        sderv2(6,3) = sderv2(6,3) + d6si*dis3 + d6sj*djs3
        sderv2(4,4) = sderv2(4,4) + 2.0_dp*(d4si*dis4 + d4sj*djs4)
        sderv2(5,4) = sderv2(5,4) + d4si*dis5 + d4sj*djs5
        sderv2(5,4) = sderv2(5,4) + d5si*dis4 + d5sj*djs4
        sderv2(6,4) = sderv2(6,4) + d4si*dis6 + d4sj*djs6
        sderv2(6,4) = sderv2(6,4) + d6si*dis4 + d6sj*djs4
        sderv2(5,5) = sderv2(5,5) + 2.0_dp*(d5si*dis5 + d5sj*djs5)
        sderv2(6,5) = sderv2(6,5) + d5si*dis6 + d5sj*djs6
        sderv2(6,5) = sderv2(6,5) + d6si*dis5 + d6sj*djs5
        sderv2(6,6) = sderv2(6,6) + 2.0_dp*(d6si*dis6 + d6sj*djs6)
!
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
        sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
        sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
        sderv2(4,1) = sderv2(4,1) + d2ijs*(dis1*djs4 + dis4*djs1)
        sderv2(5,1) = sderv2(5,1) + d2ijs*(dis1*djs5 + dis5*djs1)
        sderv2(6,1) = sderv2(6,1) + d2ijs*(dis1*djs6 + dis6*djs1)
        sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
        sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
        sderv2(4,2) = sderv2(4,2) + d2ijs*(dis2*djs4 + dis4*djs2)
        sderv2(5,2) = sderv2(5,2) + d2ijs*(dis2*djs5 + dis5*djs2)
        sderv2(6,2) = sderv2(6,2) + d2ijs*(dis2*djs6 + dis6*djs2)
        sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
        sderv2(4,3) = sderv2(4,3) + d2ijs*(dis3*djs4 + dis4*djs3)
        sderv2(5,3) = sderv2(5,3) + d2ijs*(dis3*djs5 + dis5*djs3)
        sderv2(6,3) = sderv2(6,3) + d2ijs*(dis3*djs6 + dis6*djs3)
        sderv2(4,4) = sderv2(4,4) + 2.0_dp*d2ijs*dis4*djs4
        sderv2(5,4) = sderv2(5,4) + d2ijs*(dis4*djs5 + dis5*djs4)
        sderv2(6,4) = sderv2(6,4) + d2ijs*(dis4*djs6 + dis6*djs4)
        sderv2(5,5) = sderv2(5,5) + 2.0_dp*d2ijs*dis5*djs5
        sderv2(6,5) = sderv2(6,5) + d2ijs*(dis5*djs6 + dis6*djs5)
        sderv2(6,6) = sderv2(6,6) + 2.0_dp*d2ijs*dis6*djs6
        if (abs(d2i2s).gt.1.0d-8) then
          sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
          sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
          sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
          sderv2(4,1) = sderv2(4,1) + d2i2s*dis1*dis4
          sderv2(5,1) = sderv2(5,1) + d2i2s*dis1*dis5
          sderv2(6,1) = sderv2(6,1) + d2i2s*dis1*dis6
          sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
          sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
          sderv2(4,2) = sderv2(4,2) + d2i2s*dis2*dis4
          sderv2(5,2) = sderv2(5,2) + d2i2s*dis2*dis5
          sderv2(6,2) = sderv2(6,2) + d2i2s*dis2*dis6
          sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
          sderv2(4,3) = sderv2(4,3) + d2i2s*dis3*dis4
          sderv2(5,3) = sderv2(5,3) + d2i2s*dis3*dis5
          sderv2(6,3) = sderv2(6,3) + d2i2s*dis3*dis6
          sderv2(4,4) = sderv2(4,4) + d2i2s*dis4*dis4
          sderv2(5,4) = sderv2(5,4) + d2i2s*dis4*dis5
          sderv2(6,4) = sderv2(6,4) + d2i2s*dis4*dis6
          sderv2(5,5) = sderv2(5,5) + d2i2s*dis5*dis5
          sderv2(6,5) = sderv2(6,5) + d2i2s*dis5*dis6
          sderv2(6,6) = sderv2(6,6) + d2i2s*dis6*dis6
        endif
        if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
          sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
          sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
          sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
          sderv2(4,1) = sderv2(4,1) + d2j2s*djs1*djs4
          sderv2(5,1) = sderv2(5,1) + d2j2s*djs1*djs5
          sderv2(6,1) = sderv2(6,1) + d2j2s*djs1*djs6
          sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
          sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
          sderv2(4,2) = sderv2(4,2) + d2j2s*djs2*djs4
          sderv2(5,2) = sderv2(5,2) + d2j2s*djs2*djs5
          sderv2(6,2) = sderv2(6,2) + d2j2s*djs2*djs6
          sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
          sderv2(4,3) = sderv2(4,3) + d2j2s*djs3*djs4
          sderv2(5,3) = sderv2(5,3) + d2j2s*djs3*djs5
          sderv2(6,3) = sderv2(6,3) + d2j2s*djs3*djs6
          sderv2(4,4) = sderv2(4,4) + d2j2s*djs4*djs4
          sderv2(5,4) = sderv2(5,4) + d2j2s*djs4*djs5
          sderv2(6,4) = sderv2(6,4) + d2j2s*djs4*djs6
          sderv2(5,5) = sderv2(5,5) + d2j2s*djs5*djs5
          sderv2(6,5) = sderv2(6,5) + d2j2s*djs5*djs6
          sderv2(6,6) = sderv2(6,6) + d2j2s*djs6*djs6
        endif
      elseif (ndim.eq.2) then
        d1si = dsi(1)
        d2si = dsi(2)
        d3si = dsi(3)
        d1sj = dsj(1)
        d2sj = dsj(2)
        d3sj = dsj(3)
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
        sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2
        sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1
        sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3
        sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1
        sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2 + d2sj*djs2)
        sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3
        sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2
        sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3)
!
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
        sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
        sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
        sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
        sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
        sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
        if (abs(d2i2s).gt.1.0d-8) then
          sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
          sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
          sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
          sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
          sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
          sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
        endif
        if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
          sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
          sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
          sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
          sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
          sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
          sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
        endif
      elseif (ndim.eq.1) then
        d1si = dsi(1)
        d1sj = dsj(1)
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
        sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
        if (abs(d2i2s).gt.1.0d-8) then
          sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
        endif
        if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
          sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
        endif
      endif
    endif
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charged')
#endif
!
  return
  end

  subroutine d2charge2D(i,j,nor,ix,iy,iz,jx,jy,jz,lopi,lopj,dei,dej, &
                        d1xi,d1yi,d1zi,d1xj,d1yj,d1zj,ds1i,ds1j,d2i2,d2ij,d2j2, &
                        d2self,d2trm1dij,d2trm1diz,d2trm1djz,dtrm1di,dtrm1dj)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from EEM/QEq for the 2-D case. Only
!  needed for reciprocal space part.
!
!   3/01 Created from d2charge
!  11/02 Second strain derivatives corrected
!  11/02 Error in derv3 calculation for QEq(H) case corrected
!   9/04 Order of terms in dqds switched
!   9/04 Modified to generalise to charge derivatives other than from EEM
!   9/04 Contribution from d2q/dalpha.dbeta added
!   7/05 lReverse set to true for k.le.(i.or.j)
!  11/07 Unused variables removed
!   2/18 Trace added
!   1/19 Modified to that derivative terms are passed in rather than 
!        constructed locally, as per d2charge following strain modifications
!   1/19 Correction to volume related term at finite strain
!   1/19 Error in lReverse flag setting corrected
!   2/19 Corrections to alpha/beta order for derv2
!   7/19 lReverse and alpha/beta changes reversed
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   2/21 Correction for GFNFF added 
!   6/21 lexact_d2q added
!   6/21 Exact d2q terms rearranged to speed up
!   7/21 Timing added
!   7/21 lexact_d2q replaced by loldd2q
!   7/21 New algorithm to avoid 4th power scaling added
!  10/21 lgfnff moved to control
!   9/23 Trapping of shells added
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
  use control
  use current
  use derivatives
  use element
  use m_strain,       only : strainddetds, straindet
  use optimisation
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: ix
  integer(i4), intent(in)  :: iy
  integer(i4), intent(in)  :: iz
  integer(i4), intent(in)  :: j
  integer(i4), intent(in)  :: jx
  integer(i4), intent(in)  :: jy
  integer(i4), intent(in)  :: jz
  integer(i4), intent(in)  :: nor
  logical,     intent(in)  :: lopi
  logical,     intent(in)  :: lopj
  real(dp),    intent(in)  :: d1xi          ! d2E/dqi.dx
  real(dp),    intent(in)  :: d1yi          ! d2E/dqi.dy
  real(dp),    intent(in)  :: d1zi          ! d2E/dqi.dz
  real(dp),    intent(in)  :: d1xj          ! d2E/dqj.dx
  real(dp),    intent(in)  :: d1yj          ! d2E/dqj.dy
  real(dp),    intent(in)  :: d1zj          ! d2E/dqj.dz
  real(dp),    intent(in)  :: d2i2(*)
  real(dp),    intent(in)  :: d2ij(*)
  real(dp),    intent(in)  :: d2j2(*)
  real(dp),    intent(in)  :: d2self
  real(dp),    intent(in)  :: d2trm1dij
  real(dp),    intent(in)  :: d2trm1diz
  real(dp),    intent(in)  :: d2trm1djz
  real(dp),    intent(in)  :: dei(*)
  real(dp),    intent(in)  :: dej(*)
  real(dp),    intent(in)  :: ds1i(*)       ! d2E/dqi.d(epsilon)
  real(dp),    intent(in)  :: ds1j(*)       ! d2E/dqj.d(epsilon)
  real(dp),    intent(in)  :: dtrm1di
  real(dp),    intent(in)  :: dtrm1dj
!
!  Local variables
!
  integer(i4)              :: ind
  integer(i4)              :: indk
  integer(i4)              :: indl
  integer(i4)              :: iv
  integer(i4)              :: ixl
  integer(i4)              :: iyl
  integer(i4)              :: izl
  integer(i4)              :: jxl
  integer(i4)              :: jyl
  integer(i4)              :: jzl
  integer(i4)              :: k
  integer(i4)              :: kk
  integer(i4)              :: kl
  integer(i4)              :: kx
  integer(i4)              :: ky
  integer(i4)              :: kz
  integer(i4)              :: kxl
  integer(i4)              :: kyl
  integer(i4)              :: kzl
  integer(i4)              :: l
  integer(i4)              :: lx
  integer(i4)              :: ly
  integer(i4)              :: lz
  integer(i4)              :: lxl
  integer(i4)              :: lyl
  integer(i4)              :: lzl
  integer(i4)              :: m
  integer(i4)              :: mnxx
  integer(i4)              :: mnxy
  integer(i4)              :: mnxz
  integer(i4)              :: mnyx
  integer(i4)              :: mnyy
  integer(i4)              :: mnyz
  integer(i4)              :: mnzx
  integer(i4)              :: mnzy
  integer(i4)              :: mnzz
  integer(i4)              :: mx
  integer(i4)              :: my
  integer(i4)              :: mz
  integer(i4)              :: n
  integer(i4)              :: nx
  logical                  :: lNonIJQDeriv
  logical                  :: lopk
  logical                  :: lopl
  logical                  :: lReverse
  real(dp)                 :: d1ix
  real(dp)                 :: d1iy
  real(dp)                 :: d1iz
  real(dp)                 :: d1jx
  real(dp)                 :: d1jy
  real(dp)                 :: d1jz
  real(dp)                 :: d2i2s
  real(dp)                 :: d2ijs
  real(dp)                 :: d2j2s
  real(dp)                 :: d1si   
  real(dp)                 :: d2si   
  real(dp)                 :: d3si
  real(dp)                 :: d1sj   
  real(dp)                 :: d2sj   
  real(dp)                 :: d3sj
  real(dp)                 :: deisum
  real(dp)                 :: dejsum
  real(dp)                 :: dikx   
  real(dp)                 :: diky
  real(dp)                 :: dikz
  real(dp)                 :: d2ikx
  real(dp)                 :: d2iky
  real(dp)                 :: d2ikz
  real(dp)                 :: dilx
  real(dp)                 :: dily
  real(dp)                 :: dilz 
  real(dp)                 :: dis1 
  real(dp)                 :: dis2 
  real(dp)                 :: dis3
  real(dp)                 :: djkx
  real(dp)                 :: djky
  real(dp)                 :: djkz
  real(dp)                 :: d2jkx
  real(dp)                 :: d2jky
  real(dp)                 :: d2jkz
  real(dp)                 :: djlx   
  real(dp)                 :: djly   
  real(dp)                 :: djlz  
  real(dp)                 :: djs1  
  real(dp)                 :: djs2  
  real(dp)                 :: djs3
  real(dp)                 :: dsi(3)
  real(dp)                 :: dsi2(3)
  real(dp)                 :: dsj(3)
  real(dp)                 :: dsj2(3)
  real(dp)                 :: g_cpu_time
  real(dp)                 :: time1
  real(dp)                 :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charge2D')
#endif
!
  time1 = g_cpu_time()
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d1ix = 0.0_dp
  d1iy = 0.0_dp
  d1iz = 0.0_dp
  d1jx = 0.0_dp
  d1jy = 0.0_dp
  d1jz = 0.0_dp
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d2ijs = d2ijs + d2ij(iv)
    d2i2s = d2i2s + d2i2(iv)
    d2j2s = d2j2s + d2j2(iv)
  enddo
  d1ix = d1xi
  d1iy = d1yi
  d1iz = d1zi
  d1jx = d1xj
  d1jy = d1yj
  d1jz = d1zj
!
  d1iz = d1iz + d2trm1diz
  d1jz = d1jz + d2trm1djz
  d2ijs = d2ijs + d2trm1dij
!
  if (leem) then
    if (nat(i).gt.maxele) then
      d2ijs = 0.0_dp  
      d2i2s = 0.0_dp  
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2ijs = 0.0_dp  
      d2i2s = 0.0_dp  
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp  
      d2j2s = 0.0_dp  
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp  
      d2j2s = 0.0_dp  
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    endif
  endif
!
  if (.not.loldd2q) then
!******************
!  New algorithm  *
!******************
    d2edqdq(j,i) = d2edqdq(j,i) + d2ijs
    d2edq2(i) = d2edq2(i) + d2i2s
    if (i.ne.j) then
      d2edqdq(i,j) = d2edqdq(i,j) + d2ijs
      d2edq2(j) = d2edq2(j) + d2j2s
    endif
  endif
!
  if (lstr) then
    dis1 = dqds(1,i)
    dis2 = dqds(2,i)
    dis3 = dqds(3,i)
    djs1 = dqds(1,j)
    djs2 = dqds(2,j)
    djs3 = dqds(3,j)
    do kk = 1,nstrains
      dsi(kk) = ds1i(kk)
      dsj(kk) = ds1j(kk)
      dsi2(kk) = dsi(kk)
      dsj2(kk) = dsj(kk)
    enddo
    if (lfinitestrain) then
      do kk = 1,nstrains
        do iv = 1,nor
          dsi(kk) = dsi(kk) - dei(iv)*strainddetds(kk)*straindet
          dsj(kk) = dsj(kk) - dej(iv)*strainddetds(kk)*straindet
        enddo
        dsi(kk) = dsi(kk) - dtrm1di*strainddetds(kk)*straindet
        dsj(kk) = dsj(kk) - dtrm1dj*strainddetds(kk)*straindet
      enddo
    else
      do iv = 1,nor
        dsi(1) = dsi(1) - dei(iv)
        dsj(1) = dsj(1) - dej(iv)
        dsi(2) = dsi(2) - dei(iv)
        dsj(2) = dsj(2) - dej(iv)
      enddo
      dsi(1) = dsi(1) - dtrm1di
      dsj(1) = dsj(1) - dtrm1dj
      dsi(2) = dsi(2) - dtrm1di
      dsj(2) = dsj(2) - dtrm1dj
    endif
    if (leem.or.lgfnff) then
      do kk = 1,nstrains
        dsi2(kk) = dsi(kk)
        dsj2(kk) = dsj(kk)
      enddo
    endif
  endif
!
!  The following flag causes the second derivatives of the charge distribution 
!  with respect to two atoms neither of which are the i/j of the energy term
!  to be calculated. At the moment this contribution is not used in any method
!  and so the flag is set to false for efficiency.
!
  lNonIJQDeriv = .false.
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff.and..not.lreaxFF) then
    deisum = dtrm1di
    dejsum = dtrm1dj
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(i)
      k = nqatomptr(m,i)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (i.gt.k) then
        derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,i)
        derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,i)
        derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,i)
        derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,i)
        derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,i)
        derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,i)
        derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,i)
        derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,i)
        derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,i)
      elseif (i.lt.k) then
        derv2(ix,kx) = derv2(ix,kx) + deisum*d2qdxyz2(mnxx,i)
        derv2(iy,kx) = derv2(iy,kx) + deisum*d2qdxyz2(mnxy,i)
        derv2(iz,kx) = derv2(iz,kx) + deisum*d2qdxyz2(mnxz,i)
        derv2(ix,ky) = derv2(ix,ky) + deisum*d2qdxyz2(mnxy,i)
        derv2(iy,ky) = derv2(iy,ky) + deisum*d2qdxyz2(mnyy,i)
        derv2(iz,ky) = derv2(iz,ky) + deisum*d2qdxyz2(mnyz,i)
        derv2(ix,kz) = derv2(ix,kz) + deisum*d2qdxyz2(mnxz,i)
        derv2(iy,kz) = derv2(iy,kz) + deisum*d2qdxyz2(mnyz,i)
        derv2(iz,kz) = derv2(iz,kz) + deisum*d2qdxyz2(mnzz,i)
      endif
!
!  Strain - strain contribution
!
      if (lstr.and.ndim.gt.0) then
        ind = 0
        do kk = 1,nstrains
          do kl = 1,kk-1
            ind = ind + 1
            sderv2(kk,kl) = sderv2(kk,kl) + deisum*d2qds2(ind,i)
          enddo
          ind = ind + 1
          sderv2(kk,kk) = sderv2(kk,kk) + deisum*d2qds2(ind,i)
!
          derv3(kx,kk) = derv3(kx,kk) + deisum*d2qdxyzs(kk,mx,i)
          derv3(ky,kk) = derv3(ky,kk) + deisum*d2qdxyzs(kk,my,i)
          derv3(kz,kk) = derv3(kz,kk) + deisum*d2qdxyzs(kk,mz,i)
          derv3(ix,kk) = derv3(ix,kk) - deisum*d2qdxyzs(kk,mx,i)
          derv3(iy,kk) = derv3(iy,kk) - deisum*d2qdxyzs(kk,my,i)
          derv3(iz,kk) = derv3(iz,kk) - deisum*d2qdxyzs(kk,mz,i)
        enddo
      endif
    enddo
!
    do m = 1,nqatoms(j)
      k = nqatomptr(m,j)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (j.gt.k) then
        derv2(kx,jx) = derv2(kx,jx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(ky,jx) = derv2(ky,jx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(kz,jx) = derv2(kz,jx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(kx,jy) = derv2(kx,jy) + dejsum*d2qdxyz2(mnxy,j)
        derv2(ky,jy) = derv2(ky,jy) + dejsum*d2qdxyz2(mnyy,j)
        derv2(kz,jy) = derv2(kz,jy) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kx,jz) = derv2(kx,jz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(ky,jz) = derv2(ky,jz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kz,jz) = derv2(kz,jz) + dejsum*d2qdxyz2(mnzz,j)
      elseif (j.lt.k) then
        derv2(jx,kx) = derv2(jx,kx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(jy,kx) = derv2(jy,kx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jz,kx) = derv2(jz,kx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jx,ky) = derv2(jx,ky) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jy,ky) = derv2(jy,ky) + dejsum*d2qdxyz2(mnyy,j)
        derv2(jz,ky) = derv2(jz,ky) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jx,kz) = derv2(jx,kz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jy,kz) = derv2(jy,kz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jz,kz) = derv2(jz,kz) + dejsum*d2qdxyz2(mnzz,j)
      endif
!
!  Strain - strain contribution
!
      if (lstr.and.ndim.gt.0) then
        ind = 0
        do kk = 1,nstrains
          do kl = 1,kk-1
            ind = ind + 1
            sderv2(kk,kl) = sderv2(kk,kl) + dejsum*d2qds2(ind,j)
          enddo
          ind = ind + 1 
          sderv2(kk,kk) = sderv2(kk,kk) + dejsum*d2qds2(ind,j)
!
          derv3(kx,kk) = derv3(kx,kk) + dejsum*d2qdxyzs(kk,mx,j)
          derv3(ky,kk) = derv3(ky,kk) + dejsum*d2qdxyzs(kk,my,j)
          derv3(kz,kk) = derv3(kz,kk) + dejsum*d2qdxyzs(kk,mz,j)
          derv3(jx,kk) = derv3(jx,kk) - dejsum*d2qdxyzs(kk,mx,j)
          derv3(jy,kk) = derv3(jy,kk) - dejsum*d2qdxyzs(kk,my,j)
          derv3(jz,kk) = derv3(jz,kk) - dejsum*d2qdxyzs(kk,mz,j)
        enddo
      endif
    enddo
!
!  Terms involving Q derivatives for non-i/j case
!
    if (lNonIJQDeriv) then
      do m = 2,nqatoms(i)
        k = nqatomptr(m,i)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,i)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,i)
          derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,i)
          derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,i)
          derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,i)
          derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,i)
          derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,i)
          derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,i)
          derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,i)
          derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,i)
        enddo
      enddo
!
      do m = 2,nqatoms(j)
        k = nqatomptr(m,j)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,j)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,j)
          derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,j)
          derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,j)
          derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,j)
          derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,j)
          derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,j)
          derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,j)
          derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,j)
          derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,j)
        enddo
      enddo
    endif
  endif
!
!  Loop over atoms to add ij-k contributions
!
!  Note indices must be set to cope with the fact that only
!  half of derv2 is constructed in upper triangle and the
!  effects of lopi/lopj/lopk if freezing is being used.
!
  kx = - 2
  ky = - 1
  kz =   0
  do k = 1,numat
    lopk = (.not.lfreeze.or.lopf(k))
    indk = 3*(k - 1)
    if (lopk) then
      kx = kx + 3
      ky = ky + 3
      kz = kz + 3
    endif
    dikx = dqdxyz(indk+1,i)
    diky = dqdxyz(indk+2,i)
    dikz = dqdxyz(indk+3,i)
    djkx = dqdxyz(indk+1,j)
    djky = dqdxyz(indk+2,j)
    djkz = dqdxyz(indk+3,j)
    if (i.ne.k.and.(lopi.or.lopk)) then
      if (lopi.and.lopk) then
        if (k.le.i) then
          ixl = ix
          iyl = iy
          izl = iz
          kxl = kx
          kyl = ky
          kzl = kz
          lReverse = .true.
        else
          ixl = kx
          iyl = ky
          izl = kz
          kxl = ix
          kyl = iy
          kzl = iz
          lReverse = .true.
        endif
      elseif (lopi) then
        ixl = ix
        iyl = iy
        izl = iz
        kxl = ix
        kyl = iy
        kzl = iz
        lReverse = .false.
      elseif (lopk) then
        ixl = kx
        iyl = ky
        izl = kz
        kxl = kx
        kyl = ky
        kzl = kz
        lReverse = .true.
      endif
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-k
!
      if (lReverse) then
        derv2(kxl,ixl) = derv2(kxl,ixl) - d1ix*dikx
        derv2(kyl,ixl) = derv2(kyl,ixl) - d1iy*dikx
        derv2(kzl,ixl) = derv2(kzl,ixl) - d1iz*dikx
        derv2(kxl,iyl) = derv2(kxl,iyl) - d1ix*diky
        derv2(kyl,iyl) = derv2(kyl,iyl) - d1iy*diky
        derv2(kzl,iyl) = derv2(kzl,iyl) - d1iz*diky
        derv2(kxl,izl) = derv2(kxl,izl) - d1ix*dikz
        derv2(kyl,izl) = derv2(kyl,izl) - d1iy*dikz
        derv2(kzl,izl) = derv2(kzl,izl) - d1iz*dikz
!
        derv2(kxl,ixl) = derv2(kxl,ixl) - d1jx*djkx
        derv2(kyl,ixl) = derv2(kyl,ixl) - d1jy*djkx
        derv2(kzl,ixl) = derv2(kzl,ixl) - d1jz*djkx
        derv2(kxl,iyl) = derv2(kxl,iyl) - d1jx*djky
        derv2(kyl,iyl) = derv2(kyl,iyl) - d1jy*djky
        derv2(kzl,iyl) = derv2(kzl,iyl) - d1jz*djky
        derv2(kxl,izl) = derv2(kxl,izl) - d1jx*djkz
        derv2(kyl,izl) = derv2(kyl,izl) - d1jy*djkz
        derv2(kzl,izl) = derv2(kzl,izl) - d1jz*djkz
      else
        derv2(kxl,ixl) = derv2(kxl,ixl) - d1ix*dikx
        derv2(kyl,ixl) = derv2(kyl,ixl) - d1ix*diky
        derv2(kzl,ixl) = derv2(kzl,ixl) - d1ix*dikz
        derv2(kxl,iyl) = derv2(kxl,iyl) - d1iy*dikx
        derv2(kyl,iyl) = derv2(kyl,iyl) - d1iy*diky
        derv2(kzl,iyl) = derv2(kzl,iyl) - d1iy*dikz
        derv2(kxl,izl) = derv2(kxl,izl) - d1iz*dikx
        derv2(kyl,izl) = derv2(kyl,izl) - d1iz*diky
        derv2(kzl,izl) = derv2(kzl,izl) - d1iz*dikz
!
        derv2(kxl,ixl) = derv2(kxl,ixl) - d1jx*djkx
        derv2(kyl,ixl) = derv2(kyl,ixl) - d1jx*djky
        derv2(kzl,ixl) = derv2(kzl,ixl) - d1jx*djkz
        derv2(kxl,iyl) = derv2(kxl,iyl) - d1jy*djkx
        derv2(kyl,iyl) = derv2(kyl,iyl) - d1jy*djky
        derv2(kzl,iyl) = derv2(kzl,iyl) - d1jy*djkz
        derv2(kxl,izl) = derv2(kxl,izl) - d1jz*djkx
        derv2(kyl,izl) = derv2(kyl,izl) - d1jz*djky
        derv2(kzl,izl) = derv2(kzl,izl) - d1jz*djkz
      endif
    endif
    if (j.ne.k.and.(lopj.or.lopk)) then
      if (lopj.and.lopk) then
        if (k.le.j) then
          jxl = jx
          jyl = jy
          jzl = jz
          kxl = kx
          kyl = ky
          kzl = kz
          lReverse = .true.
        else
          jxl = kx
          jyl = ky
          jzl = kz
          kxl = jx
          kyl = jy
          kzl = jz
          lReverse = .true.
        endif
      elseif (lopj) then
        jxl = jx
        jyl = jy
        jzl = jz
        kxl = jx
        kyl = jy
        kzl = jz
        lReverse = .false.
      elseif (lopk) then
        jxl = kx
        jyl = ky
        jzl = kz
        kxl = kx
        kyl = ky
        kzl = kz
        lReverse = .true.
      endif
!
!  d2E/(d(alpha).dq) x dq/d(beta) : j-k
!
      if (lReverse) then
        derv2(kxl,jxl) = derv2(kxl,jxl) + d1jx*djkx
        derv2(kyl,jxl) = derv2(kyl,jxl) + d1jy*djkx
        derv2(kzl,jxl) = derv2(kzl,jxl) + d1jz*djkx
        derv2(kxl,jyl) = derv2(kxl,jyl) + d1jx*djky
        derv2(kyl,jyl) = derv2(kyl,jyl) + d1jy*djky
        derv2(kzl,jyl) = derv2(kzl,jyl) + d1jz*djky
        derv2(kxl,jzl) = derv2(kxl,jzl) + d1jx*djkz
        derv2(kyl,jzl) = derv2(kyl,jzl) + d1jy*djkz
        derv2(kzl,jzl) = derv2(kzl,jzl) + d1jz*djkz
!
        derv2(kxl,jxl) = derv2(kxl,jxl) + d1ix*dikx
        derv2(kyl,jxl) = derv2(kyl,jxl) + d1iy*dikx
        derv2(kzl,jxl) = derv2(kzl,jxl) + d1iz*dikx
        derv2(kxl,jyl) = derv2(kxl,jyl) + d1ix*diky
        derv2(kyl,jyl) = derv2(kyl,jyl) + d1iy*diky
        derv2(kzl,jyl) = derv2(kzl,jyl) + d1iz*diky
        derv2(kxl,jzl) = derv2(kxl,jzl) + d1ix*dikz
        derv2(kyl,jzl) = derv2(kyl,jzl) + d1iy*dikz
        derv2(kzl,jzl) = derv2(kzl,jzl) + d1iz*dikz
      else
        derv2(kxl,jxl) = derv2(kxl,jxl) + d1jx*djkx
        derv2(kyl,jxl) = derv2(kyl,jxl) + d1jx*djky
        derv2(kzl,jxl) = derv2(kzl,jxl) + d1jx*djkz
        derv2(kxl,jyl) = derv2(kxl,jyl) + d1jy*djkx
        derv2(kyl,jyl) = derv2(kyl,jyl) + d1jy*djky
        derv2(kzl,jyl) = derv2(kzl,jyl) + d1jy*djkz
        derv2(kxl,jzl) = derv2(kxl,jzl) + d1jz*djkx
        derv2(kyl,jzl) = derv2(kyl,jzl) + d1jz*djky
        derv2(kzl,jzl) = derv2(kzl,jzl) + d1jz*djkz
!
        derv2(kxl,jxl) = derv2(kxl,jxl) + d1ix*dikx
        derv2(kyl,jxl) = derv2(kyl,jxl) + d1ix*diky
        derv2(kzl,jxl) = derv2(kzl,jxl) + d1ix*dikz
        derv2(kxl,jyl) = derv2(kxl,jyl) + d1iy*dikx
        derv2(kyl,jyl) = derv2(kyl,jyl) + d1iy*diky
        derv2(kzl,jyl) = derv2(kzl,jyl) + d1iy*dikz
        derv2(kxl,jzl) = derv2(kxl,jzl) + d1iz*dikx
        derv2(kyl,jzl) = derv2(kyl,jzl) + d1iz*diky
        derv2(kzl,jzl) = derv2(kzl,jzl) + d1iz*dikz
      endif
    endif
    if (lstr.and.lopk) then
!
!  Mix strain-internal contribution
!
!  d2E/d(strain).dq x dq/d(alpha)
!
      do kl = 1,nstrains
        derv3(kx,kl) = derv3(kx,kl) + dsi2(kl)*dikx + dsj2(kl)*djkx
        derv3(ky,kl) = derv3(ky,kl) + dsi2(kl)*diky + dsj2(kl)*djky
        derv3(kz,kl) = derv3(kz,kl) + dsi2(kl)*dikz + dsj2(kl)*djkz
      enddo
!
!  d2E/dq.dq x dq/d(alpha) x dq/d(strain)
!
      do kl = 1,nstrains
        derv3(kx,kl) = derv3(kx,kl) + d2ijs*(dqds(kl,i)*djkx + dqds(kl,j)*dikx)
        derv3(ky,kl) = derv3(ky,kl) + d2ijs*(dqds(kl,i)*djky + dqds(kl,j)*diky)
        derv3(kz,kl) = derv3(kz,kl) + d2ijs*(dqds(kl,i)*djkz + dqds(kl,j)*dikz)
      enddo
      if (abs(d2i2s).gt.1.0d-8) then
        do kl = 1,nstrains
          derv3(kx,kl) = derv3(kx,kl) + d2i2s*dqds(kl,i)*dikx
          derv3(ky,kl) = derv3(ky,kl) + d2i2s*dqds(kl,i)*diky
          derv3(kz,kl) = derv3(kz,kl) + d2i2s*dqds(kl,i)*dikz
          derv3(ix,kl) = derv3(ix,kl) - d2i2s*dqds(kl,i)*dikx
          derv3(iy,kl) = derv3(iy,kl) - d2i2s*dqds(kl,i)*diky
          derv3(iz,kl) = derv3(iz,kl) - d2i2s*dqds(kl,i)*dikz
        enddo
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        do kl = 1,nstrains
          derv3(kx,kl) = derv3(kx,kl) + d2j2s*dqds(kl,j)*djkx
          derv3(ky,kl) = derv3(ky,kl) + d2j2s*dqds(kl,j)*djky
          derv3(kz,kl) = derv3(kz,kl) + d2j2s*dqds(kl,j)*djkz
          derv3(ix,kl) = derv3(ix,kl) - d2j2s*dqds(kl,j)*djkx
          derv3(iy,kl) = derv3(iy,kl) - d2j2s*dqds(kl,j)*djky
          derv3(iz,kl) = derv3(iz,kl) - d2j2s*dqds(kl,j)*djkz
        enddo
      endif
    endif
    if (loldd2q) then
!******************
!  Old algorithm  *
!******************
!
!  Loop over atoms to add ij-kl correction
!
      lx = - 2
      ly = - 1
      lz =   0
!
!  Scale charge derivatives for k by d2ijs
!
      d2ikx = d2ijs*dikx
      d2iky = d2ijs*diky
      d2ikz = d2ijs*dikz
      d2jkx = d2ijs*djkx
      d2jky = d2ijs*djky
      d2jkz = d2ijs*djkz
!
      do l = 1,k-1
        lopl = (.not.lfreeze.or.lopf(l))
        indl = 3*(l-1)
        if (lopl) then
          lx = lx + 3
          ly = ly + 3
          lz = lz + 3
        endif
        if (lopk.or.lopl) then
          if (lopk.and.lopl) then
            kxl = kx
            kyl = ky
            kzl = kz
            lxl = lx
            lyl = ly
            lzl = lz
          elseif (lopk) then
            kxl = kx
            kyl = ky
            kzl = kz
            lxl = kx
            lyl = ky
            lzl = kz
          elseif (lopl) then
            kxl = lx
            kyl = ly
            kzl = lz
            lxl = lx
            lyl = ly
            lzl = lz
          endif
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
          dilx = dqdxyz(indl+1,i)
          dily = dqdxyz(indl+2,i)
          dilz = dqdxyz(indl+3,i)
          djlx = dqdxyz(indl+1,j)
          djly = dqdxyz(indl+2,j)
          djlz = dqdxyz(indl+3,j)
!
          derv2(lxl,kxl) = derv2(lxl,kxl) + dilx*d2jkx + djlx*d2ikx
          derv2(lyl,kxl) = derv2(lyl,kxl) + dilx*d2jky + djlx*d2iky
          derv2(lzl,kxl) = derv2(lzl,kxl) + dilx*d2jkz + djlx*d2ikz
          derv2(lxl,kyl) = derv2(lxl,kyl) + dily*d2jkx + djly*d2ikx
          derv2(lyl,kyl) = derv2(lyl,kyl) + dily*d2jky + djly*d2iky
          derv2(lzl,kyl) = derv2(lzl,kyl) + dily*d2jkz + djly*d2ikz
          derv2(lxl,kzl) = derv2(lxl,kzl) + dilz*d2jkx + djlz*d2ikx
          derv2(lyl,kzl) = derv2(lyl,kzl) + dilz*d2jky + djlz*d2iky
          derv2(lzl,kzl) = derv2(lzl,kzl) + dilz*d2jkz + djlz*d2ikz
!
!  d2Edi2/d2Edj2 - only non-zero for QEq/H
!
          if (abs(d2i2s).gt.1.0d-8) then
            derv2(lxl,kxl) = derv2(lxl,kxl) + d2i2s*dilx*dikx
            derv2(lyl,kxl) = derv2(lyl,kxl) + d2i2s*dilx*diky
            derv2(lzl,kxl) = derv2(lzl,kxl) + d2i2s*dilx*dikz
            derv2(lxl,kyl) = derv2(lxl,kyl) + d2i2s*dily*dikx
            derv2(lyl,kyl) = derv2(lyl,kyl) + d2i2s*dily*diky
            derv2(lzl,kyl) = derv2(lzl,kyl) + d2i2s*dily*dikz
            derv2(lxl,kzl) = derv2(lxl,kzl) + d2i2s*dilz*dikx
            derv2(lyl,kzl) = derv2(lyl,kzl) + d2i2s*dilz*diky
            derv2(lzl,kzl) = derv2(lzl,kzl) + d2i2s*dilz*dikz
          endif
          if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
            derv2(lxl,kxl) = derv2(lxl,kxl) + d2j2s*djlx*djkx
            derv2(lyl,kxl) = derv2(lyl,kxl) + d2j2s*djlx*djky
            derv2(lzl,kxl) = derv2(lzl,kxl) + d2j2s*djlx*djkz
            derv2(lxl,kyl) = derv2(lxl,kyl) + d2j2s*djly*djkx
            derv2(lyl,kyl) = derv2(lyl,kyl) + d2j2s*djly*djky
            derv2(lzl,kyl) = derv2(lzl,kyl) + d2j2s*djly*djkz
            derv2(lxl,kzl) = derv2(lxl,kzl) + d2j2s*djlz*djkx
            derv2(lyl,kzl) = derv2(lyl,kzl) + d2j2s*djlz*djky
            derv2(lzl,kzl) = derv2(lzl,kzl) + d2j2s*djlz*djkz
          endif
        endif
!
!  End of loop over l
!
      enddo
    endif
!
!  End of loop over k
!
  enddo
!
!  Strain - charge derivatives terms
!
  if (lstr) then
!
!  Mixed strain-internal term
!
!  d2E/d(alpha).dq x dq/d(strain)
!
    if (lopi) then
      do kl = 1,nstrains
        derv3(ix,kl) = derv3(ix,kl) - d1ix*dqds(kl,i)
        derv3(iy,kl) = derv3(iy,kl) - d1iy*dqds(kl,i)
        derv3(iz,kl) = derv3(iz,kl) - d1iz*dqds(kl,i)
        derv3(ix,kl) = derv3(ix,kl) - d1jx*dqds(kl,j)
        derv3(iy,kl) = derv3(iy,kl) - d1jy*dqds(kl,j)
        derv3(iz,kl) = derv3(iz,kl) - d1jz*dqds(kl,j)
      enddo
    endif
    if (lopj) then
      do kl = 1,nstrains
        derv3(jx,kl) = derv3(jx,kl) + d1jx*dqds(kl,j)
        derv3(jy,kl) = derv3(jy,kl) + d1jy*dqds(kl,j)
        derv3(jz,kl) = derv3(jz,kl) + d1jz*dqds(kl,j)
        derv3(jx,kl) = derv3(jx,kl) + d1ix*dqds(kl,i)
        derv3(jy,kl) = derv3(jy,kl) + d1iy*dqds(kl,i)
        derv3(jz,kl) = derv3(jz,kl) + d1iz*dqds(kl,i)
      enddo
    endif
!
!  Strain-strain terms for charge derivatives
!
    d1si = dsi(1)
    d2si = dsi(2)
    d3si = dsi(3)
    d1sj = dsj(1)
    d2sj = dsj(2)
    d3sj = dsj(3)
    sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
    sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2
    sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1
    sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3
    sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1
    sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2 + d2sj*djs2)
    sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3
    sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2
    sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3)
!
    sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
    sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
    sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
    sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
    sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
    sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
    if (abs(d2i2s).gt.1.0d-8) then
      sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
      sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
      sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
      sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
      sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
      sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
    endif
    if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
      sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
      sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
      sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
      sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
      sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
      sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
    endif
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charge2D')
#endif
!
  return
  end

  subroutine d2charge2Dd(iloc,i,j,nor,ix,iy,iz,jx,jy,jz,dei,dej, &
                         d1xi,d1yi,d1zi,d1xj,d1yj,d1zj,ds1i,ds1j,d2i2,d2ij,d2j2, &
                         d2self,d2trm1dij,d2trm1diz,d2trm1djz,dtrm1di,dtrm1dj)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from EEM/QEq for the 2-D case. Only
!  needed for reciprocal space part.
!  Parallel distributed memory version.
!
!   1/22 Created from d2charge2D and d2charged
!   9/23 Trapping of shells added
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
  use control
  use current
  use derivatives
  use element
  use m_strain,       only : strainddetds, straindet
  use optimisation
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: iloc
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: ix
  integer(i4), intent(in)  :: iy
  integer(i4), intent(in)  :: iz
  integer(i4), intent(in)  :: j
  integer(i4), intent(in)  :: jx
  integer(i4), intent(in)  :: jy
  integer(i4), intent(in)  :: jz
  integer(i4), intent(in)  :: nor
  real(dp),    intent(in)  :: d1xi          ! d2E/dqi.dx
  real(dp),    intent(in)  :: d1yi          ! d2E/dqi.dy
  real(dp),    intent(in)  :: d1zi          ! d2E/dqi.dz
  real(dp),    intent(in)  :: d1xj          ! d2E/dqj.dx
  real(dp),    intent(in)  :: d1yj          ! d2E/dqj.dy
  real(dp),    intent(in)  :: d1zj          ! d2E/dqj.dz
  real(dp),    intent(in)  :: d2i2(*)
  real(dp),    intent(in)  :: d2ij(*)
  real(dp),    intent(in)  :: d2j2(*)
  real(dp),    intent(in)  :: d2self
  real(dp),    intent(in)  :: d2trm1dij
  real(dp),    intent(in)  :: d2trm1diz
  real(dp),    intent(in)  :: d2trm1djz
  real(dp),    intent(in)  :: dei(*)
  real(dp),    intent(in)  :: dej(*)
  real(dp),    intent(in)  :: ds1i(*)       ! d2E/dqi.d(epsilon)
  real(dp),    intent(in)  :: ds1j(*)       ! d2E/dqj.d(epsilon)
  real(dp),    intent(in)  :: dtrm1di
  real(dp),    intent(in)  :: dtrm1dj
!
!  Local variables
!
  integer(i4)              :: ind
  integer(i4)              :: iv
  integer(i4)              :: ixf
  integer(i4)              :: iyf
  integer(i4)              :: izf
  integer(i4)              :: k
  integer(i4)              :: kk
  integer(i4)              :: kl
  integer(i4)              :: kx
  integer(i4)              :: ky
  integer(i4)              :: kz
  integer(i4)              :: m
  integer(i4)              :: mnxx
  integer(i4)              :: mnxy
  integer(i4)              :: mnxz
  integer(i4)              :: mnyy
  integer(i4)              :: mnyz
  integer(i4)              :: mnzz
  integer(i4)              :: mx
  integer(i4)              :: my
  integer(i4)              :: mz
  real(dp)                 :: d1ix
  real(dp)                 :: d1iy
  real(dp)                 :: d1iz
  real(dp)                 :: d1jx
  real(dp)                 :: d1jy
  real(dp)                 :: d1jz
  real(dp)                 :: d2i2s
  real(dp)                 :: d2ijs
  real(dp)                 :: d2j2s
  real(dp)                 :: d1si   
  real(dp)                 :: d2si   
  real(dp)                 :: d3si
  real(dp)                 :: d1sj   
  real(dp)                 :: d2sj   
  real(dp)                 :: d3sj
  real(dp)                 :: deisum
  real(dp)                 :: dejsum
  real(dp)                 :: dikx   
  real(dp)                 :: diky
  real(dp)                 :: dikz
  real(dp)                 :: dis1 
  real(dp)                 :: dis2 
  real(dp)                 :: dis3
  real(dp)                 :: djs1  
  real(dp)                 :: djs2  
  real(dp)                 :: djs3
  real(dp)                 :: dsi(3)
  real(dp)                 :: dsi2(3)
  real(dp)                 :: dsj(3)
  real(dp)                 :: dsj2(3)
  real(dp)                 :: g_cpu_time
  real(dp)                 :: time1
  real(dp)                 :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charge2Dd')
#endif
!
  time1 = g_cpu_time()
!
  ixf = 3*(i-1) + 1
  iyf = ixf + 1
  izf = ixf + 2
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d1ix = 0.0_dp
  d1iy = 0.0_dp
  d1iz = 0.0_dp
  d1jx = 0.0_dp
  d1jy = 0.0_dp
  d1jz = 0.0_dp
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d2ijs = d2ijs + d2ij(iv)
    d2i2s = d2i2s + d2i2(iv)
    d2j2s = d2j2s + d2j2(iv)
  enddo
  d1ix = d1xi
  d1iy = d1yi
  d1iz = d1zi
  d1jx = d1xj
  d1jy = d1yj
  d1jz = d1zj
!
  d1iz = d1iz + d2trm1diz
  d1jz = d1jz + d2trm1djz
  d2ijs = d2ijs + d2trm1dij
!
  if (leem) then
    if (nat(i).gt.maxele) then
      d2ijs = 0.0_dp  
      d2i2s = 0.0_dp  
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2ijs = 0.0_dp  
      d2i2s = 0.0_dp  
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp  
      d2j2s = 0.0_dp  
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp  
      d2j2s = 0.0_dp  
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    endif
  endif
!
!  Add d1ix / d1iy / d1iz to sums for globalisation later
!
  dedqc(1,i) = dedqc(1,i) + d1ix
  dedqc(2,i) = dedqc(2,i) + d1iy
  dedqc(3,i) = dedqc(3,i) + d1iz
!
!  Add d1jx / d1iy / d1jz to sums for globalisation later
!
  d2edqc(ixf,j) = d2edqc(ixf,j) + d1jx
  d2edqc(iyf,j) = d2edqc(iyf,j) + d1jy
  d2edqc(izf,j) = d2edqc(izf,j) + d1jz
!******************
!  New algorithm  *
!******************
  d2edqdq(j,iloc) = d2edqdq(j,iloc) + d2ijs
  d2edq2(i) = d2edq2(i) + d2i2s
!
  if (lstr) then
    dis1 = dqds(1,i)
    dis2 = dqds(2,i)
    dis3 = dqds(3,i)
    djs1 = dqds(1,j)
    djs2 = dqds(2,j)
    djs3 = dqds(3,j)
    do kk = 1,nstrains
      dsi(kk) = ds1i(kk)
      dsj(kk) = ds1j(kk)
      dsi2(kk) = dsi(kk)
      dsj2(kk) = dsj(kk)
    enddo
    if (lfinitestrain) then
      do kk = 1,nstrains
        do iv = 1,nor
          dsi(kk) = dsi(kk) - dei(iv)*strainddetds(kk)*straindet
          dsj(kk) = dsj(kk) - dej(iv)*strainddetds(kk)*straindet
        enddo
        dsi(kk) = dsi(kk) - dtrm1di*strainddetds(kk)*straindet
        dsj(kk) = dsj(kk) - dtrm1dj*strainddetds(kk)*straindet
      enddo
    else
      do iv = 1,nor
        dsi(1) = dsi(1) - dei(iv)
        dsj(1) = dsj(1) - dej(iv)
        dsi(2) = dsi(2) - dei(iv)
        dsj(2) = dsj(2) - dej(iv)
      enddo
      dsi(1) = dsi(1) - dtrm1di
      dsj(1) = dsj(1) - dtrm1dj
      dsi(2) = dsi(2) - dtrm1di
      dsj(2) = dsj(2) - dtrm1dj
    endif
    if (leem.or.lgfnff) then
      do kk = 1,nstrains
        dsi2(kk) = dsi(kk)
        dsj2(kk) = dsj(kk)
      enddo
    endif
!
!  Add contributions to arrays for globalisation and use later
!
    if (i.eq.j) then
      do kl = 1,nstrains
        ds2g(kl,i)   = ds2g(kl,i)   + 2.0_dp*dsi2(kl)
        ds2gs(kl,i)  = ds2gs(kl,i)  + 2.0_dp*d2ijs*dqds(kl,j)
        d2s2gs(kl,i) = d2s2gs(kl,i) + d2i2s*dqds(kl,i)
      enddo
    else
      do kl = 1,nstrains
        ds2g(kl,i)   = ds2g(kl,i)   + dsi2(kl)
        ds2gs(kl,i)  = ds2gs(kl,i)  + d2ijs*dqds(kl,j)
        d2s2gs(kl,i) = d2s2gs(kl,i) + d2i2s*dqds(kl,i)
      enddo
    endif
  endif
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff.and..not.lreaxFF) then
    deisum = dtrm1di
    dejsum = dtrm1dj
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(iloc)
      k = nqatomptr(m,iloc)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (i.gt.k) then
        derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,iloc)
        derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,iloc)
        derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,iloc)
      elseif (i.lt.k) then
        derv2(ix,kx) = derv2(ix,kx) + deisum*d2qdxyz2(mnxx,iloc)
        derv2(iy,kx) = derv2(iy,kx) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(iz,kx) = derv2(iz,kx) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(ix,ky) = derv2(ix,ky) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(iy,ky) = derv2(iy,ky) + deisum*d2qdxyz2(mnyy,iloc)
        derv2(iz,ky) = derv2(iz,ky) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(ix,kz) = derv2(ix,kz) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(iy,kz) = derv2(iy,kz) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(iz,kz) = derv2(iz,kz) + deisum*d2qdxyz2(mnzz,iloc)
      endif
!
!  Strain - strain contribution
!
      if (lstr.and.ndim.gt.0) then
        ind = 0
        do kk = 1,nstrains
          do kl = 1,kk-1
            ind = ind + 1
            sderv2(kk,kl) = sderv2(kk,kl) + deisum*d2qds2(ind,iloc)
          enddo
          ind = ind + 1
          sderv2(kk,kk) = sderv2(kk,kk) + deisum*d2qds2(ind,iloc)
!
          derv3(ix,kk) = derv3(ix,kk) - deisum*d2qdxyzs(kk,mx,iloc)
          derv3(iy,kk) = derv3(iy,kk) - deisum*d2qdxyzs(kk,my,iloc)
          derv3(iz,kk) = derv3(iz,kk) - deisum*d2qdxyzs(kk,mz,iloc)
        enddo
      endif
    enddo
  endif
!
!  Loop over atoms to add ij-k contributions
!
  kx = - 2
  ky = - 1
  kz =   0
  do k = 1,numat
    kx = kx + 3
    ky = ky + 3
    kz = kz + 3
    dikx = dqdxyz(kx,i)
    diky = dqdxyz(ky,i)
    dikz = dqdxyz(kz,i)
!
    if (lstr) then
      if (abs(d2i2s).gt.1.0d-8) then
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) - d2i2s*dqds(kl,i)*dikx
          derv3(iy,kl) = derv3(iy,kl) - d2i2s*dqds(kl,i)*diky
          derv3(iz,kl) = derv3(iz,kl) - d2i2s*dqds(kl,i)*dikz
        enddo
      endif
    endif
!
!  End of loop over k
!
  enddo
!
!  Strain - charge derivatives terms
!
  if (lstr) then
!
!  Mixed strain-internal term
!
!  d2E/d(alpha).dq x dq/d(strain)
!
    do kl = 1,nstrains
      derv3(ix,kl) = derv3(ix,kl) - d1ix*dqds(kl,i)
      derv3(iy,kl) = derv3(iy,kl) - d1iy*dqds(kl,i)
      derv3(iz,kl) = derv3(iz,kl) - d1iz*dqds(kl,i)
      derv3(ix,kl) = derv3(ix,kl) - d1jx*dqds(kl,j)
      derv3(iy,kl) = derv3(iy,kl) - d1jy*dqds(kl,j)
      derv3(iz,kl) = derv3(iz,kl) - d1jz*dqds(kl,j)
    enddo
!
!  Strain-strain terms for charge derivatives
!
!  NB: Only add if j is less than or equal to i to avoid double counting
!
    if (j.le.i) then
      d1si = dsi(1)
      d2si = dsi(2)
      d3si = dsi(3)
      d1sj = dsj(1)
      d2sj = dsj(2)
      d3sj = dsj(3)
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1)
      sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2
      sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1
      sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3
      sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2 + d2sj*djs2)
      sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3
      sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3)
!
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
      sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
      sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
      sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
!
      if (abs(d2i2s).gt.1.0d-8) then
        sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
        sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
        sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
        sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
        sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
        sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
        sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
        sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
        sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
        sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
        sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
      endif
    endif
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charge2Dd')
#endif
!
  return
  end

  subroutine d2charge3(i,j,k,nor,ix,iy,iz,jx,jy,jz,kx,ky,kz,lopi,lopj,lopk,dei,dej,dek, &
                       d1q,d1qs,d2q)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from variable charge models as a 
!  result of three body potentials.
!
!  10/04 Created from d2charge
!  10/04 Corrections for sderv2, as applied to d2charge applied here
!   7/05 lReverse set to true for k.le.(i.or.j)
!  11/07 Unused variables removed
!   2/18 Trace added
!   1/19 Finite strain corrections added
!   1/19 Error in lReverse flag setting corrected
!   2/19 Corrections to alpha/beta order for derv2
!   7/19 lReverse and alpha/beta changes reversed
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   7/21 QEq terms corrected 
!  12/21 pGFNFF modifications added
!
!  On entry:
!
!  d1q(3,3,3) = derivative of first derivative of i - j vector with respect to charge
!  d1qs(6,3)  = derivative of first strain derivative with respect to charge
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, December 2021
!
  use control
  use current
  use derivatives
  use element
  use optimisation
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: i
  integer(i4), intent(in)    :: ix
  integer(i4), intent(in)    :: iy
  integer(i4), intent(in)    :: iz
  integer(i4), intent(in)    :: j
  integer(i4), intent(in)    :: jx
  integer(i4), intent(in)    :: jy
  integer(i4), intent(in)    :: jz
  integer(i4), intent(in)    :: k
  integer(i4), intent(in)    :: kx
  integer(i4), intent(in)    :: ky
  integer(i4), intent(in)    :: kz
  integer(i4), intent(in)    :: nor
  logical,     intent(in)    :: lopi
  logical,     intent(in)    :: lopj
  logical,     intent(in)    :: lopk
  real(dp),    intent(in)    :: d1q(3,3,3)
  real(dp),    intent(in)    :: d1qs(6,3)
  real(dp),    intent(in)    :: d2q(6,*)
  real(dp),    intent(in)    :: dei(*)
  real(dp),    intent(in)    :: dej(*)
  real(dp),    intent(in)    :: dek(*)
!
!  Local variables
!
  integer(i4)                :: ind
  integer(i4)                :: indl
  integer(i4)                :: indm
  integer(i4)                :: iv
  integer(i4)                :: ixl
  integer(i4)                :: iyl
  integer(i4)                :: izl
  integer(i4)                :: jxl
  integer(i4)                :: jyl
  integer(i4)                :: jzl
  integer(i4)                :: kk
  integer(i4)                :: kl
  integer(i4)                :: kxl
  integer(i4)                :: kyl
  integer(i4)                :: kzl
  integer(i4)                :: l
  integer(i4)                :: lx
  integer(i4)                :: ly
  integer(i4)                :: lz
  integer(i4)                :: lxl
  integer(i4)                :: lyl
  integer(i4)                :: lzl
  integer(i4)                :: m
  integer(i4)                :: mnxx
  integer(i4)                :: mnxy
  integer(i4)                :: mnxz
  integer(i4)                :: mnyx
  integer(i4)                :: mnyy
  integer(i4)                :: mnyz
  integer(i4)                :: mnzx
  integer(i4)                :: mnzy
  integer(i4)                :: mnzz
  integer(i4)                :: mx
  integer(i4)                :: my
  integer(i4)                :: mz
  integer(i4)                :: mxl
  integer(i4)                :: myl
  integer(i4)                :: mzl
  integer(i4)                :: n
  integer(i4)                :: nx
  integer(i4)                :: p
  integer(i4)                :: px
  integer(i4)                :: py
  integer(i4)                :: pz
  logical                    :: lopl
  logical                    :: lopm
  logical                    :: lNonIJQDeriv
  logical                    :: lReverse
  real(dp)                   :: d1ql(3,3,3)
  real(dp)                   :: d1qsl(6,3)
  real(dp)                   :: d1si
  real(dp)                   :: d2si
  real(dp)                   :: d3si
  real(dp)                   :: d4si
  real(dp)                   :: d5si
  real(dp)                   :: d6si
  real(dp)                   :: d1sj
  real(dp)                   :: d2sj
  real(dp)                   :: d3sj
  real(dp)                   :: d4sj
  real(dp)                   :: d5sj
  real(dp)                   :: d6sj
  real(dp)                   :: d1sk
  real(dp)                   :: d2sk
  real(dp)                   :: d3sk
  real(dp)                   :: d4sk
  real(dp)                   :: d5sk
  real(dp)                   :: d6sk
  real(dp)                   :: d2i2s
  real(dp)                   :: d2ijs
  real(dp)                   :: d2iks
  real(dp)                   :: d2j2s
  real(dp)                   :: d2jks
  real(dp)                   :: d2k2s
  real(dp)                   :: deisum
  real(dp)                   :: dejsum
  real(dp)                   :: deksum
  real(dp)                   :: dilx
  real(dp)                   :: dily
  real(dp)                   :: dilz
  real(dp)                   :: dimx
  real(dp)                   :: dimy
  real(dp)                   :: dimz
  real(dp)                   :: dis1
  real(dp)                   :: dis2
  real(dp)                   :: dis3
  real(dp)                   :: dis4
  real(dp)                   :: dis5
  real(dp)                   :: dis6
  real(dp)                   :: djlx
  real(dp)                   :: djly
  real(dp)                   :: djlz
  real(dp)                   :: djmx
  real(dp)                   :: djmy
  real(dp)                   :: djmz
  real(dp)                   :: djs1
  real(dp)                   :: djs2
  real(dp)                   :: djs3
  real(dp)                   :: djs4
  real(dp)                   :: djs5
  real(dp)                   :: djs6
  real(dp)                   :: dklx
  real(dp)                   :: dkly
  real(dp)                   :: dklz
  real(dp)                   :: dkmx
  real(dp)                   :: dkmy
  real(dp)                   :: dkmz
  real(dp)                   :: dks1
  real(dp)                   :: dks2
  real(dp)                   :: dks3
  real(dp)                   :: dks4
  real(dp)                   :: dks5
  real(dp)                   :: dks6
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charge3')
#endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d2i2s = 0.0_dp
  d2ijs = 0.0_dp
  d2iks = 0.0_dp
  d2j2s = 0.0_dp
  d2jks = 0.0_dp
  d2k2s = 0.0_dp
  do iv = 1,nor
    d2i2s = d2i2s + d2q(1,iv)
    d2ijs = d2ijs + d2q(2,iv)
    d2iks = d2iks + d2q(3,iv)
    d2j2s = d2j2s + d2q(4,iv)
    d2jks = d2jks + d2q(5,iv)
    d2k2s = d2k2s + d2q(6,iv)
  enddo
  d1ql(1:3,1:3,1:3) = d1q(1:3,1:3,1:3)
  if (lstr) d1qsl(1:nstrains,1:3) = d1qs(1:nstrains,1:3)
  if (leem) then
    if (nat(i).gt.maxele) then
      d2i2s = 0.0_dp
      d2ijs = 0.0_dp
      d2iks = 0.0_dp
      d1ql(1:3,1:3,1) = 0.0_dp
      if (lstr) d1qsl(1:nstrains,1) = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2i2s = 0.0_dp
      d2ijs = 0.0_dp
      d2iks = 0.0_dp
      d1ql(1:3,1:3,1) = 0.0_dp
      if (lstr) d1qsl(1:nstrains,1) = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d2jks = 0.0_dp
      d1ql(1:3,1:3,2) = 0.0_dp
      if (lstr) d1qsl(1:nstrains,2) = 0.0_dp
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d2jks = 0.0_dp
      d1ql(1:3,1:3,2) = 0.0_dp
      if (lstr) d1qsl(1:nstrains,2) = 0.0_dp
    endif
    if (nat(k).gt.maxele) then
      d2iks = 0.0_dp
      d2jks = 0.0_dp
      d2k2s = 0.0_dp
      d1ql(1:3,1:3,3) = 0.0_dp
      if (lstr) d1qsl(1:nstrains,3) = 0.0_dp
    elseif (.not.lelementOK(nat(k)).or.nregionno(nrelf2a(k)).ne.1) then
      d2iks = 0.0_dp
      d2jks = 0.0_dp
      d2k2s = 0.0_dp
      d1ql(1:3,1:3,3) = 0.0_dp
      if (lstr) d1qsl(1:nstrains,3) = 0.0_dp
    endif
  endif
  if (lstr) then
    if (ndim.eq.3) then
      dis1 = dqds(1,i)
      dis2 = dqds(2,i)
      dis3 = dqds(3,i)
      dis4 = dqds(4,i)
      dis5 = dqds(5,i)
      dis6 = dqds(6,i)
      djs1 = dqds(1,j)
      djs2 = dqds(2,j)
      djs3 = dqds(3,j)
      djs4 = dqds(4,j)
      djs5 = dqds(5,j)
      djs6 = dqds(6,j)
      dks1 = dqds(1,k)
      dks2 = dqds(2,k)
      dks3 = dqds(3,k)
      dks4 = dqds(4,k)
      dks5 = dqds(5,k)
      dks6 = dqds(6,k)
    elseif (ndim.eq.2) then
      dis1 = dqds(1,i)
      dis2 = dqds(2,i)
      dis3 = dqds(3,i)
      djs1 = dqds(1,j)
      djs2 = dqds(2,j)
      djs3 = dqds(3,j)
      dks1 = dqds(1,k)
      dks2 = dqds(2,k)
      dks3 = dqds(3,k)
    elseif (ndim.eq.1) then
      dis1 = dqds(1,i)
      djs1 = dqds(1,j)
      dks1 = dqds(1,k)
    endif
  endif
!
!  The following flag causes the second derivatives of the charge distribution 
!  with respect to two atoms neither of which are the i/j of the energy term
!  to be calculated. At the moment this contribution is not used in any method
!  and so the flag is set to false for efficiency.
!
  lNonIJQDeriv = .false.
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lreaxFF) then
    deisum = 0.0_dp
    dejsum = 0.0_dp
    deksum = 0.0_dp
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
      deksum = deksum + dek(iv)
    enddo
!     
!  Strain - strain contribution  
!     
    if (lstr.and.ndim.gt.0) then
      ind = 0                
      do kk = 1,nstrains     
        do kl = 1,kk-1       
          ind = ind + 1      
          sderv2(kl,kk) = sderv2(kl,kk) + deisum*d2qds2(ind,i) + dejsum*d2qds2(ind,j) + deksum*d2qds2(ind,k)
          sderv2(kk,kl) = sderv2(kk,kl) + deisum*d2qds2(ind,i) + dejsum*d2qds2(ind,j) + deksum*d2qds2(ind,k)
        enddo                
        ind = ind + 1        
        sderv2(kk,kk) = sderv2(kk,kk) + deisum*d2qds2(ind,i) + dejsum*d2qds2(ind,j) + deksum*d2qds2(ind,k)
      enddo                  
    endif
!
!  Terms involving Q derivatives for i/j/k and one other atom
!
    do m = 1,nqatoms(i)
      p = nqatomptr(m,i)
      px = 3*(p - 1) + 1
      py = px + 1
      pz = py + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (i.gt.p) then
        derv2(px,ix) = derv2(px,ix) + deisum*d2qdxyz2(mnxx,i)
        derv2(py,ix) = derv2(py,ix) + deisum*d2qdxyz2(mnxy,i)
        derv2(pz,ix) = derv2(pz,ix) + deisum*d2qdxyz2(mnxz,i)
        derv2(px,iy) = derv2(px,iy) + deisum*d2qdxyz2(mnxy,i)
        derv2(py,iy) = derv2(py,iy) + deisum*d2qdxyz2(mnyy,i)
        derv2(pz,iy) = derv2(pz,iy) + deisum*d2qdxyz2(mnyz,i)
        derv2(px,iz) = derv2(px,iz) + deisum*d2qdxyz2(mnxz,i)
        derv2(py,iz) = derv2(py,iz) + deisum*d2qdxyz2(mnyz,i)
        derv2(pz,iz) = derv2(pz,iz) + deisum*d2qdxyz2(mnzz,i)
      elseif (i.lt.p) then
        derv2(ix,px) = derv2(ix,px) + deisum*d2qdxyz2(mnxx,i)
        derv2(iy,px) = derv2(iy,px) + deisum*d2qdxyz2(mnxy,i)
        derv2(iz,px) = derv2(iz,px) + deisum*d2qdxyz2(mnxz,i)
        derv2(ix,py) = derv2(ix,py) + deisum*d2qdxyz2(mnxy,i)
        derv2(iy,py) = derv2(iy,py) + deisum*d2qdxyz2(mnyy,i)
        derv2(iz,py) = derv2(iz,py) + deisum*d2qdxyz2(mnyz,i)
        derv2(ix,pz) = derv2(ix,pz) + deisum*d2qdxyz2(mnxz,i)
        derv2(iy,pz) = derv2(iy,pz) + deisum*d2qdxyz2(mnyz,i)
        derv2(iz,pz) = derv2(iz,pz) + deisum*d2qdxyz2(mnzz,i)
      endif
!
!  Mixed strain contribution
!
      if (lstr.and.ndim.gt.0) then
        do kk = 1,nstrains
          derv3(px,kk) = derv3(px,kk) + deisum*d2qdxyzs(kk,mx,i)
          derv3(py,kk) = derv3(py,kk) + deisum*d2qdxyzs(kk,my,i)
          derv3(pz,kk) = derv3(pz,kk) + deisum*d2qdxyzs(kk,mz,i)
          derv3(ix,kk) = derv3(ix,kk) - deisum*d2qdxyzs(kk,mx,i)
          derv3(iy,kk) = derv3(iy,kk) - deisum*d2qdxyzs(kk,my,i)
          derv3(iz,kk) = derv3(iz,kk) - deisum*d2qdxyzs(kk,mz,i)
        enddo
      endif
    enddo
!
    do m = 1,nqatoms(j)
      p = nqatomptr(m,j)
      px = 3*(p - 1) + 1
      py = px + 1
      pz = py + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (j.gt.p) then
        derv2(px,jx) = derv2(px,jx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(py,jx) = derv2(py,jx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(pz,jx) = derv2(pz,jx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(px,jy) = derv2(px,jy) + dejsum*d2qdxyz2(mnxy,j)
        derv2(py,jy) = derv2(py,jy) + dejsum*d2qdxyz2(mnyy,j)
        derv2(pz,jy) = derv2(pz,jy) + dejsum*d2qdxyz2(mnyz,j)
        derv2(px,jz) = derv2(px,jz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(py,jz) = derv2(py,jz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(pz,jz) = derv2(pz,jz) + dejsum*d2qdxyz2(mnzz,j)
      elseif (j.lt.p) then
        derv2(jx,px) = derv2(jx,px) + dejsum*d2qdxyz2(mnxx,j)
        derv2(jy,px) = derv2(jy,px) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jz,px) = derv2(jz,px) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jx,py) = derv2(jx,py) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jy,py) = derv2(jy,py) + dejsum*d2qdxyz2(mnyy,j)
        derv2(jz,py) = derv2(jz,py) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jx,pz) = derv2(jx,pz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jy,pz) = derv2(jy,pz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jz,pz) = derv2(jz,pz) + dejsum*d2qdxyz2(mnzz,j)
      endif
!
!  Mixed - strain contribution
!
      if (lstr.and.ndim.gt.0) then
        do kk = 1,nstrains
          derv3(px,kk) = derv3(px,kk) + dejsum*d2qdxyzs(kk,mx,j)
          derv3(py,kk) = derv3(py,kk) + dejsum*d2qdxyzs(kk,my,j)
          derv3(pz,kk) = derv3(pz,kk) + dejsum*d2qdxyzs(kk,mz,j)
          derv3(jx,kk) = derv3(jx,kk) - dejsum*d2qdxyzs(kk,mx,j)
          derv3(jy,kk) = derv3(jy,kk) - dejsum*d2qdxyzs(kk,my,j)
          derv3(jz,kk) = derv3(jz,kk) - dejsum*d2qdxyzs(kk,mz,j)
        enddo
      endif
    enddo
!
    do m = 1,nqatoms(k)
      p = nqatomptr(m,k)
      px = 3*(p - 1) + 1
      py = px + 1
      pz = py + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (k.gt.p) then
        derv2(px,kx) = derv2(px,kx) + deksum*d2qdxyz2(mnxx,k)
        derv2(py,kx) = derv2(py,kx) + deksum*d2qdxyz2(mnxy,k)
        derv2(pz,kx) = derv2(pz,kx) + deksum*d2qdxyz2(mnxz,k)
        derv2(px,ky) = derv2(px,ky) + deksum*d2qdxyz2(mnxy,k)
        derv2(py,ky) = derv2(py,ky) + deksum*d2qdxyz2(mnyy,k)
        derv2(pz,ky) = derv2(pz,ky) + deksum*d2qdxyz2(mnyz,k)
        derv2(px,kz) = derv2(px,kz) + deksum*d2qdxyz2(mnxz,k)
        derv2(py,kz) = derv2(py,kz) + deksum*d2qdxyz2(mnyz,k)
        derv2(pz,kz) = derv2(pz,kz) + deksum*d2qdxyz2(mnzz,k)
      elseif (k.lt.p) then
        derv2(kx,px) = derv2(kx,px) + deksum*d2qdxyz2(mnxx,k)
        derv2(ky,px) = derv2(ky,px) + deksum*d2qdxyz2(mnxy,k)
        derv2(kz,px) = derv2(kz,px) + deksum*d2qdxyz2(mnxz,k)
        derv2(kx,py) = derv2(kx,py) + deksum*d2qdxyz2(mnxy,k)
        derv2(ky,py) = derv2(ky,py) + deksum*d2qdxyz2(mnyy,k)
        derv2(kz,py) = derv2(kz,py) + deksum*d2qdxyz2(mnyz,k)
        derv2(kx,pz) = derv2(kx,pz) + deksum*d2qdxyz2(mnxz,k)
        derv2(ky,pz) = derv2(ky,pz) + deksum*d2qdxyz2(mnyz,k)
        derv2(kz,pz) = derv2(kz,pz) + deksum*d2qdxyz2(mnzz,k)
      endif
!
!  Mixed - strain contribution
!
      if (lstr.and.ndim.gt.0) then
        do kk = 1,nstrains
          derv3(px,kk) = derv3(px,kk) + deksum*d2qdxyzs(kk,mx,k)
          derv3(py,kk) = derv3(py,kk) + deksum*d2qdxyzs(kk,my,k)
          derv3(pz,kk) = derv3(pz,kk) + deksum*d2qdxyzs(kk,mz,k)
          derv3(kx,kk) = derv3(kx,kk) - deksum*d2qdxyzs(kk,mx,k)
          derv3(ky,kk) = derv3(ky,kk) - deksum*d2qdxyzs(kk,my,k)
          derv3(kz,kk) = derv3(kz,kk) - deksum*d2qdxyzs(kk,mz,k)
        enddo
      endif
    enddo
!
!  Terms involving Q derivatives for non-i/j case
!
    if (lNonIJQDeriv) then
      do m = 2,nqatoms(i)
        p = nqatomptr(m,i)
        px = 3*(p - 1) + 1
        py = px + 1
        pz = py + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,i)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,px) = derv2(lx,px) + deisum*d2qdxyz2(mnxx,i)
          derv2(ly,px) = derv2(ly,px) + deisum*d2qdxyz2(mnxy,i)
          derv2(lz,px) = derv2(lz,px) + deisum*d2qdxyz2(mnxz,i)
          derv2(lx,py) = derv2(lx,py) + deisum*d2qdxyz2(mnyx,i)
          derv2(ly,py) = derv2(ly,py) + deisum*d2qdxyz2(mnyy,i)
          derv2(lz,py) = derv2(lz,py) + deisum*d2qdxyz2(mnyz,i)
          derv2(lx,pz) = derv2(lx,pz) + deisum*d2qdxyz2(mnzx,i)
          derv2(ly,pz) = derv2(ly,pz) + deisum*d2qdxyz2(mnzy,i)
          derv2(lz,pz) = derv2(lz,pz) + deisum*d2qdxyz2(mnzz,i)
        enddo
      enddo
!
      do m = 2,nqatoms(j)
        p = nqatomptr(m,j)
        px = 3*(p - 1) + 1
        py = px + 1
        pz = py + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,j)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,px) = derv2(lx,px) + deisum*d2qdxyz2(mnxx,j)
          derv2(ly,px) = derv2(ly,px) + deisum*d2qdxyz2(mnxy,j)
          derv2(lz,px) = derv2(lz,px) + deisum*d2qdxyz2(mnxz,j)
          derv2(lx,py) = derv2(lx,py) + deisum*d2qdxyz2(mnyx,j)
          derv2(ly,py) = derv2(ly,py) + deisum*d2qdxyz2(mnyy,j)
          derv2(lz,py) = derv2(lz,py) + deisum*d2qdxyz2(mnyz,j)
          derv2(lx,pz) = derv2(lx,pz) + deisum*d2qdxyz2(mnzx,j)
          derv2(ly,pz) = derv2(ly,pz) + deisum*d2qdxyz2(mnzy,j)
          derv2(lz,pz) = derv2(lz,pz) + deisum*d2qdxyz2(mnzz,j)
        enddo
      enddo
!
      do m = 2,nqatoms(k)
        p = nqatomptr(m,k)
        px = 3*(p - 1) + 1
        py = px + 1
        pz = py + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,k)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,px) = derv2(lx,px) + deksum*d2qdxyz2(mnxx,k)
          derv2(ly,px) = derv2(ly,px) + deksum*d2qdxyz2(mnxy,k)
          derv2(lz,px) = derv2(lz,px) + deksum*d2qdxyz2(mnxz,k)
          derv2(lx,py) = derv2(lx,py) + deksum*d2qdxyz2(mnyx,k)
          derv2(ly,py) = derv2(ly,py) + deksum*d2qdxyz2(mnyy,k)
          derv2(lz,py) = derv2(lz,py) + deksum*d2qdxyz2(mnyz,k)
          derv2(lx,pz) = derv2(lx,pz) + deksum*d2qdxyz2(mnzx,k)
          derv2(ly,pz) = derv2(ly,pz) + deksum*d2qdxyz2(mnzy,k)
          derv2(lz,pz) = derv2(lz,pz) + deksum*d2qdxyz2(mnzz,k)
        enddo
      enddo
    endif
  endif
!
!  Loop over atoms to add ij-m contributions
!
!  Note indices must be set to cope with the fact that only
!  half of derv2 is constructed in upper triangle and the
!  effects of lopi/lopj/lopk if freezing is being used.
!
  mx = - 2
  my = - 1
  mz =   0
  do m = 1,numat
    lopm = (.not.lfreeze.or.lopf(m))
    indm = 3*(m-1)
    if (lopm) then
      mx = mx + 3
      my = my + 3
      mz = mz + 3
    endif
    dimx = dqdxyz(indm+1,i)
    dimy = dqdxyz(indm+2,i)
    dimz = dqdxyz(indm+3,i)
    djmx = dqdxyz(indm+1,j)
    djmy = dqdxyz(indm+2,j)
    djmz = dqdxyz(indm+3,j)
    dkmx = dqdxyz(indm+1,k)
    dkmy = dqdxyz(indm+2,k)
    dkmz = dqdxyz(indm+3,k)
    if (i.ne.m.and.(lopi.or.lopm)) then
      if (lopi.and.lopm) then
        if (m.le.i) then
          ixl = ix
          iyl = iy
          izl = iz
          mxl = mx
          myl = my
          mzl = mz
          lReverse = .true.
        else
          ixl = mx
          iyl = my
          izl = mz
          mxl = ix
          myl = iy
          mzl = iz
          lReverse = .true.
        endif
      elseif (lopi) then
        ixl = ix
        iyl = iy
        izl = iz
        mxl = ix
        myl = iy
        mzl = iz
        lReverse = .false.
      elseif (lopm) then
        ixl = mx
        iyl = my
        izl = mz
        mxl = mx
        myl = my
        mzl = mz
        lReverse = .true.
      endif
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-m
!
      if (lReverse) then
        derv2(mxl,ixl) = derv2(mxl,ixl) - d1ql(1,1,1)*dimx - d1ql(1,2,1)*dimx
        derv2(myl,ixl) = derv2(myl,ixl) - d1ql(2,1,1)*dimx - d1ql(2,2,1)*dimx
        derv2(mzl,ixl) = derv2(mzl,ixl) - d1ql(3,1,1)*dimx - d1ql(3,2,1)*dimx
        derv2(mxl,iyl) = derv2(mxl,iyl) - d1ql(1,1,1)*dimy - d1ql(1,2,1)*dimy
        derv2(myl,iyl) = derv2(myl,iyl) - d1ql(2,1,1)*dimy - d1ql(2,2,1)*dimy
        derv2(mzl,iyl) = derv2(mzl,iyl) - d1ql(3,1,1)*dimy - d1ql(3,2,1)*dimy
        derv2(mxl,izl) = derv2(mxl,izl) - d1ql(1,1,1)*dimz - d1ql(1,2,1)*dimz
        derv2(myl,izl) = derv2(myl,izl) - d1ql(2,1,1)*dimz - d1ql(2,2,1)*dimz
        derv2(mzl,izl) = derv2(mzl,izl) - d1ql(3,1,1)*dimz - d1ql(3,2,1)*dimz
!
        derv2(mxl,ixl) = derv2(mxl,ixl) - d1ql(1,1,2)*djmx - d1ql(1,2,2)*djmx
        derv2(myl,ixl) = derv2(myl,ixl) - d1ql(2,1,2)*djmx - d1ql(2,2,2)*djmx
        derv2(mzl,ixl) = derv2(mzl,ixl) - d1ql(3,1,2)*djmx - d1ql(3,2,2)*djmx
        derv2(mxl,iyl) = derv2(mxl,iyl) - d1ql(1,1,2)*djmy - d1ql(1,2,2)*djmy
        derv2(myl,iyl) = derv2(myl,iyl) - d1ql(2,1,2)*djmy - d1ql(2,2,2)*djmy
        derv2(mzl,iyl) = derv2(mzl,iyl) - d1ql(3,1,2)*djmy - d1ql(3,2,2)*djmy
        derv2(mxl,izl) = derv2(mxl,izl) - d1ql(1,1,2)*djmz - d1ql(1,2,2)*djmz
        derv2(myl,izl) = derv2(myl,izl) - d1ql(2,1,2)*djmz - d1ql(2,2,2)*djmz
        derv2(mzl,izl) = derv2(mzl,izl) - d1ql(3,1,2)*djmz - d1ql(3,2,2)*djmz
!
        derv2(mxl,ixl) = derv2(mxl,ixl) - d1ql(1,1,3)*dkmx - d1ql(1,2,3)*dkmx
        derv2(myl,ixl) = derv2(myl,ixl) - d1ql(2,1,3)*dkmx - d1ql(2,2,3)*dkmx
        derv2(mzl,ixl) = derv2(mzl,ixl) - d1ql(3,1,3)*dkmx - d1ql(3,2,3)*dkmx
        derv2(mxl,iyl) = derv2(mxl,iyl) - d1ql(1,1,3)*dkmy - d1ql(1,2,3)*dkmy
        derv2(myl,iyl) = derv2(myl,iyl) - d1ql(2,1,3)*dkmy - d1ql(2,2,3)*dkmy
        derv2(mzl,iyl) = derv2(mzl,iyl) - d1ql(3,1,3)*dkmy - d1ql(3,2,3)*dkmy
        derv2(mxl,izl) = derv2(mxl,izl) - d1ql(1,1,3)*dkmz - d1ql(1,2,3)*dkmz
        derv2(myl,izl) = derv2(myl,izl) - d1ql(2,1,3)*dkmz - d1ql(2,2,3)*dkmz
        derv2(mzl,izl) = derv2(mzl,izl) - d1ql(3,1,3)*dkmz - d1ql(3,2,3)*dkmz
      else
        derv2(mxl,ixl) = derv2(mxl,ixl) - d1ql(1,1,1)*dimx - d1ql(1,2,1)*dimx
        derv2(myl,ixl) = derv2(myl,ixl) - d1ql(1,1,1)*dimy - d1ql(1,2,1)*dimy
        derv2(mzl,ixl) = derv2(mzl,ixl) - d1ql(1,1,1)*dimz - d1ql(1,2,1)*dimz
        derv2(mxl,iyl) = derv2(mxl,iyl) - d1ql(2,1,1)*dimx - d1ql(2,2,1)*dimx
        derv2(myl,iyl) = derv2(myl,iyl) - d1ql(2,1,1)*dimy - d1ql(2,2,1)*dimy
        derv2(mzl,iyl) = derv2(mzl,iyl) - d1ql(2,1,1)*dimz - d1ql(2,2,1)*dimz
        derv2(mxl,izl) = derv2(mxl,izl) - d1ql(3,1,1)*dimx - d1ql(3,2,1)*dimx
        derv2(myl,izl) = derv2(myl,izl) - d1ql(3,1,1)*dimy - d1ql(3,2,1)*dimy
        derv2(mzl,izl) = derv2(mzl,izl) - d1ql(3,1,1)*dimz - d1ql(3,2,1)*dimz
!
        derv2(mxl,ixl) = derv2(mxl,ixl) - d1ql(1,1,2)*djmx - d1ql(1,2,2)*djmx
        derv2(myl,ixl) = derv2(myl,ixl) - d1ql(1,1,2)*djmy - d1ql(1,2,2)*djmy
        derv2(mzl,ixl) = derv2(mzl,ixl) - d1ql(1,1,2)*djmz - d1ql(1,2,2)*djmz
        derv2(mxl,iyl) = derv2(mxl,iyl) - d1ql(2,1,2)*djmx - d1ql(2,2,2)*djmx
        derv2(myl,iyl) = derv2(myl,iyl) - d1ql(2,1,2)*djmy - d1ql(2,2,2)*djmy
        derv2(mzl,iyl) = derv2(mzl,iyl) - d1ql(2,1,2)*djmz - d1ql(2,2,2)*djmz
        derv2(mxl,izl) = derv2(mxl,izl) - d1ql(3,1,2)*djmx - d1ql(3,2,2)*djmx
        derv2(myl,izl) = derv2(myl,izl) - d1ql(3,1,2)*djmy - d1ql(3,2,2)*djmy
        derv2(mzl,izl) = derv2(mzl,izl) - d1ql(3,1,2)*djmz - d1ql(3,2,2)*djmz
!
        derv2(mxl,ixl) = derv2(mxl,ixl) - d1ql(1,1,3)*dkmx - d1ql(1,2,3)*dkmx
        derv2(myl,ixl) = derv2(myl,ixl) - d1ql(1,1,3)*dkmy - d1ql(1,2,3)*dkmy
        derv2(mzl,ixl) = derv2(mzl,ixl) - d1ql(1,1,3)*dkmz - d1ql(1,2,3)*dkmz
        derv2(mxl,iyl) = derv2(mxl,iyl) - d1ql(2,1,3)*dkmx - d1ql(2,2,3)*dkmx
        derv2(myl,iyl) = derv2(myl,iyl) - d1ql(2,1,3)*dkmy - d1ql(2,2,3)*dkmy
        derv2(mzl,iyl) = derv2(mzl,iyl) - d1ql(2,1,3)*dkmz - d1ql(2,2,3)*dkmz
        derv2(mxl,izl) = derv2(mxl,izl) - d1ql(3,1,3)*dkmx - d1ql(3,2,3)*dkmx
        derv2(myl,izl) = derv2(myl,izl) - d1ql(3,1,3)*dkmy - d1ql(3,2,3)*dkmy
        derv2(mzl,izl) = derv2(mzl,izl) - d1ql(3,1,3)*dkmz - d1ql(3,2,3)*dkmz
      endif
    endif
    if (j.ne.m.and.(lopj.or.lopm)) then
      if (lopj.and.lopm) then
        if (m.le.j) then
          jxl = jx
          jyl = jy
          jzl = jz
          mxl = mx
          myl = my
          mzl = mz
          lReverse = .true.
        else
          jxl = mx
          jyl = my
          jzl = mz
          mxl = jx
          myl = jy
          mzl = jz
          lReverse = .true.
        endif
      elseif (lopj) then
        jxl = jx
        jyl = jy
        jzl = jz
        mxl = jx
        myl = jy
        mzl = jz
        lReverse = .false.
      elseif (lopm) then
        jxl = mx
        jyl = my
        jzl = mz
        mxl = mx
        myl = my
        mzl = mz
        lReverse = .true.
      endif
!
!  d2E/(d(alpha).dq) x dq/d(beta) : j-m
!
      if (lReverse) then
        derv2(mxl,jxl) = derv2(mxl,jxl) + d1ql(1,1,1)*dimx - d1ql(1,3,1)*dimx
        derv2(myl,jxl) = derv2(myl,jxl) + d1ql(2,1,1)*dimx - d1ql(2,3,1)*dimx
        derv2(mzl,jxl) = derv2(mzl,jxl) + d1ql(3,1,1)*dimx - d1ql(3,3,1)*dimx
        derv2(mxl,jyl) = derv2(mxl,jyl) + d1ql(1,1,1)*dimy - d1ql(1,3,1)*dimy
        derv2(myl,jyl) = derv2(myl,jyl) + d1ql(2,1,1)*dimy - d1ql(2,3,1)*dimy
        derv2(mzl,jyl) = derv2(mzl,jyl) + d1ql(3,1,1)*dimy - d1ql(3,3,1)*dimy
        derv2(mxl,jzl) = derv2(mxl,jzl) + d1ql(1,1,1)*dimz - d1ql(1,3,1)*dimz
        derv2(myl,jzl) = derv2(myl,jzl) + d1ql(2,1,1)*dimz - d1ql(2,3,1)*dimz
        derv2(mzl,jzl) = derv2(mzl,jzl) + d1ql(3,1,1)*dimz - d1ql(3,3,1)*dimz
!
        derv2(mxl,jxl) = derv2(mxl,jxl) + d1ql(1,1,2)*djmx - d1ql(1,3,2)*djmx
        derv2(myl,jxl) = derv2(myl,jxl) + d1ql(2,1,2)*djmx - d1ql(2,3,2)*djmx
        derv2(mzl,jxl) = derv2(mzl,jxl) + d1ql(3,1,2)*djmx - d1ql(3,3,2)*djmx
        derv2(mxl,jyl) = derv2(mxl,jyl) + d1ql(1,1,2)*djmy - d1ql(1,3,2)*djmy
        derv2(myl,jyl) = derv2(myl,jyl) + d1ql(2,1,2)*djmy - d1ql(2,3,2)*djmy
        derv2(mzl,jyl) = derv2(mzl,jyl) + d1ql(3,1,2)*djmy - d1ql(3,3,2)*djmy
        derv2(mxl,jzl) = derv2(mxl,jzl) + d1ql(1,1,2)*djmz - d1ql(1,3,2)*djmz
        derv2(myl,jzl) = derv2(myl,jzl) + d1ql(2,1,2)*djmz - d1ql(2,3,2)*djmz
        derv2(mzl,jzl) = derv2(mzl,jzl) + d1ql(3,1,2)*djmz - d1ql(3,3,2)*djmz
!
        derv2(mxl,jxl) = derv2(mxl,jxl) + d1ql(1,1,3)*dkmx - d1ql(1,3,3)*dkmx
        derv2(myl,jxl) = derv2(myl,jxl) + d1ql(2,1,3)*dkmx - d1ql(2,3,3)*dkmx
        derv2(mzl,jxl) = derv2(mzl,jxl) + d1ql(3,1,3)*dkmx - d1ql(3,3,3)*dkmx
        derv2(mxl,jyl) = derv2(mxl,jyl) + d1ql(1,1,3)*dkmy - d1ql(1,3,3)*dkmy
        derv2(myl,jyl) = derv2(myl,jyl) + d1ql(2,1,3)*dkmy - d1ql(2,3,3)*dkmy
        derv2(mzl,jyl) = derv2(mzl,jyl) + d1ql(3,1,3)*dkmy - d1ql(3,3,3)*dkmy
        derv2(mxl,jzl) = derv2(mxl,jzl) + d1ql(1,1,3)*dkmz - d1ql(1,3,3)*dkmz
        derv2(myl,jzl) = derv2(myl,jzl) + d1ql(2,1,3)*dkmz - d1ql(2,3,3)*dkmz
        derv2(mzl,jzl) = derv2(mzl,jzl) + d1ql(3,1,3)*dkmz - d1ql(3,3,3)*dkmz
      else
        derv2(mxl,jxl) = derv2(mxl,jxl) + d1ql(1,1,1)*dimx - d1ql(1,3,1)*dimx
        derv2(myl,jxl) = derv2(myl,jxl) + d1ql(1,1,1)*dimy - d1ql(1,3,1)*dimy
        derv2(mzl,jxl) = derv2(mzl,jxl) + d1ql(1,1,1)*dimz - d1ql(1,3,1)*dimz
        derv2(mxl,jyl) = derv2(mxl,jyl) + d1ql(2,1,1)*dimx - d1ql(2,3,1)*dimx
        derv2(myl,jyl) = derv2(myl,jyl) + d1ql(2,1,1)*dimy - d1ql(2,3,1)*dimy
        derv2(mzl,jyl) = derv2(mzl,jyl) + d1ql(2,1,1)*dimz - d1ql(2,3,1)*dimz
        derv2(mxl,jzl) = derv2(mxl,jzl) + d1ql(3,1,1)*dimx - d1ql(3,3,1)*dimx
        derv2(myl,jzl) = derv2(myl,jzl) + d1ql(3,1,1)*dimy - d1ql(3,3,1)*dimy
        derv2(mzl,jzl) = derv2(mzl,jzl) + d1ql(3,1,1)*dimz - d1ql(3,3,1)*dimz
!
        derv2(mxl,jxl) = derv2(mxl,jxl) + d1ql(1,1,2)*djmx - d1ql(1,3,2)*djmx
        derv2(myl,jxl) = derv2(myl,jxl) + d1ql(1,1,2)*djmy - d1ql(1,3,2)*djmy
        derv2(mzl,jxl) = derv2(mzl,jxl) + d1ql(1,1,2)*djmz - d1ql(1,3,2)*djmz
        derv2(mxl,jyl) = derv2(mxl,jyl) + d1ql(2,1,2)*djmx - d1ql(2,3,2)*djmx
        derv2(myl,jyl) = derv2(myl,jyl) + d1ql(2,1,2)*djmy - d1ql(2,3,2)*djmy
        derv2(mzl,jyl) = derv2(mzl,jyl) + d1ql(2,1,2)*djmz - d1ql(2,3,2)*djmz
        derv2(mxl,jzl) = derv2(mxl,jzl) + d1ql(3,1,2)*djmx - d1ql(3,3,2)*djmx
        derv2(myl,jzl) = derv2(myl,jzl) + d1ql(3,1,2)*djmy - d1ql(3,3,2)*djmy
        derv2(mzl,jzl) = derv2(mzl,jzl) + d1ql(3,1,2)*djmz - d1ql(3,3,2)*djmz
!
        derv2(mxl,jxl) = derv2(mxl,jxl) + d1ql(1,1,3)*dkmx - d1ql(1,3,3)*dkmx
        derv2(myl,jxl) = derv2(myl,jxl) + d1ql(1,1,3)*dkmy - d1ql(1,3,3)*dkmy
        derv2(mzl,jxl) = derv2(mzl,jxl) + d1ql(1,1,3)*dkmz - d1ql(1,3,3)*dkmz
        derv2(mxl,jyl) = derv2(mxl,jyl) + d1ql(2,1,3)*dkmx - d1ql(2,3,3)*dkmx
        derv2(myl,jyl) = derv2(myl,jyl) + d1ql(2,1,3)*dkmy - d1ql(2,3,3)*dkmy
        derv2(mzl,jyl) = derv2(mzl,jyl) + d1ql(2,1,3)*dkmz - d1ql(2,3,3)*dkmz
        derv2(mxl,jzl) = derv2(mxl,jzl) + d1ql(3,1,3)*dkmx - d1ql(3,3,3)*dkmx
        derv2(myl,jzl) = derv2(myl,jzl) + d1ql(3,1,3)*dkmy - d1ql(3,3,3)*dkmy
        derv2(mzl,jzl) = derv2(mzl,jzl) + d1ql(3,1,3)*dkmz - d1ql(3,3,3)*dkmz
      endif
    endif
    if (k.ne.m.and.(lopk.or.lopm)) then
      if (lopk.and.lopm) then
        if (m.le.k) then
          kxl = kx
          kyl = ky
          kzl = kz
          mxl = mx
          myl = my
          mzl = mz
          lReverse = .false.
        else
          kxl = mx
          kyl = my
          kzl = mz
          mxl = kx
          myl = ky
          mzl = kz
          lReverse = .true.
        endif
      elseif (lopk) then
        kxl = kx
        kyl = ky
        kzl = kz
        mxl = kx
        myl = ky
        mzl = kz
        lReverse = .false.
      elseif (lopm) then
        kxl = mx
        kyl = my
        kzl = mz
        mxl = mx
        myl = my
        mzl = mz
        lReverse = .true.
      endif
!
!  d2E/(d(alpha).dq) x dq/d(beta) : k-m
!
      if (lReverse) then
        derv2(mxl,kxl) = derv2(mxl,kxl) + d1ql(1,2,1)*dimx + d1ql(1,3,1)*dimx
        derv2(myl,kxl) = derv2(myl,kxl) + d1ql(2,2,1)*dimx + d1ql(2,3,1)*dimx
        derv2(mzl,kxl) = derv2(mzl,kxl) + d1ql(3,2,1)*dimx + d1ql(3,3,1)*dimx
        derv2(mxl,kyl) = derv2(mxl,kyl) + d1ql(1,2,1)*dimy + d1ql(1,3,1)*dimy
        derv2(myl,kyl) = derv2(myl,kyl) + d1ql(2,2,1)*dimy + d1ql(2,3,1)*dimy
        derv2(mzl,kyl) = derv2(mzl,kyl) + d1ql(3,2,1)*dimy + d1ql(3,3,1)*dimy
        derv2(mxl,kzl) = derv2(mxl,kzl) + d1ql(1,2,1)*dimz + d1ql(1,3,1)*dimz
        derv2(myl,kzl) = derv2(myl,kzl) + d1ql(2,2,1)*dimz + d1ql(2,3,1)*dimz
        derv2(mzl,kzl) = derv2(mzl,kzl) + d1ql(3,2,1)*dimz + d1ql(3,3,1)*dimz
!
        derv2(mxl,kxl) = derv2(mxl,kxl) + d1ql(1,2,2)*djmx + d1ql(1,3,2)*djmx
        derv2(myl,kxl) = derv2(myl,kxl) + d1ql(2,2,2)*djmx + d1ql(2,3,2)*djmx
        derv2(mzl,kxl) = derv2(mzl,kxl) + d1ql(3,2,2)*djmx + d1ql(3,3,2)*djmx
        derv2(mxl,kyl) = derv2(mxl,kyl) + d1ql(1,2,2)*djmy + d1ql(1,3,2)*djmy
        derv2(myl,kyl) = derv2(myl,kyl) + d1ql(2,2,2)*djmy + d1ql(2,3,2)*djmy
        derv2(mzl,kyl) = derv2(mzl,kyl) + d1ql(3,2,2)*djmy + d1ql(3,3,2)*djmy
        derv2(mxl,kzl) = derv2(mxl,kzl) + d1ql(1,2,2)*djmz + d1ql(1,3,2)*djmz
        derv2(myl,kzl) = derv2(myl,kzl) + d1ql(2,2,2)*djmz + d1ql(2,3,2)*djmz
        derv2(mzl,kzl) = derv2(mzl,kzl) + d1ql(3,2,2)*djmz + d1ql(3,3,2)*djmz
!
        derv2(mxl,kxl) = derv2(mxl,kxl) + d1ql(1,2,3)*dkmx + d1ql(1,3,3)*dkmx
        derv2(myl,kxl) = derv2(myl,kxl) + d1ql(2,2,3)*dkmx + d1ql(2,3,3)*dkmx
        derv2(mzl,kxl) = derv2(mzl,kxl) + d1ql(3,2,3)*dkmx + d1ql(3,3,3)*dkmx
        derv2(mxl,kyl) = derv2(mxl,kyl) + d1ql(1,2,3)*dkmy + d1ql(1,3,3)*dkmy
        derv2(myl,kyl) = derv2(myl,kyl) + d1ql(2,2,3)*dkmy + d1ql(2,3,3)*dkmy
        derv2(mzl,kyl) = derv2(mzl,kyl) + d1ql(3,2,3)*dkmy + d1ql(3,3,3)*dkmy
        derv2(mxl,kzl) = derv2(mxl,kzl) + d1ql(1,2,3)*dkmz + d1ql(1,3,3)*dkmz
        derv2(myl,kzl) = derv2(myl,kzl) + d1ql(2,2,3)*dkmz + d1ql(2,3,3)*dkmz
        derv2(mzl,kzl) = derv2(mzl,kzl) + d1ql(3,2,3)*dkmz + d1ql(3,3,3)*dkmz
      else
        derv2(mxl,kxl) = derv2(mxl,kxl) + d1ql(1,2,1)*dimx + d1ql(1,3,1)*dimx
        derv2(myl,kxl) = derv2(myl,kxl) + d1ql(1,2,1)*dimy + d1ql(1,3,1)*dimy
        derv2(mzl,kxl) = derv2(mzl,kxl) + d1ql(1,2,1)*dimz + d1ql(1,3,1)*dimz
        derv2(mxl,kyl) = derv2(mxl,kyl) + d1ql(2,2,1)*dimx + d1ql(2,3,1)*dimx
        derv2(myl,kyl) = derv2(myl,kyl) + d1ql(2,2,1)*dimy + d1ql(2,3,1)*dimy
        derv2(mzl,kyl) = derv2(mzl,kyl) + d1ql(2,2,1)*dimz + d1ql(2,3,1)*dimz
        derv2(mxl,kzl) = derv2(mxl,kzl) + d1ql(3,2,1)*dimx + d1ql(3,3,1)*dimx
        derv2(myl,kzl) = derv2(myl,kzl) + d1ql(3,2,1)*dimy + d1ql(3,3,1)*dimy
        derv2(mzl,kzl) = derv2(mzl,kzl) + d1ql(3,2,1)*dimz + d1ql(3,3,1)*dimz
!
        derv2(mxl,kxl) = derv2(mxl,kxl) + d1ql(1,2,2)*djmx + d1ql(1,3,2)*djmx
        derv2(myl,kxl) = derv2(myl,kxl) + d1ql(1,2,2)*djmy + d1ql(1,3,2)*djmy
        derv2(mzl,kxl) = derv2(mzl,kxl) + d1ql(1,2,2)*djmz + d1ql(1,3,2)*djmz
        derv2(mxl,kyl) = derv2(mxl,kyl) + d1ql(2,2,2)*djmx + d1ql(2,3,2)*djmx
        derv2(myl,kyl) = derv2(myl,kyl) + d1ql(2,2,2)*djmy + d1ql(2,3,2)*djmy
        derv2(mzl,kyl) = derv2(mzl,kyl) + d1ql(2,2,2)*djmz + d1ql(2,3,2)*djmz
        derv2(mxl,kzl) = derv2(mxl,kzl) + d1ql(3,2,2)*djmx + d1ql(3,3,2)*djmx
        derv2(myl,kzl) = derv2(myl,kzl) + d1ql(3,2,2)*djmy + d1ql(3,3,2)*djmy
        derv2(mzl,kzl) = derv2(mzl,kzl) + d1ql(3,2,2)*djmz + d1ql(3,3,2)*djmz
!
        derv2(mxl,kxl) = derv2(mxl,kxl) + d1ql(1,2,3)*dkmx + d1ql(1,3,3)*dkmx
        derv2(myl,kxl) = derv2(myl,kxl) + d1ql(1,2,3)*dkmy + d1ql(1,3,3)*dkmy
        derv2(mzl,kxl) = derv2(mzl,kxl) + d1ql(1,2,3)*dkmz + d1ql(1,3,3)*dkmz
        derv2(mxl,kyl) = derv2(mxl,kyl) + d1ql(2,2,3)*dkmx + d1ql(2,3,3)*dkmx
        derv2(myl,kyl) = derv2(myl,kyl) + d1ql(2,2,3)*dkmy + d1ql(2,3,3)*dkmy
        derv2(mzl,kyl) = derv2(mzl,kyl) + d1ql(2,2,3)*dkmz + d1ql(2,3,3)*dkmz
        derv2(mxl,kzl) = derv2(mxl,kzl) + d1ql(3,2,3)*dkmx + d1ql(3,3,3)*dkmx
        derv2(myl,kzl) = derv2(myl,kzl) + d1ql(3,2,3)*dkmy + d1ql(3,3,3)*dkmy
        derv2(mzl,kzl) = derv2(mzl,kzl) + d1ql(3,2,3)*dkmz + d1ql(3,3,3)*dkmz
      endif
    endif
    if (lstr.and.lopm) then
!
!  Mixed strain-internal contribution
!
!  d2E/d(strain).dq x dq/d(alpha)
!
      do kl = 1,nstrains
        derv3(mx,kl) = derv3(mx,kl) + d1qsl(kl,1)*dimx + d1qsl(kl,2)*djmx + d1qsl(kl,3)*dkmx
        derv3(my,kl) = derv3(my,kl) + d1qsl(kl,1)*dimy + d1qsl(kl,2)*djmy + d1qsl(kl,3)*dkmy
        derv3(mz,kl) = derv3(mz,kl) + d1qsl(kl,1)*dimz + d1qsl(kl,2)*djmz + d1qsl(kl,3)*dkmz
      enddo
!
!  d2E/dq.dq x dq/d(alpha) x dq/d(strain)
!
      do kl = 1,nstrains
        derv3(mx,kl) = derv3(mx,kl) + d2ijs*(dqds(kl,i)*djmx + dqds(kl,j)*dimx)
        derv3(my,kl) = derv3(my,kl) + d2ijs*(dqds(kl,i)*djmy + dqds(kl,j)*dimy)
        derv3(mz,kl) = derv3(mz,kl) + d2ijs*(dqds(kl,i)*djmz + dqds(kl,j)*dimz)
        derv3(mx,kl) = derv3(mx,kl) + d2iks*(dqds(kl,i)*dkmx + dqds(kl,k)*dimx)
        derv3(my,kl) = derv3(my,kl) + d2iks*(dqds(kl,i)*dkmy + dqds(kl,k)*dimy)
        derv3(mz,kl) = derv3(mz,kl) + d2iks*(dqds(kl,i)*dkmz + dqds(kl,k)*dimz)
        derv3(mx,kl) = derv3(mx,kl) + d2jks*(dqds(kl,j)*dkmx + dqds(kl,k)*djmx)
        derv3(my,kl) = derv3(my,kl) + d2jks*(dqds(kl,j)*dkmy + dqds(kl,k)*djmy)
        derv3(mz,kl) = derv3(mz,kl) + d2jks*(dqds(kl,j)*dkmz + dqds(kl,k)*djmz)
      enddo
      if (abs(d2i2s).gt.1.0d-8) then
        do kl = 1,nstrains
          derv3(mx,kl) = derv3(mx,kl) + d2i2s*dqds(kl,i)*dimx
          derv3(my,kl) = derv3(my,kl) + d2i2s*dqds(kl,i)*dimy
          derv3(mz,kl) = derv3(mz,kl) + d2i2s*dqds(kl,i)*dimz
          derv3(ix,kl) = derv3(ix,kl) - d2i2s*dqds(kl,i)*dimx
          derv3(iy,kl) = derv3(iy,kl) - d2i2s*dqds(kl,i)*dimy
          derv3(iz,kl) = derv3(iz,kl) - d2i2s*dqds(kl,i)*dimz
        enddo
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        do kl = 1,nstrains
          derv3(mx,kl) = derv3(mx,kl) + d2j2s*dqds(kl,j)*djmx
          derv3(my,kl) = derv3(my,kl) + d2j2s*dqds(kl,j)*djmy
          derv3(mz,kl) = derv3(mz,kl) + d2j2s*dqds(kl,j)*djmz
          derv3(jx,kl) = derv3(jx,kl) - d2j2s*dqds(kl,j)*djmx
          derv3(jy,kl) = derv3(jy,kl) - d2j2s*dqds(kl,j)*djmy
          derv3(jz,kl) = derv3(jz,kl) - d2j2s*dqds(kl,j)*djmz
        enddo
      endif
      if (abs(d2k2s).gt.1.0d-8.and.i.ne.k.and.j.ne.k) then
        do kl = 1,nstrains
          derv3(mx,kl) = derv3(mx,kl) + d2k2s*dqds(kl,k)*dkmx
          derv3(my,kl) = derv3(my,kl) + d2k2s*dqds(kl,k)*dkmy
          derv3(mz,kl) = derv3(mz,kl) + d2k2s*dqds(kl,k)*dkmz
          derv3(kx,kl) = derv3(kx,kl) - d2k2s*dqds(kl,k)*dkmx
          derv3(ky,kl) = derv3(ky,kl) - d2k2s*dqds(kl,k)*dkmy
          derv3(kz,kl) = derv3(kz,kl) - d2k2s*dqds(kl,k)*dkmz
        enddo
      endif
    endif
!
!  Loop over atoms to add ij-ml correction
!
    lx = - 2
    ly = - 1
    lz =   0
    do l = 1,m-1
      lopl = (.not.lfreeze.or.lopf(l))
      indl = 3*(l-1)
      if (lopl) then
        lx = lx + 3
        ly = ly + 3
        lz = lz + 3
      endif
      if (lopm.or.lopl) then
        if (lopm.and.lopl) then
          mxl = mx
          myl = my
          mzl = mz
          lxl = lx
          lyl = ly
          lzl = lz
        elseif (lopm) then
          mxl = mx
          myl = my
          mzl = mz
          lxl = mx
          lyl = my
          lzl = mz
        elseif (lopl) then
          mxl = lx
          myl = ly
          mzl = lz
          lxl = lx
          lyl = ly
          lzl = lz
        endif
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
        dilx = dqdxyz(indl+1,i)
        dily = dqdxyz(indl+2,i)
        dilz = dqdxyz(indl+3,i)
        djlx = dqdxyz(indl+1,j)
        djly = dqdxyz(indl+2,j)
        djlz = dqdxyz(indl+3,j)
        dklx = dqdxyz(indl+1,k)
        dkly = dqdxyz(indl+2,k)
        dklz = dqdxyz(indl+3,k)
!
        derv2(lxl,mxl) = derv2(lxl,mxl) + d2ijs*dilx*djmx + d2iks*dilx*dkmx + d2jks*djlx*dkmx
        derv2(lyl,mxl) = derv2(lyl,mxl) + d2ijs*dilx*djmy + d2iks*dilx*dkmy + d2jks*djlx*dkmy
        derv2(lzl,mxl) = derv2(lzl,mxl) + d2ijs*dilx*djmz + d2iks*dilx*dkmz + d2jks*djlx*dkmz
        derv2(lxl,myl) = derv2(lxl,myl) + d2ijs*dily*djmx + d2iks*dily*dkmx + d2jks*djly*dkmx
        derv2(lyl,myl) = derv2(lyl,myl) + d2ijs*dily*djmy + d2iks*dily*dkmy + d2jks*djly*dkmy
        derv2(lzl,myl) = derv2(lzl,myl) + d2ijs*dily*djmz + d2iks*dily*dkmz + d2jks*djly*dkmz
        derv2(lxl,mzl) = derv2(lxl,mzl) + d2ijs*dilz*djmx + d2iks*dilz*dkmx + d2jks*djlz*dkmx
        derv2(lyl,mzl) = derv2(lyl,mzl) + d2ijs*dilz*djmy + d2iks*dilz*dkmy + d2jks*djlz*dkmy
        derv2(lzl,mzl) = derv2(lzl,mzl) + d2ijs*dilz*djmz + d2iks*dilz*dkmz + d2jks*djlz*dkmz
!
        derv2(lxl,mxl) = derv2(lxl,mxl) + d2ijs*djlx*dimx + d2iks*dklx*dimx + d2jks*dklx*djmx
        derv2(lyl,mxl) = derv2(lyl,mxl) + d2ijs*djlx*dimy + d2iks*dklx*dimy + d2jks*dklx*djmy
        derv2(lzl,mxl) = derv2(lzl,mxl) + d2ijs*djlx*dimz + d2iks*dklx*dimz + d2jks*dklx*djmz
        derv2(lxl,myl) = derv2(lxl,myl) + d2ijs*djly*dimx + d2iks*dkly*dimx + d2jks*dkly*djmx
        derv2(lyl,myl) = derv2(lyl,myl) + d2ijs*djly*dimy + d2iks*dkly*dimy + d2jks*dkly*djmy
        derv2(lzl,myl) = derv2(lzl,myl) + d2ijs*djly*dimz + d2iks*dkly*dimz + d2jks*dkly*djmz
        derv2(lxl,mzl) = derv2(lxl,mzl) + d2ijs*djlz*dimx + d2iks*dklz*dimx + d2jks*dklz*djmx
        derv2(lyl,mzl) = derv2(lyl,mzl) + d2ijs*djlz*dimy + d2iks*dklz*dimy + d2jks*dklz*djmy
        derv2(lzl,mzl) = derv2(lzl,mzl) + d2ijs*djlz*dimz + d2iks*dklz*dimz + d2jks*dklz*djmz
!
!  d2Edi2/d2Edj2/d2Edk2
!
        if (abs(d2i2s).gt.1.0d-8) then
          derv2(lxl,mxl) = derv2(lxl,mxl) + d2i2s*dilx*dimx
          derv2(lyl,mxl) = derv2(lyl,mxl) + d2i2s*dilx*dimy
          derv2(lzl,mxl) = derv2(lzl,mxl) + d2i2s*dilx*dimz
          derv2(lxl,myl) = derv2(lxl,myl) + d2i2s*dily*dimx
          derv2(lyl,myl) = derv2(lyl,myl) + d2i2s*dily*dimy
          derv2(lzl,myl) = derv2(lzl,myl) + d2i2s*dily*dimz
          derv2(lxl,mzl) = derv2(lxl,mzl) + d2i2s*dilz*dimx
          derv2(lyl,mzl) = derv2(lyl,mzl) + d2i2s*dilz*dimy
          derv2(lzl,mzl) = derv2(lzl,mzl) + d2i2s*dilz*dimz
        endif
        if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
          derv2(lxl,mxl) = derv2(lxl,mxl) + d2j2s*djlx*djmx
          derv2(lyl,mxl) = derv2(lyl,mxl) + d2j2s*djlx*djmy
          derv2(lzl,mxl) = derv2(lzl,mxl) + d2j2s*djlx*djmz
          derv2(lxl,myl) = derv2(lxl,myl) + d2j2s*djly*djmx
          derv2(lyl,myl) = derv2(lyl,myl) + d2j2s*djly*djmy
          derv2(lzl,myl) = derv2(lzl,myl) + d2j2s*djly*djmz
          derv2(lxl,mzl) = derv2(lxl,mzl) + d2j2s*djlz*djmx
          derv2(lyl,mzl) = derv2(lyl,mzl) + d2j2s*djlz*djmy
          derv2(lzl,mzl) = derv2(lzl,mzl) + d2j2s*djlz*djmz
        endif
        if (abs(d2k2s).gt.1.0d-8.and.i.ne.k.and.j.ne.k) then
          derv2(lxl,mxl) = derv2(lxl,mxl) + d2k2s*dklx*dkmx
          derv2(lyl,mxl) = derv2(lyl,mxl) + d2k2s*dklx*dkmy
          derv2(lzl,mxl) = derv2(lzl,mxl) + d2k2s*dklx*dkmz
          derv2(lxl,myl) = derv2(lxl,myl) + d2k2s*dkly*dkmx
          derv2(lyl,myl) = derv2(lyl,myl) + d2k2s*dkly*dkmy
          derv2(lzl,myl) = derv2(lzl,myl) + d2k2s*dkly*dkmz
          derv2(lxl,mzl) = derv2(lxl,mzl) + d2k2s*dklz*dkmx
          derv2(lyl,mzl) = derv2(lyl,mzl) + d2k2s*dklz*dkmy
          derv2(lzl,mzl) = derv2(lzl,mzl) + d2k2s*dklz*dkmz
        endif
      endif
!
!  End of loop over l
!
    enddo
!
!  End of loop over m
!
  enddo
!
!  Strain - charge derivatives terms
!
  if (lstr) then
!
!  Mixed strain-internal term
!
!  d2E/d(alpha).dq x dq/d(strain)
!
    if (lopi) then
      do kl = 1,nstrains
        derv3(ix,kl) = derv3(ix,kl) - (d1ql(1,1,1) + d1ql(1,2,1))*dqds(kl,i)
        derv3(iy,kl) = derv3(iy,kl) - (d1ql(2,1,1) + d1ql(2,2,1))*dqds(kl,i)
        derv3(iz,kl) = derv3(iz,kl) - (d1ql(3,1,1) + d1ql(3,2,1))*dqds(kl,i)
        derv3(ix,kl) = derv3(ix,kl) - (d1ql(1,1,2) + d1ql(1,2,2))*dqds(kl,j)
        derv3(iy,kl) = derv3(iy,kl) - (d1ql(2,1,2) + d1ql(2,2,2))*dqds(kl,j)
        derv3(iz,kl) = derv3(iz,kl) - (d1ql(3,1,2) + d1ql(3,2,2))*dqds(kl,j)
        derv3(ix,kl) = derv3(ix,kl) - (d1ql(1,1,3) + d1ql(1,2,3))*dqds(kl,k)
        derv3(iy,kl) = derv3(iy,kl) - (d1ql(2,1,3) + d1ql(2,2,3))*dqds(kl,k)
        derv3(iz,kl) = derv3(iz,kl) - (d1ql(3,1,3) + d1ql(3,2,3))*dqds(kl,k)
      enddo
    endif
    if (lopj) then
      do kl = 1,nstrains
        derv3(jx,kl) = derv3(jx,kl) + (d1ql(1,1,1) - d1ql(1,3,1))*dqds(kl,i)
        derv3(jy,kl) = derv3(jy,kl) + (d1ql(2,1,1) - d1ql(2,3,1))*dqds(kl,i)
        derv3(jz,kl) = derv3(jz,kl) + (d1ql(3,1,1) - d1ql(3,3,1))*dqds(kl,i)
        derv3(jx,kl) = derv3(jx,kl) + (d1ql(1,1,2) - d1ql(1,3,2))*dqds(kl,j)
        derv3(jy,kl) = derv3(jy,kl) + (d1ql(2,1,2) - d1ql(2,3,2))*dqds(kl,j)
        derv3(jz,kl) = derv3(jz,kl) + (d1ql(3,1,2) - d1ql(3,3,2))*dqds(kl,j)
        derv3(jx,kl) = derv3(jx,kl) + (d1ql(1,1,3) - d1ql(1,3,3))*dqds(kl,k)
        derv3(jy,kl) = derv3(jy,kl) + (d1ql(2,1,3) - d1ql(2,3,3))*dqds(kl,k)
        derv3(jz,kl) = derv3(jz,kl) + (d1ql(3,1,3) - d1ql(3,3,3))*dqds(kl,k)
      enddo
    endif
    if (lopk) then
      do kl = 1,nstrains
        derv3(kx,kl) = derv3(kx,kl) + (d1ql(1,2,1) + d1ql(1,3,1))*dqds(kl,i)
        derv3(ky,kl) = derv3(ky,kl) + (d1ql(2,2,1) + d1ql(2,3,1))*dqds(kl,i)
        derv3(kz,kl) = derv3(kz,kl) + (d1ql(3,2,1) + d1ql(3,3,1))*dqds(kl,i)
        derv3(kx,kl) = derv3(kx,kl) + (d1ql(1,2,2) + d1ql(1,3,2))*dqds(kl,j)
        derv3(ky,kl) = derv3(ky,kl) + (d1ql(2,2,2) + d1ql(2,3,2))*dqds(kl,j)
        derv3(kz,kl) = derv3(kz,kl) + (d1ql(3,2,2) + d1ql(3,3,2))*dqds(kl,j)
        derv3(kx,kl) = derv3(kx,kl) + (d1ql(1,2,3) + d1ql(1,3,3))*dqds(kl,k)
        derv3(ky,kl) = derv3(ky,kl) + (d1ql(2,2,3) + d1ql(2,3,3))*dqds(kl,k)
        derv3(kz,kl) = derv3(kz,kl) + (d1ql(3,2,3) + d1ql(3,3,3))*dqds(kl,k)
      enddo
    endif
!
!  Strain-strain terms for charge derivatives
!
    if (ndim.eq.3) then
      d1si = d1qsl(1,1)
      d2si = d1qsl(2,1)
      d3si = d1qsl(3,1)
      d4si = d1qsl(4,1)
      d5si = d1qsl(5,1)
      d6si = d1qsl(6,1)
      d1sj = d1qsl(1,2)
      d2sj = d1qsl(2,2)
      d3sj = d1qsl(3,2)
      d4sj = d1qsl(4,2)
      d5sj = d1qsl(5,2)
      d6sj = d1qsl(6,2)
      d1sk = d1qsl(1,3)
      d2sk = d1qsl(2,3)
      d3sk = d1qsl(3,3)
      d4sk = d1qsl(4,3)
      d5sk = d1qsl(5,3)
      d6sk = d1qsl(6,3)
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1 + d1sk*dks1)
      sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2 + d1sk*dks2
      sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1 + d2sk*dks1
      sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3 + d1sk*dks3
      sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1 + d3sk*dks1
      sderv2(4,1) = sderv2(4,1) + d1si*dis4 + d1sj*djs4 + d1sk*dks4
      sderv2(4,1) = sderv2(4,1) + d4si*dis1 + d4sj*djs1 + d4sk*dks1
      sderv2(5,1) = sderv2(5,1) + d1si*dis5 + d1sj*djs5 + d1sk*dks5
      sderv2(5,1) = sderv2(5,1) + d5si*dis1 + d5sj*djs1 + d5sk*dks1
      sderv2(6,1) = sderv2(6,1) + d1si*dis6 + d1sj*djs6 + d1sk*dks6
      sderv2(6,1) = sderv2(6,1) + d6si*dis1 + d6sj*djs1 + d6sk*dks1
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2 + d2sj*djs2 + d2sk*dks2)
      sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3 + d2sk*dks3
      sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2 + d3sk*dks2
      sderv2(4,2) = sderv2(4,2) + d2si*dis4 + d2sj*djs4 + d2sk*dks4
      sderv2(4,2) = sderv2(4,2) + d4si*dis2 + d4sj*djs2 + d4sk*dks2
      sderv2(5,2) = sderv2(5,2) + d2si*dis5 + d2sj*djs5 + d2sk*dks5
      sderv2(5,2) = sderv2(5,2) + d5si*dis2 + d5sj*djs2 + d5sk*dks2
      sderv2(6,2) = sderv2(6,2) + d2si*dis6 + d2sj*djs6 + d2sk*dks6
      sderv2(6,2) = sderv2(6,2) + d6si*dis2 + d6sj*djs2 + d6sk*dks2
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3 + d3sk*dks3)
      sderv2(4,3) = sderv2(4,3) + d3si*dis4 + d3sj*djs4 + d3sk*dks4
      sderv2(4,3) = sderv2(4,3) + d4si*dis3 + d4sj*djs3 + d4sk*dks3
      sderv2(5,3) = sderv2(5,3) + d3si*dis5 + d3sj*djs5 + d3sk*dks5
      sderv2(5,3) = sderv2(5,3) + d5si*dis3 + d5sj*djs3 + d5sk*dks3
      sderv2(6,3) = sderv2(6,3) + d3si*dis6 + d3sj*djs6 + d3sk*dks6
      sderv2(6,3) = sderv2(6,3) + d6si*dis3 + d6sj*djs3 + d6sk*dks3
      sderv2(4,4) = sderv2(4,4) + 2.0_dp*(d4si*dis4 + d4sj*djs4 + d4sk*dks4)
      sderv2(5,4) = sderv2(5,4) + d4si*dis5 + d4sj*djs5 + d4sk*dks5
      sderv2(5,4) = sderv2(5,4) + d5si*dis4 + d5sj*djs4 + d5sk*dks4
      sderv2(6,4) = sderv2(6,4) + d4si*dis6 + d4sj*djs6 + d4sk*dks6
      sderv2(6,4) = sderv2(6,4) + d6si*dis4 + d6sj*djs4 + d6sk*dks4
      sderv2(5,5) = sderv2(5,5) + 2.0_dp*(d5si*dis5 + d5sj*djs5 + d5sk*dks5)
      sderv2(6,5) = sderv2(6,5) + d5si*dis6 + d5sj*djs6 + d5sk*dks6
      sderv2(6,5) = sderv2(6,5) + d6si*dis5 + d6sj*djs5 + d6sk*dks5
      sderv2(6,6) = sderv2(6,6) + 2.0_dp*(d6si*dis6 + d6sj*djs6 + d6sk*dks6)
!
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
      sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
      sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
      sderv2(4,1) = sderv2(4,1) + d2ijs*(dis1*djs4 + dis4*djs1)
      sderv2(5,1) = sderv2(5,1) + d2ijs*(dis1*djs5 + dis5*djs1)
      sderv2(6,1) = sderv2(6,1) + d2ijs*(dis1*djs6 + dis6*djs1)
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
      sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
      sderv2(4,2) = sderv2(4,2) + d2ijs*(dis2*djs4 + dis4*djs2)
      sderv2(5,2) = sderv2(5,2) + d2ijs*(dis2*djs5 + dis5*djs2)
      sderv2(6,2) = sderv2(6,2) + d2ijs*(dis2*djs6 + dis6*djs2)
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
      sderv2(4,3) = sderv2(4,3) + d2ijs*(dis3*djs4 + dis4*djs3)
      sderv2(5,3) = sderv2(5,3) + d2ijs*(dis3*djs5 + dis5*djs3)
      sderv2(6,3) = sderv2(6,3) + d2ijs*(dis3*djs6 + dis6*djs3)
      sderv2(4,4) = sderv2(4,4) + 2.0_dp*d2ijs*dis4*djs4
      sderv2(5,4) = sderv2(5,4) + d2ijs*(dis4*djs5 + dis5*djs4)
      sderv2(6,4) = sderv2(6,4) + d2ijs*(dis4*djs6 + dis6*djs4)
      sderv2(5,5) = sderv2(5,5) + 2.0_dp*d2ijs*dis5*djs5
      sderv2(6,5) = sderv2(6,5) + d2ijs*(dis5*djs6 + dis6*djs5)
      sderv2(6,6) = sderv2(6,6) + 2.0_dp*d2ijs*dis6*djs6
!
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2iks*dis1*dks1
      sderv2(2,1) = sderv2(2,1) + d2iks*(dis1*dks2 + dis2*dks1)
      sderv2(3,1) = sderv2(3,1) + d2iks*(dis1*dks3 + dis3*dks1)
      sderv2(4,1) = sderv2(4,1) + d2iks*(dis1*dks4 + dis4*dks1)
      sderv2(5,1) = sderv2(5,1) + d2iks*(dis1*dks5 + dis5*dks1)
      sderv2(6,1) = sderv2(6,1) + d2iks*(dis1*dks6 + dis6*dks1)
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2iks*dis2*dks2
      sderv2(3,2) = sderv2(3,2) + d2iks*(dis2*dks3 + dis3*dks2)
      sderv2(4,2) = sderv2(4,2) + d2iks*(dis2*dks4 + dis4*dks2)
      sderv2(5,2) = sderv2(5,2) + d2iks*(dis2*dks5 + dis5*dks2)
      sderv2(6,2) = sderv2(6,2) + d2iks*(dis2*dks6 + dis6*dks2)
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2iks*dis3*dks3
      sderv2(4,3) = sderv2(4,3) + d2iks*(dis3*dks4 + dis4*dks3)
      sderv2(5,3) = sderv2(5,3) + d2iks*(dis3*dks5 + dis5*dks3)
      sderv2(6,3) = sderv2(6,3) + d2iks*(dis3*dks6 + dis6*dks3)
      sderv2(4,4) = sderv2(4,4) + 2.0_dp*d2iks*dis4*dks4
      sderv2(5,4) = sderv2(5,4) + d2iks*(dis4*dks5 + dis5*dks4)
      sderv2(6,4) = sderv2(6,4) + d2iks*(dis4*dks6 + dis6*dks4)
      sderv2(5,5) = sderv2(5,5) + 2.0_dp*d2iks*dis5*dks5
      sderv2(6,5) = sderv2(6,5) + d2iks*(dis5*dks6 + dis6*dks5)
      sderv2(6,6) = sderv2(6,6) + 2.0_dp*d2iks*dis6*dks6
!
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2jks*djs1*dks1
      sderv2(2,1) = sderv2(2,1) + d2jks*(djs1*dks2 + djs2*dks1)
      sderv2(3,1) = sderv2(3,1) + d2jks*(djs1*dks3 + djs3*dks1)
      sderv2(4,1) = sderv2(4,1) + d2jks*(djs1*dks4 + djs4*dks1)
      sderv2(5,1) = sderv2(5,1) + d2jks*(djs1*dks5 + djs5*dks1)
      sderv2(6,1) = sderv2(6,1) + d2jks*(djs1*dks6 + djs6*dks1)
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2jks*djs2*dks2
      sderv2(3,2) = sderv2(3,2) + d2jks*(djs2*dks3 + djs3*dks2)
      sderv2(4,2) = sderv2(4,2) + d2jks*(djs2*dks4 + djs4*dks2)
      sderv2(5,2) = sderv2(5,2) + d2jks*(djs2*dks5 + djs5*dks2)
      sderv2(6,2) = sderv2(6,2) + d2jks*(djs2*dks6 + djs6*dks2)
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2jks*djs3*dks3
      sderv2(4,3) = sderv2(4,3) + d2jks*(djs3*dks4 + djs4*dks3)
      sderv2(5,3) = sderv2(5,3) + d2jks*(djs3*dks5 + djs5*dks3)
      sderv2(6,3) = sderv2(6,3) + d2jks*(djs3*dks6 + djs6*dks3)
      sderv2(4,4) = sderv2(4,4) + 2.0_dp*d2jks*djs4*dks4
      sderv2(5,4) = sderv2(5,4) + d2jks*(djs4*dks5 + djs5*dks4)
      sderv2(6,4) = sderv2(6,4) + d2jks*(djs4*dks6 + djs6*dks4)
      sderv2(5,5) = sderv2(5,5) + 2.0_dp*d2jks*djs5*dks5
      sderv2(6,5) = sderv2(6,5) + d2jks*(djs5*dks6 + djs6*dks5)
      sderv2(6,6) = sderv2(6,6) + 2.0_dp*d2jks*djs6*dks6
!
      if (abs(d2i2s).gt.1.0d-8) then
        sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
        sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
        sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
        sderv2(4,1) = sderv2(4,1) + d2i2s*dis1*dis4
        sderv2(5,1) = sderv2(5,1) + d2i2s*dis1*dis5
        sderv2(6,1) = sderv2(6,1) + d2i2s*dis1*dis6
        sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
        sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
        sderv2(4,2) = sderv2(4,2) + d2i2s*dis2*dis4
        sderv2(5,2) = sderv2(5,2) + d2i2s*dis2*dis5
        sderv2(6,2) = sderv2(6,2) + d2i2s*dis2*dis6
        sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
        sderv2(4,3) = sderv2(4,3) + d2i2s*dis3*dis4
        sderv2(5,3) = sderv2(5,3) + d2i2s*dis3*dis5
        sderv2(6,3) = sderv2(6,3) + d2i2s*dis3*dis6
        sderv2(4,4) = sderv2(4,4) + d2i2s*dis4*dis4
        sderv2(5,4) = sderv2(5,4) + d2i2s*dis4*dis5
        sderv2(6,4) = sderv2(6,4) + d2i2s*dis4*dis6
        sderv2(5,5) = sderv2(5,5) + d2i2s*dis5*dis5
        sderv2(6,5) = sderv2(6,5) + d2i2s*dis5*dis6
        sderv2(6,6) = sderv2(6,6) + d2i2s*dis6*dis6
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
        sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
        sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
        sderv2(4,1) = sderv2(4,1) + d2j2s*djs1*djs4
        sderv2(5,1) = sderv2(5,1) + d2j2s*djs1*djs5
        sderv2(6,1) = sderv2(6,1) + d2j2s*djs1*djs6
        sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
        sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
        sderv2(4,2) = sderv2(4,2) + d2j2s*djs2*djs4
        sderv2(5,2) = sderv2(5,2) + d2j2s*djs2*djs5
        sderv2(6,2) = sderv2(6,2) + d2j2s*djs2*djs6
        sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
        sderv2(4,3) = sderv2(4,3) + d2j2s*djs3*djs4
        sderv2(5,3) = sderv2(5,3) + d2j2s*djs3*djs5
        sderv2(6,3) = sderv2(6,3) + d2j2s*djs3*djs6
        sderv2(4,4) = sderv2(4,4) + d2j2s*djs4*djs4
        sderv2(5,4) = sderv2(5,4) + d2j2s*djs4*djs5
        sderv2(6,4) = sderv2(6,4) + d2j2s*djs4*djs6
        sderv2(5,5) = sderv2(5,5) + d2j2s*djs5*djs5
        sderv2(6,5) = sderv2(6,5) + d2j2s*djs5*djs6
        sderv2(6,6) = sderv2(6,6) + d2j2s*djs6*djs6
      endif
      if (abs(d2k2s).gt.1.0d-8.and.i.ne.k.and.j.ne.k) then
        sderv2(1,1) = sderv2(1,1) + d2k2s*dks1*dks1
        sderv2(2,1) = sderv2(2,1) + d2k2s*dks1*dks2
        sderv2(3,1) = sderv2(3,1) + d2k2s*dks1*dks3
        sderv2(4,1) = sderv2(4,1) + d2k2s*dks1*dks4
        sderv2(5,1) = sderv2(5,1) + d2k2s*dks1*dks5
        sderv2(6,1) = sderv2(6,1) + d2k2s*dks1*dks6
        sderv2(2,2) = sderv2(2,2) + d2k2s*dks2*dks2
        sderv2(3,2) = sderv2(3,2) + d2k2s*dks2*dks3
        sderv2(4,2) = sderv2(4,2) + d2k2s*dks2*dks4
        sderv2(5,2) = sderv2(5,2) + d2k2s*dks2*dks5
        sderv2(6,2) = sderv2(6,2) + d2k2s*dks2*dks6
        sderv2(3,3) = sderv2(3,3) + d2k2s*dks3*dks3
        sderv2(4,3) = sderv2(4,3) + d2k2s*dks3*dks4
        sderv2(5,3) = sderv2(5,3) + d2k2s*dks3*dks5
        sderv2(6,3) = sderv2(6,3) + d2k2s*dks3*dks6
        sderv2(4,4) = sderv2(4,4) + d2k2s*dks4*dks4
        sderv2(5,4) = sderv2(5,4) + d2k2s*dks4*dks5
        sderv2(6,4) = sderv2(6,4) + d2k2s*dks4*dks6
        sderv2(5,5) = sderv2(5,5) + d2k2s*dks5*dks5
        sderv2(6,5) = sderv2(6,5) + d2k2s*dks5*dks6
        sderv2(6,6) = sderv2(6,6) + d2k2s*dks6*dks6
      endif
    elseif (ndim.eq.2) then
      d1si = d1qsl(1,1)
      d2si = d1qsl(2,1)
      d3si = d1qsl(3,1)
      d1sj = d1qsl(1,2)
      d2sj = d1qsl(2,2)
      d3sj = d1qsl(3,2)
      d1sk = d1qsl(1,3)
      d2sk = d1qsl(2,3)
      d3sk = d1qsl(3,3)
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1 + d1sk*dks1)
      sderv2(2,1) = sderv2(2,1) + d1si*dis2 + d1sj*djs2 + d1sk*dks2
      sderv2(2,1) = sderv2(2,1) + d2si*dis1 + d2sj*djs1 + d2sk*dks1
      sderv2(3,1) = sderv2(3,1) + d1si*dis3 + d1sj*djs3 + d1sk*dks3
      sderv2(3,1) = sderv2(3,1) + d3si*dis1 + d3sj*djs1 + d3sk*dks1
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*(d2si*dis2 + d2sj*djs2 + d2sk*dks2)
      sderv2(3,2) = sderv2(3,2) + d2si*dis3 + d2sj*djs3 + d2sk*dks3
      sderv2(3,2) = sderv2(3,2) + d3si*dis2 + d3sj*djs2 + d3sk*dks2
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*(d3si*dis3 + d3sj*djs3 + d3sk*dks3)
!
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
      sderv2(2,1) = sderv2(2,1) + d2ijs*(dis1*djs2 + dis2*djs1)
      sderv2(3,1) = sderv2(3,1) + d2ijs*(dis1*djs3 + dis3*djs1)
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2ijs*dis2*djs2
      sderv2(3,2) = sderv2(3,2) + d2ijs*(dis2*djs3 + dis3*djs2)
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2ijs*dis3*djs3
!
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2iks*dis1*dks1
      sderv2(2,1) = sderv2(2,1) + d2iks*(dis1*dks2 + dis2*dks1)
      sderv2(3,1) = sderv2(3,1) + d2iks*(dis1*dks3 + dis3*dks1)
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2iks*dis2*dks2
      sderv2(3,2) = sderv2(3,2) + d2iks*(dis2*dks3 + dis3*dks2)
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2iks*dis3*dks3
!
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2jks*djs1*dks1
      sderv2(2,1) = sderv2(2,1) + d2jks*(djs1*dks2 + djs2*dks1)
      sderv2(3,1) = sderv2(3,1) + d2jks*(djs1*dks3 + djs3*dks1)
      sderv2(2,2) = sderv2(2,2) + 2.0_dp*d2jks*djs2*dks2
      sderv2(3,2) = sderv2(3,2) + d2jks*(djs2*dks3 + dis3*dks2)
      sderv2(3,3) = sderv2(3,3) + 2.0_dp*d2jks*djs3*dks3
      if (abs(d2i2s).gt.1.0d-8) then
        sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
        sderv2(2,1) = sderv2(2,1) + d2i2s*dis1*dis2
        sderv2(3,1) = sderv2(3,1) + d2i2s*dis1*dis3
        sderv2(2,2) = sderv2(2,2) + d2i2s*dis2*dis2
        sderv2(3,2) = sderv2(3,2) + d2i2s*dis2*dis3
        sderv2(3,3) = sderv2(3,3) + d2i2s*dis3*dis3
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
        sderv2(2,1) = sderv2(2,1) + d2j2s*djs1*djs2
        sderv2(3,1) = sderv2(3,1) + d2j2s*djs1*djs3
        sderv2(2,2) = sderv2(2,2) + d2j2s*djs2*djs2
        sderv2(3,2) = sderv2(3,2) + d2j2s*djs2*djs3
        sderv2(3,3) = sderv2(3,3) + d2j2s*djs3*djs3
      endif
      if (abs(d2k2s).gt.1.0d-8.and.i.ne.k.and.j.ne.k) then
        sderv2(1,1) = sderv2(1,1) + d2k2s*dks1*dks1
        sderv2(2,1) = sderv2(2,1) + d2k2s*dks1*dks2
        sderv2(3,1) = sderv2(3,1) + d2k2s*dks1*dks3
        sderv2(2,2) = sderv2(2,2) + d2k2s*dks2*dks2
        sderv2(3,2) = sderv2(3,2) + d2k2s*dks2*dks3
        sderv2(3,3) = sderv2(3,3) + d2k2s*dks3*dks3
      endif
    elseif (ndim.eq.1) then
      d1si = d1qsl(1,1)
      d1sj = d1qsl(1,2)
      d1sk = d1qsl(1,3)
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*(d1si*dis1 + d1sj*djs1 + d1sk*dks1)
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2ijs*dis1*djs1
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2iks*dis1*dks1
      sderv2(1,1) = sderv2(1,1) + 2.0_dp*d2jks*djs1*dks1
      if (abs(d2i2s).gt.1.0d-8) then
        sderv2(1,1) = sderv2(1,1) + d2i2s*dis1*dis1
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        sderv2(1,1) = sderv2(1,1) + d2j2s*djs1*djs1
      endif
      if (abs(d2k2s).gt.1.0d-8.and.i.ne.k.and.j.ne.k) then
        sderv2(1,1) = sderv2(1,1) + d2k2s*dks1*dks1
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('d2charge3')
#endif
!
  return
  end

  subroutine d2chargep(i,j,nor,ix,iy,iz,jx,jy,jz,xkv,ykv,zkv,dei,dej, &
                       d1ixr,d1iyr,d1izr,d1jxr,d1jyr,d1jzr,d2i2r,d2ijr, &
                       d2j2r,d2self,dei0,dej0,lreal,lDoSelf,lupper,nbosptr)
!
!  Calculates the contribution to the dynamical matrix due to charge derivatives from EEM/QEq. 
!
!  At present this is gamma point only.
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!
!   2/98 Created from d2charge
!  10/02 ReaxFF modifications added
!  11/02 K vector now passed in for phasing
!   3/03 Frozen charge modifications added
!   9/04 Modified to generalise to charge derivatives other than from EEM
!   9/04 Contribution from d2q/dalpha.dbeta added
!   9/04 dei0/dej0 arguments added
!   9/04 xtmp/ytmp/ztmp premultiplied before entry to routine
!   7/05 Streitz and Mintmire modifications added
!  11/07 Unused variables removed
!   3/08 ReaxFF charge scheme removed as this is handled elsewhere
!   9/12 Pacha added
!   2/18 Trace added
!   5/18 Multiple qranges added
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   2/21 Modifications added for GFNFF
!   3/21 Corrections to addition of terms to derv2 for some terms
!   6/21 Exact d2q terms rearranged to speed up
!   7/21 Timing added
!   7/21 New algorithm to avoid 4th power scaling added
!  10/21 lgfnff moved to control
!   1/22 lDoSelf added to arguments
!  12/22 Optional argument, nbosptr added
!   5/23 lupper added as argument
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
!  Julian Gale, CIC, Curtin University, May 2023
!
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use reaxFFdata,     only : reaxFFmu
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)           :: i
  integer(i4), intent(in)           :: ix
  integer(i4), intent(in)           :: iy
  integer(i4), intent(in)           :: iz
  integer(i4), intent(in)           :: j
  integer(i4), intent(in)           :: jx
  integer(i4), intent(in)           :: jy
  integer(i4), intent(in)           :: jz
  integer(i4), intent(in)           :: nor
  integer(i4), intent(in), optional :: nbosptr(*)       ! Only required for reaxFF calls - points to reaxFF species from atoms
  logical,     intent(in)           :: lreal
  logical,     intent(in)           :: lDoSelf
  logical,     intent(in)           :: lupper           ! Indicates whether second derivatives should be written to upper or lower triangle
  real(dp),    intent(in)           :: dei(*)
  real(dp),    intent(in)           :: dej(*)
  real(dp),    intent(in)           :: dei0
  real(dp),    intent(in)           :: dej0
  real(dp),    intent(in)           :: d1ixr
  real(dp),    intent(in)           :: d1iyr
  real(dp),    intent(in)           :: d1izr
  real(dp),    intent(in)           :: d1jxr
  real(dp),    intent(in)           :: d1jyr
  real(dp),    intent(in)           :: d1jzr
  real(dp),    intent(in)           :: d2i2r(*)
  real(dp),    intent(in)           :: d2ijr(*)
  real(dp),    intent(in)           :: d2j2r(*)
  real(dp),    intent(in)           :: d2self
  real(dp),    intent(in)           :: xkv
  real(dp),    intent(in)           :: ykv
  real(dp),    intent(in)           :: zkv
!
!  Local variables     
!
  integer(i4)                       :: iv    
  integer(i4)                       :: k
  integer(i4)                       :: kx
  integer(i4)                       :: ky 
  integer(i4)                       :: kz  
  integer(i4)                       :: l      
  integer(i4)                       :: lx     
  integer(i4)                       :: ly    
  integer(i4)                       :: lz    
  integer(i4)                       :: m
  integer(i4)                       :: mnxx
  integer(i4)                       :: mnxy
  integer(i4)                       :: mnxz
  integer(i4)                       :: mnyx
  integer(i4)                       :: mnyy
  integer(i4)                       :: mnyz
  integer(i4)                       :: mnzx
  integer(i4)                       :: mnzy
  integer(i4)                       :: mnzz
  integer(i4)                       :: mx
  integer(i4)                       :: my
  integer(i4)                       :: mz
  integer(i4)                       :: n
  integer(i4)                       :: nqr
  integer(i4)                       :: nx
  integer(i4)                       :: nk
  logical                           :: lNonIJQDeriv
  real(dp)                          :: d1ix
  real(dp)                          :: d1iy
  real(dp)                          :: d1iz
  real(dp)                          :: d1jx
  real(dp)                          :: d1jy
  real(dp)                          :: d1jz
  real(dp)                          :: d2i2s
  real(dp)                          :: d2ijs
  real(dp)                          :: d2j2s
  real(dp)                          :: d2qk 
  real(dp)                          :: deisum
  real(dp)                          :: dejsum
  real(dp)                          :: dikx
  real(dp)                          :: diky
  real(dp)                          :: dikz
  real(dp)                          :: d2ikx
  real(dp)                          :: d2iky
  real(dp)                          :: d2ikz
  real(dp)                          :: dilx
  real(dp)                          :: dily
  real(dp)                          :: dilz 
  real(dp)                          :: djkx  
  real(dp)                          :: djky  
  real(dp)                          :: djkz   
  real(dp)                          :: d2jkx
  real(dp)                          :: d2jky
  real(dp)                          :: d2jkz
  real(dp)                          :: djlx   
  real(dp)                          :: djly
  real(dp)                          :: djlz
  real(dp)                          :: dkix
  real(dp)                          :: dkiy
  real(dp)                          :: dkiz 
  real(dp)                          :: dkjx 
  real(dp)                          :: dkjy 
  real(dp)                          :: dkjz 
  real(dp)                          :: ock 
  real(dp)                          :: g_cpu_time
  real(dp)                          :: qlk 
  real(dp)                          :: time1
  real(dp)                          :: time2
  real(dp), parameter               :: tsqrt2pi = 0.797884560802866_dp
  real(dp)                          :: zetah0
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargep')
#endif
!
  time1 = g_cpu_time()
! 
!  For ReaxFF check that nbosptr is present
!   
  if (lreaxFF) then
    if (.not.present(nbosptr)) then
      call outerror('d2chargep called for reaxFF without passing nbosptr',0_i4)
      call stopnow('d2chargep')
    endif
  endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d1ix = d1ixr
  d1iy = d1iyr
  d1iz = d1izr
  d1jx = d1jxr
  d1jy = d1jyr
  d1jz = d1jzr
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d2ijs = d2ijs + d2ijr(iv)
    d2i2s = d2i2s + d2i2r(iv)
    d2j2s = d2j2s + d2j2r(iv)
  enddo
  if (leem) then
    if (nat(i).gt.maxele) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp 
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp 
    endif
  endif
!
  if (.not.loldd2q) then
!******************
!  New algorithm  *
!******************
    d2edqdq(j,i) = d2edqdq(j,i) + d2ijs
    d2edq2(i) = d2edq2(i) + d2i2s
    if (i.ne.j) then
      d2edqdq(i,j) = d2edqdq(i,j) + d2ijs
      d2edq2(j) = d2edq2(j) + d2j2s
    endif
  endif
!
!  The following flag causes the second derivatives of the charge distribution 
!  with respect to two atoms neither of which are the i/j of the energy term
!  to be calculated. At the moment this contribution is not used in any method
!  and so the flag is set to false for efficiency.
!
  lNonIJQDeriv = .false.
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff.and..not.lreaxFF) then
    deisum = dei0
    dejsum = dej0
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(i)
      k = nqatomptr(m,i)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (lupper) then
        if (i.gt.k) then
          derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,i)
          derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,i)
          derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,i)
          derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,i)
          derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,i)
          derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,i)
          derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,i)
          derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,i)
          derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,i)
        elseif (i.lt.k) then
          derv2(ix,kx) = derv2(ix,kx) + deisum*d2qdxyz2(mnxx,i)
          derv2(iy,kx) = derv2(iy,kx) + deisum*d2qdxyz2(mnxy,i)
          derv2(iz,kx) = derv2(iz,kx) + deisum*d2qdxyz2(mnxz,i)
          derv2(ix,ky) = derv2(ix,ky) + deisum*d2qdxyz2(mnxy,i)
          derv2(iy,ky) = derv2(iy,ky) + deisum*d2qdxyz2(mnyy,i)
          derv2(iz,ky) = derv2(iz,ky) + deisum*d2qdxyz2(mnyz,i)
          derv2(ix,kz) = derv2(ix,kz) + deisum*d2qdxyz2(mnxz,i)
          derv2(iy,kz) = derv2(iy,kz) + deisum*d2qdxyz2(mnyz,i)
          derv2(iz,kz) = derv2(iz,kz) + deisum*d2qdxyz2(mnzz,i)
        endif
      else
        if (i.gt.k) then
          derv2(ix,kx) = derv2(ix,kx) + deisum*d2qdxyz2(mnxx,i)
          derv2(iy,kx) = derv2(iy,kx) + deisum*d2qdxyz2(mnxy,i)
          derv2(iz,kx) = derv2(iz,kx) + deisum*d2qdxyz2(mnxz,i)
          derv2(ix,ky) = derv2(ix,ky) + deisum*d2qdxyz2(mnxy,i)
          derv2(iy,ky) = derv2(iy,ky) + deisum*d2qdxyz2(mnyy,i)
          derv2(iz,ky) = derv2(iz,ky) + deisum*d2qdxyz2(mnyz,i)
          derv2(ix,kz) = derv2(ix,kz) + deisum*d2qdxyz2(mnxz,i)
          derv2(iy,kz) = derv2(iy,kz) + deisum*d2qdxyz2(mnyz,i)
          derv2(iz,kz) = derv2(iz,kz) + deisum*d2qdxyz2(mnzz,i)
        elseif (i.lt.k) then
          derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,i)
          derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,i)
          derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,i)
          derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,i)
          derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,i)
          derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,i)
          derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,i)
          derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,i)
          derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,i)
        endif
      endif
    enddo
!
    do m = 1,nqatoms(j)
      k = nqatomptr(m,j)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (j.gt.k) then
        derv2(kx,jx) = derv2(kx,jx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(ky,jx) = derv2(ky,jx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(kz,jx) = derv2(kz,jx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(kx,jy) = derv2(kx,jy) + dejsum*d2qdxyz2(mnxy,j)
        derv2(ky,jy) = derv2(ky,jy) + dejsum*d2qdxyz2(mnyy,j)
        derv2(kz,jy) = derv2(kz,jy) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kx,jz) = derv2(kx,jz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(ky,jz) = derv2(ky,jz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kz,jz) = derv2(kz,jz) + dejsum*d2qdxyz2(mnzz,j)
      elseif (j.lt.k) then
        derv2(jx,kx) = derv2(jx,kx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(jy,kx) = derv2(jy,kx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jz,kx) = derv2(jz,kx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jx,ky) = derv2(jx,ky) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jy,ky) = derv2(jy,ky) + dejsum*d2qdxyz2(mnyy,j)
        derv2(jz,ky) = derv2(jz,ky) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jx,kz) = derv2(jx,kz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jy,kz) = derv2(jy,kz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jz,kz) = derv2(jz,kz) + dejsum*d2qdxyz2(mnzz,j)
      endif
    enddo
!
!  Terms involving Q derivatives for non-i/j case
!
    if (lNonIJQDeriv) then
      do m = 2,nqatoms(i)
        k = nqatomptr(m,i)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,i)
          if (l.ne.k) then
            lx = 3*(l - 1) + 1
            ly = lx + 1
            lz = ly + 1
            nx = 3*(n - 1) + 1
!
            mnxx = mx*(mx - 1)/2 + nx
            mnxy = mnxx + 1
            mnxz = mnxx + 2
            mnyx = my*(my - 1)/2 + nx
            mnyy = mnyx + 1
            mnyz = mnyx + 2
            mnzx = mz*(mz - 1)/2 + nx
            mnzy = mnzx + 1
            mnzz = mnzx + 2
!
            if (lupper) then
              if (k.gt.l) then
                derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,i)
                derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,i)
                derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,i)
                derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,i)
                derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,i)
                derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,i)
                derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,i)
                derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,i)
                derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,i)
              else
                derv2(kx,lx) = derv2(kx,lx) + deisum*d2qdxyz2(mnxx,i)
                derv2(ky,lx) = derv2(ky,lx) + deisum*d2qdxyz2(mnxy,i)
                derv2(kz,lx) = derv2(kz,lx) + deisum*d2qdxyz2(mnxz,i)
                derv2(kx,ly) = derv2(kx,ly) + deisum*d2qdxyz2(mnyx,i)
                derv2(ky,ly) = derv2(ky,ly) + deisum*d2qdxyz2(mnyy,i)
                derv2(kz,ly) = derv2(kz,ly) + deisum*d2qdxyz2(mnyz,i)
                derv2(kx,lz) = derv2(kx,lz) + deisum*d2qdxyz2(mnzx,i)
                derv2(ky,lz) = derv2(ky,lz) + deisum*d2qdxyz2(mnzy,i)
                derv2(kz,lz) = derv2(kz,lz) + deisum*d2qdxyz2(mnzz,i)
              endif
            else
              if (k.gt.l) then
                derv2(kx,lx) = derv2(kx,lx) + deisum*d2qdxyz2(mnxx,i)
                derv2(ky,lx) = derv2(ky,lx) + deisum*d2qdxyz2(mnxy,i)
                derv2(kz,lx) = derv2(kz,lx) + deisum*d2qdxyz2(mnxz,i)
                derv2(kx,ly) = derv2(kx,ly) + deisum*d2qdxyz2(mnyx,i)
                derv2(ky,ly) = derv2(ky,ly) + deisum*d2qdxyz2(mnyy,i)
                derv2(kz,ly) = derv2(kz,ly) + deisum*d2qdxyz2(mnyz,i)
                derv2(kx,lz) = derv2(kx,lz) + deisum*d2qdxyz2(mnzx,i)
                derv2(ky,lz) = derv2(ky,lz) + deisum*d2qdxyz2(mnzy,i)
                derv2(kz,lz) = derv2(kz,lz) + deisum*d2qdxyz2(mnzz,i)
              else
                derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,i)
                derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,i)
                derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,i)
                derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,i)
                derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,i)
                derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,i)
                derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,i)
                derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,i)
                derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,i)
              endif
            endif
          endif
        enddo
      enddo
!
      do m = 2,nqatoms(j)
        k = nqatomptr(m,j)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,j)
          if (l.ne.k) then
            lx = 3*(l - 1) + 1
            ly = lx + 1
            lz = ly + 1
            nx = 3*(n - 1) + 1
!
            mnxx = mx*(mx - 1)/2 + nx
            mnxy = mnxx + 1
            mnxz = mnxx + 2
            mnyx = my*(my - 1)/2 + nx
            mnyy = mnyx + 1
            mnyz = mnyx + 2
            mnzx = mz*(mz - 1)/2 + nx
            mnzy = mnzx + 1
            mnzz = mnzx + 2
!
            if (lupper) then
              if (k.gt.l) then
                derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,j)
                derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,j)
                derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,j)
                derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,j)
                derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,j)
                derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,j)
                derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,j)
                derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,j)
                derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,j)
              else
                derv2(kx,lx) = derv2(kx,lx) + deisum*d2qdxyz2(mnxx,j)
                derv2(ky,lx) = derv2(ky,lx) + deisum*d2qdxyz2(mnxy,j)
                derv2(kz,lx) = derv2(kz,lx) + deisum*d2qdxyz2(mnxz,j)
                derv2(kx,ly) = derv2(kx,ly) + deisum*d2qdxyz2(mnyx,j)
                derv2(ky,ly) = derv2(ky,ly) + deisum*d2qdxyz2(mnyy,j)
                derv2(kz,ly) = derv2(kz,ly) + deisum*d2qdxyz2(mnyz,j)
                derv2(kx,lz) = derv2(kx,lz) + deisum*d2qdxyz2(mnzx,j)
                derv2(ky,lz) = derv2(ky,lz) + deisum*d2qdxyz2(mnzy,j)
                derv2(kz,lz) = derv2(kz,lz) + deisum*d2qdxyz2(mnzz,j)
              endif
            else
              if (k.gt.l) then
                derv2(kx,lx) = derv2(kx,lx) + deisum*d2qdxyz2(mnxx,j)
                derv2(ky,lx) = derv2(ky,lx) + deisum*d2qdxyz2(mnxy,j)
                derv2(kz,lx) = derv2(kz,lx) + deisum*d2qdxyz2(mnxz,j)
                derv2(kx,ly) = derv2(kx,ly) + deisum*d2qdxyz2(mnyx,j)
                derv2(ky,ly) = derv2(ky,ly) + deisum*d2qdxyz2(mnyy,j)
                derv2(kz,ly) = derv2(kz,ly) + deisum*d2qdxyz2(mnyz,j)
                derv2(kx,lz) = derv2(kx,lz) + deisum*d2qdxyz2(mnzx,j)
                derv2(ky,lz) = derv2(ky,lz) + deisum*d2qdxyz2(mnzy,j)
                derv2(kz,lz) = derv2(kz,lz) + deisum*d2qdxyz2(mnzz,j)
              else
                derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,j)
                derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,j)
                derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,j)
                derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,j)
                derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,j)
                derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,j)
                derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,j)
                derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,j)
                derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,j)
              endif
            endif
          endif
        enddo
      enddo
    endif
  endif
!
!  Loop over atoms to add ij-k contributions
!
  kx = - 2
  ky = - 1
  kz =   0
  do k = 1,numat
    kx = kx + 3
    ky = ky + 3
    kz = kz + 3
    dikx = dqdxyz(kx,i)
    diky = dqdxyz(ky,i)
    dikz = dqdxyz(kz,i)
    djkx = dqdxyz(kx,j)
    djky = dqdxyz(ky,j)
    djkz = dqdxyz(kz,j)
    if (i.ne.k) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-k
!
      if (lupper) then
        if (k.le.i) then
          derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx - d1jx*djkx
          derv2(ky,ix) = derv2(ky,ix) - d1iy*dikx - d1jy*djkx
          derv2(kz,ix) = derv2(kz,ix) - d1iz*dikx - d1jz*djkx
          derv2(kx,iy) = derv2(kx,iy) - d1ix*diky - d1jx*djky
          derv2(ky,iy) = derv2(ky,iy) - d1iy*diky - d1jy*djky
          derv2(kz,iy) = derv2(kz,iy) - d1iz*diky - d1jz*djky
          derv2(kx,iz) = derv2(kx,iz) - d1ix*dikz - d1jx*djkz
          derv2(ky,iz) = derv2(ky,iz) - d1iy*dikz - d1jy*djkz
          derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz - d1jz*djkz
        else
          derv2(ix,kx) = derv2(ix,kx) - d1ix*dikx - d1jx*djkx
          derv2(iy,kx) = derv2(iy,kx) - d1iy*dikx - d1jy*djkx
          derv2(iz,kx) = derv2(iz,kx) - d1iz*dikx - d1jz*djkx
          derv2(ix,ky) = derv2(ix,ky) - d1ix*diky - d1jx*djky
          derv2(iy,ky) = derv2(iy,ky) - d1iy*diky - d1jy*djky
          derv2(iz,ky) = derv2(iz,ky) - d1iz*diky - d1jz*djky
          derv2(ix,kz) = derv2(ix,kz) - d1ix*dikz - d1jx*djkz
          derv2(iy,kz) = derv2(iy,kz) - d1iy*dikz - d1jy*djkz
          derv2(iz,kz) = derv2(iz,kz) - d1iz*dikz - d1jz*djkz
        endif
      else
        if (k.le.i) then
          derv2(ix,kx) = derv2(ix,kx) - d1ix*dikx - d1jx*djkx
          derv2(iy,kx) = derv2(iy,kx) - d1ix*diky - d1jx*djky
          derv2(iz,kx) = derv2(iz,kx) - d1ix*dikz - d1jx*djkz
          derv2(ix,ky) = derv2(ix,ky) - d1iy*dikx - d1jy*djkx
          derv2(iy,ky) = derv2(iy,ky) - d1iy*diky - d1jy*djky
          derv2(iz,ky) = derv2(iz,ky) - d1iy*dikz - d1jy*djkz
          derv2(ix,kz) = derv2(ix,kz) - d1iz*dikx - d1jz*djkx
          derv2(iy,kz) = derv2(iy,kz) - d1iz*diky - d1jz*djky
          derv2(iz,kz) = derv2(iz,kz) - d1iz*dikz - d1jz*djkz
        else
          derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx - d1jx*djkx
          derv2(ky,ix) = derv2(ky,ix) - d1ix*diky - d1jx*djky
          derv2(kz,ix) = derv2(kz,ix) - d1ix*dikz - d1jx*djkz
          derv2(kx,iy) = derv2(kx,iy) - d1iy*dikx - d1jy*djkx
          derv2(ky,iy) = derv2(ky,iy) - d1iy*diky - d1jy*djky
          derv2(kz,iy) = derv2(kz,iy) - d1iy*dikz - d1jy*djkz
          derv2(kx,iz) = derv2(kx,iz) - d1iz*dikx - d1jz*djkx
          derv2(ky,iz) = derv2(ky,iz) - d1iz*diky - d1jz*djky
          derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz - d1jz*djkz
        endif
      endif
    endif
    if (j.ne.k) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : j-k
!
      if (lupper) then
        if (k.le.j) then
          derv2(kx,jx) = derv2(kx,jx) + d1jx*djkx + d1ix*dikx
          derv2(ky,jx) = derv2(ky,jx) + d1jy*djkx + d1iy*dikx
          derv2(kz,jx) = derv2(kz,jx) + d1jz*djkx + d1iz*dikx
          derv2(kx,jy) = derv2(kx,jy) + d1jx*djky + d1ix*diky
          derv2(ky,jy) = derv2(ky,jy) + d1jy*djky + d1iy*diky
          derv2(kz,jy) = derv2(kz,jy) + d1jz*djky + d1iz*diky
          derv2(kx,jz) = derv2(kx,jz) + d1jx*djkz + d1ix*dikz
          derv2(ky,jz) = derv2(ky,jz) + d1jy*djkz + d1iy*dikz
          derv2(kz,jz) = derv2(kz,jz) + d1jz*djkz + d1iz*dikz
        else
          derv2(jx,kx) = derv2(jx,kx) + d1jx*djkx + d1ix*dikx
          derv2(jy,kx) = derv2(jy,kx) + d1jy*djkx + d1iy*dikx
          derv2(jz,kx) = derv2(jz,kx) + d1jz*djkx + d1iz*dikx
          derv2(jx,ky) = derv2(jx,ky) + d1jx*djky + d1ix*diky
          derv2(jy,ky) = derv2(jy,ky) + d1jy*djky + d1iy*diky
          derv2(jz,ky) = derv2(jz,ky) + d1jz*djky + d1iz*diky
          derv2(jx,kz) = derv2(jx,kz) + d1jx*djkz + d1ix*dikz
          derv2(jy,kz) = derv2(jy,kz) + d1jy*djkz + d1iy*dikz
          derv2(jz,kz) = derv2(jz,kz) + d1jz*djkz + d1iz*dikz
        endif
      else
        if (k.le.j) then
          derv2(jx,kx) = derv2(jx,kx) + d1jx*djkx + d1ix*dikx
          derv2(jy,kx) = derv2(jy,kx) + d1jx*djky + d1ix*diky
          derv2(jz,kx) = derv2(jz,kx) + d1jx*djkz + d1ix*dikz
          derv2(jx,ky) = derv2(jx,ky) + d1jy*djkx + d1iy*dikx
          derv2(jy,ky) = derv2(jy,ky) + d1jy*djky + d1iy*diky
          derv2(jz,ky) = derv2(jz,ky) + d1jy*djkz + d1iy*dikz
          derv2(jx,kz) = derv2(jx,kz) + d1jz*djkx + d1iz*dikx
          derv2(jy,kz) = derv2(jy,kz) + d1jz*djky + d1iz*diky
          derv2(jz,kz) = derv2(jz,kz) + d1jz*djkz + d1iz*dikz
        else
          derv2(kx,jx) = derv2(kx,jx) + d1jx*djkx + d1ix*dikx
          derv2(ky,jx) = derv2(ky,jx) + d1jx*djky + d1ix*diky
          derv2(kz,jx) = derv2(kz,jx) + d1jx*djkz + d1ix*dikz
          derv2(kx,jy) = derv2(kx,jy) + d1jy*djkx + d1iy*dikx
          derv2(ky,jy) = derv2(ky,jy) + d1jy*djky + d1iy*diky
          derv2(kz,jy) = derv2(kz,jy) + d1jy*djkz + d1iy*dikz
          derv2(kx,jz) = derv2(kx,jz) + d1jz*djkx + d1iz*dikx
          derv2(ky,jz) = derv2(ky,jz) + d1jz*djky + d1iz*diky
          derv2(kz,jz) = derv2(kz,jz) + d1jz*djkz + d1iz*dikz
        endif
      endif
    endif
    if (lreal.and.lDoSelf.and.i.ne.j) then
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
      nk = nat(k)
      ock = occuf(k)
      if (leem.and.nk.le.maxele) then
        if (lmultiqrange.and.neemrptr(k).ne.0) then
          nqr = nqrnow(neemrptr(k))
        else
          nqr = 1
        endif
        if (lqeq) then
          if (nk.ne.1) then
            d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
          else
            qlk = qf(k)
            zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
            d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp+2.0_dp*qlk/zetah0)
          endif
        else
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        endif
      elseif (lgfnff) then
        d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
      elseif (lreaxFF) then
        d2qk = 2.0_dp*ock*reaxFFmu(nbosptr(k))
      else
        d2qk = 0.0_dp
      endif
      dkix = dqdxyz(ix,k)
      dkiy = dqdxyz(iy,k)
      dkiz = dqdxyz(iz,k)
      dkjx = dqdxyz(jx,k)
      dkjy = dqdxyz(jy,k)
      dkjz = dqdxyz(jz,k)
      if (lupper) then
        if (i.ge.j) then
          derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
          derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
          derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
          derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
          derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
          derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
          derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
          derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
          derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
        else
          derv2(ix,jx) = derv2(ix,jx) + d2qk*dkix*dkjx
          derv2(iy,jx) = derv2(iy,jx) + d2qk*dkix*dkjy
          derv2(iz,jx) = derv2(iz,jx) + d2qk*dkix*dkjz
          derv2(ix,jy) = derv2(ix,jy) + d2qk*dkiy*dkjx
          derv2(iy,jy) = derv2(iy,jy) + d2qk*dkiy*dkjy
          derv2(iz,jy) = derv2(iz,jy) + d2qk*dkiy*dkjz
          derv2(ix,jz) = derv2(ix,jz) + d2qk*dkiz*dkjx
          derv2(iy,jz) = derv2(iy,jz) + d2qk*dkiz*dkjy
          derv2(iz,jz) = derv2(iz,jz) + d2qk*dkiz*dkjz
        endif
      else
        if (i.ge.j) then
          derv2(ix,jx) = derv2(ix,jx) + d2qk*dkix*dkjx
          derv2(iy,jx) = derv2(iy,jx) + d2qk*dkix*dkjy
          derv2(iz,jx) = derv2(iz,jx) + d2qk*dkix*dkjz
          derv2(ix,jy) = derv2(ix,jy) + d2qk*dkiy*dkjx
          derv2(iy,jy) = derv2(iy,jy) + d2qk*dkiy*dkjy
          derv2(iz,jy) = derv2(iz,jy) + d2qk*dkiy*dkjz
          derv2(ix,jz) = derv2(ix,jz) + d2qk*dkiz*dkjx
          derv2(iy,jz) = derv2(iy,jz) + d2qk*dkiz*dkjy
          derv2(iz,jz) = derv2(iz,jz) + d2qk*dkiz*dkjz
        else
          derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
          derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
          derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
          derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
          derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
          derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
          derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
          derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
          derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
        endif
      endif
    endif
    if (loldd2q) then
!******************
!  Old algorithm  *
!******************
!
!  Loop over atoms to add ij-kl correction
!
      lx = - 2
      ly = - 1
      lz =   0
!
!  Scale charge derivatives for k by d2ijs
!
      d2ikx = d2ijs*dikx
      d2iky = d2ijs*diky
      d2ikz = d2ijs*dikz
      d2jkx = d2ijs*djkx
      d2jky = d2ijs*djky
      d2jkz = d2ijs*djkz
!
      if (lupper) then
        do l = 1,k-1
          lx = lx + 3
          ly = ly + 3
          lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
          dilx = dqdxyz(lx,i)
          dily = dqdxyz(ly,i)
          dilz = dqdxyz(lz,i)
          djlx = dqdxyz(lx,j)
          djly = dqdxyz(ly,j)
          djlz = dqdxyz(lz,j)
!
          derv2(lx,kx) = derv2(lx,kx) + dilx*d2jkx + djlx*d2ikx
          derv2(ly,kx) = derv2(ly,kx) + dilx*d2jky + djlx*d2iky
          derv2(lz,kx) = derv2(lz,kx) + dilx*d2jkz + djlx*d2ikz
          derv2(lx,ky) = derv2(lx,ky) + dily*d2jkx + djly*d2ikx
          derv2(ly,ky) = derv2(ly,ky) + dily*d2jky + djly*d2iky
          derv2(lz,ky) = derv2(lz,ky) + dily*d2jkz + djly*d2ikz
          derv2(lx,kz) = derv2(lx,kz) + dilz*d2jkx + djlz*d2ikx
          derv2(ly,kz) = derv2(ly,kz) + dilz*d2jky + djlz*d2iky
          derv2(lz,kz) = derv2(lz,kz) + dilz*d2jkz + djlz*d2ikz
!
!  d2Edi2/d2Edj2 - only non-zero for QEq/H
!
          if (abs(d2i2s).gt.1.0d-8) then
            derv2(lx,kx) = derv2(lx,kx) + d2i2s*dilx*dikx
            derv2(ly,kx) = derv2(ly,kx) + d2i2s*dilx*diky
            derv2(lz,kx) = derv2(lz,kx) + d2i2s*dilx*dikz
            derv2(lx,ky) = derv2(lx,ky) + d2i2s*dily*dikx
            derv2(ly,ky) = derv2(ly,ky) + d2i2s*dily*diky
            derv2(lz,ky) = derv2(lz,ky) + d2i2s*dily*dikz
            derv2(lx,kz) = derv2(lx,kz) + d2i2s*dilz*dikx
            derv2(ly,kz) = derv2(ly,kz) + d2i2s*dilz*diky
            derv2(lz,kz) = derv2(lz,kz) + d2i2s*dilz*dikz
          endif
          if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
            derv2(lx,kx) = derv2(lx,kx) + d2j2s*djlx*djkx
            derv2(ly,kx) = derv2(ly,kx) + d2j2s*djlx*djky
            derv2(lz,kx) = derv2(lz,kx) + d2j2s*djlx*djkz
            derv2(lx,ky) = derv2(lx,ky) + d2j2s*djly*djkx
            derv2(ly,ky) = derv2(ly,ky) + d2j2s*djly*djky
            derv2(lz,ky) = derv2(lz,ky) + d2j2s*djly*djkz
            derv2(lx,kz) = derv2(lx,kz) + d2j2s*djlz*djkx
            derv2(ly,kz) = derv2(ly,kz) + d2j2s*djlz*djky
            derv2(lz,kz) = derv2(lz,kz) + d2j2s*djlz*djkz
          endif
!
!  End of loop over l
!
        enddo
      else
        do l = 1,k-1
          lx = lx + 3
          ly = ly + 3
          lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta)
!
          dilx = dqdxyz(lx,i)
          dily = dqdxyz(ly,i)
          dilz = dqdxyz(lz,i)
          djlx = dqdxyz(lx,j)
          djly = dqdxyz(ly,j)
          djlz = dqdxyz(lz,j)
!
          derv2(kx,lx) = derv2(kx,lx) + dilx*d2jkx + djlx*d2ikx
          derv2(ky,lx) = derv2(ky,lx) + dily*d2jkx + djly*d2ikx
          derv2(kz,lx) = derv2(kz,lx) + dilz*d2jkx + djlz*d2ikx
          derv2(kx,ly) = derv2(kx,ly) + dilx*d2jky + djlx*d2iky
          derv2(ky,ly) = derv2(ky,ly) + dily*d2jky + djly*d2iky
          derv2(kz,ly) = derv2(kz,ly) + dilz*d2jky + djlz*d2iky
          derv2(kx,lz) = derv2(kx,lz) + dilx*d2jkz + djlx*d2ikz
          derv2(ky,lz) = derv2(ky,lz) + dily*d2jkz + djly*d2ikz
          derv2(kz,lz) = derv2(kz,lz) + dilz*d2jkz + djlz*d2ikz
!
!  d2Edi2/d2Edj2 - only non-zero for QEq/H
!
          if (abs(d2i2s).gt.1.0d-8) then
            derv2(kx,lx) = derv2(kx,lx) + d2i2s*dikx*dilx
            derv2(ky,lx) = derv2(ky,lx) + d2i2s*dikx*dily
            derv2(kz,lx) = derv2(kz,lx) + d2i2s*dikx*dilz
            derv2(kx,ly) = derv2(kx,ly) + d2i2s*diky*dilx
            derv2(ky,ly) = derv2(ky,ly) + d2i2s*diky*dily
            derv2(kz,ly) = derv2(kz,ly) + d2i2s*diky*dilz
            derv2(kx,lz) = derv2(kx,lz) + d2i2s*dikz*dilx
            derv2(ky,lz) = derv2(ky,lz) + d2i2s*dikz*dily
            derv2(kz,lz) = derv2(kz,lz) + d2i2s*dikz*dilz
          endif
          if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
            derv2(kx,lx) = derv2(kx,lx) + d2j2s*djkx*djlx
            derv2(ky,lx) = derv2(ky,lx) + d2j2s*djkx*djly
            derv2(kz,lx) = derv2(kz,lx) + d2j2s*djkx*djlz
            derv2(kx,ly) = derv2(kx,ly) + d2j2s*djky*djlx
            derv2(ky,ly) = derv2(ky,ly) + d2j2s*djky*djly
            derv2(kz,ly) = derv2(kz,ly) + d2j2s*djky*djlz
            derv2(kx,lz) = derv2(kx,lz) + d2j2s*djkz*djlx
            derv2(ky,lz) = derv2(ky,lz) + d2j2s*djkz*djly
            derv2(kz,lz) = derv2(kz,lz) + d2j2s*djkz*djlz
          endif
!
!  End of loop over l
!
        enddo
      endif
    endif
!
!  End of loop over k
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargep')
#endif
!
  return
  end

  subroutine d2chargepd(iloc,i,j,nor,ix,iy,iz,jx,jy,jz,xkv,ykv,zkv,dei,dej, &
                        d1ixr,d1iyr,d1izr,d1jxr,d1jyr,d1jzr,d2i2r,d2ijr, &
                        d2j2r,d2self,dei0,dej0,lreal,lDoSelf,nbosptr)
!
!  Calculates the contribution to the dynamical matrix due to charge derivatives from EEM/QEq. 
!  This version is for calls from distributed memory routines.
!
!  At present this is gamma point only.
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!
!   1/22 Created from d2chargep and d2charged
!   1/22 lDoSelf added to d2chargepd arguments
!  12/22 Optional argument, nbosptr added
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
!  Copyright Curtin University 2022
!
!  Julian Gale, CIC, Curtin University, December 2022
!
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use reaxFFdata,     only : reaxFFmu
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)           :: i       ! Atom i
  integer(i4), intent(in)           :: iloc    ! Atom i number on local node
  integer(i4), intent(in)           :: ix
  integer(i4), intent(in)           :: iy
  integer(i4), intent(in)           :: iz
  integer(i4), intent(in)           :: j
  integer(i4), intent(in)           :: jx
  integer(i4), intent(in)           :: jy
  integer(i4), intent(in)           :: jz
  integer(i4), intent(in)           :: nor
  integer(i4), intent(in), optional :: nbosptr(*)       ! Only required for reaxFF calls - points to reaxFF species from atoms
  logical,     intent(in)           :: lreal
  logical,     intent(in)           :: lDoSelf
  real(dp),    intent(in)           :: dei(*)
  real(dp),    intent(in)           :: dej(*)
  real(dp),    intent(in)           :: dei0
  real(dp),    intent(in)           :: dej0
  real(dp),    intent(in)           :: d1ixr
  real(dp),    intent(in)           :: d1iyr
  real(dp),    intent(in)           :: d1izr
  real(dp),    intent(in)           :: d1jxr
  real(dp),    intent(in)           :: d1jyr
  real(dp),    intent(in)           :: d1jzr
  real(dp),    intent(in)           :: d2i2r(*)
  real(dp),    intent(in)           :: d2ijr(*)
  real(dp),    intent(in)           :: d2j2r(*)
  real(dp),    intent(in)           :: d2self
  real(dp),    intent(in)           :: xkv
  real(dp),    intent(in)           :: ykv
  real(dp),    intent(in)           :: zkv
!
!  Local variables     
!
  integer(i4)                       :: iv    
  integer(i4)                       :: ixf
  integer(i4)                       :: iyf
  integer(i4)                       :: izf
  integer(i4)                       :: k
  integer(i4)                       :: kx
  integer(i4)                       :: ky 
  integer(i4)                       :: kz  
  integer(i4)                       :: m
  integer(i4)                       :: mnxx
  integer(i4)                       :: mnxy
  integer(i4)                       :: mnxz
  integer(i4)                       :: mnyy
  integer(i4)                       :: mnyz
  integer(i4)                       :: mnzz
  integer(i4)                       :: mx
  integer(i4)                       :: my
  integer(i4)                       :: mz
  integer(i4)                       :: nqr
  integer(i4)                       :: nk
  real(dp)                          :: d1ix
  real(dp)                          :: d1iy
  real(dp)                          :: d1iz
  real(dp)                          :: d1jx
  real(dp)                          :: d1jy
  real(dp)                          :: d1jz
  real(dp)                          :: d2i2s
  real(dp)                          :: d2ijs
  real(dp)                          :: d2j2s
  real(dp)                          :: d2qk 
  real(dp)                          :: deisum
  real(dp)                          :: dejsum
  real(dp)                          :: dkix
  real(dp)                          :: dkiy
  real(dp)                          :: dkiz 
  real(dp)                          :: dkjx 
  real(dp)                          :: dkjy 
  real(dp)                          :: dkjz 
  real(dp)                          :: ock 
  real(dp)                          :: g_cpu_time
  real(dp)                          :: qlk 
  real(dp)                          :: time1
  real(dp)                          :: time2
  real(dp), parameter               :: tsqrt2pi = 0.797884560802866_dp
  real(dp)                          :: zetah0
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargepd')
#endif
!
  time1 = g_cpu_time()
! 
!  For ReaxFF check that nbosptr is present
!   
  if (lreaxFF) then
    if (.not.present(nbosptr)) then
      call outerror('d2chargepd called for reaxFF without passing nbosptr',0_i4)
      call stopnow('d2chargepd')
    endif
  endif
!
  ixf = 3*(i-1) + 1
  iyf = ixf + 1
  izf = ixf + 2
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d1ix = d1ixr
  d1iy = d1iyr
  d1iz = d1izr
  d1jx = d1jxr
  d1jy = d1jyr
  d1jz = d1jzr
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d2ijs = d2ijs + d2ijr(iv)
    d2i2s = d2i2s + d2i2r(iv)
    d2j2s = d2j2s + d2j2r(iv)
  enddo
  if (leem) then
    if (nat(i).gt.maxele) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp 
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp 
    endif
  endif
!
!  Add d1ix / d1iy / d1iz to sums for globalisation later
!
  dedqc(1,i) = dedqc(1,i) + d1ix
  dedqc(2,i) = dedqc(2,i) + d1iy
  dedqc(3,i) = dedqc(3,i) + d1iz
!
!  Add d1jx / d1iy / d1jz to sums for globalisation later
!
  d2edqc(ixf,j) = d2edqc(ixf,j) + d1jx
  d2edqc(iyf,j) = d2edqc(iyf,j) + d1jy
  d2edqc(izf,j) = d2edqc(izf,j) + d1jz
!******************
!  New algorithm  *
!******************
  d2edqdq(j,iloc) = d2edqdq(j,iloc) + d2ijs
  d2edq2(i) = d2edq2(i) + d2i2s
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff.and..not.lreaxFF) then
    deisum = dei0
    dejsum = dej0
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(iloc)
      k = nqatomptr(m,iloc)
      if (k.ne.i) then
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
!
        mnxx = mx*(mx + 1)/2 
        mnyy = my*(my + 1)/2 
        mnzz = mz*(mz + 1)/2 
        mnxy = mnyy - 1
        mnxz = mnzz - 2
        mnyz = mnzz - 1
!
        derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,iloc)
        derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,iloc)
        derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,iloc)
        derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,iloc)
        derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,iloc)
        derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,iloc)
      endif
    enddo
  endif
!
!  Loop over atoms to add ij-k contributions
!
  kx = - 2
  ky = - 1
  kz =   0
  do k = 1,numat
    kx = kx + 3
    ky = ky + 3
    kz = kz + 3
!
    if (lreal.and.lDoSelf.and.i.ne.j) then
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
      nk = nat(k)
      ock = occuf(k)
      if (leem.and.nk.le.maxele) then
        if (lmultiqrange.and.neemrptr(k).ne.0) then
          nqr = nqrnow(neemrptr(k))
        else
          nqr = 1
        endif
        if (lqeq) then
          if (nk.ne.1) then
            d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
          else
            qlk = qf(k)
            zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
            d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp+2.0_dp*qlk/zetah0)
          endif
        else
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        endif
      elseif (lgfnff) then
        d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
      elseif (lreaxFF) then
        d2qk = 2.0_dp*ock*reaxFFmu(nbosptr(k))
      else
        d2qk = 0.0_dp
      endif
      dkix = dqdxyz(ixf,k)
      dkiy = dqdxyz(iyf,k)
      dkiz = dqdxyz(izf,k)
      dkjx = dqdxyz(jx,k)
      dkjy = dqdxyz(jy,k)
      dkjz = dqdxyz(jz,k)
      derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
      derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
      derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
      derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
      derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
      derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
      derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
      derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
      derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
    endif
!
!  End of loop over k
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargepd')
#endif
!
  return
  end

  subroutine d2charge3p(i,j,k,nor,ix,iy,iz,jx,jy,jz,kx,ky,kz,xkv,ykv,zkv,dei,dej,dek,d1q,d2q)
!
!  Calculates the contribution to the dynamical matrix due to charge derivatives from variable 
!  charge models as a result of three body potentials.
!
!  At present this is gamma point only. NB: Only one triangle of matrix is generated since
!  this is sufficient for the current algorithm at gamma.
!
!  10/04 Created from d2charge3
!  11/07 Unused variables removed
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   7/21 QEq terms corrected 
!   9/23 Trapping of shells added
!
!  On entry:
!
!  d1q(3,3,3) = derivative of first derivative of i - j vector with respect to charge
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
  use control
  use current
  use derivatives
  use element
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: i
  integer(i4), intent(in)    :: ix
  integer(i4), intent(in)    :: iy
  integer(i4), intent(in)    :: iz
  integer(i4), intent(in)    :: j
  integer(i4), intent(in)    :: jx
  integer(i4), intent(in)    :: jy
  integer(i4), intent(in)    :: jz
  integer(i4), intent(in)    :: k
  integer(i4), intent(in)    :: kx
  integer(i4), intent(in)    :: ky
  integer(i4), intent(in)    :: kz
  integer(i4), intent(in)    :: nor
  real(dp),    intent(in)    :: d1q(3,3,3)
  real(dp),    intent(in)    :: d2q(6,*)
  real(dp),    intent(in)    :: dei(*)
  real(dp),    intent(in)    :: dej(*)
  real(dp),    intent(in)    :: dek(*)
  real(dp),    intent(in)    :: xkv
  real(dp),    intent(in)    :: ykv
  real(dp),    intent(in)    :: zkv
!
!  Local variables
!
  integer(i4)                :: iv
  integer(i4)                :: l
  integer(i4)                :: lx
  integer(i4)                :: ly
  integer(i4)                :: lz
  integer(i4)                :: m
  integer(i4)                :: mnxx
  integer(i4)                :: mnxy
  integer(i4)                :: mnxz
  integer(i4)                :: mnyx
  integer(i4)                :: mnyy
  integer(i4)                :: mnyz
  integer(i4)                :: mnzx
  integer(i4)                :: mnzy
  integer(i4)                :: mnzz
  integer(i4)                :: mx
  integer(i4)                :: my
  integer(i4)                :: mz
  integer(i4)                :: n
  integer(i4)                :: nx
  integer(i4)                :: p
  integer(i4)                :: px
  integer(i4)                :: py
  integer(i4)                :: pz
  logical                    :: lNonIJQDeriv
  real(dp)                   :: d1ql(3,3,3)
  real(dp)                   :: d2i2s
  real(dp)                   :: d2ijs
  real(dp)                   :: d2iks
  real(dp)                   :: d2j2s
  real(dp)                   :: d2jks
  real(dp)                   :: d2k2s
  real(dp)                   :: deisum
  real(dp)                   :: dejsum
  real(dp)                   :: deksum
  real(dp)                   :: dilx
  real(dp)                   :: dily
  real(dp)                   :: dilz
  real(dp)                   :: dimx
  real(dp)                   :: dimy
  real(dp)                   :: dimz
  real(dp)                   :: djlx
  real(dp)                   :: djly
  real(dp)                   :: djlz
  real(dp)                   :: djmx
  real(dp)                   :: djmy
  real(dp)                   :: djmz
  real(dp)                   :: dklx
  real(dp)                   :: dkly
  real(dp)                   :: dklz
  real(dp)                   :: dkmx
  real(dp)                   :: dkmy
  real(dp)                   :: dkmz
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charge3p')
#endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d2i2s = 0.0_dp
  d2ijs = 0.0_dp
  d2iks = 0.0_dp
  d2j2s = 0.0_dp
  d2jks = 0.0_dp
  d2k2s = 0.0_dp
  do iv = 1,nor
    d2i2s = d2i2s + d2q(1,iv)
    d2ijs = d2ijs + d2q(2,iv)
    d2iks = d2iks + d2q(3,iv)
    d2j2s = d2j2s + d2q(4,iv)
    d2jks = d2jks + d2q(5,iv)
    d2k2s = d2k2s + d2q(6,iv)
  enddo
  d1ql(1:3,1:3,1:3) = d1q(1:3,1:3,1:3)
  if (leem) then
    if (nat(i).gt.maxele) then
      d2i2s = 0.0_dp
      d2ijs = 0.0_dp
      d2iks = 0.0_dp
      d1ql(1:3,1:3,1) = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2i2s = 0.0_dp
      d2ijs = 0.0_dp
      d2iks = 0.0_dp
      d1ql(1:3,1:3,1) = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d2jks = 0.0_dp
      d1ql(1:3,1:3,2) = 0.0_dp
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d2jks = 0.0_dp
      d1ql(1:3,1:3,2) = 0.0_dp
    endif
    if (nat(k).gt.maxele) then
      d2iks = 0.0_dp
      d2jks = 0.0_dp
      d2k2s = 0.0_dp
      d1ql(1:3,1:3,3) = 0.0_dp
    elseif (.not.lelementOK(nat(k)).or.nregionno(nrelf2a(k)).ne.1) then
      d2iks = 0.0_dp
      d2jks = 0.0_dp
      d2k2s = 0.0_dp
      d1ql(1:3,1:3,3) = 0.0_dp
    endif
  endif
!
!  The following flag causes the second derivatives of the charge distribution 
!  with respect to two atoms neither of which are the i/j of the energy term
!  to be calculated. At the moment this contribution is not used in any method
!  and so the flag is set to false for efficiency.
!
  lNonIJQDeriv = .false.
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lreaxFF) then
    deisum = 0.0_dp
    dejsum = 0.0_dp
    deksum = 0.0_dp
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
      deksum = deksum + dek(iv)
    enddo
!
!  Terms involving Q derivatives for i/j/k and one other atom
!
    do m = 1,nqatoms(i)
      p = nqatomptr(m,i)
      px = 3*(p - 1) + 1
      py = px + 1
      pz = py + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (i.gt.p) then
        derv2(px,ix) = derv2(px,ix) + deisum*d2qdxyz2(mnxx,i)
        derv2(py,ix) = derv2(py,ix) + deisum*d2qdxyz2(mnxy,i)
        derv2(pz,ix) = derv2(pz,ix) + deisum*d2qdxyz2(mnxz,i)
        derv2(px,iy) = derv2(px,iy) + deisum*d2qdxyz2(mnxy,i)
        derv2(py,iy) = derv2(py,iy) + deisum*d2qdxyz2(mnyy,i)
        derv2(pz,iy) = derv2(pz,iy) + deisum*d2qdxyz2(mnyz,i)
        derv2(px,iz) = derv2(px,iz) + deisum*d2qdxyz2(mnxz,i)
        derv2(py,iz) = derv2(py,iz) + deisum*d2qdxyz2(mnyz,i)
        derv2(pz,iz) = derv2(pz,iz) + deisum*d2qdxyz2(mnzz,i)
      elseif (i.lt.p) then
        derv2(ix,px) = derv2(ix,px) + deisum*d2qdxyz2(mnxx,i)
        derv2(iy,px) = derv2(iy,px) + deisum*d2qdxyz2(mnxy,i)
        derv2(iz,px) = derv2(iz,px) + deisum*d2qdxyz2(mnxz,i)
        derv2(ix,py) = derv2(ix,py) + deisum*d2qdxyz2(mnxy,i)
        derv2(iy,py) = derv2(iy,py) + deisum*d2qdxyz2(mnyy,i)
        derv2(iz,py) = derv2(iz,py) + deisum*d2qdxyz2(mnyz,i)
        derv2(ix,pz) = derv2(ix,pz) + deisum*d2qdxyz2(mnxz,i)
        derv2(iy,pz) = derv2(iy,pz) + deisum*d2qdxyz2(mnyz,i)
        derv2(iz,pz) = derv2(iz,pz) + deisum*d2qdxyz2(mnzz,i)
      endif
    enddo
!
    do m = 1,nqatoms(j)
      p = nqatomptr(m,j)
      px = 3*(p - 1) + 1
      py = px + 1
      pz = py + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (j.gt.p) then
        derv2(px,jx) = derv2(px,jx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(py,jx) = derv2(py,jx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(pz,jx) = derv2(pz,jx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(px,jy) = derv2(px,jy) + dejsum*d2qdxyz2(mnxy,j)
        derv2(py,jy) = derv2(py,jy) + dejsum*d2qdxyz2(mnyy,j)
        derv2(pz,jy) = derv2(pz,jy) + dejsum*d2qdxyz2(mnyz,j)
        derv2(px,jz) = derv2(px,jz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(py,jz) = derv2(py,jz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(pz,jz) = derv2(pz,jz) + dejsum*d2qdxyz2(mnzz,j)
      elseif (j.lt.p) then
        derv2(jx,px) = derv2(jx,px) + dejsum*d2qdxyz2(mnxx,j)
        derv2(jy,px) = derv2(jy,px) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jz,px) = derv2(jz,px) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jx,py) = derv2(jx,py) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jy,py) = derv2(jy,py) + dejsum*d2qdxyz2(mnyy,j)
        derv2(jz,py) = derv2(jz,py) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jx,pz) = derv2(jx,pz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jy,pz) = derv2(jy,pz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jz,pz) = derv2(jz,pz) + dejsum*d2qdxyz2(mnzz,j)
      endif
    enddo
!
    do m = 1,nqatoms(k)
      p = nqatomptr(m,k)
      px = 3*(p - 1) + 1
      py = px + 1
      pz = py + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (k.gt.p) then
        derv2(px,kx) = derv2(px,kx) + deksum*d2qdxyz2(mnxx,k)
        derv2(py,kx) = derv2(py,kx) + deksum*d2qdxyz2(mnxy,k)
        derv2(pz,kx) = derv2(pz,kx) + deksum*d2qdxyz2(mnxz,k)
        derv2(px,ky) = derv2(px,ky) + deksum*d2qdxyz2(mnxy,k)
        derv2(py,ky) = derv2(py,ky) + deksum*d2qdxyz2(mnyy,k)
        derv2(pz,ky) = derv2(pz,ky) + deksum*d2qdxyz2(mnyz,k)
        derv2(px,kz) = derv2(px,kz) + deksum*d2qdxyz2(mnxz,k)
        derv2(py,kz) = derv2(py,kz) + deksum*d2qdxyz2(mnyz,k)
        derv2(pz,kz) = derv2(pz,kz) + deksum*d2qdxyz2(mnzz,k)
      elseif (k.lt.p) then
        derv2(kx,px) = derv2(kx,px) + deksum*d2qdxyz2(mnxx,k)
        derv2(ky,px) = derv2(ky,px) + deksum*d2qdxyz2(mnxy,k)
        derv2(kz,px) = derv2(kz,px) + deksum*d2qdxyz2(mnxz,k)
        derv2(kx,py) = derv2(kx,py) + deksum*d2qdxyz2(mnxy,k)
        derv2(ky,py) = derv2(ky,py) + deksum*d2qdxyz2(mnyy,k)
        derv2(kz,py) = derv2(kz,py) + deksum*d2qdxyz2(mnyz,k)
        derv2(kx,pz) = derv2(kx,pz) + deksum*d2qdxyz2(mnxz,k)
        derv2(ky,pz) = derv2(ky,pz) + deksum*d2qdxyz2(mnyz,k)
        derv2(kz,pz) = derv2(kz,pz) + deksum*d2qdxyz2(mnzz,k)
      endif
    enddo
!
!  Terms involving Q derivatives for non-i/j case
!
    if (lNonIJQDeriv) then
      do m = 2,nqatoms(i)
        p = nqatomptr(m,i)
        px = 3*(p - 1) + 1
        py = px + 1
        pz = py + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,i)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,px) = derv2(lx,px) + deisum*d2qdxyz2(mnxx,i)
          derv2(ly,px) = derv2(ly,px) + deisum*d2qdxyz2(mnxy,i)
          derv2(lz,px) = derv2(lz,px) + deisum*d2qdxyz2(mnxz,i)
          derv2(lx,py) = derv2(lx,py) + deisum*d2qdxyz2(mnyx,i)
          derv2(ly,py) = derv2(ly,py) + deisum*d2qdxyz2(mnyy,i)
          derv2(lz,py) = derv2(lz,py) + deisum*d2qdxyz2(mnyz,i)
          derv2(lx,pz) = derv2(lx,pz) + deisum*d2qdxyz2(mnzx,i)
          derv2(ly,pz) = derv2(ly,pz) + deisum*d2qdxyz2(mnzy,i)
          derv2(lz,pz) = derv2(lz,pz) + deisum*d2qdxyz2(mnzz,i)
        enddo
      enddo
!
      do m = 2,nqatoms(j)
        p = nqatomptr(m,j)
        px = 3*(p - 1) + 1
        py = px + 1
        pz = py + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,j)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,px) = derv2(lx,px) + deisum*d2qdxyz2(mnxx,j)
          derv2(ly,px) = derv2(ly,px) + deisum*d2qdxyz2(mnxy,j)
          derv2(lz,px) = derv2(lz,px) + deisum*d2qdxyz2(mnxz,j)
          derv2(lx,py) = derv2(lx,py) + deisum*d2qdxyz2(mnyx,j)
          derv2(ly,py) = derv2(ly,py) + deisum*d2qdxyz2(mnyy,j)
          derv2(lz,py) = derv2(lz,py) + deisum*d2qdxyz2(mnyz,j)
          derv2(lx,pz) = derv2(lx,pz) + deisum*d2qdxyz2(mnzx,j)
          derv2(ly,pz) = derv2(ly,pz) + deisum*d2qdxyz2(mnzy,j)
          derv2(lz,pz) = derv2(lz,pz) + deisum*d2qdxyz2(mnzz,j)
        enddo
      enddo
!
      do m = 2,nqatoms(k)
        p = nqatomptr(m,k)
        px = 3*(p - 1) + 1
        py = px + 1
        pz = py + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,k)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,px) = derv2(lx,px) + deksum*d2qdxyz2(mnxx,k)
          derv2(ly,px) = derv2(ly,px) + deksum*d2qdxyz2(mnxy,k)
          derv2(lz,px) = derv2(lz,px) + deksum*d2qdxyz2(mnxz,k)
          derv2(lx,py) = derv2(lx,py) + deksum*d2qdxyz2(mnyx,k)
          derv2(ly,py) = derv2(ly,py) + deksum*d2qdxyz2(mnyy,k)
          derv2(lz,py) = derv2(lz,py) + deksum*d2qdxyz2(mnyz,k)
          derv2(lx,pz) = derv2(lx,pz) + deksum*d2qdxyz2(mnzx,k)
          derv2(ly,pz) = derv2(ly,pz) + deksum*d2qdxyz2(mnzy,k)
          derv2(lz,pz) = derv2(lz,pz) + deksum*d2qdxyz2(mnzz,k)
        enddo
      enddo
    endif
  endif
!
!  Loop over atoms to add ij-m contributions
!
  mx = - 2
  my = - 1
  mz =   0
  do m = 1,numat
    mx = mx + 3
    my = my + 3
    mz = mz + 3
!
    dimx = dqdxyz(mx,i)
    dimy = dqdxyz(my,i)
    dimz = dqdxyz(mz,i)
    djmx = dqdxyz(mx,j)
    djmy = dqdxyz(my,j)
    djmz = dqdxyz(mz,j)
    dkmx = dqdxyz(mx,k)
    dkmy = dqdxyz(my,k)
    dkmz = dqdxyz(mz,k)
    if (i.ne.m) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-m
!
      if (m.le.i) then
        derv2(mx,ix) = derv2(mx,ix) - d1ql(1,1,1)*dimx - d1ql(1,2,1)*dimx
        derv2(my,ix) = derv2(my,ix) - d1ql(1,1,1)*dimy - d1ql(1,2,1)*dimy
        derv2(mz,ix) = derv2(mz,ix) - d1ql(1,1,1)*dimz - d1ql(1,2,1)*dimz
        derv2(mx,iy) = derv2(mx,iy) - d1ql(2,1,1)*dimx - d1ql(2,2,1)*dimx
        derv2(my,iy) = derv2(my,iy) - d1ql(2,1,1)*dimy - d1ql(2,2,1)*dimy
        derv2(mz,iy) = derv2(mz,iy) - d1ql(2,1,1)*dimz - d1ql(2,2,1)*dimz
        derv2(mx,iz) = derv2(mx,iz) - d1ql(3,1,1)*dimx - d1ql(3,2,1)*dimx
        derv2(my,iz) = derv2(my,iz) - d1ql(3,1,1)*dimy - d1ql(3,2,1)*dimy
        derv2(mz,iz) = derv2(mz,iz) - d1ql(3,1,1)*dimz - d1ql(3,2,1)*dimz
!
        derv2(mx,ix) = derv2(mx,ix) - d1ql(1,1,2)*djmx - d1ql(1,2,2)*djmx
        derv2(my,ix) = derv2(my,ix) - d1ql(1,1,2)*djmy - d1ql(1,2,2)*djmy
        derv2(mz,ix) = derv2(mz,ix) - d1ql(1,1,2)*djmz - d1ql(1,2,2)*djmz
        derv2(mx,iy) = derv2(mx,iy) - d1ql(2,1,2)*djmx - d1ql(2,2,2)*djmx
        derv2(my,iy) = derv2(my,iy) - d1ql(2,1,2)*djmy - d1ql(2,2,2)*djmy
        derv2(mz,iy) = derv2(mz,iy) - d1ql(2,1,2)*djmz - d1ql(2,2,2)*djmz
        derv2(mx,iz) = derv2(mx,iz) - d1ql(3,1,2)*djmx - d1ql(3,2,2)*djmx
        derv2(my,iz) = derv2(my,iz) - d1ql(3,1,2)*djmy - d1ql(3,2,2)*djmy
        derv2(mz,iz) = derv2(mz,iz) - d1ql(3,1,2)*djmz - d1ql(3,2,2)*djmz
!
        derv2(mx,ix) = derv2(mx,ix) - d1ql(1,1,3)*dkmx - d1ql(1,2,3)*dkmx
        derv2(my,ix) = derv2(my,ix) - d1ql(1,1,3)*dkmy - d1ql(1,2,3)*dkmy
        derv2(mz,ix) = derv2(mz,ix) - d1ql(1,1,3)*dkmz - d1ql(1,2,3)*dkmz
        derv2(mx,iy) = derv2(mx,iy) - d1ql(2,1,3)*dkmx - d1ql(2,2,3)*dkmx
        derv2(my,iy) = derv2(my,iy) - d1ql(2,1,3)*dkmy - d1ql(2,2,3)*dkmy
        derv2(mz,iy) = derv2(mz,iy) - d1ql(2,1,3)*dkmz - d1ql(2,2,3)*dkmz
        derv2(mx,iz) = derv2(mx,iz) - d1ql(3,1,3)*dkmx - d1ql(3,2,3)*dkmx
        derv2(my,iz) = derv2(my,iz) - d1ql(3,1,3)*dkmy - d1ql(3,2,3)*dkmy
        derv2(mz,iz) = derv2(mz,iz) - d1ql(3,1,3)*dkmz - d1ql(3,2,3)*dkmz
      else
        derv2(ix,mx) = derv2(ix,mx) - d1ql(1,1,1)*dimx - d1ql(1,2,1)*dimx
        derv2(iy,mx) = derv2(iy,mx) - d1ql(2,1,1)*dimx - d1ql(2,2,1)*dimx
        derv2(iz,mx) = derv2(iz,mx) - d1ql(3,1,1)*dimx - d1ql(3,2,1)*dimx
        derv2(ix,my) = derv2(ix,my) - d1ql(1,1,1)*dimy - d1ql(1,2,1)*dimy
        derv2(iy,my) = derv2(iy,my) - d1ql(2,1,1)*dimy - d1ql(2,2,1)*dimy
        derv2(iz,my) = derv2(iz,my) - d1ql(3,1,1)*dimy - d1ql(3,2,1)*dimy
        derv2(ix,mz) = derv2(ix,mz) - d1ql(1,1,1)*dimz - d1ql(1,2,1)*dimz
        derv2(iy,mz) = derv2(iy,mz) - d1ql(2,1,1)*dimz - d1ql(2,2,1)*dimz
        derv2(iz,mz) = derv2(iz,mz) - d1ql(3,1,1)*dimz - d1ql(3,2,1)*dimz
!
        derv2(ix,mx) = derv2(ix,mx) - d1ql(1,1,2)*djmx - d1ql(1,2,2)*djmx
        derv2(iy,mx) = derv2(iy,mx) - d1ql(2,1,2)*djmx - d1ql(2,2,2)*djmx
        derv2(iz,mx) = derv2(iz,mx) - d1ql(3,1,2)*djmx - d1ql(3,2,2)*djmx
        derv2(ix,my) = derv2(ix,my) - d1ql(1,1,2)*djmy - d1ql(1,2,2)*djmy
        derv2(iy,my) = derv2(iy,my) - d1ql(2,1,2)*djmy - d1ql(2,2,2)*djmy
        derv2(iz,my) = derv2(iz,my) - d1ql(3,1,2)*djmy - d1ql(3,2,2)*djmy
        derv2(ix,mz) = derv2(ix,mz) - d1ql(1,1,2)*djmz - d1ql(1,2,2)*djmz
        derv2(iy,mz) = derv2(iy,mz) - d1ql(2,1,2)*djmz - d1ql(2,2,2)*djmz
        derv2(iz,mz) = derv2(iz,mz) - d1ql(3,1,2)*djmz - d1ql(3,2,2)*djmz
!
        derv2(ix,mx) = derv2(ix,mx) - d1ql(1,1,3)*dkmx - d1ql(1,2,3)*dkmx
        derv2(iy,mx) = derv2(iy,mx) - d1ql(2,1,3)*dkmx - d1ql(2,2,3)*dkmx
        derv2(iz,mx) = derv2(iz,mx) - d1ql(3,1,3)*dkmx - d1ql(3,2,3)*dkmx
        derv2(ix,my) = derv2(ix,my) - d1ql(1,1,3)*dkmy - d1ql(1,2,3)*dkmy
        derv2(iy,my) = derv2(iy,my) - d1ql(2,1,3)*dkmy - d1ql(2,2,3)*dkmy
        derv2(iz,my) = derv2(iz,my) - d1ql(3,1,3)*dkmy - d1ql(3,2,3)*dkmy
        derv2(ix,mz) = derv2(ix,mz) - d1ql(1,1,3)*dkmz - d1ql(1,2,3)*dkmz
        derv2(iy,mz) = derv2(iy,mz) - d1ql(2,1,3)*dkmz - d1ql(2,2,3)*dkmz
        derv2(iz,mz) = derv2(iz,mz) - d1ql(3,1,3)*dkmz - d1ql(3,2,3)*dkmz
      endif
    endif
    if (j.ne.m) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : j-m
!
      if (m.le.j) then
        derv2(mx,jx) = derv2(mx,jx) + d1ql(1,1,1)*dimx - d1ql(1,3,1)*dimx
        derv2(my,jx) = derv2(my,jx) + d1ql(1,1,1)*dimy - d1ql(1,3,1)*dimy
        derv2(mz,jx) = derv2(mz,jx) + d1ql(1,1,1)*dimz - d1ql(1,3,1)*dimz
        derv2(mx,jy) = derv2(mx,jy) + d1ql(2,1,1)*dimx - d1ql(2,3,1)*dimx
        derv2(my,jy) = derv2(my,jy) + d1ql(2,1,1)*dimy - d1ql(2,3,1)*dimy
        derv2(mz,jy) = derv2(mz,jy) + d1ql(2,1,1)*dimz - d1ql(2,3,1)*dimz
        derv2(mx,jz) = derv2(mx,jz) + d1ql(3,1,1)*dimx - d1ql(3,3,1)*dimx
        derv2(my,jz) = derv2(my,jz) + d1ql(3,1,1)*dimy - d1ql(3,3,1)*dimy
        derv2(mz,jz) = derv2(mz,jz) + d1ql(3,1,1)*dimz - d1ql(3,3,1)*dimz
!
        derv2(mx,jx) = derv2(mx,jx) + d1ql(1,1,2)*djmx - d1ql(1,3,2)*djmx
        derv2(my,jx) = derv2(my,jx) + d1ql(1,1,2)*djmy - d1ql(1,3,2)*djmy
        derv2(mz,jx) = derv2(mz,jx) + d1ql(1,1,2)*djmz - d1ql(1,3,2)*djmz
        derv2(mx,jy) = derv2(mx,jy) + d1ql(2,1,2)*djmx - d1ql(2,3,2)*djmx
        derv2(my,jy) = derv2(my,jy) + d1ql(2,1,2)*djmy - d1ql(2,3,2)*djmy
        derv2(mz,jy) = derv2(mz,jy) + d1ql(2,1,2)*djmz - d1ql(2,3,2)*djmz
        derv2(mx,jz) = derv2(mx,jz) + d1ql(3,1,2)*djmx - d1ql(3,3,2)*djmx
        derv2(my,jz) = derv2(my,jz) + d1ql(3,1,2)*djmy - d1ql(3,3,2)*djmy
        derv2(mz,jz) = derv2(mz,jz) + d1ql(3,1,2)*djmz - d1ql(3,3,2)*djmz
!
        derv2(mx,jx) = derv2(mx,jx) + d1ql(1,1,3)*dkmx - d1ql(1,3,3)*dkmx
        derv2(my,jx) = derv2(my,jx) + d1ql(1,1,3)*dkmy - d1ql(1,3,3)*dkmy
        derv2(mz,jx) = derv2(mz,jx) + d1ql(1,1,3)*dkmz - d1ql(1,3,3)*dkmz
        derv2(mx,jy) = derv2(mx,jy) + d1ql(2,1,3)*dkmx - d1ql(2,3,3)*dkmx
        derv2(my,jy) = derv2(my,jy) + d1ql(2,1,3)*dkmy - d1ql(2,3,3)*dkmy
        derv2(mz,jy) = derv2(mz,jy) + d1ql(2,1,3)*dkmz - d1ql(2,3,3)*dkmz
        derv2(mx,jz) = derv2(mx,jz) + d1ql(3,1,3)*dkmx - d1ql(3,3,3)*dkmx
        derv2(my,jz) = derv2(my,jz) + d1ql(3,1,3)*dkmy - d1ql(3,3,3)*dkmy
        derv2(mz,jz) = derv2(mz,jz) + d1ql(3,1,3)*dkmz - d1ql(3,3,3)*dkmz
      else
        derv2(jx,mx) = derv2(jx,mx) + d1ql(1,1,1)*dimx - d1ql(1,3,1)*dimx
        derv2(jy,mx) = derv2(jy,mx) + d1ql(2,1,1)*dimx - d1ql(2,3,1)*dimx
        derv2(jz,mx) = derv2(jz,mx) + d1ql(3,1,1)*dimx - d1ql(3,3,1)*dimx
        derv2(jx,my) = derv2(jx,my) + d1ql(1,1,1)*dimy - d1ql(1,3,1)*dimy
        derv2(jy,my) = derv2(jy,my) + d1ql(2,1,1)*dimy - d1ql(2,3,1)*dimy
        derv2(jz,my) = derv2(jz,my) + d1ql(3,1,1)*dimy - d1ql(3,3,1)*dimy
        derv2(jx,mz) = derv2(jx,mz) + d1ql(1,1,1)*dimz - d1ql(1,3,1)*dimz
        derv2(jy,mz) = derv2(jy,mz) + d1ql(2,1,1)*dimz - d1ql(2,3,1)*dimz
        derv2(jz,mz) = derv2(jz,mz) + d1ql(3,1,1)*dimz - d1ql(3,3,1)*dimz
!
        derv2(jx,mx) = derv2(jx,mx) + d1ql(1,1,2)*djmx - d1ql(1,3,2)*djmx
        derv2(jy,mx) = derv2(jy,mx) + d1ql(2,1,2)*djmx - d1ql(2,3,2)*djmx
        derv2(jz,mx) = derv2(jz,mx) + d1ql(3,1,2)*djmx - d1ql(3,3,2)*djmx
        derv2(jx,my) = derv2(jx,my) + d1ql(1,1,2)*djmy - d1ql(1,3,2)*djmy
        derv2(jy,my) = derv2(jy,my) + d1ql(2,1,2)*djmy - d1ql(2,3,2)*djmy
        derv2(jz,my) = derv2(jz,my) + d1ql(3,1,2)*djmy - d1ql(3,3,2)*djmy
        derv2(jx,mz) = derv2(jx,mz) + d1ql(1,1,2)*djmz - d1ql(1,3,2)*djmz
        derv2(jy,mz) = derv2(jy,mz) + d1ql(2,1,2)*djmz - d1ql(2,3,2)*djmz
        derv2(jz,mz) = derv2(jz,mz) + d1ql(3,1,2)*djmz - d1ql(3,3,2)*djmz
!
        derv2(jx,mx) = derv2(jx,mx) + d1ql(1,1,3)*dkmx - d1ql(1,3,3)*dkmx
        derv2(jy,mx) = derv2(jy,mx) + d1ql(2,1,3)*dkmx - d1ql(2,3,3)*dkmx
        derv2(jz,mx) = derv2(jz,mx) + d1ql(3,1,3)*dkmx - d1ql(3,3,3)*dkmx
        derv2(jx,my) = derv2(jx,my) + d1ql(1,1,3)*dkmy - d1ql(1,3,3)*dkmy
        derv2(jy,my) = derv2(jy,my) + d1ql(2,1,3)*dkmy - d1ql(2,3,3)*dkmy
        derv2(jz,my) = derv2(jz,my) + d1ql(3,1,3)*dkmy - d1ql(3,3,3)*dkmy
        derv2(jx,mz) = derv2(jx,mz) + d1ql(1,1,3)*dkmz - d1ql(1,3,3)*dkmz
        derv2(jy,mz) = derv2(jy,mz) + d1ql(2,1,3)*dkmz - d1ql(2,3,3)*dkmz
        derv2(jz,mz) = derv2(jz,mz) + d1ql(3,1,3)*dkmz - d1ql(3,3,3)*dkmz
      endif
    endif
    if (k.ne.m) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : k-m
!
      if (m.le.k) then
        derv2(mx,kx) = derv2(mx,kx) + d1ql(1,2,1)*dimx + d1ql(1,3,1)*dimx
        derv2(my,kx) = derv2(my,kx) + d1ql(1,2,1)*dimy + d1ql(1,3,1)*dimy
        derv2(mz,kx) = derv2(mz,kx) + d1ql(1,2,1)*dimz + d1ql(1,3,1)*dimz
        derv2(mx,ky) = derv2(mx,ky) + d1ql(2,2,1)*dimx + d1ql(2,3,1)*dimx
        derv2(my,ky) = derv2(my,ky) + d1ql(2,2,1)*dimy + d1ql(2,3,1)*dimy
        derv2(mz,ky) = derv2(mz,ky) + d1ql(2,2,1)*dimz + d1ql(2,3,1)*dimz
        derv2(mx,kz) = derv2(mx,kz) + d1ql(3,2,1)*dimx + d1ql(3,3,1)*dimx
        derv2(my,kz) = derv2(my,kz) + d1ql(3,2,1)*dimy + d1ql(3,3,1)*dimy
        derv2(mz,kz) = derv2(mz,kz) + d1ql(3,2,1)*dimz + d1ql(3,3,1)*dimz
!
        derv2(mx,kx) = derv2(mx,kx) + d1ql(1,2,2)*djmx + d1ql(1,3,2)*djmx
        derv2(my,kx) = derv2(my,kx) + d1ql(1,2,2)*djmy + d1ql(1,3,2)*djmy
        derv2(mz,kx) = derv2(mz,kx) + d1ql(1,2,2)*djmz + d1ql(1,3,2)*djmz
        derv2(mx,ky) = derv2(mx,ky) + d1ql(2,2,2)*djmx + d1ql(2,3,2)*djmx
        derv2(my,ky) = derv2(my,ky) + d1ql(2,2,2)*djmy + d1ql(2,3,2)*djmy
        derv2(mz,ky) = derv2(mz,ky) + d1ql(2,2,2)*djmz + d1ql(2,3,2)*djmz
        derv2(mx,kz) = derv2(mx,kz) + d1ql(3,2,2)*djmx + d1ql(3,3,2)*djmx
        derv2(my,kz) = derv2(my,kz) + d1ql(3,2,2)*djmy + d1ql(3,3,2)*djmy
        derv2(mz,kz) = derv2(mz,kz) + d1ql(3,2,2)*djmz + d1ql(3,3,2)*djmz
!
        derv2(mx,kx) = derv2(mx,kx) + d1ql(1,2,3)*dkmx + d1ql(1,3,3)*dkmx
        derv2(my,kx) = derv2(my,kx) + d1ql(1,2,3)*dkmy + d1ql(1,3,3)*dkmy
        derv2(mz,kx) = derv2(mz,kx) + d1ql(1,2,3)*dkmz + d1ql(1,3,3)*dkmz
        derv2(mx,ky) = derv2(mx,ky) + d1ql(2,2,3)*dkmx + d1ql(2,3,3)*dkmx
        derv2(my,ky) = derv2(my,ky) + d1ql(2,2,3)*dkmy + d1ql(2,3,3)*dkmy
        derv2(mz,ky) = derv2(mz,ky) + d1ql(2,2,3)*dkmz + d1ql(2,3,3)*dkmz
        derv2(mx,kz) = derv2(mx,kz) + d1ql(3,2,3)*dkmx + d1ql(3,3,3)*dkmx
        derv2(my,kz) = derv2(my,kz) + d1ql(3,2,3)*dkmy + d1ql(3,3,3)*dkmy
        derv2(mz,kz) = derv2(mz,kz) + d1ql(3,2,3)*dkmz + d1ql(3,3,3)*dkmz
      else
        derv2(kx,mx) = derv2(kx,mx) + d1ql(1,2,1)*dimx + d1ql(1,3,1)*dimx
        derv2(ky,mx) = derv2(ky,mx) + d1ql(2,2,1)*dimx + d1ql(2,3,1)*dimx
        derv2(kz,mx) = derv2(kz,mx) + d1ql(3,2,1)*dimx + d1ql(3,3,1)*dimx
        derv2(kx,my) = derv2(kx,my) + d1ql(1,2,1)*dimy + d1ql(1,3,1)*dimy
        derv2(ky,my) = derv2(ky,my) + d1ql(2,2,1)*dimy + d1ql(2,3,1)*dimy
        derv2(kz,my) = derv2(kz,my) + d1ql(3,2,1)*dimy + d1ql(3,3,1)*dimy
        derv2(kx,mz) = derv2(kx,mz) + d1ql(1,2,1)*dimz + d1ql(1,3,1)*dimz
        derv2(ky,mz) = derv2(ky,mz) + d1ql(2,2,1)*dimz + d1ql(2,3,1)*dimz
        derv2(kz,mz) = derv2(kz,mz) + d1ql(3,2,1)*dimz + d1ql(3,3,1)*dimz
!
        derv2(kx,mx) = derv2(kx,mx) + d1ql(1,2,2)*djmx + d1ql(1,3,2)*djmx
        derv2(ky,mx) = derv2(ky,mx) + d1ql(2,2,2)*djmx + d1ql(2,3,2)*djmx
        derv2(kz,mx) = derv2(kz,mx) + d1ql(3,2,2)*djmx + d1ql(3,3,2)*djmx
        derv2(kx,my) = derv2(kx,my) + d1ql(1,2,2)*djmy + d1ql(1,3,2)*djmy
        derv2(ky,my) = derv2(ky,my) + d1ql(2,2,2)*djmy + d1ql(2,3,2)*djmy
        derv2(kz,my) = derv2(kz,my) + d1ql(3,2,2)*djmy + d1ql(3,3,2)*djmy
        derv2(kx,mz) = derv2(kx,mz) + d1ql(1,2,2)*djmz + d1ql(1,3,2)*djmz
        derv2(ky,mz) = derv2(ky,mz) + d1ql(2,2,2)*djmz + d1ql(2,3,2)*djmz
        derv2(kz,mz) = derv2(kz,mz) + d1ql(3,2,2)*djmz + d1ql(3,3,2)*djmz
!
        derv2(kx,mx) = derv2(kx,mx) + d1ql(1,2,3)*dkmx + d1ql(1,3,3)*dkmx
        derv2(ky,mx) = derv2(ky,mx) + d1ql(2,2,3)*dkmx + d1ql(2,3,3)*dkmx
        derv2(kz,mx) = derv2(kz,mx) + d1ql(3,2,3)*dkmx + d1ql(3,3,3)*dkmx
        derv2(kx,my) = derv2(kx,my) + d1ql(1,2,3)*dkmy + d1ql(1,3,3)*dkmy
        derv2(ky,my) = derv2(ky,my) + d1ql(2,2,3)*dkmy + d1ql(2,3,3)*dkmy
        derv2(kz,my) = derv2(kz,my) + d1ql(3,2,3)*dkmy + d1ql(3,3,3)*dkmy
        derv2(kx,mz) = derv2(kx,mz) + d1ql(1,2,3)*dkmz + d1ql(1,3,3)*dkmz
        derv2(ky,mz) = derv2(ky,mz) + d1ql(2,2,3)*dkmz + d1ql(2,3,3)*dkmz
        derv2(kz,mz) = derv2(kz,mz) + d1ql(3,2,3)*dkmz + d1ql(3,3,3)*dkmz
      endif
    endif
!
!  Loop over atoms to add ij-ml correction
!
    lx = - 2
    ly = - 1
    lz =   0
    do l = 1,m-1
      lx = lx + 3
      ly = ly + 3
      lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
      dilx = dqdxyz(lx,i)
      dily = dqdxyz(ly,i)
      dilz = dqdxyz(lz,i)
      djlx = dqdxyz(lx,j)
      djly = dqdxyz(ly,j)
      djlz = dqdxyz(lz,j)
      dklx = dqdxyz(lx,k)
      dkly = dqdxyz(ly,k)
      dklz = dqdxyz(lz,k)
!
      derv2(lx,mx) = derv2(lx,mx) + d2ijs*dilx*djmx + d2iks*dilx*dkmx + d2jks*djlx*dkmx
      derv2(ly,mx) = derv2(ly,mx) + d2ijs*dilx*djmy + d2iks*dilx*dkmy + d2jks*djlx*dkmy
      derv2(lz,mx) = derv2(lz,mx) + d2ijs*dilx*djmz + d2iks*dilx*dkmz + d2jks*djlx*dkmz
      derv2(lx,my) = derv2(lx,my) + d2ijs*dily*djmx + d2iks*dily*dkmx + d2jks*djly*dkmx
      derv2(ly,my) = derv2(ly,my) + d2ijs*dily*djmy + d2iks*dily*dkmy + d2jks*djly*dkmy
      derv2(lz,my) = derv2(lz,my) + d2ijs*dily*djmz + d2iks*dily*dkmz + d2jks*djly*dkmz
      derv2(lx,mz) = derv2(lx,mz) + d2ijs*dilz*djmx + d2iks*dilz*dkmx + d2jks*djlz*dkmx
      derv2(ly,mz) = derv2(ly,mz) + d2ijs*dilz*djmy + d2iks*dilz*dkmy + d2jks*djlz*dkmy
      derv2(lz,mz) = derv2(lz,mz) + d2ijs*dilz*djmz + d2iks*dilz*dkmz + d2jks*djlz*dkmz
!
      derv2(lx,mx) = derv2(lx,mx) + d2ijs*djlx*dimx + d2iks*dklx*dimx + d2jks*dklx*djmx
      derv2(ly,mx) = derv2(ly,mx) + d2ijs*djlx*dimy + d2iks*dklx*dimy + d2jks*dklx*djmy
      derv2(lz,mx) = derv2(lz,mx) + d2ijs*djlx*dimz + d2iks*dklx*dimz + d2jks*dklx*djmz
      derv2(lx,my) = derv2(lx,my) + d2ijs*djly*dimx + d2iks*dkly*dimx + d2jks*dkly*djmx
      derv2(ly,my) = derv2(ly,my) + d2ijs*djly*dimy + d2iks*dkly*dimy + d2jks*dkly*djmy
      derv2(lz,my) = derv2(lz,my) + d2ijs*djly*dimz + d2iks*dkly*dimz + d2jks*dkly*djmz
      derv2(lx,mz) = derv2(lx,mz) + d2ijs*djlz*dimx + d2iks*dklz*dimx + d2jks*dklz*djmx
      derv2(ly,mz) = derv2(ly,mz) + d2ijs*djlz*dimy + d2iks*dklz*dimy + d2jks*dklz*djmy
      derv2(lz,mz) = derv2(lz,mz) + d2ijs*djlz*dimz + d2iks*dklz*dimz + d2jks*dklz*djmz
!
!  d2Edi2/d2Edj2/d2Edk2
!
      if (abs(d2i2s).gt.1.0d-8) then
        derv2(lx,mx) = derv2(lx,mx) + d2i2s*dilx*dimx
        derv2(ly,mx) = derv2(ly,mx) + d2i2s*dilx*dimy
        derv2(lz,mx) = derv2(lz,mx) + d2i2s*dilx*dimz
        derv2(lx,my) = derv2(lx,my) + d2i2s*dily*dimx
        derv2(ly,my) = derv2(ly,my) + d2i2s*dily*dimy
        derv2(lz,my) = derv2(lz,my) + d2i2s*dily*dimz
        derv2(lx,mz) = derv2(lx,mz) + d2i2s*dilz*dimx
        derv2(ly,mz) = derv2(ly,mz) + d2i2s*dilz*dimy
        derv2(lz,mz) = derv2(lz,mz) + d2i2s*dilz*dimz
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        derv2(lx,mx) = derv2(lx,mx) + d2j2s*djlx*djmx
        derv2(ly,mx) = derv2(ly,mx) + d2j2s*djlx*djmy
        derv2(lz,mx) = derv2(lz,mx) + d2j2s*djlx*djmz
        derv2(lx,my) = derv2(lx,my) + d2j2s*djly*djmx
        derv2(ly,my) = derv2(ly,my) + d2j2s*djly*djmy
        derv2(lz,my) = derv2(lz,my) + d2j2s*djly*djmz
        derv2(lx,mz) = derv2(lx,mz) + d2j2s*djlz*djmx
        derv2(ly,mz) = derv2(ly,mz) + d2j2s*djlz*djmy
        derv2(lz,mz) = derv2(lz,mz) + d2j2s*djlz*djmz
      endif
      if (abs(d2k2s).gt.1.0d-8.and.i.ne.k.and.j.ne.k) then
        derv2(lx,mx) = derv2(lx,mx) + d2k2s*dklx*dkmx
        derv2(ly,mx) = derv2(ly,mx) + d2k2s*dklx*dkmy
        derv2(lz,mx) = derv2(lz,mx) + d2k2s*dklx*dkmz
        derv2(lx,my) = derv2(lx,my) + d2k2s*dkly*dkmx
        derv2(ly,my) = derv2(ly,my) + d2k2s*dkly*dkmy
        derv2(lz,my) = derv2(lz,my) + d2k2s*dkly*dkmz
        derv2(lx,mz) = derv2(lx,mz) + d2k2s*dklz*dkmx
        derv2(ly,mz) = derv2(ly,mz) + d2k2s*dklz*dkmy
        derv2(lz,mz) = derv2(lz,mz) + d2k2s*dklz*dkmz
      endif
!
!  End of loop over l
!
    enddo
!
!  End of loop over m
!
  enddo
#ifdef TRACE
  call trace_out('d2charge3p')
#endif
!
  return
  end

  subroutine d2charge_finish
!
!  Calculates the second contribution to the second derivative matrices
!  due to charge derivatives from variable charge models using the pre-summed terms
!
!  NB: It is assumed that all atoms are included in the 2nd derivatives
!      as lfreeze is incompatible with variable charges
!
!   7/21 Created from d2charge
!   7/21 QEq self term corrections added
!   1/22 Symmetrisation of derv2 added
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
!  Copyright Curtin University 2022
!
!  Julian Gale, CIC, Curtin University, January 2022
!
  use control
  use current
  use derivatives
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: kx
  integer(i4)                                  :: ky
  integer(i4)                                  :: kz
  integer(i4)                                  :: l
  integer(i4)                                  :: lx
  integer(i4)                                  :: ly
  integer(i4)                                  :: lz
  integer(i4)                                  :: status
  real(dp)                                     :: d2i2
  real(dp)                                     :: d2ij
  real(dp)                                     :: dikx
  real(dp)                                     :: diky
  real(dp)                                     :: dikz
  real(dp)                                     :: dilx
  real(dp)                                     :: dily
  real(dp)                                     :: dilz
  real(dp),    dimension(:), allocatable, save :: dEdxyzdq
  real(dp)                                     :: dqikx
  real(dp)                                     :: dqiky
  real(dp)                                     :: dqikz
  real(dp)                                     :: dqilx
  real(dp)                                     :: dqily
  real(dp)                                     :: dqilz
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: time1
  real(dp)                                     :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charge_finish')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(dEdxyzdq(3*numat),stat=status)
  if (status/=0) call outofmemory('d2charge_finish','dEdxyzdq')
!***************
!  First pass  *
!***************
!
!  Loop over i-j
!
  do i = 1,numat
    dEdxyzdq(1:3*numat) = 0.0_dp
    do j = 1,i
      kx = - 2
      ky = - 1
      kz =   0
      d2ij = d2edqdq(j,i)
      do k = 1,numat
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
        dEdxyzdq(kx) = dEdxyzdq(kx) + d2ij*dqdxyz(kx,j)
        dEdxyzdq(ky) = dEdxyzdq(ky) + d2ij*dqdxyz(ky,j)
        dEdxyzdq(kz) = dEdxyzdq(kz) + d2ij*dqdxyz(kz,j)
      enddo
    enddo
!****************
!  Second pass  *
!****************
    kx = - 2
    ky = - 1
    kz =   0
    do k = 1,numat
      kx = kx + 3
      ky = ky + 3
      kz = kz + 3
!
      dikx = dEdxyzdq(kx)
      diky = dEdxyzdq(ky)
      dikz = dEdxyzdq(kz)
!
      dqikx = dqdxyz(kx,i)
      dqiky = dqdxyz(ky,i)
      dqikz = dqdxyz(kz,i)
!
      lx = - 2
      ly = - 1
      lz =   0
!
      do l = 1,k-1
        lx = lx + 3
        ly = ly + 3
        lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
        dilx = dEdxyzdq(lx)
        dily = dEdxyzdq(ly)
        dilz = dEdxyzdq(lz)
!
        dqilx = dqdxyz(lx,i)
        dqily = dqdxyz(ly,i)
        dqilz = dqdxyz(lz,i)
!
        if (l.lt.k) then
          derv2(lx,kx) = derv2(lx,kx) + dqilx*dikx + dqikx*dilx
          derv2(ly,kx) = derv2(ly,kx) + dqilx*diky + dqiky*dilx
          derv2(lz,kx) = derv2(lz,kx) + dqilx*dikz + dqikz*dilx
          derv2(lx,ky) = derv2(lx,ky) + dqily*dikx + dqikx*dily
          derv2(ly,ky) = derv2(ly,ky) + dqily*diky + dqiky*dily
          derv2(lz,ky) = derv2(lz,ky) + dqily*dikz + dqikz*dily
          derv2(lx,kz) = derv2(lx,kz) + dqilz*dikx + dqikx*dilz
          derv2(ly,kz) = derv2(ly,kz) + dqilz*diky + dqiky*dilz
          derv2(lz,kz) = derv2(lz,kz) + dqilz*dikz + dqikz*dilz
        else
          derv2(kx,lx) = derv2(kx,lx) + dqilx*dikx + dqikx*dilx
          derv2(ky,lx) = derv2(ky,lx) + dqily*dikx + dqikx*dily
          derv2(kz,lx) = derv2(kz,lx) + dqilz*dikx + dqikx*dilz
          derv2(kx,ly) = derv2(kx,ly) + dqilx*diky + dqiky*dilx
          derv2(ky,ly) = derv2(ky,ly) + dqily*diky + dqiky*dily
          derv2(kz,ly) = derv2(kz,ly) + dqilz*diky + dqiky*dilz
          derv2(kx,lz) = derv2(kx,lz) + dqilx*dikz + dqikz*dilx
          derv2(ky,lz) = derv2(ky,lz) + dqily*dikz + dqikz*dily
          derv2(kz,lz) = derv2(kz,lz) + dqilz*dikz + dqikz*dilz
        endif
!
!  End of loop over l
!
      enddo
!
!  End of loop over k
!
    enddo
!*************************
!  Extra term for QEq/H  *
!*************************
    if (abs(d2edq2(i)).gt.1.0d-8) then
      d2i2 = d2edq2(i)
!
      kx = - 2
      ky = - 1
      kz =   0
      do k = 1,numat
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
!
        dqikx = dqdxyz(kx,i)*d2i2
        dqiky = dqdxyz(ky,i)*d2i2
        dqikz = dqdxyz(kz,i)*d2i2
!
        lx = - 2
        ly = - 1
        lz =   0
!
        do l = 1,k-1
          lx = lx + 3
          ly = ly + 3
          lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta)
!
          dqilx = dqdxyz(lx,i)
          dqily = dqdxyz(ly,i)
          dqilz = dqdxyz(lz,i)
!
          derv2(lx,kx) = derv2(lx,kx) + dqilx*dqikx
          derv2(ly,kx) = derv2(ly,kx) + dqilx*dqiky
          derv2(lz,kx) = derv2(lz,kx) + dqilx*dqikz
          derv2(lx,ky) = derv2(lx,ky) + dqily*dqikx
          derv2(ly,ky) = derv2(ly,ky) + dqily*dqiky
          derv2(lz,ky) = derv2(lz,ky) + dqily*dqikz
          derv2(lx,kz) = derv2(lx,kz) + dqilz*dqikx
          derv2(ly,kz) = derv2(ly,kz) + dqilz*dqiky
          derv2(lz,kz) = derv2(lz,kz) + dqilz*dqikz
!
!  End of loop over l
!
        enddo
!
!  End of loop over k
!
      enddo
    endif
  enddo
!
!  Free local memory
!
  deallocate(dEdxyzdq,stat=status)
  if (status/=0) call deallocate_error('d2charge_finish','dEdxyzdq')
!
!  Symmetrise second derivative matrix after applying charge derivatives
!
  do i = 2,3*numat
    do j = 1,i-1
      derv2(i,j) = derv2(j,i)
    enddo
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charge_finish')
#endif
!
  return
  end

  subroutine d2charged_finish
!
!  Calculates the second contribution to the second derivative matrices
!  due to charge derivatives from variable charge models using the pre-summed terms
!  This version is for distributed memory parallel second derivatives.
!
!  NB: It is assumed that all atoms are included in the 2nd derivatives
!      as lfreeze is incompatible with variable charges
!
!   1/22 Created from d2charge_finish
!   1/22 dedqc and d2edqc added for variable charge second derivatives
!   6/22 Modified to avoid gfortran warning
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
!  Copyright Curtin University 2022
!
!  Julian Gale, CIC, Curtin University, June 2022
!
  use control
  use current
  use derivatives
  use parallel
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ind
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixf
  integer(i4)                                  :: iyf
  integer(i4)                                  :: izf
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: kx
  integer(i4)                                  :: ky
  integer(i4)                                  :: kz
  integer(i4)                                  :: kxf
  integer(i4)                                  :: kyf
  integer(i4)                                  :: kzf
  integer(i4)                                  :: l
  integer(i4)                                  :: lx
  integer(i4)                                  :: ly
  integer(i4)                                  :: lz
  integer(i4)                                  :: status
  real(dp)                                     :: d1ix
  real(dp)                                     :: d1iy
  real(dp)                                     :: d1iz
  real(dp)                                     :: d1kx
  real(dp)                                     :: d1ky
  real(dp)                                     :: d1kz
  real(dp)                                     :: d2i2
  real(dp)                                     :: d2ij
  real(dp)                                     :: dikx
  real(dp)                                     :: diky
  real(dp)                                     :: dikz
  real(dp)                                     :: dilx
  real(dp)                                     :: dily
  real(dp)                                     :: dilz
  real(dp)                                     :: djix
  real(dp)                                     :: djiy
  real(dp)                                     :: djiz
  real(dp)                                     :: dkix
  real(dp)                                     :: dkiy
  real(dp)                                     :: dkiz
  real(dp),    dimension(:), allocatable, save :: dEdxyzdq
  real(dp)                                     :: dqikx
  real(dp)                                     :: dqiky
  real(dp)                                     :: dqikz
  real(dp)                                     :: dqilx
  real(dp)                                     :: dqily
  real(dp)                                     :: dqilz
  real(dp),    dimension(:), allocatable, save :: dtmp
  real(dp),    dimension(:), allocatable, save :: dtmp2
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: time1
  real(dp)                                     :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charged_finish')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(dEdxyzdq(3*numat),stat=status)
  if (status/=0) call outofmemory('d2charged_finish','dEdxyzdq')
  allocate(dtmp(3*numat),stat=status)
  if (status/=0) call outofmemory('d2charged_finish','dtmp')
  allocate(dtmp2(3*numat),stat=status)
  if (status/=0) call outofmemory('d2charged_finish','dtmp2')
!
!  Sum d2edq2 over all nodes
!
  dtmp(1:numat) = d2edq2(1:numat)
  call sumall(dtmp,d2edq2,numat,"d2charged_finish","d2edq2") 
!
!  Sum dedqc over all nodes
!
  ind = 0
  do i = 1,numat
    do j = 1,3
      ind = ind + 1
      dtmp(ind) = dedqc(j,i)
    enddo
  enddo
  call sumall(dtmp,dedqc,3_i4*numat,"d2charged_finish","dedqc") 
!
!  Sum d2edqc over all nodes
!
  do i = 1,numat
    dtmp(1:3*numat) = d2edqc(1:3*numat,i)
    call sumall(dtmp,dtmp2,3_i4*numat,"d2charged_finish","d2edqc") 
    d2edqc(1:3*numat,i) = dtmp2(1:3*numat)
  enddo
  if (lstr) then
!
!  Sum ds2g / ds2gs / d2s2gs over all nodes
!
    do kl = 1,nstrains
      dtmp(1:numat) = ds2g(kl,1:numat)
      dtmp(numat+1:2*numat) = ds2gs(kl,1:numat)
      dtmp(2*numat+1:3*numat) = d2s2gs(kl,1:numat)
!
      call sumall(dtmp,dtmp2,3_i4*numat,"d2charged_finish","ds2g") 
!
      ds2g(kl,1:numat) = dtmp2(1:numat)
      ds2gs(kl,1:numat) = dtmp2(numat+1:2*numat)
      d2s2gs(kl,1:numat) = dtmp2(2*numat+1:3*numat)
    enddo
  endif
!***************
!  First pass  *
!***************
!
!  Loop over i-j
!
  do i = 1,numat
    dEdxyzdq(1:3*numat) = 0.0_dp
    iloc = atom2local(i)
    if (iloc.gt.0) then
      dtmp(1:numat) = d2edqdq(1:numat,iloc)
    endif
    call sendall(dtmp,numat,atom2node(i),"d2charged_finish","d2edqdq")
    call mpbarrier
    do j = 1,i
      kx = - 2
      ky = - 1
      kz =   0
      d2ij = dtmp(j)
      do k = 1,numat
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
        dEdxyzdq(kx) = dEdxyzdq(kx) + d2ij*dqdxyz(kx,j)
        dEdxyzdq(ky) = dEdxyzdq(ky) + d2ij*dqdxyz(ky,j)
        dEdxyzdq(kz) = dEdxyzdq(kz) + d2ij*dqdxyz(kz,j)
      enddo
    enddo
!****************
!  Second pass  *
!****************
    kx = - 2
    ky = - 1
    kz =   0
    do kk = 1,natomsonnode
      k = node2atom(kk)
      kx = kx + 3
      ky = ky + 3
      kz = kz + 3
!
      kxf = 3*(k-1) + 1
      kyf = kxf + 1
      kzf = kxf + 2
!
      dikx = dEdxyzdq(kxf)
      diky = dEdxyzdq(kyf)
      dikz = dEdxyzdq(kzf)
!
      dqikx = dqdxyz(kxf,i)
      dqiky = dqdxyz(kyf,i)
      dqikz = dqdxyz(kzf,i)
!
      lx = - 2
      ly = - 1
      lz =   0
!
      do l = 1,numat
        lx = lx + 3
        ly = ly + 3
        lz = lz + 3
!
        if (l.eq.k) cycle
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
        dilx = dEdxyzdq(lx)
        dily = dEdxyzdq(ly)
        dilz = dEdxyzdq(lz)
!
        dqilx = dqdxyz(lx,i)
        dqily = dqdxyz(ly,i)
        dqilz = dqdxyz(lz,i)
!
        derv2(lx,kx) = derv2(lx,kx) + dqilx*dikx + dqikx*dilx
        derv2(ly,kx) = derv2(ly,kx) + dqilx*diky + dqiky*dilx
        derv2(lz,kx) = derv2(lz,kx) + dqilx*dikz + dqikz*dilx
        derv2(lx,ky) = derv2(lx,ky) + dqily*dikx + dqikx*dily
        derv2(ly,ky) = derv2(ly,ky) + dqily*diky + dqiky*dily
        derv2(lz,ky) = derv2(lz,ky) + dqily*dikz + dqikz*dily
        derv2(lx,kz) = derv2(lx,kz) + dqilz*dikx + dqikx*dilz
        derv2(ly,kz) = derv2(ly,kz) + dqilz*diky + dqiky*dilz
        derv2(lz,kz) = derv2(lz,kz) + dqilz*dikz + dqikz*dilz
!
!  End of loop over l
!
      enddo
!
!  End of loop over k
!
    enddo
!*************************
!  Extra term for QEq/H  *
!*************************
    if (abs(d2edq2(i)).gt.1.0d-8) then
      d2i2 = d2edq2(i)
!
      kx = - 2
      ky = - 1
      kz =   0
      do kk = 1,natomsonnode
        k = node2atom(kk)
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
!
        dqikx = dqdxyz(kx,i)*d2i2
        dqiky = dqdxyz(ky,i)*d2i2
        dqikz = dqdxyz(kz,i)*d2i2
!
        lx = - 2
        ly = - 1
        lz =   0
!
        do l = 1,numat
          lx = lx + 3
          ly = ly + 3
          lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta)
!
          dqilx = dqdxyz(lx,i)
          dqily = dqdxyz(ly,i)
          dqilz = dqdxyz(lz,i)
!
          derv2(lx,kx) = derv2(lx,kx) + dqilx*dqikx
          derv2(ly,kx) = derv2(ly,kx) + dqilx*dqiky
          derv2(lz,kx) = derv2(lz,kx) + dqilx*dqikz
          derv2(lx,ky) = derv2(lx,ky) + dqily*dqikx
          derv2(ly,ky) = derv2(ly,ky) + dqily*dqiky
          derv2(lz,ky) = derv2(lz,ky) + dqily*dqikz
          derv2(lx,kz) = derv2(lx,kz) + dqilz*dqikx
          derv2(ly,kz) = derv2(ly,kz) + dqilz*dqiky
          derv2(lz,kz) = derv2(lz,kz) + dqilz*dqikz
!
!  End of loop over l
!
        enddo
!
!  End of loop over k
!
      enddo
    endif
  enddo
!****************************
!  Contribution from dedqc  *
!****************************
  ix = - 2
  iy = - 1
  iz =   0
  do ii = 1,natomsonnode
    i = node2atom(ii)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    ixf = 3*(i-1) + 1
    iyf = ixf + 1
    izf = ixf + 2
!
    kx = - 2
    ky = - 1
    kz =   0
    do k = 1,numat
      kx = kx + 3
      ky = ky + 3
      kz = kz + 3
!
      dikx = dqdxyz(kx,i)
      diky = dqdxyz(ky,i)
      dikz = dqdxyz(kz,i)
      dkix = dqdxyz(ixf,k)
      dkiy = dqdxyz(iyf,k)
      dkiz = dqdxyz(izf,k)
!
      d1ix = dedqc(1,i)
      d1iy = dedqc(2,i)
      d1iz = dedqc(3,i)
      d1kx = dedqc(1,k)
      d1ky = dedqc(2,k)
      d1kz = dedqc(3,k)
!
      if (i.ne.k) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-k
!
        if (k.gt.i) then
          derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx
          derv2(ky,ix) = derv2(ky,ix) - d1ix*diky
          derv2(kz,ix) = derv2(kz,ix) - d1ix*dikz
          derv2(kx,iy) = derv2(kx,iy) - d1iy*dikx
          derv2(ky,iy) = derv2(ky,iy) - d1iy*diky
          derv2(kz,iy) = derv2(kz,iy) - d1iy*dikz
          derv2(kx,iz) = derv2(kx,iz) - d1iz*dikx
          derv2(ky,iz) = derv2(ky,iz) - d1iz*diky
          derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz
!
          derv2(kx,ix) = derv2(kx,ix) - d1kx*dkix
          derv2(ky,ix) = derv2(ky,ix) - d1kx*dkiy
          derv2(kz,ix) = derv2(kz,ix) - d1kx*dkiz
          derv2(kx,iy) = derv2(kx,iy) - d1ky*dkix
          derv2(ky,iy) = derv2(ky,iy) - d1ky*dkiy
          derv2(kz,iy) = derv2(kz,iy) - d1ky*dkiz
          derv2(kx,iz) = derv2(kx,iz) - d1kz*dkix
          derv2(ky,iz) = derv2(ky,iz) - d1kz*dkiy
          derv2(kz,iz) = derv2(kz,iz) - d1kz*dkiz
        else
          derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx
          derv2(ky,ix) = derv2(ky,ix) - d1iy*dikx
          derv2(kz,ix) = derv2(kz,ix) - d1iz*dikx
          derv2(kx,iy) = derv2(kx,iy) - d1ix*diky
          derv2(ky,iy) = derv2(ky,iy) - d1iy*diky
          derv2(kz,iy) = derv2(kz,iy) - d1iz*diky
          derv2(kx,iz) = derv2(kx,iz) - d1ix*dikz
          derv2(ky,iz) = derv2(ky,iz) - d1iy*dikz
          derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz
!
          derv2(kx,ix) = derv2(kx,ix) - d1kx*dkix
          derv2(ky,ix) = derv2(ky,ix) - d1ky*dkix
          derv2(kz,ix) = derv2(kz,ix) - d1kz*dkix
          derv2(kx,iy) = derv2(kx,iy) - d1kx*dkiy
          derv2(ky,iy) = derv2(ky,iy) - d1ky*dkiy
          derv2(kz,iy) = derv2(kz,iy) - d1kz*dkiy
          derv2(kx,iz) = derv2(kx,iz) - d1kx*dkiz
          derv2(ky,iz) = derv2(ky,iz) - d1ky*dkiz
          derv2(kz,iz) = derv2(kz,iz) - d1kz*dkiz
        endif
      endif
    enddo
  enddo
!*****************************
!  Contribution from d2edqc  *
!*****************************
  ix = - 2
  iy = - 1
  iz =   0
  do ii = 1,natomsonnode
    i = node2atom(ii)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!
    ixf = 3*(i-1) + 1
    iyf = ixf + 1
    izf = ixf + 2
    do j = 1,numat
      djix = dqdxyz(ixf,j)
      djiy = dqdxyz(iyf,j)
      djiz = dqdxyz(izf,j)
      kx = -2
      ky = -1
      kz =  0
      do k = 1,numat
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
        if (k.gt.i) then
          derv2(kx,ix) = derv2(kx,ix) - d2edqc(kx,j)*djix - d2edqc(ixf,j)*dqdxyz(kx,j)
          derv2(kx,iy) = derv2(kx,iy) - d2edqc(ky,j)*djix - d2edqc(iyf,j)*dqdxyz(kx,j)
          derv2(kx,iz) = derv2(kx,iz) - d2edqc(kz,j)*djix - d2edqc(izf,j)*dqdxyz(kx,j)
          derv2(ky,ix) = derv2(ky,ix) - d2edqc(kx,j)*djiy - d2edqc(ixf,j)*dqdxyz(ky,j)
          derv2(ky,iy) = derv2(ky,iy) - d2edqc(ky,j)*djiy - d2edqc(iyf,j)*dqdxyz(ky,j)
          derv2(ky,iz) = derv2(ky,iz) - d2edqc(kz,j)*djiy - d2edqc(izf,j)*dqdxyz(ky,j)
          derv2(kz,ix) = derv2(kz,ix) - d2edqc(kx,j)*djiz - d2edqc(ixf,j)*dqdxyz(kz,j)
          derv2(kz,iy) = derv2(kz,iy) - d2edqc(ky,j)*djiz - d2edqc(iyf,j)*dqdxyz(kz,j)
          derv2(kz,iz) = derv2(kz,iz) - d2edqc(kz,j)*djiz - d2edqc(izf,j)*dqdxyz(kz,j)
        else
          derv2(kx,ix) = derv2(kx,ix) - d2edqc(kx,j)*djix - d2edqc(ixf,j)*dqdxyz(kx,j)
          derv2(kx,iy) = derv2(kx,iy) - d2edqc(kx,j)*djiy - d2edqc(ixf,j)*dqdxyz(ky,j)
          derv2(kx,iz) = derv2(kx,iz) - d2edqc(kx,j)*djiz - d2edqc(ixf,j)*dqdxyz(kz,j)
          derv2(ky,ix) = derv2(ky,ix) - d2edqc(ky,j)*djix - d2edqc(iyf,j)*dqdxyz(kx,j)
          derv2(ky,iy) = derv2(ky,iy) - d2edqc(ky,j)*djiy - d2edqc(iyf,j)*dqdxyz(ky,j)
          derv2(ky,iz) = derv2(ky,iz) - d2edqc(ky,j)*djiz - d2edqc(iyf,j)*dqdxyz(kz,j)
          derv2(kz,ix) = derv2(kz,ix) - d2edqc(kz,j)*djix - d2edqc(izf,j)*dqdxyz(kx,j)
          derv2(kz,iy) = derv2(kz,iy) - d2edqc(kz,j)*djiy - d2edqc(izf,j)*dqdxyz(ky,j)
          derv2(kz,iz) = derv2(kz,iz) - d2edqc(kz,j)*djiz - d2edqc(izf,j)*dqdxyz(kz,j)
        endif
      enddo
    enddo
  enddo
!*****************************
!  Contributions to strains  *
!*****************************
  if (lstr) then
    ix = - 2
    iy = - 1
    iz =   0
    do ii = 1,natomsonnode
      i = node2atom(ii)
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      ixf = 3*(i-1) + 1
      iyf = ixf + 1
      izf = ixf + 2
!
      do k = 1,numat
        dkix = dqdxyz(ixf,k)
        dkiy = dqdxyz(iyf,k)
        dkiz = dqdxyz(izf,k)
!
!  Mix strain-internal contribution
!
!  d2E/d(strain).dq x dq/d(alpha)
!
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) + ds2g(kl,k)*dkix
          derv3(iy,kl) = derv3(iy,kl) + ds2g(kl,k)*dkiy
          derv3(iz,kl) = derv3(iz,kl) + ds2g(kl,k)*dkiz
        enddo
!
!  d2E/dq.dq x dq/d(alpha) x dq/d(strain)
!
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) + ds2gs(kl,k)*dkix
          derv3(iy,kl) = derv3(iy,kl) + ds2gs(kl,k)*dkiy
          derv3(iz,kl) = derv3(iz,kl) + ds2gs(kl,k)*dkiz
        enddo
!
        do kl = 1,nstrains
          derv3(ix,kl) = derv3(ix,kl) + d2s2gs(kl,k)*dkix
          derv3(iy,kl) = derv3(iy,kl) + d2s2gs(kl,k)*dkiy
          derv3(iz,kl) = derv3(iz,kl) + d2s2gs(kl,k)*dkiz
        enddo
      enddo
    enddo
  endif
!
!  Free local memory
!
  deallocate(dtmp2,stat=status)
  if (status/=0) call deallocate_error('d2charged_finish','dtmp2')
  deallocate(dtmp,stat=status)
  if (status/=0) call deallocate_error('d2charged_finish','dtmp')
  deallocate(dEdxyzdq,stat=status)
  if (status/=0) call deallocate_error('d2charged_finish','dEdxyzdq')
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2charged_finish')
#endif
!
  return
  end

  subroutine d2chargepd_finish
!
!  Calculates the second contribution to the second derivative matrices
!  due to charge derivatives from variable charge models using the pre-summed terms
!  This version is for distributed memory parallel second derivatives.
!
!  NB: It is assumed that all atoms are included in the 2nd derivatives
!      as lfreeze is incompatible with variable charges
!
!   1/22 Created from d2charge_finish
!   1/22 dedqc and d2edqc added for variable charge second derivatives
!   6/22 Modified to avoid gfortran warning
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
!  Copyright Curtin University 2022
!
!  Julian Gale, CIC, Curtin University, June 2022
!
  use control
  use current
  use derivatives
  use parallel
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ind
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixf
  integer(i4)                                  :: iyf
  integer(i4)                                  :: izf
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kx
  integer(i4)                                  :: ky
  integer(i4)                                  :: kz
  integer(i4)                                  :: kxf
  integer(i4)                                  :: kyf
  integer(i4)                                  :: kzf
  integer(i4)                                  :: l
  integer(i4)                                  :: lx
  integer(i4)                                  :: ly
  integer(i4)                                  :: lz
  integer(i4)                                  :: status
  real(dp)                                     :: d1ix
  real(dp)                                     :: d1iy
  real(dp)                                     :: d1iz
  real(dp)                                     :: d1kx
  real(dp)                                     :: d1ky
  real(dp)                                     :: d1kz
  real(dp)                                     :: d2i2
  real(dp)                                     :: d2ij
  real(dp)                                     :: dikx
  real(dp)                                     :: diky
  real(dp)                                     :: dikz
  real(dp)                                     :: dilx
  real(dp)                                     :: dily
  real(dp)                                     :: dilz
  real(dp)                                     :: djix
  real(dp)                                     :: djiy
  real(dp)                                     :: djiz
  real(dp)                                     :: dkix
  real(dp)                                     :: dkiy
  real(dp)                                     :: dkiz
  real(dp),    dimension(:), allocatable, save :: dEdxyzdq
  real(dp)                                     :: dqikx
  real(dp)                                     :: dqiky
  real(dp)                                     :: dqikz
  real(dp)                                     :: dqilx
  real(dp)                                     :: dqily
  real(dp)                                     :: dqilz
  real(dp),    dimension(:), allocatable, save :: dtmp
  real(dp),    dimension(:), allocatable, save :: dtmp2
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: time1
  real(dp)                                     :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargepd_finish')
#endif
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(dEdxyzdq(3*numat),stat=status)
  if (status/=0) call outofmemory('d2chargepd_finish','dEdxyzdq')
  allocate(dtmp(3*numat),stat=status)
  if (status/=0) call outofmemory('d2chargepd_finish','dtmp')
  allocate(dtmp2(3*numat),stat=status)
  if (status/=0) call outofmemory('d2chargepd_finish','dtmp2')
!
!  Sum d2edq2 over all nodes
!
  dtmp(1:numat) = d2edq2(1:numat)
  call sumall(dtmp,d2edq2,numat,"d2chargepd_finish","d2edq2") 
!
!  Sum dedqc over all nodes
!
  ind = 0
  do i = 1,numat
    do j = 1,3
      ind = ind + 1
      dtmp(ind) = dedqc(j,i)
    enddo
  enddo
  call sumall(dtmp,dedqc,3_i4*numat,"d2chargepd_finish","dedqc") 
!
!  Sum d2edqc over all nodes
!
  do i = 1,numat
    dtmp(1:3*numat) = d2edqc(1:3*numat,i)
    call sumall(dtmp,dtmp2,3_i4*numat,"d2chargepd_finish","d2edqc") 
    d2edqc(1:3*numat,i) = dtmp2(1:3*numat)
  enddo
!***************
!  First pass  *
!***************
!
!  Loop over i-j
!
  do i = 1,numat
    dEdxyzdq(1:3*numat) = 0.0_dp
    iloc = atom2local(i)
    if (iloc.gt.0) then
      dtmp(1:numat) = d2edqdq(1:numat,iloc)
    endif
    call sendall(dtmp,numat,atom2node(i),"d2chargepd_finish","d2edqdq")
    call mpbarrier
    do j = 1,i
      kx = - 2
      ky = - 1
      kz =   0
      d2ij = dtmp(j)
      do k = 1,numat
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
        dEdxyzdq(kx) = dEdxyzdq(kx) + d2ij*dqdxyz(kx,j)
        dEdxyzdq(ky) = dEdxyzdq(ky) + d2ij*dqdxyz(ky,j)
        dEdxyzdq(kz) = dEdxyzdq(kz) + d2ij*dqdxyz(kz,j)
      enddo
    enddo
!****************
!  Second pass  *
!****************
    kx = - 2
    ky = - 1
    kz =   0
    do kk = 1,natomsonnode
      k = node2atom(kk)
      kx = kx + 3
      ky = ky + 3
      kz = kz + 3
!
      kxf = 3*(k-1) + 1
      kyf = kxf + 1
      kzf = kxf + 2
!
      dikx = dEdxyzdq(kxf)
      diky = dEdxyzdq(kyf)
      dikz = dEdxyzdq(kzf)
!
      dqikx = dqdxyz(kxf,i)
      dqiky = dqdxyz(kyf,i)
      dqikz = dqdxyz(kzf,i)
!
      lx = - 2
      ly = - 1
      lz =   0
!
      do l = 1,numat
        lx = lx + 3
        ly = ly + 3
        lz = lz + 3
!
        if (l.eq.k) cycle
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
        dilx = dEdxyzdq(lx)
        dily = dEdxyzdq(ly)
        dilz = dEdxyzdq(lz)
!
        dqilx = dqdxyz(lx,i)
        dqily = dqdxyz(ly,i)
        dqilz = dqdxyz(lz,i)
!
        derv2(lx,kx) = derv2(lx,kx) + dqilx*dikx + dqikx*dilx
        derv2(ly,kx) = derv2(ly,kx) + dqilx*diky + dqiky*dilx
        derv2(lz,kx) = derv2(lz,kx) + dqilx*dikz + dqikz*dilx
        derv2(lx,ky) = derv2(lx,ky) + dqily*dikx + dqikx*dily
        derv2(ly,ky) = derv2(ly,ky) + dqily*diky + dqiky*dily
        derv2(lz,ky) = derv2(lz,ky) + dqily*dikz + dqikz*dily
        derv2(lx,kz) = derv2(lx,kz) + dqilz*dikx + dqikx*dilz
        derv2(ly,kz) = derv2(ly,kz) + dqilz*diky + dqiky*dilz
        derv2(lz,kz) = derv2(lz,kz) + dqilz*dikz + dqikz*dilz
!
!  End of loop over l
!
      enddo
!
!  End of loop over k
!
    enddo
!*************************
!  Extra term for QEq/H  *
!*************************
    if (abs(d2edq2(i)).gt.1.0d-8) then
      d2i2 = d2edq2(i)
!
      kx = - 2
      ky = - 1
      kz =   0
      do kk = 1,natomsonnode
        k = node2atom(kk)
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
!
        dqikx = dqdxyz(kx,i)*d2i2
        dqiky = dqdxyz(ky,i)*d2i2
        dqikz = dqdxyz(kz,i)*d2i2
!
        lx = - 2
        ly = - 1
        lz =   0
!
        do l = 1,numat
          lx = lx + 3
          ly = ly + 3
          lz = lz + 3
!
          if (l.eq.k) cycle
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta)
!
          dqilx = dqdxyz(lx,i)
          dqily = dqdxyz(ly,i)
          dqilz = dqdxyz(lz,i)
!
          derv2(lx,kx) = derv2(lx,kx) + dqilx*dqikx
          derv2(ly,kx) = derv2(ly,kx) + dqilx*dqiky
          derv2(lz,kx) = derv2(lz,kx) + dqilx*dqikz
          derv2(lx,ky) = derv2(lx,ky) + dqily*dqikx
          derv2(ly,ky) = derv2(ly,ky) + dqily*dqiky
          derv2(lz,ky) = derv2(lz,ky) + dqily*dqikz
          derv2(lx,kz) = derv2(lx,kz) + dqilz*dqikx
          derv2(ly,kz) = derv2(ly,kz) + dqilz*dqiky
          derv2(lz,kz) = derv2(lz,kz) + dqilz*dqikz
!
!  End of loop over l
!
        enddo
!
!  End of loop over k
!
      enddo
    endif
  enddo
!****************************
!  Contribution from dedqc  *
!****************************
  ix = - 2
  iy = - 1
  iz =   0
  do ii = 1,natomsonnode
    i = node2atom(ii)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
    ixf = 3*(i-1) + 1
    iyf = ixf + 1
    izf = ixf + 2
!
    kx = - 2
    ky = - 1
    kz =   0
    do k = 1,numat
      kx = kx + 3
      ky = ky + 3
      kz = kz + 3
!
      dikx = dqdxyz(kx,i)
      diky = dqdxyz(ky,i)
      dikz = dqdxyz(kz,i)
      dkix = dqdxyz(ixf,k)
      dkiy = dqdxyz(iyf,k)
      dkiz = dqdxyz(izf,k)
!
      d1ix = dedqc(1,i)
      d1iy = dedqc(2,i)
      d1iz = dedqc(3,i)
      d1kx = dedqc(1,k)
      d1ky = dedqc(2,k)
      d1kz = dedqc(3,k)
!
      if (i.ne.k) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-k
!
        if (k.gt.i) then
          derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx
          derv2(ky,ix) = derv2(ky,ix) - d1ix*diky
          derv2(kz,ix) = derv2(kz,ix) - d1ix*dikz
          derv2(kx,iy) = derv2(kx,iy) - d1iy*dikx
          derv2(ky,iy) = derv2(ky,iy) - d1iy*diky
          derv2(kz,iy) = derv2(kz,iy) - d1iy*dikz
          derv2(kx,iz) = derv2(kx,iz) - d1iz*dikx
          derv2(ky,iz) = derv2(ky,iz) - d1iz*diky
          derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz
!
          derv2(kx,ix) = derv2(kx,ix) - d1kx*dkix
          derv2(ky,ix) = derv2(ky,ix) - d1kx*dkiy
          derv2(kz,ix) = derv2(kz,ix) - d1kx*dkiz
          derv2(kx,iy) = derv2(kx,iy) - d1ky*dkix
          derv2(ky,iy) = derv2(ky,iy) - d1ky*dkiy
          derv2(kz,iy) = derv2(kz,iy) - d1ky*dkiz
          derv2(kx,iz) = derv2(kx,iz) - d1kz*dkix
          derv2(ky,iz) = derv2(ky,iz) - d1kz*dkiy
          derv2(kz,iz) = derv2(kz,iz) - d1kz*dkiz
        else
          derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx
          derv2(ky,ix) = derv2(ky,ix) - d1iy*dikx
          derv2(kz,ix) = derv2(kz,ix) - d1iz*dikx
          derv2(kx,iy) = derv2(kx,iy) - d1ix*diky
          derv2(ky,iy) = derv2(ky,iy) - d1iy*diky
          derv2(kz,iy) = derv2(kz,iy) - d1iz*diky
          derv2(kx,iz) = derv2(kx,iz) - d1ix*dikz
          derv2(ky,iz) = derv2(ky,iz) - d1iy*dikz
          derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz
!
          derv2(kx,ix) = derv2(kx,ix) - d1kx*dkix
          derv2(ky,ix) = derv2(ky,ix) - d1ky*dkix
          derv2(kz,ix) = derv2(kz,ix) - d1kz*dkix
          derv2(kx,iy) = derv2(kx,iy) - d1kx*dkiy
          derv2(ky,iy) = derv2(ky,iy) - d1ky*dkiy
          derv2(kz,iy) = derv2(kz,iy) - d1kz*dkiy
          derv2(kx,iz) = derv2(kx,iz) - d1kx*dkiz
          derv2(ky,iz) = derv2(ky,iz) - d1ky*dkiz
          derv2(kz,iz) = derv2(kz,iz) - d1kz*dkiz
        endif
      endif
    enddo
  enddo
!*****************************
!  Contribution from d2edqc  *
!*****************************
  ix = - 2
  iy = - 1
  iz =   0
  do ii = 1,natomsonnode
    i = node2atom(ii)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!
    ixf = 3*(i-1) + 1
    iyf = ixf + 1
    izf = ixf + 2
    do j = 1,numat
      djix = dqdxyz(ixf,j)
      djiy = dqdxyz(iyf,j)
      djiz = dqdxyz(izf,j)
      kx = -2
      ky = -1
      kz =  0
      do k = 1,numat
        kx = kx + 3
        ky = ky + 3
        kz = kz + 3
        if (k.gt.i) then
          derv2(kx,ix) = derv2(kx,ix) - d2edqc(kx,j)*djix - d2edqc(ixf,j)*dqdxyz(kx,j)
          derv2(kx,iy) = derv2(kx,iy) - d2edqc(ky,j)*djix - d2edqc(iyf,j)*dqdxyz(kx,j)
          derv2(kx,iz) = derv2(kx,iz) - d2edqc(kz,j)*djix - d2edqc(izf,j)*dqdxyz(kx,j)
          derv2(ky,ix) = derv2(ky,ix) - d2edqc(kx,j)*djiy - d2edqc(ixf,j)*dqdxyz(ky,j)
          derv2(ky,iy) = derv2(ky,iy) - d2edqc(ky,j)*djiy - d2edqc(iyf,j)*dqdxyz(ky,j)
          derv2(ky,iz) = derv2(ky,iz) - d2edqc(kz,j)*djiy - d2edqc(izf,j)*dqdxyz(ky,j)
          derv2(kz,ix) = derv2(kz,ix) - d2edqc(kx,j)*djiz - d2edqc(ixf,j)*dqdxyz(kz,j)
          derv2(kz,iy) = derv2(kz,iy) - d2edqc(ky,j)*djiz - d2edqc(iyf,j)*dqdxyz(kz,j)
          derv2(kz,iz) = derv2(kz,iz) - d2edqc(kz,j)*djiz - d2edqc(izf,j)*dqdxyz(kz,j)
        elseif (k.lt.i) then
          derv2(kx,ix) = derv2(kx,ix) - d2edqc(kx,j)*djix - d2edqc(ixf,j)*dqdxyz(kx,j)
          derv2(kx,iy) = derv2(kx,iy) - d2edqc(kx,j)*djiy - d2edqc(ixf,j)*dqdxyz(ky,j)
          derv2(kx,iz) = derv2(kx,iz) - d2edqc(kx,j)*djiz - d2edqc(ixf,j)*dqdxyz(kz,j)
          derv2(ky,ix) = derv2(ky,ix) - d2edqc(ky,j)*djix - d2edqc(iyf,j)*dqdxyz(kx,j)
          derv2(ky,iy) = derv2(ky,iy) - d2edqc(ky,j)*djiy - d2edqc(iyf,j)*dqdxyz(ky,j)
          derv2(ky,iz) = derv2(ky,iz) - d2edqc(ky,j)*djiz - d2edqc(iyf,j)*dqdxyz(kz,j)
          derv2(kz,ix) = derv2(kz,ix) - d2edqc(kz,j)*djix - d2edqc(izf,j)*dqdxyz(kx,j)
          derv2(kz,iy) = derv2(kz,iy) - d2edqc(kz,j)*djiy - d2edqc(izf,j)*dqdxyz(ky,j)
          derv2(kz,iz) = derv2(kz,iz) - d2edqc(kz,j)*djiz - d2edqc(izf,j)*dqdxyz(kz,j)
        endif
      enddo
    enddo
  enddo
!
!  Free local memory
!
  deallocate(dtmp2,stat=status)
  if (status/=0) call deallocate_error('d2chargepd_finish','dtmp2')
  deallocate(dtmp,stat=status)
  if (status/=0) call deallocate_error('d2chargepd_finish','dtmp')
  deallocate(dEdxyzdq,stat=status)
  if (status/=0) call deallocate_error('d2chargepd_finish','dEdxyzdq')
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargepd_finish')
#endif
!
  return
  end

  subroutine d2chargep2D(i,j,nor,xtmp,ytmp,ix,iy,iz,jx,jy,jz,dei,dej,d1ir,d1zi,d1jr,d1zj, &
                         d2i2r,d2ijr,d2j2r,d2self,d2trm1dij,d2trm1diz,d2trm1djz,dtrm1di, &
                         dtrm1dj)
!
!  Calculates the contribution to the phonon matrices
!  due to charge derivatives from EEM/QEq. 2-D version.
!
!   3/01 Created from d2chargep
!   3/03 Frozen charge atoms accommodated
!   9/04 Modified to generalise to charge derivatives other than from EEM
!  11/07 Unused variables removed
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   2/21 GFNFF modifications added
!   6/21 Exact d2q terms rearranged to speed up
!   7/21 Timing added
!   7/21 New algorithm to avoid 4th power scaling added
!  10/21 lgfnff moved to control
!   9/23 Trapping of shells added
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
  use control
  use current
  use derivatives
  use element
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  integer(i4), intent(in)    :: i
  integer(i4), intent(in)    :: ix
  integer(i4), intent(in)    :: iy
  integer(i4), intent(in)    :: iz
  integer(i4), intent(in)    :: j
  integer(i4), intent(in)    :: jx
  integer(i4), intent(in)    :: jy
  integer(i4), intent(in)    :: jz
  integer(i4), intent(in)    :: nor
  real(dp),    intent(in)    :: d1ir(*)
  real(dp),    intent(in)    :: d1jr(*)
  real(dp),    intent(in)    :: d1zi(*)
  real(dp),    intent(in)    :: d1zj(*)
  real(dp),    intent(in)    :: dtrm1di
  real(dp),    intent(in)    :: dtrm1dj
  real(dp),    intent(in)    :: d2i2r(*)
  real(dp),    intent(in)    :: d2ijr(*)
  real(dp),    intent(in)    :: d2j2r(*)
  real(dp),    intent(in)    :: d2self
  real(dp),    intent(in)    :: d2trm1dij
  real(dp),    intent(in)    :: d2trm1diz
  real(dp),    intent(in)    :: d2trm1djz
  real(dp),    intent(in)    :: dei(*)
  real(dp),    intent(in)    :: dej(*)
  real(dp),    intent(in)    :: xtmp(*)
  real(dp),    intent(in)    :: ytmp(*)
!
!  Local variables
!
  integer(i4)                :: iv
  integer(i4)                :: k
  integer(i4)                :: kx
  integer(i4)                :: ky
  integer(i4)                :: kz
  integer(i4)                :: l
  integer(i4)                :: lx
  integer(i4)                :: ly
  integer(i4)                :: lz
  integer(i4)                :: m
  integer(i4)                :: mnxx
  integer(i4)                :: mnxy
  integer(i4)                :: mnxz
  integer(i4)                :: mnyx
  integer(i4)                :: mnyy
  integer(i4)                :: mnyz
  integer(i4)                :: mnzx
  integer(i4)                :: mnzy
  integer(i4)                :: mnzz
  integer(i4)                :: mx
  integer(i4)                :: my
  integer(i4)                :: mz
  integer(i4)                :: n
  integer(i4)                :: nx
  logical                    :: lNonIJQDeriv
  real(dp)                   :: d1ix
  real(dp)                   :: d1iy
  real(dp)                   :: d1iz
  real(dp)                   :: d1jx
  real(dp)                   :: d1jy
  real(dp)                   :: d1jz
  real(dp)                   :: d2ijs
  real(dp)                   :: d2i2s
  real(dp)                   :: d2j2s
  real(dp)                   :: deisum
  real(dp)                   :: dejsum
  real(dp)                   :: dikx
  real(dp)                   :: diky
  real(dp)                   :: dikz
  real(dp)                   :: d2ikx
  real(dp)                   :: d2iky
  real(dp)                   :: d2ikz
  real(dp)                   :: dilx
  real(dp)                   :: dily
  real(dp)                   :: dilz
  real(dp)                   :: djkx
  real(dp)                   :: djky
  real(dp)                   :: djkz
  real(dp)                   :: d2jkx
  real(dp)                   :: d2jky
  real(dp)                   :: d2jkz
  real(dp)                   :: djlx
  real(dp)                   :: djly
  real(dp)                   :: djlz
  real(dp)                   :: g_cpu_time
  real(dp)                   :: time1
  real(dp)                   :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargep2D')
#endif
!
  time1 = g_cpu_time()
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d1ix = 0.0_dp
  d1iy = 0.0_dp
  d1iz = 0.0_dp
  d1jx = 0.0_dp
  d1jy = 0.0_dp
  d1jz = 0.0_dp
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d1ix = d1ix + d1ir(iv)*xtmp(iv)
    d1iy = d1iy + d1ir(iv)*ytmp(iv)
    d1iz = d1iz + d1zi(iv)
    d1jx = d1jx + d1jr(iv)*xtmp(iv)
    d1jy = d1jy + d1jr(iv)*ytmp(iv)
    d1jz = d1jz + d1zj(iv)
    d2ijs = d2ijs + d2ijr(iv)
    d2i2s = d2i2s + d2i2r(iv)
    d2j2s = d2j2s + d2j2r(iv)
  enddo
  d1iz = d1iz + d2trm1diz
  d1jz = d1jz + d2trm1djz
  d2ijs = d2ijs + d2trm1dij
  if (leem) then
    if (nat(i).gt.maxele) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    endif
  endif
!
  if (.not.loldd2q) then
!******************
!  New algorithm  *
!******************
    d2edqdq(j,i) = d2edqdq(j,i) + d2ijs
    d2edq2(i) = d2edq2(i) + d2i2s
    if (i.ne.j) then
      d2edqdq(i,j) = d2edqdq(i,j) + d2ijs
      d2edq2(j) = d2edq2(j) + d2j2s
    endif
  endif
!
!  The following flag causes the second derivatives of the charge distribution 
!  with respect to two atoms neither of which are the i/j of the energy term
!  to be calculated. At the moment this contribution is not used in any method
!  and so the flag is set to false for efficiency.
!
  lNonIJQDeriv = .false.
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff.and..not.lreaxFF) then
    deisum = dtrm1di
    dejsum = dtrm1dj
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(i)
      k = nqatomptr(m,i)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (i.gt.k) then
        derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,i)
        derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,i)
        derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,i)
        derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,i)
        derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,i)
        derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,i)
        derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,i)
        derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,i)
        derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,i)
      elseif (i.lt.k) then
        derv2(ix,kx) = derv2(ix,kx) + deisum*d2qdxyz2(mnxx,i)
        derv2(iy,kx) = derv2(iy,kx) + deisum*d2qdxyz2(mnxy,i)
        derv2(iz,kx) = derv2(iz,kx) + deisum*d2qdxyz2(mnxz,i)
        derv2(ix,ky) = derv2(ix,ky) + deisum*d2qdxyz2(mnxy,i)
        derv2(iy,ky) = derv2(iy,ky) + deisum*d2qdxyz2(mnyy,i)
        derv2(iz,ky) = derv2(iz,ky) + deisum*d2qdxyz2(mnyz,i)
        derv2(ix,kz) = derv2(ix,kz) + deisum*d2qdxyz2(mnxz,i)
        derv2(iy,kz) = derv2(iy,kz) + deisum*d2qdxyz2(mnyz,i)
        derv2(iz,kz) = derv2(iz,kz) + deisum*d2qdxyz2(mnzz,i)
      endif
    enddo
!
    do m = 1,nqatoms(j)
      k = nqatomptr(m,j)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (j.gt.k) then
        derv2(kx,jx) = derv2(kx,jx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(ky,jx) = derv2(ky,jx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(kz,jx) = derv2(kz,jx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(kx,jy) = derv2(kx,jy) + dejsum*d2qdxyz2(mnxy,j)
        derv2(ky,jy) = derv2(ky,jy) + dejsum*d2qdxyz2(mnyy,j)
        derv2(kz,jy) = derv2(kz,jy) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kx,jz) = derv2(kx,jz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(ky,jz) = derv2(ky,jz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kz,jz) = derv2(kz,jz) + dejsum*d2qdxyz2(mnzz,j)
      elseif (j.lt.k) then
        derv2(jx,kx) = derv2(jx,kx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(jy,kx) = derv2(jy,kx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jz,kx) = derv2(jz,kx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jx,ky) = derv2(jx,ky) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jy,ky) = derv2(jy,ky) + dejsum*d2qdxyz2(mnyy,j)
        derv2(jz,ky) = derv2(jz,ky) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jx,kz) = derv2(jx,kz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jy,kz) = derv2(jy,kz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jz,kz) = derv2(jz,kz) + dejsum*d2qdxyz2(mnzz,j)
      endif
    enddo
!
!  Terms involving Q derivatives for non-i/j case
!
    if (lNonIJQDeriv) then
      do m = 2,nqatoms(i)
        k = nqatomptr(m,i)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,i)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,i)
          derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,i)
          derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,i)
          derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,i)
          derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,i)
          derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,i)
          derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,i)
          derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,i)
          derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,i)
        enddo
      enddo
!
      do m = 2,nqatoms(j)
        k = nqatomptr(m,j)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,j)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,j)
          derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,j)
          derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,j)
          derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,j)
          derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,j)
          derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,j)
          derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,j)
          derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,j)
          derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,j)
        enddo
      enddo
    endif
  endif
!
!  Loop over atoms to add ij-k contributions
!
  kx = - 2
  ky = - 1
  kz =   0
  do k = 1,numat
    kx = kx + 3
    ky = ky + 3
    kz = kz + 3
    dikx = dqdxyz(kx,i)
    diky = dqdxyz(ky,i)
    dikz = dqdxyz(kz,i)
    djkx = dqdxyz(kx,j)
    djky = dqdxyz(ky,j)
    djkz = dqdxyz(kz,j)
    if (i.ne.k) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-k
!
      if (k.le.i) then
        derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx
        derv2(ky,ix) = derv2(ky,ix) - d1iy*dikx
        derv2(kz,ix) = derv2(kz,ix) - d1iz*dikx
        derv2(kx,iy) = derv2(kx,iy) - d1ix*diky
        derv2(ky,iy) = derv2(ky,iy) - d1iy*diky
        derv2(kz,iy) = derv2(kz,iy) - d1iz*diky
        derv2(kx,iz) = derv2(kx,iz) - d1ix*dikz
        derv2(ky,iz) = derv2(ky,iz) - d1iy*dikz
        derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz
!
        derv2(kx,ix) = derv2(kx,ix) - d1jx*djkx
        derv2(ky,ix) = derv2(ky,ix) - d1jy*djkx
        derv2(kz,ix) = derv2(kz,ix) - d1jz*djkx
        derv2(kx,iy) = derv2(kx,iy) - d1jx*djky
        derv2(ky,iy) = derv2(ky,iy) - d1jy*djky
        derv2(kz,iy) = derv2(kz,iy) - d1jz*djky
        derv2(kx,iz) = derv2(kx,iz) - d1jx*djkz
        derv2(ky,iz) = derv2(ky,iz) - d1jy*djkz
        derv2(kz,iz) = derv2(kz,iz) - d1jz*djkz
      else
        derv2(ix,kx) = derv2(ix,kx) - d1ix*dikx
        derv2(iy,kx) = derv2(iy,kx) - d1iy*dikx
        derv2(iz,kx) = derv2(iz,kx) - d1iz*dikx
        derv2(ix,ky) = derv2(ix,ky) - d1ix*diky
        derv2(iy,ky) = derv2(iy,ky) - d1iy*diky
        derv2(iz,ky) = derv2(iz,ky) - d1iz*diky
        derv2(ix,kz) = derv2(ix,kz) - d1ix*dikz
        derv2(iy,kz) = derv2(iy,kz) - d1iy*dikz
        derv2(iz,kz) = derv2(iz,kz) - d1iz*dikz
!
        derv2(ix,kx) = derv2(ix,kx) - d1jx*djkx
        derv2(iy,kx) = derv2(iy,kx) - d1jy*djkx
        derv2(iz,kx) = derv2(iz,kx) - d1jz*djkx
        derv2(ix,ky) = derv2(ix,ky) - d1jx*djky
        derv2(iy,ky) = derv2(iy,ky) - d1jy*djky
        derv2(iz,ky) = derv2(iz,ky) - d1jz*djky
        derv2(ix,kz) = derv2(ix,kz) - d1jx*djkz
        derv2(iy,kz) = derv2(iy,kz) - d1jy*djkz
        derv2(iz,kz) = derv2(iz,kz) - d1jz*djkz
      endif
    endif
    if (j.ne.k) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : j-k
!
      if (k.le.j) then
        derv2(kx,jx) = derv2(kx,jx) + d1jx*djkx
        derv2(ky,jx) = derv2(ky,jx) + d1jy*djkx
        derv2(kz,jx) = derv2(kz,jx) + d1jz*djkx
        derv2(kx,jy) = derv2(kx,jy) + d1jx*djky
        derv2(ky,jy) = derv2(ky,jy) + d1jy*djky
        derv2(kz,jy) = derv2(kz,jy) + d1jz*djky
        derv2(kx,jz) = derv2(kx,jz) + d1jx*djkz
        derv2(ky,jz) = derv2(ky,jz) + d1jy*djkz
        derv2(kz,jz) = derv2(kz,jz) + d1jz*djkz
!
        derv2(kx,jx) = derv2(kx,jx) + d1ix*dikx
        derv2(ky,jx) = derv2(ky,jx) + d1iy*dikx
        derv2(kz,jx) = derv2(kz,jx) + d1iz*dikx
        derv2(kx,jy) = derv2(kx,jy) + d1ix*diky
        derv2(ky,jy) = derv2(ky,jy) + d1iy*diky
        derv2(kz,jy) = derv2(kz,jy) + d1iz*diky
        derv2(kx,jz) = derv2(kx,jz) + d1ix*dikz
        derv2(ky,jz) = derv2(ky,jz) + d1iy*dikz
        derv2(kz,jz) = derv2(kz,jz) + d1iz*dikz
      else
        derv2(jx,kx) = derv2(jx,kx) + d1jx*djkx
        derv2(jy,kx) = derv2(jy,kx) + d1jy*djkx
        derv2(jz,kx) = derv2(jz,kx) + d1jz*djkx
        derv2(jx,ky) = derv2(jx,ky) + d1jx*djky
        derv2(jy,ky) = derv2(jy,ky) + d1jy*djky
        derv2(jz,ky) = derv2(jz,ky) + d1jz*djky
        derv2(jx,kz) = derv2(jx,kz) + d1jx*djkz
        derv2(jy,kz) = derv2(jy,kz) + d1jy*djkz
        derv2(jz,kz) = derv2(jz,kz) + d1jz*djkz
!
        derv2(jx,kx) = derv2(jx,kx) + d1ix*dikx
        derv2(jy,kx) = derv2(jy,kx) + d1iy*dikx
        derv2(jz,kx) = derv2(jz,kx) + d1iz*dikx
        derv2(jx,ky) = derv2(jx,ky) + d1ix*diky
        derv2(jy,ky) = derv2(jy,ky) + d1iy*diky
        derv2(jz,ky) = derv2(jz,ky) + d1iz*diky
        derv2(jx,kz) = derv2(jx,kz) + d1ix*dikz
        derv2(jy,kz) = derv2(jy,kz) + d1iy*dikz
        derv2(jz,kz) = derv2(jz,kz) + d1iz*dikz
      endif
    endif
    if (loldd2q) then
!******************
!  Old algorithm  *
!******************
!
!  Loop over atoms to add ij-kl correction
!
      lx = - 2
      ly = - 1
      lz =   0
!
!  Scale charge derivatives for k by d2ijs
!
      d2ikx = d2ijs*dikx
      d2iky = d2ijs*diky
      d2ikz = d2ijs*dikz
      d2jkx = d2ijs*djkx
      d2jky = d2ijs*djky
      d2jkz = d2ijs*djkz
!
      do l = 1,k-1
        lx = lx + 3
        ly = ly + 3
        lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
        dilx = dqdxyz(lx,i)
        dily = dqdxyz(ly,i)
        dilz = dqdxyz(lz,i)
        djlx = dqdxyz(lx,j)
        djly = dqdxyz(ly,j)
        djlz = dqdxyz(lz,j)
!
        derv2(lx,kx) = derv2(lx,kx) + dilx*d2jkx + djlx*d2ikx
        derv2(ly,kx) = derv2(ly,kx) + dilx*d2jky + djlx*d2iky
        derv2(lz,kx) = derv2(lz,kx) + dilx*d2jkz + djlx*d2ikz
        derv2(lx,ky) = derv2(lx,ky) + dily*d2jkx + djly*d2ikx
        derv2(ly,ky) = derv2(ly,ky) + dily*d2jky + djly*d2iky
        derv2(lz,ky) = derv2(lz,ky) + dily*d2jkz + djly*d2ikz
        derv2(lx,kz) = derv2(lx,kz) + dilz*d2jkx + djlz*d2ikx
        derv2(ly,kz) = derv2(ly,kz) + dilz*d2jky + djlz*d2iky
        derv2(lz,kz) = derv2(lz,kz) + dilz*d2jkz + djlz*d2ikz
!
!  d2Edi2/d2Edj2 - only non-zero for QEq/H
!
        if (abs(d2i2s).gt.1.0d-8) then
          derv2(lx,kx) = derv2(lx,kx) + d2i2s*dilx*dikx
          derv2(ly,kx) = derv2(ly,kx) + d2i2s*dilx*diky
          derv2(lz,kx) = derv2(lz,kx) + d2i2s*dilx*dikz
          derv2(lx,ky) = derv2(lx,ky) + d2i2s*dily*dikx
          derv2(ly,ky) = derv2(ly,ky) + d2i2s*dily*diky
          derv2(lz,ky) = derv2(lz,ky) + d2i2s*dily*dikz
          derv2(lx,kz) = derv2(lx,kz) + d2i2s*dilz*dikx
          derv2(ly,kz) = derv2(ly,kz) + d2i2s*dilz*diky
          derv2(lz,kz) = derv2(lz,kz) + d2i2s*dilz*dikz
        endif
        if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
          derv2(lx,kx) = derv2(lx,kx) + d2j2s*djlx*djkx
          derv2(ly,kx) = derv2(ly,kx) + d2j2s*djlx*djky
          derv2(lz,kx) = derv2(lz,kx) + d2j2s*djlx*djkz
          derv2(lx,ky) = derv2(lx,ky) + d2j2s*djly*djkx
          derv2(ly,ky) = derv2(ly,ky) + d2j2s*djly*djky
          derv2(lz,ky) = derv2(lz,ky) + d2j2s*djly*djkz
          derv2(lx,kz) = derv2(lx,kz) + d2j2s*djlz*djkx
          derv2(ly,kz) = derv2(ly,kz) + d2j2s*djlz*djky
          derv2(lz,kz) = derv2(lz,kz) + d2j2s*djlz*djkz
        endif
!
!  End of loop over l
!
      enddo
    endif
!
!  End of loop over k
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargep2D')
#endif
!
  return
  end

  subroutine d2chargep2Dd(iloc,i,j,nor,xtmp,ytmp,ix,iy,iz,jx,jy,jz,dei,dej,d1ir,d1zi,d1jr,d1zj, &
                          d2i2r,d2ijr,d2j2r,d2self,d2trm1dij,d2trm1diz,d2trm1djz,dtrm1di, &
                          dtrm1dj)
!
!  Calculates the contribution to the phonon matrices
!  due to charge derivatives from EEM/QEq. 2-D version.
!  This version is for calls from distributed memory routines.
!
!   1/22 Created from d2chargep2D and d2charged
!   9/23 Trapping of shells added
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
  use control
  use current
  use derivatives
  use element
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed arguments
!
  integer(i4), intent(in)    :: i       ! Atom i
  integer(i4), intent(in)    :: iloc    ! Atom i number on local node
  integer(i4), intent(in)    :: ix
  integer(i4), intent(in)    :: iy
  integer(i4), intent(in)    :: iz
  integer(i4), intent(in)    :: j
  integer(i4), intent(in)    :: jx
  integer(i4), intent(in)    :: jy
  integer(i4), intent(in)    :: jz
  integer(i4), intent(in)    :: nor
  real(dp),    intent(in)    :: d1ir(*)
  real(dp),    intent(in)    :: d1jr(*)
  real(dp),    intent(in)    :: d1zi(*)
  real(dp),    intent(in)    :: d1zj(*)
  real(dp),    intent(in)    :: dtrm1di
  real(dp),    intent(in)    :: dtrm1dj
  real(dp),    intent(in)    :: d2i2r(*)
  real(dp),    intent(in)    :: d2ijr(*)
  real(dp),    intent(in)    :: d2j2r(*)
  real(dp),    intent(in)    :: d2self
  real(dp),    intent(in)    :: d2trm1dij
  real(dp),    intent(in)    :: d2trm1diz
  real(dp),    intent(in)    :: d2trm1djz
  real(dp),    intent(in)    :: dei(*)
  real(dp),    intent(in)    :: dej(*)
  real(dp),    intent(in)    :: xtmp(*)
  real(dp),    intent(in)    :: ytmp(*)
!
!  Local variables
!
  integer(i4)                :: iv
  integer(i4)                :: ixf
  integer(i4)                :: iyf
  integer(i4)                :: izf
  integer(i4)                :: k
  integer(i4)                :: kx
  integer(i4)                :: ky
  integer(i4)                :: kz
  integer(i4)                :: m
  integer(i4)                :: mnxx
  integer(i4)                :: mnxy
  integer(i4)                :: mnxz
  integer(i4)                :: mnyy
  integer(i4)                :: mnyz
  integer(i4)                :: mnzz
  integer(i4)                :: mx
  integer(i4)                :: my
  integer(i4)                :: mz
  real(dp)                   :: d1ix
  real(dp)                   :: d1iy
  real(dp)                   :: d1iz
  real(dp)                   :: d1jx
  real(dp)                   :: d1jy
  real(dp)                   :: d1jz
  real(dp)                   :: d2ijs
  real(dp)                   :: d2i2s
  real(dp)                   :: d2j2s
  real(dp)                   :: deisum
  real(dp)                   :: dejsum
  real(dp)                   :: g_cpu_time
  real(dp)                   :: time1
  real(dp)                   :: time2
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargep2Dd')
#endif
!
  time1 = g_cpu_time()
!
  ixf = 3*(i-1) + 1
  iyf = ixf + 1
  izf = ixf + 2
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d1ix = 0.0_dp
  d1iy = 0.0_dp
  d1iz = 0.0_dp
  d1jx = 0.0_dp
  d1jy = 0.0_dp
  d1jz = 0.0_dp
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d1ix = d1ix + d1ir(iv)*xtmp(iv)
    d1iy = d1iy + d1ir(iv)*ytmp(iv)
    d1iz = d1iz + d1zi(iv)
    d1jx = d1jx + d1jr(iv)*xtmp(iv)
    d1jy = d1jy + d1jr(iv)*ytmp(iv)
    d1jz = d1jz + d1zj(iv)
    d2ijs = d2ijs + d2ijr(iv)
    d2i2s = d2i2s + d2i2r(iv)
    d2j2s = d2j2s + d2j2r(iv)
  enddo
  d1iz = d1iz + d2trm1diz
  d1jz = d1jz + d2trm1djz
  d2ijs = d2ijs + d2trm1dij
  if (leem) then
    if (nat(i).gt.maxele) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp
    endif
  endif
!
!  Add d1ix / d1iy / d1iz to sums for globalisation later
!
  dedqc(1,i) = dedqc(1,i) + d1ix
  dedqc(2,i) = dedqc(2,i) + d1iy
  dedqc(3,i) = dedqc(3,i) + d1iz
!
!  Add d1jx / d1iy / d1jz to sums for globalisation later
!
  d2edqc(ixf,j) = d2edqc(ixf,j) + d1jx
  d2edqc(iyf,j) = d2edqc(iyf,j) + d1jy
  d2edqc(izf,j) = d2edqc(izf,j) + d1jz
!******************
!  New algorithm  *
!******************
  d2edqdq(j,iloc) = d2edqdq(j,iloc) + d2ijs
  d2edq2(i) = d2edq2(i) + d2i2s
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff.and..not.lreaxFF) then
    deisum = dtrm1di
    dejsum = dtrm1dj
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(iloc)
      k = nqatomptr(m,iloc)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,i)
      derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,i)
      derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,i)
      derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,i)
      derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,i)
      derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,i)
      derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,i)
      derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,i)
      derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,i)
    enddo
  endif
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargep2Dd')
#endif
!
  return
  end

  subroutine d2charge3fc(i,j,k,nor,ix,iy,iz,jx,jy,jz,kx,ky,kz,dei,dej,dek,d1q,d2q)
!
!  Calculates the contribution to the dynamical matrix due to charge derivatives from variable 
!  charge models as a result of three body potentials.
!  Version for unphased second derivatives.
!
!  10/14 Created from d2charge3p
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   7/21 QEq terms corrected 
!   9/23 Trapping of shells added
!
!  On entry:
!
!  d1q(3,3,3) = derivative of first derivative of i - j vector with respect to charge
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
  use control
  use current
  use derivatives
  use element
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: i
  integer(i4), intent(in)    :: ix
  integer(i4), intent(in)    :: iy
  integer(i4), intent(in)    :: iz
  integer(i4), intent(in)    :: j
  integer(i4), intent(in)    :: jx
  integer(i4), intent(in)    :: jy
  integer(i4), intent(in)    :: jz
  integer(i4), intent(in)    :: k
  integer(i4), intent(in)    :: kx
  integer(i4), intent(in)    :: ky
  integer(i4), intent(in)    :: kz
  integer(i4), intent(in)    :: nor
  real(dp),    intent(in)    :: d1q(3,3,3)
  real(dp),    intent(in)    :: d2q(6,*)
  real(dp),    intent(in)    :: dei(*)
  real(dp),    intent(in)    :: dej(*)
  real(dp),    intent(in)    :: dek(*)
!
!  Local variables
!
  integer(i4)                :: iv
  integer(i4)                :: l
  integer(i4)                :: lx
  integer(i4)                :: ly
  integer(i4)                :: lz
  integer(i4)                :: m
  integer(i4)                :: mnxx
  integer(i4)                :: mnxy
  integer(i4)                :: mnxz
  integer(i4)                :: mnyx
  integer(i4)                :: mnyy
  integer(i4)                :: mnyz
  integer(i4)                :: mnzx
  integer(i4)                :: mnzy
  integer(i4)                :: mnzz
  integer(i4)                :: mx
  integer(i4)                :: my
  integer(i4)                :: mz
  integer(i4)                :: n
  integer(i4)                :: nx
  integer(i4)                :: p
  integer(i4)                :: px
  integer(i4)                :: py
  integer(i4)                :: pz
  logical                    :: lNonIJQDeriv
  real(dp)                   :: d1ql(3,3,3)
  real(dp)                   :: d2i2s
  real(dp)                   :: d2ijs
  real(dp)                   :: d2iks
  real(dp)                   :: d2j2s
  real(dp)                   :: d2jks
  real(dp)                   :: d2k2s
  real(dp)                   :: deisum
  real(dp)                   :: dejsum
  real(dp)                   :: deksum
  real(dp)                   :: dilx
  real(dp)                   :: dily
  real(dp)                   :: dilz
  real(dp)                   :: dimx
  real(dp)                   :: dimy
  real(dp)                   :: dimz
  real(dp)                   :: djlx
  real(dp)                   :: djly
  real(dp)                   :: djlz
  real(dp)                   :: djmx
  real(dp)                   :: djmy
  real(dp)                   :: djmz
  real(dp)                   :: dklx
  real(dp)                   :: dkly
  real(dp)                   :: dklz
  real(dp)                   :: dkmx
  real(dp)                   :: dkmy
  real(dp)                   :: dkmz
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2charge3fc')
#endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d2i2s = 0.0_dp
  d2ijs = 0.0_dp
  d2iks = 0.0_dp
  d2j2s = 0.0_dp
  d2jks = 0.0_dp
  d2k2s = 0.0_dp
  do iv = 1,nor
    d2i2s = d2i2s + d2q(1,iv)
    d2ijs = d2ijs + d2q(2,iv)
    d2iks = d2iks + d2q(3,iv)
    d2j2s = d2j2s + d2q(4,iv)
    d2jks = d2jks + d2q(5,iv)
    d2k2s = d2k2s + d2q(6,iv)
  enddo
  d1ql(1:3,1:3,1:3) = d1q(1:3,1:3,1:3)
  if (leem) then
    if (nat(i).gt.maxele) then
      d2i2s = 0.0_dp
      d2ijs = 0.0_dp
      d2iks = 0.0_dp
      d1ql(1:3,1:3,1) = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2i2s = 0.0_dp
      d2ijs = 0.0_dp
      d2iks = 0.0_dp
      d1ql(1:3,1:3,1) = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d2jks = 0.0_dp
      d1ql(1:3,1:3,2) = 0.0_dp
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d2jks = 0.0_dp
      d1ql(1:3,1:3,2) = 0.0_dp
    endif
    if (nat(k).gt.maxele) then
      d2iks = 0.0_dp
      d2jks = 0.0_dp
      d2k2s = 0.0_dp
      d1ql(1:3,1:3,3) = 0.0_dp
    elseif (.not.lelementOK(nat(k)).or.nregionno(nrelf2a(k)).ne.1) then
      d2iks = 0.0_dp
      d2jks = 0.0_dp
      d2k2s = 0.0_dp
      d1ql(1:3,1:3,3) = 0.0_dp
    endif
  endif
!
!  The following flag causes the second derivatives of the charge distribution 
!  with respect to two atoms neither of which are the i/j of the energy term
!  to be calculated. At the moment this contribution is not used in any method
!  and so the flag is set to false for efficiency.
!
  lNonIJQDeriv = .false.
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lreaxFF) then
    deisum = 0.0_dp
    dejsum = 0.0_dp
    deksum = 0.0_dp
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
      deksum = deksum + dek(iv)
    enddo
!
!  Terms involving Q derivatives for i/j/k and one other atom
!
    do m = 1,nqatoms(i)
      p = nqatomptr(m,i)
      px = 3*(p - 1) + 1
      py = px + 1
      pz = py + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (i.gt.p) then
        derv2(px,ix) = derv2(px,ix) + deisum*d2qdxyz2(mnxx,i)
        derv2(py,ix) = derv2(py,ix) + deisum*d2qdxyz2(mnxy,i)
        derv2(pz,ix) = derv2(pz,ix) + deisum*d2qdxyz2(mnxz,i)
        derv2(px,iy) = derv2(px,iy) + deisum*d2qdxyz2(mnxy,i)
        derv2(py,iy) = derv2(py,iy) + deisum*d2qdxyz2(mnyy,i)
        derv2(pz,iy) = derv2(pz,iy) + deisum*d2qdxyz2(mnyz,i)
        derv2(px,iz) = derv2(px,iz) + deisum*d2qdxyz2(mnxz,i)
        derv2(py,iz) = derv2(py,iz) + deisum*d2qdxyz2(mnyz,i)
        derv2(pz,iz) = derv2(pz,iz) + deisum*d2qdxyz2(mnzz,i)
      elseif (i.lt.p) then
        derv2(ix,px) = derv2(ix,px) + deisum*d2qdxyz2(mnxx,i)
        derv2(iy,px) = derv2(iy,px) + deisum*d2qdxyz2(mnxy,i)
        derv2(iz,px) = derv2(iz,px) + deisum*d2qdxyz2(mnxz,i)
        derv2(ix,py) = derv2(ix,py) + deisum*d2qdxyz2(mnxy,i)
        derv2(iy,py) = derv2(iy,py) + deisum*d2qdxyz2(mnyy,i)
        derv2(iz,py) = derv2(iz,py) + deisum*d2qdxyz2(mnyz,i)
        derv2(ix,pz) = derv2(ix,pz) + deisum*d2qdxyz2(mnxz,i)
        derv2(iy,pz) = derv2(iy,pz) + deisum*d2qdxyz2(mnyz,i)
        derv2(iz,pz) = derv2(iz,pz) + deisum*d2qdxyz2(mnzz,i)
      endif
    enddo
!
    do m = 1,nqatoms(j)
      p = nqatomptr(m,j)
      px = 3*(p - 1) + 1
      py = px + 1
      pz = py + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (j.gt.p) then
        derv2(px,jx) = derv2(px,jx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(py,jx) = derv2(py,jx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(pz,jx) = derv2(pz,jx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(px,jy) = derv2(px,jy) + dejsum*d2qdxyz2(mnxy,j)
        derv2(py,jy) = derv2(py,jy) + dejsum*d2qdxyz2(mnyy,j)
        derv2(pz,jy) = derv2(pz,jy) + dejsum*d2qdxyz2(mnyz,j)
        derv2(px,jz) = derv2(px,jz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(py,jz) = derv2(py,jz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(pz,jz) = derv2(pz,jz) + dejsum*d2qdxyz2(mnzz,j)
      elseif (j.lt.p) then
        derv2(jx,px) = derv2(jx,px) + dejsum*d2qdxyz2(mnxx,j)
        derv2(jy,px) = derv2(jy,px) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jz,px) = derv2(jz,px) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jx,py) = derv2(jx,py) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jy,py) = derv2(jy,py) + dejsum*d2qdxyz2(mnyy,j)
        derv2(jz,py) = derv2(jz,py) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jx,pz) = derv2(jx,pz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jy,pz) = derv2(jy,pz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jz,pz) = derv2(jz,pz) + dejsum*d2qdxyz2(mnzz,j)
      endif
    enddo
!
    do m = 1,nqatoms(k)
      p = nqatomptr(m,k)
      px = 3*(p - 1) + 1
      py = px + 1
      pz = py + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (k.gt.p) then
        derv2(px,kx) = derv2(px,kx) + deksum*d2qdxyz2(mnxx,k)
        derv2(py,kx) = derv2(py,kx) + deksum*d2qdxyz2(mnxy,k)
        derv2(pz,kx) = derv2(pz,kx) + deksum*d2qdxyz2(mnxz,k)
        derv2(px,ky) = derv2(px,ky) + deksum*d2qdxyz2(mnxy,k)
        derv2(py,ky) = derv2(py,ky) + deksum*d2qdxyz2(mnyy,k)
        derv2(pz,ky) = derv2(pz,ky) + deksum*d2qdxyz2(mnyz,k)
        derv2(px,kz) = derv2(px,kz) + deksum*d2qdxyz2(mnxz,k)
        derv2(py,kz) = derv2(py,kz) + deksum*d2qdxyz2(mnyz,k)
        derv2(pz,kz) = derv2(pz,kz) + deksum*d2qdxyz2(mnzz,k)
      elseif (k.lt.p) then
        derv2(kx,px) = derv2(kx,px) + deksum*d2qdxyz2(mnxx,k)
        derv2(ky,px) = derv2(ky,px) + deksum*d2qdxyz2(mnxy,k)
        derv2(kz,px) = derv2(kz,px) + deksum*d2qdxyz2(mnxz,k)
        derv2(kx,py) = derv2(kx,py) + deksum*d2qdxyz2(mnxy,k)
        derv2(ky,py) = derv2(ky,py) + deksum*d2qdxyz2(mnyy,k)
        derv2(kz,py) = derv2(kz,py) + deksum*d2qdxyz2(mnyz,k)
        derv2(kx,pz) = derv2(kx,pz) + deksum*d2qdxyz2(mnxz,k)
        derv2(ky,pz) = derv2(ky,pz) + deksum*d2qdxyz2(mnyz,k)
        derv2(kz,pz) = derv2(kz,pz) + deksum*d2qdxyz2(mnzz,k)
      endif
    enddo
!
!  Terms involving Q derivatives for non-i/j case
!
    if (lNonIJQDeriv) then
      do m = 2,nqatoms(i)
        p = nqatomptr(m,i)
        px = 3*(p - 1) + 1
        py = px + 1
        pz = py + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,i)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,px) = derv2(lx,px) + deisum*d2qdxyz2(mnxx,i)
          derv2(ly,px) = derv2(ly,px) + deisum*d2qdxyz2(mnxy,i)
          derv2(lz,px) = derv2(lz,px) + deisum*d2qdxyz2(mnxz,i)
          derv2(lx,py) = derv2(lx,py) + deisum*d2qdxyz2(mnyx,i)
          derv2(ly,py) = derv2(ly,py) + deisum*d2qdxyz2(mnyy,i)
          derv2(lz,py) = derv2(lz,py) + deisum*d2qdxyz2(mnyz,i)
          derv2(lx,pz) = derv2(lx,pz) + deisum*d2qdxyz2(mnzx,i)
          derv2(ly,pz) = derv2(ly,pz) + deisum*d2qdxyz2(mnzy,i)
          derv2(lz,pz) = derv2(lz,pz) + deisum*d2qdxyz2(mnzz,i)
        enddo
      enddo
!
      do m = 2,nqatoms(j)
        p = nqatomptr(m,j)
        px = 3*(p - 1) + 1
        py = px + 1
        pz = py + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,j)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,px) = derv2(lx,px) + deisum*d2qdxyz2(mnxx,j)
          derv2(ly,px) = derv2(ly,px) + deisum*d2qdxyz2(mnxy,j)
          derv2(lz,px) = derv2(lz,px) + deisum*d2qdxyz2(mnxz,j)
          derv2(lx,py) = derv2(lx,py) + deisum*d2qdxyz2(mnyx,j)
          derv2(ly,py) = derv2(ly,py) + deisum*d2qdxyz2(mnyy,j)
          derv2(lz,py) = derv2(lz,py) + deisum*d2qdxyz2(mnyz,j)
          derv2(lx,pz) = derv2(lx,pz) + deisum*d2qdxyz2(mnzx,j)
          derv2(ly,pz) = derv2(ly,pz) + deisum*d2qdxyz2(mnzy,j)
          derv2(lz,pz) = derv2(lz,pz) + deisum*d2qdxyz2(mnzz,j)
        enddo
      enddo
!
      do m = 2,nqatoms(k)
        p = nqatomptr(m,k)
        px = 3*(p - 1) + 1
        py = px + 1
        pz = py + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,k)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,px) = derv2(lx,px) + deksum*d2qdxyz2(mnxx,k)
          derv2(ly,px) = derv2(ly,px) + deksum*d2qdxyz2(mnxy,k)
          derv2(lz,px) = derv2(lz,px) + deksum*d2qdxyz2(mnxz,k)
          derv2(lx,py) = derv2(lx,py) + deksum*d2qdxyz2(mnyx,k)
          derv2(ly,py) = derv2(ly,py) + deksum*d2qdxyz2(mnyy,k)
          derv2(lz,py) = derv2(lz,py) + deksum*d2qdxyz2(mnyz,k)
          derv2(lx,pz) = derv2(lx,pz) + deksum*d2qdxyz2(mnzx,k)
          derv2(ly,pz) = derv2(ly,pz) + deksum*d2qdxyz2(mnzy,k)
          derv2(lz,pz) = derv2(lz,pz) + deksum*d2qdxyz2(mnzz,k)
        enddo
      enddo
    endif
  endif
!
!  Loop over atoms to add ij-m contributions
!
  mx = - 2
  my = - 1
  mz =   0
  do m = 1,numat
    mx = mx + 3
    my = my + 3
    mz = mz + 3
!
    dimx = dqdxyz(mx,i)
    dimy = dqdxyz(my,i)
    dimz = dqdxyz(mz,i)
    djmx = dqdxyz(mx,j)
    djmy = dqdxyz(my,j)
    djmz = dqdxyz(mz,j)
    dkmx = dqdxyz(mx,k)
    dkmy = dqdxyz(my,k)
    dkmz = dqdxyz(mz,k)
    if (i.ne.m) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-m
!
      if (m.le.i) then
        derv2(mx,ix) = derv2(mx,ix) - d1ql(1,1,1)*dimx - d1ql(1,2,1)*dimx
        derv2(my,ix) = derv2(my,ix) - d1ql(1,1,1)*dimy - d1ql(1,2,1)*dimy
        derv2(mz,ix) = derv2(mz,ix) - d1ql(1,1,1)*dimz - d1ql(1,2,1)*dimz
        derv2(mx,iy) = derv2(mx,iy) - d1ql(2,1,1)*dimx - d1ql(2,2,1)*dimx
        derv2(my,iy) = derv2(my,iy) - d1ql(2,1,1)*dimy - d1ql(2,2,1)*dimy
        derv2(mz,iy) = derv2(mz,iy) - d1ql(2,1,1)*dimz - d1ql(2,2,1)*dimz
        derv2(mx,iz) = derv2(mx,iz) - d1ql(3,1,1)*dimx - d1ql(3,2,1)*dimx
        derv2(my,iz) = derv2(my,iz) - d1ql(3,1,1)*dimy - d1ql(3,2,1)*dimy
        derv2(mz,iz) = derv2(mz,iz) - d1ql(3,1,1)*dimz - d1ql(3,2,1)*dimz
!
        derv2(mx,ix) = derv2(mx,ix) - d1ql(1,1,2)*djmx - d1ql(1,2,2)*djmx
        derv2(my,ix) = derv2(my,ix) - d1ql(1,1,2)*djmy - d1ql(1,2,2)*djmy
        derv2(mz,ix) = derv2(mz,ix) - d1ql(1,1,2)*djmz - d1ql(1,2,2)*djmz
        derv2(mx,iy) = derv2(mx,iy) - d1ql(2,1,2)*djmx - d1ql(2,2,2)*djmx
        derv2(my,iy) = derv2(my,iy) - d1ql(2,1,2)*djmy - d1ql(2,2,2)*djmy
        derv2(mz,iy) = derv2(mz,iy) - d1ql(2,1,2)*djmz - d1ql(2,2,2)*djmz
        derv2(mx,iz) = derv2(mx,iz) - d1ql(3,1,2)*djmx - d1ql(3,2,2)*djmx
        derv2(my,iz) = derv2(my,iz) - d1ql(3,1,2)*djmy - d1ql(3,2,2)*djmy
        derv2(mz,iz) = derv2(mz,iz) - d1ql(3,1,2)*djmz - d1ql(3,2,2)*djmz
!
        derv2(mx,ix) = derv2(mx,ix) - d1ql(1,1,3)*dkmx - d1ql(1,2,3)*dkmx
        derv2(my,ix) = derv2(my,ix) - d1ql(1,1,3)*dkmy - d1ql(1,2,3)*dkmy
        derv2(mz,ix) = derv2(mz,ix) - d1ql(1,1,3)*dkmz - d1ql(1,2,3)*dkmz
        derv2(mx,iy) = derv2(mx,iy) - d1ql(2,1,3)*dkmx - d1ql(2,2,3)*dkmx
        derv2(my,iy) = derv2(my,iy) - d1ql(2,1,3)*dkmy - d1ql(2,2,3)*dkmy
        derv2(mz,iy) = derv2(mz,iy) - d1ql(2,1,3)*dkmz - d1ql(2,2,3)*dkmz
        derv2(mx,iz) = derv2(mx,iz) - d1ql(3,1,3)*dkmx - d1ql(3,2,3)*dkmx
        derv2(my,iz) = derv2(my,iz) - d1ql(3,1,3)*dkmy - d1ql(3,2,3)*dkmy
        derv2(mz,iz) = derv2(mz,iz) - d1ql(3,1,3)*dkmz - d1ql(3,2,3)*dkmz
      else
        derv2(ix,mx) = derv2(ix,mx) - d1ql(1,1,1)*dimx - d1ql(1,2,1)*dimx
        derv2(iy,mx) = derv2(iy,mx) - d1ql(2,1,1)*dimx - d1ql(2,2,1)*dimx
        derv2(iz,mx) = derv2(iz,mx) - d1ql(3,1,1)*dimx - d1ql(3,2,1)*dimx
        derv2(ix,my) = derv2(ix,my) - d1ql(1,1,1)*dimy - d1ql(1,2,1)*dimy
        derv2(iy,my) = derv2(iy,my) - d1ql(2,1,1)*dimy - d1ql(2,2,1)*dimy
        derv2(iz,my) = derv2(iz,my) - d1ql(3,1,1)*dimy - d1ql(3,2,1)*dimy
        derv2(ix,mz) = derv2(ix,mz) - d1ql(1,1,1)*dimz - d1ql(1,2,1)*dimz
        derv2(iy,mz) = derv2(iy,mz) - d1ql(2,1,1)*dimz - d1ql(2,2,1)*dimz
        derv2(iz,mz) = derv2(iz,mz) - d1ql(3,1,1)*dimz - d1ql(3,2,1)*dimz
!
        derv2(ix,mx) = derv2(ix,mx) - d1ql(1,1,2)*djmx - d1ql(1,2,2)*djmx
        derv2(iy,mx) = derv2(iy,mx) - d1ql(2,1,2)*djmx - d1ql(2,2,2)*djmx
        derv2(iz,mx) = derv2(iz,mx) - d1ql(3,1,2)*djmx - d1ql(3,2,2)*djmx
        derv2(ix,my) = derv2(ix,my) - d1ql(1,1,2)*djmy - d1ql(1,2,2)*djmy
        derv2(iy,my) = derv2(iy,my) - d1ql(2,1,2)*djmy - d1ql(2,2,2)*djmy
        derv2(iz,my) = derv2(iz,my) - d1ql(3,1,2)*djmy - d1ql(3,2,2)*djmy
        derv2(ix,mz) = derv2(ix,mz) - d1ql(1,1,2)*djmz - d1ql(1,2,2)*djmz
        derv2(iy,mz) = derv2(iy,mz) - d1ql(2,1,2)*djmz - d1ql(2,2,2)*djmz
        derv2(iz,mz) = derv2(iz,mz) - d1ql(3,1,2)*djmz - d1ql(3,2,2)*djmz
!
        derv2(ix,mx) = derv2(ix,mx) - d1ql(1,1,3)*dkmx - d1ql(1,2,3)*dkmx
        derv2(iy,mx) = derv2(iy,mx) - d1ql(2,1,3)*dkmx - d1ql(2,2,3)*dkmx
        derv2(iz,mx) = derv2(iz,mx) - d1ql(3,1,3)*dkmx - d1ql(3,2,3)*dkmx
        derv2(ix,my) = derv2(ix,my) - d1ql(1,1,3)*dkmy - d1ql(1,2,3)*dkmy
        derv2(iy,my) = derv2(iy,my) - d1ql(2,1,3)*dkmy - d1ql(2,2,3)*dkmy
        derv2(iz,my) = derv2(iz,my) - d1ql(3,1,3)*dkmy - d1ql(3,2,3)*dkmy
        derv2(ix,mz) = derv2(ix,mz) - d1ql(1,1,3)*dkmz - d1ql(1,2,3)*dkmz
        derv2(iy,mz) = derv2(iy,mz) - d1ql(2,1,3)*dkmz - d1ql(2,2,3)*dkmz
        derv2(iz,mz) = derv2(iz,mz) - d1ql(3,1,3)*dkmz - d1ql(3,2,3)*dkmz
      endif
    endif
    if (j.ne.m) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : j-m
!
      if (m.le.j) then
        derv2(mx,jx) = derv2(mx,jx) + d1ql(1,1,1)*dimx - d1ql(1,3,1)*dimx
        derv2(my,jx) = derv2(my,jx) + d1ql(1,1,1)*dimy - d1ql(1,3,1)*dimy
        derv2(mz,jx) = derv2(mz,jx) + d1ql(1,1,1)*dimz - d1ql(1,3,1)*dimz
        derv2(mx,jy) = derv2(mx,jy) + d1ql(2,1,1)*dimx - d1ql(2,3,1)*dimx
        derv2(my,jy) = derv2(my,jy) + d1ql(2,1,1)*dimy - d1ql(2,3,1)*dimy
        derv2(mz,jy) = derv2(mz,jy) + d1ql(2,1,1)*dimz - d1ql(2,3,1)*dimz
        derv2(mx,jz) = derv2(mx,jz) + d1ql(3,1,1)*dimx - d1ql(3,3,1)*dimx
        derv2(my,jz) = derv2(my,jz) + d1ql(3,1,1)*dimy - d1ql(3,3,1)*dimy
        derv2(mz,jz) = derv2(mz,jz) + d1ql(3,1,1)*dimz - d1ql(3,3,1)*dimz
!
        derv2(mx,jx) = derv2(mx,jx) + d1ql(1,1,2)*djmx - d1ql(1,3,2)*djmx
        derv2(my,jx) = derv2(my,jx) + d1ql(1,1,2)*djmy - d1ql(1,3,2)*djmy
        derv2(mz,jx) = derv2(mz,jx) + d1ql(1,1,2)*djmz - d1ql(1,3,2)*djmz
        derv2(mx,jy) = derv2(mx,jy) + d1ql(2,1,2)*djmx - d1ql(2,3,2)*djmx
        derv2(my,jy) = derv2(my,jy) + d1ql(2,1,2)*djmy - d1ql(2,3,2)*djmy
        derv2(mz,jy) = derv2(mz,jy) + d1ql(2,1,2)*djmz - d1ql(2,3,2)*djmz
        derv2(mx,jz) = derv2(mx,jz) + d1ql(3,1,2)*djmx - d1ql(3,3,2)*djmx
        derv2(my,jz) = derv2(my,jz) + d1ql(3,1,2)*djmy - d1ql(3,3,2)*djmy
        derv2(mz,jz) = derv2(mz,jz) + d1ql(3,1,2)*djmz - d1ql(3,3,2)*djmz
!
        derv2(mx,jx) = derv2(mx,jx) + d1ql(1,1,3)*dkmx - d1ql(1,3,3)*dkmx
        derv2(my,jx) = derv2(my,jx) + d1ql(1,1,3)*dkmy - d1ql(1,3,3)*dkmy
        derv2(mz,jx) = derv2(mz,jx) + d1ql(1,1,3)*dkmz - d1ql(1,3,3)*dkmz
        derv2(mx,jy) = derv2(mx,jy) + d1ql(2,1,3)*dkmx - d1ql(2,3,3)*dkmx
        derv2(my,jy) = derv2(my,jy) + d1ql(2,1,3)*dkmy - d1ql(2,3,3)*dkmy
        derv2(mz,jy) = derv2(mz,jy) + d1ql(2,1,3)*dkmz - d1ql(2,3,3)*dkmz
        derv2(mx,jz) = derv2(mx,jz) + d1ql(3,1,3)*dkmx - d1ql(3,3,3)*dkmx
        derv2(my,jz) = derv2(my,jz) + d1ql(3,1,3)*dkmy - d1ql(3,3,3)*dkmy
        derv2(mz,jz) = derv2(mz,jz) + d1ql(3,1,3)*dkmz - d1ql(3,3,3)*dkmz
      else
        derv2(jx,mx) = derv2(jx,mx) + d1ql(1,1,1)*dimx - d1ql(1,3,1)*dimx
        derv2(jy,mx) = derv2(jy,mx) + d1ql(2,1,1)*dimx - d1ql(2,3,1)*dimx
        derv2(jz,mx) = derv2(jz,mx) + d1ql(3,1,1)*dimx - d1ql(3,3,1)*dimx
        derv2(jx,my) = derv2(jx,my) + d1ql(1,1,1)*dimy - d1ql(1,3,1)*dimy
        derv2(jy,my) = derv2(jy,my) + d1ql(2,1,1)*dimy - d1ql(2,3,1)*dimy
        derv2(jz,my) = derv2(jz,my) + d1ql(3,1,1)*dimy - d1ql(3,3,1)*dimy
        derv2(jx,mz) = derv2(jx,mz) + d1ql(1,1,1)*dimz - d1ql(1,3,1)*dimz
        derv2(jy,mz) = derv2(jy,mz) + d1ql(2,1,1)*dimz - d1ql(2,3,1)*dimz
        derv2(jz,mz) = derv2(jz,mz) + d1ql(3,1,1)*dimz - d1ql(3,3,1)*dimz
!
        derv2(jx,mx) = derv2(jx,mx) + d1ql(1,1,2)*djmx - d1ql(1,3,2)*djmx
        derv2(jy,mx) = derv2(jy,mx) + d1ql(2,1,2)*djmx - d1ql(2,3,2)*djmx
        derv2(jz,mx) = derv2(jz,mx) + d1ql(3,1,2)*djmx - d1ql(3,3,2)*djmx
        derv2(jx,my) = derv2(jx,my) + d1ql(1,1,2)*djmy - d1ql(1,3,2)*djmy
        derv2(jy,my) = derv2(jy,my) + d1ql(2,1,2)*djmy - d1ql(2,3,2)*djmy
        derv2(jz,my) = derv2(jz,my) + d1ql(3,1,2)*djmy - d1ql(3,3,2)*djmy
        derv2(jx,mz) = derv2(jx,mz) + d1ql(1,1,2)*djmz - d1ql(1,3,2)*djmz
        derv2(jy,mz) = derv2(jy,mz) + d1ql(2,1,2)*djmz - d1ql(2,3,2)*djmz
        derv2(jz,mz) = derv2(jz,mz) + d1ql(3,1,2)*djmz - d1ql(3,3,2)*djmz
!
        derv2(jx,mx) = derv2(jx,mx) + d1ql(1,1,3)*dkmx - d1ql(1,3,3)*dkmx
        derv2(jy,mx) = derv2(jy,mx) + d1ql(2,1,3)*dkmx - d1ql(2,3,3)*dkmx
        derv2(jz,mx) = derv2(jz,mx) + d1ql(3,1,3)*dkmx - d1ql(3,3,3)*dkmx
        derv2(jx,my) = derv2(jx,my) + d1ql(1,1,3)*dkmy - d1ql(1,3,3)*dkmy
        derv2(jy,my) = derv2(jy,my) + d1ql(2,1,3)*dkmy - d1ql(2,3,3)*dkmy
        derv2(jz,my) = derv2(jz,my) + d1ql(3,1,3)*dkmy - d1ql(3,3,3)*dkmy
        derv2(jx,mz) = derv2(jx,mz) + d1ql(1,1,3)*dkmz - d1ql(1,3,3)*dkmz
        derv2(jy,mz) = derv2(jy,mz) + d1ql(2,1,3)*dkmz - d1ql(2,3,3)*dkmz
        derv2(jz,mz) = derv2(jz,mz) + d1ql(3,1,3)*dkmz - d1ql(3,3,3)*dkmz
      endif
    endif
    if (k.ne.m) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : k-m
!
      if (m.le.k) then
        derv2(mx,kx) = derv2(mx,kx) + d1ql(1,2,1)*dimx + d1ql(1,3,1)*dimx
        derv2(my,kx) = derv2(my,kx) + d1ql(1,2,1)*dimy + d1ql(1,3,1)*dimy
        derv2(mz,kx) = derv2(mz,kx) + d1ql(1,2,1)*dimz + d1ql(1,3,1)*dimz
        derv2(mx,ky) = derv2(mx,ky) + d1ql(2,2,1)*dimx + d1ql(2,3,1)*dimx
        derv2(my,ky) = derv2(my,ky) + d1ql(2,2,1)*dimy + d1ql(2,3,1)*dimy
        derv2(mz,ky) = derv2(mz,ky) + d1ql(2,2,1)*dimz + d1ql(2,3,1)*dimz
        derv2(mx,kz) = derv2(mx,kz) + d1ql(3,2,1)*dimx + d1ql(3,3,1)*dimx
        derv2(my,kz) = derv2(my,kz) + d1ql(3,2,1)*dimy + d1ql(3,3,1)*dimy
        derv2(mz,kz) = derv2(mz,kz) + d1ql(3,2,1)*dimz + d1ql(3,3,1)*dimz
!
        derv2(mx,kx) = derv2(mx,kx) + d1ql(1,2,2)*djmx + d1ql(1,3,2)*djmx
        derv2(my,kx) = derv2(my,kx) + d1ql(1,2,2)*djmy + d1ql(1,3,2)*djmy
        derv2(mz,kx) = derv2(mz,kx) + d1ql(1,2,2)*djmz + d1ql(1,3,2)*djmz
        derv2(mx,ky) = derv2(mx,ky) + d1ql(2,2,2)*djmx + d1ql(2,3,2)*djmx
        derv2(my,ky) = derv2(my,ky) + d1ql(2,2,2)*djmy + d1ql(2,3,2)*djmy
        derv2(mz,ky) = derv2(mz,ky) + d1ql(2,2,2)*djmz + d1ql(2,3,2)*djmz
        derv2(mx,kz) = derv2(mx,kz) + d1ql(3,2,2)*djmx + d1ql(3,3,2)*djmx
        derv2(my,kz) = derv2(my,kz) + d1ql(3,2,2)*djmy + d1ql(3,3,2)*djmy
        derv2(mz,kz) = derv2(mz,kz) + d1ql(3,2,2)*djmz + d1ql(3,3,2)*djmz
!
        derv2(mx,kx) = derv2(mx,kx) + d1ql(1,2,3)*dkmx + d1ql(1,3,3)*dkmx
        derv2(my,kx) = derv2(my,kx) + d1ql(1,2,3)*dkmy + d1ql(1,3,3)*dkmy
        derv2(mz,kx) = derv2(mz,kx) + d1ql(1,2,3)*dkmz + d1ql(1,3,3)*dkmz
        derv2(mx,ky) = derv2(mx,ky) + d1ql(2,2,3)*dkmx + d1ql(2,3,3)*dkmx
        derv2(my,ky) = derv2(my,ky) + d1ql(2,2,3)*dkmy + d1ql(2,3,3)*dkmy
        derv2(mz,ky) = derv2(mz,ky) + d1ql(2,2,3)*dkmz + d1ql(2,3,3)*dkmz
        derv2(mx,kz) = derv2(mx,kz) + d1ql(3,2,3)*dkmx + d1ql(3,3,3)*dkmx
        derv2(my,kz) = derv2(my,kz) + d1ql(3,2,3)*dkmy + d1ql(3,3,3)*dkmy
        derv2(mz,kz) = derv2(mz,kz) + d1ql(3,2,3)*dkmz + d1ql(3,3,3)*dkmz
      else
        derv2(kx,mx) = derv2(kx,mx) + d1ql(1,2,1)*dimx + d1ql(1,3,1)*dimx
        derv2(ky,mx) = derv2(ky,mx) + d1ql(2,2,1)*dimx + d1ql(2,3,1)*dimx
        derv2(kz,mx) = derv2(kz,mx) + d1ql(3,2,1)*dimx + d1ql(3,3,1)*dimx
        derv2(kx,my) = derv2(kx,my) + d1ql(1,2,1)*dimy + d1ql(1,3,1)*dimy
        derv2(ky,my) = derv2(ky,my) + d1ql(2,2,1)*dimy + d1ql(2,3,1)*dimy
        derv2(kz,my) = derv2(kz,my) + d1ql(3,2,1)*dimy + d1ql(3,3,1)*dimy
        derv2(kx,mz) = derv2(kx,mz) + d1ql(1,2,1)*dimz + d1ql(1,3,1)*dimz
        derv2(ky,mz) = derv2(ky,mz) + d1ql(2,2,1)*dimz + d1ql(2,3,1)*dimz
        derv2(kz,mz) = derv2(kz,mz) + d1ql(3,2,1)*dimz + d1ql(3,3,1)*dimz
!
        derv2(kx,mx) = derv2(kx,mx) + d1ql(1,2,2)*djmx + d1ql(1,3,2)*djmx
        derv2(ky,mx) = derv2(ky,mx) + d1ql(2,2,2)*djmx + d1ql(2,3,2)*djmx
        derv2(kz,mx) = derv2(kz,mx) + d1ql(3,2,2)*djmx + d1ql(3,3,2)*djmx
        derv2(kx,my) = derv2(kx,my) + d1ql(1,2,2)*djmy + d1ql(1,3,2)*djmy
        derv2(ky,my) = derv2(ky,my) + d1ql(2,2,2)*djmy + d1ql(2,3,2)*djmy
        derv2(kz,my) = derv2(kz,my) + d1ql(3,2,2)*djmy + d1ql(3,3,2)*djmy
        derv2(kx,mz) = derv2(kx,mz) + d1ql(1,2,2)*djmz + d1ql(1,3,2)*djmz
        derv2(ky,mz) = derv2(ky,mz) + d1ql(2,2,2)*djmz + d1ql(2,3,2)*djmz
        derv2(kz,mz) = derv2(kz,mz) + d1ql(3,2,2)*djmz + d1ql(3,3,2)*djmz
!
        derv2(kx,mx) = derv2(kx,mx) + d1ql(1,2,3)*dkmx + d1ql(1,3,3)*dkmx
        derv2(ky,mx) = derv2(ky,mx) + d1ql(2,2,3)*dkmx + d1ql(2,3,3)*dkmx
        derv2(kz,mx) = derv2(kz,mx) + d1ql(3,2,3)*dkmx + d1ql(3,3,3)*dkmx
        derv2(kx,my) = derv2(kx,my) + d1ql(1,2,3)*dkmy + d1ql(1,3,3)*dkmy
        derv2(ky,my) = derv2(ky,my) + d1ql(2,2,3)*dkmy + d1ql(2,3,3)*dkmy
        derv2(kz,my) = derv2(kz,my) + d1ql(3,2,3)*dkmy + d1ql(3,3,3)*dkmy
        derv2(kx,mz) = derv2(kx,mz) + d1ql(1,2,3)*dkmz + d1ql(1,3,3)*dkmz
        derv2(ky,mz) = derv2(ky,mz) + d1ql(2,2,3)*dkmz + d1ql(2,3,3)*dkmz
        derv2(kz,mz) = derv2(kz,mz) + d1ql(3,2,3)*dkmz + d1ql(3,3,3)*dkmz
      endif
    endif
!
!  Loop over atoms to add ij-ml correction
!
    lx = - 2
    ly = - 1
    lz =   0
    do l = 1,m-1
      lx = lx + 3
      ly = ly + 3
      lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
      dilx = dqdxyz(lx,i)
      dily = dqdxyz(ly,i)
      dilz = dqdxyz(lz,i)
      djlx = dqdxyz(lx,j)
      djly = dqdxyz(ly,j)
      djlz = dqdxyz(lz,j)
      dklx = dqdxyz(lx,k)
      dkly = dqdxyz(ly,k)
      dklz = dqdxyz(lz,k)
!
      derv2(lx,mx) = derv2(lx,mx) + d2ijs*dilx*djmx + d2iks*dilx*dkmx + d2jks*djlx*dkmx
      derv2(ly,mx) = derv2(ly,mx) + d2ijs*dilx*djmy + d2iks*dilx*dkmy + d2jks*djlx*dkmy
      derv2(lz,mx) = derv2(lz,mx) + d2ijs*dilx*djmz + d2iks*dilx*dkmz + d2jks*djlx*dkmz
      derv2(lx,my) = derv2(lx,my) + d2ijs*dily*djmx + d2iks*dily*dkmx + d2jks*djly*dkmx
      derv2(ly,my) = derv2(ly,my) + d2ijs*dily*djmy + d2iks*dily*dkmy + d2jks*djly*dkmy
      derv2(lz,my) = derv2(lz,my) + d2ijs*dily*djmz + d2iks*dily*dkmz + d2jks*djly*dkmz
      derv2(lx,mz) = derv2(lx,mz) + d2ijs*dilz*djmx + d2iks*dilz*dkmx + d2jks*djlz*dkmx
      derv2(ly,mz) = derv2(ly,mz) + d2ijs*dilz*djmy + d2iks*dilz*dkmy + d2jks*djlz*dkmy
      derv2(lz,mz) = derv2(lz,mz) + d2ijs*dilz*djmz + d2iks*dilz*dkmz + d2jks*djlz*dkmz
!
      derv2(lx,mx) = derv2(lx,mx) + d2ijs*djlx*dimx + d2iks*dklx*dimx + d2jks*dklx*djmx
      derv2(ly,mx) = derv2(ly,mx) + d2ijs*djlx*dimy + d2iks*dklx*dimy + d2jks*dklx*djmy
      derv2(lz,mx) = derv2(lz,mx) + d2ijs*djlx*dimz + d2iks*dklx*dimz + d2jks*dklx*djmz
      derv2(lx,my) = derv2(lx,my) + d2ijs*djly*dimx + d2iks*dkly*dimx + d2jks*dkly*djmx
      derv2(ly,my) = derv2(ly,my) + d2ijs*djly*dimy + d2iks*dkly*dimy + d2jks*dkly*djmy
      derv2(lz,my) = derv2(lz,my) + d2ijs*djly*dimz + d2iks*dkly*dimz + d2jks*dkly*djmz
      derv2(lx,mz) = derv2(lx,mz) + d2ijs*djlz*dimx + d2iks*dklz*dimx + d2jks*dklz*djmx
      derv2(ly,mz) = derv2(ly,mz) + d2ijs*djlz*dimy + d2iks*dklz*dimy + d2jks*dklz*djmy
      derv2(lz,mz) = derv2(lz,mz) + d2ijs*djlz*dimz + d2iks*dklz*dimz + d2jks*dklz*djmz
!
!  d2Edi2/d2Edj2/d2Edk2
!
      if (abs(d2i2s).gt.1.0d-8) then
        derv2(lx,mx) = derv2(lx,mx) + d2i2s*dilx*dimx
        derv2(ly,mx) = derv2(ly,mx) + d2i2s*dilx*dimy
        derv2(lz,mx) = derv2(lz,mx) + d2i2s*dilx*dimz
        derv2(lx,my) = derv2(lx,my) + d2i2s*dily*dimx
        derv2(ly,my) = derv2(ly,my) + d2i2s*dily*dimy
        derv2(lz,my) = derv2(lz,my) + d2i2s*dily*dimz
        derv2(lx,mz) = derv2(lx,mz) + d2i2s*dilz*dimx
        derv2(ly,mz) = derv2(ly,mz) + d2i2s*dilz*dimy
        derv2(lz,mz) = derv2(lz,mz) + d2i2s*dilz*dimz
      endif
      if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
        derv2(lx,mx) = derv2(lx,mx) + d2j2s*djlx*djmx
        derv2(ly,mx) = derv2(ly,mx) + d2j2s*djlx*djmy
        derv2(lz,mx) = derv2(lz,mx) + d2j2s*djlx*djmz
        derv2(lx,my) = derv2(lx,my) + d2j2s*djly*djmx
        derv2(ly,my) = derv2(ly,my) + d2j2s*djly*djmy
        derv2(lz,my) = derv2(lz,my) + d2j2s*djly*djmz
        derv2(lx,mz) = derv2(lx,mz) + d2j2s*djlz*djmx
        derv2(ly,mz) = derv2(ly,mz) + d2j2s*djlz*djmy
        derv2(lz,mz) = derv2(lz,mz) + d2j2s*djlz*djmz
      endif
      if (abs(d2k2s).gt.1.0d-8.and.i.ne.k.and.j.ne.k) then
        derv2(lx,mx) = derv2(lx,mx) + d2k2s*dklx*dkmx
        derv2(ly,mx) = derv2(ly,mx) + d2k2s*dklx*dkmy
        derv2(lz,mx) = derv2(lz,mx) + d2k2s*dklx*dkmz
        derv2(lx,my) = derv2(lx,my) + d2k2s*dkly*dkmx
        derv2(ly,my) = derv2(ly,my) + d2k2s*dkly*dkmy
        derv2(lz,my) = derv2(lz,my) + d2k2s*dkly*dkmz
        derv2(lx,mz) = derv2(lx,mz) + d2k2s*dklz*dkmx
        derv2(ly,mz) = derv2(ly,mz) + d2k2s*dklz*dkmy
        derv2(lz,mz) = derv2(lz,mz) + d2k2s*dklz*dkmz
      endif
!
!  End of loop over l
!
    enddo
!
!  End of loop over m
!
  enddo
#ifdef TRACE
  call trace_out('d2charge3fc')
#endif
!
  return
  end

  subroutine d2chargefc(i,j,nor,ix,iy,iz,jx,jy,jz,dei,dej, &
                       d1ixr,d1iyr,d1izr,d1jxr,d1jyr,d1jzr,d2i2r,d2ijr, &
                       d2j2r,d2self,dei0,dej0,lreal,nbosptr)
!
!  Calculates the contribution to the dynamical matrix due to charge derivatives from EEM/QEq. 
!  Non-phased fc_supercell version.
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!
!  10/14 Created from d2chargep
!   2/18 Trace added
!   5/18 Multiple qranges added
!   8/19 Trap for neemrptr being zero added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   2/21 Modifications added for GFNFF
!   6/21 Exact d2q terms rearranged to speed up
!   7/21 Timing added
!   7/21 New algorithm to avoid 4th power scaling added
!  10/21 lgfnff moved to control
!  12/22 Optional argument, nbosptr added
!   9/23 Trapping of shells added
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
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use reaxFFdata,     only : reaxFFmu
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)           :: i
  integer(i4), intent(in)           :: ix
  integer(i4), intent(in)           :: iy
  integer(i4), intent(in)           :: iz
  integer(i4), intent(in)           :: j
  integer(i4), intent(in)           :: jx
  integer(i4), intent(in)           :: jy
  integer(i4), intent(in)           :: jz
  integer(i4), intent(in)           :: nor
  integer(i4), intent(in), optional :: nbosptr(*)       ! Only required for reaxFF calls - points to reaxFF species from atoms
  logical,     intent(in)           :: lreal
  real(dp),    intent(in)           :: dei(*)
  real(dp),    intent(in)           :: dej(*)
  real(dp),    intent(in)           :: dei0
  real(dp),    intent(in)           :: dej0
  real(dp),    intent(in)           :: d1ixr
  real(dp),    intent(in)           :: d1iyr
  real(dp),    intent(in)           :: d1izr
  real(dp),    intent(in)           :: d1jxr
  real(dp),    intent(in)           :: d1jyr
  real(dp),    intent(in)           :: d1jzr
  real(dp),    intent(in)           :: d2i2r(*)
  real(dp),    intent(in)           :: d2ijr(*)
  real(dp),    intent(in)           :: d2j2r(*)
  real(dp),    intent(in)           :: d2self
!
!  Local variables     
!
  integer(i4)                       :: iv    
  integer(i4)                       :: k
  integer(i4)                       :: kx
  integer(i4)                       :: ky 
  integer(i4)                       :: kz  
  integer(i4)                       :: l      
  integer(i4)                       :: lx     
  integer(i4)                       :: ly    
  integer(i4)                       :: lz    
  integer(i4)                       :: m
  integer(i4)                       :: mnxx
  integer(i4)                       :: mnxy
  integer(i4)                       :: mnxz
  integer(i4)                       :: mnyx
  integer(i4)                       :: mnyy
  integer(i4)                       :: mnyz
  integer(i4)                       :: mnzx
  integer(i4)                       :: mnzy
  integer(i4)                       :: mnzz
  integer(i4)                       :: mx
  integer(i4)                       :: my
  integer(i4)                       :: mz
  integer(i4)                       :: n
  integer(i4)                       :: nqr
  integer(i4)                       :: nx
  integer(i4)                       :: nk
  logical                           :: lNonIJQDeriv
  real(dp)                          :: d1ix
  real(dp)                          :: d1iy
  real(dp)                          :: d1iz
  real(dp)                          :: d1jx
  real(dp)                          :: d1jy
  real(dp)                          :: d1jz
  real(dp)                          :: d2i2s
  real(dp)                          :: d2ijs
  real(dp)                          :: d2j2s
  real(dp)                          :: d2qk 
  real(dp)                          :: deisum
  real(dp)                          :: dejsum
  real(dp)                          :: dikx
  real(dp)                          :: diky
  real(dp)                          :: dikz
  real(dp)                          :: d2ikx
  real(dp)                          :: d2iky
  real(dp)                          :: d2ikz
  real(dp)                          :: dilx
  real(dp)                          :: dily
  real(dp)                          :: dilz 
  real(dp)                          :: djkx  
  real(dp)                          :: djky  
  real(dp)                          :: djkz   
  real(dp)                          :: d2jkx
  real(dp)                          :: d2jky
  real(dp)                          :: d2jkz
  real(dp)                          :: djlx   
  real(dp)                          :: djly
  real(dp)                          :: djlz
  real(dp)                          :: dkix
  real(dp)                          :: dkiy
  real(dp)                          :: dkiz 
  real(dp)                          :: dkjx 
  real(dp)                          :: dkjy 
  real(dp)                          :: dkjz 
  real(dp)                          :: ock 
  real(dp)                          :: qlk 
  real(dp)                          :: g_cpu_time
  real(dp)                          :: time1
  real(dp)                          :: time2
  real(dp), parameter               :: tsqrt2pi = 0.797884560802866_dp
  real(dp)                          :: zetah0
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargefc')
#endif
!
  time1 = g_cpu_time()
!
!  For ReaxFF check that nbosptr is present
!
  if (lreaxFF) then
    if (.not.present(nbosptr)) then
      call outerror('d2chargefc called for reaxFF without passing nbosptr',0_i4)
      call stopnow('d2chargefc')
    endif
  endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  d1ix = d1ixr
  d1iy = d1iyr
  d1iz = d1izr
  d1jx = d1jxr
  d1jy = d1jyr
  d1jz = d1jzr
  d2ijs = d2self
  d2i2s = 0.0_dp
  d2j2s = 0.0_dp
  do iv = 1,nor
    d2ijs = d2ijs + d2ijr(iv)
    d2i2s = d2i2s + d2i2r(iv)
    d2j2s = d2j2s + d2j2r(iv)
  enddo
  if (leem) then
    if (nat(i).gt.maxele) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    elseif (.not.lelementOK(nat(i)).or.nregionno(nrelf2a(i)).ne.1) then
      d2ijs = 0.0_dp
      d2i2s = 0.0_dp
      d1ix = 0.0_dp
      d1iy = 0.0_dp
      d1iz = 0.0_dp
    endif
    if (nat(j).gt.maxele) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp 
    elseif (.not.lelementOK(nat(j)).or.nregionno(nrelf2a(j)).ne.1) then
      d2ijs = 0.0_dp
      d2j2s = 0.0_dp
      d1jx = 0.0_dp
      d1jy = 0.0_dp
      d1jz = 0.0_dp 
    endif
  endif
  if (.not.loldd2q) then
!******************
!  New algorithm  *
!******************
    d2edqdq(j,i) = d2edqdq(j,i) + d2ijs
    d2edq2(i) = d2edq2(i) + d2i2s
    if (i.ne.j) then
      d2edqdq(i,j) = d2edqdq(i,j) + d2ijs
      d2edq2(j) = d2edq2(j) + d2j2s
    endif
  endif
!
!  The following flag causes the second derivatives of the charge distribution 
!  with respect to two atoms neither of which are the i/j of the energy term
!  to be calculated. At the moment this contribution is not used in any method
!  and so the flag is set to false for efficiency.
!
  lNonIJQDeriv = .false.
!
!  Second derivative of charge with respect to coordinate contributions
!
!  This term is zero for Hellman-Feynman case
!
  if (.not.leem.and..not.lgfnff.and..not.lreaxFF) then
    deisum = dei0
    dejsum = dej0
    do iv = 1,nor
      deisum = deisum + dei(iv)
      dejsum = dejsum + dej(iv)
    enddo
!
!  Terms involving Q derivatives for i/j and one other atom
!
    do m = 1,nqatoms(i)
      k = nqatomptr(m,i)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (i.gt.k) then
        derv2(kx,ix) = derv2(kx,ix) + deisum*d2qdxyz2(mnxx,i)
        derv2(ky,ix) = derv2(ky,ix) + deisum*d2qdxyz2(mnxy,i)
        derv2(kz,ix) = derv2(kz,ix) + deisum*d2qdxyz2(mnxz,i)
        derv2(kx,iy) = derv2(kx,iy) + deisum*d2qdxyz2(mnxy,i)
        derv2(ky,iy) = derv2(ky,iy) + deisum*d2qdxyz2(mnyy,i)
        derv2(kz,iy) = derv2(kz,iy) + deisum*d2qdxyz2(mnyz,i)
        derv2(kx,iz) = derv2(kx,iz) + deisum*d2qdxyz2(mnxz,i)
        derv2(ky,iz) = derv2(ky,iz) + deisum*d2qdxyz2(mnyz,i)
        derv2(kz,iz) = derv2(kz,iz) + deisum*d2qdxyz2(mnzz,i)
      elseif (i.lt.k) then
        derv2(ix,kx) = derv2(ix,kx) + deisum*d2qdxyz2(mnxx,i)
        derv2(iy,kx) = derv2(iy,kx) + deisum*d2qdxyz2(mnxy,i)
        derv2(iz,kx) = derv2(iz,kx) + deisum*d2qdxyz2(mnxz,i)
        derv2(ix,ky) = derv2(ix,ky) + deisum*d2qdxyz2(mnxy,i)
        derv2(iy,ky) = derv2(iy,ky) + deisum*d2qdxyz2(mnyy,i)
        derv2(iz,ky) = derv2(iz,ky) + deisum*d2qdxyz2(mnyz,i)
        derv2(ix,kz) = derv2(ix,kz) + deisum*d2qdxyz2(mnxz,i)
        derv2(iy,kz) = derv2(iy,kz) + deisum*d2qdxyz2(mnyz,i)
        derv2(iz,kz) = derv2(iz,kz) + deisum*d2qdxyz2(mnzz,i)
      endif
    enddo
!
    do m = 1,nqatoms(j)
      k = nqatomptr(m,j)
      kx = 3*(k - 1) + 1
      ky = kx + 1
      kz = ky + 1
      mx = 3*(m - 1) + 1
      my = mx + 1
      mz = my + 1
!
      mnxx = mx*(mx + 1)/2 
      mnyy = my*(my + 1)/2 
      mnzz = mz*(mz + 1)/2 
      mnxy = mnyy - 1
      mnxz = mnzz - 2
      mnyz = mnzz - 1
!
      if (j.gt.k) then
        derv2(kx,jx) = derv2(kx,jx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(ky,jx) = derv2(ky,jx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(kz,jx) = derv2(kz,jx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(kx,jy) = derv2(kx,jy) + dejsum*d2qdxyz2(mnxy,j)
        derv2(ky,jy) = derv2(ky,jy) + dejsum*d2qdxyz2(mnyy,j)
        derv2(kz,jy) = derv2(kz,jy) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kx,jz) = derv2(kx,jz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(ky,jz) = derv2(ky,jz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(kz,jz) = derv2(kz,jz) + dejsum*d2qdxyz2(mnzz,j)
      elseif (j.lt.k) then
        derv2(jx,kx) = derv2(jx,kx) + dejsum*d2qdxyz2(mnxx,j)
        derv2(jy,kx) = derv2(jy,kx) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jz,kx) = derv2(jz,kx) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jx,ky) = derv2(jx,ky) + dejsum*d2qdxyz2(mnxy,j)
        derv2(jy,ky) = derv2(jy,ky) + dejsum*d2qdxyz2(mnyy,j)
        derv2(jz,ky) = derv2(jz,ky) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jx,kz) = derv2(jx,kz) + dejsum*d2qdxyz2(mnxz,j)
        derv2(jy,kz) = derv2(jy,kz) + dejsum*d2qdxyz2(mnyz,j)
        derv2(jz,kz) = derv2(jz,kz) + dejsum*d2qdxyz2(mnzz,j)
      endif
    enddo
!
!  Terms involving Q derivatives for non-i/j case
!
    if (lNonIJQDeriv) then
      do m = 2,nqatoms(i)
        k = nqatomptr(m,i)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,i)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,i)
          derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,i)
          derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,i)
          derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,i)
          derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,i)
          derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,i)
          derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,i)
          derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,i)
          derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,i)
        enddo
      enddo
!
      do m = 2,nqatoms(j)
        k = nqatomptr(m,j)
        kx = 3*(k - 1) + 1
        ky = kx + 1
        kz = ky + 1
        mx = 3*(m - 1) + 1
        my = mx + 1
        mz = my + 1
        do n = 1,m-1
          l = nqatomptr(n,j)
          lx = 3*(l - 1) + 1
          ly = lx + 1
          lz = ly + 1
          nx = 3*(n - 1) + 1
!
          mnxx = mx*(mx - 1)/2 + nx
          mnxy = mnxx + 1
          mnxz = mnxx + 2
          mnyx = my*(my - 1)/2 + nx
          mnyy = mnyx + 1
          mnyz = mnyx + 2
          mnzx = mz*(mz - 1)/2 + nx
          mnzy = mnzx + 1
          mnzz = mnzx + 2
!
          derv2(lx,kx) = derv2(lx,kx) + deisum*d2qdxyz2(mnxx,j)
          derv2(ly,kx) = derv2(ly,kx) + deisum*d2qdxyz2(mnxy,j)
          derv2(lz,kx) = derv2(lz,kx) + deisum*d2qdxyz2(mnxz,j)
          derv2(lx,ky) = derv2(lx,ky) + deisum*d2qdxyz2(mnyx,j)
          derv2(ly,ky) = derv2(ly,ky) + deisum*d2qdxyz2(mnyy,j)
          derv2(lz,ky) = derv2(lz,ky) + deisum*d2qdxyz2(mnyz,j)
          derv2(lx,kz) = derv2(lx,kz) + deisum*d2qdxyz2(mnzx,j)
          derv2(ly,kz) = derv2(ly,kz) + deisum*d2qdxyz2(mnzy,j)
          derv2(lz,kz) = derv2(lz,kz) + deisum*d2qdxyz2(mnzz,j)
        enddo
      enddo
    endif
  endif
!
!  Loop over atoms to add ij-k contributions
!
  kx = - 2
  ky = - 1
  kz =   0
  do k = 1,numat
    kx = kx + 3
    ky = ky + 3
    kz = kz + 3
    dikx = dqdxyz(kx,i)
    diky = dqdxyz(ky,i)
    dikz = dqdxyz(kz,i)
    djkx = dqdxyz(kx,j)
    djky = dqdxyz(ky,j)
    djkz = dqdxyz(kz,j)
    if (i.ne.k) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : i-k
!
      if (k.le.i) then
        derv2(kx,ix) = derv2(kx,ix) - d1ix*dikx
        derv2(ky,ix) = derv2(ky,ix) - d1ix*diky
        derv2(kz,ix) = derv2(kz,ix) - d1ix*dikz
        derv2(kx,iy) = derv2(kx,iy) - d1iy*dikx
        derv2(ky,iy) = derv2(ky,iy) - d1iy*diky
        derv2(kz,iy) = derv2(kz,iy) - d1iy*dikz
        derv2(kx,iz) = derv2(kx,iz) - d1iz*dikx
        derv2(ky,iz) = derv2(ky,iz) - d1iz*diky
        derv2(kz,iz) = derv2(kz,iz) - d1iz*dikz
!
        derv2(kx,ix) = derv2(kx,ix) - d1jx*djkx
        derv2(ky,ix) = derv2(ky,ix) - d1jx*djky
        derv2(kz,ix) = derv2(kz,ix) - d1jx*djkz
        derv2(kx,iy) = derv2(kx,iy) - d1jy*djkx
        derv2(ky,iy) = derv2(ky,iy) - d1jy*djky
        derv2(kz,iy) = derv2(kz,iy) - d1jy*djkz
        derv2(kx,iz) = derv2(kx,iz) - d1jz*djkx
        derv2(ky,iz) = derv2(ky,iz) - d1jz*djky
        derv2(kz,iz) = derv2(kz,iz) - d1jz*djkz
      else
        derv2(ix,kx) = derv2(ix,kx) - d1ix*dikx
        derv2(iy,kx) = derv2(iy,kx) - d1iy*dikx
        derv2(iz,kx) = derv2(iz,kx) - d1iz*dikx
        derv2(ix,ky) = derv2(ix,ky) - d1ix*diky
        derv2(iy,ky) = derv2(iy,ky) - d1iy*diky
        derv2(iz,ky) = derv2(iz,ky) - d1iz*diky
        derv2(ix,kz) = derv2(ix,kz) - d1ix*dikz
        derv2(iy,kz) = derv2(iy,kz) - d1iy*dikz
        derv2(iz,kz) = derv2(iz,kz) - d1iz*dikz
!
        derv2(ix,kx) = derv2(ix,kx) - d1jx*djkx
        derv2(iy,kx) = derv2(iy,kx) - d1jy*djkx
        derv2(iz,kx) = derv2(iz,kx) - d1jz*djkx
        derv2(ix,ky) = derv2(ix,ky) - d1jx*djky
        derv2(iy,ky) = derv2(iy,ky) - d1jy*djky
        derv2(iz,ky) = derv2(iz,ky) - d1jz*djky
        derv2(ix,kz) = derv2(ix,kz) - d1jx*djkz
        derv2(iy,kz) = derv2(iy,kz) - d1jy*djkz
        derv2(iz,kz) = derv2(iz,kz) - d1jz*djkz
      endif
    endif
    if (j.ne.k) then
!
!  d2E/(d(alpha).dq) x dq/d(beta) : j-k
!
      if (k.le.j) then
        derv2(kx,jx) = derv2(kx,jx) + d1jx*djkx
        derv2(ky,jx) = derv2(ky,jx) + d1jx*djky
        derv2(kz,jx) = derv2(kz,jx) + d1jx*djkz
        derv2(kx,jy) = derv2(kx,jy) + d1jy*djkx
        derv2(ky,jy) = derv2(ky,jy) + d1jy*djky
        derv2(kz,jy) = derv2(kz,jy) + d1jy*djkz
        derv2(kx,jz) = derv2(kx,jz) + d1jz*djkx
        derv2(ky,jz) = derv2(ky,jz) + d1jz*djky
        derv2(kz,jz) = derv2(kz,jz) + d1jz*djkz
!
        derv2(kx,jx) = derv2(kx,jx) + d1ix*dikx
        derv2(ky,jx) = derv2(ky,jx) + d1ix*diky
        derv2(kz,jx) = derv2(kz,jx) + d1ix*dikz
        derv2(kx,jy) = derv2(kx,jy) + d1iy*dikx
        derv2(ky,jy) = derv2(ky,jy) + d1iy*diky
        derv2(kz,jy) = derv2(kz,jy) + d1iy*dikz
        derv2(kx,jz) = derv2(kx,jz) + d1iz*dikx
        derv2(ky,jz) = derv2(ky,jz) + d1iz*diky
        derv2(kz,jz) = derv2(kz,jz) + d1iz*dikz
      else
        derv2(jx,kx) = derv2(jx,kx) + d1jx*djkx
        derv2(jy,kx) = derv2(jy,kx) + d1jy*djkx
        derv2(jz,kx) = derv2(jz,kx) + d1jz*djkx
        derv2(jx,ky) = derv2(jx,ky) + d1jx*djky
        derv2(jy,ky) = derv2(jy,ky) + d1jy*djky
        derv2(jz,ky) = derv2(jz,ky) + d1jz*djky
        derv2(jx,kz) = derv2(jx,kz) + d1jx*djkz
        derv2(jy,kz) = derv2(jy,kz) + d1jy*djkz
        derv2(jz,kz) = derv2(jz,kz) + d1jz*djkz
!
        derv2(jx,kx) = derv2(jx,kx) + d1ix*dikx
        derv2(jy,kx) = derv2(jy,kx) + d1iy*dikx
        derv2(jz,kx) = derv2(jz,kx) + d1iz*dikx
        derv2(jx,ky) = derv2(jx,ky) + d1ix*diky
        derv2(jy,ky) = derv2(jy,ky) + d1iy*diky
        derv2(jz,ky) = derv2(jz,ky) + d1iz*diky
        derv2(jx,kz) = derv2(jx,kz) + d1ix*dikz
        derv2(jy,kz) = derv2(jy,kz) + d1iy*dikz
        derv2(jz,kz) = derv2(jz,kz) + d1iz*dikz
      endif
    endif
    if (lreal.and.i.ne.j) then
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
      nk = nat(k)
      ock = occuf(k)
      if (leem.and.nk.le.maxele) then
        if (lmultiqrange.and.neemrptr(k).ne.0) then
          nqr = nqrnow(neemrptr(k))
        else
          nqr = 1
        endif
        if (lqeq) then
          if (nk.ne.1) then
            d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
          else
            qlk = qf(k)
            zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
            d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp+2.0_dp*qlk/zetah0)
          endif
        else
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        endif
      elseif (lgfnff) then
        d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
      elseif (lreaxFF) then
        d2qk = 2.0_dp*ock*reaxFFmu(nbosptr(k))
      else
        d2qk = 0.0_dp
      endif
      dkix = dqdxyz(ix,k)
      dkiy = dqdxyz(iy,k)
      dkiz = dqdxyz(iz,k)
      dkjx = dqdxyz(jx,k)
      dkjy = dqdxyz(jy,k)
      dkjz = dqdxyz(jz,k)
      derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
      derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
      derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
      derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
      derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
      derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
      derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
      derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
      derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
    endif
    if (loldd2q) then
!******************
!  Old algorithm  *
!******************
!
!  Loop over atoms to add ij-kl correction
!
      lx = - 2
      ly = - 1
      lz =   0
!
!  Scale charge derivatives for k by d2ijs
!
      d2ikx = d2ijs*dikx
      d2iky = d2ijs*diky
      d2ikz = d2ijs*dikz
      d2jkx = d2ijs*djkx
      d2jky = d2ijs*djky
      d2jkz = d2ijs*djkz
!
      do l = 1,k-1
        lx = lx + 3
        ly = ly + 3
        lz = lz + 3
!
!  d2E/(dq.dq) x dq/d(alpha) x dq/d(beta) 
!
        dilx = dqdxyz(lx,i)
        dily = dqdxyz(ly,i)
        dilz = dqdxyz(lz,i)
        djlx = dqdxyz(lx,j)
        djly = dqdxyz(ly,j)
        djlz = dqdxyz(lz,j)
!
        derv2(lx,kx) = derv2(lx,kx) + dilx*d2jkx + djlx*d2ikx
        derv2(ly,kx) = derv2(ly,kx) + dilx*d2jky + djlx*d2iky
        derv2(lz,kx) = derv2(lz,kx) + dilx*d2jkz + djlx*d2ikz
        derv2(lx,ky) = derv2(lx,ky) + dily*d2jkx + djly*d2ikx
        derv2(ly,ky) = derv2(ly,ky) + dily*d2jky + djly*d2iky
        derv2(lz,ky) = derv2(lz,ky) + dily*d2jkz + djly*d2ikz
        derv2(lx,kz) = derv2(lx,kz) + dilz*d2jkx + djlz*d2ikx
        derv2(ly,kz) = derv2(ly,kz) + dilz*d2jky + djlz*d2iky
        derv2(lz,kz) = derv2(lz,kz) + dilz*d2jkz + djlz*d2ikz
!
!  d2Edi2/d2Edj2 - only non-zero for QEq/H
!
        if (abs(d2i2s).gt.1.0d-8) then
          derv2(lx,kx) = derv2(lx,kx) + d2i2s*dilx*dikx
          derv2(ly,kx) = derv2(ly,kx) + d2i2s*dilx*diky
          derv2(lz,kx) = derv2(lz,kx) + d2i2s*dilx*dikz
          derv2(lx,ky) = derv2(lx,ky) + d2i2s*dily*dikx
          derv2(ly,ky) = derv2(ly,ky) + d2i2s*dily*diky
          derv2(lz,ky) = derv2(lz,ky) + d2i2s*dily*dikz
          derv2(lx,kz) = derv2(lx,kz) + d2i2s*dilz*dikx
          derv2(ly,kz) = derv2(ly,kz) + d2i2s*dilz*diky
          derv2(lz,kz) = derv2(lz,kz) + d2i2s*dilz*dikz
        endif
        if (abs(d2j2s).gt.1.0d-8.and.i.ne.j) then
          derv2(lx,kx) = derv2(lx,kx) + d2j2s*djlx*djkx
          derv2(ly,kx) = derv2(ly,kx) + d2j2s*djlx*djky
          derv2(lz,kx) = derv2(lz,kx) + d2j2s*djlx*djkz
          derv2(lx,ky) = derv2(lx,ky) + d2j2s*djly*djkx
          derv2(ly,ky) = derv2(ly,ky) + d2j2s*djly*djky
          derv2(lz,ky) = derv2(lz,ky) + d2j2s*djly*djkz
          derv2(lx,kz) = derv2(lx,kz) + d2j2s*djlz*djkx
          derv2(ly,kz) = derv2(ly,kz) + d2j2s*djlz*djky
          derv2(lz,kz) = derv2(lz,kz) + d2j2s*djlz*djkz
        endif
!
!  End of loop over l
!
      enddo
    endif
!
!  End of loop over k
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargefc')
#endif
!
  return
  end
!
  subroutine d2chargeself(i,j,ix,iy,iz,jx,jy,jz,lopi,lopj,nbosptr)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from variable charge models for the self term
!  Should only be called from real space and once per atom pair!
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!
!   4/23 Created from d2charge
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
!  Julian Gale, CIC, Curtin University, April 2023
!
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use optimisation
  use reaxFFdata,     only : reaxFFmu
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)           :: i
  integer(i4), intent(in)           :: ix
  integer(i4), intent(in)           :: iy
  integer(i4), intent(in)           :: iz
  integer(i4), intent(in)           :: j
  integer(i4), intent(in)           :: jx
  integer(i4), intent(in)           :: jy
  integer(i4), intent(in)           :: jz
  integer(i4), intent(in), optional :: nbosptr(*)       ! Only required for reaxFF calls - points to reaxFF species from atoms
  logical,     intent(in)           :: lopi
  logical,     intent(in)           :: lopj
!
!  Local variables
!
  integer(i4)                       :: iix
  integer(i4)                       :: iiy
  integer(i4)                       :: iiz
  integer(i4)                       :: indi
  integer(i4)                       :: indj
  integer(i4)                       :: jjx
  integer(i4)                       :: jjy
  integer(i4)                       :: jjz
  integer(i4)                       :: k
  integer(i4)                       :: kl
  integer(i4)                       :: nqr
  integer(i4)                       :: nk
  real(dp)                          :: d2qk
  real(dp)                          :: dkix
  real(dp)                          :: dkiy
  real(dp)                          :: dkiz
  real(dp)                          :: dkjx
  real(dp)                          :: dkjy
  real(dp)                          :: dkjz
  real(dp)                          :: dk1
  real(dp)                          :: dk2
  real(dp)                          :: dk3
  real(dp)                          :: dk4
  real(dp)                          :: dk5
  real(dp)                          :: dk6
  real(dp)                          :: g_cpu_time
  real(dp)                          :: time1
  real(dp)                          :: time2
  real(dp)                          :: ock
  real(dp)                          :: qlk
  real(dp)                          :: zetah0
  real(dp), parameter               :: tsqrt2pi = 0.797884560802866_dp
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
!
  if (.not.lopi.and..not.lopj) return
#ifdef TRACE
  call trace_in('d2chargeself')
#endif
!
  time1 = g_cpu_time()
!
!  For ReaxFF check that nbosptr is present
!
  if (lreaxFF) then
    if (.not.present(nbosptr)) then
      call outerror('d2chargeself called for reaxFF without passing nbosptr',0_i4)
      call stopnow('d2chargeself')
    endif
  endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  if (lopi.and.lopj) then
    iix = ix
    iiy = iy
    iiz = iz
    jjx = jx
    jjy = jy
    jjz = jz
  elseif (lopi) then
    iix = ix
    iiy = iy
    iiz = iz
    jjx = ix
    jjy = iy
    jjz = iz
  elseif (lopj) then
    iix = jx
    iiy = jy
    iiz = jz
    jjx = jx
    jjy = jy
    jjz = jz
  endif
  indi = 3*(i - 1)
  indj = 3*(j - 1)
!
!  Loop over atoms to add ij-k contributions
!
  do k = 1,numat
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
    nk = nat(k)
    ock = occuf(k)
    if (leem.and.nk.le.maxele) then
      if (lmultiqrange.and.neemrptr(k).ne.0) then
        nqr = nqrnow(neemrptr(k))
      else
        nqr = 1
      endif
      if (lqeq) then
        if (nk.ne.1) then
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        else
          qlk = qf(k)
          zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
          d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp + 2.0_dp*qlk/zetah0)
        endif
      else
        d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
      endif
    elseif (lgfnff) then
      d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
    elseif (lreaxFF) then
      d2qk = 2.0_dp*ock*reaxFFmu(nbosptr(k))
    else
      d2qk = 0.0_dp
    endif
    dkix = dqdxyz(indi+1,k)
    dkiy = dqdxyz(indi+2,k)
    dkiz = dqdxyz(indi+3,k)
    dkjx = dqdxyz(indj+1,k)
    dkjy = dqdxyz(indj+2,k)
    dkjz = dqdxyz(indj+3,k)
    derv2(jjx,iix) = derv2(jjx,iix) + d2qk*dkjx*dkix
    derv2(jjy,iix) = derv2(jjy,iix) + d2qk*dkjx*dkiy
    derv2(jjz,iix) = derv2(jjz,iix) + d2qk*dkjx*dkiz
    derv2(jjx,iiy) = derv2(jjx,iiy) + d2qk*dkjy*dkix
    derv2(jjy,iiy) = derv2(jjy,iiy) + d2qk*dkjy*dkiy
    derv2(jjz,iiy) = derv2(jjz,iiy) + d2qk*dkjy*dkiz
    derv2(jjx,iiz) = derv2(jjx,iiz) + d2qk*dkjz*dkix
    derv2(jjy,iiz) = derv2(jjy,iiz) + d2qk*dkjz*dkiy
    derv2(jjz,iiz) = derv2(jjz,iiz) + d2qk*dkjz*dkiz
    if (lstr.and.i.eq.j) then
!
!  Mixed strain terms from self energy - only need to do once for each atom - hence check on i = j
!
      if (lopi) then
        do kl = 1,nstrains
          derv3(iix,kl) = derv3(iix,kl) + d2qk*dqds(kl,k)*dkix
          derv3(iiy,kl) = derv3(iiy,kl) + d2qk*dqds(kl,k)*dkiy
          derv3(iiz,kl) = derv3(iiz,kl) + d2qk*dqds(kl,k)*dkiz
        enddo
      endif
!
!  Strain-strain derivatives - only need to do if i=j=k
!
      if (i.eq.k) then
        if (ndim.eq.3) then
          dk1 = dqds(1,k)
          dk2 = dqds(2,k)
          dk3 = dqds(3,k)
          dk4 = dqds(4,k)
          dk5 = dqds(5,k)
          dk6 = dqds(6,k)
          sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
          sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
          sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
          sderv2(4,1) = sderv2(4,1) + d2qk*dk4*dk1
          sderv2(5,1) = sderv2(5,1) + d2qk*dk5*dk1
          sderv2(6,1) = sderv2(6,1) + d2qk*dk6*dk1
          sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
          sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
          sderv2(4,2) = sderv2(4,2) + d2qk*dk4*dk2
          sderv2(5,2) = sderv2(5,2) + d2qk*dk5*dk2
          sderv2(6,2) = sderv2(6,2) + d2qk*dk6*dk2
          sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
          sderv2(4,3) = sderv2(4,3) + d2qk*dk4*dk3
          sderv2(5,3) = sderv2(5,3) + d2qk*dk5*dk3
          sderv2(6,3) = sderv2(6,3) + d2qk*dk6*dk3
          sderv2(4,4) = sderv2(4,4) + d2qk*dk4*dk4
          sderv2(5,4) = sderv2(5,4) + d2qk*dk5*dk4
          sderv2(6,4) = sderv2(6,4) + d2qk*dk6*dk4
          sderv2(5,5) = sderv2(5,5) + d2qk*dk5*dk5
          sderv2(6,5) = sderv2(6,5) + d2qk*dk6*dk5
          sderv2(6,6) = sderv2(6,6) + d2qk*dk6*dk6
        elseif (ndim.eq.2) then
          dk1 = dqds(1,k)
          dk2 = dqds(2,k)
          dk3 = dqds(3,k)
          sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
          sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
          sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
          sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
          sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
          sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
        elseif (ndim.eq.1) then
          dk1 = dqds(1,k)
          sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
        endif
      endif
    endif
!
!  End of loop over k
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargeself')
#endif
!
  return
  end
!
  subroutine d2chargeselfp(i,j,ix,iy,iz,jx,jy,jz,nbosptr)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from variable charge models for the self term
!  Should only be called from real space and once per atom pair!
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!
!   5/23 Created from d2chargeself
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
!  Julian Gale, CIC, Curtin University, May 2023
!
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use optimisation
  use reaxFFdata,     only : reaxFFmu
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)           :: i
  integer(i4), intent(in)           :: ix
  integer(i4), intent(in)           :: iy
  integer(i4), intent(in)           :: iz
  integer(i4), intent(in)           :: j
  integer(i4), intent(in)           :: jx
  integer(i4), intent(in)           :: jy
  integer(i4), intent(in)           :: jz
  integer(i4), intent(in), optional :: nbosptr(*)       ! Only required for reaxFF calls - points to reaxFF species from atoms
!
!  Local variables
!
  integer(i4)                       :: indi
  integer(i4)                       :: indj
  integer(i4)                       :: k
  integer(i4)                       :: nqr
  integer(i4)                       :: nk
  real(dp)                          :: d2qk
  real(dp)                          :: dkix
  real(dp)                          :: dkiy
  real(dp)                          :: dkiz
  real(dp)                          :: dkjx
  real(dp)                          :: dkjy
  real(dp)                          :: dkjz
  real(dp)                          :: g_cpu_time
  real(dp)                          :: time1
  real(dp)                          :: time2
  real(dp)                          :: ock
  real(dp)                          :: qlk
  real(dp)                          :: zetah0
  real(dp), parameter               :: tsqrt2pi = 0.797884560802866_dp
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargeselfp')
#endif
!
  time1 = g_cpu_time()
!
!  For ReaxFF check that nbosptr is present
!
  if (lreaxFF) then
    if (.not.present(nbosptr)) then
      call outerror('d2chargeselfp called for reaxFF without passing nbosptr',0_i4)
      call stopnow('d2chargeselfp')
    endif
  endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  indi = 3*(i - 1)
  indj = 3*(j - 1)
!
!  Loop over atoms to add ij-k contributions
!
  do k = 1,numat
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
    nk = nat(k)
    ock = occuf(k)
    if (leem.and.nk.le.maxele) then
      if (lmultiqrange.and.neemrptr(k).ne.0) then
        nqr = nqrnow(neemrptr(k))
      else
        nqr = 1
      endif
      if (lqeq) then
        if (nk.ne.1) then
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        else
          qlk = qf(k)
          zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
          d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp + 2.0_dp*qlk/zetah0)
        endif
      else
        d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
      endif
    elseif (lgfnff) then
      d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
    elseif (lreaxFF) then
      d2qk = 2.0_dp*ock*reaxFFmu(nbosptr(k))
    else
      d2qk = 0.0_dp
    endif
    dkix = dqdxyz(indi+1,k)
    dkiy = dqdxyz(indi+2,k)
    dkiz = dqdxyz(indi+3,k)
    dkjx = dqdxyz(indj+1,k)
    dkjy = dqdxyz(indj+2,k)
    dkjz = dqdxyz(indj+3,k)
    if (jx.ge.ix) then
      derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
      derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
      derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
      derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
      derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
      derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
      derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
      derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
      derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
    else
      derv2(ix,jx) = derv2(ix,jx) + d2qk*dkix*dkjx
      derv2(iy,jx) = derv2(iy,jx) + d2qk*dkix*dkjy
      derv2(iz,jx) = derv2(iz,jx) + d2qk*dkix*dkjz
      derv2(ix,jy) = derv2(ix,jy) + d2qk*dkiy*dkjx
      derv2(iy,jy) = derv2(iy,jy) + d2qk*dkiy*dkjy
      derv2(iz,jy) = derv2(iz,jy) + d2qk*dkiy*dkjz
      derv2(ix,jz) = derv2(ix,jz) + d2qk*dkiz*dkjx
      derv2(iy,jz) = derv2(iy,jz) + d2qk*dkiz*dkjy
      derv2(iz,jz) = derv2(iz,jz) + d2qk*dkiz*dkjz
    endif
!
!  End of loop over k
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargeselfp')
#endif
!
  return
  end
!
  subroutine d2chargeselfd(iloc,i,j,ix,iy,iz,jx,jy,jz,nbosptr)
!
!  Calculates the contribution to the second derivative matrices
!  due to charge derivatives from variable charge models for the self term
!  Should only be called from real space and once per atom pair!
!  Parallel version.
!
!  lreal = .true.  => call from real space routine
!        = .false. => call from reciprocal space routine
!
!   7/23 Created from d2charged
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
!  Julian Gale, CIC, Curtin University, July 2023
!
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use optimisation
  use parallel
  use reaxFFdata,     only : reaxFFmu
  use symmetry
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)           :: i       ! Atom i
  integer(i4), intent(in)           :: iloc    ! Atom i number on local node
  integer(i4), intent(in)           :: ix
  integer(i4), intent(in)           :: iy
  integer(i4), intent(in)           :: iz
  integer(i4), intent(in)           :: j
  integer(i4), intent(in)           :: jx
  integer(i4), intent(in)           :: jy
  integer(i4), intent(in)           :: jz
  integer(i4), intent(in), optional :: nbosptr(*)       ! Only required for reaxFF calls - points to reaxFF species from atoms
!
!  Local variables
!
  integer(i4)                       :: indi
  integer(i4)                       :: indj
  integer(i4)                       :: ixf
  integer(i4)                       :: iyf
  integer(i4)                       :: izf
  integer(i4)                       :: k
  integer(i4)                       :: kl
  integer(i4)                       :: nqr
  integer(i4)                       :: nk
  real(dp)                          :: d2qk
  real(dp)                          :: dkix
  real(dp)                          :: dkiy
  real(dp)                          :: dkiz
  real(dp)                          :: dkjx
  real(dp)                          :: dkjy
  real(dp)                          :: dkjz
  real(dp)                          :: dk1
  real(dp)                          :: dk2
  real(dp)                          :: dk3
  real(dp)                          :: dk4
  real(dp)                          :: dk5
  real(dp)                          :: dk6
  real(dp)                          :: g_cpu_time
  real(dp)                          :: time1
  real(dp)                          :: time2
  real(dp)                          :: ock
  real(dp)                          :: qlk
  real(dp)                          :: zetah0
  real(dp), parameter               :: tsqrt2pi = 0.797884560802866_dp
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargeselfd')
#endif
!
  time1 = g_cpu_time()
! 
!  For ReaxFF check that nbosptr is present
!   
  if (lreaxFF) then
    if (.not.present(nbosptr)) then
      call outerror('d2chargeselfd called for reaxFF without passing nbosptr',0_i4)
      call stopnow('d2chargeselfd')
    endif
  endif
!
!  Set pointers to i-j second derivative block for eself derivatives
!
  indi = 3*(i - 1)
  indj = 3*(j - 1)
  ixf = indi + 1
  iyf = indi + 2
  izf = indi + 3
!
!  Loop over atoms to add ij-k contributions
!
  do k = 1,numat
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
    nk = nat(k)
    ock = occuf(k)
    if (leem.and.nk.le.maxele) then
      if (lmultiqrange.and.neemrptr(k).ne.0) then
        nqr = nqrnow(neemrptr(k))
      else
        nqr = 1
      endif
      if (lqeq) then
        if (nk.ne.1) then
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        else
          qlk = qf(k)
          zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
          d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp + 2.0_dp*qlk/zetah0)
        endif
      else
        d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
      endif
    elseif (lgfnff) then
      d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
    elseif (lreaxFF) then
      d2qk = 2.0_dp*ock*reaxFFmu(nbosptr(k))
    else
      d2qk = 0.0_dp
    endif
    dkix = dqdxyz(indi+1,k)
    dkiy = dqdxyz(indi+2,k)
    dkiz = dqdxyz(indi+3,k)
    dkjx = dqdxyz(indj+1,k)
    dkjy = dqdxyz(indj+2,k)
    dkjz = dqdxyz(indj+3,k)
    derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
    derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
    derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
    derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
    derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
    derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
    derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
    derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
    derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
    if (lstr.and.i.eq.j) then
!
!  Mixed strain terms from self energy - only need to do once for each atom - hence check on i = j
!
      do kl = 1,nstrains
        derv3(ix,kl) = derv3(ix,kl) + d2qk*dqds(kl,k)*dkix
        derv3(iy,kl) = derv3(iy,kl) + d2qk*dqds(kl,k)*dkiy
        derv3(iz,kl) = derv3(iz,kl) + d2qk*dqds(kl,k)*dkiz
      enddo
!
!  Strain-strain derivatives - only need to do if i=j=k
!
      if (i.eq.k) then
        if (ndim.eq.3) then
          dk1 = dqds(1,k)
          dk2 = dqds(2,k)
          dk3 = dqds(3,k)
          dk4 = dqds(4,k)
          dk5 = dqds(5,k)
          dk6 = dqds(6,k)
          sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
          sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
          sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
          sderv2(4,1) = sderv2(4,1) + d2qk*dk4*dk1
          sderv2(5,1) = sderv2(5,1) + d2qk*dk5*dk1
          sderv2(6,1) = sderv2(6,1) + d2qk*dk6*dk1
          sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
          sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
          sderv2(4,2) = sderv2(4,2) + d2qk*dk4*dk2
          sderv2(5,2) = sderv2(5,2) + d2qk*dk5*dk2
          sderv2(6,2) = sderv2(6,2) + d2qk*dk6*dk2
          sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
          sderv2(4,3) = sderv2(4,3) + d2qk*dk4*dk3
          sderv2(5,3) = sderv2(5,3) + d2qk*dk5*dk3
          sderv2(6,3) = sderv2(6,3) + d2qk*dk6*dk3
          sderv2(4,4) = sderv2(4,4) + d2qk*dk4*dk4
          sderv2(5,4) = sderv2(5,4) + d2qk*dk5*dk4
          sderv2(6,4) = sderv2(6,4) + d2qk*dk6*dk4
          sderv2(5,5) = sderv2(5,5) + d2qk*dk5*dk5
          sderv2(6,5) = sderv2(6,5) + d2qk*dk6*dk5
          sderv2(6,6) = sderv2(6,6) + d2qk*dk6*dk6
        elseif (ndim.eq.2) then
          dk1 = dqds(1,k)
          dk2 = dqds(2,k)
          dk3 = dqds(3,k)
          sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
          sderv2(2,1) = sderv2(2,1) + d2qk*dk2*dk1
          sderv2(3,1) = sderv2(3,1) + d2qk*dk3*dk1
          sderv2(2,2) = sderv2(2,2) + d2qk*dk2*dk2
          sderv2(3,2) = sderv2(3,2) + d2qk*dk3*dk2
          sderv2(3,3) = sderv2(3,3) + d2qk*dk3*dk3
        elseif (ndim.eq.1) then
          dk1 = dqds(1,k)
          sderv2(1,1) = sderv2(1,1) + d2qk*dk1*dk1
        endif
      endif
    endif
!
!  End of loop over k
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargeselfd')
#endif
!
  return
  end
!
  subroutine d2chargeselfpd(iloc,i,j,ix,iy,iz,jx,jy,jz,xkv,ykv,zkv,nbosptr)
!
!  Calculates the contribution to the dynamical matrix due to charge derivatives from EEM/QEq. 
!  This version is for calls from distributed memory routines.
!
!  At present this is gamma point only.
!
!   7/23 Created from d2chargeselfp and d2charged
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
!  Julian Gale, CIC, Curtin University, July 2023
!
  use control
  use current
  use derivatives
  use eemdata
  use element
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_gam
  use reaxFFdata,     only : reaxFFmu
  use times,          only : td2charge
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)           :: i       ! Atom i
  integer(i4), intent(in)           :: iloc    ! Atom i number on local node
  integer(i4), intent(in)           :: ix
  integer(i4), intent(in)           :: iy
  integer(i4), intent(in)           :: iz
  integer(i4), intent(in)           :: j
  integer(i4), intent(in)           :: jx
  integer(i4), intent(in)           :: jy
  integer(i4), intent(in)           :: jz
  integer(i4), intent(in), optional :: nbosptr(*)       ! Only required for reaxFF calls - points to reaxFF species from atoms
  real(dp),    intent(in)           :: xkv
  real(dp),    intent(in)           :: ykv
  real(dp),    intent(in)           :: zkv
!
!  Local variables     
!
  integer(i4)                       :: ixf
  integer(i4)                       :: iyf
  integer(i4)                       :: izf
  integer(i4)                       :: k
  integer(i4)                       :: nqr
  integer(i4)                       :: nk
  real(dp)                          :: d2qk 
  real(dp)                          :: dkix
  real(dp)                          :: dkiy
  real(dp)                          :: dkiz 
  real(dp)                          :: dkjx 
  real(dp)                          :: dkjy 
  real(dp)                          :: dkjz 
  real(dp)                          :: ock 
  real(dp)                          :: g_cpu_time
  real(dp)                          :: qlk 
  real(dp)                          :: time1
  real(dp)                          :: time2
  real(dp), parameter               :: tsqrt2pi = 0.797884560802866_dp
  real(dp)                          :: zetah0
!
!  Check that charge derivatives are in use
!
  if (.not.lDoQDeriv2) return
#ifdef TRACE
  call trace_in('d2chargeselfpd')
#endif
!
  time1 = g_cpu_time()
! 
!  For ReaxFF check that nbosptr is present
!   
  if (lreaxFF) then
    if (.not.present(nbosptr)) then
      call outerror('d2chargeselfpd called for reaxFF without passing nbosptr',0_i4)
      call stopnow('d2chargeselfpd')
    endif
  endif
!
  ixf = 3*(i-1) + 1
  iyf = ixf + 1
  izf = ixf + 2
!
!  Loop over atoms to add ij-k contributions
!
  do k = 1,numat
    if (i.ne.j) then
!**************************************************
!  Contributions from derivatives of self energy  *
!**************************************************
      nk = nat(k)
      ock = occuf(k)
      if (leem.and.nk.le.maxele) then
        if (lmultiqrange.and.neemrptr(k).ne.0) then
          nqr = nqrnow(neemrptr(k))
        else
          nqr = 1
        endif
        if (lqeq) then
          if (nk.ne.1) then
            d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
          else
            qlk = qf(k)
            zetah0 = 0.529177_dp*0.75_dp/radrange(nqr,1,neemtype)
            d2qk = 2.0_dp*ock*murange(nqr,1,neemtype)*(1.0_dp+2.0_dp*qlk/zetah0)
          endif
        else
          d2qk = 2.0_dp*ock*murange(nqr,nk,neemtype)
        endif
      elseif (lgfnff) then
        d2qk = ock*(tsqrt2pi/sqrt(gfnff_eeq_alp(k)) + gfnff_eeq_gam(k))*angstoev
      elseif (lreaxFF) then
        d2qk = 2.0_dp*ock*reaxFFmu(nbosptr(k))
      else
        d2qk = 0.0_dp
      endif
      dkix = dqdxyz(ixf,k)
      dkiy = dqdxyz(iyf,k)
      dkiz = dqdxyz(izf,k)
      dkjx = dqdxyz(jx,k)
      dkjy = dqdxyz(jy,k)
      dkjz = dqdxyz(jz,k)
      derv2(jx,ix) = derv2(jx,ix) + d2qk*dkjx*dkix
      derv2(jy,ix) = derv2(jy,ix) + d2qk*dkjx*dkiy
      derv2(jz,ix) = derv2(jz,ix) + d2qk*dkjx*dkiz
      derv2(jx,iy) = derv2(jx,iy) + d2qk*dkjy*dkix
      derv2(jy,iy) = derv2(jy,iy) + d2qk*dkjy*dkiy
      derv2(jz,iy) = derv2(jz,iy) + d2qk*dkjy*dkiz
      derv2(jx,iz) = derv2(jx,iz) + d2qk*dkjz*dkix
      derv2(jy,iz) = derv2(jy,iz) + d2qk*dkjz*dkiy
      derv2(jz,iz) = derv2(jz,iz) + d2qk*dkjz*dkiz
    endif
!
!  End of loop over k
!
  enddo
!
!  Timing
!
  time2 = g_cpu_time()
  td2charge = td2charge + time2 - time1
#ifdef TRACE
  call trace_out('d2chargeselfpd')
#endif
!
  return
  end

end module m_qderiv
