!*******************************************************************************************************
!  Routines that add the coordination number derivatives from GFNFF to the appropriate arrays in GULP  *
!*******************************************************************************************************
  subroutine gfnff_drv_dcn(i,dEdlogcn,logcn,dlogcndcn)
!
!  Computes the first derivatives of the log coordination number for the neighbour list in GFNFF
!
!  On entry : 
!
!  dEdlogcn        = derivatives of energy w.r.t. log of the coordination number
!  logcn           = log coordination numbers
!
!  On exit :
!
!  dlogcn          = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns         = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!  10/20 Created
!   2/21 Use of pointer to non-zero terms added
!   2/21 cut no longer passed in as not needed
!   2/21 Derivatives of coordination number now precomputed
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
!  Julian Gale, CIC, Curtin University, February 2021
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: dEdlogcn
  real(dp),    intent(in)                          :: logcn(*)
  real(dp),    intent(in)                          :: dlogcndcn(*)
!
!  Local variables
!
  integer(i4)                                      :: j
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  real(dp)                                         :: d1trm
  real(dp)                                         :: dcnidrij
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dr2ds(6)
  real(dp)                                         :: d2r2dx2(3,3)
  real(dp)                                         :: d2r2ds2(6,6)
  real(dp)                                         :: d2r2dsdx(6,3)
  real(dp)                                         :: dxyz(3)
#ifdef TRACE
  call trace_in('gfnff_drv_dcn')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
!
!  Compute derivatives of total coordination number for i
!
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    j = nbrno(ni,i)
    dcnidrij = d1cndr_cn(n,i)
    d1trm = dlogcnidcni*dcnidrij*dEdlogcn
    dxyz(1) = d1trm*xnbr(ni,i)
    dxyz(2) = d1trm*ynbr(ni,i)
    dxyz(3) = d1trm*znbr(ni,i)
!
    xdrv(i) = xdrv(i) - dxyz(1)
    ydrv(i) = ydrv(i) - dxyz(2)
    zdrv(i) = zdrv(i) - dxyz(3)
    xdrv(j) = xdrv(j) + dxyz(1)
    ydrv(j) = ydrv(j) + dxyz(2)
    zdrv(j) = zdrv(j) + dxyz(3)
!
    if (lstr) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
      do kl = 1,nstrains
        ks = nstrptr(kl)
        rstrd(kl) = rstrd(kl) + d1trm*dr2ds(ks)
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv_dcn')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn(i,j,xij,yij,zij,dEdlogcni,dEdlogcnj,d2Edlogcni2,d2Edlogcnj2,d2Edlogcnidlogcnj, &
                            d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn,d2logcndcn2,lgrad2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!
!  On entry : 
!
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  dEdlogcnj         = derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcnj2       = second derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcnidlogcnj = second derivatives of energy w.r.t. log of the coordination number for i and j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!  11/20 Created from gfnff_drv_dcn with second derivatives added
!   2/21 Use of pointer to non-zero terms added
!   2/21 cut no longer passed in as not needed
!   2/21 Derivatives of coordination number now precomputed
!   2/21 Tolerance introduced to reduce second derivative cost
!   4/21 Use of drop tolerance extended to further improve speed
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
!  Julian Gale, CIC, Curtin University, March 2021
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: dEdlogcnj
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: d2Edlogcnidlogcnj
  real(dp),    intent(in)                          :: d2Edlogcnj2
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: d3trm1
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dcnjdrjl
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2logcnjdcnj2
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2ilds(6)
  real(dp)                                         :: d2r2ildx2(6)
  real(dp)                                         :: d2r2ilds2(6,6)
  real(dp)                                         :: d2r2ildsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
  real(dp)                                         :: dr2jlds(6)
  real(dp)                                         :: d2r2jldx2(6)
  real(dp)                                         :: d2r2jlds2(6,6)
  real(dp)                                         :: d2r2jldsdx(6,3)
  real(dp)                                         :: dxyz(3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  if (lgrad2) then
    d2logcnidcni2 = d2logcndcn2(i)
    d2logcnjdcnj2 = d2logcndcn2(j)
  endif
  if (lgrad2.and.lstr) then
    call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2,d2r2ijdsdx,d2r2ijds2,.false.)
  endif
!
!  Set up for second derivatives
!
  if (lgrad2) then
    indi = 3*(i-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
!
    indj = 3*(j-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
  endif
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
    dxyz(1) = d1trm*xnbr(ni,i)
    dxyz(2) = d1trm*ynbr(ni,i)
    dxyz(3) = d1trm*znbr(ni,i)
!
    xdrv(i) = xdrv(i) - dxyz(1)
    ydrv(i) = ydrv(i) - dxyz(2)
    zdrv(i) = zdrv(i) - dxyz(3)
    xdrv(k) = xdrv(k) + dxyz(1)
    ydrv(k) = ydrv(k) + dxyz(2)
    zdrv(k) = zdrv(k) + dxyz(3)
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                        d2r2ikdsdx,d2r2ikds2,lgrad2)
    endif
    if (lstr) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        rstrd(kl) = rstrd(kl) + d1trm*dr2ikds(ks)
      enddo
    endif
    if (lgrad2) then
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
!
      d2cnidrik2 = d2cndr_cn(n,i)
      d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
             (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!
      if (abs(d2trm).gt.gfnff_cnc6tol) then
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
        if (k.ge.i) then
          derv2(kx,ix) = derv2(kx,ix) - d2trm*d2r2ikdx2(1)
          derv2(ky,ix) = derv2(ky,ix) - d2trm*d2r2ikdx2(6)
          derv2(kz,ix) = derv2(kz,ix) - d2trm*d2r2ikdx2(5)
          derv2(kx,iy) = derv2(kx,iy) - d2trm*d2r2ikdx2(6)
          derv2(ky,iy) = derv2(ky,iy) - d2trm*d2r2ikdx2(2)
          derv2(kz,iy) = derv2(kz,iy) - d2trm*d2r2ikdx2(4)
          derv2(kx,iz) = derv2(kx,iz) - d2trm*d2r2ikdx2(5)
          derv2(ky,iz) = derv2(ky,iz) - d2trm*d2r2ikdx2(4)
          derv2(kz,iz) = derv2(kz,iz) - d2trm*d2r2ikdx2(3)
          derv2(kx,ix) = derv2(kx,ix) - d1trm
          derv2(ky,iy) = derv2(ky,iy) - d1trm
          derv2(kz,iz) = derv2(kz,iz) - d1trm
        else
          derv2(ix,kx) = derv2(ix,kx) - d2trm*d2r2ikdx2(1)
          derv2(iy,kx) = derv2(iy,kx) - d2trm*d2r2ikdx2(6)
          derv2(iz,kx) = derv2(iz,kx) - d2trm*d2r2ikdx2(5)
          derv2(ix,ky) = derv2(ix,ky) - d2trm*d2r2ikdx2(6)
          derv2(iy,ky) = derv2(iy,ky) - d2trm*d2r2ikdx2(2)
          derv2(iz,ky) = derv2(iz,ky) - d2trm*d2r2ikdx2(4)
          derv2(ix,kz) = derv2(ix,kz) - d2trm*d2r2ikdx2(5)
          derv2(iy,kz) = derv2(iy,kz) - d2trm*d2r2ikdx2(4)
          derv2(iz,kz) = derv2(iz,kz) - d2trm*d2r2ikdx2(3)
          derv2(ix,kx) = derv2(ix,kx) - d1trm
          derv2(iy,ky) = derv2(iy,ky) - d1trm
          derv2(iz,kz) = derv2(iz,kz) - d1trm
        endif
!
        if (lstr) then
!
!  Mixed derivatives
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ikdsdx(ks,1) - xnbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ikdsdx(ks,2) - ynbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ikdsdx(ks,3) - znbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2ikdsdx(ks,1) + xnbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2ikdsdx(ks,2) + ynbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2ikdsdx(ks,3) + znbr(ni,i)*d2trm*dr2ikds(ks)
          enddo
!
!  Strain-strain
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2ikds(kt)*dr2ikds(ks) + d1trm*d2r2ikds2(kt,ks)
            enddo
          enddo
        endif
      endif
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
      d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
      if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  I-K
!
        if (k.ge.i) then
          derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xij
          derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xij
          derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xij
          derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*yij
          derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*yij
          derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*yij
          derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*zij
          derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*zij
          derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*zij
        else
          derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xij
          derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*yij
          derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*zij
          derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xij
          derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*yij
          derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*zij
          derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xij
          derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*yij
          derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*zij
        endif
!
!  J-K
!
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) + d3trm*xnbr(ni,i)*xij
          derv2(ky,jx) = derv2(ky,jx) + d3trm*ynbr(ni,i)*xij
          derv2(kz,jx) = derv2(kz,jx) + d3trm*znbr(ni,i)*xij
          derv2(kx,jy) = derv2(kx,jy) + d3trm*xnbr(ni,i)*yij
          derv2(ky,jy) = derv2(ky,jy) + d3trm*ynbr(ni,i)*yij
          derv2(kz,jy) = derv2(kz,jy) + d3trm*znbr(ni,i)*yij
          derv2(kx,jz) = derv2(kx,jz) + d3trm*xnbr(ni,i)*zij
          derv2(ky,jz) = derv2(ky,jz) + d3trm*ynbr(ni,i)*zij
          derv2(kz,jz) = derv2(kz,jz) + d3trm*znbr(ni,i)*zij
        else
          derv2(jx,kx) = derv2(jx,kx) + d3trm*xnbr(ni,i)*xij
          derv2(jy,kx) = derv2(jy,kx) + d3trm*xnbr(ni,i)*yij
          derv2(jz,kx) = derv2(jz,kx) + d3trm*xnbr(ni,i)*zij
          derv2(jx,ky) = derv2(jx,ky) + d3trm*ynbr(ni,i)*xij
          derv2(jy,ky) = derv2(jy,ky) + d3trm*ynbr(ni,i)*yij
          derv2(jz,ky) = derv2(jz,ky) + d3trm*ynbr(ni,i)*zij
          derv2(jx,kz) = derv2(jx,kz) + d3trm*znbr(ni,i)*xij
          derv2(jy,kz) = derv2(jy,kz) + d3trm*znbr(ni,i)*yij
          derv2(jz,kz) = derv2(jz,kz) + d3trm*znbr(ni,i)*zij
        endif
!
!  I-J
!
        if (j.ge.i) then
          derv2(jx,ix) = derv2(jx,ix) - d3trm*xnbr(ni,i)*xij
          derv2(jy,ix) = derv2(jy,ix) - d3trm*xnbr(ni,i)*yij
          derv2(jz,ix) = derv2(jz,ix) - d3trm*xnbr(ni,i)*zij
          derv2(jx,iy) = derv2(jx,iy) - d3trm*ynbr(ni,i)*xij
          derv2(jy,iy) = derv2(jy,iy) - d3trm*ynbr(ni,i)*yij
          derv2(jz,iy) = derv2(jz,iy) - d3trm*ynbr(ni,i)*zij
          derv2(jx,iz) = derv2(jx,iz) - d3trm*znbr(ni,i)*xij
          derv2(jy,iz) = derv2(jy,iz) - d3trm*znbr(ni,i)*yij
          derv2(jz,iz) = derv2(jz,iz) - d3trm*znbr(ni,i)*zij
        else
          derv2(ix,jx) = derv2(ix,jx) - d3trm*xnbr(ni,i)*xij
          derv2(iy,jx) = derv2(iy,jx) - d3trm*ynbr(ni,i)*xij
          derv2(iz,jx) = derv2(iz,jx) - d3trm*znbr(ni,i)*xij
          derv2(ix,jy) = derv2(ix,jy) - d3trm*xnbr(ni,i)*yij
          derv2(iy,jy) = derv2(iy,jy) - d3trm*ynbr(ni,i)*yij
          derv2(iz,jy) = derv2(iz,jy) - d3trm*znbr(ni,i)*yij
          derv2(ix,jz) = derv2(ix,jz) - d3trm*xnbr(ni,i)*zij
          derv2(iy,jz) = derv2(iy,jz) - d3trm*ynbr(ni,i)*zij
          derv2(iz,jz) = derv2(iz,jz) - d3trm*znbr(ni,i)*zij
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ijds(ks)
          enddo
!
         do kl = 1,nstrains
           ks = nstrptr(kl)
           derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2ikds(ks)
           derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2ikds(ks)
           derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2ikds(ks)
           derv3(jx,kl) = derv3(jx,kl) + xij*d3trm*dr2ikds(ks)
           derv3(jy,kl) = derv3(jy,kl) + yij*d3trm*dr2ikds(ks)
           derv3(jz,kl) = derv3(jz,kl) + zij*d3trm*dr2ikds(ks)
         enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2ikds(ks)) 
            enddo
          enddo
        endif
      endif
!
      d3trm1 = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
      do n2 = 1,n-1
        dcnidril = d1cndr_cn(n2,i)
        d3trm = d3trm1*dcnidril
!
        if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
        ni2 = nbrno_cn(n2,i)
        l = nbrno(ni2,i)
        indl = 3*(l-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
!
        if (lstr) then
          call real1strterm(ndim,xnbr(ni2,i),ynbr(ni2,i),znbr(ni2,i),0.0_dp,0.0_dp,0.0_dp,dr2ilds,d2r2ildx2, &
                            d2r2ildsdx,d2r2ilds2,.false.)
        endif
!
!  I-K
!
        if (k.ge.i) then
          derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        else
          derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  I-L
!
        if (l.ge.i) then
          derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ly,ix) = derv2(ly,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(lz,ix) = derv2(lz,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(lx,iy) = derv2(lx,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(lz,iy) = derv2(lz,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(lx,iz) = derv2(lx,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ly,iz) = derv2(ly,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        else
          derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(iy,lx) = derv2(iy,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(iz,lx) = derv2(iz,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ix,ly) = derv2(ix,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(iz,ly) = derv2(iz,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(ix,lz) = derv2(ix,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(iy,lz) = derv2(iy,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  K-L
!
        if (l.ge.k) then
          derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
        else
          derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ilds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(ni2,i)*d3trm*dr2ikds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ilds(ks) + dr2ilds(kt)*dr2ikds(ks))
            enddo
          enddo
        endif
      enddo
!
      d3trm1 = d2Edlogcnidlogcnj*dlogcnidcni*dlogcnjdcnj*dcnidrik
!------------------------------------------------------------------------------------------
!  Loop over neighbours of j for second derivatives for 2 different coordination numbers  |
!------------------------------------------------------------------------------------------
      do n2 = 1,nnbr_cn(j)
        dcnjdrjl = d1cndr_cn(n2,j)
        d3trm = d3trm1*dcnjdrjl
!
        if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
        nj = nbrno_cn(n2,j)
        l = nbrno(nj,j)
        indl = 3*(l-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
!
        if (lstr) then
          call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jlds,d2r2jldx2, &
                            d2r2jldsdx,d2r2jlds2,.false.)
        endif
!
!  I-J
!
        if (j.ge.i) then
          derv2(jx,ix) = derv2(jx,ix) + d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(jy,ix) = derv2(jy,ix) + d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(jz,ix) = derv2(jz,ix) + d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(jx,iy) = derv2(jx,iy) + d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(jy,iy) = derv2(jy,iy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(jz,iy) = derv2(jz,iy) + d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(jx,iz) = derv2(jx,iz) + d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(jy,iz) = derv2(jy,iz) + d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(jz,iz) = derv2(jz,iz) + d3trm*znbr(nj,j)*znbr(ni,i)
        else
          derv2(ix,jx) = derv2(ix,jx) + d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(iy,jx) = derv2(iy,jx) + d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(iz,jx) = derv2(iz,jx) + d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(ix,jy) = derv2(ix,jy) + d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(iy,jy) = derv2(iy,jy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(iz,jy) = derv2(iz,jy) + d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(ix,jz) = derv2(ix,jz) + d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(iy,jz) = derv2(iy,jz) + d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(iz,jz) = derv2(iz,jz) + d3trm*znbr(nj,j)*znbr(ni,i)
        endif
!
!  I-L
!
        if (l.ge.i) then
          derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(ly,ix) = derv2(ly,ix) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(lz,ix) = derv2(lz,ix) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(lx,iy) = derv2(lx,iy) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(lz,iy) = derv2(lz,iy) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(lx,iz) = derv2(lx,iz) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(ly,iz) = derv2(ly,iz) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(nj,j)*znbr(ni,i)
        else
          derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(iy,lx) = derv2(iy,lx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(iz,lx) = derv2(iz,lx) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(ix,ly) = derv2(ix,ly) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(iz,ly) = derv2(iz,ly) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(ix,lz) = derv2(ix,lz) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(iy,lz) = derv2(iy,lz) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(nj,j)*znbr(ni,i)
        endif
!
!  J-K
!
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(ky,jx) = derv2(ky,jx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(kz,jx) = derv2(kz,jx) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(kx,jy) = derv2(kx,jy) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(ky,jy) = derv2(ky,jy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(kz,jy) = derv2(kz,jy) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(kx,jz) = derv2(kx,jz) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(ky,jz) = derv2(ky,jz) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(kz,jz) = derv2(kz,jz) - d3trm*znbr(nj,j)*znbr(ni,i)
        else
          derv2(jx,kx) = derv2(jx,kx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(jy,kx) = derv2(jy,kx) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(jz,kx) = derv2(jz,kx) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(jx,ky) = derv2(jx,ky) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(jy,ky) = derv2(jy,ky) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(jz,ky) = derv2(jz,ky) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(jx,kz) = derv2(jx,kz) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(jy,kz) = derv2(jy,kz) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(jz,kz) = derv2(jz,kz) - d3trm*znbr(nj,j)*znbr(ni,i)
        endif
!
!  K-L
!     
        if (l.ge.k) then
          derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
          derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(nj,j)
          derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(nj,j)
          derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(nj,j)
          derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(nj,j)
          derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(nj,j)
          derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(nj,j)
          derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(nj,j)
          derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(nj,j)
        else
          derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
          derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(nj,j)
          derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(nj,j)
          derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(nj,j)
          derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(nj,j)
          derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(nj,j)
          derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(nj,j)
          derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(nj,j)
          derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(nj,j)
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2jlds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(nj,j)*d3trm*dr2ikds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2ikds(ks))
            enddo
          enddo
        endif
      enddo
    endif
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    dcnjdrjk = d1cndr_cn(n,j)
    d1trm = dlogcnjdcnj*dcnjdrjk*dEdlogcnj
    dxyz(1) = d1trm*xnbr(nj,j)
    dxyz(2) = d1trm*ynbr(nj,j)
    dxyz(3) = d1trm*znbr(nj,j)
!
    xdrv(j) = xdrv(j) - dxyz(1)
    ydrv(j) = ydrv(j) - dxyz(2)
    zdrv(j) = zdrv(j) - dxyz(3)
    xdrv(k) = xdrv(k) + dxyz(1)
    ydrv(k) = ydrv(k) + dxyz(2)
    zdrv(k) = zdrv(k) + dxyz(3)
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jkds,d2r2jkdx2, &
                        d2r2jkdsdx,d2r2jkds2,lgrad2)
    endif
    if (lstr) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        rstrd(kl) = rstrd(kl) + d1trm*dr2jkds(ks)
      enddo
    endif
    if (lgrad2) then
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
!
      d2cnjdrjk2 = d2cndr_cn(n,j)
      d2trm = dEdlogcnj*dlogcnjdcnj*d2cnjdrjk2 + &
             (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjk
!
      if (abs(d2trm).gt.gfnff_cnc6tol) then
!----------------------------------------
!  Second derivatives of CN w.r.t. rjk  |
!----------------------------------------
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) - d2trm*d2r2jkdx2(1)
          derv2(ky,jx) = derv2(ky,jx) - d2trm*d2r2jkdx2(6)
          derv2(kz,jx) = derv2(kz,jx) - d2trm*d2r2jkdx2(5)
          derv2(kx,jy) = derv2(kx,jy) - d2trm*d2r2jkdx2(6)
          derv2(ky,jy) = derv2(ky,jy) - d2trm*d2r2jkdx2(2)
          derv2(kz,jy) = derv2(kz,jy) - d2trm*d2r2jkdx2(4)
          derv2(kx,jz) = derv2(kx,jz) - d2trm*d2r2jkdx2(5)
          derv2(ky,jz) = derv2(ky,jz) - d2trm*d2r2jkdx2(4)
          derv2(kz,jz) = derv2(kz,jz) - d2trm*d2r2jkdx2(3)
          derv2(kx,jx) = derv2(kx,jx) - d1trm
          derv2(ky,jy) = derv2(ky,jy) - d1trm
          derv2(kz,jz) = derv2(kz,jz) - d1trm
        else
          derv2(jx,kx) = derv2(jx,kx) - d2trm*d2r2jkdx2(1)
          derv2(jy,kx) = derv2(jy,kx) - d2trm*d2r2jkdx2(6)
          derv2(jz,kx) = derv2(jz,kx) - d2trm*d2r2jkdx2(5)
          derv2(jx,ky) = derv2(jx,ky) - d2trm*d2r2jkdx2(6)
          derv2(jy,ky) = derv2(jy,ky) - d2trm*d2r2jkdx2(2)
          derv2(jz,ky) = derv2(jz,ky) - d2trm*d2r2jkdx2(4)
          derv2(jx,kz) = derv2(jx,kz) - d2trm*d2r2jkdx2(5)
          derv2(jy,kz) = derv2(jy,kz) - d2trm*d2r2jkdx2(4)
          derv2(jz,kz) = derv2(jz,kz) - d2trm*d2r2jkdx2(3)
          derv2(jx,kx) = derv2(jx,kx) - d1trm
          derv2(jy,ky) = derv2(jy,ky) - d1trm
          derv2(jz,kz) = derv2(jz,kz) - d1trm
        endif
!
        if (lstr) then
!
!  Mixed derivatives
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - d1trm*d2r2jkdsdx(ks,1) - xnbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(jy,kl) = derv3(jy,kl) - d1trm*d2r2jkdsdx(ks,2) - ynbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(jz,kl) = derv3(jz,kl) - d1trm*d2r2jkdsdx(ks,3) - znbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2jkdsdx(ks,1) + xnbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2jkdsdx(ks,2) + ynbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2jkdsdx(ks,3) + znbr(nj,j)*d2trm*dr2jkds(ks)
          enddo
!
!  Strain-strain
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2jkds(kt)*dr2jkds(ks) + d1trm*d2r2jkds2(kt,ks)
            enddo
          enddo
        endif
      endif
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
      d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
      if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  J-K
!
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) + d3trm*xnbr(nj,j)*xij
          derv2(ky,jx) = derv2(ky,jx) + d3trm*ynbr(nj,j)*xij
          derv2(kz,jx) = derv2(kz,jx) + d3trm*znbr(nj,j)*xij
          derv2(kx,jy) = derv2(kx,jy) + d3trm*xnbr(nj,j)*yij
          derv2(ky,jy) = derv2(ky,jy) + d3trm*ynbr(nj,j)*yij
          derv2(kz,jy) = derv2(kz,jy) + d3trm*znbr(nj,j)*yij
          derv2(kx,jz) = derv2(kx,jz) + d3trm*xnbr(nj,j)*zij
          derv2(ky,jz) = derv2(ky,jz) + d3trm*ynbr(nj,j)*zij
          derv2(kz,jz) = derv2(kz,jz) + d3trm*znbr(nj,j)*zij
        else
          derv2(jx,kx) = derv2(jx,kx) + d3trm*xnbr(nj,j)*xij
          derv2(jy,kx) = derv2(jy,kx) + d3trm*xnbr(nj,j)*yij
          derv2(jz,kx) = derv2(jz,kx) + d3trm*xnbr(nj,j)*zij
          derv2(jx,ky) = derv2(jx,ky) + d3trm*ynbr(nj,j)*xij
          derv2(jy,ky) = derv2(jy,ky) + d3trm*ynbr(nj,j)*yij
          derv2(jz,ky) = derv2(jz,ky) + d3trm*ynbr(nj,j)*zij
          derv2(jx,kz) = derv2(jx,kz) + d3trm*znbr(nj,j)*xij
          derv2(jy,kz) = derv2(jy,kz) + d3trm*znbr(nj,j)*yij
          derv2(jz,kz) = derv2(jz,kz) + d3trm*znbr(nj,j)*zij
        endif
!
!  I-K
!
        if (k.ge.i) then
          derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(nj,j)*xij
          derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(nj,j)*xij
          derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(nj,j)*xij
          derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(nj,j)*yij
          derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(nj,j)*yij
          derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(nj,j)*yij
          derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(nj,j)*zij
          derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(nj,j)*zij
          derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(nj,j)*zij
        else
          derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(nj,j)*xij
          derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(nj,j)*yij
          derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(nj,j)*zij
          derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(nj,j)*xij
          derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(nj,j)*yij
          derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(nj,j)*zij
          derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(nj,j)*xij
          derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(nj,j)*yij
          derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(nj,j)*zij
        endif
!
!  J-I
!
        if (i.ge.j) then
          derv2(ix,jx) = derv2(ix,jx) + d3trm*xnbr(nj,j)*xij
          derv2(iy,jx) = derv2(iy,jx) + d3trm*xnbr(nj,j)*yij
          derv2(iz,jx) = derv2(iz,jx) + d3trm*xnbr(nj,j)*zij
          derv2(ix,jy) = derv2(ix,jy) + d3trm*ynbr(nj,j)*xij
          derv2(iy,jy) = derv2(iy,jy) + d3trm*ynbr(nj,j)*yij
          derv2(iz,jy) = derv2(iz,jy) + d3trm*ynbr(nj,j)*zij
          derv2(ix,jz) = derv2(ix,jz) + d3trm*znbr(nj,j)*xij
          derv2(iy,jz) = derv2(iy,jz) + d3trm*znbr(nj,j)*yij
          derv2(iz,jz) = derv2(iz,jz) + d3trm*znbr(nj,j)*zij
        else
          derv2(jx,ix) = derv2(jx,ix) + d3trm*xnbr(nj,j)*xij
          derv2(jy,ix) = derv2(jy,ix) + d3trm*ynbr(nj,j)*xij
          derv2(jz,ix) = derv2(jz,ix) + d3trm*znbr(nj,j)*xij
          derv2(jx,iy) = derv2(jx,iy) + d3trm*xnbr(nj,j)*yij
          derv2(jy,iy) = derv2(jy,iy) + d3trm*ynbr(nj,j)*yij
          derv2(jz,iy) = derv2(jz,iy) + d3trm*znbr(nj,j)*yij
          derv2(jx,iz) = derv2(jx,iz) + d3trm*xnbr(nj,j)*zij
          derv2(jy,iz) = derv2(jy,iz) + d3trm*ynbr(nj,j)*zij
          derv2(jz,iz) = derv2(jz,iz) + d3trm*znbr(nj,j)*zij
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(nj,j)*d3trm*dr2ijds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) + xij*d3trm*dr2jkds(ks)
            derv3(jy,kl) = derv3(jy,kl) + yij*d3trm*dr2jkds(ks)
            derv3(jz,kl) = derv3(jz,kl) + zij*d3trm*dr2jkds(ks)
            derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2jkds(ks)
            derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2jkds(ks)
            derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2jkds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2jkds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2jkds(ks))
            enddo
          enddo
        endif
      endif
!
      d3trm1 = (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of j for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
      do n2 = 1,n-1
        dcnjdrjl = d1cndr_cn(n2,j)
        d3trm = d3trm1*dcnjdrjl
!
        if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
        nj2 = nbrno_cn(n2,j)
        l = nbrno(nj2,j)
        indl = 3*(l-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
!
        if (lstr) then
          call real1strterm(ndim,xnbr(nj2,j),ynbr(nj2,j),znbr(nj2,j),0.0_dp,0.0_dp,0.0_dp,dr2jlds,d2r2jldx2, &
                            d2r2jldsdx,d2r2jlds2,.false.)
        endif
!
!  J-K
!
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(ky,jx) = derv2(ky,jx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(kz,jx) = derv2(kz,jx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(kx,jy) = derv2(kx,jy) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(ky,jy) = derv2(ky,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(kz,jy) = derv2(kz,jy) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(kx,jz) = derv2(kx,jz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(ky,jz) = derv2(ky,jz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(kz,jz) = derv2(kz,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        else
          derv2(jx,kx) = derv2(jx,kx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(jy,kx) = derv2(jy,kx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(jz,kx) = derv2(jz,kx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(jx,ky) = derv2(jx,ky) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(jy,ky) = derv2(jy,ky) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(jz,ky) = derv2(jz,ky) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(jx,kz) = derv2(jx,kz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(jy,kz) = derv2(jy,kz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(jz,kz) = derv2(jz,kz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
!
!  J-L
!
        if (l.ge.j) then
          derv2(lx,jx) = derv2(lx,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(ly,jx) = derv2(ly,jx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(lz,jx) = derv2(lz,jx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(lx,jy) = derv2(lx,jy) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(ly,jy) = derv2(ly,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(lz,jy) = derv2(lz,jy) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(lx,jz) = derv2(lx,jz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(ly,jz) = derv2(ly,jz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(lz,jz) = derv2(lz,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        else
          derv2(jx,lx) = derv2(jx,lx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(jy,lx) = derv2(jy,lx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(jz,lx) = derv2(jz,lx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(jx,ly) = derv2(jx,ly) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(jy,ly) = derv2(jy,ly) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(jz,ly) = derv2(jz,ly) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(jx,lz) = derv2(jx,lz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(jy,lz) = derv2(jy,lz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(jz,lz) = derv2(jz,lz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
!
!  K-L
!     
        if (l.ge.k) then
          derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(nj,j)*znbr(nj2,j)
        else
          derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(nj,j)*d3trm*dr2jlds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(nj2,j)*d3trm*dr2jkds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2jkds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2jkds(ks))
            enddo
          enddo
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_x(i,j,xij,yij,zij,d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn,lgrad2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  This version is based on gfnff_drv2_dcn but only computes the cross terms between distance/CN
!
!  On entry : 
!
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!   5/22 Created from gfnff_drv2_dcn
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_x')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
!
  if (lgrad2.and.lstr) then
    call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2,d2r2ijdsdx,d2r2ijds2,.false.)
  endif
!
!  Set up for second derivatives
!
  if (lgrad2) then
    indi = 3*(i-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
!
    indj = 3*(j-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
  endif
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                        d2r2ikdsdx,d2r2ikds2,lgrad2)
    endif
    if (lgrad2) then
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
      d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
      if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  I-K
!
        if (k.ge.i) then
          derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xij
          derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xij
          derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xij
          derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*yij
          derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*yij
          derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*yij
          derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*zij
          derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*zij
          derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*zij
        else
          derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xij
          derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*yij
          derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*zij
          derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xij
          derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*yij
          derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*zij
          derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xij
          derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*yij
          derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*zij
        endif
!
!  J-K
!
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) + d3trm*xnbr(ni,i)*xij
          derv2(ky,jx) = derv2(ky,jx) + d3trm*ynbr(ni,i)*xij
          derv2(kz,jx) = derv2(kz,jx) + d3trm*znbr(ni,i)*xij
          derv2(kx,jy) = derv2(kx,jy) + d3trm*xnbr(ni,i)*yij
          derv2(ky,jy) = derv2(ky,jy) + d3trm*ynbr(ni,i)*yij
          derv2(kz,jy) = derv2(kz,jy) + d3trm*znbr(ni,i)*yij
          derv2(kx,jz) = derv2(kx,jz) + d3trm*xnbr(ni,i)*zij
          derv2(ky,jz) = derv2(ky,jz) + d3trm*ynbr(ni,i)*zij
          derv2(kz,jz) = derv2(kz,jz) + d3trm*znbr(ni,i)*zij
        else
          derv2(jx,kx) = derv2(jx,kx) + d3trm*xnbr(ni,i)*xij
          derv2(jy,kx) = derv2(jy,kx) + d3trm*xnbr(ni,i)*yij
          derv2(jz,kx) = derv2(jz,kx) + d3trm*xnbr(ni,i)*zij
          derv2(jx,ky) = derv2(jx,ky) + d3trm*ynbr(ni,i)*xij
          derv2(jy,ky) = derv2(jy,ky) + d3trm*ynbr(ni,i)*yij
          derv2(jz,ky) = derv2(jz,ky) + d3trm*ynbr(ni,i)*zij
          derv2(jx,kz) = derv2(jx,kz) + d3trm*znbr(ni,i)*xij
          derv2(jy,kz) = derv2(jy,kz) + d3trm*znbr(ni,i)*yij
          derv2(jz,kz) = derv2(jz,kz) + d3trm*znbr(ni,i)*zij
        endif
!
!  I-J
!
        if (j.ge.i) then
          derv2(jx,ix) = derv2(jx,ix) - d3trm*xnbr(ni,i)*xij
          derv2(jy,ix) = derv2(jy,ix) - d3trm*xnbr(ni,i)*yij
          derv2(jz,ix) = derv2(jz,ix) - d3trm*xnbr(ni,i)*zij
          derv2(jx,iy) = derv2(jx,iy) - d3trm*ynbr(ni,i)*xij
          derv2(jy,iy) = derv2(jy,iy) - d3trm*ynbr(ni,i)*yij
          derv2(jz,iy) = derv2(jz,iy) - d3trm*ynbr(ni,i)*zij
          derv2(jx,iz) = derv2(jx,iz) - d3trm*znbr(ni,i)*xij
          derv2(jy,iz) = derv2(jy,iz) - d3trm*znbr(ni,i)*yij
          derv2(jz,iz) = derv2(jz,iz) - d3trm*znbr(ni,i)*zij
        else
          derv2(ix,jx) = derv2(ix,jx) - d3trm*xnbr(ni,i)*xij
          derv2(iy,jx) = derv2(iy,jx) - d3trm*ynbr(ni,i)*xij
          derv2(iz,jx) = derv2(iz,jx) - d3trm*znbr(ni,i)*xij
          derv2(ix,jy) = derv2(ix,jy) - d3trm*xnbr(ni,i)*yij
          derv2(iy,jy) = derv2(iy,jy) - d3trm*ynbr(ni,i)*yij
          derv2(iz,jy) = derv2(iz,jy) - d3trm*znbr(ni,i)*yij
          derv2(ix,jz) = derv2(ix,jz) - d3trm*xnbr(ni,i)*zij
          derv2(iy,jz) = derv2(iy,jz) - d3trm*ynbr(ni,i)*zij
          derv2(iz,jz) = derv2(iz,jz) - d3trm*znbr(ni,i)*zij
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ijds(ks)
          enddo
!
         do kl = 1,nstrains
           ks = nstrptr(kl)
           derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2ikds(ks)
           derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2ikds(ks)
           derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2ikds(ks)
           derv3(jx,kl) = derv3(jx,kl) + xij*d3trm*dr2ikds(ks)
           derv3(jy,kl) = derv3(jy,kl) + yij*d3trm*dr2ikds(ks)
           derv3(jz,kl) = derv3(jz,kl) + zij*d3trm*dr2ikds(ks)
         enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2ikds(ks)) 
            enddo
          enddo
        endif
      endif
    endif
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    dcnjdrjk = d1cndr_cn(n,j)
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jkds,d2r2jkdx2, &
                        d2r2jkdsdx,d2r2jkds2,lgrad2)
    endif
    if (lgrad2) then
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
      d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
      if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  J-K
!
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) + d3trm*xnbr(nj,j)*xij
          derv2(ky,jx) = derv2(ky,jx) + d3trm*ynbr(nj,j)*xij
          derv2(kz,jx) = derv2(kz,jx) + d3trm*znbr(nj,j)*xij
          derv2(kx,jy) = derv2(kx,jy) + d3trm*xnbr(nj,j)*yij
          derv2(ky,jy) = derv2(ky,jy) + d3trm*ynbr(nj,j)*yij
          derv2(kz,jy) = derv2(kz,jy) + d3trm*znbr(nj,j)*yij
          derv2(kx,jz) = derv2(kx,jz) + d3trm*xnbr(nj,j)*zij
          derv2(ky,jz) = derv2(ky,jz) + d3trm*ynbr(nj,j)*zij
          derv2(kz,jz) = derv2(kz,jz) + d3trm*znbr(nj,j)*zij
        else
          derv2(jx,kx) = derv2(jx,kx) + d3trm*xnbr(nj,j)*xij
          derv2(jy,kx) = derv2(jy,kx) + d3trm*xnbr(nj,j)*yij
          derv2(jz,kx) = derv2(jz,kx) + d3trm*xnbr(nj,j)*zij
          derv2(jx,ky) = derv2(jx,ky) + d3trm*ynbr(nj,j)*xij
          derv2(jy,ky) = derv2(jy,ky) + d3trm*ynbr(nj,j)*yij
          derv2(jz,ky) = derv2(jz,ky) + d3trm*ynbr(nj,j)*zij
          derv2(jx,kz) = derv2(jx,kz) + d3trm*znbr(nj,j)*xij
          derv2(jy,kz) = derv2(jy,kz) + d3trm*znbr(nj,j)*yij
          derv2(jz,kz) = derv2(jz,kz) + d3trm*znbr(nj,j)*zij
        endif
!
!  I-K
!
        if (k.ge.i) then
          derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(nj,j)*xij
          derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(nj,j)*xij
          derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(nj,j)*xij
          derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(nj,j)*yij
          derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(nj,j)*yij
          derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(nj,j)*yij
          derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(nj,j)*zij
          derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(nj,j)*zij
          derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(nj,j)*zij
        else
          derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(nj,j)*xij
          derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(nj,j)*yij
          derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(nj,j)*zij
          derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(nj,j)*xij
          derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(nj,j)*yij
          derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(nj,j)*zij
          derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(nj,j)*xij
          derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(nj,j)*yij
          derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(nj,j)*zij
        endif
!
!  J-I
!
        if (i.ge.j) then
          derv2(ix,jx) = derv2(ix,jx) + d3trm*xnbr(nj,j)*xij
          derv2(iy,jx) = derv2(iy,jx) + d3trm*xnbr(nj,j)*yij
          derv2(iz,jx) = derv2(iz,jx) + d3trm*xnbr(nj,j)*zij
          derv2(ix,jy) = derv2(ix,jy) + d3trm*ynbr(nj,j)*xij
          derv2(iy,jy) = derv2(iy,jy) + d3trm*ynbr(nj,j)*yij
          derv2(iz,jy) = derv2(iz,jy) + d3trm*ynbr(nj,j)*zij
          derv2(ix,jz) = derv2(ix,jz) + d3trm*znbr(nj,j)*xij
          derv2(iy,jz) = derv2(iy,jz) + d3trm*znbr(nj,j)*yij
          derv2(iz,jz) = derv2(iz,jz) + d3trm*znbr(nj,j)*zij
        else
          derv2(jx,ix) = derv2(jx,ix) + d3trm*xnbr(nj,j)*xij
          derv2(jy,ix) = derv2(jy,ix) + d3trm*ynbr(nj,j)*xij
          derv2(jz,ix) = derv2(jz,ix) + d3trm*znbr(nj,j)*xij
          derv2(jx,iy) = derv2(jx,iy) + d3trm*xnbr(nj,j)*yij
          derv2(jy,iy) = derv2(jy,iy) + d3trm*ynbr(nj,j)*yij
          derv2(jz,iy) = derv2(jz,iy) + d3trm*znbr(nj,j)*yij
          derv2(jx,iz) = derv2(jx,iz) + d3trm*xnbr(nj,j)*zij
          derv2(jy,iz) = derv2(jy,iz) + d3trm*ynbr(nj,j)*zij
          derv2(jz,iz) = derv2(jz,iz) + d3trm*znbr(nj,j)*zij
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(nj,j)*d3trm*dr2ijds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) + xij*d3trm*dr2jkds(ks)
            derv3(jy,kl) = derv3(jy,kl) + yij*d3trm*dr2jkds(ks)
            derv3(jz,kl) = derv3(jz,kl) + zij*d3trm*dr2jkds(ks)
            derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2jkds(ks)
            derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2jkds(ks)
            derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2jkds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2jkds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2jkds(ks))
            enddo
          enddo
        endif
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_x')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_nox(i,j,dEdlogcni,dEdlogcnj,d2Edlogcni2,d2Edlogcnj2,d2Edlogcnidlogcnj, &
                                dlogcndcn,d2logcndcn2,lgrad2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  This version is like gfnff_drv2_dcn except that no cross terms are added for distance/CN
!
!  On entry : 
!
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  dEdlogcnj         = derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcnj2       = second derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcnidlogcnj = second derivatives of energy w.r.t. log of the coordination number for i and j
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!   5/22 Created from gfnff_drv2_dcn
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: dEdlogcnj
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: d2Edlogcnidlogcnj
  real(dp),    intent(in)                          :: d2Edlogcnj2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: d3trm1
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dcnjdrjl
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2logcnjdcnj2
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2ilds(6)
  real(dp)                                         :: d2r2ildx2(6)
  real(dp)                                         :: d2r2ilds2(6,6)
  real(dp)                                         :: d2r2ildsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
  real(dp)                                         :: dr2jlds(6)
  real(dp)                                         :: d2r2jldx2(6)
  real(dp)                                         :: d2r2jlds2(6,6)
  real(dp)                                         :: d2r2jldsdx(6,3)
  real(dp)                                         :: dxyz(3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_nox')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  if (lgrad2) then
    d2logcnidcni2 = d2logcndcn2(i)
    d2logcnjdcnj2 = d2logcndcn2(j)
  endif
!
!  Set up for second derivatives
!
  if (lgrad2) then
    indi = 3*(i-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
!
    indj = 3*(j-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
  endif
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
    dxyz(1) = d1trm*xnbr(ni,i)
    dxyz(2) = d1trm*ynbr(ni,i)
    dxyz(3) = d1trm*znbr(ni,i)
!
    xdrv(i) = xdrv(i) - dxyz(1)
    ydrv(i) = ydrv(i) - dxyz(2)
    zdrv(i) = zdrv(i) - dxyz(3)
    xdrv(k) = xdrv(k) + dxyz(1)
    ydrv(k) = ydrv(k) + dxyz(2)
    zdrv(k) = zdrv(k) + dxyz(3)
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                        d2r2ikdsdx,d2r2ikds2,lgrad2)
    endif
    if (lstr) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        rstrd(kl) = rstrd(kl) + d1trm*dr2ikds(ks)
      enddo
    endif
    if (lgrad2) then
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
!
      d2cnidrik2 = d2cndr_cn(n,i)
      d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
             (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!
      if (abs(d2trm).gt.gfnff_cnc6tol) then
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
        if (k.ge.i) then
          derv2(kx,ix) = derv2(kx,ix) - d2trm*d2r2ikdx2(1)
          derv2(ky,ix) = derv2(ky,ix) - d2trm*d2r2ikdx2(6)
          derv2(kz,ix) = derv2(kz,ix) - d2trm*d2r2ikdx2(5)
          derv2(kx,iy) = derv2(kx,iy) - d2trm*d2r2ikdx2(6)
          derv2(ky,iy) = derv2(ky,iy) - d2trm*d2r2ikdx2(2)
          derv2(kz,iy) = derv2(kz,iy) - d2trm*d2r2ikdx2(4)
          derv2(kx,iz) = derv2(kx,iz) - d2trm*d2r2ikdx2(5)
          derv2(ky,iz) = derv2(ky,iz) - d2trm*d2r2ikdx2(4)
          derv2(kz,iz) = derv2(kz,iz) - d2trm*d2r2ikdx2(3)
          derv2(kx,ix) = derv2(kx,ix) - d1trm
          derv2(ky,iy) = derv2(ky,iy) - d1trm
          derv2(kz,iz) = derv2(kz,iz) - d1trm
        else
          derv2(ix,kx) = derv2(ix,kx) - d2trm*d2r2ikdx2(1)
          derv2(iy,kx) = derv2(iy,kx) - d2trm*d2r2ikdx2(6)
          derv2(iz,kx) = derv2(iz,kx) - d2trm*d2r2ikdx2(5)
          derv2(ix,ky) = derv2(ix,ky) - d2trm*d2r2ikdx2(6)
          derv2(iy,ky) = derv2(iy,ky) - d2trm*d2r2ikdx2(2)
          derv2(iz,ky) = derv2(iz,ky) - d2trm*d2r2ikdx2(4)
          derv2(ix,kz) = derv2(ix,kz) - d2trm*d2r2ikdx2(5)
          derv2(iy,kz) = derv2(iy,kz) - d2trm*d2r2ikdx2(4)
          derv2(iz,kz) = derv2(iz,kz) - d2trm*d2r2ikdx2(3)
          derv2(ix,kx) = derv2(ix,kx) - d1trm
          derv2(iy,ky) = derv2(iy,ky) - d1trm
          derv2(iz,kz) = derv2(iz,kz) - d1trm
        endif
!
        if (lstr) then
!
!  Mixed derivatives
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ikdsdx(ks,1) - xnbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ikdsdx(ks,2) - ynbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ikdsdx(ks,3) - znbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2ikdsdx(ks,1) + xnbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2ikdsdx(ks,2) + ynbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2ikdsdx(ks,3) + znbr(ni,i)*d2trm*dr2ikds(ks)
          enddo
!
!  Strain-strain
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2ikds(kt)*dr2ikds(ks) + d1trm*d2r2ikds2(kt,ks)
            enddo
          enddo
        endif
      endif
!
      d3trm1 = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
      do n2 = 1,n-1
        dcnidril = d1cndr_cn(n2,i)
        d3trm = d3trm1*dcnidril
!
        if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
        ni2 = nbrno_cn(n2,i)
        l = nbrno(ni2,i)
        indl = 3*(l-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
!
        if (lstr) then
          call real1strterm(ndim,xnbr(ni2,i),ynbr(ni2,i),znbr(ni2,i),0.0_dp,0.0_dp,0.0_dp,dr2ilds,d2r2ildx2, &
                            d2r2ildsdx,d2r2ilds2,.false.)
        endif
!
!  I-K
!
        if (k.ge.i) then
          derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        else
          derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  I-L
!
        if (l.ge.i) then
          derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ly,ix) = derv2(ly,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(lz,ix) = derv2(lz,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(lx,iy) = derv2(lx,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(lz,iy) = derv2(lz,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(lx,iz) = derv2(lx,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ly,iz) = derv2(ly,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        else
          derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(iy,lx) = derv2(iy,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(iz,lx) = derv2(iz,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ix,ly) = derv2(ix,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(iz,ly) = derv2(iz,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(ix,lz) = derv2(ix,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(iy,lz) = derv2(iy,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  K-L
!
        if (l.ge.k) then
          derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
        else
          derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ilds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(ni2,i)*d3trm*dr2ikds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ilds(ks) + dr2ilds(kt)*dr2ikds(ks))
            enddo
          enddo
        endif
      enddo
!
      d3trm1 = d2Edlogcnidlogcnj*dlogcnidcni*dlogcnjdcnj*dcnidrik
!------------------------------------------------------------------------------------------
!  Loop over neighbours of j for second derivatives for 2 different coordination numbers  |
!------------------------------------------------------------------------------------------
      do n2 = 1,nnbr_cn(j)
        dcnjdrjl = d1cndr_cn(n2,j)
        d3trm = d3trm1*dcnjdrjl
!
        if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
        nj = nbrno_cn(n2,j)
        l = nbrno(nj,j)
        indl = 3*(l-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
!
        if (lstr) then
          call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jlds,d2r2jldx2, &
                            d2r2jldsdx,d2r2jlds2,.false.)
        endif
!
!  I-J
!
        if (j.ge.i) then
          derv2(jx,ix) = derv2(jx,ix) + d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(jy,ix) = derv2(jy,ix) + d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(jz,ix) = derv2(jz,ix) + d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(jx,iy) = derv2(jx,iy) + d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(jy,iy) = derv2(jy,iy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(jz,iy) = derv2(jz,iy) + d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(jx,iz) = derv2(jx,iz) + d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(jy,iz) = derv2(jy,iz) + d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(jz,iz) = derv2(jz,iz) + d3trm*znbr(nj,j)*znbr(ni,i)
        else
          derv2(ix,jx) = derv2(ix,jx) + d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(iy,jx) = derv2(iy,jx) + d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(iz,jx) = derv2(iz,jx) + d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(ix,jy) = derv2(ix,jy) + d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(iy,jy) = derv2(iy,jy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(iz,jy) = derv2(iz,jy) + d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(ix,jz) = derv2(ix,jz) + d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(iy,jz) = derv2(iy,jz) + d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(iz,jz) = derv2(iz,jz) + d3trm*znbr(nj,j)*znbr(ni,i)
        endif
!
!  I-L
!
        if (l.ge.i) then
          derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(ly,ix) = derv2(ly,ix) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(lz,ix) = derv2(lz,ix) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(lx,iy) = derv2(lx,iy) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(lz,iy) = derv2(lz,iy) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(lx,iz) = derv2(lx,iz) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(ly,iz) = derv2(ly,iz) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(nj,j)*znbr(ni,i)
        else
          derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(iy,lx) = derv2(iy,lx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(iz,lx) = derv2(iz,lx) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(ix,ly) = derv2(ix,ly) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(iz,ly) = derv2(iz,ly) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(ix,lz) = derv2(ix,lz) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(iy,lz) = derv2(iy,lz) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(nj,j)*znbr(ni,i)
        endif
!
!  J-K
!
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(ky,jx) = derv2(ky,jx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(kz,jx) = derv2(kz,jx) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(kx,jy) = derv2(kx,jy) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(ky,jy) = derv2(ky,jy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(kz,jy) = derv2(kz,jy) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(kx,jz) = derv2(kx,jz) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(ky,jz) = derv2(ky,jz) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(kz,jz) = derv2(kz,jz) - d3trm*znbr(nj,j)*znbr(ni,i)
        else
          derv2(jx,kx) = derv2(jx,kx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(jy,kx) = derv2(jy,kx) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(jz,kx) = derv2(jz,kx) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(jx,ky) = derv2(jx,ky) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(jy,ky) = derv2(jy,ky) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(jz,ky) = derv2(jz,ky) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(jx,kz) = derv2(jx,kz) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(jy,kz) = derv2(jy,kz) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(jz,kz) = derv2(jz,kz) - d3trm*znbr(nj,j)*znbr(ni,i)
        endif
!
!  K-L
!     
        if (l.ge.k) then
          derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
          derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(nj,j)
          derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(nj,j)
          derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(nj,j)
          derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(nj,j)
          derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(nj,j)
          derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(nj,j)
          derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(nj,j)
          derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(nj,j)
        else
          derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
          derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(nj,j)
          derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(nj,j)
          derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(nj,j)
          derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(nj,j)
          derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(nj,j)
          derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(nj,j)
          derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(nj,j)
          derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(nj,j)
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2jlds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(nj,j)*d3trm*dr2ikds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2ikds(ks))
            enddo
          enddo
        endif
      enddo
    endif
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    dcnjdrjk = d1cndr_cn(n,j)
    d1trm = dlogcnjdcnj*dcnjdrjk*dEdlogcnj
    dxyz(1) = d1trm*xnbr(nj,j)
    dxyz(2) = d1trm*ynbr(nj,j)
    dxyz(3) = d1trm*znbr(nj,j)
!
    xdrv(j) = xdrv(j) - dxyz(1)
    ydrv(j) = ydrv(j) - dxyz(2)
    zdrv(j) = zdrv(j) - dxyz(3)
    xdrv(k) = xdrv(k) + dxyz(1)
    ydrv(k) = ydrv(k) + dxyz(2)
    zdrv(k) = zdrv(k) + dxyz(3)
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jkds,d2r2jkdx2, &
                        d2r2jkdsdx,d2r2jkds2,lgrad2)
    endif
    if (lstr) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        rstrd(kl) = rstrd(kl) + d1trm*dr2jkds(ks)
      enddo
    endif
    if (lgrad2) then
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
!
      d2cnjdrjk2 = d2cndr_cn(n,j)
      d2trm = dEdlogcnj*dlogcnjdcnj*d2cnjdrjk2 + &
             (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjk
!
      if (abs(d2trm).gt.gfnff_cnc6tol) then
!----------------------------------------
!  Second derivatives of CN w.r.t. rjk  |
!----------------------------------------
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) - d2trm*d2r2jkdx2(1)
          derv2(ky,jx) = derv2(ky,jx) - d2trm*d2r2jkdx2(6)
          derv2(kz,jx) = derv2(kz,jx) - d2trm*d2r2jkdx2(5)
          derv2(kx,jy) = derv2(kx,jy) - d2trm*d2r2jkdx2(6)
          derv2(ky,jy) = derv2(ky,jy) - d2trm*d2r2jkdx2(2)
          derv2(kz,jy) = derv2(kz,jy) - d2trm*d2r2jkdx2(4)
          derv2(kx,jz) = derv2(kx,jz) - d2trm*d2r2jkdx2(5)
          derv2(ky,jz) = derv2(ky,jz) - d2trm*d2r2jkdx2(4)
          derv2(kz,jz) = derv2(kz,jz) - d2trm*d2r2jkdx2(3)
          derv2(kx,jx) = derv2(kx,jx) - d1trm
          derv2(ky,jy) = derv2(ky,jy) - d1trm
          derv2(kz,jz) = derv2(kz,jz) - d1trm
        else
          derv2(jx,kx) = derv2(jx,kx) - d2trm*d2r2jkdx2(1)
          derv2(jy,kx) = derv2(jy,kx) - d2trm*d2r2jkdx2(6)
          derv2(jz,kx) = derv2(jz,kx) - d2trm*d2r2jkdx2(5)
          derv2(jx,ky) = derv2(jx,ky) - d2trm*d2r2jkdx2(6)
          derv2(jy,ky) = derv2(jy,ky) - d2trm*d2r2jkdx2(2)
          derv2(jz,ky) = derv2(jz,ky) - d2trm*d2r2jkdx2(4)
          derv2(jx,kz) = derv2(jx,kz) - d2trm*d2r2jkdx2(5)
          derv2(jy,kz) = derv2(jy,kz) - d2trm*d2r2jkdx2(4)
          derv2(jz,kz) = derv2(jz,kz) - d2trm*d2r2jkdx2(3)
          derv2(jx,kx) = derv2(jx,kx) - d1trm
          derv2(jy,ky) = derv2(jy,ky) - d1trm
          derv2(jz,kz) = derv2(jz,kz) - d1trm
        endif
!
        if (lstr) then
!
!  Mixed derivatives
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - d1trm*d2r2jkdsdx(ks,1) - xnbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(jy,kl) = derv3(jy,kl) - d1trm*d2r2jkdsdx(ks,2) - ynbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(jz,kl) = derv3(jz,kl) - d1trm*d2r2jkdsdx(ks,3) - znbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2jkdsdx(ks,1) + xnbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2jkdsdx(ks,2) + ynbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2jkdsdx(ks,3) + znbr(nj,j)*d2trm*dr2jkds(ks)
          enddo
!
!  Strain-strain
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2jkds(kt)*dr2jkds(ks) + d1trm*d2r2jkds2(kt,ks)
            enddo
          enddo
        endif
      endif
!
      d3trm1 = (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of j for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
      do n2 = 1,n-1
        dcnjdrjl = d1cndr_cn(n2,j)
        d3trm = d3trm1*dcnjdrjl
!
        if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
        nj2 = nbrno_cn(n2,j)
        l = nbrno(nj2,j)
        indl = 3*(l-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
!
        if (lstr) then
          call real1strterm(ndim,xnbr(nj2,j),ynbr(nj2,j),znbr(nj2,j),0.0_dp,0.0_dp,0.0_dp,dr2jlds,d2r2jldx2, &
                            d2r2jldsdx,d2r2jlds2,.false.)
        endif
!
!  J-K
!
        if (k.ge.j) then
          derv2(kx,jx) = derv2(kx,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(ky,jx) = derv2(ky,jx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(kz,jx) = derv2(kz,jx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(kx,jy) = derv2(kx,jy) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(ky,jy) = derv2(ky,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(kz,jy) = derv2(kz,jy) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(kx,jz) = derv2(kx,jz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(ky,jz) = derv2(ky,jz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(kz,jz) = derv2(kz,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        else
          derv2(jx,kx) = derv2(jx,kx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(jy,kx) = derv2(jy,kx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(jz,kx) = derv2(jz,kx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(jx,ky) = derv2(jx,ky) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(jy,ky) = derv2(jy,ky) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(jz,ky) = derv2(jz,ky) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(jx,kz) = derv2(jx,kz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(jy,kz) = derv2(jy,kz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(jz,kz) = derv2(jz,kz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
!
!  J-L
!
        if (l.ge.j) then
          derv2(lx,jx) = derv2(lx,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(ly,jx) = derv2(ly,jx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(lz,jx) = derv2(lz,jx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(lx,jy) = derv2(lx,jy) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(ly,jy) = derv2(ly,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(lz,jy) = derv2(lz,jy) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(lx,jz) = derv2(lx,jz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(ly,jz) = derv2(ly,jz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(lz,jz) = derv2(lz,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        else
          derv2(jx,lx) = derv2(jx,lx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(jy,lx) = derv2(jy,lx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(jz,lx) = derv2(jz,lx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(jx,ly) = derv2(jx,ly) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(jy,ly) = derv2(jy,ly) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(jz,ly) = derv2(jz,ly) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(jx,lz) = derv2(jx,lz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(jy,lz) = derv2(jy,lz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(jz,lz) = derv2(jz,lz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
!
!  K-L
!     
        if (l.ge.k) then
          derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(nj,j)*znbr(nj2,j)
        else
          derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(nj,j)*d3trm*dr2jlds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(nj2,j)*d3trm*dr2jkds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2jkds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2jkds(ks))
            enddo
          enddo
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_nox')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnd(lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,j,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                             dEdlogcni,dEdlogcnj,d2Edlogcni2,d2Edlogcnj2,d2Edlogcnidlogcnj, &
                             d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn,d2logcndcn2,lgrad2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Distributed memory parallel version
!
!  On entry : 
!
!  lilocal           = if true then i is local to this processor
!  ljlocal           = if true then j is local to this processor
!  ix,  iy,  iz      = second derivative Cartesian elements for i (local)
!  ixf, iyf, izf     = second derivative Cartesian elements for i (global)
!  jx,  jy,  jz      = second derivative Cartesian elements for j (local)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  dEdlogcnj         = derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcnj2       = second derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcnidlogcnj = second derivatives of energy w.r.t. log of the coordination number for i and j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!   3/22 Created from gfnff_drv2_dcn 
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
!  Julian Gale, CIC, Curtin University, March 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use parallel
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: dEdlogcnj
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: d2Edlogcnidlogcnj
  real(dp),    intent(in)                          :: d2Edlogcnj2
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
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
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  logical                                          :: liklocal
  logical                                          :: lijklocal
  logical                                          :: lklocal
  logical                                          :: lllocal
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: d3trm1
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dcnjdrjl
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2logcnjdcnj2
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2ilds(6)
  real(dp)                                         :: d2r2ildx2(6)
  real(dp)                                         :: d2r2ilds2(6,6)
  real(dp)                                         :: d2r2ildsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
  real(dp)                                         :: dr2jlds(6)
  real(dp)                                         :: d2r2jldx2(6)
  real(dp)                                         :: d2r2jlds2(6,6)
  real(dp)                                         :: d2r2jldsdx(6,3)
  real(dp)                                         :: dxyz(3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnd')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  if (lgrad2) then
    d2logcnidcni2 = d2logcndcn2(i)
    d2logcnjdcnj2 = d2logcndcn2(j)
  endif
  if (lgrad2.and.lstr) then
    call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2,d2r2ijdsdx,d2r2ijds2,.false.)
  endif
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
    dxyz(1) = d1trm*xnbr(ni,i)
    dxyz(2) = d1trm*ynbr(ni,i)
    dxyz(3) = d1trm*znbr(ni,i)
!
    if (lilocal) then
      xdrv(i) = xdrv(i) - dxyz(1)
      ydrv(i) = ydrv(i) - dxyz(2)
      zdrv(i) = zdrv(i) - dxyz(3)
      xdrv(k) = xdrv(k) + dxyz(1)
      ydrv(k) = ydrv(k) + dxyz(2)
      zdrv(k) = zdrv(k) + dxyz(3)
    endif
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                        d2r2ikdsdx,d2r2ikds2,lgrad2)
    endif
    if (lstr) then
      if (lilocal) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + d1trm*dr2ikds(ks)
        enddo
      endif
    endif
    if (lgrad2) then
      indk = 3*(k-1)
      kxf = indk + 1
      kyf = indk + 2
      kzf = indk + 3
!
      kloc = atom2local(k)
      lklocal = (kloc.ne.0)
      liklocal = (lilocal.or.lklocal)
      lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
      if (lklocal) then
        indk = 3*(kloc-1)
        kx = indk + 1
        ky = indk + 2
        kz = indk + 3
      endif
!
      d2cnidrik2 = d2cndr_cn(n,i)
      d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
             (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!
      if (abs(d2trm).gt.gfnff_cnc6tol) then
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
        if (lilocal) then
          derv2(kxf,ix) = derv2(kxf,ix) - d2trm*d2r2ikdx2(1)
          derv2(kyf,ix) = derv2(kyf,ix) - d2trm*d2r2ikdx2(6)
          derv2(kzf,ix) = derv2(kzf,ix) - d2trm*d2r2ikdx2(5)
          derv2(kxf,iy) = derv2(kxf,iy) - d2trm*d2r2ikdx2(6)
          derv2(kyf,iy) = derv2(kyf,iy) - d2trm*d2r2ikdx2(2)
          derv2(kzf,iy) = derv2(kzf,iy) - d2trm*d2r2ikdx2(4)
          derv2(kxf,iz) = derv2(kxf,iz) - d2trm*d2r2ikdx2(5)
          derv2(kyf,iz) = derv2(kyf,iz) - d2trm*d2r2ikdx2(4)
          derv2(kzf,iz) = derv2(kzf,iz) - d2trm*d2r2ikdx2(3)
          derv2(kxf,ix) = derv2(kxf,ix) - d1trm
          derv2(kyf,iy) = derv2(kyf,iy) - d1trm
          derv2(kzf,iz) = derv2(kzf,iz) - d1trm
        endif
        if (lklocal) then
          derv2(ixf,kx) = derv2(ixf,kx) - d2trm*d2r2ikdx2(1)
          derv2(iyf,kx) = derv2(iyf,kx) - d2trm*d2r2ikdx2(6)
          derv2(izf,kx) = derv2(izf,kx) - d2trm*d2r2ikdx2(5)
          derv2(ixf,ky) = derv2(ixf,ky) - d2trm*d2r2ikdx2(6)
          derv2(iyf,ky) = derv2(iyf,ky) - d2trm*d2r2ikdx2(2)
          derv2(izf,ky) = derv2(izf,ky) - d2trm*d2r2ikdx2(4)
          derv2(ixf,kz) = derv2(ixf,kz) - d2trm*d2r2ikdx2(5)
          derv2(iyf,kz) = derv2(iyf,kz) - d2trm*d2r2ikdx2(4)
          derv2(izf,kz) = derv2(izf,kz) - d2trm*d2r2ikdx2(3)
          derv2(ixf,kx) = derv2(ixf,kx) - d1trm
          derv2(iyf,ky) = derv2(iyf,ky) - d1trm
          derv2(izf,kz) = derv2(izf,kz) - d1trm
        endif
!
        if (lstr) then
!
!  Mixed derivatives
!
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ikdsdx(ks,1) - xnbr(ni,i)*d2trm*dr2ikds(ks)
              derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ikdsdx(ks,2) - ynbr(ni,i)*d2trm*dr2ikds(ks)
              derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ikdsdx(ks,3) - znbr(ni,i)*d2trm*dr2ikds(ks)
            enddo
          endif
          if (lklocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2ikdsdx(ks,1) + xnbr(ni,i)*d2trm*dr2ikds(ks)
              derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2ikdsdx(ks,2) + ynbr(ni,i)*d2trm*dr2ikds(ks)
              derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2ikdsdx(ks,3) + znbr(ni,i)*d2trm*dr2ikds(ks)
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
                sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2ikds(kt)*dr2ikds(ks) + d1trm*d2r2ikds2(kt,ks)
              enddo
            enddo
          endif
        endif
      endif
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
      d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
      if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  I-K
!
        if (lilocal) then
          derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xij
          derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xij
          derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xij
          derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*yij
          derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*yij
          derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*yij
          derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*zij
          derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*zij
          derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*zij
        endif
        if (lklocal) then
          derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xij
          derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*yij
          derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*zij
          derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xij
          derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*yij
          derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*zij
          derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xij
          derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*yij
          derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*zij
        endif
!
!  J-K
!
        if (ljlocal) then
          derv2(kxf,jx) = derv2(kxf,jx) + d3trm*xnbr(ni,i)*xij
          derv2(kyf,jx) = derv2(kyf,jx) + d3trm*ynbr(ni,i)*xij
          derv2(kzf,jx) = derv2(kzf,jx) + d3trm*znbr(ni,i)*xij
          derv2(kxf,jy) = derv2(kxf,jy) + d3trm*xnbr(ni,i)*yij
          derv2(kyf,jy) = derv2(kyf,jy) + d3trm*ynbr(ni,i)*yij
          derv2(kzf,jy) = derv2(kzf,jy) + d3trm*znbr(ni,i)*yij
          derv2(kxf,jz) = derv2(kxf,jz) + d3trm*xnbr(ni,i)*zij
          derv2(kyf,jz) = derv2(kyf,jz) + d3trm*ynbr(ni,i)*zij
          derv2(kzf,jz) = derv2(kzf,jz) + d3trm*znbr(ni,i)*zij
        endif
        if (lklocal) then
          derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(ni,i)*xij
          derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(ni,i)*yij
          derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(ni,i)*zij
          derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(ni,i)*xij
          derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(ni,i)*yij
          derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(ni,i)*zij
          derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(ni,i)*xij
          derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(ni,i)*yij
          derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(ni,i)*zij
        endif
!
!  I-J
!
        if (lilocal) then
          derv2(jxf,ix) = derv2(jxf,ix) - d3trm*xnbr(ni,i)*xij
          derv2(jyf,ix) = derv2(jyf,ix) - d3trm*xnbr(ni,i)*yij
          derv2(jzf,ix) = derv2(jzf,ix) - d3trm*xnbr(ni,i)*zij
          derv2(jxf,iy) = derv2(jxf,iy) - d3trm*ynbr(ni,i)*xij
          derv2(jyf,iy) = derv2(jyf,iy) - d3trm*ynbr(ni,i)*yij
          derv2(jzf,iy) = derv2(jzf,iy) - d3trm*ynbr(ni,i)*zij
          derv2(jxf,iz) = derv2(jxf,iz) - d3trm*znbr(ni,i)*xij
          derv2(jyf,iz) = derv2(jyf,iz) - d3trm*znbr(ni,i)*yij
          derv2(jzf,iz) = derv2(jzf,iz) - d3trm*znbr(ni,i)*zij
        endif
        if (ljlocal) then
          derv2(ixf,jx) = derv2(ixf,jx) - d3trm*xnbr(ni,i)*xij
          derv2(iyf,jx) = derv2(iyf,jx) - d3trm*ynbr(ni,i)*xij
          derv2(izf,jx) = derv2(izf,jx) - d3trm*znbr(ni,i)*xij
          derv2(ixf,jy) = derv2(ixf,jy) - d3trm*xnbr(ni,i)*yij
          derv2(iyf,jy) = derv2(iyf,jy) - d3trm*ynbr(ni,i)*yij
          derv2(izf,jy) = derv2(izf,jy) - d3trm*znbr(ni,i)*yij
          derv2(ixf,jz) = derv2(ixf,jz) - d3trm*xnbr(ni,i)*zij
          derv2(iyf,jz) = derv2(iyf,jz) - d3trm*ynbr(ni,i)*zij
          derv2(izf,jz) = derv2(izf,jz) - d3trm*znbr(ni,i)*zij
        endif
!
!  Strain terms
!
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ijds(ks)
              derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ijds(ks)
              derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ijds(ks)
              derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2ikds(ks)
              derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2ikds(ks)
              derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2ikds(ks)
            enddo
!
            do kk = 1,nstrains
              ks = nstrptr(kk)
              do kl = 1,nstrains
                kt = nstrptr(kl)
                sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2ikds(ks)) 
              enddo
            enddo
          endif
!
          if (ljlocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(jx,kl) = derv3(jx,kl) + xij*d3trm*dr2ikds(ks)
              derv3(jy,kl) = derv3(jy,kl) + yij*d3trm*dr2ikds(ks)
              derv3(jz,kl) = derv3(jz,kl) + zij*d3trm*dr2ikds(ks)
            enddo
          endif
!
          if (lklocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ijds(ks)
              derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ijds(ks)
              derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ijds(ks)
            enddo
          endif
        endif
      endif
!
      d3trm1 = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
      do n2 = 1,n-1
        dcnidril = d1cndr_cn(n2,i)
        d3trm = d3trm1*dcnidril
!
        if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
        ni2 = nbrno_cn(n2,i)
        l = nbrno(ni2,i)
!
        lloc = atom2local(l)
        lllocal = (lloc.ne.0)
!
!  If none of the atoms are local then cycle
!
        if (.not.liklocal.and..not.lllocal) cycle
!
        indl = 3*(l-1)
        lxf = indl + 1
        lyf = indl + 2
        lzf = indl + 3
!
        if (lllocal) then
          indl = 3*(lloc-1)
          lx = indl + 1
          ly = indl + 2
          lz = indl + 3
        endif
!
        if (lstr) then
          call real1strterm(ndim,xnbr(ni2,i),ynbr(ni2,i),znbr(ni2,i),0.0_dp,0.0_dp,0.0_dp,dr2ilds,d2r2ildx2, &
                            d2r2ildsdx,d2r2ilds2,.false.)
        endif
!
!  I-K
!
        if (lilocal) then
          derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
        if (lklocal) then
          derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  I-L
!
        if (lilocal) then
          derv2(lxf,ix) = derv2(lxf,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(lyf,ix) = derv2(lyf,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(lzf,ix) = derv2(lzf,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(lxf,iy) = derv2(lxf,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(lyf,iy) = derv2(lyf,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(lzf,iy) = derv2(lzf,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(lxf,iz) = derv2(lxf,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(lyf,iz) = derv2(lyf,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(lzf,iz) = derv2(lzf,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
        if (lllocal) then
          derv2(ixf,lx) = derv2(ixf,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(iyf,lx) = derv2(iyf,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(izf,lx) = derv2(izf,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ixf,ly) = derv2(ixf,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(iyf,ly) = derv2(iyf,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(izf,ly) = derv2(izf,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(ixf,lz) = derv2(ixf,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(iyf,lz) = derv2(iyf,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(izf,lz) = derv2(izf,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  K-L
!
        if (lklocal) then
          derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
        if (lllocal) then
          derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  Strain terms
!
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ilds(ks)
              derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ilds(ks)
              derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ilds(ks)
              derv3(ix,kl) = derv3(ix,kl) - xnbr(ni2,i)*d3trm*dr2ikds(ks)
              derv3(iy,kl) = derv3(iy,kl) - ynbr(ni2,i)*d3trm*dr2ikds(ks)
              derv3(iz,kl) = derv3(iz,kl) - znbr(ni2,i)*d3trm*dr2ikds(ks)
            enddo
!
            do kk = 1,nstrains
              ks = nstrptr(kk)
              do kl = 1,nstrains
                kt = nstrptr(kl)
                sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ilds(ks) + dr2ilds(kt)*dr2ikds(ks))
              enddo
            enddo
          endif
!
          if (lklocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ilds(ks)
              derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ilds(ks)
              derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ilds(ks)
            enddo
          endif
!
          if (lllocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(lx,kl) = derv3(lx,kl) + xnbr(ni2,i)*d3trm*dr2ikds(ks)
              derv3(ly,kl) = derv3(ly,kl) + ynbr(ni2,i)*d3trm*dr2ikds(ks)
              derv3(lz,kl) = derv3(lz,kl) + znbr(ni2,i)*d3trm*dr2ikds(ks)
            enddo
          endif
        endif
      enddo
!
      d3trm1 = d2Edlogcnidlogcnj*dlogcnidcni*dlogcnjdcnj*dcnidrik
!------------------------------------------------------------------------------------------
!  Loop over neighbours of j for second derivatives for 2 different coordination numbers  |
!------------------------------------------------------------------------------------------
      do n2 = 1,nnbr_cn(j)
        dcnjdrjl = d1cndr_cn(n2,j)
        d3trm = d3trm1*dcnjdrjl
!
        if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
        nj = nbrno_cn(n2,j)
        l = nbrno(nj,j)
        lloc = atom2local(l)
        lllocal = (lloc.ne.0)
!
!  If none of the atoms are local then cycle
!
        if (.not.lijklocal.and..not.lllocal) cycle
!
        indl = 3*(l-1)
        lxf = indl + 1
        lyf = indl + 2
        lzf = indl + 3
!
        if (lllocal) then
          indl = 3*(lloc-1)
          lx = indl + 1
          ly = indl + 2
          lz = indl + 3
        endif
!
        if (lstr) then
          call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jlds,d2r2jldx2, &
                            d2r2jldsdx,d2r2jlds2,.false.)
        endif
!
!  I-J
!
        if (lilocal) then
          derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*znbr(ni,i)
        endif
        if (ljlocal) then
          derv2(ixf,jx) = derv2(ixf,jx) + d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(iyf,jx) = derv2(iyf,jx) + d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(izf,jx) = derv2(izf,jx) + d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(ixf,jy) = derv2(ixf,jy) + d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(iyf,jy) = derv2(iyf,jy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(izf,jy) = derv2(izf,jy) + d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(ixf,jz) = derv2(ixf,jz) + d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(iyf,jz) = derv2(iyf,jz) + d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(izf,jz) = derv2(izf,jz) + d3trm*znbr(nj,j)*znbr(ni,i)
        endif
!
!  I-L
!
        if (lilocal) then
          derv2(lxf,ix) = derv2(lxf,ix) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(lyf,ix) = derv2(lyf,ix) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(lzf,ix) = derv2(lzf,ix) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(lxf,iy) = derv2(lxf,iy) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(lyf,iy) = derv2(lyf,iy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(lzf,iy) = derv2(lzf,iy) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(lxf,iz) = derv2(lxf,iz) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(lyf,iz) = derv2(lyf,iz) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(lzf,iz) = derv2(lzf,iz) - d3trm*znbr(nj,j)*znbr(ni,i)
        endif
        if (lllocal) then
          derv2(ixf,lx) = derv2(ixf,lx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(iyf,lx) = derv2(iyf,lx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(izf,lx) = derv2(izf,lx) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(ixf,ly) = derv2(ixf,ly) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(iyf,ly) = derv2(iyf,ly) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(izf,ly) = derv2(izf,ly) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(ixf,lz) = derv2(ixf,lz) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(iyf,lz) = derv2(iyf,lz) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(izf,lz) = derv2(izf,lz) - d3trm*znbr(nj,j)*znbr(ni,i)
        endif
!
!  J-K
!
        if (ljlocal) then
          derv2(kxf,jx) = derv2(kxf,jx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(kyf,jx) = derv2(kyf,jx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(kzf,jx) = derv2(kzf,jx) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(kxf,jy) = derv2(kxf,jy) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(kyf,jy) = derv2(kyf,jy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(kzf,jy) = derv2(kzf,jy) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(kxf,jz) = derv2(kxf,jz) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(kyf,jz) = derv2(kyf,jz) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(kzf,jz) = derv2(kzf,jz) - d3trm*znbr(nj,j)*znbr(ni,i)
        endif
        if (lklocal) then
          derv2(jxf,kx) = derv2(jxf,kx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
          derv2(jyf,kx) = derv2(jyf,kx) - d3trm*ynbr(nj,j)*xnbr(ni,i)
          derv2(jzf,kx) = derv2(jzf,kx) - d3trm*znbr(nj,j)*xnbr(ni,i)
          derv2(jxf,ky) = derv2(jxf,ky) - d3trm*xnbr(nj,j)*ynbr(ni,i)
          derv2(jyf,ky) = derv2(jyf,ky) - d3trm*ynbr(nj,j)*ynbr(ni,i)
          derv2(jzf,ky) = derv2(jzf,ky) - d3trm*znbr(nj,j)*ynbr(ni,i)
          derv2(jxf,kz) = derv2(jxf,kz) - d3trm*xnbr(nj,j)*znbr(ni,i)
          derv2(jyf,kz) = derv2(jyf,kz) - d3trm*ynbr(nj,j)*znbr(ni,i)
          derv2(jzf,kz) = derv2(jzf,kz) - d3trm*znbr(nj,j)*znbr(ni,i)
        endif
!
!  K-L
!     
        if (lklocal) then
          derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
          derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(ni,i)*ynbr(nj,j)
          derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(ni,i)*znbr(nj,j)
          derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(ni,i)*xnbr(nj,j)
          derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(ni,i)*ynbr(nj,j)
          derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(ni,i)*znbr(nj,j)
          derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(ni,i)*xnbr(nj,j)
          derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(ni,i)*ynbr(nj,j)
          derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(ni,i)*znbr(nj,j)
        endif
        if (lllocal) then
          derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
          derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(ni,i)*xnbr(nj,j)
          derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(ni,i)*xnbr(nj,j)
          derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(ni,i)*ynbr(nj,j)
          derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(ni,i)*ynbr(nj,j)
          derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(ni,i)*ynbr(nj,j)
          derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(ni,i)*znbr(nj,j)
          derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(ni,i)*znbr(nj,j)
          derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(ni,i)*znbr(nj,j)
        endif
!
!  Strain terms
!
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2jlds(ks)
              derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2jlds(ks)
              derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2jlds(ks)
            enddo
!
            do kk = 1,nstrains
              ks = nstrptr(kk)
              do kl = 1,nstrains
                kt = nstrptr(kl)
                sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2ikds(ks))
              enddo
            enddo
          endif
!
          if (lklocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2jlds(ks)
              derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2jlds(ks)
              derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2jlds(ks)
            enddo
          endif
!
          if (ljlocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2ikds(ks)
              derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2ikds(ks)
              derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2ikds(ks)
            enddo
          endif
!
          if (lllocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(lx,kl) = derv3(lx,kl) + xnbr(nj,j)*d3trm*dr2ikds(ks)
              derv3(ly,kl) = derv3(ly,kl) + ynbr(nj,j)*d3trm*dr2ikds(ks)
              derv3(lz,kl) = derv3(lz,kl) + znbr(nj,j)*d3trm*dr2ikds(ks)
            enddo
          endif
        endif
      enddo
    endif
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    dcnjdrjk = d1cndr_cn(n,j)
    d1trm = dlogcnjdcnj*dcnjdrjk*dEdlogcnj
    dxyz(1) = d1trm*xnbr(nj,j)
    dxyz(2) = d1trm*ynbr(nj,j)
    dxyz(3) = d1trm*znbr(nj,j)
!
    if (lilocal) then
      xdrv(j) = xdrv(j) - dxyz(1)
      ydrv(j) = ydrv(j) - dxyz(2)
      zdrv(j) = zdrv(j) - dxyz(3)
      xdrv(k) = xdrv(k) + dxyz(1)
      ydrv(k) = ydrv(k) + dxyz(2)
      zdrv(k) = zdrv(k) + dxyz(3)
    endif
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jkds,d2r2jkdx2, &
                        d2r2jkdsdx,d2r2jkds2,lgrad2)
    endif
    if (lstr) then
      if (lilocal) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + d1trm*dr2jkds(ks)
        enddo
      endif
    endif
    if (lgrad2) then
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
      d2cnjdrjk2 = d2cndr_cn(n,j)
      d2trm = dEdlogcnj*dlogcnjdcnj*d2cnjdrjk2 + &
             (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjk
!
      if (abs(d2trm).gt.gfnff_cnc6tol) then
!----------------------------------------
!  Second derivatives of CN w.r.t. rjk  |
!----------------------------------------
        if (ljlocal) then
          derv2(kxf,jx) = derv2(kxf,jx) - d2trm*d2r2jkdx2(1)
          derv2(kyf,jx) = derv2(kyf,jx) - d2trm*d2r2jkdx2(6)
          derv2(kzf,jx) = derv2(kzf,jx) - d2trm*d2r2jkdx2(5)
          derv2(kxf,jy) = derv2(kxf,jy) - d2trm*d2r2jkdx2(6)
          derv2(kyf,jy) = derv2(kyf,jy) - d2trm*d2r2jkdx2(2)
          derv2(kzf,jy) = derv2(kzf,jy) - d2trm*d2r2jkdx2(4)
          derv2(kxf,jz) = derv2(kxf,jz) - d2trm*d2r2jkdx2(5)
          derv2(kyf,jz) = derv2(kyf,jz) - d2trm*d2r2jkdx2(4)
          derv2(kzf,jz) = derv2(kzf,jz) - d2trm*d2r2jkdx2(3)
          derv2(kxf,jx) = derv2(kxf,jx) - d1trm
          derv2(kyf,jy) = derv2(kyf,jy) - d1trm
          derv2(kzf,jz) = derv2(kzf,jz) - d1trm
        endif
        if (lklocal) then
          derv2(jxf,kx) = derv2(jxf,kx) - d2trm*d2r2jkdx2(1)
          derv2(jyf,kx) = derv2(jyf,kx) - d2trm*d2r2jkdx2(6)
          derv2(jzf,kx) = derv2(jzf,kx) - d2trm*d2r2jkdx2(5)
          derv2(jxf,ky) = derv2(jxf,ky) - d2trm*d2r2jkdx2(6)
          derv2(jyf,ky) = derv2(jyf,ky) - d2trm*d2r2jkdx2(2)
          derv2(jzf,ky) = derv2(jzf,ky) - d2trm*d2r2jkdx2(4)
          derv2(jxf,kz) = derv2(jxf,kz) - d2trm*d2r2jkdx2(5)
          derv2(jyf,kz) = derv2(jyf,kz) - d2trm*d2r2jkdx2(4)
          derv2(jzf,kz) = derv2(jzf,kz) - d2trm*d2r2jkdx2(3)
          derv2(jxf,kx) = derv2(jxf,kx) - d1trm
          derv2(jyf,ky) = derv2(jyf,ky) - d1trm
          derv2(jzf,kz) = derv2(jzf,kz) - d1trm
        endif
!
        if (lstr) then
!
!  Mixed derivatives
!
          if (ljlocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(jx,kl) = derv3(jx,kl) - d1trm*d2r2jkdsdx(ks,1) - xnbr(nj,j)*d2trm*dr2jkds(ks)
              derv3(jy,kl) = derv3(jy,kl) - d1trm*d2r2jkdsdx(ks,2) - ynbr(nj,j)*d2trm*dr2jkds(ks)
              derv3(jz,kl) = derv3(jz,kl) - d1trm*d2r2jkdsdx(ks,3) - znbr(nj,j)*d2trm*dr2jkds(ks)
            enddo
          endif
          if (lklocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2jkdsdx(ks,1) + xnbr(nj,j)*d2trm*dr2jkds(ks)
              derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2jkdsdx(ks,2) + ynbr(nj,j)*d2trm*dr2jkds(ks)
              derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2jkdsdx(ks,3) + znbr(nj,j)*d2trm*dr2jkds(ks)
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
                sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2jkds(kt)*dr2jkds(ks) + d1trm*d2r2jkds2(kt,ks)
              enddo
            enddo
          endif
        endif
      endif
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
      d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
      if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  J-K
!
        if (ljlocal) then
          derv2(kxf,jx) = derv2(kxf,jx) + d3trm*xnbr(nj,j)*xij
          derv2(kyf,jx) = derv2(kyf,jx) + d3trm*ynbr(nj,j)*xij
          derv2(kzf,jx) = derv2(kzf,jx) + d3trm*znbr(nj,j)*xij
          derv2(kxf,jy) = derv2(kxf,jy) + d3trm*xnbr(nj,j)*yij
          derv2(kyf,jy) = derv2(kyf,jy) + d3trm*ynbr(nj,j)*yij
          derv2(kzf,jy) = derv2(kzf,jy) + d3trm*znbr(nj,j)*yij
          derv2(kxf,jz) = derv2(kxf,jz) + d3trm*xnbr(nj,j)*zij
          derv2(kyf,jz) = derv2(kyf,jz) + d3trm*ynbr(nj,j)*zij
          derv2(kzf,jz) = derv2(kzf,jz) + d3trm*znbr(nj,j)*zij
        endif
        if (lklocal) then
          derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(nj,j)*xij
          derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(nj,j)*yij
          derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(nj,j)*zij
          derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(nj,j)*xij
          derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(nj,j)*yij
          derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(nj,j)*zij
          derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(nj,j)*xij
          derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(nj,j)*yij
          derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(nj,j)*zij
        endif
!
!  I-K
!
        if (lilocal) then
          derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(nj,j)*xij
          derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(nj,j)*xij
          derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(nj,j)*xij
          derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(nj,j)*yij
          derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(nj,j)*yij
          derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(nj,j)*yij
          derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(nj,j)*zij
          derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(nj,j)*zij
          derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(nj,j)*zij
        endif
        if (lklocal) then
          derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(nj,j)*xij
          derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(nj,j)*yij
          derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(nj,j)*zij
          derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(nj,j)*xij
          derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(nj,j)*yij
          derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(nj,j)*zij
          derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(nj,j)*xij
          derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(nj,j)*yij
          derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(nj,j)*zij
        endif
!
!  J-I
!
        if (ljlocal) then
          derv2(ixf,jx) = derv2(ixf,jx) + d3trm*xnbr(nj,j)*xij
          derv2(iyf,jx) = derv2(iyf,jx) + d3trm*xnbr(nj,j)*yij
          derv2(izf,jx) = derv2(izf,jx) + d3trm*xnbr(nj,j)*zij
          derv2(ixf,jy) = derv2(ixf,jy) + d3trm*ynbr(nj,j)*xij
          derv2(iyf,jy) = derv2(iyf,jy) + d3trm*ynbr(nj,j)*yij
          derv2(izf,jy) = derv2(izf,jy) + d3trm*ynbr(nj,j)*zij
          derv2(ixf,jz) = derv2(ixf,jz) + d3trm*znbr(nj,j)*xij
          derv2(iyf,jz) = derv2(iyf,jz) + d3trm*znbr(nj,j)*yij
          derv2(izf,jz) = derv2(izf,jz) + d3trm*znbr(nj,j)*zij
        endif
        if (lilocal) then
          derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xij
          derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xij
          derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xij
          derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*yij
          derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*yij
          derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*yij
          derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*zij
          derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*zij
          derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*zij
        endif
!
!  Strain terms
!
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2jkds(ks)
              derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2jkds(ks)
              derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2jkds(ks)
            enddo
!
            do kk = 1,nstrains
              ks = nstrptr(kk)
              do kl = 1,nstrains
                kt = nstrptr(kl)
                sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2jkds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2jkds(ks))
              enddo
            enddo
          endif
!
          if (ljlocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2ijds(ks)
              derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2ijds(ks)
              derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2ijds(ks)
              derv3(jx,kl) = derv3(jx,kl) + xij*d3trm*dr2jkds(ks)
              derv3(jy,kl) = derv3(jy,kl) + yij*d3trm*dr2jkds(ks)
              derv3(jz,kl) = derv3(jz,kl) + zij*d3trm*dr2jkds(ks)
            enddo
          endif
!
          if (lklocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(kx,kl) = derv3(kx,kl) + xnbr(nj,j)*d3trm*dr2ijds(ks)
              derv3(ky,kl) = derv3(ky,kl) + ynbr(nj,j)*d3trm*dr2ijds(ks)
              derv3(kz,kl) = derv3(kz,kl) + znbr(nj,j)*d3trm*dr2ijds(ks)
            enddo
          endif
        endif
      endif
!
      d3trm1 = (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of j for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
      do n2 = 1,n-1
        dcnjdrjl = d1cndr_cn(n2,j)
        d3trm = d3trm1*dcnjdrjl
!
        if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
        nj2 = nbrno_cn(n2,j)
        l = nbrno(nj2,j)
        lloc = atom2local(l)
        lllocal = (lloc.ne.0)
!
!  If none of the atoms are local then cycle
!
        if (.not.lijklocal.and..not.lllocal) cycle
!
        indl = 3*(l-1)
        lxf = indl + 1
        lyf = indl + 2
        lzf = indl + 3
!
        if (lllocal) then
          indl = 3*(lloc-1)
          lx = indl + 1
          ly = indl + 2
          lz = indl + 3
        endif
!
        if (lstr) then
          call real1strterm(ndim,xnbr(nj2,j),ynbr(nj2,j),znbr(nj2,j),0.0_dp,0.0_dp,0.0_dp,dr2jlds,d2r2jldx2, &
                            d2r2jldsdx,d2r2jlds2,.false.)
        endif
!
!  J-K
!
        if (ljlocal) then
          derv2(kxf,jx) = derv2(kxf,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(kyf,jx) = derv2(kyf,jx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(kzf,jx) = derv2(kzf,jx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(kxf,jy) = derv2(kxf,jy) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(kyf,jy) = derv2(kyf,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(kzf,jy) = derv2(kzf,jy) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(kxf,jz) = derv2(kxf,jz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(kyf,jz) = derv2(kyf,jz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(kzf,jz) = derv2(kzf,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
        if (lklocal) then
          derv2(jxf,kx) = derv2(jxf,kx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(jyf,kx) = derv2(jyf,kx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(jzf,kx) = derv2(jzf,kx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(jxf,ky) = derv2(jxf,ky) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(jyf,ky) = derv2(jyf,ky) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(jzf,ky) = derv2(jzf,ky) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(jxf,kz) = derv2(jxf,kz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(jyf,kz) = derv2(jyf,kz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(jzf,kz) = derv2(jzf,kz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
!
!  J-L
!
        if (ljlocal) then
          derv2(lxf,jx) = derv2(lxf,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(lyf,jx) = derv2(lyf,jx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(lzf,jx) = derv2(lzf,jx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(lxf,jy) = derv2(lxf,jy) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(lyf,jy) = derv2(lyf,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(lzf,jy) = derv2(lzf,jy) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(lxf,jz) = derv2(lxf,jz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(lyf,jz) = derv2(lyf,jz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(lzf,jz) = derv2(lzf,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
        if (lllocal) then
          derv2(jxf,lx) = derv2(jxf,lx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(jyf,lx) = derv2(jyf,lx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(jzf,lx) = derv2(jzf,lx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(jxf,ly) = derv2(jxf,ly) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(jyf,ly) = derv2(jyf,ly) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(jzf,ly) = derv2(jzf,ly) - d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(jxf,lz) = derv2(jxf,lz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(jyf,lz) = derv2(jyf,lz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(jzf,lz) = derv2(jzf,lz) - d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
!
!  K-L
!     
        if (lklocal) then
          derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
        if (lllocal) then
          derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
          derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
          derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(nj,j)*xnbr(nj2,j)
          derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
          derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
          derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(nj,j)*ynbr(nj2,j)
          derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(nj,j)*znbr(nj2,j)
          derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(nj,j)*znbr(nj2,j)
          derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(nj,j)*znbr(nj2,j)
        endif
!
!  Strain terms
!
        if (lstr) then
          if (ljlocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2jlds(ks)
              derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2jlds(ks)
              derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2jlds(ks)
              derv3(jx,kl) = derv3(jx,kl) - xnbr(nj2,j)*d3trm*dr2jkds(ks)
              derv3(jy,kl) = derv3(jy,kl) - ynbr(nj2,j)*d3trm*dr2jkds(ks)
              derv3(jz,kl) = derv3(jz,kl) - znbr(nj2,j)*d3trm*dr2jkds(ks)
            enddo
          endif
!
          if (lklocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(kx,kl) = derv3(kx,kl) + xnbr(nj,j)*d3trm*dr2jlds(ks)
              derv3(ky,kl) = derv3(ky,kl) + ynbr(nj,j)*d3trm*dr2jlds(ks)
              derv3(kz,kl) = derv3(kz,kl) + znbr(nj,j)*d3trm*dr2jlds(ks)
            enddo
          endif
!
          if (lllocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(lx,kl) = derv3(lx,kl) + xnbr(nj2,j)*d3trm*dr2jkds(ks)
              derv3(ly,kl) = derv3(ly,kl) + ynbr(nj2,j)*d3trm*dr2jkds(ks)
              derv3(lz,kl) = derv3(lz,kl) + znbr(nj2,j)*d3trm*dr2jkds(ks)
            enddo
          endif
!
          if (lilocal) then
            do kk = 1,nstrains
              ks = nstrptr(kk)
              do kl = 1,nstrains
                kt = nstrptr(kl)
                sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2jkds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2jkds(ks))
              enddo
            enddo
          endif
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnd')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnd_x1(i,ix,iy,iz,j,jxf,jyf,jzf,xij,yij,zij, &
                                d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn,ofctij)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Distributed memory parallel version. This version only does crossterms between distance/coordination 
!  Version where i is local. Only called if lgrad2 is true
!
!  On entry : 
!
!  ix,  iy,  iz      = second derivative Cartesian elements for i (local)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  ofctij            = scale factor certain terms
!
!   5/22 Created from gfnff_drv2_dcnd_x
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use parallel
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: ofctij
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnd_x1')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  if (lstr) then
    call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2,d2r2ijdsdx,d2r2ijds2,.false.)
  endif
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
!
    call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                      d2r2ikdsdx,d2r2ikds2,.true.)
    indk = 3*(k-1)
    kxf = indk + 1
    kyf = indk + 2
    kzf = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
    d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
    if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  I-K
!
      derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xij
      derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xij
      derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xij
      derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*yij
      derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*yij
      derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*yij
      derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*zij
      derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*zij
      derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*zij
!
!  I-J
!
      derv2(jxf,ix) = derv2(jxf,ix) - d3trm*xnbr(ni,i)*xij
      derv2(jyf,ix) = derv2(jyf,ix) - d3trm*xnbr(ni,i)*yij
      derv2(jzf,ix) = derv2(jzf,ix) - d3trm*xnbr(ni,i)*zij
      derv2(jxf,iy) = derv2(jxf,iy) - d3trm*ynbr(ni,i)*xij
      derv2(jyf,iy) = derv2(jyf,iy) - d3trm*ynbr(ni,i)*yij
      derv2(jzf,iy) = derv2(jzf,iy) - d3trm*ynbr(ni,i)*zij
      derv2(jxf,iz) = derv2(jxf,iz) - d3trm*znbr(ni,i)*xij
      derv2(jyf,iz) = derv2(jyf,iz) - d3trm*znbr(ni,i)*yij
      derv2(jzf,iz) = derv2(jzf,iz) - d3trm*znbr(ni,i)*zij
!
!  Strain terms
!
      if (lstr) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ijds(ks)
          derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ijds(ks)
          derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ijds(ks)
          derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2ikds(ks)
          derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2ikds(ks)
          derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2ikds(ks)
        enddo
!
        do kk = 1,nstrains
          ks = nstrptr(kk)
          do kl = 1,nstrains
            kt = nstrptr(kl)
            sderv2(kl,kk) = sderv2(kl,kk) + ofctij*d3trm*(dr2ikds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2ikds(ks)) 
          enddo
        enddo
      endif
    endif
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
!
    dcnjdrjk = d1cndr_cn(n,j)
!
    call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jkds,d2r2jkdx2, &
                      d2r2jkdsdx,d2r2jkds2,.true.)
!
    indk = 3*(k-1)
    kxf = indk + 1
    kyf = indk + 2
    kzf = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
    d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
    if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  I-K
!
      derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(nj,j)*xij
      derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(nj,j)*xij
      derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(nj,j)*xij
      derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(nj,j)*yij
      derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(nj,j)*yij
      derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(nj,j)*yij
      derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(nj,j)*zij
      derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(nj,j)*zij
      derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(nj,j)*zij
!
!  J-I
!
      derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xij
      derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xij
      derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xij
      derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*yij
      derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*yij
      derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*yij
      derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*zij
      derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*zij
      derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*zij
!
!  Strain terms
!
      if (lstr) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2jkds(ks)
          derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2jkds(ks)
          derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2jkds(ks)
        enddo
!
        do kk = 1,nstrains
          ks = nstrptr(kk)
          do kl = 1,nstrains
            kt = nstrptr(kl)
            sderv2(kl,kk) = sderv2(kl,kk) + ofctij*d3trm*(dr2jkds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2jkds(ks))
          enddo
        enddo
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnd_x1')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnd_x2(i,ixf,iyf,izf,j,jxf,jyf,jzf,xij,yij,zij, &
                                d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Distributed memory parallel version. This version only does crossterms between distance/coordination 
!  Version where i and j are global and k is local. Only called when lgrad2 is true.
!
!  On entry : 
!
!  ixf, iyf, izf     = second derivative Cartesian elements for i (global)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!
!   5/22 Created from gfnff_drv2_dcnd_x
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use parallel
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: k
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: ks
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  logical                                          :: lklocal
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnd_x2')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  if (lstr) then
    call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2,d2r2ijdsdx,d2r2ijds2,.false.)
  endif
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
!
    if (.not.lklocal) cycle
!
    dcnidrik = d1cndr_cn(n,i)
!
    call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                      d2r2ikdsdx,d2r2ikds2,.true.)
!
    indk = 3*(kloc-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
    d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
    if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  I-K
!
      derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xij
      derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*yij
      derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*zij
      derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xij
      derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*yij
      derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*zij
      derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xij
      derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*yij
      derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*zij
!
!  J-K
!
      derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(ni,i)*xij
      derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(ni,i)*yij
      derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(ni,i)*zij
      derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(ni,i)*xij
      derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(ni,i)*yij
      derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(ni,i)*zij
      derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(ni,i)*xij
      derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(ni,i)*yij
      derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(ni,i)*zij
!
!  Strain terms
!
      if (lstr) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ijds(ks)
          derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ijds(ks)
          derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ijds(ks)
        enddo
      endif
    endif
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
!
    if (.not.lklocal) cycle
!
    dcnjdrjk = d1cndr_cn(n,j)
!
    call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jkds,d2r2jkdx2, &
                      d2r2jkdsdx,d2r2jkds2,.true.)
!
    indk = 3*(kloc-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
    d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
    if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  J-K
!
      derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(nj,j)*xij
      derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(nj,j)*yij
      derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(nj,j)*zij
      derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(nj,j)*xij
      derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(nj,j)*yij
      derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(nj,j)*zij
      derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(nj,j)*xij
      derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(nj,j)*yij
      derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(nj,j)*zij
!
!  I-K
!
      derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(nj,j)*xij
      derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(nj,j)*yij
      derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(nj,j)*zij
      derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(nj,j)*xij
      derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(nj,j)*yij
      derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(nj,j)*zij
      derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(nj,j)*xij
      derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(nj,j)*yij
      derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(nj,j)*zij
!
!  Strain terms
!
      if (lstr) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(kx,kl) = derv3(kx,kl) + xnbr(nj,j)*d3trm*dr2ijds(ks)
          derv3(ky,kl) = derv3(ky,kl) + ynbr(nj,j)*d3trm*dr2ijds(ks)
          derv3(kz,kl) = derv3(kz,kl) + znbr(nj,j)*d3trm*dr2ijds(ks)
        enddo
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnd_x2')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnd_x(lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,j,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                               d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Distributed memory parallel version. This version only does crossterms between distance/coordination 
!  Only called if lgrad2 is true
!
!  On entry : 
!
!  lilocal           = if true then i is local to this processor
!  ljlocal           = if true then j is local to this processor
!  ix,  iy,  iz      = second derivative Cartesian elements for i (local)
!  ixf, iyf, izf     = second derivative Cartesian elements for i (global)
!  jx,  jy,  jz      = second derivative Cartesian elements for j (local)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!
!   5/22 Created from gfnff_drv2_dcnd 
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use parallel
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: indk
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
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  logical                                          :: liklocal
  logical                                          :: lijklocal
  logical                                          :: lklocal
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnd_x')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  if (lstr) then
    call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2,d2r2ijdsdx,d2r2ijds2,.false.)
  endif
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
!
    call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                      d2r2ikdsdx,d2r2ikds2,.true.)
!
    indk = 3*(k-1)
    kxf = indk + 1
    kyf = indk + 2
    kzf = indk + 3
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    liklocal = (lilocal.or.lklocal)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    if (lklocal) then
      indk = 3*(kloc-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
    endif
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
    d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
    if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  I-K
!
      if (lilocal) then
        derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xij
        derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xij
        derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xij
        derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*yij
        derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*yij
        derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*yij
        derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*zij
        derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*zij
        derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*zij
      endif
      if (lklocal) then
        derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xij
        derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*yij
        derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*zij
        derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xij
        derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*yij
        derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*zij
        derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xij
        derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*yij
        derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*zij
      endif
!
!  J-K
!
      if (ljlocal) then
        derv2(kxf,jx) = derv2(kxf,jx) + d3trm*xnbr(ni,i)*xij
        derv2(kyf,jx) = derv2(kyf,jx) + d3trm*ynbr(ni,i)*xij
        derv2(kzf,jx) = derv2(kzf,jx) + d3trm*znbr(ni,i)*xij
        derv2(kxf,jy) = derv2(kxf,jy) + d3trm*xnbr(ni,i)*yij
        derv2(kyf,jy) = derv2(kyf,jy) + d3trm*ynbr(ni,i)*yij
        derv2(kzf,jy) = derv2(kzf,jy) + d3trm*znbr(ni,i)*yij
        derv2(kxf,jz) = derv2(kxf,jz) + d3trm*xnbr(ni,i)*zij
        derv2(kyf,jz) = derv2(kyf,jz) + d3trm*ynbr(ni,i)*zij
        derv2(kzf,jz) = derv2(kzf,jz) + d3trm*znbr(ni,i)*zij
      endif
      if (lklocal) then
        derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(ni,i)*xij
        derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(ni,i)*yij
        derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(ni,i)*zij
        derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(ni,i)*xij
        derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(ni,i)*yij
        derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(ni,i)*zij
        derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(ni,i)*xij
        derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(ni,i)*yij
        derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(ni,i)*zij
      endif
!
!  I-J
!
      if (lilocal) then
        derv2(jxf,ix) = derv2(jxf,ix) - d3trm*xnbr(ni,i)*xij
        derv2(jyf,ix) = derv2(jyf,ix) - d3trm*xnbr(ni,i)*yij
        derv2(jzf,ix) = derv2(jzf,ix) - d3trm*xnbr(ni,i)*zij
        derv2(jxf,iy) = derv2(jxf,iy) - d3trm*ynbr(ni,i)*xij
        derv2(jyf,iy) = derv2(jyf,iy) - d3trm*ynbr(ni,i)*yij
        derv2(jzf,iy) = derv2(jzf,iy) - d3trm*ynbr(ni,i)*zij
        derv2(jxf,iz) = derv2(jxf,iz) - d3trm*znbr(ni,i)*xij
        derv2(jyf,iz) = derv2(jyf,iz) - d3trm*znbr(ni,i)*yij
        derv2(jzf,iz) = derv2(jzf,iz) - d3trm*znbr(ni,i)*zij
      endif
      if (ljlocal) then
        derv2(ixf,jx) = derv2(ixf,jx) - d3trm*xnbr(ni,i)*xij
        derv2(iyf,jx) = derv2(iyf,jx) - d3trm*ynbr(ni,i)*xij
        derv2(izf,jx) = derv2(izf,jx) - d3trm*znbr(ni,i)*xij
        derv2(ixf,jy) = derv2(ixf,jy) - d3trm*xnbr(ni,i)*yij
        derv2(iyf,jy) = derv2(iyf,jy) - d3trm*ynbr(ni,i)*yij
        derv2(izf,jy) = derv2(izf,jy) - d3trm*znbr(ni,i)*yij
        derv2(ixf,jz) = derv2(ixf,jz) - d3trm*xnbr(ni,i)*zij
        derv2(iyf,jz) = derv2(iyf,jz) - d3trm*ynbr(ni,i)*zij
        derv2(izf,jz) = derv2(izf,jz) - d3trm*znbr(ni,i)*zij
      endif
!
!  Strain terms
!
      if (lstr) then
        if (lilocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2ikds(ks)
            derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2ikds(ks)
            derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2ikds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2ikds(ks)) 
            enddo
          enddo
        endif
!
        if (ljlocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) + xij*d3trm*dr2ikds(ks)
            derv3(jy,kl) = derv3(jy,kl) + yij*d3trm*dr2ikds(ks)
            derv3(jz,kl) = derv3(jz,kl) + zij*d3trm*dr2ikds(ks)
          enddo
        endif
!
        if (lklocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ijds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ijds(ks)
          enddo
        endif
      endif
    endif
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    dcnjdrjk = d1cndr_cn(n,j)
!
    call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jkds,d2r2jkdx2, &
                      d2r2jkdsdx,d2r2jkds2,.true.)
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
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
    d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
    if (abs(d3trm).gt.gfnff_cnc6tol) then
!
!  J-K
!
      if (ljlocal) then
        derv2(kxf,jx) = derv2(kxf,jx) + d3trm*xnbr(nj,j)*xij
        derv2(kyf,jx) = derv2(kyf,jx) + d3trm*ynbr(nj,j)*xij
        derv2(kzf,jx) = derv2(kzf,jx) + d3trm*znbr(nj,j)*xij
        derv2(kxf,jy) = derv2(kxf,jy) + d3trm*xnbr(nj,j)*yij
        derv2(kyf,jy) = derv2(kyf,jy) + d3trm*ynbr(nj,j)*yij
        derv2(kzf,jy) = derv2(kzf,jy) + d3trm*znbr(nj,j)*yij
        derv2(kxf,jz) = derv2(kxf,jz) + d3trm*xnbr(nj,j)*zij
        derv2(kyf,jz) = derv2(kyf,jz) + d3trm*ynbr(nj,j)*zij
        derv2(kzf,jz) = derv2(kzf,jz) + d3trm*znbr(nj,j)*zij
      endif
      if (lklocal) then
        derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(nj,j)*xij
        derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(nj,j)*yij
        derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(nj,j)*zij
        derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(nj,j)*xij
        derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(nj,j)*yij
        derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(nj,j)*zij
        derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(nj,j)*xij
        derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(nj,j)*yij
        derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(nj,j)*zij
      endif
!
!  I-K
!
      if (lilocal) then
        derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(nj,j)*xij
        derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(nj,j)*xij
        derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(nj,j)*xij
        derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(nj,j)*yij
        derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(nj,j)*yij
        derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(nj,j)*yij
        derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(nj,j)*zij
        derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(nj,j)*zij
        derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(nj,j)*zij
      endif
      if (lklocal) then
        derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(nj,j)*xij
        derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(nj,j)*yij
        derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(nj,j)*zij
        derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(nj,j)*xij
        derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(nj,j)*yij
        derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(nj,j)*zij
        derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(nj,j)*xij
        derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(nj,j)*yij
        derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(nj,j)*zij
      endif
!
!  J-I
!
      if (ljlocal) then
        derv2(ixf,jx) = derv2(ixf,jx) + d3trm*xnbr(nj,j)*xij
        derv2(iyf,jx) = derv2(iyf,jx) + d3trm*xnbr(nj,j)*yij
        derv2(izf,jx) = derv2(izf,jx) + d3trm*xnbr(nj,j)*zij
        derv2(ixf,jy) = derv2(ixf,jy) + d3trm*ynbr(nj,j)*xij
        derv2(iyf,jy) = derv2(iyf,jy) + d3trm*ynbr(nj,j)*yij
        derv2(izf,jy) = derv2(izf,jy) + d3trm*ynbr(nj,j)*zij
        derv2(ixf,jz) = derv2(ixf,jz) + d3trm*znbr(nj,j)*xij
        derv2(iyf,jz) = derv2(iyf,jz) + d3trm*znbr(nj,j)*yij
        derv2(izf,jz) = derv2(izf,jz) + d3trm*znbr(nj,j)*zij
      endif
      if (lilocal) then
        derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xij
        derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xij
        derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xij
        derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*yij
        derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*yij
        derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*yij
        derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*zij
        derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*zij
        derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*zij
      endif
!
!  Strain terms
!
      if (lstr) then
        if (lilocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xij*d3trm*dr2jkds(ks)
            derv3(iy,kl) = derv3(iy,kl) - yij*d3trm*dr2jkds(ks)
            derv3(iz,kl) = derv3(iz,kl) - zij*d3trm*dr2jkds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2jkds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2jkds(ks))
            enddo
          enddo
        endif
!
        if (ljlocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(jx,kl) = derv3(jx,kl) + xij*d3trm*dr2jkds(ks)
            derv3(jy,kl) = derv3(jy,kl) + yij*d3trm*dr2jkds(ks)
            derv3(jz,kl) = derv3(jz,kl) + zij*d3trm*dr2jkds(ks)
          enddo
        endif
!
        if (lklocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(nj,j)*d3trm*dr2ijds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(nj,j)*d3trm*dr2ijds(ks)
          enddo
        endif
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnd_x')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnd_nox(lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,j,jx,jy,jz,jxf,jyf,jzf, &
                                 dEdlogcni,dEdlogcnj,d2Edlogcni2,d2Edlogcnj2,d2Edlogcnidlogcnj, &
                                 dlogcndcn,d2logcndcn2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Distributed memory parallel version. This version excludes crossterms between distance/coordination 
!  Only called if lgrad2 is true
!
!  On entry : 
!
!  lilocal           = if true then i is local to this processor
!  ljlocal           = if true then j is local to this processor
!  ix,  iy,  iz      = second derivative Cartesian elements for i (local)
!  ixf, iyf, izf     = second derivative Cartesian elements for i (global)
!  jx,  jy,  jz      = second derivative Cartesian elements for j (local)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  dEdlogcnj         = derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcnj2       = second derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcnidlogcnj = second derivatives of energy w.r.t. log of the coordination number for i and j
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!
!   5/22 Created from gfnff_drv2_dcnd 
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use parallel
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: dEdlogcnj
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: d2Edlogcnidlogcnj
  real(dp),    intent(in)                          :: d2Edlogcnj2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
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
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  logical                                          :: liklocal
  logical                                          :: lijklocal
  logical                                          :: lklocal
  logical                                          :: lllocal
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: d3trm1
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dcnjdrjl
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2logcnjdcnj2
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2ilds(6)
  real(dp)                                         :: d2r2ildx2(6)
  real(dp)                                         :: d2r2ilds2(6,6)
  real(dp)                                         :: d2r2ildsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
  real(dp)                                         :: dr2jlds(6)
  real(dp)                                         :: d2r2jldx2(6)
  real(dp)                                         :: d2r2jlds2(6,6)
  real(dp)                                         :: d2r2jldsdx(6,3)
  real(dp)                                         :: dxyz(3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnd_nox')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  d2logcnidcni2 = d2logcndcn2(i)
  d2logcnjdcnj2 = d2logcndcn2(j)
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
    dxyz(1) = d1trm*xnbr(ni,i)
    dxyz(2) = d1trm*ynbr(ni,i)
    dxyz(3) = d1trm*znbr(ni,i)
!
    if (lilocal) then
      xdrv(i) = xdrv(i) - dxyz(1)
      ydrv(i) = ydrv(i) - dxyz(2)
      zdrv(i) = zdrv(i) - dxyz(3)
      xdrv(k) = xdrv(k) + dxyz(1)
      ydrv(k) = ydrv(k) + dxyz(2)
      zdrv(k) = zdrv(k) + dxyz(3)
    endif
!
    call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                      d2r2ikdsdx,d2r2ikds2,.true.)
    if (lstr) then
      if (lilocal) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + d1trm*dr2ikds(ks)
        enddo
      endif
    endif
!
    indk = 3*(k-1)
    kxf = indk + 1
    kyf = indk + 2
    kzf = indk + 3
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    liklocal = (lilocal.or.lklocal)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    if (lklocal) then
      indk = 3*(kloc-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
    endif
!
    d2cnidrik2 = d2cndr_cn(n,i)
    d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
           (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!
    if (abs(d2trm).gt.gfnff_cnc6tol) then
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
      if (lilocal) then
        derv2(kxf,ix) = derv2(kxf,ix) - d2trm*d2r2ikdx2(1)
        derv2(kyf,ix) = derv2(kyf,ix) - d2trm*d2r2ikdx2(6)
        derv2(kzf,ix) = derv2(kzf,ix) - d2trm*d2r2ikdx2(5)
        derv2(kxf,iy) = derv2(kxf,iy) - d2trm*d2r2ikdx2(6)
        derv2(kyf,iy) = derv2(kyf,iy) - d2trm*d2r2ikdx2(2)
        derv2(kzf,iy) = derv2(kzf,iy) - d2trm*d2r2ikdx2(4)
        derv2(kxf,iz) = derv2(kxf,iz) - d2trm*d2r2ikdx2(5)
        derv2(kyf,iz) = derv2(kyf,iz) - d2trm*d2r2ikdx2(4)
        derv2(kzf,iz) = derv2(kzf,iz) - d2trm*d2r2ikdx2(3)
        derv2(kxf,ix) = derv2(kxf,ix) - d1trm
        derv2(kyf,iy) = derv2(kyf,iy) - d1trm
        derv2(kzf,iz) = derv2(kzf,iz) - d1trm
      endif
      if (lklocal) then
        derv2(ixf,kx) = derv2(ixf,kx) - d2trm*d2r2ikdx2(1)
        derv2(iyf,kx) = derv2(iyf,kx) - d2trm*d2r2ikdx2(6)
        derv2(izf,kx) = derv2(izf,kx) - d2trm*d2r2ikdx2(5)
        derv2(ixf,ky) = derv2(ixf,ky) - d2trm*d2r2ikdx2(6)
        derv2(iyf,ky) = derv2(iyf,ky) - d2trm*d2r2ikdx2(2)
        derv2(izf,ky) = derv2(izf,ky) - d2trm*d2r2ikdx2(4)
        derv2(ixf,kz) = derv2(ixf,kz) - d2trm*d2r2ikdx2(5)
        derv2(iyf,kz) = derv2(iyf,kz) - d2trm*d2r2ikdx2(4)
        derv2(izf,kz) = derv2(izf,kz) - d2trm*d2r2ikdx2(3)
        derv2(ixf,kx) = derv2(ixf,kx) - d1trm
        derv2(iyf,ky) = derv2(iyf,ky) - d1trm
        derv2(izf,kz) = derv2(izf,kz) - d1trm
      endif
!
      if (lstr) then
!
!  Mixed derivatives
!
        if (lilocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ikdsdx(ks,1) - xnbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ikdsdx(ks,2) - ynbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ikdsdx(ks,3) - znbr(ni,i)*d2trm*dr2ikds(ks)
          enddo
        endif
        if (lklocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2ikdsdx(ks,1) + xnbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2ikdsdx(ks,2) + ynbr(ni,i)*d2trm*dr2ikds(ks)
            derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2ikdsdx(ks,3) + znbr(ni,i)*d2trm*dr2ikds(ks)
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
              sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2ikds(kt)*dr2ikds(ks) + d1trm*d2r2ikds2(kt,ks)
            enddo
          enddo
        endif
      endif
    endif
!
    d3trm1 = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      dcnidril = d1cndr_cn(n2,i)
      d3trm = d3trm1*dcnidril
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
      ni2 = nbrno_cn(n2,i)
      l = nbrno(ni2,i)
!
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
!  If none of the atoms are local then cycle
!
      if (.not.liklocal.and..not.lllocal) cycle
!
      indl = 3*(l-1)
      lxf = indl + 1
      lyf = indl + 2
      lzf = indl + 3
!
      if (lllocal) then
        indl = 3*(lloc-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
      endif
!
      if (lstr) then
        call real1strterm(ndim,xnbr(ni2,i),ynbr(ni2,i),znbr(ni2,i),0.0_dp,0.0_dp,0.0_dp,dr2ilds,d2r2ildx2, &
                          d2r2ildsdx,d2r2ilds2,.false.)
      endif
!
!  I-K
!
      if (lilocal) then
        derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
      if (lklocal) then
        derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  I-L
!
      if (lilocal) then
        derv2(lxf,ix) = derv2(lxf,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,ix) = derv2(lyf,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,ix) = derv2(lzf,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lxf,iy) = derv2(lxf,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,iy) = derv2(lyf,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,iy) = derv2(lzf,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lxf,iz) = derv2(lxf,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,iz) = derv2(lyf,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,iz) = derv2(lzf,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
      if (lllocal) then
        derv2(ixf,lx) = derv2(ixf,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,lx) = derv2(iyf,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(izf,lx) = derv2(izf,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ixf,ly) = derv2(ixf,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(iyf,ly) = derv2(iyf,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(izf,ly) = derv2(izf,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(ixf,lz) = derv2(ixf,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(iyf,lz) = derv2(iyf,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(izf,lz) = derv2(izf,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  K-L
!
      if (lklocal) then
        derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
      if (lllocal) then
        derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  Strain terms
!
      if (lstr) then
        if (lilocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni2,i)*d3trm*dr2ikds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ilds(ks) + dr2ilds(kt)*dr2ikds(ks))
            enddo
          enddo
        endif
!
        if (lklocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ilds(ks)
          enddo
        endif
!
        if (lllocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(ni2,i)*d3trm*dr2ikds(ks)
          enddo
        endif
      endif
    enddo
!
    d3trm1 = d2Edlogcnidlogcnj*dlogcnidcni*dlogcnjdcnj*dcnidrik
!------------------------------------------------------------------------------------------
!  Loop over neighbours of j for second derivatives for 2 different coordination numbers  |
!------------------------------------------------------------------------------------------
    do n2 = 1,nnbr_cn(j)
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = d3trm1*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
      nj = nbrno_cn(n2,j)
      l = nbrno(nj,j)
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
!  If none of the atoms are local then cycle
!
      if (.not.lijklocal.and..not.lllocal) cycle
!
      indl = 3*(l-1)
      lxf = indl + 1
      lyf = indl + 2
      lzf = indl + 3
!
      if (lllocal) then
        indl = 3*(lloc-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
      endif
!
      if (lstr) then
        call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jlds,d2r2jldx2, &
                          d2r2jldsdx,d2r2jlds2,.false.)
      endif
!
!  I-J
!
      if (lilocal) then
        derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*znbr(ni,i)
      endif
      if (ljlocal) then
        derv2(ixf,jx) = derv2(ixf,jx) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iyf,jx) = derv2(iyf,jx) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(izf,jx) = derv2(izf,jx) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ixf,jy) = derv2(ixf,jy) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iyf,jy) = derv2(iyf,jy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(izf,jy) = derv2(izf,jy) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ixf,jz) = derv2(ixf,jz) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iyf,jz) = derv2(iyf,jz) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(izf,jz) = derv2(izf,jz) + d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  I-L
!
      if (lilocal) then
        derv2(lxf,ix) = derv2(lxf,ix) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(lyf,ix) = derv2(lyf,ix) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(lzf,ix) = derv2(lzf,ix) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(lxf,iy) = derv2(lxf,iy) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(lyf,iy) = derv2(lyf,iy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(lzf,iy) = derv2(lzf,iy) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(lxf,iz) = derv2(lxf,iz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(lyf,iz) = derv2(lyf,iz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(lzf,iz) = derv2(lzf,iz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
      if (lllocal) then
        derv2(ixf,lx) = derv2(ixf,lx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iyf,lx) = derv2(iyf,lx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(izf,lx) = derv2(izf,lx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ixf,ly) = derv2(ixf,ly) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iyf,ly) = derv2(iyf,ly) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(izf,ly) = derv2(izf,ly) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ixf,lz) = derv2(ixf,lz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iyf,lz) = derv2(iyf,lz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(izf,lz) = derv2(izf,lz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  J-K
!
      if (ljlocal) then
        derv2(kxf,jx) = derv2(kxf,jx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(kyf,jx) = derv2(kyf,jx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(kzf,jx) = derv2(kzf,jx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(kxf,jy) = derv2(kxf,jy) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(kyf,jy) = derv2(kyf,jy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(kzf,jy) = derv2(kzf,jy) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(kxf,jz) = derv2(kxf,jz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(kyf,jz) = derv2(kyf,jz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(kzf,jz) = derv2(kzf,jz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
      if (lklocal) then
        derv2(jxf,kx) = derv2(jxf,kx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jyf,kx) = derv2(jyf,kx) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jzf,kx) = derv2(jzf,kx) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jxf,ky) = derv2(jxf,ky) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jyf,ky) = derv2(jyf,ky) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jzf,ky) = derv2(jzf,ky) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jxf,kz) = derv2(jxf,kz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jyf,kz) = derv2(jyf,kz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jzf,kz) = derv2(jzf,kz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  K-L
!     
      if (lklocal) then
        derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(ni,i)*znbr(nj,j)
      endif
      if (lllocal) then
        derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(ni,i)*znbr(nj,j)
      endif
!
!  Strain terms
!
      if (lstr) then
        if (lilocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2jlds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2ikds(ks))
            enddo
          enddo
        endif
!
        if (lklocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2jlds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2jlds(ks)
          enddo
        endif
!
        if (ljlocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2ikds(ks)
          enddo
        endif
!
        if (lllocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(nj,j)*d3trm*dr2ikds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(nj,j)*d3trm*dr2ikds(ks)
          enddo
        endif
      endif
    enddo
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    dcnjdrjk = d1cndr_cn(n,j)
    d1trm = dlogcnjdcnj*dcnjdrjk*dEdlogcnj
    dxyz(1) = d1trm*xnbr(nj,j)
    dxyz(2) = d1trm*ynbr(nj,j)
    dxyz(3) = d1trm*znbr(nj,j)
!
    if (lilocal) then
      xdrv(j) = xdrv(j) - dxyz(1)
      ydrv(j) = ydrv(j) - dxyz(2)
      zdrv(j) = zdrv(j) - dxyz(3)
      xdrv(k) = xdrv(k) + dxyz(1)
      ydrv(k) = ydrv(k) + dxyz(2)
      zdrv(k) = zdrv(k) + dxyz(3)
    endif
!
    call real1strterm(ndim,xnbr(nj,j),ynbr(nj,j),znbr(nj,j),0.0_dp,0.0_dp,0.0_dp,dr2jkds,d2r2jkdx2, &
                      d2r2jkdsdx,d2r2jkds2,.true.)
!
    if (lstr) then
      if (lilocal) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + d1trm*dr2jkds(ks)
        enddo
      endif
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
    d2cnjdrjk2 = d2cndr_cn(n,j)
    d2trm = dEdlogcnj*dlogcnjdcnj*d2cnjdrjk2 + &
           (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjk
!
    if (abs(d2trm).gt.gfnff_cnc6tol) then
!----------------------------------------
!  Second derivatives of CN w.r.t. rjk  |
!----------------------------------------
      if (ljlocal) then
        derv2(kxf,jx) = derv2(kxf,jx) - d2trm*d2r2jkdx2(1)
        derv2(kyf,jx) = derv2(kyf,jx) - d2trm*d2r2jkdx2(6)
        derv2(kzf,jx) = derv2(kzf,jx) - d2trm*d2r2jkdx2(5)
        derv2(kxf,jy) = derv2(kxf,jy) - d2trm*d2r2jkdx2(6)
        derv2(kyf,jy) = derv2(kyf,jy) - d2trm*d2r2jkdx2(2)
        derv2(kzf,jy) = derv2(kzf,jy) - d2trm*d2r2jkdx2(4)
        derv2(kxf,jz) = derv2(kxf,jz) - d2trm*d2r2jkdx2(5)
        derv2(kyf,jz) = derv2(kyf,jz) - d2trm*d2r2jkdx2(4)
        derv2(kzf,jz) = derv2(kzf,jz) - d2trm*d2r2jkdx2(3)
        derv2(kxf,jx) = derv2(kxf,jx) - d1trm
        derv2(kyf,jy) = derv2(kyf,jy) - d1trm
        derv2(kzf,jz) = derv2(kzf,jz) - d1trm
      endif
      if (lklocal) then
        derv2(jxf,kx) = derv2(jxf,kx) - d2trm*d2r2jkdx2(1)
        derv2(jyf,kx) = derv2(jyf,kx) - d2trm*d2r2jkdx2(6)
        derv2(jzf,kx) = derv2(jzf,kx) - d2trm*d2r2jkdx2(5)
        derv2(jxf,ky) = derv2(jxf,ky) - d2trm*d2r2jkdx2(6)
        derv2(jyf,ky) = derv2(jyf,ky) - d2trm*d2r2jkdx2(2)
        derv2(jzf,ky) = derv2(jzf,ky) - d2trm*d2r2jkdx2(4)
        derv2(jxf,kz) = derv2(jxf,kz) - d2trm*d2r2jkdx2(5)
        derv2(jyf,kz) = derv2(jyf,kz) - d2trm*d2r2jkdx2(4)
        derv2(jzf,kz) = derv2(jzf,kz) - d2trm*d2r2jkdx2(3)
        derv2(jxf,kx) = derv2(jxf,kx) - d1trm
        derv2(jyf,ky) = derv2(jyf,ky) - d1trm
        derv2(jzf,kz) = derv2(jzf,kz) - d1trm
      endif
!
      if (lstr) then
!
!  Mixed derivatives
!
        if (ljlocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - d1trm*d2r2jkdsdx(ks,1) - xnbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(jy,kl) = derv3(jy,kl) - d1trm*d2r2jkdsdx(ks,2) - ynbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(jz,kl) = derv3(jz,kl) - d1trm*d2r2jkdsdx(ks,3) - znbr(nj,j)*d2trm*dr2jkds(ks)
          enddo
        endif
        if (lklocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2jkdsdx(ks,1) + xnbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2jkdsdx(ks,2) + ynbr(nj,j)*d2trm*dr2jkds(ks)
            derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2jkdsdx(ks,3) + znbr(nj,j)*d2trm*dr2jkds(ks)
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
              sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2jkds(kt)*dr2jkds(ks) + d1trm*d2r2jkds2(kt,ks)
            enddo
          enddo
        endif
      endif
    endif
!
    d3trm1 = (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of j for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = d3trm1*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
      nj2 = nbrno_cn(n2,j)
      l = nbrno(nj2,j)
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
!  If none of the atoms are local then cycle
!
      if (.not.lijklocal.and..not.lllocal) cycle
!
      indl = 3*(l-1)
      lxf = indl + 1
      lyf = indl + 2
      lzf = indl + 3
!
      if (lllocal) then
        indl = 3*(lloc-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
      endif
!
      if (lstr) then
        call real1strterm(ndim,xnbr(nj2,j),ynbr(nj2,j),znbr(nj2,j),0.0_dp,0.0_dp,0.0_dp,dr2jlds,d2r2jldx2, &
                          d2r2jldsdx,d2r2jlds2,.false.)
      endif
!
!  J-K
!
      if (ljlocal) then
        derv2(kxf,jx) = derv2(kxf,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(kyf,jx) = derv2(kyf,jx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kzf,jx) = derv2(kzf,jx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kxf,jy) = derv2(kxf,jy) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(kyf,jy) = derv2(kyf,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kzf,jy) = derv2(kzf,jy) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kxf,jz) = derv2(kxf,jz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(kyf,jz) = derv2(kyf,jz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kzf,jz) = derv2(kzf,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
      if (lklocal) then
        derv2(jxf,kx) = derv2(jxf,kx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,kx) = derv2(jyf,kx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,kx) = derv2(jzf,kx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jxf,ky) = derv2(jxf,ky) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,ky) = derv2(jyf,ky) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,ky) = derv2(jzf,ky) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jxf,kz) = derv2(jxf,kz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,kz) = derv2(jyf,kz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,kz) = derv2(jzf,kz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  J-L
!
      if (ljlocal) then
        derv2(lxf,jx) = derv2(lxf,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,jx) = derv2(lyf,jx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,jx) = derv2(lzf,jx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lxf,jy) = derv2(lxf,jy) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,jy) = derv2(lyf,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,jy) = derv2(lzf,jy) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lxf,jz) = derv2(lxf,jz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,jz) = derv2(lyf,jz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,jz) = derv2(lzf,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
      if (lllocal) then
        derv2(jxf,lx) = derv2(jxf,lx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,lx) = derv2(jyf,lx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jzf,lx) = derv2(jzf,lx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jxf,ly) = derv2(jxf,ly) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jyf,ly) = derv2(jyf,ly) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,ly) = derv2(jzf,ly) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jxf,lz) = derv2(jxf,lz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jyf,lz) = derv2(jyf,lz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jzf,lz) = derv2(jzf,lz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  K-L
!     
      if (lklocal) then
        derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
      if (lllocal) then
        derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  Strain terms
!
      if (lstr) then
        if (ljlocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(jx,kl) = derv3(jx,kl) - xnbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(jy,kl) = derv3(jy,kl) - ynbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(jz,kl) = derv3(jz,kl) - znbr(nj2,j)*d3trm*dr2jkds(ks)
          enddo
        endif
!
        if (lklocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(nj,j)*d3trm*dr2jlds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(nj,j)*d3trm*dr2jlds(ks)
          enddo
        endif
!
        if (lllocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(nj2,j)*d3trm*dr2jkds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(nj2,j)*d3trm*dr2jkds(ks)
          enddo
        endif
!
        if (lilocal) then
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2jkds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2jkds(ks))
            enddo
          enddo
        endif
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnd_nox')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_self(i,dEdlogcni,d2Edlogcni2,dlogcndcn,d2logcndcn2,lgrad2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Special case version of gfnff_drv2_dcn for i only self terms
!
!  On entry : 
!
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   2/21 Created from gfnff_drv2_dcn 
!   2/21 Use of pointer to non-zero terms added
!   2/21 cut no longer passed in as not needed
!   2/21 Derivatives of coordination number now precomputed
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
!  Julian Gale, CIC, Curtin University, February 2021
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2ilds(6)
  real(dp)                                         :: d2r2ildx2(6)
  real(dp)                                         :: d2r2ilds2(6,6)
  real(dp)                                         :: d2r2ildsdx(6,3)
  real(dp)                                         :: dxyz(3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_self')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  if (lgrad2) then
    d2logcnidcni2 = d2logcndcn2(i)
  endif
!
!  Set up for second derivatives
!
  if (lgrad2) then
    indi = 3*(i-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
  endif
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
    dxyz(1) = d1trm*xnbr(ni,i)
    dxyz(2) = d1trm*ynbr(ni,i)
    dxyz(3) = d1trm*znbr(ni,i)
!
    xdrv(i) = xdrv(i) - dxyz(1)
    ydrv(i) = ydrv(i) - dxyz(2)
    zdrv(i) = zdrv(i) - dxyz(3)
    xdrv(k) = xdrv(k) + dxyz(1)
    ydrv(k) = ydrv(k) + dxyz(2)
    zdrv(k) = zdrv(k) + dxyz(3)
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                        d2r2ikdsdx,d2r2ikds2,lgrad2)
    endif
    if (lstr) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        rstrd(kl) = rstrd(kl) + d1trm*dr2ikds(ks)
      enddo
    endif
    if (lgrad2) then
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
!
      d2cnidrik2 = d2cndr_cn(n,i)
      d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
             (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
      if (k.ge.i) then
        derv2(kx,ix) = derv2(kx,ix) - d2trm*d2r2ikdx2(1)
        derv2(ky,ix) = derv2(ky,ix) - d2trm*d2r2ikdx2(6)
        derv2(kz,ix) = derv2(kz,ix) - d2trm*d2r2ikdx2(5)
        derv2(kx,iy) = derv2(kx,iy) - d2trm*d2r2ikdx2(6)
        derv2(ky,iy) = derv2(ky,iy) - d2trm*d2r2ikdx2(2)
        derv2(kz,iy) = derv2(kz,iy) - d2trm*d2r2ikdx2(4)
        derv2(kx,iz) = derv2(kx,iz) - d2trm*d2r2ikdx2(5)
        derv2(ky,iz) = derv2(ky,iz) - d2trm*d2r2ikdx2(4)
        derv2(kz,iz) = derv2(kz,iz) - d2trm*d2r2ikdx2(3)
        derv2(kx,ix) = derv2(kx,ix) - d1trm
        derv2(ky,iy) = derv2(ky,iy) - d1trm
        derv2(kz,iz) = derv2(kz,iz) - d1trm
      else
        derv2(ix,kx) = derv2(ix,kx) - d2trm*d2r2ikdx2(1)
        derv2(iy,kx) = derv2(iy,kx) - d2trm*d2r2ikdx2(6)
        derv2(iz,kx) = derv2(iz,kx) - d2trm*d2r2ikdx2(5)
        derv2(ix,ky) = derv2(ix,ky) - d2trm*d2r2ikdx2(6)
        derv2(iy,ky) = derv2(iy,ky) - d2trm*d2r2ikdx2(2)
        derv2(iz,ky) = derv2(iz,ky) - d2trm*d2r2ikdx2(4)
        derv2(ix,kz) = derv2(ix,kz) - d2trm*d2r2ikdx2(5)
        derv2(iy,kz) = derv2(iy,kz) - d2trm*d2r2ikdx2(4)
        derv2(iz,kz) = derv2(iz,kz) - d2trm*d2r2ikdx2(3)
        derv2(ix,kx) = derv2(ix,kx) - d1trm
        derv2(iy,ky) = derv2(iy,ky) - d1trm
        derv2(iz,kz) = derv2(iz,kz) - d1trm
      endif
!
      if (lstr) then
!
!  Mixed derivatives
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ikdsdx(ks,1) - xnbr(ni,i)*d2trm*dr2ikds(ks)
          derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ikdsdx(ks,2) - ynbr(ni,i)*d2trm*dr2ikds(ks)
          derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ikdsdx(ks,3) - znbr(ni,i)*d2trm*dr2ikds(ks)
          derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2ikdsdx(ks,1) + xnbr(ni,i)*d2trm*dr2ikds(ks)
          derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2ikdsdx(ks,2) + ynbr(ni,i)*d2trm*dr2ikds(ks)
          derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2ikdsdx(ks,3) + znbr(ni,i)*d2trm*dr2ikds(ks)
        enddo
!
!  Strain-strain
!
        do kk = 1,nstrains
          ks = nstrptr(kk)
          do kl = 1,nstrains
            kt = nstrptr(kl)
            sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2ikds(kt)*dr2ikds(ks) + d1trm*d2r2ikds2(kt,ks)
          enddo
        enddo
      endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
      do n2 = 1,n-1
        ni2 = nbrno_cn(n2,i)
        l = nbrno(ni2,i)
        indl = 3*(l-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
!
        dcnidril = d1cndr_cn(n2,i)
        d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
        if (lstr) then
          call real1strterm(ndim,xnbr(ni2,i),ynbr(ni2,i),znbr(ni2,i),0.0_dp,0.0_dp,0.0_dp,dr2ilds,d2r2ildx2, &
                            d2r2ildsdx,d2r2ilds2,.false.)
        endif
!
!  I-K
!
        if (k.ge.i) then
          derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        else
          derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  I-L
!
        if (l.ge.i) then
          derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ly,ix) = derv2(ly,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(lz,ix) = derv2(lz,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(lx,iy) = derv2(lx,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(lz,iy) = derv2(lz,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(lx,iz) = derv2(lx,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ly,iz) = derv2(ly,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        else
          derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(iy,lx) = derv2(iy,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(iz,lx) = derv2(iz,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ix,ly) = derv2(ix,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(iz,ly) = derv2(iz,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(ix,lz) = derv2(ix,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(iy,lz) = derv2(iy,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  K-L
!
        if (l.ge.k) then
          derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
        else
          derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
          derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
          derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
          derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
          derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
          derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
          derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
          derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
          derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
        endif
!
!  I-I - because this is on-diagonal it will be set by the translational invariance condition
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ilds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(lx,kl) = derv3(lx,kl) + xnbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(ly,kl) = derv3(ly,kl) + ynbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(lz,kl) = derv3(lz,kl) + znbr(ni2,i)*d3trm*dr2ikds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ilds(ks) + dr2ilds(kt)*dr2ikds(ks))
            enddo
          enddo
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_self')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_selfd_p1(iloc,i,dEdlogcni,d2Edlogcni2,dlogcndcn,d2logcndcn2,lgrad2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Special case version of gfnff_drv2_dcn for i only self terms
!  Distributed memory parallel version.
!  Part 1 where i is local to the processor.
!
!  On entry : 
!
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   2/22 Created from gfnff_drv2_dcn_self
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
!  Julian Gale, CIC, Curtin University, February 2022
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use parallel
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: iloc
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2ilds(6)
  real(dp)                                         :: d2r2ildx2(6)
  real(dp)                                         :: d2r2ilds2(6,6)
  real(dp)                                         :: d2r2ildsdx(6,3)
  real(dp)                                         :: dxyz(3)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_selfd_p1')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  if (lgrad2) then
    d2logcnidcni2 = d2logcndcn2(i)
  endif
!
!  Set up for second derivatives
!
  if (lgrad2) then
    indi = 3*(iloc-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
  endif
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
    dxyz(1) = d1trm*xnbr(ni,i)
    dxyz(2) = d1trm*ynbr(ni,i)
    dxyz(3) = d1trm*znbr(ni,i)
!
    xdrv(i) = xdrv(i) - dxyz(1)
    ydrv(i) = ydrv(i) - dxyz(2)
    zdrv(i) = zdrv(i) - dxyz(3)
    xdrv(k) = xdrv(k) + dxyz(1)
    ydrv(k) = ydrv(k) + dxyz(2)
    zdrv(k) = zdrv(k) + dxyz(3)
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                        d2r2ikdsdx,d2r2ikds2,lgrad2)
    endif
    if (lstr) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        rstrd(kl) = rstrd(kl) + d1trm*dr2ikds(ks)
      enddo
    endif
    if (lgrad2) then
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
!
      d2cnidrik2 = d2cndr_cn(n,i)
      d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
             (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
      derv2(kx,ix) = derv2(kx,ix) - d2trm*d2r2ikdx2(1)
      derv2(ky,ix) = derv2(ky,ix) - d2trm*d2r2ikdx2(6)
      derv2(kz,ix) = derv2(kz,ix) - d2trm*d2r2ikdx2(5)
      derv2(kx,iy) = derv2(kx,iy) - d2trm*d2r2ikdx2(6)
      derv2(ky,iy) = derv2(ky,iy) - d2trm*d2r2ikdx2(2)
      derv2(kz,iy) = derv2(kz,iy) - d2trm*d2r2ikdx2(4)
      derv2(kx,iz) = derv2(kx,iz) - d2trm*d2r2ikdx2(5)
      derv2(ky,iz) = derv2(ky,iz) - d2trm*d2r2ikdx2(4)
      derv2(kz,iz) = derv2(kz,iz) - d2trm*d2r2ikdx2(3)
      derv2(kx,ix) = derv2(kx,ix) - d1trm
      derv2(ky,iy) = derv2(ky,iy) - d1trm
      derv2(kz,iz) = derv2(kz,iz) - d1trm
!
      if (lstr) then
!
!  Mixed derivatives
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ikdsdx(ks,1) - xnbr(ni,i)*d2trm*dr2ikds(ks)
          derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ikdsdx(ks,2) - ynbr(ni,i)*d2trm*dr2ikds(ks)
          derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ikdsdx(ks,3) - znbr(ni,i)*d2trm*dr2ikds(ks)
        enddo
!
!  Strain-strain
!
        do kk = 1,nstrains
          ks = nstrptr(kk)
          do kl = 1,nstrains
            kt = nstrptr(kl)
            sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2ikds(kt)*dr2ikds(ks) + d1trm*d2r2ikds2(kt,ks)
          enddo
        enddo
      endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
      do n2 = 1,n-1
        ni2 = nbrno_cn(n2,i)
        l = nbrno(ni2,i)
        indl = 3*(l-1)
        lx = indl + 1
        ly = indl + 2
        lz = indl + 3
!
        dcnidril = d1cndr_cn(n2,i)
        d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
        if (lstr) then
          call real1strterm(ndim,xnbr(ni2,i),ynbr(ni2,i),znbr(ni2,i),0.0_dp,0.0_dp,0.0_dp,dr2ilds,d2r2ildx2, &
                            d2r2ildsdx,d2r2ilds2,.false.)
        endif
!
!  I-K
!
        derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
!
!  I-L
!
        derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ly,ix) = derv2(ly,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lz,ix) = derv2(lz,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lx,iy) = derv2(lx,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lz,iy) = derv2(lz,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lx,iz) = derv2(lx,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ly,iz) = derv2(ly,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
!
!  I-I - because this is on-diagonal it will be set by the translational invariance condition
!
!  Strain terms
!
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni,i)*d3trm*dr2ilds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni,i)*d3trm*dr2ilds(ks)
          enddo
!
          do kl = 1,nstrains
            ks = nstrptr(kl)
            derv3(ix,kl) = derv3(ix,kl) - xnbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iy,kl) = derv3(iy,kl) - ynbr(ni2,i)*d3trm*dr2ikds(ks)
            derv3(iz,kl) = derv3(iz,kl) - znbr(ni2,i)*d3trm*dr2ikds(ks)
          enddo
!
          do kk = 1,nstrains
            ks = nstrptr(kk)
            do kl = 1,nstrains
              kt = nstrptr(kl)
              sderv2(kl,kk) = sderv2(kl,kk) + d3trm*(dr2ikds(kt)*dr2ilds(ks) + dr2ilds(kt)*dr2ikds(ks))
            enddo
          enddo
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_selfd_p1')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_selfd_p2(i,dEdlogcni,d2Edlogcni2,dlogcndcn,d2logcndcn2,lgrad2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Special case version of gfnff_drv2_dcn for i only self terms
!  Distributed memory parallel version.
!  Part 2 where i is global
!
!  On entry : 
!
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   2/22 Created from gfnff_drv2_dcn_self
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
!  Julian Gale, CIC, Curtin University, February 2022
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use parallel
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: k
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: ks
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2ilds(6)
  real(dp)                                         :: d2r2ildx2(6)
  real(dp)                                         :: d2r2ilds2(6,6)
  real(dp)                                         :: d2r2ildsdx(6,3)
!
!  Only needed for second derivatives
!
  if (.not.lgrad2) return
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_selfd_p2')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d2logcnidcni2 = d2logcndcn2(i)
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    kloc = atom2local(k)
!
!  k is not local then cycle
!
    if (kloc.eq.0) cycle
!
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
!
    call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ikds,d2r2ikdx2, &
                      d2r2ikdsdx,d2r2ikds2,lgrad2)
!
    indk = 3*(kloc-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!
    d2cnidrik2 = d2cndr_cn(n,i)
    d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
           (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
    derv2(ix,kx) = derv2(ix,kx) - d2trm*d2r2ikdx2(1)
    derv2(iy,kx) = derv2(iy,kx) - d2trm*d2r2ikdx2(6)
    derv2(iz,kx) = derv2(iz,kx) - d2trm*d2r2ikdx2(5)
    derv2(ix,ky) = derv2(ix,ky) - d2trm*d2r2ikdx2(6)
    derv2(iy,ky) = derv2(iy,ky) - d2trm*d2r2ikdx2(2)
    derv2(iz,ky) = derv2(iz,ky) - d2trm*d2r2ikdx2(4)
    derv2(ix,kz) = derv2(ix,kz) - d2trm*d2r2ikdx2(5)
    derv2(iy,kz) = derv2(iy,kz) - d2trm*d2r2ikdx2(4)
    derv2(iz,kz) = derv2(iz,kz) - d2trm*d2r2ikdx2(3)
    derv2(ix,kx) = derv2(ix,kx) - d1trm
    derv2(iy,ky) = derv2(iy,ky) - d1trm
    derv2(iz,kz) = derv2(iz,kz) - d1trm
!
    if (lstr) then
!
!  Mixed derivatives
!
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(kx,kl) = derv3(kx,kl) + d1trm*d2r2ikdsdx(ks,1) + xnbr(ni,i)*d2trm*dr2ikds(ks)
        derv3(ky,kl) = derv3(ky,kl) + d1trm*d2r2ikdsdx(ks,2) + ynbr(ni,i)*d2trm*dr2ikds(ks)
        derv3(kz,kl) = derv3(kz,kl) + d1trm*d2r2ikdsdx(ks,3) + znbr(ni,i)*d2trm*dr2ikds(ks)
      enddo
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,nnbr_cn(i)
!
!  Skip n2 = n
!
      if (n2.eq.n) cycle
!
      ni2 = nbrno_cn(n2,i)
      l = nbrno(ni2,i)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnidril = d1cndr_cn(n2,i)
      d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
      if (lstr) then
        call real1strterm(ndim,xnbr(ni2,i),ynbr(ni2,i),znbr(ni2,i),0.0_dp,0.0_dp,0.0_dp,dr2ilds,d2r2ildx2, &
                          d2r2ildsdx,d2r2ilds2,.false.)
      endif
!
!  I-K
!
      derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
      derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
      derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
      derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
      derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
      derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
      derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
      derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
      derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
!
!  K-L
!
      derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
      derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
      derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
      derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
      derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
      derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
      derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
      derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
      derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
!
!  I-I - because this is on-diagonal it will be set by the translational invariance condition
!
!  Strain terms
!
      if (lstr) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(kx,kl) = derv3(kx,kl) + xnbr(ni,i)*d3trm*dr2ilds(ks)
          derv3(ky,kl) = derv3(ky,kl) + ynbr(ni,i)*d3trm*dr2ilds(ks)
          derv3(kz,kl) = derv3(kz,kl) + znbr(ni,i)*d3trm*dr2ilds(ks)
        enddo
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_selfd_p2')
#endif
!
  return
  end
!
!
  subroutine gfnff_drv2_dcn_dq(i,iloc,d2Edlogcnidq,dlogcndcn)
!
!  Computes the second derivatives of the log coordination number and charge for GFNFF
!
!  On entry : 
!
!  d2Edlogcnidq      = second derivatives of energy w.r.t. log of the coordination number for i and charge i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!   2/21 Created from gfnff_drv2_dcn_self
!   2/21 Use of pointer to non-zero terms added
!   2/21 cut no longer passed in as not needed
!   2/21 Derivatives of coordination number now precomputed
!   1/22 iloc argument added to point to position of dqdxyz terms in parallel
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
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: iloc
  real(dp),    intent(in)                          :: d2Edlogcnidq
  real(dp),    intent(in)                          :: dlogcndcn(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  real(dp)                                         :: d1trm
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dcnidrij 
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dijx
  real(dp)                                         :: dijy
  real(dp)                                         :: dijz
  real(dp)                                         :: dikx
  real(dp)                                         :: diky
  real(dp)                                         :: dikz
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_dq')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d1trm = d2Edlogcnidq*dlogcnidcni
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
!  Loop over coordination number pairs
!
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    j = nbrno(ni,i)
    indj = 3*(j-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
    dcnidrij = d1cndr_cn(n,i)
    d1trm = d2Edlogcnidq*dlogcnidcni*dcnidrij
    dijx = d1trm*xnbr(ni,i)
    dijy = d1trm*ynbr(ni,i)
    dijz = d1trm*znbr(ni,i)
!
    if (lstr) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2, &
                        d2r2ijdsdx,d2r2ijds2,.false.)
    endif
!
!  Loop over atoms for charge derivatives
!
    do k = 1,numat
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
      dikx = dqdxyz(indk+1,iloc)
      diky = dqdxyz(indk+2,iloc)
      dikz = dqdxyz(indk+3,iloc)
!
!  NB: If order of indices of derv2 are changed then swap the j/k for the term that is added
!
      if (kx.gt.ix) then
        derv2(kx,ix) = derv2(kx,ix) - dijx*dikx
        derv2(ky,ix) = derv2(ky,ix) - dijx*diky
        derv2(kz,ix) = derv2(kz,ix) - dijx*dikz
        derv2(kx,iy) = derv2(kx,iy) - dijy*dikx
        derv2(ky,iy) = derv2(ky,iy) - dijy*diky
        derv2(kz,iy) = derv2(kz,iy) - dijy*dikz
        derv2(kx,iz) = derv2(kx,iz) - dijz*dikx
        derv2(ky,iz) = derv2(ky,iz) - dijz*diky
        derv2(kz,iz) = derv2(kz,iz) - dijz*dikz
      else
        derv2(ix,kx) = derv2(ix,kx) - dikx*dijx
        derv2(iy,kx) = derv2(iy,kx) - diky*dijx
        derv2(iz,kx) = derv2(iz,kx) - dikz*dijx
        derv2(ix,ky) = derv2(ix,ky) - dikx*dijy
        derv2(iy,ky) = derv2(iy,ky) - diky*dijy
        derv2(iz,ky) = derv2(iz,ky) - dikz*dijy
        derv2(ix,kz) = derv2(ix,kz) - dikx*dijz
        derv2(iy,kz) = derv2(iy,kz) - diky*dijz
        derv2(iz,kz) = derv2(iz,kz) - dikz*dijz
      endif
!
      if (kx.gt.jx) then
        derv2(kx,jx) = derv2(kx,jx) + dijx*dikx
        derv2(ky,jx) = derv2(ky,jx) + dijx*diky
        derv2(kz,jx) = derv2(kz,jx) + dijx*dikz
        derv2(kx,jy) = derv2(kx,jy) + dijy*dikx
        derv2(ky,jy) = derv2(ky,jy) + dijy*diky
        derv2(kz,jy) = derv2(kz,jy) + dijy*dikz
        derv2(kx,jz) = derv2(kx,jz) + dijz*dikx
        derv2(ky,jz) = derv2(ky,jz) + dijz*diky
        derv2(kz,jz) = derv2(kz,jz) + dijz*dikz
      else
        derv2(jx,kx) = derv2(jx,kx) + dikx*dijx
        derv2(jy,kx) = derv2(jy,kx) + diky*dijx
        derv2(jz,kx) = derv2(jz,kx) + dikz*dijx
        derv2(jx,ky) = derv2(jx,ky) + dikx*dijy
        derv2(jy,ky) = derv2(jy,ky) + diky*dijy
        derv2(jz,ky) = derv2(jz,ky) + dikz*dijy
        derv2(jx,kz) = derv2(jx,kz) + dikx*dijz
        derv2(jy,kz) = derv2(jy,kz) + diky*dijz
        derv2(jz,kz) = derv2(jz,kz) + dikz*dijz
      endif
!
      if (ix.gt.jx) then
        derv2(ix,jx) = derv2(ix,jx) - dijx*dikx
        derv2(iy,jx) = derv2(iy,jx) - dijx*diky
        derv2(iz,jx) = derv2(iz,jx) - dijx*dikz
        derv2(ix,jy) = derv2(ix,jy) - dijy*dikx
        derv2(iy,jy) = derv2(iy,jy) - dijy*diky
        derv2(iz,jy) = derv2(iz,jy) - dijy*dikz
        derv2(ix,jz) = derv2(ix,jz) - dijz*dikx
        derv2(iy,jz) = derv2(iy,jz) - dijz*diky
        derv2(iz,jz) = derv2(iz,jz) - dijz*dikz
      else
        derv2(jx,ix) = derv2(jx,ix) - dikx*dijx
        derv2(jy,ix) = derv2(jy,ix) - diky*dijx
        derv2(jz,ix) = derv2(jz,ix) - dikz*dijx
        derv2(jx,iy) = derv2(jx,iy) - dikx*dijy
        derv2(jy,iy) = derv2(jy,iy) - diky*dijy
        derv2(jz,iy) = derv2(jz,iy) - dikz*dijy
        derv2(jx,iz) = derv2(jx,iz) - dikx*dijz
        derv2(jy,iz) = derv2(jy,iz) - diky*dijz
        derv2(jz,iz) = derv2(jz,iz) - dikz*dijz
      endif
!
      if (lstr) then
!
!  Mixed derivatives
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(ix,kl) = derv3(ix,kl) - d1trm*dr2ijds(ks)*dikx
          derv3(iy,kl) = derv3(iy,kl) - d1trm*dr2ijds(ks)*diky
          derv3(iz,kl) = derv3(iz,kl) - d1trm*dr2ijds(ks)*dikz
          derv3(kx,kl) = derv3(kx,kl) + d1trm*dr2ijds(ks)*dikx
          derv3(ky,kl) = derv3(ky,kl) + d1trm*dr2ijds(ks)*diky
          derv3(kz,kl) = derv3(kz,kl) + d1trm*dr2ijds(ks)*dikz
        enddo
      endif
!
!  End of loop over k
!
    enddo
!
!  Strain contributions
!
    if (lstr) then
!
!  Mixed derivatives
!
      do kl = 1,nstrains
        derv3(ix,kl) = derv3(ix,kl) - dijx*dqds(kl,iloc)
        derv3(iy,kl) = derv3(iy,kl) - dijy*dqds(kl,iloc)
        derv3(iz,kl) = derv3(iz,kl) - dijz*dqds(kl,iloc)
        derv3(jx,kl) = derv3(jx,kl) + dijx*dqds(kl,iloc)
        derv3(jy,kl) = derv3(jy,kl) + dijy*dqds(kl,iloc)
        derv3(jz,kl) = derv3(jz,kl) + dijz*dqds(kl,iloc)
      enddo
!
!  Strain-strain
!
      do kk = 1,nstrains
        ks = nstrptr(kk)
        do kl = 1,nstrains
          kt = nstrptr(kl)
          sderv2(kl,kk) = sderv2(kl,kk) + d1trm*(dr2ijds(kt)*dqds(kk,iloc) + dr2ijds(ks)*dqds(kl,iloc))
        enddo
      enddo
    endif
!
!  End of loop over j
!
  enddo

#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_dq')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_dqd_p1(iloc,i,d2Edlogcnidq,dlogcndcn)
!
!  Computes the second derivatives of the log coordination number and charge for GFNFF
!  Distributed memory parallel version.
!  Part 1 where i is loc to processor.
!
!  On entry : 
!
!  d2Edlogcnidq      = second derivatives of energy w.r.t. log of the coordination number for i and charge i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!   2/22 Created from gfnff_drv2_dcn_dq
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
!  Julian Gale, CIC, Curtin University, February 2022
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use parallel
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: iloc
  real(dp),    intent(in)                          :: d2Edlogcnidq
  real(dp),    intent(in)                          :: dlogcndcn(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  real(dp)                                         :: d1trm
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dcnidrij 
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dijx
  real(dp)                                         :: dijy
  real(dp)                                         :: dijz
  real(dp)                                         :: dikx
  real(dp)                                         :: diky
  real(dp)                                         :: dikz
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_dqd_p1')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d1trm = d2Edlogcnidq*dlogcnidcni
!
!  Set up for second derivatives
!
  indi = 3*(iloc-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
!  Loop over coordination number pairs
!
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    j = nbrno(ni,i)
    indj = 3*(j-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
    dcnidrij = d1cndr_cn(n,i)
    d1trm = d2Edlogcnidq*dlogcnidcni*dcnidrij
    dijx = d1trm*xnbr(ni,i)
    dijy = d1trm*ynbr(ni,i)
    dijz = d1trm*znbr(ni,i)
!
    if (lstr) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2, &
                        d2r2ijdsdx,d2r2ijds2,.false.)
    endif
!
!  Loop over atoms for charge derivatives
!
    do k = 1,numat
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
      dikx = dqdxyz(indk+1,i)
      diky = dqdxyz(indk+2,i)
      dikz = dqdxyz(indk+3,i)
!
      if (k.gt.i) then
        derv2(kx,ix) = derv2(kx,ix) - dijx*dikx
        derv2(ky,ix) = derv2(ky,ix) - dijx*diky
        derv2(kz,ix) = derv2(kz,ix) - dijx*dikz
        derv2(kx,iy) = derv2(kx,iy) - dijy*dikx
        derv2(ky,iy) = derv2(ky,iy) - dijy*diky
        derv2(kz,iy) = derv2(kz,iy) - dijy*dikz
        derv2(kx,iz) = derv2(kx,iz) - dijz*dikx
        derv2(ky,iz) = derv2(ky,iz) - dijz*diky
        derv2(kz,iz) = derv2(kz,iz) - dijz*dikz
      else
        derv2(kx,ix) = derv2(kx,ix) - dijx*dikx
        derv2(ky,ix) = derv2(ky,ix) - dijy*dikx
        derv2(kz,ix) = derv2(kz,ix) - dijz*dikx
        derv2(kx,iy) = derv2(kx,iy) - dijx*diky
        derv2(ky,iy) = derv2(ky,iy) - dijy*diky
        derv2(kz,iy) = derv2(kz,iy) - dijz*diky
        derv2(kx,iz) = derv2(kx,iz) - dijx*dikz
        derv2(ky,iz) = derv2(ky,iz) - dijy*dikz
        derv2(kz,iz) = derv2(kz,iz) - dijz*dikz
      endif
!
      if (j.gt.i) then
        derv2(jx,ix) = derv2(jx,ix) - dikx*dijx
        derv2(jy,ix) = derv2(jy,ix) - diky*dijx
        derv2(jz,ix) = derv2(jz,ix) - dikz*dijx
        derv2(jx,iy) = derv2(jx,iy) - dikx*dijy
        derv2(jy,iy) = derv2(jy,iy) - diky*dijy
        derv2(jz,iy) = derv2(jz,iy) - dikz*dijy
        derv2(jx,iz) = derv2(jx,iz) - dikx*dijz
        derv2(jy,iz) = derv2(jy,iz) - diky*dijz
        derv2(jz,iz) = derv2(jz,iz) - dikz*dijz
      else
        derv2(jx,ix) = derv2(jx,ix) - dikx*dijx
        derv2(jy,ix) = derv2(jy,ix) - dikx*dijy
        derv2(jz,ix) = derv2(jz,ix) - dikx*dijz
        derv2(jx,iy) = derv2(jx,iy) - diky*dijx
        derv2(jy,iy) = derv2(jy,iy) - diky*dijy
        derv2(jz,iy) = derv2(jz,iy) - diky*dijz
        derv2(jx,iz) = derv2(jx,iz) - dikz*dijx
        derv2(jy,iz) = derv2(jy,iz) - dikz*dijy
        derv2(jz,iz) = derv2(jz,iz) - dikz*dijz
      endif
!
      if (lstr) then
!
!  Mixed derivatives
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(ix,kl) = derv3(ix,kl) - d1trm*dr2ijds(ks)*dikx
          derv3(iy,kl) = derv3(iy,kl) - d1trm*dr2ijds(ks)*diky
          derv3(iz,kl) = derv3(iz,kl) - d1trm*dr2ijds(ks)*dikz
        enddo
      endif
!
!  End of loop over k
!
    enddo
!
!  Strain contributions
!
    if (lstr) then
!
!  Mixed derivatives
!
      do kl = 1,nstrains
        derv3(ix,kl) = derv3(ix,kl) - dijx*dqds(kl,i)
        derv3(iy,kl) = derv3(iy,kl) - dijy*dqds(kl,i)
        derv3(iz,kl) = derv3(iz,kl) - dijz*dqds(kl,i)
      enddo
!
!  Strain-strain
!
      do kk = 1,nstrains
        ks = nstrptr(kk)
        do kl = 1,nstrains
          kt = nstrptr(kl)
          sderv2(kl,kk) = sderv2(kl,kk) + d1trm*(dr2ijds(kt)*dqds(kk,i) + dr2ijds(ks)*dqds(kl,i))
        enddo
      enddo
    endif
!
!  End of loop over j
!
  enddo

#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_dqd_p1')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_dqd_p2(i,d2Edlogcnidq,dlogcndcn)
!
!  Computes the second derivatives of the log coordination number and charge for GFNFF
!  Distributed memory parallel version.
!  Part 2 where i is global.
!
!  On entry : 
!
!  d2Edlogcnidq      = second derivatives of energy w.r.t. log of the coordination number for i and charge i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2            = if .true. then compute the second derivatives
!
!   2/22 Created from gfnff_drv2_dcn_dq
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
!  Julian Gale, CIC, Curtin University, February 2022
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use parallel
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: d2Edlogcnidq
  real(dp),    intent(in)                          :: dlogcndcn(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jloc
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: jxf
  integer(i4)                                      :: jyf
  integer(i4)                                      :: jzf
  integer(i4)                                      :: k
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: ks
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  real(dp)                                         :: d1trm
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dcnidrij 
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dijx
  real(dp)                                         :: dijy
  real(dp)                                         :: dijz
  real(dp)                                         :: dikx
  real(dp)                                         :: diky
  real(dp)                                         :: dikz
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_dqd_p2')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d1trm = d2Edlogcnidq*dlogcnidcni
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
!  Loop over coordination number pairs
!
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    j = nbrno(ni,i)
    jloc = atom2local(j)
    if (jloc.ne.0) then
      indj = 3*(jloc-1)
      jx = indj + 1
      jy = indj + 2
      jz = indj + 3
    endif
!
    indj = 3*(j-1)
    jxf = indj + 1
    jyf = indj + 2
    jzf = indj + 3
!
    dcnidrij = d1cndr_cn(n,i)
    d1trm = d2Edlogcnidq*dlogcnidcni*dcnidrij
    dijx = d1trm*xnbr(ni,i)
    dijy = d1trm*ynbr(ni,i)
    dijz = d1trm*znbr(ni,i)
!
    if (lstr) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2, &
                        d2r2ijdsdx,d2r2ijds2,.false.)
    endif
!
!  Loop over atoms for charge derivatives
!
    do k = 1,numat
      kloc = atom2local(k)
!
!  If neither j or k are local then cycle
!
      if (jloc+kloc.eq.0) cycle
!
      if (kloc.ne.0) then
        indk = 3*(kloc-1)
        kx = indk + 1
        ky = indk + 2
        kz = indk + 3
      endif
!
      indk = 3*(k-1)
      kxf = indk + 1
      kyf = indk + 2
      kzf = indk + 3
!
      dikx = dqdxyz(indk+1,i)
      diky = dqdxyz(indk+2,i)
      dikz = dqdxyz(indk+3,i)
!
!  NB: If order of indices of derv2 are changed then swap the j/k for the term that is added
!
      if (kloc.ne.0) then
        if (k.gt.i) then
          derv2(ix,kx) = derv2(ix,kx) - dikx*dijx
          derv2(iy,kx) = derv2(iy,kx) - dikx*dijy
          derv2(iz,kx) = derv2(iz,kx) - dikx*dijz
          derv2(ix,ky) = derv2(ix,ky) - diky*dijx
          derv2(iy,ky) = derv2(iy,ky) - diky*dijy
          derv2(iz,ky) = derv2(iz,ky) - diky*dijz
          derv2(ix,kz) = derv2(ix,kz) - dikz*dijx
          derv2(iy,kz) = derv2(iy,kz) - dikz*dijy
          derv2(iz,kz) = derv2(iz,kz) - dikz*dijz
        else
          derv2(ix,kx) = derv2(ix,kx) - dikx*dijx
          derv2(iy,kx) = derv2(iy,kx) - diky*dijx
          derv2(iz,kx) = derv2(iz,kx) - dikz*dijx
          derv2(ix,ky) = derv2(ix,ky) - dikx*dijy
          derv2(iy,ky) = derv2(iy,ky) - diky*dijy
          derv2(iz,ky) = derv2(iz,ky) - dikz*dijy
          derv2(ix,kz) = derv2(ix,kz) - dikx*dijz
          derv2(iy,kz) = derv2(iy,kz) - diky*dijz
          derv2(iz,kz) = derv2(iz,kz) - dikz*dijz
        endif
      endif
!
      if (jloc.ne.0) then
        if (k.gt.j) then
          derv2(kxf,jx) = derv2(kxf,jx) + dijx*dikx
          derv2(kyf,jx) = derv2(kyf,jx) + dijx*diky
          derv2(kzf,jx) = derv2(kzf,jx) + dijx*dikz
          derv2(kxf,jy) = derv2(kxf,jy) + dijy*dikx
          derv2(kyf,jy) = derv2(kyf,jy) + dijy*diky
          derv2(kzf,jy) = derv2(kzf,jy) + dijy*dikz
          derv2(kxf,jz) = derv2(kxf,jz) + dijz*dikx
          derv2(kyf,jz) = derv2(kyf,jz) + dijz*diky
          derv2(kzf,jz) = derv2(kzf,jz) + dijz*dikz
        else
          derv2(kxf,jx) = derv2(kxf,jx) + dijx*dikx
          derv2(kyf,jx) = derv2(kyf,jx) + dijy*dikx
          derv2(kzf,jx) = derv2(kzf,jx) + dijz*dikx
          derv2(kxf,jy) = derv2(kxf,jy) + dijx*diky
          derv2(kyf,jy) = derv2(kyf,jy) + dijy*diky
          derv2(kzf,jy) = derv2(kzf,jy) + dijz*diky
          derv2(kxf,jz) = derv2(kxf,jz) + dijx*dikz
          derv2(kyf,jz) = derv2(kyf,jz) + dijy*dikz
          derv2(kzf,jz) = derv2(kzf,jz) + dijz*dikz
        endif
      endif
!
      if (kloc.ne.0) then
        if (j.gt.k) then
          derv2(jxf,kx) = derv2(jxf,kx) + dikx*dijx
          derv2(jyf,kx) = derv2(jyf,kx) + diky*dijx
          derv2(jzf,kx) = derv2(jzf,kx) + dikz*dijx
          derv2(jxf,ky) = derv2(jxf,ky) + dikx*dijy
          derv2(jyf,ky) = derv2(jyf,ky) + diky*dijy
          derv2(jzf,ky) = derv2(jzf,ky) + dikz*dijy
          derv2(jxf,kz) = derv2(jxf,kz) + dikx*dijz
          derv2(jyf,kz) = derv2(jyf,kz) + diky*dijz
          derv2(jzf,kz) = derv2(jzf,kz) + dikz*dijz
        else
          derv2(jxf,kx) = derv2(jxf,kx) + dikx*dijx
          derv2(jyf,kx) = derv2(jyf,kx) + dikx*dijy
          derv2(jzf,kx) = derv2(jzf,kx) + dikx*dijz
          derv2(jxf,ky) = derv2(jxf,ky) + diky*dijx
          derv2(jyf,ky) = derv2(jyf,ky) + diky*dijy
          derv2(jzf,ky) = derv2(jzf,ky) + diky*dijz
          derv2(jxf,kz) = derv2(jxf,kz) + dikz*dijx
          derv2(jyf,kz) = derv2(jyf,kz) + dikz*dijy
          derv2(jzf,kz) = derv2(jzf,kz) + dikz*dijz
        endif
      endif
!
      if (jloc.ne.0) then
        if (j.gt.i) then
          derv2(ix,jx) = derv2(ix,jx) - dijx*dikx
          derv2(iy,jx) = derv2(iy,jx) - dijy*dikx
          derv2(iz,jx) = derv2(iz,jx) - dijz*dikx
          derv2(ix,jy) = derv2(ix,jy) - dijx*diky
          derv2(iy,jy) = derv2(iy,jy) - dijy*diky
          derv2(iz,jy) = derv2(iz,jy) - dijz*diky
          derv2(ix,jz) = derv2(ix,jz) - dijx*dikz
          derv2(iy,jz) = derv2(iy,jz) - dijy*dikz
          derv2(iz,jz) = derv2(iz,jz) - dijz*dikz
        else
          derv2(ix,jx) = derv2(ix,jx) - dijx*dikx
          derv2(iy,jx) = derv2(iy,jx) - dijx*diky
          derv2(iz,jx) = derv2(iz,jx) - dijx*dikz
          derv2(ix,jy) = derv2(ix,jy) - dijy*dikx
          derv2(iy,jy) = derv2(iy,jy) - dijy*diky
          derv2(iz,jy) = derv2(iz,jy) - dijy*dikz
          derv2(ix,jz) = derv2(ix,jz) - dijz*dikx
          derv2(iy,jz) = derv2(iy,jz) - dijz*diky
          derv2(iz,jz) = derv2(iz,jz) - dijz*dikz
        endif
      endif
!
      if (lstr.and.kloc.ne.0) then
!
!  Mixed derivatives
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          derv3(kx,kl) = derv3(kx,kl) + d1trm*dr2ijds(ks)*dikx
          derv3(ky,kl) = derv3(ky,kl) + d1trm*dr2ijds(ks)*diky
          derv3(kz,kl) = derv3(kz,kl) + d1trm*dr2ijds(ks)*dikz
        enddo
      endif
!
!  End of loop over k
!
    enddo
!
!  Strain contributions
!
    if (lstr.and.jloc.ne.0) then
!
!  Mixed derivatives
!
      do kl = 1,nstrains
        derv3(jx,kl) = derv3(jx,kl) + dijx*dqds(kl,i)
        derv3(jy,kl) = derv3(jy,kl) + dijy*dqds(kl,i)
        derv3(jz,kl) = derv3(jz,kl) + dijz*dqds(kl,i)
      enddo
    endif
!
!  End of loop over j
!
  enddo

#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_dqd_p2')
#endif
!
  return
  end
!
  subroutine gfnff_drv_dcn_q(i,cn,dcn,dqx,dqy,dqz,dqs,maxlhs)
!
!  Computes the derivatives of the coordination number contribution to Z in EEM
!  and adds to the charge first derivatives
!  Has to compute the terms for the derivatives of atom i
!
!  On entry : 
!
!  cn              = coordination number
!  dcn             = derivatives of coordination number
!
!  On exit :
!
!  dqx             = contribution of dZ/dx due to CN added
!  dqy             = contribution of dZ/dy due to CN added
!  dqz             = contribution of dZ/dz due to CN added
!  dqs             = contribution of dZ/ds due to CN added
!  maxlhs          = first dimension of dqs
!
!   2/21 Created
!   2/21 Use of pointer to non-zero terms added
!   2/21 cut no longer passed in as not needed
!   2/21 Derivatives of coordination number now precomputed
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
!  Julian Gale, CIC, Curtin University, February 2021
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: maxlhs
  real(dp),    intent(inout)                       :: cn(*)
  real(dp),    intent(inout)                       :: dcn(*)
  real(dp),    intent(inout)                       :: dqx(*)
  real(dp),    intent(inout)                       :: dqy(*)
  real(dp),    intent(inout)                       :: dqz(*)
  real(dp),    intent(inout)                       :: dqs(maxlhs,6)
!
!  Local variables
!
  integer(i4)                                      :: j
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  real(dp)                                         :: d1trm
  real(dp)                                         :: dcnidrij
  real(dp)                                         :: dr2ds(6)
  real(dp)                                         :: d2r2dx2(3,3)
  real(dp)                                         :: d2r2ds2(6,6)
  real(dp)                                         :: d2r2dsdx(6,3)
  real(dp)                                         :: dxyz(3)
  real(dp)                                         :: dZdcn
  real(dp)                                         :: sr_cn
#ifdef TRACE
  call trace_in('gfnff_drv_dcn_q')
#endif
!
!  Compute derivatives of i based on it's own coordination number
!
  sr_cn = sqrt(cn(i))
  dZdcn = - 0.5_dp*gfnff_eeq_cnf(i)*dcn(i)/sr_cn
!
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    j = nbrno(ni,i)
    dcnidrij = d1cndr_cn(n,i)
    d1trm = dZdcn*dcnidrij
    dxyz(1) = d1trm*xnbr(ni,i)
    dxyz(2) = d1trm*ynbr(ni,i)
    dxyz(3) = d1trm*znbr(ni,i)
!
    dqx(i) = dqx(i) - dxyz(1)
    dqy(i) = dqy(i) - dxyz(2)
    dqz(i) = dqz(i) - dxyz(3)
!
    if (lstr) then
      call real1strterm(ndim,xnbr(ni,i),ynbr(ni,i),znbr(ni,i),0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
      do kl = 1,nstrains
        ks = nstrptr(kl)
        dqs(i,kl) = dqs(i,kl) + d1trm*dr2ds(ks)
      enddo
    endif
  enddo
!
!  Loop over all atoms to look for derivatives with respect to i
!
  do j = 1,numat
    if (i.eq.j) cycle 
    sr_cn = sqrt(cn(j))
    dZdcn = - 0.5_dp*gfnff_eeq_cnf(j)*dcn(j)/sr_cn
    do n = 1,nnbr_cn(j)
      ni = nbrno_cn(n,j)
      if (nbrno(ni,j).eq.i) then
        dcnidrij = d1cndr_cn(n,j)
        d1trm = dZdcn*dcnidrij
        dxyz(1) = d1trm*xnbr(ni,j)
        dxyz(2) = d1trm*ynbr(ni,j)
        dxyz(3) = d1trm*znbr(ni,j)
!
        dqx(j) = dqx(j) + dxyz(1)
        dqy(j) = dqy(j) + dxyz(2)
        dqz(j) = dqz(j) + dxyz(3)
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv_dcn_q')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnp(i,j,xij,yij,zij,dEdlogcni,dEdlogcnj,d2Edlogcni2,d2Edlogcnj2,d2Edlogcnidlogcnj, &
                             d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn,d2logcndcn2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Phonon version.
!
!  On entry : 
!
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  dEdlogcnj         = derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcnj2       = second derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcnidlogcnj = second derivatives of energy w.r.t. log of the coordination number for i and j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   3/21 Created from gfnff_drv2_dcn
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
!  Julian Gale, CIC, Curtin University, March 2021
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: dEdlogcnj
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: d2Edlogcnidlogcnj
  real(dp),    intent(in)                          :: d2Edlogcnj2
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dcnjdrjl
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2logcnjdcnj2
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2jkdx2(6)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnp')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  d2logcnidcni2 = d2logcndcn2(i)
  d2logcnjdcnj2 = d2logcndcn2(j)
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
  indj = 3*(j-1)
  jx = indj + 1
  jy = indj + 2
  jz = indj + 3
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
!
    d2r2ikdx2(1) = xnbr(ni,i)*xnbr(ni,i)
    d2r2ikdx2(2) = ynbr(ni,i)*ynbr(ni,i)
    d2r2ikdx2(3) = znbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(4) = ynbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(5) = xnbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(6) = xnbr(ni,i)*ynbr(ni,i)
!
    indk = 3*(k-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!
    d2cnidrik2 = d2cndr_cn(n,i)
    d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
           (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
    if (k.ge.i) then
      derv2(kx,ix) = derv2(kx,ix) - d2trm*d2r2ikdx2(1)
      derv2(ky,ix) = derv2(ky,ix) - d2trm*d2r2ikdx2(6)
      derv2(kz,ix) = derv2(kz,ix) - d2trm*d2r2ikdx2(5)
      derv2(kx,iy) = derv2(kx,iy) - d2trm*d2r2ikdx2(6)
      derv2(ky,iy) = derv2(ky,iy) - d2trm*d2r2ikdx2(2)
      derv2(kz,iy) = derv2(kz,iy) - d2trm*d2r2ikdx2(4)
      derv2(kx,iz) = derv2(kx,iz) - d2trm*d2r2ikdx2(5)
      derv2(ky,iz) = derv2(ky,iz) - d2trm*d2r2ikdx2(4)
      derv2(kz,iz) = derv2(kz,iz) - d2trm*d2r2ikdx2(3)
      derv2(kx,ix) = derv2(kx,ix) - d1trm
      derv2(ky,iy) = derv2(ky,iy) - d1trm
      derv2(kz,iz) = derv2(kz,iz) - d1trm
    else
      derv2(ix,kx) = derv2(ix,kx) - d2trm*d2r2ikdx2(1)
      derv2(iy,kx) = derv2(iy,kx) - d2trm*d2r2ikdx2(6)
      derv2(iz,kx) = derv2(iz,kx) - d2trm*d2r2ikdx2(5)
      derv2(ix,ky) = derv2(ix,ky) - d2trm*d2r2ikdx2(6)
      derv2(iy,ky) = derv2(iy,ky) - d2trm*d2r2ikdx2(2)
      derv2(iz,ky) = derv2(iz,ky) - d2trm*d2r2ikdx2(4)
      derv2(ix,kz) = derv2(ix,kz) - d2trm*d2r2ikdx2(5)
      derv2(iy,kz) = derv2(iy,kz) - d2trm*d2r2ikdx2(4)
      derv2(iz,kz) = derv2(iz,kz) - d2trm*d2r2ikdx2(3)
      derv2(ix,kx) = derv2(ix,kx) - d1trm
      derv2(iy,ky) = derv2(iy,ky) - d1trm
      derv2(iz,kz) = derv2(iz,kz) - d1trm
    endif
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
    d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
!
!  I-K
!
    if (k.ge.i) then
      derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xij
      derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xij
      derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xij
      derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*yij
      derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*yij
      derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*yij
      derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*zij
      derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*zij
      derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*zij
    else
      derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xij
      derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*yij
      derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*zij
      derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xij
      derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*yij
      derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*zij
      derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xij
      derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*yij
      derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*zij
    endif
!
!  J-K
!
    if (k.ge.j) then
      derv2(kx,jx) = derv2(kx,jx) + d3trm*xnbr(ni,i)*xij
      derv2(ky,jx) = derv2(ky,jx) + d3trm*ynbr(ni,i)*xij
      derv2(kz,jx) = derv2(kz,jx) + d3trm*znbr(ni,i)*xij
      derv2(kx,jy) = derv2(kx,jy) + d3trm*xnbr(ni,i)*yij
      derv2(ky,jy) = derv2(ky,jy) + d3trm*ynbr(ni,i)*yij
      derv2(kz,jy) = derv2(kz,jy) + d3trm*znbr(ni,i)*yij
      derv2(kx,jz) = derv2(kx,jz) + d3trm*xnbr(ni,i)*zij
      derv2(ky,jz) = derv2(ky,jz) + d3trm*ynbr(ni,i)*zij
      derv2(kz,jz) = derv2(kz,jz) + d3trm*znbr(ni,i)*zij
    else
      derv2(jx,kx) = derv2(jx,kx) + d3trm*xnbr(ni,i)*xij
      derv2(jy,kx) = derv2(jy,kx) + d3trm*xnbr(ni,i)*yij
      derv2(jz,kx) = derv2(jz,kx) + d3trm*xnbr(ni,i)*zij
      derv2(jx,ky) = derv2(jx,ky) + d3trm*ynbr(ni,i)*xij
      derv2(jy,ky) = derv2(jy,ky) + d3trm*ynbr(ni,i)*yij
      derv2(jz,ky) = derv2(jz,ky) + d3trm*ynbr(ni,i)*zij
      derv2(jx,kz) = derv2(jx,kz) + d3trm*znbr(ni,i)*xij
      derv2(jy,kz) = derv2(jy,kz) + d3trm*znbr(ni,i)*yij
      derv2(jz,kz) = derv2(jz,kz) + d3trm*znbr(ni,i)*zij
    endif
!
!  I-J
!
    if (j.ge.i) then
      derv2(jx,ix) = derv2(jx,ix) - d3trm*xnbr(ni,i)*xij
      derv2(jy,ix) = derv2(jy,ix) - d3trm*xnbr(ni,i)*yij
      derv2(jz,ix) = derv2(jz,ix) - d3trm*xnbr(ni,i)*zij
      derv2(jx,iy) = derv2(jx,iy) - d3trm*ynbr(ni,i)*xij
      derv2(jy,iy) = derv2(jy,iy) - d3trm*ynbr(ni,i)*yij
      derv2(jz,iy) = derv2(jz,iy) - d3trm*ynbr(ni,i)*zij
      derv2(jx,iz) = derv2(jx,iz) - d3trm*znbr(ni,i)*xij
      derv2(jy,iz) = derv2(jy,iz) - d3trm*znbr(ni,i)*yij
      derv2(jz,iz) = derv2(jz,iz) - d3trm*znbr(ni,i)*zij
    else
      derv2(ix,jx) = derv2(ix,jx) - d3trm*xnbr(ni,i)*xij
      derv2(iy,jx) = derv2(iy,jx) - d3trm*ynbr(ni,i)*xij
      derv2(iz,jx) = derv2(iz,jx) - d3trm*znbr(ni,i)*xij
      derv2(ix,jy) = derv2(ix,jy) - d3trm*xnbr(ni,i)*yij
      derv2(iy,jy) = derv2(iy,jy) - d3trm*ynbr(ni,i)*yij
      derv2(iz,jy) = derv2(iz,jy) - d3trm*znbr(ni,i)*yij
      derv2(ix,jz) = derv2(ix,jz) - d3trm*xnbr(ni,i)*zij
      derv2(iy,jz) = derv2(iy,jz) - d3trm*ynbr(ni,i)*zij
      derv2(iz,jz) = derv2(iz,jz) - d3trm*znbr(ni,i)*zij
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      ni2 = nbrno_cn(n2,i)
      l = nbrno(ni2,i)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnidril = d1cndr_cn(n2,i)
      d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  I-K
!
      if (k.ge.i) then
        derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      else
        derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  I-L
!
      if (l.ge.i) then
        derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ly,ix) = derv2(ly,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lz,ix) = derv2(lz,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lx,iy) = derv2(lx,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lz,iy) = derv2(lz,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lx,iz) = derv2(lx,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ly,iz) = derv2(ly,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      else
        derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iy,lx) = derv2(iy,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(iz,lx) = derv2(iz,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ix,ly) = derv2(ix,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(iz,ly) = derv2(iz,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(ix,lz) = derv2(ix,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(iy,lz) = derv2(iy,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  K-L
!
      if (l.ge.k) then
        derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      else
        derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
    enddo
!------------------------------------------------------------------------------------------
!  Loop over neighbours of j for second derivatives for 2 different coordination numbers  |
!------------------------------------------------------------------------------------------
    do n2 = 1,nnbr_cn(j)
      nj = nbrno_cn(n2,j)
      l = nbrno(nj,j)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = d2Edlogcnidlogcnj*dlogcnidcni*dlogcnjdcnj*dcnidrik*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  I-J
!
      if (j.ge.i) then
        derv2(jx,ix) = derv2(jx,ix) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jy,ix) = derv2(jy,ix) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jz,ix) = derv2(jz,ix) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jx,iy) = derv2(jx,iy) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jy,iy) = derv2(jy,iy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jz,iy) = derv2(jz,iy) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jx,iz) = derv2(jx,iz) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jy,iz) = derv2(jy,iz) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jz,iz) = derv2(jz,iz) + d3trm*znbr(nj,j)*znbr(ni,i)
      else
        derv2(ix,jx) = derv2(ix,jx) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iy,jx) = derv2(iy,jx) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(iz,jx) = derv2(iz,jx) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ix,jy) = derv2(ix,jy) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iy,jy) = derv2(iy,jy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(iz,jy) = derv2(iz,jy) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ix,jz) = derv2(ix,jz) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iy,jz) = derv2(iy,jz) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(iz,jz) = derv2(iz,jz) + d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  I-L
!
      if (l.ge.i) then
        derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(ly,ix) = derv2(ly,ix) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(lz,ix) = derv2(lz,ix) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(lx,iy) = derv2(lx,iy) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(lz,iy) = derv2(lz,iy) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(lx,iz) = derv2(lx,iz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ly,iz) = derv2(ly,iz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(nj,j)*znbr(ni,i)
      else
        derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iy,lx) = derv2(iy,lx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(iz,lx) = derv2(iz,lx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ix,ly) = derv2(ix,ly) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(iz,ly) = derv2(iz,ly) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ix,lz) = derv2(ix,lz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iy,lz) = derv2(iy,lz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  J-K
!
      if (k.ge.j) then
        derv2(kx,jx) = derv2(kx,jx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(ky,jx) = derv2(ky,jx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(kz,jx) = derv2(kz,jx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(kx,jy) = derv2(kx,jy) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(ky,jy) = derv2(ky,jy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(kz,jy) = derv2(kz,jy) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(kx,jz) = derv2(kx,jz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(ky,jz) = derv2(ky,jz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(kz,jz) = derv2(kz,jz) - d3trm*znbr(nj,j)*znbr(ni,i)
      else
        derv2(jx,kx) = derv2(jx,kx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jy,kx) = derv2(jy,kx) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jz,kx) = derv2(jz,kx) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jx,ky) = derv2(jx,ky) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jy,ky) = derv2(jy,ky) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jz,ky) = derv2(jz,ky) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jx,kz) = derv2(jx,kz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jy,kz) = derv2(jy,kz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jz,kz) = derv2(jz,kz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  K-L
!     
      if (l.ge.k) then
        derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(nj,j)
      else
        derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(nj,j)
      endif
    enddo
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    dcnjdrjk = d1cndr_cn(n,j)
    d1trm = dlogcnjdcnj*dcnjdrjk*dEdlogcnj
!
    d2r2jkdx2(1) = xnbr(nj,j)*xnbr(nj,j)
    d2r2jkdx2(2) = ynbr(nj,j)*ynbr(nj,j)
    d2r2jkdx2(3) = znbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(4) = ynbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(5) = xnbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(6) = xnbr(nj,j)*ynbr(nj,j)
!
    indk = 3*(k-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!
    d2cnjdrjk2 = d2cndr_cn(n,j)
    d2trm = dEdlogcnj*dlogcnjdcnj*d2cnjdrjk2 + &
           (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjk
!----------------------------------------
!  Second derivatives of CN w.r.t. rjk  |
!----------------------------------------
    if (k.ge.j) then
      derv2(kx,jx) = derv2(kx,jx) - d2trm*d2r2jkdx2(1)
      derv2(ky,jx) = derv2(ky,jx) - d2trm*d2r2jkdx2(6)
      derv2(kz,jx) = derv2(kz,jx) - d2trm*d2r2jkdx2(5)
      derv2(kx,jy) = derv2(kx,jy) - d2trm*d2r2jkdx2(6)
      derv2(ky,jy) = derv2(ky,jy) - d2trm*d2r2jkdx2(2)
      derv2(kz,jy) = derv2(kz,jy) - d2trm*d2r2jkdx2(4)
      derv2(kx,jz) = derv2(kx,jz) - d2trm*d2r2jkdx2(5)
      derv2(ky,jz) = derv2(ky,jz) - d2trm*d2r2jkdx2(4)
      derv2(kz,jz) = derv2(kz,jz) - d2trm*d2r2jkdx2(3)
      derv2(kx,jx) = derv2(kx,jx) - d1trm
      derv2(ky,jy) = derv2(ky,jy) - d1trm
      derv2(kz,jz) = derv2(kz,jz) - d1trm
    else
      derv2(jx,kx) = derv2(jx,kx) - d2trm*d2r2jkdx2(1)
      derv2(jy,kx) = derv2(jy,kx) - d2trm*d2r2jkdx2(6)
      derv2(jz,kx) = derv2(jz,kx) - d2trm*d2r2jkdx2(5)
      derv2(jx,ky) = derv2(jx,ky) - d2trm*d2r2jkdx2(6)
      derv2(jy,ky) = derv2(jy,ky) - d2trm*d2r2jkdx2(2)
      derv2(jz,ky) = derv2(jz,ky) - d2trm*d2r2jkdx2(4)
      derv2(jx,kz) = derv2(jx,kz) - d2trm*d2r2jkdx2(5)
      derv2(jy,kz) = derv2(jy,kz) - d2trm*d2r2jkdx2(4)
      derv2(jz,kz) = derv2(jz,kz) - d2trm*d2r2jkdx2(3)
      derv2(jx,kx) = derv2(jx,kx) - d1trm
      derv2(jy,ky) = derv2(jy,ky) - d1trm
      derv2(jz,kz) = derv2(jz,kz) - d1trm
    endif
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
    d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
!
!  J-K
!
    if (k.ge.j) then
      derv2(kx,jx) = derv2(kx,jx) + d3trm*xnbr(nj,j)*xij
      derv2(ky,jx) = derv2(ky,jx) + d3trm*ynbr(nj,j)*xij
      derv2(kz,jx) = derv2(kz,jx) + d3trm*znbr(nj,j)*xij
      derv2(kx,jy) = derv2(kx,jy) + d3trm*xnbr(nj,j)*yij
      derv2(ky,jy) = derv2(ky,jy) + d3trm*ynbr(nj,j)*yij
      derv2(kz,jy) = derv2(kz,jy) + d3trm*znbr(nj,j)*yij
      derv2(kx,jz) = derv2(kx,jz) + d3trm*xnbr(nj,j)*zij
      derv2(ky,jz) = derv2(ky,jz) + d3trm*ynbr(nj,j)*zij
      derv2(kz,jz) = derv2(kz,jz) + d3trm*znbr(nj,j)*zij
    else
      derv2(jx,kx) = derv2(jx,kx) + d3trm*xnbr(nj,j)*xij
      derv2(jy,kx) = derv2(jy,kx) + d3trm*xnbr(nj,j)*yij
      derv2(jz,kx) = derv2(jz,kx) + d3trm*xnbr(nj,j)*zij
      derv2(jx,ky) = derv2(jx,ky) + d3trm*ynbr(nj,j)*xij
      derv2(jy,ky) = derv2(jy,ky) + d3trm*ynbr(nj,j)*yij
      derv2(jz,ky) = derv2(jz,ky) + d3trm*ynbr(nj,j)*zij
      derv2(jx,kz) = derv2(jx,kz) + d3trm*znbr(nj,j)*xij
      derv2(jy,kz) = derv2(jy,kz) + d3trm*znbr(nj,j)*yij
      derv2(jz,kz) = derv2(jz,kz) + d3trm*znbr(nj,j)*zij
    endif
!
!  I-K
!
    if (k.ge.i) then
      derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(nj,j)*xij
      derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(nj,j)*xij
      derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(nj,j)*xij
      derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(nj,j)*yij
      derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(nj,j)*yij
      derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(nj,j)*yij
      derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(nj,j)*zij
      derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(nj,j)*zij
      derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(nj,j)*zij
    else
      derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(nj,j)*xij
      derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(nj,j)*yij
      derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(nj,j)*zij
      derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(nj,j)*xij
      derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(nj,j)*yij
      derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(nj,j)*zij
      derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(nj,j)*xij
      derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(nj,j)*yij
      derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(nj,j)*zij
    endif
!
!  J-I
!
    if (i.ge.j) then
      derv2(ix,jx) = derv2(ix,jx) + d3trm*xnbr(nj,j)*xij
      derv2(iy,jx) = derv2(iy,jx) + d3trm*xnbr(nj,j)*yij
      derv2(iz,jx) = derv2(iz,jx) + d3trm*xnbr(nj,j)*zij
      derv2(ix,jy) = derv2(ix,jy) + d3trm*ynbr(nj,j)*xij
      derv2(iy,jy) = derv2(iy,jy) + d3trm*ynbr(nj,j)*yij
      derv2(iz,jy) = derv2(iz,jy) + d3trm*ynbr(nj,j)*zij
      derv2(ix,jz) = derv2(ix,jz) + d3trm*znbr(nj,j)*xij
      derv2(iy,jz) = derv2(iy,jz) + d3trm*znbr(nj,j)*yij
      derv2(iz,jz) = derv2(iz,jz) + d3trm*znbr(nj,j)*zij
    else
      derv2(jx,ix) = derv2(jx,ix) + d3trm*xnbr(nj,j)*xij
      derv2(jy,ix) = derv2(jy,ix) + d3trm*ynbr(nj,j)*xij
      derv2(jz,ix) = derv2(jz,ix) + d3trm*znbr(nj,j)*xij
      derv2(jx,iy) = derv2(jx,iy) + d3trm*xnbr(nj,j)*yij
      derv2(jy,iy) = derv2(jy,iy) + d3trm*ynbr(nj,j)*yij
      derv2(jz,iy) = derv2(jz,iy) + d3trm*znbr(nj,j)*yij
      derv2(jx,iz) = derv2(jx,iz) + d3trm*xnbr(nj,j)*zij
      derv2(jy,iz) = derv2(jy,iz) + d3trm*ynbr(nj,j)*zij
      derv2(jz,iz) = derv2(jz,iz) + d3trm*znbr(nj,j)*zij
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of j for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      nj2 = nbrno_cn(n2,j)
      l = nbrno(nj2,j)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  J-K
!
      if (k.ge.j) then
        derv2(kx,jx) = derv2(kx,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(ky,jx) = derv2(ky,jx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kz,jx) = derv2(kz,jx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kx,jy) = derv2(kx,jy) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(ky,jy) = derv2(ky,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kz,jy) = derv2(kz,jy) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kx,jz) = derv2(kx,jz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(ky,jz) = derv2(ky,jz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kz,jz) = derv2(kz,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      else
        derv2(jx,kx) = derv2(jx,kx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jy,kx) = derv2(jy,kx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jz,kx) = derv2(jz,kx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jx,ky) = derv2(jx,ky) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jy,ky) = derv2(jy,ky) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jz,ky) = derv2(jz,ky) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jx,kz) = derv2(jx,kz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jy,kz) = derv2(jy,kz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jz,kz) = derv2(jz,kz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  J-L
!
      if (l.ge.j) then
        derv2(lx,jx) = derv2(lx,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(ly,jx) = derv2(ly,jx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lz,jx) = derv2(lz,jx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lx,jy) = derv2(lx,jy) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(ly,jy) = derv2(ly,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lz,jy) = derv2(lz,jy) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lx,jz) = derv2(lx,jz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(ly,jz) = derv2(ly,jz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lz,jz) = derv2(lz,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      else
        derv2(jx,lx) = derv2(jx,lx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jy,lx) = derv2(jy,lx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jz,lx) = derv2(jz,lx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jx,ly) = derv2(jx,ly) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jy,ly) = derv2(jy,ly) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jz,ly) = derv2(jz,ly) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jx,lz) = derv2(jx,lz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jy,lz) = derv2(jy,lz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jz,lz) = derv2(jz,lz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  K-L
!     
      if (l.ge.k) then
        derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      else
        derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnp')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnp_nox(i,j,dEdlogcni,dEdlogcnj,d2Edlogcni2,d2Edlogcnj2,d2Edlogcnidlogcnj, &
                                 dlogcndcn,d2logcndcn2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Phonon version. Based on gfnff_drv2_dcnp but excludes cross terms for distance/CN
!
!  On entry : 
!
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  dEdlogcnj         = derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcnj2       = second derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcnidlogcnj = second derivatives of energy w.r.t. log of the coordination number for i and j
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   5/22 Created from gfnff_drv2_dcnp
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: dEdlogcnj
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: d2Edlogcnidlogcnj
  real(dp),    intent(in)                          :: d2Edlogcnj2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dcnjdrjl
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2logcnjdcnj2
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2jkdx2(6)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnp')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  d2logcnidcni2 = d2logcndcn2(i)
  d2logcnjdcnj2 = d2logcndcn2(j)
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
  indj = 3*(j-1)
  jx = indj + 1
  jy = indj + 2
  jz = indj + 3
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
!
    d2r2ikdx2(1) = xnbr(ni,i)*xnbr(ni,i)
    d2r2ikdx2(2) = ynbr(ni,i)*ynbr(ni,i)
    d2r2ikdx2(3) = znbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(4) = ynbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(5) = xnbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(6) = xnbr(ni,i)*ynbr(ni,i)
!
    indk = 3*(k-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!
    d2cnidrik2 = d2cndr_cn(n,i)
    d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
           (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
    if (k.ge.i) then
      derv2(kx,ix) = derv2(kx,ix) - d2trm*d2r2ikdx2(1)
      derv2(ky,ix) = derv2(ky,ix) - d2trm*d2r2ikdx2(6)
      derv2(kz,ix) = derv2(kz,ix) - d2trm*d2r2ikdx2(5)
      derv2(kx,iy) = derv2(kx,iy) - d2trm*d2r2ikdx2(6)
      derv2(ky,iy) = derv2(ky,iy) - d2trm*d2r2ikdx2(2)
      derv2(kz,iy) = derv2(kz,iy) - d2trm*d2r2ikdx2(4)
      derv2(kx,iz) = derv2(kx,iz) - d2trm*d2r2ikdx2(5)
      derv2(ky,iz) = derv2(ky,iz) - d2trm*d2r2ikdx2(4)
      derv2(kz,iz) = derv2(kz,iz) - d2trm*d2r2ikdx2(3)
      derv2(kx,ix) = derv2(kx,ix) - d1trm
      derv2(ky,iy) = derv2(ky,iy) - d1trm
      derv2(kz,iz) = derv2(kz,iz) - d1trm
    else
      derv2(ix,kx) = derv2(ix,kx) - d2trm*d2r2ikdx2(1)
      derv2(iy,kx) = derv2(iy,kx) - d2trm*d2r2ikdx2(6)
      derv2(iz,kx) = derv2(iz,kx) - d2trm*d2r2ikdx2(5)
      derv2(ix,ky) = derv2(ix,ky) - d2trm*d2r2ikdx2(6)
      derv2(iy,ky) = derv2(iy,ky) - d2trm*d2r2ikdx2(2)
      derv2(iz,ky) = derv2(iz,ky) - d2trm*d2r2ikdx2(4)
      derv2(ix,kz) = derv2(ix,kz) - d2trm*d2r2ikdx2(5)
      derv2(iy,kz) = derv2(iy,kz) - d2trm*d2r2ikdx2(4)
      derv2(iz,kz) = derv2(iz,kz) - d2trm*d2r2ikdx2(3)
      derv2(ix,kx) = derv2(ix,kx) - d1trm
      derv2(iy,ky) = derv2(iy,ky) - d1trm
      derv2(iz,kz) = derv2(iz,kz) - d1trm
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      ni2 = nbrno_cn(n2,i)
      l = nbrno(ni2,i)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnidril = d1cndr_cn(n2,i)
      d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  I-K
!
      if (k.ge.i) then
        derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      else
        derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  I-L
!
      if (l.ge.i) then
        derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ly,ix) = derv2(ly,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lz,ix) = derv2(lz,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lx,iy) = derv2(lx,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lz,iy) = derv2(lz,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lx,iz) = derv2(lx,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ly,iz) = derv2(ly,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      else
        derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iy,lx) = derv2(iy,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(iz,lx) = derv2(iz,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ix,ly) = derv2(ix,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(iz,ly) = derv2(iz,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(ix,lz) = derv2(ix,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(iy,lz) = derv2(iy,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  K-L
!
      if (l.ge.k) then
        derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      else
        derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
    enddo
!------------------------------------------------------------------------------------------
!  Loop over neighbours of j for second derivatives for 2 different coordination numbers  |
!------------------------------------------------------------------------------------------
    do n2 = 1,nnbr_cn(j)
      nj = nbrno_cn(n2,j)
      l = nbrno(nj,j)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = d2Edlogcnidlogcnj*dlogcnidcni*dlogcnjdcnj*dcnidrik*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  I-J
!
      if (j.ge.i) then
        derv2(jx,ix) = derv2(jx,ix) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jy,ix) = derv2(jy,ix) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jz,ix) = derv2(jz,ix) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jx,iy) = derv2(jx,iy) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jy,iy) = derv2(jy,iy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jz,iy) = derv2(jz,iy) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jx,iz) = derv2(jx,iz) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jy,iz) = derv2(jy,iz) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jz,iz) = derv2(jz,iz) + d3trm*znbr(nj,j)*znbr(ni,i)
      else
        derv2(ix,jx) = derv2(ix,jx) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iy,jx) = derv2(iy,jx) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(iz,jx) = derv2(iz,jx) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ix,jy) = derv2(ix,jy) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iy,jy) = derv2(iy,jy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(iz,jy) = derv2(iz,jy) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ix,jz) = derv2(ix,jz) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iy,jz) = derv2(iy,jz) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(iz,jz) = derv2(iz,jz) + d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  I-L
!
      if (l.ge.i) then
        derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(ly,ix) = derv2(ly,ix) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(lz,ix) = derv2(lz,ix) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(lx,iy) = derv2(lx,iy) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(lz,iy) = derv2(lz,iy) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(lx,iz) = derv2(lx,iz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ly,iz) = derv2(ly,iz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(nj,j)*znbr(ni,i)
      else
        derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iy,lx) = derv2(iy,lx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(iz,lx) = derv2(iz,lx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ix,ly) = derv2(ix,ly) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(iz,ly) = derv2(iz,ly) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ix,lz) = derv2(ix,lz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iy,lz) = derv2(iy,lz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  J-K
!
      if (k.ge.j) then
        derv2(kx,jx) = derv2(kx,jx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(ky,jx) = derv2(ky,jx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(kz,jx) = derv2(kz,jx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(kx,jy) = derv2(kx,jy) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(ky,jy) = derv2(ky,jy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(kz,jy) = derv2(kz,jy) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(kx,jz) = derv2(kx,jz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(ky,jz) = derv2(ky,jz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(kz,jz) = derv2(kz,jz) - d3trm*znbr(nj,j)*znbr(ni,i)
      else
        derv2(jx,kx) = derv2(jx,kx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jy,kx) = derv2(jy,kx) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jz,kx) = derv2(jz,kx) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jx,ky) = derv2(jx,ky) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jy,ky) = derv2(jy,ky) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jz,ky) = derv2(jz,ky) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jx,kz) = derv2(jx,kz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jy,kz) = derv2(jy,kz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jz,kz) = derv2(jz,kz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  K-L
!     
      if (l.ge.k) then
        derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(nj,j)
      else
        derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(nj,j)
      endif
    enddo
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    dcnjdrjk = d1cndr_cn(n,j)
    d1trm = dlogcnjdcnj*dcnjdrjk*dEdlogcnj
!
    d2r2jkdx2(1) = xnbr(nj,j)*xnbr(nj,j)
    d2r2jkdx2(2) = ynbr(nj,j)*ynbr(nj,j)
    d2r2jkdx2(3) = znbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(4) = ynbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(5) = xnbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(6) = xnbr(nj,j)*ynbr(nj,j)
!
    indk = 3*(k-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!
    d2cnjdrjk2 = d2cndr_cn(n,j)
    d2trm = dEdlogcnj*dlogcnjdcnj*d2cnjdrjk2 + &
           (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjk
!----------------------------------------
!  Second derivatives of CN w.r.t. rjk  |
!----------------------------------------
    if (k.ge.j) then
      derv2(kx,jx) = derv2(kx,jx) - d2trm*d2r2jkdx2(1)
      derv2(ky,jx) = derv2(ky,jx) - d2trm*d2r2jkdx2(6)
      derv2(kz,jx) = derv2(kz,jx) - d2trm*d2r2jkdx2(5)
      derv2(kx,jy) = derv2(kx,jy) - d2trm*d2r2jkdx2(6)
      derv2(ky,jy) = derv2(ky,jy) - d2trm*d2r2jkdx2(2)
      derv2(kz,jy) = derv2(kz,jy) - d2trm*d2r2jkdx2(4)
      derv2(kx,jz) = derv2(kx,jz) - d2trm*d2r2jkdx2(5)
      derv2(ky,jz) = derv2(ky,jz) - d2trm*d2r2jkdx2(4)
      derv2(kz,jz) = derv2(kz,jz) - d2trm*d2r2jkdx2(3)
      derv2(kx,jx) = derv2(kx,jx) - d1trm
      derv2(ky,jy) = derv2(ky,jy) - d1trm
      derv2(kz,jz) = derv2(kz,jz) - d1trm
    else
      derv2(jx,kx) = derv2(jx,kx) - d2trm*d2r2jkdx2(1)
      derv2(jy,kx) = derv2(jy,kx) - d2trm*d2r2jkdx2(6)
      derv2(jz,kx) = derv2(jz,kx) - d2trm*d2r2jkdx2(5)
      derv2(jx,ky) = derv2(jx,ky) - d2trm*d2r2jkdx2(6)
      derv2(jy,ky) = derv2(jy,ky) - d2trm*d2r2jkdx2(2)
      derv2(jz,ky) = derv2(jz,ky) - d2trm*d2r2jkdx2(4)
      derv2(jx,kz) = derv2(jx,kz) - d2trm*d2r2jkdx2(5)
      derv2(jy,kz) = derv2(jy,kz) - d2trm*d2r2jkdx2(4)
      derv2(jz,kz) = derv2(jz,kz) - d2trm*d2r2jkdx2(3)
      derv2(jx,kx) = derv2(jx,kx) - d1trm
      derv2(jy,ky) = derv2(jy,ky) - d1trm
      derv2(jz,kz) = derv2(jz,kz) - d1trm
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of j for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      nj2 = nbrno_cn(n2,j)
      l = nbrno(nj2,j)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  J-K
!
      if (k.ge.j) then
        derv2(kx,jx) = derv2(kx,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(ky,jx) = derv2(ky,jx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kz,jx) = derv2(kz,jx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kx,jy) = derv2(kx,jy) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(ky,jy) = derv2(ky,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kz,jy) = derv2(kz,jy) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kx,jz) = derv2(kx,jz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(ky,jz) = derv2(ky,jz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kz,jz) = derv2(kz,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      else
        derv2(jx,kx) = derv2(jx,kx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jy,kx) = derv2(jy,kx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jz,kx) = derv2(jz,kx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jx,ky) = derv2(jx,ky) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jy,ky) = derv2(jy,ky) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jz,ky) = derv2(jz,ky) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jx,kz) = derv2(jx,kz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jy,kz) = derv2(jy,kz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jz,kz) = derv2(jz,kz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  J-L
!
      if (l.ge.j) then
        derv2(lx,jx) = derv2(lx,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(ly,jx) = derv2(ly,jx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lz,jx) = derv2(lz,jx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lx,jy) = derv2(lx,jy) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(ly,jy) = derv2(ly,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lz,jy) = derv2(lz,jy) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lx,jz) = derv2(lx,jz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(ly,jz) = derv2(ly,jz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lz,jz) = derv2(lz,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      else
        derv2(jx,lx) = derv2(jx,lx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jy,lx) = derv2(jy,lx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jz,lx) = derv2(jz,lx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jx,ly) = derv2(jx,ly) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jy,ly) = derv2(jy,ly) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jz,ly) = derv2(jz,ly) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jx,lz) = derv2(jx,lz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jy,lz) = derv2(jy,lz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jz,lz) = derv2(jz,lz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  K-L
!     
      if (l.ge.k) then
        derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      else
        derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnp_nox')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnp_x(i,j,xij,yij,zij,d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Phonon version. Based on gfnff_drv2_dcnp but only does the cross terms for distance/CN
!
!  On entry : 
!
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   5/22 Created from gfnff_drv2_dcnp
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnp_x')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
  indj = 3*(j-1)
  jx = indj + 1
  jy = indj + 2
  jz = indj + 3
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
!
    indk = 3*(k-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
    d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
!
!  I-K
!
    if (k.ge.i) then
      derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xij
      derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xij
      derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xij
      derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*yij
      derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*yij
      derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*yij
      derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*zij
      derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*zij
      derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*zij
    else
      derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xij
      derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*yij
      derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*zij
      derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xij
      derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*yij
      derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*zij
      derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xij
      derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*yij
      derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*zij
    endif
!
!  J-K
!
    if (k.ge.j) then
      derv2(kx,jx) = derv2(kx,jx) + d3trm*xnbr(ni,i)*xij
      derv2(ky,jx) = derv2(ky,jx) + d3trm*ynbr(ni,i)*xij
      derv2(kz,jx) = derv2(kz,jx) + d3trm*znbr(ni,i)*xij
      derv2(kx,jy) = derv2(kx,jy) + d3trm*xnbr(ni,i)*yij
      derv2(ky,jy) = derv2(ky,jy) + d3trm*ynbr(ni,i)*yij
      derv2(kz,jy) = derv2(kz,jy) + d3trm*znbr(ni,i)*yij
      derv2(kx,jz) = derv2(kx,jz) + d3trm*xnbr(ni,i)*zij
      derv2(ky,jz) = derv2(ky,jz) + d3trm*ynbr(ni,i)*zij
      derv2(kz,jz) = derv2(kz,jz) + d3trm*znbr(ni,i)*zij
    else
      derv2(jx,kx) = derv2(jx,kx) + d3trm*xnbr(ni,i)*xij
      derv2(jy,kx) = derv2(jy,kx) + d3trm*xnbr(ni,i)*yij
      derv2(jz,kx) = derv2(jz,kx) + d3trm*xnbr(ni,i)*zij
      derv2(jx,ky) = derv2(jx,ky) + d3trm*ynbr(ni,i)*xij
      derv2(jy,ky) = derv2(jy,ky) + d3trm*ynbr(ni,i)*yij
      derv2(jz,ky) = derv2(jz,ky) + d3trm*ynbr(ni,i)*zij
      derv2(jx,kz) = derv2(jx,kz) + d3trm*znbr(ni,i)*xij
      derv2(jy,kz) = derv2(jy,kz) + d3trm*znbr(ni,i)*yij
      derv2(jz,kz) = derv2(jz,kz) + d3trm*znbr(ni,i)*zij
    endif
!
!  I-J
!
    if (j.ge.i) then
      derv2(jx,ix) = derv2(jx,ix) - d3trm*xnbr(ni,i)*xij
      derv2(jy,ix) = derv2(jy,ix) - d3trm*xnbr(ni,i)*yij
      derv2(jz,ix) = derv2(jz,ix) - d3trm*xnbr(ni,i)*zij
      derv2(jx,iy) = derv2(jx,iy) - d3trm*ynbr(ni,i)*xij
      derv2(jy,iy) = derv2(jy,iy) - d3trm*ynbr(ni,i)*yij
      derv2(jz,iy) = derv2(jz,iy) - d3trm*ynbr(ni,i)*zij
      derv2(jx,iz) = derv2(jx,iz) - d3trm*znbr(ni,i)*xij
      derv2(jy,iz) = derv2(jy,iz) - d3trm*znbr(ni,i)*yij
      derv2(jz,iz) = derv2(jz,iz) - d3trm*znbr(ni,i)*zij
    else
      derv2(ix,jx) = derv2(ix,jx) - d3trm*xnbr(ni,i)*xij
      derv2(iy,jx) = derv2(iy,jx) - d3trm*ynbr(ni,i)*xij
      derv2(iz,jx) = derv2(iz,jx) - d3trm*znbr(ni,i)*xij
      derv2(ix,jy) = derv2(ix,jy) - d3trm*xnbr(ni,i)*yij
      derv2(iy,jy) = derv2(iy,jy) - d3trm*ynbr(ni,i)*yij
      derv2(iz,jy) = derv2(iz,jy) - d3trm*znbr(ni,i)*yij
      derv2(ix,jz) = derv2(ix,jz) - d3trm*xnbr(ni,i)*zij
      derv2(iy,jz) = derv2(iy,jz) - d3trm*ynbr(ni,i)*zij
      derv2(iz,jz) = derv2(iz,jz) - d3trm*znbr(ni,i)*zij
    endif
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
    dcnjdrjk = d1cndr_cn(n,j)
!
    indk = 3*(k-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
    d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
!
!  J-K
!
    if (k.ge.j) then
      derv2(kx,jx) = derv2(kx,jx) + d3trm*xnbr(nj,j)*xij
      derv2(ky,jx) = derv2(ky,jx) + d3trm*ynbr(nj,j)*xij
      derv2(kz,jx) = derv2(kz,jx) + d3trm*znbr(nj,j)*xij
      derv2(kx,jy) = derv2(kx,jy) + d3trm*xnbr(nj,j)*yij
      derv2(ky,jy) = derv2(ky,jy) + d3trm*ynbr(nj,j)*yij
      derv2(kz,jy) = derv2(kz,jy) + d3trm*znbr(nj,j)*yij
      derv2(kx,jz) = derv2(kx,jz) + d3trm*xnbr(nj,j)*zij
      derv2(ky,jz) = derv2(ky,jz) + d3trm*ynbr(nj,j)*zij
      derv2(kz,jz) = derv2(kz,jz) + d3trm*znbr(nj,j)*zij
    else
      derv2(jx,kx) = derv2(jx,kx) + d3trm*xnbr(nj,j)*xij
      derv2(jy,kx) = derv2(jy,kx) + d3trm*xnbr(nj,j)*yij
      derv2(jz,kx) = derv2(jz,kx) + d3trm*xnbr(nj,j)*zij
      derv2(jx,ky) = derv2(jx,ky) + d3trm*ynbr(nj,j)*xij
      derv2(jy,ky) = derv2(jy,ky) + d3trm*ynbr(nj,j)*yij
      derv2(jz,ky) = derv2(jz,ky) + d3trm*ynbr(nj,j)*zij
      derv2(jx,kz) = derv2(jx,kz) + d3trm*znbr(nj,j)*xij
      derv2(jy,kz) = derv2(jy,kz) + d3trm*znbr(nj,j)*yij
      derv2(jz,kz) = derv2(jz,kz) + d3trm*znbr(nj,j)*zij
    endif
!
!  I-K
!
    if (k.ge.i) then
      derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(nj,j)*xij
      derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(nj,j)*xij
      derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(nj,j)*xij
      derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(nj,j)*yij
      derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(nj,j)*yij
      derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(nj,j)*yij
      derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(nj,j)*zij
      derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(nj,j)*zij
      derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(nj,j)*zij
    else
      derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(nj,j)*xij
      derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(nj,j)*yij
      derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(nj,j)*zij
      derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(nj,j)*xij
      derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(nj,j)*yij
      derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(nj,j)*zij
      derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(nj,j)*xij
      derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(nj,j)*yij
      derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(nj,j)*zij
    endif
!
!  J-I
!
    if (i.ge.j) then
      derv2(ix,jx) = derv2(ix,jx) + d3trm*xnbr(nj,j)*xij
      derv2(iy,jx) = derv2(iy,jx) + d3trm*xnbr(nj,j)*yij
      derv2(iz,jx) = derv2(iz,jx) + d3trm*xnbr(nj,j)*zij
      derv2(ix,jy) = derv2(ix,jy) + d3trm*ynbr(nj,j)*xij
      derv2(iy,jy) = derv2(iy,jy) + d3trm*ynbr(nj,j)*yij
      derv2(iz,jy) = derv2(iz,jy) + d3trm*ynbr(nj,j)*zij
      derv2(ix,jz) = derv2(ix,jz) + d3trm*znbr(nj,j)*xij
      derv2(iy,jz) = derv2(iy,jz) + d3trm*znbr(nj,j)*yij
      derv2(iz,jz) = derv2(iz,jz) + d3trm*znbr(nj,j)*zij
    else
      derv2(jx,ix) = derv2(jx,ix) + d3trm*xnbr(nj,j)*xij
      derv2(jy,ix) = derv2(jy,ix) + d3trm*ynbr(nj,j)*xij
      derv2(jz,ix) = derv2(jz,ix) + d3trm*znbr(nj,j)*xij
      derv2(jx,iy) = derv2(jx,iy) + d3trm*xnbr(nj,j)*yij
      derv2(jy,iy) = derv2(jy,iy) + d3trm*ynbr(nj,j)*yij
      derv2(jz,iy) = derv2(jz,iy) + d3trm*znbr(nj,j)*yij
      derv2(jx,iz) = derv2(jx,iz) + d3trm*xnbr(nj,j)*zij
      derv2(jy,iz) = derv2(jy,iz) + d3trm*ynbr(nj,j)*zij
      derv2(jz,iz) = derv2(jz,iz) + d3trm*znbr(nj,j)*zij
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnp_x')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnpd(lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,j,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                              dEdlogcni,dEdlogcnj,d2Edlogcni2,d2Edlogcnj2,d2Edlogcnidlogcnj, &
                              d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn,d2logcndcn2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Phonon version. Distributed memory parallel version.
!
!  On entry : 
!
!  ix,  iy,  iz      = second derivative Cartesian elements for i (local)
!  ixf, iyf, izf     = second derivative Cartesian elements for i (global)
!  jx,  jy,  jz      = second derivative Cartesian elements for j (local)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  dEdlogcnj         = derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcnj2       = second derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcnidlogcnj = second derivatives of energy w.r.t. log of the coordination number for i and j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   4/22 Created from gfnff_drv2_dcnp
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
!  Julian Gale, CIC, Curtin University, April 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: dEdlogcnj
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: d2Edlogcnidlogcnj
  real(dp),    intent(in)                          :: d2Edlogcnj2
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
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
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  logical                                          :: lijklocal
  logical                                          :: lklocal
  logical                                          :: lllocal
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dcnjdrjl
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2logcnjdcnj2
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2jkdx2(6)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnpd')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  d2logcnidcni2 = d2logcndcn2(i)
  d2logcnjdcnj2 = d2logcndcn2(j)
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
!
    d2r2ikdx2(1) = xnbr(ni,i)*xnbr(ni,i)
    d2r2ikdx2(2) = ynbr(ni,i)*ynbr(ni,i)
    d2r2ikdx2(3) = znbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(4) = ynbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(5) = xnbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(6) = xnbr(ni,i)*ynbr(ni,i)
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
    d2cnidrik2 = d2cndr_cn(n,i)
    d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
           (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
    if (lilocal) then
      derv2(kxf,ix) = derv2(kxf,ix) - d2trm*d2r2ikdx2(1)
      derv2(kyf,ix) = derv2(kyf,ix) - d2trm*d2r2ikdx2(6)
      derv2(kzf,ix) = derv2(kzf,ix) - d2trm*d2r2ikdx2(5)
      derv2(kxf,iy) = derv2(kxf,iy) - d2trm*d2r2ikdx2(6)
      derv2(kyf,iy) = derv2(kyf,iy) - d2trm*d2r2ikdx2(2)
      derv2(kzf,iy) = derv2(kzf,iy) - d2trm*d2r2ikdx2(4)
      derv2(kxf,iz) = derv2(kxf,iz) - d2trm*d2r2ikdx2(5)
      derv2(kyf,iz) = derv2(kyf,iz) - d2trm*d2r2ikdx2(4)
      derv2(kzf,iz) = derv2(kzf,iz) - d2trm*d2r2ikdx2(3)
      derv2(kxf,ix) = derv2(kxf,ix) - d1trm
      derv2(kyf,iy) = derv2(kyf,iy) - d1trm
      derv2(kzf,iz) = derv2(kzf,iz) - d1trm
    endif
!
    if (lklocal) then
      derv2(ixf,kx) = derv2(ixf,kx) - d2trm*d2r2ikdx2(1)
      derv2(iyf,kx) = derv2(iyf,kx) - d2trm*d2r2ikdx2(6)
      derv2(izf,kx) = derv2(izf,kx) - d2trm*d2r2ikdx2(5)
      derv2(ixf,ky) = derv2(ixf,ky) - d2trm*d2r2ikdx2(6)
      derv2(iyf,ky) = derv2(iyf,ky) - d2trm*d2r2ikdx2(2)
      derv2(izf,ky) = derv2(izf,ky) - d2trm*d2r2ikdx2(4)
      derv2(ixf,kz) = derv2(ixf,kz) - d2trm*d2r2ikdx2(5)
      derv2(iyf,kz) = derv2(iyf,kz) - d2trm*d2r2ikdx2(4)
      derv2(izf,kz) = derv2(izf,kz) - d2trm*d2r2ikdx2(3)
      derv2(ixf,kx) = derv2(ixf,kx) - d1trm
      derv2(iyf,ky) = derv2(iyf,ky) - d1trm
      derv2(izf,kz) = derv2(izf,kz) - d1trm
    endif
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
    d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
!
!  I-K
!
    if (lilocal) then
      derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xij
      derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xij
      derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xij
      derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*yij
      derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*yij
      derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*yij
      derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*zij
      derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*zij
      derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*zij
    endif
!
    if (lklocal) then
      derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xij
      derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*yij
      derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*zij
      derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xij
      derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*yij
      derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*zij
      derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xij
      derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*yij
      derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*zij
    endif
!
!  J-K
!
    if (ljlocal) then
      derv2(kxf,jx) = derv2(kxf,jx) + d3trm*xnbr(ni,i)*xij
      derv2(kyf,jx) = derv2(kyf,jx) + d3trm*ynbr(ni,i)*xij
      derv2(kzf,jx) = derv2(kzf,jx) + d3trm*znbr(ni,i)*xij
      derv2(kxf,jy) = derv2(kxf,jy) + d3trm*xnbr(ni,i)*yij
      derv2(kyf,jy) = derv2(kyf,jy) + d3trm*ynbr(ni,i)*yij
      derv2(kzf,jy) = derv2(kzf,jy) + d3trm*znbr(ni,i)*yij
      derv2(kxf,jz) = derv2(kxf,jz) + d3trm*xnbr(ni,i)*zij
      derv2(kyf,jz) = derv2(kyf,jz) + d3trm*ynbr(ni,i)*zij
      derv2(kzf,jz) = derv2(kzf,jz) + d3trm*znbr(ni,i)*zij
    endif
!
    if (lklocal) then
      derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(ni,i)*xij
      derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(ni,i)*yij
      derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(ni,i)*zij
      derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(ni,i)*xij
      derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(ni,i)*yij
      derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(ni,i)*zij
      derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(ni,i)*xij
      derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(ni,i)*yij
      derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(ni,i)*zij
    endif
!
!  I-J
!
    if (lilocal) then
      derv2(jxf,ix) = derv2(jxf,ix) - d3trm*xnbr(ni,i)*xij
      derv2(jyf,ix) = derv2(jyf,ix) - d3trm*xnbr(ni,i)*yij
      derv2(jzf,ix) = derv2(jzf,ix) - d3trm*xnbr(ni,i)*zij
      derv2(jxf,iy) = derv2(jxf,iy) - d3trm*ynbr(ni,i)*xij
      derv2(jyf,iy) = derv2(jyf,iy) - d3trm*ynbr(ni,i)*yij
      derv2(jzf,iy) = derv2(jzf,iy) - d3trm*ynbr(ni,i)*zij
      derv2(jxf,iz) = derv2(jxf,iz) - d3trm*znbr(ni,i)*xij
      derv2(jyf,iz) = derv2(jyf,iz) - d3trm*znbr(ni,i)*yij
      derv2(jzf,iz) = derv2(jzf,iz) - d3trm*znbr(ni,i)*zij
    endif
!
    if (ljlocal) then
      derv2(ixf,jx) = derv2(ixf,jx) - d3trm*xnbr(ni,i)*xij
      derv2(iyf,jx) = derv2(iyf,jx) - d3trm*ynbr(ni,i)*xij
      derv2(izf,jx) = derv2(izf,jx) - d3trm*znbr(ni,i)*xij
      derv2(ixf,jy) = derv2(ixf,jy) - d3trm*xnbr(ni,i)*yij
      derv2(iyf,jy) = derv2(iyf,jy) - d3trm*ynbr(ni,i)*yij
      derv2(izf,jy) = derv2(izf,jy) - d3trm*znbr(ni,i)*yij
      derv2(ixf,jz) = derv2(ixf,jz) - d3trm*xnbr(ni,i)*zij
      derv2(iyf,jz) = derv2(iyf,jz) - d3trm*ynbr(ni,i)*zij
      derv2(izf,jz) = derv2(izf,jz) - d3trm*znbr(ni,i)*zij
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      ni2 = nbrno_cn(n2,i)
      l = nbrno(ni2,i)
!
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
      if (.not.lijklocal.and..not.lllocal) cycle
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
      dcnidril = d1cndr_cn(n2,i)
      d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  I-K
!
      if (lilocal) then
        derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
      if (lklocal) then
        derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  I-L
!
      if (lilocal) then
        derv2(lxf,ix) = derv2(lxf,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,ix) = derv2(lyf,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,ix) = derv2(lzf,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lxf,iy) = derv2(lxf,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,iy) = derv2(lyf,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,iy) = derv2(lzf,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lxf,iz) = derv2(lxf,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,iz) = derv2(lyf,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,iz) = derv2(lzf,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
      if (lllocal) then
        derv2(ixf,lx) = derv2(ixf,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,lx) = derv2(iyf,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(izf,lx) = derv2(izf,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ixf,ly) = derv2(ixf,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(iyf,ly) = derv2(iyf,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(izf,ly) = derv2(izf,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(ixf,lz) = derv2(ixf,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(iyf,lz) = derv2(iyf,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(izf,lz) = derv2(izf,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  K-L
!
      if (lklocal) then
        derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
      if (lllocal) then
        derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
    enddo
!------------------------------------------------------------------------------------------
!  Loop over neighbours of j for second derivatives for 2 different coordination numbers  |
!------------------------------------------------------------------------------------------
    do n2 = 1,nnbr_cn(j)
      nj = nbrno_cn(n2,j)
      l = nbrno(nj,j)
!
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
      if (.not.lijklocal.and..not.lllocal) cycle
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
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = d2Edlogcnidlogcnj*dlogcnidcni*dlogcnjdcnj*dcnidrik*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  I-J
!
      if (lilocal) then
        derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
      if (ljlocal) then
        derv2(ixf,jx) = derv2(ixf,jx) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iyf,jx) = derv2(iyf,jx) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(izf,jx) = derv2(izf,jx) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ixf,jy) = derv2(ixf,jy) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iyf,jy) = derv2(iyf,jy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(izf,jy) = derv2(izf,jy) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ixf,jz) = derv2(ixf,jz) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iyf,jz) = derv2(iyf,jz) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(izf,jz) = derv2(izf,jz) + d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  I-L
!
      if (lilocal) then
        derv2(lxf,ix) = derv2(lxf,ix) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(lyf,ix) = derv2(lyf,ix) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(lzf,ix) = derv2(lzf,ix) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(lxf,iy) = derv2(lxf,iy) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(lyf,iy) = derv2(lyf,iy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(lzf,iy) = derv2(lzf,iy) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(lxf,iz) = derv2(lxf,iz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(lyf,iz) = derv2(lyf,iz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(lzf,iz) = derv2(lzf,iz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
      if (lllocal) then
        derv2(ixf,lx) = derv2(ixf,lx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iyf,lx) = derv2(iyf,lx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(izf,lx) = derv2(izf,lx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ixf,ly) = derv2(ixf,ly) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iyf,ly) = derv2(iyf,ly) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(izf,ly) = derv2(izf,ly) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ixf,lz) = derv2(ixf,lz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iyf,lz) = derv2(iyf,lz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(izf,lz) = derv2(izf,lz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  J-K
!
      if (ljlocal) then
        derv2(kxf,jx) = derv2(kxf,jx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(kyf,jx) = derv2(kyf,jx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(kzf,jx) = derv2(kzf,jx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(kxf,jy) = derv2(kxf,jy) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(kyf,jy) = derv2(kyf,jy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(kzf,jy) = derv2(kzf,jy) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(kxf,jz) = derv2(kxf,jz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(kyf,jz) = derv2(kyf,jz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(kzf,jz) = derv2(kzf,jz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
      if (lklocal) then
        derv2(jxf,kx) = derv2(jxf,kx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jyf,kx) = derv2(jyf,kx) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jzf,kx) = derv2(jzf,kx) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jxf,ky) = derv2(jxf,ky) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jyf,ky) = derv2(jyf,ky) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jzf,ky) = derv2(jzf,ky) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jxf,kz) = derv2(jxf,kz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jyf,kz) = derv2(jyf,kz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jzf,kz) = derv2(jzf,kz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  K-L
!     
      if (lklocal) then
        derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(ni,i)*znbr(nj,j)
      endif
!
      if (lllocal) then
        derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(ni,i)*znbr(nj,j)
      endif
    enddo
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    dcnjdrjk = d1cndr_cn(n,j)
    d1trm = dlogcnjdcnj*dcnjdrjk*dEdlogcnj
!
    d2r2jkdx2(1) = xnbr(nj,j)*xnbr(nj,j)
    d2r2jkdx2(2) = ynbr(nj,j)*ynbr(nj,j)
    d2r2jkdx2(3) = znbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(4) = ynbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(5) = xnbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(6) = xnbr(nj,j)*ynbr(nj,j)
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
    d2cnjdrjk2 = d2cndr_cn(n,j)
    d2trm = dEdlogcnj*dlogcnjdcnj*d2cnjdrjk2 + &
           (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjk
!----------------------------------------
!  Second derivatives of CN w.r.t. rjk  |
!----------------------------------------
    if (ljlocal) then
      derv2(kxf,jx) = derv2(kxf,jx) - d2trm*d2r2jkdx2(1)
      derv2(kyf,jx) = derv2(kyf,jx) - d2trm*d2r2jkdx2(6)
      derv2(kzf,jx) = derv2(kzf,jx) - d2trm*d2r2jkdx2(5)
      derv2(kxf,jy) = derv2(kxf,jy) - d2trm*d2r2jkdx2(6)
      derv2(kyf,jy) = derv2(kyf,jy) - d2trm*d2r2jkdx2(2)
      derv2(kzf,jy) = derv2(kzf,jy) - d2trm*d2r2jkdx2(4)
      derv2(kxf,jz) = derv2(kxf,jz) - d2trm*d2r2jkdx2(5)
      derv2(kyf,jz) = derv2(kyf,jz) - d2trm*d2r2jkdx2(4)
      derv2(kzf,jz) = derv2(kzf,jz) - d2trm*d2r2jkdx2(3)
      derv2(kxf,jx) = derv2(kxf,jx) - d1trm
      derv2(kyf,jy) = derv2(kyf,jy) - d1trm
      derv2(kzf,jz) = derv2(kzf,jz) - d1trm
    endif
!
    if (lklocal) then
      derv2(jxf,kx) = derv2(jxf,kx) - d2trm*d2r2jkdx2(1)
      derv2(jyf,kx) = derv2(jyf,kx) - d2trm*d2r2jkdx2(6)
      derv2(jzf,kx) = derv2(jzf,kx) - d2trm*d2r2jkdx2(5)
      derv2(jxf,ky) = derv2(jxf,ky) - d2trm*d2r2jkdx2(6)
      derv2(jyf,ky) = derv2(jyf,ky) - d2trm*d2r2jkdx2(2)
      derv2(jzf,ky) = derv2(jzf,ky) - d2trm*d2r2jkdx2(4)
      derv2(jxf,kz) = derv2(jxf,kz) - d2trm*d2r2jkdx2(5)
      derv2(jyf,kz) = derv2(jyf,kz) - d2trm*d2r2jkdx2(4)
      derv2(jzf,kz) = derv2(jzf,kz) - d2trm*d2r2jkdx2(3)
      derv2(jxf,kx) = derv2(jxf,kx) - d1trm
      derv2(jyf,ky) = derv2(jyf,ky) - d1trm
      derv2(jzf,kz) = derv2(jzf,kz) - d1trm
    endif
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
    d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
!
!  J-K
!
    if (ljlocal) then
      derv2(kxf,jx) = derv2(kxf,jx) + d3trm*xnbr(nj,j)*xij
      derv2(kyf,jx) = derv2(kyf,jx) + d3trm*ynbr(nj,j)*xij
      derv2(kzf,jx) = derv2(kzf,jx) + d3trm*znbr(nj,j)*xij
      derv2(kxf,jy) = derv2(kxf,jy) + d3trm*xnbr(nj,j)*yij
      derv2(kyf,jy) = derv2(kyf,jy) + d3trm*ynbr(nj,j)*yij
      derv2(kzf,jy) = derv2(kzf,jy) + d3trm*znbr(nj,j)*yij
      derv2(kxf,jz) = derv2(kxf,jz) + d3trm*xnbr(nj,j)*zij
      derv2(kyf,jz) = derv2(kyf,jz) + d3trm*ynbr(nj,j)*zij
      derv2(kzf,jz) = derv2(kzf,jz) + d3trm*znbr(nj,j)*zij
    endif
!
    if (lklocal) then
      derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(nj,j)*xij
      derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(nj,j)*yij
      derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(nj,j)*zij
      derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(nj,j)*xij
      derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(nj,j)*yij
      derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(nj,j)*zij
      derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(nj,j)*xij
      derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(nj,j)*yij
      derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(nj,j)*zij
    endif
!
!  I-K
!
    if (lilocal) then
      derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(nj,j)*xij
      derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(nj,j)*xij
      derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(nj,j)*xij
      derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(nj,j)*yij
      derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(nj,j)*yij
      derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(nj,j)*yij
      derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(nj,j)*zij
      derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(nj,j)*zij
      derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(nj,j)*zij
    endif
!
    if (lklocal) then
      derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(nj,j)*xij
      derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(nj,j)*yij
      derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(nj,j)*zij
      derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(nj,j)*xij
      derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(nj,j)*yij
      derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(nj,j)*zij
      derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(nj,j)*xij
      derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(nj,j)*yij
      derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(nj,j)*zij
    endif
!
!  J-I
!
    if (ljlocal) then
      derv2(ixf,jx) = derv2(ixf,jx) + d3trm*xnbr(nj,j)*xij
      derv2(iyf,jx) = derv2(iyf,jx) + d3trm*xnbr(nj,j)*yij
      derv2(izf,jx) = derv2(izf,jx) + d3trm*xnbr(nj,j)*zij
      derv2(ixf,jy) = derv2(ixf,jy) + d3trm*ynbr(nj,j)*xij
      derv2(iyf,jy) = derv2(iyf,jy) + d3trm*ynbr(nj,j)*yij
      derv2(izf,jy) = derv2(izf,jy) + d3trm*ynbr(nj,j)*zij
      derv2(ixf,jz) = derv2(ixf,jz) + d3trm*znbr(nj,j)*xij
      derv2(iyf,jz) = derv2(iyf,jz) + d3trm*znbr(nj,j)*yij
      derv2(izf,jz) = derv2(izf,jz) + d3trm*znbr(nj,j)*zij
    endif
!
    if (lilocal) then
      derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xij
      derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xij
      derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xij
      derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*yij
      derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*yij
      derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*yij
      derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*zij
      derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*zij
      derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*zij
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of j for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      nj2 = nbrno_cn(n2,j)
      l = nbrno(nj2,j)
!
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
      if (.not.lijklocal.and..not.lllocal) cycle
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
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  J-K
!
      if (ljlocal) then
        derv2(kxf,jx) = derv2(kxf,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(kyf,jx) = derv2(kyf,jx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kzf,jx) = derv2(kzf,jx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kxf,jy) = derv2(kxf,jy) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(kyf,jy) = derv2(kyf,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kzf,jy) = derv2(kzf,jy) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kxf,jz) = derv2(kxf,jz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(kyf,jz) = derv2(kyf,jz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kzf,jz) = derv2(kzf,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
      if (lklocal) then
        derv2(jxf,kx) = derv2(jxf,kx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,kx) = derv2(jyf,kx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,kx) = derv2(jzf,kx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jxf,ky) = derv2(jxf,ky) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,ky) = derv2(jyf,ky) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,ky) = derv2(jzf,ky) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jxf,kz) = derv2(jxf,kz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,kz) = derv2(jyf,kz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,kz) = derv2(jzf,kz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  J-L
!
      if (ljlocal) then
        derv2(lxf,jx) = derv2(lxf,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,jx) = derv2(lyf,jx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,jx) = derv2(lzf,jx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lxf,jy) = derv2(lxf,jy) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,jy) = derv2(lyf,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,jy) = derv2(lzf,jy) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lxf,jz) = derv2(lxf,jz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,jz) = derv2(lyf,jz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,jz) = derv2(lzf,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
      if (lllocal) then
        derv2(jxf,lx) = derv2(jxf,lx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,lx) = derv2(jyf,lx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jzf,lx) = derv2(jzf,lx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jxf,ly) = derv2(jxf,ly) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jyf,ly) = derv2(jyf,ly) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,ly) = derv2(jzf,ly) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jxf,lz) = derv2(jxf,lz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jyf,lz) = derv2(jyf,lz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jzf,lz) = derv2(jzf,lz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  K-L
!     
      if (lklocal) then
        derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
      if (lllocal) then
        derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnpd')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnpd_nox(lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,j,jx,jy,jz,jxf,jyf,jzf, &
                                  dEdlogcni,dEdlogcnj,d2Edlogcni2,d2Edlogcnj2,d2Edlogcnidlogcnj, &
                                  dlogcndcn,d2logcndcn2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Phonon version. Distributed memory parallel version. Version that excludes crossterms for distance/CN
!
!  On entry : 
!
!  ix,  iy,  iz      = second derivative Cartesian elements for i (local)
!  ixf, iyf, izf     = second derivative Cartesian elements for i (global)
!  jx,  jy,  jz      = second derivative Cartesian elements for j (local)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  dEdlogcnj         = derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcnj2       = second derivatives of energy w.r.t. log of the coordination number for j
!  d2Edlogcnidlogcnj = second derivatives of energy w.r.t. log of the coordination number for i and j
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   5/22 Created from gfnff_drv2_dcnpd
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: dEdlogcnj
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: d2Edlogcnidlogcnj
  real(dp),    intent(in)                          :: d2Edlogcnj2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
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
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  logical                                          :: lijklocal
  logical                                          :: lklocal
  logical                                          :: lllocal
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dcnjdrjl
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2logcnjdcnj2
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2jkdx2(6)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnpd_nox')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
  d2logcnidcni2 = d2logcndcn2(i)
  d2logcnjdcnj2 = d2logcndcn2(j)
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
!
    d2r2ikdx2(1) = xnbr(ni,i)*xnbr(ni,i)
    d2r2ikdx2(2) = ynbr(ni,i)*ynbr(ni,i)
    d2r2ikdx2(3) = znbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(4) = ynbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(5) = xnbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(6) = xnbr(ni,i)*ynbr(ni,i)
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
    d2cnidrik2 = d2cndr_cn(n,i)
    d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
           (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
    if (lilocal) then
      derv2(kxf,ix) = derv2(kxf,ix) - d2trm*d2r2ikdx2(1)
      derv2(kyf,ix) = derv2(kyf,ix) - d2trm*d2r2ikdx2(6)
      derv2(kzf,ix) = derv2(kzf,ix) - d2trm*d2r2ikdx2(5)
      derv2(kxf,iy) = derv2(kxf,iy) - d2trm*d2r2ikdx2(6)
      derv2(kyf,iy) = derv2(kyf,iy) - d2trm*d2r2ikdx2(2)
      derv2(kzf,iy) = derv2(kzf,iy) - d2trm*d2r2ikdx2(4)
      derv2(kxf,iz) = derv2(kxf,iz) - d2trm*d2r2ikdx2(5)
      derv2(kyf,iz) = derv2(kyf,iz) - d2trm*d2r2ikdx2(4)
      derv2(kzf,iz) = derv2(kzf,iz) - d2trm*d2r2ikdx2(3)
      derv2(kxf,ix) = derv2(kxf,ix) - d1trm
      derv2(kyf,iy) = derv2(kyf,iy) - d1trm
      derv2(kzf,iz) = derv2(kzf,iz) - d1trm
    endif
!
    if (lklocal) then
      derv2(ixf,kx) = derv2(ixf,kx) - d2trm*d2r2ikdx2(1)
      derv2(iyf,kx) = derv2(iyf,kx) - d2trm*d2r2ikdx2(6)
      derv2(izf,kx) = derv2(izf,kx) - d2trm*d2r2ikdx2(5)
      derv2(ixf,ky) = derv2(ixf,ky) - d2trm*d2r2ikdx2(6)
      derv2(iyf,ky) = derv2(iyf,ky) - d2trm*d2r2ikdx2(2)
      derv2(izf,ky) = derv2(izf,ky) - d2trm*d2r2ikdx2(4)
      derv2(ixf,kz) = derv2(ixf,kz) - d2trm*d2r2ikdx2(5)
      derv2(iyf,kz) = derv2(iyf,kz) - d2trm*d2r2ikdx2(4)
      derv2(izf,kz) = derv2(izf,kz) - d2trm*d2r2ikdx2(3)
      derv2(ixf,kx) = derv2(ixf,kx) - d1trm
      derv2(iyf,ky) = derv2(iyf,ky) - d1trm
      derv2(izf,kz) = derv2(izf,kz) - d1trm
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      ni2 = nbrno_cn(n2,i)
      l = nbrno(ni2,i)
!
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
      if (.not.lijklocal.and..not.lllocal) cycle
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
      dcnidril = d1cndr_cn(n2,i)
      d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  I-K
!
      if (lilocal) then
        derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
      if (lklocal) then
        derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  I-L
!
      if (lilocal) then
        derv2(lxf,ix) = derv2(lxf,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,ix) = derv2(lyf,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,ix) = derv2(lzf,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lxf,iy) = derv2(lxf,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,iy) = derv2(lyf,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,iy) = derv2(lzf,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lxf,iz) = derv2(lxf,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,iz) = derv2(lyf,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,iz) = derv2(lzf,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
      if (lllocal) then
        derv2(ixf,lx) = derv2(ixf,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iyf,lx) = derv2(iyf,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(izf,lx) = derv2(izf,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ixf,ly) = derv2(ixf,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(iyf,ly) = derv2(iyf,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(izf,ly) = derv2(izf,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(ixf,lz) = derv2(ixf,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(iyf,lz) = derv2(iyf,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(izf,lz) = derv2(izf,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  K-L
!
      if (lklocal) then
        derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
      if (lllocal) then
        derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
    enddo
!------------------------------------------------------------------------------------------
!  Loop over neighbours of j for second derivatives for 2 different coordination numbers  |
!------------------------------------------------------------------------------------------
    do n2 = 1,nnbr_cn(j)
      nj = nbrno_cn(n2,j)
      l = nbrno(nj,j)
!
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
      if (.not.lijklocal.and..not.lllocal) cycle
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
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = d2Edlogcnidlogcnj*dlogcnidcni*dlogcnjdcnj*dcnidrik*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  I-J
!
      if (lilocal) then
        derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
      if (ljlocal) then
        derv2(ixf,jx) = derv2(ixf,jx) + d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iyf,jx) = derv2(iyf,jx) + d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(izf,jx) = derv2(izf,jx) + d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ixf,jy) = derv2(ixf,jy) + d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iyf,jy) = derv2(iyf,jy) + d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(izf,jy) = derv2(izf,jy) + d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ixf,jz) = derv2(ixf,jz) + d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iyf,jz) = derv2(iyf,jz) + d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(izf,jz) = derv2(izf,jz) + d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  I-L
!
      if (lilocal) then
        derv2(lxf,ix) = derv2(lxf,ix) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(lyf,ix) = derv2(lyf,ix) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(lzf,ix) = derv2(lzf,ix) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(lxf,iy) = derv2(lxf,iy) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(lyf,iy) = derv2(lyf,iy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(lzf,iy) = derv2(lzf,iy) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(lxf,iz) = derv2(lxf,iz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(lyf,iz) = derv2(lyf,iz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(lzf,iz) = derv2(lzf,iz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
      if (lllocal) then
        derv2(ixf,lx) = derv2(ixf,lx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(iyf,lx) = derv2(iyf,lx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(izf,lx) = derv2(izf,lx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(ixf,ly) = derv2(ixf,ly) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(iyf,ly) = derv2(iyf,ly) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(izf,ly) = derv2(izf,ly) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(ixf,lz) = derv2(ixf,lz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(iyf,lz) = derv2(iyf,lz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(izf,lz) = derv2(izf,lz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  J-K
!
      if (ljlocal) then
        derv2(kxf,jx) = derv2(kxf,jx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(kyf,jx) = derv2(kyf,jx) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(kzf,jx) = derv2(kzf,jx) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(kxf,jy) = derv2(kxf,jy) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(kyf,jy) = derv2(kyf,jy) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(kzf,jy) = derv2(kzf,jy) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(kxf,jz) = derv2(kxf,jz) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(kyf,jz) = derv2(kyf,jz) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(kzf,jz) = derv2(kzf,jz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
      if (lklocal) then
        derv2(jxf,kx) = derv2(jxf,kx) - d3trm*xnbr(nj,j)*xnbr(ni,i)
        derv2(jyf,kx) = derv2(jyf,kx) - d3trm*ynbr(nj,j)*xnbr(ni,i)
        derv2(jzf,kx) = derv2(jzf,kx) - d3trm*znbr(nj,j)*xnbr(ni,i)
        derv2(jxf,ky) = derv2(jxf,ky) - d3trm*xnbr(nj,j)*ynbr(ni,i)
        derv2(jyf,ky) = derv2(jyf,ky) - d3trm*ynbr(nj,j)*ynbr(ni,i)
        derv2(jzf,ky) = derv2(jzf,ky) - d3trm*znbr(nj,j)*ynbr(ni,i)
        derv2(jxf,kz) = derv2(jxf,kz) - d3trm*xnbr(nj,j)*znbr(ni,i)
        derv2(jyf,kz) = derv2(jyf,kz) - d3trm*ynbr(nj,j)*znbr(ni,i)
        derv2(jzf,kz) = derv2(jzf,kz) - d3trm*znbr(nj,j)*znbr(ni,i)
      endif
!
!  K-L
!     
      if (lklocal) then
        derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(ni,i)*znbr(nj,j)
      endif
!
      if (lllocal) then
        derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(ni,i)*xnbr(nj,j)
        derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(ni,i)*xnbr(nj,j)
        derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(ni,i)*xnbr(nj,j)
        derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(ni,i)*ynbr(nj,j)
        derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(ni,i)*ynbr(nj,j)
        derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(ni,i)*ynbr(nj,j)
        derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(ni,i)*znbr(nj,j)
        derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(ni,i)*znbr(nj,j)
        derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(ni,i)*znbr(nj,j)
      endif
    enddo
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    dcnjdrjk = d1cndr_cn(n,j)
    d1trm = dlogcnjdcnj*dcnjdrjk*dEdlogcnj
!
    d2r2jkdx2(1) = xnbr(nj,j)*xnbr(nj,j)
    d2r2jkdx2(2) = ynbr(nj,j)*ynbr(nj,j)
    d2r2jkdx2(3) = znbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(4) = ynbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(5) = xnbr(nj,j)*znbr(nj,j)
    d2r2jkdx2(6) = xnbr(nj,j)*ynbr(nj,j)
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
    d2cnjdrjk2 = d2cndr_cn(n,j)
    d2trm = dEdlogcnj*dlogcnjdcnj*d2cnjdrjk2 + &
           (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjk
!----------------------------------------
!  Second derivatives of CN w.r.t. rjk  |
!----------------------------------------
    if (ljlocal) then
      derv2(kxf,jx) = derv2(kxf,jx) - d2trm*d2r2jkdx2(1)
      derv2(kyf,jx) = derv2(kyf,jx) - d2trm*d2r2jkdx2(6)
      derv2(kzf,jx) = derv2(kzf,jx) - d2trm*d2r2jkdx2(5)
      derv2(kxf,jy) = derv2(kxf,jy) - d2trm*d2r2jkdx2(6)
      derv2(kyf,jy) = derv2(kyf,jy) - d2trm*d2r2jkdx2(2)
      derv2(kzf,jy) = derv2(kzf,jy) - d2trm*d2r2jkdx2(4)
      derv2(kxf,jz) = derv2(kxf,jz) - d2trm*d2r2jkdx2(5)
      derv2(kyf,jz) = derv2(kyf,jz) - d2trm*d2r2jkdx2(4)
      derv2(kzf,jz) = derv2(kzf,jz) - d2trm*d2r2jkdx2(3)
      derv2(kxf,jx) = derv2(kxf,jx) - d1trm
      derv2(kyf,jy) = derv2(kyf,jy) - d1trm
      derv2(kzf,jz) = derv2(kzf,jz) - d1trm
    endif
!
    if (lklocal) then
      derv2(jxf,kx) = derv2(jxf,kx) - d2trm*d2r2jkdx2(1)
      derv2(jyf,kx) = derv2(jyf,kx) - d2trm*d2r2jkdx2(6)
      derv2(jzf,kx) = derv2(jzf,kx) - d2trm*d2r2jkdx2(5)
      derv2(jxf,ky) = derv2(jxf,ky) - d2trm*d2r2jkdx2(6)
      derv2(jyf,ky) = derv2(jyf,ky) - d2trm*d2r2jkdx2(2)
      derv2(jzf,ky) = derv2(jzf,ky) - d2trm*d2r2jkdx2(4)
      derv2(jxf,kz) = derv2(jxf,kz) - d2trm*d2r2jkdx2(5)
      derv2(jyf,kz) = derv2(jyf,kz) - d2trm*d2r2jkdx2(4)
      derv2(jzf,kz) = derv2(jzf,kz) - d2trm*d2r2jkdx2(3)
      derv2(jxf,kx) = derv2(jxf,kx) - d1trm
      derv2(jyf,ky) = derv2(jyf,ky) - d1trm
      derv2(jzf,kz) = derv2(jzf,kz) - d1trm
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of j for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      nj2 = nbrno_cn(n2,j)
      l = nbrno(nj2,j)
!
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
      if (.not.lijklocal.and..not.lllocal) cycle
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
      dcnjdrjl = d1cndr_cn(n2,j)
      d3trm = (d2Edlogcnj2*dlogcnjdcnj*dlogcnjdcnj + dEdlogcnj*d2logcnjdcnj2)*dcnjdrjk*dcnjdrjl
!
      if (abs(d3trm).lt.gfnff_cnc6tol) cycle
!
!  J-K
!
      if (ljlocal) then
        derv2(kxf,jx) = derv2(kxf,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(kyf,jx) = derv2(kyf,jx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kzf,jx) = derv2(kzf,jx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kxf,jy) = derv2(kxf,jy) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(kyf,jy) = derv2(kyf,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kzf,jy) = derv2(kzf,jy) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kxf,jz) = derv2(kxf,jz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(kyf,jz) = derv2(kyf,jz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kzf,jz) = derv2(kzf,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
      if (lklocal) then
        derv2(jxf,kx) = derv2(jxf,kx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,kx) = derv2(jyf,kx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,kx) = derv2(jzf,kx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jxf,ky) = derv2(jxf,ky) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,ky) = derv2(jyf,ky) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,ky) = derv2(jzf,ky) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jxf,kz) = derv2(jxf,kz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,kz) = derv2(jyf,kz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,kz) = derv2(jzf,kz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  J-L
!
      if (ljlocal) then
        derv2(lxf,jx) = derv2(lxf,jx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,jx) = derv2(lyf,jx) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,jx) = derv2(lzf,jx) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lxf,jy) = derv2(lxf,jy) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,jy) = derv2(lyf,jy) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,jy) = derv2(lzf,jy) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lxf,jz) = derv2(lxf,jz) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,jz) = derv2(lyf,jz) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,jz) = derv2(lzf,jz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
      if (lllocal) then
        derv2(jxf,lx) = derv2(jxf,lx) - d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(jyf,lx) = derv2(jyf,lx) - d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(jzf,lx) = derv2(jzf,lx) - d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(jxf,ly) = derv2(jxf,ly) - d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(jyf,ly) = derv2(jyf,ly) - d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(jzf,ly) = derv2(jzf,ly) - d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(jxf,lz) = derv2(jxf,lz) - d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(jyf,lz) = derv2(jyf,lz) - d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(jzf,lz) = derv2(jzf,lz) - d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
!  K-L
!     
      if (lklocal) then
        derv2(lxf,kx) = derv2(lxf,kx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,kx) = derv2(lyf,kx) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,kx) = derv2(lzf,kx) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(lxf,ky) = derv2(lxf,ky) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,ky) = derv2(lyf,ky) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,ky) = derv2(lzf,ky) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(lxf,kz) = derv2(lxf,kz) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(lyf,kz) = derv2(lyf,kz) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(lzf,kz) = derv2(lzf,kz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
!
      if (lllocal) then
        derv2(kxf,lx) = derv2(kxf,lx) + d3trm*xnbr(nj,j)*xnbr(nj2,j)
        derv2(kyf,lx) = derv2(kyf,lx) + d3trm*ynbr(nj,j)*xnbr(nj2,j)
        derv2(kzf,lx) = derv2(kzf,lx) + d3trm*znbr(nj,j)*xnbr(nj2,j)
        derv2(kxf,ly) = derv2(kxf,ly) + d3trm*xnbr(nj,j)*ynbr(nj2,j)
        derv2(kyf,ly) = derv2(kyf,ly) + d3trm*ynbr(nj,j)*ynbr(nj2,j)
        derv2(kzf,ly) = derv2(kzf,ly) + d3trm*znbr(nj,j)*ynbr(nj2,j)
        derv2(kxf,lz) = derv2(kxf,lz) + d3trm*xnbr(nj,j)*znbr(nj2,j)
        derv2(kyf,lz) = derv2(kyf,lz) + d3trm*ynbr(nj,j)*znbr(nj2,j)
        derv2(kzf,lz) = derv2(kzf,lz) + d3trm*znbr(nj,j)*znbr(nj2,j)
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnpd_nox')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnpd_x(lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,j,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                                d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Phonon version. Distributed memory parallel version. Version that only does crossterms for distance/CN
!
!  On entry : 
!
!  ix,  iy,  iz      = second derivative Cartesian elements for i (local)
!  ixf, iyf, izf     = second derivative Cartesian elements for i (global)
!  jx,  jy,  jz      = second derivative Cartesian elements for j (local)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   5/22 Created from gfnff_drv2_dcnpd
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  logical                                          :: lijklocal
  logical                                          :: lklocal
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnpd_x')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    dcnidrik = d1cndr_cn(n,i)
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
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
    d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
!
!  I-K
!
    if (lilocal) then
      derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xij
      derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xij
      derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xij
      derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*yij
      derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*yij
      derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*yij
      derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*zij
      derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*zij
      derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*zij
    endif
!
    if (lklocal) then
      derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xij
      derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*yij
      derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*zij
      derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xij
      derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*yij
      derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*zij
      derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xij
      derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*yij
      derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*zij
    endif
!
!  J-K
!
    if (ljlocal) then
      derv2(kxf,jx) = derv2(kxf,jx) + d3trm*xnbr(ni,i)*xij
      derv2(kyf,jx) = derv2(kyf,jx) + d3trm*ynbr(ni,i)*xij
      derv2(kzf,jx) = derv2(kzf,jx) + d3trm*znbr(ni,i)*xij
      derv2(kxf,jy) = derv2(kxf,jy) + d3trm*xnbr(ni,i)*yij
      derv2(kyf,jy) = derv2(kyf,jy) + d3trm*ynbr(ni,i)*yij
      derv2(kzf,jy) = derv2(kzf,jy) + d3trm*znbr(ni,i)*yij
      derv2(kxf,jz) = derv2(kxf,jz) + d3trm*xnbr(ni,i)*zij
      derv2(kyf,jz) = derv2(kyf,jz) + d3trm*ynbr(ni,i)*zij
      derv2(kzf,jz) = derv2(kzf,jz) + d3trm*znbr(ni,i)*zij
    endif
!
    if (lklocal) then
      derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(ni,i)*xij
      derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(ni,i)*yij
      derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(ni,i)*zij
      derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(ni,i)*xij
      derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(ni,i)*yij
      derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(ni,i)*zij
      derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(ni,i)*xij
      derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(ni,i)*yij
      derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(ni,i)*zij
    endif
!
!  I-J
!
    if (lilocal) then
      derv2(jxf,ix) = derv2(jxf,ix) - d3trm*xnbr(ni,i)*xij
      derv2(jyf,ix) = derv2(jyf,ix) - d3trm*xnbr(ni,i)*yij
      derv2(jzf,ix) = derv2(jzf,ix) - d3trm*xnbr(ni,i)*zij
      derv2(jxf,iy) = derv2(jxf,iy) - d3trm*ynbr(ni,i)*xij
      derv2(jyf,iy) = derv2(jyf,iy) - d3trm*ynbr(ni,i)*yij
      derv2(jzf,iy) = derv2(jzf,iy) - d3trm*ynbr(ni,i)*zij
      derv2(jxf,iz) = derv2(jxf,iz) - d3trm*znbr(ni,i)*xij
      derv2(jyf,iz) = derv2(jyf,iz) - d3trm*znbr(ni,i)*yij
      derv2(jzf,iz) = derv2(jzf,iz) - d3trm*znbr(ni,i)*zij
    endif
!
    if (ljlocal) then
      derv2(ixf,jx) = derv2(ixf,jx) - d3trm*xnbr(ni,i)*xij
      derv2(iyf,jx) = derv2(iyf,jx) - d3trm*ynbr(ni,i)*xij
      derv2(izf,jx) = derv2(izf,jx) - d3trm*znbr(ni,i)*xij
      derv2(ixf,jy) = derv2(ixf,jy) - d3trm*xnbr(ni,i)*yij
      derv2(iyf,jy) = derv2(iyf,jy) - d3trm*ynbr(ni,i)*yij
      derv2(izf,jy) = derv2(izf,jy) - d3trm*znbr(ni,i)*yij
      derv2(ixf,jz) = derv2(ixf,jz) - d3trm*xnbr(ni,i)*zij
      derv2(iyf,jz) = derv2(iyf,jz) - d3trm*ynbr(ni,i)*zij
      derv2(izf,jz) = derv2(izf,jz) - d3trm*znbr(ni,i)*zij
    endif
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
    lijklocal = (lilocal.or.ljlocal.or.lklocal)
!
    dcnjdrjk = d1cndr_cn(n,j)
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
    d2cnjdrjk2 = d2cndr_cn(n,j)
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
    d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
!
!  J-K
!
    if (ljlocal) then
      derv2(kxf,jx) = derv2(kxf,jx) + d3trm*xnbr(nj,j)*xij
      derv2(kyf,jx) = derv2(kyf,jx) + d3trm*ynbr(nj,j)*xij
      derv2(kzf,jx) = derv2(kzf,jx) + d3trm*znbr(nj,j)*xij
      derv2(kxf,jy) = derv2(kxf,jy) + d3trm*xnbr(nj,j)*yij
      derv2(kyf,jy) = derv2(kyf,jy) + d3trm*ynbr(nj,j)*yij
      derv2(kzf,jy) = derv2(kzf,jy) + d3trm*znbr(nj,j)*yij
      derv2(kxf,jz) = derv2(kxf,jz) + d3trm*xnbr(nj,j)*zij
      derv2(kyf,jz) = derv2(kyf,jz) + d3trm*ynbr(nj,j)*zij
      derv2(kzf,jz) = derv2(kzf,jz) + d3trm*znbr(nj,j)*zij
    endif
!
    if (lklocal) then
      derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(nj,j)*xij
      derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(nj,j)*yij
      derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(nj,j)*zij
      derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(nj,j)*xij
      derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(nj,j)*yij
      derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(nj,j)*zij
      derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(nj,j)*xij
      derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(nj,j)*yij
      derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(nj,j)*zij
    endif
!
!  I-K
!
    if (lilocal) then
      derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(nj,j)*xij
      derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(nj,j)*xij
      derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(nj,j)*xij
      derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(nj,j)*yij
      derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(nj,j)*yij
      derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(nj,j)*yij
      derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(nj,j)*zij
      derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(nj,j)*zij
      derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(nj,j)*zij
    endif
!
    if (lklocal) then
      derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(nj,j)*xij
      derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(nj,j)*yij
      derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(nj,j)*zij
      derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(nj,j)*xij
      derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(nj,j)*yij
      derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(nj,j)*zij
      derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(nj,j)*xij
      derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(nj,j)*yij
      derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(nj,j)*zij
    endif
!
!  J-I
!
    if (ljlocal) then
      derv2(ixf,jx) = derv2(ixf,jx) + d3trm*xnbr(nj,j)*xij
      derv2(iyf,jx) = derv2(iyf,jx) + d3trm*xnbr(nj,j)*yij
      derv2(izf,jx) = derv2(izf,jx) + d3trm*xnbr(nj,j)*zij
      derv2(ixf,jy) = derv2(ixf,jy) + d3trm*ynbr(nj,j)*xij
      derv2(iyf,jy) = derv2(iyf,jy) + d3trm*ynbr(nj,j)*yij
      derv2(izf,jy) = derv2(izf,jy) + d3trm*ynbr(nj,j)*zij
      derv2(ixf,jz) = derv2(ixf,jz) + d3trm*znbr(nj,j)*xij
      derv2(iyf,jz) = derv2(iyf,jz) + d3trm*znbr(nj,j)*yij
      derv2(izf,jz) = derv2(izf,jz) + d3trm*znbr(nj,j)*zij
    endif
!
    if (lilocal) then
      derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xij
      derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xij
      derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xij
      derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*yij
      derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*yij
      derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*yij
      derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*zij
      derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*zij
      derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*zij
    endif
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnpd_x')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnpd_x1(i,ix,iy,iz,j,jxf,jyf,jzf,xij,yij,zij, &
                                 d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Phonon version. Distributed memory parallel version. Version that only does crossterms for distance/CN
!  Version where i is local.
!
!  On entry : 
!
!  ix,  iy,  iz      = second derivative Cartesian elements for i (local)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   5/22 Created from gfnff_drv2_dcnpd_x
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: k
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnpd_x1')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
!
    dcnidrik = d1cndr_cn(n,i)
!
    indk = 3*(k-1)
    kxf = indk + 1
    kyf = indk + 2
    kzf = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
    d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
!
!  I-K
!
    derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(ni,i)*xij
    derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(ni,i)*xij
    derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(ni,i)*xij
    derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(ni,i)*yij
    derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(ni,i)*yij
    derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(ni,i)*yij
    derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(ni,i)*zij
    derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(ni,i)*zij
    derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(ni,i)*zij
!
!  I-J
!
    derv2(jxf,ix) = derv2(jxf,ix) - d3trm*xnbr(ni,i)*xij
    derv2(jyf,ix) = derv2(jyf,ix) - d3trm*xnbr(ni,i)*yij
    derv2(jzf,ix) = derv2(jzf,ix) - d3trm*xnbr(ni,i)*zij
    derv2(jxf,iy) = derv2(jxf,iy) - d3trm*ynbr(ni,i)*xij
    derv2(jyf,iy) = derv2(jyf,iy) - d3trm*ynbr(ni,i)*yij
    derv2(jzf,iy) = derv2(jzf,iy) - d3trm*ynbr(ni,i)*zij
    derv2(jxf,iz) = derv2(jxf,iz) - d3trm*znbr(ni,i)*xij
    derv2(jyf,iz) = derv2(jyf,iz) - d3trm*znbr(ni,i)*yij
    derv2(jzf,iz) = derv2(jzf,iz) - d3trm*znbr(ni,i)*zij
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
!
    dcnjdrjk = d1cndr_cn(n,j)
!
    indk = 3*(k-1)
    kxf = indk + 1
    kyf = indk + 2
    kzf = indk + 3
!
    d2cnjdrjk2 = d2cndr_cn(n,j)
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
    d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
!
!  I-K
!
    derv2(kxf,ix) = derv2(kxf,ix) - d3trm*xnbr(nj,j)*xij
    derv2(kyf,ix) = derv2(kyf,ix) - d3trm*ynbr(nj,j)*xij
    derv2(kzf,ix) = derv2(kzf,ix) - d3trm*znbr(nj,j)*xij
    derv2(kxf,iy) = derv2(kxf,iy) - d3trm*xnbr(nj,j)*yij
    derv2(kyf,iy) = derv2(kyf,iy) - d3trm*ynbr(nj,j)*yij
    derv2(kzf,iy) = derv2(kzf,iy) - d3trm*znbr(nj,j)*yij
    derv2(kxf,iz) = derv2(kxf,iz) - d3trm*xnbr(nj,j)*zij
    derv2(kyf,iz) = derv2(kyf,iz) - d3trm*ynbr(nj,j)*zij
    derv2(kzf,iz) = derv2(kzf,iz) - d3trm*znbr(nj,j)*zij
!
!  J-I
!
    derv2(jxf,ix) = derv2(jxf,ix) + d3trm*xnbr(nj,j)*xij
    derv2(jyf,ix) = derv2(jyf,ix) + d3trm*ynbr(nj,j)*xij
    derv2(jzf,ix) = derv2(jzf,ix) + d3trm*znbr(nj,j)*xij
    derv2(jxf,iy) = derv2(jxf,iy) + d3trm*xnbr(nj,j)*yij
    derv2(jyf,iy) = derv2(jyf,iy) + d3trm*ynbr(nj,j)*yij
    derv2(jzf,iy) = derv2(jzf,iy) + d3trm*znbr(nj,j)*yij
    derv2(jxf,iz) = derv2(jxf,iz) + d3trm*xnbr(nj,j)*zij
    derv2(jyf,iz) = derv2(jyf,iz) + d3trm*ynbr(nj,j)*zij
    derv2(jzf,iz) = derv2(jzf,iz) + d3trm*znbr(nj,j)*zij
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnpd_x1')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnpd_x2(i,ixf,iyf,izf,j,jxf,jyf,jzf,xij,yij,zij, &
                                 d2Edlogcnidrij,d2Edlogcnjdrij,dlogcndcn)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Phonon version. Distributed memory parallel version. Version that only does crossterms for distance/CN
!  Version where i and j are global and k is local.
!
!  On entry : 
!
!  ixf, iyf, izf     = second derivative Cartesian elements for i (global)
!  jxf, jyf, jzf     = second derivative Cartesian elements for j (global)
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d2Edlogcnidrij    = second derivatives of energy w.r.t. log of the coordination number for i and rij
!  d2Edlogcnjdrij    = second derivatives of energy w.r.t. log of the coordination number for j and rij
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   5/22 Created from gfnff_drv2_dcnpd_x
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
!  Julian Gale, CIC, Curtin University, May 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  real(dp),    intent(in)                          :: d2Edlogcnidrij
  real(dp),    intent(in)                          :: d2Edlogcnjdrij
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  logical                                          :: lklocal
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: d2cnjdrjk2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnpd_x2')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  dlogcnjdcnj = dlogcndcn(j)
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
!
    if (.not.lklocal) cycle
!
    dcnidrik = d1cndr_cn(n,i)
!
    indk = 3*(kloc-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!-------------------------------------------
!  Second derivatives between rij and rik  |
!-------------------------------------------
    d3trm = d2Edlogcnidrij*dlogcnidcni*dcnidrik
!
!  I-K
!
    derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(ni,i)*xij
    derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(ni,i)*yij
    derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(ni,i)*zij
    derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(ni,i)*xij
    derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(ni,i)*yij
    derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(ni,i)*zij
    derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(ni,i)*xij
    derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(ni,i)*yij
    derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(ni,i)*zij
!
!  J-K
!
    derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(ni,i)*xij
    derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(ni,i)*yij
    derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(ni,i)*zij
    derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(ni,i)*xij
    derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(ni,i)*yij
    derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(ni,i)*zij
    derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(ni,i)*xij
    derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(ni,i)*yij
    derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(ni,i)*zij
  enddo
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for j  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(j)
    nj = nbrno_cn(n,j)
    k = nbrno(nj,j)
!
    kloc = atom2local(k)
    lklocal = (kloc.ne.0)
!
    if (.not.lklocal) cycle
!
    dcnjdrjk = d1cndr_cn(n,j)
!
    indk = 3*(kloc-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!
    d2cnjdrjk2 = d2cndr_cn(n,j)
!-------------------------------------------
!  Second derivatives between rij and rjk  |
!-------------------------------------------
    d3trm = d2Edlogcnjdrij*dlogcnjdcnj*dcnjdrjk
!
!  J-K
!
    derv2(jxf,kx) = derv2(jxf,kx) + d3trm*xnbr(nj,j)*xij
    derv2(jyf,kx) = derv2(jyf,kx) + d3trm*xnbr(nj,j)*yij
    derv2(jzf,kx) = derv2(jzf,kx) + d3trm*xnbr(nj,j)*zij
    derv2(jxf,ky) = derv2(jxf,ky) + d3trm*ynbr(nj,j)*xij
    derv2(jyf,ky) = derv2(jyf,ky) + d3trm*ynbr(nj,j)*yij
    derv2(jzf,ky) = derv2(jzf,ky) + d3trm*ynbr(nj,j)*zij
    derv2(jxf,kz) = derv2(jxf,kz) + d3trm*znbr(nj,j)*xij
    derv2(jyf,kz) = derv2(jyf,kz) + d3trm*znbr(nj,j)*yij
    derv2(jzf,kz) = derv2(jzf,kz) + d3trm*znbr(nj,j)*zij
!
!  I-K
!
    derv2(ixf,kx) = derv2(ixf,kx) - d3trm*xnbr(nj,j)*xij
    derv2(iyf,kx) = derv2(iyf,kx) - d3trm*xnbr(nj,j)*yij
    derv2(izf,kx) = derv2(izf,kx) - d3trm*xnbr(nj,j)*zij
    derv2(ixf,ky) = derv2(ixf,ky) - d3trm*ynbr(nj,j)*xij
    derv2(iyf,ky) = derv2(iyf,ky) - d3trm*ynbr(nj,j)*yij
    derv2(izf,ky) = derv2(izf,ky) - d3trm*ynbr(nj,j)*zij
    derv2(ixf,kz) = derv2(ixf,kz) - d3trm*znbr(nj,j)*xij
    derv2(iyf,kz) = derv2(iyf,kz) - d3trm*znbr(nj,j)*yij
    derv2(izf,kz) = derv2(izf,kz) - d3trm*znbr(nj,j)*zij
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnpd_x2')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_selfp(i,dEdlogcni,d2Edlogcni2,dlogcndcn,d2logcndcn2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Special case version of gfnff_drv2_dcn for i only self terms. Phonon version.
!
!  On entry : 
!
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   3/21 Created from gfnff_drv2_dcn_self
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
!  Julian Gale, CIC, Curtin University, March 2021
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: k
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2r2ikdx2(6)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_selfp')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d2logcnidcni2 = d2logcndcn2(i)
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
!
    d2r2ikdx2(1) = xnbr(ni,i)*xnbr(ni,i)
    d2r2ikdx2(2) = ynbr(ni,i)*ynbr(ni,i)
    d2r2ikdx2(3) = znbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(4) = ynbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(5) = xnbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(6) = xnbr(ni,i)*ynbr(ni,i)
!
    indk = 3*(k-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!
    d2cnidrik2 = d2cndr_cn(n,i)
    d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
           (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
    if (k.ge.i) then
      derv2(kx,ix) = derv2(kx,ix) - d2trm*d2r2ikdx2(1)
      derv2(ky,ix) = derv2(ky,ix) - d2trm*d2r2ikdx2(6)
      derv2(kz,ix) = derv2(kz,ix) - d2trm*d2r2ikdx2(5)
      derv2(kx,iy) = derv2(kx,iy) - d2trm*d2r2ikdx2(6)
      derv2(ky,iy) = derv2(ky,iy) - d2trm*d2r2ikdx2(2)
      derv2(kz,iy) = derv2(kz,iy) - d2trm*d2r2ikdx2(4)
      derv2(kx,iz) = derv2(kx,iz) - d2trm*d2r2ikdx2(5)
      derv2(ky,iz) = derv2(ky,iz) - d2trm*d2r2ikdx2(4)
      derv2(kz,iz) = derv2(kz,iz) - d2trm*d2r2ikdx2(3)
      derv2(kx,ix) = derv2(kx,ix) - d1trm
      derv2(ky,iy) = derv2(ky,iy) - d1trm
      derv2(kz,iz) = derv2(kz,iz) - d1trm
    else
      derv2(ix,kx) = derv2(ix,kx) - d2trm*d2r2ikdx2(1)
      derv2(iy,kx) = derv2(iy,kx) - d2trm*d2r2ikdx2(6)
      derv2(iz,kx) = derv2(iz,kx) - d2trm*d2r2ikdx2(5)
      derv2(ix,ky) = derv2(ix,ky) - d2trm*d2r2ikdx2(6)
      derv2(iy,ky) = derv2(iy,ky) - d2trm*d2r2ikdx2(2)
      derv2(iz,ky) = derv2(iz,ky) - d2trm*d2r2ikdx2(4)
      derv2(ix,kz) = derv2(ix,kz) - d2trm*d2r2ikdx2(5)
      derv2(iy,kz) = derv2(iy,kz) - d2trm*d2r2ikdx2(4)
      derv2(iz,kz) = derv2(iz,kz) - d2trm*d2r2ikdx2(3)
      derv2(ix,kx) = derv2(ix,kx) - d1trm
      derv2(iy,ky) = derv2(iy,ky) - d1trm
      derv2(iz,kz) = derv2(iz,kz) - d1trm
    endif
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      ni2 = nbrno_cn(n2,i)
      l = nbrno(ni2,i)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnidril = d1cndr_cn(n2,i)
      d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
!  I-K
!
      if (k.ge.i) then
        derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      else
        derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  I-L
!
      if (l.ge.i) then
        derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ly,ix) = derv2(ly,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lz,ix) = derv2(lz,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lx,iy) = derv2(lx,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lz,iy) = derv2(lz,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lx,iz) = derv2(lx,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ly,iz) = derv2(ly,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      else
        derv2(ix,lx) = derv2(ix,lx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(iy,lx) = derv2(iy,lx) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(iz,lx) = derv2(iz,lx) - d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ix,ly) = derv2(ix,ly) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(iy,ly) = derv2(iy,ly) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(iz,ly) = derv2(iz,ly) - d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(ix,lz) = derv2(ix,lz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(iy,lz) = derv2(iy,lz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(iz,lz) = derv2(iz,lz) - d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
!
!  K-L
!
      if (l.ge.k) then
        derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      else
        derv2(kx,lx) = derv2(kx,lx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
        derv2(ky,lx) = derv2(ky,lx) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
        derv2(kz,lx) = derv2(kz,lx) + d3trm*znbr(ni,i)*xnbr(ni2,i)
        derv2(kx,ly) = derv2(kx,ly) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
        derv2(ky,ly) = derv2(ky,ly) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
        derv2(kz,ly) = derv2(kz,ly) + d3trm*znbr(ni,i)*ynbr(ni2,i)
        derv2(kx,lz) = derv2(kx,lz) + d3trm*xnbr(ni,i)*znbr(ni2,i)
        derv2(ky,lz) = derv2(ky,lz) + d3trm*ynbr(ni,i)*znbr(ni2,i)
        derv2(kz,lz) = derv2(kz,lz) + d3trm*znbr(ni,i)*znbr(ni2,i)
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_selfp')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_dqp(i,iloc,d2Edlogcnidq,dlogcndcn)
!
!  Computes the second derivatives of the log coordination number and charge for GFNFF
!  Phonon version.
!
!  On entry : 
!
!  d2Edlogcnidq      = second derivatives of energy w.r.t. log of the coordination number for i and charge i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!
!   3/21 Created from gfnff_drv2_dcn_dq
!   1/22 iloc argument added for parallel case
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
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: iloc
  real(dp),    intent(in)                          :: d2Edlogcnidq
  real(dp),    intent(in)                          :: dlogcndcn(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  real(dp)                                         :: d1trm
  real(dp)                                         :: dcnidrij 
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dijx
  real(dp)                                         :: dijy
  real(dp)                                         :: dijz
  real(dp)                                         :: dikx
  real(dp)                                         :: diky
  real(dp)                                         :: dikz
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_dqp')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d1trm = d2Edlogcnidq*dlogcnidcni
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
!  Loop over coordination number pairs
!
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    j = nbrno(ni,i)
    indj = 3*(j-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
    dcnidrij = d1cndr_cn(n,i)
    d1trm = d2Edlogcnidq*dlogcnidcni*dcnidrij
    dijx = d1trm*xnbr(ni,i)
    dijy = d1trm*ynbr(ni,i)
    dijz = d1trm*znbr(ni,i)
!
!  Loop over atoms for charge derivatives
!
    do k = 1,numat
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
      dikx = dqdxyz(indk+1,iloc)
      diky = dqdxyz(indk+2,iloc)
      dikz = dqdxyz(indk+3,iloc)
!
!  NB: If order of indices of derv2 are changed then swap the j/k for the term that is added
!
      if (kx.gt.ix) then
        derv2(kx,ix) = derv2(kx,ix) - dijx*dikx
        derv2(ky,ix) = derv2(ky,ix) - dijx*diky
        derv2(kz,ix) = derv2(kz,ix) - dijx*dikz
        derv2(kx,iy) = derv2(kx,iy) - dijy*dikx
        derv2(ky,iy) = derv2(ky,iy) - dijy*diky
        derv2(kz,iy) = derv2(kz,iy) - dijy*dikz
        derv2(kx,iz) = derv2(kx,iz) - dijz*dikx
        derv2(ky,iz) = derv2(ky,iz) - dijz*diky
        derv2(kz,iz) = derv2(kz,iz) - dijz*dikz
      else
        derv2(ix,kx) = derv2(ix,kx) - dikx*dijx
        derv2(iy,kx) = derv2(iy,kx) - diky*dijx
        derv2(iz,kx) = derv2(iz,kx) - dikz*dijx
        derv2(ix,ky) = derv2(ix,ky) - dikx*dijy
        derv2(iy,ky) = derv2(iy,ky) - diky*dijy
        derv2(iz,ky) = derv2(iz,ky) - dikz*dijy
        derv2(ix,kz) = derv2(ix,kz) - dikx*dijz
        derv2(iy,kz) = derv2(iy,kz) - diky*dijz
        derv2(iz,kz) = derv2(iz,kz) - dikz*dijz
      endif
!
      if (kx.gt.jx) then
        derv2(kx,jx) = derv2(kx,jx) + dijx*dikx
        derv2(ky,jx) = derv2(ky,jx) + dijx*diky
        derv2(kz,jx) = derv2(kz,jx) + dijx*dikz
        derv2(kx,jy) = derv2(kx,jy) + dijy*dikx
        derv2(ky,jy) = derv2(ky,jy) + dijy*diky
        derv2(kz,jy) = derv2(kz,jy) + dijy*dikz
        derv2(kx,jz) = derv2(kx,jz) + dijz*dikx
        derv2(ky,jz) = derv2(ky,jz) + dijz*diky
        derv2(kz,jz) = derv2(kz,jz) + dijz*dikz
      else
        derv2(jx,kx) = derv2(jx,kx) + dikx*dijx
        derv2(jy,kx) = derv2(jy,kx) + diky*dijx
        derv2(jz,kx) = derv2(jz,kx) + dikz*dijx
        derv2(jx,ky) = derv2(jx,ky) + dikx*dijy
        derv2(jy,ky) = derv2(jy,ky) + diky*dijy
        derv2(jz,ky) = derv2(jz,ky) + dikz*dijy
        derv2(jx,kz) = derv2(jx,kz) + dikx*dijz
        derv2(jy,kz) = derv2(jy,kz) + diky*dijz
        derv2(jz,kz) = derv2(jz,kz) + dikz*dijz
      endif
!
      if (ix.gt.jx) then
        derv2(ix,jx) = derv2(ix,jx) - dijx*dikx
        derv2(iy,jx) = derv2(iy,jx) - dijx*diky
        derv2(iz,jx) = derv2(iz,jx) - dijx*dikz
        derv2(ix,jy) = derv2(ix,jy) - dijy*dikx
        derv2(iy,jy) = derv2(iy,jy) - dijy*diky
        derv2(iz,jy) = derv2(iz,jy) - dijy*dikz
        derv2(ix,jz) = derv2(ix,jz) - dijz*dikx
        derv2(iy,jz) = derv2(iy,jz) - dijz*diky
        derv2(iz,jz) = derv2(iz,jz) - dijz*dikz
      else
        derv2(jx,ix) = derv2(jx,ix) - dikx*dijx
        derv2(jy,ix) = derv2(jy,ix) - diky*dijx
        derv2(jz,ix) = derv2(jz,ix) - dikz*dijx
        derv2(jx,iy) = derv2(jx,iy) - dikx*dijy
        derv2(jy,iy) = derv2(jy,iy) - diky*dijy
        derv2(jz,iy) = derv2(jz,iy) - dikz*dijy
        derv2(jx,iz) = derv2(jx,iz) - dikx*dijz
        derv2(jy,iz) = derv2(jy,iz) - diky*dijz
        derv2(jz,iz) = derv2(jz,iz) - dikz*dijz
      endif
!
!  End of loop over k
!
    enddo
!
!  End of loop over j
!
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_dqp')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_selfpd_p1(iloc,i,dEdlogcni,d2Edlogcni2,dlogcndcn,d2logcndcn2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Special case version of gfnff_drv2_dcn for i only self terms. 
!  Phonon version.
!  Distributed memory parallel version.
!  Part 1 where i is local to the processor.
!
!  On entry : 
!
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   2/22 Created from gfnff_drv2_dcn_selfp
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
!  Julian Gale, CIC, Curtin University, February 2022
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: iloc
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: k
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2r2ikdx2(6)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_selfpd_p1')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d2logcnidcni2 = d2logcndcn2(i)
!
!  Set up for second derivatives
!
  indi = 3*(iloc-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
!
    d2r2ikdx2(1) = xnbr(ni,i)*xnbr(ni,i)
    d2r2ikdx2(2) = ynbr(ni,i)*ynbr(ni,i)
    d2r2ikdx2(3) = znbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(4) = ynbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(5) = xnbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(6) = xnbr(ni,i)*ynbr(ni,i)
!
    indk = 3*(k-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!
    d2cnidrik2 = d2cndr_cn(n,i)
    d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
           (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
    derv2(kx,ix) = derv2(kx,ix) - d2trm*d2r2ikdx2(1)
    derv2(ky,ix) = derv2(ky,ix) - d2trm*d2r2ikdx2(6)
    derv2(kz,ix) = derv2(kz,ix) - d2trm*d2r2ikdx2(5)
    derv2(kx,iy) = derv2(kx,iy) - d2trm*d2r2ikdx2(6)
    derv2(ky,iy) = derv2(ky,iy) - d2trm*d2r2ikdx2(2)
    derv2(kz,iy) = derv2(kz,iy) - d2trm*d2r2ikdx2(4)
    derv2(kx,iz) = derv2(kx,iz) - d2trm*d2r2ikdx2(5)
    derv2(ky,iz) = derv2(ky,iz) - d2trm*d2r2ikdx2(4)
    derv2(kz,iz) = derv2(kz,iz) - d2trm*d2r2ikdx2(3)
    derv2(kx,ix) = derv2(kx,ix) - d1trm
    derv2(ky,iy) = derv2(ky,iy) - d1trm
    derv2(kz,iz) = derv2(kz,iz) - d1trm
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,n-1
      ni2 = nbrno_cn(n2,i)
      l = nbrno(ni2,i)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnidril = d1cndr_cn(n2,i)
      d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
!  I-K
!
      derv2(kx,ix) = derv2(kx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
      derv2(ky,ix) = derv2(ky,ix) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
      derv2(kz,ix) = derv2(kz,ix) - d3trm*znbr(ni,i)*xnbr(ni2,i)
      derv2(kx,iy) = derv2(kx,iy) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
      derv2(ky,iy) = derv2(ky,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
      derv2(kz,iy) = derv2(kz,iy) - d3trm*znbr(ni,i)*ynbr(ni2,i)
      derv2(kx,iz) = derv2(kx,iz) - d3trm*xnbr(ni,i)*znbr(ni2,i)
      derv2(ky,iz) = derv2(ky,iz) - d3trm*ynbr(ni,i)*znbr(ni2,i)
      derv2(kz,iz) = derv2(kz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
!
!  I-L
!
      derv2(lx,ix) = derv2(lx,ix) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
      derv2(ly,ix) = derv2(ly,ix) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
      derv2(lz,ix) = derv2(lz,ix) - d3trm*xnbr(ni,i)*znbr(ni2,i)
      derv2(lx,iy) = derv2(lx,iy) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
      derv2(ly,iy) = derv2(ly,iy) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
      derv2(lz,iy) = derv2(lz,iy) - d3trm*ynbr(ni,i)*znbr(ni2,i)
      derv2(lx,iz) = derv2(lx,iz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
      derv2(ly,iz) = derv2(ly,iz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
      derv2(lz,iz) = derv2(lz,iz) - d3trm*znbr(ni,i)*znbr(ni2,i)
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_selfpd_p1')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_selfpd_p2(i,dEdlogcni,d2Edlogcni2,dlogcndcn,d2logcndcn2)
!
!  Computes the derivatives of the log coordination number for the neighbour list in GFNFF
!  Special case version of gfnff_drv2_dcn for i only self terms. 
!  Phonon version.
!  Distributed memory parallel version.
!  Part 2 where i is global
!
!  On entry : 
!
!  dEdlogcni         = derivatives of energy w.r.t. log of the coordination number for i
!  d2Edlogcni2       = second derivatives of energy w.r.t. log of the coordination number for i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!  d2logcndcn2       = second derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  dlogcn            = log coordination number derivatives w.r.t. coordinates if lgrad1 is true
!  dlogcns           = log coordination number derivatives w.r.t. strain if lgrad1 is true
!
!   2/22 Created from gfnff_drv2_dcn_selfp
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
!  Julian Gale, CIC, Curtin University, February 2022
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: dEdlogcni
  real(dp),    intent(in)                          :: d2Edlogcni2
  real(dp),    intent(in)                          :: dlogcndcn(*)
  real(dp),    intent(in)                          :: d2logcndcn2(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: l
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: n
  integer(i4)                                      :: n2
  integer(i4)                                      :: ni
  integer(i4)                                      :: ni2
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d3trm
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnidril
  real(dp)                                         :: d2cnidrik2
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: d2logcnidcni2
  real(dp)                                         :: d2r2ikdx2(6)
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_selfpd_p2')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d2logcnidcni2 = d2logcndcn2(i)
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!-----------------------------------------------------------
!  Compute derivatives of total coordination number for i  |
!-----------------------------------------------------------
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    k = nbrno(ni,i)
    kloc = atom2local(k)
!
!  k is not local then cycle
!
    if (kloc.eq.0) cycle
!
    dcnidrik = d1cndr_cn(n,i)
    d1trm = dEdlogcni*dlogcnidcni*dcnidrik
!
    d2r2ikdx2(1) = xnbr(ni,i)*xnbr(ni,i)
    d2r2ikdx2(2) = ynbr(ni,i)*ynbr(ni,i)
    d2r2ikdx2(3) = znbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(4) = ynbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(5) = xnbr(ni,i)*znbr(ni,i)
    d2r2ikdx2(6) = xnbr(ni,i)*ynbr(ni,i)
!
    indk = 3*(kloc-1)
    kx = indk + 1
    ky = indk + 2
    kz = indk + 3
!
    d2cnidrik2 = d2cndr_cn(n,i)
    d2trm = dEdlogcni*dlogcnidcni*d2cnidrik2 + &
           (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidrik
!----------------------------------------
!  Second derivatives of CN w.r.t. rik  |
!----------------------------------------
    derv2(ix,kx) = derv2(ix,kx) - d2trm*d2r2ikdx2(1)
    derv2(iy,kx) = derv2(iy,kx) - d2trm*d2r2ikdx2(6)
    derv2(iz,kx) = derv2(iz,kx) - d2trm*d2r2ikdx2(5)
    derv2(ix,ky) = derv2(ix,ky) - d2trm*d2r2ikdx2(6)
    derv2(iy,ky) = derv2(iy,ky) - d2trm*d2r2ikdx2(2)
    derv2(iz,ky) = derv2(iz,ky) - d2trm*d2r2ikdx2(4)
    derv2(ix,kz) = derv2(ix,kz) - d2trm*d2r2ikdx2(5)
    derv2(iy,kz) = derv2(iy,kz) - d2trm*d2r2ikdx2(4)
    derv2(iz,kz) = derv2(iz,kz) - d2trm*d2r2ikdx2(3)
    derv2(ix,kx) = derv2(ix,kx) - d1trm
    derv2(iy,ky) = derv2(iy,ky) - d1trm
    derv2(iz,kz) = derv2(iz,kz) - d1trm
!-------------------------------------------------------------------------------------
!  Loop over second neighbour of i for second derivatives for 2 different distances  |
!-------------------------------------------------------------------------------------
    do n2 = 1,nnbr_cn(i)
!
!  Skip n2 = n
!
      if (n2.eq.n) cycle
!
      ni2 = nbrno_cn(n2,i)
      l = nbrno(ni2,i)
      indl = 3*(l-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
!
      dcnidril = d1cndr_cn(n2,i)
      d3trm = (d2Edlogcni2*dlogcnidcni*dlogcnidcni + dEdlogcni*d2logcnidcni2)*dcnidrik*dcnidril
!
!  I-K
!
      derv2(ix,kx) = derv2(ix,kx) - d3trm*xnbr(ni,i)*xnbr(ni2,i)
      derv2(iy,kx) = derv2(iy,kx) - d3trm*xnbr(ni,i)*ynbr(ni2,i)
      derv2(iz,kx) = derv2(iz,kx) - d3trm*xnbr(ni,i)*znbr(ni2,i)
      derv2(ix,ky) = derv2(ix,ky) - d3trm*ynbr(ni,i)*xnbr(ni2,i)
      derv2(iy,ky) = derv2(iy,ky) - d3trm*ynbr(ni,i)*ynbr(ni2,i)
      derv2(iz,ky) = derv2(iz,ky) - d3trm*ynbr(ni,i)*znbr(ni2,i)
      derv2(ix,kz) = derv2(ix,kz) - d3trm*znbr(ni,i)*xnbr(ni2,i)
      derv2(iy,kz) = derv2(iy,kz) - d3trm*znbr(ni,i)*ynbr(ni2,i)
      derv2(iz,kz) = derv2(iz,kz) - d3trm*znbr(ni,i)*znbr(ni2,i)
!
!  K-L
!
      derv2(lx,kx) = derv2(lx,kx) + d3trm*xnbr(ni,i)*xnbr(ni2,i)
      derv2(ly,kx) = derv2(ly,kx) + d3trm*xnbr(ni,i)*ynbr(ni2,i)
      derv2(lz,kx) = derv2(lz,kx) + d3trm*xnbr(ni,i)*znbr(ni2,i)
      derv2(lx,ky) = derv2(lx,ky) + d3trm*ynbr(ni,i)*xnbr(ni2,i)
      derv2(ly,ky) = derv2(ly,ky) + d3trm*ynbr(ni,i)*ynbr(ni2,i)
      derv2(lz,ky) = derv2(lz,ky) + d3trm*ynbr(ni,i)*znbr(ni2,i)
      derv2(lx,kz) = derv2(lx,kz) + d3trm*znbr(ni,i)*xnbr(ni2,i)
      derv2(ly,kz) = derv2(ly,kz) + d3trm*znbr(ni,i)*ynbr(ni2,i)
      derv2(lz,kz) = derv2(lz,kz) + d3trm*znbr(ni,i)*znbr(ni2,i)
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_selfpd_p2')
#endif
!
  return
  end
!
!
  subroutine gfnff_drv2_dcn_dqpd_p1(iloc,i,d2Edlogcnidq,dlogcndcn)
!
!  Computes the second derivatives of the log coordination number and charge for GFNFF
!  Phonon version.
!  Distributed memory parallel version.
!  Part 1 where i is local to the processor.
!
!  On entry : 
!
!  d2Edlogcnidq      = second derivatives of energy w.r.t. log of the coordination number for i and charge i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!
!   2/22 Created from gfnff_drv2_dcn_dqp
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
!  Julian Gale, CIC, Curtin University, February 2022
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: iloc
  real(dp),    intent(in)                          :: d2Edlogcnidq
  real(dp),    intent(in)                          :: dlogcndcn(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  real(dp)                                         :: d1trm
  real(dp)                                         :: dcnidrij 
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dijx
  real(dp)                                         :: dijy
  real(dp)                                         :: dijz
  real(dp)                                         :: dikx
  real(dp)                                         :: diky
  real(dp)                                         :: dikz
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_dqpd_p1')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d1trm = d2Edlogcnidq*dlogcnidcni
!
!  Set up for second derivatives
!
  indi = 3*(iloc-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
!  Loop over coordination number pairs
!
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    j = nbrno(ni,i)
    indj = 3*(j-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
    dcnidrij = d1cndr_cn(n,i)
    d1trm = d2Edlogcnidq*dlogcnidcni*dcnidrij
    dijx = d1trm*xnbr(ni,i)
    dijy = d1trm*ynbr(ni,i)
    dijz = d1trm*znbr(ni,i)
!
!  Loop over atoms for charge derivatives
!
    do k = 1,numat
      indk = 3*(k-1)
      kx = indk + 1
      ky = indk + 2
      kz = indk + 3
      dikx = dqdxyz(indk+1,i)
      diky = dqdxyz(indk+2,i)
      dikz = dqdxyz(indk+3,i)
!
      if (k.gt.i) then
        derv2(kx,ix) = derv2(kx,ix) - dijx*dikx
        derv2(ky,ix) = derv2(ky,ix) - dijx*diky
        derv2(kz,ix) = derv2(kz,ix) - dijx*dikz
        derv2(kx,iy) = derv2(kx,iy) - dijy*dikx
        derv2(ky,iy) = derv2(ky,iy) - dijy*diky
        derv2(kz,iy) = derv2(kz,iy) - dijy*dikz
        derv2(kx,iz) = derv2(kx,iz) - dijz*dikx
        derv2(ky,iz) = derv2(ky,iz) - dijz*diky
        derv2(kz,iz) = derv2(kz,iz) - dijz*dikz
      else
        derv2(kx,ix) = derv2(kx,ix) - dijx*dikx
        derv2(ky,ix) = derv2(ky,ix) - dijy*dikx
        derv2(kz,ix) = derv2(kz,ix) - dijz*dikx
        derv2(kx,iy) = derv2(kx,iy) - dijx*diky
        derv2(ky,iy) = derv2(ky,iy) - dijy*diky
        derv2(kz,iy) = derv2(kz,iy) - dijz*diky
        derv2(kx,iz) = derv2(kx,iz) - dijx*dikz
        derv2(ky,iz) = derv2(ky,iz) - dijy*dikz
        derv2(kz,iz) = derv2(kz,iz) - dijz*dikz
      endif
!
      if (j.gt.i) then
        derv2(jx,ix) = derv2(jx,ix) - dikx*dijx
        derv2(jy,ix) = derv2(jy,ix) - diky*dijx
        derv2(jz,ix) = derv2(jz,ix) - dikz*dijx
        derv2(jx,iy) = derv2(jx,iy) - dikx*dijy
        derv2(jy,iy) = derv2(jy,iy) - diky*dijy
        derv2(jz,iy) = derv2(jz,iy) - dikz*dijy
        derv2(jx,iz) = derv2(jx,iz) - dikx*dijz
        derv2(jy,iz) = derv2(jy,iz) - diky*dijz
        derv2(jz,iz) = derv2(jz,iz) - dikz*dijz
      else
        derv2(jx,ix) = derv2(jx,ix) - dikx*dijx
        derv2(jy,ix) = derv2(jy,ix) - dikx*dijy
        derv2(jz,ix) = derv2(jz,ix) - dikx*dijz
        derv2(jx,iy) = derv2(jx,iy) - diky*dijx
        derv2(jy,iy) = derv2(jy,iy) - diky*dijy
        derv2(jz,iy) = derv2(jz,iy) - diky*dijz
        derv2(jx,iz) = derv2(jx,iz) - dikz*dijx
        derv2(jy,iz) = derv2(jy,iz) - dikz*dijy
        derv2(jz,iz) = derv2(jz,iz) - dikz*dijz
      endif
!
!  End of loop over k
!
    enddo
!
!  End of loop over j
!
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_dqpd_p1')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_dqpd_p2(i,d2Edlogcnidq,dlogcndcn)
!
!  Computes the second derivatives of the log coordination number and charge for GFNFF
!  Phonon version.
!  Distributed memory parallel version.
!  Part 2 where i is global
!
!  On entry : 
!
!  d2Edlogcnidq      = second derivatives of energy w.r.t. log of the coordination number for i and charge i
!  dlogcndcn         = derivative of log coordination numbers w.r.t. coordination numbers
!
!   2/22 Created from gfnff_drv2_dcn_dqp
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
!  Julian Gale, CIC, Curtin University, February 2022
!
  use datatypes
  use current
  use derivatives
  use m_gfnff_nbr
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: d2Edlogcnidq
  real(dp),    intent(in)                          :: dlogcndcn(*)
!
!  Local variables
!
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jloc
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: jxf
  integer(i4)                                      :: jyf
  integer(i4)                                      :: jzf
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: n
  integer(i4)                                      :: ni
  real(dp)                                         :: d1trm
  real(dp)                                         :: dcnidrij 
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dijx
  real(dp)                                         :: dijy
  real(dp)                                         :: dijz
  real(dp)                                         :: dikx
  real(dp)                                         :: diky
  real(dp)                                         :: dikz
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_dqpd_p2')
#endif
!******************************************************************************
!  Loop over pairs of atoms to compute coordination number derivatives for i  *
!******************************************************************************
  dlogcnidcni = dlogcndcn(i)
  d1trm = d2Edlogcnidq*dlogcnidcni
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
!  Loop over coordination number pairs
!
  do n = 1,nnbr_cn(i)
    ni = nbrno_cn(n,i)
    j = nbrno(ni,i)
    jloc = atom2local(j)
    if (jloc.ne.0) then
      indj = 3*(jloc-1)
      jx = indj + 1
      jy = indj + 2
      jz = indj + 3
    endif
!
    indj = 3*(j-1)
    jxf = indj + 1
    jyf = indj + 2
    jzf = indj + 3
!
    dcnidrij = d1cndr_cn(n,i)
    d1trm = d2Edlogcnidq*dlogcnidcni*dcnidrij
    dijx = d1trm*xnbr(ni,i)
    dijy = d1trm*ynbr(ni,i)
    dijz = d1trm*znbr(ni,i)
!
!  Loop over atoms for charge derivatives
!
    do k = 1,numat
      kloc = atom2local(k)
!
!  If neither j or k are local then cycle
!
      if (jloc+kloc.eq.0) cycle
!
      if (kloc.ne.0) then
        indk = 3*(kloc-1)
        kx = indk + 1
        ky = indk + 2
        kz = indk + 3
      endif
!
      indk = 3*(k-1)
      kxf = indk + 1
      kyf = indk + 2
      kzf = indk + 3
!
      dikx = dqdxyz(indk+1,i)
      diky = dqdxyz(indk+2,i)
      dikz = dqdxyz(indk+3,i)
!
!  NB: If order of indices of derv2 are changed then swap the j/k for the term that is added
!
      if (kloc.ne.0) then
        if (k.gt.i) then
          derv2(ix,kx) = derv2(ix,kx) - dikx*dijx
          derv2(iy,kx) = derv2(iy,kx) - dikx*dijy
          derv2(iz,kx) = derv2(iz,kx) - dikx*dijz
          derv2(ix,ky) = derv2(ix,ky) - diky*dijx
          derv2(iy,ky) = derv2(iy,ky) - diky*dijy
          derv2(iz,ky) = derv2(iz,ky) - diky*dijz
          derv2(ix,kz) = derv2(ix,kz) - dikz*dijx
          derv2(iy,kz) = derv2(iy,kz) - dikz*dijy
          derv2(iz,kz) = derv2(iz,kz) - dikz*dijz
        else
          derv2(ix,kx) = derv2(ix,kx) - dikx*dijx
          derv2(iy,kx) = derv2(iy,kx) - diky*dijx
          derv2(iz,kx) = derv2(iz,kx) - dikz*dijx
          derv2(ix,ky) = derv2(ix,ky) - dikx*dijy
          derv2(iy,ky) = derv2(iy,ky) - diky*dijy
          derv2(iz,ky) = derv2(iz,ky) - dikz*dijy
          derv2(ix,kz) = derv2(ix,kz) - dikx*dijz
          derv2(iy,kz) = derv2(iy,kz) - diky*dijz
          derv2(iz,kz) = derv2(iz,kz) - dikz*dijz
        endif
      endif
!
      if (jloc.ne.0) then
        if (k.gt.j) then
          derv2(kxf,jx) = derv2(kxf,jx) + dijx*dikx
          derv2(kyf,jx) = derv2(kyf,jx) + dijx*diky
          derv2(kzf,jx) = derv2(kzf,jx) + dijx*dikz
          derv2(kxf,jy) = derv2(kxf,jy) + dijy*dikx
          derv2(kyf,jy) = derv2(kyf,jy) + dijy*diky
          derv2(kzf,jy) = derv2(kzf,jy) + dijy*dikz
          derv2(kxf,jz) = derv2(kxf,jz) + dijz*dikx
          derv2(kyf,jz) = derv2(kyf,jz) + dijz*diky
          derv2(kzf,jz) = derv2(kzf,jz) + dijz*dikz
        else
          derv2(kxf,jx) = derv2(kxf,jx) + dijx*dikx
          derv2(kyf,jx) = derv2(kyf,jx) + dijy*dikx
          derv2(kzf,jx) = derv2(kzf,jx) + dijz*dikx
          derv2(kxf,jy) = derv2(kxf,jy) + dijx*diky
          derv2(kyf,jy) = derv2(kyf,jy) + dijy*diky
          derv2(kzf,jy) = derv2(kzf,jy) + dijz*diky
          derv2(kxf,jz) = derv2(kxf,jz) + dijx*dikz
          derv2(kyf,jz) = derv2(kyf,jz) + dijy*dikz
          derv2(kzf,jz) = derv2(kzf,jz) + dijz*dikz
        endif
      endif
!
      if (kloc.ne.0) then
        if (j.gt.k) then
          derv2(jxf,kx) = derv2(jxf,kx) + dikx*dijx
          derv2(jyf,kx) = derv2(jyf,kx) + diky*dijx
          derv2(jzf,kx) = derv2(jzf,kx) + dikz*dijx
          derv2(jxf,ky) = derv2(jxf,ky) + dikx*dijy
          derv2(jyf,ky) = derv2(jyf,ky) + diky*dijy
          derv2(jzf,ky) = derv2(jzf,ky) + dikz*dijy
          derv2(jxf,kz) = derv2(jxf,kz) + dikx*dijz
          derv2(jyf,kz) = derv2(jyf,kz) + diky*dijz
          derv2(jzf,kz) = derv2(jzf,kz) + dikz*dijz
        else
          derv2(jxf,kx) = derv2(jxf,kx) + dikx*dijx
          derv2(jyf,kx) = derv2(jyf,kx) + dikx*dijy
          derv2(jzf,kx) = derv2(jzf,kx) + dikx*dijz
          derv2(jxf,ky) = derv2(jxf,ky) + diky*dijx
          derv2(jyf,ky) = derv2(jyf,ky) + diky*dijy
          derv2(jzf,ky) = derv2(jzf,ky) + diky*dijz
          derv2(jxf,kz) = derv2(jxf,kz) + dikz*dijx
          derv2(jyf,kz) = derv2(jyf,kz) + dikz*dijy
          derv2(jzf,kz) = derv2(jzf,kz) + dikz*dijz
        endif
      endif
!
      if (jloc.ne.0) then
        if (j.gt.i) then
          derv2(ix,jx) = derv2(ix,jx) - dijx*dikx
          derv2(iy,jx) = derv2(iy,jx) - dijy*dikx
          derv2(iz,jx) = derv2(iz,jx) - dijz*dikx
          derv2(ix,jy) = derv2(ix,jy) - dijx*diky
          derv2(iy,jy) = derv2(iy,jy) - dijy*diky
          derv2(iz,jy) = derv2(iz,jy) - dijz*diky
          derv2(ix,jz) = derv2(ix,jz) - dijx*dikz
          derv2(iy,jz) = derv2(iy,jz) - dijy*dikz
          derv2(iz,jz) = derv2(iz,jz) - dijz*dikz
        else
          derv2(ix,jx) = derv2(ix,jx) - dijx*dikx
          derv2(iy,jx) = derv2(iy,jx) - dijx*diky
          derv2(iz,jx) = derv2(iz,jx) - dijx*dikz
          derv2(ix,jy) = derv2(ix,jy) - dijy*dikx
          derv2(iy,jy) = derv2(iy,jy) - dijy*diky
          derv2(iz,jy) = derv2(iz,jy) - dijy*dikz
          derv2(ix,jz) = derv2(ix,jz) - dijz*dikx
          derv2(iy,jz) = derv2(iy,jz) - dijz*diky
          derv2(iz,jz) = derv2(iz,jz) - dijz*dikz
        endif
      endif
!
!  End of loop over k
!
    enddo
!
!  End of loop over j
!
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_dqpd_p2')
#endif
!
  return
  end
