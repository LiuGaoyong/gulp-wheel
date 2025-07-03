  subroutine valencebondpd(xkv,ykv,zkv)
!
!  Calculates the second derivatives for the Valance Bond potentials.
!  Phonon version. Parallel distributed memory.
!
!  On entry : 
!
!  xkv             = x kvector component
!  ykv             = y kvector component
!  zkv             = z kvector component
!
!  12/21 Created from valencebondp
!   2/22 Molecule centre of mass handling added
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
  use bondvalence
  use control,        only : lrigid
  use current
  use derivatives
  use element,        only : maxele
  use iochannels
  use m_strain,       only : real1strterm, cartstrterm
  use m_vb_nbr
  use molecule
  use parallel
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(in)                          :: xkv
  real(dp),    intent(in)                          :: ykv
  real(dp),    intent(in)                          :: zkv
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: iloc
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: ixf
  integer(i4)                                      :: iyf
  integer(i4)                                      :: izf
  integer(i4)                                      :: j
  integer(i4)                                      :: jloc
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nmi
  integer(i4)                                      :: nmj
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: ntypj
  integer(i4)                                      :: nv
  integer(i4)                                      :: nvb
  integer(i4)                                      :: nvbij
  integer(i4), dimension(:),     allocatable, save :: nvbijptr
  integer(i4)                                      :: nve
  integer(i4)                                      :: nvi
  integer(i4)                                      :: status
  logical                                          :: lanylocal
  real(dp)                                         :: bvsum
  real(dp)                                         :: bvvsum(3)
  real(dp)                                         :: g_cpu_time
  real(dp),    dimension(:),     allocatable, save :: d1i
  real(dp),    dimension(:),     allocatable, save :: d2i
  real(dp)                                         :: d1wji(3,3)
  real(dp)                                         :: d1wki(3,3)
  real(dp)                                         :: d2wji(3,3,3)
  real(dp)                                         :: dwji2(3)
  real(dp)                                         :: dwki2(3)
  real(dp)                                         :: d2wji2(3,3)
  real(dp)                                         :: d2wkj2(3,3)
  real(dp)                                         :: deidbvsum
  real(dp)                                         :: d2eidbvsum2
  real(dp)                                         :: deidwi2
  real(dp)                                         :: d2eidwi22
  real(dp)                                         :: deijdr
  real(dp)                                         :: d2eidr2
  real(dp)                                         :: drij2ds(6)
  real(dp)                                         :: d2rij2dx2(6)
  real(dp)                                         :: d2rij2ds2(6,6)
  real(dp)                                         :: d2rij2dsdx(6,3)
  real(dp)                                         :: dxyzijds(6,3)
  real(dp)                                         :: d2xyzijdsdx(6,3,3)
  real(dp)                                         :: d2xyzijds2(6,6,3)
  real(dp)                                         :: cosij
  real(dp)                                         :: cosik
  real(dp)                                         :: cosjk
  real(dp)                                         :: sinij
  real(dp)                                         :: sinik
  real(dp)                                         :: sinjk
  real(dp)                                         :: oneij
  real(dp)                                         :: oneik
  real(dp)                                         :: onejk
  real(dp)                                         :: rij
  real(dp)                                         :: rrij
  real(dp)                                         :: rrij3
  real(dp)                                         :: rrij5
  real(dp)                                         :: rik
  real(dp)                                         :: rrik
  real(dp)                                         :: rrik3
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp)                                         :: vij
  real(dp)                                         :: d1vijdr
  real(dp)                                         :: d2vijdr
  real(dp)                                         :: wi2
  real(dp)                                         :: xcom
  real(dp)                                         :: ycom
  real(dp)                                         :: zcom
  real(dp)                                         :: xcomi
  real(dp)                                         :: ycomi
  real(dp)                                         :: zcomi
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xki
  real(dp)                                         :: yki
  real(dp)                                         :: zki
  real(dp),    dimension(:),     allocatable, save :: bv
  real(dp),    dimension(:,:),   allocatable, save :: bvv
#ifdef TRACE
  call trace_in('valencebondpd')
#endif
!
  t1 = g_cpu_time()
!***************************************
!  Generate neighbour lists for atoms  *
!***************************************
  call vb_getnbr
!
!  Allocate local memory 
!
  allocate(nvbijptr(maxvalbond),stat=status)
  if (status/=0) call outofmemory('valencebondpd','nvbijptr')
  allocate(bv(maxvbnbr),stat=status)
  if (status/=0) call outofmemory('valencebondpd','vb')
  allocate(bvv(3,maxvbnbr),stat=status)
  if (status/=0) call outofmemory('valencebondpd','vbv')
!
  allocate(d1i(maxvbnbr),stat=status)
  if (status/=0) call outofmemory('valencebondpd','d1i')
  allocate(d2i(maxvbnbr),stat=status)   
  if (status/=0) call outofmemory('valencebondpd','d2i')
!***************************************************
!  Loop over atoms to compute valence bond energy  *
!***************************************************
  ix = - 2
  iy = - 1
  iz =   0
  do iloc = 1,natomsonnode
    i = node2atom(iloc)
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!
    ixf = 3*(i-1) + 1
    iyf = ixf + 1
    izf = ixf + 2
!
!  Set variables relating to i
!
    nati = nat(i)
    ntypi = nftype(i)
    nmi = natmol(i)
!
!  Set COM coordinates
!
    if (lrigid.and.nmi.gt.0) then
      xcomi = molxyz(1,natinmol(i),nmi)
      ycomi = molxyz(2,natinmol(i),nmi)
      zcomi = molxyz(3,natinmol(i),nmi)
    else
      xcomi = 0.0_dp
      ycomi = 0.0_dp
      zcomi = 0.0_dp
    endif
!
!  Initialise bond valence and bond valence vector
!
    bv(1:nvbnbr(i)) = 0.0_dp
    bvv(1:3,1:nvbnbr(i)) = 0.0_dp
!
!  Set total number of distances for neighbours of i
!
    nneighi2 = nvbnbr(i) + nvbnbr(i)*(nvbnbr(i) + 1)/2 
!
!  Initialise derivative storage for neighbours of i
!
    d1i(1:nvbnbr(i)) = 0.0_dp
    d2i(1:nvbnbr(i)) = 0.0_dp
!
!  Loop over neighbours of i
!
    do ni = 1,nvbnbr(i)
      j = nvbnbrno(ni,i)
!
!  Set variables relating to j
!
      natj = nat(j)
      ntypj = nftype(j)
!
!  Set up i-j quantities
!
      rij = rvbnbr(ni,i)
      xji = xvbnbr(ni,i)
      yji = yvbnbr(ni,i)
      zji = zvbnbr(ni,i)
      rrij = 1.0_dp/rij
!
!  Find valence bond potential between i and j
!
      nvbij = 0
      do nvi = 1,nvalbondatm(i)
        nv = nvalbondatmptr(nvi,i)
        if (nVBspecB2(nv).eq.maxele.or.(nVBspecB2(nv).eq.natj.and.(nVBtypeB2(nv).eq.ntypj.or.nVBtypeB2(nv).eq.0))) then
          nvbij = nvbij + 1
          nvbijptr(nvbij) = nv
        endif
      enddo
      if (nvbij.gt.0) then
!*********************************
!  Calculate valence bond terms  *
!*********************************
        do nvb = 1,nvbij
          nv = nvbijptr(nvb)
          if (rij.ge.rVBmin(nv).and.rij.le.rVBmax(nv)) then
            if (nvalbondtype(nv).eq.1) then
!
!  Exponential form
!
              vij = VBwgtB(nv)*exp((VBparB(1,nv) - rij)/VBparB(3,nv))
              d1vijdr = - vij*rrij/VBparB(3,nv)
              d2vijdr = - d1vijdr*rrij*(1.0_dp/VBparB(3,nv) + rrij)
            else
!
!  Power law form
!
              vij = VBwgtB(nv)*(VBparB(1,nv)/rij)**VBparB(2,nv)
              d1vijdr = - VBparB(2,nv)*vij*rrij*rrij
              d2vijdr = - (VBparB(2,nv) + 2.0_dp)*d1vijdr*rrij*rrij
            endif
            bv(ni) = bv(ni) + vij
            bvv(1,ni) = bvv(1,ni) + vij*xji*rrij
            bvv(2,ni) = bvv(2,ni) + vij*yji*rrij
            bvv(3,ni) = bvv(3,ni) + vij*zji*rrij
!
            d1i(ni) = d1i(ni) + d1vijdr
            d2i(ni) = d2i(ni) + d2vijdr
          endif
        enddo
      endif
    enddo
!**********************************
!  Calculate valence bond energy  *
!**********************************
    bvsum = 0.0_dp
    bvvsum(1:3) = 0.0_dp
    do ni = 1,nvbnbr(i)
      bvsum = bvsum + bv(ni)
      bvvsum(1:3) = bvvsum(1:3) + bvv(1:3,ni)
    enddo
    wi2 = bvvsum(1)**2 + bvvsum(2)**2 + bvvsum(3)**2
!***************************************
!  Derivatives of valence bond energy  *
!***************************************
!
!  Compute derivatives of energy with respect to bvsum / wi2
!
    deidbvsum   = 0.0_dp
    d2eidbvsum2 = 0.0_dp
    deidwi2     = 0.0_dp
    d2eidwi22   = 0.0_dp
    do nve = 1,nvalener
      if (nVBspecE(nve).eq.nati) then
        if (nVBtypeE(nve).eq.ntypi.or.nVBtypeE(nve).eq.0) then
          deidbvsum   = deidbvsum   + 2.0_dp*VBparE(1,nve)*(bvsum - VBparE(2,nve))
          deidwi2     = deidwi2     + 2.0_dp*VBparE(3,nve)*(wi2 - VBparE(4,nve)**2)
          d2eidbvsum2 = d2eidbvsum2 + 2.0_dp*VBparE(1,nve)
          d2eidwi22   = d2eidwi22   + 2.0_dp*VBparE(3,nve)
        endif
      endif
    enddo
!
!  Loop over neighbours of i
!
    do ni = 1,nvbnbr(i)
      j = nvbnbrno(ni,i)
!
      jx = 3*(j-1) + 1
      jy = jx + 1
      jz = jx + 2
!
!  Set up i-j quantities
!
      rij = rvbnbr(ni,i)
      rrij = 1.0_dp/rij
      xji = xvbnbr(ni,i)
      yji = yvbnbr(ni,i)
      zji = zvbnbr(ni,i)
!
!  Set COM coordinates
!
      nmj = natmol(j)
      if (lrigid.and.nmj.gt.0) then
        xcom = molxyz(1,natinmol(j),nmj) - xcomi
        ycom = molxyz(2,natinmol(j),nmj) - ycomi
        zcom = molxyz(3,natinmol(j),nmj) - zcomi
      else
        xcom = - xcomi
        ycom = - ycomi
        zcom = - zcomi
      endif
!
!  Set diagonal blocks to be difference from gamma point
!
      if (i.eq.j) then
        oneij = 1.0_dp
      else
        oneij = 0.0_dp
      endif
!
!  Compute phase factor for i-j
!
      cosij = xkv*xji + ykv*yji + zkv*zji
      sinij = sin(cosij)
      cosij = cos(cosij) - oneij
!
!  Derivatives of E_bv
!
      deijdr = deidbvsum*d1i(ni)
!
      call real1strterm(ndim,xji,yji,zji,xcom,ycom,zcom,drij2ds,d2rij2dx2,d2rij2dsdx,d2rij2ds2,.true.)
!
!  Derivatives of E_bvv
!
      rrij3 = rrij**3
      d1wji(1,1) = d1i(ni)*xji*xji*rrij - bv(ni)*(rrij3*xji*xji - rrij)
      d1wji(2,1) = d1i(ni)*yji*xji*rrij - bv(ni)*rrij3*yji*xji
      d1wji(3,1) = d1i(ni)*zji*xji*rrij - bv(ni)*rrij3*zji*xji
!
      d1wji(1,2) = d1i(ni)*xji*yji*rrij - bv(ni)*rrij3*xji*yji
      d1wji(2,2) = d1i(ni)*yji*yji*rrij - bv(ni)*(rrij3*yji*yji - rrij)
      d1wji(3,2) = d1i(ni)*zji*yji*rrij - bv(ni)*rrij3*zji*yji
!
      d1wji(1,3) = d1i(ni)*xji*zji*rrij - bv(ni)*rrij3*xji*zji
      d1wji(2,3) = d1i(ni)*yji*zji*rrij - bv(ni)*rrij3*yji*zji
      d1wji(3,3) = d1i(ni)*zji*zji*rrij - bv(ni)*(rrij3*zji*zji - rrij)
!
      dwji2(1) = 2.0_dp*(bvvsum(1)*d1wji(1,1) + bvvsum(2)*d1wji(1,2) + bvvsum(3)*d1wji(1,3))
      dwji2(2) = 2.0_dp*(bvvsum(1)*d1wji(2,1) + bvvsum(2)*d1wji(2,2) + bvvsum(3)*d1wji(2,3))
      dwji2(3) = 2.0_dp*(bvvsum(1)*d1wji(3,1) + bvvsum(2)*d1wji(3,2) + bvvsum(3)*d1wji(3,3))
!
      call cartstrterm(ndim,xji,yji,zji,xcom,ycom,zcom,dxyzijds,d2xyzijdsdx,d2xyzijds2,.true.)
!
!  Second derivatives
!
      jx = 3*(j-1) + 1
      jy = jx + 1
      jz = jx + 2
!
!  Second derivatives of a single Vij
!
      d2eidr2 = deidbvsum*d2i(ni) + d2eidbvsum2*d1i(ni)**2
!
      derv2(jx,ix) = derv2(jx,ix) - d2eidr2*d2rij2dx2(1)*cosij
      derv2(jy,ix) = derv2(jy,ix) - d2eidr2*d2rij2dx2(6)*cosij
      derv2(jz,ix) = derv2(jz,ix) - d2eidr2*d2rij2dx2(5)*cosij
      derv2(jx,iy) = derv2(jx,iy) - d2eidr2*d2rij2dx2(6)*cosij
      derv2(jy,iy) = derv2(jy,iy) - d2eidr2*d2rij2dx2(2)*cosij
      derv2(jz,iy) = derv2(jz,iy) - d2eidr2*d2rij2dx2(4)*cosij
      derv2(jx,iz) = derv2(jx,iz) - d2eidr2*d2rij2dx2(5)*cosij
      derv2(jy,iz) = derv2(jy,iz) - d2eidr2*d2rij2dx2(4)*cosij
      derv2(jz,iz) = derv2(jz,iz) - d2eidr2*d2rij2dx2(3)*cosij
      derv2(jx,ix) = derv2(jx,ix) - deijdr*cosij
      derv2(jy,iy) = derv2(jy,iy) - deijdr*cosij
      derv2(jz,iz) = derv2(jz,iz) - deijdr*cosij
!
      dervi(jx,ix) = dervi(jx,ix) - d2eidr2*d2rij2dx2(1)*sinij
      dervi(jy,ix) = dervi(jy,ix) - d2eidr2*d2rij2dx2(6)*sinij
      dervi(jz,ix) = dervi(jz,ix) - d2eidr2*d2rij2dx2(5)*sinij
      dervi(jx,iy) = dervi(jx,iy) - d2eidr2*d2rij2dx2(6)*sinij
      dervi(jy,iy) = dervi(jy,iy) - d2eidr2*d2rij2dx2(2)*sinij
      dervi(jz,iy) = dervi(jz,iy) - d2eidr2*d2rij2dx2(4)*sinij
      dervi(jx,iz) = dervi(jx,iz) - d2eidr2*d2rij2dx2(5)*sinij
      dervi(jy,iz) = dervi(jy,iz) - d2eidr2*d2rij2dx2(4)*sinij
      dervi(jz,iz) = dervi(jz,iz) - d2eidr2*d2rij2dx2(3)*sinij
      dervi(jx,ix) = dervi(jx,ix) - deijdr*sinij
      dervi(jy,iy) = dervi(jy,iy) - deijdr*sinij
      dervi(jz,iz) = dervi(jz,iz) - deijdr*sinij
!
!  Second derivatives of bvv
!
      rrij5 = rrij**5
      d2wji(1,1,1) = xji*xji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     3.0_dp*xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,1,1) = yji*xji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,1,1) = zji*xji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,2,1) = xji*yji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,2,1) = yji*yji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,2,1) = zji*yji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(1,3,1) = xji*zji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,3,1) = yji*zji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(3,3,1) = zji*zji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,1,2) = xji*xji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,1,2) = yji*xji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,1,2) = zji*xji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(1,2,2) = xji*yji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,2,2) = yji*yji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     3.0_dp*yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,2,2) = zji*yji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,3,2) = xji*zji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(2,3,2) = yji*zji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,3,2) = zji*zji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,1,3) = xji*xji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,1,3) = yji*xji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(3,1,3) = zji*xji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,2,3) = xji*yji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(2,2,3) = yji*yji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,2,3) = zji*yji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,3,3) = xji*zji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,3,3) = yji*zji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,3,3) = zji*zji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     3.0_dp*zji*(d1i(ni)*rrij - bv(ni)*rrij3)
!
      d2wji2(1,1) = 2.0_dp*(bvvsum(1)*d2wji(1,1,1) + bvvsum(2)*d2wji(1,1,2) + bvvsum(3)*d2wji(1,1,3) + &
                            d1wji(1,1)*d1wji(1,1) + d1wji(1,2)*d1wji(1,2) + d1wji(1,3)*d1wji(1,3))
      d2wji2(2,1) = 2.0_dp*(bvvsum(1)*d2wji(2,1,1) + bvvsum(2)*d2wji(2,1,2) + bvvsum(3)*d2wji(2,1,3) + &
                            d1wji(2,1)*d1wji(1,1) + d1wji(2,2)*d1wji(1,2) + d1wji(2,3)*d1wji(1,3))
      d2wji2(3,1) = 2.0_dp*(bvvsum(1)*d2wji(3,1,1) + bvvsum(2)*d2wji(3,1,2) + bvvsum(3)*d2wji(3,1,3) + &
                            d1wji(3,1)*d1wji(1,1) + d1wji(3,2)*d1wji(1,2) + d1wji(3,3)*d1wji(1,3))
      d2wji2(1,2) = 2.0_dp*(bvvsum(1)*d2wji(1,2,1) + bvvsum(2)*d2wji(1,2,2) + bvvsum(3)*d2wji(1,2,3) + &
                            d1wji(1,1)*d1wji(2,1) + d1wji(1,2)*d1wji(2,2) + d1wji(1,3)*d1wji(2,3))
      d2wji2(2,2) = 2.0_dp*(bvvsum(1)*d2wji(2,2,1) + bvvsum(2)*d2wji(2,2,2) + bvvsum(3)*d2wji(2,2,3) + &
                            d1wji(2,1)*d1wji(2,1) + d1wji(2,2)*d1wji(2,2) + d1wji(2,3)*d1wji(2,3))
      d2wji2(3,2) = 2.0_dp*(bvvsum(1)*d2wji(3,2,1) + bvvsum(2)*d2wji(3,2,2) + bvvsum(3)*d2wji(3,2,3) + &
                            d1wji(3,1)*d1wji(2,1) + d1wji(3,2)*d1wji(2,2) + d1wji(3,3)*d1wji(2,3))
      d2wji2(1,3) = 2.0_dp*(bvvsum(1)*d2wji(1,3,1) + bvvsum(2)*d2wji(1,3,2) + bvvsum(3)*d2wji(1,3,3) + &
                            d1wji(1,1)*d1wji(3,1) + d1wji(1,2)*d1wji(3,2) + d1wji(1,3)*d1wji(3,3))
      d2wji2(2,3) = 2.0_dp*(bvvsum(1)*d2wji(2,3,1) + bvvsum(2)*d2wji(2,3,2) + bvvsum(3)*d2wji(2,3,3) + &
                            d1wji(2,1)*d1wji(3,1) + d1wji(2,2)*d1wji(3,2) + d1wji(2,3)*d1wji(3,3))
      d2wji2(3,3) = 2.0_dp*(bvvsum(1)*d2wji(3,3,1) + bvvsum(2)*d2wji(3,3,2) + bvvsum(3)*d2wji(3,3,3) + &
                            d1wji(3,1)*d1wji(3,1) + d1wji(3,2)*d1wji(3,2) + d1wji(3,3)*d1wji(3,3))
!
      derv2(jx,ix) = derv2(jx,ix) - (d2eidwi22*dwji2(1)*dwji2(1) + deidwi2*d2wji2(1,1))*cosij
      derv2(jy,ix) = derv2(jy,ix) - (d2eidwi22*dwji2(2)*dwji2(1) + deidwi2*d2wji2(2,1))*cosij
      derv2(jz,ix) = derv2(jz,ix) - (d2eidwi22*dwji2(3)*dwji2(1) + deidwi2*d2wji2(3,1))*cosij
      derv2(jx,iy) = derv2(jx,iy) - (d2eidwi22*dwji2(1)*dwji2(2) + deidwi2*d2wji2(1,2))*cosij
      derv2(jy,iy) = derv2(jy,iy) - (d2eidwi22*dwji2(2)*dwji2(2) + deidwi2*d2wji2(2,2))*cosij
      derv2(jz,iy) = derv2(jz,iy) - (d2eidwi22*dwji2(3)*dwji2(2) + deidwi2*d2wji2(3,2))*cosij
      derv2(jx,iz) = derv2(jx,iz) - (d2eidwi22*dwji2(1)*dwji2(3) + deidwi2*d2wji2(1,3))*cosij
      derv2(jy,iz) = derv2(jy,iz) - (d2eidwi22*dwji2(2)*dwji2(3) + deidwi2*d2wji2(2,3))*cosij
      derv2(jz,iz) = derv2(jz,iz) - (d2eidwi22*dwji2(3)*dwji2(3) + deidwi2*d2wji2(3,3))*cosij
!
      dervi(jx,ix) = dervi(jx,ix) - (d2eidwi22*dwji2(1)*dwji2(1) + deidwi2*d2wji2(1,1))*sinij
      dervi(jy,ix) = dervi(jy,ix) - (d2eidwi22*dwji2(2)*dwji2(1) + deidwi2*d2wji2(2,1))*sinij
      dervi(jz,ix) = dervi(jz,ix) - (d2eidwi22*dwji2(3)*dwji2(1) + deidwi2*d2wji2(3,1))*sinij
      dervi(jx,iy) = dervi(jx,iy) - (d2eidwi22*dwji2(1)*dwji2(2) + deidwi2*d2wji2(1,2))*sinij
      dervi(jy,iy) = dervi(jy,iy) - (d2eidwi22*dwji2(2)*dwji2(2) + deidwi2*d2wji2(2,2))*sinij
      dervi(jz,iy) = dervi(jz,iy) - (d2eidwi22*dwji2(3)*dwji2(2) + deidwi2*d2wji2(3,2))*sinij
      dervi(jx,iz) = dervi(jx,iz) - (d2eidwi22*dwji2(1)*dwji2(3) + deidwi2*d2wji2(1,3))*sinij
      dervi(jy,iz) = dervi(jy,iz) - (d2eidwi22*dwji2(2)*dwji2(3) + deidwi2*d2wji2(2,3))*sinij
      dervi(jz,iz) = dervi(jz,iz) - (d2eidwi22*dwji2(3)*dwji2(3) + deidwi2*d2wji2(3,3))*sinij
!
!  Loop over second neighbour of i
!
      do nj = 1,ni-1
        k = nvbnbrno(nj,i)
        kx = 3*(k-1) + 1
        ky = kx + 1
        kz = kx + 2
!
!  Set up i-k quantities
!
        rik = rvbnbr(nj,i)
        rrik = 1.0_dp/rik
        xki = xvbnbr(nj,i)
        yki = yvbnbr(nj,i)
        zki = zvbnbr(nj,i)
!
        if (i.eq.k) then
          oneik = 1.0_dp
        else
          oneik = 0.0_dp
        endif
        if (j.eq.k) then
          onejk = 1.0_dp
        else
          onejk = 0.0_dp
        endif
!
        cosik = xkv*xki + ykv*yki + zkv*zki
        sinik = sin(cosik)
        cosik = cos(cosik) - oneik
        cosjk = xkv*(xki-xji) + ykv*(yki-yji) + zkv*(zki-zji)
        sinjk = sin(cosjk)
        cosjk = cos(cosjk) - onejk
!
        d2eidr2 = d2eidbvsum2*d1i(ni)*d1i(nj)
!
!  Derivatives that would have been added in add_drv2_2
!
!  I-K
!
        derv2(kx,ix) = derv2(kx,ix) - d2eidr2*xki*xji*cosik
        derv2(ky,ix) = derv2(ky,ix) - d2eidr2*yki*xji*cosik
        derv2(kz,ix) = derv2(kz,ix) - d2eidr2*zki*xji*cosik
        derv2(kx,iy) = derv2(kx,iy) - d2eidr2*xki*yji*cosik
        derv2(ky,iy) = derv2(ky,iy) - d2eidr2*yki*yji*cosik
        derv2(kz,iy) = derv2(kz,iy) - d2eidr2*zki*yji*cosik
        derv2(kx,iz) = derv2(kx,iz) - d2eidr2*xki*zji*cosik
        derv2(ky,iz) = derv2(ky,iz) - d2eidr2*yki*zji*cosik
        derv2(kz,iz) = derv2(kz,iz) - d2eidr2*zki*zji*cosik
!
        dervi(kx,ix) = dervi(kx,ix) - d2eidr2*xki*xji*sinik
        dervi(ky,ix) = dervi(ky,ix) - d2eidr2*yki*xji*sinik
        dervi(kz,ix) = dervi(kz,ix) - d2eidr2*zki*xji*sinik
        dervi(kx,iy) = dervi(kx,iy) - d2eidr2*xki*yji*sinik
        dervi(ky,iy) = dervi(ky,iy) - d2eidr2*yki*yji*sinik
        dervi(kz,iy) = dervi(kz,iy) - d2eidr2*zki*yji*sinik
        dervi(kx,iz) = dervi(kx,iz) - d2eidr2*xki*zji*sinik
        dervi(ky,iz) = dervi(ky,iz) - d2eidr2*yki*zji*sinik
        dervi(kz,iz) = dervi(kz,iz) - d2eidr2*zki*zji*sinik
!
!  I-J
!
        derv2(jx,ix) = derv2(jx,ix) - d2eidr2*xki*xji*cosij
        derv2(jy,ix) = derv2(jy,ix) - d2eidr2*xki*yji*cosij
        derv2(jz,ix) = derv2(jz,ix) - d2eidr2*xki*zji*cosij
        derv2(jx,iy) = derv2(jx,iy) - d2eidr2*yki*xji*cosij
        derv2(jy,iy) = derv2(jy,iy) - d2eidr2*yki*yji*cosij
        derv2(jz,iy) = derv2(jz,iy) - d2eidr2*yki*zji*cosij
        derv2(jx,iz) = derv2(jx,iz) - d2eidr2*zki*xji*cosij
        derv2(jy,iz) = derv2(jy,iz) - d2eidr2*zki*yji*cosij
        derv2(jz,iz) = derv2(jz,iz) - d2eidr2*zki*zji*cosij
!
        dervi(jx,ix) = dervi(jx,ix) - d2eidr2*xki*xji*sinij
        dervi(jy,ix) = dervi(jy,ix) - d2eidr2*xki*yji*sinij
        dervi(jz,ix) = dervi(jz,ix) - d2eidr2*xki*zji*sinij
        dervi(jx,iy) = dervi(jx,iy) - d2eidr2*yki*xji*sinij
        dervi(jy,iy) = dervi(jy,iy) - d2eidr2*yki*yji*sinij
        dervi(jz,iy) = dervi(jz,iy) - d2eidr2*yki*zji*sinij
        dervi(jx,iz) = dervi(jx,iz) - d2eidr2*zki*xji*sinij
        dervi(jy,iz) = dervi(jy,iz) - d2eidr2*zki*yji*sinij
        dervi(jz,iz) = dervi(jz,iz) - d2eidr2*zki*zji*sinij
!
        rrik3 = rrik**3
        d1wki(1,1) = d1i(nj)*xki*xki*rrik - bv(nj)*(rrik3*xki*xki - rrik)
        d1wki(2,1) = d1i(nj)*yki*xki*rrik - bv(nj)*rrik3*yki*xki
        d1wki(3,1) = d1i(nj)*zki*xki*rrik - bv(nj)*rrik3*zki*xki
!
        d1wki(1,2) = d1i(nj)*xki*yki*rrik - bv(nj)*rrik3*xki*yki
        d1wki(2,2) = d1i(nj)*yki*yki*rrik - bv(nj)*(rrik3*yki*yki - rrik)
        d1wki(3,2) = d1i(nj)*zki*yki*rrik - bv(nj)*rrik3*zki*yki
!
        d1wki(1,3) = d1i(nj)*xki*zki*rrik - bv(nj)*rrik3*xki*zki
        d1wki(2,3) = d1i(nj)*yki*zki*rrik - bv(nj)*rrik3*yki*zki
        d1wki(3,3) = d1i(nj)*zki*zki*rrik - bv(nj)*(rrik3*zki*zki - rrik)
!
        dwki2(1) = 2.0_dp*(bvvsum(1)*d1wki(1,1) + bvvsum(2)*d1wki(1,2) + bvvsum(3)*d1wki(1,3))
        dwki2(2) = 2.0_dp*(bvvsum(1)*d1wki(2,1) + bvvsum(2)*d1wki(2,2) + bvvsum(3)*d1wki(2,3))
        dwki2(3) = 2.0_dp*(bvvsum(1)*d1wki(3,1) + bvvsum(2)*d1wki(3,2) + bvvsum(3)*d1wki(3,3))
!
        d2wkj2(1,1) = 2.0_dp*(d1wki(1,1)*d1wji(1,1) + d1wki(1,2)*d1wji(1,2) + d1wki(1,3)*d1wji(1,3))
        d2wkj2(2,1) = 2.0_dp*(d1wki(2,1)*d1wji(1,1) + d1wki(2,2)*d1wji(1,2) + d1wki(2,3)*d1wji(1,3))
        d2wkj2(3,1) = 2.0_dp*(d1wki(3,1)*d1wji(1,1) + d1wki(3,2)*d1wji(1,2) + d1wki(3,3)*d1wji(1,3))
        d2wkj2(1,2) = 2.0_dp*(d1wki(1,1)*d1wji(2,1) + d1wki(1,2)*d1wji(2,2) + d1wki(1,3)*d1wji(2,3))
        d2wkj2(2,2) = 2.0_dp*(d1wki(2,1)*d1wji(2,1) + d1wki(2,2)*d1wji(2,2) + d1wki(2,3)*d1wji(2,3))
        d2wkj2(3,2) = 2.0_dp*(d1wki(3,1)*d1wji(2,1) + d1wki(3,2)*d1wji(2,2) + d1wki(3,3)*d1wji(2,3))
        d2wkj2(1,3) = 2.0_dp*(d1wki(1,1)*d1wji(3,1) + d1wki(1,2)*d1wji(3,2) + d1wki(1,3)*d1wji(3,3))
        d2wkj2(2,3) = 2.0_dp*(d1wki(2,1)*d1wji(3,1) + d1wki(2,2)*d1wji(3,2) + d1wki(2,3)*d1wji(3,3))
        d2wkj2(3,3) = 2.0_dp*(d1wki(3,1)*d1wji(3,1) + d1wki(3,2)*d1wji(3,2) + d1wki(3,3)*d1wji(3,3))
!
        derv2(jx,ix) = derv2(jx,ix) - (deidwi2*d2wkj2(1,1) + d2eidwi22*dwki2(1)*dwji2(1))*cosij
        derv2(jy,ix) = derv2(jy,ix) - (deidwi2*d2wkj2(1,2) + d2eidwi22*dwki2(1)*dwji2(2))*cosij
        derv2(jz,ix) = derv2(jz,ix) - (deidwi2*d2wkj2(1,3) + d2eidwi22*dwki2(1)*dwji2(3))*cosij
        derv2(jx,iy) = derv2(jx,iy) - (deidwi2*d2wkj2(2,1) + d2eidwi22*dwki2(2)*dwji2(1))*cosij
        derv2(jy,iy) = derv2(jy,iy) - (deidwi2*d2wkj2(2,2) + d2eidwi22*dwki2(2)*dwji2(2))*cosij
        derv2(jz,iy) = derv2(jz,iy) - (deidwi2*d2wkj2(2,3) + d2eidwi22*dwki2(2)*dwji2(3))*cosij
        derv2(jx,iz) = derv2(jx,iz) - (deidwi2*d2wkj2(3,1) + d2eidwi22*dwki2(3)*dwji2(1))*cosij
        derv2(jy,iz) = derv2(jy,iz) - (deidwi2*d2wkj2(3,2) + d2eidwi22*dwki2(3)*dwji2(2))*cosij
        derv2(jz,iz) = derv2(jz,iz) - (deidwi2*d2wkj2(3,3) + d2eidwi22*dwki2(3)*dwji2(3))*cosij
!
        dervi(jx,ix) = dervi(jx,ix) - (deidwi2*d2wkj2(1,1) + d2eidwi22*dwki2(1)*dwji2(1))*sinij
        dervi(jy,ix) = dervi(jy,ix) - (deidwi2*d2wkj2(1,2) + d2eidwi22*dwki2(1)*dwji2(2))*sinij
        dervi(jz,ix) = dervi(jz,ix) - (deidwi2*d2wkj2(1,3) + d2eidwi22*dwki2(1)*dwji2(3))*sinij
        dervi(jx,iy) = dervi(jx,iy) - (deidwi2*d2wkj2(2,1) + d2eidwi22*dwki2(2)*dwji2(1))*sinij
        dervi(jy,iy) = dervi(jy,iy) - (deidwi2*d2wkj2(2,2) + d2eidwi22*dwki2(2)*dwji2(2))*sinij
        dervi(jz,iy) = dervi(jz,iy) - (deidwi2*d2wkj2(2,3) + d2eidwi22*dwki2(2)*dwji2(3))*sinij
        dervi(jx,iz) = dervi(jx,iz) - (deidwi2*d2wkj2(3,1) + d2eidwi22*dwki2(3)*dwji2(1))*sinij
        dervi(jy,iz) = dervi(jy,iz) - (deidwi2*d2wkj2(3,2) + d2eidwi22*dwki2(3)*dwji2(2))*sinij
        dervi(jz,iz) = dervi(jz,iz) - (deidwi2*d2wkj2(3,3) + d2eidwi22*dwki2(3)*dwji2(3))*sinij
!
        derv2(kx,ix) = derv2(kx,ix) - (deidwi2*d2wkj2(1,1) + d2eidwi22*dwji2(1)*dwki2(1))*cosik
        derv2(ky,ix) = derv2(ky,ix) - (deidwi2*d2wkj2(2,1) + d2eidwi22*dwji2(1)*dwki2(2))*cosik
        derv2(kz,ix) = derv2(kz,ix) - (deidwi2*d2wkj2(3,1) + d2eidwi22*dwji2(1)*dwki2(3))*cosik
        derv2(kx,iy) = derv2(kx,iy) - (deidwi2*d2wkj2(1,2) + d2eidwi22*dwji2(2)*dwki2(1))*cosik
        derv2(ky,iy) = derv2(ky,iy) - (deidwi2*d2wkj2(2,2) + d2eidwi22*dwji2(2)*dwki2(2))*cosik
        derv2(kz,iy) = derv2(kz,iy) - (deidwi2*d2wkj2(3,2) + d2eidwi22*dwji2(2)*dwki2(3))*cosik
        derv2(kx,iz) = derv2(kx,iz) - (deidwi2*d2wkj2(1,3) + d2eidwi22*dwji2(3)*dwki2(1))*cosik
        derv2(ky,iz) = derv2(ky,iz) - (deidwi2*d2wkj2(2,3) + d2eidwi22*dwji2(3)*dwki2(2))*cosik
        derv2(kz,iz) = derv2(kz,iz) - (deidwi2*d2wkj2(3,3) + d2eidwi22*dwji2(3)*dwki2(3))*cosik
!
        dervi(kx,ix) = dervi(kx,ix) - (deidwi2*d2wkj2(1,1) + d2eidwi22*dwji2(1)*dwki2(1))*sinik
        dervi(ky,ix) = dervi(ky,ix) - (deidwi2*d2wkj2(2,1) + d2eidwi22*dwji2(1)*dwki2(2))*sinik
        dervi(kz,ix) = dervi(kz,ix) - (deidwi2*d2wkj2(3,1) + d2eidwi22*dwji2(1)*dwki2(3))*sinik
        dervi(kx,iy) = dervi(kx,iy) - (deidwi2*d2wkj2(1,2) + d2eidwi22*dwji2(2)*dwki2(1))*sinik
        dervi(ky,iy) = dervi(ky,iy) - (deidwi2*d2wkj2(2,2) + d2eidwi22*dwji2(2)*dwki2(2))*sinik
        dervi(kz,iy) = dervi(kz,iy) - (deidwi2*d2wkj2(3,2) + d2eidwi22*dwji2(2)*dwki2(3))*sinik
        dervi(kx,iz) = dervi(kx,iz) - (deidwi2*d2wkj2(1,3) + d2eidwi22*dwji2(3)*dwki2(1))*sinik
        dervi(ky,iz) = dervi(ky,iz) - (deidwi2*d2wkj2(2,3) + d2eidwi22*dwji2(3)*dwki2(2))*sinik
        dervi(kz,iz) = dervi(kz,iz) - (deidwi2*d2wkj2(3,3) + d2eidwi22*dwji2(3)*dwki2(3))*sinik
      enddo
    enddo
  enddo
!***********************************************************************************
!  Second pass through derivatives to capture terms where j is local instead of i  *
!***********************************************************************************
  ix = - 2
  iy = - 1
  iz =   0
  iloop: do i = 1,numat
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!
!  Check whether any neighbours of i are local to this node
!
    lanylocal = .false.
    do ni = 1,nvbnbr(i)
      j = nvbnbrno(ni,i)
!
!  Set flag as to whether any j atoms are local to this node
!
      if (atom2local(j).gt.0) lanylocal = .true.
    enddo
!
    if (.not.lanylocal) cycle iloop
!
!  Set variables relating to i
!
    nati = nat(i)
    ntypi = nftype(i)
    nmi = natmol(i)
!
!  Set COM coordinates
!
    if (lrigid.and.nmi.gt.0) then
      xcomi = molxyz(1,natinmol(i),nmi)
      ycomi = molxyz(2,natinmol(i),nmi)
      zcomi = molxyz(3,natinmol(i),nmi)
    else
      xcomi = 0.0_dp
      ycomi = 0.0_dp
      zcomi = 0.0_dp
    endif
!
!  Initialise bond valence and bond valence vector
!
    bv(1:nvbnbr(i)) = 0.0_dp
    bvv(1:3,1:nvbnbr(i)) = 0.0_dp
!
!  Set total number of distances for neighbours of i
!
    nneighi2 = nvbnbr(i) + nvbnbr(i)*(nvbnbr(i) + 1)/2 
!
!  Initialise derivative storage for neighbours of i
!
    d1i(1:nvbnbr(i)) = 0.0_dp
    d2i(1:nvbnbr(i)) = 0.0_dp
!
!  Loop over neighbours of i
!
    do ni = 1,nvbnbr(i)
      j = nvbnbrno(ni,i)
!
!  Set variables relating to j
!
      natj = nat(j)
      ntypj = nftype(j)
!
!  Set up i-j quantities
!
      rij = rvbnbr(ni,i)
      xji = xvbnbr(ni,i)
      yji = yvbnbr(ni,i)
      zji = zvbnbr(ni,i)
      rrij = 1.0_dp/rij
!
!  Find valence bond potential between i and j
!
      nvbij = 0
      do nvi = 1,nvalbondatm(i)
        nv = nvalbondatmptr(nvi,i)
        if (nVBspecB2(nv).eq.maxele.or.(nVBspecB2(nv).eq.natj.and.(nVBtypeB2(nv).eq.ntypj.or.nVBtypeB2(nv).eq.0))) then
          nvbij = nvbij + 1
          nvbijptr(nvbij) = nv
        endif
      enddo
      if (nvbij.gt.0) then
!*********************************
!  Calculate valence bond terms  *
!*********************************
        do nvb = 1,nvbij
          nv = nvbijptr(nvb)
          if (rij.ge.rVBmin(nv).and.rij.le.rVBmax(nv)) then
            if (nvalbondtype(nv).eq.1) then
!
!  Exponential form
!
              vij = VBwgtB(nv)*exp((VBparB(1,nv) - rij)/VBparB(3,nv))
              d1vijdr = - vij*rrij/VBparB(3,nv)
              d2vijdr = - d1vijdr*rrij*(1.0_dp/VBparB(3,nv) + rrij)
            else
!
!  Power law form
!
              vij = VBwgtB(nv)*(VBparB(1,nv)/rij)**VBparB(2,nv)
              d1vijdr = - VBparB(2,nv)*vij*rrij*rrij
              d2vijdr = - (VBparB(2,nv) + 2.0_dp)*d1vijdr*rrij*rrij
            endif
            bv(ni) = bv(ni) + vij
            bvv(1,ni) = bvv(1,ni) + vij*xji*rrij
            bvv(2,ni) = bvv(2,ni) + vij*yji*rrij
            bvv(3,ni) = bvv(3,ni) + vij*zji*rrij
!
            d1i(ni) = d1i(ni) + d1vijdr
            d2i(ni) = d2i(ni) + d2vijdr
          endif
        enddo
      endif
    enddo
!**********************************
!  Calculate valence bond energy  *
!**********************************
    bvsum = 0.0_dp
    bvvsum(1:3) = 0.0_dp
    do ni = 1,nvbnbr(i)
      bvsum = bvsum + bv(ni)
      bvvsum(1:3) = bvvsum(1:3) + bvv(1:3,ni)
    enddo
    wi2 = bvvsum(1)**2 + bvvsum(2)**2 + bvvsum(3)**2
!***************************************
!  Derivatives of valence bond energy  *
!***************************************
!
!  Compute derivatives of energy with respect to bvsum / wi2
!
    deidbvsum   = 0.0_dp
    d2eidbvsum2 = 0.0_dp
    deidwi2     = 0.0_dp
    d2eidwi22   = 0.0_dp
    do nve = 1,nvalener
      if (nVBspecE(nve).eq.nati) then
        if (nVBtypeE(nve).eq.ntypi.or.nVBtypeE(nve).eq.0) then
          deidbvsum   = deidbvsum   + 2.0_dp*VBparE(1,nve)*(bvsum - VBparE(2,nve))
          deidwi2     = deidwi2     + 2.0_dp*VBparE(3,nve)*(wi2 - VBparE(4,nve)**2)
          d2eidbvsum2 = d2eidbvsum2 + 2.0_dp*VBparE(1,nve)
          d2eidwi22   = d2eidwi22   + 2.0_dp*VBparE(3,nve)
        endif
      endif
    enddo
!
!  Loop over neighbours of i
!
    jloop: do ni = 1,nvbnbr(i)
      j = nvbnbrno(ni,i)
      jloc = atom2local(j)
!
!  Skip if j is not local to this node
!
      if (jloc.eq.0) cycle jloop
!
      jx = 3*(jloc-1) + 1
      jy = jx + 1
      jz = jx + 2
!
!  Set up i-j quantities
!
      rij = rvbnbr(ni,i)
      rrij = 1.0_dp/rij
      xji = xvbnbr(ni,i)
      yji = yvbnbr(ni,i)
      zji = zvbnbr(ni,i)
!
!  Set COM coordinates
!
      nmj = natmol(j)
      if (lrigid.and.nmj.gt.0) then
        xcom = molxyz(1,natinmol(j),nmj) - xcomi
        ycom = molxyz(2,natinmol(j),nmj) - ycomi
        zcom = molxyz(3,natinmol(j),nmj) - zcomi
      else
        xcom = - xcomi
        ycom = - ycomi
        zcom = - zcomi
      endif
!
!  Set diagonal blocks to be difference from gamma point
!
      if (i.eq.j) then
        oneij = 1.0_dp
      else
        oneij = 0.0_dp
      endif
!
!  Compute phase factor for i-j
!
      cosij = xkv*xji + ykv*yji + zkv*zji
      sinij = sin(cosij)
      cosij = cos(cosij) - oneij
!
!  Derivatives of E_bv
!
      deijdr = deidbvsum*d1i(ni)
!
      call real1strterm(ndim,xji,yji,zji,xcom,ycom,zcom,drij2ds,d2rij2dx2,d2rij2dsdx,d2rij2ds2,.true.)
!
!  Derivatives of E_bvv
!
      rrij3 = rrij**3
      d1wji(1,1) = d1i(ni)*xji*xji*rrij - bv(ni)*(rrij3*xji*xji - rrij)
      d1wji(2,1) = d1i(ni)*yji*xji*rrij - bv(ni)*rrij3*yji*xji
      d1wji(3,1) = d1i(ni)*zji*xji*rrij - bv(ni)*rrij3*zji*xji
!
      d1wji(1,2) = d1i(ni)*xji*yji*rrij - bv(ni)*rrij3*xji*yji
      d1wji(2,2) = d1i(ni)*yji*yji*rrij - bv(ni)*(rrij3*yji*yji - rrij)
      d1wji(3,2) = d1i(ni)*zji*yji*rrij - bv(ni)*rrij3*zji*yji
!
      d1wji(1,3) = d1i(ni)*xji*zji*rrij - bv(ni)*rrij3*xji*zji
      d1wji(2,3) = d1i(ni)*yji*zji*rrij - bv(ni)*rrij3*yji*zji
      d1wji(3,3) = d1i(ni)*zji*zji*rrij - bv(ni)*(rrij3*zji*zji - rrij)
!
      dwji2(1) = 2.0_dp*(bvvsum(1)*d1wji(1,1) + bvvsum(2)*d1wji(1,2) + bvvsum(3)*d1wji(1,3))
      dwji2(2) = 2.0_dp*(bvvsum(1)*d1wji(2,1) + bvvsum(2)*d1wji(2,2) + bvvsum(3)*d1wji(2,3))
      dwji2(3) = 2.0_dp*(bvvsum(1)*d1wji(3,1) + bvvsum(2)*d1wji(3,2) + bvvsum(3)*d1wji(3,3))
!
      call cartstrterm(ndim,xji,yji,zji,xcom,ycom,zcom,dxyzijds,d2xyzijdsdx,d2xyzijds2,.true.)
!
!  Second derivatives
!
      jx = 3*(jloc-1) + 1
      jy = jx + 1
      jz = jx + 2
!
!  Second derivatives of a single Vij
!
      d2eidr2 = deidbvsum*d2i(ni) + d2eidbvsum2*d1i(ni)**2
!
      derv2(ix,jx) = derv2(ix,jx) - d2eidr2*d2rij2dx2(1)*cosij
      derv2(iy,jx) = derv2(iy,jx) - d2eidr2*d2rij2dx2(6)*cosij
      derv2(iz,jx) = derv2(iz,jx) - d2eidr2*d2rij2dx2(5)*cosij
      derv2(ix,jy) = derv2(ix,jy) - d2eidr2*d2rij2dx2(6)*cosij
      derv2(iy,jy) = derv2(iy,jy) - d2eidr2*d2rij2dx2(2)*cosij
      derv2(iz,jy) = derv2(iz,jy) - d2eidr2*d2rij2dx2(4)*cosij
      derv2(ix,jz) = derv2(ix,jz) - d2eidr2*d2rij2dx2(5)*cosij
      derv2(iy,jz) = derv2(iy,jz) - d2eidr2*d2rij2dx2(4)*cosij
      derv2(iz,jz) = derv2(iz,jz) - d2eidr2*d2rij2dx2(3)*cosij
      derv2(ix,jx) = derv2(ix,jx) - deijdr*cosij
      derv2(iy,jy) = derv2(iy,jy) - deijdr*cosij
      derv2(iz,jz) = derv2(iz,jz) - deijdr*cosij
!
      dervi(ix,jx) = dervi(ix,jx) + d2eidr2*d2rij2dx2(1)*sinij
      dervi(iy,jx) = dervi(iy,jx) + d2eidr2*d2rij2dx2(6)*sinij
      dervi(iz,jx) = dervi(iz,jx) + d2eidr2*d2rij2dx2(5)*sinij
      dervi(ix,jy) = dervi(ix,jy) + d2eidr2*d2rij2dx2(6)*sinij
      dervi(iy,jy) = dervi(iy,jy) + d2eidr2*d2rij2dx2(2)*sinij
      dervi(iz,jy) = dervi(iz,jy) + d2eidr2*d2rij2dx2(4)*sinij
      dervi(ix,jz) = dervi(ix,jz) + d2eidr2*d2rij2dx2(5)*sinij
      dervi(iy,jz) = dervi(iy,jz) + d2eidr2*d2rij2dx2(4)*sinij
      dervi(iz,jz) = dervi(iz,jz) + d2eidr2*d2rij2dx2(3)*sinij
      dervi(ix,jx) = dervi(ix,jx) + deijdr*sinij
      dervi(iy,jy) = dervi(iy,jy) + deijdr*sinij
      dervi(iz,jz) = dervi(iz,jz) + deijdr*sinij
!
!  Second derivatives of bvv
!
      rrij5 = rrij**5
      d2wji(1,1,1) = xji*xji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     3.0_dp*xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,1,1) = yji*xji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,1,1) = zji*xji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,2,1) = xji*yji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,2,1) = yji*yji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,2,1) = zji*yji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(1,3,1) = xji*zji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,3,1) = yji*zji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(3,3,1) = zji*zji*xji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,1,2) = xji*xji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,1,2) = yji*xji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,1,2) = zji*xji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(1,2,2) = xji*yji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,2,2) = yji*yji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     3.0_dp*yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,2,2) = zji*yji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,3,2) = xji*zji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(2,3,2) = yji*zji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,3,2) = zji*zji*yji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,1,3) = xji*xji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,1,3) = yji*xji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(3,1,3) = zji*xji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,2,3) = xji*yji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5)
      d2wji(2,2,3) = yji*yji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     zji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,2,3) = zji*yji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(1,3,3) = xji*zji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     xji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(2,3,3) = yji*zji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     yji*(d1i(ni)*rrij - bv(ni)*rrij3)
      d2wji(3,3,3) = zji*zji*zji*(d2i(ni)*rrij - 2.0_dp*d1i(ni)*rrij3 + 3.0_dp*bv(ni)*rrij5) + &
                     3.0_dp*zji*(d1i(ni)*rrij - bv(ni)*rrij3)
!
      d2wji2(1,1) = 2.0_dp*(bvvsum(1)*d2wji(1,1,1) + bvvsum(2)*d2wji(1,1,2) + bvvsum(3)*d2wji(1,1,3) + &
                            d1wji(1,1)*d1wji(1,1) + d1wji(1,2)*d1wji(1,2) + d1wji(1,3)*d1wji(1,3))
      d2wji2(2,1) = 2.0_dp*(bvvsum(1)*d2wji(2,1,1) + bvvsum(2)*d2wji(2,1,2) + bvvsum(3)*d2wji(2,1,3) + &
                            d1wji(2,1)*d1wji(1,1) + d1wji(2,2)*d1wji(1,2) + d1wji(2,3)*d1wji(1,3))
      d2wji2(3,1) = 2.0_dp*(bvvsum(1)*d2wji(3,1,1) + bvvsum(2)*d2wji(3,1,2) + bvvsum(3)*d2wji(3,1,3) + &
                            d1wji(3,1)*d1wji(1,1) + d1wji(3,2)*d1wji(1,2) + d1wji(3,3)*d1wji(1,3))
      d2wji2(1,2) = 2.0_dp*(bvvsum(1)*d2wji(1,2,1) + bvvsum(2)*d2wji(1,2,2) + bvvsum(3)*d2wji(1,2,3) + &
                            d1wji(1,1)*d1wji(2,1) + d1wji(1,2)*d1wji(2,2) + d1wji(1,3)*d1wji(2,3))
      d2wji2(2,2) = 2.0_dp*(bvvsum(1)*d2wji(2,2,1) + bvvsum(2)*d2wji(2,2,2) + bvvsum(3)*d2wji(2,2,3) + &
                            d1wji(2,1)*d1wji(2,1) + d1wji(2,2)*d1wji(2,2) + d1wji(2,3)*d1wji(2,3))
      d2wji2(3,2) = 2.0_dp*(bvvsum(1)*d2wji(3,2,1) + bvvsum(2)*d2wji(3,2,2) + bvvsum(3)*d2wji(3,2,3) + &
                            d1wji(3,1)*d1wji(2,1) + d1wji(3,2)*d1wji(2,2) + d1wji(3,3)*d1wji(2,3))
      d2wji2(1,3) = 2.0_dp*(bvvsum(1)*d2wji(1,3,1) + bvvsum(2)*d2wji(1,3,2) + bvvsum(3)*d2wji(1,3,3) + &
                            d1wji(1,1)*d1wji(3,1) + d1wji(1,2)*d1wji(3,2) + d1wji(1,3)*d1wji(3,3))
      d2wji2(2,3) = 2.0_dp*(bvvsum(1)*d2wji(2,3,1) + bvvsum(2)*d2wji(2,3,2) + bvvsum(3)*d2wji(2,3,3) + &
                            d1wji(2,1)*d1wji(3,1) + d1wji(2,2)*d1wji(3,2) + d1wji(2,3)*d1wji(3,3))
      d2wji2(3,3) = 2.0_dp*(bvvsum(1)*d2wji(3,3,1) + bvvsum(2)*d2wji(3,3,2) + bvvsum(3)*d2wji(3,3,3) + &
                            d1wji(3,1)*d1wji(3,1) + d1wji(3,2)*d1wji(3,2) + d1wji(3,3)*d1wji(3,3))
!
      derv2(ix,jx) = derv2(ix,jx) - (d2eidwi22*dwji2(1)*dwji2(1) + deidwi2*d2wji2(1,1))*cosij
      derv2(iy,jx) = derv2(iy,jx) - (d2eidwi22*dwji2(2)*dwji2(1) + deidwi2*d2wji2(2,1))*cosij
      derv2(iz,jx) = derv2(iz,jx) - (d2eidwi22*dwji2(3)*dwji2(1) + deidwi2*d2wji2(3,1))*cosij
      derv2(ix,jy) = derv2(ix,jy) - (d2eidwi22*dwji2(1)*dwji2(2) + deidwi2*d2wji2(1,2))*cosij
      derv2(iy,jy) = derv2(iy,jy) - (d2eidwi22*dwji2(2)*dwji2(2) + deidwi2*d2wji2(2,2))*cosij
      derv2(iz,jy) = derv2(iz,jy) - (d2eidwi22*dwji2(3)*dwji2(2) + deidwi2*d2wji2(3,2))*cosij
      derv2(ix,jz) = derv2(ix,jz) - (d2eidwi22*dwji2(1)*dwji2(3) + deidwi2*d2wji2(1,3))*cosij
      derv2(iy,jz) = derv2(iy,jz) - (d2eidwi22*dwji2(2)*dwji2(3) + deidwi2*d2wji2(2,3))*cosij
      derv2(iz,jz) = derv2(iz,jz) - (d2eidwi22*dwji2(3)*dwji2(3) + deidwi2*d2wji2(3,3))*cosij
!
      dervi(ix,jx) = dervi(ix,jx) + (d2eidwi22*dwji2(1)*dwji2(1) + deidwi2*d2wji2(1,1))*sinij
      dervi(iy,jx) = dervi(iy,jx) + (d2eidwi22*dwji2(2)*dwji2(1) + deidwi2*d2wji2(2,1))*sinij
      dervi(iz,jx) = dervi(iz,jx) + (d2eidwi22*dwji2(3)*dwji2(1) + deidwi2*d2wji2(3,1))*sinij
      dervi(ix,jy) = dervi(ix,jy) + (d2eidwi22*dwji2(1)*dwji2(2) + deidwi2*d2wji2(1,2))*sinij
      dervi(iy,jy) = dervi(iy,jy) + (d2eidwi22*dwji2(2)*dwji2(2) + deidwi2*d2wji2(2,2))*sinij
      dervi(iz,jy) = dervi(iz,jy) + (d2eidwi22*dwji2(3)*dwji2(2) + deidwi2*d2wji2(3,2))*sinij
      dervi(ix,jz) = dervi(ix,jz) + (d2eidwi22*dwji2(1)*dwji2(3) + deidwi2*d2wji2(1,3))*sinij
      dervi(iy,jz) = dervi(iy,jz) + (d2eidwi22*dwji2(2)*dwji2(3) + deidwi2*d2wji2(2,3))*sinij
      dervi(iz,jz) = dervi(iz,jz) + (d2eidwi22*dwji2(3)*dwji2(3) + deidwi2*d2wji2(3,3))*sinij
!
!  Loop over second neighbour of i
!
      do nj = 1,nvbnbr(i)
!
!  Exclude self term
!
        if (ni.eq.nj) cycle
!
        k = nvbnbrno(nj,i)
        kx = 3*(k-1) + 1
        ky = kx + 1
        kz = kx + 2
!
!  Set up i-k quantities
!
        rik = rvbnbr(nj,i)
        rrik = 1.0_dp/rik
        xki = xvbnbr(nj,i)
        yki = yvbnbr(nj,i)
        zki = zvbnbr(nj,i)
!
        if (i.eq.k) then
          oneik = 1.0_dp
        else
          oneik = 0.0_dp
        endif
        if (j.eq.k) then
          onejk = 1.0_dp
        else
          onejk = 0.0_dp
        endif
!
        cosik = xkv*xki + ykv*yki + zkv*zki
        sinik = sin(cosik)
        cosik = cos(cosik) - oneik
        cosjk = xkv*(xki-xji) + ykv*(yki-yji) + zkv*(zki-zji)
        sinjk = sin(cosjk)
        cosjk = cos(cosjk) - onejk
!
        d2eidr2 = d2eidbvsum2*d1i(ni)*d1i(nj)
!
!  Derivatives that would have been added in add_drv2_2
!
        derv2(ix,jx) = derv2(ix,jx) - d2eidr2*xki*xji*cosij
        derv2(iy,jx) = derv2(iy,jx) - d2eidr2*yki*xji*cosij
        derv2(iz,jx) = derv2(iz,jx) - d2eidr2*zki*xji*cosij
        derv2(ix,jy) = derv2(ix,jy) - d2eidr2*xki*yji*cosij
        derv2(iy,jy) = derv2(iy,jy) - d2eidr2*yki*yji*cosij
        derv2(iz,jy) = derv2(iz,jy) - d2eidr2*zki*yji*cosij
        derv2(ix,jz) = derv2(ix,jz) - d2eidr2*xki*zji*cosij
        derv2(iy,jz) = derv2(iy,jz) - d2eidr2*yki*zji*cosij
        derv2(iz,jz) = derv2(iz,jz) - d2eidr2*zki*zji*cosij
!
        dervi(ix,jx) = dervi(ix,jx) + d2eidr2*xki*xji*sinij
        dervi(iy,jx) = dervi(iy,jx) + d2eidr2*yki*xji*sinij
        dervi(iz,jx) = dervi(iz,jx) + d2eidr2*zki*xji*sinij
        dervi(ix,jy) = dervi(ix,jy) + d2eidr2*xki*yji*sinij
        dervi(iy,jy) = dervi(iy,jy) + d2eidr2*yki*yji*sinij
        dervi(iz,jy) = dervi(iz,jy) + d2eidr2*zki*yji*sinij
        dervi(ix,jz) = dervi(ix,jz) + d2eidr2*xki*zji*sinij
        dervi(iy,jz) = dervi(iy,jz) + d2eidr2*yki*zji*sinij
        dervi(iz,jz) = dervi(iz,jz) + d2eidr2*zki*zji*sinij
!
        derv2(kx,jx) = derv2(kx,jx) + d2eidr2*xji*xki*cosjk
        derv2(ky,jx) = derv2(ky,jx) + d2eidr2*xji*yki*cosjk
        derv2(kz,jx) = derv2(kz,jx) + d2eidr2*xji*zki*cosjk
        derv2(kx,jy) = derv2(kx,jy) + d2eidr2*yji*xki*cosjk
        derv2(ky,jy) = derv2(ky,jy) + d2eidr2*yji*yki*cosjk
        derv2(kz,jy) = derv2(kz,jy) + d2eidr2*yji*zki*cosjk
        derv2(kx,jz) = derv2(kx,jz) + d2eidr2*zji*xki*cosjk
        derv2(ky,jz) = derv2(ky,jz) + d2eidr2*zji*yki*cosjk
        derv2(kz,jz) = derv2(kz,jz) + d2eidr2*zji*zki*cosjk
!
        dervi(kx,jx) = dervi(kx,jx) + d2eidr2*xji*xki*sinjk
        dervi(ky,jx) = dervi(ky,jx) + d2eidr2*xji*yki*sinjk
        dervi(kz,jx) = dervi(kz,jx) + d2eidr2*xji*zki*sinjk
        dervi(kx,jy) = dervi(kx,jy) + d2eidr2*yji*xki*sinjk
        dervi(ky,jy) = dervi(ky,jy) + d2eidr2*yji*yki*sinjk
        dervi(kz,jy) = dervi(kz,jy) + d2eidr2*yji*zki*sinjk
        dervi(kx,jz) = dervi(kx,jz) + d2eidr2*zji*xki*sinjk
        dervi(ky,jz) = dervi(ky,jz) + d2eidr2*zji*yki*sinjk
        dervi(kz,jz) = dervi(kz,jz) + d2eidr2*zji*zki*sinjk
!
        rrik3 = rrik**3
        d1wki(1,1) = d1i(nj)*xki*xki*rrik - bv(nj)*(rrik3*xki*xki - rrik)
        d1wki(2,1) = d1i(nj)*yki*xki*rrik - bv(nj)*rrik3*yki*xki
        d1wki(3,1) = d1i(nj)*zki*xki*rrik - bv(nj)*rrik3*zki*xki
!
        d1wki(1,2) = d1i(nj)*xki*yki*rrik - bv(nj)*rrik3*xki*yki
        d1wki(2,2) = d1i(nj)*yki*yki*rrik - bv(nj)*(rrik3*yki*yki - rrik)
        d1wki(3,2) = d1i(nj)*zki*yki*rrik - bv(nj)*rrik3*zki*yki
!
        d1wki(1,3) = d1i(nj)*xki*zki*rrik - bv(nj)*rrik3*xki*zki
        d1wki(2,3) = d1i(nj)*yki*zki*rrik - bv(nj)*rrik3*yki*zki
        d1wki(3,3) = d1i(nj)*zki*zki*rrik - bv(nj)*(rrik3*zki*zki - rrik)
!
        dwki2(1) = 2.0_dp*(bvvsum(1)*d1wki(1,1) + bvvsum(2)*d1wki(1,2) + bvvsum(3)*d1wki(1,3))
        dwki2(2) = 2.0_dp*(bvvsum(1)*d1wki(2,1) + bvvsum(2)*d1wki(2,2) + bvvsum(3)*d1wki(2,3))
        dwki2(3) = 2.0_dp*(bvvsum(1)*d1wki(3,1) + bvvsum(2)*d1wki(3,2) + bvvsum(3)*d1wki(3,3))
!
        d2wkj2(1,1) = 2.0_dp*(d1wki(1,1)*d1wji(1,1) + d1wki(1,2)*d1wji(1,2) + d1wki(1,3)*d1wji(1,3))
        d2wkj2(2,1) = 2.0_dp*(d1wki(2,1)*d1wji(1,1) + d1wki(2,2)*d1wji(1,2) + d1wki(2,3)*d1wji(1,3))
        d2wkj2(3,1) = 2.0_dp*(d1wki(3,1)*d1wji(1,1) + d1wki(3,2)*d1wji(1,2) + d1wki(3,3)*d1wji(1,3))
        d2wkj2(1,2) = 2.0_dp*(d1wki(1,1)*d1wji(2,1) + d1wki(1,2)*d1wji(2,2) + d1wki(1,3)*d1wji(2,3))
        d2wkj2(2,2) = 2.0_dp*(d1wki(2,1)*d1wji(2,1) + d1wki(2,2)*d1wji(2,2) + d1wki(2,3)*d1wji(2,3))
        d2wkj2(3,2) = 2.0_dp*(d1wki(3,1)*d1wji(2,1) + d1wki(3,2)*d1wji(2,2) + d1wki(3,3)*d1wji(2,3))
        d2wkj2(1,3) = 2.0_dp*(d1wki(1,1)*d1wji(3,1) + d1wki(1,2)*d1wji(3,2) + d1wki(1,3)*d1wji(3,3))
        d2wkj2(2,3) = 2.0_dp*(d1wki(2,1)*d1wji(3,1) + d1wki(2,2)*d1wji(3,2) + d1wki(2,3)*d1wji(3,3))
        d2wkj2(3,3) = 2.0_dp*(d1wki(3,1)*d1wji(3,1) + d1wki(3,2)*d1wji(3,2) + d1wki(3,3)*d1wji(3,3))
!
        derv2(ix,jx) = derv2(ix,jx) - (deidwi2*d2wkj2(1,1) + d2eidwi22*dwki2(1)*dwji2(1))*cosij
        derv2(iy,jx) = derv2(iy,jx) - (deidwi2*d2wkj2(2,1) + d2eidwi22*dwki2(2)*dwji2(1))*cosij
        derv2(iz,jx) = derv2(iz,jx) - (deidwi2*d2wkj2(3,1) + d2eidwi22*dwki2(3)*dwji2(1))*cosij
        derv2(ix,jy) = derv2(ix,jy) - (deidwi2*d2wkj2(1,2) + d2eidwi22*dwki2(1)*dwji2(2))*cosij
        derv2(iy,jy) = derv2(iy,jy) - (deidwi2*d2wkj2(2,2) + d2eidwi22*dwki2(2)*dwji2(2))*cosij
        derv2(iz,jy) = derv2(iz,jy) - (deidwi2*d2wkj2(3,2) + d2eidwi22*dwki2(3)*dwji2(2))*cosij
        derv2(ix,jz) = derv2(ix,jz) - (deidwi2*d2wkj2(1,3) + d2eidwi22*dwki2(1)*dwji2(3))*cosij
        derv2(iy,jz) = derv2(iy,jz) - (deidwi2*d2wkj2(2,3) + d2eidwi22*dwki2(2)*dwji2(3))*cosij
        derv2(iz,jz) = derv2(iz,jz) - (deidwi2*d2wkj2(3,3) + d2eidwi22*dwki2(3)*dwji2(3))*cosij
!
        dervi(ix,jx) = dervi(ix,jx) + (deidwi2*d2wkj2(1,1) + d2eidwi22*dwki2(1)*dwji2(1))*sinij
        dervi(iy,jx) = dervi(iy,jx) + (deidwi2*d2wkj2(2,1) + d2eidwi22*dwki2(2)*dwji2(1))*sinij
        dervi(iz,jx) = dervi(iz,jx) + (deidwi2*d2wkj2(3,1) + d2eidwi22*dwki2(3)*dwji2(1))*sinij
        dervi(ix,jy) = dervi(ix,jy) + (deidwi2*d2wkj2(1,2) + d2eidwi22*dwki2(1)*dwji2(2))*sinij
        dervi(iy,jy) = dervi(iy,jy) + (deidwi2*d2wkj2(2,2) + d2eidwi22*dwki2(2)*dwji2(2))*sinij
        dervi(iz,jy) = dervi(iz,jy) + (deidwi2*d2wkj2(3,2) + d2eidwi22*dwki2(3)*dwji2(2))*sinij
        dervi(ix,jz) = dervi(ix,jz) + (deidwi2*d2wkj2(1,3) + d2eidwi22*dwki2(1)*dwji2(3))*sinij
        dervi(iy,jz) = dervi(iy,jz) + (deidwi2*d2wkj2(2,3) + d2eidwi22*dwki2(2)*dwji2(3))*sinij
        dervi(iz,jz) = dervi(iz,jz) + (deidwi2*d2wkj2(3,3) + d2eidwi22*dwki2(3)*dwji2(3))*sinij
!
        derv2(kx,jx) = derv2(kx,jx) + (deidwi2*d2wkj2(1,1) + d2eidwi22*dwji2(1)*dwki2(1))*cosjk
        derv2(ky,jx) = derv2(ky,jx) + (deidwi2*d2wkj2(2,1) + d2eidwi22*dwji2(1)*dwki2(2))*cosjk
        derv2(kz,jx) = derv2(kz,jx) + (deidwi2*d2wkj2(3,1) + d2eidwi22*dwji2(1)*dwki2(3))*cosjk
        derv2(kx,jy) = derv2(kx,jy) + (deidwi2*d2wkj2(1,2) + d2eidwi22*dwji2(2)*dwki2(1))*cosjk
        derv2(ky,jy) = derv2(ky,jy) + (deidwi2*d2wkj2(2,2) + d2eidwi22*dwji2(2)*dwki2(2))*cosjk
        derv2(kz,jy) = derv2(kz,jy) + (deidwi2*d2wkj2(3,2) + d2eidwi22*dwji2(2)*dwki2(3))*cosjk
        derv2(kx,jz) = derv2(kx,jz) + (deidwi2*d2wkj2(1,3) + d2eidwi22*dwji2(3)*dwki2(1))*cosjk
        derv2(ky,jz) = derv2(ky,jz) + (deidwi2*d2wkj2(2,3) + d2eidwi22*dwji2(3)*dwki2(2))*cosjk
        derv2(kz,jz) = derv2(kz,jz) + (deidwi2*d2wkj2(3,3) + d2eidwi22*dwji2(3)*dwki2(3))*cosjk
!
        dervi(kx,jx) = dervi(kx,jx) + (deidwi2*d2wkj2(1,1) + d2eidwi22*dwji2(1)*dwki2(1))*sinjk
        dervi(ky,jx) = dervi(ky,jx) + (deidwi2*d2wkj2(2,1) + d2eidwi22*dwji2(1)*dwki2(2))*sinjk
        dervi(kz,jx) = dervi(kz,jx) + (deidwi2*d2wkj2(3,1) + d2eidwi22*dwji2(1)*dwki2(3))*sinjk
        dervi(kx,jy) = dervi(kx,jy) + (deidwi2*d2wkj2(1,2) + d2eidwi22*dwji2(2)*dwki2(1))*sinjk
        dervi(ky,jy) = dervi(ky,jy) + (deidwi2*d2wkj2(2,2) + d2eidwi22*dwji2(2)*dwki2(2))*sinjk
        dervi(kz,jy) = dervi(kz,jy) + (deidwi2*d2wkj2(3,2) + d2eidwi22*dwji2(2)*dwki2(3))*sinjk
        dervi(kx,jz) = dervi(kx,jz) + (deidwi2*d2wkj2(1,3) + d2eidwi22*dwji2(3)*dwki2(1))*sinjk
        dervi(ky,jz) = dervi(ky,jz) + (deidwi2*d2wkj2(2,3) + d2eidwi22*dwji2(3)*dwki2(2))*sinjk
        dervi(kz,jz) = dervi(kz,jz) + (deidwi2*d2wkj2(3,3) + d2eidwi22*dwji2(3)*dwki2(3))*sinjk
      enddo
    enddo jloop
  enddo iloop
!
!  Free local memory
!
  deallocate(d2i,stat=status)
  if (status/=0) call deallocate_error('valencebondpd','d2i')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('valencebondpd','d1i')
  deallocate(bvv,stat=status)
  if (status/=0) call deallocate_error('valencebondpd','vbv')
  deallocate(bv,stat=status)
  if (status/=0) call deallocate_error('valencebondpd','vb')
  deallocate(nvbijptr,stat=status)
  if (status/=0) call deallocate_error('valencebondpd','nvbijptr')
!
  t2 = g_cpu_time()
  tval = tval + t2 - t1
#ifdef TRACE
  call trace_out('valencebondpd')
#endif
!
  return
  end
