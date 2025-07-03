  subroutine valencebondd(eval,lgrad1,lgrad2)
!
!  Calculates the energy and derivatives for the Valance Bond potentials.
!  Distributed memory version.
!
!  NB: Does not currently support lfreeze
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!
!  On exit :
!
!  eval            = the value of the energy contribution
!
!  12/21 Created from valencebond
!   2/22 Molecule centre of mass handling added
!   7/23 nregionnocfg added
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
  use datatypes
  use bondvalence
  use configurations, only : nregions, lsliceatom, nregiontype
  use control,        only : keyword, lseok, lrigid
  use current
  use derivatives
  use element,        only : maxele
  use energies,       only : eattach, esregion2
  use energies,       only : eregion2region, siteenergy
  use iochannels
  use m_strain,       only : real1strterm, cartstrterm
  use m_vb_nbr
  use molecule
  use parallel
  use symmetry,       only : lstr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),    intent(out)                         :: eval
  logical,     intent(in)                          :: lgrad1
  logical,     intent(in)                          :: lgrad2
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
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nmi
  integer(i4)                                      :: nmj
  integer(i4)                                      :: nmk
  integer(i4)                                      :: nn
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregiontypi
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
  logical                                          :: lreg2ok
  logical                                          :: lslicei
  real(dp)                                         :: bvsum
  real(dp)                                         :: bvvsum(3)
  real(dp)                                         :: g_cpu_time
  real(dp),    dimension(:),     allocatable, save :: d1i
  real(dp),    dimension(:),     allocatable, save :: d2i
  real(dp)                                         :: d1wji(3,3)
  real(dp)                                         :: d1wki(3,3)
  real(dp)                                         :: d1wjis(6,3)
  real(dp)                                         :: d1wkis(6,3)
  real(dp)                                         :: d2wji(3,3,3)
  real(dp)                                         :: dwji2(3)
  real(dp)                                         :: dwki2(3)
  real(dp)                                         :: d2wji2(3,3)
  real(dp)                                         :: d2wji2sx(3)
  real(dp)                                         :: d2wjis2(3)
  real(dp)                                         :: d2wjisx(3,3)
  real(dp)                                         :: d2wki2sx(3)
  real(dp)                                         :: d2wkj2(3,3)
  real(dp)                                         :: dwji2s(6)
  real(dp)                                         :: dwki2s(6)
  real(dp)                                         :: d2wji2s2
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
  real(dp)                                         :: drik2ds(6)
  real(dp)                                         :: d2rik2dx2(6)
  real(dp)                                         :: d2rik2ds2(6,6)
  real(dp)                                         :: d2rik2dsdx(6,3)
  real(dp)                                         :: dxyzijds(6,3)
  real(dp)                                         :: d2xyzijdsdx(6,3,3)
  real(dp)                                         :: d2xyzijds2(6,6,3)
  real(dp)                                         :: dxyzikds(6,3)
  real(dp)                                         :: d2xyzikdsdx(6,3,3)
  real(dp)                                         :: d2xyzikds2(6,6,3)
  real(dp)                                         :: eij
  real(dp)                                         :: eijbv
  real(dp)                                         :: eijbvv
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
  real(dp)                                         :: xcomk
  real(dp)                                         :: ycomk
  real(dp)                                         :: zcomk
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xki
  real(dp)                                         :: yki
  real(dp)                                         :: zki
  real(dp),    dimension(:),     allocatable, save :: bv
  real(dp),    dimension(:,:),   allocatable, save :: bvv
#ifdef TRACE
  call trace_in('valencebondd')
#endif
!
  t1 = g_cpu_time()
!
  lreg2ok = (lseok.and.nregions(ncf).gt.1)
!***************************************
!  Generate neighbour lists for atoms  *
!***************************************
  call vb_getnbr
!
  if (index(keyword,'debu').ne.0) then
    write(ioout,'(/,''  Number of Valence Bond neighbours for atoms :'',/)')
    write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    do i = 1,numat
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nvbnbr(i)
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do i = 1,numat
      write(ioout,'(i4,8(1x,i4))') i,(nvbnbrno(nn,i),nn=1,nvbnbr(i))
    enddo
  endif
!
!  Allocate local memory 
!
  allocate(nvbijptr(maxvalbond),stat=status)
  if (status/=0) call outofmemory('valencebondd','nvbijptr')
  allocate(bv(maxvbnbr),stat=status)
  if (status/=0) call outofmemory('valencebondd','vb')
  allocate(bvv(3,maxvbnbr),stat=status)
  if (status/=0) call outofmemory('valencebondd','vbv')
!
  if (lgrad1) then
    allocate(d1i(maxvbnbr),stat=status)
    if (status/=0) call outofmemory('valencebondd','d1i')
  else
    allocate(d1i(1),stat=status)
    if (status/=0) call outofmemory('valencebondd','d1i')
  endif
  if (lgrad2) then
    allocate(d2i(maxvbnbr),stat=status)   
    if (status/=0) call outofmemory('valencebondd','d2i')
  else
    allocate(d2i(1),stat=status)   
    if (status/=0) call outofmemory('valencebondd','d2i')
  endif
!
!  Initialise Bond Order energy
!
  eval = 0.0_dp
!
!  Opening banner for energy decomposition
!
  if (lPrintVB) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  VB  : Atom No.            Energy (eV)       Valence / vector '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
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
    nregioni = nregionno(nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    lslicei = lsliceatom(nsft + nrelf2a(i))
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
    if (lgrad1) then
      d1i(1:nvbnbr(i)) = 0.0_dp
      if (lgrad2) then
        d2i(1:nvbnbr(i)) = 0.0_dp
      endif
    endif
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
              if (lgrad1) then
                d1vijdr = - vij*rrij/VBparB(3,nv)
                if (lgrad2) then
                  d2vijdr = - d1vijdr*rrij*(1.0_dp/VBparB(3,nv) + rrij)
                endif
              endif
            else
!
!  Power law form
!
              vij = VBwgtB(nv)*(VBparB(1,nv)/rij)**VBparB(2,nv)
              if (lgrad1) then
                d1vijdr = - VBparB(2,nv)*vij*rrij*rrij
                if (lgrad2) then
                  d2vijdr = - (VBparB(2,nv) + 2.0_dp)*d1vijdr*rrij*rrij
                endif
              endif
            endif
            bv(ni) = bv(ni) + vij
            bvv(1,ni) = bvv(1,ni) + vij*xji*rrij
            bvv(2,ni) = bvv(2,ni) + vij*yji*rrij
            bvv(3,ni) = bvv(3,ni) + vij*zji*rrij
            if (lgrad1) then
              d1i(ni) = d1i(ni) + d1vijdr
              if (lgrad2) then
                d2i(ni) = d2i(ni) + d2vijdr
              endif
            endif
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
!
    eijbv = 0.0_dp
    eijbvv = 0.0_dp
    do nve = 1,nvalener
      if (nVBspecE(nve).eq.nati) then
        if (nVBtypeE(nve).eq.ntypi.or.nVBtypeE(nve).eq.0) then
          eijbv  = eijbv  + VBparE(1,nve)*(bvsum - VBparE(2,nve))**2 
          eijbvv = eijbvv + VBparE(3,nve)*(wi2 - VBparE(4,nve)**2)**2
        endif
      endif
    enddo
    eij = eijbv + eijbvv
!
    if (lPrintVB) then
      write(ioout,'(4x,i12,2x,'' BV : '',f16.10,1x,f12.8)') i,eijbv,bvsum
      write(ioout,'(16x,2x,'' BVV: '',f16.10,1x,3f12.8)') eijbvv,bvvsum(1:3)
    endif
!
!  Add to surface energy totals if appropriate
!
    if (lreg2ok.and.nregioni.gt.1) then
      esregion2 = esregion2 + eij
    else
      eval = eval + eij
    endif
    if (lslicei) eattach = eattach + eij
!
    eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eij
!
    siteenergy(i) = siteenergy(i) + eij
!***************************************
!  Derivatives of valence bond energy  *
!***************************************
    if (lgrad1) then
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
!  Derivatives of E_bv
!
        deijdr = deidbvsum*d1i(ni)
!
        xdrv(i) = xdrv(i) - deijdr*xji
        ydrv(i) = ydrv(i) - deijdr*yji
        zdrv(i) = zdrv(i) - deijdr*zji
!
        xdrv(j) = xdrv(j) + deijdr*xji
        ydrv(j) = ydrv(j) + deijdr*yji
        zdrv(j) = zdrv(j) + deijdr*zji
!
        if (lstr.or.lgrad2) then
          call real1strterm(ndim,xji,yji,zji,xcom,ycom,zcom,drij2ds,d2rij2dx2,d2rij2dsdx,d2rij2ds2,lgrad2)
        endif
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijdr*drij2ds(ks)
          enddo
        endif
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
        xdrv(i) = xdrv(i) - deidwi2*dwji2(1)
        ydrv(i) = ydrv(i) - deidwi2*dwji2(2)
        zdrv(i) = zdrv(i) - deidwi2*dwji2(3)
!
        xdrv(j) = xdrv(j) + deidwi2*dwji2(1)
        ydrv(j) = ydrv(j) + deidwi2*dwji2(2)
        zdrv(j) = zdrv(j) + deidwi2*dwji2(3)
!
        if (lstr.or.lgrad2) then
          call cartstrterm(ndim,xji,yji,zji,xcom,ycom,zcom,dxyzijds,d2xyzijdsdx,d2xyzijds2,lgrad2)
        endif
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            d1wjis(kl,1) = drij2ds(ks)*xji*rrij*(d1i(ni) - bv(ni)*rrij*rrij) + bv(ni)*rrij*dxyzijds(kl,1)
            d1wjis(kl,2) = drij2ds(ks)*yji*rrij*(d1i(ni) - bv(ni)*rrij*rrij) + bv(ni)*rrij*dxyzijds(kl,2)
            d1wjis(kl,3) = drij2ds(ks)*zji*rrij*(d1i(ni) - bv(ni)*rrij*rrij) + bv(ni)*rrij*dxyzijds(kl,3)
!
            dwji2s(kl) = 2.0_dp*(bvvsum(1)*d1wjis(kl,1) + bvvsum(2)*d1wjis(kl,2) + bvvsum(3)*d1wjis(kl,3))
!
            rstrd(kl) = rstrd(kl) + deidwi2*dwji2s(kl)
          enddo
        endif
!
!  Second derivatives
!
        if (lgrad2) then
          jx = 3*(j-1) + 1
          jy = jx + 1
          jz = jx + 2
!
!  Second derivatives of a single Vij
!
          d2eidr2 = deidbvsum*d2i(ni) + d2eidbvsum2*d1i(ni)**2
!
          derv2(jx,ix) = derv2(jx,ix) - d2eidr2*d2rij2dx2(1)
          derv2(jy,ix) = derv2(jy,ix) - d2eidr2*d2rij2dx2(6)
          derv2(jz,ix) = derv2(jz,ix) - d2eidr2*d2rij2dx2(5)
          derv2(jx,iy) = derv2(jx,iy) - d2eidr2*d2rij2dx2(6)
          derv2(jy,iy) = derv2(jy,iy) - d2eidr2*d2rij2dx2(2)
          derv2(jz,iy) = derv2(jz,iy) - d2eidr2*d2rij2dx2(4)
          derv2(jx,iz) = derv2(jx,iz) - d2eidr2*d2rij2dx2(5)
          derv2(jy,iz) = derv2(jy,iz) - d2eidr2*d2rij2dx2(4)
          derv2(jz,iz) = derv2(jz,iz) - d2eidr2*d2rij2dx2(3)
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
                derv3(ix,kl) = derv3(ix,kl) - deijdr*d2rij2dsdx(ks,1) - xji*d2eidr2*drij2ds(ks)
                derv3(iy,kl) = derv3(iy,kl) - deijdr*d2rij2dsdx(ks,2) - yji*d2eidr2*drij2ds(ks)
                derv3(iz,kl) = derv3(iz,kl) - deijdr*d2rij2dsdx(ks,3) - zji*d2eidr2*drij2ds(ks)
              enddo
            endif
!
!  Strain-strain
!
            do kk = 1,nstrains
              ks = nstrptr(kk)
              do kl = 1,nstrains
                kt = nstrptr(kl)
                sderv2(kl,kk) = sderv2(kl,kk) + d2eidr2*drij2ds(kt)*drij2ds(ks) + deijdr*d2rij2ds2(kt,ks)
              enddo
            enddo
          endif
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
          derv2(jx,ix) = derv2(jx,ix) - d2eidwi22*dwji2(1)*dwji2(1) - deidwi2*d2wji2(1,1)
          derv2(jy,ix) = derv2(jy,ix) - d2eidwi22*dwji2(2)*dwji2(1) - deidwi2*d2wji2(2,1)
          derv2(jz,ix) = derv2(jz,ix) - d2eidwi22*dwji2(3)*dwji2(1) - deidwi2*d2wji2(3,1)
          derv2(jx,iy) = derv2(jx,iy) - d2eidwi22*dwji2(1)*dwji2(2) - deidwi2*d2wji2(1,2)
          derv2(jy,iy) = derv2(jy,iy) - d2eidwi22*dwji2(2)*dwji2(2) - deidwi2*d2wji2(2,2)
          derv2(jz,iy) = derv2(jz,iy) - d2eidwi22*dwji2(3)*dwji2(2) - deidwi2*d2wji2(3,2)
          derv2(jx,iz) = derv2(jx,iz) - d2eidwi22*dwji2(1)*dwji2(3) - deidwi2*d2wji2(1,3)
          derv2(jy,iz) = derv2(jy,iz) - d2eidwi22*dwji2(2)*dwji2(3) - deidwi2*d2wji2(2,3)
          derv2(jz,iz) = derv2(jz,iz) - d2eidwi22*dwji2(3)*dwji2(3) - deidwi2*d2wji2(3,3)
!
          if (lstr) then
            do kk = 1,nstrains
              ks = nstrptr(kk)
!
              d2wjisx(1,1) = d2rij2dsdx(ks,1)*xji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,1,1) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(xji*dxyzijds(kk,1) + drij2ds(ks)) + &
                             drij2ds(ks)*xji*xji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
              d2wjisx(2,1) = d2rij2dsdx(ks,2)*xji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,1,2) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(yji*dxyzijds(kk,1)) + &
                             drij2ds(ks)*yji*xji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
              d2wjisx(3,1) = d2rij2dsdx(ks,3)*xji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,1,3) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(zji*dxyzijds(kk,1)) + &
                             drij2ds(ks)*zji*xji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
!
              d2wjisx(1,2) = d2rij2dsdx(ks,1)*yji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,2,1) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(xji*dxyzijds(kk,2)) + &
                             drij2ds(ks)*xji*yji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
              d2wjisx(2,2) = d2rij2dsdx(ks,2)*yji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,2,2) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(yji*dxyzijds(kk,2) + drij2ds(ks)) + &
                             drij2ds(ks)*yji*yji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
              d2wjisx(3,2) = d2rij2dsdx(ks,3)*yji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,2,3) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(zji*dxyzijds(kk,2)) + &
                             drij2ds(ks)*zji*yji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
!
              d2wjisx(1,3) = d2rij2dsdx(ks,1)*zji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,3,1) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(xji*dxyzijds(kk,3)) + &
                             drij2ds(ks)*xji*zji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
              d2wjisx(2,3) = d2rij2dsdx(ks,2)*zji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,3,2) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(yji*dxyzijds(kk,3)) + &
                             drij2ds(ks)*yji*zji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
              d2wjisx(3,3) = d2rij2dsdx(ks,3)*zji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,3,3) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(zji*dxyzijds(kk,3) + drij2ds(ks)) + &
                             drij2ds(ks)*zji*zji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
!
              d2wji2sx(1) = 2.0_dp*(bvvsum(1)*d2wjisx(1,1) + bvvsum(2)*d2wjisx(1,2) + bvvsum(3)*d2wjisx(1,3) + &
                                    d1wjis(kk,1)*d1wji(1,1) + d1wjis(kk,2)*d1wji(1,2) + d1wjis(kk,3)*d1wji(1,3))
              d2wji2sx(2) = 2.0_dp*(bvvsum(1)*d2wjisx(2,1) + bvvsum(2)*d2wjisx(2,2) + bvvsum(3)*d2wjisx(2,3) + &
                                    d1wjis(kk,1)*d1wji(2,1) + d1wjis(kk,2)*d1wji(2,2) + d1wjis(kk,3)*d1wji(2,3))
              d2wji2sx(3) = 2.0_dp*(bvvsum(1)*d2wjisx(3,1) + bvvsum(2)*d2wjisx(3,2) + bvvsum(3)*d2wjisx(3,3) + &
                                    d1wjis(kk,1)*d1wji(3,1) + d1wjis(kk,2)*d1wji(3,2) + d1wjis(kk,3)*d1wji(3,3))
!
              derv3(ix,kk) = derv3(ix,kk) - d2eidwi22*dwji2s(kk)*dwji2(1) - deidwi2*d2wji2sx(1)
              derv3(iy,kk) = derv3(iy,kk) - d2eidwi22*dwji2s(kk)*dwji2(2) - deidwi2*d2wji2sx(2)
              derv3(iz,kk) = derv3(iz,kk) - d2eidwi22*dwji2s(kk)*dwji2(3) - deidwi2*d2wji2sx(3)
!
              do kl = 1,nstrains
                kt = nstrptr(kl)
                d2wjis2(1) = d2rij2ds2(kt,ks)*xji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijds2(kl,kk,1) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(drij2ds(ks)*dxyzijds(kl,1) + drij2ds(kt)*dxyzijds(kk,1)) + &
                             drij2ds(ks)*drij2ds(kt)*xji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
                d2wjis2(2) = d2rij2ds2(kt,ks)*yji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijds2(kl,kk,2) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(drij2ds(ks)*dxyzijds(kl,2) + drij2ds(kt)*dxyzijds(kk,2)) + &
                             drij2ds(ks)*drij2ds(kt)*yji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
                d2wjis2(3) = d2rij2ds2(kt,ks)*zji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijds2(kl,kk,3) + &
                             (rrij*d1i(ni) - bv(ni)*rrij3)*(drij2ds(ks)*dxyzijds(kl,3) + drij2ds(kt)*dxyzijds(kk,3)) + &
                             drij2ds(ks)*drij2ds(kt)*zji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
!
                d2wji2s2 = 2.0_dp*(bvvsum(1)*d2wjis2(1) + bvvsum(2)*d2wjis2(2) + bvvsum(3)*d2wjis2(3) + &
                                   d1wjis(kl,1)*d1wjis(kk,1) + d1wjis(kl,2)*d1wjis(kk,2) + d1wjis(kl,3)*d1wjis(kk,3))
                sderv2(kl,kk) = sderv2(kl,kk) + d2eidwi22*dwji2s(kl)*dwji2s(kk) + deidwi2*d2wji2s2
              enddo
            enddo
          endif
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
!  Set COM coordinates
!
            nmk = natmol(k)
            if (lrigid.and.nmk.gt.0) then
              xcomk = molxyz(1,natinmol(k),nmk) - xcomi
              ycomk = molxyz(2,natinmol(k),nmk) - ycomi
              zcomk = molxyz(3,natinmol(k),nmk) - zcomi
            else
              xcomk = - xcomi
              ycomk = - ycomi
              zcomk = - zcomi
            endif
!
            d2eidr2 = d2eidbvsum2*d1i(ni)*d1i(nj)
!
            if (lstr) then
              call real1strterm(ndim,xki,yki,zki,xcomk,ycomk,zcomk,drik2ds,d2rik2dx2,d2rik2dsdx,d2rik2ds2,.false.)
              call cartstrterm(ndim,xki,yki,zki,xcomk,ycomk,zcomk,dxyzikds,d2xyzikdsdx,d2xyzikds2,.false.)
            endif
!
            derv2(ixf,ix) = derv2(ixf,ix) + d2eidr2*xki*xji
            derv2(iyf,ix) = derv2(iyf,ix) + d2eidr2*yki*xji
            derv2(izf,ix) = derv2(izf,ix) + d2eidr2*zki*xji
            derv2(ixf,iy) = derv2(ixf,iy) + d2eidr2*xki*yji
            derv2(iyf,iy) = derv2(iyf,iy) + d2eidr2*yki*yji
            derv2(izf,iy) = derv2(izf,iy) + d2eidr2*zki*yji
            derv2(ixf,iz) = derv2(ixf,iz) + d2eidr2*xki*zji
            derv2(iyf,iz) = derv2(iyf,iz) + d2eidr2*yki*zji
            derv2(izf,iz) = derv2(izf,iz) + d2eidr2*zki*zji
!
            derv2(kx,ix) = derv2(kx,ix) - d2eidr2*xki*xji
            derv2(ky,ix) = derv2(ky,ix) - d2eidr2*yki*xji
            derv2(kz,ix) = derv2(kz,ix) - d2eidr2*zki*xji
            derv2(kx,iy) = derv2(kx,iy) - d2eidr2*xki*yji
            derv2(ky,iy) = derv2(ky,iy) - d2eidr2*yki*yji
            derv2(kz,iy) = derv2(kz,iy) - d2eidr2*zki*yji
            derv2(kx,iz) = derv2(kx,iz) - d2eidr2*xki*zji
            derv2(ky,iz) = derv2(ky,iz) - d2eidr2*yki*zji
            derv2(kz,iz) = derv2(kz,iz) - d2eidr2*zki*zji
!
            derv2(jx,ix) = derv2(jx,ix) - d2eidr2*xki*xji
            derv2(jy,ix) = derv2(jy,ix) - d2eidr2*xki*yji
            derv2(jz,ix) = derv2(jz,ix) - d2eidr2*xki*zji
            derv2(jx,iy) = derv2(jx,iy) - d2eidr2*yki*xji
            derv2(jy,iy) = derv2(jy,iy) - d2eidr2*yki*yji
            derv2(jz,iy) = derv2(jz,iy) - d2eidr2*yki*zji
            derv2(jx,iz) = derv2(jx,iz) - d2eidr2*zki*xji
            derv2(jy,iz) = derv2(jy,iz) - d2eidr2*zki*yji
            derv2(jz,iz) = derv2(jz,iz) - d2eidr2*zki*zji
!
            if (lstr) then
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(ix,kl) = derv3(ix,kl) - xji*d2eidr2*drik2ds(ks)
                derv3(iy,kl) = derv3(iy,kl) - yji*d2eidr2*drik2ds(ks)
                derv3(iz,kl) = derv3(iz,kl) - zji*d2eidr2*drik2ds(ks)
              enddo
!
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(ix,kl) = derv3(ix,kl) - xki*d2eidr2*drij2ds(ks)
                derv3(iy,kl) = derv3(iy,kl) - yki*d2eidr2*drij2ds(ks)
                derv3(iz,kl) = derv3(iz,kl) - zki*d2eidr2*drij2ds(ks)
              enddo
!
              do kk = 1,nstrains
                ks = nstrptr(kk)
                do kl = 1,nstrains
                  kt = nstrptr(kl)
                  sderv2(kl,kk) = sderv2(kl,kk) + d2eidr2*(drij2ds(kt)*drik2ds(ks) + drik2ds(kt)*drij2ds(ks))
                enddo
              enddo
            endif
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
            derv2(jx,ix) = derv2(jx,ix) - deidwi2*d2wkj2(1,1) - d2eidwi22*dwki2(1)*dwji2(1)
            derv2(jy,ix) = derv2(jy,ix) - deidwi2*d2wkj2(1,2) - d2eidwi22*dwki2(1)*dwji2(2)
            derv2(jz,ix) = derv2(jz,ix) - deidwi2*d2wkj2(1,3) - d2eidwi22*dwki2(1)*dwji2(3)
            derv2(jx,iy) = derv2(jx,iy) - deidwi2*d2wkj2(2,1) - d2eidwi22*dwki2(2)*dwji2(1)
            derv2(jy,iy) = derv2(jy,iy) - deidwi2*d2wkj2(2,2) - d2eidwi22*dwki2(2)*dwji2(2)
            derv2(jz,iy) = derv2(jz,iy) - deidwi2*d2wkj2(2,3) - d2eidwi22*dwki2(2)*dwji2(3)
            derv2(jx,iz) = derv2(jx,iz) - deidwi2*d2wkj2(3,1) - d2eidwi22*dwki2(3)*dwji2(1)
            derv2(jy,iz) = derv2(jy,iz) - deidwi2*d2wkj2(3,2) - d2eidwi22*dwki2(3)*dwji2(2)
            derv2(jz,iz) = derv2(jz,iz) - deidwi2*d2wkj2(3,3) - d2eidwi22*dwki2(3)*dwji2(3)
!
            derv2(kx,ix) = derv2(kx,ix) - deidwi2*d2wkj2(1,1) - d2eidwi22*dwji2(1)*dwki2(1)
            derv2(ky,ix) = derv2(ky,ix) - deidwi2*d2wkj2(2,1) - d2eidwi22*dwji2(1)*dwki2(2)
            derv2(kz,ix) = derv2(kz,ix) - deidwi2*d2wkj2(3,1) - d2eidwi22*dwji2(1)*dwki2(3)
            derv2(kx,iy) = derv2(kx,iy) - deidwi2*d2wkj2(1,2) - d2eidwi22*dwji2(2)*dwki2(1)
            derv2(ky,iy) = derv2(ky,iy) - deidwi2*d2wkj2(2,2) - d2eidwi22*dwji2(2)*dwki2(2)
            derv2(kz,iy) = derv2(kz,iy) - deidwi2*d2wkj2(3,2) - d2eidwi22*dwji2(2)*dwki2(3)
            derv2(kx,iz) = derv2(kx,iz) - deidwi2*d2wkj2(1,3) - d2eidwi22*dwji2(3)*dwki2(1)
            derv2(ky,iz) = derv2(ky,iz) - deidwi2*d2wkj2(2,3) - d2eidwi22*dwji2(3)*dwki2(2)
            derv2(kz,iz) = derv2(kz,iz) - deidwi2*d2wkj2(3,3) - d2eidwi22*dwji2(3)*dwki2(3)
!
!  Strain second derivatives
!
            if (lstr) then
              do kk = 1,nstrains
                ks = nstrptr(kk)
                d1wkis(kk,1) = drik2ds(ks)*xki*rrik*(d1i(nj) - bv(nj)*rrik*rrik) + bv(nj)*rrik*dxyzikds(kk,1)
                d1wkis(kk,2) = drik2ds(ks)*yki*rrik*(d1i(nj) - bv(nj)*rrik*rrik) + bv(nj)*rrik*dxyzikds(kk,2)
                d1wkis(kk,3) = drik2ds(ks)*zki*rrik*(d1i(nj) - bv(nj)*rrik*rrik) + bv(nj)*rrik*dxyzikds(kk,3)
!
                dwki2s(kk) = 2.0_dp*(bvvsum(1)*d1wkis(kk,1) + bvvsum(2)*d1wkis(kk,2) + bvvsum(3)*d1wkis(kk,3))
!
                d2wki2sx(1) = 2.0_dp*(d1wkis(kk,1)*d1wji(1,1) + d1wkis(kk,2)*d1wji(1,2) + d1wkis(kk,3)*d1wji(1,3))
                d2wki2sx(2) = 2.0_dp*(d1wkis(kk,1)*d1wji(2,1) + d1wkis(kk,2)*d1wji(2,2) + d1wkis(kk,3)*d1wji(2,3))
                d2wki2sx(3) = 2.0_dp*(d1wkis(kk,1)*d1wji(3,1) + d1wkis(kk,2)*d1wji(3,2) + d1wkis(kk,3)*d1wji(3,3))
!
                derv3(ix,kk) = derv3(ix,kk) - d2eidwi22*dwki2s(kk)*dwji2(1) - deidwi2*d2wki2sx(1)
                derv3(iy,kk) = derv3(iy,kk) - d2eidwi22*dwki2s(kk)*dwji2(2) - deidwi2*d2wki2sx(2)
                derv3(iz,kk) = derv3(iz,kk) - d2eidwi22*dwki2s(kk)*dwji2(3) - deidwi2*d2wki2sx(3)
!
                d2wki2sx(1) = 2.0_dp*(d1wjis(kk,1)*d1wki(1,1) + d1wjis(kk,2)*d1wki(1,2) + d1wjis(kk,3)*d1wki(1,3))
                d2wki2sx(2) = 2.0_dp*(d1wjis(kk,1)*d1wki(2,1) + d1wjis(kk,2)*d1wki(2,2) + d1wjis(kk,3)*d1wki(2,3))
                d2wki2sx(3) = 2.0_dp*(d1wjis(kk,1)*d1wki(3,1) + d1wjis(kk,2)*d1wki(3,2) + d1wjis(kk,3)*d1wki(3,3))
!
                derv3(ix,kk) = derv3(ix,kk) - d2eidwi22*dwji2s(kk)*dwki2(1) - deidwi2*d2wki2sx(1)
                derv3(iy,kk) = derv3(iy,kk) - d2eidwi22*dwji2s(kk)*dwki2(2) - deidwi2*d2wki2sx(2)
                derv3(iz,kk) = derv3(iz,kk) - d2eidwi22*dwji2s(kk)*dwki2(3) - deidwi2*d2wki2sx(3)
              enddo
!
              do kk = 1,nstrains
                do kl = 1,nstrains
                  d2wji2s2 = 2.0_dp*(d1wkis(kk,1)*d1wjis(kl,1) + d1wkis(kk,2)*d1wjis(kl,2) + d1wkis(kk,3)*d1wjis(kl,3))
                  sderv2(kk,kl) = sderv2(kk,kl) + d2eidwi22*dwki2s(kk)*dwji2s(kl) + deidwi2*d2wji2s2
                  d2wji2s2 = 2.0_dp*(d1wkis(kl,1)*d1wjis(kk,1) + d1wkis(kl,2)*d1wjis(kk,2) + d1wkis(kl,3)*d1wjis(kk,3))
                  sderv2(kk,kl) = sderv2(kk,kl) + d2eidwi22*dwki2s(kl)*dwji2s(kk) + deidwi2*d2wji2s2
                enddo
              enddo
            endif
          enddo
        endif
      enddo
    endif
  enddo
!
  if (lgrad2) then
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
      nregioni = nregionno(nrelf2a(i))
      nregiontypi = nregiontype(nregioni,ncf)
      lslicei = lsliceatom(nsft + nrelf2a(i))
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
                if (lgrad1) then
                  d1vijdr = - vij*rrij/VBparB(3,nv)
                  if (lgrad2) then
                    d2vijdr = - d1vijdr*rrij*(1.0_dp/VBparB(3,nv) + rrij)
                  endif
                endif
              else
!
!  Power law form
!
                vij = VBwgtB(nv)*(VBparB(1,nv)/rij)**VBparB(2,nv)
                if (lgrad1) then
                  d1vijdr = - VBparB(2,nv)*vij*rrij*rrij
                  if (lgrad2) then
                    d2vijdr = - (VBparB(2,nv) + 2.0_dp)*d1vijdr*rrij*rrij
                  endif
                endif
              endif
              bv(ni) = bv(ni) + vij
              bvv(1,ni) = bvv(1,ni) + vij*xji*rrij
              bvv(2,ni) = bvv(2,ni) + vij*yji*rrij
              bvv(3,ni) = bvv(3,ni) + vij*zji*rrij
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
!  Derivatives of E_bv
!
        deijdr = deidbvsum*d1i(ni)
!
        if (lstr) then
          call real1strterm(ndim,xji,yji,zji,xcom,ycom,zcom,drij2ds,d2rij2dx2,d2rij2dsdx,d2rij2ds2,.true.)
        endif
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
        if (lstr) then
          call cartstrterm(ndim,xji,yji,zji,xcom,ycom,zcom,dxyzijds,d2xyzijdsdx,d2xyzijds2,.true.)
        endif
        if (lstr) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            d1wjis(kl,1) = drij2ds(ks)*xji*rrij*(d1i(ni) - bv(ni)*rrij*rrij) + bv(ni)*rrij*dxyzijds(kl,1)
            d1wjis(kl,2) = drij2ds(ks)*yji*rrij*(d1i(ni) - bv(ni)*rrij*rrij) + bv(ni)*rrij*dxyzijds(kl,2)
            d1wjis(kl,3) = drij2ds(ks)*zji*rrij*(d1i(ni) - bv(ni)*rrij*rrij) + bv(ni)*rrij*dxyzijds(kl,3)
!
            dwji2s(kl) = 2.0_dp*(bvvsum(1)*d1wjis(kl,1) + bvvsum(2)*d1wjis(kl,2) + bvvsum(3)*d1wjis(kl,3))
          enddo
        endif
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
        derv2(ix,jx) = derv2(ix,jx) - d2eidr2*d2rij2dx2(1)
        derv2(iy,jx) = derv2(iy,jx) - d2eidr2*d2rij2dx2(6)
        derv2(iz,jx) = derv2(iz,jx) - d2eidr2*d2rij2dx2(5)
        derv2(ix,jy) = derv2(ix,jy) - d2eidr2*d2rij2dx2(6)
        derv2(iy,jy) = derv2(iy,jy) - d2eidr2*d2rij2dx2(2)
        derv2(iz,jy) = derv2(iz,jy) - d2eidr2*d2rij2dx2(4)
        derv2(ix,jz) = derv2(ix,jz) - d2eidr2*d2rij2dx2(5)
        derv2(iy,jz) = derv2(iy,jz) - d2eidr2*d2rij2dx2(4)
        derv2(iz,jz) = derv2(iz,jz) - d2eidr2*d2rij2dx2(3)
        derv2(ix,jx) = derv2(ix,jx) - deijdr
        derv2(iy,jy) = derv2(iy,jy) - deijdr
        derv2(iz,jz) = derv2(iz,jz) - deijdr
!
        if (lstr) then
!
!  Mixed derivatives
!
          if (i.ne.j) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(jx,kl) = derv3(jx,kl) + deijdr*d2rij2dsdx(ks,1) + xji*d2eidr2*drij2ds(ks)
              derv3(jy,kl) = derv3(jy,kl) + deijdr*d2rij2dsdx(ks,2) + yji*d2eidr2*drij2ds(ks)
              derv3(jz,kl) = derv3(jz,kl) + deijdr*d2rij2dsdx(ks,3) + zji*d2eidr2*drij2ds(ks)
            enddo
          endif
        endif
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
        derv2(ix,jx) = derv2(ix,jx) - d2eidwi22*dwji2(1)*dwji2(1) - deidwi2*d2wji2(1,1)
        derv2(iy,jx) = derv2(iy,jx) - d2eidwi22*dwji2(2)*dwji2(1) - deidwi2*d2wji2(2,1)
        derv2(iz,jx) = derv2(iz,jx) - d2eidwi22*dwji2(3)*dwji2(1) - deidwi2*d2wji2(3,1)
        derv2(ix,jy) = derv2(ix,jy) - d2eidwi22*dwji2(1)*dwji2(2) - deidwi2*d2wji2(1,2)
        derv2(iy,jy) = derv2(iy,jy) - d2eidwi22*dwji2(2)*dwji2(2) - deidwi2*d2wji2(2,2)
        derv2(iz,jy) = derv2(iz,jy) - d2eidwi22*dwji2(3)*dwji2(2) - deidwi2*d2wji2(3,2)
        derv2(ix,jz) = derv2(ix,jz) - d2eidwi22*dwji2(1)*dwji2(3) - deidwi2*d2wji2(1,3)
        derv2(iy,jz) = derv2(iy,jz) - d2eidwi22*dwji2(2)*dwji2(3) - deidwi2*d2wji2(2,3)
        derv2(iz,jz) = derv2(iz,jz) - d2eidwi22*dwji2(3)*dwji2(3) - deidwi2*d2wji2(3,3)
!
        if (lstr) then
          do kk = 1,nstrains
            ks = nstrptr(kk)
!
            d2wjisx(1,1) = d2rij2dsdx(ks,1)*xji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,1,1) + &
                           (rrij*d1i(ni) - bv(ni)*rrij3)*(xji*dxyzijds(kk,1) + drij2ds(ks)) + &
                           drij2ds(ks)*xji*xji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
            d2wjisx(2,1) = d2rij2dsdx(ks,2)*xji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,1,2) + &
                           (rrij*d1i(ni) - bv(ni)*rrij3)*(yji*dxyzijds(kk,1)) + &
                           drij2ds(ks)*yji*xji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
            d2wjisx(3,1) = d2rij2dsdx(ks,3)*xji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,1,3) + &
                           (rrij*d1i(ni) - bv(ni)*rrij3)*(zji*dxyzijds(kk,1)) + &
                           drij2ds(ks)*zji*xji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
!
            d2wjisx(1,2) = d2rij2dsdx(ks,1)*yji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,2,1) + &
                           (rrij*d1i(ni) - bv(ni)*rrij3)*(xji*dxyzijds(kk,2)) + &
                           drij2ds(ks)*xji*yji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
            d2wjisx(2,2) = d2rij2dsdx(ks,2)*yji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,2,2) + &
                           (rrij*d1i(ni) - bv(ni)*rrij3)*(yji*dxyzijds(kk,2) + drij2ds(ks)) + &
                           drij2ds(ks)*yji*yji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
            d2wjisx(3,2) = d2rij2dsdx(ks,3)*yji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,2,3) + &
                           (rrij*d1i(ni) - bv(ni)*rrij3)*(zji*dxyzijds(kk,2)) + &
                           drij2ds(ks)*zji*yji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
!
            d2wjisx(1,3) = d2rij2dsdx(ks,1)*zji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,3,1) + &
                           (rrij*d1i(ni) - bv(ni)*rrij3)*(xji*dxyzijds(kk,3)) + &
                           drij2ds(ks)*xji*zji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
            d2wjisx(2,3) = d2rij2dsdx(ks,2)*zji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,3,2) + &
                           (rrij*d1i(ni) - bv(ni)*rrij3)*(yji*dxyzijds(kk,3)) + &
                           drij2ds(ks)*yji*zji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
            d2wjisx(3,3) = d2rij2dsdx(ks,3)*zji*(rrij*d1i(ni) - bv(ni)*rrij3) + bv(ni)*rrij*d2xyzijdsdx(kk,3,3) + &
                           (rrij*d1i(ni) - bv(ni)*rrij3)*(zji*dxyzijds(kk,3) + drij2ds(ks)) + &
                           drij2ds(ks)*zji*zji*(rrij*d2i(ni) - 2.0_dp*rrij3*d1i(ni) + 3.0_dp*bv(ni)*rrij5)
!
            d2wji2sx(1) = 2.0_dp*(bvvsum(1)*d2wjisx(1,1) + bvvsum(2)*d2wjisx(1,2) + bvvsum(3)*d2wjisx(1,3) + &
                                  d1wjis(kk,1)*d1wji(1,1) + d1wjis(kk,2)*d1wji(1,2) + d1wjis(kk,3)*d1wji(1,3))
            d2wji2sx(2) = 2.0_dp*(bvvsum(1)*d2wjisx(2,1) + bvvsum(2)*d2wjisx(2,2) + bvvsum(3)*d2wjisx(2,3) + &
                                  d1wjis(kk,1)*d1wji(2,1) + d1wjis(kk,2)*d1wji(2,2) + d1wjis(kk,3)*d1wji(2,3))
            d2wji2sx(3) = 2.0_dp*(bvvsum(1)*d2wjisx(3,1) + bvvsum(2)*d2wjisx(3,2) + bvvsum(3)*d2wjisx(3,3) + &
                                  d1wjis(kk,1)*d1wji(3,1) + d1wjis(kk,2)*d1wji(3,2) + d1wjis(kk,3)*d1wji(3,3))
!
            derv3(jx,kk) = derv3(jx,kk) + d2eidwi22*dwji2s(kk)*dwji2(1) + deidwi2*d2wji2sx(1)
            derv3(jy,kk) = derv3(jy,kk) + d2eidwi22*dwji2s(kk)*dwji2(2) + deidwi2*d2wji2sx(2)
            derv3(jz,kk) = derv3(jz,kk) + d2eidwi22*dwji2s(kk)*dwji2(3) + deidwi2*d2wji2sx(3)
          enddo
        endif
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
!  Set COM coordinates
!
          nmk = natmol(k)
          if (lrigid.and.nmk.gt.0) then
            xcomk = molxyz(1,natinmol(k),nmk) - xcomi
            ycomk = molxyz(2,natinmol(k),nmk) - ycomi
            zcomk = molxyz(3,natinmol(k),nmk) - zcomi
          else
            xcomk = - xcomi
            ycomk = - ycomi
            zcomk = - zcomi
          endif
!
          d2eidr2 = d2eidbvsum2*d1i(ni)*d1i(nj)
!
          if (lstr) then
            call real1strterm(ndim,xki,yki,zki,xcomk,ycomk,zcomk,drik2ds,d2rik2dx2,d2rik2dsdx,d2rik2ds2,.false.)
            call cartstrterm(ndim,xki,yki,zki,xcomk,ycomk,zcomk,dxyzikds,d2xyzikdsdx,d2xyzikds2,.false.)
          endif
!
          derv2(kx,jx) = derv2(kx,jx) + d2eidr2*xji*xki
          derv2(ky,jx) = derv2(ky,jx) + d2eidr2*xji*yki
          derv2(kz,jx) = derv2(kz,jx) + d2eidr2*xji*zki
          derv2(kx,jy) = derv2(kx,jy) + d2eidr2*yji*xki
          derv2(ky,jy) = derv2(ky,jy) + d2eidr2*yji*yki
          derv2(kz,jy) = derv2(kz,jy) + d2eidr2*yji*zki
          derv2(kx,jz) = derv2(kx,jz) + d2eidr2*zji*xki
          derv2(ky,jz) = derv2(ky,jz) + d2eidr2*zji*yki
          derv2(kz,jz) = derv2(kz,jz) + d2eidr2*zji*zki
!
          derv2(ix,jx) = derv2(ix,jx) - d2eidr2*xki*xji
          derv2(iy,jx) = derv2(iy,jx) - d2eidr2*yki*xji
          derv2(iz,jx) = derv2(iz,jx) - d2eidr2*zki*xji
          derv2(ix,jy) = derv2(ix,jy) - d2eidr2*xki*yji
          derv2(iy,jy) = derv2(iy,jy) - d2eidr2*yki*yji
          derv2(iz,jy) = derv2(iz,jy) - d2eidr2*zki*yji
          derv2(ix,jz) = derv2(ix,jz) - d2eidr2*xki*zji
          derv2(iy,jz) = derv2(iy,jz) - d2eidr2*yki*zji
          derv2(iz,jz) = derv2(iz,jz) - d2eidr2*zki*zji
!
          if (lstr) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              derv3(jx,kl) = derv3(jx,kl) + xji*d2eidr2*drik2ds(ks)
              derv3(jy,kl) = derv3(jy,kl) + yji*d2eidr2*drik2ds(ks)
              derv3(jz,kl) = derv3(jz,kl) + zji*d2eidr2*drik2ds(ks)
            enddo
          endif
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
          derv2(ix,jx) = derv2(ix,jx) - deidwi2*d2wkj2(1,1) - d2eidwi22*dwki2(1)*dwji2(1)
          derv2(iy,jx) = derv2(iy,jx) - deidwi2*d2wkj2(2,1) - d2eidwi22*dwki2(2)*dwji2(1)
          derv2(iz,jx) = derv2(iz,jx) - deidwi2*d2wkj2(3,1) - d2eidwi22*dwki2(3)*dwji2(1)
          derv2(ix,jy) = derv2(ix,jy) - deidwi2*d2wkj2(1,2) - d2eidwi22*dwki2(1)*dwji2(2)
          derv2(iy,jy) = derv2(iy,jy) - deidwi2*d2wkj2(2,2) - d2eidwi22*dwki2(2)*dwji2(2)
          derv2(iz,jy) = derv2(iz,jy) - deidwi2*d2wkj2(3,2) - d2eidwi22*dwki2(3)*dwji2(2)
          derv2(ix,jz) = derv2(ix,jz) - deidwi2*d2wkj2(1,3) - d2eidwi22*dwki2(1)*dwji2(3)
          derv2(iy,jz) = derv2(iy,jz) - deidwi2*d2wkj2(2,3) - d2eidwi22*dwki2(2)*dwji2(3)
          derv2(iz,jz) = derv2(iz,jz) - deidwi2*d2wkj2(3,3) - d2eidwi22*dwki2(3)*dwji2(3)
!
          derv2(kx,jx) = derv2(kx,jx) + deidwi2*d2wkj2(1,1) + d2eidwi22*dwji2(1)*dwki2(1)
          derv2(ky,jx) = derv2(ky,jx) + deidwi2*d2wkj2(2,1) + d2eidwi22*dwji2(1)*dwki2(2)
          derv2(kz,jx) = derv2(kz,jx) + deidwi2*d2wkj2(3,1) + d2eidwi22*dwji2(1)*dwki2(3)
          derv2(kx,jy) = derv2(kx,jy) + deidwi2*d2wkj2(1,2) + d2eidwi22*dwji2(2)*dwki2(1)
          derv2(ky,jy) = derv2(ky,jy) + deidwi2*d2wkj2(2,2) + d2eidwi22*dwji2(2)*dwki2(2)
          derv2(kz,jy) = derv2(kz,jy) + deidwi2*d2wkj2(3,2) + d2eidwi22*dwji2(2)*dwki2(3)
          derv2(kx,jz) = derv2(kx,jz) + deidwi2*d2wkj2(1,3) + d2eidwi22*dwji2(3)*dwki2(1)
          derv2(ky,jz) = derv2(ky,jz) + deidwi2*d2wkj2(2,3) + d2eidwi22*dwji2(3)*dwki2(2)
          derv2(kz,jz) = derv2(kz,jz) + deidwi2*d2wkj2(3,3) + d2eidwi22*dwji2(3)*dwki2(3)
!
!  Strain second derivatives
!
          if (lstr) then
            do kk = 1,nstrains
              ks = nstrptr(kk)
              d1wkis(kk,1) = drik2ds(ks)*xki*rrik*(d1i(nj) - bv(nj)*rrik*rrik) + bv(nj)*rrik*dxyzikds(kk,1)
              d1wkis(kk,2) = drik2ds(ks)*yki*rrik*(d1i(nj) - bv(nj)*rrik*rrik) + bv(nj)*rrik*dxyzikds(kk,2)
              d1wkis(kk,3) = drik2ds(ks)*zki*rrik*(d1i(nj) - bv(nj)*rrik*rrik) + bv(nj)*rrik*dxyzikds(kk,3)
!
              dwki2s(kk) = 2.0_dp*(bvvsum(1)*d1wkis(kk,1) + bvvsum(2)*d1wkis(kk,2) + bvvsum(3)*d1wkis(kk,3))
!
              d2wki2sx(1) = 2.0_dp*(d1wkis(kk,1)*d1wji(1,1) + d1wkis(kk,2)*d1wji(1,2) + d1wkis(kk,3)*d1wji(1,3))
              d2wki2sx(2) = 2.0_dp*(d1wkis(kk,1)*d1wji(2,1) + d1wkis(kk,2)*d1wji(2,2) + d1wkis(kk,3)*d1wji(2,3))
              d2wki2sx(3) = 2.0_dp*(d1wkis(kk,1)*d1wji(3,1) + d1wkis(kk,2)*d1wji(3,2) + d1wkis(kk,3)*d1wji(3,3))
!
              derv3(jx,kk) = derv3(jx,kk) + d2eidwi22*dwki2s(kk)*dwji2(1) + deidwi2*d2wki2sx(1)
              derv3(jy,kk) = derv3(jy,kk) + d2eidwi22*dwki2s(kk)*dwji2(2) + deidwi2*d2wki2sx(2)
              derv3(jz,kk) = derv3(jz,kk) + d2eidwi22*dwki2s(kk)*dwji2(3) + deidwi2*d2wki2sx(3)
            enddo
          endif
        enddo
      enddo jloop
    enddo iloop
  endif
!
!  Closing banner for energy decomposition
!
  if (lPrintVB) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  Free local memory
!
  deallocate(d2i,stat=status)
  if (status/=0) call deallocate_error('valencebondd','d2i')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('valencebondd','d1i')
  deallocate(bvv,stat=status)
  if (status/=0) call deallocate_error('valencebondd','vbv')
  deallocate(bv,stat=status)
  if (status/=0) call deallocate_error('valencebondd','vb')
  deallocate(nvbijptr,stat=status)
  if (status/=0) call deallocate_error('valencebondd','nvbijptr')
!
  t2 = g_cpu_time()
  tval = tval + t2 - t1
#ifdef TRACE
  call trace_out('valencebondd')
#endif
!
  return
  end
