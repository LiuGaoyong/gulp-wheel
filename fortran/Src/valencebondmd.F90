  subroutine valencebondmd(eval,lgrad1)
!
!  Calculates the energy and first derivatives for the Valance Bond potentials.
!
!  NB: Does not currently support lfreeze
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
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
  use spatial
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
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ic
  integer(i4)                                      :: j
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: ni
  integer(i4)                                      :: nmi
  integer(i4)                                      :: nmj
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
  logical                                          :: lreg2ok
  logical                                          :: lslicei
  real(dp)                                         :: bvsum
  real(dp)                                         :: bvvsum(3)
  real(dp)                                         :: g_cpu_time
  real(dp),    dimension(:),     allocatable, save :: d1i
  real(dp)                                         :: d1wji(3,3)
  real(dp)                                         :: d1wjis(6,3)
  real(dp)                                         :: dwji2(3)
  real(dp)                                         :: dwji2s(6)
  real(dp)                                         :: deidbvsum
  real(dp)                                         :: deidwi2
  real(dp)                                         :: deijdr
  real(dp)                                         :: drij2ds(6)
  real(dp)                                         :: d2rij2dx2(6)
  real(dp)                                         :: d2rij2ds2(6,6)
  real(dp)                                         :: d2rij2dsdx(6,3)
  real(dp)                                         :: dxyzijds(6,3)
  real(dp)                                         :: d2xyzijdsdx(6,3,3)
  real(dp)                                         :: d2xyzijds2(6,6,3)
  real(dp)                                         :: eij
  real(dp)                                         :: eijbv
  real(dp)                                         :: eijbvv
  real(dp)                                         :: rij
  real(dp)                                         :: rrij
  real(dp)                                         :: rrij3
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp)                                         :: vij
  real(dp)                                         :: d1vijdr
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
  real(dp),    dimension(:),     allocatable, save :: bv
  real(dp),    dimension(:,:),   allocatable, save :: bvv
#ifdef TRACE
  call trace_in('valencebondmd')
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
  if (ioproc.and.index(keyword,'debu').ne.0) then
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
  if (status/=0) call outofmemory('valencebondmd','nvbijptr')
  allocate(bv(maxvbnbr),stat=status)
  if (status/=0) call outofmemory('valencebondmd','vb')
  allocate(bvv(3,maxvbnbr),stat=status)
  if (status/=0) call outofmemory('valencebondmd','vbv')
!
  if (lgrad1) then
    allocate(d1i(maxvbnbr),stat=status)
    if (status/=0) call outofmemory('valencebondmd','d1i')
  else
    allocate(d1i(1),stat=status)
    if (status/=0) call outofmemory('valencebondmd','d1i')
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
  do ic = 1,natompernode
    i = natomnodeptr(ic)
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
              endif
            else
!
!  Power law form
!
              vij = VBwgtB(nv)*(VBparB(1,nv)/rij)**VBparB(2,nv)
              if (lgrad1) then
                d1vijdr = - VBparB(2,nv)*vij*rrij*rrij
              endif
            endif
            bv(ni) = bv(ni) + vij
            bvv(1,ni) = bvv(1,ni) + vij*xji*rrij
            bvv(2,ni) = bvv(2,ni) + vij*yji*rrij
            bvv(3,ni) = bvv(3,ni) + vij*zji*rrij
            if (lgrad1) then
              d1i(ni) = d1i(ni) + d1vijdr
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
      deidwi2     = 0.0_dp
      do nve = 1,nvalener
        if (nVBspecE(nve).eq.nati) then
          if (nVBtypeE(nve).eq.ntypi.or.nVBtypeE(nve).eq.0) then
            deidbvsum   = deidbvsum   + 2.0_dp*VBparE(1,nve)*(bvsum - VBparE(2,nve))
            deidwi2     = deidwi2     + 2.0_dp*VBparE(3,nve)*(wi2 - VBparE(4,nve)**2)
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
        if (lstr) then
          call real1strterm(ndim,xji,yji,zji,xcom,ycom,zcom,drij2ds,d2rij2dx2,d2rij2dsdx,d2rij2ds2,.false.)
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
        if (lstr) then
          call cartstrterm(ndim,xji,yji,zji,xcom,ycom,zcom,dxyzijds,d2xyzijdsdx,d2xyzijds2,.false.)
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
      enddo
    endif
  enddo
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
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('valencebondmd','d1i')
  deallocate(bvv,stat=status)
  if (status/=0) call deallocate_error('valencebondmd','vbv')
  deallocate(bv,stat=status)
  if (status/=0) call deallocate_error('valencebondmd','vb')
  deallocate(nvbijptr,stat=status)
  if (status/=0) call deallocate_error('valencebondmd','nvbijptr')
!
  t2 = g_cpu_time()
  tval = tval + t2 - t1
#ifdef TRACE
  call trace_out('valencebondmd')
#endif
!
  return
  end
