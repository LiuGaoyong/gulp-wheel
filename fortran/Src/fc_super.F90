  subroutine build_fc_supercell(nfcsuper)
!
!  Generate supercell for force constant calculation
!  NB: Assumes that the supercell will only be used for a phonon calculation and then destroyed.
!      This means the setup of things not required for this will be incomplete!
!
!   7/23 Created from build_supercell
!   7/23 extpot added
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
  use bondorderdata,  only : nboQ
  use bondvalence,    only : nvalbond
  use configurations
  use control
  use current
  use eemdata,        only : nqrnow
  use element,        only : lqeq, maxele, lSandM, lpacha, lgasteiger
  use m_vb_nbr,       only : vb_getcut
  use molecule
  use parallel
  use partial 
  use polarise
  use spatial,        only : lspatialok
  use spatialbo,      only : lspatialBOok
  use shells
  use species
  use sutton,         only : lsuttonc
  use thresholds,     only : thresh_q, thresh_c6
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: nfcsuper(3)   ! Force constant supercell dimensions
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: ind
  integer(i4)                                  :: isfct
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: kk
  integer(i4)                                  :: nnew
  integer(i4)                                  :: ntotal
  integer(i4), dimension(:), allocatable       :: ntemp
  integer(i4), dimension(:), allocatable       :: ntemp2
  logical,     dimension(:), allocatable       :: ltemp
  logical                                      :: lvariablecharge
  integer(i4)                                  :: status
  real(dp)                                     :: c6trm
  real(dp)                                     :: qtrm
  real(dp)                                     :: rix
  real(dp)                                     :: riy
  real(dp)                                     :: riz
  real(dp),    dimension(:), allocatable       :: rtemp
  real(dp)                                     :: xadd
  real(dp)                                     :: yadd
  real(dp)                                     :: zadd
  real(dp),    dimension(:), allocatable       :: xtmp
  real(dp),    dimension(:), allocatable       :: ytmp
  real(dp),    dimension(:), allocatable       :: ztmp
#ifdef TRACE
  call trace_in('build_fc_super')
#endif
!
!  Check that system is periodic
!
  if (ndim.lt.1) then
    call outerror('fc_supercell cannot be specified for a cluster',0_i4)
    call stopnow('build_fc_supercell')
  endif
!
!  Convert supercell index into components
!
  ix = nfcsuper(1)
  iy = nfcsuper(2)
  iz = nfcsuper(3)
!
!  Check to see whether the atomic arrays need re-allocating
!
  isfct = ix*iy*iz
  ntotal = numat*isfct
  nnew = ntotal - numat
  if (ntotal.gt.maxat) then
    maxat = ntotal
    call changemaxat
  endif
  if (nnew+nasum.gt.maxatot) then
    maxatot = nnew + nasum
    call changemaxatot
  endif
!
!  Expand ion attributes for supercell without symmetry
!
  ind = nsft - isfct
  allocate(ltemp(numat),stat=status)
  if (status/=0) call outofmemory('build_fc_supercell','ltemp')
  allocate(ntemp(numat),stat=status)
  if (status/=0) call outofmemory('build_fc_supercell','ntemp')
  allocate(ntemp2(numat),stat=status)
  if (status/=0) call outofmemory('build_fc_supercell','ntemp2')
  allocate(rtemp(numat),stat=status)
  if (status/=0) call outofmemory('build_fc_supercell','rtemp')
  do i = 1,numat
    ltemp(i) = lbsmat(nrelf2a(i))
    ntemp(i) = nregionno(nrelf2a(i))
    ntemp2(i) = nspecptr(nrelf2a(i))
    rtemp(i) = extpot(nrelf2a(i))
  enddo
!
  do i = numat,1,-1
    ind = isfct*(i-1)
    do j = 1,isfct
      ind = ind + 1
      nat(ind) = nat(i)
      iatn(ind) = nat(i)
      neqv(ind) = 1
      nrelf2a(ind) = i
      nrela2f(ind) = i
      natype(ind) = nftype(i)
      nftype(ind) = nftype(i)
      nspecptr(ind) = ntemp2(i)
      occua(ind) = occuf(i)
      occuf(ind) = occuf(i)
      qa(ind) = qf(i)
      qf(ind) = qf(i)
      cna(ind) = cnf(i)
      cnf(ind) = cnf(i)
      oxa(ind) = oxf(i)
      oxf(ind) = oxf(i)
      rada(ind) = radf(i)
      radf(ind) = radf(i)
      nregionno(ind) = ntemp(i)
      lbsmat(ind) = ltemp(i)
      extpot(ind) = rtemp(i)
      if (lc6one) then
        c6a(ind) = c6f(i)
        c6f(ind) = c6f(i)
      endif
    enddo
  enddo
!
  deallocate(rtemp,stat=status)
  if (status/=0) call deallocate_error('build_fc_supercell','rtemp')
  deallocate(ntemp2,stat=status)
  if (status/=0) call deallocate_error('build_fc_supercell','ntemp2')
  deallocate(ntemp,stat=status)
  if (status/=0) call deallocate_error('build_fc_supercell','ntemp')
  deallocate(ltemp,stat=status)
  if (status/=0) call deallocate_error('build_fc_supercell','ltemp')
!
  allocate(xtmp(numat),stat=status)
  if (status/=0) call outofmemory('build_fc_supercell','xtmp')
  allocate(ytmp(numat),stat=status)
  if (status/=0) call outofmemory('build_fc_supercell','ytmp')
  allocate(ztmp(numat),stat=status)
  if (status/=0) call outofmemory('build_fc_supercell','ztmp')
!****************
!  Coordinates  *
!****************
!
!  Correct fractional coordinates for new cell size
!
  xtmp(1:numat) = xfrac(1:numat)
  ytmp(1:numat) = yfrac(1:numat)
  ztmp(1:numat) = zfrac(1:numat)
!
  ind = 0
  rix = 1.0_dp/dble(ix)
  riy = 1.0_dp/dble(iy)
  riz = 1.0_dp/dble(iz)
  do i = 1,numat
    do ii = 1,ix
      xadd = rix*(ii-1)
      do jj = 1,iy
        do kk = 1,iz
          ind = ind + 1
          xfrac(ind) = rix*xtmp(i) + xadd
        enddo
      enddo
    enddo
  enddo
  if (ndim.ge.2) then
    ind = 0
    do i = 1,numat
      do ii = 1,ix
        do jj = 1,iy
          yadd = riy*(jj-1)
          do kk = 1,iz
            ind = ind + 1
            yfrac(ind) = riy*ytmp(i) + yadd
          enddo
        enddo
      enddo
    enddo
  else
    ind = 0
    do i = 1,numat
      do ii = 1,ix
        do jj = 1,iy
          do kk = 1,iz
            ind = ind + 1
            yfrac(ind) = ytmp(i)
          enddo
        enddo
      enddo
    enddo
  endif
  if (ndim.eq.3) then
    ind = 0
    do i = 1,numat
      do ii = 1,ix
        do jj = 1,iy
          do kk = 1,iz
            zadd = riz*(kk-1)
            ind = ind + 1
            zfrac(ind) = riz*ztmp(i) + zadd
          enddo
        enddo
      enddo
    enddo
  else
    ind = 0
    do i = 1,numat
      do ii = 1,ix
        do jj = 1,iy
          do kk = 1,iz
            ind = ind + 1
            zfrac(ind) = ztmp(i)
          enddo
        enddo
      enddo
    enddo
  endif
!
  xafrac(1:numat) = xfrac(1:numat)
  yafrac(1:numat) = yfrac(1:numat)
  zafrac(1:numat) = zfrac(1:numat)
!
  deallocate(ztmp,stat=status)
  if (status/=0) call deallocate_error('build_fc_supercell','ztmp')
  deallocate(ytmp,stat=status)
  if (status/=0) call deallocate_error('build_fc_supercell','ytmp')
  deallocate(xtmp,stat=status)
  if (status/=0) call deallocate_error('build_fc_supercell','xtmp')
!********************
!  Unit cell setup  *
!********************
!
!  Change cell to supercell
!
  do j = 1,3
    rv(j,1) = dble(ix)*rv(j,1)
    rv(j,2) = dble(iy)*rv(j,2)
    rv(j,3) = dble(iz)*rv(j,3)
  enddo
!
!  Correct Cartesian coordinates for new cell size
!
  do i = 1,ntotal
    xclat(i) = rv(1,1)*xfrac(i) + rv(1,2)*yfrac(i) + rv(1,3)*zfrac(i)
    yclat(i) = rv(2,1)*xfrac(i) + rv(2,2)*yfrac(i) + rv(2,3)*zfrac(i)
    zclat(i) = rv(3,1)*xfrac(i) + rv(3,2)*yfrac(i) + rv(3,3)*zfrac(i)
  enddo
!
  xalat(1:ntotal) = xclat(ntotal)
  yalat(1:ntotal) = yclat(ntotal)
  zalat(1:ntotal) = zclat(ntotal)
!
  if (ndim.eq.3) then
    call uncell3D(rv,a,b,c,alpha,beta,gamma)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    if (c.gt.1.0d-12) then
      recipc = 1.0_dp/c
    else
      recipc = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = rv(2,1)
    r1z = rv(3,1)
    r2x = rv(1,2)
    r2y = rv(2,2)
    r2z = rv(3,2)
    r3x = rv(1,3)
    r3y = rv(2,3)
    r3z = rv(3,3)
  elseif (ndim.eq.2) then
    call uncell2D(rv,a,b,alpha)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    if (b.gt.1.0d-12) then
      recipb = 1.0_dp/b
    else
      recipb = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = rv(2,1)
    r1z = 0.0_dp
    r2x = rv(1,2)
    r2y = rv(2,2)
    r2z = 0.0_dp
    r3x = 0.0_dp
    r3y = 0.0_dp
    r3z = 0.0_dp
  elseif (ndim.eq.1) then
    call uncell1D(rv,a)
    if (a.gt.1.0d-12) then
      recipa = 1.0_dp/a
    else
      recipa = 0.0_dp
    endif
    r1x = rv(1,1)
    r1y = 0.0_dp
    r1z = 0.0_dp
    r2x = 0.0_dp
    r2y = 0.0_dp
    r2z = 0.0_dp
    r3x = 0.0_dp
    r3y = 0.0_dp
    r3z = 0.0_dp
  endif
!
!  Set up cell lists
!
  call rlist
  call rlist2
!
!  Change number of atoms 
!
  numat = ntotal
  nasym = ntotal
!*********************
!  Misc setup tasks  *
!*********************
  totalcharge = totalcharge*dble(isfct)
  nbsmat = nbsmat*isfct
  nqrnow(1:ntotal) = 1
  numofspec(1:nspec) = 0
  do i = 1,ntotal
    ii = nspecptr(i)
    numofspec(ii) = numofspec(ii) + 1
  enddo
  if (lewald) call setewald
  lspatialok = .false.     ! Turn off spatial for now
  lspatialBOok = .false.   ! Turn off spatial for now
!
  lvariablecharge = (((leem.or.lqeq.or.lSandM.or.lpacha.or.lgfnff.or.leht).and.(.not.lnoqeem)).or.(nboQ.gt.0)) 
  ncharge = 0
  do i = 1,ntotal
    qtrm = qf(i)*occuf(i)
    if (abs(qtrm).gt.thresh_q.or.lvariablecharge) then
      ncharge = ncharge + 1
      nchargeptr(ncharge) = i
    endif
  enddo
  if (lc6one) then
    nchargec6 = 0
    do i = 1,ntotal
      qtrm = abs(qf(i)*occuf(i))
      c6trm = abs(c6f(i)*occuf(i))
      if (qtrm.gt.thresh_q.or.c6trm.gt.thresh_c6) then
        nchargec6 = nchargec6 + 1
        nchargec6ptr(nchargec6) = i
      endif
    enddo
  else
    nchargec6 = ncharge
    nchargec6ptr(1:ncharge) = nchargeptr(1:ncharge)
  endif
!
  ncore  = 0
  nshell = 0
  do i = 1,ntotal
    if (nat(i).gt.maxele) then
      nshell = nshell + 1
      nshptr(nshell) = i
      ncoshptr(i) = nshell
    else
      ncore = ncore + 1
      ncoptr(ncore) = i
      ncoshptr(i) = ncore
    endif
  enddo
!*************************************
!  Molecule and bonding information  *
!*************************************
  if (index(keyword,'mol').ne.0) then 
    if (lmolatom.or.lmolrigid) then
      call setmola
    else
      call setmol
    endif
  else
! 
!  If molecule is not being called then set nasymnomol/numatnomol/ncorenomol
! 
    nasymnomol = nasym
    numatnomol = numat
    ncorenomol = ncore
    nbsmatnomol = nbsmat
    do i = 1,nasymnomol
      nasymnomolptr(i) = i
      nasymnomolrptr(i) = i
    enddo
    do i = 1,numatnomol          
      numatnomolptr(i) = i
      numatnomolrptr(i) = i
    enddo
    do i = 1,ncorenomol                        
      ncorenomolptr(i) = i
      ncorenomolrptr(i) = i
    enddo
  endif
!*****************************************************************
!  Final tasks that require numat to be set before being called  *
!*****************************************************************
!
!  EAM set up if needed
!
  if (lsuttonc) call seteam
!
!  Polarisability set up if needed
!
  if (lpolar) call setpolar
!
!  Valence bond set up if needed
!
  if (nvalbond.gt.0) call vb_getcut
!
!  Setup Gasteiger charges if needed since they are not geometry dependent
!
  if (lgasteiger) then
    call gasteiger(.false.)
  endif
!
!  Setup AA-CLP one-centre force field parameters
!
  if (laaclp) then
    call setaaclp
  endif
!
!  GFNFF setup
!
  if (lgfnff) then
    call setup_gfnff
  endif
!
!  Set breathing shell pointer
!
  call setbsmptr(nbs,nbss,nbsptr,(ndim.eq.2))
!
!  Set partial occupancy pointer
!
  call setoccshptr(nsfoc,nbsfoc,iocshptr,ibocshptr,(ndim.eq.2))
  call setoccptr(ncfoc,nsfoc,nbfoc,iocptr,ibocptr,(ndim.eq.2),lpocc)
  ncsfoc = ncfoc + nsfoc
!
!  Calculate parallel division of work :
!
!  Spatial - divide cells over processors
!  Non-spatial - divide atoms over processors
!
  call setatomnoded2
  call setatomdistribution('a')
  call setatomnodes(numat,nprocs,procid,lspatialok)
  call setatomnodesbo(numat,nprocs,procid,lspatialBOok)
#ifdef TRACE
  call trace_out('build_fc_super')
#endif
!
  return
  end
!
  subroutine unbuild_fc_supercell(nfcsuper)
!
!  Revert supercell for force constant calculation back to regular cell
!
!   7/23 Created from build_fc_supercell
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
  use configurations
  use control
  use current
  use molecule
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                   intent(in)    :: nfcsuper(3)   ! Force constant supercell dimensions
!
!  Local variables
!
#ifdef TRACE
  call trace_in('build_fc_super')
#endif
!
!  Check that system is periodic
!
  if (ndim.lt.1) then
    call outerror('fc_supercell cannot be specified for a cluster',0_i4)
    call stopnow('unbuild_fc_supercell')
  endif
!
!  Call setup to revert to the original configuration arrays
!
  call setup(.true.)
#ifdef TRACE
  call trace_out('unbuild_fc_super')
#endif
!
  return
  end
