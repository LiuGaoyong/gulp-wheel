  subroutine dcharged_gfnff(lprint,emat,maxe,lstrainin,lgrad2,cn,dcn)
!
!  Calculates the first derivative of the charge with respect to the coordinates of the atoms/strains for GFNFF.
!  NB: If only first derivatives are needed then dqdxyz/dqds are distributed across all nodes, but
!      if second derivatives are required then these arrays need to be global at present
!
!  NB: Might be problems for fixed charge case!!
!
!  Distributed memory parallel version.
!
!   1/22 Created from dcharged
!   1/22 lgrad2 added to arguments
!   9/23 Electric field added
!   9/23 dqdVxyz added
!   9/23 Strain derivatives in the presence of a field corrected
!  10/23 Electric field charge derivatives corrected for symmetry-adapted case
!
!  On entry:
!
!    emat      = inverse matrix calculated during EEM/QEq
!    lprint    = logical indicating whether printed output is wanted
!    lstrainin = if .true. then strain derivatives are calculated
!                if dimensionality is greater than 0
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
!  Julian Gale, CIC, Curtin University, October 2023
!
  use g_constants
  use control
  use current
  use derivatives
  use eemdata
  use field,          only : lfieldcfg, ntdfieldcfg
  use element
#ifdef MPI
  use general,        only : cutw, etaw
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_rad
#endif
  use iochannels
  use kspace
#ifdef MPI
  use m_strain,       only : strainddetds, straindet
#endif
  use m_strain,       only : gstrterms, real1strterm, cartstrterm
  use parallel
  use qmedata
  use shells
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  integer(i4), intent(in)                      :: maxe
  real(dp),    intent(inout)                   :: emat(maxe,*)
  logical,     intent(in)                      :: lprint
  logical,     intent(in)                      :: lstrainin
  logical,     intent(in)                      :: lgrad2
  real(dp),                     intent(in)     :: cn(numat)
  real(dp),                     intent(in)     :: dcn(numat)
#ifdef MPI
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: igroup
  integer(i4)                                  :: ii
  integer(i4)                                  :: iloc
  integer(i4)                                  :: ind
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4)                                  :: iresid
  integer(i4)                                  :: is
  integer(i4)                                  :: iv
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jloc
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: ml1
  integer(i4)                                  :: ml2
  integer(i4)                                  :: ml3
  integer(i4)                                  :: ns
  integer(i4)                                  :: status
!
  integer                                      :: MPIerror
  integer                                      :: idesd(9)
  integer                                      :: idesq(9)
  integer                                      :: idest(9)
  integer                                      :: ifails
  integer                                      :: ld
  integer                                      :: n3a
  integer                                      :: nb
  integer                                      :: ncs
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer                                      :: Request
  integer(i4), dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
!
  logical                                      :: lstrain
  real(dp)                                     :: accf2
  real(dp)                                     :: arg
  real(dp)                                     :: argtest
  real(dp)                                     :: cosa
  real(dp)                                     :: costrm
  real(dp)                                     :: cuts2
  real(dp)                                     :: d2gamr2
  real(dp)                                     :: darg1
  real(dp)                                     :: darg2
  real(dp),    dimension(:),   allocatable     :: dtmp
  real(dp),    dimension(:),   allocatable     :: dtmp2
  real(dp)                                     :: g_derf
  real(dp)                                     :: g_derfc
  real(dp)                                     :: derfc1
  real(dp)                                     :: derfc2
  real(dp)                                     :: derfez
  real(dp)                                     :: dexp1
  real(dp)                                     :: dexp2
  real(dp)                                     :: dexp3
  real(dp)                                     :: dexp4
  real(dp)                                     :: dexpz
  real(dp)                                     :: dqme(3)
  real(dp)                                     :: d2qme(6)
  real(dp)                                     :: dgam
  real(dp),    dimension(:,:), allocatable     :: dqtmp
  real(dp),    dimension(:,:), allocatable     :: dqs
  real(dp),    dimension(:,:), allocatable     :: dqx
  real(dp),    dimension(:,:), allocatable     :: dqy
  real(dp),    dimension(:,:), allocatable     :: dqz
  real(dp)                                     :: dtrm1
  real(dp)                                     :: dr2ds(6)
  real(dp)                                     :: d2r2dx2(3,3)
  real(dp)                                     :: d2r2ds2(6,6)
  real(dp)                                     :: d2r2dsdx(6,3)
  real(dp)                                     :: dxyzds(6,3)
  real(dp)                                     :: d2xyzdsdx(6,3,3)
  real(dp)                                     :: d2xyzds2(6,6,3)
  real(dp)                                     :: errfcn
  real(dp)                                     :: etaloc
  real(dp)                                     :: etaz
  real(dp)                                     :: etaz2
  real(dp)                                     :: etrm
  real(dp)                                     :: fieldx
  real(dp)                                     :: fieldy 
  real(dp)                                     :: fieldz
  real(dp)                                     :: gam
  real(dp)                                     :: Gmax
  real(dp)                                     :: kexperfc
  real(dp)                                     :: kvec
  real(dp)                                     :: qme
  real(dp)                                     :: qli
  real(dp)                                     :: qlii
  real(dp)                                     :: qlj
  real(dp)                                     :: qljj
  real(dp)                                     :: rconv
  real(dp)                                     :: rexp
  real(dp)                                     :: rkvec
  real(dp)                                     :: rl
  real(dp)                                     :: rmax
  real(dp)                                     :: rq
  real(dp)                                     :: rq2
  real(dp)                                     :: rr2
  real(dp)                                     :: rrr
  real(dp)                                     :: rv2
  real(dp)                                     :: rx
  real(dp)                                     :: rxi
  real(dp)                                     :: rxj
  real(dp)                                     :: rxk
  real(dp)                                     :: ry
  real(dp)                                     :: ry2
  real(dp)                                     :: ryi
  real(dp)                                     :: ryj
  real(dp)                                     :: ryk
  real(dp)                                     :: rz
  real(dp)                                     :: rz2
  real(dp)                                     :: rzi
  real(dp)                                     :: rzj
  real(dp)                                     :: rzk
  real(dp)                                     :: setaloc
  real(dp)                                     :: sina
  real(dp)                                     :: sineq
  real(dp)                                     :: smallestG
  real(dp)                                     :: strm
  real(dp)                                     :: strm1
  real(dp)                                     :: strm2
  real(dp)                                     :: sum
  real(dp)                                     :: trmi
  real(dp)                                     :: xci
  real(dp)                                     :: yci
  real(dp)                                     :: zci
  real(dp)                                     :: xd
  real(dp)                                     :: yd
  real(dp)                                     :: zd
  real(dp)                                     :: zetaij
  real(dp)                                     :: ztrm1
#ifdef TRACE
  call trace_in('dcharged_gfnff')
#endif
!
!  Set local strain flag
!
  lstrain = (lstrainin.and.nstrains.gt.0)
!
!  Check the memory for the charge derivatives
!
!  NB: If using parallel second derivatives then dqdxyz/dqds must be globalised
!      otherwise they can be distributed
!
  if (lgrad2) then
    if (numat.gt.maxd2qu) then
      maxd2qu = numat
      call changemaxd2q
    endif
  else
    if (natomsonnode.gt.maxd2qu) then
      maxd2qu = natomsonnode
      call changemaxd2q
    endif
  endif
!
  n3a = 3*numat
!
  if (n3a.gt.maxd2q) then
    maxd2q = n3a
    call changemaxd2q
  endif
!
!  Zero dq/dxyz and dq/ds
!
  do i = 1,natomsonnode
    do j = 1,n3a
      dqdxyz(j,i) = 0.0_dp
    enddo
    dqdVxyz(1:3,i) = 0.0_dp
  enddo
  if (ndim.gt.0) then
    do i = 1,natomsonnode
      do j = 1,nstrains
        dqds(j,i) = 0.0_dp
      enddo
    enddo
  endif
!
!  If electrostatics have been turned off then there is no point in proceeding.
!
  if (.not.lDoElectrostatics) then
#ifdef TRACE
    call trace_out('dcharged_gfnff')
#endif
    return
  endif
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
! 
!  Get electric field if needed                
! 
    call electricfieldparts(fieldx,fieldy,fieldz) 
  endif
!
  cuts2 = cuts*cuts
  rconv = 1.0_dp/autoangs
!
!  Set radius for Coulomb correction
!
  rq = gfnff_eeq_rad
  rq2 = rq**2
!
!  Allocate local memory
!
  allocate(dqtmp(max(6_i4,numat),max(2_i4,natomsonnode)),stat=status)
  if (status/=0) call outofmemory('dcharged_gfnff','dqtmp')
  allocate(dqx(numat,natomsonnode),stat=status)
  if (status/=0) call outofmemory('dcharged_gfnff','dqx')
  allocate(dqy(numat,natomsonnode),stat=status)
  if (status/=0) call outofmemory('dcharged_gfnff','dqy')
  allocate(dqz(numat,natomsonnode),stat=status)
  if (status/=0) call outofmemory('dcharged_gfnff','dqz')
!
!  Set up Blacs descriptors for matrices
!
  nb = nblocksize
  ifails = 0
  ncs = numat
  ld = maxe
  call descinit( idesd, ncs, ncs, nb, nb, 0, 0, iBlacsContext, ld, ifails )
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('dcharged_gfnff')
  endif
!
  ifails = 0
  ncs = numat
  ld = numat
  call descinit( idesq, ncs, ncs, nb, nb, 0, 0, iBlacsContext, ld, ifails )
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('dcharged_gfnff')
  endif
!
  ifails = 0
  ncs = numat
  ld = max(6,numat)
  call descinit( idest, ncs, ncs, nb, nb, 0, 0, iBlacsContext, ld, ifails )
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('dcharged_gfnff')
  endif
!
  if (ifails.ne.0) then
    call outerror('initialisation in descinit failed',0_i4)
    call stopnow('dcharged_gfnff')
  endif
!
  if (lstrain) then
    allocate(dqs(numat,nstrains),stat=status)
    if (status/=0) call outofmemory('dcharged_gfnff','dqs')
!
!  Zero temporary strain derivative arrays
!
    do i = 1,nstrains
      do j = 1,numat
        dqs(j,i) = 0.0_dp
      enddo
    enddo
  endif
  if (ndim.gt.1.or.lwolf) then
!******************
!  Periodic case  *
!******************
!
!  Define constants
!
    if (lwolf) then
      radmax = cutw
      etaloc = etaw*etaw
      setaloc = etaw
    elseif (lewald) then
      if (ndim.eq.2) then
        rpieta = 1.0_dp / sqrt(pi * eta)
        rhseta = 0.5_dp / seta
        accf2 = accf*accf
        argtest = sqrt(3.0+0.5*accf2) - sqrt(3.0)
        smallestG = min(kv(1,1),kv(2,2))
      endif
      radmax = accf/seta
      eta4 = 0.25_dp/eta
      etaloc = eta
      setaloc = seta
    else
      radmax = 0.0_dp
      eta4 = 0.25_dp
      etaloc = eta
      setaloc = seta
    endif
!
    rmax = max(radmax,rq)
    rmax2 = rmax*rmax
!
!  Estimate upper limits for looping
!
    if (ndim.eq.3) then
      rv2 = rv(1,1)**2 + rv(2,1)**2 + rv(3,1)**2
      rv2 = sqrt(rv2)
      ml1 = rmax/rv2 + 1
      rv2 = rv(1,2)**2 + rv(2,2)**2 + rv(3,2)**2
      rv2 = sqrt(rv2)
      ml2 = rmax/rv2 + 1
      rv2 = rv(1,3)**2 + rv(2,3)**2 + rv(3,3)**2
      rv2 = sqrt(rv2)
      ml3 = rmax/rv2 + 1
    elseif (ndim.eq.2) then
      rv2 = rv(1,1)**2 + rv(2,1)**2
      rv2 = sqrt(rv2)
      ml1 = rmax/rv2 + 1
      rv2 = rv(1,2)**2 + rv(2,2)**2
      rv2 = sqrt(rv2)
      ml2 = rmax/rv2 + 1
      ml3 = 0
    elseif (ndim.eq.1) then
      rv2 = rv(1,1)
      ml1 = rmax/rv2 + 1
      ml2 = 0
      ml3 = 0
    elseif (ndim.eq.0) then
      ml1 = 0
      ml2 = 0
      ml3 = 0
    endif
!**********************************
!  Reciprocal space contribution  *
!**********************************
    if (lnorecip.or.lwolf) goto 5
    if (lstrain) then
      call gstrterms(ndim,maxkvec,nkvec,xrk,yrk,zrk,dg2ds,d2g2dx2,d2g2ds2,.false.)
    endif
!
!  Start loop over atoms for coordinate differentiation
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      qli = qf(i)
      ind = 3*(i - 1)
!
!  Zero temporary derivative arrays
!
      do j = 1,numat
        dqx(j,iloc) = 0.0_dp
        dqy(j,iloc) = 0.0_dp
        dqz(j,iloc) = 0.0_dp
      enddo
      do j = 1,numat
        xd = xclat(j) - xci
        yd = yclat(j) - yci
        zd = zclat(j) - zci
        qlj = qf(j)
!
!  Evaluate derivative of reciprocal space elements of A
!
        if (ndim.eq.3) then
          do iv = 1,nkvec
            argc(iv) = xrk(iv)*xd + yrk(iv)*yd + zrk(iv)*zd
            sine(iv) = sin(argc(iv))*ktrm(iv)
!
            dqx(i,iloc) = dqx(i,iloc) + sine(iv)*xrk(iv)*qlj
            dqy(i,iloc) = dqy(i,iloc) + sine(iv)*yrk(iv)*qlj
            dqz(i,iloc) = dqz(i,iloc) + sine(iv)*zrk(iv)*qlj
!
            dqx(j,iloc) = dqx(j,iloc) + sine(iv)*xrk(iv)*qli
            dqy(j,iloc) = dqy(j,iloc) + sine(iv)*yrk(iv)*qli
            dqz(j,iloc) = dqz(j,iloc) + sine(iv)*zrk(iv)*qli
!
            if (lstrain) then
              costrm = cos(argc(iv))*angstoev*qlj
              strm1 = costrm*ktrms(iv)
              strm2 = costrm*ktrm(iv)
              if (lfinitestrain) then
                do is = 1,nstrains
                  dqs(i,is) = dqs(i,is) + strm1*dg2ds(iv,is) - strm2*strainddetds(is)*straindet
                enddo
              else
                do is = 1,nstrains
                  dqs(i,is) = dqs(i,is) + strm1*dg2ds(iv,is)
                enddo
                dqs(i,1) = dqs(i,1) - strm2
                dqs(i,2) = dqs(i,2) - strm2
                dqs(i,3) = dqs(i,3) - strm2
              endif
            endif
          enddo
        elseif (ndim.eq.2) then
!
!  First term - K vector independent
!
          etaz = seta*zd
          etaz2 = etaz*etaz
          derfez = g_derf(etaz)
          dexpz  = exp(-etaz2)
          etrm   = - vol4pi*(zd*derfez + dexpz*rpieta)*angstoev
          dtrm1  = - vol4pi*derfez
!
          dqz(i,iloc) = dqz(i,iloc) - dtrm1*qlj
          dqz(j,iloc) = dqz(j,iloc) - dtrm1*qli
!
          if (lstrain) then
            dqs(i,1) = dqs(i,1) - etrm*qlj
            dqs(i,2) = dqs(i,2) - etrm*qlj
          endif
!
!  Find local kvector cut-off
!
          if (abs(etaz).gt.argtest) then
            Gmax = abs(accf2/zd)
          else
            Gmax = sqrt(4.0*eta*(accf2-etaz2))
          endif
          if (Gmax.ge.smallestG) then
            do iv = 1,nkvec
              kvec = kmod(iv)
              if (kvec.le.Gmax) then
                arg = xrk(iv)*xd + yrk(iv)*yd
                sina = sin(arg)*ktrm(iv)
                cosa = cos(arg)*ktrm(iv)
                dexp1 = exp(kvec*zd)
                dexp2 = 1.0_dp/dexp1
                darg1 = kvec*rhseta + etaz
                darg2 = kvec*rhseta - etaz
                dexp3 = exp(-(darg1)**2)
                dexp4 = exp(-(darg2)**2)
                derfc1 = g_derfc(darg1)
                derfc2 = g_derfc(darg2)
                kexperfc = dexp1*derfc1 + dexp2*derfc2
                sineq = sina*kexperfc
                ztrm1 = (kvec*(dexp1*derfc1-dexp2*derfc2) - tweatpi*(dexp1*dexp3-dexp2*dexp4))
!
                dqx(i,iloc) = dqx(i,iloc) + sineq*xrk(iv)*qlj
                dqy(i,iloc) = dqy(i,iloc) + sineq*yrk(iv)*qlj
                dqz(i,iloc) = dqz(i,iloc) - cosa*ztrm1*qlj
!
                dqx(j,iloc) = dqx(j,iloc) + sineq*xrk(iv)*qli
                dqy(j,iloc) = dqy(j,iloc) + sineq*yrk(iv)*qli
                dqz(j,iloc) = dqz(j,iloc) - cosa*ztrm1*qli
!
                if (lstrain) then
                  costrm = cosa*angstoev*qlj
                  rkvec = 1.0_dp/kvec
                  strm = rkvec*(-rkvec*kexperfc + zd*(dexp1*derfc1-dexp2*derfc2) - &
                         rpieta*(dexp1*dexp3+dexp2*dexp4))
                  strm1 = costrm*strm
                  strm2 = costrm*kexperfc
                  if (lfinitestrain) then
                    do is = 1,nstrains
                      dqs(i,is) = dqs(i,is) + strm1*dg2ds(iv,is) - strm2*strainddetds(is)*straindet
                    enddo
                  else
                    do is = 1,nstrains
                      dqs(i,is) = dqs(i,is) + strm1*dg2ds(iv,is)
                    enddo
                    dqs(i,1) = dqs(i,1) - strm2
                    dqs(i,2) = dqs(i,2) - strm2
                  endif
                endif
              endif
            enddo
          endif
        endif
!
!  End inner atom loop
!
      enddo
!
!  End outer atom loop
!
    enddo
!
!  Global matrix multiplies
!
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqx,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqx(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqy,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqy(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqz,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqz(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
!
!  Multiply dqx/dqy/dqz by inverse matrix
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      do j = 1,numat
        ind = 3*(j - 1)
        dqdxyz(ind+1,iloc) = dqdxyz(ind+1,iloc) - dqx(j,iloc)*angstoev
        dqdxyz(ind+2,iloc) = dqdxyz(ind+2,iloc) - dqy(j,iloc)*angstoev
        dqdxyz(ind+3,iloc) = dqdxyz(ind+3,iloc) - dqz(j,iloc)*angstoev
      enddo
    enddo
5   continue
!*************************
!  Real space summation  *
!*************************
    if (lnoreal.and.lstrain) then
!
!  Globalise dqs
!
      do kl = 1,nstrains
        do j = 1,numat
          dqtmp(j,1) = dqs(j,kl)
        enddo
        call sumall(dqtmp(1,1),dqtmp(1,2),numat,"dcharged_gfnff","dqs")
        do j = 1,numat
          dqs(j,kl) = dqtmp(j,2)
        enddo
      enddo
!
!  Strain terms
!
      do jloc = 1,natomsonnode
        j = node2atom(jloc)
        do kl = 1,nstrains
          sum = 0.0_dp
          do k = 1,numat
            sum = sum + dqs(k,kl)*emat(k,jloc)
          enddo
          dqds(kl,jloc) = dqds(kl,jloc) - sum
        enddo
      enddo
      goto 135
    endif
!
!  Start loop over atoms for coordinate differentiation
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      ind = 3*(i-1)
!
!  Zero temporary derivative arrays
!
      do j = 1,numat
        dqx(j,iloc) = 0.0_dp
        dqy(j,iloc) = 0.0_dp
        dqz(j,iloc) = 0.0_dp
      enddo
      do j = 1,numat
        qljj = qf(j)
        qlj = qljj*occuf(j)
        rx = xclat(j) - xci
        ry = yclat(j) - yci
        rz = zclat(j) - zci
!
!  Loop over cell vectors
!
        rxi = rx - (ml1 + 1)*r1x
        ryi = ry - (ml1 + 1)*r1y
        rzi = rz - (ml1 + 1)*r1z
        do ii = -ml1,ml1
          rxi = rxi + r1x
          ryi = ryi + r1y
          rzi = rzi + r1z
          rxj = rxi - (ml2+1)*r2x
          ryj = ryi - (ml2+1)*r2y
          rzj = rzi - (ml2+1)*r2z
          do jj = -ml2,ml2
            rxj = rxj + r2x
            ryj = ryj + r2y
            rzj = rzj + r2z
            rxk = rxj - (ml3+1)*r3x
            ryk = ryj - (ml3+1)*r3y
            rzk = rzj - (ml3+1)*r3z
            do 120 kk = -ml3,ml3
              rxk = rxk + r3x
              ryk = ryk + r3y
              rzk = rzk + r3z
!
!  Calculate distance squared
!
              rr2 = rxk*rxk + ryk*ryk + rzk*rzk
!
!  Exclude distances outside maximum cutoff
!
              if (rr2.gt.rmax2) goto 120
!
!  Trap self term
!
              if (rr2.lt.1.0d-15) goto 120
              rl = sqrt(rr2)
              rrr = 1.0_dp/rl
              if (rr2.lt.cuts2) then
!
!  Core-shell interaction
!
                dtrm1 = rrr*rrr*rrr
              elseif (rr2.lt.rq2) then
                zetaij = 1.0_dp/sqrt(gfnff_eeq_alp(i)+gfnff_eeq_alp(j))
                call gfnff_gamma(zetaij,rl,gam,dgam,d2gamr2,.true.,.false.)
                dtrm1 = (rrr*rrr + dgam)*rrr
              else
                dtrm1 = 0.0_dp
              endif
!
!  Complementary error function
!
              errfcn = g_derfc(setaloc*rl)
              trmi = errfcn/rl
              rexp = tweatpi*exp(-etaloc*rr2)
              dtrm1 = dtrm1 - (trmi+rexp)*rrr*rrr
!
!  First derivatives of matrix elements
!
              dtrm1 = dtrm1*angstoev
!
              dqx(i,iloc) = dqx(i,iloc) - qljj*dtrm1*rxk
              dqy(i,iloc) = dqy(i,iloc) - qljj*dtrm1*ryk
              dqz(i,iloc) = dqz(i,iloc) - qljj*dtrm1*rzk
!
              dqx(j,iloc) = dqx(j,iloc) - qlii*dtrm1*rxk
              dqy(j,iloc) = dqy(j,iloc) - qlii*dtrm1*ryk
              dqz(j,iloc) = dqz(j,iloc) - qlii*dtrm1*rzk
!
              if (lstrain) then
!
!  Strain
!
                call real1strterm(ndim,rxk,ryk,rzk,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
                do is = 1,nstrains
                  ns = nstrptr(is)
                  dqs(i,is) = dqs(i,is) + qljj*dtrm1*dr2ds(ns)
                enddo
              endif
!
!  End of loops over lattice vectors
!
120         continue
          enddo
        enddo
!
!  End of loop over inner atom 
!
      enddo
!
!  Add contributions from coordination number to z
!
      call gfnff_drv_dcn_q(i,cn,dcn,dqx(1,iloc),dqy(1,iloc),dqz(1,iloc),dqs,numat)
!
!  End of loop over outer atom
!
    enddo
!
!  Global matrix multiplies
!
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqx,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqx(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqy,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqy(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqz,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqz(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
!
!  Multiply dqx/dqy/dqz by inverse matrix
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      do j = 1,numat
        ind = 3*(j - 1)
        dqdxyz(ind+1,iloc) = dqdxyz(ind+1,iloc) - dqx(j,iloc)
        dqdxyz(ind+2,iloc) = dqdxyz(ind+2,iloc) - dqy(j,iloc)
        dqdxyz(ind+3,iloc) = dqdxyz(ind+3,iloc) - dqz(j,iloc)
      enddo
    enddo
!
    if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!      
!  Electric field
!    
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        do j = 1,numat
          indj = 3*(j - 1)
          dqdxyz(indj+1,iloc) = dqdxyz(indj+1,iloc) + fieldx*emat(j,iloc)
          dqdxyz(indj+2,iloc) = dqdxyz(indj+2,iloc) + fieldy*emat(j,iloc)
          dqdxyz(indj+3,iloc) = dqdxyz(indj+3,iloc) + fieldz*emat(j,iloc)
        enddo
      enddo 
    endif
!
!  Compute electric field derivatives of charges
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      do j = 1,numat
        dqdVxyz(1,iloc) = dqdVxyz(1,iloc) + xclat(j)*emat(j,iloc)
        dqdVxyz(2,iloc) = dqdVxyz(2,iloc) + yclat(j)*emat(j,iloc)
        dqdVxyz(3,iloc) = dqdVxyz(3,iloc) + zclat(j)*emat(j,iloc)
      enddo
    enddo 
!
    if (lstrain) then
!
!  Globalise dqs
!
      do kl = 1,nstrains
        do j = 1,numat
          dqtmp(j,1) = dqs(j,kl)
        enddo
        call sumall(dqtmp(1,1),dqtmp(1,2),numat,"dcharged_gfnff","dqs")
        do j = 1,numat
          dqs(j,kl) = dqtmp(j,2) 
        enddo
      enddo
!
!  Strain terms
!
      do jloc = 1,natomsonnode
        j = node2atom(jloc)
        do kl = 1,nstrains
          sum = 0.0_dp
          do k = 1,numat
            sum = sum + dqs(k,kl)*emat(k,jloc)
          enddo
          dqds(kl,jloc) = dqds(kl,jloc) - sum
        enddo
      enddo
!
      if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  Electric field contribution to strains
!
        do j = 1,numat
          call cartstrterm(ndim,xalat(j),yalat(j),zalat(j),0.0_dp,0.0_dp,0.0_dp,dxyzds, &
                           d2xyzdsdx,d2xyzds2,.false.)
          do iloc = 1,natomsonnode
            do kl = 1,nstrains
              dqds(kl,iloc) = dqds(kl,iloc) + fieldx*dxyzds(kl,1)*emat(j,iloc) + &
                                              fieldy*dxyzds(kl,2)*emat(j,iloc) + &
                                              fieldz*dxyzds(kl,3)*emat(j,iloc)
            enddo
          enddo
        enddo
      endif
    endif
!**********************
!  End periodic case  *
!**********************
  else
!***********************
!  Cluster case / 1-D  *
!***********************
    if (lnoreal) goto 135
    if (ndim.eq.1) then
      call setmaxcell1D(maxloop(1))
      ml1 = maxloop(1)
    else
      ml1 = 0
      r1x = 0.0_dp
    endif
!
!  Start loop over cluster atoms
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      ind = 3*(i-1)
!
!  Zero temporary derivative arrays
!
      dqx(1:numat,iloc) = 0.0_dp
      dqy(1:numat,iloc) = 0.0_dp
      dqz(1:numat,iloc) = 0.0_dp
!
!  Loop over other atoms and build daij/d(alpha)
!
      do j = 1,numat
        if (i.ne.j.or.lstrain) then
          qljj = qf(j)
          qlj = qljj*occuf(j)
!
!  Find relative vector between atoms
!
          rx = xclat(j) - xci
          ry = yclat(j) - yci
          rz = zclat(j) - zci
!
!  Calculate Euler-Maclaurin correction to 1-D sum
!
          if (ndim.eq.1) then
            qme = 0.0_dp
            dqme(1) = 0.0_dp
            dqme(2) = 0.0_dp
            dqme(3) = 0.0_dp
            call qmatrix1D(rx,ry,rz,.true.,.false.,qme,dqme,d2qme)
            dqme(1) = dqme(1)*angstoev
            dqme(2) = dqme(2)*angstoev
            dqme(3) = dqme(3)*angstoev
!
            dqx(i,iloc) = dqx(i,iloc) - qljj*dqme(1)
            dqy(i,iloc) = dqy(i,iloc) - qljj*dqme(2)
            dqz(i,iloc) = dqz(i,iloc) - qljj*dqme(3)
!
            dqx(j,iloc) = dqx(j,iloc) - qlii*dqme(1)
            dqy(j,iloc) = dqy(j,iloc) - qlii*dqme(2)
            dqz(j,iloc) = dqz(j,iloc) - qlii*dqme(3)
          endif
!
!  Calculate distances for search
!
          rx = rx - (ml1+1)*r1x
          ry2 = ry*ry
          rz2 = rz*rz
!
!  Loop over lattice vectors
!
          do ii = -ml1,ml1
            rx = rx + r1x
            rr2 = rx*rx + ry2 + rz2
            rl = sqrt(rr2)
            if (rr2.gt.cuts2) then
              rrr = 1.0_dp/rl
!*********************
!  GFNFF-EEM scheme  *
!*********************
              if (rr2.lt.rq2) then
                zetaij = 1.0_dp/sqrt(gfnff_eeq_alp(i)+gfnff_eeq_alp(j))
                call gfnff_gamma(zetaij,rl,gam,dgam,d2gamr2,.true.,.false.)
                dtrm1 = dgam*rrr
              else
                dtrm1 = - rrr*rrr*rrr
              endif
!
!  First derivatives of matrix elements
!
              dtrm1 = dtrm1*angstoev
!
              dqx(i,iloc) = dqx(i,iloc) - qljj*dtrm1*rx
              dqy(i,iloc) = dqy(i,iloc) - qljj*dtrm1*ry
              dqz(i,iloc) = dqz(i,iloc) - qljj*dtrm1*rz
!
              dqx(j,iloc) = dqx(j,iloc) - qlii*dtrm1*rx
              dqy(j,iloc) = dqy(j,iloc) - qlii*dtrm1*ry
              dqz(j,iloc) = dqz(j,iloc) - qlii*dtrm1*rz
            endif
            if (lstrain) then
!             
!  Strain 
!             
              call real1strterm(ndim,rx,ry,rz,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
              dqs(i,1) = dqs(i,1) + qljj*dtrm1*dr2ds(1)
            endif
          enddo
!       
!  End of loop over lattice vectors
!
        endif
      enddo
!
!  Add contributions from coordination number to z
!
      call gfnff_drv_dcn_q(i,cn,dcn,dqx(1,iloc),dqy(1,iloc),dqz(1,iloc),dqs,numat)
!
!  End of outer loop over atoms
!
    enddo
!
!  Global matrix multiplies
!
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqx,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqx(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqy,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqy(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
    call pdgemm('t','n',numat,numat,numat,1.0d0,dqz,1,1,idesq,emat,1,1,idesd,0.0d0,dqtmp,1,1,idest)
    dqz(1:numat,1:natomsonnode) = dqtmp(1:numat,1:natomsonnode)
!
!  Multiply dqx/dqy/dqz by inverse matrix
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      do j = 1,numat
        ind  = 3*(j - 1)
        dqdxyz(ind+1,iloc) = dqdxyz(ind+1,iloc) - dqx(j,iloc)
        dqdxyz(ind+2,iloc) = dqdxyz(ind+2,iloc) - dqy(j,iloc)
        dqdxyz(ind+3,iloc) = dqdxyz(ind+3,iloc) - dqz(j,iloc)
      enddo
    enddo
!
    if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  Electric field
!
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        do j = 1,numat
          indj = 3*(j - 1)
          dqdxyz(indj+1,iloc) = dqdxyz(indj+1,iloc) + fieldx*emat(j,iloc)
          dqdxyz(indj+2,iloc) = dqdxyz(indj+2,iloc) + fieldy*emat(j,iloc)
          dqdxyz(indj+3,iloc) = dqdxyz(indj+3,iloc) + fieldz*emat(j,iloc)
        enddo
      enddo
    endif
!
!  Compute electric field derivatives of charges
!
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      do j = 1,numat
        dqdVxyz(1,iloc) = dqdVxyz(1,iloc) + xclat(j)*emat(j,iloc)
        dqdVxyz(2,iloc) = dqdVxyz(2,iloc) + yclat(j)*emat(j,iloc)
        dqdVxyz(3,iloc) = dqdVxyz(3,iloc) + zclat(j)*emat(j,iloc)
      enddo
    enddo
!
    if (lstrain) then
!
!  Globalise dqs
!
      do j = 1,numat
        dqtmp(j,1) = dqs(j,1)
      enddo
      call sumall(dqtmp(1,1),dqtmp(1,2),numat,"dcharged_gfnff","dqs")
      do j = 1,numat
        dqs(j,1) = dqtmp(j,2) 
      enddo
!
!  Strain terms
!
      do jloc = 1,natomsonnode
        j = node2atom(jloc)
        sum = 0.0_dp
        do k = 1,numat
          sum = sum + dqs(k,1)*emat(k,jloc)
        enddo
        dqds(1,jloc) = dqds(1,jloc) - sum
      enddo
    endif
!******************************
!  End of cluster / 1-D case  *
!******************************
  endif
135 continue
!***********************************************************************************
!  Enforce sum rules that total charge derivative equals zero for each coordinate  *
!***********************************************************************************
!
!  Sum elements and count number of active atoms
!
  do i = 1,numat
    dqx(i,1) = 0.0_dp
    dqy(i,1) = 0.0_dp
    dqz(i,1) = 0.0_dp
  enddo
  do iloc = 1,natomsonnode
    jx = - 2
    jy = - 1
    jz =   0
    do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
      dqx(j,1) = dqx(j,1) + dqdxyz(jx,iloc)
      dqy(j,1) = dqy(j,1) + dqdxyz(jy,iloc)
      dqz(j,1) = dqz(j,1) + dqdxyz(jz,iloc)
    enddo
  enddo
!
!  Globalise dqx / dqy / dqz 
!
  call sumall(dqx,dqtmp(1,1),numat,"dcharged_gfnff","dqx")
  dqx(1:numat,1) = dqtmp(1:numat,1)
  call sumall(dqy,dqtmp(1,1),numat,"dcharged_gfnff","dqx")
  dqy(1:numat,1) = dqtmp(1:numat,1)
  call sumall(dqz,dqtmp(1,1),numat,"dcharged_gfnff","dqx")
  dqz(1:numat,1) = dqtmp(1:numat,1)
!
!  Average error and subtract from active elements
!
  do i = 1,numat
    dqx(i,1) = - dqx(i,1)/dble(numat)
    dqy(i,1) = - dqy(i,1)/dble(numat)
    dqz(i,1) = - dqz(i,1)/dble(numat)
  enddo
  do iloc = 1,natomsonnode
    i = node2atom(iloc)
    jx = - 2
    jy = - 1
    jz =   0
    do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
      dqdxyz(jx,iloc) = dqdxyz(jx,iloc) + dqx(j,1)
      dqdxyz(jy,iloc) = dqdxyz(jy,iloc) + dqy(j,1)
      dqdxyz(jz,iloc) = dqdxyz(jz,iloc) + dqz(j,1)
    enddo
  enddo
!
!  Free local memory
!
  if (lstrain) then
    deallocate(dqs,stat=status)
    if (status/=0) call deallocate_error('dcharged_gfnff','dqs')
  endif
  deallocate(dqz,stat=status)
  if (status/=0) call deallocate_error('dcharged_gfnff','dqz')
  deallocate(dqy,stat=status)
  if (status/=0) call deallocate_error('dcharged_gfnff','dqy')
  deallocate(dqx,stat=status)
  if (status/=0) call deallocate_error('dcharged_gfnff','dqx')
  deallocate(dqtmp,stat=status)
  if (status/=0) call deallocate_error('dcharged_gfnff','dqtmp')
!
  if (lgrad2) then
!**********************************************************************
!  For 2nd derivative case globalise the derivatives and then output  *
!**********************************************************************
    if (lstrain) then
      ntmp = n3a + nstrains + 3
    else
      ntmp = n3a + 3
    endif
!
    allocate(dtmp(ntmp),stat=status)
    if (status/=0) call outofmemory('dcharged_gfnff','dtmp')
!
!
    do i = numat,1,-1
      if (procid.eq.atom2node(i)) then
        iloc = atom2local(i)
        dtmp(1:n3a) = dqdxyz(1:n3a,iloc)
        dtmp(n3a+1:n3a+3) = dqdVxyz(1:3,iloc)
        if (lstrain) then
          dtmp(n3a+4:n3a+3+nstrains) = dqds(1:nstrains,iloc)
        endif
      endif
      call sendall(dtmp,ntmp,atom2node(i),"dcharged_gfnff","dtmp")
      dqdxyz(1:n3a,i) = dtmp(1:n3a)
      dqdVxyz(1:3,i) = dtmp(n3a+1:n3a+3)
      if (lstrain) then
        dqds(1:nstrains,i) = dtmp(n3a+4:n3a+3+nstrains)
      endif
    enddo
!
    deallocate(dtmp,stat=status)
    if (status/=0) call deallocate_error('dcharged_gfnff','dtmp')
!
    if (lprint.and.ioproc) then
      if (ndim.eq.0) then
        write(ioout,'(/,''  First derivatives of charge distribution : (Angstroms**-1)'',/)')
      else
        write(ioout,'(/,''  First derivatives of charge distribution : '',/)')
        write(ioout,'(''  Strain :'',/)')
        write(ioout,'(''  Atom   '',6(i10))') (j,j=1,nstrains)
        do i = 1,numat
          write(ioout,'(i6,4x,6f10.6)') i,(dqds(j,i),j=1,nstrains)
        enddo
        write(ioout,'(/,''  Coordinate (Angstroms**-1) :'',/)')
      endif
      igroup = numat/6
      iresid = numat - igroup*6
      indi = 0
      if (igroup.gt.0) then
        do i = 1,igroup
          write(ioout,'(''  Atom   '',6(i10))') (indi+j,j=1,6)
          indj = 0
          do j = 1,numat
            write(ioout,'(i6,'' x '',4x,6f10.6)') j,(dqdxyz(indj+1,indi+k),k=1,6)
            write(ioout,'(i6,'' y '',4x,6f10.6)') j,(dqdxyz(indj+2,indi+k),k=1,6)
            write(ioout,'(i6,'' z '',4x,6f10.6)') j,(dqdxyz(indj+3,indi+k),k=1,6)
            indj = indj + 3
          enddo
          indi = indi + 6
          write(ioout,'(/)')
        enddo
      endif
      if (iresid.gt.0) then
        write(ioout,'(''  Atom   '',6(i10))') (indi+j,j=1,iresid)
        indj = 0
        do j = 1,numat
          write(ioout,'(i6,'' x '',4x,6f10.6)') j,(dqdxyz(indj+1,indi+k),k=1,iresid)
          write(ioout,'(i6,'' y '',4x,6f10.6)') j,(dqdxyz(indj+2,indi+k),k=1,iresid)
          write(ioout,'(i6,'' z '',4x,6f10.6)') j,(dqdxyz(indj+3,indi+k),k=1,iresid)
          indj = indj + 3
        enddo
      endif
      write(ioout,'(/)')
      write(ioout,'(''  First derivatives of charge distribution for electric field : (Angstroms/V)'',/)')
      write(ioout,'(''  Atom         X             Y            Z'')')
      do i = 1,numat
        write(ioout,'(i6,3(4x,f10.6))') i,(dqdVxyz(k,i),k=1,3)
      enddo
      write(ioout,'(/)')
    endif
  else
!*************************************************************
!  For 1st derivative case only communicate only for output  *
!*************************************************************
    if (lprint) then
      call mpbarrier
      if (ndim.eq.0) then
        if (ioproc) then
          write(ioout,'(/,''  First derivatives of charge distribution : (Angstroms**-1)'',/)')
        endif
      else
        if (ioproc) then
          write(ioout,'(/,''  First derivatives of charge distribution : '',/)')
          write(ioout,'(''  Strain :'',/)')
          write(ioout,'(''  Atom   '',6(i10))') (j,j=1,nstrains)
        endif
        if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
          ntmp = nstrains
          ntag = 1
          allocate(dtmp2(ntmp),stat=status)
          if (status/=0) call outofmemory('dcharged_gfnff','dtmp2')
          allocate(StatMPI(MPI_Status_Size),stat=status)
          if (status/=0) call outofmemory('dcharged_gfnff','StatMPI')
        endif
        call mpbarrier
        do i = 1,numat
          iloc = atom2local(i)
          if (lioproconly.and.atom2node(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = atom2node(i)
              call MPI_IRecv(dtmp2,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (iloc.gt.0) then
              dtmp2(1:nstrains) = dqds(1:nstrains,iloc)
!
!  Post send
!
              call MPI_ISend(dtmp2,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_World,Request,MPIerror)
            endif
            if (ioproc.or.iloc.gt.0) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
            if (ioproc) then
!
!  Write on I/O node
!
              write(ioout,'(i6,4x,6f10.6)') i,(dtmp2(j),j=1,nstrains)
            endif
          else
            if (iloc.gt.0) then
              write(ioout,'(i6,4x,6f10.6)') i,(dqds(j,atom2local(i)),j=1,nstrains)
            endif
          endif
          call mpbarrier
        enddo
        if (lioproconly) then
          deallocate(StatMPI,stat=status)
          if (status/=0) call deallocate_error('dcharged_gfnff','StatMPI')
          deallocate(dtmp2,stat=status)
          if (status/=0) call deallocate_error('dcharged_gfnff','dtmp2')
        endif
        if (ioproc) then
          write(ioout,'(/,''  Coordinate (Angstroms**-1) :'',/)')
        endif
      endif
!
      if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
        ntmp = n3a
        ntag = 1
        allocate(dtmp2(ntmp),stat=status)
        if (status/=0) call outofmemory('dcharged_gfnff','dtmp2')
        allocate(StatMPI(MPI_Status_Size),stat=status)
        if (status/=0) call outofmemory('dcharged_gfnff','StatMPI')
      endif
      call mpbarrier
      do i = 1,numat
        iloc = atom2local(i)
        if (lioproconly.and.atom2node(i).ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = atom2node(i)
            call MPI_IRecv(dtmp2,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (iloc.gt.0) then
            dtmp2(1:n3a) = dqdxyz(1:n3a,iloc)
!
!  Post send
!
            call MPI_ISend(dtmp2,ntmp,MPI_double_precision,0, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
          if (ioproc.or.iloc.gt.0) then
            call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          endif
          if (ioproc) then
!
!  Write on I/O node
!
            write(ioout,'(''  Atom   '',i10)') i
            indj = 0
            do j = 1,numat
              write(ioout,'(i6,'' x '',4x,f10.6)') j,dtmp2(indj+1)
              write(ioout,'(i6,'' y '',4x,f10.6)') j,dtmp2(indj+2)
              write(ioout,'(i6,'' z '',4x,f10.6)') j,dtmp2(indj+3)
              indj = indj + 3
            enddo
            write(ioout,'(/)')
          endif
        else
          if (iloc.gt.0) then
            write(ioout,'(''  Atom   '',i10)') i
            indj = 0
            do j = 1,numat
              write(ioout,'(i6,'' x '',4x,f10.6)') j,dqdxyz(indj+1,iloc)
              write(ioout,'(i6,'' y '',4x,f10.6)') j,dqdxyz(indj+2,iloc)
              write(ioout,'(i6,'' z '',4x,f10.6)') j,dqdxyz(indj+3,iloc)
              indj = indj + 3
            enddo
            write(ioout,'(/)')
          endif
        endif
        call mpbarrier
      enddo
      if (lioproconly) then
        deallocate(StatMPI,stat=status)
        if (status/=0) call deallocate_error('dcharged_gfnff','StatMPI')
        deallocate(dtmp2,stat=status)
        if (status/=0) call deallocate_error('dcharged_gfnff','dtmp2')
      endif
!
      call mpbarrier
      if (ioproc) then
        write(ioout,'(''  First derivatives of charge distribution for electric field : (Angstroms/V)'',/)')
        write(ioout,'(''  Atom         X             Y            Z'')')
      endif
!
      allocate(dtmp(numat),stat=status)
      if (status/=0) call outofmemory('dcharged_gfnff','dtmp')
      allocate(dtmp2(numat),stat=status)
      if (status/=0) call outofmemory('dcharged_gfnff','dtmp2')
!
      dtmp(1:numat) = 0.0_dp
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        dtmp(i) = dqdVxyz(1,iloc)
      enddo
!
      call sumall(dtmp,dtmp2,numat,"dcharged_gfnff","dqdVxyz")
!
      do i = 1,numat
        dqdVxyz(1,i) = dtmp2(i)
      enddo
!
      dtmp(1:numat) = 0.0_dp
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        dtmp(i) = dqdVxyz(2,iloc)
      enddo
!
      call sumall(dtmp,dtmp2,numat,"dcharged_gfnff","dqdVxyz")
!
      do i = 1,numat
        dqdVxyz(2,i) = dtmp2(i)
      enddo
!
      dtmp(1:numat) = 0.0_dp
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        dtmp(i) = dqdVxyz(3,iloc)
      enddo
!
      call sumall(dtmp,dtmp2,numat,"dcharged_gfnff","dqdVxyz")
!
      do i = 1,numat
        dqdVxyz(3,i) = dtmp2(i)
      enddo
!
      deallocate(dtmp2,stat=status)
      if (status/=0) call deallocate_error('dcharged_gfnff','dtmp2')
      deallocate(dtmp,stat=status)
      if (status/=0) call deallocate_error('dcharged_gfnff','dtmp')
!
      if (ioproc) then
        do i = 1,numat
          write(ioout,'(i6,3(4x,f10.6))') i,(dqdVxyz(k,i),k=1,3)
        enddo
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('dcharged_gfnff')
#endif
#else
  call outerror('dcharged_gfnff called when not compiled with MPI',0_i4)
  call stopnow('dcharged_gfnff')
#endif
!
  return
  end
