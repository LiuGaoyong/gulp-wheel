  subroutine dcharge_gfnff(lprint,emat,maxe,lstrainin,cn,dcn)
!
!  Calculates the first derivative of the charge with respect
!  to the coordinates of the atoms for GFNFF.
!
!   2/21 Created from dcharge
!   2/21 cut_cn no longer passed in as not needed
!   6/21 Use of xclat/yclat/zclat changed to xalat/yalat/zalat for benefit of MD
!   7/21 In non-MD case xclat/yclat/zclat are now used to allow for symmetry
!   7/21 Coordinates changed to xclat/yclat/zclat for all cases as shouldn't be 
!        called during MD
!  10/21 Modules rearranged
!   1/22 Field terms removed as not supported for GFNFF at present
!   3/22 dqs now always allocated
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
  use element
  use eemdata
  use field,          only : lfieldcfg, ntdfieldcfg
  use general,        only : cutw, etaw
  use gulp_gfnff,     only : gfnff_eeq_alp, gfnff_eeq_rad
  use iochannels
  use kspace
  use m_strain,       only : gstrterms, strainddetds, straindet
  use m_strain,       only : real1strterm, cartstrterm
  use parallel
  use qmedata
  use shells
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),                  intent(in)     :: maxe
  real(dp),                     intent(inout)  :: emat(maxe,*)
  logical,                      intent(in)     :: lprint
  logical,                      intent(in)     :: lstrainin
  real(dp),                     intent(in)     :: cn(numat)
  real(dp),                     intent(in)     :: dcn(numat)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ieem
  integer(i4)                                  :: igroup
  integer(i4)                                  :: ii
  integer(i4)                                  :: ind
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4)                                  :: iresid
  integer(i4)                                  :: is
  integer(i4)                                  :: iv
  integer(i4)                                  :: j
  integer(i4)                                  :: jeem
  integer(i4)                                  :: jj
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: kl
  integer(i4)                                  :: ml1
  integer(i4)                                  :: ml2
  integer(i4)                                  :: ml3
  integer(i4)                                  :: ndqs
  integer(i4)                                  :: ns
  integer(i4)                                  :: status
  logical                                      :: lstrain
  real(dp)                                     :: accf2
  real(dp)                                     :: arg
  real(dp)                                     :: argtest
  real(dp)                                     :: cosa
  real(dp)                                     :: costrm
  real(dp)                                     :: cuts2
  real(dp)                                     :: darg1
  real(dp)                                     :: darg2
  real(dp)                                     :: dgam
  real(dp)                                     :: d2gamr2
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
  real(dp),    dimension(:,:), allocatable     :: dqs
  real(dp),    dimension(:),   allocatable     :: dqx
  real(dp),    dimension(:),   allocatable     :: dqy
  real(dp),    dimension(:),   allocatable     :: dqz
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
  real(dp)                                     :: sumx
  real(dp)                                     :: sumy
  real(dp)                                     :: sumz
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
  call trace_in('dcharge_gfnff')
#endif
!
!  Set local strain flag
!
  lstrain = (lstrainin.and.nstrains.gt.0)
!
!  Check the memory for the charge derivatives
!
  if (numat.gt.maxd2qu) then
    maxd2qu = numat
    call changemaxd2q
  endif
  if (3*numat.gt.maxd2q) then
    maxd2q = 3*numat
    call changemaxd2q
  endif
!
!  Zero dq/dxyz and dq/ds
!
  do i = 1,numat
    do j = 1,3*numat
      dqdxyz(j,i) = 0.0_dp
    enddo
    dqdVxyz(1:3,i) = 0.0_dp
  enddo
  if (ndim.gt.0) then
    do i = 1,numat
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
    call trace_out('dcharge_gfnff')
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
  allocate(dqx(numat),stat=status)
  if (status/=0) call outofmemory('dcharge_gfnff','dqx')
  allocate(dqy(numat),stat=status)
  if (status/=0) call outofmemory('dcharge_gfnff','dqy')
  allocate(dqz(numat),stat=status)
  if (status/=0) call outofmemory('dcharge_gfnff','dqz')
!
  if (lstrain) then
    ndqs = numat
    allocate(dqs(numat,nstrains),stat=status)
    if (status/=0) call outofmemory('dcharge_gfnff','dqs')
!
!  Zero temporary strain derivative arrays
!
    do i = 1,nstrains
      do j = 1,numat
        dqs(j,i) = 0.0_dp
      enddo
    enddo
  else
    ndqs = 1
    allocate(dqs(ndqs,1),stat=status)
    if (status/=0) call outofmemory('dcharge_gfnff','dqs')
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
    do ieem = 1,neem
      i = neemptr(ieem)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      qli = qf(i)
      ind = 3*(i - 1)
!
!  Zero temporary derivative arrays
!
      do j = 1,neem
        dqx(j) = 0.0_dp
        dqy(j) = 0.0_dp
        dqz(j) = 0.0_dp
      enddo
      do jeem = 1,numat
        j = neemptr(jeem)
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
            dqx(ieem) = dqx(ieem) + sine(iv)*xrk(iv)*qlj
            dqy(ieem) = dqy(ieem) + sine(iv)*yrk(iv)*qlj
            dqz(ieem) = dqz(ieem) + sine(iv)*zrk(iv)*qlj
!
            dqx(jeem) = dqx(jeem) + sine(iv)*xrk(iv)*qli
            dqy(jeem) = dqy(jeem) + sine(iv)*yrk(iv)*qli
            dqz(jeem) = dqz(jeem) + sine(iv)*zrk(iv)*qli
!
            if (lstrain) then
              costrm = cos(argc(iv))*angstoev*qlj
              strm1 = costrm*ktrms(iv)
              strm2 = costrm*ktrm(iv)
              if (lfinitestrain) then
                do is = 1,nstrains
                  dqs(ieem,is) = dqs(ieem,is) + strm1*dg2ds(iv,is) - strm2*strainddetds(is)*straindet
                enddo
              else
                do is = 1,nstrains
                  dqs(ieem,is) = dqs(ieem,is) + strm1*dg2ds(iv,is)
                enddo
                dqs(ieem,1) = dqs(ieem,1) - strm2
                dqs(ieem,2) = dqs(ieem,2) - strm2
                dqs(ieem,3) = dqs(ieem,3) - strm2
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
          dqz(ieem) = dqz(ieem) - dtrm1*qlj
!
          dqz(jeem) = dqz(jeem) - dtrm1*qli
!
          if (lstrain) then
            dqs(ieem,1) = dqs(ieem,1) - etrm*qlj
            dqs(ieem,2) = dqs(ieem,2) - etrm*qlj
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
                dqx(ieem) = dqx(ieem) + sineq*xrk(iv)*qlj
                dqy(ieem) = dqy(ieem) + sineq*yrk(iv)*qlj
                dqz(ieem) = dqz(ieem) - cosa*ztrm1*qlj
!
                dqx(jeem) = dqx(jeem) + sineq*xrk(iv)*qli
                dqy(jeem) = dqy(jeem) + sineq*yrk(iv)*qli
                dqz(jeem) = dqz(jeem) - cosa*ztrm1*qli
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
                      dqs(ieem,is) = dqs(ieem,is) + strm1*dg2ds(iv,is) - strm2*strainddetds(is)*straindet
                    enddo
                  else
                    do is = 1,nstrains
                      dqs(ieem,is) = dqs(ieem,is) + strm1*dg2ds(iv,is) 
                    enddo
                    dqs(ieem,1) = dqs(ieem,1) - strm2
                    dqs(ieem,2) = dqs(ieem,2) - strm2
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
!  Multiply dqx/dqy/dqz by inverse matrix
!
      do jeem = 1,neem
        j = neemptr(jeem)
        sumx = 0.0_dp
        sumy = 0.0_dp
        sumz = 0.0_dp
        do k = 1,neem
          sumx = sumx + dqx(k)*emat(k,jeem)
          sumy = sumy + dqy(k)*emat(k,jeem)
          sumz = sumz + dqz(k)*emat(k,jeem)
        enddo
        sumx = sumx*angstoev
        sumy = sumy*angstoev
        sumz = sumz*angstoev
        dqdxyz(ind+1,j) = dqdxyz(ind+1,j) - sumx
        dqdxyz(ind+2,j) = dqdxyz(ind+2,j) - sumy
        dqdxyz(ind+3,j) = dqdxyz(ind+3,j) - sumz
      enddo
    enddo
5   continue
!*************************
!  Real space summation  *
!*************************
    if (lnoreal.and.lstrain) then
!
!  Strain terms
!
      do jeem = 1,neem
        j = neemptr(jeem)
        do kl = 1,nstrains
          sum = 0.0_dp
          do k = 1,neem
            sum = sum + dqs(k,kl)*emat(k,jeem)
          enddo
          dqds(kl,j) = dqds(kl,j) - sum
        enddo
      enddo
      goto 135
    endif
!
!  Start loop over atoms for coordinate differentiation
!
    do ieem = 1,neem
      i = neemptr(ieem)
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      ind = 3*(i-1)
!
!  Zero temporary derivative arrays
!
      do j = 1,neem
        dqx(j) = 0.0_dp
        dqy(j) = 0.0_dp
        dqz(j) = 0.0_dp
      enddo
      jeem = 0
      do jeem = 1,neem
        j = neemptr(jeem)
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
              dqx(ieem) = dqx(ieem) - qljj*dtrm1*rxk
              dqy(ieem) = dqy(ieem) - qljj*dtrm1*ryk
              dqz(ieem) = dqz(ieem) - qljj*dtrm1*rzk
!
              dqx(jeem) = dqx(jeem) - qlii*dtrm1*rxk
              dqy(jeem) = dqy(jeem) - qlii*dtrm1*ryk
              dqz(jeem) = dqz(jeem) - qlii*dtrm1*rzk
!
              if (lstrain) then
!
!  Strain
!
                call real1strterm(ndim,rxk,ryk,rzk,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
                do is = 1,nstrains
                  ns = nstrptr(is)
                  dqs(ieem,is) = dqs(ieem,is) + qljj*dtrm1*dr2ds(ns)
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
      call gfnff_drv_dcn_q(i,cn,dcn,dqx,dqy,dqz,dqs,numat)
!
!  Multiply dqx/dqy/dqz by inverse matrix
!
      do jeem = 1,neem
        j = neemptr(jeem)
        sumx = 0.0_dp
        sumy = 0.0_dp
        sumz = 0.0_dp
        do k = 1,neem
          sumx = sumx + dqx(k)*emat(k,jeem)
          sumy = sumy + dqy(k)*emat(k,jeem)
          sumz = sumz + dqz(k)*emat(k,jeem)
        enddo
        dqdxyz(ind+1,j) = dqdxyz(ind+1,j) - sumx
        dqdxyz(ind+2,j) = dqdxyz(ind+2,j) - sumy
        dqdxyz(ind+3,j) = dqdxyz(ind+3,j) - sumz
      enddo
!
      if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!       
!  Electric field contribution
!   
        do jeem = 1,neem
          j = neemptr(jeem)
          indj = 3*(j - 1)
          dqdxyz(indj+1,i) = dqdxyz(indj+1,i) + fieldx*emat(jeem,ieem)
          dqdxyz(indj+2,i) = dqdxyz(indj+2,i) + fieldy*emat(jeem,ieem)
          dqdxyz(indj+3,i) = dqdxyz(indj+3,i) + fieldz*emat(jeem,ieem)
        enddo
      endif
!
!  Compute electric field derivatives of charges
!
      do jeem = 1,neem
        j = neemptr(jeem)
        dqdVxyz(1,i) = dqdVxyz(1,i) + xclat(j)*emat(jeem,ieem)
        dqdVxyz(2,i) = dqdVxyz(2,i) + yclat(j)*emat(jeem,ieem)
        dqdVxyz(3,i) = dqdVxyz(3,i) + zclat(j)*emat(jeem,ieem)
      enddo
!
!  End loop over i
!
    enddo
!
    if (lstrain) then
!
!  Strain terms
!
      do jeem = 1,neem
        j = neemptr(jeem)
        do kl = 1,nstrains
          sum = 0.0_dp
          do k = 1,neem
            sum = sum + dqs(k,kl)*emat(k,jeem)
          enddo
          dqds(kl,j) = dqds(kl,j) - sum
        enddo
      enddo
!
      if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  Electric field contribution to strains
!
        do jeem = 1,neem
          j = neemptr(jeem)
          call cartstrterm(ndim,xalat(j),yalat(j),zalat(j),0.0_dp,0.0_dp,0.0_dp,dxyzds, &
                           d2xyzdsdx,d2xyzds2,.false.)
          do ieem = 1,neem
            i = neemptr(ieem)
            do kl = 1,nstrains
              dqds(kl,i) = dqds(kl,i) + fieldx*dxyzds(kl,1)*emat(jeem,ieem) + &
                                        fieldy*dxyzds(kl,2)*emat(jeem,ieem) + &
                                        fieldz*dxyzds(kl,3)*emat(jeem,ieem)
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
    do ieem = 1,neem
      i = neemptr(ieem)
      qlii = qf(i)
      xci = xclat(i)
      yci = yclat(i)
      zci = zclat(i)
      ind = 3*(i-1)
!
!  Zero temporary derivative arrays
!
      dqx(1:numat) = 0.0_dp
      dqy(1:numat) = 0.0_dp
      dqz(1:numat) = 0.0_dp
!
!  Loop over other atoms and build daij/d(alpha)
!
      do jeem = 1,neem
        if (ieem.ne.jeem.or.lstrain) then
          j = neemptr(jeem)
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
            dqx(ieem) = dqx(ieem) - qljj*dqme(1)
            dqy(ieem) = dqy(ieem) - qljj*dqme(2)
            dqz(ieem) = dqz(ieem) - qljj*dqme(3)
!
            dqx(jeem) = dqx(jeem) - qlii*dqme(1)
            dqy(jeem) = dqy(jeem) - qlii*dqme(2)
            dqz(jeem) = dqz(jeem) - qlii*dqme(3)
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
            if (rr2.gt.cuts2) then
              rl = sqrt(rr2)
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
              dqx(ieem) = dqx(ieem) - qljj*dtrm1*rx
              dqy(ieem) = dqy(ieem) - qljj*dtrm1*ry
              dqz(ieem) = dqz(ieem) - qljj*dtrm1*rz
!
              dqx(jeem) = dqx(jeem) - qlii*dtrm1*rx
              dqy(jeem) = dqy(jeem) - qlii*dtrm1*ry
              dqz(jeem) = dqz(jeem) - qlii*dtrm1*rz
            endif
            if (lstrain) then
!             
!  Strain 
!             
              call real1strterm(ndim,rx,ry,rz,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
              dqs(ieem,1) = dqs(ieem,1) + qljj*dtrm1*dr2ds(1)
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
      call gfnff_drv_dcn_q(i,cn,dcn,dqx,dqy,dqz,dqs,numat)
!
!  Multiply dqx/dqy/dqz by inverse matrix 
!
      do jeem = 1,neem
        j = neemptr(jeem)
        sumx = 0.0_dp
        sumy = 0.0_dp
        sumz = 0.0_dp
        do k = 1,neem
          sumx = sumx + dqx(k)*emat(k,jeem)
          sumy = sumy + dqy(k)*emat(k,jeem)
          sumz = sumz + dqz(k)*emat(k,jeem)
        enddo
        dqdxyz(ind+1,j) = dqdxyz(ind+1,j) - sumx
        dqdxyz(ind+2,j) = dqdxyz(ind+2,j) - sumy
        dqdxyz(ind+3,j) = dqdxyz(ind+3,j) - sumz
      enddo
!
      if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  Electric field contribution
!  
        do jeem = 1,neem
          j = neemptr(jeem)
          indj = 3*(j - 1)
          dqdxyz(indj+1,i) = dqdxyz(indj+1,i) + fieldx*emat(jeem,ieem)
          dqdxyz(indj+2,i) = dqdxyz(indj+2,i) + fieldy*emat(jeem,ieem)
          dqdxyz(indj+3,i) = dqdxyz(indj+3,i) + fieldz*emat(jeem,ieem)
        enddo
      endif
!
!  Compute electric field derivatives of charges
!
      do jeem = 1,neem
        j = neemptr(jeem)
        dqdVxyz(1,i) = dqdVxyz(1,i) + xclat(j)*emat(jeem,ieem)
        dqdVxyz(2,i) = dqdVxyz(2,i) + yclat(j)*emat(jeem,ieem)
        dqdVxyz(3,i) = dqdVxyz(3,i) + zclat(j)*emat(jeem,ieem)
      enddo
    enddo
    if (lstrain) then
!         
!  Strain terms
!           
      do jeem = 1,neem
        j = neemptr(jeem)
        do kl = 1,nstrains
          sum = 0.0_dp
          do k = 1,neem
            sum = sum + dqs(k,kl)*emat(k,jeem)
          enddo
          dqds(kl,j) = dqds(kl,j) - sum
        enddo
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
    dqx(i) = 0.0_dp
    dqy(i) = 0.0_dp
    dqz(i) = 0.0_dp
  enddo
  do i = 1,numat
    jx = - 2
    jy = - 1
    jz =   0
    do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
      dqx(j) = dqx(j) + dqdxyz(jx,i)
      dqy(j) = dqy(j) + dqdxyz(jy,i)
      dqz(j) = dqz(j) + dqdxyz(jz,i)
    enddo
  enddo
!
!  Average error and subtract from active elements
!
  if (neem.gt.0) then
    do i = 1,numat
      dqx(i) = - dqx(i)/dble(neem)
      dqy(i) = - dqy(i)/dble(neem)
      dqz(i) = - dqz(i)/dble(neem)
    enddo
    do ieem = 1,neem
      i = neemptr(ieem)
      jx = - 2
      jy = - 1
      jz =   0
      do j = 1,numat
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3
        dqdxyz(jx,i) = dqdxyz(jx,i) + dqx(j)
        dqdxyz(jy,i) = dqdxyz(jy,i) + dqy(j)
        dqdxyz(jz,i) = dqdxyz(jz,i) + dqz(j)
      enddo
    enddo
  endif
!
!  Free local memory
!
  deallocate(dqs,stat=status)
  if (status/=0) call deallocate_error('dcharge_gfnff','dqs')
  deallocate(dqz,stat=status)
  if (status/=0) call deallocate_error('dcharge_gfnff','dqz')
  deallocate(dqy,stat=status)
  if (status/=0) call deallocate_error('dcharge_gfnff','dqy')
  deallocate(dqx,stat=status)
  if (status/=0) call deallocate_error('dcharge_gfnff','dqx')
!
  if (lprint.and.ioproc) then
    if (ndim.eq.0) then
      write(ioout,'(/,''  First derivatives of GFNFF charge distribution : (Angstroms**-1)'',/)')
    else
      write(ioout,'(/,''  First derivatives of GFNFF charge distribution : '',/)')
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
    write(ioout,'(''  First derivatives of GFNFF charge distribution for electric field : (Angstroms/V)'',/)')
    write(ioout,'(''  Atom         X             Y            Z'')')
    do i = 1,numat
      write(ioout,'(i6,3(4x,f10.6))') i,(dqdVxyz(k,i),k=1,3)
    enddo
    write(ioout,'(/)')
  endif
#ifdef TRACE
  call trace_out('dcharge_gfnff')
#endif
!
  return
  end
