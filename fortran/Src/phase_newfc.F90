  subroutine phase_newfc(xk,yk,zk,nfccells,nfcsuper,fccells,fcderv2,maxd2fc,nd2fc)
!
!  Calculates the phased dynamical matrix for a k point from
!  the second derivatives stored for a supercell.
!  New version that directly uses the second derivative matrix for a supercell.
!
!  NB: Assumes that configuration has been reset to original one after supercell calculation.
!
!  Note that at present the reciprocal space and 1-D real space
!  Coulomb sums are computed using the standard algorithm.
!
!   7/23 Created from phase_fc
!   7/23 lbsmatcfg added
!   8/23 Breathing shells removed
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
!  Julian Gale, CIC, Curtin University, August 2023
!
  use g_constants
  use control
  use current
  use datatypes
  use derivatives
  use element
  use kspace
  use molecule
  use parallel
  use shells
  use symmetry
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),    intent(in)                   :: nfccells
  integer(i4),    intent(in)                   :: nfcsuper(3)
  integer(i4),    intent(in)                   :: nd2fc
  integer(i4),    intent(in)                   :: maxd2fc
  real(dp),       intent(in)                   :: fccells(3,nfccells)
  real(dp),       intent(in)                   :: xk
  real(dp),       intent(in)                   :: yk
  real(dp),       intent(in)                   :: zk
  real(dp),       intent(in)                   :: fcderv2(maxd2fc,nd2fc)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  integer(i4)                                  :: iim
  integer(i4)                                  :: iimn
  integer(i4)                                  :: iimx
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixs
  integer(i4)                                  :: iys
  integer(i4)                                  :: izs
  integer(i4)                                  :: j
  integer(i4)                                  :: jcell
  integer(i4)                                  :: jj
  integer(i4)                                  :: jjm
  integer(i4)                                  :: jjmn
  integer(i4)                                  :: jjmx
  integer(i4)                                  :: jx
  integer(i4)                                  :: jy
  integer(i4)                                  :: jz
  integer(i4)                                  :: jxs
  integer(i4)                                  :: jys
  integer(i4)                                  :: jzs
  integer(i4)                                  :: kk
  integer(i4)                                  :: kkm
  integer(i4)                                  :: kkmn
  integer(i4)                                  :: kkmx
  integer(i4)                                  :: maxlim
  real(dp)                                     :: cosk
  real(dp)                                     :: oneij
  real(dp)                                     :: sink
  real(dp)                                     :: kvf(3,3)
  real(dp)                                     :: r1xs
  real(dp)                                     :: r1ys
  real(dp)                                     :: r1zs
  real(dp)                                     :: r2xs
  real(dp)                                     :: r2ys
  real(dp)                                     :: r2zs
  real(dp)                                     :: r3xs
  real(dp)                                     :: r3ys
  real(dp)                                     :: r3zs
  real(dp)                                     :: r2
  real(dp)                                     :: r2min
  real(dp)                                     :: xcd
  real(dp)                                     :: ycd
  real(dp)                                     :: zcd
  real(dp)                                     :: xcdi
  real(dp)                                     :: ycdi
  real(dp)                                     :: zcdi
  real(dp)                                     :: xcdj
  real(dp)                                     :: ycdj
  real(dp)                                     :: zcdj
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xj
  real(dp)                                     :: yj
  real(dp)                                     :: zj
  real(dp)                                     :: xji
  real(dp)                                     :: yji
  real(dp)                                     :: zji
  real(dp)                                     :: xmin
  real(dp)                                     :: ymin
  real(dp)                                     :: zmin
  real(dp)                                     :: xkv
  real(dp)                                     :: ykv
  real(dp)                                     :: zkv
#ifdef TRACE
  call trace_in('phase_newfc')
#endif
!****************************
!  Zero second derivatives  *
!****************************
  maxlim = 3*numat
  if (maxlim.gt.maxd2u) then
    maxd2u = maxlim
    call changemaxd2
  endif
  if (maxlim.gt.maxd2) then
    maxd2 = maxlim
    call changemaxd2
  endif
  do i = 1,maxlim
    do j = 1,maxlim
      derv2(j,i) = 0.0_dp
      dervi(j,i) = 0.0_dp
    enddo
  enddo
!
!  If group velocities are to be computed then zero k derivatives of dynamical matrix
!
  if (lgroupvelocity) then
    derv2dk(1:3,1:maxlim,1:maxlim) = 0.0_dpc
  endif
!
!  Select appropriate K vectors
!
  if (lkfull.and.ndim.eq.3) then
    call kvector3Df(kvf)
  else
    kvf(1:3,1:3) = kv(1:3,1:3)
  endif
!***************************
!  Calculate phase factor  *
!***************************
  if (ndim.eq.3) then
    xkv = xk*kvf(1,1) + yk*kvf(1,2) + zk*kvf(1,3)
    ykv = xk*kvf(2,1) + yk*kvf(2,2) + zk*kvf(2,3)
    zkv = xk*kvf(3,1) + yk*kvf(3,2) + zk*kvf(3,3)
  elseif (ndim.eq.2) then
    xkv = xk*kvf(1,1) + yk*kvf(1,2)
    ykv = xk*kvf(2,1) + yk*kvf(2,2)
    zkv = 0.0_dp
  elseif (ndim.eq.1) then
    xkv = xk*kvf(1,1)
    ykv = 0.0_dp
    zkv = 0.0_dp
  endif
!   
!  Set up values for minimum image search
!   
  if (ndim.eq.3) then
    iimn = -1
    iimx =  1
    jjmn = -1
    jjmx =  1
    kkmn = -1
    kkmx =  1
  elseif (ndim.eq.2) then
    iimn = -1
    iimx =  1
    jjmn = -1
    jjmx =  1
    kkmn =  0
    kkmx =  0
  elseif (ndim.eq.1) then
    iimn = -1
    iimx =  1
    jjmn =  0
    jjmx =  0
    kkmn =  0
    kkmx =  0
  endif
!
  r1xs = rv(1,1)*dble(nfcsuper(1))
  r1ys = rv(2,1)*dble(nfcsuper(1))
  r1zs = rv(3,1)*dble(nfcsuper(1))
  r2xs = rv(1,2)*dble(nfcsuper(2))
  r2ys = rv(2,2)*dble(nfcsuper(2))
  r2zs = rv(3,2)*dble(nfcsuper(2))
  r3xs = rv(1,3)*dble(nfcsuper(3))
  r3ys = rv(2,3)*dble(nfcsuper(3))
  r3zs = rv(3,3)*dble(nfcsuper(3))
!*************************
!  Loop over first atom  *
!*************************
  ix = -2
  iy = -1
  iz =  0
  do i = 1,numat
    ix = ix + 3
    iy = iy + 3
    iz = iz + 3
!
    ixs = 3*nfccells*(i-1) + 1
    iys = ixs + 1
    izs = ixs + 2
!
    xi = xclat(i)
    yi = yclat(i)
    zi = zclat(i)
!
!  Loop over second atom
!
    jx = -2
    jy = -1
    jz =  0
    jxs = -2
    jys = -1
    jzs =  0
    do j = 1,numat
      jx = jx + 3
      jy = jy + 3
      jz = jz + 3
!
!  Loop over cells for j
!
      do jcell = 1,nfccells
        jxs = jxs + 3
        jys = jys + 3
        jzs = jzs + 3
!
        xj = xclat(j) + fccells(1,jcell)
        yj = yclat(j) + fccells(2,jcell)
        zj = zclat(j) + fccells(3,jcell)
!
        xji = xj - xi
        yji = yj - yi
        zji = zj - zi
!******************************************
!  Find minimum image based on supercell  *
!******************************************
        if (lra) then
          r2min = 1.0d10
          xcd = xji + (iimn-1)*r1xs
          do ii = iimn,iimx
            xcd = xcd + r1xs
            ycd = yji + (jjmn-1)*r2ys
            do jj = jjmn,jjmx
              ycd = ycd + r2ys
              zcd = zji + (kkmn-1)*r3zs
              do kk = kkmn,kkmx
                zcd = zcd + r3zs
                r2 = xcd*xcd + ycd*ycd + zcd*zcd
                if (r2.lt.r2min) then
                  r2min = r2
                  xmin = xcd
                  ymin = ycd
                  zmin = zcd
                  iim = ii
                  jjm = jj
                  kkm = kk
                endif
              enddo
            enddo
          enddo
          r2 = r2min
          xji = xmin
          yji = ymin
          zji = zmin
        else
          r2min = 1.0d10
          xcdi = xji + (iimn-1)*r1xs
          ycdi = yji + (iimn-1)*r1ys
          zcdi = zji + (iimn-1)*r1zs
          do ii = iimn,iimx
            xcdi = xcdi + r1xs
            ycdi = ycdi + r1ys
            zcdi = zcdi + r1zs
            xcdj = xcdi + (jjmn-1)*r2xs
            ycdj = ycdi + (jjmn-1)*r2ys
            zcdj = zcdi + (jjmn-1)*r2zs
            do jj = jjmn,jjmx
              xcdj = xcdj + r2xs
              ycdj = ycdj + r2ys
              zcdj = zcdj + r2zs
              xcd = xcdj + (kkmn-1)*r3xs
              ycd = ycdj + (kkmn-1)*r3ys
              zcd = zcdj + (kkmn-1)*r3zs
              do kk = kkmn,kkmx
                xcd = xcd + r3xs
                ycd = ycd + r3ys
                zcd = zcd + r3zs
                r2 = xcd*xcd + ycd*ycd + zcd*zcd
                if (r2.lt.r2min) then
                  r2min = r2
                  xmin = xcd
                  ymin = ycd
                  zmin = zcd
                  iim = ii
                  jjm = jj
                  kkm = kk
                endif
              enddo
            enddo
          enddo
          r2 = r2min
          xji = xmin
          yji = ymin
          zji = zmin
        endif
!
!  Compute pair-wise phase factor
!
        if (i.eq.j) then
          oneij = 1.0_dp
        else
          oneij = 0.0_dp
        endif
        cosk = xkv*xji + ykv*yji + zkv*zji
        sink = sin(cosk)
        cosk = cos(cosk) - oneij
!
!  Add phased contribution to the main arrays
!
        derv2(jx,ix) = derv2(jx,ix) + fcderv2(jxs,ixs)*cosk
        derv2(jy,ix) = derv2(jy,ix) + fcderv2(jys,ixs)*cosk
        derv2(jz,ix) = derv2(jz,ix) + fcderv2(jzs,ixs)*cosk
        derv2(jx,iy) = derv2(jx,iy) + fcderv2(jxs,iys)*cosk
        derv2(jy,iy) = derv2(jy,iy) + fcderv2(jys,iys)*cosk
        derv2(jz,iy) = derv2(jz,iy) + fcderv2(jzs,iys)*cosk
        derv2(jx,iz) = derv2(jx,iz) + fcderv2(jxs,izs)*cosk
        derv2(jy,iz) = derv2(jy,iz) + fcderv2(jys,izs)*cosk
        derv2(jz,iz) = derv2(jz,iz) + fcderv2(jzs,izs)*cosk
!
        dervi(jx,ix) = dervi(jx,ix) + fcderv2(jxs,ixs)*sink
        dervi(jy,ix) = dervi(jy,ix) + fcderv2(jys,ixs)*sink
        dervi(jz,ix) = dervi(jz,ix) + fcderv2(jzs,ixs)*sink
        dervi(jx,iy) = dervi(jx,iy) + fcderv2(jxs,iys)*sink
        dervi(jy,iy) = dervi(jy,iy) + fcderv2(jys,iys)*sink
        dervi(jz,iy) = dervi(jz,iy) + fcderv2(jzs,iys)*sink
        dervi(jx,iz) = dervi(jx,iz) + fcderv2(jxs,izs)*sink
        dervi(jy,iz) = dervi(jy,iz) + fcderv2(jys,izs)*sink
        dervi(jz,iz) = dervi(jz,iz) + fcderv2(jzs,izs)*sink
!
!  End of loop over cells for second atom
!
      enddo
!
!  End of loop over second atom
!
    enddo
!
!  End loop over first atom
!
  enddo
#ifdef TRACE
  call trace_out('phase_newfc')
#endif
!
  return
  end
