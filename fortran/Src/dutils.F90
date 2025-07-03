!*********************************************
!  Utility routines for derivative handling  *
!*********************************************
!
!  Structure of d1 is as follows:
!
!  1 -> nneigh(i) => derivatives between i and neighbours of i
!  nneigh(i) + 1 -> nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 => derivatives between neighbours of i and each other
!
!  if (ltorderv) is true:
!
!  nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + 1 -> nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh => i to neighbours of neighbours of i
!  nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh + 1 -> 
!    nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + nneigh(i)*(mneigh+1)*mneigh => neighbours of i to neighbours of neighbours.
!
!  if (lextended) is true:
!
!  nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + nneigh(i)*(mneigh+1)*mneigh + 1 -> 
!    nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + nneigh(i)*(mneigh+1)*mneigh + (nj-1)*mneigh*mneigh + (nl-1)*mneigh + nm => nearest neighbours in third shell
!
  subroutine d2addpdmsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nREBOatomRptr,d1,d2,xkv,ykv,zkv, &
                        nonzero,nonzeroptr,ltorderv,lextended)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian phased second derivative arrays.
!
!  Distributed memory parallel version. 
!  Sparse storage of derivatives version.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  xkv             = x component of K vector
!  ykv             = y component of K vector
!  zkv             = z component of K vector
!  nonzero         = number of nonzero elements in d1
!  nonzeroptr(*)   = pointer to nonzero elements in d1
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!
!   8/23 Created from d2addpdm
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
  use brennerdata
  use control,        only : lgroupvelocity
  use derivatives
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nonzero
  integer(i4), intent(in)             :: nonzeroptr(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  real(dp),    intent(in)             :: xkv
  real(dp),    intent(in)             :: ykv
  real(dp),    intent(in)             :: zkv
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ind0
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: ixf,iyf,izf
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: jxf,jyf,jzf
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: kxf,kyf,kzf
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: lxf,lyf,lzf
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nf1
  integer(i4)                         :: nf2
  integer(i4)                         :: nji
  integer(i4)                         :: nkj
  integer(i4)                         :: nlk
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: ntmp
  complex(dpc)                        :: cdk(3)
  real(dp)                            :: cosik
  real(dp)                            :: cosil
  real(dp)                            :: cosjk
  real(dp)                            :: cosjl
  real(dp)                            :: d2k
  real(dp)                            :: d2ks
  real(dp)                            :: oneik
  real(dp)                            :: oneil
  real(dp)                            :: onejk
  real(dp)                            :: onejl
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: sinik
  real(dp)                            :: sinil
  real(dp)                            :: sinjk
  real(dp)                            :: sinjl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i
  real(dp)                            :: z2i
  logical                             :: lilocal
  logical                             :: lijlocal
  logical                             :: ljlocal
  logical                             :: lklocal
  logical                             :: lkllocal
  logical                             :: lllocal
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lsameijkl
#ifdef TRACE
  call trace_in('d2addpdmsp')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndisttor = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
    if (lextended) then
      ndist = ndisttor + nneigh(nri)*mneigh*mneigh
    else
      ndist = ndisttor
    endif
  else
    ndist = nneigh2
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = nneigh2 + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
!  Loop over pairs of interatomic distances
!
  do nf1 = 1,nonzero
    n1 = nonzeroptr(nf1)
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!     
!  Check that neighbour numbers are valid
!         
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))       
!
      if (liok.and.ljok) then
!         
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    elseif (n1.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!         
!  Check that atoms are valid
!       
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
!  
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    else
!
!  Extended derivatives
!
      ntmp = n1 - ndisttor
      nji = (ntmp - 1)/(mneigh*mneigh) + 1
      ntmp = ntmp - (nji - 1)*mneigh*mneigh
      nkj = (ntmp - 1)/mneigh + 1
      ntmp = ntmp - (nkj - 1)*mneigh
      nlk = ntmp
!
      nrj = neighno(nji,nri)
      if (nrj.eq.0) cycle
      if (nkj.gt.nneigh(nrj)) cycle
      nrk = neighno(nkj,nrj)
      if (nrk.eq.0) cycle
      if (nlk.gt.nneigh(nrk)) cycle
      n1i = nrk
      n1j = neighno(nlk,nrk)
      if (n1j.eq.0) cycle
! 
!  Check that atoms are valid
!
      liok = (nkj.le.nneigh(nrj))
      ljok = (nlk.le.nneigh(nrk))
!
      if (liok.and.ljok) then
        x2i = xneigh(nji,nri) + xneigh(nkj,nrj)
        y2i = yneigh(nji,nri) + yneigh(nkj,nrj)
        z2i = zneigh(nji,nri) + zneigh(nkj,nrj)
!
        xd1 = xneigh(nlk,nrk)
        yd1 = yneigh(nlk,nrk)
        zd1 = zneigh(nlk,nrk)
      endif
    endif
!  
!  Skip the following section as it is irrelevant if i and j are not valid
!           
    if (liok.and.ljok) then
!
      lilocal = (atom2local(n1i).gt.0)
      ljlocal = (atom2local(n1j).gt.0)
      lijlocal = (lilocal.or.ljlocal)
!
!  Set second derivative location pointers
!
      ixf = 3*(n1i - 1) + 1
      iyf = ixf + 1
      izf = iyf + 1
!
      if (lilocal) then
        ix = 3*(atom2local(n1i) - 1) + 1
        iy = ix + 1
        iz = iy + 1
      endif
!
      jxf = 3*(n1j - 1) + 1
      jyf = jxf + 1
      jzf = jyf + 1
!
      if (ljlocal) then
        jx = 3*(atom2local(n1j) - 1) + 1
        jy = jx + 1
        jz = jy + 1
      endif
!
      ind0 = nf1*(nf1-1)/2
      do nf2 = 1,nf1
        n2 = nonzeroptr(nf2)
        ind = ind0 + nf2
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(nf2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!  
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!               
!  Calculate vector between atoms
!             
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          elseif (n2.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!               
!  Check that atoms are valid
!             
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(neighno(nji,nri))) 
!               
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          else
!
!  Extended derivatives
!
            ntmp = n2 - ndisttor
            nji = (ntmp - 1)/(mneigh*mneigh) + 1
            ntmp = ntmp - (nji - 1)*mneigh*mneigh
            nkj = (ntmp - 1)/mneigh + 1
            ntmp = ntmp - (nkj - 1)*mneigh
            nlk = ntmp
!
            nrj = neighno(nji,nri)
            if (nrj.eq.0) cycle
            if (nkj.gt.nneigh(nrj)) cycle
            nrk = neighno(nkj,nrj)
            if (nrk.eq.0) cycle
            if (nlk.gt.nneigh(nrk)) cycle
            n2k = nrk
            n2l = neighno(nlk,nrk)
            if (n2l.eq.0) cycle
!
!  Check that atoms are valid
!
            lkok = (nkj.le.nneigh(nrj))
            llok = (nlk.le.nneigh(nrk))
!
            if (liok.and.ljok) then
              xik = xneigh(nji,nri) + xneigh(nkj,nrj) - x2i
              yik = yneigh(nji,nri) + yneigh(nkj,nrj) - y2i
              zik = zneigh(nji,nri) + zneigh(nkj,nrj) - z2i
!
              xd2 = xneigh(nlk,nrk)
              yd2 = yneigh(nlk,nrk)
              zd2 = zneigh(nlk,nrk)
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
            lklocal = (atom2local(n2k).gt.0)
            lllocal = (atom2local(n2l).gt.0)
            lkllocal = (lklocal.or.lllocal)
!
!  Check whether any of i, j, k or l are local before going further
!
            if (lijlocal.or.lkllocal) then
!
!  Complete remaining vectors
!
              xil = xik + xd2
              yil = yik + yd2
              zil = zik + zd2
              xjk = xik - xd1
              yjk = yik - yd1
              zjk = zik - zd1
              xjl = xjk + xd2
              yjl = yjk + yd2
              zjl = zjk + zd2
!
!  Compute distances squared
!
              r2ik = xik*xik + yik*yik + zik*zik
              r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
              kxf = 3*(n2k - 1) + 1
              kyf = kxf + 1
              kzf = kyf + 1
!
              if (lklocal) then
                kx = 3*(atom2local(n2k) - 1) + 1
                ky = kx + 1
                kz = ky + 1
              endif
!
              lxf = 3*(n2l - 1) + 1
              lyf = lxf + 1
              lzf = lyf + 1
!
              if (lllocal) then
                lx = 3*(atom2local(n2l) - 1) + 1
                ly = lx + 1
                lz = ly + 1
              endif
!
              lsameijkl = .false.
              if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
!  Second derivatives : i - k
!
              if (lilocal.or.lklocal) then
                if (n1i.eq.n2k) then
                  oneik = 1.0_dp
                else
                  oneik = 0.0_dp
                endif
                cosik = xkv*xik + ykv*yik + zkv*zik
                sinik = sin(cosik)
                cosik = cos(cosik) - oneik
!
                if (lsameijkl) then
                  cosik = 0.5_dp*cosik
                  sinik = 0.5_dp*sinik
                endif
              endif
!
              if (lilocal) then
                derv2(kxf,ix) = derv2(kxf,ix) + xd1*xd2*d2(ind)*cosik
                derv2(kyf,ix) = derv2(kyf,ix) + xd1*yd2*d2(ind)*cosik
                derv2(kzf,ix) = derv2(kzf,ix) + xd1*zd2*d2(ind)*cosik
                derv2(kxf,iy) = derv2(kxf,iy) + yd1*xd2*d2(ind)*cosik
                derv2(kyf,iy) = derv2(kyf,iy) + yd1*yd2*d2(ind)*cosik
                derv2(kzf,iy) = derv2(kzf,iy) + yd1*zd2*d2(ind)*cosik
                derv2(kxf,iz) = derv2(kxf,iz) + zd1*xd2*d2(ind)*cosik
                derv2(kyf,iz) = derv2(kyf,iz) + zd1*yd2*d2(ind)*cosik
                derv2(kzf,iz) = derv2(kzf,iz) + zd1*zd2*d2(ind)*cosik
!
                dervi(kxf,ix) = dervi(kxf,ix) + xd1*xd2*d2(ind)*sinik
                dervi(kyf,ix) = dervi(kyf,ix) + xd1*yd2*d2(ind)*sinik
                dervi(kzf,ix) = dervi(kzf,ix) + xd1*zd2*d2(ind)*sinik
                dervi(kxf,iy) = dervi(kxf,iy) + yd1*xd2*d2(ind)*sinik
                dervi(kyf,iy) = dervi(kyf,iy) + yd1*yd2*d2(ind)*sinik
                dervi(kzf,iy) = dervi(kzf,iy) + yd1*zd2*d2(ind)*sinik
                dervi(kxf,iz) = dervi(kxf,iz) + zd1*xd2*d2(ind)*sinik
                dervi(kyf,iz) = dervi(kyf,iz) + zd1*yd2*d2(ind)*sinik
                dervi(kzf,iz) = dervi(kzf,iz) + zd1*zd2*d2(ind)*sinik
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosik
                  d2ks = d2(ind)*sinik
                  cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,kxf,ix) = derv2dk(1:3,kxf,ix) + xd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,ix) = derv2dk(1:3,kyf,ix) + xd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,ix) = derv2dk(1:3,kzf,ix) + xd1*zd2*cdk(1:3)
                  derv2dk(1:3,kxf,iy) = derv2dk(1:3,kxf,iy) + yd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,iy) = derv2dk(1:3,kyf,iy) + yd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,iy) = derv2dk(1:3,kzf,iy) + yd1*zd2*cdk(1:3)
                  derv2dk(1:3,kxf,iz) = derv2dk(1:3,kxf,iz) + zd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,iz) = derv2dk(1:3,kyf,iz) + zd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,iz) = derv2dk(1:3,kzf,iz) + zd1*zd2*cdk(1:3)
                endif
!
                if (n1.eq.n2) then
                  derv2(kxf,ix) = derv2(kxf,ix) + d1(nf2)*cosik
                  derv2(kyf,iy) = derv2(kyf,iy) + d1(nf2)*cosik
                  derv2(kzf,iz) = derv2(kzf,iz) + d1(nf2)*cosik
                  dervi(kxf,ix) = dervi(kxf,ix) + d1(nf2)*sinik
                  dervi(kyf,iy) = dervi(kyf,iy) + d1(nf2)*sinik
                  dervi(kzf,iz) = dervi(kzf,iz) + d1(nf2)*sinik
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(nf2)*cosik
                    d2ks = d1(nf2)*sinik
                    cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,kxf,ix) = derv2dk(1:3,kxf,ix) + cdk(1:3)
                    derv2dk(1:3,kyf,iy) = derv2dk(1:3,kyf,iy) + cdk(1:3)
                    derv2dk(1:3,kzf,iz) = derv2dk(1:3,kzf,iz) + cdk(1:3)
                  endif
                endif
              endif
!
              if (lklocal) then
                derv2(ixf,kx) = derv2(ixf,kx) + xd1*xd2*d2(ind)*cosik
                derv2(ixf,ky) = derv2(ixf,ky) + xd1*yd2*d2(ind)*cosik
                derv2(ixf,kz) = derv2(ixf,kz) + xd1*zd2*d2(ind)*cosik
                derv2(iyf,kx) = derv2(iyf,kx) + yd1*xd2*d2(ind)*cosik
                derv2(iyf,ky) = derv2(iyf,ky) + yd1*yd2*d2(ind)*cosik
                derv2(iyf,kz) = derv2(iyf,kz) + yd1*zd2*d2(ind)*cosik
                derv2(izf,kx) = derv2(izf,kx) + zd1*xd2*d2(ind)*cosik
                derv2(izf,ky) = derv2(izf,ky) + zd1*yd2*d2(ind)*cosik
                derv2(izf,kz) = derv2(izf,kz) + zd1*zd2*d2(ind)*cosik
!
                dervi(ixf,kx) = dervi(ixf,kx) - xd1*xd2*d2(ind)*sinik
                dervi(ixf,ky) = dervi(ixf,ky) - xd1*yd2*d2(ind)*sinik
                dervi(ixf,kz) = dervi(ixf,kz) - xd1*zd2*d2(ind)*sinik
                dervi(iyf,kx) = dervi(iyf,kx) - yd1*xd2*d2(ind)*sinik
                dervi(iyf,ky) = dervi(iyf,ky) - yd1*yd2*d2(ind)*sinik
                dervi(iyf,kz) = dervi(iyf,kz) - yd1*zd2*d2(ind)*sinik
                dervi(izf,kx) = dervi(izf,kx) - zd1*xd2*d2(ind)*sinik
                dervi(izf,ky) = dervi(izf,ky) - zd1*yd2*d2(ind)*sinik
                dervi(izf,kz) = dervi(izf,kz) - zd1*zd2*d2(ind)*sinik
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosik
                  d2ks = d2(ind)*sinik
                  cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,ixf,kx) = derv2dk(1:3,ixf,kx) + xd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,ixf,ky) = derv2dk(1:3,ixf,ky) + xd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,ixf,kz) = derv2dk(1:3,ixf,kz) + xd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,kx) = derv2dk(1:3,iyf,kx) + yd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,ky) = derv2dk(1:3,iyf,ky) + yd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,kz) = derv2dk(1:3,iyf,kz) + yd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,kx) = derv2dk(1:3,izf,kx) + zd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,ky) = derv2dk(1:3,izf,ky) + zd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,kz) = derv2dk(1:3,izf,kz) + zd1*zd2*conjg(cdk(1:3))
                endif
!
                if (n1.eq.n2) then
                  derv2(ixf,kx) = derv2(ixf,kx) + d1(nf2)*cosik
                  derv2(iyf,ky) = derv2(iyf,ky) + d1(nf2)*cosik
                  derv2(izf,kz) = derv2(izf,kz) + d1(nf2)*cosik
                  dervi(ixf,kx) = dervi(ixf,kx) - d1(nf2)*sinik
                  dervi(iyf,ky) = dervi(iyf,ky) - d1(nf2)*sinik
                  dervi(izf,kz) = dervi(izf,kz) - d1(nf2)*sinik
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(nf2)*cosik
                    d2ks = d1(nf2)*sinik
                    cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,ixf,kx) = derv2dk(1:3,ixf,kx) + conjg(cdk(1:3))
                    derv2dk(1:3,iyf,ky) = derv2dk(1:3,iyf,ky) + conjg(cdk(1:3))
                    derv2dk(1:3,izf,kz) = derv2dk(1:3,izf,kz) + conjg(cdk(1:3))
                  endif
                endif
              endif
!
!  Second derivatives : i - l
!
              if (lilocal.or.lllocal) then
                if (n1i.eq.n2l) then
                  oneil = 1.0_dp
                else
                  oneil = 0.0_dp
                endif
                cosil = xkv*xil + ykv*yil + zkv*zil
                sinil = sin(cosil)
                cosil = cos(cosil) - oneil
!
                if (lsameijkl) then
                  cosil = 0.5_dp*cosil
                  sinil = 0.5_dp*sinil
                endif
              endif
!
              if (lilocal) then
                derv2(lxf,ix) = derv2(lxf,ix) - xd1*xd2*d2(ind)*cosil
                derv2(lyf,ix) = derv2(lyf,ix) - xd1*yd2*d2(ind)*cosil
                derv2(lzf,ix) = derv2(lzf,ix) - xd1*zd2*d2(ind)*cosil
                derv2(lxf,iy) = derv2(lxf,iy) - yd1*xd2*d2(ind)*cosil
                derv2(lyf,iy) = derv2(lyf,iy) - yd1*yd2*d2(ind)*cosil
                derv2(lzf,iy) = derv2(lzf,iy) - yd1*zd2*d2(ind)*cosil
                derv2(lxf,iz) = derv2(lxf,iz) - zd1*xd2*d2(ind)*cosil
                derv2(lyf,iz) = derv2(lyf,iz) - zd1*yd2*d2(ind)*cosil
                derv2(lzf,iz) = derv2(lzf,iz) - zd1*zd2*d2(ind)*cosil
!
                dervi(lxf,ix) = dervi(lxf,ix) - xd1*xd2*d2(ind)*sinil
                dervi(lyf,ix) = dervi(lyf,ix) - xd1*yd2*d2(ind)*sinil
                dervi(lzf,ix) = dervi(lzf,ix) - xd1*zd2*d2(ind)*sinil
                dervi(lxf,iy) = dervi(lxf,iy) - yd1*xd2*d2(ind)*sinil
                dervi(lyf,iy) = dervi(lyf,iy) - yd1*yd2*d2(ind)*sinil
                dervi(lzf,iy) = dervi(lzf,iy) - yd1*zd2*d2(ind)*sinil
                dervi(lxf,iz) = dervi(lxf,iz) - zd1*xd2*d2(ind)*sinil
                dervi(lyf,iz) = dervi(lyf,iz) - zd1*yd2*d2(ind)*sinil
                dervi(lzf,iz) = dervi(lzf,iz) - zd1*zd2*d2(ind)*sinil
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosil
                  d2ks = d2(ind)*sinil
                  cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,lxf,ix) = derv2dk(1:3,lxf,ix) - xd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,ix) = derv2dk(1:3,lyf,ix) - xd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,ix) = derv2dk(1:3,lzf,ix) - xd1*zd2*cdk(1:3)
                  derv2dk(1:3,lxf,iy) = derv2dk(1:3,lxf,iy) - yd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,iy) = derv2dk(1:3,lyf,iy) - yd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,iy) = derv2dk(1:3,lzf,iy) - yd1*zd2*cdk(1:3)
                  derv2dk(1:3,lxf,iz) = derv2dk(1:3,lxf,iz) - zd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,iz) = derv2dk(1:3,lyf,iz) - zd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,iz) = derv2dk(1:3,lzf,iz) - zd1*zd2*cdk(1:3)
                endif
!
                if (n1.eq.n2) then
                  derv2(lxf,ix) = derv2(lxf,ix) - d1(nf2)*cosil
                  derv2(lyf,iy) = derv2(lyf,iy) - d1(nf2)*cosil
                  derv2(lzf,iz) = derv2(lzf,iz) - d1(nf2)*cosil
                  dervi(lxf,ix) = dervi(lxf,ix) - d1(nf2)*sinil
                  dervi(lyf,iy) = dervi(lyf,iy) - d1(nf2)*sinil
                  dervi(lzf,iz) = dervi(lzf,iz) - d1(nf2)*sinil
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(nf2)*cosil
                    d2ks = d1(nf2)*sinil
                    cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,lxf,ix) = derv2dk(1:3,lxf,ix) - cdk(1:3)
                    derv2dk(1:3,lyf,iy) = derv2dk(1:3,lyf,iy) - cdk(1:3)
                    derv2dk(1:3,lzf,iz) = derv2dk(1:3,lzf,iz) - cdk(1:3)
                  endif
                endif
              endif
!
              if (lllocal) then
                derv2(ixf,lx) = derv2(ixf,lx) - xd1*xd2*d2(ind)*cosil
                derv2(ixf,ly) = derv2(ixf,ly) - xd1*yd2*d2(ind)*cosil
                derv2(ixf,lz) = derv2(ixf,lz) - xd1*zd2*d2(ind)*cosil
                derv2(iyf,lx) = derv2(iyf,lx) - yd1*xd2*d2(ind)*cosil
                derv2(iyf,ly) = derv2(iyf,ly) - yd1*yd2*d2(ind)*cosil
                derv2(iyf,lz) = derv2(iyf,lz) - yd1*zd2*d2(ind)*cosil
                derv2(izf,lx) = derv2(izf,lx) - zd1*xd2*d2(ind)*cosil
                derv2(izf,ly) = derv2(izf,ly) - zd1*yd2*d2(ind)*cosil
                derv2(izf,lz) = derv2(izf,lz) - zd1*zd2*d2(ind)*cosil
!
                dervi(ixf,lx) = dervi(ixf,lx) + xd1*xd2*d2(ind)*sinil
                dervi(ixf,ly) = dervi(ixf,ly) + xd1*yd2*d2(ind)*sinil
                dervi(ixf,lz) = dervi(ixf,lz) + xd1*zd2*d2(ind)*sinil
                dervi(iyf,lx) = dervi(iyf,lx) + yd1*xd2*d2(ind)*sinil
                dervi(iyf,ly) = dervi(iyf,ly) + yd1*yd2*d2(ind)*sinil
                dervi(iyf,lz) = dervi(iyf,lz) + yd1*zd2*d2(ind)*sinil
                dervi(izf,lx) = dervi(izf,lx) + zd1*xd2*d2(ind)*sinil
                dervi(izf,ly) = dervi(izf,ly) + zd1*yd2*d2(ind)*sinil
                dervi(izf,lz) = dervi(izf,lz) + zd1*zd2*d2(ind)*sinil
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosil
                  d2ks = d2(ind)*sinil
                  cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,ixf,lx) = derv2dk(1:3,ixf,lx) - xd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,ixf,ly) = derv2dk(1:3,ixf,ly) - xd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,ixf,lz) = derv2dk(1:3,ixf,lz) - xd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,lx) = derv2dk(1:3,iyf,lx) - yd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,ly) = derv2dk(1:3,iyf,ly) - yd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,lz) = derv2dk(1:3,iyf,lz) - yd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,lx) = derv2dk(1:3,izf,lx) - zd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,ly) = derv2dk(1:3,izf,ly) - zd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,lz) = derv2dk(1:3,izf,lz) - zd1*zd2*conjg(cdk(1:3))
                endif
!
                if (n1.eq.n2) then
                  derv2(ixf,lx) = derv2(ixf,lx) - d1(nf2)*cosil
                  derv2(iyf,ly) = derv2(iyf,ly) - d1(nf2)*cosil
                  derv2(izf,lz) = derv2(izf,lz) - d1(nf2)*cosil
                  dervi(ixf,lx) = dervi(ixf,lx) + d1(nf2)*sinil
                  dervi(iyf,ly) = dervi(iyf,ly) + d1(nf2)*sinil
                  dervi(izf,lz) = dervi(izf,lz) + d1(nf2)*sinil
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(nf2)*cosil
                    d2ks = d1(nf2)*sinil
                    cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,ixf,lx) = derv2dk(1:3,ixf,lx) - conjg(cdk(1:3))
                    derv2dk(1:3,iyf,ly) = derv2dk(1:3,iyf,ly) - conjg(cdk(1:3))
                    derv2dk(1:3,izf,lz) = derv2dk(1:3,izf,lz) - conjg(cdk(1:3))
                  endif
                endif
              endif
!
!  Second derivatives : j - k
!
              if (ljlocal.or.lklocal) then
                if (n1j.eq.n2k) then
                  onejk = 1.0_dp
                else
                  onejk = 0.0_dp
                endif
                cosjk = xkv*xjk + ykv*yjk + zkv*zjk
                sinjk = sin(cosjk)
                cosjk = cos(cosjk) - onejk
!
                if (lsameijkl) then
                  cosjk = 0.5_dp*cosjk
                  sinjk = 0.5_dp*sinjk
                endif
              endif
!
              if (ljlocal) then
                derv2(kxf,jx) = derv2(kxf,jx) - xd1*xd2*d2(ind)*cosjk
                derv2(kyf,jx) = derv2(kyf,jx) - xd1*yd2*d2(ind)*cosjk
                derv2(kzf,jx) = derv2(kzf,jx) - xd1*zd2*d2(ind)*cosjk
                derv2(kxf,jy) = derv2(kxf,jy) - yd1*xd2*d2(ind)*cosjk
                derv2(kyf,jy) = derv2(kyf,jy) - yd1*yd2*d2(ind)*cosjk
                derv2(kzf,jy) = derv2(kzf,jy) - yd1*zd2*d2(ind)*cosjk
                derv2(kxf,jz) = derv2(kxf,jz) - zd1*xd2*d2(ind)*cosjk
                derv2(kyf,jz) = derv2(kyf,jz) - zd1*yd2*d2(ind)*cosjk
                derv2(kzf,jz) = derv2(kzf,jz) - zd1*zd2*d2(ind)*cosjk
!
                dervi(kxf,jx) = dervi(kxf,jx) - xd1*xd2*d2(ind)*sinjk
                dervi(kyf,jx) = dervi(kyf,jx) - xd1*yd2*d2(ind)*sinjk
                dervi(kzf,jx) = dervi(kzf,jx) - xd1*zd2*d2(ind)*sinjk
                dervi(kxf,jy) = dervi(kxf,jy) - yd1*xd2*d2(ind)*sinjk
                dervi(kyf,jy) = dervi(kyf,jy) - yd1*yd2*d2(ind)*sinjk
                dervi(kzf,jy) = dervi(kzf,jy) - yd1*zd2*d2(ind)*sinjk
                dervi(kxf,jz) = dervi(kxf,jz) - zd1*xd2*d2(ind)*sinjk
                dervi(kyf,jz) = dervi(kyf,jz) - zd1*yd2*d2(ind)*sinjk
                dervi(kzf,jz) = dervi(kzf,jz) - zd1*zd2*d2(ind)*sinjk
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosjk
                  d2ks = d2(ind)*sinjk
                  cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,kxf,jx) = derv2dk(1:3,kxf,jx) - xd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,jx) = derv2dk(1:3,kyf,jx) - xd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,jx) = derv2dk(1:3,kzf,jx) - xd1*zd2*cdk(1:3)
                  derv2dk(1:3,kxf,jy) = derv2dk(1:3,kxf,jy) - yd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,jy) = derv2dk(1:3,kyf,jy) - yd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,jy) = derv2dk(1:3,kzf,jy) - yd1*zd2*cdk(1:3)
                  derv2dk(1:3,kxf,jz) = derv2dk(1:3,kxf,jz) - zd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,jz) = derv2dk(1:3,kyf,jz) - zd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,jz) = derv2dk(1:3,kzf,jz) - zd1*zd2*cdk(1:3)
                endif
!
                if (n1.eq.n2) then
                  derv2(kxf,jx) = derv2(kxf,jx) - d1(nf2)*cosjk
                  derv2(kyf,jy) = derv2(kyf,jy) - d1(nf2)*cosjk
                  derv2(kzf,jz) = derv2(kzf,jz) - d1(nf2)*cosjk
                  dervi(kxf,jx) = dervi(kxf,jx) - d1(nf2)*sinjk
                  dervi(kyf,jy) = dervi(kyf,jy) - d1(nf2)*sinjk
                  dervi(kzf,jz) = dervi(kzf,jz) - d1(nf2)*sinjk
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(nf2)*cosjk
                    d2ks = d1(nf2)*sinjk
                    cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,kxf,jx) = derv2dk(1:3,kxf,jx) - cdk(1:3)
                    derv2dk(1:3,kyf,jy) = derv2dk(1:3,kyf,jy) - cdk(1:3)
                    derv2dk(1:3,kzf,jz) = derv2dk(1:3,kzf,jz) - cdk(1:3)
                  endif
                endif
              endif
!
              if (lklocal) then
                derv2(jxf,kx) = derv2(jxf,kx) - xd1*xd2*d2(ind)*cosjk
                derv2(jxf,ky) = derv2(jxf,ky) - xd1*yd2*d2(ind)*cosjk
                derv2(jxf,kz) = derv2(jxf,kz) - xd1*zd2*d2(ind)*cosjk
                derv2(jyf,kx) = derv2(jyf,kx) - yd1*xd2*d2(ind)*cosjk
                derv2(jyf,ky) = derv2(jyf,ky) - yd1*yd2*d2(ind)*cosjk
                derv2(jyf,kz) = derv2(jyf,kz) - yd1*zd2*d2(ind)*cosjk
                derv2(jzf,kx) = derv2(jzf,kx) - zd1*xd2*d2(ind)*cosjk
                derv2(jzf,ky) = derv2(jzf,ky) - zd1*yd2*d2(ind)*cosjk
                derv2(jzf,kz) = derv2(jzf,kz) - zd1*zd2*d2(ind)*cosjk
!
                dervi(jxf,kx) = dervi(jxf,kx) + xd1*xd2*d2(ind)*sinjk
                dervi(jxf,ky) = dervi(jxf,ky) + xd1*yd2*d2(ind)*sinjk
                dervi(jxf,kz) = dervi(jxf,kz) + xd1*zd2*d2(ind)*sinjk
                dervi(jyf,kx) = dervi(jyf,kx) + yd1*xd2*d2(ind)*sinjk
                dervi(jyf,ky) = dervi(jyf,ky) + yd1*yd2*d2(ind)*sinjk
                dervi(jyf,kz) = dervi(jyf,kz) + yd1*zd2*d2(ind)*sinjk
                dervi(jzf,kx) = dervi(jzf,kx) + zd1*xd2*d2(ind)*sinjk
                dervi(jzf,ky) = dervi(jzf,ky) + zd1*yd2*d2(ind)*sinjk
                dervi(jzf,kz) = dervi(jzf,kz) + zd1*zd2*d2(ind)*sinjk
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosjk
                  d2ks = d2(ind)*sinjk
                  cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,jxf,kx) = derv2dk(1:3,jxf,kx) - xd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jxf,ky) = derv2dk(1:3,jxf,ky) - xd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jxf,kz) = derv2dk(1:3,jxf,kz) - xd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,kx) = derv2dk(1:3,jyf,kx) - yd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,ky) = derv2dk(1:3,jyf,ky) - yd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,kz) = derv2dk(1:3,jyf,kz) - yd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,kx) = derv2dk(1:3,jzf,kx) - zd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,ky) = derv2dk(1:3,jzf,ky) - zd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,kz) = derv2dk(1:3,jzf,kz) - zd1*zd2*conjg(cdk(1:3))
                endif
!
                if (n1.eq.n2) then
                  derv2(jxf,kx) = derv2(jxf,kx) - d1(nf2)*cosjk
                  derv2(jyf,ky) = derv2(jyf,ky) - d1(nf2)*cosjk
                  derv2(jzf,kz) = derv2(jzf,kz) - d1(nf2)*cosjk
                  dervi(jxf,kx) = dervi(jxf,kx) + d1(nf2)*sinjk
                  dervi(jyf,ky) = dervi(jyf,ky) + d1(nf2)*sinjk
                  dervi(jzf,kz) = dervi(jzf,kz) + d1(nf2)*sinjk
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(nf2)*cosjk
                    d2ks = d1(nf2)*sinjk
                    cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,jxf,kx) = derv2dk(1:3,jxf,kx) - conjg(cdk(1:3))
                    derv2dk(1:3,jyf,ky) = derv2dk(1:3,jyf,ky) - conjg(cdk(1:3))
                    derv2dk(1:3,jzf,kz) = derv2dk(1:3,jzf,kz) - conjg(cdk(1:3))
                  endif
                endif
              endif
!
!  Second derivatives : j - l
!
              if (ljlocal.or.lllocal) then
                if (n1j.eq.n2l) then
                  onejl = 1.0_dp
                else
                  onejl = 0.0_dp
                endif
                cosjl = xkv*xjl + ykv*yjl + zkv*zjl
                sinjl = sin(cosjl)
                cosjl = cos(cosjl) - onejl
!
                if (lsameijkl) then
                  cosjl = 0.5_dp*cosjl
                  sinjl = 0.5_dp*sinjl
                endif
              endif
!
              if (ljlocal) then
                derv2(lxf,jx) = derv2(lxf,jx) + xd1*xd2*d2(ind)*cosjl
                derv2(lyf,jx) = derv2(lyf,jx) + xd1*yd2*d2(ind)*cosjl
                derv2(lzf,jx) = derv2(lzf,jx) + xd1*zd2*d2(ind)*cosjl
                derv2(lxf,jy) = derv2(lxf,jy) + yd1*xd2*d2(ind)*cosjl
                derv2(lyf,jy) = derv2(lyf,jy) + yd1*yd2*d2(ind)*cosjl
                derv2(lzf,jy) = derv2(lzf,jy) + yd1*zd2*d2(ind)*cosjl
                derv2(lxf,jz) = derv2(lxf,jz) + zd1*xd2*d2(ind)*cosjl
                derv2(lyf,jz) = derv2(lyf,jz) + zd1*yd2*d2(ind)*cosjl
                derv2(lzf,jz) = derv2(lzf,jz) + zd1*zd2*d2(ind)*cosjl
!
                dervi(lxf,jx) = dervi(lxf,jx) + xd1*xd2*d2(ind)*sinjl
                dervi(lyf,jx) = dervi(lyf,jx) + xd1*yd2*d2(ind)*sinjl
                dervi(lzf,jx) = dervi(lzf,jx) + xd1*zd2*d2(ind)*sinjl
                dervi(lxf,jy) = dervi(lxf,jy) + yd1*xd2*d2(ind)*sinjl
                dervi(lyf,jy) = dervi(lyf,jy) + yd1*yd2*d2(ind)*sinjl
                dervi(lzf,jy) = dervi(lzf,jy) + yd1*zd2*d2(ind)*sinjl
                dervi(lxf,jz) = dervi(lxf,jz) + zd1*xd2*d2(ind)*sinjl
                dervi(lyf,jz) = dervi(lyf,jz) + zd1*yd2*d2(ind)*sinjl
                dervi(lzf,jz) = dervi(lzf,jz) + zd1*zd2*d2(ind)*sinjl
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosjl
                  d2ks = d2(ind)*sinjl
                  cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,lxf,jx) = derv2dk(1:3,lxf,jx) + xd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,jx) = derv2dk(1:3,lyf,jx) + xd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,jx) = derv2dk(1:3,lzf,jx) + xd1*zd2*cdk(1:3)
                  derv2dk(1:3,lxf,jy) = derv2dk(1:3,lxf,jy) + yd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,jy) = derv2dk(1:3,lyf,jy) + yd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,jy) = derv2dk(1:3,lzf,jy) + yd1*zd2*cdk(1:3)
                  derv2dk(1:3,lxf,jz) = derv2dk(1:3,lxf,jz) + zd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,jz) = derv2dk(1:3,lyf,jz) + zd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,jz) = derv2dk(1:3,lzf,jz) + zd1*zd2*cdk(1:3)
                endif
!
                if (n1.eq.n2) then
                  derv2(lxf,jx) = derv2(lxf,jx) + d1(nf2)*cosjl
                  derv2(lyf,jy) = derv2(lyf,jy) + d1(nf2)*cosjl
                  derv2(lzf,jz) = derv2(lzf,jz) + d1(nf2)*cosjl
                  dervi(lxf,jx) = dervi(lxf,jx) + d1(nf2)*sinjl
                  dervi(lyf,jy) = dervi(lyf,jy) + d1(nf2)*sinjl
                  dervi(lzf,jz) = dervi(lzf,jz) + d1(nf2)*sinjl
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(nf2)*cosjl
                    d2ks = d1(nf2)*sinjl
                    cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,lxf,jx) = derv2dk(1:3,lxf,jx) + cdk(1:3)
                    derv2dk(1:3,lyf,jy) = derv2dk(1:3,lyf,jy) + cdk(1:3)
                    derv2dk(1:3,lzf,jz) = derv2dk(1:3,lzf,jz) + cdk(1:3)
                  endif
                endif
              endif
!
              if (lllocal) then
                derv2(jxf,lx) = derv2(jxf,lx) + xd1*xd2*d2(ind)*cosjl
                derv2(jxf,ly) = derv2(jxf,ly) + xd1*yd2*d2(ind)*cosjl
                derv2(jxf,lz) = derv2(jxf,lz) + xd1*zd2*d2(ind)*cosjl
                derv2(jyf,lx) = derv2(jyf,lx) + yd1*xd2*d2(ind)*cosjl
                derv2(jyf,ly) = derv2(jyf,ly) + yd1*yd2*d2(ind)*cosjl
                derv2(jyf,lz) = derv2(jyf,lz) + yd1*zd2*d2(ind)*cosjl
                derv2(jzf,lx) = derv2(jzf,lx) + zd1*xd2*d2(ind)*cosjl
                derv2(jzf,ly) = derv2(jzf,ly) + zd1*yd2*d2(ind)*cosjl
                derv2(jzf,lz) = derv2(jzf,lz) + zd1*zd2*d2(ind)*cosjl
!
                dervi(jxf,lx) = dervi(jxf,lx) - xd1*xd2*d2(ind)*sinjl
                dervi(jxf,ly) = dervi(jxf,ly) - xd1*yd2*d2(ind)*sinjl
                dervi(jxf,lz) = dervi(jxf,lz) - xd1*zd2*d2(ind)*sinjl
                dervi(jyf,lx) = dervi(jyf,lx) - yd1*xd2*d2(ind)*sinjl
                dervi(jyf,ly) = dervi(jyf,ly) - yd1*yd2*d2(ind)*sinjl
                dervi(jyf,lz) = dervi(jyf,lz) - yd1*zd2*d2(ind)*sinjl
                dervi(jzf,lx) = dervi(jzf,lx) - zd1*xd2*d2(ind)*sinjl
                dervi(jzf,ly) = dervi(jzf,ly) - zd1*yd2*d2(ind)*sinjl
                dervi(jzf,lz) = dervi(jzf,lz) - zd1*zd2*d2(ind)*sinjl
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosjl
                  d2ks = d2(ind)*sinjl
                  cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,jxf,lx) = derv2dk(1:3,jxf,lx) + xd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jxf,ly) = derv2dk(1:3,jxf,ly) + xd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jxf,lz) = derv2dk(1:3,jxf,lz) + xd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,lx) = derv2dk(1:3,jyf,lx) + yd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,ly) = derv2dk(1:3,jyf,ly) + yd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,lz) = derv2dk(1:3,jyf,lz) + yd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,lx) = derv2dk(1:3,jzf,lx) + zd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,ly) = derv2dk(1:3,jzf,ly) + zd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,lz) = derv2dk(1:3,jzf,lz) + zd1*zd2*conjg(cdk(1:3))
                endif
!
                if (n1.eq.n2) then
                  derv2(jxf,lx) = derv2(jxf,lx) + d1(nf2)*cosjl
                  derv2(jyf,ly) = derv2(jyf,ly) + d1(nf2)*cosjl
                  derv2(jzf,lz) = derv2(jzf,lz) + d1(nf2)*cosjl
                  dervi(jxf,lx) = dervi(jxf,lx) - d1(nf2)*sinjl
                  dervi(jyf,ly) = dervi(jyf,ly) - d1(nf2)*sinjl
                  dervi(jzf,lz) = dervi(jzf,lz) - d1(nf2)*sinjl
!
                  if (lgroupvelocity) then
!
!  Group velocities
!           
                    d2k  = d1(nf2)*cosjl
                    d2ks = d1(nf2)*sinjl
                    cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!           
                    derv2dk(1:3,jxf,lx) = derv2dk(1:3,jxf,lx) + conjg(cdk(1:3))
                    derv2dk(1:3,jyf,ly) = derv2dk(1:3,jyf,ly) + conjg(cdk(1:3))
                    derv2dk(1:3,jzf,lz) = derv2dk(1:3,jzf,lz) + conjg(cdk(1:3))
                  endif
                endif
              endif
            endif
          endif
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addpdmsp')
#endif
!
  return
  end
!
  subroutine d2adddmsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1,d2, &
                       nonzero,nonzeroptr,ltorderv,lextended)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array.
!
!  Distributed memory parallel version.
!  Sparse storage of derivatives version.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic distances
!  d2              = array of first derivatives w.r.t. interatomic distances
!  nonzero         = number of nonzero elements in d1
!  nonzeroptr(*)   = pointer to nonzero elements in d1
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!
!   8/23 Created from d2adddm and d2addsp
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
  use brennerdata
  use current,        only : ndim, nstrains
  use current,        only : nstrains, nstrptr
  use derivatives
  use parallel
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nonzero
  integer(i4), intent(in)             :: nonzeroptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ind0
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: is1,is2
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: ixf,iyf,izf
  integer(i4)                         :: jxf,jyf,jzf
  integer(i4)                         :: kxf,kyf,kzf
  integer(i4)                         :: lxf,lyf,lzf
  integer(i4)                         :: ix1,iy1,iz1
  integer(i4)                         :: ix2,iy2,iz2
  integer(i4)                         :: jx1,jy1,jz1
  integer(i4)                         :: jx2,jy2,jz2
  integer(i4)                         :: kx1,ky1,kz1
  integer(i4)                         :: kx2,ky2,kz2
  integer(i4)                         :: lx1,ly1,lz1
  integer(i4)                         :: lx2,ly2,lz2
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nf1
  integer(i4)                         :: nf2
  integer(i4)                         :: nji
  integer(i4)                         :: nkj
  integer(i4)                         :: nlk
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: ns1
  integer(i4)                         :: ns2
  integer(i4)                         :: ntmp
  real(dp)                            :: one1
  real(dp)                            :: one2
  real(dp)                            :: dr2ds1(6)
  real(dp)                            :: dr2ds2(6)
  real(dp)                            :: d2r2ds21(6,6)
  real(dp)                            :: d2r2ds22(6,6)
  real(dp)                            :: d2r2dsdx1(6,3)
  real(dp)                            :: d2r2dsdx2(6,3)
  real(dp)                            :: d2r2dx21(3,3)
  real(dp)                            :: d2r2dx22(3,3)
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: ldo12
  logical                             :: ldo21
  logical                             :: lsameijkl
  logical                             :: lilocal
  logical                             :: ljlocal
  logical                             :: lklocal
  logical                             :: lllocal
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lopi
  logical                             :: lopj
  logical                             :: lopk
  logical                             :: lopl
#ifdef TRACE
  call trace_in('d2adddm')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndisttor = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
    if (lextended) then
      ndist = ndisttor + nneigh(nri)*mneigh*mneigh
    else
      ndist = ndisttor
    endif
  else
    ndist = nneigh2 
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = nneigh2 + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
!  Loop over pairs of interatomic distances
!
  do nf1 = 1,nonzero
    n1 = nonzeroptr(nf1)
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Check that neighbour numbers are valid
!
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    elseif (n1.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
!  Check that atoms are valid
!
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
      else
        ljok = .false.
      endif
!
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!         
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    else
!
!  Extended derivatives
!
      ntmp = n1 - ndisttor
      nji = (ntmp - 1)/(mneigh*mneigh) + 1
      ntmp = ntmp - (nji - 1)*mneigh*mneigh
      nkj = (ntmp - 1)/mneigh + 1
      ntmp = ntmp - (nkj - 1)*mneigh
      nlk = ntmp
!
      nrj = neighno(nji,nri)
      if (nrj.eq.0) cycle
      if (nkj.gt.nneigh(nrj)) cycle
      nrk = neighno(nkj,nrj)
      if (nrk.eq.0) cycle
      if (nlk.gt.nneigh(nrk)) cycle
      n1i = nrk
      n1j = neighno(nlk,nrk)
      if (n1j.eq.0) cycle
!
!  Check that atoms are valid
!
      liok = (nkj.le.nneigh(nrj))
      ljok = (nlk.le.nneigh(nrk))
!
      if (liok.and.ljok) then
        x2i = xneigh(nji,nri) + xneigh(nkj,nrj)
        y2i = yneigh(nji,nri) + yneigh(nkj,nrj)
        z2i = zneigh(nji,nri) + zneigh(nkj,nrj)
!
        xd1 = xneigh(nlk,nrk)
        yd1 = yneigh(nlk,nrk)
        zd1 = zneigh(nlk,nrk)
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
!
!  Set up products for strain
!
      if (lstr) then
        call real1strterm(ndim,xd1,yd1,zd1,0.0_dp,0.0_dp,0.0_dp,dr2ds1,d2r2dx21,d2r2dsdx1,d2r2ds21,.true.)
      endif
!
!  Set second derivative location pointers
!
      lopi = (nfreeatom(n1i).gt.0)
      lopj = (nfreeatom(n1j).gt.0)
!
      lilocal = (atom2local(n1i).gt.0)
      if (lilocal) then
        ix = 3*(atom2local(nfreeatom(n1i)) - 1) + 1
        iy = ix + 1
        iz = iy + 1
      endif
      ixf = 3*(nfreeatom(n1i) - 1) + 1
      iyf = ixf + 1
      izf = iyf + 1
!
      ljlocal = (atom2local(n1j).gt.0)
      if (ljlocal) then
        jx = 3*(atom2local(nfreeatom(n1j)) - 1) + 1
        jy = jx + 1
        jz = jy + 1
      endif
      jxf = 3*(nfreeatom(n1j) - 1) + 1
      jyf = jxf + 1
      jzf = jyf + 1
    endif
!
!  If neither i nor j are being optimised and strain derivatives are not needed then atoms are not valid
!
    if (.not.lopi.and..not.lopj.and..not.lstr) then
      liok = .false.
      ljok = .false.
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
      ind0 = nf1*(nf1-1)/2
      do nf2 = 1,nf1
        n2 = nonzeroptr(nf2)
        ind = ind0 + nf2
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(nf2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!  
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          elseif (n2.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
            else
              llok = .false.
            endif
!
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!               
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          else
!
!  Extended derivatives
!
            ntmp = n2 - ndisttor
            nji = (ntmp - 1)/(mneigh*mneigh) + 1
            ntmp = ntmp - (nji - 1)*mneigh*mneigh
            nkj = (ntmp - 1)/mneigh + 1
            ntmp = ntmp - (nkj - 1)*mneigh
            nlk = ntmp
!
            nrj = neighno(nji,nri)
            if (nrj.eq.0) cycle
            if (nkj.gt.nneigh(nrj)) cycle
            nrk = neighno(nkj,nrj)
            if (nrk.eq.0) cycle
            if (nlk.gt.nneigh(nrk)) cycle
            n2k = nrk
            n2l = neighno(nlk,nrk)
            if (n2l.eq.0) cycle
!
!  Check that atoms are valid
!
            lkok = (nkj.le.nneigh(nrj))
            llok = (nlk.le.nneigh(nrk))
!
            if (liok.and.ljok) then
              xik = xneigh(nji,nri) + xneigh(nkj,nrj) - x2i
              yik = yneigh(nji,nri) + yneigh(nkj,nrj) - y2i
              zik = zneigh(nji,nri) + zneigh(nkj,nrj) - z2i
!
              xd2 = xneigh(nlk,nrk)
              yd2 = yneigh(nlk,nrk)
              zd2 = zneigh(nlk,nrk)
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
!  Set up products for strain
!
            if (lstr) then
              call real1strterm(ndim,xd2,yd2,zd2,0.0_dp,0.0_dp,0.0_dp,dr2ds2,d2r2dx22,d2r2dsdx2,d2r2ds22,.true.)
            endif
!
            lopk = (nfreeatom(n2k).gt.0)
            lopl = (nfreeatom(n2l).gt.0)
!
            lklocal = (atom2local(n2k).gt.0)
            if (lklocal) then
              kx = 3*(atom2local(nfreeatom(n2k)) - 1) + 1
              ky = kx + 1
              kz = ky + 1
            endif
            kxf = 3*(nfreeatom(n2k) - 1) + 1
            kyf = kxf + 1
            kzf = kyf + 1
!
            lllocal = (atom2local(n2l).gt.0)
            if (lllocal) then
              lx = 3*(atom2local(nfreeatom(n2l)) - 1) + 1
              ly = lx + 1
              lz = ly + 1
            endif
            lxf = 3*(nfreeatom(n2l) - 1) + 1
            lyf = lxf + 1
            lzf = lyf + 1
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if ((lopi.or.lopk).and.n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              if (lopi.and.lopk) then
                ldo21 = lilocal
                ldo12 = lklocal
                if (ldo21) then
                  ix1 = ix
                  iy1 = iy
                  iz1 = iz
                  kx2 = kxf
                  ky2 = kyf
                  kz2 = kzf
                endif
                if (ldo12) then
                  ix2 = ixf
                  iy2 = iyf
                  iz2 = izf
                  kx1 = kx
                  ky1 = ky
                  kz1 = kz
                endif
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ldo21 = lilocal
                ldo12 = lilocal
                if (ldo21) then
                  ix1 = ix
                  iy1 = iy
                  iz1 = iz
                  kx2 = ixf
                  ky2 = iyf
                  kz2 = izf
                endif
                if (ldo12) then
                  ix2 = ixf
                  iy2 = iyf
                  iz2 = izf
                  kx1 = ix
                  ky1 = iy
                  kz1 = iz
                endif
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                ldo21 = lklocal
                ldo12 = lklocal
                if (ldo21) then
                  ix1 = kx
                  iy1 = ky
                  iz1 = kz
                  kx2 = kxf
                  ky2 = kyf
                  kz2 = kzf
                endif
                if (ldo12) then
                  ix2 = kxf
                  iy2 = kyf
                  iz2 = kzf
                  kx1 = kx
                  ky1 = ky
                  kz1 = kz
                endif
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              if (ldo21) then
                derv2(kx2,ix1) = derv2(kx2,ix1) + xd1*xd2*d2(ind)*one1
                derv2(ky2,ix1) = derv2(ky2,ix1) + xd1*yd2*d2(ind)*one1
                derv2(kz2,ix1) = derv2(kz2,ix1) + xd1*zd2*d2(ind)*one1
                derv2(kx2,iy1) = derv2(kx2,iy1) + yd1*xd2*d2(ind)*one1
                derv2(ky2,iy1) = derv2(ky2,iy1) + yd1*yd2*d2(ind)*one1
                derv2(kz2,iy1) = derv2(kz2,iy1) + yd1*zd2*d2(ind)*one1
                derv2(kx2,iz1) = derv2(kx2,iz1) + zd1*xd2*d2(ind)*one1
                derv2(ky2,iz1) = derv2(ky2,iz1) + zd1*yd2*d2(ind)*one1
                derv2(kz2,iz1) = derv2(kz2,iz1) + zd1*zd2*d2(ind)*one1
                if (n1.eq.n2) then
                  derv2(kx2,ix1) = derv2(kx2,ix1) + d1(nf2)*one1
                  derv2(ky2,iy1) = derv2(ky2,iy1) + d1(nf2)*one1
                  derv2(kz2,iz1) = derv2(kz2,iz1) + d1(nf2)*one1
                endif
              endif
              if (.not.lsameijkl.and.ldo12) then
                derv2(ix2,kx1) = derv2(ix2,kx1) + xd1*xd2*d2(ind)*one2
                derv2(ix2,ky1) = derv2(ix2,ky1) + xd1*yd2*d2(ind)*one2
                derv2(ix2,kz1) = derv2(ix2,kz1) + xd1*zd2*d2(ind)*one2
                derv2(iy2,kx1) = derv2(iy2,kx1) + yd1*xd2*d2(ind)*one2
                derv2(iy2,ky1) = derv2(iy2,ky1) + yd1*yd2*d2(ind)*one2
                derv2(iy2,kz1) = derv2(iy2,kz1) + yd1*zd2*d2(ind)*one2
                derv2(iz2,kx1) = derv2(iz2,kx1) + zd1*xd2*d2(ind)*one2
                derv2(iz2,ky1) = derv2(iz2,ky1) + zd1*yd2*d2(ind)*one2
                derv2(iz2,kz1) = derv2(iz2,kz1) + zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(ix2,kx1) = derv2(ix2,kx1) + d1(nf2)*one2
                  derv2(iy2,ky1) = derv2(iy2,ky1) + d1(nf2)*one2
                  derv2(iz2,kz1) = derv2(iz2,kz1) + d1(nf2)*one2
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if ((lopi.or.lopl).and.n1i.ne.n2l) then
              if (lopi.and.lopl) then
                ldo21 = lilocal
                ldo12 = lllocal
                if (ldo21) then
                  ix1 = ix
                  iy1 = iy
                  iz1 = iz
                  lx2 = lxf
                  ly2 = lyf
                  lz2 = lzf
                endif
                if (ldo12) then
                  ix2 = ixf
                  iy2 = iyf
                  iz2 = izf
                  lx1 = lx
                  ly1 = ly
                  lz1 = lz
                endif
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ldo21 = lilocal
                ldo12 = lilocal
                if (ldo21) then
                  ix1 = ix
                  iy1 = iy
                  iz1 = iz
                  lx2 = ixf
                  ly2 = iyf
                  lz2 = izf
                endif
                if (ldo12) then
                  ix2 = ixf
                  iy2 = iyf
                  iz2 = izf
                  lx1 = ix
                  ly1 = iy
                  lz1 = iz
                endif
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                ldo21 = lllocal
                ldo12 = lllocal
                if (ldo21) then
                  ix1 = lx
                  iy1 = ly
                  iz1 = lz
                  lx2 = lxf
                  ly2 = lyf
                  lz2 = lzf
                endif
                if (ldo12) then
                  ix2 = lxf
                  iy2 = lyf
                  iz2 = lzf
                  lx1 = lx
                  ly1 = ly
                  lz1 = lz
                endif
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              if (ldo21) then
                derv2(lx2,ix1) = derv2(lx2,ix1) - xd1*xd2*d2(ind)*one1
                derv2(ly2,ix1) = derv2(ly2,ix1) - xd1*yd2*d2(ind)*one1
                derv2(lz2,ix1) = derv2(lz2,ix1) - xd1*zd2*d2(ind)*one1
                derv2(lx2,iy1) = derv2(lx2,iy1) - yd1*xd2*d2(ind)*one1
                derv2(ly2,iy1) = derv2(ly2,iy1) - yd1*yd2*d2(ind)*one1
                derv2(lz2,iy1) = derv2(lz2,iy1) - yd1*zd2*d2(ind)*one1
                derv2(lx2,iz1) = derv2(lx2,iz1) - zd1*xd2*d2(ind)*one1
                derv2(ly2,iz1) = derv2(ly2,iz1) - zd1*yd2*d2(ind)*one1
                derv2(lz2,iz1) = derv2(lz2,iz1) - zd1*zd2*d2(ind)*one1
                if (n1.eq.n2) then
                  derv2(lx2,ix1) = derv2(lx2,ix1) - d1(nf2)*one1
                  derv2(ly2,iy1) = derv2(ly2,iy1) - d1(nf2)*one1
                  derv2(lz2,iz1) = derv2(lz2,iz1) - d1(nf2)*one1
                endif
              endif
              if (.not.lsameijkl.and.ldo12) then
                derv2(ix2,lx1) = derv2(ix2,lx1) - xd1*xd2*d2(ind)*one2
                derv2(ix2,ly1) = derv2(ix2,ly1) - xd1*yd2*d2(ind)*one2
                derv2(ix2,lz1) = derv2(ix2,lz1) - xd1*zd2*d2(ind)*one2
                derv2(iy2,lx1) = derv2(iy2,lx1) - yd1*xd2*d2(ind)*one2
                derv2(iy2,ly1) = derv2(iy2,ly1) - yd1*yd2*d2(ind)*one2
                derv2(iy2,lz1) = derv2(iy2,lz1) - yd1*zd2*d2(ind)*one2
                derv2(iz2,lx1) = derv2(iz2,lx1) - zd1*xd2*d2(ind)*one2
                derv2(iz2,ly1) = derv2(iz2,ly1) - zd1*yd2*d2(ind)*one2
                derv2(iz2,lz1) = derv2(iz2,lz1) - zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(ix2,lx1) = derv2(ix2,lx1) - d1(nf2)*one2
                  derv2(iy2,ly1) = derv2(iy2,ly1) - d1(nf2)*one2
                  derv2(iz2,lz1) = derv2(iz2,lz1) - d1(nf2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if ((lopj.or.lopk).and.n1j.ne.n2k) then
              if (lopj.and.lopk) then
                ldo21 = ljlocal
                ldo12 = lklocal
                if (ldo21) then
                  jx1 = jx
                  jy1 = jy
                  jz1 = jz
                  kx2 = kxf
                  ky2 = kyf
                  kz2 = kzf
                endif
                if (ldo12) then
                  jx2 = jxf
                  jy2 = jyf
                  jz2 = jzf
                  kx1 = kx
                  ky1 = ky
                  kz1 = kz
                endif
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                ldo21 = ljlocal
                ldo12 = ljlocal
                if (ldo21) then
                  jx1 = jx
                  jy1 = jy
                  jz1 = jz
                  kx2 = jxf
                  ky2 = jyf
                  kz2 = jzf
                endif
                if (ldo12) then
                  jx2 = jxf
                  jy2 = jyf
                  jz2 = jzf
                  kx1 = jx
                  ky1 = jy
                  kz1 = jz
                endif
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                ldo21 = lklocal
                ldo12 = lklocal
                if (ldo21) then
                  jx1 = kx
                  jy1 = ky
                  jz1 = kz
                  kx2 = kxf
                  ky2 = kyf
                  kz2 = kzf
                endif
                if (ldo12) then
                  jx2 = kxf
                  iy2 = kyf
                  iz2 = kzf
                  kx1 = kx
                  ky1 = ky
                  kz1 = kz
                endif
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              if (ldo21) then
                derv2(kx2,jx1) = derv2(kx2,jx1) - xd1*xd2*d2(ind)*one1
                derv2(ky2,jx1) = derv2(ky2,jx1) - xd1*yd2*d2(ind)*one1
                derv2(kz2,jx1) = derv2(kz2,jx1) - xd1*zd2*d2(ind)*one1
                derv2(kx2,jy1) = derv2(kx2,jy1) - yd1*xd2*d2(ind)*one1
                derv2(ky2,jy1) = derv2(ky2,jy1) - yd1*yd2*d2(ind)*one1
                derv2(kz2,jy1) = derv2(kz2,jy1) - yd1*zd2*d2(ind)*one1
                derv2(kx2,jz1) = derv2(kx2,jz1) - zd1*xd2*d2(ind)*one1
                derv2(ky2,jz1) = derv2(ky2,jz1) - zd1*yd2*d2(ind)*one1
                derv2(kz2,jz1) = derv2(kz2,jz1) - zd1*zd2*d2(ind)*one1
                if (n1.eq.n2) then
                  derv2(kx2,jx1) = derv2(kx2,jx1) - d1(nf2)*one1
                  derv2(ky2,jy1) = derv2(ky2,jy1) - d1(nf2)*one1
                  derv2(kz2,jz1) = derv2(kz2,jz1) - d1(nf2)*one1
                endif
              endif
              if (.not.lsameijkl.and.ldo12) then
                derv2(jx2,kx1) = derv2(jx2,kx1) - xd1*xd2*d2(ind)*one2
                derv2(jx2,ky1) = derv2(jx2,ky1) - xd1*yd2*d2(ind)*one2
                derv2(jx2,kz1) = derv2(jx2,kz1) - xd1*zd2*d2(ind)*one2
                derv2(jy2,kx1) = derv2(jy2,kx1) - yd1*xd2*d2(ind)*one2
                derv2(jy2,ky1) = derv2(jy2,ky1) - yd1*yd2*d2(ind)*one2
                derv2(jy2,kz1) = derv2(jy2,kz1) - yd1*zd2*d2(ind)*one2
                derv2(jz2,kx1) = derv2(jz2,kx1) - zd1*xd2*d2(ind)*one2
                derv2(jz2,ky1) = derv2(jz2,ky1) - zd1*yd2*d2(ind)*one2
                derv2(jz2,kz1) = derv2(jz2,kz1) - zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(jx2,kx1) = derv2(jx2,kx1) - d1(nf2)*one2
                  derv2(jy2,ky1) = derv2(jy2,ky1) - d1(nf2)*one2
                  derv2(jz2,kz1) = derv2(jz2,kz1) - d1(nf2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if ((lopj.or.lopl).and.n1j.ne.n2l) then
              if (lopj.and.lopl) then
                ldo21 = ljlocal
                ldo12 = lllocal
                if (ldo21) then
                  jx1 = jx
                  jy1 = jy
                  jz1 = jz
                  lx2 = lxf
                  ly2 = lyf
                  lz2 = lzf
                endif
                if (ldo12) then
                  jx2 = jxf
                  jy2 = jyf
                  jz2 = jzf
                  lx1 = lx
                  ly1 = ly
                  lz1 = lz
                endif
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                ldo21 = ljlocal
                ldo12 = ljlocal
                if (ldo21) then
                  jx1 = jx
                  jy1 = jy
                  jz1 = jz
                  lx2 = jxf
                  ly2 = jyf
                  lz2 = jzf
                endif
                if (ldo12) then
                  jx2 = jxf
                  jy2 = jyf
                  jz2 = jzf
                  lx1 = jx
                  ly1 = jy
                  lz1 = jz
                endif
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                ldo21 = lllocal
                ldo12 = lllocal
                if (ldo21) then
                  jx1 = lx
                  jy1 = ly
                  jz1 = lz
                  lx2 = lxf
                  ly2 = lyf
                  lz2 = lzf
                endif
                if (ldo12) then
                  jx2 = lxf
                  jy2 = lyf
                  jz2 = lzf
                  lx1 = lx
                  ly1 = ly
                  lz1 = lz
                endif
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              if (ldo21) then
                derv2(lx2,jx1) = derv2(lx2,jx1) + xd1*xd2*d2(ind)*one1
                derv2(ly2,jx1) = derv2(ly2,jx1) + xd1*yd2*d2(ind)*one1
                derv2(lz2,jx1) = derv2(lz2,jx1) + xd1*zd2*d2(ind)*one1
                derv2(lx2,jy1) = derv2(lx2,jy1) + yd1*xd2*d2(ind)*one1
                derv2(ly2,jy1) = derv2(ly2,jy1) + yd1*yd2*d2(ind)*one1
                derv2(lz2,jy1) = derv2(lz2,jy1) + yd1*zd2*d2(ind)*one1
                derv2(lx2,jz1) = derv2(lx2,jz1) + zd1*xd2*d2(ind)*one1
                derv2(ly2,jz1) = derv2(ly2,jz1) + zd1*yd2*d2(ind)*one1
                derv2(lz2,jz1) = derv2(lz2,jz1) + zd1*zd2*d2(ind)*one1
                if (n1.eq.n2) then
                  derv2(lx2,jx1) = derv2(lx2,jx1) + d1(nf2)*one1
                  derv2(ly2,jy1) = derv2(ly2,jy1) + d1(nf2)*one1
                  derv2(lz2,jz1) = derv2(lz2,jz1) + d1(nf2)*one1
                endif
              endif
              if (.not.lsameijkl.and.ldo12) then
                derv2(jx2,lx1) = derv2(jx2,lx1) + xd1*xd2*d2(ind)*one2
                derv2(jx2,ly1) = derv2(jx2,ly1) + xd1*yd2*d2(ind)*one2
                derv2(jx2,lz1) = derv2(jx2,lz1) + xd1*zd2*d2(ind)*one2
                derv2(jy2,lx1) = derv2(jy2,lx1) + yd1*xd2*d2(ind)*one2
                derv2(jy2,ly1) = derv2(jy2,ly1) + yd1*yd2*d2(ind)*one2
                derv2(jy2,lz1) = derv2(jy2,lz1) + yd1*zd2*d2(ind)*one2
                derv2(jz2,lx1) = derv2(jz2,lx1) + zd1*xd2*d2(ind)*one2
                derv2(jz2,ly1) = derv2(jz2,ly1) + zd1*yd2*d2(ind)*one2
                derv2(jz2,lz1) = derv2(jz2,lz1) + zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(jx2,lx1) = derv2(jx2,lx1) + d1(nf2)*one2
                  derv2(jy2,ly1) = derv2(jy2,ly1) + d1(nf2)*one2
                  derv2(jz2,lz1) = derv2(jz2,lz1) + d1(nf2)*one2
                endif
              endif
            endif
            if (lstr) then
              if (lilocal) then
!
!  Strain - strain second derivatives
!
                do is1 = 1,nstrains
                  ns1 = nstrptr(is1)
                  do is2 = 1,nstrains
                    ns2 = nstrptr(is2)
                    sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds2(ns2)*dr2ds1(ns1)
                    if (n1.ne.n2) then
                      sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds1(ns2)*dr2ds2(ns1)
                    else
                      sderv2(is2,is1) = sderv2(is2,is1) + d1(nf2)*d2r2ds22(ns2,ns1)
                    endif
                  enddo
                enddo
              endif
!
!  Internal - strain second derivatives
!
              if (lilocal) then
                do is1 = 1,nstrains
                  ns1 = nstrptr(is1)
                  derv3(ix,is1) = derv3(ix,is1) - xd1*d2(ind)*dr2ds2(ns1)
                  derv3(iy,is1) = derv3(iy,is1) - yd1*d2(ind)*dr2ds2(ns1)
                  derv3(iz,is1) = derv3(iz,is1) - zd1*d2(ind)*dr2ds2(ns1)
                enddo
                if (n1.eq.n2) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(ix,is1) = derv3(ix,is1) - d1(nf2)*d2r2dsdx2(ns1,1)
                    derv3(iy,is1) = derv3(iy,is1) - d1(nf2)*d2r2dsdx2(ns1,2)
                    derv3(iz,is1) = derv3(iz,is1) - d1(nf2)*d2r2dsdx2(ns1,3)
                  enddo
                endif
              endif
              if (ljlocal) then
                do is1 = 1,nstrains
                  ns1 = nstrptr(is1)
                  derv3(jx,is1) = derv3(jx,is1) + xd1*d2(ind)*dr2ds2(ns1)
                  derv3(jy,is1) = derv3(jy,is1) + yd1*d2(ind)*dr2ds2(ns1)
                  derv3(jz,is1) = derv3(jz,is1) + zd1*d2(ind)*dr2ds2(ns1)
                enddo
                if (n1.eq.n2) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(jx,is1) = derv3(jx,is1) + d1(nf2)*d2r2dsdx2(ns1,1)
                    derv3(jy,is1) = derv3(jy,is1) + d1(nf2)*d2r2dsdx2(ns1,2)
                    derv3(jz,is1) = derv3(jz,is1) + d1(nf2)*d2r2dsdx2(ns1,3)
                  enddo
                endif
              endif
              if (n1.ne.n2) then
                if (lklocal) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(kx,is1) = derv3(kx,is1) - xd2*d2(ind)*dr2ds1(ns1)
                    derv3(ky,is1) = derv3(ky,is1) - yd2*d2(ind)*dr2ds1(ns1)
                    derv3(kz,is1) = derv3(kz,is1) - zd2*d2(ind)*dr2ds1(ns1)
                  enddo
                endif
                if (lllocal) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(lx,is1) = derv3(lx,is1) + xd2*d2(ind)*dr2ds1(ns1)
                    derv3(ly,is1) = derv3(ly,is1) + yd2*d2(ind)*dr2ds1(ns1)
                    derv3(lz,is1) = derv3(lz,is1) + zd2*d2(ind)*dr2ds1(ns1)
                  enddo
                endif
              endif
            endif
          endif
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('d2adddmsp')
#endif
!
  return
  end
!
  subroutine d1addsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1, &
                     nonzero,nonzeroptr,ltorderv,lextended,ldoregions)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian first derivative arrays.
!  Sparse storage of derivatives version.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  d1              = array of first derivatives w.r.t. interatomic distances
!  nonzero         = number of nonzero elements in d1
!  nonzeroptr(*)   = pointer to nonzero elements in d1
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!  ldoregions      = if true then region derivatives should be computed
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!
!   8/23 Created from d1add and d2addsp
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
  use brennerdata
  use control,        only : latomicstress
  use current,        only : ndim, nrelf2a, nregionno
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nonzero
  integer(i4), intent(in)             :: nonzeroptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1(*)
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
  logical,     intent(in)             :: ldoregions
!
!  Local variables
!
  integer(i4)                         :: io
  integer(i4)                         :: is
  integer(i4)                         :: jo
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nf1
  integer(i4)                         :: nji
  integer(i4)                         :: nkj
  integer(i4)                         :: nlk
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: nregi
  integer(i4)                         :: nregj
  integer(i4)                         :: ns
  integer(i4)                         :: ntmp
  real(dp)                            :: dr2ds(6)
  real(dp)                            :: d2r2dx2(3,3)
  real(dp)                            :: d2r2ds2(6,6)
  real(dp)                            :: d2r2dsdx(6,3)
  real(dp)                            :: xd
  real(dp)                            :: yd
  real(dp)                            :: zd
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
#ifdef TRACE
  call trace_in('d1addsp')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndisttor = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
    if (lextended) then
      ndist = ndisttor + nneigh(nri)*mneigh*mneigh
    else
      ndist = ndisttor
    endif
  else
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = nneigh2 + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
!  Loop over pairs of interatomic distances
!
  do nf1 = 1,nonzero
    n1 = nonzeroptr(nf1)
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Get vector between atoms
!
      xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
      yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
      zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
!  Convert neighbour numbers to real atoms
!
      n1i = neighno(n1i,nri)
      n1j = neighno(n1j,nri)
    elseif (n1.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
      if (n1i.eq.0) then
        xi = 0.0_dp
        yi = 0.0_dp
        zi = 0.0_dp
        n1i = i
      else
        xi = xneigh(n1i,nri)
        yi = yneigh(n1i,nri)
        zi = zneigh(n1i,nri)
        n1i = neighno(n1i,nri)
      endif
!
      xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
      yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
      zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
      n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
    else
!
!  Extended derivatives
!
      ntmp = n1 - ndisttor
      nji = (ntmp - 1)/(mneigh*mneigh) + 1
      ntmp = ntmp - (nji - 1)*mneigh*mneigh
      nkj = (ntmp - 1)/mneigh + 1
      ntmp = ntmp - (nkj - 1)*mneigh
      nlk = ntmp
!
      nrj = neighno(nji,nri)
      if (nkj.gt.nneigh(nrj)) cycle
      nrk = neighno(nkj,nrj)
      if (nlk.gt.nneigh(nrk)) cycle
      n1i = nrk
      n1j = neighno(nlk,nrk)
!
      xd1 = xneigh(nlk,nrk)
      yd1 = yneigh(nlk,nrk)
      zd1 = zneigh(nlk,nrk)
    endif
!
!  Scale coordinate differences by the first derivative
!
    xd = xd1*d1(nf1)
    yd = yd1*d1(nf1)
    zd = zd1*d1(nf1)
!
    io = nfreeatom(n1i)
    jo = nfreeatom(n1j)
    if (io.gt.0) then
      xdrv(n1i) = xdrv(n1i) - xd
      ydrv(n1i) = ydrv(n1i) - yd
      zdrv(n1i) = zdrv(n1i) - zd
    endif 
    if (jo.gt.0) then
      xdrv(n1j) = xdrv(n1j) + xd
      ydrv(n1j) = ydrv(n1j) + yd
      zdrv(n1j) = zdrv(n1j) + zd
    endif 
    if (ldoregions) then
      nregi = nregionno(nrelf2a(n1i))
      nregj = nregionno(nrelf2a(n1j))
      if (nregi.ne.nregj) then
        xregdrv(nregi) = xregdrv(nregi) - xd
        yregdrv(nregi) = yregdrv(nregi) - yd
        zregdrv(nregi) = zregdrv(nregi) - zd
        xregdrv(nregj) = xregdrv(nregj) + xd
        yregdrv(nregj) = yregdrv(nregj) + yd
        zregdrv(nregj) = zregdrv(nregj) + zd
      endif 
    endif   
    if (lstr) then
      call real1strterm(ndim,xd1,yd1,zd1,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
      do is = 1,nstrains
        ns = nstrptr(is)
        rstrd(is) = rstrd(is) + d1(nf1)*dr2ds(ns)
      enddo 
      if (latomicstress) then
        do is = 1,nstrains
          ns = nstrptr(is)
          atomicstress(is,n1i) = atomicstress(is,n1i) + 0.5_dp*d1(nf1)*dr2ds(ns)
          atomicstress(is,n1j) = atomicstress(is,n1j) + 0.5_dp*d1(nf1)*dr2ds(ns)
        enddo 
      endif 
    endif
  enddo
#ifdef TRACE
  call trace_out('d1addsp')
#endif
!
  return
  end
!
  subroutine d2addspp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nREBOatomRptr,d1,d2, &
                      xkv,ykv,zkv,nonzero,nonzeroptr,ltorderv,lextended)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array.
!  Phonon version.
!  Sparse storage of derivatives version.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic distances in sparse form
!  d2              = array of first derivatives w.r.t. interatomic distances in sparse form
!  nonzero         = number of nonzero elements in d1
!  nonzeroptr(*)   = pointer to nonzero elements in d1
!  xkv             = x component of K vector
!  ykv             = y component of K vector
!  zkv             = z component of K vector
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!
!   5/23 Created from d2addsp and d2addp
!   8/23 d2 changed to sparse format and lextended added
!   8/23 Order of arguments changed to be more consistent with other routines
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
  use brennerdata
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nonzero
  integer(i4), intent(in)             :: nonzeroptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  real(dp),    intent(in)             :: xkv
  real(dp),    intent(in)             :: ykv
  real(dp),    intent(in)             :: zkv
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ind0
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nf1
  integer(i4)                         :: nf2
  integer(i4)                         :: nji
  integer(i4)                         :: nkj
  integer(i4)                         :: nlk
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: ntmp
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: lsameijkl
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
#ifdef TRACE
  call trace_in('d2addspp')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndisttor = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
  else
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = nneigh2 + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
!  Loop over pairs of interatomic distances
!
  do nf1 = 1,nonzero
    n1 = nonzeroptr(nf1)
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Check that neighbour numbers are valid
!
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    elseif (n1.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
!  Check that atoms are valid
!
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        if (nREBOatomRptr(neighno(nji,nri)).eq.0) then
          ljok = .false.
        else
          ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
        endif
      else
        ljok = .false.
      endif
!
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!         
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    else
!
!  Extended derivatives
!
      ntmp = n1 - ndisttor
      nji = (ntmp - 1)/(mneigh*mneigh) + 1
      ntmp = ntmp - (nji - 1)*mneigh*mneigh
      nkj = (ntmp - 1)/mneigh + 1
      ntmp = ntmp - (nkj - 1)*mneigh
      nlk = ntmp
!
      nrj = neighno(nji,nri)
      if (nrj.eq.0) cycle
      if (nkj.gt.nneigh(nrj)) cycle
      nrk = neighno(nkj,nrj)
      if (nrk.eq.0) cycle
      if (nlk.gt.nneigh(nrk)) cycle
      n1i = nrk
      n1j = neighno(nlk,nrk)
      if (n1j.eq.0) cycle
!   
!  Check that atoms are valid
!  
      liok = (nkj.le.nneigh(nrj))
      ljok = (nlk.le.nneigh(nrk))
!
      if (liok.and.ljok) then
        x2i = xneigh(nji,nri) + xneigh(nkj,nrj)
        y2i = yneigh(nji,nri) + yneigh(nkj,nrj)
        z2i = zneigh(nji,nri) + zneigh(nkj,nrj)
!
        xd1 = xneigh(nlk,nrk)
        yd1 = yneigh(nlk,nrk)
        zd1 = zneigh(nlk,nrk)
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
!
!  Set second derivative location pointers
!
      ix = 3*(n1i - 1) + 1
      iy = ix + 1
      iz = iy + 1
      jx = 3*(n1j - 1) + 1
      jy = jx + 1
      jz = jy + 1
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
      ind0 = nf1*(nf1-1)/2
      do nf2 = 1,nf1
        n2 = nonzeroptr(nf2)
        ind = ind0 + nf2
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(nf2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!  
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          elseif (n2.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              if (nREBOatomRptr(neighno(nji,nri)).eq.0) then
                llok = .false.
              else
                llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
              endif
            else
              llok = .false.
            endif
!
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!               
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          else
! 
!  Extended derivatives               
! 
            ntmp = n2 - ndisttor
            nji = (ntmp - 1)/(mneigh*mneigh) + 1
            ntmp = ntmp - (nji - 1)*mneigh*mneigh
            nkj = (ntmp - 1)/mneigh + 1
            ntmp = ntmp - (nkj - 1)*mneigh
            nlk = ntmp
! 
            nrj = neighno(nji,nri)
            if (nrj.eq.0) cycle
            if (nkj.gt.nneigh(nrj)) cycle
            nrk = neighno(nkj,nrj)
            if (nrk.eq.0) cycle
            if (nlk.gt.nneigh(nrk)) cycle
            n2k = nrk
            n2l = neighno(nlk,nrk)
            if (n2l.eq.0) cycle
!   
!  Check that atoms are valid         
!  
            lkok = (nkj.le.nneigh(nrj))
            llok = (nlk.le.nneigh(nrk))
! 
            if (liok.and.ljok) then         
              xik = xneigh(nji,nri) + xneigh(nkj,nrj) - x2i
              yik = yneigh(nji,nri) + yneigh(nkj,nrj) - y2i
              zik = zneigh(nji,nri) + zneigh(nkj,nrj) - z2i
! 
              xd2 = xneigh(nlk,nrk)
              yd2 = yneigh(nlk,nrk)
              zd2 = zneigh(nlk,nrk)
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
            kx = 3*(n2k - 1) + 1
            ky = kx + 1
            kz = ky + 1
            lx = 3*(n2l - 1) + 1
            ly = lx + 1
            lz = ly + 1
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if (n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              derv2(kx,ix) = derv2(kx,ix) + xd1*xd2*d2(ind)
              derv2(ky,ix) = derv2(ky,ix) + xd1*yd2*d2(ind)
              derv2(kz,ix) = derv2(kz,ix) + xd1*zd2*d2(ind)
              derv2(kx,iy) = derv2(kx,iy) + yd1*xd2*d2(ind)
              derv2(ky,iy) = derv2(ky,iy) + yd1*yd2*d2(ind)
              derv2(kz,iy) = derv2(kz,iy) + yd1*zd2*d2(ind)
              derv2(kx,iz) = derv2(kx,iz) + zd1*xd2*d2(ind)
              derv2(ky,iz) = derv2(ky,iz) + zd1*yd2*d2(ind)
              derv2(kz,iz) = derv2(kz,iz) + zd1*zd2*d2(ind)
              if (n1.eq.n2) then
                derv2(kx,ix) = derv2(kx,ix) + d1(nf2)
                derv2(ky,iy) = derv2(ky,iy) + d1(nf2)
                derv2(kz,iz) = derv2(kz,iz) + d1(nf2)
              endif
              if (.not.lsameijkl) then
                derv2(ix,kx) = derv2(ix,kx) + xd1*xd2*d2(ind)
                derv2(ix,ky) = derv2(ix,ky) + xd1*yd2*d2(ind)
                derv2(ix,kz) = derv2(ix,kz) + xd1*zd2*d2(ind)
                derv2(iy,kx) = derv2(iy,kx) + yd1*xd2*d2(ind)
                derv2(iy,ky) = derv2(iy,ky) + yd1*yd2*d2(ind)
                derv2(iy,kz) = derv2(iy,kz) + yd1*zd2*d2(ind)
                derv2(iz,kx) = derv2(iz,kx) + zd1*xd2*d2(ind)
                derv2(iz,ky) = derv2(iz,ky) + zd1*yd2*d2(ind)
                derv2(iz,kz) = derv2(iz,kz) + zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(ix,kx) = derv2(ix,kx) + d1(nf2)
                  derv2(iy,ky) = derv2(iy,ky) + d1(nf2)
                  derv2(iz,kz) = derv2(iz,kz) + d1(nf2)
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if (n1i.ne.n2l) then
              derv2(lx,ix) = derv2(lx,ix) - xd1*xd2*d2(ind)
              derv2(ly,ix) = derv2(ly,ix) - xd1*yd2*d2(ind)
              derv2(lz,ix) = derv2(lz,ix) - xd1*zd2*d2(ind)
              derv2(lx,iy) = derv2(lx,iy) - yd1*xd2*d2(ind)
              derv2(ly,iy) = derv2(ly,iy) - yd1*yd2*d2(ind)
              derv2(lz,iy) = derv2(lz,iy) - yd1*zd2*d2(ind)
              derv2(lx,iz) = derv2(lx,iz) - zd1*xd2*d2(ind)
              derv2(ly,iz) = derv2(ly,iz) - zd1*yd2*d2(ind)
              derv2(lz,iz) = derv2(lz,iz) - zd1*zd2*d2(ind)
              if (n1.eq.n2) then
                derv2(lx,ix) = derv2(lx,ix) - d1(nf2)
                derv2(ly,iy) = derv2(ly,iy) - d1(nf2)
                derv2(lz,iz) = derv2(lz,iz) - d1(nf2)
              endif
              if (.not.lsameijkl) then
                derv2(ix,lx) = derv2(ix,lx) - xd1*xd2*d2(ind)
                derv2(ix,ly) = derv2(ix,ly) - xd1*yd2*d2(ind)
                derv2(ix,lz) = derv2(ix,lz) - xd1*zd2*d2(ind)
                derv2(iy,lx) = derv2(iy,lx) - yd1*xd2*d2(ind)
                derv2(iy,ly) = derv2(iy,ly) - yd1*yd2*d2(ind)
                derv2(iy,lz) = derv2(iy,lz) - yd1*zd2*d2(ind)
                derv2(iz,lx) = derv2(iz,lx) - zd1*xd2*d2(ind)
                derv2(iz,ly) = derv2(iz,ly) - zd1*yd2*d2(ind)
                derv2(iz,lz) = derv2(iz,lz) - zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(ix,lx) = derv2(ix,lx) - d1(nf2)
                  derv2(iy,ly) = derv2(iy,ly) - d1(nf2)
                  derv2(iz,lz) = derv2(iz,lz) - d1(nf2)
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if (n1j.ne.n2k) then
              derv2(kx,jx) = derv2(kx,jx) - xd1*xd2*d2(ind)
              derv2(ky,jx) = derv2(ky,jx) - xd1*yd2*d2(ind)
              derv2(kz,jx) = derv2(kz,jx) - xd1*zd2*d2(ind)
              derv2(kx,jy) = derv2(kx,jy) - yd1*xd2*d2(ind)
              derv2(ky,jy) = derv2(ky,jy) - yd1*yd2*d2(ind)
              derv2(kz,jy) = derv2(kz,jy) - yd1*zd2*d2(ind)
              derv2(kx,jz) = derv2(kx,jz) - zd1*xd2*d2(ind)
              derv2(ky,jz) = derv2(ky,jz) - zd1*yd2*d2(ind)
              derv2(kz,jz) = derv2(kz,jz) - zd1*zd2*d2(ind)
              if (n1.eq.n2) then
                derv2(kx,jx) = derv2(kx,jx) - d1(nf2)
                derv2(ky,jy) = derv2(ky,jy) - d1(nf2)
                derv2(kz,jz) = derv2(kz,jz) - d1(nf2)
              endif
              if (.not.lsameijkl) then
                derv2(jx,kx) = derv2(jx,kx) - xd1*xd2*d2(ind)
                derv2(jx,ky) = derv2(jx,ky) - xd1*yd2*d2(ind)
                derv2(jx,kz) = derv2(jx,kz) - xd1*zd2*d2(ind)
                derv2(jy,kx) = derv2(jy,kx) - yd1*xd2*d2(ind)
                derv2(jy,ky) = derv2(jy,ky) - yd1*yd2*d2(ind)
                derv2(jy,kz) = derv2(jy,kz) - yd1*zd2*d2(ind)
                derv2(jz,kx) = derv2(jz,kx) - zd1*xd2*d2(ind)
                derv2(jz,ky) = derv2(jz,ky) - zd1*yd2*d2(ind)
                derv2(jz,kz) = derv2(jz,kz) - zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(jx,kx) = derv2(jx,kx) - d1(nf2)
                  derv2(jy,ky) = derv2(jy,ky) - d1(nf2)
                  derv2(jz,kz) = derv2(jz,kz) - d1(nf2)
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if (n1j.ne.n2l) then
              derv2(lx,jx) = derv2(lx,jx) + xd1*xd2*d2(ind)
              derv2(ly,jx) = derv2(ly,jx) + xd1*yd2*d2(ind)
              derv2(lz,jx) = derv2(lz,jx) + xd1*zd2*d2(ind)
              derv2(lx,jy) = derv2(lx,jy) + yd1*xd2*d2(ind)
              derv2(ly,jy) = derv2(ly,jy) + yd1*yd2*d2(ind)
              derv2(lz,jy) = derv2(lz,jy) + yd1*zd2*d2(ind)
              derv2(lx,jz) = derv2(lx,jz) + zd1*xd2*d2(ind)
              derv2(ly,jz) = derv2(ly,jz) + zd1*yd2*d2(ind)
              derv2(lz,jz) = derv2(lz,jz) + zd1*zd2*d2(ind)
              if (n1.eq.n2) then
                derv2(lx,jx) = derv2(lx,jx) + d1(nf2)
                derv2(ly,jy) = derv2(ly,jy) + d1(nf2)
                derv2(lz,jz) = derv2(lz,jz) + d1(nf2)
              endif
              if (.not.lsameijkl) then
                derv2(jx,lx) = derv2(jx,lx) + xd1*xd2*d2(ind)
                derv2(jx,ly) = derv2(jx,ly) + xd1*yd2*d2(ind)
                derv2(jx,lz) = derv2(jx,lz) + xd1*zd2*d2(ind)
                derv2(jy,lx) = derv2(jy,lx) + yd1*xd2*d2(ind)
                derv2(jy,ly) = derv2(jy,ly) + yd1*yd2*d2(ind)
                derv2(jy,lz) = derv2(jy,lz) + yd1*zd2*d2(ind)
                derv2(jz,lx) = derv2(jz,lx) + zd1*xd2*d2(ind)
                derv2(jz,ly) = derv2(jz,ly) + zd1*yd2*d2(ind)
                derv2(jz,lz) = derv2(jz,lz) + zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(jx,lx) = derv2(jx,lx) + d1(nf2)
                  derv2(jy,ly) = derv2(jy,ly) + d1(nf2)
                  derv2(jz,lz) = derv2(jz,lz) + d1(nf2)
                endif
              endif
            endif
          endif
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addspp')
#endif
!
  return
  end
!
  subroutine d2addmixdp(lilocal,lmlocal,lnlocal,i,ix,iy,iz,ixf,iyf,izf,maxneigh,mneigh,nneigh,neighno, &
                        xneigh,yneigh,zneigh,nREBOatomRptr,d2,mx,my,mz,mxf,myf,mzf, &
                        nx,ny,nz,nxf,nyf,nzf,xmn,ymn,zmn,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivatives. Phonon version.
!  This routines handles a combination of the derivatives of a single vector 
!  and derivatives stored in d1 format (but named d2 since these are 2nd
!  derivatives between the two distances).
!
!  On entry : 
!
!  lilocal         = if .true. then i is local to this processor
!  lmlocal         = if .true. then m is local to this processor
!  lnlocal         = if .true. then n is local to this processor
!  i               = atom number of i to whom others are connected
!  ix(f)           = position of x coordinate for i in second derivative matrix (f for full; otherwise local)
!  iy(f)           = position of y coordinate for i in second derivative matrix (f for full; otherwise local)
!  iz(f)           = position of z coordinate for i in second derivative matrix (f for full; otherwise local)
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  d2              = sparse array of second derivatives w.r.t. two interatomic distances
!  ltorderv        = if true then include torsional derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  mx(f)           = position of x coordinate for m in second derivative matrix (f for full; otherwise local)
!  my(f)           = position of y coordinate for m in second derivative matrix (f for full; otherwise local)
!  mz(f)           = position of z coordinate for m in second derivative matrix (f for full; otherwise local)
!  nx(f)           = position of x coordinate for n in second derivative matrix (f for full; otherwise local)
!  ny(f)           = position of y coordinate for n in second derivative matrix (f for full; otherwise local)
!  nz(f)           = position of z coordinate for n in second derivative matrix (f for full; otherwise local)
!  xmn, ymn, zmn   = Cartesian components of vector from m to n
!
!   5/23 Created from d2addmixdspp
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use derivatives
  use parallel,       only : atom2local
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: ix
  integer(i4), intent(in)             :: iy
  integer(i4), intent(in)             :: iz
  integer(i4), intent(in)             :: ixf
  integer(i4), intent(in)             :: iyf
  integer(i4), intent(in)             :: izf
  integer(i4), intent(in)             :: mx
  integer(i4), intent(in)             :: my
  integer(i4), intent(in)             :: mz
  integer(i4), intent(in)             :: mxf
  integer(i4), intent(in)             :: myf
  integer(i4), intent(in)             :: mzf
  integer(i4), intent(in)             :: nx
  integer(i4), intent(in)             :: ny
  integer(i4), intent(in)             :: nz
  integer(i4), intent(in)             :: nxf
  integer(i4), intent(in)             :: nyf
  integer(i4), intent(in)             :: nzf
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  logical,     intent(in)             :: lilocal
  logical,     intent(in)             :: lmlocal
  logical,     intent(in)             :: lnlocal
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  real(dp),    intent(in)             :: xmn
  real(dp),    intent(in)             :: ymn
  real(dp),    intent(in)             :: zmn
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: indn2n
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: kloc
  integer(i4)                         :: kx
  integer(i4)                         :: ky
  integer(i4)                         :: kz
  integer(i4)                         :: kxf
  integer(i4)                         :: kyf
  integer(i4)                         :: kzf
  integer(i4)                         :: l
  integer(i4)                         :: lloc
  integer(i4)                         :: lx
  integer(i4)                         :: ly
  integer(i4)                         :: lz
  integer(i4)                         :: lxf
  integer(i4)                         :: lyf
  integer(i4)                         :: lzf
  integer(i4)                         :: n1j2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: ntmp
  logical                             :: lklocal
  logical                             :: lllocal
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xkl
  real(dp)                            :: ykl
  real(dp)                            :: zkl
#ifdef TRACE
  call trace_in('d2addmixdspp')
#endif
!
!  Set up terms for i
!
  nri = nREBOatomRptr(i)
!
!  Set up indices for boundaries
!
  indn2n = nneigh(nri) + nneigh(nri)*(nneigh(nri)+1)/2
  if (ltorderv) then
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh
    ndist = ndisttor
  else
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
    ndist = indn2n
  endif
!
  do ind = 1,ndist
    if (abs(d2(ind)).gt.d2tol) then
      l = ind
      if (l.le.nneigh(nri)) then
!
!  i -> neighbours
!
        k = neighno(l,nri)
        xik = xneigh(l,nri)
        yik = yneigh(l,nri)
        zik = zneigh(l,nri)
!
        kloc = atom2local(k)
        lklocal = (kloc.gt.0)
!
        kx = 3*(kloc-1) + 1
        ky = kx + 1
        kz = kx + 2
!
        kxf = 3*(k-1) + 1
        kyf = kxf + 1
        kzf = kxf + 2
!
        call add_drv2_2pd(lilocal,lklocal,lmlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xik,yik,zik,xmn,ymn,zmn,d2(ind))
      elseif (l.le.indn2n) then
!
!  Loop over neighbours -> neighbours
!
        ntmp = l - nneigh(nri)
        nk = 1
        n1j2 = 1
        do while (n1j2.lt.ntmp)
          nk = nk + 1
          n1j2 = n1j2 + nk
        enddo
        n1j2 = n1j2 - nk
        nl = ntmp - n1j2
!
        k = neighno(nk,nri)
        l = neighno(nl,nri)
!
        kloc = atom2local(k)
        lloc = atom2local(l)
        lklocal = (kloc.gt.0)
        lllocal = (lloc.gt.0)
!
        xkl = (xneigh(nk,nri) - xneigh(nl,nri))
        ykl = (yneigh(nk,nri) - yneigh(nl,nri))
        zkl = (zneigh(nk,nri) - zneigh(nl,nri))
!
        kx = 3*(kloc-1) + 1
        ky = kx + 1
        kz = kx + 2
!
        kxf = 3*(k-1) + 1
        kyf = kxf + 1
        kzf = kxf + 2
!
        lx = 3*(lloc-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        lxf = 3*(l-1) + 1
        lyf = lxf + 1
        lzf = lxf + 2
!
        call add_drv2_2pd(lklocal,lllocal,lmlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind))
      elseif (l.le.ndisttor) then
        ntmp = l - indn2n
        nj = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
        ntmp = ntmp - (nj - 1)*(mneigh + 1)*mneigh
        nk = (ntmp - 1)/mneigh
        nl = ntmp - nk*mneigh
!
        if (nk.eq.0) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
          j = neighno(nj,nri)
          nrj = nREBOatomRptr(j)
          l = neighno(nl,nrj)
!
          lloc = atom2local(l)
          lllocal = (lloc.gt.0)
!
          xil = xneigh(nl,nrj) + xneigh(nj,nri)
          yil = yneigh(nl,nrj) + yneigh(nj,nri)
          zil = zneigh(nl,nrj) + zneigh(nj,nri)
!
          lx = 3*(lloc-1) + 1
          ly = lx + 1
          lz = lx + 2
!
          lxf = 3*(l-1) + 1
          lyf = lxf + 1
          lzf = lxf + 2
!
          call add_drv2_2pd(lilocal,lllocal,lmlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                            mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xil,yil,zil,xmn,ymn,zmn,d2(ind))
        else
!
!  Loop over neighbours -> neighbours of neighbours
!
          j = neighno(nj,nri)
          nrj = nREBOatomRptr(j)
          k = neighno(nk,nri)
          l = neighno(nl,nrj)
!
          kloc = atom2local(k)
          lloc = atom2local(l)
          lklocal = (kloc.gt.0)
          lllocal = (lloc.gt.0)
!
          xkl = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
          ykl = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
          zkl = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
!
          kx = 3*(kloc-1) + 1
          ky = kx + 1
          kz = kx + 2
!
          kxf = 3*(k-1) + 1
          kyf = kxf + 1
          kzf = kxf + 2
!
          lx = 3*(lloc-1) + 1
          ly = lx + 1
          lz = lx + 2
!
          lxf = 3*(l-1) + 1
          lyf = lxf + 1
          lzf = lxf + 2
!
          call add_drv2_2pd(lklocal,lllocal,lmlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                            mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind))
        endif
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addmixdp')
#endif
!
  return
  end
!
  subroutine d2addmixdspp(lilocal,lmlocal,lnlocal,i,ix,iy,iz,ixf,iyf,izf,maxneigh,mneigh,nneigh,neighno, &
                          xneigh,yneigh,zneigh,nREBOatomRptr,d2,nonzero,nonzeroptr,mx,my,mz,mxf,myf,mzf, &
                          nx,ny,nz,nxf,nyf,nzf,xmn,ymn,zmn,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivatives. Phonon version.
!  This routines handles a combination of the derivatives of a single vector 
!  and derivatives stored in d1 format (but named d2 since these are 2nd
!  derivatives between the two distances).
!  Version that uses sparsity in d2 with distributed parallel algorithm
!
!  On entry : 
!
!  lilocal         = if .true. then i is local to this processor
!  lmlocal         = if .true. then m is local to this processor
!  lnlocal         = if .true. then n is local to this processor
!  i               = atom number of i to whom others are connected
!  ix(f)           = position of x coordinate for i in second derivative matrix (f for full; otherwise local)
!  iy(f)           = position of y coordinate for i in second derivative matrix (f for full; otherwise local)
!  iz(f)           = position of z coordinate for i in second derivative matrix (f for full; otherwise local)
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  d2              = sparse array of second derivatives w.r.t. two interatomic distances
!  nonzero         = number of nonzero elements in d2
!  nonzeroptr(*)   = pointer to nonzero elements in d2
!  ltorderv        = if true then include torsional derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  mx(f)           = position of x coordinate for m in second derivative matrix (f for full; otherwise local)
!  my(f)           = position of y coordinate for m in second derivative matrix (f for full; otherwise local)
!  mz(f)           = position of z coordinate for m in second derivative matrix (f for full; otherwise local)
!  nx(f)           = position of x coordinate for n in second derivative matrix (f for full; otherwise local)
!  ny(f)           = position of y coordinate for n in second derivative matrix (f for full; otherwise local)
!  nz(f)           = position of z coordinate for n in second derivative matrix (f for full; otherwise local)
!  xmn, ymn, zmn   = Cartesian components of vector from m to n
!
!   5/23 Created from d2addmixdsp and d2addmixspp
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
  use brennerdata
  use derivatives
  use parallel,       only : atom2local
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: ix
  integer(i4), intent(in)             :: iy
  integer(i4), intent(in)             :: iz
  integer(i4), intent(in)             :: ixf
  integer(i4), intent(in)             :: iyf
  integer(i4), intent(in)             :: izf
  integer(i4), intent(in)             :: mx
  integer(i4), intent(in)             :: my
  integer(i4), intent(in)             :: mz
  integer(i4), intent(in)             :: mxf
  integer(i4), intent(in)             :: myf
  integer(i4), intent(in)             :: mzf
  integer(i4), intent(in)             :: nx
  integer(i4), intent(in)             :: ny
  integer(i4), intent(in)             :: nz
  integer(i4), intent(in)             :: nxf
  integer(i4), intent(in)             :: nyf
  integer(i4), intent(in)             :: nzf
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nonzero
  integer(i4), intent(in)             :: nonzeroptr(*)
  logical,     intent(in)             :: lilocal
  logical,     intent(in)             :: lmlocal
  logical,     intent(in)             :: lnlocal
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  real(dp),    intent(in)             :: xmn
  real(dp),    intent(in)             :: ymn
  real(dp),    intent(in)             :: zmn
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: indn2n
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: kloc
  integer(i4)                         :: kx
  integer(i4)                         :: ky
  integer(i4)                         :: kz
  integer(i4)                         :: kxf
  integer(i4)                         :: kyf
  integer(i4)                         :: kzf
  integer(i4)                         :: l
  integer(i4)                         :: lloc
  integer(i4)                         :: lx
  integer(i4)                         :: ly
  integer(i4)                         :: lz
  integer(i4)                         :: lxf
  integer(i4)                         :: lyf
  integer(i4)                         :: lzf
  integer(i4)                         :: n1j2
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: ntmp
  logical                             :: lklocal
  logical                             :: lllocal
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xkl
  real(dp)                            :: ykl
  real(dp)                            :: zkl
#ifdef TRACE
  call trace_in('d2addmixdspp')
#endif
!
!  Set up terms for i
!
  nri = nREBOatomRptr(i)
!
!  Set up indices for boundaries
!
  indn2n = nneigh(nri) + nneigh(nri)*(nneigh(nri)+1)/2
  if (ltorderv) then
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh
  else
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
  do ind = 1,nonzero
    l = nonzeroptr(ind)
    if (l.le.nneigh(nri)) then
!
!  i -> neighbours
!
      k = neighno(l,nri)
      xik = xneigh(l,nri)
      yik = yneigh(l,nri)
      zik = zneigh(l,nri)
!
      kloc = atom2local(k)
      lklocal = (kloc.gt.0)
!
      kx = 3*(kloc-1) + 1
      ky = kx + 1
      kz = kx + 2
!
      kxf = 3*(k-1) + 1
      kyf = kxf + 1
      kzf = kxf + 2
!
      call add_drv2_2pd(lilocal,lklocal,lmlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xik,yik,zik,xmn,ymn,zmn,d2(ind))
    elseif (l.le.indn2n) then
!
!  Loop over neighbours -> neighbours
!
      ntmp = l - nneigh(nri)
      nk = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        nk = nk + 1
        n1j2 = n1j2 + nk
      enddo
      n1j2 = n1j2 - nk
      nl = ntmp - n1j2
!
      k = neighno(nk,nri)
      l = neighno(nl,nri)
!
      kloc = atom2local(k)
      lloc = atom2local(l)
      lklocal = (kloc.gt.0)
      lllocal = (lloc.gt.0)
!
      xkl = (xneigh(nk,nri) - xneigh(nl,nri))
      ykl = (yneigh(nk,nri) - yneigh(nl,nri))
      zkl = (zneigh(nk,nri) - zneigh(nl,nri))
!
      kx = 3*(kloc-1) + 1
      ky = kx + 1
      kz = kx + 2
!
      kxf = 3*(k-1) + 1
      kyf = kxf + 1
      kzf = kxf + 2
!
      lx = 3*(lloc-1) + 1
      ly = lx + 1
      lz = lx + 2
!
      lxf = 3*(l-1) + 1
      lyf = lxf + 1
      lzf = lxf + 2
!
      call add_drv2_2pd(lklocal,lllocal,lmlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                        mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind))
    elseif (l.le.ndisttor) then
      ntmp = l - indn2n
      nj = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nj - 1)*(mneigh + 1)*mneigh
      nk = (ntmp - 1)/mneigh
      nl = ntmp - nk*mneigh
!
      if (nk.eq.0) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        l = neighno(nl,nrj)
!
        lloc = atom2local(l)
        lllocal = (lloc.gt.0)
!
        xil = xneigh(nl,nrj) + xneigh(nj,nri)
        yil = yneigh(nl,nrj) + yneigh(nj,nri)
        zil = zneigh(nl,nrj) + zneigh(nj,nri)
!
        lx = 3*(lloc-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        lxf = 3*(l-1) + 1
        lyf = lxf + 1
        lzf = lxf + 2
!
        call add_drv2_2pd(lilocal,lllocal,lmlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                          mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xil,yil,zil,xmn,ymn,zmn,d2(ind))
      else
!
!  Loop over neighbours -> neighbours of neighbours
!
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        k = neighno(nk,nri)
        l = neighno(nl,nrj)
!
        kloc = atom2local(k)
        lloc = atom2local(l)
        lklocal = (kloc.gt.0)
        lllocal = (lloc.gt.0)
!
        xkl = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
        ykl = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
        zkl = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
!
        kx = 3*(kloc-1) + 1
        ky = kx + 1
        kz = kx + 2
!
        kxf = 3*(k-1) + 1
        kyf = kxf + 1
        kzf = kxf + 2
!
        lx = 3*(lloc-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        lxf = 3*(l-1) + 1
        lyf = lxf + 1
        lzf = lxf + 2
!
        call add_drv2_2pd(lklocal,lllocal,lmlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind))
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addmixdspp')
#endif
!
  return
  end
!
  subroutine d2addsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1,d2, &
                     nonzero,nonzeroptr,ltorderv,lextended)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array.
!  Sparse storage of derivatives version.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic distances in sparse form
!  d2              = array of first derivatives w.r.t. interatomic distances in sparse form
!  nonzero         = number of nonzero elements in d1
!  nonzeroptr(*)   = pointer to nonzero elements in d1
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!
!   5/23 Created from d2add
!   8/23 d2 changed to sparse format and lextended added
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
  use brennerdata
  use current,        only : ndim, nstrains
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nonzero
  integer(i4), intent(in)             :: nonzeroptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ind0
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: is1,is2
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: ixc,iyc,izc
  integer(i4)                         :: jxc,jyc,jzc
  integer(i4)                         :: kxc,kyc,kzc
  integer(i4)                         :: lxc,lyc,lzc
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nf1
  integer(i4)                         :: nf2
  integer(i4)                         :: nji
  integer(i4)                         :: nkj
  integer(i4)                         :: nlk
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: ns1
  integer(i4)                         :: ns2
  integer(i4)                         :: ntmp
  real(dp)                            :: one1
  real(dp)                            :: one2
  real(dp)                            :: dr2ds1(6)
  real(dp)                            :: dr2ds2(6)
  real(dp)                            :: d2r2dx21(3,3)
  real(dp)                            :: d2r2dx22(3,3)
  real(dp)                            :: d2r2ds21(6,6)
  real(dp)                            :: d2r2ds22(6,6)
  real(dp)                            :: d2r2dsdx1(6,3)
  real(dp)                            :: d2r2dsdx2(6,3)
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: lsameijkl
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lopi
  logical                             :: lopj
  logical                             :: lopk
  logical                             :: lopl
#ifdef TRACE
  call trace_in('d2addsp')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndisttor = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
    if (lextended) then
      ndist = ndisttor + nneigh(nri)*mneigh*mneigh
    else
      ndist = ndisttor
    endif
  else
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = nneigh2 + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
!  Loop over pairs of interatomic distances
!
  do nf1 = 1,nonzero
    n1 = nonzeroptr(nf1)
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Check that neighbour numbers are valid
!
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    elseif (n1.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
!  Check that atoms are valid
!
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        if (nREBOatomRptr(neighno(nji,nri)).eq.0) then
          ljok = .false.
        else
          ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
        endif
      else
        ljok = .false.
      endif
!
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!         
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    else
!
!  Extended derivatives
!
      ntmp = n1 - ndisttor
      nji = (ntmp - 1)/(mneigh*mneigh) + 1
      ntmp = ntmp - (nji - 1)*mneigh*mneigh
      nkj = (ntmp - 1)/mneigh + 1
      ntmp = ntmp - (nkj - 1)*mneigh
      nlk = ntmp
!
      nrj = neighno(nji,nri)
      if (nrj.eq.0) cycle
      if (nkj.gt.nneigh(nrj)) cycle
      nrk = neighno(nkj,nrj)
      if (nrk.eq.0) cycle
      if (nlk.gt.nneigh(nrk)) cycle
      n1i = nrk
      n1j = neighno(nlk,nrk)
      if (n1j.eq.0) cycle
!   
!  Check that atoms are valid
!  
      liok = (nkj.le.nneigh(nrj))
      ljok = (nlk.le.nneigh(nrk))
!
      if (liok.and.ljok) then
        x2i = xneigh(nji,nri) + xneigh(nkj,nrj)
        y2i = yneigh(nji,nri) + yneigh(nkj,nrj)
        z2i = zneigh(nji,nri) + zneigh(nkj,nrj)
!
        xd1 = xneigh(nlk,nrk)
        yd1 = yneigh(nlk,nrk)
        zd1 = zneigh(nlk,nrk)
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
!
!  Set up products for strain
!
      if (lstr) then
        call real1strterm(ndim,xd1,yd1,zd1,0.0_dp,0.0_dp,0.0_dp,dr2ds1,d2r2dx21,d2r2dsdx1,d2r2ds21,.true.)
      endif
!
!  Set second derivative location pointers
!
      lopi = (nfreeatom(n1i).gt.0)
      lopj = (nfreeatom(n1j).gt.0)
      ix = 3*(nfreeatom(n1i) - 1) + 1
      iy = ix + 1
      iz = iy + 1
      jx = 3*(nfreeatom(n1j) - 1) + 1
      jy = jx + 1
      jz = jy + 1
    endif
!
!  If neither i nor j are being optimised and strain derivatives are not needed then atoms are not valid
!
    if (.not.lopi.and..not.lopj.and..not.lstr) then
      liok = .false.
      ljok = .false.
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
      ind0 = nf1*(nf1-1)/2
      do nf2 = 1,nf1
        n2 = nonzeroptr(nf2)
        ind = ind0 + nf2
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(nf2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!  
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          elseif (n2.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              if (nREBOatomRptr(neighno(nji,nri)).eq.0) then
                llok = .false.
              else
                llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
              endif
            else
              llok = .false.
            endif
!
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!               
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          else
! 
!  Extended derivatives               
! 
            ntmp = n2 - ndisttor
            nji = (ntmp - 1)/(mneigh*mneigh) + 1
            ntmp = ntmp - (nji - 1)*mneigh*mneigh
            nkj = (ntmp - 1)/mneigh + 1
            ntmp = ntmp - (nkj - 1)*mneigh
            nlk = ntmp
! 
            nrj = neighno(nji,nri)
            if (nrj.eq.0) cycle
            if (nkj.gt.nneigh(nrj)) cycle
            nrk = neighno(nkj,nrj)
            if (nrk.eq.0) cycle
            if (nlk.gt.nneigh(nrk)) cycle
            n2k = nrk
            n2l = neighno(nlk,nrk)
            if (n2l.eq.0) cycle
!   
!  Check that atoms are valid         
!  
            lkok = (nkj.le.nneigh(nrj))
            llok = (nlk.le.nneigh(nrk))
! 
            if (liok.and.ljok) then         
              xik = xneigh(nji,nri) + xneigh(nkj,nrj) - x2i
              yik = yneigh(nji,nri) + yneigh(nkj,nrj) - y2i
              zik = zneigh(nji,nri) + zneigh(nkj,nrj) - z2i
! 
              xd2 = xneigh(nlk,nrk)
              yd2 = yneigh(nlk,nrk)
              zd2 = zneigh(nlk,nrk)
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
!  Set up products for strain
!
            if (lstr) then
              call real1strterm(ndim,xd2,yd2,zd2,0.0_dp,0.0_dp,0.0_dp,dr2ds2,d2r2dx22,d2r2dsdx2,d2r2ds22,.true.)
            endif
!
            lopk = (nfreeatom(n2k).gt.0)
            lopl = (nfreeatom(n2l).gt.0)
            kx = 3*(nfreeatom(n2k) - 1) + 1
            ky = kx + 1
            kz = ky + 1
            lx = 3*(nfreeatom(n2l) - 1) + 1
            ly = lx + 1
            lz = ly + 1
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if ((lopi.or.lopk).and.n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              if (lopi.and.lopk) then
                ixc = ix
                iyc = iy
                izc = iz
                kxc = kx
                kyc = ky
                kzc = kz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ixc = ix
                iyc = iy
                izc = iz
                kxc = ix
                kyc = iy
                kzc = iz
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                ixc = kx
                iyc = ky
                izc = kz
                kxc = kx
                kyc = ky
                kzc = kz
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(kxc,ixc) = derv2(kxc,ixc) + xd1*xd2*d2(ind)*one1
              derv2(kyc,ixc) = derv2(kyc,ixc) + xd1*yd2*d2(ind)*one1
              derv2(kzc,ixc) = derv2(kzc,ixc) + xd1*zd2*d2(ind)*one1
              derv2(kxc,iyc) = derv2(kxc,iyc) + yd1*xd2*d2(ind)*one1
              derv2(kyc,iyc) = derv2(kyc,iyc) + yd1*yd2*d2(ind)*one1
              derv2(kzc,iyc) = derv2(kzc,iyc) + yd1*zd2*d2(ind)*one1
              derv2(kxc,izc) = derv2(kxc,izc) + zd1*xd2*d2(ind)*one1
              derv2(kyc,izc) = derv2(kyc,izc) + zd1*yd2*d2(ind)*one1
              derv2(kzc,izc) = derv2(kzc,izc) + zd1*zd2*d2(ind)*one1
              if (n1.eq.n2) then
                derv2(kxc,ixc) = derv2(kxc,ixc) + d1(nf2)*one1
                derv2(kyc,iyc) = derv2(kyc,iyc) + d1(nf2)*one1
                derv2(kzc,izc) = derv2(kzc,izc) + d1(nf2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(ixc,kxc) = derv2(ixc,kxc) + xd1*xd2*d2(ind)*one2
                derv2(ixc,kyc) = derv2(ixc,kyc) + xd1*yd2*d2(ind)*one2
                derv2(ixc,kzc) = derv2(ixc,kzc) + xd1*zd2*d2(ind)*one2
                derv2(iyc,kxc) = derv2(iyc,kxc) + yd1*xd2*d2(ind)*one2
                derv2(iyc,kyc) = derv2(iyc,kyc) + yd1*yd2*d2(ind)*one2
                derv2(iyc,kzc) = derv2(iyc,kzc) + yd1*zd2*d2(ind)*one2
                derv2(izc,kxc) = derv2(izc,kxc) + zd1*xd2*d2(ind)*one2
                derv2(izc,kyc) = derv2(izc,kyc) + zd1*yd2*d2(ind)*one2
                derv2(izc,kzc) = derv2(izc,kzc) + zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(ixc,kxc) = derv2(ixc,kxc) + d1(nf2)*one2
                  derv2(iyc,kyc) = derv2(iyc,kyc) + d1(nf2)*one2
                  derv2(izc,kzc) = derv2(izc,kzc) + d1(nf2)*one2
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if ((lopi.or.lopl).and.n1i.ne.n2l) then
              if (lopi.and.lopl) then
                ixc = ix
                iyc = iy
                izc = iz
                lxc = lx
                lyc = ly
                lzc = lz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ixc = ix
                iyc = iy
                izc = iz
                lxc = ix
                lyc = iy
                lzc = iz
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                ixc = lx
                iyc = ly
                izc = lz
                lxc = lx
                lyc = ly
                lzc = lz
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(lxc,ixc) = derv2(lxc,ixc) - xd1*xd2*d2(ind)*one1
              derv2(lyc,ixc) = derv2(lyc,ixc) - xd1*yd2*d2(ind)*one1
              derv2(lzc,ixc) = derv2(lzc,ixc) - xd1*zd2*d2(ind)*one1
              derv2(lxc,iyc) = derv2(lxc,iyc) - yd1*xd2*d2(ind)*one1
              derv2(lyc,iyc) = derv2(lyc,iyc) - yd1*yd2*d2(ind)*one1
              derv2(lzc,iyc) = derv2(lzc,iyc) - yd1*zd2*d2(ind)*one1
              derv2(lxc,izc) = derv2(lxc,izc) - zd1*xd2*d2(ind)*one1
              derv2(lyc,izc) = derv2(lyc,izc) - zd1*yd2*d2(ind)*one1
              derv2(lzc,izc) = derv2(lzc,izc) - zd1*zd2*d2(ind)*one1
              if (n1.eq.n2) then
                derv2(lxc,ixc) = derv2(lxc,ixc) - d1(nf2)*one1
                derv2(lyc,iyc) = derv2(lyc,iyc) - d1(nf2)*one1
                derv2(lzc,izc) = derv2(lzc,izc) - d1(nf2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(ixc,lxc) = derv2(ixc,lxc) - xd1*xd2*d2(ind)*one2
                derv2(ixc,lyc) = derv2(ixc,lyc) - xd1*yd2*d2(ind)*one2
                derv2(ixc,lzc) = derv2(ixc,lzc) - xd1*zd2*d2(ind)*one2
                derv2(iyc,lxc) = derv2(iyc,lxc) - yd1*xd2*d2(ind)*one2
                derv2(iyc,lyc) = derv2(iyc,lyc) - yd1*yd2*d2(ind)*one2
                derv2(iyc,lzc) = derv2(iyc,lzc) - yd1*zd2*d2(ind)*one2
                derv2(izc,lxc) = derv2(izc,lxc) - zd1*xd2*d2(ind)*one2
                derv2(izc,lyc) = derv2(izc,lyc) - zd1*yd2*d2(ind)*one2
                derv2(izc,lzc) = derv2(izc,lzc) - zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(ixc,lxc) = derv2(ixc,lxc) - d1(nf2)*one2
                  derv2(iyc,lyc) = derv2(iyc,lyc) - d1(nf2)*one2
                  derv2(izc,lzc) = derv2(izc,lzc) - d1(nf2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if ((lopj.or.lopk).and.n1j.ne.n2k) then
              if (lopj.and.lopk) then
                jxc = jx
                jyc = jy
                jzc = jz
                kxc = kx
                kyc = ky
                kzc = kz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                jxc = jx
                jyc = jy
                jzc = jz
                kxc = jx
                kyc = jy
                kzc = jz
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                jxc = kx
                jyc = ky
                jzc = kz
                kxc = kx
                kyc = ky
                kzc = kz
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(kxc,jxc) = derv2(kxc,jxc) - xd1*xd2*d2(ind)*one1
              derv2(kyc,jxc) = derv2(kyc,jxc) - xd1*yd2*d2(ind)*one1
              derv2(kzc,jxc) = derv2(kzc,jxc) - xd1*zd2*d2(ind)*one1
              derv2(kxc,jyc) = derv2(kxc,jyc) - yd1*xd2*d2(ind)*one1
              derv2(kyc,jyc) = derv2(kyc,jyc) - yd1*yd2*d2(ind)*one1
              derv2(kzc,jyc) = derv2(kzc,jyc) - yd1*zd2*d2(ind)*one1
              derv2(kxc,jzc) = derv2(kxc,jzc) - zd1*xd2*d2(ind)*one1
              derv2(kyc,jzc) = derv2(kyc,jzc) - zd1*yd2*d2(ind)*one1
              derv2(kzc,jzc) = derv2(kzc,jzc) - zd1*zd2*d2(ind)*one1
              if (n1.eq.n2) then
                derv2(kxc,jxc) = derv2(kxc,jxc) - d1(nf2)*one1
                derv2(kyc,jyc) = derv2(kyc,jyc) - d1(nf2)*one1
                derv2(kzc,jzc) = derv2(kzc,jzc) - d1(nf2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(jxc,kxc) = derv2(jxc,kxc) - xd1*xd2*d2(ind)*one2
                derv2(jxc,kyc) = derv2(jxc,kyc) - xd1*yd2*d2(ind)*one2
                derv2(jxc,kzc) = derv2(jxc,kzc) - xd1*zd2*d2(ind)*one2
                derv2(jyc,kxc) = derv2(jyc,kxc) - yd1*xd2*d2(ind)*one2
                derv2(jyc,kyc) = derv2(jyc,kyc) - yd1*yd2*d2(ind)*one2
                derv2(jyc,kzc) = derv2(jyc,kzc) - yd1*zd2*d2(ind)*one2
                derv2(jzc,kxc) = derv2(jzc,kxc) - zd1*xd2*d2(ind)*one2
                derv2(jzc,kyc) = derv2(jzc,kyc) - zd1*yd2*d2(ind)*one2
                derv2(jzc,kzc) = derv2(jzc,kzc) - zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(jxc,kxc) = derv2(jxc,kxc) - d1(nf2)*one2
                  derv2(jyc,kyc) = derv2(jyc,kyc) - d1(nf2)*one2
                  derv2(jzc,kzc) = derv2(jzc,kzc) - d1(nf2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if ((lopj.or.lopl).and.n1j.ne.n2l) then
              if (lopj.and.lopl) then
                jxc = jx
                jyc = jy
                jzc = jz
                lxc = lx
                lyc = ly
                lzc = lz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                jxc = jx
                jyc = jy
                jzc = jz
                lxc = jx
                lyc = jy
                lzc = jz
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                jxc = lx
                jyc = ly
                jzc = lz
                lxc = lx
                lyc = ly
                lzc = lz
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(lxc,jxc) = derv2(lxc,jxc) + xd1*xd2*d2(ind)*one1
              derv2(lyc,jxc) = derv2(lyc,jxc) + xd1*yd2*d2(ind)*one1
              derv2(lzc,jxc) = derv2(lzc,jxc) + xd1*zd2*d2(ind)*one1
              derv2(lxc,jyc) = derv2(lxc,jyc) + yd1*xd2*d2(ind)*one1
              derv2(lyc,jyc) = derv2(lyc,jyc) + yd1*yd2*d2(ind)*one1
              derv2(lzc,jyc) = derv2(lzc,jyc) + yd1*zd2*d2(ind)*one1
              derv2(lxc,jzc) = derv2(lxc,jzc) + zd1*xd2*d2(ind)*one1
              derv2(lyc,jzc) = derv2(lyc,jzc) + zd1*yd2*d2(ind)*one1
              derv2(lzc,jzc) = derv2(lzc,jzc) + zd1*zd2*d2(ind)*one1
              if (n1.eq.n2) then
                derv2(lxc,jxc) = derv2(lxc,jxc) + d1(nf2)*one1
                derv2(lyc,jyc) = derv2(lyc,jyc) + d1(nf2)*one1
                derv2(lzc,jzc) = derv2(lzc,jzc) + d1(nf2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(jxc,lxc) = derv2(jxc,lxc) + xd1*xd2*d2(ind)*one2
                derv2(jxc,lyc) = derv2(jxc,lyc) + xd1*yd2*d2(ind)*one2
                derv2(jxc,lzc) = derv2(jxc,lzc) + xd1*zd2*d2(ind)*one2
                derv2(jyc,lxc) = derv2(jyc,lxc) + yd1*xd2*d2(ind)*one2
                derv2(jyc,lyc) = derv2(jyc,lyc) + yd1*yd2*d2(ind)*one2
                derv2(jyc,lzc) = derv2(jyc,lzc) + yd1*zd2*d2(ind)*one2
                derv2(jzc,lxc) = derv2(jzc,lxc) + zd1*xd2*d2(ind)*one2
                derv2(jzc,lyc) = derv2(jzc,lyc) + zd1*yd2*d2(ind)*one2
                derv2(jzc,lzc) = derv2(jzc,lzc) + zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(jxc,lxc) = derv2(jxc,lxc) + d1(nf2)*one2
                  derv2(jyc,lyc) = derv2(jyc,lyc) + d1(nf2)*one2
                  derv2(jzc,lzc) = derv2(jzc,lzc) + d1(nf2)*one2
                endif
              endif
            endif
            if (lstr) then
!
!  Strain - strain second derivatives
!
              do is1 = 1,nstrains
                ns1 = nstrptr(is1)
                do is2 = 1,nstrains
                  ns2 = nstrptr(is2)
                  sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds2(ns2)*dr2ds1(ns1)
                  if (n1.ne.n2) then
                    sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds1(ns2)*dr2ds2(ns1)
                  else
                    sderv2(is2,is1) = sderv2(is2,is1) + d1(nf2)*d2r2ds22(ns2,ns1)
                  endif
                enddo
              enddo
!
!  Internal - strain second derivatives
!
              do is1 = 1,nstrains
                ns1 = nstrptr(is1)
                derv3(ix,is1) = derv3(ix,is1) - xd1*d2(ind)*dr2ds2(ns1)
                derv3(iy,is1) = derv3(iy,is1) - yd1*d2(ind)*dr2ds2(ns1)
                derv3(iz,is1) = derv3(iz,is1) - zd1*d2(ind)*dr2ds2(ns1)
                derv3(jx,is1) = derv3(jx,is1) + xd1*d2(ind)*dr2ds2(ns1)
                derv3(jy,is1) = derv3(jy,is1) + yd1*d2(ind)*dr2ds2(ns1)
                derv3(jz,is1) = derv3(jz,is1) + zd1*d2(ind)*dr2ds2(ns1)
                if (n1.ne.n2) then
                  derv3(kx,is1) = derv3(kx,is1) - xd2*d2(ind)*dr2ds1(ns1)
                  derv3(ky,is1) = derv3(ky,is1) - yd2*d2(ind)*dr2ds1(ns1)
                  derv3(kz,is1) = derv3(kz,is1) - zd2*d2(ind)*dr2ds1(ns1)
                  derv3(lx,is1) = derv3(lx,is1) + xd2*d2(ind)*dr2ds1(ns1)
                  derv3(ly,is1) = derv3(ly,is1) + yd2*d2(ind)*dr2ds1(ns1)
                  derv3(lz,is1) = derv3(lz,is1) + zd2*d2(ind)*dr2ds1(ns1)
                else
                  derv3(ix,is1) = derv3(ix,is1) - d1(nf2)*d2r2dsdx2(ns1,1)
                  derv3(iy,is1) = derv3(iy,is1) - d1(nf2)*d2r2dsdx2(ns1,2)
                  derv3(iz,is1) = derv3(iz,is1) - d1(nf2)*d2r2dsdx2(ns1,3)
                  derv3(jx,is1) = derv3(jx,is1) + d1(nf2)*d2r2dsdx2(ns1,1)
                  derv3(jy,is1) = derv3(jy,is1) + d1(nf2)*d2r2dsdx2(ns1,2)
                  derv3(jz,is1) = derv3(jz,is1) + d1(nf2)*d2r2dsdx2(ns1,3)
                endif
              enddo
            endif
          endif
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addsp')
#endif
!
  return
  end
!
  subroutine d2addmixd(lilocal,lmlocal,lnlocal,i,ix,iy,iz,ixf,iyf,izf,maxneigh,mneigh,nneigh,neighno, &
                       xneigh,yneigh,zneigh,nREBOatomRptr,d2,mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf, &
                       xmn,ymn,zmn,dr2dsmn,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivatives.
!  This routines handles a combination of the derivatives of a single vector 
!  and derivatives stored in d1 format (but named d2 since these are 2nd
!  derivatives between the two distances).
!
!  On entry : 
!
!  lilocal         = if .true. then i is local to this processor
!  lmlocal         = if .true. then m is local to this processor
!  lnlocal         = if .true. then n is local to this processor
!  i               = atom number of i to whom others are connected
!  ix(f)           = position of x coordinate for i in second derivative matrix (f for full; otherwise local)
!  iy(f)           = position of y coordinate for i in second derivative matrix (f for full; otherwise local)
!  iz(f)           = position of z coordinate for i in second derivative matrix (f for full; otherwise local)
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  d2              = sparse array of second derivatives w.r.t. two interatomic distances
!  ltorderv        = if true then include torsional derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  mx(f)           = position of x coordinate for m in second derivative matrix (f for full; otherwise local)
!  my(f)           = position of y coordinate for m in second derivative matrix (f for full; otherwise local)
!  mz(f)           = position of z coordinate for m in second derivative matrix (f for full; otherwise local)
!  nx(f)           = position of x coordinate for n in second derivative matrix (f for full; otherwise local)
!  ny(f)           = position of y coordinate for n in second derivative matrix (f for full; otherwise local)
!  nz(f)           = position of z coordinate for n in second derivative matrix (f for full; otherwise local)
!  xmn, ymn, zmn   = Cartesian components of vector from m to n
!  dr2dsmn         = strain derivatives from m-n vector
!
!   5/23 Created from d2addmixdsp
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use current,        only : ndim
  use derivatives
  use m_strain,       only : real1strterm
  use parallel,       only : atom2local
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: ix
  integer(i4), intent(in)             :: iy
  integer(i4), intent(in)             :: iz
  integer(i4), intent(in)             :: ixf
  integer(i4), intent(in)             :: iyf
  integer(i4), intent(in)             :: izf
  integer(i4), intent(in)             :: mx
  integer(i4), intent(in)             :: my
  integer(i4), intent(in)             :: mz
  integer(i4), intent(in)             :: mxf
  integer(i4), intent(in)             :: myf
  integer(i4), intent(in)             :: mzf
  integer(i4), intent(in)             :: nx
  integer(i4), intent(in)             :: ny
  integer(i4), intent(in)             :: nz
  integer(i4), intent(in)             :: nxf
  integer(i4), intent(in)             :: nyf
  integer(i4), intent(in)             :: nzf
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  logical,     intent(in)             :: lilocal
  logical,     intent(in)             :: lmlocal
  logical,     intent(in)             :: lnlocal
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  real(dp),    intent(in)             :: xmn
  real(dp),    intent(in)             :: ymn
  real(dp),    intent(in)             :: zmn
  real(dp),    intent(in)             :: dr2dsmn(6)
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: indn2n
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: kloc
  integer(i4)                         :: kx
  integer(i4)                         :: ky
  integer(i4)                         :: kz
  integer(i4)                         :: kxf
  integer(i4)                         :: kyf
  integer(i4)                         :: kzf
  integer(i4)                         :: l
  integer(i4)                         :: lloc
  integer(i4)                         :: lx
  integer(i4)                         :: ly
  integer(i4)                         :: lz
  integer(i4)                         :: lxf
  integer(i4)                         :: lyf
  integer(i4)                         :: lzf
  integer(i4)                         :: n1j2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: ntmp
  logical                             :: lklocal
  logical                             :: lllocal
  real(dp)                            :: dr2ds(6)
  real(dp)                            :: d2r2dx2(3,3)
  real(dp)                            :: d2r2ds2(6,6)
  real(dp)                            :: d2r2dsdx(6,3)
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xkl
  real(dp)                            :: ykl
  real(dp)                            :: zkl
#ifdef TRACE
  call trace_in('d2addmixd')
#endif
!
!  Set up terms for i
!
  nri = nREBOatomRptr(i)
!
!  Set up indices for boundaries
!
  indn2n = nneigh(nri) + nneigh(nri)*(nneigh(nri)+1)/2
  if (ltorderv) then
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh
    ndist = ndisttor
  else
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
    ndist = indn2n
  endif
!
  do ind = 1,ndist
    if (abs(d2(ind)).gt.d2tol) then
      l = ind
      if (l.le.nneigh(nri)) then
!
!  i -> neighbours
!
        k = neighno(l,nri)
        xik = xneigh(l,nri)
        yik = yneigh(l,nri)
        zik = zneigh(l,nri)
!
        kloc = atom2local(k)
        lklocal = (kloc.gt.0)
!
        kx = 3*(kloc-1) + 1
        ky = kx + 1
        kz = kx + 2
!
        kxf = 3*(k-1) + 1
        kyf = kxf + 1
        kzf = kxf + 2
!
        call real1strterm(ndim,xik,yik,zik,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
        call add_drv2_2dm(lilocal,lklocal,lmlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xik,yik,zik,xmn,ymn,zmn, &
                          d2(ind),dr2ds,dr2dsmn,lilocal)
      elseif (l.le.indn2n) then
!
!  Loop over neighbours -> neighbours
!
        ntmp = l - nneigh(nri)
        nk = 1
        n1j2 = 1
        do while (n1j2.lt.ntmp)
          nk = nk + 1
          n1j2 = n1j2 + nk
        enddo
        n1j2 = n1j2 - nk
        nl = ntmp - n1j2
!
        k = neighno(nk,nri)
        l = neighno(nl,nri)
!
        kloc = atom2local(k)
        lloc = atom2local(l)
        lklocal = (kloc.gt.0)
        lllocal = (lloc.gt.0)
!
        xkl = (xneigh(nk,nri) - xneigh(nl,nri))
        ykl = (yneigh(nk,nri) - yneigh(nl,nri))
        zkl = (zneigh(nk,nri) - zneigh(nl,nri))
!
        kx = 3*(kloc-1) + 1
        ky = kx + 1
        kz = kx + 2
!
        kxf = 3*(k-1) + 1
        kyf = kxf + 1
        kzf = kxf + 2
!
        lx = 3*(lloc-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        lxf = 3*(l-1) + 1
        lyf = lxf + 1
        lzf = lxf + 2
!
        call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
        call add_drv2_2dm(lklocal,lllocal,lmlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xmn,ymn,zmn, &
                          d2(ind),dr2ds,dr2dsmn,lilocal)
      elseif (l.le.ndisttor) then
        ntmp = l - indn2n
        nj = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
        ntmp = ntmp - (nj - 1)*(mneigh + 1)*mneigh
        nk = (ntmp - 1)/mneigh
        nl = ntmp - nk*mneigh
!
        if (nk.eq.0) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
          j = neighno(nj,nri)
          nrj = nREBOatomRptr(j)
          l = neighno(nl,nrj)
!
          lloc = atom2local(l)
          lllocal = (lloc.gt.0)
!
          xil = xneigh(nl,nrj) + xneigh(nj,nri)
          yil = yneigh(nl,nrj) + yneigh(nj,nri)
          zil = zneigh(nl,nrj) + zneigh(nj,nri)
!
          lx = 3*(lloc-1) + 1
          ly = lx + 1
          lz = lx + 2
!
          lxf = 3*(l-1) + 1
          lyf = lxf + 1
          lzf = lxf + 2
!
          call real1strterm(ndim,xil,yil,zil,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
          call add_drv2_2dm(lilocal,lllocal,lmlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                            mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xil,yil,zil,xmn,ymn,zmn, &
                            d2(ind),dr2ds,dr2dsmn,lilocal)
        else
!
!  Loop over neighbours -> neighbours of neighbours
!
          j = neighno(nj,nri)
          nrj = nREBOatomRptr(j)
          k = neighno(nk,nri)
          l = neighno(nl,nrj)
!
          kloc = atom2local(k)
          lloc = atom2local(l)
          lklocal = (kloc.gt.0)
          lllocal = (lloc.gt.0)
!
          xkl = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
          ykl = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
          zkl = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
!
          kx = 3*(kloc-1) + 1
          ky = kx + 1
          kz = kx + 2
!
          kxf = 3*(k-1) + 1
          kyf = kxf + 1
          kzf = kxf + 2
!
          lx = 3*(lloc-1) + 1
          ly = lx + 1
          lz = lx + 2
!
          lxf = 3*(l-1) + 1
          lyf = lxf + 1
          lzf = lxf + 2
!
          call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
          call add_drv2_2dm(lklocal,lllocal,lmlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                            mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xmn,ymn,zmn, &
                            d2(ind),dr2ds,dr2dsmn,lilocal)
        endif
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addmixd')
#endif
!
  return
  end
!
  subroutine d2addmixspp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nREBOatomRptr,d2, &
                         nonzero,nonzeroptr,mx,my,mz,nx,ny,nz,xmn,ymn,zmn,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivatives. Phonon version.
!  This routines handles a combination of the derivatives of a single vector 
!  and derivatives stored in d1 format (but named d2 since these are 2nd
!  derivatives between the two distances).
!  Version that uses sparsity in d2
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  d2              = sparse array of second derivatives w.r.t. two interatomic distances
!  nonzero         = number of nonzero elements in d2
!  nonzeroptr(*)   = pointer to nonzero elements in d2
!  ltorderv        = if true then include torsional derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  mx, my, mz      = indices of the atom m in the 2nd derivative matrix
!  nx, ny, nz      = indices of the atom n in the 2nd derivative matrix
!  xmn, ymn, zmn   = Cartesian components of vector from m to n
!
!   5/23 Created from d2addmixsp
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
  use brennerdata
  use derivatives
  use m_strain,       only : real1strterm
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: mx
  integer(i4), intent(in)             :: my
  integer(i4), intent(in)             :: mz
  integer(i4), intent(in)             :: nx
  integer(i4), intent(in)             :: ny
  integer(i4), intent(in)             :: nz
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nonzero
  integer(i4), intent(in)             :: nonzeroptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  real(dp),    intent(in)             :: xmn
  real(dp),    intent(in)             :: ymn
  real(dp),    intent(in)             :: zmn
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: indn2n
  integer(i4)                         :: ix
  integer(i4)                         :: iy
  integer(i4)                         :: iz
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: kx
  integer(i4)                         :: ky
  integer(i4)                         :: kz
  integer(i4)                         :: l
  integer(i4)                         :: lx
  integer(i4)                         :: ly
  integer(i4)                         :: lz
  integer(i4)                         :: n1j2
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: ntmp
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xkl
  real(dp)                            :: ykl
  real(dp)                            :: zkl
#ifdef TRACE
  call trace_in('d2addmixspp')
#endif
!
!  Set up terms for i
!
  ix = 3*(i-1) + 1
  iy = ix + 1
  iz = ix + 2
  nri = nREBOatomRptr(i)
!
!  Set up indices for boundaries
!
  indn2n = nneigh(nri) + nneigh(nri)*(nneigh(nri)+1)/2
  if (ltorderv) then
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh
  else
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
  do ind = 1,nonzero
    l = nonzeroptr(ind)
    if (l.le.nneigh(nri)) then
!
!  i -> neighbours
!
      k = neighno(l,nri)
      xik = xneigh(l,nri)
      yik = yneigh(l,nri)
      zik = zneigh(l,nri)
!
      kx = 3*(k-1) + 1
      ky = kx + 1
      kz = kx + 2
!
      call add_drv2_2p(ix,iy,iz,kx,ky,kz,mx,my,mz,nx,ny,nz,xik,yik,zik,xmn,ymn,zmn,d2(ind))
    elseif (l.le.indn2n) then
!
!  Loop over neighbours -> neighbours
!
      ntmp = l - nneigh(nri)
      nk = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        nk = nk + 1
        n1j2 = n1j2 + nk
      enddo
      n1j2 = n1j2 - nk
      nl = ntmp - n1j2
!
      k = neighno(nk,nri)
      l = neighno(nl,nri)
      xkl = (xneigh(nk,nri) - xneigh(nl,nri))
      ykl = (yneigh(nk,nri) - yneigh(nl,nri))
      zkl = (zneigh(nk,nri) - zneigh(nl,nri))
!
      kx = 3*(k-1) + 1
      ky = kx + 1
      kz = kx + 2
      lx = 3*(l-1) + 1
      ly = lx + 1
      lz = lx + 2
!
      call add_drv2_2p(kx,ky,kz,lx,ly,lz,mx,my,mz,nx,ny,nz,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind))
    elseif (l.le.ndisttor) then
      ntmp = l - indn2n
      nj = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nj - 1)*(mneigh + 1)*mneigh
      nk = (ntmp - 1)/mneigh
      nl = ntmp - nk*mneigh
!
      if (nk.eq.0) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        l = neighno(nl,nrj)
        xil = xneigh(nl,nrj) + xneigh(nj,nri)
        yil = yneigh(nl,nrj) + yneigh(nj,nri)
        zil = zneigh(nl,nrj) + zneigh(nj,nri)
!
        lx = 3*(l-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        call add_drv2_2p(ix,iy,iz,lx,ly,lz,mx,my,mz,nx,ny,nz,xil,yil,zil,xmn,ymn,zmn,d2(ind))
      else
!
!  Loop over neighbours -> neighbours of neighbours
!
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        k = neighno(nk,nri)
        kx = 3*(k-1) + 1
        ky = kx + 1
        kz = kx + 2
        l = neighno(nl,nrj)
        xkl = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
        ykl = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
        zkl = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
        lx = 3*(l-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        call add_drv2_2p(kx,ky,kz,lx,ly,lz,mx,my,mz,nx,ny,nz,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind))
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addmixspp')
#endif
!
  return
  end
!
  subroutine d2addmixdsp(lilocal,lmlocal,lnlocal,i,ix,iy,iz,ixf,iyf,izf,maxneigh,mneigh,nneigh,neighno, &
                         xneigh,yneigh,zneigh,nREBOatomRptr,d2,nonzero,nonzeroptr,mx,my,mz,mxf,myf,mzf, &
                         nx,ny,nz,nxf,nyf,nzf,xmn,ymn,zmn,dr2dsmn,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivatives.
!  This routines handles a combination of the derivatives of a single vector 
!  and derivatives stored in d1 format (but named d2 since these are 2nd
!  derivatives between the two distances).
!  Version that uses sparsity in d2 with distributed parallel algorithm
!
!  On entry : 
!
!  lilocal         = if .true. then i is local to this processor
!  lmlocal         = if .true. then m is local to this processor
!  lnlocal         = if .true. then n is local to this processor
!  i               = atom number of i to whom others are connected
!  ix(f)           = position of x coordinate for i in second derivative matrix (f for full; otherwise local)
!  iy(f)           = position of y coordinate for i in second derivative matrix (f for full; otherwise local)
!  iz(f)           = position of z coordinate for i in second derivative matrix (f for full; otherwise local)
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  d2              = sparse array of second derivatives w.r.t. two interatomic distances
!  nonzero         = number of nonzero elements in d2
!  nonzeroptr(*)   = pointer to nonzero elements in d2
!  ltorderv        = if true then include torsional derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  mx(f)           = position of x coordinate for m in second derivative matrix (f for full; otherwise local)
!  my(f)           = position of y coordinate for m in second derivative matrix (f for full; otherwise local)
!  mz(f)           = position of z coordinate for m in second derivative matrix (f for full; otherwise local)
!  nx(f)           = position of x coordinate for n in second derivative matrix (f for full; otherwise local)
!  ny(f)           = position of y coordinate for n in second derivative matrix (f for full; otherwise local)
!  nz(f)           = position of z coordinate for n in second derivative matrix (f for full; otherwise local)
!  xmn, ymn, zmn   = Cartesian components of vector from m to n
!  dr2dsmn         = strain derivatives from m-n vector
!
!   5/23 Created from d2addmixsp
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
  use brennerdata
  use current,        only : ndim
  use derivatives
  use m_strain,       only : real1strterm
  use parallel,       only : atom2local
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: ix
  integer(i4), intent(in)             :: iy
  integer(i4), intent(in)             :: iz
  integer(i4), intent(in)             :: ixf
  integer(i4), intent(in)             :: iyf
  integer(i4), intent(in)             :: izf
  integer(i4), intent(in)             :: mx
  integer(i4), intent(in)             :: my
  integer(i4), intent(in)             :: mz
  integer(i4), intent(in)             :: mxf
  integer(i4), intent(in)             :: myf
  integer(i4), intent(in)             :: mzf
  integer(i4), intent(in)             :: nx
  integer(i4), intent(in)             :: ny
  integer(i4), intent(in)             :: nz
  integer(i4), intent(in)             :: nxf
  integer(i4), intent(in)             :: nyf
  integer(i4), intent(in)             :: nzf
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nonzero
  integer(i4), intent(in)             :: nonzeroptr(*)
  logical,     intent(in)             :: lilocal
  logical,     intent(in)             :: lmlocal
  logical,     intent(in)             :: lnlocal
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  real(dp),    intent(in)             :: xmn
  real(dp),    intent(in)             :: ymn
  real(dp),    intent(in)             :: zmn
  real(dp),    intent(in)             :: dr2dsmn(6)
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: indn2n
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: kloc
  integer(i4)                         :: kx
  integer(i4)                         :: ky
  integer(i4)                         :: kz
  integer(i4)                         :: kxf
  integer(i4)                         :: kyf
  integer(i4)                         :: kzf
  integer(i4)                         :: l
  integer(i4)                         :: lloc
  integer(i4)                         :: lx
  integer(i4)                         :: ly
  integer(i4)                         :: lz
  integer(i4)                         :: lxf
  integer(i4)                         :: lyf
  integer(i4)                         :: lzf
  integer(i4)                         :: n1j2
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: ntmp
  logical                             :: lklocal
  logical                             :: lllocal
  real(dp)                            :: dr2ds(6)
  real(dp)                            :: d2r2dx2(3,3)
  real(dp)                            :: d2r2ds2(6,6)
  real(dp)                            :: d2r2dsdx(6,3)
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xkl
  real(dp)                            :: ykl
  real(dp)                            :: zkl
#ifdef TRACE
  call trace_in('d2addmixdsp')
#endif
!
!  Set up terms for i
!
  nri = nREBOatomRptr(i)
!
!  Set up indices for boundaries
!
  indn2n = nneigh(nri) + nneigh(nri)*(nneigh(nri)+1)/2
  if (ltorderv) then
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh
  else
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
  do ind = 1,nonzero
    l = nonzeroptr(ind)
    if (l.le.nneigh(nri)) then
!
!  i -> neighbours
!
      k = neighno(l,nri)
      xik = xneigh(l,nri)
      yik = yneigh(l,nri)
      zik = zneigh(l,nri)
!
      kloc = atom2local(k)
      lklocal = (kloc.gt.0)
!
      kx = 3*(kloc-1) + 1
      ky = kx + 1
      kz = kx + 2
!
      kxf = 3*(k-1) + 1
      kyf = kxf + 1
      kzf = kxf + 2
!
      call real1strterm(ndim,xik,yik,zik,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
      call add_drv2_2dm(lilocal,lklocal,lmlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xik,yik,zik,xmn,ymn,zmn, &
                        d2(ind),dr2ds,dr2dsmn,lilocal)
    elseif (l.le.indn2n) then
!
!  Loop over neighbours -> neighbours
!
      ntmp = l - nneigh(nri)
      nk = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        nk = nk + 1
        n1j2 = n1j2 + nk
      enddo
      n1j2 = n1j2 - nk
      nl = ntmp - n1j2
!
      k = neighno(nk,nri)
      l = neighno(nl,nri)
!
      kloc = atom2local(k)
      lloc = atom2local(l)
      lklocal = (kloc.gt.0)
      lllocal = (lloc.gt.0)
!
      xkl = (xneigh(nk,nri) - xneigh(nl,nri))
      ykl = (yneigh(nk,nri) - yneigh(nl,nri))
      zkl = (zneigh(nk,nri) - zneigh(nl,nri))
!
      kx = 3*(kloc-1) + 1
      ky = kx + 1
      kz = kx + 2
!
      kxf = 3*(k-1) + 1
      kyf = kxf + 1
      kzf = kxf + 2
!
      lx = 3*(lloc-1) + 1
      ly = lx + 1
      lz = lx + 2
!
      lxf = 3*(l-1) + 1
      lyf = lxf + 1
      lzf = lxf + 2
!
      call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
      call add_drv2_2dm(lklocal,lllocal,lmlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                        mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xmn,ymn,zmn, &
                        d2(ind),dr2ds,dr2dsmn,lilocal)
    elseif (l.le.ndisttor) then
      ntmp = l - indn2n
      nj = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nj - 1)*(mneigh + 1)*mneigh
      nk = (ntmp - 1)/mneigh
      nl = ntmp - nk*mneigh
!
      if (nk.eq.0) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        l = neighno(nl,nrj)
!
        lloc = atom2local(l)
        lllocal = (lloc.gt.0)
!
        xil = xneigh(nl,nrj) + xneigh(nj,nri)
        yil = yneigh(nl,nrj) + yneigh(nj,nri)
        zil = zneigh(nl,nrj) + zneigh(nj,nri)
!
        lx = 3*(lloc-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        lxf = 3*(l-1) + 1
        lyf = lxf + 1
        lzf = lxf + 2
!
        call real1strterm(ndim,xil,yil,zil,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
        call add_drv2_2dm(lilocal,lllocal,lmlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                          mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xil,yil,zil,xmn,ymn,zmn, &
                          d2(ind),dr2ds,dr2dsmn,lilocal)
      else
!
!  Loop over neighbours -> neighbours of neighbours
!
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        k = neighno(nk,nri)
        l = neighno(nl,nrj)
!
        kloc = atom2local(k)
        lloc = atom2local(l)
        lklocal = (kloc.gt.0)
        lllocal = (lloc.gt.0)
!
        xkl = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
        ykl = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
        zkl = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
!
        kx = 3*(kloc-1) + 1
        ky = kx + 1
        kz = kx + 2
!
        kxf = 3*(k-1) + 1
        kyf = kxf + 1
        kzf = kxf + 2
!
        lx = 3*(lloc-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        lxf = 3*(l-1) + 1
        lyf = lxf + 1
        lzf = lxf + 2
!
        call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
        call add_drv2_2dm(lklocal,lllocal,lmlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          mx,my,mz,mxf,myf,mzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xmn,ymn,zmn, &
                          d2(ind),dr2ds,dr2dsmn,lilocal)
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addmixdsp')
#endif
!
  return
  end
!
  subroutine d2addmixsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d2, &
                        nonzero,nonzeroptr,mx,my,mz,nx,ny,nz,xmn,ymn,zmn,dr2dsmn,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivatives.
!  This routines handles a combination of the derivatives of a single vector 
!  and derivatives stored in d1 format (but named d2 since these are 2nd
!  derivatives between the two distances).
!  Version that uses sparsity in d2
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  d2              = sparse array of second derivatives w.r.t. two interatomic distances
!  nonzero         = number of nonzero elements in d2
!  nonzeroptr(*)   = pointer to nonzero elements in d2
!  ltorderv        = if true then include torsional derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  mx, my, mz      = indices of the atom m in the 2nd derivative matrix
!  nx, ny, nz      = indices of the atom n in the 2nd derivative matrix
!  xmn, ymn, zmn   = Cartesian components of vector from m to n
!  dr2dsmn         = strain derivatives from m-n vector
!
!   5/23 Created from d2addmix
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
  use brennerdata
  use current,        only : ndim
  use derivatives
  use m_strain,       only : real1strterm
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: mx
  integer(i4), intent(in)             :: my
  integer(i4), intent(in)             :: mz
  integer(i4), intent(in)             :: nx
  integer(i4), intent(in)             :: ny
  integer(i4), intent(in)             :: nz
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nonzero
  integer(i4), intent(in)             :: nonzeroptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  real(dp),    intent(in)             :: xmn
  real(dp),    intent(in)             :: ymn
  real(dp),    intent(in)             :: zmn
  real(dp),    intent(in)             :: dr2dsmn(6)
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: indn2n
  integer(i4)                         :: io
  integer(i4)                         :: ix
  integer(i4)                         :: iy
  integer(i4)                         :: iz
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: ko
  integer(i4)                         :: kx
  integer(i4)                         :: ky
  integer(i4)                         :: kz
  integer(i4)                         :: l
  integer(i4)                         :: lo
  integer(i4)                         :: lx
  integer(i4)                         :: ly
  integer(i4)                         :: lz
  integer(i4)                         :: n1j2
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: ntmp
  real(dp)                            :: dr2ds(6)
  real(dp)                            :: d2r2dx2(3,3)
  real(dp)                            :: d2r2ds2(6,6)
  real(dp)                            :: d2r2dsdx(6,3)
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xkl
  real(dp)                            :: ykl
  real(dp)                            :: zkl
#ifdef TRACE
  call trace_in('d2addmixsp')
#endif
!
!  Set up terms for i
!
  io = nfreeatom(i)
  ix = 3*(io-1) + 1
  iy = ix + 1
  iz = ix + 2
  nri = nREBOatomRptr(i)
!
!  Set up indices for boundaries
!
  indn2n = nneigh(nri) + nneigh(nri)*(nneigh(nri)+1)/2
  if (ltorderv) then
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh
  else
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = indn2n + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
  do ind = 1,nonzero
    l = nonzeroptr(ind)
    if (l.le.nneigh(nri)) then
!
!  i -> neighbours
!
      k = neighno(l,nri)
      ko = nfreeatom(k)
      if (io.gt.0.and.ko.gt.0) then
        xik = xneigh(l,nri)
        yik = yneigh(l,nri)
        zik = zneigh(l,nri)
!
        kx = 3*(ko-1) + 1
        ky = kx + 1
        kz = kx + 2
!
        call real1strterm(ndim,xik,yik,zik,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
        call add_drv2_2(ix,iy,iz,kx,ky,kz,mx,my,mz,nx,ny,nz,xik,yik,zik,xmn,ymn,zmn,d2(ind),dr2ds,dr2dsmn)
      endif
    elseif (l.le.indn2n) then
!
!  Loop over neighbours -> neighbours
!
      ntmp = l - nneigh(nri)
      nk = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        nk = nk + 1
        n1j2 = n1j2 + nk
      enddo
      n1j2 = n1j2 - nk
      nl = ntmp - n1j2
!
      k = neighno(nk,nri)
      ko = nfreeatom(k)
      l = neighno(nl,nri)
      lo = nfreeatom(l)
      if (ko.gt.0.and.lo.gt.0) then
        xkl = (xneigh(nk,nri) - xneigh(nl,nri))
        ykl = (yneigh(nk,nri) - yneigh(nl,nri))
        zkl = (zneigh(nk,nri) - zneigh(nl,nri))
!
        kx = 3*(ko-1) + 1
        ky = kx + 1
        kz = kx + 2
        lx = 3*(lo-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
        call add_drv2_2(kx,ky,kz,lx,ly,lz,mx,my,mz,nx,ny,nz,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind),dr2ds,dr2dsmn)
      endif
    elseif (l.le.ndisttor) then
      ntmp = l - indn2n
      nj = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nj - 1)*(mneigh + 1)*mneigh
      nk = (ntmp - 1)/mneigh
      nl = ntmp - nk*mneigh
!
      if (nk.eq.0) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        l = neighno(nl,nrj)
        lo = nfreeatom(l)
        if (io.gt.0.and.lo.gt.0) then
          xil = xneigh(nl,nrj) + xneigh(nj,nri)
          yil = yneigh(nl,nrj) + yneigh(nj,nri)
          zil = zneigh(nl,nrj) + zneigh(nj,nri)
!
          lx = 3*(lo-1) + 1
          ly = lx + 1
          lz = lx + 2
!
          call real1strterm(ndim,xil,yil,zil,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
          call add_drv2_2(ix,iy,iz,lx,ly,lz,mx,my,mz,nx,ny,nz,xil,yil,zil,xmn,ymn,zmn,d2(ind),dr2ds,dr2dsmn)
        endif
      else
!
!  Loop over neighbours -> neighbours of neighbours
!
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        k = neighno(nk,nri)
        ko = nfreeatom(k)
        kx = 3*(ko-1) + 1
        ky = kx + 1
        kz = kx + 2
        l = neighno(nl,nrj)
        lo = nfreeatom(l)
        if (ko.gt.0.and.lo.gt.0) then
          xkl = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
          ykl = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
          zkl = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
          lx = 3*(lo-1) + 1
          ly = lx + 1
          lz = lx + 2
!
          call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
          call add_drv2_2(kx,ky,kz,lx,ly,lz,mx,my,mz,nx,ny,nz,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind),dr2ds,dr2dsmn)
        endif
      endif
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addmixsp')
#endif
!
  return
  end
!
  subroutine d2addmix(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d2, &
                      mx,my,mz,nx,ny,nz,xmn,ymn,zmn,dr2dsmn,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivatives.
!  This routines handles a combination of the derivatives of a single vector 
!  and derivatives stored in d1 format (but named d2 since these are 2nd
!  derivatives between the two distances).
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  d2              = array of second derivatives w.r.t. two interatomic distances
!  ltorderv        = if true then include torsional derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  mx, my, mz      = indices of the atom m in the 2nd derivative matrix
!  nx, ny, nz      = indices of the atom n in the 2nd derivative matrix
!  xmn, ymn, zmn   = Cartesian components of vector from m to n
!  dr2dsmn         = strain derivatives from m-n vector
!
!  12/22 Created from d1add
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use current,        only : ndim
  use derivatives
  use m_strain,       only : real1strterm
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: mx
  integer(i4), intent(in)             :: my
  integer(i4), intent(in)             :: mz
  integer(i4), intent(in)             :: nx
  integer(i4), intent(in)             :: ny
  integer(i4), intent(in)             :: nz
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  real(dp),    intent(in)             :: xmn
  real(dp),    intent(in)             :: ymn
  real(dp),    intent(in)             :: zmn
  real(dp),    intent(in)             :: dr2dsmn(6)
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: io
  integer(i4)                         :: ix
  integer(i4)                         :: iy
  integer(i4)                         :: iz
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: ko
  integer(i4)                         :: kx
  integer(i4)                         :: ky
  integer(i4)                         :: kz
  integer(i4)                         :: l
  integer(i4)                         :: lo
  integer(i4)                         :: lx
  integer(i4)                         :: ly
  integer(i4)                         :: lz
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  real(dp)                            :: dr2ds(6)
  real(dp)                            :: d2r2dx2(3,3)
  real(dp)                            :: d2r2ds2(6,6)
  real(dp)                            :: d2r2dsdx(6,3)
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xkl
  real(dp)                            :: ykl
  real(dp)                            :: zkl
#ifdef TRACE
  call trace_in('d2addmix')
#endif
!
!  Set up terms for i
!
  io = nfreeatom(i)
  ix = 3*(io-1) + 1
  iy = ix + 1
  iz = ix + 2
!
!  Loop over i -> neighbours
!
  ind = 0
  nri = nREBOatomRptr(i)
  do nk = 1,nneigh(nri)
    ind = ind + 1
    if (abs(d2(ind)).gt.d2tol) then
      k = neighno(nk,nri)
      ko = nfreeatom(k)
      if (io.gt.0.and.ko.gt.0) then
        xik = xneigh(nk,nri)
        yik = yneigh(nk,nri)
        zik = zneigh(nk,nri)
!
        kx = 3*(ko-1) + 1
        ky = kx + 1
        kz = kx + 2
!
        call real1strterm(ndim,xik,yik,zik,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
        call add_drv2_2(ix,iy,iz,kx,ky,kz,mx,my,mz,nx,ny,nz,xik,yik,zik,xmn,ymn,zmn,d2(ind),dr2ds,dr2dsmn)
      endif
    endif
  enddo
!
!  Loop over neighbours -> neighbours
!
  do nk = 2,nneigh(nri)
    k = neighno(nk,nri)
    ko = nfreeatom(k)
    do nl = 1,nk - 1
      ind = nneigh(nri) + nk*(nk-1)/2 + nl
      if (abs(d2(ind)).gt.d2tol) then
        l = neighno(nl,nri)
        lo = nfreeatom(l)
        if (ko.gt.0.and.lo.gt.0) then
          xkl = (xneigh(nk,nri) - xneigh(nl,nri))
          ykl = (yneigh(nk,nri) - yneigh(nl,nri))
          zkl = (zneigh(nk,nri) - zneigh(nl,nri))
!
          kx = 3*(ko-1) + 1
          ky = kx + 1
          kz = kx + 2
          lx = 3*(lo-1) + 1
          ly = lx + 1
          lz = lx + 2
!
          call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
          call add_drv2_2(kx,ky,kz,lx,ly,lz,mx,my,mz,nx,ny,nz,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind),dr2ds,dr2dsmn)
        endif
      endif
    enddo
  enddo
!
!  Loop over i -> neighbours of neighbours
!
  if (ltorderv) then
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh
      do nl = 1,nneigh(nrj)
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol) then
          l = neighno(nl,nrj)
          lo = nfreeatom(l)
          if (io.gt.0.and.lo.gt.0) then
            xil = xneigh(nl,nrj) + xneigh(nj,nri)
            yil = yneigh(nl,nrj) + yneigh(nj,nri)
            zil = zneigh(nl,nrj) + zneigh(nj,nri)
!
            lx = 3*(lo-1) + 1
            ly = lx + 1
            lz = lx + 2
!
            call real1strterm(ndim,xil,yil,zil,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
            call add_drv2_2(ix,iy,iz,lx,ly,lz,mx,my,mz,nx,ny,nz,xil,yil,zil,xmn,ymn,zmn,d2(ind),dr2ds,dr2dsmn)
          endif
        endif
      enddo
    enddo
!
!  Loop over neighbours -> neighbours of neighbours
!
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      do nk = 1,nneigh(nri)
        k = neighno(nk,nri)
        ko = nfreeatom(k)
        kx = 3*(ko-1) + 1
        ky = kx + 1
        kz = kx + 2
        ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh + nk*mneigh
        do nl = 1,nneigh(nrj)
          ind = ind + 1
          if (abs(d2(ind)).gt.d2tol) then
            l = neighno(nl,nrj)
            lo = nfreeatom(l)
            if (ko.gt.0.and.lo.gt.0) then
              xkl = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
              ykl = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
              zkl = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
              lx = 3*(lo-1) + 1
              ly = lx + 1
              lz = lx + 2
!
              call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
              call add_drv2_2(kx,ky,kz,lx,ly,lz,mx,my,mz,nx,ny,nz,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind),dr2ds,dr2dsmn)
            endif
          endif
        enddo
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('d2addmix')
#endif
!
  return
  end
!
  subroutine d2addmixp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nREBOatomRptr,d2, &
                       mx,my,mz,nx,ny,nz,xmn,ymn,zmn,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivatives. Phonon version.
!  This routines handles a combination of the derivatives of a single vector 
!  and derivatives stored in d1 format (but named d2 since these are 2nd
!  derivatives between the two distances).
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  d2              = array of second derivatives w.r.t. two interatomic distances
!  ltorderv        = if true then include torsional derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  mx, my, mz      = indices of the atom m in the 2nd derivative matrix
!  nx, ny, nz      = indices of the atom n in the 2nd derivative matrix
!  xmn, ymn, zmn   = Cartesian components of vector from m to n
!
!   4/23 Created from d2addmixp
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: mx
  integer(i4), intent(in)             :: my
  integer(i4), intent(in)             :: mz
  integer(i4), intent(in)             :: nx
  integer(i4), intent(in)             :: ny
  integer(i4), intent(in)             :: nz
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  real(dp),    intent(in)             :: xmn
  real(dp),    intent(in)             :: ymn
  real(dp),    intent(in)             :: zmn
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix
  integer(i4)                         :: iy
  integer(i4)                         :: iz
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: kx
  integer(i4)                         :: ky
  integer(i4)                         :: kz
  integer(i4)                         :: l
  integer(i4)                         :: lx
  integer(i4)                         :: ly
  integer(i4)                         :: lz
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xkl
  real(dp)                            :: ykl
  real(dp)                            :: zkl
#ifdef TRACE
  call trace_in('d2addmixp')
#endif
!
!  Set up terms for i
!
  ix = 3*(i-1) + 1
  iy = ix + 1
  iz = ix + 2
!
!  Loop over i -> neighbours
!
  ind = 0
  nri = nREBOatomRptr(i)
  do nk = 1,nneigh(nri)
    ind = ind + 1
    if (abs(d2(ind)).gt.d2tol) then
      k = neighno(nk,nri)
      xik = xneigh(nk,nri)
      yik = yneigh(nk,nri)
      zik = zneigh(nk,nri)
!
      kx = 3*(k-1) + 1
      ky = kx + 1
      kz = kx + 2
!
      call add_drv2_2p(ix,iy,iz,kx,ky,kz,mx,my,mz,nx,ny,nz,xik,yik,zik,xmn,ymn,zmn,d2(ind))
    endif
  enddo
!
!  Loop over neighbours -> neighbours
!
  do nk = 2,nneigh(nri)
    k = neighno(nk,nri)
    do nl = 1,nk - 1
      ind = nneigh(nri) + nk*(nk-1)/2 + nl
      if (abs(d2(ind)).gt.d2tol) then
        l = neighno(nl,nri)
        xkl = (xneigh(nk,nri) - xneigh(nl,nri))
        ykl = (yneigh(nk,nri) - yneigh(nl,nri))
        zkl = (zneigh(nk,nri) - zneigh(nl,nri))
!
        kx = 3*(k-1) + 1
        ky = kx + 1
        kz = kx + 2
        lx = 3*(l-1) + 1
        ly = lx + 1
        lz = lx + 2
!
        call add_drv2_2p(kx,ky,kz,lx,ly,lz,mx,my,mz,nx,ny,nz,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind))
      endif
    enddo
  enddo
!
!  Loop over i -> neighbours of neighbours
!
  if (ltorderv) then
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh
      do nl = 1,nneigh(nrj)
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol) then
          l = neighno(nl,nrj)
          xil = xneigh(nl,nrj) + xneigh(nj,nri)
          yil = yneigh(nl,nrj) + yneigh(nj,nri)
          zil = zneigh(nl,nrj) + zneigh(nj,nri)
!
          lx = 3*(l-1) + 1
          ly = lx + 1
          lz = lx + 2
!
          call add_drv2_2p(ix,iy,iz,lx,ly,lz,mx,my,mz,nx,ny,nz,xil,yil,zil,xmn,ymn,zmn,d2(ind))
        endif
      enddo
    enddo
!
!  Loop over neighbours -> neighbours of neighbours
!
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      do nk = 1,nneigh(nri)
        k = neighno(nk,nri)
        kx = 3*(k-1) + 1
        ky = kx + 1
        kz = kx + 2
        ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh + nk*mneigh
        do nl = 1,nneigh(nrj)
          ind = ind + 1
          if (abs(d2(ind)).gt.d2tol) then
            l = neighno(nl,nrj)
            xkl = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
            ykl = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
            zkl = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
            lx = 3*(l-1) + 1
            ly = lx + 1
            lz = lx + 2
            call add_drv2_2p(kx,ky,kz,lx,ly,lz,mx,my,mz,nx,ny,nz,xkl,ykl,zkl,xmn,ymn,zmn,d2(ind))
          endif
        enddo
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('d2addmixp')
#endif
!
  return
  end
!
  subroutine d1add(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1,ltorderv,lextended,ldoregions)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian first derivative arrays.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!  ldoregions      = if true then region derivatives should be computed
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!
!   5/02 Created
!   6/02 Strains added
!   6/02 Virial added
!   6/02 Modified to allow for torsional terms
!   7/02 Indexing further modified for torsions
!   8/02 Freezing added
!  12/03 Torsional derivatives made optional
!   6/07 nREBOatomRptr added as input
!   4/07 mneigh added as a distinct number from maxneigh
!  11/09 Region derivatives added
!   4/12 Explicit virial calculation removed as no longer needed
!   5/12 Atomic stresses added
!   4/16 ldoregions flag added
!   2/18 Trace added
!   9/18 Strain module added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!  12/22 Extended derivatives added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use control,        only : latomicstress
  use current,        only : ndim, nrelf2a, nregionno
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1(*)
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
  logical,     intent(in)             :: ldoregions
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: io
  integer(i4)                         :: is
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: ko
  integer(i4)                         :: l
  integer(i4)                         :: lo
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nregi
  integer(i4)                         :: nregk
  integer(i4)                         :: nregl
  integer(i4)                         :: ns
  real(dp)                            :: dr2ds(6)
  real(dp)                            :: d2r2dx2(3,3)
  real(dp)                            :: d2r2ds2(6,6)
  real(dp)                            :: d2r2dsdx(6,3)
  real(dp)                            :: xd
  real(dp)                            :: yd
  real(dp)                            :: zd
  real(dp)                            :: xdiff
  real(dp)                            :: ydiff
  real(dp)                            :: zdiff
#ifdef TRACE
  call trace_in('d1add')
#endif
!
!  Loop over i -> neighbours
!
  ind = 0
  nri = nREBOatomRptr(i)
  do nk = 1,nneigh(nri)
    ind = ind + 1
    k = neighno(nk,nri)
    xd = xneigh(nk,nri)*d1(ind)
    yd = yneigh(nk,nri)*d1(ind)
    zd = zneigh(nk,nri)*d1(ind)
    io = nfreeatom(i)
    ko = nfreeatom(k)
    if (io.gt.0) then
      xdrv(i) = xdrv(i) - xd
      ydrv(i) = ydrv(i) - yd
      zdrv(i) = zdrv(i) - zd
    endif
    if (ko.gt.0) then
      xdrv(k) = xdrv(k) + xd
      ydrv(k) = ydrv(k) + yd
      zdrv(k) = zdrv(k) + zd
    endif
    if (ldoregions) then
      nregi = nregionno(nrelf2a(i))
      nregk = nregionno(nrelf2a(k))
      if (nregi.ne.nregk) then
        xregdrv(nregi) = xregdrv(nregi) - xd
        yregdrv(nregi) = yregdrv(nregi) - yd
        zregdrv(nregi) = zregdrv(nregi) - zd
        xregdrv(nregk) = xregdrv(nregk) + xd
        yregdrv(nregk) = yregdrv(nregk) + yd
        zregdrv(nregk) = zregdrv(nregk) + zd
      endif
    endif
    if (lstr) then
      call real1strterm(ndim,xneigh(nk,nri),yneigh(nk,nri),zneigh(nk,nri),0.0_dp,0.0_dp,0.0_dp, &
        dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
      do is = 1,nstrains
        ns = nstrptr(is)
        rstrd(is) = rstrd(is) + d1(ind)*dr2ds(ns)
      enddo
      if (latomicstress) then
        do is = 1,nstrains
          ns = nstrptr(is)
          atomicstress(is,i) = atomicstress(is,i) + 0.5_dp*d1(ind)*dr2ds(ns)
          atomicstress(is,k) = atomicstress(is,k) + 0.5_dp*d1(ind)*dr2ds(ns)
        enddo
      endif
    endif
  enddo
!
!  Loop over neighbours -> neighbours
!
  do nk = 2,nneigh(nri)
    k = neighno(nk,nri)
    do nl = 1,nk - 1
      ind = nneigh(nri) + nk*(nk-1)/2 + nl
      if (abs(d1(ind)).gt.d1tol) then
        l = neighno(nl,nri)
        xdiff = (xneigh(nk,nri) - xneigh(nl,nri))
        ydiff = (yneigh(nk,nri) - yneigh(nl,nri))
        zdiff = (zneigh(nk,nri) - zneigh(nl,nri))
        xd = xdiff*d1(ind)
        yd = ydiff*d1(ind)
        zd = zdiff*d1(ind)
        lo = nfreeatom(l)
        ko = nfreeatom(k)
        if (lo.gt.0) then
          xdrv(l) = xdrv(l) - xd
          ydrv(l) = ydrv(l) - yd
          zdrv(l) = zdrv(l) - zd
        endif
        if (ko.gt.0) then
          xdrv(k) = xdrv(k) + xd
          ydrv(k) = ydrv(k) + yd
          zdrv(k) = zdrv(k) + zd
        endif
!
        if (ldoregions) then
          nregl = nregionno(nrelf2a(l))
          nregk = nregionno(nrelf2a(k))
          if (nregl.ne.nregk) then
            xregdrv(nregl) = xregdrv(nregl) - xd
            yregdrv(nregl) = yregdrv(nregl) - yd
            zregdrv(nregl) = zregdrv(nregl) - zd
            xregdrv(nregk) = xregdrv(nregk) + xd
            yregdrv(nregk) = yregdrv(nregk) + yd
            zregdrv(nregk) = zregdrv(nregk) + zd
          endif
        endif
        if (lstr) then
          call real1strterm(ndim,xdiff,ydiff,zdiff,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do is = 1,nstrains
            ns = nstrptr(is)
            rstrd(is) = rstrd(is) + d1(ind)*dr2ds(ns)
          enddo
          if (latomicstress) then
            do is = 1,nstrains
              ns = nstrptr(is)
              atomicstress(is,l) = atomicstress(is,l) + 0.5_dp*d1(ind)*dr2ds(ns)
              atomicstress(is,k) = atomicstress(is,k) + 0.5_dp*d1(ind)*dr2ds(ns)
            enddo
          endif
        endif
      endif
    enddo
  enddo
!
!  Loop over i -> neighbours of neighbours
!
  if (ltorderv) then
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh
      do nl = 1,nneigh(nrj)
        ind = ind + 1
        if (abs(d1(ind)).gt.d1tol) then
          l = neighno(nl,nrj)
          xdiff = xneigh(nl,nrj) + xneigh(nj,nri)
          ydiff = yneigh(nl,nrj) + yneigh(nj,nri)
          zdiff = zneigh(nl,nrj) + zneigh(nj,nri)
          xd = xdiff*d1(ind)
          yd = ydiff*d1(ind)
          zd = zdiff*d1(ind)
          io = nfreeatom(i)
          lo = nfreeatom(l)
          if (io.gt.0) then
            xdrv(i) = xdrv(i) - xd
            ydrv(i) = ydrv(i) - yd
            zdrv(i) = zdrv(i) - zd
          endif
          if (lo.gt.0) then
            xdrv(l) = xdrv(l) + xd
            ydrv(l) = ydrv(l) + yd
            zdrv(l) = zdrv(l) + zd
          endif
!
          if (ldoregions) then
            nregi = nregionno(nrelf2a(i))
            nregl = nregionno(nrelf2a(l))
            if (nregi.ne.nregl) then
              xregdrv(nregi) = xregdrv(nregi) - xd
              yregdrv(nregi) = yregdrv(nregi) - yd
              zregdrv(nregi) = zregdrv(nregi) - zd
              xregdrv(nregl) = xregdrv(nregl) + xd
              yregdrv(nregl) = yregdrv(nregl) + yd
              zregdrv(nregl) = zregdrv(nregl) + zd
            endif
          endif
!
          if (lstr) then
            call real1strterm(ndim,xdiff,ydiff,zdiff,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
            do is = 1,nstrains
              ns = nstrptr(is)
              rstrd(is) = rstrd(is) + d1(ind)*dr2ds(ns)
            enddo
            if (latomicstress) then
              do is = 1,nstrains
                ns = nstrptr(is)
                atomicstress(is,i) = atomicstress(is,i) + 0.5_dp*d1(ind)*dr2ds(ns)
                atomicstress(is,l) = atomicstress(is,l) + 0.5_dp*d1(ind)*dr2ds(ns)
              enddo
            endif
          endif
        endif
      enddo
    enddo
!
!  Loop over neighbours -> neighbours of neighbours
!
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      do nk = 1,nneigh(nri)
        k = neighno(nk,nri)
        ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh + nk*mneigh
        do nl = 1,nneigh(nrj)
          ind = ind + 1
          if (abs(d1(ind)).gt.d1tol) then
            l = neighno(nl,nrj)
            xdiff = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
            ydiff = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
            zdiff = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
            xd = xdiff*d1(ind)
            yd = ydiff*d1(ind)
            zd = zdiff*d1(ind)
            ko = nfreeatom(k)
            lo = nfreeatom(l)
            if (ko.gt.0) then
              xdrv(k) = xdrv(k) - xd
              ydrv(k) = ydrv(k) - yd
              zdrv(k) = zdrv(k) - zd
            endif
            if (lo.gt.0) then
              xdrv(l) = xdrv(l) + xd
              ydrv(l) = ydrv(l) + yd
              zdrv(l) = zdrv(l) + zd
            endif
!
            if (ldoregions) then
              nregk = nregionno(nrelf2a(k))
              nregl = nregionno(nrelf2a(l))
              if (nregk.ne.nregl) then
                xregdrv(nregk) = xregdrv(nregk) - xd
                yregdrv(nregk) = yregdrv(nregk) - yd
                zregdrv(nregk) = zregdrv(nregk) - zd
                xregdrv(nregl) = xregdrv(nregl) + xd
                yregdrv(nregl) = yregdrv(nregl) + yd
                zregdrv(nregl) = zregdrv(nregl) + zd
              endif
            endif
!
            if (lstr) then
              call real1strterm(ndim,xdiff,ydiff,zdiff,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
              do is = 1,nstrains
                ns = nstrptr(is)
                rstrd(is) = rstrd(is) + d1(ind)*dr2ds(ns)
              enddo
              if (latomicstress) then
                do is = 1,nstrains
                  ns = nstrptr(is)
                  atomicstress(is,k) = atomicstress(is,k) + 0.5_dp*d1(ind)*dr2ds(ns)
                  atomicstress(is,l) = atomicstress(is,l) + 0.5_dp*d1(ind)*dr2ds(ns)
                enddo
              endif
            endif
          endif
        enddo
      enddo
    enddo
!
    if (lextended) then
!
!  Extended derivatives over third shell
!
      do nj = 1,nneigh(nri)
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        do nk = 1,nneigh(nrj)
          k = neighno(nk,nrj)
          nrk = nREBOatomRptr(k)
          ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + nneigh(nri)*(mneigh+1)*mneigh + (nj-1)*mneigh*mneigh + (nk-1)*mneigh
          do nl = 1,nneigh(nrk)
            ind = ind + 1
            if (abs(d1(ind)).gt.d1tol) then
              l = neighno(nl,nrk)
              xdiff = xneigh(nl,nrk)
              ydiff = yneigh(nl,nrk)
              zdiff = zneigh(nl,nrk)
              xd = xdiff*d1(ind)
              yd = ydiff*d1(ind)
              zd = zdiff*d1(ind)
              ko = nfreeatom(k)
              lo = nfreeatom(l)
              if (ko.gt.0) then
                xdrv(k) = xdrv(k) - xd
                ydrv(k) = ydrv(k) - yd
                zdrv(k) = zdrv(k) - zd
              endif
              if (lo.gt.0) then
                xdrv(l) = xdrv(l) + xd
                ydrv(l) = ydrv(l) + yd
                zdrv(l) = zdrv(l) + zd
              endif
!
              if (ldoregions) then
                nregk = nregionno(nrelf2a(k))
                nregl = nregionno(nrelf2a(l))
                if (nregk.ne.nregl) then
                  xregdrv(nregk) = xregdrv(nregk) - xd
                  yregdrv(nregk) = yregdrv(nregk) - yd
                  zregdrv(nregk) = zregdrv(nregk) - zd
                  xregdrv(nregl) = xregdrv(nregl) + xd
                  yregdrv(nregl) = yregdrv(nregl) + yd
                  zregdrv(nregl) = zregdrv(nregl) + zd
                endif
              endif
!
              if (lstr) then
                call real1strterm(ndim,xdiff,ydiff,zdiff,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
                do is = 1,nstrains
                  ns = nstrptr(is)
                  rstrd(is) = rstrd(is) + d1(ind)*dr2ds(ns)
                enddo
                if (latomicstress) then
                  do is = 1,nstrains
                    ns = nstrptr(is)
                    atomicstress(is,k) = atomicstress(is,k) + 0.5_dp*d1(ind)*dr2ds(ns)
                    atomicstress(is,l) = atomicstress(is,l) + 0.5_dp*d1(ind)*dr2ds(ns)
                  enddo
                endif
              endif
            endif
          enddo
        enddo
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('d1add')
#endif
!
  return
  end
!
  subroutine d1addfc(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,maxlhs,d1cell,matom,nREBOatomRptr, &
                     d1,ltorderv,lextended)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian first derivative arrays.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  ineigh          = cell indices for vector to neighbour of atoms
!  matom           = atom being finite differenced
!  maxlhs          = second dimension of d1cell
!  d1cell          = array for storing forces separated by cell index
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!
!   1/15 Created from d1add and d2addfc
!   1/15 imi1/2/3 values corrected
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   5/23 Derivative tolerances as global variables introduced
!   7/23 Extended derivatives added
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
  use brennerdata
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: ineigh(3_i4,maxneigh,*)
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: maxlhs
  integer(i4), intent(in)             :: matom
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(out)            :: d1cell(4,maxlhs,*)
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: imi1
  integer(i4)                         :: imi2
  integer(i4)                         :: imi3
  integer(i4)                         :: imk1
  integer(i4)                         :: imk2
  integer(i4)                         :: imk3
  integer(i4)                         :: iml1
  integer(i4)                         :: iml2
  integer(i4)                         :: iml3
  integer(i4)                         :: inm1
  integer(i4)                         :: inm2
  integer(i4)                         :: inm3
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: l
  integer(i4)                         :: ncindm
  integer(i4)                         :: ncindp
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  logical                             :: lfound
  real(dp)                            :: xd
  real(dp)                            :: yd
  real(dp)                            :: zd
  real(dp)                            :: xdiff
  real(dp)                            :: ydiff
  real(dp)                            :: zdiff
#ifdef TRACE
  call trace_in('d1addfc')
#endif
!
!  Map atom i to relevant atom in arrays
!
  nri = nREBOatomRptr(i)
!
!  Find cell indices for matom - first search i and neighbours of i
!
  inm1 = 0
  inm2 = 0
  inm3 = 0
  if (i.ne.matom) then
    lfound = .false.
    nj = 0
    do while (.not.lfound.and.nj.lt.nneigh(nri))
      nj = nj + 1
      lfound = (neighno(nj,nri).eq.matom)
    enddo
    if (lfound) then
      inm1 = ineigh(1,nj,nri) 
      inm2 = ineigh(2,nj,nri) 
      inm3 = ineigh(3,nj,nri) 
    else
!
!  If matom hasn't been found then search neighbours of neighbours of i
!
      nj = 0
      do while (.not.lfound.and.nj.lt.nneigh(nri))
        nj = nj + 1
        nrj = nREBOatomRptr(neighno(nj,nri))
        nk = 0
        do while (.not.lfound.and.nk.lt.nneigh(nrj))
          nk = nk + 1
          lfound = (neighno(nk,nrj).eq.matom)
        enddo
      enddo
      if (lfound) then
        inm1 = ineigh(1,nk,nrj) + ineigh(1,nj,nri)
        inm2 = ineigh(2,nk,nrj) + ineigh(2,nj,nri)
        inm3 = ineigh(3,nk,nrj) + ineigh(3,nj,nri)
      endif
    endif
  endif
!
!  Loop over i -> neighbours
!
  ind = 0
  do nk = 1,nneigh(nri)
    ind = ind + 1
    k = neighno(nk,nri)
    xd = xneigh(nk,nri)*d1(ind)
    yd = yneigh(nk,nri)*d1(ind)
    zd = zneigh(nk,nri)*d1(ind)
!
!  Compute cell numbers from matom to i and k
!
    imi1 = - inm1
    imi2 = - inm2
    imi3 = - inm3
    imk1 = ineigh(1,nk,nri) - inm1
    imk2 = ineigh(2,nk,nri) - inm2
    imk3 = ineigh(3,nk,nri) - inm3
!
!  Compute cell indices for i
!
    if (abs(imi1).gt.nd2cell(1).or. &
        abs(imi2).gt.nd2cell(2).or. &
        abs(imi3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
      ncindm = nd2central
    else
!
!  Compute index
!
      ncindm = nd2cellptr(nd2cell(1)+1+imi1, &
                          nd2cell(2)+1+imi2, &
                          nd2cell(3)+1+imi3)
    endif
!
!  Compute cell indices for k
!
    if (abs(imk1).gt.nd2cell(1).or. &
        abs(imk2).gt.nd2cell(2).or. &
        abs(imk3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
      ncindp = nd2central
    else
!
!  Compute index
!
      ncindp = nd2cellptr(nd2cell(1)+1+imk1, &
                          nd2cell(2)+1+imk2, &
                          nd2cell(3)+1+imk3)
    endif
!
    d1cell(1,i,ncindm) = d1cell(1,i,ncindm) - xd
    d1cell(2,i,ncindm) = d1cell(2,i,ncindm) - yd
    d1cell(3,i,ncindm) = d1cell(3,i,ncindm) - zd
!
    d1cell(1,k,ncindp) = d1cell(1,k,ncindp) + xd
    d1cell(2,k,ncindp) = d1cell(2,k,ncindp) + yd
    d1cell(3,k,ncindp) = d1cell(3,k,ncindp) + zd
  enddo
!
!  Loop over neighbours -> neighbours
!
  do nk = 2,nneigh(nri)
    k = neighno(nk,nri)
    do nl = 1,nk - 1
      ind = nneigh(nri) + nk*(nk-1)/2 + nl
      if (abs(d1(ind)).gt.d1tol) then
        l = neighno(nl,nri)
        xdiff = (xneigh(nk,nri) - xneigh(nl,nri))
        ydiff = (yneigh(nk,nri) - yneigh(nl,nri))
        zdiff = (zneigh(nk,nri) - zneigh(nl,nri))
        xd = xdiff*d1(ind)
        yd = ydiff*d1(ind)
        zd = zdiff*d1(ind)
!
!  Compute cell numbers from matom to l and k
!
        imk1 = ineigh(1,nk,nri) - inm1
        imk2 = ineigh(2,nk,nri) - inm2
        imk3 = ineigh(3,nk,nri) - inm3
        iml1 = ineigh(1,nl,nri) - inm1
        iml2 = ineigh(2,nl,nri) - inm2
        iml3 = ineigh(3,nl,nri) - inm3
!
!  Compute cell indices for k
!
        if (abs(imk1).gt.nd2cell(1).or. &
            abs(imk2).gt.nd2cell(2).or. &
            abs(imk3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
          ncindp = nd2central
        else
!
!  Compute index
!
          ncindp = nd2cellptr(nd2cell(1)+1+imk1, &
                              nd2cell(2)+1+imk2, &
                              nd2cell(3)+1+imk3)
        endif
!
!  Compute cell indices for l
!
        if (abs(iml1).gt.nd2cell(1).or. &
            abs(iml2).gt.nd2cell(2).or. &
            abs(iml3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
          ncindm = nd2central
        else
!
!  Compute index
!
          ncindm = nd2cellptr(nd2cell(1)+1+iml1, &
                              nd2cell(2)+1+iml2, &
                              nd2cell(3)+1+iml3)
        endif
!
        d1cell(1,l,ncindm) = d1cell(1,l,ncindm) - xd
        d1cell(2,l,ncindm) = d1cell(2,l,ncindm) - yd
        d1cell(3,l,ncindm) = d1cell(3,l,ncindm) - zd
!
        d1cell(1,k,ncindp) = d1cell(1,k,ncindp) + xd
        d1cell(2,k,ncindp) = d1cell(2,k,ncindp) + yd
        d1cell(3,k,ncindp) = d1cell(3,k,ncindp) + zd
      endif
    enddo
  enddo
!
!  Loop over i -> neighbours of neighbours
!
  if (ltorderv) then
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh
      do nl = 1,nneigh(nrj)
        ind = ind + 1
        if (abs(d1(ind)).gt.d1tol) then
          l = neighno(nl,nrj)
          xdiff = xneigh(nl,nrj) + xneigh(nj,nri)
          ydiff = yneigh(nl,nrj) + yneigh(nj,nri)
          zdiff = zneigh(nl,nrj) + zneigh(nj,nri)
          xd = xdiff*d1(ind)
          yd = ydiff*d1(ind)
          zd = zdiff*d1(ind)
!
!  Compute cell numbers from matom to l and k
!
          imi1 = - inm1
          imi2 = - inm2
          imi3 = - inm3
          iml1 = ineigh(1,nl,nrj) + ineigh(1,nj,nri) - inm1
          iml2 = ineigh(2,nl,nrj) + ineigh(2,nj,nri) - inm2
          iml3 = ineigh(3,nl,nrj) + ineigh(3,nj,nri) - inm3
!
!  Compute cell indices for i
!         
          if (abs(imi1).gt.nd2cell(1).or. &
              abs(imi2).gt.nd2cell(2).or. &
              abs(imi3).gt.nd2cell(3)) then
!         
!  Find cell index - if outside user range then assign to central cell
!         
            ncindm = nd2central
          else
!   
!  Compute index
!  
            ncindm = nd2cellptr(nd2cell(1)+1+imi1, &
                                nd2cell(2)+1+imi2, &
                                nd2cell(3)+1+imi3)
          endif
!
!  Compute cell indices for l
!
          if (abs(iml1).gt.nd2cell(1).or. &
              abs(iml2).gt.nd2cell(2).or. &
              abs(iml3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
            ncindp = nd2central
          else
!
!  Compute index
!
            ncindp = nd2cellptr(nd2cell(1)+1+iml1, &
                                nd2cell(2)+1+iml2, &
                                nd2cell(3)+1+iml3)
          endif
!
          d1cell(1,i,ncindm) = d1cell(1,i,ncindm) - xd
          d1cell(2,i,ncindm) = d1cell(2,i,ncindm) - yd
          d1cell(3,i,ncindm) = d1cell(3,i,ncindm) - zd
!
          d1cell(1,l,ncindp) = d1cell(1,l,ncindp) + xd
          d1cell(2,l,ncindp) = d1cell(2,l,ncindp) + yd
          d1cell(3,l,ncindp) = d1cell(3,l,ncindp) + zd
        endif
      enddo
    enddo
!
!  Loop over neighbours -> neighbours of neighbours
!
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      do nk = 1,nneigh(nri)
        k = neighno(nk,nri)
        ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh + nk*mneigh
        do nl = 1,nneigh(nrj)
          ind = ind + 1
          if (abs(d1(ind)).gt.d1tol) then
            l = neighno(nl,nrj)
            xdiff = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
            ydiff = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
            zdiff = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
            xd = xdiff*d1(ind)
            yd = ydiff*d1(ind)
            zd = zdiff*d1(ind)
!
!  Compute cell numbers from matom to k and l
!
            imk1 = ineigh(1,nl,nrj) + ineigh(1,nj,nri) - ineigh(1,nk,nri) - inm1
            imk2 = ineigh(2,nl,nrj) + ineigh(2,nj,nri) - ineigh(2,nk,nri) - inm2
            imk3 = ineigh(3,nl,nrj) + ineigh(3,nj,nri) - ineigh(3,nk,nri) - inm3
            iml1 = ineigh(1,nl,nrj) + ineigh(1,nj,nri) - inm1
            iml2 = ineigh(2,nl,nrj) + ineigh(2,nj,nri) - inm2
            iml3 = ineigh(3,nl,nrj) + ineigh(3,nj,nri) - inm3
!
!  Compute cell indices for k
!           
            if (abs(imk1).gt.nd2cell(1).or. &
                abs(imk2).gt.nd2cell(2).or. &
                abs(imk3).gt.nd2cell(3)) then
!           
!  Find cell index - if outside user range then assign to central cell
!           
              ncindm = nd2central
            else
!
!  Compute index
!
              ncindm = nd2cellptr(nd2cell(1)+1+imk1, &
                                  nd2cell(2)+1+imk2, &
                                  nd2cell(3)+1+imk3)
            endif
!
!  Compute cell indices for l
!
            if (abs(iml1).gt.nd2cell(1).or. &
                abs(iml2).gt.nd2cell(2).or. &
                abs(iml3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
              ncindp = nd2central
            else
!
!  Compute index
!
              ncindp = nd2cellptr(nd2cell(1)+1+iml1, &
                                  nd2cell(2)+1+iml2, &
                                  nd2cell(3)+1+iml3)
            endif
!
            d1cell(1,k,ncindm) = d1cell(1,k,ncindm) - xd
            d1cell(2,k,ncindm) = d1cell(2,k,ncindm) - yd
            d1cell(3,k,ncindm) = d1cell(3,k,ncindm) - zd
!
            d1cell(1,l,ncindp) = d1cell(1,l,ncindp) + xd
            d1cell(2,l,ncindp) = d1cell(2,l,ncindp) + yd
            d1cell(3,l,ncindp) = d1cell(3,l,ncindp) + zd
          endif
        enddo
      enddo
    enddo
!
    if (lextended) then
!
!  Extended derivatives over third shell
!
      do nj = 1,nneigh(nri)
        j = neighno(nj,nri)
        nrj = nREBOatomRptr(j)
        do nk = 1,nneigh(nrj)
          k = neighno(nk,nrj)
          nrk = nREBOatomRptr(k)
          ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + nneigh(nri)*(mneigh+1)*mneigh + (nj-1)*mneigh*mneigh + (nk-1)*mneigh
          do nl = 1,nneigh(nrk)
            ind = ind + 1
            if (abs(d1(ind)).gt.d1tol) then
              l = neighno(nl,nrk)
              xdiff = xneigh(nl,nrk)
              ydiff = yneigh(nl,nrk)
              zdiff = zneigh(nl,nrk)
              xd = xdiff*d1(ind)
              yd = ydiff*d1(ind)
              zd = zdiff*d1(ind)
!
!  Compute cell numbers from matom to l and k
!
              imk1 = ineigh(1,nk,nrj) + ineigh(1,nj,nri) - inm1
              imk2 = ineigh(2,nk,nrj) + ineigh(2,nj,nri) - inm2
              imk3 = ineigh(3,nk,nrj) + ineigh(3,nj,nri) - inm3
              iml1 = ineigh(1,nl,nrk) + ineigh(1,nk,nrj) + ineigh(1,nj,nri) - inm1
              iml2 = ineigh(2,nl,nrk) + ineigh(2,nk,nrj) + ineigh(2,nj,nri) - inm2
              iml3 = ineigh(3,nl,nrk) + ineigh(3,nk,nrj) + ineigh(3,nj,nri) - inm3
!
!  Compute cell indices for k
!
              if (abs(imk1).gt.nd2cell(1).or. &
                  abs(imk2).gt.nd2cell(2).or. &
                  abs(imk3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncindp = nd2central
              else
!
!  Compute index
!
                ncindp = nd2cellptr(nd2cell(1)+1+imk1, &
                                    nd2cell(2)+1+imk2, &
                                    nd2cell(3)+1+imk3)
              endif
!
!  Compute cell indices for l
!
              if (abs(iml1).gt.nd2cell(1).or. &
                  abs(iml2).gt.nd2cell(2).or. &
                  abs(iml3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncindm = nd2central
              else
!
!  Compute index
!
                ncindm = nd2cellptr(nd2cell(1)+1+iml1, &
                                    nd2cell(2)+1+iml2, &
                                    nd2cell(3)+1+iml3)
              endif
!
              d1cell(1,k,ncindm) = d1cell(1,k,ncindm) - xd
              d1cell(2,k,ncindm) = d1cell(2,k,ncindm) - yd
              d1cell(3,k,ncindm) = d1cell(3,k,ncindm) - zd
!
              d1cell(1,l,ncindp) = d1cell(1,l,ncindp) + xd
              d1cell(2,l,ncindp) = d1cell(2,l,ncindp) + yd
              d1cell(3,l,ncindp) = d1cell(3,l,ncindp) + zd
            endif
          enddo
        enddo
      enddo
    endif
  endif
!
#ifdef TRACE
  call trace_out('d1addfc')
#endif
!
  return
  end
!
  subroutine d1addfcc(i,maxneigh,mneigh,nneigh,neighno,ineigh,maxlhs,d1cell,matom,nREBOatomRptr,d1,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian first derivative arrays. This version
!  is as per d1addfc except that Cartesian derivative components are
!  passed in rather than computed here.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  ineigh          = cell indices for vector to neighbour of atoms
!  matom           = atom being finite differenced
!  maxlhs          = second dimension of d1cell
!  d1cell          = array for storing forces separated by cell index
!  d1              = array of Cartesian first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!
!  11/17 Created from d1addfc
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: ineigh(3_i4,maxneigh,*)
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: maxlhs
  integer(i4), intent(in)             :: matom
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1(3,*)
  real(dp),    intent(out)            :: d1cell(4,maxlhs,*)
  logical,     intent(in)             :: ltorderv
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: imi1
  integer(i4)                         :: imi2
  integer(i4)                         :: imi3
  integer(i4)                         :: imk1
  integer(i4)                         :: imk2
  integer(i4)                         :: imk3
  integer(i4)                         :: iml1
  integer(i4)                         :: iml2
  integer(i4)                         :: iml3
  integer(i4)                         :: inm1
  integer(i4)                         :: inm2
  integer(i4)                         :: inm3
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: l
  integer(i4)                         :: ncindm
  integer(i4)                         :: ncindp
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  logical                             :: lfound
  real(dp)                            :: xd
  real(dp)                            :: yd
  real(dp)                            :: zd
#ifdef TRACE
  call trace_in('d1addfcc')
#endif
!
!  Map atom i to relevant atom in arrays
!
  nri = nREBOatomRptr(i)
!
!  Find cell indices for matom - first search i and neighbours of i
!
  inm1 = 0
  inm2 = 0
  inm3 = 0
  if (i.ne.matom) then
    lfound = .false.
    nj = 0
    do while (.not.lfound.and.nj.lt.nneigh(nri))
      nj = nj + 1
      lfound = (neighno(nj,nri).eq.matom)
    enddo
    if (lfound) then
      inm1 = ineigh(1,nj,nri) 
      inm2 = ineigh(2,nj,nri) 
      inm3 = ineigh(3,nj,nri) 
    else
!
!  If matom hasn't been found then search neighbours of neighbours of i
!
      nj = 0
      do while (.not.lfound.and.nj.lt.nneigh(nri))
        nj = nj + 1
        nrj = nREBOatomRptr(neighno(nj,nri))
        nk = 0
        do while (.not.lfound.and.nk.lt.nneigh(nrj))
          nk = nk + 1
          lfound = (neighno(nk,nrj).eq.matom)
        enddo
      enddo
      if (lfound) then
        inm1 = ineigh(1,nk,nrj) + ineigh(1,nj,nri)
        inm2 = ineigh(2,nk,nrj) + ineigh(2,nj,nri)
        inm3 = ineigh(3,nk,nrj) + ineigh(3,nj,nri)
      endif
    endif
  endif
!
!  Loop over i -> neighbours
!
  ind = 0
  do nk = 1,nneigh(nri)
    ind = ind + 1
    k = neighno(nk,nri)
    xd = d1(1,ind)
    yd = d1(2,ind)
    zd = d1(3,ind)
!
!  Compute cell numbers from matom to i and k
!
    imi1 = - inm1
    imi2 = - inm2
    imi3 = - inm3
    imk1 = ineigh(1,nk,nri) - inm1
    imk2 = ineigh(2,nk,nri) - inm2
    imk3 = ineigh(3,nk,nri) - inm3
!
!  Compute cell indices for i
!
    if (abs(imi1).gt.nd2cell(1).or. &
        abs(imi2).gt.nd2cell(2).or. &
        abs(imi3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
      ncindm = nd2central
    else
!
!  Compute index
!
      ncindm = nd2cellptr(nd2cell(1)+1+imi1, &
                          nd2cell(2)+1+imi2, &
                          nd2cell(3)+1+imi3)
    endif
!
!  Compute cell indices for k
!
    if (abs(imk1).gt.nd2cell(1).or. &
        abs(imk2).gt.nd2cell(2).or. &
        abs(imk3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
      ncindp = nd2central
    else
!
!  Compute index
!
      ncindp = nd2cellptr(nd2cell(1)+1+imk1, &
                          nd2cell(2)+1+imk2, &
                          nd2cell(3)+1+imk3)
    endif
!
    d1cell(1,i,ncindm) = d1cell(1,i,ncindm) - xd
    d1cell(2,i,ncindm) = d1cell(2,i,ncindm) - yd
    d1cell(3,i,ncindm) = d1cell(3,i,ncindm) - zd
!
    d1cell(1,k,ncindp) = d1cell(1,k,ncindp) + xd
    d1cell(2,k,ncindp) = d1cell(2,k,ncindp) + yd
    d1cell(3,k,ncindp) = d1cell(3,k,ncindp) + zd
  enddo
!
!  Loop over neighbours -> neighbours
!
  do nk = 2,nneigh(nri)
    k = neighno(nk,nri)
    do nl = 1,nk - 1
      ind = nneigh(nri) + nk*(nk-1)/2 + nl
      if (abs(d1(1,ind)).gt.d1tol.or.abs(d1(2,ind)).gt.d1tol.or.abs(d1(3,ind)).gt.d1tol) then
        l = neighno(nl,nri)
        xd = d1(1,ind)
        yd = d1(2,ind)
        zd = d1(3,ind)
!
!  Compute cell numbers from matom to l and k
!
        imk1 = ineigh(1,nk,nri) - inm1
        imk2 = ineigh(2,nk,nri) - inm2
        imk3 = ineigh(3,nk,nri) - inm3
        iml1 = ineigh(1,nl,nri) - inm1
        iml2 = ineigh(2,nl,nri) - inm2
        iml3 = ineigh(3,nl,nri) - inm3
!
!  Compute cell indices for k
!
        if (abs(imk1).gt.nd2cell(1).or. &
            abs(imk2).gt.nd2cell(2).or. &
            abs(imk3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
          ncindp = nd2central
        else
!
!  Compute index
!
          ncindp = nd2cellptr(nd2cell(1)+1+imk1, &
                              nd2cell(2)+1+imk2, &
                              nd2cell(3)+1+imk3)
        endif
!
!  Compute cell indices for l
!
        if (abs(iml1).gt.nd2cell(1).or. &
            abs(iml2).gt.nd2cell(2).or. &
            abs(iml3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
          ncindm = nd2central
        else
!
!  Compute index
!
          ncindm = nd2cellptr(nd2cell(1)+1+iml1, &
                              nd2cell(2)+1+iml2, &
                              nd2cell(3)+1+iml3)
        endif
!
        d1cell(1,l,ncindm) = d1cell(1,l,ncindm) - xd
        d1cell(2,l,ncindm) = d1cell(2,l,ncindm) - yd
        d1cell(3,l,ncindm) = d1cell(3,l,ncindm) - zd
!
        d1cell(1,k,ncindp) = d1cell(1,k,ncindp) + xd
        d1cell(2,k,ncindp) = d1cell(2,k,ncindp) + yd
        d1cell(3,k,ncindp) = d1cell(3,k,ncindp) + zd
      endif
    enddo
  enddo
!
!  Loop over i -> neighbours of neighbours
!
  if (ltorderv) then
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh
      do nl = 1,nneigh(nrj)
        ind = ind + 1
        if (abs(d1(1,ind)).gt.d1tol.or.abs(d1(2,ind)).gt.d1tol.or.abs(d1(3,ind)).gt.d1tol) then
          l = neighno(nl,nrj)
          xd = d1(1,ind)
          yd = d1(2,ind)
          zd = d1(3,ind)
!
!  Compute cell numbers from matom to l and k
!
          imi1 = - inm1
          imi2 = - inm2
          imi3 = - inm3
          iml1 = ineigh(1,nl,nrj) + ineigh(1,nj,nri) - inm1
          iml2 = ineigh(2,nl,nrj) + ineigh(2,nj,nri) - inm2
          iml3 = ineigh(3,nl,nrj) + ineigh(3,nj,nri) - inm3
!
!  Compute cell indices for i
!         
          if (abs(imi1).gt.nd2cell(1).or. &
              abs(imi2).gt.nd2cell(2).or. &
              abs(imi3).gt.nd2cell(3)) then
!         
!  Find cell index - if outside user range then assign to central cell
!         
            ncindm = nd2central
          else
!   
!  Compute index
!  
            ncindm = nd2cellptr(nd2cell(1)+1+imi1, &
                                nd2cell(2)+1+imi2, &
                                nd2cell(3)+1+imi3)
          endif
!
!  Compute cell indices for l
!
          if (abs(iml1).gt.nd2cell(1).or. &
              abs(iml2).gt.nd2cell(2).or. &
              abs(iml3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
            ncindp = nd2central
          else
!
!  Compute index
!
            ncindp = nd2cellptr(nd2cell(1)+1+iml1, &
                                nd2cell(2)+1+iml2, &
                                nd2cell(3)+1+iml3)
          endif
!
          d1cell(1,i,ncindm) = d1cell(1,i,ncindm) - xd
          d1cell(2,i,ncindm) = d1cell(2,i,ncindm) - yd
          d1cell(3,i,ncindm) = d1cell(3,i,ncindm) - zd
!
          d1cell(1,l,ncindp) = d1cell(1,l,ncindp) + xd
          d1cell(2,l,ncindp) = d1cell(2,l,ncindp) + yd
          d1cell(3,l,ncindp) = d1cell(3,l,ncindp) + zd
        endif
      enddo
    enddo
!
!  Loop over neighbours -> neighbours of neighbours
!
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      do nk = 1,nneigh(nri)
        k = neighno(nk,nri)
        ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh + nk*mneigh
        do nl = 1,nneigh(nrj)
          ind = ind + 1
          if (abs(d1(1,ind)).gt.d1tol.or.abs(d1(2,ind)).gt.d1tol.or.abs(d1(3,ind)).gt.d1tol) then
            l = neighno(nl,nrj)
            xd = d1(1,ind)
            yd = d1(2,ind)
            zd = d1(3,ind)
!
!  Compute cell numbers from matom to k and l
!
            imk1 = ineigh(1,nl,nrj) + ineigh(1,nj,nri) - ineigh(1,nk,nri) - inm1
            imk2 = ineigh(2,nl,nrj) + ineigh(2,nj,nri) - ineigh(2,nk,nri) - inm2
            imk3 = ineigh(3,nl,nrj) + ineigh(3,nj,nri) - ineigh(3,nk,nri) - inm3
            iml1 = ineigh(1,nl,nrj) + ineigh(1,nj,nri) - inm1
            iml2 = ineigh(2,nl,nrj) + ineigh(2,nj,nri) - inm2
            iml3 = ineigh(3,nl,nrj) + ineigh(3,nj,nri) - inm3
!
!  Compute cell indices for k
!           
            if (abs(imk1).gt.nd2cell(1).or. &
                abs(imk2).gt.nd2cell(2).or. &
                abs(imk3).gt.nd2cell(3)) then
!           
!  Find cell index - if outside user range then assign to central cell
!           
              ncindm = nd2central
            else
!
!  Compute index
!
              ncindm = nd2cellptr(nd2cell(1)+1+imk1, &
                                  nd2cell(2)+1+imk2, &
                                  nd2cell(3)+1+imk3)
            endif
!
!  Compute cell indices for l
!
            if (abs(iml1).gt.nd2cell(1).or. &
                abs(iml2).gt.nd2cell(2).or. &
                abs(iml3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
              ncindp = nd2central
            else
!
!  Compute index
!
              ncindp = nd2cellptr(nd2cell(1)+1+iml1, &
                                  nd2cell(2)+1+iml2, &
                                  nd2cell(3)+1+iml3)
            endif
!
            d1cell(1,k,ncindm) = d1cell(1,k,ncindm) - xd
            d1cell(2,k,ncindm) = d1cell(2,k,ncindm) - yd
            d1cell(3,k,ncindm) = d1cell(3,k,ncindm) - zd
!
            d1cell(1,l,ncindp) = d1cell(1,l,ncindp) + xd
            d1cell(2,l,ncindp) = d1cell(2,l,ncindp) + yd
            d1cell(3,l,ncindp) = d1cell(3,l,ncindp) + zd
          endif
        enddo
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('d1addfcc')
#endif
!
  return
  end
!
  subroutine d1addc(i,maxneigh,mneigh,nneigh,neighno,nfreeatom,nREBOatomRptr,d1,d1s,ltorderv,ldoregions)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian first derivative arrays. This version
!  is as per d1add except that Cartesian derivative components are
!  passed in rather than computed here.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  nfreeatom       = pointer from atom to free atom to move
!  d1              = array of first derivatives w.r.t. Cartesian components
!  d1s             = array of first derivatives w.r.t. strains
!  ltorderv        = if true then include torsional derivatives
!  ldoregions      = if true then region derivatives should be computed
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!
!  10/10 Created from d1add
!   5/11 Check on d1 being > 0 reinstated using sum of components of d1
!   4/12 Explicit virial calculation removed as no longer needed
!   5/12 Atomic stresses added
!   4/16 ldoregions flag added
!   2/18 Trace added
!   9/18 Strain module added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   5/23 Derivative tolerances as global variables introduced
!   7/23 nregi now defined at the start
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
  use brennerdata
  use control,        only : latomicstress
  use current,        only : ndim, nrelf2a, nregionno
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1(3,*)
  real(dp),    intent(in)             :: d1s(6,*)
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: ldoregions
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: io
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: ko
  integer(i4)                         :: l
  integer(i4)                         :: lo
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nregi
  integer(i4)                         :: nregk
  integer(i4)                         :: nregl
  real(dp)                            :: d1sum
  real(dp)                            :: xd
  real(dp)                            :: yd
  real(dp)                            :: zd
#ifdef TRACE
  call trace_in('d1addc')
#endif
!
!  Loop over i -> neighbours
!
  ind = 0
  nri = nREBOatomRptr(i)
  nregi = nregionno(nrelf2a(i))
  do nk = 1,nneigh(nri)
    ind = ind + 1
    k = neighno(nk,nri)
    xd = d1(1,ind)
    yd = d1(2,ind)
    zd = d1(3,ind)
    io = nfreeatom(i)
    ko = nfreeatom(k)
    if (io.gt.0) then
      xdrv(i) = xdrv(i) - xd
      ydrv(i) = ydrv(i) - yd
      zdrv(i) = zdrv(i) - zd
    endif
    if (ko.gt.0) then
      xdrv(k) = xdrv(k) + xd
      ydrv(k) = ydrv(k) + yd
      zdrv(k) = zdrv(k) + zd
    endif
    if (ldoregions) then
      nregk = nregionno(nrelf2a(k))
      if (nregi.ne.nregk) then
        xregdrv(nregi) = xregdrv(nregi) - xd
        yregdrv(nregi) = yregdrv(nregi) - yd
        zregdrv(nregi) = zregdrv(nregi) - zd
        xregdrv(nregk) = xregdrv(nregk) + xd
        yregdrv(nregk) = yregdrv(nregk) + yd
        zregdrv(nregk) = zregdrv(nregk) + zd
      endif
    endif
    if (lstr) then
      select case(ndim)
        case(1)
          rstrd(1) = rstrd(1) + d1s(1,ind)
        case(2)
          rstrd(1) = rstrd(1) + d1s(1,ind)
          rstrd(2) = rstrd(2) + d1s(2,ind)
          rstrd(3) = rstrd(3) + d1s(3,ind)
        case(3)
          rstrd(1) = rstrd(1) + d1s(1,ind)
          rstrd(2) = rstrd(2) + d1s(2,ind)
          rstrd(3) = rstrd(3) + d1s(3,ind)
          rstrd(4) = rstrd(4) + d1s(4,ind)
          rstrd(5) = rstrd(5) + d1s(5,ind)
          rstrd(6) = rstrd(6) + d1s(6,ind)
      end select
      if (latomicstress) then
        select case(ndim)
          case(1)
            atomicstress(1,i) = atomicstress(1,i) + 0.5_dp*d1s(1,ind)
            atomicstress(1,k) = atomicstress(1,k) + 0.5_dp*d1s(1,ind)
          case(2)
            atomicstress(1,i) = atomicstress(1,i) + 0.5_dp*d1s(1,ind)
            atomicstress(2,i) = atomicstress(2,i) + 0.5_dp*d1s(2,ind)
            atomicstress(3,i) = atomicstress(3,i) + 0.5_dp*d1s(3,ind)
            atomicstress(1,k) = atomicstress(1,k) + 0.5_dp*d1s(1,ind)
            atomicstress(2,k) = atomicstress(2,k) + 0.5_dp*d1s(2,ind)
            atomicstress(3,k) = atomicstress(3,k) + 0.5_dp*d1s(3,ind)
          case(3)
            atomicstress(1,i) = atomicstress(1,i) + 0.5_dp*d1s(1,ind)
            atomicstress(2,i) = atomicstress(2,i) + 0.5_dp*d1s(2,ind)
            atomicstress(3,i) = atomicstress(3,i) + 0.5_dp*d1s(3,ind)
            atomicstress(4,i) = atomicstress(4,i) + 0.5_dp*d1s(4,ind)
            atomicstress(5,i) = atomicstress(5,i) + 0.5_dp*d1s(5,ind)
            atomicstress(6,i) = atomicstress(6,i) + 0.5_dp*d1s(6,ind)
            atomicstress(1,k) = atomicstress(1,k) + 0.5_dp*d1s(1,ind)
            atomicstress(2,k) = atomicstress(2,k) + 0.5_dp*d1s(2,ind)
            atomicstress(3,k) = atomicstress(3,k) + 0.5_dp*d1s(3,ind)
            atomicstress(4,k) = atomicstress(4,k) + 0.5_dp*d1s(4,ind)
            atomicstress(5,k) = atomicstress(5,k) + 0.5_dp*d1s(5,ind)
            atomicstress(6,k) = atomicstress(6,k) + 0.5_dp*d1s(6,ind)
        end select
      endif
    endif
  enddo
!
!  Loop over neighbours -> neighbours
!
  do nk = 2,nneigh(nri)
    k = neighno(nk,nri)
    do nl = 1,nk - 1
      ind = nneigh(nri) + nk*(nk-1)/2 + nl
      d1sum = abs(d1(1,ind)) + abs(d1(2,ind)) + abs(d1(3,ind))
      if (d1sum.gt.d1tol) then
        l = neighno(nl,nri)
        xd = d1(1,ind)
        yd = d1(2,ind)
        zd = d1(3,ind)
        lo = nfreeatom(l)
        ko = nfreeatom(k)
        if (lo.gt.0) then
          xdrv(l) = xdrv(l) - xd
          ydrv(l) = ydrv(l) - yd
          zdrv(l) = zdrv(l) - zd
        endif
        if (ko.gt.0) then
          xdrv(k) = xdrv(k) + xd
          ydrv(k) = ydrv(k) + yd
          zdrv(k) = zdrv(k) + zd
        endif
!
        if (ldoregions) then
          nregl = nregionno(nrelf2a(l))
          nregk = nregionno(nrelf2a(k))
          if (nregl.ne.nregk) then
            xregdrv(nregl) = xregdrv(nregl) - xd
            yregdrv(nregl) = yregdrv(nregl) - yd
            zregdrv(nregl) = zregdrv(nregl) - zd
            xregdrv(nregk) = xregdrv(nregk) + xd
            yregdrv(nregk) = yregdrv(nregk) + yd
            zregdrv(nregk) = zregdrv(nregk) + zd
          endif
        endif
        if (lstr) then
          select case(ndim)
            case(1)
              rstrd(1) = rstrd(1) + d1s(1,ind)
            case(2)
              rstrd(1) = rstrd(1) + d1s(1,ind)
              rstrd(2) = rstrd(2) + d1s(2,ind)
              rstrd(3) = rstrd(3) + d1s(3,ind)
            case(3)
              rstrd(1) = rstrd(1) + d1s(1,ind)
              rstrd(2) = rstrd(2) + d1s(2,ind)
              rstrd(3) = rstrd(3) + d1s(3,ind)
              rstrd(4) = rstrd(4) + d1s(4,ind)
              rstrd(5) = rstrd(5) + d1s(5,ind)
              rstrd(6) = rstrd(6) + d1s(6,ind)
          end select
          if (latomicstress) then
            select case(ndim)
              case(1)
                atomicstress(1,l) = atomicstress(1,l) + 0.5_dp*d1s(1,ind)
                atomicstress(1,k) = atomicstress(1,k) + 0.5_dp*d1s(1,ind)
              case(2)
                atomicstress(1,l) = atomicstress(1,l) + 0.5_dp*d1s(1,ind)
                atomicstress(2,l) = atomicstress(2,l) + 0.5_dp*d1s(2,ind)
                atomicstress(3,l) = atomicstress(3,l) + 0.5_dp*d1s(3,ind)
                atomicstress(1,k) = atomicstress(1,k) + 0.5_dp*d1s(1,ind)
                atomicstress(2,k) = atomicstress(2,k) + 0.5_dp*d1s(2,ind)
                atomicstress(3,k) = atomicstress(3,k) + 0.5_dp*d1s(3,ind)
              case(3)
                atomicstress(1,l) = atomicstress(1,l) + 0.5_dp*d1s(1,ind)
                atomicstress(2,l) = atomicstress(2,l) + 0.5_dp*d1s(2,ind)
                atomicstress(3,l) = atomicstress(3,l) + 0.5_dp*d1s(3,ind)
                atomicstress(4,l) = atomicstress(4,l) + 0.5_dp*d1s(4,ind)
                atomicstress(5,l) = atomicstress(5,l) + 0.5_dp*d1s(5,ind)
                atomicstress(6,l) = atomicstress(6,l) + 0.5_dp*d1s(6,ind)
                atomicstress(1,k) = atomicstress(1,k) + 0.5_dp*d1s(1,ind)
                atomicstress(2,k) = atomicstress(2,k) + 0.5_dp*d1s(2,ind)
                atomicstress(3,k) = atomicstress(3,k) + 0.5_dp*d1s(3,ind)
                atomicstress(4,k) = atomicstress(4,k) + 0.5_dp*d1s(4,ind)
                atomicstress(5,k) = atomicstress(5,k) + 0.5_dp*d1s(5,ind)
                atomicstress(6,k) = atomicstress(6,k) + 0.5_dp*d1s(6,ind)
            end select
          endif
        endif
      endif
    enddo
  enddo
!
!  Loop over i -> neighbours of neighbours
!
  if (ltorderv) then
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh
      do nl = 1,nneigh(nrj)
        ind = ind + 1
        d1sum = abs(d1(1,ind)) + abs(d1(2,ind)) + abs(d1(3,ind))
        if (d1sum.gt.d1tol) then
          l = neighno(nl,nrj)
          xd = d1(1,ind)
          yd = d1(2,ind)
          zd = d1(3,ind)
          io = nfreeatom(i)
          lo = nfreeatom(l)
          if (io.gt.0) then
            xdrv(i) = xdrv(i) - xd
            ydrv(i) = ydrv(i) - yd
            zdrv(i) = zdrv(i) - zd
          endif
          if (lo.gt.0) then
            xdrv(l) = xdrv(l) + xd
            ydrv(l) = ydrv(l) + yd
            zdrv(l) = zdrv(l) + zd
          endif
!
          if (ldoregions) then
            nregi = nregionno(nrelf2a(i))
            nregl = nregionno(nrelf2a(l))
            if (nregi.ne.nregl) then
              xregdrv(nregi) = xregdrv(nregi) - xd
              yregdrv(nregi) = yregdrv(nregi) - yd
              zregdrv(nregi) = zregdrv(nregi) - zd
              xregdrv(nregl) = xregdrv(nregl) + xd
              yregdrv(nregl) = yregdrv(nregl) + yd
              zregdrv(nregl) = zregdrv(nregl) + zd
            endif
          endif
!
          if (lstr) then
            select case(ndim)
              case(1)
                rstrd(1) = rstrd(1) + d1s(1,ind)
              case(2)
                rstrd(1) = rstrd(1) + d1s(1,ind)
                rstrd(2) = rstrd(2) + d1s(2,ind)
                rstrd(3) = rstrd(3) + d1s(3,ind)
              case(3)
                rstrd(1) = rstrd(1) + d1s(1,ind)
                rstrd(2) = rstrd(2) + d1s(2,ind)
                rstrd(3) = rstrd(3) + d1s(3,ind)
                rstrd(4) = rstrd(4) + d1s(4,ind)
                rstrd(5) = rstrd(5) + d1s(5,ind)
                rstrd(6) = rstrd(6) + d1s(6,ind)
            end select
            if (latomicstress) then
              select case(ndim)
                case(1)
                  atomicstress(1,i) = atomicstress(1,i) + 0.5_dp*d1s(1,ind)
                  atomicstress(1,l) = atomicstress(1,l) + 0.5_dp*d1s(1,ind)
                case(2)
                  atomicstress(1,i) = atomicstress(1,i) + 0.5_dp*d1s(1,ind)
                  atomicstress(2,i) = atomicstress(2,i) + 0.5_dp*d1s(2,ind)
                  atomicstress(3,i) = atomicstress(3,i) + 0.5_dp*d1s(3,ind)
                  atomicstress(1,l) = atomicstress(1,l) + 0.5_dp*d1s(1,ind)
                  atomicstress(2,l) = atomicstress(2,l) + 0.5_dp*d1s(2,ind)
                  atomicstress(3,l) = atomicstress(3,l) + 0.5_dp*d1s(3,ind)
                case(3)
                  atomicstress(1,i) = atomicstress(1,i) + 0.5_dp*d1s(1,ind)
                  atomicstress(2,i) = atomicstress(2,i) + 0.5_dp*d1s(2,ind)
                  atomicstress(3,i) = atomicstress(3,i) + 0.5_dp*d1s(3,ind)
                  atomicstress(4,i) = atomicstress(4,i) + 0.5_dp*d1s(4,ind)
                  atomicstress(5,i) = atomicstress(5,i) + 0.5_dp*d1s(5,ind)
                  atomicstress(6,i) = atomicstress(6,i) + 0.5_dp*d1s(6,ind)
                  atomicstress(1,l) = atomicstress(1,l) + 0.5_dp*d1s(1,ind)
                  atomicstress(2,l) = atomicstress(2,l) + 0.5_dp*d1s(2,ind)
                  atomicstress(3,l) = atomicstress(3,l) + 0.5_dp*d1s(3,ind)
                  atomicstress(4,l) = atomicstress(4,l) + 0.5_dp*d1s(4,ind)
                  atomicstress(5,l) = atomicstress(5,l) + 0.5_dp*d1s(5,ind)
                  atomicstress(6,l) = atomicstress(6,l) + 0.5_dp*d1s(6,ind)
              end select
            endif
          endif
        endif
      enddo
    enddo
!
!  Loop over neighbours -> neighbours of neighbours
!
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      do nk = 1,nneigh(nri)
        k = neighno(nk,nri)
        ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh + nk*mneigh
        do nl = 1,nneigh(nrj)
          ind = ind + 1
          d1sum = abs(d1(1,ind)) + abs(d1(2,ind)) + abs(d1(3,ind))
          if (d1sum.gt.d1tol) then
            l = neighno(nl,nrj)
            xd = d1(1,ind)
            yd = d1(2,ind)
            zd = d1(3,ind)
            ko = nfreeatom(k)
            lo = nfreeatom(l)
            if (ko.gt.0) then
              xdrv(k) = xdrv(k) - xd
              ydrv(k) = ydrv(k) - yd
              zdrv(k) = zdrv(k) - zd
            endif
            if (lo.gt.0) then
              xdrv(l) = xdrv(l) + xd
              ydrv(l) = ydrv(l) + yd
              zdrv(l) = zdrv(l) + zd
            endif
!
            if (ldoregions) then
              nregk = nregionno(nrelf2a(k))
              nregl = nregionno(nrelf2a(l))
              if (nregk.ne.nregl) then
                xregdrv(nregk) = xregdrv(nregk) - xd
                yregdrv(nregk) = yregdrv(nregk) - yd
                zregdrv(nregk) = zregdrv(nregk) - zd
                xregdrv(nregl) = xregdrv(nregl) + xd
                yregdrv(nregl) = yregdrv(nregl) + yd
                zregdrv(nregl) = zregdrv(nregl) + zd
              endif
            endif
!
            if (lstr) then
              select case(ndim)
                case(1)
                  rstrd(1) = rstrd(1) + d1s(1,ind)
                case(2)
                  rstrd(1) = rstrd(1) + d1s(1,ind)
                  rstrd(2) = rstrd(2) + d1s(2,ind)
                  rstrd(3) = rstrd(3) + d1s(3,ind)
                case(3)
                  rstrd(1) = rstrd(1) + d1s(1,ind)
                  rstrd(2) = rstrd(2) + d1s(2,ind)
                  rstrd(3) = rstrd(3) + d1s(3,ind)
                  rstrd(4) = rstrd(4) + d1s(4,ind)
                  rstrd(5) = rstrd(5) + d1s(5,ind)
                  rstrd(6) = rstrd(6) + d1s(6,ind)
              end select
              if (latomicstress) then
                select case(ndim)
                  case(1)
                    atomicstress(1,k) = atomicstress(1,k) + 0.5_dp*d1s(1,ind)
                    atomicstress(1,l) = atomicstress(1,l) + 0.5_dp*d1s(1,ind)
                  case(2)
                    atomicstress(1,k) = atomicstress(1,k) + 0.5_dp*d1s(1,ind)
                    atomicstress(2,k) = atomicstress(2,k) + 0.5_dp*d1s(2,ind)
                    atomicstress(3,k) = atomicstress(3,k) + 0.5_dp*d1s(3,ind)
                    atomicstress(1,l) = atomicstress(1,l) + 0.5_dp*d1s(1,ind)
                    atomicstress(2,l) = atomicstress(2,l) + 0.5_dp*d1s(2,ind)
                    atomicstress(3,l) = atomicstress(3,l) + 0.5_dp*d1s(3,ind)
                  case(3)
                    atomicstress(1,k) = atomicstress(1,k) + 0.5_dp*d1s(1,ind)
                    atomicstress(2,k) = atomicstress(2,k) + 0.5_dp*d1s(2,ind)
                    atomicstress(3,k) = atomicstress(3,k) + 0.5_dp*d1s(3,ind)
                    atomicstress(4,k) = atomicstress(4,k) + 0.5_dp*d1s(4,ind)
                    atomicstress(5,k) = atomicstress(5,k) + 0.5_dp*d1s(5,ind)
                    atomicstress(6,k) = atomicstress(6,k) + 0.5_dp*d1s(6,ind)
                    atomicstress(1,l) = atomicstress(1,l) + 0.5_dp*d1s(1,ind)
                    atomicstress(2,l) = atomicstress(2,l) + 0.5_dp*d1s(2,ind)
                    atomicstress(3,l) = atomicstress(3,l) + 0.5_dp*d1s(3,ind)
                    atomicstress(4,l) = atomicstress(4,l) + 0.5_dp*d1s(4,ind)
                    atomicstress(5,l) = atomicstress(5,l) + 0.5_dp*d1s(5,ind)
                    atomicstress(6,l) = atomicstress(6,l) + 0.5_dp*d1s(6,ind)
                end select
              endif
            endif
          endif
        enddo
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('d1addc')
#endif
!
  return
  end
!
  subroutine d1addnn(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1nn,ldoregions)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian first derivative arrays.
!  This version has a different structure to d1add suited for interactions that are for 
!  nearest neighbours only.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  d1nn            = array of first derivatives w.r.t. interatomic distances
!  ldoregions      = if true then region derivatives should be computed
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!
!   8/22 Created
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use control,        only : latomicstress
  use current,        only : ndim, nrelf2a, nregionno
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1nn(*)
  logical,     intent(in)             :: ldoregions
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: io
  integer(i4)                         :: is
  integer(i4)                         :: k
  integer(i4)                         :: ko
  integer(i4)                         :: l
  integer(i4)                         :: lo
  integer(i4)                         :: nri
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nregi
  integer(i4)                         :: nregk
  integer(i4)                         :: nregl
  integer(i4)                         :: ns
  real(dp)                            :: dr2ds(6)
  real(dp)                            :: d2r2dx2(3,3)
  real(dp)                            :: d2r2ds2(6,6)
  real(dp)                            :: d2r2dsdx(6,3)
  real(dp)                            :: xd
  real(dp)                            :: yd
  real(dp)                            :: zd
  real(dp)                            :: xdiff
  real(dp)                            :: ydiff
  real(dp)                            :: zdiff
#ifdef TRACE
  call trace_in('d1addnn')
#endif
!
!  Loop over i -> neighbours
!
  ind = 0
  nri = nREBOatomRptr(i)
  do nk = 1,nneigh(nri)
    ind = ind + 1
    k = neighno(nk,nri)
    xd = xneigh(nk,nri)*d1nn(ind)
    yd = yneigh(nk,nri)*d1nn(ind)
    zd = zneigh(nk,nri)*d1nn(ind)
    io = nfreeatom(i)
    ko = nfreeatom(k)
    if (io.gt.0) then
      xdrv(i) = xdrv(i) - xd
      ydrv(i) = ydrv(i) - yd
      zdrv(i) = zdrv(i) - zd
    endif
    if (ko.gt.0) then
      xdrv(k) = xdrv(k) + xd
      ydrv(k) = ydrv(k) + yd
      zdrv(k) = zdrv(k) + zd
    endif
    if (ldoregions) then
      nregi = nregionno(nrelf2a(i))
      nregk = nregionno(nrelf2a(k))
      if (nregi.ne.nregk) then
        xregdrv(nregi) = xregdrv(nregi) - xd
        yregdrv(nregi) = yregdrv(nregi) - yd
        zregdrv(nregi) = zregdrv(nregi) - zd
        xregdrv(nregk) = xregdrv(nregk) + xd
        yregdrv(nregk) = yregdrv(nregk) + yd
        zregdrv(nregk) = zregdrv(nregk) + zd
      endif
    endif
    if (lstr) then
      call real1strterm(ndim,xneigh(nk,nri),yneigh(nk,nri),zneigh(nk,nri),0.0_dp,0.0_dp,0.0_dp, &
        dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
      do is = 1,nstrains
        ns = nstrptr(is)
        rstrd(is) = rstrd(is) + d1nn(ind)*dr2ds(ns)
      enddo
      if (latomicstress) then
        do is = 1,nstrains
          ns = nstrptr(is)
          atomicstress(is,i) = atomicstress(is,i) + 0.5_dp*d1nn(ind)*dr2ds(ns)
          atomicstress(is,k) = atomicstress(is,k) + 0.5_dp*d1nn(ind)*dr2ds(ns)
        enddo
      endif
    endif
  enddo
!
!  Loop over neighbours of i -> neighbours of neighbours of i
!
  do nk = 1,nneigh(nri)
    k = neighno(nk,nri)
    do nl = 1,nneigh(k)
      ind = nneigh(nri) + (nk-1)*maxneigh + nl
      if (abs(d1nn(ind)).gt.d1tol) then
        l = neighno(nl,k)
        xdiff = - xneigh(nl,k)
        ydiff = - yneigh(nl,k)
        zdiff = - zneigh(nl,k)
        xd = xdiff*d1nn(ind)
        yd = ydiff*d1nn(ind)
        zd = zdiff*d1nn(ind)
        lo = nfreeatom(l)
        ko = nfreeatom(k)
        if (lo.gt.0) then
          xdrv(l) = xdrv(l) - xd
          ydrv(l) = ydrv(l) - yd
          zdrv(l) = zdrv(l) - zd
        endif
        if (ko.gt.0) then
          xdrv(k) = xdrv(k) + xd
          ydrv(k) = ydrv(k) + yd
          zdrv(k) = zdrv(k) + zd
        endif
!
        if (ldoregions) then
          nregl = nregionno(nrelf2a(l))
          nregk = nregionno(nrelf2a(k))
          if (nregl.ne.nregk) then
            xregdrv(nregl) = xregdrv(nregl) - xd
            yregdrv(nregl) = yregdrv(nregl) - yd
            zregdrv(nregl) = zregdrv(nregl) - zd
            xregdrv(nregk) = xregdrv(nregk) + xd
            yregdrv(nregk) = yregdrv(nregk) + yd
            zregdrv(nregk) = zregdrv(nregk) + zd
          endif
        endif
        if (lstr) then
          call real1strterm(ndim,xdiff,ydiff,zdiff,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do is = 1,nstrains
            ns = nstrptr(is)
            rstrd(is) = rstrd(is) + d1nn(ind)*dr2ds(ns)
          enddo
          if (latomicstress) then
            do is = 1,nstrains
              ns = nstrptr(is)
              atomicstress(is,l) = atomicstress(is,l) + 0.5_dp*d1nn(ind)*dr2ds(ns)
              atomicstress(is,k) = atomicstress(is,k) + 0.5_dp*d1nn(ind)*dr2ds(ns)
            enddo
          endif
        endif
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('d1addnn')
#endif
!
  return
  end
!
  subroutine d1adds(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,nauatom,neqvatom,d1, &
                    ltorderv,ldoregions)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian first derivative arrays.
!
!  Symmetry adapted version of d1add
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  nauatom         = pointer from full cell atom to asymmetric unit
!  neqvatom        = number of equivalent positions for asymmetric unit atom
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!  ldoregions      = if true then region derivatives should be computed
!
!   9/02 Created from d1add
!  12/03 Torsional derivatives made optional
!   6/07 nREBOatomRptr added as input
!   4/08 mneigh added as a distinct number from maxneigh
!  11/09 Region derivatives added
!   3/10 Bug in handling of region derivatives when atom number is 0 fixed
!   4/12 Explicit virial calculation removed as no longer needed
!   5/12 Atomic stresses added
!   5/12 Atomic stresses removed for routines involving symmetry
!   4/16 ldoregions flag added
!   2/18 Trace added
!   9/18 Strain module added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use current,        only : ndim, nregionno
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  integer(i4), intent(in)             :: nauatom(*)
  integer(i4), intent(in)             :: neqvatom(*)
  real(dp),    intent(in)             :: d1(*)
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: ldoregions
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ia
  integer(i4)                         :: io
  integer(i4)                         :: is
  integer(i4)                         :: j
  integer(i4)                         :: k
  integer(i4)                         :: ka
  integer(i4)                         :: ko
  integer(i4)                         :: l
  integer(i4)                         :: la
  integer(i4)                         :: lo
  integer(i4)                         :: nj
  integer(i4)                         :: nk
  integer(i4)                         :: nl
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nregi
  integer(i4)                         :: nregk
  integer(i4)                         :: nregl
  integer(i4)                         :: ns
  real(dp)                            :: dr2ds(6)
  real(dp)                            :: d2r2dx2(3,3)
  real(dp)                            :: d2r2ds2(6,6)
  real(dp)                            :: d2r2dsdx(6,3)
  real(dp)                            :: xd
  real(dp)                            :: yd
  real(dp)                            :: zd
  real(dp)                            :: xdiff
  real(dp)                            :: ydiff
  real(dp)                            :: zdiff
#ifdef TRACE
  call trace_in('d1adds')
#endif
!
!  Loop over i -> neighbours
!
  ind = 0
  nri = nREBOatomRptr(i)
  do nk = 1,nneigh(nri)
    ind = ind + 1
    k = neighno(nk,nri)
    xd = xneigh(nk,nri)*d1(ind)
    yd = yneigh(nk,nri)*d1(ind)
    zd = zneigh(nk,nri)*d1(ind)
    io = nfreeatom(i)
    ko = nfreeatom(k)
    ia = nauatom(i)
    ka = nauatom(k)
    if (io.gt.0.and.ia.gt.0) then
      xdrv(ia) = xdrv(ia) - xd*dble(neqvatom(ia))
      ydrv(ia) = ydrv(ia) - yd*dble(neqvatom(ia))
      zdrv(ia) = zdrv(ia) - zd*dble(neqvatom(ia))
      if (ldoregions) then
        nregi = nregionno(ia)
        xregdrv(nregi) = xregdrv(nregi) - xd*dble(neqvatom(ia))
        yregdrv(nregi) = yregdrv(nregi) - yd*dble(neqvatom(ia))
        zregdrv(nregi) = zregdrv(nregi) - zd*dble(neqvatom(ia))
      endif
    endif
    if (ko.gt.0.and.ka.gt.0) then
      xdrv(ka) = xdrv(ka) + xd*dble(neqvatom(ka))
      ydrv(ka) = ydrv(ka) + yd*dble(neqvatom(ka))
      zdrv(ka) = zdrv(ka) + zd*dble(neqvatom(ka))
      if (ldoregions) then
        nregk = nregionno(ka)
        xregdrv(nregk) = xregdrv(nregk) + xd*dble(neqvatom(ka))
        yregdrv(nregk) = yregdrv(nregk) + yd*dble(neqvatom(ka))
        zregdrv(nregk) = zregdrv(nregk) + zd*dble(neqvatom(ka))
      endif
    endif
!
    if (lstr) then
      call real1strterm(ndim,xneigh(nk,nri),yneigh(nk,nri),zneigh(nk,nri),0.0_dp,0.0_dp,0.0_dp, &
        dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
      do is = 1,nstrains
        ns = nstrptr(is)
        rstrd(is) = rstrd(is) + d1(ind)*dr2ds(ns)
      enddo
    endif
  enddo
!
!  Loop over neighbours -> neighbours
!
  do nk = 2,nneigh(nri)
    k = neighno(nk,nri)
    do nl = 1,nk - 1
      ind = nneigh(nri) + nk*(nk-1)/2 + nl
      if (abs(d1(ind)).gt.d1tol) then
        l = neighno(nl,nri)
        xdiff = (xneigh(nk,nri) - xneigh(nl,nri))
        ydiff = (yneigh(nk,nri) - yneigh(nl,nri))
        zdiff = (zneigh(nk,nri) - zneigh(nl,nri))
        xd = xdiff*d1(ind)
        yd = ydiff*d1(ind)
        zd = zdiff*d1(ind)
        lo = nfreeatom(l)
        ko = nfreeatom(k)
        la = nauatom(l)
        ka = nauatom(k)
        if (lo.gt.0.and.la.gt.0) then
          xdrv(la) = xdrv(la) - xd*dble(neqvatom(la))
          ydrv(la) = ydrv(la) - yd*dble(neqvatom(la))
          zdrv(la) = zdrv(la) - zd*dble(neqvatom(la))
          if (ldoregions) then
            nregl = nregionno(la)
            xregdrv(nregl) = xregdrv(nregl) - xd*dble(neqvatom(la))
            yregdrv(nregl) = yregdrv(nregl) - yd*dble(neqvatom(la))
            zregdrv(nregl) = zregdrv(nregl) - zd*dble(neqvatom(la))
          endif
        endif
        if (ko.gt.0.and.ka.gt.0) then
          xdrv(ka) = xdrv(ka) + xd*dble(neqvatom(ka))
          ydrv(ka) = ydrv(ka) + yd*dble(neqvatom(ka))
          zdrv(ka) = zdrv(ka) + zd*dble(neqvatom(ka))
          if (ldoregions) then
            nregk = nregionno(ka)
            xregdrv(nregk) = xregdrv(nregk) + xd*dble(neqvatom(ka))
            yregdrv(nregk) = yregdrv(nregk) + yd*dble(neqvatom(ka))
            zregdrv(nregk) = zregdrv(nregk) + zd*dble(neqvatom(ka))
          endif
        endif
!
        if (lstr) then
          call real1strterm(ndim,xdiff,ydiff,zdiff,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          do is = 1,nstrains
            ns = nstrptr(is)
            rstrd(is) = rstrd(is) + d1(ind)*dr2ds(ns)
          enddo
        endif
      endif
    enddo
  enddo
!
!  Loop over i -> neighbours of neighbours
!
  if (ltorderv) then
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh
      do nl = 1,nneigh(nrj)
        ind = ind + 1
        if (abs(d1(ind)).gt.d1tol) then
          l = neighno(nl,nrj)
          xdiff = xneigh(nl,nrj) + xneigh(nj,nri)
          ydiff = yneigh(nl,nrj) + yneigh(nj,nri)
          zdiff = zneigh(nl,nrj) + zneigh(nj,nri)
          xd = xdiff*d1(ind)
          yd = ydiff*d1(ind)
          zd = zdiff*d1(ind)
          io = nfreeatom(i)
          lo = nfreeatom(l)
          ia = nauatom(i)
          la = nauatom(l)
          if (io.gt.0.and.ia.gt.0) then
            xdrv(ia) = xdrv(ia) - xd*dble(neqvatom(ia))
            ydrv(ia) = ydrv(ia) - yd*dble(neqvatom(ia))
            zdrv(ia) = zdrv(ia) - zd*dble(neqvatom(ia))
            if (ldoregions) then
              nregi = nregionno(ia)
              xregdrv(nregi) = xregdrv(nregi) - xd*dble(neqvatom(ia))
              yregdrv(nregi) = yregdrv(nregi) - yd*dble(neqvatom(ia))
              zregdrv(nregi) = zregdrv(nregi) - zd*dble(neqvatom(ia))
            endif
          endif
          if (lo.gt.0.and.la.gt.0) then
            xdrv(la) = xdrv(la) + xd*dble(neqvatom(la))
            ydrv(la) = ydrv(la) + yd*dble(neqvatom(la))
            zdrv(la) = zdrv(la) + zd*dble(neqvatom(la))
            if (ldoregions) then
              nregl = nregionno(la)
              xregdrv(nregl) = xregdrv(nregl) + xd*dble(neqvatom(la))
              yregdrv(nregl) = yregdrv(nregl) + yd*dble(neqvatom(la))
              zregdrv(nregl) = zregdrv(nregl) + zd*dble(neqvatom(la))
            endif
          endif
!
          if (lstr) then
            call real1strterm(ndim,xdiff,ydiff,zdiff,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
            do is = 1,nstrains
              ns = nstrptr(is)
              rstrd(is) = rstrd(is) + d1(ind)*dr2ds(ns)
            enddo
          endif
        endif
      enddo
    enddo
!
!  Loop over neighbours -> neighbours of neighbours
!
    do nj = 1,nneigh(nri)
      j = neighno(nj,nri)
      nrj = nREBOatomRptr(j)
      do nk = 1,nneigh(nri)
        k = neighno(nk,nri)
        ind = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh + nk*mneigh
        do nl = 1,nneigh(nrj)
          ind = ind + 1
          if (abs(d1(ind)).gt.d1tol) then
            l = neighno(nl,nrj)
            xdiff = xneigh(nl,nrj) + xneigh(nj,nri) - xneigh(nk,nri)
            ydiff = yneigh(nl,nrj) + yneigh(nj,nri) - yneigh(nk,nri)
            zdiff = zneigh(nl,nrj) + zneigh(nj,nri) - zneigh(nk,nri)
            xd = xdiff*d1(ind)
            yd = ydiff*d1(ind)
            zd = zdiff*d1(ind)
            ko = nfreeatom(k)
            lo = nfreeatom(l)
            ka = nauatom(k)
            la = nauatom(l)
            if (ko.gt.0.and.ka.gt.0) then
              xdrv(ka) = xdrv(ka) - xd*dble(neqvatom(ka))
              ydrv(ka) = ydrv(ka) - yd*dble(neqvatom(ka))
              zdrv(ka) = zdrv(ka) - zd*dble(neqvatom(ka))
              if (ldoregions) then
                nregk = nregionno(ka)
                xregdrv(nregk) = xregdrv(nregk) - xd*dble(neqvatom(ka))
                yregdrv(nregk) = yregdrv(nregk) - yd*dble(neqvatom(ka))
                zregdrv(nregk) = zregdrv(nregk) - zd*dble(neqvatom(ka))
              endif
            endif
            if (lo.gt.0.and.la.gt.0) then
              xdrv(la) = xdrv(la) + xd*dble(neqvatom(la))
              ydrv(la) = ydrv(la) + yd*dble(neqvatom(la))
              zdrv(la) = zdrv(la) + zd*dble(neqvatom(la))
              if (ldoregions) then
                nregl = nregionno(la)
                xregdrv(nregl) = xregdrv(nregl) + xd*dble(neqvatom(la))
                yregdrv(nregl) = yregdrv(nregl) + yd*dble(neqvatom(la))
                zregdrv(nregl) = zregdrv(nregl) + zd*dble(neqvatom(la))
              endif
            endif
!
            if (lstr) then
              call real1strterm(ndim,xdiff,ydiff,zdiff,0.0_dp,0.0_dp,0.0_dp,dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
              do is = 1,nstrains
                ns = nstrptr(is)
                rstrd(is) = rstrd(is) + d1(ind)*dr2ds(ns)
              enddo
            endif
          endif
        enddo
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('d1adds')
#endif
!
  return
  end
!
  subroutine d2add(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1,d2,ltorderv,lextended)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!
!   6/02 Created from d1add
!   6/02 Modified to allow for torsional terms
!   7/02 Indexing further modified for torsions
!   8/02 Strain derivatives added
!   8/02 Freezing added
!   8/02 Removal of uninvalid terms added
!   9/02 Correction to lsameijkl added - need to check distances too
!  12/03 Torsional derivatives made optional
!   6/07 nREBOatomRptr added as input
!  10/07 Correction made to setting of n1i/n2k for torsional case
!  11/07 Second criterion for lsameijkl removed as it gives wrong second
!        derivatives for reaxff.
!  11/07 Unused variables removed
!   1/10 Bug fixed for neighbour of i -> neighbour of neighbour of i case
!   2/18 Trace added
!   9/18 Strain module added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!  12/22 Extended derivatives added
!   2/23 ind now set explicitly before n2 loop in case there is a cycle during extended derivatives
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use current,        only : ndim, nstrains
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: is1,is2
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: ixc,iyc,izc
  integer(i4)                         :: jxc,jyc,jzc
  integer(i4)                         :: kxc,kyc,kzc
  integer(i4)                         :: lxc,lyc,lzc
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nji
  integer(i4)                         :: nkj
  integer(i4)                         :: nlk
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: ns1
  integer(i4)                         :: ns2
  integer(i4)                         :: ntmp
  real(dp)                            :: one1
  real(dp)                            :: one2
  real(dp)                            :: dr2ds1(6)
  real(dp)                            :: dr2ds2(6)
  real(dp)                            :: d2r2dx21(3,3)
  real(dp)                            :: d2r2dx22(3,3)
  real(dp)                            :: d2r2ds21(6,6)
  real(dp)                            :: d2r2ds22(6,6)
  real(dp)                            :: d2r2dsdx1(6,3)
  real(dp)                            :: d2r2dsdx2(6,3)
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: lsameijkl
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lopi
  logical                             :: lopj
  logical                             :: lopk
  logical                             :: lopl
#ifdef TRACE
  call trace_in('d2add')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndisttor = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
    if (lextended) then
      ndist = ndisttor + nneigh(nri)*mneigh*mneigh
    else
      ndist = ndisttor
    endif
  else
    ndist = nneigh2 
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = nneigh2 + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
!  Loop over pairs of interatomic distances
!
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Check that neighbour numbers are valid
!
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    elseif (n1.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
!  Check that atoms are valid
!
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        if (nREBOatomRptr(neighno(nji,nri)).eq.0) then
          ljok = .false.
        else
          ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
        endif
      else
        ljok = .false.
      endif
!
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!         
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    else
!
!  Extended derivatives
!
      ntmp = n1 - ndisttor
      nji = (ntmp - 1)/(mneigh*mneigh) + 1
      ntmp = ntmp - (nji - 1)*mneigh*mneigh
      nkj = (ntmp - 1)/mneigh + 1
      ntmp = ntmp - (nkj - 1)*mneigh
      nlk = ntmp
!
      nrj = neighno(nji,nri)
      if (nrj.eq.0) cycle
      if (nkj.gt.nneigh(nrj)) cycle
      nrk = neighno(nkj,nrj)
      if (nrk.eq.0) cycle
      if (nlk.gt.nneigh(nrk)) cycle
      n1i = nrk
      n1j = neighno(nlk,nrk)
      if (n1j.eq.0) cycle
!   
!  Check that atoms are valid
!  
      liok = (nkj.le.nneigh(nrj))
      ljok = (nlk.le.nneigh(nrk))
!
      if (liok.and.ljok) then
        x2i = xneigh(nji,nri) + xneigh(nkj,nrj)
        y2i = yneigh(nji,nri) + yneigh(nkj,nrj)
        z2i = zneigh(nji,nri) + zneigh(nkj,nrj)
!
        xd1 = xneigh(nlk,nrk)
        yd1 = yneigh(nlk,nrk)
        zd1 = zneigh(nlk,nrk)
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
!
!  Set up products for strain
!
      if (lstr) then
        call real1strterm(ndim,xd1,yd1,zd1,0.0_dp,0.0_dp,0.0_dp,dr2ds1,d2r2dx21,d2r2dsdx1,d2r2ds21,.true.)
      endif
!
!  Set second derivative location pointers
!
      lopi = (nfreeatom(n1i).gt.0)
      lopj = (nfreeatom(n1j).gt.0)
      ix = 3*(nfreeatom(n1i) - 1) + 1
      iy = ix + 1
      iz = iy + 1
      jx = 3*(nfreeatom(n1j) - 1) + 1
      jy = jx + 1
      jz = jy + 1
    endif
!
!  If neither i nor j are being optimised and strain derivatives are not needed then atoms are not valid
!
    if (.not.lopi.and..not.lopj.and..not.lstr) then
      liok = .false.
      ljok = .false.
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
      ind = n1*(n1-1)/2
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!  
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          elseif (n2.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              if (nREBOatomRptr(neighno(nji,nri)).eq.0) then
                llok = .false.
              else
                llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
              endif
            else
              llok = .false.
            endif
!
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!               
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          else
! 
!  Extended derivatives               
! 
            ntmp = n2 - ndisttor
            nji = (ntmp - 1)/(mneigh*mneigh) + 1
            ntmp = ntmp - (nji - 1)*mneigh*mneigh
            nkj = (ntmp - 1)/mneigh + 1
            ntmp = ntmp - (nkj - 1)*mneigh
            nlk = ntmp
! 
            nrj = neighno(nji,nri)
            if (nrj.eq.0) cycle
            if (nkj.gt.nneigh(nrj)) cycle
            nrk = neighno(nkj,nrj)
            if (nrk.eq.0) cycle
            if (nlk.gt.nneigh(nrk)) cycle
            n2k = nrk
            n2l = neighno(nlk,nrk)
            if (n2l.eq.0) cycle
!   
!  Check that atoms are valid         
!  
            lkok = (nkj.le.nneigh(nrj))
            llok = (nlk.le.nneigh(nrk))
! 
            if (liok.and.ljok) then         
              xik = xneigh(nji,nri) + xneigh(nkj,nrj) - x2i
              yik = yneigh(nji,nri) + yneigh(nkj,nrj) - y2i
              zik = zneigh(nji,nri) + zneigh(nkj,nrj) - z2i
! 
              xd2 = xneigh(nlk,nrk)
              yd2 = yneigh(nlk,nrk)
              zd2 = zneigh(nlk,nrk)
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
!  Set up products for strain
!
            if (lstr) then
              call real1strterm(ndim,xd2,yd2,zd2,0.0_dp,0.0_dp,0.0_dp,dr2ds2,d2r2dx22,d2r2dsdx2,d2r2ds22,.true.)
            endif
!
            lopk = (nfreeatom(n2k).gt.0)
            lopl = (nfreeatom(n2l).gt.0)
            kx = 3*(nfreeatom(n2k) - 1) + 1
            ky = kx + 1
            kz = ky + 1
            lx = 3*(nfreeatom(n2l) - 1) + 1
            ly = lx + 1
            lz = ly + 1
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if ((lopi.or.lopk).and.n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              if (lopi.and.lopk) then
                ixc = ix
                iyc = iy
                izc = iz
                kxc = kx
                kyc = ky
                kzc = kz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ixc = ix
                iyc = iy
                izc = iz
                kxc = ix
                kyc = iy
                kzc = iz
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                ixc = kx
                iyc = ky
                izc = kz
                kxc = kx
                kyc = ky
                kzc = kz
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(kxc,ixc) = derv2(kxc,ixc) + xd1*xd2*d2(ind)*one1
              derv2(kyc,ixc) = derv2(kyc,ixc) + xd1*yd2*d2(ind)*one1
              derv2(kzc,ixc) = derv2(kzc,ixc) + xd1*zd2*d2(ind)*one1
              derv2(kxc,iyc) = derv2(kxc,iyc) + yd1*xd2*d2(ind)*one1
              derv2(kyc,iyc) = derv2(kyc,iyc) + yd1*yd2*d2(ind)*one1
              derv2(kzc,iyc) = derv2(kzc,iyc) + yd1*zd2*d2(ind)*one1
              derv2(kxc,izc) = derv2(kxc,izc) + zd1*xd2*d2(ind)*one1
              derv2(kyc,izc) = derv2(kyc,izc) + zd1*yd2*d2(ind)*one1
              derv2(kzc,izc) = derv2(kzc,izc) + zd1*zd2*d2(ind)*one1
              if (n1.eq.n2) then
                derv2(kxc,ixc) = derv2(kxc,ixc) + d1(n2)*one1
                derv2(kyc,iyc) = derv2(kyc,iyc) + d1(n2)*one1
                derv2(kzc,izc) = derv2(kzc,izc) + d1(n2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(ixc,kxc) = derv2(ixc,kxc) + xd1*xd2*d2(ind)*one2
                derv2(ixc,kyc) = derv2(ixc,kyc) + xd1*yd2*d2(ind)*one2
                derv2(ixc,kzc) = derv2(ixc,kzc) + xd1*zd2*d2(ind)*one2
                derv2(iyc,kxc) = derv2(iyc,kxc) + yd1*xd2*d2(ind)*one2
                derv2(iyc,kyc) = derv2(iyc,kyc) + yd1*yd2*d2(ind)*one2
                derv2(iyc,kzc) = derv2(iyc,kzc) + yd1*zd2*d2(ind)*one2
                derv2(izc,kxc) = derv2(izc,kxc) + zd1*xd2*d2(ind)*one2
                derv2(izc,kyc) = derv2(izc,kyc) + zd1*yd2*d2(ind)*one2
                derv2(izc,kzc) = derv2(izc,kzc) + zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(ixc,kxc) = derv2(ixc,kxc) + d1(n2)*one2
                  derv2(iyc,kyc) = derv2(iyc,kyc) + d1(n2)*one2
                  derv2(izc,kzc) = derv2(izc,kzc) + d1(n2)*one2
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if ((lopi.or.lopl).and.n1i.ne.n2l) then
              if (lopi.and.lopl) then
                ixc = ix
                iyc = iy
                izc = iz
                lxc = lx
                lyc = ly
                lzc = lz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ixc = ix
                iyc = iy
                izc = iz
                lxc = ix
                lyc = iy
                lzc = iz
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                ixc = lx
                iyc = ly
                izc = lz
                lxc = lx
                lyc = ly
                lzc = lz
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(lxc,ixc) = derv2(lxc,ixc) - xd1*xd2*d2(ind)*one1
              derv2(lyc,ixc) = derv2(lyc,ixc) - xd1*yd2*d2(ind)*one1
              derv2(lzc,ixc) = derv2(lzc,ixc) - xd1*zd2*d2(ind)*one1
              derv2(lxc,iyc) = derv2(lxc,iyc) - yd1*xd2*d2(ind)*one1
              derv2(lyc,iyc) = derv2(lyc,iyc) - yd1*yd2*d2(ind)*one1
              derv2(lzc,iyc) = derv2(lzc,iyc) - yd1*zd2*d2(ind)*one1
              derv2(lxc,izc) = derv2(lxc,izc) - zd1*xd2*d2(ind)*one1
              derv2(lyc,izc) = derv2(lyc,izc) - zd1*yd2*d2(ind)*one1
              derv2(lzc,izc) = derv2(lzc,izc) - zd1*zd2*d2(ind)*one1
              if (n1.eq.n2) then
                derv2(lxc,ixc) = derv2(lxc,ixc) - d1(n2)*one1
                derv2(lyc,iyc) = derv2(lyc,iyc) - d1(n2)*one1
                derv2(lzc,izc) = derv2(lzc,izc) - d1(n2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(ixc,lxc) = derv2(ixc,lxc) - xd1*xd2*d2(ind)*one2
                derv2(ixc,lyc) = derv2(ixc,lyc) - xd1*yd2*d2(ind)*one2
                derv2(ixc,lzc) = derv2(ixc,lzc) - xd1*zd2*d2(ind)*one2
                derv2(iyc,lxc) = derv2(iyc,lxc) - yd1*xd2*d2(ind)*one2
                derv2(iyc,lyc) = derv2(iyc,lyc) - yd1*yd2*d2(ind)*one2
                derv2(iyc,lzc) = derv2(iyc,lzc) - yd1*zd2*d2(ind)*one2
                derv2(izc,lxc) = derv2(izc,lxc) - zd1*xd2*d2(ind)*one2
                derv2(izc,lyc) = derv2(izc,lyc) - zd1*yd2*d2(ind)*one2
                derv2(izc,lzc) = derv2(izc,lzc) - zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(ixc,lxc) = derv2(ixc,lxc) - d1(n2)*one2
                  derv2(iyc,lyc) = derv2(iyc,lyc) - d1(n2)*one2
                  derv2(izc,lzc) = derv2(izc,lzc) - d1(n2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if ((lopj.or.lopk).and.n1j.ne.n2k) then
              if (lopj.and.lopk) then
                jxc = jx
                jyc = jy
                jzc = jz
                kxc = kx
                kyc = ky
                kzc = kz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                jxc = jx
                jyc = jy
                jzc = jz
                kxc = jx
                kyc = jy
                kzc = jz
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                jxc = kx
                jyc = ky
                jzc = kz
                kxc = kx
                kyc = ky
                kzc = kz
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(kxc,jxc) = derv2(kxc,jxc) - xd1*xd2*d2(ind)*one1
              derv2(kyc,jxc) = derv2(kyc,jxc) - xd1*yd2*d2(ind)*one1
              derv2(kzc,jxc) = derv2(kzc,jxc) - xd1*zd2*d2(ind)*one1
              derv2(kxc,jyc) = derv2(kxc,jyc) - yd1*xd2*d2(ind)*one1
              derv2(kyc,jyc) = derv2(kyc,jyc) - yd1*yd2*d2(ind)*one1
              derv2(kzc,jyc) = derv2(kzc,jyc) - yd1*zd2*d2(ind)*one1
              derv2(kxc,jzc) = derv2(kxc,jzc) - zd1*xd2*d2(ind)*one1
              derv2(kyc,jzc) = derv2(kyc,jzc) - zd1*yd2*d2(ind)*one1
              derv2(kzc,jzc) = derv2(kzc,jzc) - zd1*zd2*d2(ind)*one1
              if (n1.eq.n2) then
                derv2(kxc,jxc) = derv2(kxc,jxc) - d1(n2)*one1
                derv2(kyc,jyc) = derv2(kyc,jyc) - d1(n2)*one1
                derv2(kzc,jzc) = derv2(kzc,jzc) - d1(n2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(jxc,kxc) = derv2(jxc,kxc) - xd1*xd2*d2(ind)*one2
                derv2(jxc,kyc) = derv2(jxc,kyc) - xd1*yd2*d2(ind)*one2
                derv2(jxc,kzc) = derv2(jxc,kzc) - xd1*zd2*d2(ind)*one2
                derv2(jyc,kxc) = derv2(jyc,kxc) - yd1*xd2*d2(ind)*one2
                derv2(jyc,kyc) = derv2(jyc,kyc) - yd1*yd2*d2(ind)*one2
                derv2(jyc,kzc) = derv2(jyc,kzc) - yd1*zd2*d2(ind)*one2
                derv2(jzc,kxc) = derv2(jzc,kxc) - zd1*xd2*d2(ind)*one2
                derv2(jzc,kyc) = derv2(jzc,kyc) - zd1*yd2*d2(ind)*one2
                derv2(jzc,kzc) = derv2(jzc,kzc) - zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(jxc,kxc) = derv2(jxc,kxc) - d1(n2)*one2
                  derv2(jyc,kyc) = derv2(jyc,kyc) - d1(n2)*one2
                  derv2(jzc,kzc) = derv2(jzc,kzc) - d1(n2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if ((lopj.or.lopl).and.n1j.ne.n2l) then
              if (lopj.and.lopl) then
                jxc = jx
                jyc = jy
                jzc = jz
                lxc = lx
                lyc = ly
                lzc = lz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                jxc = jx
                jyc = jy
                jzc = jz
                lxc = jx
                lyc = jy
                lzc = jz
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                jxc = lx
                jyc = ly
                jzc = lz
                lxc = lx
                lyc = ly
                lzc = lz
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(lxc,jxc) = derv2(lxc,jxc) + xd1*xd2*d2(ind)*one1
              derv2(lyc,jxc) = derv2(lyc,jxc) + xd1*yd2*d2(ind)*one1
              derv2(lzc,jxc) = derv2(lzc,jxc) + xd1*zd2*d2(ind)*one1
              derv2(lxc,jyc) = derv2(lxc,jyc) + yd1*xd2*d2(ind)*one1
              derv2(lyc,jyc) = derv2(lyc,jyc) + yd1*yd2*d2(ind)*one1
              derv2(lzc,jyc) = derv2(lzc,jyc) + yd1*zd2*d2(ind)*one1
              derv2(lxc,jzc) = derv2(lxc,jzc) + zd1*xd2*d2(ind)*one1
              derv2(lyc,jzc) = derv2(lyc,jzc) + zd1*yd2*d2(ind)*one1
              derv2(lzc,jzc) = derv2(lzc,jzc) + zd1*zd2*d2(ind)*one1
              if (n1.eq.n2) then
                derv2(lxc,jxc) = derv2(lxc,jxc) + d1(n2)*one1
                derv2(lyc,jyc) = derv2(lyc,jyc) + d1(n2)*one1
                derv2(lzc,jzc) = derv2(lzc,jzc) + d1(n2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(jxc,lxc) = derv2(jxc,lxc) + xd1*xd2*d2(ind)*one2
                derv2(jxc,lyc) = derv2(jxc,lyc) + xd1*yd2*d2(ind)*one2
                derv2(jxc,lzc) = derv2(jxc,lzc) + xd1*zd2*d2(ind)*one2
                derv2(jyc,lxc) = derv2(jyc,lxc) + yd1*xd2*d2(ind)*one2
                derv2(jyc,lyc) = derv2(jyc,lyc) + yd1*yd2*d2(ind)*one2
                derv2(jyc,lzc) = derv2(jyc,lzc) + yd1*zd2*d2(ind)*one2
                derv2(jzc,lxc) = derv2(jzc,lxc) + zd1*xd2*d2(ind)*one2
                derv2(jzc,lyc) = derv2(jzc,lyc) + zd1*yd2*d2(ind)*one2
                derv2(jzc,lzc) = derv2(jzc,lzc) + zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(jxc,lxc) = derv2(jxc,lxc) + d1(n2)*one2
                  derv2(jyc,lyc) = derv2(jyc,lyc) + d1(n2)*one2
                  derv2(jzc,lzc) = derv2(jzc,lzc) + d1(n2)*one2
                endif
              endif
            endif
            if (lstr) then
!
!  Strain - strain second derivatives
!
              do is1 = 1,nstrains
                ns1 = nstrptr(is1)
                do is2 = 1,nstrains
                  ns2 = nstrptr(is2)
                  sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds2(ns2)*dr2ds1(ns1)
                  if (n1.ne.n2) then
                    sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds1(ns2)*dr2ds2(ns1)
                  else
                    sderv2(is2,is1) = sderv2(is2,is1) + d1(n2)*d2r2ds22(ns2,ns1)
                  endif
                enddo
              enddo
!
!  Internal - strain second derivatives
!
              do is1 = 1,nstrains
                ns1 = nstrptr(is1)
                derv3(ix,is1) = derv3(ix,is1) - xd1*d2(ind)*dr2ds2(ns1)
                derv3(iy,is1) = derv3(iy,is1) - yd1*d2(ind)*dr2ds2(ns1)
                derv3(iz,is1) = derv3(iz,is1) - zd1*d2(ind)*dr2ds2(ns1)
                derv3(jx,is1) = derv3(jx,is1) + xd1*d2(ind)*dr2ds2(ns1)
                derv3(jy,is1) = derv3(jy,is1) + yd1*d2(ind)*dr2ds2(ns1)
                derv3(jz,is1) = derv3(jz,is1) + zd1*d2(ind)*dr2ds2(ns1)
                if (n1.ne.n2) then
                  derv3(kx,is1) = derv3(kx,is1) - xd2*d2(ind)*dr2ds1(ns1)
                  derv3(ky,is1) = derv3(ky,is1) - yd2*d2(ind)*dr2ds1(ns1)
                  derv3(kz,is1) = derv3(kz,is1) - zd2*d2(ind)*dr2ds1(ns1)
                  derv3(lx,is1) = derv3(lx,is1) + xd2*d2(ind)*dr2ds1(ns1)
                  derv3(ly,is1) = derv3(ly,is1) + yd2*d2(ind)*dr2ds1(ns1)
                  derv3(lz,is1) = derv3(lz,is1) + zd2*d2(ind)*dr2ds1(ns1)
                else
                  derv3(ix,is1) = derv3(ix,is1) - d1(n2)*d2r2dsdx2(ns1,1)
                  derv3(iy,is1) = derv3(iy,is1) - d1(n2)*d2r2dsdx2(ns1,2)
                  derv3(iz,is1) = derv3(iz,is1) - d1(n2)*d2r2dsdx2(ns1,3)
                  derv3(jx,is1) = derv3(jx,is1) + d1(n2)*d2r2dsdx2(ns1,1)
                  derv3(jy,is1) = derv3(jy,is1) + d1(n2)*d2r2dsdx2(ns1,2)
                  derv3(jz,is1) = derv3(jz,is1) + d1(n2)*d2r2dsdx2(ns1,3)
                endif
              enddo
            endif
          endif
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('d2add')
#endif
!
  return
  end
!
  subroutine d2addnn(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1nn,d2nn)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array.
!  Nearest-neighbour only variant of d2add.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1nn            = array of first derivatives w.r.t. interatomic distances
!  d2nn            = array of first derivatives w.r.t. interatomic distances
!
!   8/22 Created from d2add
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use current,        only : ndim, nstrains
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1nn(*)
  real(dp),    intent(in)             :: d2nn(*)
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: is1,is2
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: ixc,iyc,izc
  integer(i4)                         :: jxc,jyc,jzc
  integer(i4)                         :: kxc,kyc,kzc
  integer(i4)                         :: lxc,lyc,lzc
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: ndist
  integer(i4)                         :: nn1i
  integer(i4)                         :: nn1j
  integer(i4)                         :: nn2k
  integer(i4)                         :: nn2l
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: ns1
  integer(i4)                         :: ns2
  real(dp)                            :: one1
  real(dp)                            :: one2
  real(dp)                            :: dr2ds1(6)
  real(dp)                            :: dr2ds2(6)
  real(dp)                            :: d2r2dx21(3,3)
  real(dp)                            :: d2r2dx22(3,3)
  real(dp)                            :: d2r2ds21(6,6)
  real(dp)                            :: d2r2ds22(6,6)
  real(dp)                            :: d2r2dsdx1(6,3)
  real(dp)                            :: d2r2dsdx2(6,3)
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: lsameijkl
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lopi
  logical                             :: lopj
  logical                             :: lopk
  logical                             :: lopl
#ifdef TRACE
  call trace_in('d2addnn')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*mneigh
  ndist = nneigh2 
!
!  Loop over pairs of interatomic distances
!
  ind = 0
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    else
!
!  Neighbours of i -> neighbours of neighbours of i
!
!  Find indices and check that neighbour numbers are valid
!
      nn1i = (n1-nneigh(nri)-1)/mneigh + 1
      liok = (nn1i.le.nneigh(nri))
      if (liok) then
        n1i = neighno(nn1i,nri)
        nn1j = n1 - nneigh(nri) - (nn1i-1)*mneigh
        ljok = (nn1j.le.nneigh(n1i))
        if (ljok) then
          n1j = neighno(nn1j,n1i)
        endif
      else
        ljok = .false.
      endif
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = - xneigh(nn1j,n1i)
        yd1 = - yneigh(nn1j,n1i)
        zd1 = - zneigh(nn1j,n1i)
!  
        x2i = xneigh(nn1i,nri)
        y2i = yneigh(nn1i,nri)
        z2i = zneigh(nn1i,nri)
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
!
!  Set up products for strain
!
      if (lstr) then
        call real1strterm(ndim,xd1,yd1,zd1,0.0_dp,0.0_dp,0.0_dp,dr2ds1,d2r2dx21,d2r2dsdx1,d2r2ds21,.true.)
      endif
!
!  Set second derivative location pointers
!
      lopi = (nfreeatom(n1i).gt.0)
      lopj = (nfreeatom(n1j).gt.0)
      ix = 3*(nfreeatom(n1i) - 1) + 1
      iy = ix + 1
      iz = iy + 1
      jx = 3*(nfreeatom(n1j) - 1) + 1
      jy = jx + 1
      jz = jy + 1
    endif
!
!  If neither i nor j are being optimised and strain derivatives are not needed then atoms are not valid
!
    if (.not.lopi.and..not.lopj.and..not.lstr) then
      liok = .false.
      ljok = .false.
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2nn(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1nn(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          else
!
!  Neighbours of i -> neighbours of neighbours of i
!
!  Find indices and check that neighbours are valid
!
            nn2k = (n2-nneigh(nri)-1)/mneigh + 1
            lkok = (nn2k.le.nneigh(nri))
            if (lkok) then
              n2k = neighno(nn2k,nri)
              nn2l = n2 - nneigh(nri) - (nn2k-1)*mneigh
              llok = (nn2l.le.nneigh(n2k))
              if (llok) then
                n2l = neighno(nn2l,n2k)
              endif
            else
              llok = .false.
            endif
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = - xneigh(nn2l,n2k)
              yd2 = - yneigh(nn2l,n2k)
              zd2 = - zneigh(nn2l,n2k)
!  
              xik = xneigh(nn2k,nri) - x2i
              yik = yneigh(nn2k,nri) - y2i
              zik = zneigh(nn2k,nri) - z2i
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
!  Set up products for strain
!
            if (lstr) then
              call real1strterm(ndim,xd2,yd2,zd2,0.0_dp,0.0_dp,0.0_dp,dr2ds2,d2r2dx22,d2r2dsdx2,d2r2ds22,.true.)
            endif
!
            lopk = (nfreeatom(n2k).gt.0)
            lopl = (nfreeatom(n2l).gt.0)
            kx = 3*(nfreeatom(n2k) - 1) + 1
            ky = kx + 1
            kz = ky + 1
            lx = 3*(nfreeatom(n2l) - 1) + 1
            ly = lx + 1
            lz = ly + 1
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if ((lopi.or.lopk).and.n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              if (lopi.and.lopk) then
                ixc = ix
                iyc = iy
                izc = iz
                kxc = kx
                kyc = ky
                kzc = kz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ixc = ix
                iyc = iy
                izc = iz
                kxc = ix
                kyc = iy
                kzc = iz
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                ixc = kx
                iyc = ky
                izc = kz
                kxc = kx
                kyc = ky
                kzc = kz
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(kxc,ixc) = derv2(kxc,ixc) + xd1*xd2*d2nn(ind)*one1
              derv2(kyc,ixc) = derv2(kyc,ixc) + xd1*yd2*d2nn(ind)*one1
              derv2(kzc,ixc) = derv2(kzc,ixc) + xd1*zd2*d2nn(ind)*one1
              derv2(kxc,iyc) = derv2(kxc,iyc) + yd1*xd2*d2nn(ind)*one1
              derv2(kyc,iyc) = derv2(kyc,iyc) + yd1*yd2*d2nn(ind)*one1
              derv2(kzc,iyc) = derv2(kzc,iyc) + yd1*zd2*d2nn(ind)*one1
              derv2(kxc,izc) = derv2(kxc,izc) + zd1*xd2*d2nn(ind)*one1
              derv2(kyc,izc) = derv2(kyc,izc) + zd1*yd2*d2nn(ind)*one1
              derv2(kzc,izc) = derv2(kzc,izc) + zd1*zd2*d2nn(ind)*one1
              if (n1.eq.n2) then
                derv2(kxc,ixc) = derv2(kxc,ixc) + d1nn(n2)*one1
                derv2(kyc,iyc) = derv2(kyc,iyc) + d1nn(n2)*one1
                derv2(kzc,izc) = derv2(kzc,izc) + d1nn(n2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(ixc,kxc) = derv2(ixc,kxc) + xd1*xd2*d2nn(ind)*one2
                derv2(ixc,kyc) = derv2(ixc,kyc) + xd1*yd2*d2nn(ind)*one2
                derv2(ixc,kzc) = derv2(ixc,kzc) + xd1*zd2*d2nn(ind)*one2
                derv2(iyc,kxc) = derv2(iyc,kxc) + yd1*xd2*d2nn(ind)*one2
                derv2(iyc,kyc) = derv2(iyc,kyc) + yd1*yd2*d2nn(ind)*one2
                derv2(iyc,kzc) = derv2(iyc,kzc) + yd1*zd2*d2nn(ind)*one2
                derv2(izc,kxc) = derv2(izc,kxc) + zd1*xd2*d2nn(ind)*one2
                derv2(izc,kyc) = derv2(izc,kyc) + zd1*yd2*d2nn(ind)*one2
                derv2(izc,kzc) = derv2(izc,kzc) + zd1*zd2*d2nn(ind)*one2
                if (n1.eq.n2) then
                  derv2(ixc,kxc) = derv2(ixc,kxc) + d1nn(n2)*one2
                  derv2(iyc,kyc) = derv2(iyc,kyc) + d1nn(n2)*one2
                  derv2(izc,kzc) = derv2(izc,kzc) + d1nn(n2)*one2
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if ((lopi.or.lopl).and.n1i.ne.n2l) then
              if (lopi.and.lopl) then
                ixc = ix
                iyc = iy
                izc = iz
                lxc = lx
                lyc = ly
                lzc = lz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ixc = ix
                iyc = iy
                izc = iz
                lxc = ix
                lyc = iy
                lzc = iz
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                ixc = lx
                iyc = ly
                izc = lz
                lxc = lx
                lyc = ly
                lzc = lz
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(lxc,ixc) = derv2(lxc,ixc) - xd1*xd2*d2nn(ind)*one1
              derv2(lyc,ixc) = derv2(lyc,ixc) - xd1*yd2*d2nn(ind)*one1
              derv2(lzc,ixc) = derv2(lzc,ixc) - xd1*zd2*d2nn(ind)*one1
              derv2(lxc,iyc) = derv2(lxc,iyc) - yd1*xd2*d2nn(ind)*one1
              derv2(lyc,iyc) = derv2(lyc,iyc) - yd1*yd2*d2nn(ind)*one1
              derv2(lzc,iyc) = derv2(lzc,iyc) - yd1*zd2*d2nn(ind)*one1
              derv2(lxc,izc) = derv2(lxc,izc) - zd1*xd2*d2nn(ind)*one1
              derv2(lyc,izc) = derv2(lyc,izc) - zd1*yd2*d2nn(ind)*one1
              derv2(lzc,izc) = derv2(lzc,izc) - zd1*zd2*d2nn(ind)*one1
              if (n1.eq.n2) then
                derv2(lxc,ixc) = derv2(lxc,ixc) - d1nn(n2)*one1
                derv2(lyc,iyc) = derv2(lyc,iyc) - d1nn(n2)*one1
                derv2(lzc,izc) = derv2(lzc,izc) - d1nn(n2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(ixc,lxc) = derv2(ixc,lxc) - xd1*xd2*d2nn(ind)*one2
                derv2(ixc,lyc) = derv2(ixc,lyc) - xd1*yd2*d2nn(ind)*one2
                derv2(ixc,lzc) = derv2(ixc,lzc) - xd1*zd2*d2nn(ind)*one2
                derv2(iyc,lxc) = derv2(iyc,lxc) - yd1*xd2*d2nn(ind)*one2
                derv2(iyc,lyc) = derv2(iyc,lyc) - yd1*yd2*d2nn(ind)*one2
                derv2(iyc,lzc) = derv2(iyc,lzc) - yd1*zd2*d2nn(ind)*one2
                derv2(izc,lxc) = derv2(izc,lxc) - zd1*xd2*d2nn(ind)*one2
                derv2(izc,lyc) = derv2(izc,lyc) - zd1*yd2*d2nn(ind)*one2
                derv2(izc,lzc) = derv2(izc,lzc) - zd1*zd2*d2nn(ind)*one2
                if (n1.eq.n2) then
                  derv2(ixc,lxc) = derv2(ixc,lxc) - d1nn(n2)*one2
                  derv2(iyc,lyc) = derv2(iyc,lyc) - d1nn(n2)*one2
                  derv2(izc,lzc) = derv2(izc,lzc) - d1nn(n2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if ((lopj.or.lopk).and.n1j.ne.n2k) then
              if (lopj.and.lopk) then
                jxc = jx
                jyc = jy
                jzc = jz
                kxc = kx
                kyc = ky
                kzc = kz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                jxc = jx
                jyc = jy
                jzc = jz
                kxc = jx
                kyc = jy
                kzc = jz
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                jxc = kx
                jyc = ky
                jzc = kz
                kxc = kx
                kyc = ky
                kzc = kz
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(kxc,jxc) = derv2(kxc,jxc) - xd1*xd2*d2nn(ind)*one1
              derv2(kyc,jxc) = derv2(kyc,jxc) - xd1*yd2*d2nn(ind)*one1
              derv2(kzc,jxc) = derv2(kzc,jxc) - xd1*zd2*d2nn(ind)*one1
              derv2(kxc,jyc) = derv2(kxc,jyc) - yd1*xd2*d2nn(ind)*one1
              derv2(kyc,jyc) = derv2(kyc,jyc) - yd1*yd2*d2nn(ind)*one1
              derv2(kzc,jyc) = derv2(kzc,jyc) - yd1*zd2*d2nn(ind)*one1
              derv2(kxc,jzc) = derv2(kxc,jzc) - zd1*xd2*d2nn(ind)*one1
              derv2(kyc,jzc) = derv2(kyc,jzc) - zd1*yd2*d2nn(ind)*one1
              derv2(kzc,jzc) = derv2(kzc,jzc) - zd1*zd2*d2nn(ind)*one1
              if (n1.eq.n2) then
                derv2(kxc,jxc) = derv2(kxc,jxc) - d1nn(n2)*one1
                derv2(kyc,jyc) = derv2(kyc,jyc) - d1nn(n2)*one1
                derv2(kzc,jzc) = derv2(kzc,jzc) - d1nn(n2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(jxc,kxc) = derv2(jxc,kxc) - xd1*xd2*d2nn(ind)*one2
                derv2(jxc,kyc) = derv2(jxc,kyc) - xd1*yd2*d2nn(ind)*one2
                derv2(jxc,kzc) = derv2(jxc,kzc) - xd1*zd2*d2nn(ind)*one2
                derv2(jyc,kxc) = derv2(jyc,kxc) - yd1*xd2*d2nn(ind)*one2
                derv2(jyc,kyc) = derv2(jyc,kyc) - yd1*yd2*d2nn(ind)*one2
                derv2(jyc,kzc) = derv2(jyc,kzc) - yd1*zd2*d2nn(ind)*one2
                derv2(jzc,kxc) = derv2(jzc,kxc) - zd1*xd2*d2nn(ind)*one2
                derv2(jzc,kyc) = derv2(jzc,kyc) - zd1*yd2*d2nn(ind)*one2
                derv2(jzc,kzc) = derv2(jzc,kzc) - zd1*zd2*d2nn(ind)*one2
                if (n1.eq.n2) then
                  derv2(jxc,kxc) = derv2(jxc,kxc) - d1nn(n2)*one2
                  derv2(jyc,kyc) = derv2(jyc,kyc) - d1nn(n2)*one2
                  derv2(jzc,kzc) = derv2(jzc,kzc) - d1nn(n2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if ((lopj.or.lopl).and.n1j.ne.n2l) then
              if (lopj.and.lopl) then
                jxc = jx
                jyc = jy
                jzc = jz
                lxc = lx
                lyc = ly
                lzc = lz
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                jxc = jx
                jyc = jy
                jzc = jz
                lxc = jx
                lyc = jy
                lzc = jz
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                jxc = lx
                jyc = ly
                jzc = lz
                lxc = lx
                lyc = ly
                lzc = lz
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              derv2(lxc,jxc) = derv2(lxc,jxc) + xd1*xd2*d2nn(ind)*one1
              derv2(lyc,jxc) = derv2(lyc,jxc) + xd1*yd2*d2nn(ind)*one1
              derv2(lzc,jxc) = derv2(lzc,jxc) + xd1*zd2*d2nn(ind)*one1
              derv2(lxc,jyc) = derv2(lxc,jyc) + yd1*xd2*d2nn(ind)*one1
              derv2(lyc,jyc) = derv2(lyc,jyc) + yd1*yd2*d2nn(ind)*one1
              derv2(lzc,jyc) = derv2(lzc,jyc) + yd1*zd2*d2nn(ind)*one1
              derv2(lxc,jzc) = derv2(lxc,jzc) + zd1*xd2*d2nn(ind)*one1
              derv2(lyc,jzc) = derv2(lyc,jzc) + zd1*yd2*d2nn(ind)*one1
              derv2(lzc,jzc) = derv2(lzc,jzc) + zd1*zd2*d2nn(ind)*one1
              if (n1.eq.n2) then
                derv2(lxc,jxc) = derv2(lxc,jxc) + d1nn(n2)*one1
                derv2(lyc,jyc) = derv2(lyc,jyc) + d1nn(n2)*one1
                derv2(lzc,jzc) = derv2(lzc,jzc) + d1nn(n2)*one1
              endif
              if (.not.lsameijkl) then
                derv2(jxc,lxc) = derv2(jxc,lxc) + xd1*xd2*d2nn(ind)*one2
                derv2(jxc,lyc) = derv2(jxc,lyc) + xd1*yd2*d2nn(ind)*one2
                derv2(jxc,lzc) = derv2(jxc,lzc) + xd1*zd2*d2nn(ind)*one2
                derv2(jyc,lxc) = derv2(jyc,lxc) + yd1*xd2*d2nn(ind)*one2
                derv2(jyc,lyc) = derv2(jyc,lyc) + yd1*yd2*d2nn(ind)*one2
                derv2(jyc,lzc) = derv2(jyc,lzc) + yd1*zd2*d2nn(ind)*one2
                derv2(jzc,lxc) = derv2(jzc,lxc) + zd1*xd2*d2nn(ind)*one2
                derv2(jzc,lyc) = derv2(jzc,lyc) + zd1*yd2*d2nn(ind)*one2
                derv2(jzc,lzc) = derv2(jzc,lzc) + zd1*zd2*d2nn(ind)*one2
                if (n1.eq.n2) then
                  derv2(jxc,lxc) = derv2(jxc,lxc) + d1nn(n2)*one2
                  derv2(jyc,lyc) = derv2(jyc,lyc) + d1nn(n2)*one2
                  derv2(jzc,lzc) = derv2(jzc,lzc) + d1nn(n2)*one2
                endif
              endif
            endif
            if (lstr) then
!
!  Strain - strain second derivatives
!
              do is1 = 1,nstrains
                ns1 = nstrptr(is1)
                do is2 = 1,nstrains
                  ns2 = nstrptr(is2)
                  sderv2(is2,is1) = sderv2(is2,is1) + d2nn(ind)*dr2ds2(ns2)*dr2ds1(ns1)
                  if (n1.ne.n2) then
                    sderv2(is2,is1) = sderv2(is2,is1) + d2nn(ind)*dr2ds1(ns2)*dr2ds2(ns1)
                  else
                    sderv2(is2,is1) = sderv2(is2,is1) + d1nn(n2)*d2r2ds22(ns2,ns1)
                  endif
                enddo
              enddo
!
!  Internal - strain second derivatives
!
              do is1 = 1,nstrains
                ns1 = nstrptr(is1)
                derv3(ix,is1) = derv3(ix,is1) - xd1*d2nn(ind)*dr2ds2(ns1)
                derv3(iy,is1) = derv3(iy,is1) - yd1*d2nn(ind)*dr2ds2(ns1)
                derv3(iz,is1) = derv3(iz,is1) - zd1*d2nn(ind)*dr2ds2(ns1)
                derv3(jx,is1) = derv3(jx,is1) + xd1*d2nn(ind)*dr2ds2(ns1)
                derv3(jy,is1) = derv3(jy,is1) + yd1*d2nn(ind)*dr2ds2(ns1)
                derv3(jz,is1) = derv3(jz,is1) + zd1*d2nn(ind)*dr2ds2(ns1)
                if (n1.ne.n2) then
                  derv3(kx,is1) = derv3(kx,is1) - xd2*d2nn(ind)*dr2ds1(ns1)
                  derv3(ky,is1) = derv3(ky,is1) - yd2*d2nn(ind)*dr2ds1(ns1)
                  derv3(kz,is1) = derv3(kz,is1) - zd2*d2nn(ind)*dr2ds1(ns1)
                  derv3(lx,is1) = derv3(lx,is1) + xd2*d2nn(ind)*dr2ds1(ns1)
                  derv3(ly,is1) = derv3(ly,is1) + yd2*d2nn(ind)*dr2ds1(ns1)
                  derv3(lz,is1) = derv3(lz,is1) + zd2*d2nn(ind)*dr2ds1(ns1)
                else
                  derv3(ix,is1) = derv3(ix,is1) - d1nn(n2)*d2r2dsdx2(ns1,1)
                  derv3(iy,is1) = derv3(iy,is1) - d1nn(n2)*d2r2dsdx2(ns1,2)
                  derv3(iz,is1) = derv3(iz,is1) - d1nn(n2)*d2r2dsdx2(ns1,3)
                  derv3(jx,is1) = derv3(jx,is1) + d1nn(n2)*d2r2dsdx2(ns1,1)
                  derv3(jy,is1) = derv3(jy,is1) + d1nn(n2)*d2r2dsdx2(ns1,2)
                  derv3(jz,is1) = derv3(jz,is1) + d1nn(n2)*d2r2dsdx2(ns1,3)
                endif
              enddo
            endif
          endif
        endif
      enddo
    else
!
!  Increment ind pointer by number of missed terms
!
      ind = ind + n1
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addnn')
#endif
!
  return
  end
!
  subroutine d2adds(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nfreeatomau, &
                    nREBOatomRptr,nauatom,neqvatom,d1,d2,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array.
!
!  Symmetry adapted version of d2add
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move in full cell
!  nfreeatomau     = pointer from atom to free atom to move in asymmetric unit
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  nauatom         = pointer from atom to asymmetric unit atom
!  neqvatom        = number of equivalent positions for asymmetric unit atom
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!
!   9/02 Created from d2add
!  12/03 Torsional derivatives made optional
!   6/07 nREBOatomRptr added as input
!  10/07 Correction made to setting of n1i/n2k for torsional case
!  11/07 Unused variables removed
!   1/10 Bug fixed for neighbour of i -> neighbour of neighbour of i case
!   2/18 Trace added
!   9/18 Strain module added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use current,        only : ndim, nstrains
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nfreeatomau(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  integer(i4), intent(in)             :: nauatom(*)
  integer(i4), intent(in)             :: neqvatom(*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: is1,is2
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: ixa,iya,iza
  integer(i4)                         :: jxa,jya,jza
  integer(i4)                         :: kxa,kya,kza
  integer(i4)                         :: lxa,lya,lza
  integer(i4)                         :: ixc,iyc,izc
  integer(i4)                         :: jxc,jyc,jzc
  integer(i4)                         :: kxc,kyc,kzc
  integer(i4)                         :: lxc,lyc,lzc
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: nji
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: ns1
  integer(i4)                         :: ns2
  integer(i4)                         :: ntmp
  real(dp)                            :: dr2ds1(6)
  real(dp)                            :: dr2ds2(6)
  real(dp)                            :: d2r2ds21(6,6)
  real(dp)                            :: d2r2ds22(6,6)
  real(dp)                            :: d2r2dsdx1(6,3)
  real(dp)                            :: d2r2dsdx2(6,3)
  real(dp)                            :: d2r2dx21(3,3)
  real(dp)                            :: d2r2dx22(3,3)
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: rneq
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: lsameijkl
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lopi
  logical                             :: lopj
  logical                             :: lopk
  logical                             :: lopl
#ifdef TRACE
  call trace_in('d2adds')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndist = nneigh2 + nneigh(nri)*(maxneigh + 1)*mneigh
  else
    ndist = nneigh2
  endif
!
!  Loop over pairs of interatomic distances
!
  ind = 0
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Check that neighbour numbers are valid
!
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    else
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
!  Check that atoms are valid
!
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
      else
        ljok = .false.
      endif
!
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!         
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
!
!  Set up products for strain
!
      if (lstr) then
        call real1strterm(ndim,xd1,yd1,zd1,0.0_dp,0.0_dp,0.0_dp,dr2ds1,d2r2dx21,d2r2dsdx1,d2r2ds21,.true.)
      endif
!
!  Set second derivative location pointers
!
      lopi = (nfreeatom(n1i).gt.0)
      lopj = (nfreeatom(n1j).gt.0)
      ix = 3*(nfreeatom(n1i) - 1) + 1
      iy = ix + 1
      iz = iy + 1
      jx = 3*(nfreeatom(n1j) - 1) + 1
      jy = jx + 1
      jz = jy + 1
      if (nauatom(n1i).gt.0) then
        ixa = 3*(nfreeatomau(nauatom(n1i)) - 1) + 1
        iya = ixa + 1
        iza = iya + 1
      endif
      if (nauatom(n1j).gt.0) then
        jxa = 3*(nfreeatomau(nauatom(n1j)) - 1) + 1
        jya = jxa + 1
        jza = jya + 1
      endif
    endif
!
!  If neither i nor j are being optimised and strain derivatives are not needed then atoms are not valid
!
    if (.not.lopi.and..not.lopj.and..not.lstr) then
      liok = .false.
      ljok = .false.
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!  
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          else
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
            else
              llok = .false.
            endif
!
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!               
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
!  Set up products for strain
!
            if (lstr) then
              call real1strterm(ndim,xd2,yd2,zd2,0.0_dp,0.0_dp,0.0_dp,dr2ds2,d2r2dx22,d2r2dsdx2,d2r2ds22,.true.)
            endif
!
            lopk = (nfreeatom(n2k).gt.0)
            lopl = (nfreeatom(n2l).gt.0)
            kx = 3*(nfreeatom(n2k) - 1) + 1
            ky = kx + 1
            kz = ky + 1
            lx = 3*(nfreeatom(n2l) - 1) + 1
            ly = lx + 1
            lz = ly + 1
            if (nauatom(n2k).gt.0) then
              kxa = 3*(nfreeatomau(nauatom(n2k)) - 1) + 1
              kya = kxa + 1
              kza = kya + 1
            endif
            if (nauatom(n2l).gt.0) then
              lxa = 3*(nfreeatomau(nauatom(n2l)) - 1) + 1
              lya = lxa + 1
              lza = lya + 1
            endif
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if ((lopi.or.lopk).and.n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              if (lopi.and.nauatom(n1i).gt.0) then
                rneq = dble(neqvatom(nauatom(n1i)))
                if (lopk) then
                  kxc = kx
                  kyc = ky
                  kzc = kz
                else
                  kxc = ix
                  kyc = iy
                  kzc = iz
                endif
                derv2(kxc,ixa) = derv2(kxc,ixa) + xd1*xd2*d2(ind)*rneq
                derv2(kyc,ixa) = derv2(kyc,ixa) + xd1*yd2*d2(ind)*rneq
                derv2(kzc,ixa) = derv2(kzc,ixa) + xd1*zd2*d2(ind)*rneq
                derv2(kxc,iya) = derv2(kxc,iya) + yd1*xd2*d2(ind)*rneq
                derv2(kyc,iya) = derv2(kyc,iya) + yd1*yd2*d2(ind)*rneq
                derv2(kzc,iya) = derv2(kzc,iya) + yd1*zd2*d2(ind)*rneq
                derv2(kxc,iza) = derv2(kxc,iza) + zd1*xd2*d2(ind)*rneq
                derv2(kyc,iza) = derv2(kyc,iza) + zd1*yd2*d2(ind)*rneq
                derv2(kzc,iza) = derv2(kzc,iza) + zd1*zd2*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(kxc,ixa) = derv2(kxc,ixa) + d1(n2)*rneq
                  derv2(kyc,iya) = derv2(kyc,iya) + d1(n2)*rneq
                  derv2(kzc,iza) = derv2(kzc,iza) + d1(n2)*rneq
                endif
              endif
              if (lopk.and.nauatom(n2k).gt.0.and..not.lsameijkl) then
                rneq = dble(neqvatom(nauatom(n2k)))
                if (lopi) then
                  ixc = ix
                  iyc = iy
                  izc = iz
                else
                  ixc = kx 
                  iyc = ky 
                  izc = kz
                endif
                derv2(ixc,kxa) = derv2(ixc,kxa) + xd2*xd1*d2(ind)*rneq
                derv2(iyc,kxa) = derv2(iyc,kxa) + xd2*yd1*d2(ind)*rneq
                derv2(izc,kxa) = derv2(izc,kxa) + xd2*zd1*d2(ind)*rneq
                derv2(ixc,kya) = derv2(ixc,kya) + yd2*xd1*d2(ind)*rneq
                derv2(iyc,kya) = derv2(iyc,kya) + yd2*yd1*d2(ind)*rneq
                derv2(izc,kya) = derv2(izc,kya) + yd2*zd1*d2(ind)*rneq
                derv2(ixc,kza) = derv2(ixc,kza) + zd2*xd1*d2(ind)*rneq
                derv2(iyc,kza) = derv2(iyc,kza) + zd2*yd1*d2(ind)*rneq
                derv2(izc,kza) = derv2(izc,kza) + zd2*zd1*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(ixc,kxa) = derv2(ixc,kxa) + d1(n2)*rneq
                  derv2(iyc,kya) = derv2(iyc,kya) + d1(n2)*rneq
                  derv2(izc,kza) = derv2(izc,kza) + d1(n2)*rneq
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if ((lopi.or.lopl).and.n1i.ne.n2l) then
              if (lopi.and.nauatom(n1i).gt.0) then
                rneq = dble(neqvatom(nauatom(n1i)))
                if (lopl) then
                  lxc = lx
                  lyc = ly
                  lzc = lz
                else
                  lxc = ix
                  lyc = iy
                  lzc = iz
                endif
                derv2(lxc,ixa) = derv2(lxc,ixa) - xd1*xd2*d2(ind)*rneq
                derv2(lyc,ixa) = derv2(lyc,ixa) - xd1*yd2*d2(ind)*rneq
                derv2(lzc,ixa) = derv2(lzc,ixa) - xd1*zd2*d2(ind)*rneq
                derv2(lxc,iya) = derv2(lxc,iya) - yd1*xd2*d2(ind)*rneq
                derv2(lyc,iya) = derv2(lyc,iya) - yd1*yd2*d2(ind)*rneq
                derv2(lzc,iya) = derv2(lzc,iya) - yd1*zd2*d2(ind)*rneq
                derv2(lxc,iza) = derv2(lxc,iza) - zd1*xd2*d2(ind)*rneq
                derv2(lyc,iza) = derv2(lyc,iza) - zd1*yd2*d2(ind)*rneq
                derv2(lzc,iza) = derv2(lzc,iza) - zd1*zd2*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(lxc,ixa) = derv2(lxc,ixa) - d1(n2)*rneq
                  derv2(lyc,iya) = derv2(lyc,iya) - d1(n2)*rneq
                  derv2(lzc,iza) = derv2(lzc,iza) - d1(n2)*rneq
                endif
              endif
              if (lopl.and.nauatom(n2l).gt.0.and..not.lsameijkl) then
                rneq = dble(neqvatom(nauatom(n2l)))
                if (lopi) then
                  ixc = ix
                  iyc = iy
                  izc = iz
                else
                  ixc = lx 
                  iyc = ly 
                  izc = lz
                endif
                derv2(ixc,lxa) = derv2(ixc,lxa) - xd2*xd1*d2(ind)*rneq
                derv2(iyc,lxa) = derv2(iyc,lxa) - xd2*yd1*d2(ind)*rneq
                derv2(izc,lxa) = derv2(izc,lxa) - xd2*zd1*d2(ind)*rneq
                derv2(ixc,lya) = derv2(ixc,lya) - yd2*xd1*d2(ind)*rneq
                derv2(iyc,lya) = derv2(iyc,lya) - yd2*yd1*d2(ind)*rneq
                derv2(izc,lya) = derv2(izc,lya) - yd2*zd1*d2(ind)*rneq
                derv2(ixc,lza) = derv2(ixc,lza) - zd2*xd1*d2(ind)*rneq
                derv2(iyc,lza) = derv2(iyc,lza) - zd2*yd1*d2(ind)*rneq
                derv2(izc,lza) = derv2(izc,lza) - zd2*zd1*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(ixc,lxa) = derv2(ixc,lxa) - d1(n2)*rneq
                  derv2(iyc,lya) = derv2(iyc,lya) - d1(n2)*rneq
                  derv2(izc,lza) = derv2(izc,lza) - d1(n2)*rneq
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if ((lopj.or.lopk).and.n1j.ne.n2k) then
              if (lopj.and.nauatom(n1j).gt.0) then
                rneq = dble(neqvatom(nauatom(n1j)))
                if (lopk) then
                  kxc = kx
                  kyc = ky
                  kzc = kz
                else
                  kxc = jx
                  kyc = jy
                  kzc = jz
                endif
                derv2(kxc,jxa) = derv2(kxc,jxa) - xd1*xd2*d2(ind)*rneq
                derv2(kyc,jxa) = derv2(kyc,jxa) - xd1*yd2*d2(ind)*rneq
                derv2(kzc,jxa) = derv2(kzc,jxa) - xd1*zd2*d2(ind)*rneq
                derv2(kxc,jya) = derv2(kxc,jya) - yd1*xd2*d2(ind)*rneq
                derv2(kyc,jya) = derv2(kyc,jya) - yd1*yd2*d2(ind)*rneq
                derv2(kzc,jya) = derv2(kzc,jya) - yd1*zd2*d2(ind)*rneq
                derv2(kxc,jza) = derv2(kxc,jza) - zd1*xd2*d2(ind)*rneq
                derv2(kyc,jza) = derv2(kyc,jza) - zd1*yd2*d2(ind)*rneq
                derv2(kzc,jza) = derv2(kzc,jza) - zd1*zd2*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(kxc,jxa) = derv2(kxc,jxa) - d1(n2)*rneq
                  derv2(kyc,jya) = derv2(kyc,jya) - d1(n2)*rneq
                  derv2(kzc,jza) = derv2(kzc,jza) - d1(n2)*rneq
                endif
              endif
              if (lopk.and.nauatom(n2k).gt.0.and..not.lsameijkl) then
                rneq = dble(neqvatom(nauatom(n2k)))
                if (lopj) then
                  jxc = jx
                  jyc = jy
                  jzc = jz
                else
                  jxc = kx 
                  jyc = ky 
                  jzc = kz
                endif
                derv2(jxc,kxa) = derv2(jxc,kxa) - xd2*xd1*d2(ind)*rneq
                derv2(jyc,kxa) = derv2(jyc,kxa) - xd2*yd1*d2(ind)*rneq
                derv2(jzc,kxa) = derv2(jzc,kxa) - xd2*zd1*d2(ind)*rneq
                derv2(jxc,kya) = derv2(jxc,kya) - yd2*xd1*d2(ind)*rneq
                derv2(jyc,kya) = derv2(jyc,kya) - yd2*yd1*d2(ind)*rneq
                derv2(jzc,kya) = derv2(jzc,kya) - yd2*zd1*d2(ind)*rneq
                derv2(jxc,kza) = derv2(jxc,kza) - zd2*xd1*d2(ind)*rneq
                derv2(jyc,kza) = derv2(jyc,kza) - zd2*yd1*d2(ind)*rneq
                derv2(jzc,kza) = derv2(jzc,kza) - zd2*zd1*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(jxc,kxa) = derv2(jxc,kxa) - d1(n2)*rneq
                  derv2(jyc,kya) = derv2(jyc,kya) - d1(n2)*rneq
                  derv2(jzc,kza) = derv2(jzc,kza) - d1(n2)*rneq
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if ((lopj.or.lopl).and.n1j.ne.n2l) then
              if (lopj.and.nauatom(n1j).gt.0) then
                rneq = dble(neqvatom(nauatom(n1j)))
                if (lopl) then
                  lxc = lx
                  lyc = ly
                  lzc = lz
                else
                  lxc = jx
                  lyc = jy
                  lzc = jz
                endif
                derv2(lxc,jxa) = derv2(lxc,jxa) + xd1*xd2*d2(ind)*rneq
                derv2(lyc,jxa) = derv2(lyc,jxa) + xd1*yd2*d2(ind)*rneq
                derv2(lzc,jxa) = derv2(lzc,jxa) + xd1*zd2*d2(ind)*rneq
                derv2(lxc,jya) = derv2(lxc,jya) + yd1*xd2*d2(ind)*rneq
                derv2(lyc,jya) = derv2(lyc,jya) + yd1*yd2*d2(ind)*rneq
                derv2(lzc,jya) = derv2(lzc,jya) + yd1*zd2*d2(ind)*rneq
                derv2(lxc,jza) = derv2(lxc,jza) + zd1*xd2*d2(ind)*rneq
                derv2(lyc,jza) = derv2(lyc,jza) + zd1*yd2*d2(ind)*rneq
                derv2(lzc,jza) = derv2(lzc,jza) + zd1*zd2*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(lxc,jxa) = derv2(lxc,jxa) + d1(n2)*rneq
                  derv2(lyc,jya) = derv2(lyc,jya) + d1(n2)*rneq
                  derv2(lzc,jza) = derv2(lzc,jza) + d1(n2)*rneq
                endif
              endif
              if (lopl.and.nauatom(n2l).gt.0.and..not.lsameijkl) then
                rneq = dble(neqvatom(nauatom(n2l)))
                if (lopj) then
                  jxc = jx
                  jyc = jy
                  jzc = jz
                else
                  jxc = lx 
                  jyc = ly 
                  jzc = lz
                endif
                derv2(jxc,lxa) = derv2(jxc,lxa) + xd2*xd1*d2(ind)*rneq
                derv2(jyc,lxa) = derv2(jyc,lxa) + xd2*yd1*d2(ind)*rneq
                derv2(jzc,lxa) = derv2(jzc,lxa) + xd2*zd1*d2(ind)*rneq
                derv2(jxc,lya) = derv2(jxc,lya) + yd2*xd1*d2(ind)*rneq
                derv2(jyc,lya) = derv2(jyc,lya) + yd2*yd1*d2(ind)*rneq
                derv2(jzc,lya) = derv2(jzc,lya) + yd2*zd1*d2(ind)*rneq
                derv2(jxc,lza) = derv2(jxc,lza) + zd2*xd1*d2(ind)*rneq
                derv2(jyc,lza) = derv2(jyc,lza) + zd2*yd1*d2(ind)*rneq
                derv2(jzc,lza) = derv2(jzc,lza) + zd2*zd1*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(jxc,lxa) = derv2(jxc,lxa) + d1(n2)*rneq
                  derv2(jyc,lya) = derv2(jyc,lya) + d1(n2)*rneq
                  derv2(jzc,lza) = derv2(jzc,lza) + d1(n2)*rneq
                endif
              endif
            endif
            if (lstr) then
!
!  Strain - strain second derivatives
!
              do is1 = 1,nstrains
                ns1 = nstrptr(is1)
                do is2 = 1,nstrains
                  ns2 = nstrptr(is2)
                  sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds2(ns2)*dr2ds1(ns1)
                  if (n1.ne.n2) then
                    sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds1(ns2)*dr2ds2(ns1)
                  else
                    sderv2(is2,is1) = sderv2(is2,is1) + d1(n2)*d2r2ds22(ns2,ns1)
                  endif
                enddo
              enddo
!
!  Internal - strain second derivatives
!
              if (lopi.and.nauatom(n1i).gt.0) then
                rneq = dble(neqvatom(nauatom(n1i)))
                do is1 = 1,nstrains
                  ns1 = nstrptr(is1)
                  derv3(ixa,is1) = derv3(ixa,is1) - xd1*d2(ind)*dr2ds2(ns1)*rneq 
                  derv3(iya,is1) = derv3(iya,is1) - yd1*d2(ind)*dr2ds2(ns1)*rneq
                  derv3(iza,is1) = derv3(iza,is1) - zd1*d2(ind)*dr2ds2(ns1)*rneq
                enddo
                if (n1.eq.n2) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(ixa,is1) = derv3(ixa,is1) - d1(n2)*d2r2dsdx2(ns1,1)*rneq 
                    derv3(iya,is1) = derv3(iya,is1) - d1(n2)*d2r2dsdx2(ns1,2)*rneq
                    derv3(iza,is1) = derv3(iza,is1) - d1(n2)*d2r2dsdx2(ns1,3)*rneq
                  enddo
                endif
              endif
              if (lopj.and.nauatom(n1j).gt.0) then
                rneq = dble(neqvatom(nauatom(n1j)))
                do is1 = 1,nstrains
                  ns1 = nstrptr(is1)
                  derv3(jxa,is1) = derv3(jxa,is1) + xd1*d2(ind)*dr2ds2(ns1)*rneq
                  derv3(jya,is1) = derv3(jya,is1) + yd1*d2(ind)*dr2ds2(ns1)*rneq
                  derv3(jza,is1) = derv3(jza,is1) + zd1*d2(ind)*dr2ds2(ns1)*rneq
                enddo
                if (n1.eq.n2) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(jxa,is1) = derv3(jxa,is1) + d1(n2)*d2r2dsdx2(ns1,1)*rneq 
                    derv3(jya,is1) = derv3(jya,is1) + d1(n2)*d2r2dsdx2(ns1,2)*rneq
                    derv3(jza,is1) = derv3(jza,is1) + d1(n2)*d2r2dsdx2(ns1,3)*rneq
                  enddo
                endif
              endif
              if (n1.ne.n2) then
                if (lopk.and.nauatom(n2k).gt.0) then
                  rneq = dble(neqvatom(nauatom(n2k)))
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(kxa,is1) = derv3(kxa,is1) - xd2*d2(ind)*dr2ds1(ns1)*rneq
                    derv3(kya,is1) = derv3(kya,is1) - yd2*d2(ind)*dr2ds1(ns1)*rneq
                    derv3(kza,is1) = derv3(kza,is1) - zd2*d2(ind)*dr2ds1(ns1)*rneq
                  enddo
                endif
                if (lopl.and.nauatom(n2l).gt.0) then
                  rneq = dble(neqvatom(nauatom(n2l)))
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(lxa,is1) = derv3(lxa,is1) + xd2*d2(ind)*dr2ds1(ns1)*rneq
                    derv3(lya,is1) = derv3(lya,is1) + yd2*d2(ind)*dr2ds1(ns1)*rneq
                    derv3(lza,is1) = derv3(lza,is1) + zd2*d2(ind)*dr2ds1(ns1)*rneq
                  enddo
                endif
              endif
            endif
          endif
        endif
      enddo
    else
!
!  Increment ind pointer by number of missed terms
!
      ind = ind + n1
    endif
  enddo
#ifdef TRACE
  call trace_out('d2adds')
#endif
!
  return
  end
!
  subroutine d2addd(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr, &
                    nregion1,d1,d2,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array. Defect
!  calculation version.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  nregion1        = number of ions in region 1
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!
!   9/02 Created from d2add
!  12/03 Torsional derivatives made optional
!   6/07 nREBOatomRptr added as input
!  10/07 Correction made to setting of n1i/n2k for torsional case
!  11/07 Unused variables removed
!   1/10 Bug fixed for neighbour of i -> neighbour of neighbour of i case
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  integer(i4), intent(in)             :: nregion1
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: nji
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: ntmp
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: lsameijkl
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lopi
  logical                             :: lopj
  logical                             :: lopk
  logical                             :: lopl
#ifdef TRACE
  call trace_in('d2addd')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndist = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
  else
    ndist = nneigh2
  endif
!
!  Loop over pairs of interatomic distances
!
  ind = 0
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Check that neighbour numbers are valid
!
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    else
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
!  Check that atoms are valid
!
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        ljok = (n1j.le.nneigh(neighno(nji,nri)))
      else
        ljok = .false.
      endif
!
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!         
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
!
!  Set second derivative location pointers
!
      lopi = (nfreeatom(n1i).gt.0)
      lopj = (nfreeatom(n1j).gt.0)
      if (n1i.le.nregion1) then
        ix = 3*(nfreeatom(n1i) - 1) + 1
        iy = ix + 1
        iz = iy + 1
      else
        ix = 3*nregion1 + 1
        iy = ix + 1
        iz = iy + 1
      endif
      if (n1j.le.nregion1) then
        jx = 3*(nfreeatom(n1j) - 1) + 1
        jy = jx + 1
        jz = jy + 1
      else
        jx = 3*nregion1 + 1
        jy = jx + 1
        jz = jy + 1
      endif
    endif
!
!  If neither i nor j are being optimised then atoms are not valid
!
    if (.not.lopi.and..not.lopj) then
      liok = .false.
      ljok = .false.
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!  
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          else
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
            else
              llok = .false.
            endif
!
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!               
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
            lopk = (nfreeatom(n2k).gt.0)
            lopl = (nfreeatom(n2l).gt.0)
            if (n2k.le.nregion1) then
              kx = 3*(nfreeatom(n2k) - 1) + 1
              ky = kx + 1
              kz = ky + 1
            else
              kx = 3*nregion1 + 1
              ky = kx + 1
              kz = ky + 1
            endif
            if (n2l.le.nregion1) then
              lx = 3*(nfreeatom(n2l) - 1) + 1
              ly = lx + 1
              lz = ly + 1
            else
              lx = 3*nregion1 + 1
              ly = lx + 1
              lz = ly + 1
            endif
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if ((lopi.or.lopk).and.n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              derv2(kx,ix) = derv2(kx,ix) + xd1*xd2*d2(ind)
              derv2(ky,ix) = derv2(ky,ix) + xd1*yd2*d2(ind)
              derv2(kz,ix) = derv2(kz,ix) + xd1*zd2*d2(ind)
              derv2(kx,iy) = derv2(kx,iy) + yd1*xd2*d2(ind)
              derv2(ky,iy) = derv2(ky,iy) + yd1*yd2*d2(ind)
              derv2(kz,iy) = derv2(kz,iy) + yd1*zd2*d2(ind)
              derv2(kx,iz) = derv2(kx,iz) + zd1*xd2*d2(ind)
              derv2(ky,iz) = derv2(ky,iz) + zd1*yd2*d2(ind)
              derv2(kz,iz) = derv2(kz,iz) + zd1*zd2*d2(ind)
              if (n1.eq.n2) then
                derv2(kx,ix) = derv2(kx,ix) + d1(n2)
                derv2(ky,iy) = derv2(ky,iy) + d1(n2)
                derv2(kz,iz) = derv2(kz,iz) + d1(n2)
              endif
              if (.not.lsameijkl) then
                derv2(ix,kx) = derv2(ix,kx) + xd1*xd2*d2(ind)
                derv2(ix,ky) = derv2(ix,ky) + xd1*yd2*d2(ind)
                derv2(ix,kz) = derv2(ix,kz) + xd1*zd2*d2(ind)
                derv2(iy,kx) = derv2(iy,kx) + yd1*xd2*d2(ind)
                derv2(iy,ky) = derv2(iy,ky) + yd1*yd2*d2(ind)
                derv2(iy,kz) = derv2(iy,kz) + yd1*zd2*d2(ind)
                derv2(iz,kx) = derv2(iz,kx) + zd1*xd2*d2(ind)
                derv2(iz,ky) = derv2(iz,ky) + zd1*yd2*d2(ind)
                derv2(iz,kz) = derv2(iz,kz) + zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(ix,kx) = derv2(ix,kx) + d1(n2)
                  derv2(iy,ky) = derv2(iy,ky) + d1(n2)
                  derv2(iz,kz) = derv2(iz,kz) + d1(n2)
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if ((lopi.or.lopl).and.n1i.ne.n2l) then
              derv2(lx,ix) = derv2(lx,ix) - xd1*xd2*d2(ind)
              derv2(ly,ix) = derv2(ly,ix) - xd1*yd2*d2(ind)
              derv2(lz,ix) = derv2(lz,ix) - xd1*zd2*d2(ind)
              derv2(lx,iy) = derv2(lx,iy) - yd1*xd2*d2(ind)
              derv2(ly,iy) = derv2(ly,iy) - yd1*yd2*d2(ind)
              derv2(lz,iy) = derv2(lz,iy) - yd1*zd2*d2(ind)
              derv2(lx,iz) = derv2(lx,iz) - zd1*xd2*d2(ind)
              derv2(ly,iz) = derv2(ly,iz) - zd1*yd2*d2(ind)
              derv2(lz,iz) = derv2(lz,iz) - zd1*zd2*d2(ind)
              if (n1.eq.n2) then
                derv2(lx,ix) = derv2(lx,ix) - d1(n2)
                derv2(ly,iy) = derv2(ly,iy) - d1(n2)
                derv2(lz,iz) = derv2(lz,iz) - d1(n2)
              endif
              if (.not.lsameijkl) then
                derv2(ix,lx) = derv2(ix,lx) - xd1*xd2*d2(ind)
                derv2(ix,ly) = derv2(ix,ly) - xd1*yd2*d2(ind)
                derv2(ix,lz) = derv2(ix,lz) - xd1*zd2*d2(ind)
                derv2(iy,lx) = derv2(iy,lx) - yd1*xd2*d2(ind)
                derv2(iy,ly) = derv2(iy,ly) - yd1*yd2*d2(ind)
                derv2(iy,lz) = derv2(iy,lz) - yd1*zd2*d2(ind)
                derv2(iz,lx) = derv2(iz,lx) - zd1*xd2*d2(ind)
                derv2(iz,ly) = derv2(iz,ly) - zd1*yd2*d2(ind)
                derv2(iz,lz) = derv2(iz,lz) - zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(ix,lx) = derv2(ix,lx) - d1(n2)
                  derv2(iy,ly) = derv2(iy,ly) - d1(n2)
                  derv2(iz,lz) = derv2(iz,lz) - d1(n2)
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if ((lopj.or.lopk).and.n1j.ne.n2k) then
              derv2(kx,jx) = derv2(kx,jx) - xd1*xd2*d2(ind)
              derv2(ky,jx) = derv2(ky,jx) - xd1*yd2*d2(ind)
              derv2(kz,jx) = derv2(kz,jx) - xd1*zd2*d2(ind)
              derv2(kx,jy) = derv2(kx,jy) - yd1*xd2*d2(ind)
              derv2(ky,jy) = derv2(ky,jy) - yd1*yd2*d2(ind)
              derv2(kz,jy) = derv2(kz,jy) - yd1*zd2*d2(ind)
              derv2(kx,jz) = derv2(kx,jz) - zd1*xd2*d2(ind)
              derv2(ky,jz) = derv2(ky,jz) - zd1*yd2*d2(ind)
              derv2(kz,jz) = derv2(kz,jz) - zd1*zd2*d2(ind)
              if (n1.eq.n2) then
                derv2(kx,jx) = derv2(kx,jx) - d1(n2)
                derv2(ky,jy) = derv2(ky,jy) - d1(n2)
                derv2(kz,jz) = derv2(kz,jz) - d1(n2)
              endif
              if (.not.lsameijkl) then
                derv2(jx,kx) = derv2(jx,kx) - xd1*xd2*d2(ind)
                derv2(jx,ky) = derv2(jx,ky) - xd1*yd2*d2(ind)
                derv2(jx,kz) = derv2(jx,kz) - xd1*zd2*d2(ind)
                derv2(jy,kx) = derv2(jy,kx) - yd1*xd2*d2(ind)
                derv2(jy,ky) = derv2(jy,ky) - yd1*yd2*d2(ind)
                derv2(jy,kz) = derv2(jy,kz) - yd1*zd2*d2(ind)
                derv2(jz,kx) = derv2(jz,kx) - zd1*xd2*d2(ind)
                derv2(jz,ky) = derv2(jz,ky) - zd1*yd2*d2(ind)
                derv2(jz,kz) = derv2(jz,kz) - zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(jx,kx) = derv2(jx,kx) - d1(n2)
                  derv2(jy,ky) = derv2(jy,ky) - d1(n2)
                  derv2(jz,kz) = derv2(jz,kz) - d1(n2)
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if ((lopj.or.lopl).and.n1j.ne.n2l) then
              derv2(lx,jx) = derv2(lx,jx) + xd1*xd2*d2(ind)
              derv2(ly,jx) = derv2(ly,jx) + xd1*yd2*d2(ind)
              derv2(lz,jx) = derv2(lz,jx) + xd1*zd2*d2(ind)
              derv2(lx,jy) = derv2(lx,jy) + yd1*xd2*d2(ind)
              derv2(ly,jy) = derv2(ly,jy) + yd1*yd2*d2(ind)
              derv2(lz,jy) = derv2(lz,jy) + yd1*zd2*d2(ind)
              derv2(lx,jz) = derv2(lx,jz) + zd1*xd2*d2(ind)
              derv2(ly,jz) = derv2(ly,jz) + zd1*yd2*d2(ind)
              derv2(lz,jz) = derv2(lz,jz) + zd1*zd2*d2(ind)
              if (n1.eq.n2) then
                derv2(lx,jx) = derv2(lx,jx) + d1(n2)
                derv2(ly,jy) = derv2(ly,jy) + d1(n2)
                derv2(lz,jz) = derv2(lz,jz) + d1(n2)
              endif
              if (.not.lsameijkl) then
                derv2(jx,lx) = derv2(jx,lx) + xd1*xd2*d2(ind)
                derv2(jx,ly) = derv2(jx,ly) + xd1*yd2*d2(ind)
                derv2(jx,lz) = derv2(jx,lz) + xd1*zd2*d2(ind)
                derv2(jy,lx) = derv2(jy,lx) + yd1*xd2*d2(ind)
                derv2(jy,ly) = derv2(jy,ly) + yd1*yd2*d2(ind)
                derv2(jy,lz) = derv2(jy,lz) + yd1*zd2*d2(ind)
                derv2(jz,lx) = derv2(jz,lx) + zd1*xd2*d2(ind)
                derv2(jz,ly) = derv2(jz,ly) + zd1*yd2*d2(ind)
                derv2(jz,lz) = derv2(jz,lz) + zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(jx,lx) = derv2(jx,lx) + d1(n2)
                  derv2(jy,ly) = derv2(jy,ly) + d1(n2)
                  derv2(jz,lz) = derv2(jz,lz) + d1(n2)
                endif
              endif
            endif
          endif
        endif
      enddo
    else
!
!  Increment ind pointer by number of missed terms
!
      ind = ind + n1
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addd')
#endif
!
  return
  end
!
  subroutine d2addds(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nfreeatomau,nREBOatomRptr, &
                     nauatom,neqvatom,nregion1,d1,d2,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array. Defect
!  calculation version
!
!  Symmetry adapted version of d2addd
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move in full cell
!  nfreeatomau     = pointer from atom to free atom to move in asymmetric unit
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  nauatom         = pointer from atom to asymmetric unit atom
!  neqvatom        = number of equivalent positions for asymmetric unit atom
!  nregion1        = number of atoms in region 1
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!
!   9/02 Created from d2adds
!  12/03 Torsional derivatives made optional
!   6/07 nREBOatomRptr added as input
!  10/07 Correction made to setting of n1i/n2k for torsional case
!  11/07 Unused variables removed
!   1/10 Bug fixed for neighbour of i -> neighbour of neighbour of i case
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nfreeatomau(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  integer(i4), intent(in)             :: nauatom(*)
  integer(i4), intent(in)             :: neqvatom(*)
  integer(i4), intent(in)             :: nregion1
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: ixa,iya,iza
  integer(i4)                         :: jxa,jya,jza
  integer(i4)                         :: kxa,kya,kza
  integer(i4)                         :: lxa,lya,lza
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: nji
  integer(i4)                         :: nri
  integer(i4)                         :: nneigh2
  integer(i4)                         :: ntmp
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: rneq
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: lsameijkl
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lopi
  logical                             :: lopj
  logical                             :: lopk
  logical                             :: lopl
#ifdef TRACE
  call trace_in('d2addds')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndist = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
  else
    ndist = nneigh2
  endif
!
!  Loop over pairs of interatomic distances
!
  ind = 0
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Check that neighbour numbers are valid
!
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    else
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
!  Check that atoms are valid
!
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
      else
        ljok = .false.
      endif
!
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!         
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
!
!  Set second derivative location pointers
!
      lopi = (nfreeatom(n1i).gt.0)
      lopj = (nfreeatom(n1j).gt.0)
      if (n1i.le.nregion1) then
        ix = 3*(nfreeatom(n1i) - 1) + 1
        iy = ix + 1
        iz = iy + 1
      else
        ix = 3*nregion1 + 1
        iy = ix + 1
        iz = iy + 1
      endif
      if (n1j.le.nregion1) then
        jx = 3*(nfreeatom(n1j) - 1) + 1
        jy = jx + 1
        jz = jy + 1
      else
        jx = 3*nregion1 + 1
        jy = jx + 1
        jz = jy + 1
      endif
      if (nauatom(n1i).gt.0) then
        ixa = 3*(nfreeatomau(nauatom(n1i)) - 1) + 1
        iya = ixa + 1
        iza = iya + 1
      endif
      if (nauatom(n1j).gt.0) then
        jxa = 3*(nfreeatomau(nauatom(n1j)) - 1) + 1
        jya = jxa + 1
        jza = jya + 1
      endif
    endif
!
!  If neither i nor j are being optimised then atoms are not valid
!
    if (.not.lopi.and..not.lopj) then
      liok = .false.
      ljok = .false.
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!  
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          else
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
            else
              llok = .false.
            endif
!
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!               
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
            lopk = (nfreeatom(n2k).gt.0)
            lopl = (nfreeatom(n2l).gt.0)
            if (n2k.le.nregion1) then
              kx = 3*(nfreeatom(n2k) - 1) + 1
              ky = kx + 1
              kz = ky + 1
            else
              kx = 3*nregion1 + 1
              ky = kx + 1
              kz = ky + 1
            endif
            if (n2l.le.nregion1) then
              lx = 3*(nfreeatom(n2l) - 1) + 1
              ly = lx + 1
              lz = ly + 1
            else
              lx = 3*nregion1 + 1
              ly = lx + 1
              lz = ly + 1
            endif
            if (nauatom(n2k).gt.0) then
              kxa = 3*(nfreeatomau(nauatom(n2k)) - 1) + 1
              kya = kxa + 1
              kza = kya + 1
            endif
            if (nauatom(n2l).gt.0) then
              lxa = 3*(nfreeatomau(nauatom(n2l)) - 1) + 1
              lya = lxa + 1
              lza = lya + 1
            endif
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if ((lopi.or.lopk).and.n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              if (lopi.and.nauatom(n1i).gt.0) then
                rneq = dble(neqvatom(nauatom(n1i)))
                derv2(kx,ixa) = derv2(kx,ixa) + xd1*xd2*d2(ind)*rneq
                derv2(ky,ixa) = derv2(ky,ixa) + xd1*yd2*d2(ind)*rneq
                derv2(kz,ixa) = derv2(kz,ixa) + xd1*zd2*d2(ind)*rneq
                derv2(kx,iya) = derv2(kx,iya) + yd1*xd2*d2(ind)*rneq
                derv2(ky,iya) = derv2(ky,iya) + yd1*yd2*d2(ind)*rneq
                derv2(kz,iya) = derv2(kz,iya) + yd1*zd2*d2(ind)*rneq
                derv2(kx,iza) = derv2(kx,iza) + zd1*xd2*d2(ind)*rneq
                derv2(ky,iza) = derv2(ky,iza) + zd1*yd2*d2(ind)*rneq
                derv2(kz,iza) = derv2(kz,iza) + zd1*zd2*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(kx,ixa) = derv2(kx,ixa) + d1(n2)*rneq
                  derv2(ky,iya) = derv2(ky,iya) + d1(n2)*rneq
                  derv2(kz,iza) = derv2(kz,iza) + d1(n2)*rneq
                endif
              endif
              if (lopk.and.nauatom(n2k).gt.0.and..not.lsameijkl) then
                rneq = dble(neqvatom(nauatom(n2k)))
                derv2(ix,kxa) = derv2(ix,kxa) + xd2*xd1*d2(ind)*rneq
                derv2(iy,kxa) = derv2(iy,kxa) + xd2*yd1*d2(ind)*rneq
                derv2(iz,kxa) = derv2(iz,kxa) + xd2*zd1*d2(ind)*rneq
                derv2(ix,kya) = derv2(ix,kya) + yd2*xd1*d2(ind)*rneq
                derv2(iy,kya) = derv2(iy,kya) + yd2*yd1*d2(ind)*rneq
                derv2(iz,kya) = derv2(iz,kya) + yd2*zd1*d2(ind)*rneq
                derv2(ix,kza) = derv2(ix,kza) + zd2*xd1*d2(ind)*rneq
                derv2(iy,kza) = derv2(iy,kza) + zd2*yd1*d2(ind)*rneq
                derv2(iz,kza) = derv2(iz,kza) + zd2*zd1*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(ix,kxa) = derv2(ix,kxa) + d1(n2)*rneq
                  derv2(iy,kya) = derv2(iy,kya) + d1(n2)*rneq
                  derv2(iz,kza) = derv2(iz,kza) + d1(n2)*rneq
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if ((lopi.or.lopl).and.n1i.ne.n2l) then
              if (lopi.and.nauatom(n1i).gt.0) then
                rneq = dble(neqvatom(nauatom(n1i)))
                derv2(lx,ixa) = derv2(lx,ixa) - xd1*xd2*d2(ind)*rneq
                derv2(ly,ixa) = derv2(ly,ixa) - xd1*yd2*d2(ind)*rneq
                derv2(lz,ixa) = derv2(lz,ixa) - xd1*zd2*d2(ind)*rneq
                derv2(lx,iya) = derv2(lx,iya) - yd1*xd2*d2(ind)*rneq
                derv2(ly,iya) = derv2(ly,iya) - yd1*yd2*d2(ind)*rneq
                derv2(lz,iya) = derv2(lz,iya) - yd1*zd2*d2(ind)*rneq
                derv2(lx,iza) = derv2(lx,iza) - zd1*xd2*d2(ind)*rneq
                derv2(ly,iza) = derv2(ly,iza) - zd1*yd2*d2(ind)*rneq
                derv2(lz,iza) = derv2(lz,iza) - zd1*zd2*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(lx,ixa) = derv2(lx,ixa) - d1(n2)*rneq
                  derv2(ly,iya) = derv2(ly,iya) - d1(n2)*rneq
                  derv2(lz,iza) = derv2(lz,iza) - d1(n2)*rneq
                endif
              endif
              if (lopl.and.nauatom(n2l).gt.0.and..not.lsameijkl) then
                rneq = dble(neqvatom(nauatom(n2l)))
                derv2(ix,lxa) = derv2(ix,lxa) - xd2*xd1*d2(ind)*rneq
                derv2(iy,lxa) = derv2(iy,lxa) - xd2*yd1*d2(ind)*rneq
                derv2(iz,lxa) = derv2(iz,lxa) - xd2*zd1*d2(ind)*rneq
                derv2(ix,lya) = derv2(ix,lya) - yd2*xd1*d2(ind)*rneq
                derv2(iy,lya) = derv2(iy,lya) - yd2*yd1*d2(ind)*rneq
                derv2(iz,lya) = derv2(iz,lya) - yd2*zd1*d2(ind)*rneq
                derv2(ix,lza) = derv2(ix,lza) - zd2*xd1*d2(ind)*rneq
                derv2(iy,lza) = derv2(iy,lza) - zd2*yd1*d2(ind)*rneq
                derv2(iz,lza) = derv2(iz,lza) - zd2*zd1*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(ix,lxa) = derv2(ix,lxa) - d1(n2)*rneq
                  derv2(iy,lya) = derv2(iy,lya) - d1(n2)*rneq
                  derv2(iz,lza) = derv2(iz,lza) - d1(n2)*rneq
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if ((lopj.or.lopk).and.n1j.ne.n2k) then
              if (lopj.and.nauatom(n1j).gt.0) then
                rneq = dble(neqvatom(nauatom(n1j)))
                derv2(kx,jxa) = derv2(kx,jxa) - xd1*xd2*d2(ind)*rneq
                derv2(ky,jxa) = derv2(ky,jxa) - xd1*yd2*d2(ind)*rneq
                derv2(kz,jxa) = derv2(kz,jxa) - xd1*zd2*d2(ind)*rneq
                derv2(kx,jya) = derv2(kx,jya) - yd1*xd2*d2(ind)*rneq
                derv2(ky,jya) = derv2(ky,jya) - yd1*yd2*d2(ind)*rneq
                derv2(kz,jya) = derv2(kz,jya) - yd1*zd2*d2(ind)*rneq
                derv2(kx,jza) = derv2(kx,jza) - zd1*xd2*d2(ind)*rneq
                derv2(ky,jza) = derv2(ky,jza) - zd1*yd2*d2(ind)*rneq
                derv2(kz,jza) = derv2(kz,jza) - zd1*zd2*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(kx,jxa) = derv2(kx,jxa) - d1(n2)*rneq
                  derv2(ky,jya) = derv2(ky,jya) - d1(n2)*rneq
                  derv2(kz,jza) = derv2(kz,jza) - d1(n2)*rneq
                endif
              endif
              if (lopk.and.nauatom(n2k).gt.0.and..not.lsameijkl) then
                rneq = dble(neqvatom(nauatom(n2k)))
                derv2(jx,kxa) = derv2(jx,kxa) - xd2*xd1*d2(ind)*rneq
                derv2(jy,kxa) = derv2(jy,kxa) - xd2*yd1*d2(ind)*rneq
                derv2(jz,kxa) = derv2(jz,kxa) - xd2*zd1*d2(ind)*rneq
                derv2(jx,kya) = derv2(jx,kya) - yd2*xd1*d2(ind)*rneq
                derv2(jy,kya) = derv2(jy,kya) - yd2*yd1*d2(ind)*rneq
                derv2(jz,kya) = derv2(jz,kya) - yd2*zd1*d2(ind)*rneq
                derv2(jx,kza) = derv2(jx,kza) - zd2*xd1*d2(ind)*rneq
                derv2(jy,kza) = derv2(jy,kza) - zd2*yd1*d2(ind)*rneq
                derv2(jz,kza) = derv2(jz,kza) - zd2*zd1*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(jx,kxa) = derv2(jx,kxa) - d1(n2)*rneq
                  derv2(jy,kya) = derv2(jy,kya) - d1(n2)*rneq
                  derv2(jz,kza) = derv2(jz,kza) - d1(n2)*rneq
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if ((lopj.or.lopl).and.n1j.ne.n2l) then
              if (lopj.and.nauatom(n1j).gt.0) then
                rneq = dble(neqvatom(nauatom(n1j)))
                derv2(lx,jxa) = derv2(lx,jxa) + xd1*xd2*d2(ind)*rneq
                derv2(ly,jxa) = derv2(ly,jxa) + xd1*yd2*d2(ind)*rneq
                derv2(lz,jxa) = derv2(lz,jxa) + xd1*zd2*d2(ind)*rneq
                derv2(lx,jya) = derv2(lx,jya) + yd1*xd2*d2(ind)*rneq
                derv2(ly,jya) = derv2(ly,jya) + yd1*yd2*d2(ind)*rneq
                derv2(lz,jya) = derv2(lz,jya) + yd1*zd2*d2(ind)*rneq
                derv2(lx,jza) = derv2(lx,jza) + zd1*xd2*d2(ind)*rneq
                derv2(ly,jza) = derv2(ly,jza) + zd1*yd2*d2(ind)*rneq
                derv2(lz,jza) = derv2(lz,jza) + zd1*zd2*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(lx,jxa) = derv2(lx,jxa) + d1(n2)*rneq
                  derv2(ly,jya) = derv2(ly,jya) + d1(n2)*rneq
                  derv2(lz,jza) = derv2(lz,jza) + d1(n2)*rneq
                endif
              endif
              if (lopl.and.nauatom(n2l).gt.0.and..not.lsameijkl) then
                rneq = dble(neqvatom(nauatom(n2l)))
                derv2(jx,lxa) = derv2(jx,lxa) + xd2*xd1*d2(ind)*rneq
                derv2(jy,lxa) = derv2(jy,lxa) + xd2*yd1*d2(ind)*rneq
                derv2(jz,lxa) = derv2(jz,lxa) + xd2*zd1*d2(ind)*rneq
                derv2(jx,lya) = derv2(jx,lya) + yd2*xd1*d2(ind)*rneq
                derv2(jy,lya) = derv2(jy,lya) + yd2*yd1*d2(ind)*rneq
                derv2(jz,lya) = derv2(jz,lya) + yd2*zd1*d2(ind)*rneq
                derv2(jx,lza) = derv2(jx,lza) + zd2*xd1*d2(ind)*rneq
                derv2(jy,lza) = derv2(jy,lza) + zd2*yd1*d2(ind)*rneq
                derv2(jz,lza) = derv2(jz,lza) + zd2*zd1*d2(ind)*rneq
                if (n1.eq.n2) then
                  derv2(jx,lxa) = derv2(jx,lxa) + d1(n2)*rneq
                  derv2(jy,lya) = derv2(jy,lya) + d1(n2)*rneq
                  derv2(jz,lza) = derv2(jz,lza) + d1(n2)*rneq
                endif
              endif
            endif
          endif
        endif
      enddo
    else
!
!  Increment ind pointer by number of missed terms
!
      ind = ind + n1
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addds')
#endif
!
  return
  end
!
  subroutine d2addp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nREBOatomRptr,d1,d2,xkv,ykv,zkv,ltorderv,lextended)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian phased second derivative arrays.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  xkv             = x component of K vector
!  ykv             = y component of K vector
!  zkv             = z component of K vector
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!
!   8/02 Created from d2add
!   8/02 Removal of uninvalid terms added
!   9/02 Correction to lsameijkl added - need to check distances too
!  12/03 Torsional derivatives made optional
!   6/07 nREBOatomRptr added as input
!  10/07 Correction made to setting of n1i/n2k for torsional case
!  11/07 Unused variables removed
!   1/10 Bug fixed for neighbour of i -> neighbour of neighbour of i case
!  12/13 Condition for multiplication by half changed to fix error in 
!        self-image phasing during phonon dispersion. 
!  10/14 Group velocities added as an option
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   4/23 Extended derivatives added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use control,        only : lgroupvelocity
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  real(dp),    intent(in)             :: xkv
  real(dp),    intent(in)             :: ykv
  real(dp),    intent(in)             :: zkv
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nji
  integer(i4)                         :: nkj
  integer(i4)                         :: nlk
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: ntmp
  complex(dpc)                        :: cdk(3)
  real(dp)                            :: cosik
  real(dp)                            :: cosil
  real(dp)                            :: cosjk
  real(dp)                            :: cosjl
  real(dp)                            :: d2k
  real(dp)                            :: d2ks
  real(dp)                            :: oneik
  real(dp)                            :: oneil
  real(dp)                            :: onejk
  real(dp)                            :: onejl
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: sinik
  real(dp)                            :: sinil
  real(dp)                            :: sinjk
  real(dp)                            :: sinjl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i
  real(dp)                            :: z2i
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lsameijkl
#ifdef TRACE
  call trace_in('d2addp')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndisttor = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
    if (lextended) then
      ndist = ndisttor + nneigh(nri)*mneigh*mneigh
    else
      ndist = ndisttor
    endif
  else
    ndist = nneigh2
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = nneigh2 + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
!  Loop over pairs of interatomic distances
!
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!     
!  Check that neighbour numbers are valid
!         
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))       
      if (liok.and.ljok) then
!         
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    elseif (n1.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!         
!  Check that atoms are valid
!       
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        if (nREBOatomRptr(neighno(nji,nri)).eq.0) then
          ljok = .false.
        else
          ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
        endif
      else
        ljok = .false.
      endif
!  
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    else
!     
!  Extended derivatives
!     
      ntmp = n1 - ndisttor
      nji = (ntmp - 1)/(mneigh*mneigh) + 1
      ntmp = ntmp - (nji - 1)*mneigh*mneigh
      nkj = (ntmp - 1)/mneigh + 1
      ntmp = ntmp - (nkj - 1)*mneigh
      nlk = ntmp
!     
      nrj = neighno(nji,nri)
      if (nrj.eq.0) cycle
      if (nkj.gt.nneigh(nrj)) cycle
      nrk = neighno(nkj,nrj)
      if (nrk.eq.0) cycle 
      if (nlk.gt.nneigh(nrk)) cycle
      n1i = nrk
      n1j = neighno(nlk,nrk)
      if (n1j.eq.0) cycle
!  
!  Check that atoms are valid
! 
      liok = (nkj.le.nneigh(nrj))
      ljok = (nlk.le.nneigh(nrk))
!
      if (liok.and.ljok) then
        x2i = xneigh(nji,nri) + xneigh(nkj,nrj)
        y2i = yneigh(nji,nri) + yneigh(nkj,nrj)
        z2i = zneigh(nji,nri) + zneigh(nkj,nrj)
!
        xd1 = xneigh(nlk,nrk)
        yd1 = yneigh(nlk,nrk)
        zd1 = zneigh(nlk,nrk)
      endif
    endif
!  
!  Skip the following section as it is irrelevant if i and j are not valid
!           
    if (liok.and.ljok) then
!
!  Set second derivative location pointers
!
      ix = 3*(n1i - 1) + 1
      iy = ix + 1
      iz = iy + 1
      jx = 3*(n1j - 1) + 1
      jy = jx + 1
      jz = jy + 1
!
      ind = n1*(n1-1)/2
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!  
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!               
!  Calculate vector between atoms
!             
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          elseif (n2.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!               
!  Check that atoms are valid
!             
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              if (nREBOatomRptr(neighno(nji,nri)).eq.0) then
                llok = .false.
              else
                llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
              endif
            else
              llok = .false.
            endif
!               
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          else
!
!  Extended derivatives
!
            ntmp = n2 - ndisttor
            nji = (ntmp - 1)/(mneigh*mneigh) + 1
            ntmp = ntmp - (nji - 1)*mneigh*mneigh
            nkj = (ntmp - 1)/mneigh + 1
            ntmp = ntmp - (nkj - 1)*mneigh
            nlk = ntmp
!
            nrj = neighno(nji,nri)
            if (nrj.eq.0) cycle
            if (nkj.gt.nneigh(nrj)) cycle
            nrk = neighno(nkj,nrj)
            if (nrk.eq.0) cycle
            if (nlk.gt.nneigh(nrk)) cycle
            n2k = nrk
            n2l = neighno(nlk,nrk)
            if (n2l.eq.0) cycle
!
!  Check that atoms are valid
!
            lkok = (nkj.le.nneigh(nrj))
            llok = (nlk.le.nneigh(nrk))
!
            if (liok.and.ljok) then
              xik = xneigh(nji,nri) + xneigh(nkj,nrj) - x2i
              yik = yneigh(nji,nri) + yneigh(nkj,nrj) - y2i
              zik = zneigh(nji,nri) + zneigh(nkj,nrj) - z2i
!
              xd2 = xneigh(nlk,nrk)
              yd2 = yneigh(nlk,nrk)
              zd2 = zneigh(nlk,nrk)
            endif
          endif
!               
!  No point continuing beyond here unless k and l are valid
!                 
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xil = xik + xd2
            yil = yik + yd2
            zil = zik + zd2
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
            kx = 3*(n2k - 1) + 1
            ky = kx + 1
            kz = ky + 1
            lx = 3*(n2l - 1) + 1
            ly = lx + 1
            lz = ly + 1
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
!  Second derivatives : i - k
!
            if (n1i.eq.n2k) then
              oneik = 1.0_dp
            else
              oneik = 0.0_dp
            endif
            cosik = xkv*xik + ykv*yik + zkv*zik
            sinik = sin(cosik)
            cosik = cos(cosik) - oneik
!
            if (lsameijkl) then
              cosik = 0.5_dp*cosik
              sinik = 0.5_dp*sinik
            endif
            derv2(kx,ix) = derv2(kx,ix) + xd1*xd2*d2(ind)*cosik
            derv2(ky,ix) = derv2(ky,ix) + xd1*yd2*d2(ind)*cosik
            derv2(kz,ix) = derv2(kz,ix) + xd1*zd2*d2(ind)*cosik
            derv2(kx,iy) = derv2(kx,iy) + yd1*xd2*d2(ind)*cosik
            derv2(ky,iy) = derv2(ky,iy) + yd1*yd2*d2(ind)*cosik
            derv2(kz,iy) = derv2(kz,iy) + yd1*zd2*d2(ind)*cosik
            derv2(kx,iz) = derv2(kx,iz) + zd1*xd2*d2(ind)*cosik
            derv2(ky,iz) = derv2(ky,iz) + zd1*yd2*d2(ind)*cosik
            derv2(kz,iz) = derv2(kz,iz) + zd1*zd2*d2(ind)*cosik
!
            dervi(kx,ix) = dervi(kx,ix) + xd1*xd2*d2(ind)*sinik
            dervi(ky,ix) = dervi(ky,ix) + xd1*yd2*d2(ind)*sinik
            dervi(kz,ix) = dervi(kz,ix) + xd1*zd2*d2(ind)*sinik
            dervi(kx,iy) = dervi(kx,iy) + yd1*xd2*d2(ind)*sinik
            dervi(ky,iy) = dervi(ky,iy) + yd1*yd2*d2(ind)*sinik
            dervi(kz,iy) = dervi(kz,iy) + yd1*zd2*d2(ind)*sinik
            dervi(kx,iz) = dervi(kx,iz) + zd1*xd2*d2(ind)*sinik
            dervi(ky,iz) = dervi(ky,iz) + zd1*yd2*d2(ind)*sinik
            dervi(kz,iz) = dervi(kz,iz) + zd1*zd2*d2(ind)*sinik
!
            if (lgroupvelocity) then
!
!  Group velocities
!
              d2k  = d2(ind)*cosik
              d2ks = d2(ind)*sinik
              cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
              derv2dk(1:3,kx,ix) = derv2dk(1:3,kx,ix) + xd1*xd2*cdk(1:3)
              derv2dk(1:3,ky,ix) = derv2dk(1:3,ky,ix) + xd1*yd2*cdk(1:3)
              derv2dk(1:3,kz,ix) = derv2dk(1:3,kz,ix) + xd1*zd2*cdk(1:3)
              derv2dk(1:3,kx,iy) = derv2dk(1:3,kx,iy) + yd1*xd2*cdk(1:3)
              derv2dk(1:3,ky,iy) = derv2dk(1:3,ky,iy) + yd1*yd2*cdk(1:3)
              derv2dk(1:3,kz,iy) = derv2dk(1:3,kz,iy) + yd1*zd2*cdk(1:3)
              derv2dk(1:3,kx,iz) = derv2dk(1:3,kx,iz) + zd1*xd2*cdk(1:3)
              derv2dk(1:3,ky,iz) = derv2dk(1:3,ky,iz) + zd1*yd2*cdk(1:3)
              derv2dk(1:3,kz,iz) = derv2dk(1:3,kz,iz) + zd1*zd2*cdk(1:3)
            endif
!
            if (n1.eq.n2) then
              derv2(kx,ix) = derv2(kx,ix) + d1(n2)*cosik
              derv2(ky,iy) = derv2(ky,iy) + d1(n2)*cosik
              derv2(kz,iz) = derv2(kz,iz) + d1(n2)*cosik
              dervi(kx,ix) = dervi(kx,ix) + d1(n2)*sinik
              dervi(ky,iy) = dervi(ky,iy) + d1(n2)*sinik
              dervi(kz,iz) = dervi(kz,iz) + d1(n2)*sinik
!
              if (lgroupvelocity) then
!
!  Group velocities
!
                d2k  = d1(n2)*cosik
                d2ks = d1(n2)*sinik
                cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                derv2dk(1:3,kx,ix) = derv2dk(1:3,kx,ix) + cdk(1:3)
                derv2dk(1:3,ky,iy) = derv2dk(1:3,ky,iy) + cdk(1:3)
                derv2dk(1:3,kz,iz) = derv2dk(1:3,kz,iz) + cdk(1:3)
              endif
            endif
!
            derv2(ix,kx) = derv2(ix,kx) + xd1*xd2*d2(ind)*cosik
            derv2(ix,ky) = derv2(ix,ky) + xd1*yd2*d2(ind)*cosik
            derv2(ix,kz) = derv2(ix,kz) + xd1*zd2*d2(ind)*cosik
            derv2(iy,kx) = derv2(iy,kx) + yd1*xd2*d2(ind)*cosik
            derv2(iy,ky) = derv2(iy,ky) + yd1*yd2*d2(ind)*cosik
            derv2(iy,kz) = derv2(iy,kz) + yd1*zd2*d2(ind)*cosik
            derv2(iz,kx) = derv2(iz,kx) + zd1*xd2*d2(ind)*cosik
            derv2(iz,ky) = derv2(iz,ky) + zd1*yd2*d2(ind)*cosik
            derv2(iz,kz) = derv2(iz,kz) + zd1*zd2*d2(ind)*cosik
!
            dervi(ix,kx) = dervi(ix,kx) - xd1*xd2*d2(ind)*sinik
            dervi(ix,ky) = dervi(ix,ky) - xd1*yd2*d2(ind)*sinik
            dervi(ix,kz) = dervi(ix,kz) - xd1*zd2*d2(ind)*sinik
            dervi(iy,kx) = dervi(iy,kx) - yd1*xd2*d2(ind)*sinik
            dervi(iy,ky) = dervi(iy,ky) - yd1*yd2*d2(ind)*sinik
            dervi(iy,kz) = dervi(iy,kz) - yd1*zd2*d2(ind)*sinik
            dervi(iz,kx) = dervi(iz,kx) - zd1*xd2*d2(ind)*sinik
            dervi(iz,ky) = dervi(iz,ky) - zd1*yd2*d2(ind)*sinik
            dervi(iz,kz) = dervi(iz,kz) - zd1*zd2*d2(ind)*sinik
!
            if (lgroupvelocity) then
!
!  Group velocities
!
              d2k  = d2(ind)*cosik
              d2ks = d2(ind)*sinik
              cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
              derv2dk(1:3,ix,kx) = derv2dk(1:3,ix,kx) + xd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,ix,ky) = derv2dk(1:3,ix,ky) + xd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,ix,kz) = derv2dk(1:3,ix,kz) + xd1*zd2*conjg(cdk(1:3))
              derv2dk(1:3,iy,kx) = derv2dk(1:3,iy,kx) + yd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,iy,ky) = derv2dk(1:3,iy,ky) + yd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,iy,kz) = derv2dk(1:3,iy,kz) + yd1*zd2*conjg(cdk(1:3))
              derv2dk(1:3,iz,kx) = derv2dk(1:3,iz,kx) + zd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,iz,ky) = derv2dk(1:3,iz,ky) + zd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,iz,kz) = derv2dk(1:3,iz,kz) + zd1*zd2*conjg(cdk(1:3))
            endif
!
            if (n1.eq.n2) then
              derv2(ix,kx) = derv2(ix,kx) + d1(n2)*cosik
              derv2(iy,ky) = derv2(iy,ky) + d1(n2)*cosik
              derv2(iz,kz) = derv2(iz,kz) + d1(n2)*cosik
              dervi(ix,kx) = dervi(ix,kx) - d1(n2)*sinik
              dervi(iy,ky) = dervi(iy,ky) - d1(n2)*sinik
              dervi(iz,kz) = dervi(iz,kz) - d1(n2)*sinik
!
              if (lgroupvelocity) then
!
!  Group velocities
!
                d2k  = d1(n2)*cosik
                d2ks = d1(n2)*sinik
                cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                derv2dk(1:3,ix,kx) = derv2dk(1:3,ix,kx) + conjg(cdk(1:3))
                derv2dk(1:3,iy,ky) = derv2dk(1:3,iy,ky) + conjg(cdk(1:3))
                derv2dk(1:3,iz,kz) = derv2dk(1:3,iz,kz) + conjg(cdk(1:3))
              endif
            endif
!
!  Second derivatives : i - l
!
            if (n1i.eq.n2l) then
              oneil = 1.0_dp
            else
              oneil = 0.0_dp
            endif
            cosil = xkv*xil + ykv*yil + zkv*zil
            sinil = sin(cosil)
            cosil = cos(cosil) - oneil
!
            if (lsameijkl) then
              cosil = 0.5_dp*cosil
              sinil = 0.5_dp*sinil
            endif
!
            derv2(lx,ix) = derv2(lx,ix) - xd1*xd2*d2(ind)*cosil
            derv2(ly,ix) = derv2(ly,ix) - xd1*yd2*d2(ind)*cosil
            derv2(lz,ix) = derv2(lz,ix) - xd1*zd2*d2(ind)*cosil
            derv2(lx,iy) = derv2(lx,iy) - yd1*xd2*d2(ind)*cosil
            derv2(ly,iy) = derv2(ly,iy) - yd1*yd2*d2(ind)*cosil
            derv2(lz,iy) = derv2(lz,iy) - yd1*zd2*d2(ind)*cosil
            derv2(lx,iz) = derv2(lx,iz) - zd1*xd2*d2(ind)*cosil
            derv2(ly,iz) = derv2(ly,iz) - zd1*yd2*d2(ind)*cosil
            derv2(lz,iz) = derv2(lz,iz) - zd1*zd2*d2(ind)*cosil
!
            dervi(lx,ix) = dervi(lx,ix) - xd1*xd2*d2(ind)*sinil
            dervi(ly,ix) = dervi(ly,ix) - xd1*yd2*d2(ind)*sinil
            dervi(lz,ix) = dervi(lz,ix) - xd1*zd2*d2(ind)*sinil
            dervi(lx,iy) = dervi(lx,iy) - yd1*xd2*d2(ind)*sinil
            dervi(ly,iy) = dervi(ly,iy) - yd1*yd2*d2(ind)*sinil
            dervi(lz,iy) = dervi(lz,iy) - yd1*zd2*d2(ind)*sinil
            dervi(lx,iz) = dervi(lx,iz) - zd1*xd2*d2(ind)*sinil
            dervi(ly,iz) = dervi(ly,iz) - zd1*yd2*d2(ind)*sinil
            dervi(lz,iz) = dervi(lz,iz) - zd1*zd2*d2(ind)*sinil
!
            if (lgroupvelocity) then
!
!  Group velocities
!
              d2k  = d2(ind)*cosil
              d2ks = d2(ind)*sinil
              cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
              derv2dk(1:3,lx,ix) = derv2dk(1:3,lx,ix) - xd1*xd2*cdk(1:3)
              derv2dk(1:3,ly,ix) = derv2dk(1:3,ly,ix) - xd1*yd2*cdk(1:3)
              derv2dk(1:3,lz,ix) = derv2dk(1:3,lz,ix) - xd1*zd2*cdk(1:3)
              derv2dk(1:3,lx,iy) = derv2dk(1:3,lx,iy) - yd1*xd2*cdk(1:3)
              derv2dk(1:3,ly,iy) = derv2dk(1:3,ly,iy) - yd1*yd2*cdk(1:3)
              derv2dk(1:3,lz,iy) = derv2dk(1:3,lz,iy) - yd1*zd2*cdk(1:3)
              derv2dk(1:3,lx,iz) = derv2dk(1:3,lx,iz) - zd1*xd2*cdk(1:3)
              derv2dk(1:3,ly,iz) = derv2dk(1:3,ly,iz) - zd1*yd2*cdk(1:3)
              derv2dk(1:3,lz,iz) = derv2dk(1:3,lz,iz) - zd1*zd2*cdk(1:3)
            endif
!
            if (n1.eq.n2) then
              derv2(lx,ix) = derv2(lx,ix) - d1(n2)*cosil
              derv2(ly,iy) = derv2(ly,iy) - d1(n2)*cosil
              derv2(lz,iz) = derv2(lz,iz) - d1(n2)*cosil
              dervi(lx,ix) = dervi(lx,ix) - d1(n2)*sinil
              dervi(ly,iy) = dervi(ly,iy) - d1(n2)*sinil
              dervi(lz,iz) = dervi(lz,iz) - d1(n2)*sinil
!
              if (lgroupvelocity) then
!
!  Group velocities
!
                d2k  = d1(n2)*cosil
                d2ks = d1(n2)*sinil
                cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                derv2dk(1:3,lx,ix) = derv2dk(1:3,lx,ix) - cdk(1:3)
                derv2dk(1:3,ly,iy) = derv2dk(1:3,ly,iy) - cdk(1:3)
                derv2dk(1:3,lz,iz) = derv2dk(1:3,lz,iz) - cdk(1:3)
              endif
            endif
!
            derv2(ix,lx) = derv2(ix,lx) - xd1*xd2*d2(ind)*cosil
            derv2(ix,ly) = derv2(ix,ly) - xd1*yd2*d2(ind)*cosil
            derv2(ix,lz) = derv2(ix,lz) - xd1*zd2*d2(ind)*cosil
            derv2(iy,lx) = derv2(iy,lx) - yd1*xd2*d2(ind)*cosil
            derv2(iy,ly) = derv2(iy,ly) - yd1*yd2*d2(ind)*cosil
            derv2(iy,lz) = derv2(iy,lz) - yd1*zd2*d2(ind)*cosil
            derv2(iz,lx) = derv2(iz,lx) - zd1*xd2*d2(ind)*cosil
            derv2(iz,ly) = derv2(iz,ly) - zd1*yd2*d2(ind)*cosil
            derv2(iz,lz) = derv2(iz,lz) - zd1*zd2*d2(ind)*cosil
!
            dervi(ix,lx) = dervi(ix,lx) + xd1*xd2*d2(ind)*sinil
            dervi(ix,ly) = dervi(ix,ly) + xd1*yd2*d2(ind)*sinil
            dervi(ix,lz) = dervi(ix,lz) + xd1*zd2*d2(ind)*sinil
            dervi(iy,lx) = dervi(iy,lx) + yd1*xd2*d2(ind)*sinil
            dervi(iy,ly) = dervi(iy,ly) + yd1*yd2*d2(ind)*sinil
            dervi(iy,lz) = dervi(iy,lz) + yd1*zd2*d2(ind)*sinil
            dervi(iz,lx) = dervi(iz,lx) + zd1*xd2*d2(ind)*sinil
            dervi(iz,ly) = dervi(iz,ly) + zd1*yd2*d2(ind)*sinil
            dervi(iz,lz) = dervi(iz,lz) + zd1*zd2*d2(ind)*sinil
!
            if (lgroupvelocity) then
!
!  Group velocities
!
              d2k  = d2(ind)*cosil
              d2ks = d2(ind)*sinil
              cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
              derv2dk(1:3,ix,lx) = derv2dk(1:3,ix,lx) - xd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,ix,ly) = derv2dk(1:3,ix,ly) - xd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,ix,lz) = derv2dk(1:3,ix,lz) - xd1*zd2*conjg(cdk(1:3))
              derv2dk(1:3,iy,lx) = derv2dk(1:3,iy,lx) - yd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,iy,ly) = derv2dk(1:3,iy,ly) - yd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,iy,lz) = derv2dk(1:3,iy,lz) - yd1*zd2*conjg(cdk(1:3))
              derv2dk(1:3,iz,lx) = derv2dk(1:3,iz,lx) - zd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,iz,ly) = derv2dk(1:3,iz,ly) - zd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,iz,lz) = derv2dk(1:3,iz,lz) - zd1*zd2*conjg(cdk(1:3))
            endif
!
            if (n1.eq.n2) then
              derv2(ix,lx) = derv2(ix,lx) - d1(n2)*cosil
              derv2(iy,ly) = derv2(iy,ly) - d1(n2)*cosil
              derv2(iz,lz) = derv2(iz,lz) - d1(n2)*cosil
              dervi(ix,lx) = dervi(ix,lx) + d1(n2)*sinil
              dervi(iy,ly) = dervi(iy,ly) + d1(n2)*sinil
              dervi(iz,lz) = dervi(iz,lz) + d1(n2)*sinil
!
              if (lgroupvelocity) then
!
!  Group velocities
!
                d2k  = d1(n2)*cosil
                d2ks = d1(n2)*sinil
                cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                derv2dk(1:3,ix,lx) = derv2dk(1:3,ix,lx) - conjg(cdk(1:3))
                derv2dk(1:3,iy,ly) = derv2dk(1:3,iy,ly) - conjg(cdk(1:3))
                derv2dk(1:3,iz,lz) = derv2dk(1:3,iz,lz) - conjg(cdk(1:3))
              endif
            endif
!
!  Second derivatives : j - k
!
            if (n1j.eq.n2k) then
              onejk = 1.0_dp
            else
              onejk = 0.0_dp
            endif
            cosjk = xkv*xjk + ykv*yjk + zkv*zjk
            sinjk = sin(cosjk)
            cosjk = cos(cosjk) - onejk
!
            if (lsameijkl) then
              cosjk = 0.5_dp*cosjk
              sinjk = 0.5_dp*sinjk
            endif
!
            derv2(kx,jx) = derv2(kx,jx) - xd1*xd2*d2(ind)*cosjk
            derv2(ky,jx) = derv2(ky,jx) - xd1*yd2*d2(ind)*cosjk
            derv2(kz,jx) = derv2(kz,jx) - xd1*zd2*d2(ind)*cosjk
            derv2(kx,jy) = derv2(kx,jy) - yd1*xd2*d2(ind)*cosjk
            derv2(ky,jy) = derv2(ky,jy) - yd1*yd2*d2(ind)*cosjk
            derv2(kz,jy) = derv2(kz,jy) - yd1*zd2*d2(ind)*cosjk
            derv2(kx,jz) = derv2(kx,jz) - zd1*xd2*d2(ind)*cosjk
            derv2(ky,jz) = derv2(ky,jz) - zd1*yd2*d2(ind)*cosjk
            derv2(kz,jz) = derv2(kz,jz) - zd1*zd2*d2(ind)*cosjk
!
            dervi(kx,jx) = dervi(kx,jx) - xd1*xd2*d2(ind)*sinjk
            dervi(ky,jx) = dervi(ky,jx) - xd1*yd2*d2(ind)*sinjk
            dervi(kz,jx) = dervi(kz,jx) - xd1*zd2*d2(ind)*sinjk
            dervi(kx,jy) = dervi(kx,jy) - yd1*xd2*d2(ind)*sinjk
            dervi(ky,jy) = dervi(ky,jy) - yd1*yd2*d2(ind)*sinjk
            dervi(kz,jy) = dervi(kz,jy) - yd1*zd2*d2(ind)*sinjk
            dervi(kx,jz) = dervi(kx,jz) - zd1*xd2*d2(ind)*sinjk
            dervi(ky,jz) = dervi(ky,jz) - zd1*yd2*d2(ind)*sinjk
            dervi(kz,jz) = dervi(kz,jz) - zd1*zd2*d2(ind)*sinjk
!
            if (lgroupvelocity) then
!
!  Group velocities
!
              d2k  = d2(ind)*cosjk
              d2ks = d2(ind)*sinjk
              cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
              derv2dk(1:3,kx,jx) = derv2dk(1:3,kx,jx) - xd1*xd2*cdk(1:3)
              derv2dk(1:3,ky,jx) = derv2dk(1:3,ky,jx) - xd1*yd2*cdk(1:3)
              derv2dk(1:3,kz,jx) = derv2dk(1:3,kz,jx) - xd1*zd2*cdk(1:3)
              derv2dk(1:3,kx,jy) = derv2dk(1:3,kx,jy) - yd1*xd2*cdk(1:3)
              derv2dk(1:3,ky,jy) = derv2dk(1:3,ky,jy) - yd1*yd2*cdk(1:3)
              derv2dk(1:3,kz,jy) = derv2dk(1:3,kz,jy) - yd1*zd2*cdk(1:3)
              derv2dk(1:3,kx,jz) = derv2dk(1:3,kx,jz) - zd1*xd2*cdk(1:3)
              derv2dk(1:3,ky,jz) = derv2dk(1:3,ky,jz) - zd1*yd2*cdk(1:3)
              derv2dk(1:3,kz,jz) = derv2dk(1:3,kz,jz) - zd1*zd2*cdk(1:3)
            endif
!
            if (n1.eq.n2) then
              derv2(kx,jx) = derv2(kx,jx) - d1(n2)*cosjk
              derv2(ky,jy) = derv2(ky,jy) - d1(n2)*cosjk
              derv2(kz,jz) = derv2(kz,jz) - d1(n2)*cosjk
              dervi(kx,jx) = dervi(kx,jx) - d1(n2)*sinjk
              dervi(ky,jy) = dervi(ky,jy) - d1(n2)*sinjk
              dervi(kz,jz) = dervi(kz,jz) - d1(n2)*sinjk
!
              if (lgroupvelocity) then
!
!  Group velocities
!
                d2k  = d1(n2)*cosjk
                d2ks = d1(n2)*sinjk
                cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                derv2dk(1:3,kx,jx) = derv2dk(1:3,kx,jx) - cdk(1:3)
                derv2dk(1:3,ky,jy) = derv2dk(1:3,ky,jy) - cdk(1:3)
                derv2dk(1:3,kz,jz) = derv2dk(1:3,kz,jz) - cdk(1:3)
              endif
            endif
!
            derv2(jx,kx) = derv2(jx,kx) - xd1*xd2*d2(ind)*cosjk
            derv2(jx,ky) = derv2(jx,ky) - xd1*yd2*d2(ind)*cosjk
            derv2(jx,kz) = derv2(jx,kz) - xd1*zd2*d2(ind)*cosjk
            derv2(jy,kx) = derv2(jy,kx) - yd1*xd2*d2(ind)*cosjk
            derv2(jy,ky) = derv2(jy,ky) - yd1*yd2*d2(ind)*cosjk
            derv2(jy,kz) = derv2(jy,kz) - yd1*zd2*d2(ind)*cosjk
            derv2(jz,kx) = derv2(jz,kx) - zd1*xd2*d2(ind)*cosjk
            derv2(jz,ky) = derv2(jz,ky) - zd1*yd2*d2(ind)*cosjk
            derv2(jz,kz) = derv2(jz,kz) - zd1*zd2*d2(ind)*cosjk
!
            dervi(jx,kx) = dervi(jx,kx) + xd1*xd2*d2(ind)*sinjk
            dervi(jx,ky) = dervi(jx,ky) + xd1*yd2*d2(ind)*sinjk
            dervi(jx,kz) = dervi(jx,kz) + xd1*zd2*d2(ind)*sinjk
            dervi(jy,kx) = dervi(jy,kx) + yd1*xd2*d2(ind)*sinjk
            dervi(jy,ky) = dervi(jy,ky) + yd1*yd2*d2(ind)*sinjk
            dervi(jy,kz) = dervi(jy,kz) + yd1*zd2*d2(ind)*sinjk
            dervi(jz,kx) = dervi(jz,kx) + zd1*xd2*d2(ind)*sinjk
            dervi(jz,ky) = dervi(jz,ky) + zd1*yd2*d2(ind)*sinjk
            dervi(jz,kz) = dervi(jz,kz) + zd1*zd2*d2(ind)*sinjk
!
            if (lgroupvelocity) then
!
!  Group velocities
!
              d2k  = d2(ind)*cosjk
              d2ks = d2(ind)*sinjk
              cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
              derv2dk(1:3,jx,kx) = derv2dk(1:3,jx,kx) - xd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,jx,ky) = derv2dk(1:3,jx,ky) - xd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,jx,kz) = derv2dk(1:3,jx,kz) - xd1*zd2*conjg(cdk(1:3))
              derv2dk(1:3,jy,kx) = derv2dk(1:3,jy,kx) - yd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,jy,ky) = derv2dk(1:3,jy,ky) - yd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,jy,kz) = derv2dk(1:3,jy,kz) - yd1*zd2*conjg(cdk(1:3))
              derv2dk(1:3,jz,kx) = derv2dk(1:3,jz,kx) - zd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,jz,ky) = derv2dk(1:3,jz,ky) - zd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,jz,kz) = derv2dk(1:3,jz,kz) - zd1*zd2*conjg(cdk(1:3))
            endif
!
            if (n1.eq.n2) then
              derv2(jx,kx) = derv2(jx,kx) - d1(n2)*cosjk
              derv2(jy,ky) = derv2(jy,ky) - d1(n2)*cosjk
              derv2(jz,kz) = derv2(jz,kz) - d1(n2)*cosjk
              dervi(jx,kx) = dervi(jx,kx) + d1(n2)*sinjk
              dervi(jy,ky) = dervi(jy,ky) + d1(n2)*sinjk
              dervi(jz,kz) = dervi(jz,kz) + d1(n2)*sinjk
!
              if (lgroupvelocity) then
!
!  Group velocities
!
                d2k  = d1(n2)*cosjk
                d2ks = d1(n2)*sinjk
                cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                derv2dk(1:3,jx,kx) = derv2dk(1:3,jx,kx) - conjg(cdk(1:3))
                derv2dk(1:3,jy,ky) = derv2dk(1:3,jy,ky) - conjg(cdk(1:3))
                derv2dk(1:3,jz,kz) = derv2dk(1:3,jz,kz) - conjg(cdk(1:3))
              endif
            endif
!
!  Second derivatives : j - l
!
            if (n1j.eq.n2l) then
              onejl = 1.0_dp
            else
              onejl = 0.0_dp
            endif
            cosjl = xkv*xjl + ykv*yjl + zkv*zjl
            sinjl = sin(cosjl)
            cosjl = cos(cosjl) - onejl
!
            if (lsameijkl) then
              cosjl = 0.5_dp*cosjl
              sinjl = 0.5_dp*sinjl
            endif
!
            derv2(lx,jx) = derv2(lx,jx) + xd1*xd2*d2(ind)*cosjl
            derv2(ly,jx) = derv2(ly,jx) + xd1*yd2*d2(ind)*cosjl
            derv2(lz,jx) = derv2(lz,jx) + xd1*zd2*d2(ind)*cosjl
            derv2(lx,jy) = derv2(lx,jy) + yd1*xd2*d2(ind)*cosjl
            derv2(ly,jy) = derv2(ly,jy) + yd1*yd2*d2(ind)*cosjl
            derv2(lz,jy) = derv2(lz,jy) + yd1*zd2*d2(ind)*cosjl
            derv2(lx,jz) = derv2(lx,jz) + zd1*xd2*d2(ind)*cosjl
            derv2(ly,jz) = derv2(ly,jz) + zd1*yd2*d2(ind)*cosjl
            derv2(lz,jz) = derv2(lz,jz) + zd1*zd2*d2(ind)*cosjl
!
            dervi(lx,jx) = dervi(lx,jx) + xd1*xd2*d2(ind)*sinjl
            dervi(ly,jx) = dervi(ly,jx) + xd1*yd2*d2(ind)*sinjl
            dervi(lz,jx) = dervi(lz,jx) + xd1*zd2*d2(ind)*sinjl
            dervi(lx,jy) = dervi(lx,jy) + yd1*xd2*d2(ind)*sinjl
            dervi(ly,jy) = dervi(ly,jy) + yd1*yd2*d2(ind)*sinjl
            dervi(lz,jy) = dervi(lz,jy) + yd1*zd2*d2(ind)*sinjl
            dervi(lx,jz) = dervi(lx,jz) + zd1*xd2*d2(ind)*sinjl
            dervi(ly,jz) = dervi(ly,jz) + zd1*yd2*d2(ind)*sinjl
            dervi(lz,jz) = dervi(lz,jz) + zd1*zd2*d2(ind)*sinjl
!
            if (lgroupvelocity) then
!
!  Group velocities
!
              d2k  = d2(ind)*cosjl
              d2ks = d2(ind)*sinjl
              cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!
              derv2dk(1:3,lx,jx) = derv2dk(1:3,lx,jx) + xd1*xd2*cdk(1:3)
              derv2dk(1:3,ly,jx) = derv2dk(1:3,ly,jx) + xd1*yd2*cdk(1:3)
              derv2dk(1:3,lz,jx) = derv2dk(1:3,lz,jx) + xd1*zd2*cdk(1:3)
              derv2dk(1:3,lx,jy) = derv2dk(1:3,lx,jy) + yd1*xd2*cdk(1:3)
              derv2dk(1:3,ly,jy) = derv2dk(1:3,ly,jy) + yd1*yd2*cdk(1:3)
              derv2dk(1:3,lz,jy) = derv2dk(1:3,lz,jy) + yd1*zd2*cdk(1:3)
              derv2dk(1:3,lx,jz) = derv2dk(1:3,lx,jz) + zd1*xd2*cdk(1:3)
              derv2dk(1:3,ly,jz) = derv2dk(1:3,ly,jz) + zd1*yd2*cdk(1:3)
              derv2dk(1:3,lz,jz) = derv2dk(1:3,lz,jz) + zd1*zd2*cdk(1:3)
            endif
!
            if (n1.eq.n2) then
              derv2(lx,jx) = derv2(lx,jx) + d1(n2)*cosjl
              derv2(ly,jy) = derv2(ly,jy) + d1(n2)*cosjl
              derv2(lz,jz) = derv2(lz,jz) + d1(n2)*cosjl
              dervi(lx,jx) = dervi(lx,jx) + d1(n2)*sinjl
              dervi(ly,jy) = dervi(ly,jy) + d1(n2)*sinjl
              dervi(lz,jz) = dervi(lz,jz) + d1(n2)*sinjl
!
              if (lgroupvelocity) then
!
!  Group velocities
!
                d2k  = d1(n2)*cosjl
                d2ks = d1(n2)*sinjl
                cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!
                derv2dk(1:3,lx,jx) = derv2dk(1:3,lx,jx) + cdk(1:3)
                derv2dk(1:3,ly,jy) = derv2dk(1:3,ly,jy) + cdk(1:3)
                derv2dk(1:3,lz,jz) = derv2dk(1:3,lz,jz) + cdk(1:3)
              endif
            endif
!
            derv2(jx,lx) = derv2(jx,lx) + xd1*xd2*d2(ind)*cosjl
            derv2(jx,ly) = derv2(jx,ly) + xd1*yd2*d2(ind)*cosjl
            derv2(jx,lz) = derv2(jx,lz) + xd1*zd2*d2(ind)*cosjl
            derv2(jy,lx) = derv2(jy,lx) + yd1*xd2*d2(ind)*cosjl
            derv2(jy,ly) = derv2(jy,ly) + yd1*yd2*d2(ind)*cosjl
            derv2(jy,lz) = derv2(jy,lz) + yd1*zd2*d2(ind)*cosjl
            derv2(jz,lx) = derv2(jz,lx) + zd1*xd2*d2(ind)*cosjl
            derv2(jz,ly) = derv2(jz,ly) + zd1*yd2*d2(ind)*cosjl
            derv2(jz,lz) = derv2(jz,lz) + zd1*zd2*d2(ind)*cosjl
!
            dervi(jx,lx) = dervi(jx,lx) - xd1*xd2*d2(ind)*sinjl
            dervi(jx,ly) = dervi(jx,ly) - xd1*yd2*d2(ind)*sinjl
            dervi(jx,lz) = dervi(jx,lz) - xd1*zd2*d2(ind)*sinjl
            dervi(jy,lx) = dervi(jy,lx) - yd1*xd2*d2(ind)*sinjl
            dervi(jy,ly) = dervi(jy,ly) - yd1*yd2*d2(ind)*sinjl
            dervi(jy,lz) = dervi(jy,lz) - yd1*zd2*d2(ind)*sinjl
            dervi(jz,lx) = dervi(jz,lx) - zd1*xd2*d2(ind)*sinjl
            dervi(jz,ly) = dervi(jz,ly) - zd1*yd2*d2(ind)*sinjl
            dervi(jz,lz) = dervi(jz,lz) - zd1*zd2*d2(ind)*sinjl
!
            if (lgroupvelocity) then
!
!  Group velocities
!
              d2k  = d2(ind)*cosjl
              d2ks = d2(ind)*sinjl
              cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
              cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
              cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!
              derv2dk(1:3,jx,lx) = derv2dk(1:3,jx,lx) + xd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,jx,ly) = derv2dk(1:3,jx,ly) + xd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,jx,lz) = derv2dk(1:3,jx,lz) + xd1*zd2*conjg(cdk(1:3))
              derv2dk(1:3,jy,lx) = derv2dk(1:3,jy,lx) + yd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,jy,ly) = derv2dk(1:3,jy,ly) + yd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,jy,lz) = derv2dk(1:3,jy,lz) + yd1*zd2*conjg(cdk(1:3))
              derv2dk(1:3,jz,lx) = derv2dk(1:3,jz,lx) + zd1*xd2*conjg(cdk(1:3))
              derv2dk(1:3,jz,ly) = derv2dk(1:3,jz,ly) + zd1*yd2*conjg(cdk(1:3))
              derv2dk(1:3,jz,lz) = derv2dk(1:3,jz,lz) + zd1*zd2*conjg(cdk(1:3))
            endif
!
            if (n1.eq.n2) then
              derv2(jx,lx) = derv2(jx,lx) + d1(n2)*cosjl
              derv2(jy,ly) = derv2(jy,ly) + d1(n2)*cosjl
              derv2(jz,lz) = derv2(jz,lz) + d1(n2)*cosjl
              dervi(jx,lx) = dervi(jx,lx) - d1(n2)*sinjl
              dervi(jy,ly) = dervi(jy,ly) - d1(n2)*sinjl
              dervi(jz,lz) = dervi(jz,lz) - d1(n2)*sinjl
!
              if (lgroupvelocity) then
!
!  Group velocities
!           
                d2k  = d1(n2)*cosjl
                d2ks = d1(n2)*sinjl
                cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!           
                derv2dk(1:3,jx,lx) = derv2dk(1:3,jx,lx) + conjg(cdk(1:3))
                derv2dk(1:3,jy,ly) = derv2dk(1:3,jy,ly) + conjg(cdk(1:3))
                derv2dk(1:3,jz,lz) = derv2dk(1:3,jz,lz) + conjg(cdk(1:3))
              endif
            endif
          endif
        endif
      enddo
    else    
!               
!  Increment ind pointer by number of missed terms
!  
      ind = ind + n1
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addp')
#endif
!
  return
  end
!
  subroutine d2addfc(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,nREBOatomRptr,d1,d2,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian unphased second derivative arrays.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  ineigh          = cell indices for vector to neighbour of atoms
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!
!  11/14 Created from d2addp
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: ineigh(3_i4,maxneigh,*)
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
!
!  Local variables
!
  integer(i4)                         :: id11
  integer(i4)                         :: id12
  integer(i4)                         :: id13
  integer(i4)                         :: id21
  integer(i4)                         :: id22
  integer(i4)                         :: id23
  integer(i4)                         :: i21
  integer(i4)                         :: i22
  integer(i4)                         :: i23
  integer(i4)                         :: iik1
  integer(i4)                         :: iik2
  integer(i4)                         :: iik3
  integer(i4)                         :: iil1
  integer(i4)                         :: iil2
  integer(i4)                         :: iil3
  integer(i4)                         :: ijk1
  integer(i4)                         :: ijk2
  integer(i4)                         :: ijk3
  integer(i4)                         :: ijl1
  integer(i4)                         :: ijl2
  integer(i4)                         :: ijl3
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ncindikm
  integer(i4)                         :: ncindikp
  integer(i4)                         :: ncindilm
  integer(i4)                         :: ncindilp
  integer(i4)                         :: ncindjkm
  integer(i4)                         :: ncindjkp
  integer(i4)                         :: ncindjlm
  integer(i4)                         :: ncindjlp
  integer(i4)                         :: ndist
  integer(i4)                         :: nji
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: ntmp
  real(dp)                            :: oneik
  real(dp)                            :: oneil
  real(dp)                            :: onejk
  real(dp)                            :: onejl
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i
  real(dp)                            :: z2i
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lsameijkl
#ifdef TRACE
  call trace_in('d2addfc')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndist = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
  else
    ndist = nneigh2
  endif
!
!  Loop over pairs of interatomic distances
!
  ind = 0
!
!  id11 / id12 / id13 are the cell indices between i and j
!
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      id11 = ineigh(1,n1,nri)
      id12 = ineigh(2,n1,nri)
      id13 = ineigh(3,n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      i21 = 0
      i22 = 0
      i23 = 0
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!     
!  Check that neighbour numbers are valid
!         
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))       
      if (liok.and.ljok) then
!         
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
        id11 = (ineigh(1,n1j,nri) - ineigh(1,n1i,nri))
        id12 = (ineigh(2,n1j,nri) - ineigh(2,n1i,nri))
        id13 = (ineigh(3,n1j,nri) - ineigh(3,n1i,nri))
!
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
        i21 = ineigh(1,n1i,nri)
        i22 = ineigh(2,n1i,nri)
        i23 = ineigh(3,n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    else
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!         
!  Check that atoms are valid
!       
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
!  
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
          id11 = (ineigh(1,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(1,nji,nri))
          id12 = (ineigh(2,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(2,nji,nri))
          id13 = (ineigh(3,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(3,nji,nri))
          i21 = 0
          i22 = 0
          i23 = 0
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          id11 = (ineigh(1,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(1,nji,nri) - ineigh(1,n1i,nri))
          id12 = (ineigh(2,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(2,nji,nri) - ineigh(2,n1i,nri))
          id13 = (ineigh(3,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(3,nji,nri) - ineigh(3,n1i,nri))
          i21 = ineigh(1,n1i,nri)
          i22 = ineigh(2,n1i,nri)
          i23 = ineigh(3,n1i,nri)
!
!  Convert neighbour number to real atom
!
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    endif
!  
!  Skip the following section as it is irrelevant if i and j are not valid
!           
    if (liok.and.ljok) then
!
!  Set second derivative location pointers
!
      ix = 3*(n1i - 1) + 1
      iy = ix + 1
      iz = iy + 1
      jx = 3*(n1j - 1) + 1
      jy = jx + 1
      jz = jy + 1
    endif
!  
!  If i and j are not valid then there is no point continuing for this value of n1
!           
    if (liok.and.ljok) then
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
            id21 = ineigh(1,n2,nri)
            id22 = ineigh(2,n2,nri)
            id23 = ineigh(3,n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            iik1 = - i21
            iik2 = - i22
            iik3 = - i23
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!  
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!               
!  Calculate vector between atoms
!             
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
              id21 = (ineigh(1,n2l,nri) - ineigh(1,n2k,nri))
              id22 = (ineigh(2,n2l,nri) - ineigh(2,n2k,nri))
              id23 = (ineigh(3,n2l,nri) - ineigh(3,n2k,nri))
!
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
              iik1 = (ineigh(1,n2k,nri) - i21)
              iik2 = (ineigh(2,n2k,nri) - i22)
              iik3 = (ineigh(3,n2k,nri) - i23)
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          else
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!               
!  Check that atoms are valid
!             
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(neighno(nji,nri))) 
!               
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
                iik1 = - i21
                iik2 = - i22
                iik3 = - i23
                id21 = (ineigh(1,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(1,nji,nri))
                id22 = (ineigh(2,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(2,nji,nri))
                id23 = (ineigh(3,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(3,nji,nri))
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                iik1 = (ineigh(1,n2k,nri) - i21)
                iik2 = (ineigh(2,n2k,nri) - i22)
                iik3 = (ineigh(3,n2k,nri) - i23)
                id21 = (ineigh(1,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(1,nji,nri) - ineigh(1,n2k,nri))
                id22 = (ineigh(2,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(2,nji,nri) - ineigh(2,n2k,nri))
                id23 = (ineigh(3,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(3,nji,nri) - ineigh(3,n2k,nri))
!
!  Convert neighbour number to real atom
!
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          endif
!               
!  No point continuing beyond here unless k and l are valid
!                 
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xil = xik + xd2
            yil = yik + yd2
            zil = zik + zd2
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Find cell index for i-k
!
            if (abs(iik1).gt.nd2cell(1).or. &
                abs(iik2).gt.nd2cell(2).or. &
                abs(iik3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
              ncindikp = nd2central
              ncindikm = nd2central
            else
!
!  Compute index
!
              ncindikp = nd2cellptr(nd2cell(1)+1+iik1,nd2cell(2)+1+iik2,nd2cell(3)+1+iik3)
              ncindikm = nd2cellptr(nd2cell(1)+1-iik1,nd2cell(2)+1-iik2,nd2cell(3)+1-iik3)
            endif
!
!  Find cell index for i-l
!
            iil1 = iik1 + id21
            iil2 = iik2 + id22
            iil3 = iik3 + id23
            if (abs(iil1).gt.nd2cell(1).or. &
                abs(iil2).gt.nd2cell(2).or. &
                abs(iil3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
              ncindilp = nd2central
              ncindilm = nd2central
            else
!
!  Compute index
!
              ncindilp = nd2cellptr(nd2cell(1)+1+iil1,nd2cell(2)+1+iil2,nd2cell(3)+1+iil3)
              ncindilm = nd2cellptr(nd2cell(1)+1-iil1,nd2cell(2)+1-iil2,nd2cell(3)+1-iil3)
            endif
!
!  Find cell index for j-k
!
            ijk1 = iik1 - id11
            ijk2 = iik2 - id12
            ijk3 = iik3 - id13
            if (abs(ijk1).gt.nd2cell(1).or. &
                abs(ijk2).gt.nd2cell(2).or. &
                abs(ijk3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
              ncindjkp = nd2central
              ncindjkm = nd2central
            else
!
!  Compute index
!
              ncindjkp = nd2cellptr(nd2cell(1)+1+ijk1,nd2cell(2)+1+ijk2,nd2cell(3)+1+ijk3)
              ncindjkm = nd2cellptr(nd2cell(1)+1-ijk1,nd2cell(2)+1-ijk2,nd2cell(3)+1-ijk3)
            endif
!
!  Find cell index for j-l
!
            ijl1 = ijk1 + id21
            ijl2 = ijk2 + id22
            ijl3 = ijk3 + id23
            if (abs(ijl1).gt.nd2cell(1).or. & 
                abs(ijl2).gt.nd2cell(2).or. & 
                abs(ijl3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
              ncindjlp = nd2central
              ncindjlm = nd2central
            else
!
!  Compute index
!
              ncindjlp = nd2cellptr(nd2cell(1)+1+ijl1,nd2cell(2)+1+ijl2,nd2cell(3)+1+ijl3)
              ncindjlm = nd2cellptr(nd2cell(1)+1-ijl1,nd2cell(2)+1-ijl2,nd2cell(3)+1-ijl3)
            endif
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
            kx = 3*(n2k - 1) + 1
            ky = kx + 1
            kz = ky + 1
            lx = 3*(n2l - 1) + 1
            ly = lx + 1
            lz = ly + 1
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
!  Second derivatives : i - k
!
            if (lsameijkl) then
              oneik = 0.5_dp
            else
              oneik = 1.0_dp
            endif
            d2cell(kx,ix,ncindikp) = d2cell(kx,ix,ncindikp) + xd1*xd2*d2(ind)*oneik
            d2cell(ky,ix,ncindikp) = d2cell(ky,ix,ncindikp) + xd1*yd2*d2(ind)*oneik
            d2cell(kz,ix,ncindikp) = d2cell(kz,ix,ncindikp) + xd1*zd2*d2(ind)*oneik
            d2cell(kx,iy,ncindikp) = d2cell(kx,iy,ncindikp) + yd1*xd2*d2(ind)*oneik
            d2cell(ky,iy,ncindikp) = d2cell(ky,iy,ncindikp) + yd1*yd2*d2(ind)*oneik
            d2cell(kz,iy,ncindikp) = d2cell(kz,iy,ncindikp) + yd1*zd2*d2(ind)*oneik
            d2cell(kx,iz,ncindikp) = d2cell(kx,iz,ncindikp) + zd1*xd2*d2(ind)*oneik
            d2cell(ky,iz,ncindikp) = d2cell(ky,iz,ncindikp) + zd1*yd2*d2(ind)*oneik
            d2cell(kz,iz,ncindikp) = d2cell(kz,iz,ncindikp) + zd1*zd2*d2(ind)*oneik
!
            if (n1.eq.n2) then
              d2cell(kx,ix,ncindikp) = d2cell(kx,ix,ncindikp) + d1(n2)*oneik
              d2cell(ky,iy,ncindikp) = d2cell(ky,iy,ncindikp) + d1(n2)*oneik
              d2cell(kz,iz,ncindikp) = d2cell(kz,iz,ncindikp) + d1(n2)*oneik
            endif
!
            d2cell(ix,kx,ncindikm) = d2cell(ix,kx,ncindikm) + xd1*xd2*d2(ind)*oneik
            d2cell(ix,ky,ncindikm) = d2cell(ix,ky,ncindikm) + xd1*yd2*d2(ind)*oneik
            d2cell(ix,kz,ncindikm) = d2cell(ix,kz,ncindikm) + xd1*zd2*d2(ind)*oneik
            d2cell(iy,kx,ncindikm) = d2cell(iy,kx,ncindikm) + yd1*xd2*d2(ind)*oneik
            d2cell(iy,ky,ncindikm) = d2cell(iy,ky,ncindikm) + yd1*yd2*d2(ind)*oneik
            d2cell(iy,kz,ncindikm) = d2cell(iy,kz,ncindikm) + yd1*zd2*d2(ind)*oneik
            d2cell(iz,kx,ncindikm) = d2cell(iz,kx,ncindikm) + zd1*xd2*d2(ind)*oneik
            d2cell(iz,ky,ncindikm) = d2cell(iz,ky,ncindikm) + zd1*yd2*d2(ind)*oneik
            d2cell(iz,kz,ncindikm) = d2cell(iz,kz,ncindikm) + zd1*zd2*d2(ind)*oneik
!
            if (n1.eq.n2) then
              d2cell(ix,kx,ncindikm) = d2cell(ix,kx,ncindikm) + d1(n2)*oneik
              d2cell(iy,ky,ncindikm) = d2cell(iy,ky,ncindikm) + d1(n2)*oneik
              d2cell(iz,kz,ncindikm) = d2cell(iz,kz,ncindikm) + d1(n2)*oneik
            endif
!
!  Second derivatives : i - l
!
            if (lsameijkl) then
              oneil = 0.5_dp
            else
              oneil = 1.0_dp
            endif
!
            d2cell(lx,ix,ncindilp) = d2cell(lx,ix,ncindilp) - xd1*xd2*d2(ind)*oneil
            d2cell(ly,ix,ncindilp) = d2cell(ly,ix,ncindilp) - xd1*yd2*d2(ind)*oneil
            d2cell(lz,ix,ncindilp) = d2cell(lz,ix,ncindilp) - xd1*zd2*d2(ind)*oneil
            d2cell(lx,iy,ncindilp) = d2cell(lx,iy,ncindilp) - yd1*xd2*d2(ind)*oneil
            d2cell(ly,iy,ncindilp) = d2cell(ly,iy,ncindilp) - yd1*yd2*d2(ind)*oneil
            d2cell(lz,iy,ncindilp) = d2cell(lz,iy,ncindilp) - yd1*zd2*d2(ind)*oneil
            d2cell(lx,iz,ncindilp) = d2cell(lx,iz,ncindilp) - zd1*xd2*d2(ind)*oneil
            d2cell(ly,iz,ncindilp) = d2cell(ly,iz,ncindilp) - zd1*yd2*d2(ind)*oneil
            d2cell(lz,iz,ncindilp) = d2cell(lz,iz,ncindilp) - zd1*zd2*d2(ind)*oneil
!
            if (n1.eq.n2) then
              d2cell(lx,ix,ncindilp) = d2cell(lx,ix,ncindilp) - d1(n2)*oneil
              d2cell(ly,iy,ncindilp) = d2cell(ly,iy,ncindilp) - d1(n2)*oneil
              d2cell(lz,iz,ncindilp) = d2cell(lz,iz,ncindilp) - d1(n2)*oneil
            endif
!
            d2cell(ix,lx,ncindilm) = d2cell(ix,lx,ncindilm) - xd1*xd2*d2(ind)*oneil
            d2cell(ix,ly,ncindilm) = d2cell(ix,ly,ncindilm) - xd1*yd2*d2(ind)*oneil
            d2cell(ix,lz,ncindilm) = d2cell(ix,lz,ncindilm) - xd1*zd2*d2(ind)*oneil
            d2cell(iy,lx,ncindilm) = d2cell(iy,lx,ncindilm) - yd1*xd2*d2(ind)*oneil
            d2cell(iy,ly,ncindilm) = d2cell(iy,ly,ncindilm) - yd1*yd2*d2(ind)*oneil
            d2cell(iy,lz,ncindilm) = d2cell(iy,lz,ncindilm) - yd1*zd2*d2(ind)*oneil
            d2cell(iz,lx,ncindilm) = d2cell(iz,lx,ncindilm) - zd1*xd2*d2(ind)*oneil
            d2cell(iz,ly,ncindilm) = d2cell(iz,ly,ncindilm) - zd1*yd2*d2(ind)*oneil
            d2cell(iz,lz,ncindilm) = d2cell(iz,lz,ncindilm) - zd1*zd2*d2(ind)*oneil
!
            if (n1.eq.n2) then
              d2cell(ix,lx,ncindilm) = d2cell(ix,lx,ncindilm) - d1(n2)*oneil
              d2cell(iy,ly,ncindilm) = d2cell(iy,ly,ncindilm) - d1(n2)*oneil
              d2cell(iz,lz,ncindilm) = d2cell(iz,lz,ncindilm) - d1(n2)*oneil
            endif
!
!  Second derivatives : j - k
!
            if (lsameijkl) then
              onejk = 0.5_dp
            else
              onejk = 1.0_dp
            endif
!
            d2cell(kx,jx,ncindjkp) = d2cell(kx,jx,ncindjkp) - xd1*xd2*d2(ind)*onejk
            d2cell(ky,jx,ncindjkp) = d2cell(ky,jx,ncindjkp) - xd1*yd2*d2(ind)*onejk
            d2cell(kz,jx,ncindjkp) = d2cell(kz,jx,ncindjkp) - xd1*zd2*d2(ind)*onejk
            d2cell(kx,jy,ncindjkp) = d2cell(kx,jy,ncindjkp) - yd1*xd2*d2(ind)*onejk
            d2cell(ky,jy,ncindjkp) = d2cell(ky,jy,ncindjkp) - yd1*yd2*d2(ind)*onejk
            d2cell(kz,jy,ncindjkp) = d2cell(kz,jy,ncindjkp) - yd1*zd2*d2(ind)*onejk
            d2cell(kx,jz,ncindjkp) = d2cell(kx,jz,ncindjkp) - zd1*xd2*d2(ind)*onejk
            d2cell(ky,jz,ncindjkp) = d2cell(ky,jz,ncindjkp) - zd1*yd2*d2(ind)*onejk
            d2cell(kz,jz,ncindjkp) = d2cell(kz,jz,ncindjkp) - zd1*zd2*d2(ind)*onejk
!
            if (n1.eq.n2) then
              d2cell(kx,jx,ncindjkp) = d2cell(kx,jx,ncindjkp) - d1(n2)*onejk
              d2cell(ky,jy,ncindjkp) = d2cell(ky,jy,ncindjkp) - d1(n2)*onejk
              d2cell(kz,jz,ncindjkp) = d2cell(kz,jz,ncindjkp) - d1(n2)*onejk
            endif
!
            d2cell(jx,kx,ncindjkm) = d2cell(jx,kx,ncindjkm) - xd1*xd2*d2(ind)*onejk
            d2cell(jx,ky,ncindjkm) = d2cell(jx,ky,ncindjkm) - xd1*yd2*d2(ind)*onejk
            d2cell(jx,kz,ncindjkm) = d2cell(jx,kz,ncindjkm) - xd1*zd2*d2(ind)*onejk
            d2cell(jy,kx,ncindjkm) = d2cell(jy,kx,ncindjkm) - yd1*xd2*d2(ind)*onejk
            d2cell(jy,ky,ncindjkm) = d2cell(jy,ky,ncindjkm) - yd1*yd2*d2(ind)*onejk
            d2cell(jy,kz,ncindjkm) = d2cell(jy,kz,ncindjkm) - yd1*zd2*d2(ind)*onejk
            d2cell(jz,kx,ncindjkm) = d2cell(jz,kx,ncindjkm) - zd1*xd2*d2(ind)*onejk
            d2cell(jz,ky,ncindjkm) = d2cell(jz,ky,ncindjkm) - zd1*yd2*d2(ind)*onejk
            d2cell(jz,kz,ncindjkm) = d2cell(jz,kz,ncindjkm) - zd1*zd2*d2(ind)*onejk
!
            if (n1.eq.n2) then
              d2cell(jx,kx,ncindjkm) = d2cell(jx,kx,ncindjkm) - d1(n2)*onejk
              d2cell(jy,ky,ncindjkm) = d2cell(jy,ky,ncindjkm) - d1(n2)*onejk
              d2cell(jz,kz,ncindjkm) = d2cell(jz,kz,ncindjkm) - d1(n2)*onejk
            endif
!
!  Second derivatives : j - l
!
            if (lsameijkl) then
              onejl = 0.5_dp
            else
              onejl = 1.0_dp
            endif
!
            d2cell(lx,jx,ncindjlp) = d2cell(lx,jx,ncindjlp) + xd1*xd2*d2(ind)*onejl
            d2cell(ly,jx,ncindjlp) = d2cell(ly,jx,ncindjlp) + xd1*yd2*d2(ind)*onejl
            d2cell(lz,jx,ncindjlp) = d2cell(lz,jx,ncindjlp) + xd1*zd2*d2(ind)*onejl
            d2cell(lx,jy,ncindjlp) = d2cell(lx,jy,ncindjlp) + yd1*xd2*d2(ind)*onejl
            d2cell(ly,jy,ncindjlp) = d2cell(ly,jy,ncindjlp) + yd1*yd2*d2(ind)*onejl
            d2cell(lz,jy,ncindjlp) = d2cell(lz,jy,ncindjlp) + yd1*zd2*d2(ind)*onejl
            d2cell(lx,jz,ncindjlp) = d2cell(lx,jz,ncindjlp) + zd1*xd2*d2(ind)*onejl
            d2cell(ly,jz,ncindjlp) = d2cell(ly,jz,ncindjlp) + zd1*yd2*d2(ind)*onejl
            d2cell(lz,jz,ncindjlp) = d2cell(lz,jz,ncindjlp) + zd1*zd2*d2(ind)*onejl
!
            if (n1.eq.n2) then
              d2cell(lx,jx,ncindjlp) = d2cell(lx,jx,ncindjlp) + d1(n2)*onejl
              d2cell(ly,jy,ncindjlp) = d2cell(ly,jy,ncindjlp) + d1(n2)*onejl
              d2cell(lz,jz,ncindjlp) = d2cell(lz,jz,ncindjlp) + d1(n2)*onejl
            endif
!
            d2cell(jx,lx,ncindjlm) = d2cell(jx,lx,ncindjlm) + xd1*xd2*d2(ind)*onejl
            d2cell(jx,ly,ncindjlm) = d2cell(jx,ly,ncindjlm) + xd1*yd2*d2(ind)*onejl
            d2cell(jx,lz,ncindjlm) = d2cell(jx,lz,ncindjlm) + xd1*zd2*d2(ind)*onejl
            d2cell(jy,lx,ncindjlm) = d2cell(jy,lx,ncindjlm) + yd1*xd2*d2(ind)*onejl
            d2cell(jy,ly,ncindjlm) = d2cell(jy,ly,ncindjlm) + yd1*yd2*d2(ind)*onejl
            d2cell(jy,lz,ncindjlm) = d2cell(jy,lz,ncindjlm) + yd1*zd2*d2(ind)*onejl
            d2cell(jz,lx,ncindjlm) = d2cell(jz,lx,ncindjlm) + zd1*xd2*d2(ind)*onejl
            d2cell(jz,ly,ncindjlm) = d2cell(jz,ly,ncindjlm) + zd1*yd2*d2(ind)*onejl
            d2cell(jz,lz,ncindjlm) = d2cell(jz,lz,ncindjlm) + zd1*zd2*d2(ind)*onejl
!
            if (n1.eq.n2) then
              d2cell(jx,lx,ncindjlm) = d2cell(jx,lx,ncindjlm) + d1(n2)*onejl
              d2cell(jy,ly,ncindjlm) = d2cell(jy,ly,ncindjlm) + d1(n2)*onejl
              d2cell(jz,lz,ncindjlm) = d2cell(jz,lz,ncindjlm) + d1(n2)*onejl
            endif
          endif
        endif
      enddo
    else    
!               
!  Increment ind pointer by number of missed terms
!  
      ind = ind + n1
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addfc')
#endif
!
  return
  end
!
  subroutine d2addfcd(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,ineigh,nREBOatomRptr,d1,d2,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian unphased second derivative arrays.
!
!  Distributed memory parallel version
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  ineigh          = cell indices for vector to neighbour of atoms
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!
!   4/17 Created from d2addfc
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use derivatives
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: ineigh(3_i4,maxneigh,*)
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
!
!  Local variables
!
  integer(i4)                         :: id11
  integer(i4)                         :: id12
  integer(i4)                         :: id13
  integer(i4)                         :: id21
  integer(i4)                         :: id22
  integer(i4)                         :: id23
  integer(i4)                         :: i21
  integer(i4)                         :: i22
  integer(i4)                         :: i23
  integer(i4)                         :: iik1
  integer(i4)                         :: iik2
  integer(i4)                         :: iik3
  integer(i4)                         :: iil1
  integer(i4)                         :: iil2
  integer(i4)                         :: iil3
  integer(i4)                         :: ijk1
  integer(i4)                         :: ijk2
  integer(i4)                         :: ijk3
  integer(i4)                         :: ijl1
  integer(i4)                         :: ijl2
  integer(i4)                         :: ijl3
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: ixf,iyf,izf
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: jxf,jyf,jzf
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: kxf,kyf,kzf
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: lxf,lyf,lzf
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ncindikm
  integer(i4)                         :: ncindikp
  integer(i4)                         :: ncindilm
  integer(i4)                         :: ncindilp
  integer(i4)                         :: ncindjkm
  integer(i4)                         :: ncindjkp
  integer(i4)                         :: ncindjlm
  integer(i4)                         :: ncindjlp
  integer(i4)                         :: ndist
  integer(i4)                         :: nji
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: ntmp
  logical                             :: lilocal
  logical                             :: lijlocal
  logical                             :: ljlocal
  logical                             :: lklocal
  logical                             :: lkllocal
  logical                             :: lllocal
  real(dp)                            :: oneik
  real(dp)                            :: oneil
  real(dp)                            :: onejk
  real(dp)                            :: onejl
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i
  real(dp)                            :: z2i
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lsameijkl
#ifdef TRACE
  call trace_in('d2addfcd')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndist = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
  else
    ndist = nneigh2
  endif
!
!  Loop over pairs of interatomic distances
!
  ind = 0
!
!  id11 / id12 / id13 are the cell indices between i and j
!
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      id11 = ineigh(1,n1,nri)
      id12 = ineigh(2,n1,nri)
      id13 = ineigh(3,n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      i21 = 0
      i22 = 0
      i23 = 0
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!     
!  Check that neighbour numbers are valid
!         
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))       
      if (liok.and.ljok) then
!         
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
        id11 = (ineigh(1,n1j,nri) - ineigh(1,n1i,nri))
        id12 = (ineigh(2,n1j,nri) - ineigh(2,n1i,nri))
        id13 = (ineigh(3,n1j,nri) - ineigh(3,n1i,nri))
!
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
        i21 = ineigh(1,n1i,nri)
        i22 = ineigh(2,n1i,nri)
        i23 = ineigh(3,n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    else
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!         
!  Check that atoms are valid
!       
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
!  
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
          id11 = (ineigh(1,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(1,nji,nri))
          id12 = (ineigh(2,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(2,nji,nri))
          id13 = (ineigh(3,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(3,nji,nri))
          i21 = 0
          i22 = 0
          i23 = 0
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          id11 = (ineigh(1,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(1,nji,nri) - ineigh(1,n1i,nri))
          id12 = (ineigh(2,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(2,nji,nri) - ineigh(2,n1i,nri))
          id13 = (ineigh(3,n1j,nREBOatomRptr(neighno(nji,nri))) + ineigh(3,nji,nri) - ineigh(3,n1i,nri))
          i21 = ineigh(1,n1i,nri)
          i22 = ineigh(2,n1i,nri)
          i23 = ineigh(3,n1i,nri)
!
!  Convert neighbour number to real atom
!
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    endif
!  
!  Skip the following section as it is irrelevant if i and j are not valid
!           
    if (liok.and.ljok) then
!
      lilocal = (atom2local(n1i).gt.0)
      ljlocal = (atom2local(n1j).gt.0)
      lijlocal = (lilocal.or.ljlocal)
!
!  Set second derivative location pointers
!
      ixf = 3*(n1i - 1) + 1
      iyf = ixf + 1
      izf = iyf + 1
!
      if (lilocal) then
        ix = 3*(atom2local(n1i) - 1) + 1
        iy = ix + 1
        iz = iy + 1
      endif
!
      jxf = 3*(n1j - 1) + 1
      jyf = jxf + 1
      jzf = jyf + 1
!
      if (ljlocal) then
        jx = 3*(atom2local(n1j) - 1) + 1
        jy = jx + 1
        jz = jy + 1
      endif
!
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
            id21 = ineigh(1,n2,nri)
            id22 = ineigh(2,n2,nri)
            id23 = ineigh(3,n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            iik1 = - i21
            iik2 = - i22
            iik3 = - i23
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!  
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!               
!  Calculate vector between atoms
!             
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
              id21 = (ineigh(1,n2l,nri) - ineigh(1,n2k,nri))
              id22 = (ineigh(2,n2l,nri) - ineigh(2,n2k,nri))
              id23 = (ineigh(3,n2l,nri) - ineigh(3,n2k,nri))
!
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
              iik1 = (ineigh(1,n2k,nri) - i21)
              iik2 = (ineigh(2,n2k,nri) - i22)
              iik3 = (ineigh(3,n2k,nri) - i23)
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          else
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!               
!  Check that atoms are valid
!             
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(neighno(nji,nri))) 
!               
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
                iik1 = - i21
                iik2 = - i22
                iik3 = - i23
                id21 = (ineigh(1,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(1,nji,nri))
                id22 = (ineigh(2,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(2,nji,nri))
                id23 = (ineigh(3,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(3,nji,nri))
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                iik1 = (ineigh(1,n2k,nri) - i21)
                iik2 = (ineigh(2,n2k,nri) - i22)
                iik3 = (ineigh(3,n2k,nri) - i23)
                id21 = (ineigh(1,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(1,nji,nri) - ineigh(1,n2k,nri))
                id22 = (ineigh(2,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(2,nji,nri) - ineigh(2,n2k,nri))
                id23 = (ineigh(3,n2l,nREBOatomRptr(neighno(nji,nri))) + ineigh(3,nji,nri) - ineigh(3,n2k,nri))
!
!  Convert neighbour number to real atom
!
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          endif
!               
!  No point continuing beyond here unless k and l are valid
!                 
          if (lkok.and.llok) then
            lklocal = (atom2local(n2k).gt.0)
            lllocal = (atom2local(n2l).gt.0)
            lkllocal = (lklocal.or.lllocal)
!
!  Check whether any of i, j, k or l are local before going further
!
            if (lijlocal.or.lkllocal) then
!
!  Complete remaining vectors
!
              xil = xik + xd2
              yil = yik + yd2
              zil = zik + zd2
              xjk = xik - xd1
              yjk = yik - yd1
              zjk = zik - zd1
              xjl = xjk + xd2
              yjl = yjk + yd2
              zjl = zjk + zd2
!
!  Find cell index for i-k
!
              if (abs(iik1).gt.nd2cell(1).or. &
                  abs(iik2).gt.nd2cell(2).or. &
                  abs(iik3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncindikp = nd2central
                ncindikm = nd2central
              else
!
!  Compute index
!
                ncindikp = nd2cellptr(nd2cell(1)+1+iik1,nd2cell(2)+1+iik2,nd2cell(3)+1+iik3)
                ncindikm = nd2cellptr(nd2cell(1)+1-iik1,nd2cell(2)+1-iik2,nd2cell(3)+1-iik3)
              endif
!
!  Find cell index for i-l
!
              iil1 = iik1 + id21
              iil2 = iik2 + id22
              iil3 = iik3 + id23
              if (abs(iil1).gt.nd2cell(1).or. &
                  abs(iil2).gt.nd2cell(2).or. &
                  abs(iil3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncindilp = nd2central
                ncindilm = nd2central
              else
!
!  Compute index
!
                ncindilp = nd2cellptr(nd2cell(1)+1+iil1,nd2cell(2)+1+iil2,nd2cell(3)+1+iil3)
                ncindilm = nd2cellptr(nd2cell(1)+1-iil1,nd2cell(2)+1-iil2,nd2cell(3)+1-iil3)
              endif
!
!  Find cell index for j-k
!
              ijk1 = iik1 - id11
              ijk2 = iik2 - id12
              ijk3 = iik3 - id13
              if (abs(ijk1).gt.nd2cell(1).or. &
                  abs(ijk2).gt.nd2cell(2).or. &
                  abs(ijk3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncindjkp = nd2central
                ncindjkm = nd2central
              else
!
!  Compute index
!
                ncindjkp = nd2cellptr(nd2cell(1)+1+ijk1,nd2cell(2)+1+ijk2,nd2cell(3)+1+ijk3)
                ncindjkm = nd2cellptr(nd2cell(1)+1-ijk1,nd2cell(2)+1-ijk2,nd2cell(3)+1-ijk3)
              endif
!
!  Find cell index for j-l
!
              ijl1 = ijk1 + id21
              ijl2 = ijk2 + id22
              ijl3 = ijk3 + id23
              if (abs(ijl1).gt.nd2cell(1).or. & 
                  abs(ijl2).gt.nd2cell(2).or. & 
                  abs(ijl3).gt.nd2cell(3)) then
!
!  Find cell index - if outside user range then assign to central cell
!
                ncindjlp = nd2central
                ncindjlm = nd2central
              else
!
!  Compute index
!
                ncindjlp = nd2cellptr(nd2cell(1)+1+ijl1,nd2cell(2)+1+ijl2,nd2cell(3)+1+ijl3)
                ncindjlm = nd2cellptr(nd2cell(1)+1-ijl1,nd2cell(2)+1-ijl2,nd2cell(3)+1-ijl3)
              endif
!
!  Compute distances squared
!
              r2ik = xik*xik + yik*yik + zik*zik
              r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
              kxf = 3*(n2k - 1) + 1
              kyf = kxf + 1
              kzf = kyf + 1
!
              if (lklocal) then
                kx = 3*(atom2local(n2k) - 1) + 1
                ky = kx + 1
                kz = ky + 1
              endif
!
              lxf = 3*(n2l - 1) + 1
              lyf = lxf + 1
              lzf = lyf + 1
!
              if (lllocal) then
                lx = 3*(atom2local(n2l) - 1) + 1
                ly = lx + 1
                lz = ly + 1
              endif
!
              lsameijkl = .false.
              if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
!  Second derivatives : i - k
!
              if (lsameijkl) then
                oneik = 0.5_dp
              else
                oneik = 1.0_dp
              endif
              if (lilocal) then
                d2cell(kxf,ix,ncindikp) = d2cell(kxf,ix,ncindikp) + xd1*xd2*d2(ind)*oneik
                d2cell(kyf,ix,ncindikp) = d2cell(kyf,ix,ncindikp) + xd1*yd2*d2(ind)*oneik
                d2cell(kzf,ix,ncindikp) = d2cell(kzf,ix,ncindikp) + xd1*zd2*d2(ind)*oneik
                d2cell(kxf,iy,ncindikp) = d2cell(kxf,iy,ncindikp) + yd1*xd2*d2(ind)*oneik
                d2cell(kyf,iy,ncindikp) = d2cell(kyf,iy,ncindikp) + yd1*yd2*d2(ind)*oneik
                d2cell(kzf,iy,ncindikp) = d2cell(kzf,iy,ncindikp) + yd1*zd2*d2(ind)*oneik
                d2cell(kxf,iz,ncindikp) = d2cell(kxf,iz,ncindikp) + zd1*xd2*d2(ind)*oneik
                d2cell(kyf,iz,ncindikp) = d2cell(kyf,iz,ncindikp) + zd1*yd2*d2(ind)*oneik
                d2cell(kzf,iz,ncindikp) = d2cell(kzf,iz,ncindikp) + zd1*zd2*d2(ind)*oneik
!
                if (n1.eq.n2) then
                  d2cell(kxf,ix,ncindikp) = d2cell(kxf,ix,ncindikp) + d1(n2)*oneik
                  d2cell(kyf,iy,ncindikp) = d2cell(kyf,iy,ncindikp) + d1(n2)*oneik
                  d2cell(kzf,iz,ncindikp) = d2cell(kzf,iz,ncindikp) + d1(n2)*oneik
                endif
              endif
!
              if (lklocal) then
                d2cell(ixf,kx,ncindikm) = d2cell(ixf,kx,ncindikm) + xd1*xd2*d2(ind)*oneik
                d2cell(ixf,ky,ncindikm) = d2cell(ixf,ky,ncindikm) + xd1*yd2*d2(ind)*oneik
                d2cell(ixf,kz,ncindikm) = d2cell(ixf,kz,ncindikm) + xd1*zd2*d2(ind)*oneik
                d2cell(iyf,kx,ncindikm) = d2cell(iyf,kx,ncindikm) + yd1*xd2*d2(ind)*oneik
                d2cell(iyf,ky,ncindikm) = d2cell(iyf,ky,ncindikm) + yd1*yd2*d2(ind)*oneik
                d2cell(iyf,kz,ncindikm) = d2cell(iyf,kz,ncindikm) + yd1*zd2*d2(ind)*oneik
                d2cell(izf,kx,ncindikm) = d2cell(izf,kx,ncindikm) + zd1*xd2*d2(ind)*oneik
                d2cell(izf,ky,ncindikm) = d2cell(izf,ky,ncindikm) + zd1*yd2*d2(ind)*oneik
                d2cell(izf,kz,ncindikm) = d2cell(izf,kz,ncindikm) + zd1*zd2*d2(ind)*oneik
!
                if (n1.eq.n2) then
                  d2cell(ixf,kx,ncindikm) = d2cell(ixf,kx,ncindikm) + d1(n2)*oneik
                  d2cell(iyf,ky,ncindikm) = d2cell(iyf,ky,ncindikm) + d1(n2)*oneik
                  d2cell(izf,kz,ncindikm) = d2cell(izf,kz,ncindikm) + d1(n2)*oneik
                endif
              endif
!
!  Second derivatives : i - l
!
              if (lsameijkl) then
                oneil = 0.5_dp
              else
                oneil = 1.0_dp
              endif
!
              if (lilocal) then
                d2cell(lxf,ix,ncindilp) = d2cell(lxf,ix,ncindilp) - xd1*xd2*d2(ind)*oneil
                d2cell(lyf,ix,ncindilp) = d2cell(lyf,ix,ncindilp) - xd1*yd2*d2(ind)*oneil
                d2cell(lzf,ix,ncindilp) = d2cell(lzf,ix,ncindilp) - xd1*zd2*d2(ind)*oneil
                d2cell(lxf,iy,ncindilp) = d2cell(lxf,iy,ncindilp) - yd1*xd2*d2(ind)*oneil
                d2cell(lyf,iy,ncindilp) = d2cell(lyf,iy,ncindilp) - yd1*yd2*d2(ind)*oneil
                d2cell(lzf,iy,ncindilp) = d2cell(lzf,iy,ncindilp) - yd1*zd2*d2(ind)*oneil
                d2cell(lxf,iz,ncindilp) = d2cell(lxf,iz,ncindilp) - zd1*xd2*d2(ind)*oneil
                d2cell(lyf,iz,ncindilp) = d2cell(lyf,iz,ncindilp) - zd1*yd2*d2(ind)*oneil
                d2cell(lzf,iz,ncindilp) = d2cell(lzf,iz,ncindilp) - zd1*zd2*d2(ind)*oneil
!
                if (n1.eq.n2) then
                  d2cell(lxf,ix,ncindilp) = d2cell(lxf,ix,ncindilp) - d1(n2)*oneil
                  d2cell(lyf,iy,ncindilp) = d2cell(lyf,iy,ncindilp) - d1(n2)*oneil
                  d2cell(lzf,iz,ncindilp) = d2cell(lzf,iz,ncindilp) - d1(n2)*oneil
                endif
              endif
!
              if (lllocal) then
                d2cell(ixf,lx,ncindilm) = d2cell(ixf,lx,ncindilm) - xd1*xd2*d2(ind)*oneil
                d2cell(ixf,ly,ncindilm) = d2cell(ixf,ly,ncindilm) - xd1*yd2*d2(ind)*oneil
                d2cell(ixf,lz,ncindilm) = d2cell(ixf,lz,ncindilm) - xd1*zd2*d2(ind)*oneil
                d2cell(iyf,lx,ncindilm) = d2cell(iyf,lx,ncindilm) - yd1*xd2*d2(ind)*oneil
                d2cell(iyf,ly,ncindilm) = d2cell(iyf,ly,ncindilm) - yd1*yd2*d2(ind)*oneil
                d2cell(iyf,lz,ncindilm) = d2cell(iyf,lz,ncindilm) - yd1*zd2*d2(ind)*oneil
                d2cell(izf,lx,ncindilm) = d2cell(izf,lx,ncindilm) - zd1*xd2*d2(ind)*oneil
                d2cell(izf,ly,ncindilm) = d2cell(izf,ly,ncindilm) - zd1*yd2*d2(ind)*oneil
                d2cell(izf,lz,ncindilm) = d2cell(izf,lz,ncindilm) - zd1*zd2*d2(ind)*oneil
!
                if (n1.eq.n2) then
                  d2cell(ixf,lx,ncindilm) = d2cell(ixf,lx,ncindilm) - d1(n2)*oneil
                  d2cell(iyf,ly,ncindilm) = d2cell(iyf,ly,ncindilm) - d1(n2)*oneil
                  d2cell(izf,lz,ncindilm) = d2cell(izf,lz,ncindilm) - d1(n2)*oneil
                endif
              endif
!
!  Second derivatives : j - k
!
              if (lsameijkl) then
                onejk = 0.5_dp
              else
                onejk = 1.0_dp
              endif
!
              if (ljlocal) then
                d2cell(kxf,jx,ncindjkp) = d2cell(kxf,jx,ncindjkp) - xd1*xd2*d2(ind)*onejk
                d2cell(kyf,jx,ncindjkp) = d2cell(kyf,jx,ncindjkp) - xd1*yd2*d2(ind)*onejk
                d2cell(kzf,jx,ncindjkp) = d2cell(kzf,jx,ncindjkp) - xd1*zd2*d2(ind)*onejk
                d2cell(kxf,jy,ncindjkp) = d2cell(kxf,jy,ncindjkp) - yd1*xd2*d2(ind)*onejk
                d2cell(kyf,jy,ncindjkp) = d2cell(kyf,jy,ncindjkp) - yd1*yd2*d2(ind)*onejk
                d2cell(kzf,jy,ncindjkp) = d2cell(kzf,jy,ncindjkp) - yd1*zd2*d2(ind)*onejk
                d2cell(kxf,jz,ncindjkp) = d2cell(kxf,jz,ncindjkp) - zd1*xd2*d2(ind)*onejk
                d2cell(kyf,jz,ncindjkp) = d2cell(kyf,jz,ncindjkp) - zd1*yd2*d2(ind)*onejk
                d2cell(kzf,jz,ncindjkp) = d2cell(kzf,jz,ncindjkp) - zd1*zd2*d2(ind)*onejk
!
                if (n1.eq.n2) then
                  d2cell(kxf,jx,ncindjkp) = d2cell(kxf,jx,ncindjkp) - d1(n2)*onejk
                  d2cell(kyf,jy,ncindjkp) = d2cell(kyf,jy,ncindjkp) - d1(n2)*onejk
                  d2cell(kzf,jz,ncindjkp) = d2cell(kzf,jz,ncindjkp) - d1(n2)*onejk
                endif
              endif
!
              if (lklocal) then
                d2cell(jxf,kx,ncindjkm) = d2cell(jxf,kx,ncindjkm) - xd1*xd2*d2(ind)*onejk
                d2cell(jxf,ky,ncindjkm) = d2cell(jxf,ky,ncindjkm) - xd1*yd2*d2(ind)*onejk
                d2cell(jxf,kz,ncindjkm) = d2cell(jxf,kz,ncindjkm) - xd1*zd2*d2(ind)*onejk
                d2cell(jyf,kx,ncindjkm) = d2cell(jyf,kx,ncindjkm) - yd1*xd2*d2(ind)*onejk
                d2cell(jyf,ky,ncindjkm) = d2cell(jyf,ky,ncindjkm) - yd1*yd2*d2(ind)*onejk
                d2cell(jyf,kz,ncindjkm) = d2cell(jyf,kz,ncindjkm) - yd1*zd2*d2(ind)*onejk
                d2cell(jzf,kx,ncindjkm) = d2cell(jzf,kx,ncindjkm) - zd1*xd2*d2(ind)*onejk
                d2cell(jzf,ky,ncindjkm) = d2cell(jzf,ky,ncindjkm) - zd1*yd2*d2(ind)*onejk
                d2cell(jzf,kz,ncindjkm) = d2cell(jzf,kz,ncindjkm) - zd1*zd2*d2(ind)*onejk
!
                if (n1.eq.n2) then
                  d2cell(jxf,kx,ncindjkm) = d2cell(jxf,kx,ncindjkm) - d1(n2)*onejk
                  d2cell(jyf,ky,ncindjkm) = d2cell(jyf,ky,ncindjkm) - d1(n2)*onejk
                  d2cell(jzf,kz,ncindjkm) = d2cell(jzf,kz,ncindjkm) - d1(n2)*onejk
                endif
              endif
!
!  Second derivatives : j - l
!
              if (lsameijkl) then
                onejl = 0.5_dp
              else
                onejl = 1.0_dp
              endif
!
              if (ljlocal) then
                d2cell(lxf,jx,ncindjlp) = d2cell(lxf,jx,ncindjlp) + xd1*xd2*d2(ind)*onejl
                d2cell(lyf,jx,ncindjlp) = d2cell(lyf,jx,ncindjlp) + xd1*yd2*d2(ind)*onejl
                d2cell(lzf,jx,ncindjlp) = d2cell(lzf,jx,ncindjlp) + xd1*zd2*d2(ind)*onejl
                d2cell(lxf,jy,ncindjlp) = d2cell(lxf,jy,ncindjlp) + yd1*xd2*d2(ind)*onejl
                d2cell(lyf,jy,ncindjlp) = d2cell(lyf,jy,ncindjlp) + yd1*yd2*d2(ind)*onejl
                d2cell(lzf,jy,ncindjlp) = d2cell(lzf,jy,ncindjlp) + yd1*zd2*d2(ind)*onejl
                d2cell(lxf,jz,ncindjlp) = d2cell(lxf,jz,ncindjlp) + zd1*xd2*d2(ind)*onejl
                d2cell(lyf,jz,ncindjlp) = d2cell(lyf,jz,ncindjlp) + zd1*yd2*d2(ind)*onejl
                d2cell(lzf,jz,ncindjlp) = d2cell(lzf,jz,ncindjlp) + zd1*zd2*d2(ind)*onejl
!
                if (n1.eq.n2) then
                  d2cell(lxf,jx,ncindjlp) = d2cell(lxf,jx,ncindjlp) + d1(n2)*onejl
                  d2cell(lyf,jy,ncindjlp) = d2cell(lyf,jy,ncindjlp) + d1(n2)*onejl
                  d2cell(lzf,jz,ncindjlp) = d2cell(lzf,jz,ncindjlp) + d1(n2)*onejl
                endif
              endif
!
              if (lllocal) then
                d2cell(jxf,lx,ncindjlm) = d2cell(jxf,lx,ncindjlm) + xd1*xd2*d2(ind)*onejl
                d2cell(jxf,ly,ncindjlm) = d2cell(jxf,ly,ncindjlm) + xd1*yd2*d2(ind)*onejl
                d2cell(jxf,lz,ncindjlm) = d2cell(jxf,lz,ncindjlm) + xd1*zd2*d2(ind)*onejl
                d2cell(jyf,lx,ncindjlm) = d2cell(jyf,lx,ncindjlm) + yd1*xd2*d2(ind)*onejl
                d2cell(jyf,ly,ncindjlm) = d2cell(jyf,ly,ncindjlm) + yd1*yd2*d2(ind)*onejl
                d2cell(jyf,lz,ncindjlm) = d2cell(jyf,lz,ncindjlm) + yd1*zd2*d2(ind)*onejl
                d2cell(jzf,lx,ncindjlm) = d2cell(jzf,lx,ncindjlm) + zd1*xd2*d2(ind)*onejl
                d2cell(jzf,ly,ncindjlm) = d2cell(jzf,ly,ncindjlm) + zd1*yd2*d2(ind)*onejl
                d2cell(jzf,lz,ncindjlm) = d2cell(jzf,lz,ncindjlm) + zd1*zd2*d2(ind)*onejl
!
                if (n1.eq.n2) then
                  d2cell(jxf,lx,ncindjlm) = d2cell(jxf,lx,ncindjlm) + d1(n2)*onejl
                  d2cell(jyf,ly,ncindjlm) = d2cell(jyf,ly,ncindjlm) + d1(n2)*onejl
                  d2cell(jzf,lz,ncindjlm) = d2cell(jzf,lz,ncindjlm) + d1(n2)*onejl
                endif
              endif
            endif
          endif
        endif
      enddo
    else    
!               
!  Increment ind pointer by number of missed terms
!  
      ind = ind + n1
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addfcd')
#endif
!
  return
  end
!
  subroutine d2adddm(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr,d1,d2, &
                     ltorderv,lextended)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array.
!
!  Distributed memory parallel version.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!  lextended       = if true then include extended derivatives
!
!   1/17 Created from d2add
!   2/18 Trace added
!   9/18 Strain module added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  12/19 Rigid molecule modifications added
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   5/23 ind now set explicitly before n2 loop in case there is a cycle during extended derivatives
!   5/23 Extended derivatives added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use current,        only : ndim, nstrains
  use current,        only : nstrains, nstrptr
  use derivatives
  use parallel
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: is1,is2
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: ixf,iyf,izf
  integer(i4)                         :: jxf,jyf,jzf
  integer(i4)                         :: kxf,kyf,kzf
  integer(i4)                         :: lxf,lyf,lzf
  integer(i4)                         :: ix1,iy1,iz1
  integer(i4)                         :: ix2,iy2,iz2
  integer(i4)                         :: jx1,jy1,jz1
  integer(i4)                         :: jx2,jy2,jz2
  integer(i4)                         :: kx1,ky1,kz1
  integer(i4)                         :: kx2,ky2,kz2
  integer(i4)                         :: lx1,ly1,lz1
  integer(i4)                         :: lx2,ly2,lz2
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nji
  integer(i4)                         :: nkj
  integer(i4)                         :: nlk
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: ns1
  integer(i4)                         :: ns2
  integer(i4)                         :: ntmp
  real(dp)                            :: one1
  real(dp)                            :: one2
  real(dp)                            :: dr2ds1(6)
  real(dp)                            :: dr2ds2(6)
  real(dp)                            :: d2r2ds21(6,6)
  real(dp)                            :: d2r2ds22(6,6)
  real(dp)                            :: d2r2dsdx1(6,3)
  real(dp)                            :: d2r2dsdx2(6,3)
  real(dp)                            :: d2r2dx21(3,3)
  real(dp)                            :: d2r2dx22(3,3)
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: ldo12
  logical                             :: ldo21
  logical                             :: lsameijkl
  logical                             :: lilocal
  logical                             :: ljlocal
  logical                             :: lklocal
  logical                             :: lllocal
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lopi
  logical                             :: lopj
  logical                             :: lopk
  logical                             :: lopl
#ifdef TRACE
  call trace_in('d2adddm')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndisttor = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
    if (lextended) then
      ndist = ndisttor + nneigh(nri)*mneigh*mneigh
    else
      ndist = ndisttor
    endif
  else
    ndist = nneigh2 
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = nneigh2 + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
!  Loop over pairs of interatomic distances
!
  ind = 0
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Check that neighbour numbers are valid
!
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    elseif (n1.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
!  Check that atoms are valid
!
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
      else
        ljok = .false.
      endif
!
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!         
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    else
!
!  Extended derivatives
!
      ntmp = n1 - ndisttor
      nji = (ntmp - 1)/(mneigh*mneigh) + 1
      ntmp = ntmp - (nji - 1)*mneigh*mneigh
      nkj = (ntmp - 1)/mneigh + 1
      ntmp = ntmp - (nkj - 1)*mneigh
      nlk = ntmp
!
      nrj = neighno(nji,nri)
      if (nrj.eq.0) cycle
      if (nkj.gt.nneigh(nrj)) cycle
      nrk = neighno(nkj,nrj)
      if (nrk.eq.0) cycle
      if (nlk.gt.nneigh(nrk)) cycle
      n1i = nrk
      n1j = neighno(nlk,nrk)
      if (n1j.eq.0) cycle
!
!  Check that atoms are valid
!
      liok = (nkj.le.nneigh(nrj))
      ljok = (nlk.le.nneigh(nrk))
!
      if (liok.and.ljok) then
        x2i = xneigh(nji,nri) + xneigh(nkj,nrj)
        y2i = yneigh(nji,nri) + yneigh(nkj,nrj)
        z2i = zneigh(nji,nri) + zneigh(nkj,nrj)
!
        xd1 = xneigh(nlk,nrk)
        yd1 = yneigh(nlk,nrk)
        zd1 = zneigh(nlk,nrk)
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
!
!  Set up products for strain
!
      if (lstr) then
        call real1strterm(ndim,xd1,yd1,zd1,0.0_dp,0.0_dp,0.0_dp,dr2ds1,d2r2dx21,d2r2dsdx1,d2r2ds21,.true.)
      endif
!
!  Set second derivative location pointers
!
      lopi = (nfreeatom(n1i).gt.0)
      lopj = (nfreeatom(n1j).gt.0)
!
      lilocal = (atom2local(n1i).gt.0)
      if (lilocal) then
        ix = 3*(atom2local(nfreeatom(n1i)) - 1) + 1
        iy = ix + 1
        iz = iy + 1
      endif
      ixf = 3*(nfreeatom(n1i) - 1) + 1
      iyf = ixf + 1
      izf = iyf + 1
!
      ljlocal = (atom2local(n1j).gt.0)
      if (ljlocal) then
        jx = 3*(atom2local(nfreeatom(n1j)) - 1) + 1
        jy = jx + 1
        jz = jy + 1
      endif
      jxf = 3*(nfreeatom(n1j) - 1) + 1
      jyf = jxf + 1
      jzf = jyf + 1
    endif
!
!  If neither i nor j are being optimised and strain derivatives are not needed then atoms are not valid
!
    if (.not.lopi.and..not.lopj.and..not.lstr) then
      liok = .false.
      ljok = .false.
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
      ind = n1*(n1-1)/2
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!  
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          elseif (n2.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
            else
              llok = .false.
            endif
!
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!               
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          else
!
!  Extended derivatives
!
            ntmp = n2 - ndisttor
            nji = (ntmp - 1)/(mneigh*mneigh) + 1
            ntmp = ntmp - (nji - 1)*mneigh*mneigh
            nkj = (ntmp - 1)/mneigh + 1
            ntmp = ntmp - (nkj - 1)*mneigh
            nlk = ntmp
!
            nrj = neighno(nji,nri)
            if (nrj.eq.0) cycle
            if (nkj.gt.nneigh(nrj)) cycle
            nrk = neighno(nkj,nrj)
            if (nrk.eq.0) cycle
            if (nlk.gt.nneigh(nrk)) cycle
            n2k = nrk
            n2l = neighno(nlk,nrk)
            if (n2l.eq.0) cycle
!
!  Check that atoms are valid
!
            lkok = (nkj.le.nneigh(nrj))
            llok = (nlk.le.nneigh(nrk))
!
            if (liok.and.ljok) then
              xik = xneigh(nji,nri) + xneigh(nkj,nrj) - x2i
              yik = yneigh(nji,nri) + yneigh(nkj,nrj) - y2i
              zik = zneigh(nji,nri) + zneigh(nkj,nrj) - z2i
!
              xd2 = xneigh(nlk,nrk)
              yd2 = yneigh(nlk,nrk)
              zd2 = zneigh(nlk,nrk)
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
!  Set up products for strain
!
            if (lstr) then
              call real1strterm(ndim,xd2,yd2,zd2,0.0_dp,0.0_dp,0.0_dp,dr2ds2,d2r2dx22,d2r2dsdx2,d2r2ds22,.true.)
            endif
!
            lopk = (nfreeatom(n2k).gt.0)
            lopl = (nfreeatom(n2l).gt.0)
!
            lklocal = (atom2local(n2k).gt.0)
            if (lklocal) then
              kx = 3*(atom2local(nfreeatom(n2k)) - 1) + 1
              ky = kx + 1
              kz = ky + 1
            endif
            kxf = 3*(nfreeatom(n2k) - 1) + 1
            kyf = kxf + 1
            kzf = kyf + 1
!
            lllocal = (atom2local(n2l).gt.0)
            if (lllocal) then
              lx = 3*(atom2local(nfreeatom(n2l)) - 1) + 1
              ly = lx + 1
              lz = ly + 1
            endif
            lxf = 3*(nfreeatom(n2l) - 1) + 1
            lyf = lxf + 1
            lzf = lyf + 1
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if ((lopi.or.lopk).and.n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              if (lopi.and.lopk) then
                ldo21 = lilocal
                ldo12 = lklocal
                if (ldo21) then
                  ix1 = ix
                  iy1 = iy
                  iz1 = iz
                  kx2 = kxf
                  ky2 = kyf
                  kz2 = kzf
                endif
                if (ldo12) then
                  ix2 = ixf
                  iy2 = iyf
                  iz2 = izf
                  kx1 = kx
                  ky1 = ky
                  kz1 = kz
                endif
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ldo21 = lilocal
                ldo12 = lilocal
                if (ldo21) then
                  ix1 = ix
                  iy1 = iy
                  iz1 = iz
                  kx2 = ixf
                  ky2 = iyf
                  kz2 = izf
                endif
                if (ldo12) then
                  ix2 = ixf
                  iy2 = iyf
                  iz2 = izf
                  kx1 = ix
                  ky1 = iy
                  kz1 = iz
                endif
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                ldo21 = lklocal
                ldo12 = lklocal
                if (ldo21) then
                  ix1 = kx
                  iy1 = ky
                  iz1 = kz
                  kx2 = kxf
                  ky2 = kyf
                  kz2 = kzf
                endif
                if (ldo12) then
                  ix2 = kxf
                  iy2 = kyf
                  iz2 = kzf
                  kx1 = kx
                  ky1 = ky
                  kz1 = kz
                endif
                if (n1i.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              if (ldo21) then
                derv2(kx2,ix1) = derv2(kx2,ix1) + xd1*xd2*d2(ind)*one1
                derv2(ky2,ix1) = derv2(ky2,ix1) + xd1*yd2*d2(ind)*one1
                derv2(kz2,ix1) = derv2(kz2,ix1) + xd1*zd2*d2(ind)*one1
                derv2(kx2,iy1) = derv2(kx2,iy1) + yd1*xd2*d2(ind)*one1
                derv2(ky2,iy1) = derv2(ky2,iy1) + yd1*yd2*d2(ind)*one1
                derv2(kz2,iy1) = derv2(kz2,iy1) + yd1*zd2*d2(ind)*one1
                derv2(kx2,iz1) = derv2(kx2,iz1) + zd1*xd2*d2(ind)*one1
                derv2(ky2,iz1) = derv2(ky2,iz1) + zd1*yd2*d2(ind)*one1
                derv2(kz2,iz1) = derv2(kz2,iz1) + zd1*zd2*d2(ind)*one1
                if (n1.eq.n2) then
                  derv2(kx2,ix1) = derv2(kx2,ix1) + d1(n2)*one1
                  derv2(ky2,iy1) = derv2(ky2,iy1) + d1(n2)*one1
                  derv2(kz2,iz1) = derv2(kz2,iz1) + d1(n2)*one1
                endif
              endif
              if (.not.lsameijkl.and.ldo12) then
                derv2(ix2,kx1) = derv2(ix2,kx1) + xd1*xd2*d2(ind)*one2
                derv2(ix2,ky1) = derv2(ix2,ky1) + xd1*yd2*d2(ind)*one2
                derv2(ix2,kz1) = derv2(ix2,kz1) + xd1*zd2*d2(ind)*one2
                derv2(iy2,kx1) = derv2(iy2,kx1) + yd1*xd2*d2(ind)*one2
                derv2(iy2,ky1) = derv2(iy2,ky1) + yd1*yd2*d2(ind)*one2
                derv2(iy2,kz1) = derv2(iy2,kz1) + yd1*zd2*d2(ind)*one2
                derv2(iz2,kx1) = derv2(iz2,kx1) + zd1*xd2*d2(ind)*one2
                derv2(iz2,ky1) = derv2(iz2,ky1) + zd1*yd2*d2(ind)*one2
                derv2(iz2,kz1) = derv2(iz2,kz1) + zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(ix2,kx1) = derv2(ix2,kx1) + d1(n2)*one2
                  derv2(iy2,ky1) = derv2(iy2,ky1) + d1(n2)*one2
                  derv2(iz2,kz1) = derv2(iz2,kz1) + d1(n2)*one2
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if ((lopi.or.lopl).and.n1i.ne.n2l) then
              if (lopi.and.lopl) then
                ldo21 = lilocal
                ldo12 = lllocal
                if (ldo21) then
                  ix1 = ix
                  iy1 = iy
                  iz1 = iz
                  lx2 = lxf
                  ly2 = lyf
                  lz2 = lzf
                endif
                if (ldo12) then
                  ix2 = ixf
                  iy2 = iyf
                  iz2 = izf
                  lx1 = lx
                  ly1 = ly
                  lz1 = lz
                endif
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopi) then
                ldo21 = lilocal
                ldo12 = lilocal
                if (ldo21) then
                  ix1 = ix
                  iy1 = iy
                  iz1 = iz
                  lx2 = ixf
                  ly2 = iyf
                  lz2 = izf
                endif
                if (ldo12) then
                  ix2 = ixf
                  iy2 = iyf
                  iz2 = izf
                  lx1 = ix
                  ly1 = iy
                  lz1 = iz
                endif
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                ldo21 = lllocal
                ldo12 = lllocal
                if (ldo21) then
                  ix1 = lx
                  iy1 = ly
                  iz1 = lz
                  lx2 = lxf
                  ly2 = lyf
                  lz2 = lzf
                endif
                if (ldo12) then
                  ix2 = lxf
                  iy2 = lyf
                  iz2 = lzf
                  lx1 = lx
                  ly1 = ly
                  lz1 = lz
                endif
                if (n1i.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              if (ldo21) then
                derv2(lx2,ix1) = derv2(lx2,ix1) - xd1*xd2*d2(ind)*one1
                derv2(ly2,ix1) = derv2(ly2,ix1) - xd1*yd2*d2(ind)*one1
                derv2(lz2,ix1) = derv2(lz2,ix1) - xd1*zd2*d2(ind)*one1
                derv2(lx2,iy1) = derv2(lx2,iy1) - yd1*xd2*d2(ind)*one1
                derv2(ly2,iy1) = derv2(ly2,iy1) - yd1*yd2*d2(ind)*one1
                derv2(lz2,iy1) = derv2(lz2,iy1) - yd1*zd2*d2(ind)*one1
                derv2(lx2,iz1) = derv2(lx2,iz1) - zd1*xd2*d2(ind)*one1
                derv2(ly2,iz1) = derv2(ly2,iz1) - zd1*yd2*d2(ind)*one1
                derv2(lz2,iz1) = derv2(lz2,iz1) - zd1*zd2*d2(ind)*one1
                if (n1.eq.n2) then
                  derv2(lx2,ix1) = derv2(lx2,ix1) - d1(n2)*one1
                  derv2(ly2,iy1) = derv2(ly2,iy1) - d1(n2)*one1
                  derv2(lz2,iz1) = derv2(lz2,iz1) - d1(n2)*one1
                endif
              endif
              if (.not.lsameijkl.and.ldo12) then
                derv2(ix2,lx1) = derv2(ix2,lx1) - xd1*xd2*d2(ind)*one2
                derv2(ix2,ly1) = derv2(ix2,ly1) - xd1*yd2*d2(ind)*one2
                derv2(ix2,lz1) = derv2(ix2,lz1) - xd1*zd2*d2(ind)*one2
                derv2(iy2,lx1) = derv2(iy2,lx1) - yd1*xd2*d2(ind)*one2
                derv2(iy2,ly1) = derv2(iy2,ly1) - yd1*yd2*d2(ind)*one2
                derv2(iy2,lz1) = derv2(iy2,lz1) - yd1*zd2*d2(ind)*one2
                derv2(iz2,lx1) = derv2(iz2,lx1) - zd1*xd2*d2(ind)*one2
                derv2(iz2,ly1) = derv2(iz2,ly1) - zd1*yd2*d2(ind)*one2
                derv2(iz2,lz1) = derv2(iz2,lz1) - zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(ix2,lx1) = derv2(ix2,lx1) - d1(n2)*one2
                  derv2(iy2,ly1) = derv2(iy2,ly1) - d1(n2)*one2
                  derv2(iz2,lz1) = derv2(iz2,lz1) - d1(n2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if ((lopj.or.lopk).and.n1j.ne.n2k) then
              if (lopj.and.lopk) then
                ldo21 = ljlocal
                ldo12 = lklocal
                if (ldo21) then
                  jx1 = jx
                  jy1 = jy
                  jz1 = jz
                  kx2 = kxf
                  ky2 = kyf
                  kz2 = kzf
                endif
                if (ldo12) then
                  jx2 = jxf
                  jy2 = jyf
                  jz2 = jzf
                  kx1 = kx
                  ky1 = ky
                  kz1 = kz
                endif
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                ldo21 = ljlocal
                ldo12 = ljlocal
                if (ldo21) then
                  jx1 = jx
                  jy1 = jy
                  jz1 = jz
                  kx2 = jxf
                  ky2 = jyf
                  kz2 = jzf
                endif
                if (ldo12) then
                  jx2 = jxf
                  jy2 = jyf
                  jz2 = jzf
                  kx1 = jx
                  ky1 = jy
                  kz1 = jz
                endif
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopk) then
                ldo21 = lklocal
                ldo12 = lklocal
                if (ldo21) then
                  jx1 = kx
                  jy1 = ky
                  jz1 = kz
                  kx2 = kxf
                  ky2 = kyf
                  kz2 = kzf
                endif
                if (ldo12) then
                  jx2 = kxf
                  iy2 = kyf
                  iz2 = kzf
                  kx1 = kx
                  ky1 = ky
                  kz1 = kz
                endif
                if (n1j.gt.n2k) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              if (ldo21) then
                derv2(kx2,jx1) = derv2(kx2,jx1) - xd1*xd2*d2(ind)*one1
                derv2(ky2,jx1) = derv2(ky2,jx1) - xd1*yd2*d2(ind)*one1
                derv2(kz2,jx1) = derv2(kz2,jx1) - xd1*zd2*d2(ind)*one1
                derv2(kx2,jy1) = derv2(kx2,jy1) - yd1*xd2*d2(ind)*one1
                derv2(ky2,jy1) = derv2(ky2,jy1) - yd1*yd2*d2(ind)*one1
                derv2(kz2,jy1) = derv2(kz2,jy1) - yd1*zd2*d2(ind)*one1
                derv2(kx2,jz1) = derv2(kx2,jz1) - zd1*xd2*d2(ind)*one1
                derv2(ky2,jz1) = derv2(ky2,jz1) - zd1*yd2*d2(ind)*one1
                derv2(kz2,jz1) = derv2(kz2,jz1) - zd1*zd2*d2(ind)*one1
                if (n1.eq.n2) then
                  derv2(kx2,jx1) = derv2(kx2,jx1) - d1(n2)*one1
                  derv2(ky2,jy1) = derv2(ky2,jy1) - d1(n2)*one1
                  derv2(kz2,jz1) = derv2(kz2,jz1) - d1(n2)*one1
                endif
              endif
              if (.not.lsameijkl.and.ldo12) then
                derv2(jx2,kx1) = derv2(jx2,kx1) - xd1*xd2*d2(ind)*one2
                derv2(jx2,ky1) = derv2(jx2,ky1) - xd1*yd2*d2(ind)*one2
                derv2(jx2,kz1) = derv2(jx2,kz1) - xd1*zd2*d2(ind)*one2
                derv2(jy2,kx1) = derv2(jy2,kx1) - yd1*xd2*d2(ind)*one2
                derv2(jy2,ky1) = derv2(jy2,ky1) - yd1*yd2*d2(ind)*one2
                derv2(jy2,kz1) = derv2(jy2,kz1) - yd1*zd2*d2(ind)*one2
                derv2(jz2,kx1) = derv2(jz2,kx1) - zd1*xd2*d2(ind)*one2
                derv2(jz2,ky1) = derv2(jz2,ky1) - zd1*yd2*d2(ind)*one2
                derv2(jz2,kz1) = derv2(jz2,kz1) - zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(jx2,kx1) = derv2(jx2,kx1) - d1(n2)*one2
                  derv2(jy2,ky1) = derv2(jy2,ky1) - d1(n2)*one2
                  derv2(jz2,kz1) = derv2(jz2,kz1) - d1(n2)*one2
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if ((lopj.or.lopl).and.n1j.ne.n2l) then
              if (lopj.and.lopl) then
                ldo21 = ljlocal
                ldo12 = lllocal
                if (ldo21) then
                  jx1 = jx
                  jy1 = jy
                  jz1 = jz
                  lx2 = lxf
                  ly2 = lyf
                  lz2 = lzf
                endif
                if (ldo12) then
                  jx2 = jxf
                  jy2 = jyf
                  jz2 = jzf
                  lx1 = lx
                  ly1 = ly
                  lz1 = lz
                endif
                one1 = 1.0_dp
                one2 = 1.0_dp
              elseif (lopj) then
                ldo21 = ljlocal
                ldo12 = ljlocal
                if (ldo21) then
                  jx1 = jx
                  jy1 = jy
                  jz1 = jz
                  lx2 = jxf
                  ly2 = jyf
                  lz2 = jzf
                endif
                if (ldo12) then
                  jx2 = jxf
                  jy2 = jyf
                  jz2 = jzf
                  lx1 = jx
                  ly1 = jy
                  lz1 = jz
                endif
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              elseif (lopl) then
                ldo21 = lllocal
                ldo12 = lllocal
                if (ldo21) then
                  jx1 = lx
                  jy1 = ly
                  jz1 = lz
                  lx2 = lxf
                  ly2 = lyf
                  lz2 = lzf
                endif
                if (ldo12) then
                  jx2 = lxf
                  jy2 = lyf
                  jz2 = lzf
                  lx1 = lx
                  ly1 = ly
                  lz1 = lz
                endif
                if (n1j.gt.n2l) then
                  one1 = 1.0_dp
                  one2 = 0.0_dp
                else
                  one1 = 0.0_dp
                  one2 = 1.0_dp
                endif
              endif
!
              if (ldo21) then
                derv2(lx2,jx1) = derv2(lx2,jx1) + xd1*xd2*d2(ind)*one1
                derv2(ly2,jx1) = derv2(ly2,jx1) + xd1*yd2*d2(ind)*one1
                derv2(lz2,jx1) = derv2(lz2,jx1) + xd1*zd2*d2(ind)*one1
                derv2(lx2,jy1) = derv2(lx2,jy1) + yd1*xd2*d2(ind)*one1
                derv2(ly2,jy1) = derv2(ly2,jy1) + yd1*yd2*d2(ind)*one1
                derv2(lz2,jy1) = derv2(lz2,jy1) + yd1*zd2*d2(ind)*one1
                derv2(lx2,jz1) = derv2(lx2,jz1) + zd1*xd2*d2(ind)*one1
                derv2(ly2,jz1) = derv2(ly2,jz1) + zd1*yd2*d2(ind)*one1
                derv2(lz2,jz1) = derv2(lz2,jz1) + zd1*zd2*d2(ind)*one1
                if (n1.eq.n2) then
                  derv2(lx2,jx1) = derv2(lx2,jx1) + d1(n2)*one1
                  derv2(ly2,jy1) = derv2(ly2,jy1) + d1(n2)*one1
                  derv2(lz2,jz1) = derv2(lz2,jz1) + d1(n2)*one1
                endif
              endif
              if (.not.lsameijkl.and.ldo12) then
                derv2(jx2,lx1) = derv2(jx2,lx1) + xd1*xd2*d2(ind)*one2
                derv2(jx2,ly1) = derv2(jx2,ly1) + xd1*yd2*d2(ind)*one2
                derv2(jx2,lz1) = derv2(jx2,lz1) + xd1*zd2*d2(ind)*one2
                derv2(jy2,lx1) = derv2(jy2,lx1) + yd1*xd2*d2(ind)*one2
                derv2(jy2,ly1) = derv2(jy2,ly1) + yd1*yd2*d2(ind)*one2
                derv2(jy2,lz1) = derv2(jy2,lz1) + yd1*zd2*d2(ind)*one2
                derv2(jz2,lx1) = derv2(jz2,lx1) + zd1*xd2*d2(ind)*one2
                derv2(jz2,ly1) = derv2(jz2,ly1) + zd1*yd2*d2(ind)*one2
                derv2(jz2,lz1) = derv2(jz2,lz1) + zd1*zd2*d2(ind)*one2
                if (n1.eq.n2) then
                  derv2(jx2,lx1) = derv2(jx2,lx1) + d1(n2)*one2
                  derv2(jy2,ly1) = derv2(jy2,ly1) + d1(n2)*one2
                  derv2(jz2,lz1) = derv2(jz2,lz1) + d1(n2)*one2
                endif
              endif
            endif
            if (lstr) then
              if (lilocal) then
!
!  Strain - strain second derivatives
!
                do is1 = 1,nstrains
                  ns1 = nstrptr(is1)
                  do is2 = 1,nstrains
                    ns2 = nstrptr(is2)
                    sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds2(ns2)*dr2ds1(ns1)
                    if (n1.ne.n2) then
                      sderv2(is2,is1) = sderv2(is2,is1) + d2(ind)*dr2ds1(ns2)*dr2ds2(ns1)
                    else
                      sderv2(is2,is1) = sderv2(is2,is1) + d1(n2)*d2r2ds22(ns2,ns1)
                    endif
                  enddo
                enddo
              endif
!
!  Internal - strain second derivatives
!
              if (lilocal) then
                do is1 = 1,nstrains
                  ns1 = nstrptr(is1)
                  derv3(ix,is1) = derv3(ix,is1) - xd1*d2(ind)*dr2ds2(ns1)
                  derv3(iy,is1) = derv3(iy,is1) - yd1*d2(ind)*dr2ds2(ns1)
                  derv3(iz,is1) = derv3(iz,is1) - zd1*d2(ind)*dr2ds2(ns1)
                enddo
                if (n1.eq.n2) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(ix,is1) = derv3(ix,is1) - d1(n2)*d2r2dsdx2(ns1,1)
                    derv3(iy,is1) = derv3(iy,is1) - d1(n2)*d2r2dsdx2(ns1,2)
                    derv3(iz,is1) = derv3(iz,is1) - d1(n2)*d2r2dsdx2(ns1,3)
                  enddo
                endif
              endif
              if (ljlocal) then
                do is1 = 1,nstrains
                  ns1 = nstrptr(is1)
                  derv3(jx,is1) = derv3(jx,is1) + xd1*d2(ind)*dr2ds2(ns1)
                  derv3(jy,is1) = derv3(jy,is1) + yd1*d2(ind)*dr2ds2(ns1)
                  derv3(jz,is1) = derv3(jz,is1) + zd1*d2(ind)*dr2ds2(ns1)
                enddo
                if (n1.eq.n2) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(jx,is1) = derv3(jx,is1) + d1(n2)*d2r2dsdx2(ns1,1)
                    derv3(jy,is1) = derv3(jy,is1) + d1(n2)*d2r2dsdx2(ns1,2)
                    derv3(jz,is1) = derv3(jz,is1) + d1(n2)*d2r2dsdx2(ns1,3)
                  enddo
                endif
              endif
              if (n1.ne.n2) then
                if (lklocal) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(kx,is1) = derv3(kx,is1) - xd2*d2(ind)*dr2ds1(ns1)
                    derv3(ky,is1) = derv3(ky,is1) - yd2*d2(ind)*dr2ds1(ns1)
                    derv3(kz,is1) = derv3(kz,is1) - zd2*d2(ind)*dr2ds1(ns1)
                  enddo
                endif
                if (lllocal) then
                  do is1 = 1,nstrains
                    ns1 = nstrptr(is1)
                    derv3(lx,is1) = derv3(lx,is1) + xd2*d2(ind)*dr2ds1(ns1)
                    derv3(ly,is1) = derv3(ly,is1) + yd2*d2(ind)*dr2ds1(ns1)
                    derv3(lz,is1) = derv3(lz,is1) + zd2*d2(ind)*dr2ds1(ns1)
                  enddo
                endif
              endif
            endif
          endif
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('d2adddm')
#endif
!
  return
  end
!
  subroutine d2addpdm(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nREBOatomRptr,d1,d2,xkv,ykv,zkv, &
                      ltorderv,lextended)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian phased second derivative arrays.
!
!  Distributed memory parallel version. 
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  xkv             = x component of K vector
!  ykv             = y component of K vector
!  zkv             = z component of K vector
!  ltorderv        = if true then include torsional derivatives
!
!   1/17 Created from d2addp
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   5/23 Extended derivatives added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use control,        only : lgroupvelocity
  use derivatives
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  real(dp),    intent(in)             :: xkv
  real(dp),    intent(in)             :: ykv
  real(dp),    intent(in)             :: zkv
  logical,     intent(in)             :: ltorderv
  logical,     intent(in)             :: lextended
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: ixf,iyf,izf
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: jxf,jyf,jzf
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: kxf,kyf,kzf
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: lxf,lyf,lzf
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: ndisttor
  integer(i4)                         :: nji
  integer(i4)                         :: nkj
  integer(i4)                         :: nlk
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: nrj
  integer(i4)                         :: nrk
  integer(i4)                         :: ntmp
  complex(dpc)                        :: cdk(3)
  real(dp)                            :: cosik
  real(dp)                            :: cosil
  real(dp)                            :: cosjk
  real(dp)                            :: cosjl
  real(dp)                            :: d2k
  real(dp)                            :: d2ks
  real(dp)                            :: oneik
  real(dp)                            :: oneil
  real(dp)                            :: onejk
  real(dp)                            :: onejl
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: sinik
  real(dp)                            :: sinil
  real(dp)                            :: sinjk
  real(dp)                            :: sinjl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xil
  real(dp)                            :: yil
  real(dp)                            :: zil
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i
  real(dp)                            :: z2i
  logical                             :: lilocal
  logical                             :: lijlocal
  logical                             :: ljlocal
  logical                             :: lklocal
  logical                             :: lkllocal
  logical                             :: lllocal
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lsameijkl
#ifdef TRACE
  call trace_in('d2addpdm')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndisttor = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
    if (lextended) then
      ndist = ndisttor + nneigh(nri)*mneigh*mneigh
    else
      ndist = ndisttor
    endif
  else
    ndist = nneigh2
!
!  If not ltorderv set to a large dummy value that should never be reached
!
    ndisttor = nneigh2 + nneigh(nri)*(mneigh+1)*mneigh + 10_i4*mneigh*mneigh*mneigh
  endif
!
!  Loop over pairs of interatomic distances
!
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!     
!  Check that neighbour numbers are valid
!         
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))       
!
      if (liok.and.ljok) then
!         
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    elseif (n1.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!         
!  Check that atoms are valid
!       
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
!  
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    else
!
!  Extended derivatives
!
      ntmp = n1 - ndisttor
      nji = (ntmp - 1)/(mneigh*mneigh) + 1
      ntmp = ntmp - (nji - 1)*mneigh*mneigh
      nkj = (ntmp - 1)/mneigh + 1
      ntmp = ntmp - (nkj - 1)*mneigh
      nlk = ntmp
!
      nrj = neighno(nji,nri)
      if (nrj.eq.0) cycle
      if (nkj.gt.nneigh(nrj)) cycle
      nrk = neighno(nkj,nrj)
      if (nrk.eq.0) cycle
      if (nlk.gt.nneigh(nrk)) cycle
      n1i = nrk
      n1j = neighno(nlk,nrk)
      if (n1j.eq.0) cycle
! 
!  Check that atoms are valid
!
      liok = (nkj.le.nneigh(nrj))
      ljok = (nlk.le.nneigh(nrk))
!
      if (liok.and.ljok) then
        x2i = xneigh(nji,nri) + xneigh(nkj,nrj)
        y2i = yneigh(nji,nri) + yneigh(nkj,nrj)
        z2i = zneigh(nji,nri) + zneigh(nkj,nrj)
!
        xd1 = xneigh(nlk,nrk)
        yd1 = yneigh(nlk,nrk)
        zd1 = zneigh(nlk,nrk)
      endif
    endif
!  
!  Skip the following section as it is irrelevant if i and j are not valid
!           
    if (liok.and.ljok) then
!
      lilocal = (atom2local(n1i).gt.0)
      ljlocal = (atom2local(n1j).gt.0)
      lijlocal = (lilocal.or.ljlocal)
!
!  Set second derivative location pointers
!
      ixf = 3*(n1i - 1) + 1
      iyf = ixf + 1
      izf = iyf + 1
!
      if (lilocal) then
        ix = 3*(atom2local(n1i) - 1) + 1
        iy = ix + 1
        iz = iy + 1
      endif
!
      jxf = 3*(n1j - 1) + 1
      jyf = jxf + 1
      jzf = jyf + 1
!
      if (ljlocal) then
        jx = 3*(atom2local(n1j) - 1) + 1
        jy = jx + 1
        jz = jy + 1
      endif
!
      ind = n1*(n1-1)/2
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!  
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!               
!  Calculate vector between atoms
!             
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          elseif (n2.le.ndisttor) then
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!               
!  Check that atoms are valid
!             
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(neighno(nji,nri))) 
!               
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          else
!
!  Extended derivatives
!
            ntmp = n2 - ndisttor
            nji = (ntmp - 1)/(mneigh*mneigh) + 1
            ntmp = ntmp - (nji - 1)*mneigh*mneigh
            nkj = (ntmp - 1)/mneigh + 1
            ntmp = ntmp - (nkj - 1)*mneigh
            nlk = ntmp
!
            nrj = neighno(nji,nri)
            if (nrj.eq.0) cycle
            if (nkj.gt.nneigh(nrj)) cycle
            nrk = neighno(nkj,nrj)
            if (nrk.eq.0) cycle
            if (nlk.gt.nneigh(nrk)) cycle
            n2k = nrk
            n2l = neighno(nlk,nrk)
            if (n2l.eq.0) cycle
!
!  Check that atoms are valid
!
            lkok = (nkj.le.nneigh(nrj))
            llok = (nlk.le.nneigh(nrk))
!
            if (liok.and.ljok) then
              xik = xneigh(nji,nri) + xneigh(nkj,nrj) - x2i
              yik = yneigh(nji,nri) + yneigh(nkj,nrj) - y2i
              zik = zneigh(nji,nri) + zneigh(nkj,nrj) - z2i
!
              xd2 = xneigh(nlk,nrk)
              yd2 = yneigh(nlk,nrk)
              zd2 = zneigh(nlk,nrk)
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
            lklocal = (atom2local(n2k).gt.0)
            lllocal = (atom2local(n2l).gt.0)
            lkllocal = (lklocal.or.lllocal)
!
!  Check whether any of i, j, k or l are local before going further
!
            if (lijlocal.or.lkllocal) then
!
!  Complete remaining vectors
!
              xil = xik + xd2
              yil = yik + yd2
              zil = zik + zd2
              xjk = xik - xd1
              yjk = yik - yd1
              zjk = zik - zd1
              xjl = xjk + xd2
              yjl = yjk + yd2
              zjl = zjk + zd2
!
!  Compute distances squared
!
              r2ik = xik*xik + yik*yik + zik*zik
              r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
              kxf = 3*(n2k - 1) + 1
              kyf = kxf + 1
              kzf = kyf + 1
!
              if (lklocal) then
                kx = 3*(atom2local(n2k) - 1) + 1
                ky = kx + 1
                kz = ky + 1
              endif
!
              lxf = 3*(n2l - 1) + 1
              lyf = lxf + 1
              lzf = lyf + 1
!
              if (lllocal) then
                lx = 3*(atom2local(n2l) - 1) + 1
                ly = lx + 1
                lz = ly + 1
              endif
!
              lsameijkl = .false.
              if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
!  Second derivatives : i - k
!
              if (lilocal.or.lklocal) then
                if (n1i.eq.n2k) then
                  oneik = 1.0_dp
                else
                  oneik = 0.0_dp
                endif
                cosik = xkv*xik + ykv*yik + zkv*zik
                sinik = sin(cosik)
                cosik = cos(cosik) - oneik
!
                if (lsameijkl) then
                  cosik = 0.5_dp*cosik
                  sinik = 0.5_dp*sinik
                endif
              endif
!
              if (lilocal) then
                derv2(kxf,ix) = derv2(kxf,ix) + xd1*xd2*d2(ind)*cosik
                derv2(kyf,ix) = derv2(kyf,ix) + xd1*yd2*d2(ind)*cosik
                derv2(kzf,ix) = derv2(kzf,ix) + xd1*zd2*d2(ind)*cosik
                derv2(kxf,iy) = derv2(kxf,iy) + yd1*xd2*d2(ind)*cosik
                derv2(kyf,iy) = derv2(kyf,iy) + yd1*yd2*d2(ind)*cosik
                derv2(kzf,iy) = derv2(kzf,iy) + yd1*zd2*d2(ind)*cosik
                derv2(kxf,iz) = derv2(kxf,iz) + zd1*xd2*d2(ind)*cosik
                derv2(kyf,iz) = derv2(kyf,iz) + zd1*yd2*d2(ind)*cosik
                derv2(kzf,iz) = derv2(kzf,iz) + zd1*zd2*d2(ind)*cosik
!
                dervi(kxf,ix) = dervi(kxf,ix) + xd1*xd2*d2(ind)*sinik
                dervi(kyf,ix) = dervi(kyf,ix) + xd1*yd2*d2(ind)*sinik
                dervi(kzf,ix) = dervi(kzf,ix) + xd1*zd2*d2(ind)*sinik
                dervi(kxf,iy) = dervi(kxf,iy) + yd1*xd2*d2(ind)*sinik
                dervi(kyf,iy) = dervi(kyf,iy) + yd1*yd2*d2(ind)*sinik
                dervi(kzf,iy) = dervi(kzf,iy) + yd1*zd2*d2(ind)*sinik
                dervi(kxf,iz) = dervi(kxf,iz) + zd1*xd2*d2(ind)*sinik
                dervi(kyf,iz) = dervi(kyf,iz) + zd1*yd2*d2(ind)*sinik
                dervi(kzf,iz) = dervi(kzf,iz) + zd1*zd2*d2(ind)*sinik
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosik
                  d2ks = d2(ind)*sinik
                  cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,kxf,ix) = derv2dk(1:3,kxf,ix) + xd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,ix) = derv2dk(1:3,kyf,ix) + xd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,ix) = derv2dk(1:3,kzf,ix) + xd1*zd2*cdk(1:3)
                  derv2dk(1:3,kxf,iy) = derv2dk(1:3,kxf,iy) + yd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,iy) = derv2dk(1:3,kyf,iy) + yd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,iy) = derv2dk(1:3,kzf,iy) + yd1*zd2*cdk(1:3)
                  derv2dk(1:3,kxf,iz) = derv2dk(1:3,kxf,iz) + zd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,iz) = derv2dk(1:3,kyf,iz) + zd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,iz) = derv2dk(1:3,kzf,iz) + zd1*zd2*cdk(1:3)
                endif
!
                if (n1.eq.n2) then
                  derv2(kxf,ix) = derv2(kxf,ix) + d1(n2)*cosik
                  derv2(kyf,iy) = derv2(kyf,iy) + d1(n2)*cosik
                  derv2(kzf,iz) = derv2(kzf,iz) + d1(n2)*cosik
                  dervi(kxf,ix) = dervi(kxf,ix) + d1(n2)*sinik
                  dervi(kyf,iy) = dervi(kyf,iy) + d1(n2)*sinik
                  dervi(kzf,iz) = dervi(kzf,iz) + d1(n2)*sinik
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(n2)*cosik
                    d2ks = d1(n2)*sinik
                    cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,kxf,ix) = derv2dk(1:3,kxf,ix) + cdk(1:3)
                    derv2dk(1:3,kyf,iy) = derv2dk(1:3,kyf,iy) + cdk(1:3)
                    derv2dk(1:3,kzf,iz) = derv2dk(1:3,kzf,iz) + cdk(1:3)
                  endif
                endif
              endif
!
              if (lklocal) then
                derv2(ixf,kx) = derv2(ixf,kx) + xd1*xd2*d2(ind)*cosik
                derv2(ixf,ky) = derv2(ixf,ky) + xd1*yd2*d2(ind)*cosik
                derv2(ixf,kz) = derv2(ixf,kz) + xd1*zd2*d2(ind)*cosik
                derv2(iyf,kx) = derv2(iyf,kx) + yd1*xd2*d2(ind)*cosik
                derv2(iyf,ky) = derv2(iyf,ky) + yd1*yd2*d2(ind)*cosik
                derv2(iyf,kz) = derv2(iyf,kz) + yd1*zd2*d2(ind)*cosik
                derv2(izf,kx) = derv2(izf,kx) + zd1*xd2*d2(ind)*cosik
                derv2(izf,ky) = derv2(izf,ky) + zd1*yd2*d2(ind)*cosik
                derv2(izf,kz) = derv2(izf,kz) + zd1*zd2*d2(ind)*cosik
!
                dervi(ixf,kx) = dervi(ixf,kx) - xd1*xd2*d2(ind)*sinik
                dervi(ixf,ky) = dervi(ixf,ky) - xd1*yd2*d2(ind)*sinik
                dervi(ixf,kz) = dervi(ixf,kz) - xd1*zd2*d2(ind)*sinik
                dervi(iyf,kx) = dervi(iyf,kx) - yd1*xd2*d2(ind)*sinik
                dervi(iyf,ky) = dervi(iyf,ky) - yd1*yd2*d2(ind)*sinik
                dervi(iyf,kz) = dervi(iyf,kz) - yd1*zd2*d2(ind)*sinik
                dervi(izf,kx) = dervi(izf,kx) - zd1*xd2*d2(ind)*sinik
                dervi(izf,ky) = dervi(izf,ky) - zd1*yd2*d2(ind)*sinik
                dervi(izf,kz) = dervi(izf,kz) - zd1*zd2*d2(ind)*sinik
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosik
                  d2ks = d2(ind)*sinik
                  cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,ixf,kx) = derv2dk(1:3,ixf,kx) + xd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,ixf,ky) = derv2dk(1:3,ixf,ky) + xd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,ixf,kz) = derv2dk(1:3,ixf,kz) + xd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,kx) = derv2dk(1:3,iyf,kx) + yd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,ky) = derv2dk(1:3,iyf,ky) + yd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,kz) = derv2dk(1:3,iyf,kz) + yd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,kx) = derv2dk(1:3,izf,kx) + zd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,ky) = derv2dk(1:3,izf,ky) + zd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,kz) = derv2dk(1:3,izf,kz) + zd1*zd2*conjg(cdk(1:3))
                endif
!
                if (n1.eq.n2) then
                  derv2(ixf,kx) = derv2(ixf,kx) + d1(n2)*cosik
                  derv2(iyf,ky) = derv2(iyf,ky) + d1(n2)*cosik
                  derv2(izf,kz) = derv2(izf,kz) + d1(n2)*cosik
                  dervi(ixf,kx) = dervi(ixf,kx) - d1(n2)*sinik
                  dervi(iyf,ky) = dervi(iyf,ky) - d1(n2)*sinik
                  dervi(izf,kz) = dervi(izf,kz) - d1(n2)*sinik
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(n2)*cosik
                    d2ks = d1(n2)*sinik
                    cdk(1) = dcmplx(d2k*xik,d2ks*xik)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yik,d2ks*yik)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zik,d2ks*zik)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,ixf,kx) = derv2dk(1:3,ixf,kx) + conjg(cdk(1:3))
                    derv2dk(1:3,iyf,ky) = derv2dk(1:3,iyf,ky) + conjg(cdk(1:3))
                    derv2dk(1:3,izf,kz) = derv2dk(1:3,izf,kz) + conjg(cdk(1:3))
                  endif
                endif
              endif
!
!  Second derivatives : i - l
!
              if (lilocal.or.lllocal) then
                if (n1i.eq.n2l) then
                  oneil = 1.0_dp
                else
                  oneil = 0.0_dp
                endif
                cosil = xkv*xil + ykv*yil + zkv*zil
                sinil = sin(cosil)
                cosil = cos(cosil) - oneil
!
                if (lsameijkl) then
                  cosil = 0.5_dp*cosil
                  sinil = 0.5_dp*sinil
                endif
              endif
!
              if (lilocal) then
                derv2(lxf,ix) = derv2(lxf,ix) - xd1*xd2*d2(ind)*cosil
                derv2(lyf,ix) = derv2(lyf,ix) - xd1*yd2*d2(ind)*cosil
                derv2(lzf,ix) = derv2(lzf,ix) - xd1*zd2*d2(ind)*cosil
                derv2(lxf,iy) = derv2(lxf,iy) - yd1*xd2*d2(ind)*cosil
                derv2(lyf,iy) = derv2(lyf,iy) - yd1*yd2*d2(ind)*cosil
                derv2(lzf,iy) = derv2(lzf,iy) - yd1*zd2*d2(ind)*cosil
                derv2(lxf,iz) = derv2(lxf,iz) - zd1*xd2*d2(ind)*cosil
                derv2(lyf,iz) = derv2(lyf,iz) - zd1*yd2*d2(ind)*cosil
                derv2(lzf,iz) = derv2(lzf,iz) - zd1*zd2*d2(ind)*cosil
!
                dervi(lxf,ix) = dervi(lxf,ix) - xd1*xd2*d2(ind)*sinil
                dervi(lyf,ix) = dervi(lyf,ix) - xd1*yd2*d2(ind)*sinil
                dervi(lzf,ix) = dervi(lzf,ix) - xd1*zd2*d2(ind)*sinil
                dervi(lxf,iy) = dervi(lxf,iy) - yd1*xd2*d2(ind)*sinil
                dervi(lyf,iy) = dervi(lyf,iy) - yd1*yd2*d2(ind)*sinil
                dervi(lzf,iy) = dervi(lzf,iy) - yd1*zd2*d2(ind)*sinil
                dervi(lxf,iz) = dervi(lxf,iz) - zd1*xd2*d2(ind)*sinil
                dervi(lyf,iz) = dervi(lyf,iz) - zd1*yd2*d2(ind)*sinil
                dervi(lzf,iz) = dervi(lzf,iz) - zd1*zd2*d2(ind)*sinil
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosil
                  d2ks = d2(ind)*sinil
                  cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,lxf,ix) = derv2dk(1:3,lxf,ix) - xd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,ix) = derv2dk(1:3,lyf,ix) - xd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,ix) = derv2dk(1:3,lzf,ix) - xd1*zd2*cdk(1:3)
                  derv2dk(1:3,lxf,iy) = derv2dk(1:3,lxf,iy) - yd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,iy) = derv2dk(1:3,lyf,iy) - yd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,iy) = derv2dk(1:3,lzf,iy) - yd1*zd2*cdk(1:3)
                  derv2dk(1:3,lxf,iz) = derv2dk(1:3,lxf,iz) - zd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,iz) = derv2dk(1:3,lyf,iz) - zd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,iz) = derv2dk(1:3,lzf,iz) - zd1*zd2*cdk(1:3)
                endif
!
                if (n1.eq.n2) then
                  derv2(lxf,ix) = derv2(lxf,ix) - d1(n2)*cosil
                  derv2(lyf,iy) = derv2(lyf,iy) - d1(n2)*cosil
                  derv2(lzf,iz) = derv2(lzf,iz) - d1(n2)*cosil
                  dervi(lxf,ix) = dervi(lxf,ix) - d1(n2)*sinil
                  dervi(lyf,iy) = dervi(lyf,iy) - d1(n2)*sinil
                  dervi(lzf,iz) = dervi(lzf,iz) - d1(n2)*sinil
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(n2)*cosil
                    d2ks = d1(n2)*sinil
                    cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,lxf,ix) = derv2dk(1:3,lxf,ix) - cdk(1:3)
                    derv2dk(1:3,lyf,iy) = derv2dk(1:3,lyf,iy) - cdk(1:3)
                    derv2dk(1:3,lzf,iz) = derv2dk(1:3,lzf,iz) - cdk(1:3)
                  endif
                endif
              endif
!
              if (lllocal) then
                derv2(ixf,lx) = derv2(ixf,lx) - xd1*xd2*d2(ind)*cosil
                derv2(ixf,ly) = derv2(ixf,ly) - xd1*yd2*d2(ind)*cosil
                derv2(ixf,lz) = derv2(ixf,lz) - xd1*zd2*d2(ind)*cosil
                derv2(iyf,lx) = derv2(iyf,lx) - yd1*xd2*d2(ind)*cosil
                derv2(iyf,ly) = derv2(iyf,ly) - yd1*yd2*d2(ind)*cosil
                derv2(iyf,lz) = derv2(iyf,lz) - yd1*zd2*d2(ind)*cosil
                derv2(izf,lx) = derv2(izf,lx) - zd1*xd2*d2(ind)*cosil
                derv2(izf,ly) = derv2(izf,ly) - zd1*yd2*d2(ind)*cosil
                derv2(izf,lz) = derv2(izf,lz) - zd1*zd2*d2(ind)*cosil
!
                dervi(ixf,lx) = dervi(ixf,lx) + xd1*xd2*d2(ind)*sinil
                dervi(ixf,ly) = dervi(ixf,ly) + xd1*yd2*d2(ind)*sinil
                dervi(ixf,lz) = dervi(ixf,lz) + xd1*zd2*d2(ind)*sinil
                dervi(iyf,lx) = dervi(iyf,lx) + yd1*xd2*d2(ind)*sinil
                dervi(iyf,ly) = dervi(iyf,ly) + yd1*yd2*d2(ind)*sinil
                dervi(iyf,lz) = dervi(iyf,lz) + yd1*zd2*d2(ind)*sinil
                dervi(izf,lx) = dervi(izf,lx) + zd1*xd2*d2(ind)*sinil
                dervi(izf,ly) = dervi(izf,ly) + zd1*yd2*d2(ind)*sinil
                dervi(izf,lz) = dervi(izf,lz) + zd1*zd2*d2(ind)*sinil
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosil
                  d2ks = d2(ind)*sinil
                  cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,ixf,lx) = derv2dk(1:3,ixf,lx) - xd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,ixf,ly) = derv2dk(1:3,ixf,ly) - xd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,ixf,lz) = derv2dk(1:3,ixf,lz) - xd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,lx) = derv2dk(1:3,iyf,lx) - yd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,ly) = derv2dk(1:3,iyf,ly) - yd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,iyf,lz) = derv2dk(1:3,iyf,lz) - yd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,lx) = derv2dk(1:3,izf,lx) - zd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,ly) = derv2dk(1:3,izf,ly) - zd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,izf,lz) = derv2dk(1:3,izf,lz) - zd1*zd2*conjg(cdk(1:3))
                endif
!
                if (n1.eq.n2) then
                  derv2(ixf,lx) = derv2(ixf,lx) - d1(n2)*cosil
                  derv2(iyf,ly) = derv2(iyf,ly) - d1(n2)*cosil
                  derv2(izf,lz) = derv2(izf,lz) - d1(n2)*cosil
                  dervi(ixf,lx) = dervi(ixf,lx) + d1(n2)*sinil
                  dervi(iyf,ly) = dervi(iyf,ly) + d1(n2)*sinil
                  dervi(izf,lz) = dervi(izf,lz) + d1(n2)*sinil
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(n2)*cosil
                    d2ks = d1(n2)*sinil
                    cdk(1) = dcmplx(d2k*xil,d2ks*xil)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yil,d2ks*yil)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zil,d2ks*zil)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,ixf,lx) = derv2dk(1:3,ixf,lx) - conjg(cdk(1:3))
                    derv2dk(1:3,iyf,ly) = derv2dk(1:3,iyf,ly) - conjg(cdk(1:3))
                    derv2dk(1:3,izf,lz) = derv2dk(1:3,izf,lz) - conjg(cdk(1:3))
                  endif
                endif
              endif
!
!  Second derivatives : j - k
!
              if (ljlocal.or.lklocal) then
                if (n1j.eq.n2k) then
                  onejk = 1.0_dp
                else
                  onejk = 0.0_dp
                endif
                cosjk = xkv*xjk + ykv*yjk + zkv*zjk
                sinjk = sin(cosjk)
                cosjk = cos(cosjk) - onejk
!
                if (lsameijkl) then
                  cosjk = 0.5_dp*cosjk
                  sinjk = 0.5_dp*sinjk
                endif
              endif
!
              if (ljlocal) then
                derv2(kxf,jx) = derv2(kxf,jx) - xd1*xd2*d2(ind)*cosjk
                derv2(kyf,jx) = derv2(kyf,jx) - xd1*yd2*d2(ind)*cosjk
                derv2(kzf,jx) = derv2(kzf,jx) - xd1*zd2*d2(ind)*cosjk
                derv2(kxf,jy) = derv2(kxf,jy) - yd1*xd2*d2(ind)*cosjk
                derv2(kyf,jy) = derv2(kyf,jy) - yd1*yd2*d2(ind)*cosjk
                derv2(kzf,jy) = derv2(kzf,jy) - yd1*zd2*d2(ind)*cosjk
                derv2(kxf,jz) = derv2(kxf,jz) - zd1*xd2*d2(ind)*cosjk
                derv2(kyf,jz) = derv2(kyf,jz) - zd1*yd2*d2(ind)*cosjk
                derv2(kzf,jz) = derv2(kzf,jz) - zd1*zd2*d2(ind)*cosjk
!
                dervi(kxf,jx) = dervi(kxf,jx) - xd1*xd2*d2(ind)*sinjk
                dervi(kyf,jx) = dervi(kyf,jx) - xd1*yd2*d2(ind)*sinjk
                dervi(kzf,jx) = dervi(kzf,jx) - xd1*zd2*d2(ind)*sinjk
                dervi(kxf,jy) = dervi(kxf,jy) - yd1*xd2*d2(ind)*sinjk
                dervi(kyf,jy) = dervi(kyf,jy) - yd1*yd2*d2(ind)*sinjk
                dervi(kzf,jy) = dervi(kzf,jy) - yd1*zd2*d2(ind)*sinjk
                dervi(kxf,jz) = dervi(kxf,jz) - zd1*xd2*d2(ind)*sinjk
                dervi(kyf,jz) = dervi(kyf,jz) - zd1*yd2*d2(ind)*sinjk
                dervi(kzf,jz) = dervi(kzf,jz) - zd1*zd2*d2(ind)*sinjk
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosjk
                  d2ks = d2(ind)*sinjk
                  cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,kxf,jx) = derv2dk(1:3,kxf,jx) - xd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,jx) = derv2dk(1:3,kyf,jx) - xd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,jx) = derv2dk(1:3,kzf,jx) - xd1*zd2*cdk(1:3)
                  derv2dk(1:3,kxf,jy) = derv2dk(1:3,kxf,jy) - yd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,jy) = derv2dk(1:3,kyf,jy) - yd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,jy) = derv2dk(1:3,kzf,jy) - yd1*zd2*cdk(1:3)
                  derv2dk(1:3,kxf,jz) = derv2dk(1:3,kxf,jz) - zd1*xd2*cdk(1:3)
                  derv2dk(1:3,kyf,jz) = derv2dk(1:3,kyf,jz) - zd1*yd2*cdk(1:3)
                  derv2dk(1:3,kzf,jz) = derv2dk(1:3,kzf,jz) - zd1*zd2*cdk(1:3)
                endif
!
                if (n1.eq.n2) then
                  derv2(kxf,jx) = derv2(kxf,jx) - d1(n2)*cosjk
                  derv2(kyf,jy) = derv2(kyf,jy) - d1(n2)*cosjk
                  derv2(kzf,jz) = derv2(kzf,jz) - d1(n2)*cosjk
                  dervi(kxf,jx) = dervi(kxf,jx) - d1(n2)*sinjk
                  dervi(kyf,jy) = dervi(kyf,jy) - d1(n2)*sinjk
                  dervi(kzf,jz) = dervi(kzf,jz) - d1(n2)*sinjk
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(n2)*cosjk
                    d2ks = d1(n2)*sinjk
                    cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,kxf,jx) = derv2dk(1:3,kxf,jx) - cdk(1:3)
                    derv2dk(1:3,kyf,jy) = derv2dk(1:3,kyf,jy) - cdk(1:3)
                    derv2dk(1:3,kzf,jz) = derv2dk(1:3,kzf,jz) - cdk(1:3)
                  endif
                endif
              endif
!
              if (lklocal) then
                derv2(jxf,kx) = derv2(jxf,kx) - xd1*xd2*d2(ind)*cosjk
                derv2(jxf,ky) = derv2(jxf,ky) - xd1*yd2*d2(ind)*cosjk
                derv2(jxf,kz) = derv2(jxf,kz) - xd1*zd2*d2(ind)*cosjk
                derv2(jyf,kx) = derv2(jyf,kx) - yd1*xd2*d2(ind)*cosjk
                derv2(jyf,ky) = derv2(jyf,ky) - yd1*yd2*d2(ind)*cosjk
                derv2(jyf,kz) = derv2(jyf,kz) - yd1*zd2*d2(ind)*cosjk
                derv2(jzf,kx) = derv2(jzf,kx) - zd1*xd2*d2(ind)*cosjk
                derv2(jzf,ky) = derv2(jzf,ky) - zd1*yd2*d2(ind)*cosjk
                derv2(jzf,kz) = derv2(jzf,kz) - zd1*zd2*d2(ind)*cosjk
!
                dervi(jxf,kx) = dervi(jxf,kx) + xd1*xd2*d2(ind)*sinjk
                dervi(jxf,ky) = dervi(jxf,ky) + xd1*yd2*d2(ind)*sinjk
                dervi(jxf,kz) = dervi(jxf,kz) + xd1*zd2*d2(ind)*sinjk
                dervi(jyf,kx) = dervi(jyf,kx) + yd1*xd2*d2(ind)*sinjk
                dervi(jyf,ky) = dervi(jyf,ky) + yd1*yd2*d2(ind)*sinjk
                dervi(jyf,kz) = dervi(jyf,kz) + yd1*zd2*d2(ind)*sinjk
                dervi(jzf,kx) = dervi(jzf,kx) + zd1*xd2*d2(ind)*sinjk
                dervi(jzf,ky) = dervi(jzf,ky) + zd1*yd2*d2(ind)*sinjk
                dervi(jzf,kz) = dervi(jzf,kz) + zd1*zd2*d2(ind)*sinjk
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosjk
                  d2ks = d2(ind)*sinjk
                  cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,jxf,kx) = derv2dk(1:3,jxf,kx) - xd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jxf,ky) = derv2dk(1:3,jxf,ky) - xd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jxf,kz) = derv2dk(1:3,jxf,kz) - xd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,kx) = derv2dk(1:3,jyf,kx) - yd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,ky) = derv2dk(1:3,jyf,ky) - yd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,kz) = derv2dk(1:3,jyf,kz) - yd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,kx) = derv2dk(1:3,jzf,kx) - zd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,ky) = derv2dk(1:3,jzf,ky) - zd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,kz) = derv2dk(1:3,jzf,kz) - zd1*zd2*conjg(cdk(1:3))
                endif
!
                if (n1.eq.n2) then
                  derv2(jxf,kx) = derv2(jxf,kx) - d1(n2)*cosjk
                  derv2(jyf,ky) = derv2(jyf,ky) - d1(n2)*cosjk
                  derv2(jzf,kz) = derv2(jzf,kz) - d1(n2)*cosjk
                  dervi(jxf,kx) = dervi(jxf,kx) + d1(n2)*sinjk
                  dervi(jyf,ky) = dervi(jyf,ky) + d1(n2)*sinjk
                  dervi(jzf,kz) = dervi(jzf,kz) + d1(n2)*sinjk
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(n2)*cosjk
                    d2ks = d1(n2)*sinjk
                    cdk(1) = dcmplx(d2k*xjk,d2ks*xjk)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yjk,d2ks*yjk)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zjk,d2ks*zjk)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,jxf,kx) = derv2dk(1:3,jxf,kx) - conjg(cdk(1:3))
                    derv2dk(1:3,jyf,ky) = derv2dk(1:3,jyf,ky) - conjg(cdk(1:3))
                    derv2dk(1:3,jzf,kz) = derv2dk(1:3,jzf,kz) - conjg(cdk(1:3))
                  endif
                endif
              endif
!
!  Second derivatives : j - l
!
              if (ljlocal.or.lllocal) then
                if (n1j.eq.n2l) then
                  onejl = 1.0_dp
                else
                  onejl = 0.0_dp
                endif
                cosjl = xkv*xjl + ykv*yjl + zkv*zjl
                sinjl = sin(cosjl)
                cosjl = cos(cosjl) - onejl
!
                if (lsameijkl) then
                  cosjl = 0.5_dp*cosjl
                  sinjl = 0.5_dp*sinjl
                endif
              endif
!
              if (ljlocal) then
                derv2(lxf,jx) = derv2(lxf,jx) + xd1*xd2*d2(ind)*cosjl
                derv2(lyf,jx) = derv2(lyf,jx) + xd1*yd2*d2(ind)*cosjl
                derv2(lzf,jx) = derv2(lzf,jx) + xd1*zd2*d2(ind)*cosjl
                derv2(lxf,jy) = derv2(lxf,jy) + yd1*xd2*d2(ind)*cosjl
                derv2(lyf,jy) = derv2(lyf,jy) + yd1*yd2*d2(ind)*cosjl
                derv2(lzf,jy) = derv2(lzf,jy) + yd1*zd2*d2(ind)*cosjl
                derv2(lxf,jz) = derv2(lxf,jz) + zd1*xd2*d2(ind)*cosjl
                derv2(lyf,jz) = derv2(lyf,jz) + zd1*yd2*d2(ind)*cosjl
                derv2(lzf,jz) = derv2(lzf,jz) + zd1*zd2*d2(ind)*cosjl
!
                dervi(lxf,jx) = dervi(lxf,jx) + xd1*xd2*d2(ind)*sinjl
                dervi(lyf,jx) = dervi(lyf,jx) + xd1*yd2*d2(ind)*sinjl
                dervi(lzf,jx) = dervi(lzf,jx) + xd1*zd2*d2(ind)*sinjl
                dervi(lxf,jy) = dervi(lxf,jy) + yd1*xd2*d2(ind)*sinjl
                dervi(lyf,jy) = dervi(lyf,jy) + yd1*yd2*d2(ind)*sinjl
                dervi(lzf,jy) = dervi(lzf,jy) + yd1*zd2*d2(ind)*sinjl
                dervi(lxf,jz) = dervi(lxf,jz) + zd1*xd2*d2(ind)*sinjl
                dervi(lyf,jz) = dervi(lyf,jz) + zd1*yd2*d2(ind)*sinjl
                dervi(lzf,jz) = dervi(lzf,jz) + zd1*zd2*d2(ind)*sinjl
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosjl
                  d2ks = d2(ind)*sinjl
                  cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,lxf,jx) = derv2dk(1:3,lxf,jx) + xd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,jx) = derv2dk(1:3,lyf,jx) + xd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,jx) = derv2dk(1:3,lzf,jx) + xd1*zd2*cdk(1:3)
                  derv2dk(1:3,lxf,jy) = derv2dk(1:3,lxf,jy) + yd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,jy) = derv2dk(1:3,lyf,jy) + yd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,jy) = derv2dk(1:3,lzf,jy) + yd1*zd2*cdk(1:3)
                  derv2dk(1:3,lxf,jz) = derv2dk(1:3,lxf,jz) + zd1*xd2*cdk(1:3)
                  derv2dk(1:3,lyf,jz) = derv2dk(1:3,lyf,jz) + zd1*yd2*cdk(1:3)
                  derv2dk(1:3,lzf,jz) = derv2dk(1:3,lzf,jz) + zd1*zd2*cdk(1:3)
                endif
!
                if (n1.eq.n2) then
                  derv2(lxf,jx) = derv2(lxf,jx) + d1(n2)*cosjl
                  derv2(lyf,jy) = derv2(lyf,jy) + d1(n2)*cosjl
                  derv2(lzf,jz) = derv2(lzf,jz) + d1(n2)*cosjl
                  dervi(lxf,jx) = dervi(lxf,jx) + d1(n2)*sinjl
                  dervi(lyf,jy) = dervi(lyf,jy) + d1(n2)*sinjl
                  dervi(lzf,jz) = dervi(lzf,jz) + d1(n2)*sinjl
!
                  if (lgroupvelocity) then
!
!  Group velocities
!
                    d2k  = d1(n2)*cosjl
                    d2ks = d1(n2)*sinjl
                    cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!
                    derv2dk(1:3,lxf,jx) = derv2dk(1:3,lxf,jx) + cdk(1:3)
                    derv2dk(1:3,lyf,jy) = derv2dk(1:3,lyf,jy) + cdk(1:3)
                    derv2dk(1:3,lzf,jz) = derv2dk(1:3,lzf,jz) + cdk(1:3)
                  endif
                endif
              endif
!
              if (lllocal) then
                derv2(jxf,lx) = derv2(jxf,lx) + xd1*xd2*d2(ind)*cosjl
                derv2(jxf,ly) = derv2(jxf,ly) + xd1*yd2*d2(ind)*cosjl
                derv2(jxf,lz) = derv2(jxf,lz) + xd1*zd2*d2(ind)*cosjl
                derv2(jyf,lx) = derv2(jyf,lx) + yd1*xd2*d2(ind)*cosjl
                derv2(jyf,ly) = derv2(jyf,ly) + yd1*yd2*d2(ind)*cosjl
                derv2(jyf,lz) = derv2(jyf,lz) + yd1*zd2*d2(ind)*cosjl
                derv2(jzf,lx) = derv2(jzf,lx) + zd1*xd2*d2(ind)*cosjl
                derv2(jzf,ly) = derv2(jzf,ly) + zd1*yd2*d2(ind)*cosjl
                derv2(jzf,lz) = derv2(jzf,lz) + zd1*zd2*d2(ind)*cosjl
!
                dervi(jxf,lx) = dervi(jxf,lx) - xd1*xd2*d2(ind)*sinjl
                dervi(jxf,ly) = dervi(jxf,ly) - xd1*yd2*d2(ind)*sinjl
                dervi(jxf,lz) = dervi(jxf,lz) - xd1*zd2*d2(ind)*sinjl
                dervi(jyf,lx) = dervi(jyf,lx) - yd1*xd2*d2(ind)*sinjl
                dervi(jyf,ly) = dervi(jyf,ly) - yd1*yd2*d2(ind)*sinjl
                dervi(jyf,lz) = dervi(jyf,lz) - yd1*zd2*d2(ind)*sinjl
                dervi(jzf,lx) = dervi(jzf,lx) - zd1*xd2*d2(ind)*sinjl
                dervi(jzf,ly) = dervi(jzf,ly) - zd1*yd2*d2(ind)*sinjl
                dervi(jzf,lz) = dervi(jzf,lz) - zd1*zd2*d2(ind)*sinjl
!
                if (lgroupvelocity) then
!
!  Group velocities
!
                  d2k  = d2(ind)*cosjl
                  d2ks = d2(ind)*sinjl
                  cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                  cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                  cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!
                  derv2dk(1:3,jxf,lx) = derv2dk(1:3,jxf,lx) + xd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jxf,ly) = derv2dk(1:3,jxf,ly) + xd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jxf,lz) = derv2dk(1:3,jxf,lz) + xd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,lx) = derv2dk(1:3,jyf,lx) + yd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,ly) = derv2dk(1:3,jyf,ly) + yd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jyf,lz) = derv2dk(1:3,jyf,lz) + yd1*zd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,lx) = derv2dk(1:3,jzf,lx) + zd1*xd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,ly) = derv2dk(1:3,jzf,ly) + zd1*yd2*conjg(cdk(1:3))
                  derv2dk(1:3,jzf,lz) = derv2dk(1:3,jzf,lz) + zd1*zd2*conjg(cdk(1:3))
                endif
!
                if (n1.eq.n2) then
                  derv2(jxf,lx) = derv2(jxf,lx) + d1(n2)*cosjl
                  derv2(jyf,ly) = derv2(jyf,ly) + d1(n2)*cosjl
                  derv2(jzf,lz) = derv2(jzf,lz) + d1(n2)*cosjl
                  dervi(jxf,lx) = dervi(jxf,lx) - d1(n2)*sinjl
                  dervi(jyf,ly) = dervi(jyf,ly) - d1(n2)*sinjl
                  dervi(jzf,lz) = dervi(jzf,lz) - d1(n2)*sinjl
!
                  if (lgroupvelocity) then
!
!  Group velocities
!           
                    d2k  = d1(n2)*cosjl
                    d2ks = d1(n2)*sinjl
                    cdk(1) = dcmplx(d2k*xjl,d2ks*xjl)*dcmplx(0.0_dp,1.0_dp)
                    cdk(2) = dcmplx(d2k*yjl,d2ks*yjl)*dcmplx(0.0_dp,1.0_dp)
                    cdk(3) = dcmplx(d2k*zjl,d2ks*zjl)*dcmplx(0.0_dp,1.0_dp)
!           
                    derv2dk(1:3,jxf,lx) = derv2dk(1:3,jxf,lx) + conjg(cdk(1:3))
                    derv2dk(1:3,jyf,ly) = derv2dk(1:3,jyf,ly) + conjg(cdk(1:3))
                    derv2dk(1:3,jzf,lz) = derv2dk(1:3,jzf,lz) + conjg(cdk(1:3))
                  endif
                endif
              endif
            endif
          endif
        endif
      enddo
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addpdm')
#endif
!
  return
  end
!
  subroutine d2addddm(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nREBOatomRptr, &
                      nregion1,d1,d2,ltorderv)
!
!  Utility to add the pairwise derivatives for a group of atoms
!  on to the global Cartesian second derivative array. Defect
!  calculation version.
!
!  Distributed memory parallel version.
!
!  On entry : 
!
!  i               = atom number of i to whom others are connected
!  maxneigh        = left-hand dimension of 2-D arrays
!  mneigh          = value of maxneigh actually used for building index numbers
!  nneigh          = no. of neighbours of atoms
!  neighno         = atom number of each neighbour
!  xneigh          = x component of vector to neighbour of atoms
!  yneigh          = y component of vector to neighbour of atoms
!  zneigh          = z component of vector to neighbour of atoms
!  nfreeatom       = pointer from atom to free atom to move
!  nREBOatomRptr   = pointer from atoms to position in REBO set
!  nregion1        = number of ions in region 1
!  d1              = array of first derivatives w.r.t. interatomic
!                    distances
!  d2              = array of first derivatives w.r.t. interatomic
!                    distances
!  ltorderv        = if true then include torsional derivatives
!
!   5/17 Created from d2addd
!   2/18 Trace added
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!  10/22 mneigh added to arguments and used for indices
!  12/22 lsameijkl now has condition on n1 = n2 added
!   5/23 Derivative tolerances as global variables introduced
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
  use brennerdata
  use derivatives
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)             :: i
  integer(i4), intent(in)             :: maxneigh
  integer(i4), intent(in)             :: mneigh
  integer(i4), intent(in)             :: nneigh(*)
  integer(i4), intent(in)             :: neighno(maxneigh,*)
  real(dp),    intent(in)             :: xneigh(maxneigh,*)
  real(dp),    intent(in)             :: yneigh(maxneigh,*)
  real(dp),    intent(in)             :: zneigh(maxneigh,*)
  integer(i4), intent(in)             :: nfreeatom(*)
  integer(i4), intent(in)             :: nREBOatomRptr(*)
  integer(i4), intent(in)             :: nregion1
  real(dp),    intent(in)             :: d1(*)
  real(dp),    intent(in)             :: d2(*)
  logical,     intent(in)             :: ltorderv
!
!  Local variables
!
  integer(i4)                         :: ind
  integer(i4)                         :: ix,iy,iz
  integer(i4)                         :: ixf,iyf,izf
  integer(i4)                         :: jx,jy,jz
  integer(i4)                         :: jxf,jyf,jzf
  integer(i4)                         :: kx,ky,kz
  integer(i4)                         :: kxf,kyf,kzf
  integer(i4)                         :: lx,ly,lz
  integer(i4)                         :: lxf,lyf,lzf
  integer(i4)                         :: n1
  integer(i4)                         :: n1i
  integer(i4)                         :: n1j
  integer(i4)                         :: n1j2
  integer(i4)                         :: n2
  integer(i4)                         :: n2k
  integer(i4)                         :: n2l
  integer(i4)                         :: n2l2
  integer(i4)                         :: ndist
  integer(i4)                         :: nji
  integer(i4)                         :: nneigh2
  integer(i4)                         :: nri
  integer(i4)                         :: ntmp
  real(dp)                            :: r2ik
  real(dp)                            :: r2jl
  real(dp)                            :: xd1
  real(dp)                            :: yd1
  real(dp)                            :: zd1
  real(dp)                            :: xd2
  real(dp)                            :: yd2
  real(dp)                            :: zd2
  real(dp)                            :: xi
  real(dp)                            :: yi
  real(dp)                            :: zi
  real(dp)                            :: xik
  real(dp)                            :: yik
  real(dp)                            :: zik
  real(dp)                            :: xjk
  real(dp)                            :: yjk
  real(dp)                            :: zjk
  real(dp)                            :: xjl
  real(dp)                            :: yjl
  real(dp)                            :: zjl
  real(dp)                            :: x2i
  real(dp)                            :: y2i  
  real(dp)                            :: z2i
  logical                             :: lsameijkl
  logical                             :: lilocal
  logical                             :: ljlocal
  logical                             :: lklocal
  logical                             :: lllocal
  logical                             :: liok
  logical                             :: ljok
  logical                             :: lkok
  logical                             :: llok
  logical                             :: lopi
  logical                             :: lopj
  logical                             :: lopk
  logical                             :: lopl
#ifdef TRACE
  call trace_in('d2addddm')
#endif
!
!  Calculate total number of distances
!
  nri = nREBOatomRptr(i)
  nneigh2 = nneigh(nri) + nneigh(nri)*(nneigh(nri) + 1)/2
  if (ltorderv) then
    ndist = nneigh2 + nneigh(nri)*(mneigh + 1)*mneigh
  else
    ndist = nneigh2
  endif
!
!  Loop over pairs of interatomic distances
!
  ind = 0
  do n1 = 1,ndist
!******************************
!  Find pair of atoms for n1  *
!******************************
    if (n1.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
      n1i = i
      n1j = neighno(n1,nri)
      xd1 = xneigh(n1,nri)
      yd1 = yneigh(n1,nri)
      zd1 = zneigh(n1,nri)
      x2i = 0.0_dp
      y2i = 0.0_dp
      z2i = 0.0_dp
      liok = .true.
      ljok = .true.
    elseif (n1.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh(nri)
      n1j = 1
      n1j2 = 1
      do while (n1j2.lt.ntmp)
        n1j = n1j + 1
        n1j2 = n1j2 + n1j
      enddo
      n1j2 = n1j2 - n1j
      n1i = ntmp - n1j2
!
!  Check that neighbour numbers are valid
!
      liok = (n1i.le.nneigh(nri))
      ljok = (n1j.le.nneigh(nri))
      if (liok.and.ljok) then
!
!  Get vector between atoms
!
        xd1 = (xneigh(n1j,nri) - xneigh(n1i,nri))
        yd1 = (yneigh(n1j,nri) - yneigh(n1i,nri))
        zd1 = (zneigh(n1j,nri) - zneigh(n1i,nri))
!  
        x2i = xneigh(n1i,nri)
        y2i = yneigh(n1i,nri)
        z2i = zneigh(n1i,nri)
!
!  Convert neighbour numbers to real atoms
!
        n1i = neighno(n1i,nri)
        n1j = neighno(n1j,nri)
      endif
    else
!
!  Neighbours of neighbours of i -> neighbours of i
!
      ntmp = n1 - nneigh2
      nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
      ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
      n1i = (ntmp - 1)/mneigh
      n1j = ntmp - n1i*mneigh
!
!  Check that atoms are valid
!
      liok = (n1i.le.nneigh(nri))
      if (nji.le.nneigh(nri)) then
        ljok = (n1j.le.nneigh(neighno(nji,nri)))
      else
        ljok = .false.
      endif
!
      if (liok.and.ljok) then
        if (n1i.eq.0) then
          xi = 0.0_dp
          yi = 0.0_dp
          zi = 0.0_dp
          n1i = i
        else
          xi = xneigh(n1i,nri)
          yi = yneigh(n1i,nri)
          zi = zneigh(n1i,nri)
          n1i = neighno(n1i,nri)
        endif
!
        xd1 = (xneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
        yd1 = (yneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
        zd1 = (zneigh(n1j,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!         
        x2i = xi
        y2i = yi
        z2i = zi
!
        n1j = neighno(n1j,nREBOatomRptr(neighno(nji,nri)))
      endif
    endif
!
!  Skip the following section as it is irrelevant if i and j are not valid
!
    if (liok.and.ljok) then
      lopi = (nfreeatom(n1i).gt.0)
      lopj = (nfreeatom(n1j).gt.0)
    endif
!
!  If neither i nor j are being optimised then atoms are not valid
!
    if (.not.lopi.and..not.lopj) then
      liok = .false.
      ljok = .false.
    endif
!
!  If i and j are not valid then there is no point continuing for this value of n1
!
    if (liok.and.ljok) then
!
!  Set second derivative location pointers
!
      if (n1i.le.nregion1) then
        lilocal = (reg12local(n1i).gt.0)
        if (lilocal) then
          ix = 3*(reg12local(nfreeatom(n1i)) - 1) + 1
          iy = ix + 1
          iz = iy + 1
        endif
        ixf = 3*(nfreeatom(n1i) - 1) + 1
        iyf = ixf + 1
        izf = iyf + 1
      else
        lilocal = .false.
        ixf = 3*nregion1 + 1
        iyf = ixf + 1
        izf = iyf + 1
      endif
      if (n1j.le.nregion1) then
        ljlocal = (reg12local(n1j).gt.0)
        if (ljlocal) then
          jx = 3*(reg12local(nfreeatom(n1j)) - 1) + 1
          jy = jx + 1
          jz = jy + 1
        endif
        jxf = 3*(nfreeatom(n1j) - 1) + 1
        jyf = jxf + 1
        jzf = jyf + 1
      else
        ljlocal = .false.
        jxf = 3*nregion1 + 1
        jyf = jxf + 1
        jzf = jyf + 1
      endif
      do n2 = 1,n1
        ind = ind + 1
        if (abs(d2(ind)).gt.d2tol.or.(n1.eq.n2.and.abs(d1(n2)).gt.d1tol)) then
!******************************
!  Find pair of atoms for n2  *
!******************************
          if (n2.le.nneigh(nri)) then
!
!  i -> neighbours of i
!
            n2k = i
            n2l = neighno(n2,nri)
            xd2 = xneigh(n2,nri)
            yd2 = yneigh(n2,nri)
            zd2 = zneigh(n2,nri)
!
            xik = - x2i
            yik = - y2i
            zik = - z2i
!
            lkok = .true.
            llok = .true.
          elseif (n2.le.nneigh2) then
!
!  Neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh(nri)
            n2l = 1
            n2l2 = 1
            do while (n2l2.lt.ntmp)
              n2l = n2l + 1
              n2l2 = n2l2 + n2l
            enddo
            n2l2 = n2l2 - n2l
            n2k = ntmp - n2l2
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            llok = (n2l.le.nneigh(nri))
            if (lkok.and.llok) then
!
!  Calculate vector between atoms
!
              xd2 = (xneigh(n2l,nri) - xneigh(n2k,nri))
              yd2 = (yneigh(n2l,nri) - yneigh(n2k,nri))
              zd2 = (zneigh(n2l,nri) - zneigh(n2k,nri))
!  
              xik = xneigh(n2k,nri) - x2i
              yik = yneigh(n2k,nri) - y2i
              zik = zneigh(n2k,nri) - z2i
!
!  Convert neighbour numbers to real atoms
!
              n2k = neighno(n2k,nri)
              n2l = neighno(n2l,nri)
            endif
          else
!
!  Neighbours of neighbours of i -> neighbours of i
!
            ntmp = n2 - nneigh2
            nji = (ntmp - 1)/((mneigh + 1)*mneigh) + 1
            ntmp = ntmp - (nji - 1)*(mneigh + 1)*mneigh
            n2k = (ntmp - 1)/mneigh
            n2l = ntmp - n2k*mneigh
!
!  Check that atoms are valid
!
            lkok = (n2k.le.nneigh(nri))
            if (nji.le.nneigh(nri)) then
              llok = (n2l.le.nneigh(nREBOatomRptr(neighno(nji,nri))))
            else
              llok = .false.
            endif
!
            if (lkok.and.llok) then
              if (n2k.eq.0) then
                xi = 0.0_dp
                yi = 0.0_dp
                zi = 0.0_dp
                n2k = i
              else
                xi = xneigh(n2k,nri)
                yi = yneigh(n2k,nri)
                zi = zneigh(n2k,nri)
                n2k = neighno(n2k,nri)
              endif
!         
              xd2 = (xneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + xneigh(nji,nri) - xi)
              yd2 = (yneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + yneigh(nji,nri) - yi)
              zd2 = (zneigh(n2l,nREBOatomRptr(neighno(nji,nri))) + zneigh(nji,nri) - zi)
!               
              xik = xi - x2i
              yik = yi - y2i
              zik = zi - z2i
!       
              n2l = neighno(n2l,nREBOatomRptr(neighno(nji,nri)))
            endif
          endif
!
!  No point continuing beyond here unless k and l are valid
!
          if (lkok.and.llok) then
!
!  Complete remaining vectors
!
            xjk = xik - xd1
            yjk = yik - yd1
            zjk = zik - zd1
            xjl = xjk + xd2
            yjl = yjk + yd2
            zjl = zjk + zd2
!
!  Compute distances squared
!
            r2ik = xik*xik + yik*yik + zik*zik
            r2jl = xjl*xjl + yjl*yjl + zjl*zjl
!
            lopk = (nfreeatom(n2k).gt.0)
            lopl = (nfreeatom(n2l).gt.0)
!
            if (n2k.le.nregion1) then
              lklocal = (reg12local(n2k).gt.0)
              if (lklocal) then
                kx = 3*(reg12local(nfreeatom(n2k)) - 1) + 1
                ky = kx + 1
                kz = ky + 1
              endif
              kxf = 3*(nfreeatom(n2k) - 1) + 1
              kyf = kxf + 1
              kzf = kyf + 1
            else
              lklocal = .false.
              kxf = 3*nregion1 + 1
              kyf = kxf + 1
              kzf = kyf + 1
            endif
!
            if (n2l.le.nregion1) then
              lllocal = (reg12local(n2l).gt.0)
              if (lllocal) then
                lx = 3*(reg12local(nfreeatom(n2l)) - 1) + 1
                ly = lx + 1
                lz = ly + 1
              endif
              lxf = 3*(nfreeatom(n2l) - 1) + 1
              lyf = lxf + 1
              lzf = lyf + 1
            else
              lllocal = .false.
              lxf = 3*nregion1 + 1
              lyf = lxf + 1
              lzf = lyf + 1
            endif
!
            lsameijkl = .false.
            if (n1.eq.n2.and.n1i.eq.n2k.and.n1j.eq.n2l.and.r2ik.lt.1.0d-8.and.r2jl.lt.1.0d-8) lsameijkl = .true.
!
            if ((lopi.or.lopk).and.n1i.ne.n2k) then
!
!  Second derivatives : i - k
!
              if (lilocal) then
                derv2(kxf,ix) = derv2(kxf,ix) + xd1*xd2*d2(ind)
                derv2(kyf,ix) = derv2(kyf,ix) + xd1*yd2*d2(ind)
                derv2(kzf,ix) = derv2(kzf,ix) + xd1*zd2*d2(ind)
                derv2(kxf,iy) = derv2(kxf,iy) + yd1*xd2*d2(ind)
                derv2(kyf,iy) = derv2(kyf,iy) + yd1*yd2*d2(ind)
                derv2(kzf,iy) = derv2(kzf,iy) + yd1*zd2*d2(ind)
                derv2(kxf,iz) = derv2(kxf,iz) + zd1*xd2*d2(ind)
                derv2(kyf,iz) = derv2(kyf,iz) + zd1*yd2*d2(ind)
                derv2(kzf,iz) = derv2(kzf,iz) + zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(kxf,ix) = derv2(kxf,ix) + d1(n2)
                  derv2(kyf,iy) = derv2(kyf,iy) + d1(n2)
                  derv2(kzf,iz) = derv2(kzf,iz) + d1(n2)
                endif
              endif
              if (.not.lsameijkl) then
                if (lklocal) then
                  derv2(ixf,kx) = derv2(ixf,kx) + xd1*xd2*d2(ind)
                  derv2(ixf,ky) = derv2(ixf,ky) + xd1*yd2*d2(ind)
                  derv2(ixf,kz) = derv2(ixf,kz) + xd1*zd2*d2(ind)
                  derv2(iyf,kx) = derv2(iyf,kx) + yd1*xd2*d2(ind)
                  derv2(iyf,ky) = derv2(iyf,ky) + yd1*yd2*d2(ind)
                  derv2(iyf,kz) = derv2(iyf,kz) + yd1*zd2*d2(ind)
                  derv2(izf,kx) = derv2(izf,kx) + zd1*xd2*d2(ind)
                  derv2(izf,ky) = derv2(izf,ky) + zd1*yd2*d2(ind)
                  derv2(izf,kz) = derv2(izf,kz) + zd1*zd2*d2(ind)
                  if (n1.eq.n2) then
                    derv2(ixf,kx) = derv2(ixf,kx) + d1(n2)
                    derv2(iyf,ky) = derv2(iyf,ky) + d1(n2)
                    derv2(izf,kz) = derv2(izf,kz) + d1(n2)
                  endif
                endif
              endif
            endif
!
!  Second derivatives : i - l
!
            if ((lopi.or.lopl).and.n1i.ne.n2l) then
              if (lilocal) then
                derv2(lxf,ix) = derv2(lxf,ix) - xd1*xd2*d2(ind)
                derv2(lyf,ix) = derv2(lyf,ix) - xd1*yd2*d2(ind)
                derv2(lzf,ix) = derv2(lzf,ix) - xd1*zd2*d2(ind)
                derv2(lxf,iy) = derv2(lxf,iy) - yd1*xd2*d2(ind)
                derv2(lyf,iy) = derv2(lyf,iy) - yd1*yd2*d2(ind)
                derv2(lzf,iy) = derv2(lzf,iy) - yd1*zd2*d2(ind)
                derv2(lxf,iz) = derv2(lxf,iz) - zd1*xd2*d2(ind)
                derv2(lyf,iz) = derv2(lyf,iz) - zd1*yd2*d2(ind)
                derv2(lzf,iz) = derv2(lzf,iz) - zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(lxf,ix) = derv2(lxf,ix) - d1(n2)
                  derv2(lyf,iy) = derv2(lyf,iy) - d1(n2)
                  derv2(lzf,iz) = derv2(lzf,iz) - d1(n2)
                endif
              endif
              if (.not.lsameijkl) then
                if (lllocal) then
                  derv2(ixf,lx) = derv2(ixf,lx) - xd1*xd2*d2(ind)
                  derv2(ixf,ly) = derv2(ixf,ly) - xd1*yd2*d2(ind)
                  derv2(ixf,lz) = derv2(ixf,lz) - xd1*zd2*d2(ind)
                  derv2(iyf,lx) = derv2(iyf,lx) - yd1*xd2*d2(ind)
                  derv2(iyf,ly) = derv2(iyf,ly) - yd1*yd2*d2(ind)
                  derv2(iyf,lz) = derv2(iyf,lz) - yd1*zd2*d2(ind)
                  derv2(izf,lx) = derv2(izf,lx) - zd1*xd2*d2(ind)
                  derv2(izf,ly) = derv2(izf,ly) - zd1*yd2*d2(ind)
                  derv2(izf,lz) = derv2(izf,lz) - zd1*zd2*d2(ind)
                  if (n1.eq.n2) then
                    derv2(ixf,lx) = derv2(ixf,lx) - d1(n2)
                    derv2(iyf,ly) = derv2(iyf,ly) - d1(n2)
                    derv2(izf,lz) = derv2(izf,lz) - d1(n2)
                  endif
                endif
              endif
            endif
!
!  Second derivatives : j - k
!
            if ((lopj.or.lopk).and.n1j.ne.n2k) then
              if (ljlocal) then
                derv2(kxf,jx) = derv2(kxf,jx) - xd1*xd2*d2(ind)
                derv2(kyf,jx) = derv2(kyf,jx) - xd1*yd2*d2(ind)
                derv2(kzf,jx) = derv2(kzf,jx) - xd1*zd2*d2(ind)
                derv2(kxf,jy) = derv2(kxf,jy) - yd1*xd2*d2(ind)
                derv2(kyf,jy) = derv2(kyf,jy) - yd1*yd2*d2(ind)
                derv2(kzf,jy) = derv2(kzf,jy) - yd1*zd2*d2(ind)
                derv2(kxf,jz) = derv2(kxf,jz) - zd1*xd2*d2(ind)
                derv2(kyf,jz) = derv2(kyf,jz) - zd1*yd2*d2(ind)
                derv2(kzf,jz) = derv2(kzf,jz) - zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(kxf,jx) = derv2(kxf,jx) - d1(n2)
                  derv2(kyf,jy) = derv2(kyf,jy) - d1(n2)
                  derv2(kzf,jz) = derv2(kzf,jz) - d1(n2)
                endif
              endif
              if (.not.lsameijkl) then
                if (lklocal) then
                  derv2(jxf,kx) = derv2(jxf,kx) - xd1*xd2*d2(ind)
                  derv2(jxf,ky) = derv2(jxf,ky) - xd1*yd2*d2(ind)
                  derv2(jxf,kz) = derv2(jxf,kz) - xd1*zd2*d2(ind)
                  derv2(jyf,kx) = derv2(jyf,kx) - yd1*xd2*d2(ind)
                  derv2(jyf,ky) = derv2(jyf,ky) - yd1*yd2*d2(ind)
                  derv2(jyf,kz) = derv2(jyf,kz) - yd1*zd2*d2(ind)
                  derv2(jzf,kx) = derv2(jzf,kx) - zd1*xd2*d2(ind)
                  derv2(jzf,ky) = derv2(jzf,ky) - zd1*yd2*d2(ind)
                  derv2(jzf,kz) = derv2(jzf,kz) - zd1*zd2*d2(ind)
                  if (n1.eq.n2) then
                    derv2(jxf,kx) = derv2(jxf,kx) - d1(n2)
                    derv2(jyf,ky) = derv2(jyf,ky) - d1(n2)
                    derv2(jzf,kz) = derv2(jzf,kz) - d1(n2)
                  endif
                endif
              endif
            endif
!
!  Second derivatives : j - l
!
            if ((lopj.or.lopl).and.n1j.ne.n2l) then
              if (ljlocal) then
                derv2(lxf,jx) = derv2(lxf,jx) + xd1*xd2*d2(ind)
                derv2(lyf,jx) = derv2(lyf,jx) + xd1*yd2*d2(ind)
                derv2(lzf,jx) = derv2(lzf,jx) + xd1*zd2*d2(ind)
                derv2(lxf,jy) = derv2(lxf,jy) + yd1*xd2*d2(ind)
                derv2(lyf,jy) = derv2(lyf,jy) + yd1*yd2*d2(ind)
                derv2(lzf,jy) = derv2(lzf,jy) + yd1*zd2*d2(ind)
                derv2(lxf,jz) = derv2(lxf,jz) + zd1*xd2*d2(ind)
                derv2(lyf,jz) = derv2(lyf,jz) + zd1*yd2*d2(ind)
                derv2(lzf,jz) = derv2(lzf,jz) + zd1*zd2*d2(ind)
                if (n1.eq.n2) then
                  derv2(lxf,jx) = derv2(lxf,jx) + d1(n2)
                  derv2(lyf,jy) = derv2(lyf,jy) + d1(n2)
                  derv2(lzf,jz) = derv2(lzf,jz) + d1(n2)
                endif
              endif
              if (.not.lsameijkl) then
                if (lllocal) then
                  derv2(jxf,lx) = derv2(jxf,lx) + xd1*xd2*d2(ind)
                  derv2(jxf,ly) = derv2(jxf,ly) + xd1*yd2*d2(ind)
                  derv2(jxf,lz) = derv2(jxf,lz) + xd1*zd2*d2(ind)
                  derv2(jyf,lx) = derv2(jyf,lx) + yd1*xd2*d2(ind)
                  derv2(jyf,ly) = derv2(jyf,ly) + yd1*yd2*d2(ind)
                  derv2(jyf,lz) = derv2(jyf,lz) + yd1*zd2*d2(ind)
                  derv2(jzf,lx) = derv2(jzf,lx) + zd1*xd2*d2(ind)
                  derv2(jzf,ly) = derv2(jzf,ly) + zd1*yd2*d2(ind)
                  derv2(jzf,lz) = derv2(jzf,lz) + zd1*zd2*d2(ind)
                  if (n1.eq.n2) then
                    derv2(jxf,lx) = derv2(jxf,lx) + d1(n2)
                    derv2(jyf,ly) = derv2(jyf,ly) + d1(n2)
                    derv2(jzf,lz) = derv2(jzf,lz) + d1(n2)
                  endif
                endif
              endif
            endif
          endif
        endif
      enddo
    else
!
!  Increment ind pointer by number of missed terms
!
      ind = ind + n1
    endif
  enddo
#ifdef TRACE
  call trace_out('d2addddm')
#endif
!
  return
  end
!
  subroutine add_drv2_1(ix,iy,iz,jx,jy,jz,xij,yij,zij,d1trm,d2trm,dr2ijds,d2r2ijdx2, &
                        d2r2ijds2,d2r2ijdsdx)
!
!  Adds the second derivatives for a single distance
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d1trm             = first derivative of energy w.r.t. rij divided by rij
!  d2trm             = derivative of d1trm w.r.t. rij divided by rij
!
!  On exit :
!
!  Second derivatives have been added to
!
!  11/20 Created
!  12/22 Second derivatives now added to upper and lower half of matrix
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
  use datatypes
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: dr2ijds(6)
  real(dp),    intent(in)                          :: d2r2ijdx2(6)
  real(dp),    intent(in)                          :: d2r2ijds2(6,6)
  real(dp),    intent(in)                          :: d2r2ijdsdx(6,3)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_1')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  derv2(jx,ix) = derv2(jx,ix) - d2trm*d2r2ijdx2(1)
  derv2(jy,ix) = derv2(jy,ix) - d2trm*d2r2ijdx2(6)
  derv2(jz,ix) = derv2(jz,ix) - d2trm*d2r2ijdx2(5)
  derv2(jx,iy) = derv2(jx,iy) - d2trm*d2r2ijdx2(6)
  derv2(jy,iy) = derv2(jy,iy) - d2trm*d2r2ijdx2(2)
  derv2(jz,iy) = derv2(jz,iy) - d2trm*d2r2ijdx2(4)
  derv2(jx,iz) = derv2(jx,iz) - d2trm*d2r2ijdx2(5)
  derv2(jy,iz) = derv2(jy,iz) - d2trm*d2r2ijdx2(4)
  derv2(jz,iz) = derv2(jz,iz) - d2trm*d2r2ijdx2(3)
  derv2(jx,ix) = derv2(jx,ix) - d1trm
  derv2(jy,iy) = derv2(jy,iy) - d1trm
  derv2(jz,iz) = derv2(jz,iz) - d1trm
!
  if (jx.ne.ix) then
    derv2(ix,jx) = derv2(ix,jx) - d2trm*d2r2ijdx2(1)
    derv2(iy,jx) = derv2(iy,jx) - d2trm*d2r2ijdx2(6)
    derv2(iz,jx) = derv2(iz,jx) - d2trm*d2r2ijdx2(5)
    derv2(ix,jy) = derv2(ix,jy) - d2trm*d2r2ijdx2(6)
    derv2(iy,jy) = derv2(iy,jy) - d2trm*d2r2ijdx2(2)
    derv2(iz,jy) = derv2(iz,jy) - d2trm*d2r2ijdx2(4)
    derv2(ix,jz) = derv2(ix,jz) - d2trm*d2r2ijdx2(5)
    derv2(iy,jz) = derv2(iy,jz) - d2trm*d2r2ijdx2(4)
    derv2(iz,jz) = derv2(iz,jz) - d2trm*d2r2ijdx2(3)
    derv2(ix,jx) = derv2(ix,jx) - d1trm
    derv2(iy,jy) = derv2(iy,jy) - d1trm
    derv2(iz,jz) = derv2(iz,jz) - d1trm
  endif
!
  if (lstr) then
!
!  Mixed derivatives
!
    if (ix.ne.jx) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ijdsdx(ks,1) - xij*d2trm*dr2ijds(ks)
        derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ijdsdx(ks,2) - yij*d2trm*dr2ijds(ks)
        derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ijdsdx(ks,3) - zij*d2trm*dr2ijds(ks)
        derv3(jx,kl) = derv3(jx,kl) + d1trm*d2r2ijdsdx(ks,1) + xij*d2trm*dr2ijds(ks)
        derv3(jy,kl) = derv3(jy,kl) + d1trm*d2r2ijdsdx(ks,2) + yij*d2trm*dr2ijds(ks)
        derv3(jz,kl) = derv3(jz,kl) + d1trm*d2r2ijdsdx(ks,3) + zij*d2trm*dr2ijds(ks)
      enddo
    endif
!
!  Strain-strain
!
    do kk = 1,nstrains
      ks = nstrptr(kk)
      do kl = 1,nstrains
        kt = nstrptr(kl)
        sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2ijds(kt)*dr2ijds(ks) + d1trm*d2r2ijds2(kt,ks)
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('add_drv2_1')
#endif
!
  return
  end
!
  subroutine add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          xij,yij,zij,d1trm,d2trm,dr2ijds,d2r2ijdx2,d2r2ijds2,d2r2ijdsdx, &
                          ladd2sderv2)
!
!  Adds the second derivatives for a single distance
!
!  On entry : 
!
!  ix(f)             = position of x coordinate for i in second derivative matrix (f for full; otherwise local)
!  iy(f)             = position of y coordinate for i in second derivative matrix (f for full; otherwise local)
!  iz(f)             = position of z coordinate for i in second derivative matrix (f for full; otherwise local)
!  jx(f)             = position of x coordinate for j in second derivative matrix (f for full; otherwise local)
!  jy(f)             = position of y coordinate for j in second derivative matrix (f for full; otherwise local)
!  jz(f)             = position of z coordinate for j in second derivative matrix (f for full; otherwise local)
!  lilocal           = if .true. then i is local to this processor
!  ljlocal           = if .true. then j is local to this processor
!  ladd2sderv2       = if .true. then add terms to sderv2
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d1trm             = first derivative of energy w.r.t. rij divided by rij
!  d2trm             = derivative of d1trm w.r.t. rij divided by rij
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/22 Created from add_drv2_1
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
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: ladd2sderv2
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: dr2ijds(6)
  real(dp),    intent(in)                          :: d2r2ijdx2(6)
  real(dp),    intent(in)                          :: d2r2ijds2(6,6)
  real(dp),    intent(in)                          :: d2r2ijdsdx(6,3)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_1dm')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) - d2trm*d2r2ijdx2(1)
    derv2(jyf,ix) = derv2(jyf,ix) - d2trm*d2r2ijdx2(6)
    derv2(jzf,ix) = derv2(jzf,ix) - d2trm*d2r2ijdx2(5)
    derv2(jxf,iy) = derv2(jxf,iy) - d2trm*d2r2ijdx2(6)
    derv2(jyf,iy) = derv2(jyf,iy) - d2trm*d2r2ijdx2(2)
    derv2(jzf,iy) = derv2(jzf,iy) - d2trm*d2r2ijdx2(4)
    derv2(jxf,iz) = derv2(jxf,iz) - d2trm*d2r2ijdx2(5)
    derv2(jyf,iz) = derv2(jyf,iz) - d2trm*d2r2ijdx2(4)
    derv2(jzf,iz) = derv2(jzf,iz) - d2trm*d2r2ijdx2(3)
    derv2(jxf,ix) = derv2(jxf,ix) - d1trm
    derv2(jyf,iy) = derv2(jyf,iy) - d1trm
    derv2(jzf,iz) = derv2(jzf,iz) - d1trm
  endif
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) - d2trm*d2r2ijdx2(1)
    derv2(iyf,jx) = derv2(iyf,jx) - d2trm*d2r2ijdx2(6)
    derv2(izf,jx) = derv2(izf,jx) - d2trm*d2r2ijdx2(5)
    derv2(ixf,jy) = derv2(ixf,jy) - d2trm*d2r2ijdx2(6)
    derv2(iyf,jy) = derv2(iyf,jy) - d2trm*d2r2ijdx2(2)
    derv2(izf,jy) = derv2(izf,jy) - d2trm*d2r2ijdx2(4)
    derv2(ixf,jz) = derv2(ixf,jz) - d2trm*d2r2ijdx2(5)
    derv2(iyf,jz) = derv2(iyf,jz) - d2trm*d2r2ijdx2(4)
    derv2(izf,jz) = derv2(izf,jz) - d2trm*d2r2ijdx2(3)
    derv2(ixf,jx) = derv2(ixf,jx) - d1trm
    derv2(iyf,jy) = derv2(iyf,jy) - d1trm
    derv2(izf,jz) = derv2(izf,jz) - d1trm
  endif
!
  if (lstr) then
!
!  Mixed derivatives
!
    if (lilocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ijdsdx(ks,1) - xij*d2trm*dr2ijds(ks)
        derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ijdsdx(ks,2) - yij*d2trm*dr2ijds(ks)
        derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ijdsdx(ks,3) - zij*d2trm*dr2ijds(ks)
      enddo
    endif
    if (ljlocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(jx,kl) = derv3(jx,kl) + d1trm*d2r2ijdsdx(ks,1) + xij*d2trm*dr2ijds(ks)
        derv3(jy,kl) = derv3(jy,kl) + d1trm*d2r2ijdsdx(ks,2) + yij*d2trm*dr2ijds(ks)
        derv3(jz,kl) = derv3(jz,kl) + d1trm*d2r2ijdsdx(ks,3) + zij*d2trm*dr2ijds(ks)
      enddo
    endif
!
!  Strain-strain
!
    if (ladd2sderv2) then
      do kk = 1,nstrains
        ks = nstrptr(kk)
        do kl = 1,nstrains
          kt = nstrptr(kl)
          sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2ijds(kt)*dr2ijds(ks) + d1trm*d2r2ijds2(kt,ks)
        enddo
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_1dm')
#endif
!
  return
  end
!
  subroutine add_drv2_2(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,xik,yik,zik,xjl,yjl,zjl,d2trm,dr2ikds,dr2jlds)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  xik               = x component of vector from i to k
!  yik               = y component of vector from i to k
!  zik               = z component of vector from i to k
!  xjl               = x component of vector from i to k
!  yjl               = y component of vector from i to k
!  zjl               = z component of vector from i to k
!  d2trm             = second derivative of energy w.r.t. rik and rjl with division by distances
!
!  On exit :
!
!  Second derivatives have been added to
!
!  11/20 Created
!  12/22 Second derivatives now added to upper and lower half of matrix
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
  use datatypes
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xik
  real(dp),    intent(in)                          :: yik
  real(dp),    intent(in)                          :: zik
  real(dp),    intent(in)                          :: xjl
  real(dp),    intent(in)                          :: yjl
  real(dp),    intent(in)                          :: zjl
  real(dp),    intent(in)                          :: dr2ikds(6)
  real(dp),    intent(in)                          :: dr2jlds(6)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_2')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  derv2(jx,ix) = derv2(jx,ix) + d2trm*xjl*xik
  derv2(jy,ix) = derv2(jy,ix) + d2trm*yjl*xik
  derv2(jz,ix) = derv2(jz,ix) + d2trm*zjl*xik
  derv2(jx,iy) = derv2(jx,iy) + d2trm*xjl*yik
  derv2(jy,iy) = derv2(jy,iy) + d2trm*yjl*yik
  derv2(jz,iy) = derv2(jz,iy) + d2trm*zjl*yik
  derv2(jx,iz) = derv2(jx,iz) + d2trm*xjl*zik
  derv2(jy,iz) = derv2(jy,iz) + d2trm*yjl*zik
  derv2(jz,iz) = derv2(jz,iz) + d2trm*zjl*zik
!
  if (jx.ne.ix) then
    derv2(ix,jx) = derv2(ix,jx) + d2trm*xjl*xik
    derv2(iy,jx) = derv2(iy,jx) + d2trm*xjl*yik
    derv2(iz,jx) = derv2(iz,jx) + d2trm*xjl*zik
    derv2(ix,jy) = derv2(ix,jy) + d2trm*yjl*xik
    derv2(iy,jy) = derv2(iy,jy) + d2trm*yjl*yik
    derv2(iz,jy) = derv2(iz,jy) + d2trm*yjl*zik
    derv2(ix,jz) = derv2(ix,jz) + d2trm*zjl*xik
    derv2(iy,jz) = derv2(iy,jz) + d2trm*zjl*yik
    derv2(iz,jz) = derv2(iz,jz) + d2trm*zjl*zik
  endif
!
!  I-L
!
  derv2(lx,ix) = derv2(lx,ix) - d2trm*xjl*xik
  derv2(ly,ix) = derv2(ly,ix) - d2trm*yjl*xik
  derv2(lz,ix) = derv2(lz,ix) - d2trm*zjl*xik
  derv2(lx,iy) = derv2(lx,iy) - d2trm*xjl*yik
  derv2(ly,iy) = derv2(ly,iy) - d2trm*yjl*yik
  derv2(lz,iy) = derv2(lz,iy) - d2trm*zjl*yik
  derv2(lx,iz) = derv2(lx,iz) - d2trm*xjl*zik
  derv2(ly,iz) = derv2(ly,iz) - d2trm*yjl*zik
  derv2(lz,iz) = derv2(lz,iz) - d2trm*zjl*zik
!
  if (lx.ne.ix) then
    derv2(ix,lx) = derv2(ix,lx) - d2trm*xjl*xik
    derv2(iy,lx) = derv2(iy,lx) - d2trm*xjl*yik
    derv2(iz,lx) = derv2(iz,lx) - d2trm*xjl*zik
    derv2(ix,ly) = derv2(ix,ly) - d2trm*yjl*xik
    derv2(iy,ly) = derv2(iy,ly) - d2trm*yjl*yik
    derv2(iz,ly) = derv2(iz,ly) - d2trm*yjl*zik
    derv2(ix,lz) = derv2(ix,lz) - d2trm*zjl*xik
    derv2(iy,lz) = derv2(iy,lz) - d2trm*zjl*yik
    derv2(iz,lz) = derv2(iz,lz) - d2trm*zjl*zik
  endif
!
!  J-K
!
  derv2(kx,jx) = derv2(kx,jx) - d2trm*xjl*xik
  derv2(ky,jx) = derv2(ky,jx) - d2trm*xjl*yik
  derv2(kz,jx) = derv2(kz,jx) - d2trm*xjl*zik
  derv2(kx,jy) = derv2(kx,jy) - d2trm*yjl*xik
  derv2(ky,jy) = derv2(ky,jy) - d2trm*yjl*yik
  derv2(kz,jy) = derv2(kz,jy) - d2trm*yjl*zik
  derv2(kx,jz) = derv2(kx,jz) - d2trm*zjl*xik
  derv2(ky,jz) = derv2(ky,jz) - d2trm*zjl*yik
  derv2(kz,jz) = derv2(kz,jz) - d2trm*zjl*zik
!
  if (kx.ne.jx) then
    derv2(jx,kx) = derv2(jx,kx) - d2trm*xjl*xik
    derv2(jy,kx) = derv2(jy,kx) - d2trm*yjl*xik
    derv2(jz,kx) = derv2(jz,kx) - d2trm*zjl*xik
    derv2(jx,ky) = derv2(jx,ky) - d2trm*xjl*yik
    derv2(jy,ky) = derv2(jy,ky) - d2trm*yjl*yik
    derv2(jz,ky) = derv2(jz,ky) - d2trm*zjl*yik
    derv2(jx,kz) = derv2(jx,kz) - d2trm*xjl*zik
    derv2(jy,kz) = derv2(jy,kz) - d2trm*yjl*zik
    derv2(jz,kz) = derv2(jz,kz) - d2trm*zjl*zik
  endif
!
!  K-L
!     
  derv2(lx,kx) = derv2(lx,kx) + d2trm*xik*xjl
  derv2(ly,kx) = derv2(ly,kx) + d2trm*xik*yjl
  derv2(lz,kx) = derv2(lz,kx) + d2trm*xik*zjl
  derv2(lx,ky) = derv2(lx,ky) + d2trm*yik*xjl
  derv2(ly,ky) = derv2(ly,ky) + d2trm*yik*yjl
  derv2(lz,ky) = derv2(lz,ky) + d2trm*yik*zjl
  derv2(lx,kz) = derv2(lx,kz) + d2trm*zik*xjl
  derv2(ly,kz) = derv2(ly,kz) + d2trm*zik*yjl
  derv2(lz,kz) = derv2(lz,kz) + d2trm*zik*zjl
!
  if (lx.ne.kx) then
    derv2(kx,lx) = derv2(kx,lx) + d2trm*xik*xjl
    derv2(ky,lx) = derv2(ky,lx) + d2trm*yik*xjl
    derv2(kz,lx) = derv2(kz,lx) + d2trm*zik*xjl
    derv2(kx,ly) = derv2(kx,ly) + d2trm*xik*yjl
    derv2(ky,ly) = derv2(ky,ly) + d2trm*yik*yjl
    derv2(kz,ly) = derv2(kz,ly) + d2trm*zik*yjl
    derv2(kx,lz) = derv2(kx,lz) + d2trm*xik*zjl
    derv2(ky,lz) = derv2(ky,lz) + d2trm*yik*zjl
    derv2(kz,lz) = derv2(kz,lz) + d2trm*zik*zjl
  endif
!
!  Strain terms
!
  if (lstr) then
    do kl = 1,nstrains
      ks = nstrptr(kl)
      derv3(ix,kl) = derv3(ix,kl) - xik*d2trm*dr2jlds(ks)
      derv3(iy,kl) = derv3(iy,kl) - yik*d2trm*dr2jlds(ks)
      derv3(iz,kl) = derv3(iz,kl) - zik*d2trm*dr2jlds(ks)
      derv3(kx,kl) = derv3(kx,kl) + xik*d2trm*dr2jlds(ks)
      derv3(ky,kl) = derv3(ky,kl) + yik*d2trm*dr2jlds(ks)
      derv3(kz,kl) = derv3(kz,kl) + zik*d2trm*dr2jlds(ks)
    enddo
!
    do kl = 1,nstrains
      ks = nstrptr(kl)
      derv3(jx,kl) = derv3(jx,kl) - xjl*d2trm*dr2ikds(ks)
      derv3(jy,kl) = derv3(jy,kl) - yjl*d2trm*dr2ikds(ks)
      derv3(jz,kl) = derv3(jz,kl) - zjl*d2trm*dr2ikds(ks)
      derv3(lx,kl) = derv3(lx,kl) + xjl*d2trm*dr2ikds(ks)
      derv3(ly,kl) = derv3(ly,kl) + yjl*d2trm*dr2ikds(ks)
      derv3(lz,kl) = derv3(lz,kl) + zjl*d2trm*dr2ikds(ks)
    enddo
!
    do kk = 1,nstrains
      ks = nstrptr(kk)
      do kl = 1,nstrains
        kt = nstrptr(kl)
        sderv2(kl,kk) = sderv2(kl,kk) + d2trm*(dr2ikds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2ikds(ks))
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('add_drv2_2')
#endif
!
  return
  end
!
  subroutine add_drv2_2dm(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xjl,yjl,zjl, &
                          d2trm,dr2ikds,dr2jlds,ladd2sderv2)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l
!  Distributed memory parallel version.
!
!  On entry : 
!
!  ix(f)             = position of x coordinate for i in second derivative matrix (with f is in full matrix; without is local)
!  iy(f)             = position of y coordinate for i in second derivative matrix (with f is in full matrix; without is local)
!  iz(f)             = position of z coordinate for i in second derivative matrix (with f is in full matrix; without is local)
!  jx(f)             = position of x coordinate for j in second derivative matrix (with f is in full matrix; without is local)
!  jy(f)             = position of y coordinate for j in second derivative matrix (with f is in full matrix; without is local)
!  jz(f)             = position of z coordinate for j in second derivative matrix (with f is in full matrix; without is local)
!  kx(f)             = position of x coordinate for k in second derivative matrix (with f is in full matrix; without is local)
!  ky(f)             = position of y coordinate for k in second derivative matrix (with f is in full matrix; without is local)
!  kz(f)             = position of z coordinate for k in second derivative matrix (with f is in full matrix; without is local)
!  lx(f)             = position of x coordinate for l in second derivative matrix (with f is in full matrix; without is local)
!  ly(f)             = position of y coordinate for l in second derivative matrix (with f is in full matrix; without is local)
!  lz(f)             = position of z coordinate for l in second derivative matrix (with f is in full matrix; without is local)
!  lilocal           = if .true. then i is local to this processor
!  ljlocal           = if .true. then j is local to this processor
!  lklocal           = if .true. then k is local to this processor
!  lllocal           = if .true. then l is local to this processor
!  ladd2sderv2       = if .true. then add terms to sderv2
!  xik               = x component of vector from i to k
!  yik               = y component of vector from i to k
!  zik               = z component of vector from i to k
!  xjl               = x component of vector from i to k
!  yjl               = y component of vector from i to k
!  zjl               = z component of vector from i to k
!  d2trm             = second derivative of energy w.r.t. rik and rjl with division by distances
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/22 Created from add_drv2_2
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
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  logical,     intent(in)                          :: ladd2sderv2
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xik
  real(dp),    intent(in)                          :: yik
  real(dp),    intent(in)                          :: zik
  real(dp),    intent(in)                          :: xjl
  real(dp),    intent(in)                          :: yjl
  real(dp),    intent(in)                          :: zjl
  real(dp),    intent(in)                          :: dr2ikds(6)
  real(dp),    intent(in)                          :: dr2jlds(6)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_2dm')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) + d2trm*xjl*xik
    derv2(jyf,ix) = derv2(jyf,ix) + d2trm*yjl*xik
    derv2(jzf,ix) = derv2(jzf,ix) + d2trm*zjl*xik
    derv2(jxf,iy) = derv2(jxf,iy) + d2trm*xjl*yik
    derv2(jyf,iy) = derv2(jyf,iy) + d2trm*yjl*yik
    derv2(jzf,iy) = derv2(jzf,iy) + d2trm*zjl*yik
    derv2(jxf,iz) = derv2(jxf,iz) + d2trm*xjl*zik
    derv2(jyf,iz) = derv2(jyf,iz) + d2trm*yjl*zik
    derv2(jzf,iz) = derv2(jzf,iz) + d2trm*zjl*zik
  endif
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) + d2trm*xjl*xik
    derv2(iyf,jx) = derv2(iyf,jx) + d2trm*xjl*yik
    derv2(izf,jx) = derv2(izf,jx) + d2trm*xjl*zik
    derv2(ixf,jy) = derv2(ixf,jy) + d2trm*yjl*xik
    derv2(iyf,jy) = derv2(iyf,jy) + d2trm*yjl*yik
    derv2(izf,jy) = derv2(izf,jy) + d2trm*yjl*zik
    derv2(ixf,jz) = derv2(ixf,jz) + d2trm*zjl*xik
    derv2(iyf,jz) = derv2(iyf,jz) + d2trm*zjl*yik
    derv2(izf,jz) = derv2(izf,jz) + d2trm*zjl*zik
  endif
!
!  I-L
!
  if (lilocal) then
    derv2(lxf,ix) = derv2(lxf,ix) - d2trm*xjl*xik
    derv2(lyf,ix) = derv2(lyf,ix) - d2trm*yjl*xik
    derv2(lzf,ix) = derv2(lzf,ix) - d2trm*zjl*xik
    derv2(lxf,iy) = derv2(lxf,iy) - d2trm*xjl*yik
    derv2(lyf,iy) = derv2(lyf,iy) - d2trm*yjl*yik
    derv2(lzf,iy) = derv2(lzf,iy) - d2trm*zjl*yik
    derv2(lxf,iz) = derv2(lxf,iz) - d2trm*xjl*zik
    derv2(lyf,iz) = derv2(lyf,iz) - d2trm*yjl*zik
    derv2(lzf,iz) = derv2(lzf,iz) - d2trm*zjl*zik
  endif
  if (lllocal) then
    derv2(ixf,lx) = derv2(ixf,lx) - d2trm*xjl*xik
    derv2(iyf,lx) = derv2(iyf,lx) - d2trm*xjl*yik
    derv2(izf,lx) = derv2(izf,lx) - d2trm*xjl*zik
    derv2(ixf,ly) = derv2(ixf,ly) - d2trm*yjl*xik
    derv2(iyf,ly) = derv2(iyf,ly) - d2trm*yjl*yik
    derv2(izf,ly) = derv2(izf,ly) - d2trm*yjl*zik
    derv2(ixf,lz) = derv2(ixf,lz) - d2trm*zjl*xik
    derv2(iyf,lz) = derv2(iyf,lz) - d2trm*zjl*yik
    derv2(izf,lz) = derv2(izf,lz) - d2trm*zjl*zik
  endif
!
!  J-K
!
  if (ljlocal) then
    derv2(kxf,jx) = derv2(kxf,jx) - d2trm*xjl*xik
    derv2(kyf,jx) = derv2(kyf,jx) - d2trm*xjl*yik
    derv2(kzf,jx) = derv2(kzf,jx) - d2trm*xjl*zik
    derv2(kxf,jy) = derv2(kxf,jy) - d2trm*yjl*xik
    derv2(kyf,jy) = derv2(kyf,jy) - d2trm*yjl*yik
    derv2(kzf,jy) = derv2(kzf,jy) - d2trm*yjl*zik
    derv2(kxf,jz) = derv2(kxf,jz) - d2trm*zjl*xik
    derv2(kyf,jz) = derv2(kyf,jz) - d2trm*zjl*yik
    derv2(kzf,jz) = derv2(kzf,jz) - d2trm*zjl*zik
  endif
  if (lklocal) then
    derv2(jxf,kx) = derv2(jxf,kx) - d2trm*xjl*xik
    derv2(jyf,kx) = derv2(jyf,kx) - d2trm*yjl*xik
    derv2(jzf,kx) = derv2(jzf,kx) - d2trm*zjl*xik
    derv2(jxf,ky) = derv2(jxf,ky) - d2trm*xjl*yik
    derv2(jyf,ky) = derv2(jyf,ky) - d2trm*yjl*yik
    derv2(jzf,ky) = derv2(jzf,ky) - d2trm*zjl*yik
    derv2(jxf,kz) = derv2(jxf,kz) - d2trm*xjl*zik
    derv2(jyf,kz) = derv2(jyf,kz) - d2trm*yjl*zik
    derv2(jzf,kz) = derv2(jzf,kz) - d2trm*zjl*zik
  endif
!
!  K-L
!     
  if (lklocal) then
    derv2(lxf,kx) = derv2(lxf,kx) + d2trm*xik*xjl
    derv2(lyf,kx) = derv2(lyf,kx) + d2trm*xik*yjl
    derv2(lzf,kx) = derv2(lzf,kx) + d2trm*xik*zjl
    derv2(lxf,ky) = derv2(lxf,ky) + d2trm*yik*xjl
    derv2(lyf,ky) = derv2(lyf,ky) + d2trm*yik*yjl
    derv2(lzf,ky) = derv2(lzf,ky) + d2trm*yik*zjl
    derv2(lxf,kz) = derv2(lxf,kz) + d2trm*zik*xjl
    derv2(lyf,kz) = derv2(lyf,kz) + d2trm*zik*yjl
    derv2(lzf,kz) = derv2(lzf,kz) + d2trm*zik*zjl
  endif
  if (lllocal) then
    derv2(kxf,lx) = derv2(kxf,lx) + d2trm*xik*xjl
    derv2(kyf,lx) = derv2(kyf,lx) + d2trm*yik*xjl
    derv2(kzf,lx) = derv2(kzf,lx) + d2trm*zik*xjl
    derv2(kxf,ly) = derv2(kxf,ly) + d2trm*xik*yjl
    derv2(kyf,ly) = derv2(kyf,ly) + d2trm*yik*yjl
    derv2(kzf,ly) = derv2(kzf,ly) + d2trm*zik*yjl
    derv2(kxf,lz) = derv2(kxf,lz) + d2trm*xik*zjl
    derv2(kyf,lz) = derv2(kyf,lz) + d2trm*yik*zjl
    derv2(kzf,lz) = derv2(kzf,lz) + d2trm*zik*zjl
  endif
!
!  Strain terms
!
  if (lstr) then
    if (lilocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(ix,kl) = derv3(ix,kl) - xik*d2trm*dr2jlds(ks)
        derv3(iy,kl) = derv3(iy,kl) - yik*d2trm*dr2jlds(ks)
        derv3(iz,kl) = derv3(iz,kl) - zik*d2trm*dr2jlds(ks)
      enddo
    endif
    if (lklocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(kx,kl) = derv3(kx,kl) + xik*d2trm*dr2jlds(ks)
        derv3(ky,kl) = derv3(ky,kl) + yik*d2trm*dr2jlds(ks)
        derv3(kz,kl) = derv3(kz,kl) + zik*d2trm*dr2jlds(ks)
      enddo
    endif
!
    if (ljlocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(jx,kl) = derv3(jx,kl) - xjl*d2trm*dr2ikds(ks)
        derv3(jy,kl) = derv3(jy,kl) - yjl*d2trm*dr2ikds(ks)
        derv3(jz,kl) = derv3(jz,kl) - zjl*d2trm*dr2ikds(ks)
      enddo
    endif
    if (lllocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(lx,kl) = derv3(lx,kl) + xjl*d2trm*dr2ikds(ks)
        derv3(ly,kl) = derv3(ly,kl) + yjl*d2trm*dr2ikds(ks)
        derv3(lz,kl) = derv3(lz,kl) + zjl*d2trm*dr2ikds(ks)
      enddo
    endif
!
    if (ladd2sderv2) then
      do kk = 1,nstrains
        ks = nstrptr(kk)
        do kl = 1,nstrains
          kt = nstrptr(kl)
          sderv2(kl,kk) = sderv2(kl,kk) + d2trm*(dr2ikds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2ikds(ks))
        enddo
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_2dm')
#endif
!
  return
  end
!
  subroutine add_drv2_1b(ix,iy,iz,jx,jy,jz,d1trm,d2trm,d1c,d1s,d2c,d2sc)
!
!  Adds the second derivatives for a single distance
!  Version that uses pre-computed block
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  d1trm             = multiplier for second derivatives
!  d2trm             = multiplier for first derivative products
!  d1c               = first derivatives - Cartesian
!  d1s               = first derivatives - strain
!  d2c               = second derivatives - Cartesian/Cartesian
!  d2sc              = second derivatives - strain/Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!  12/20 Created
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, December 2020
!
  use datatypes
  use current,        only : nstrains
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1c(3)
  real(dp),    intent(in)                          :: d1s(6)
  real(dp),    intent(in)                          :: d2c(3,3)
  real(dp),    intent(in)                          :: d2sc(6,3)
!
!  Local variables
!
  integer(i4)                                      :: ks
#ifdef TRACE
  call trace_in('add_drv2_1b')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) - d1trm*d2c(1,1) - d2trm*d1c(1)*d1c(1)
    derv2(jy,ix) = derv2(jy,ix) - d1trm*d2c(2,1) - d2trm*d1c(2)*d1c(1)
    derv2(jz,ix) = derv2(jz,ix) - d1trm*d2c(3,1) - d2trm*d1c(3)*d1c(1)
    derv2(jx,iy) = derv2(jx,iy) - d1trm*d2c(1,2) - d2trm*d1c(1)*d1c(2)
    derv2(jy,iy) = derv2(jy,iy) - d1trm*d2c(2,2) - d2trm*d1c(2)*d1c(2)
    derv2(jz,iy) = derv2(jz,iy) - d1trm*d2c(3,2) - d2trm*d1c(3)*d1c(2)
    derv2(jx,iz) = derv2(jx,iz) - d1trm*d2c(1,3) - d2trm*d1c(1)*d1c(3)
    derv2(jy,iz) = derv2(jy,iz) - d1trm*d2c(2,3) - d2trm*d1c(2)*d1c(3)
    derv2(jz,iz) = derv2(jz,iz) - d1trm*d2c(3,3) - d2trm*d1c(3)*d1c(3)
  else
    derv2(ix,jx) = derv2(ix,jx) - d1trm*d2c(1,1) - d2trm*d1c(1)*d1c(1)
    derv2(iy,jx) = derv2(iy,jx) - d1trm*d2c(2,1) - d2trm*d1c(2)*d1c(1)
    derv2(iz,jx) = derv2(iz,jx) - d1trm*d2c(3,1) - d2trm*d1c(3)*d1c(1)
    derv2(ix,jy) = derv2(ix,jy) - d1trm*d2c(1,2) - d2trm*d1c(1)*d1c(2)
    derv2(iy,jy) = derv2(iy,jy) - d1trm*d2c(2,2) - d2trm*d1c(2)*d1c(2)
    derv2(iz,jy) = derv2(iz,jy) - d1trm*d2c(3,2) - d2trm*d1c(3)*d1c(2)
    derv2(ix,jz) = derv2(ix,jz) - d1trm*d2c(1,3) - d2trm*d1c(1)*d1c(3)
    derv2(iy,jz) = derv2(iy,jz) - d1trm*d2c(2,3) - d2trm*d1c(2)*d1c(3)
    derv2(iz,jz) = derv2(iz,jz) - d1trm*d2c(3,3) - d2trm*d1c(3)*d1c(3)
  endif
!
  if (lstr) then
!
!  Mixed derivatives
!
    if (ix.ne.jx) then
      do ks = 1,nstrains
        derv3(ix,ks) = derv3(ix,ks) - d1trm*d2sc(ks,1) - d2trm*d1s(ks)*d1c(1)
        derv3(iy,ks) = derv3(iy,ks) - d1trm*d2sc(ks,2) - d2trm*d1s(ks)*d1c(2)
        derv3(iz,ks) = derv3(iz,ks) - d1trm*d2sc(ks,3) - d2trm*d1s(ks)*d1c(3)
        derv3(jx,ks) = derv3(jx,ks) + d1trm*d2sc(ks,1) + d2trm*d1s(ks)*d1c(1)
        derv3(jy,ks) = derv3(jy,ks) + d1trm*d2sc(ks,2) + d2trm*d1s(ks)*d1c(2)
        derv3(jz,ks) = derv3(jz,ks) + d1trm*d2sc(ks,3) + d2trm*d1s(ks)*d1c(3)
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_1b')
#endif
!
  return
  end
!
  subroutine add_drv2_2b(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,d1trm,d2trm,d1ikc,d1jlc,d2c)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l
!  Version that uses pre-computed block
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d1trm             = multiplier for second derivatives
!  d2trm             = multiplier for first derivative products
!  d1ikc             = first derivatives for i-k - Cartesian
!  d1jlc             = first derivatives for j-l - Cartesian
!  d2c               = second derivatives - Cartesian/Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!  12/20 Created
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, December 2020
!
  use datatypes
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1ikc(3)
  real(dp),    intent(in)                          :: d1jlc(3)
  real(dp),    intent(in)                          :: d2c(3,3)
#ifdef TRACE
  call trace_in('add_drv2_2b')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(jy,ix) = derv2(jy,ix) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(jz,ix) = derv2(jz,ix) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(jx,iy) = derv2(jx,iy) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(jy,iy) = derv2(jy,iy) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(jz,iy) = derv2(jz,iy) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(jx,iz) = derv2(jx,iz) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(jy,iz) = derv2(jy,iz) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(jz,iz) = derv2(jz,iz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(ix,jx) = derv2(ix,jx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(iy,jx) = derv2(iy,jx) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(iz,jx) = derv2(iz,jx) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ix,jy) = derv2(ix,jy) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(iy,jy) = derv2(iy,jy) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(iz,jy) = derv2(iz,jy) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(ix,jz) = derv2(ix,jz) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(iy,jz) = derv2(iy,jz) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(iz,jz) = derv2(iz,jz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  I-L
!
  if (lx.ge.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(ly,ix) = derv2(ly,ix) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(lz,ix) = derv2(lz,ix) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(lx,iy) = derv2(lx,iy) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(ly,iy) = derv2(ly,iy) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(lz,iy) = derv2(lz,iy) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(lx,iz) = derv2(lx,iz) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ly,iz) = derv2(ly,iz) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(lz,iz) = derv2(lz,iz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(ix,lx) = derv2(ix,lx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(iy,lx) = derv2(iy,lx) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(iz,lx) = derv2(iz,lx) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ix,ly) = derv2(ix,ly) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(iy,ly) = derv2(iy,ly) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(iz,ly) = derv2(iz,ly) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(ix,lz) = derv2(ix,lz) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(iy,lz) = derv2(iy,lz) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(iz,lz) = derv2(iz,lz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  J-K
!
  if (kx.ge.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(ky,jx) = derv2(ky,jx) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(kz,jx) = derv2(kz,jx) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(kx,jy) = derv2(kx,jy) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(ky,jy) = derv2(ky,jy) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(kz,jy) = derv2(kz,jy) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(kx,jz) = derv2(kx,jz) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(ky,jz) = derv2(ky,jz) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(kz,jz) = derv2(kz,jz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(jx,kx) = derv2(jx,kx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(jy,kx) = derv2(jy,kx) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(jz,kx) = derv2(jz,kx) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(jx,ky) = derv2(jx,ky) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(jy,ky) = derv2(jy,ky) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(jz,ky) = derv2(jz,ky) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(jx,kz) = derv2(jx,kz) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(jy,kz) = derv2(jy,kz) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(jz,kz) = derv2(jz,kz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  K-L
!     
  if (lx.ge.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(ly,kx) = derv2(ly,kx) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(lz,kx) = derv2(lz,kx) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(lx,ky) = derv2(lx,ky) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(ly,ky) = derv2(ly,ky) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(lz,ky) = derv2(lz,ky) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(lx,kz) = derv2(lx,kz) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ly,kz) = derv2(ly,kz) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(lz,kz) = derv2(lz,kz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(kx,lx) = derv2(kx,lx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(ky,lx) = derv2(ky,lx) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(kz,lx) = derv2(kz,lx) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(kx,ly) = derv2(kx,ly) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(ky,ly) = derv2(ky,ly) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(kz,ly) = derv2(kz,ly) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(kx,lz) = derv2(kx,lz) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(ky,lz) = derv2(ky,lz) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(kz,lz) = derv2(kz,lz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2b')
#endif
!
  return
  end
!
  subroutine add_drv2_1bx(ix,iy,iz,jx,jy,jz,d2trm,d1c1,d1c2,d1s1,d1s2)
!
!  Adds the second derivatives for a single distance based only on products
!  of first derivatives.
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  d2trm             = multiplier for first derivative products
!  d1c1              = first derivatives - Cartesian for first term
!  d1c2              = first derivatives - Cartesian for second term
!  d1s1              = first derivatives - strain for first term
!  d1s2              = first derivatives - strain for second term
!
!  On exit :
!
!  Second derivatives have been added to
!
!  12/20 Created
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, December 2020
!
  use datatypes
  use current,        only : nstrains
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1c1(3)
  real(dp),    intent(in)                          :: d1c2(3)
  real(dp),    intent(in)                          :: d1s1(6)
  real(dp),    intent(in)                          :: d1s2(6)
!
!  Local variables
!
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_1bx')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) - d2trm*d1c1(1)*d1c2(1) - d2trm*d1c2(1)*d1c1(1)
    derv2(jy,ix) = derv2(jy,ix) - d2trm*d1c1(2)*d1c2(1) - d2trm*d1c2(2)*d1c1(1)
    derv2(jz,ix) = derv2(jz,ix) - d2trm*d1c1(3)*d1c2(1) - d2trm*d1c2(3)*d1c1(1)
    derv2(jx,iy) = derv2(jx,iy) - d2trm*d1c1(1)*d1c2(2) - d2trm*d1c2(1)*d1c1(2)
    derv2(jy,iy) = derv2(jy,iy) - d2trm*d1c1(2)*d1c2(2) - d2trm*d1c2(2)*d1c1(2)
    derv2(jz,iy) = derv2(jz,iy) - d2trm*d1c1(3)*d1c2(2) - d2trm*d1c2(3)*d1c1(2)
    derv2(jx,iz) = derv2(jx,iz) - d2trm*d1c1(1)*d1c2(3) - d2trm*d1c2(1)*d1c1(3)
    derv2(jy,iz) = derv2(jy,iz) - d2trm*d1c1(2)*d1c2(3) - d2trm*d1c2(2)*d1c1(3)
    derv2(jz,iz) = derv2(jz,iz) - d2trm*d1c1(3)*d1c2(3) - d2trm*d1c2(3)*d1c1(3)
  else
    derv2(ix,jx) = derv2(ix,jx) - d2trm*d1c1(1)*d1c2(1) - d2trm*d1c2(1)*d1c1(1)
    derv2(iy,jx) = derv2(iy,jx) - d2trm*d1c1(2)*d1c2(1) - d2trm*d1c2(2)*d1c1(1)
    derv2(iz,jx) = derv2(iz,jx) - d2trm*d1c1(3)*d1c2(1) - d2trm*d1c2(3)*d1c1(1)
    derv2(ix,jy) = derv2(ix,jy) - d2trm*d1c1(1)*d1c2(2) - d2trm*d1c2(1)*d1c1(2)
    derv2(iy,jy) = derv2(iy,jy) - d2trm*d1c1(2)*d1c2(2) - d2trm*d1c2(2)*d1c1(2)
    derv2(iz,jy) = derv2(iz,jy) - d2trm*d1c1(3)*d1c2(2) - d2trm*d1c2(3)*d1c1(2)
    derv2(ix,jz) = derv2(ix,jz) - d2trm*d1c1(1)*d1c2(3) - d2trm*d1c2(1)*d1c1(3)
    derv2(iy,jz) = derv2(iy,jz) - d2trm*d1c1(2)*d1c2(3) - d2trm*d1c2(2)*d1c1(3)
    derv2(iz,jz) = derv2(iz,jz) - d2trm*d1c1(3)*d1c2(3) - d2trm*d1c2(3)*d1c1(3)
  endif
!
  if (lstr) then
!
!  Mixed derivatives
!
    if (ix.ne.jx) then
      do ks = 1,nstrains
        derv3(ix,ks) = derv3(ix,ks) - d2trm*d1s1(ks)*d1c2(1) - d2trm*d1s2(ks)*d1c1(1)
        derv3(iy,ks) = derv3(iy,ks) - d2trm*d1s1(ks)*d1c2(2) - d2trm*d1s2(ks)*d1c1(2)
        derv3(iz,ks) = derv3(iz,ks) - d2trm*d1s1(ks)*d1c2(3) - d2trm*d1s2(ks)*d1c1(3)
        derv3(jx,ks) = derv3(jx,ks) + d2trm*d1s1(ks)*d1c2(1) + d2trm*d1s2(ks)*d1c1(1)
        derv3(jy,ks) = derv3(jy,ks) + d2trm*d1s1(ks)*d1c2(2) + d2trm*d1s2(ks)*d1c1(2)
        derv3(jz,ks) = derv3(jz,ks) + d2trm*d1s1(ks)*d1c2(3) + d2trm*d1s2(ks)*d1c1(3)
      enddo
    endif
!
!  Strain-strain derivatives
!
    do ks = 1,nstrains
      do kt = 1,nstrains
        sderv2(kt,ks) = sderv2(kt,ks) + d2trm*(d1s1(ks)*d1s2(kt) + d1s1(kt)*d1s2(ks))
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('add_drv2_1bx')
#endif
!
  return
  end
!
  subroutine add_drv2_2bx(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,d2trm,d1ikc,d1jlc,d1iks,d1jls)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l
!  Version that only adds the product of first derivatives
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d2trm             = multiplier for first derivative products
!  d1ikc             = first derivatives for i-k - Cartesian
!  d1jlc             = first derivatives for j-l - Cartesian
!  d1iks             = first derivatives for i-k - strain
!  d1jls             = first derivatives for j-l - strain
!
!  On exit :
!
!  Second derivatives have been added to
!
!  12/20 Created
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, December 2020
!
  use datatypes
  use current,        only : nstrains
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1ikc(3)
  real(dp),    intent(in)                          :: d1jlc(3)
  real(dp),    intent(in)                          :: d1iks(6)
  real(dp),    intent(in)                          :: d1jls(6)
!
!  Local variables
!
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_2bx')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(jy,ix) = derv2(jy,ix) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(jz,ix) = derv2(jz,ix) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(jx,iy) = derv2(jx,iy) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(jy,iy) = derv2(jy,iy) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(jz,iy) = derv2(jz,iy) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(jx,iz) = derv2(jx,iz) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(jy,iz) = derv2(jy,iz) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(jz,iz) = derv2(jz,iz) + d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(ix,jx) = derv2(ix,jx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(iy,jx) = derv2(iy,jx) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(iz,jx) = derv2(iz,jx) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ix,jy) = derv2(ix,jy) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(iy,jy) = derv2(iy,jy) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(iz,jy) = derv2(iz,jy) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(ix,jz) = derv2(ix,jz) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(iy,jz) = derv2(iy,jz) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(iz,jz) = derv2(iz,jz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  I-L
!
  if (lx.ge.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(ly,ix) = derv2(ly,ix) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(lz,ix) = derv2(lz,ix) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(lx,iy) = derv2(lx,iy) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(ly,iy) = derv2(ly,iy) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(lz,iy) = derv2(lz,iy) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(lx,iz) = derv2(lx,iz) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ly,iz) = derv2(ly,iz) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(lz,iz) = derv2(lz,iz) - d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(ix,lx) = derv2(ix,lx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(iy,lx) = derv2(iy,lx) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(iz,lx) = derv2(iz,lx) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ix,ly) = derv2(ix,ly) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(iy,ly) = derv2(iy,ly) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(iz,ly) = derv2(iz,ly) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(ix,lz) = derv2(ix,lz) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(iy,lz) = derv2(iy,lz) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(iz,lz) = derv2(iz,lz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  J-K
!
  if (kx.ge.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(ky,jx) = derv2(ky,jx) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(kz,jx) = derv2(kz,jx) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(kx,jy) = derv2(kx,jy) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(ky,jy) = derv2(ky,jy) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(kz,jy) = derv2(kz,jy) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(kx,jz) = derv2(kx,jz) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(ky,jz) = derv2(ky,jz) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(kz,jz) = derv2(kz,jz) - d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(jx,kx) = derv2(jx,kx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(jy,kx) = derv2(jy,kx) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(jz,kx) = derv2(jz,kx) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(jx,ky) = derv2(jx,ky) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(jy,ky) = derv2(jy,ky) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(jz,ky) = derv2(jz,ky) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(jx,kz) = derv2(jx,kz) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(jy,kz) = derv2(jy,kz) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(jz,kz) = derv2(jz,kz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  K-L
!     
  if (lx.ge.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(ly,kx) = derv2(ly,kx) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(lz,kx) = derv2(lz,kx) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(lx,ky) = derv2(lx,ky) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(ly,ky) = derv2(ly,ky) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(lz,ky) = derv2(lz,ky) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(lx,kz) = derv2(lx,kz) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ly,kz) = derv2(ly,kz) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(lz,kz) = derv2(lz,kz) + d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(kx,lx) = derv2(kx,lx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(ky,lx) = derv2(ky,lx) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(kz,lx) = derv2(kz,lx) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(kx,ly) = derv2(kx,ly) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(ky,ly) = derv2(ky,ly) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(kz,ly) = derv2(kz,ly) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(kx,lz) = derv2(kx,lz) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(ky,lz) = derv2(ky,lz) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(kz,lz) = derv2(kz,lz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (lstr) then
!
!  Mixed derivatives
!
    if (ix.ne.kx) then
      do ks = 1,nstrains
        derv3(ix,ks) = derv3(ix,ks) - d2trm*d1jls(ks)*d1ikc(1)
        derv3(iy,ks) = derv3(iy,ks) - d2trm*d1jls(ks)*d1ikc(2)
        derv3(iz,ks) = derv3(iz,ks) - d2trm*d1jls(ks)*d1ikc(3)
        derv3(kx,ks) = derv3(kx,ks) + d2trm*d1jls(ks)*d1ikc(1)
        derv3(ky,ks) = derv3(ky,ks) + d2trm*d1jls(ks)*d1ikc(2)
        derv3(kz,ks) = derv3(kz,ks) + d2trm*d1jls(ks)*d1ikc(3)
      enddo
    endif
!
    if (jx.ne.lx) then
      do ks = 1,nstrains
        derv3(jx,ks) = derv3(jx,ks) - d2trm*d1iks(ks)*d1jlc(1)
        derv3(jy,ks) = derv3(jy,ks) - d2trm*d1iks(ks)*d1jlc(2)
        derv3(jz,ks) = derv3(jz,ks) - d2trm*d1iks(ks)*d1jlc(3)
        derv3(lx,ks) = derv3(lx,ks) + d2trm*d1iks(ks)*d1jlc(1)
        derv3(ly,ks) = derv3(ly,ks) + d2trm*d1iks(ks)*d1jlc(2)
        derv3(lz,ks) = derv3(lz,ks) + d2trm*d1iks(ks)*d1jlc(3)
      enddo
    endif
!
!  Strain-strain derivatives
!
    do ks = 1,nstrains
      do kt = 1,nstrains
        sderv2(kt,ks) = sderv2(kt,ks) + d2trm*(d1iks(ks)*d1jls(kt) + d1iks(kt)*d1jls(ks))
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('add_drv2_2bx')
#endif
!
  return
  end
!
  subroutine add_drv2_1bdm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                           d1trm,d2trm,d1c,d1s,d2c,d2sc)
!
!  Adds the second derivatives for a single distance
!  Version that uses pre-computed block
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!  Distributed memory parallel version
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  d1trm             = multiplier for second derivatives
!  d2trm             = multiplier for first derivative products
!  d1c               = first derivatives - Cartesian
!  d1s               = first derivatives - strain
!  d2c               = second derivatives - Cartesian/Cartesian
!  d2sc              = second derivatives - strain/Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/22 Created from add_drv2_1b
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
  use current,        only : nstrains
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1c(3)
  real(dp),    intent(in)                          :: d1s(6)
  real(dp),    intent(in)                          :: d2c(3,3)
  real(dp),    intent(in)                          :: d2sc(6,3)
!
!  Local variables
!
  integer(i4)                                      :: ks
#ifdef TRACE
  call trace_in('add_drv2_1bdm')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) - d1trm*d2c(1,1) - d2trm*d1c(1)*d1c(1)
    derv2(jyf,ix) = derv2(jyf,ix) - d1trm*d2c(2,1) - d2trm*d1c(2)*d1c(1)
    derv2(jzf,ix) = derv2(jzf,ix) - d1trm*d2c(3,1) - d2trm*d1c(3)*d1c(1)
    derv2(jxf,iy) = derv2(jxf,iy) - d1trm*d2c(1,2) - d2trm*d1c(1)*d1c(2)
    derv2(jyf,iy) = derv2(jyf,iy) - d1trm*d2c(2,2) - d2trm*d1c(2)*d1c(2)
    derv2(jzf,iy) = derv2(jzf,iy) - d1trm*d2c(3,2) - d2trm*d1c(3)*d1c(2)
    derv2(jxf,iz) = derv2(jxf,iz) - d1trm*d2c(1,3) - d2trm*d1c(1)*d1c(3)
    derv2(jyf,iz) = derv2(jyf,iz) - d1trm*d2c(2,3) - d2trm*d1c(2)*d1c(3)
    derv2(jzf,iz) = derv2(jzf,iz) - d1trm*d2c(3,3) - d2trm*d1c(3)*d1c(3)
  endif
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) - d1trm*d2c(1,1) - d2trm*d1c(1)*d1c(1)
    derv2(iyf,jx) = derv2(iyf,jx) - d1trm*d2c(2,1) - d2trm*d1c(2)*d1c(1)
    derv2(izf,jx) = derv2(izf,jx) - d1trm*d2c(3,1) - d2trm*d1c(3)*d1c(1)
    derv2(ixf,jy) = derv2(ixf,jy) - d1trm*d2c(1,2) - d2trm*d1c(1)*d1c(2)
    derv2(iyf,jy) = derv2(iyf,jy) - d1trm*d2c(2,2) - d2trm*d1c(2)*d1c(2)
    derv2(izf,jy) = derv2(izf,jy) - d1trm*d2c(3,2) - d2trm*d1c(3)*d1c(2)
    derv2(ixf,jz) = derv2(ixf,jz) - d1trm*d2c(1,3) - d2trm*d1c(1)*d1c(3)
    derv2(iyf,jz) = derv2(iyf,jz) - d1trm*d2c(2,3) - d2trm*d1c(2)*d1c(3)
    derv2(izf,jz) = derv2(izf,jz) - d1trm*d2c(3,3) - d2trm*d1c(3)*d1c(3)
  endif
!
  if (lstr) then
!
!  Mixed derivatives
!
    if (lilocal) then
      do ks = 1,nstrains
        derv3(ix,ks) = derv3(ix,ks) - d1trm*d2sc(ks,1) - d2trm*d1s(ks)*d1c(1)
        derv3(iy,ks) = derv3(iy,ks) - d1trm*d2sc(ks,2) - d2trm*d1s(ks)*d1c(2)
        derv3(iz,ks) = derv3(iz,ks) - d1trm*d2sc(ks,3) - d2trm*d1s(ks)*d1c(3)
      enddo
    endif
    if (ljlocal) then
      do ks = 1,nstrains
        derv3(jx,ks) = derv3(jx,ks) + d1trm*d2sc(ks,1) + d2trm*d1s(ks)*d1c(1)
        derv3(jy,ks) = derv3(jy,ks) + d1trm*d2sc(ks,2) + d2trm*d1s(ks)*d1c(2)
        derv3(jz,ks) = derv3(jz,ks) + d1trm*d2sc(ks,3) + d2trm*d1s(ks)*d1c(3)
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_1bdm')
#endif
!
  return
  end
!
  subroutine add_drv2_2bdm(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                           jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d1trm,d2trm,d1ikc,d1jlc,d2c)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l
!  Version that uses pre-computed block
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!  Distributed memory parallel version
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d1trm             = multiplier for second derivatives
!  d2trm             = multiplier for first derivative products
!  d1ikc             = first derivatives for i-k - Cartesian
!  d1jlc             = first derivatives for j-l - Cartesian
!  d2c               = second derivatives - Cartesian/Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/22 Created from add_drv2_2b
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1ikc(3)
  real(dp),    intent(in)                          :: d1jlc(3)
  real(dp),    intent(in)                          :: d2c(3,3)
#ifdef TRACE
  call trace_in('add_drv2_2bdm')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(jyf,ix) = derv2(jyf,ix) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(jzf,ix) = derv2(jzf,ix) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(jxf,iy) = derv2(jxf,iy) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(jyf,iy) = derv2(jyf,iy) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(jzf,iy) = derv2(jzf,iy) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(jxf,iz) = derv2(jxf,iz) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(jyf,iz) = derv2(jyf,iz) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(jzf,iz) = derv2(jzf,iz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(iyf,jx) = derv2(iyf,jx) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(izf,jx) = derv2(izf,jx) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ixf,jy) = derv2(ixf,jy) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(iyf,jy) = derv2(iyf,jy) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(izf,jy) = derv2(izf,jy) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(ixf,jz) = derv2(ixf,jz) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(iyf,jz) = derv2(iyf,jz) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(izf,jz) = derv2(izf,jz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  I-L
!
  if (lilocal) then
    derv2(lxf,ix) = derv2(lxf,ix) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(lyf,ix) = derv2(lyf,ix) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(lzf,ix) = derv2(lzf,ix) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(lxf,iy) = derv2(lxf,iy) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(lyf,iy) = derv2(lyf,iy) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(lzf,iy) = derv2(lzf,iy) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(lxf,iz) = derv2(lxf,iz) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(lyf,iz) = derv2(lyf,iz) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(lzf,iz) = derv2(lzf,iz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (lllocal) then
    derv2(ixf,lx) = derv2(ixf,lx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(iyf,lx) = derv2(iyf,lx) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(izf,lx) = derv2(izf,lx) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ixf,ly) = derv2(ixf,ly) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(iyf,ly) = derv2(iyf,ly) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(izf,ly) = derv2(izf,ly) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(ixf,lz) = derv2(ixf,lz) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(iyf,lz) = derv2(iyf,lz) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(izf,lz) = derv2(izf,lz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  J-K
!
  if (ljlocal) then
    derv2(kxf,jx) = derv2(kxf,jx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(kyf,jx) = derv2(kyf,jx) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(kzf,jx) = derv2(kzf,jx) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(kxf,jy) = derv2(kxf,jy) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(kyf,jy) = derv2(kyf,jy) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(kzf,jy) = derv2(kzf,jy) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(kxf,jz) = derv2(kxf,jz) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(kyf,jz) = derv2(kyf,jz) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(kzf,jz) = derv2(kzf,jz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (lklocal) then
    derv2(jxf,kx) = derv2(jxf,kx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(jyf,kx) = derv2(jyf,kx) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(jzf,kx) = derv2(jzf,kx) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(jxf,ky) = derv2(jxf,ky) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(jyf,ky) = derv2(jyf,ky) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(jzf,ky) = derv2(jzf,ky) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(jxf,kz) = derv2(jxf,kz) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(jyf,kz) = derv2(jyf,kz) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(jzf,kz) = derv2(jzf,kz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  K-L
!     
  if (lklocal) then
    derv2(lxf,kx) = derv2(lxf,kx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(lyf,kx) = derv2(lyf,kx) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(lzf,kx) = derv2(lzf,kx) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(lxf,ky) = derv2(lxf,ky) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(lyf,ky) = derv2(lyf,ky) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(lzf,ky) = derv2(lzf,ky) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(lxf,kz) = derv2(lxf,kz) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(lyf,kz) = derv2(lyf,kz) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(lzf,kz) = derv2(lzf,kz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (lllocal) then
    derv2(kxf,lx) = derv2(kxf,lx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(kyf,lx) = derv2(kyf,lx) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(kzf,lx) = derv2(kzf,lx) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(kxf,ly) = derv2(kxf,ly) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(kyf,ly) = derv2(kyf,ly) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(kzf,ly) = derv2(kzf,ly) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(kxf,lz) = derv2(kxf,lz) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(kyf,lz) = derv2(kyf,lz) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(kzf,lz) = derv2(kzf,lz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2bdm')
#endif
!
  return
  end
!
  subroutine add_drv2_1bxdm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            d2trm,d1c1,d1c2,d1s1,d1s2,ladd2sderv2)
!
!  Adds the second derivatives for a single distance based only on products
!  of first derivatives.
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!  Distributed memory parallel version
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  d2trm             = multiplier for first derivative products
!  d1c1              = first derivatives - Cartesian for first term
!  d1c2              = first derivatives - Cartesian for second term
!  d1s1              = first derivatives - strain for first term
!  d1s2              = first derivatives - strain for second term
!  ladd2sderv2       = if .true. then add to sderv2
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/22 Created from add_drv2_1bx
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
  use current,        only : nstrains
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: ladd2sderv2
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1c1(3)
  real(dp),    intent(in)                          :: d1c2(3)
  real(dp),    intent(in)                          :: d1s1(6)
  real(dp),    intent(in)                          :: d1s2(6)
!
!  Local variables
!
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_1bxdm')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) - d2trm*d1c1(1)*d1c2(1) - d2trm*d1c2(1)*d1c1(1)
    derv2(jyf,ix) = derv2(jyf,ix) - d2trm*d1c1(2)*d1c2(1) - d2trm*d1c2(2)*d1c1(1)
    derv2(jzf,ix) = derv2(jzf,ix) - d2trm*d1c1(3)*d1c2(1) - d2trm*d1c2(3)*d1c1(1)
    derv2(jxf,iy) = derv2(jxf,iy) - d2trm*d1c1(1)*d1c2(2) - d2trm*d1c2(1)*d1c1(2)
    derv2(jyf,iy) = derv2(jyf,iy) - d2trm*d1c1(2)*d1c2(2) - d2trm*d1c2(2)*d1c1(2)
    derv2(jzf,iy) = derv2(jzf,iy) - d2trm*d1c1(3)*d1c2(2) - d2trm*d1c2(3)*d1c1(2)
    derv2(jxf,iz) = derv2(jxf,iz) - d2trm*d1c1(1)*d1c2(3) - d2trm*d1c2(1)*d1c1(3)
    derv2(jyf,iz) = derv2(jyf,iz) - d2trm*d1c1(2)*d1c2(3) - d2trm*d1c2(2)*d1c1(3)
    derv2(jzf,iz) = derv2(jzf,iz) - d2trm*d1c1(3)*d1c2(3) - d2trm*d1c2(3)*d1c1(3)
  endif
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) - d2trm*d1c1(1)*d1c2(1) - d2trm*d1c2(1)*d1c1(1)
    derv2(iyf,jx) = derv2(iyf,jx) - d2trm*d1c1(2)*d1c2(1) - d2trm*d1c2(2)*d1c1(1)
    derv2(izf,jx) = derv2(izf,jx) - d2trm*d1c1(3)*d1c2(1) - d2trm*d1c2(3)*d1c1(1)
    derv2(ixf,jy) = derv2(ixf,jy) - d2trm*d1c1(1)*d1c2(2) - d2trm*d1c2(1)*d1c1(2)
    derv2(iyf,jy) = derv2(iyf,jy) - d2trm*d1c1(2)*d1c2(2) - d2trm*d1c2(2)*d1c1(2)
    derv2(izf,jy) = derv2(izf,jy) - d2trm*d1c1(3)*d1c2(2) - d2trm*d1c2(3)*d1c1(2)
    derv2(ixf,jz) = derv2(ixf,jz) - d2trm*d1c1(1)*d1c2(3) - d2trm*d1c2(1)*d1c1(3)
    derv2(iyf,jz) = derv2(iyf,jz) - d2trm*d1c1(2)*d1c2(3) - d2trm*d1c2(2)*d1c1(3)
    derv2(izf,jz) = derv2(izf,jz) - d2trm*d1c1(3)*d1c2(3) - d2trm*d1c2(3)*d1c1(3)
  endif
!
  if (lstr) then
!
!  Mixed derivatives
!
    if (lilocal) then
      do ks = 1,nstrains
        derv3(ix,ks) = derv3(ix,ks) - d2trm*d1s1(ks)*d1c2(1) - d2trm*d1s2(ks)*d1c1(1)
        derv3(iy,ks) = derv3(iy,ks) - d2trm*d1s1(ks)*d1c2(2) - d2trm*d1s2(ks)*d1c1(2)
        derv3(iz,ks) = derv3(iz,ks) - d2trm*d1s1(ks)*d1c2(3) - d2trm*d1s2(ks)*d1c1(3)
      enddo
    endif
    if (ljlocal) then
      do ks = 1,nstrains
        derv3(jx,ks) = derv3(jx,ks) + d2trm*d1s1(ks)*d1c2(1) + d2trm*d1s2(ks)*d1c1(1)
        derv3(jy,ks) = derv3(jy,ks) + d2trm*d1s1(ks)*d1c2(2) + d2trm*d1s2(ks)*d1c1(2)
        derv3(jz,ks) = derv3(jz,ks) + d2trm*d1s1(ks)*d1c2(3) + d2trm*d1s2(ks)*d1c1(3)
      enddo
    endif
    if (ladd2sderv2) then
!
!  Strain-strain derivatives
!
      do ks = 1,nstrains
        do kt = 1,nstrains
          sderv2(kt,ks) = sderv2(kt,ks) + d2trm*(d1s1(ks)*d1s2(kt) + d1s1(kt)*d1s2(ks))
        enddo
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_1bxdm')
#endif
!
  return
  end
!
  subroutine add_drv2_2bxdm(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d2trm,d1ikc,d1jlc,d1iks,d1jls, &
                            ladd2sderv2)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l
!  Version that only adds the product of first derivatives
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!  Distributed memory parallel version
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d2trm             = multiplier for first derivative products
!  d1ikc             = first derivatives for i-k - Cartesian
!  d1jlc             = first derivatives for j-l - Cartesian
!  d1iks             = first derivatives for i-k - strain
!  d1jls             = first derivatives for j-l - strain
!  ladd2sderv2       = if .true. then add to sderv2
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/22 Created from add_drv2_2bx
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
  use current,        only : nstrains
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  logical,     intent(in)                          :: ladd2sderv2
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1ikc(3)
  real(dp),    intent(in)                          :: d1jlc(3)
  real(dp),    intent(in)                          :: d1iks(6)
  real(dp),    intent(in)                          :: d1jls(6)
!
!  Local variables
!
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_2bxdm')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(jyf,ix) = derv2(jyf,ix) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(jzf,ix) = derv2(jzf,ix) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(jxf,iy) = derv2(jxf,iy) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(jyf,iy) = derv2(jyf,iy) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(jzf,iy) = derv2(jzf,iy) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(jxf,iz) = derv2(jxf,iz) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(jyf,iz) = derv2(jyf,iz) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(jzf,iz) = derv2(jzf,iz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(iyf,jx) = derv2(iyf,jx) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(izf,jx) = derv2(izf,jx) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ixf,jy) = derv2(ixf,jy) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(iyf,jy) = derv2(iyf,jy) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(izf,jy) = derv2(izf,jy) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(ixf,jz) = derv2(ixf,jz) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(iyf,jz) = derv2(iyf,jz) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(izf,jz) = derv2(izf,jz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  I-L
!
  if (lilocal) then
    derv2(lxf,ix) = derv2(lxf,ix) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(lyf,ix) = derv2(lyf,ix) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(lzf,ix) = derv2(lzf,ix) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(lxf,iy) = derv2(lxf,iy) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(lyf,iy) = derv2(lyf,iy) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(lzf,iy) = derv2(lzf,iy) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(lxf,iz) = derv2(lxf,iz) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(lyf,iz) = derv2(lyf,iz) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(lzf,iz) = derv2(lzf,iz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (lllocal) then
    derv2(ixf,lx) = derv2(ixf,lx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(iyf,lx) = derv2(iyf,lx) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(izf,lx) = derv2(izf,lx) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ixf,ly) = derv2(ixf,ly) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(iyf,ly) = derv2(iyf,ly) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(izf,ly) = derv2(izf,ly) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(ixf,lz) = derv2(ixf,lz) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(iyf,lz) = derv2(iyf,lz) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(izf,lz) = derv2(izf,lz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  J-K
!
  if (ljlocal) then
    derv2(kxf,jx) = derv2(kxf,jx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(kyf,jx) = derv2(kyf,jx) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(kzf,jx) = derv2(kzf,jx) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(kxf,jy) = derv2(kxf,jy) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(kyf,jy) = derv2(kyf,jy) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(kzf,jy) = derv2(kzf,jy) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(kxf,jz) = derv2(kxf,jz) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(kyf,jz) = derv2(kyf,jz) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(kzf,jz) = derv2(kzf,jz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (lklocal) then
    derv2(jxf,kx) = derv2(jxf,kx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(jyf,kx) = derv2(jyf,kx) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(jzf,kx) = derv2(jzf,kx) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(jxf,ky) = derv2(jxf,ky) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(jyf,ky) = derv2(jyf,ky) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(jzf,ky) = derv2(jzf,ky) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(jxf,kz) = derv2(jxf,kz) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(jyf,kz) = derv2(jyf,kz) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(jzf,kz) = derv2(jzf,kz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  K-L
!     
  if (lklocal) then
    derv2(lxf,kx) = derv2(lxf,kx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(lyf,kx) = derv2(lyf,kx) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(lzf,kx) = derv2(lzf,kx) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(lxf,ky) = derv2(lxf,ky) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(lyf,ky) = derv2(lyf,ky) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(lzf,ky) = derv2(lzf,ky) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(lxf,kz) = derv2(lxf,kz) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(lyf,kz) = derv2(lyf,kz) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(lzf,kz) = derv2(lzf,kz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (lllocal) then
    derv2(kxf,lx) = derv2(kxf,lx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(kyf,lx) = derv2(kyf,lx) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(kzf,lx) = derv2(kzf,lx) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(kxf,ly) = derv2(kxf,ly) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(kyf,ly) = derv2(kyf,ly) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(kzf,ly) = derv2(kzf,ly) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(kxf,lz) = derv2(kxf,lz) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(kyf,lz) = derv2(kyf,lz) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(kzf,lz) = derv2(kzf,lz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
  if (lstr) then
!
!  Mixed derivatives
!
    if (lilocal) then
      do ks = 1,nstrains
        derv3(ix,ks) = derv3(ix,ks) - d2trm*d1jls(ks)*d1ikc(1)
        derv3(iy,ks) = derv3(iy,ks) - d2trm*d1jls(ks)*d1ikc(2)
        derv3(iz,ks) = derv3(iz,ks) - d2trm*d1jls(ks)*d1ikc(3)
      enddo
    endif
    if (lklocal) then
      do ks = 1,nstrains
        derv3(kx,ks) = derv3(kx,ks) + d2trm*d1jls(ks)*d1ikc(1)
        derv3(ky,ks) = derv3(ky,ks) + d2trm*d1jls(ks)*d1ikc(2)
        derv3(kz,ks) = derv3(kz,ks) + d2trm*d1jls(ks)*d1ikc(3)
      enddo
    endif
!
    if (ljlocal) then
      do ks = 1,nstrains
        derv3(jx,ks) = derv3(jx,ks) - d2trm*d1iks(ks)*d1jlc(1)
        derv3(jy,ks) = derv3(jy,ks) - d2trm*d1iks(ks)*d1jlc(2)
        derv3(jz,ks) = derv3(jz,ks) - d2trm*d1iks(ks)*d1jlc(3)
      enddo
    endif
    if (lllocal) then
      do ks = 1,nstrains
        derv3(lx,ks) = derv3(lx,ks) + d2trm*d1iks(ks)*d1jlc(1)
        derv3(ly,ks) = derv3(ly,ks) + d2trm*d1iks(ks)*d1jlc(2)
        derv3(lz,ks) = derv3(lz,ks) + d2trm*d1iks(ks)*d1jlc(3)
      enddo
    endif
    if (ladd2sderv2) then
!
!  Strain-strain derivatives
!
      do ks = 1,nstrains
        do kt = 1,nstrains
          sderv2(kt,ks) = sderv2(kt,ks) + d2trm*(d1iks(ks)*d1jls(kt) + d1iks(kt)*d1jls(ks))
        enddo
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_2bxdm')
#endif
!
  return
  end
!
  subroutine add_drv2_1g(ix,iy,iz,jx,jy,jz,xij,yij,zij,d1trm,d2trm,dr2ijds,d2r2ijdx2, &
                         d2r2ijds2,d2r2ijdsdx)
!
!  Adds the second derivatives for a single distance
!  Generalised version that replaces assumption of delta function for 1/2 (d2r/dadb)
!  and uses products of x/y/z
!  NB: Strain products must be handled elsewhere!
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d1trm             = first derivative of energy w.r.t. rij divided by rij
!  d2trm             = derivative of d1trm w.r.t. rij divided by rij
!
!  On exit :
!
!  Second derivatives have been added to
!
!   1/21 Created
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
!  Julian Gale, CIC, Curtin University, January 2021
!
  use datatypes
  use current,        only : nstrains, nstrptr
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: dr2ijds(6)
  real(dp),    intent(in)                          :: d2r2ijdx2(6)
  real(dp),    intent(in)                          :: d2r2ijds2(6,6)
  real(dp),    intent(in)                          :: d2r2ijdsdx(6,3)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_1g')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) - d1trm*d2r2ijdx2(1) - d2trm*xij*xij
    derv2(jy,ix) = derv2(jy,ix) - d1trm*d2r2ijdx2(6) - d2trm*yij*xij
    derv2(jz,ix) = derv2(jz,ix) - d1trm*d2r2ijdx2(5) - d2trm*zij*xij
    derv2(jx,iy) = derv2(jx,iy) - d1trm*d2r2ijdx2(6) - d2trm*xij*yij
    derv2(jy,iy) = derv2(jy,iy) - d1trm*d2r2ijdx2(2) - d2trm*yij*yij
    derv2(jz,iy) = derv2(jz,iy) - d1trm*d2r2ijdx2(4) - d2trm*zij*yij
    derv2(jx,iz) = derv2(jx,iz) - d1trm*d2r2ijdx2(5) - d2trm*xij*zij
    derv2(jy,iz) = derv2(jy,iz) - d1trm*d2r2ijdx2(4) - d2trm*yij*zij
    derv2(jz,iz) = derv2(jz,iz) - d1trm*d2r2ijdx2(3) - d2trm*zij*zij
  else
    derv2(ix,jx) = derv2(ix,jx) - d1trm*d2r2ijdx2(1) - d2trm*xij*xij
    derv2(iy,jx) = derv2(iy,jx) - d1trm*d2r2ijdx2(6) - d2trm*yij*xij
    derv2(iz,jx) = derv2(iz,jx) - d1trm*d2r2ijdx2(5) - d2trm*zij*xij
    derv2(ix,jy) = derv2(ix,jy) - d1trm*d2r2ijdx2(6) - d2trm*xij*yij
    derv2(iy,jy) = derv2(iy,jy) - d1trm*d2r2ijdx2(2) - d2trm*yij*yij
    derv2(iz,jy) = derv2(iz,jy) - d1trm*d2r2ijdx2(4) - d2trm*zij*yij
    derv2(ix,jz) = derv2(ix,jz) - d1trm*d2r2ijdx2(5) - d2trm*xij*zij
    derv2(iy,jz) = derv2(iy,jz) - d1trm*d2r2ijdx2(4) - d2trm*yij*zij
    derv2(iz,jz) = derv2(iz,jz) - d1trm*d2r2ijdx2(3) - d2trm*zij*zij
  endif
!
  if (lstr) then
!
!  Mixed derivatives
!
    if (ix.ne.jx) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ijdsdx(ks,1) - xij*d2trm*dr2ijds(ks)
        derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ijdsdx(ks,2) - yij*d2trm*dr2ijds(ks)
        derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ijdsdx(ks,3) - zij*d2trm*dr2ijds(ks)
        derv3(jx,kl) = derv3(jx,kl) + d1trm*d2r2ijdsdx(ks,1) + xij*d2trm*dr2ijds(ks)
        derv3(jy,kl) = derv3(jy,kl) + d1trm*d2r2ijdsdx(ks,2) + yij*d2trm*dr2ijds(ks)
        derv3(jz,kl) = derv3(jz,kl) + d1trm*d2r2ijdsdx(ks,3) + zij*d2trm*dr2ijds(ks)
      enddo
    endif
!
!  Strain-strain
!
    do kk = 1,nstrains
      ks = nstrptr(kk)
      do kl = 1,nstrains
        kt = nstrptr(kl)
        sderv2(kl,kk) = sderv2(kl,kk) + d1trm*d2r2ijds2(kt,ks)
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('add_drv2_1g')
#endif
!
  return
  end
!
  subroutine add_drv2_2sg(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,d1trm,dr2ik,dr2jl,dr2ikds,dr2jlds, &
                          d2trm,d2r2ikjl,d2r2ikjlds2,d2r2dsdik,d2r2dsdjl)
!
!  Adds the second derivatives for a single distance that depends on two components, i-k and j-l
!  Generalised version where matrices with the terms are passed in
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d1trm             = (1/r dE/dr)
!  d2trm             = (1/r d/dr (1/r dE/dr))
!
!  On exit :
!
!  Second derivatives have been added to
!
!   1/21 Created
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
!  Julian Gale, CIC, Curtin University, January 2021
!
  use datatypes
  use current,        only : nstrains, nstrptr
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: dr2ik(3)
  real(dp),    intent(in)                          :: d2r2ikjl(3,3)
  real(dp),    intent(in)                          :: d2r2ikjlds2(6,6)
  real(dp),    intent(in)                          :: d2r2dsdik(6,3)
  real(dp),    intent(in)                          :: d2r2dsdjl(6,3)
  real(dp),    intent(in)                          :: dr2jl(3)
  real(dp),    intent(in)                          :: dr2ikds(6)
  real(dp),    intent(in)                          :: dr2jlds(6)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_2sg')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(jy,ix) = derv2(jy,ix) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(jz,ix) = derv2(jz,ix) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(jx,iy) = derv2(jx,iy) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(jy,iy) = derv2(jy,iy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(jz,iy) = derv2(jz,iy) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(jx,iz) = derv2(jx,iz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(jy,iz) = derv2(jy,iz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(jz,iz) = derv2(jz,iz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(ix,jx) = derv2(ix,jx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(iy,jx) = derv2(iy,jx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(iz,jx) = derv2(iz,jx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ix,jy) = derv2(ix,jy) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(iy,jy) = derv2(iy,jy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(iz,jy) = derv2(iz,jy) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(ix,jz) = derv2(ix,jz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(iy,jz) = derv2(iy,jz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(iz,jz) = derv2(iz,jz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  I-L
!
  if (lx.ge.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(ly,ix) = derv2(ly,ix) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(lz,ix) = derv2(lz,ix) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(lx,iy) = derv2(lx,iy) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(ly,iy) = derv2(ly,iy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(lz,iy) = derv2(lz,iy) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(lx,iz) = derv2(lx,iz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ly,iz) = derv2(ly,iz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(lz,iz) = derv2(lz,iz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(ix,lx) = derv2(ix,lx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(iy,lx) = derv2(iy,lx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(iz,lx) = derv2(iz,lx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ix,ly) = derv2(ix,ly) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(iy,ly) = derv2(iy,ly) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(iz,ly) = derv2(iz,ly) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(ix,lz) = derv2(ix,lz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(iy,lz) = derv2(iy,lz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(iz,lz) = derv2(iz,lz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  J-K
!
  if (kx.ge.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(ky,jx) = derv2(ky,jx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(kz,jx) = derv2(kz,jx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(kx,jy) = derv2(kx,jy) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(ky,jy) = derv2(ky,jy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(kz,jy) = derv2(kz,jy) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(kx,jz) = derv2(kx,jz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(ky,jz) = derv2(ky,jz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(kz,jz) = derv2(kz,jz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(jx,kx) = derv2(jx,kx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(jy,kx) = derv2(jy,kx) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(jz,kx) = derv2(jz,kx) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(jx,ky) = derv2(jx,ky) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(jy,ky) = derv2(jy,ky) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(jz,ky) = derv2(jz,ky) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(jx,kz) = derv2(jx,kz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(jy,kz) = derv2(jy,kz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(jz,kz) = derv2(jz,kz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  K-L
!     
  if (lx.ge.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(ly,kx) = derv2(ly,kx) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(lz,kx) = derv2(lz,kx) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(lx,ky) = derv2(lx,ky) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(ly,ky) = derv2(ly,ky) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(lz,ky) = derv2(lz,ky) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(lx,kz) = derv2(lx,kz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ly,kz) = derv2(ly,kz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(lz,kz) = derv2(lz,kz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(kx,lx) = derv2(kx,lx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(ky,lx) = derv2(ky,lx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(kz,lx) = derv2(kz,lx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(kx,ly) = derv2(kx,ly) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(ky,ly) = derv2(ky,ly) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(kz,ly) = derv2(kz,ly) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(kx,lz) = derv2(kx,lz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(ky,lz) = derv2(ky,lz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(kz,lz) = derv2(kz,lz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  Strain terms
!
  if (lstr) then
    do kl = 1,nstrains
      ks = nstrptr(kl)
      derv3(ix,kl) = derv3(ix,kl) - d2trm*dr2jlds(ks)*dr2ik(1) - d1trm*d2r2dsdik(ks,1)
      derv3(iy,kl) = derv3(iy,kl) - d2trm*dr2jlds(ks)*dr2ik(2) - d1trm*d2r2dsdik(ks,2)
      derv3(iz,kl) = derv3(iz,kl) - d2trm*dr2jlds(ks)*dr2ik(3) - d1trm*d2r2dsdik(ks,3)
      derv3(kx,kl) = derv3(kx,kl) + d2trm*dr2jlds(ks)*dr2ik(1) + d1trm*d2r2dsdik(ks,1)
      derv3(ky,kl) = derv3(ky,kl) + d2trm*dr2jlds(ks)*dr2ik(2) + d1trm*d2r2dsdik(ks,2)
      derv3(kz,kl) = derv3(kz,kl) + d2trm*dr2jlds(ks)*dr2ik(3) + d1trm*d2r2dsdik(ks,3)
    enddo
!
    do kl = 1,nstrains
      ks = nstrptr(kl)
      derv3(jx,kl) = derv3(jx,kl) - d2trm*dr2ikds(ks)*dr2jl(1) - d1trm*d2r2dsdjl(ks,1)
      derv3(jy,kl) = derv3(jy,kl) - d2trm*dr2ikds(ks)*dr2jl(2) - d1trm*d2r2dsdjl(ks,2)
      derv3(jz,kl) = derv3(jz,kl) - d2trm*dr2ikds(ks)*dr2jl(3) - d1trm*d2r2dsdjl(ks,3)
      derv3(lx,kl) = derv3(lx,kl) + d2trm*dr2ikds(ks)*dr2jl(1) + d1trm*d2r2dsdjl(ks,1)
      derv3(ly,kl) = derv3(ly,kl) + d2trm*dr2ikds(ks)*dr2jl(2) + d1trm*d2r2dsdjl(ks,2)
      derv3(lz,kl) = derv3(lz,kl) + d2trm*dr2ikds(ks)*dr2jl(3) + d1trm*d2r2dsdjl(ks,3)
    enddo
!
    do kk = 1,nstrains
      ks = nstrptr(kk)
      do kl = 1,nstrains
        kt = nstrptr(kl)
        sderv2(kl,kk) = sderv2(kl,kk) + d2trm*(dr2ikds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2ikds(ks)) + &
                                        d1trm*d2r2ikjlds2(kt,ks)
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('add_drv2_2sg')
#endif
!
  return
  end
!
  subroutine add_drv2_2sgnos(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,d1trm,dr2ik,dr2jl,d2trm,d2r2ikjl)
!
!  Adds the second derivatives for a single distance that depends on two components, i-k and j-l
!  Generalised version where matrices with the terms are passed in
!  NB: This version has no strain contributions
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d1trm             = (1/r dE/dr)
!  d2trm             = (1/r d/dr (1/r dE/dr))
!
!  On exit :
!
!  Second derivatives have been added to
!
!   1/21 Created
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
!  Julian Gale, CIC, Curtin University, January 2021
!
  use datatypes
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: dr2ik(3)
  real(dp),    intent(in)                          :: d2r2ikjl(3,3)
  real(dp),    intent(in)                          :: dr2jl(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2sgnos')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(jy,ix) = derv2(jy,ix) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(jz,ix) = derv2(jz,ix) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(jx,iy) = derv2(jx,iy) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(jy,iy) = derv2(jy,iy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(jz,iy) = derv2(jz,iy) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(jx,iz) = derv2(jx,iz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(jy,iz) = derv2(jy,iz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(jz,iz) = derv2(jz,iz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(ix,jx) = derv2(ix,jx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(iy,jx) = derv2(iy,jx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(iz,jx) = derv2(iz,jx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ix,jy) = derv2(ix,jy) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(iy,jy) = derv2(iy,jy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(iz,jy) = derv2(iz,jy) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(ix,jz) = derv2(ix,jz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(iy,jz) = derv2(iy,jz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(iz,jz) = derv2(iz,jz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  I-L
!
  if (lx.ge.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(ly,ix) = derv2(ly,ix) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(lz,ix) = derv2(lz,ix) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(lx,iy) = derv2(lx,iy) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(ly,iy) = derv2(ly,iy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(lz,iy) = derv2(lz,iy) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(lx,iz) = derv2(lx,iz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ly,iz) = derv2(ly,iz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(lz,iz) = derv2(lz,iz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(ix,lx) = derv2(ix,lx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(iy,lx) = derv2(iy,lx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(iz,lx) = derv2(iz,lx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ix,ly) = derv2(ix,ly) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(iy,ly) = derv2(iy,ly) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(iz,ly) = derv2(iz,ly) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(ix,lz) = derv2(ix,lz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(iy,lz) = derv2(iy,lz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(iz,lz) = derv2(iz,lz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  J-K
!
  if (kx.ge.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(ky,jx) = derv2(ky,jx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(kz,jx) = derv2(kz,jx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(kx,jy) = derv2(kx,jy) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(ky,jy) = derv2(ky,jy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(kz,jy) = derv2(kz,jy) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(kx,jz) = derv2(kx,jz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(ky,jz) = derv2(ky,jz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(kz,jz) = derv2(kz,jz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(jx,kx) = derv2(jx,kx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(jy,kx) = derv2(jy,kx) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(jz,kx) = derv2(jz,kx) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(jx,ky) = derv2(jx,ky) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(jy,ky) = derv2(jy,ky) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(jz,ky) = derv2(jz,ky) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(jx,kz) = derv2(jx,kz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(jy,kz) = derv2(jy,kz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(jz,kz) = derv2(jz,kz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  K-L
!     
  if (lx.ge.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(ly,kx) = derv2(ly,kx) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(lz,kx) = derv2(lz,kx) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(lx,ky) = derv2(lx,ky) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(ly,ky) = derv2(ly,ky) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(lz,ky) = derv2(lz,ky) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(lx,kz) = derv2(lx,kz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ly,kz) = derv2(ly,kz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(lz,kz) = derv2(lz,kz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(kx,lx) = derv2(kx,lx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(ky,lx) = derv2(ky,lx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(kz,lx) = derv2(kz,lx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(kx,ly) = derv2(kx,ly) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(ky,ly) = derv2(ky,ly) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(kz,ly) = derv2(kz,ly) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(kx,lz) = derv2(kx,lz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(ky,lz) = derv2(ky,lz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(kz,lz) = derv2(kz,lz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2sgnos')
#endif
!
  return
  end
!
  subroutine add_drv2_2sgnosdm(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                               jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d1trm,dr2ik,dr2jl,d2trm,d2r2ikjl)
!
!  Adds the second derivatives for a single distance that depends on two components, i-k and j-l
!  Generalised version where matrices with the terms are passed in
!  NB: This version has no strain contributions
!  Distributed memory parallel version.
!
!  On entry : 
!
!  ix(f)             = position of x coordinate for i in second derivative matrix (with f is in full matrix; without is local)
!  iy(f)             = position of y coordinate for i in second derivative matrix (with f is in full matrix; without is local)
!  iz(f)             = position of z coordinate for i in second derivative matrix (with f is in full matrix; without is local)
!  jx(f)             = position of x coordinate for j in second derivative matrix (with f is in full matrix; without is local)
!  jy(f)             = position of y coordinate for j in second derivative matrix (with f is in full matrix; without is local)
!  jz(f)             = position of z coordinate for j in second derivative matrix (with f is in full matrix; without is local)
!  kx(f)             = position of x coordinate for k in second derivative matrix (with f is in full matrix; without is local)
!  ky(f)             = position of y coordinate for k in second derivative matrix (with f is in full matrix; without is local)
!  kz(f)             = position of z coordinate for k in second derivative matrix (with f is in full matrix; without is local)
!  lx(f)             = position of x coordinate for l in second derivative matrix (with f is in full matrix; without is local)
!  ly(f)             = position of y coordinate for l in second derivative matrix (with f is in full matrix; without is local)
!  lz(f)             = position of z coordinate for l in second derivative matrix (with f is in full matrix; without is local)
!  lilocal           = if .true. then i is local to this processor
!  ljlocal           = if .true. then j is local to this processor
!  lklocal           = if .true. then k is local to this processor
!  lllocal           = if .true. then l is local to this processor
!  d1trm             = (1/r dE/dr)
!  d2trm             = (1/r d/dr (1/r dE/dr))
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_2sgnos
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: dr2ik(3)
  real(dp),    intent(in)                          :: d2r2ikjl(3,3)
  real(dp),    intent(in)                          :: dr2jl(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2sgnosdm')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(jyf,ix) = derv2(jyf,ix) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(jzf,ix) = derv2(jzf,ix) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(jxf,iy) = derv2(jxf,iy) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(jyf,iy) = derv2(jyf,iy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(jzf,iy) = derv2(jzf,iy) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(jxf,iz) = derv2(jxf,iz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(jyf,iz) = derv2(jyf,iz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(jzf,iz) = derv2(jzf,iz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(iyf,jx) = derv2(iyf,jx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(izf,jx) = derv2(izf,jx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ixf,jy) = derv2(ixf,jy) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(iyf,jy) = derv2(iyf,jy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(izf,jy) = derv2(izf,jy) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(ixf,jz) = derv2(ixf,jz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(iyf,jz) = derv2(iyf,jz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(izf,jz) = derv2(izf,jz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  I-L
!
  if (lilocal) then
    derv2(lxf,ix) = derv2(lxf,ix) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(lyf,ix) = derv2(lyf,ix) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(lzf,ix) = derv2(lzf,ix) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(lxf,iy) = derv2(lxf,iy) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(lyf,iy) = derv2(lyf,iy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(lzf,iy) = derv2(lzf,iy) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(lxf,iz) = derv2(lxf,iz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(lyf,iz) = derv2(lyf,iz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(lzf,iz) = derv2(lzf,iz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lllocal) then
    derv2(ixf,lx) = derv2(ixf,lx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(iyf,lx) = derv2(iyf,lx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(izf,lx) = derv2(izf,lx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ixf,ly) = derv2(ixf,ly) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(iyf,ly) = derv2(iyf,ly) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(izf,ly) = derv2(izf,ly) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(ixf,lz) = derv2(ixf,lz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(iyf,lz) = derv2(iyf,lz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(izf,lz) = derv2(izf,lz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  J-K
!
  if (ljlocal) then
    derv2(kxf,jx) = derv2(kxf,jx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(kyf,jx) = derv2(kyf,jx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(kzf,jx) = derv2(kzf,jx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(kxf,jy) = derv2(kxf,jy) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(kyf,jy) = derv2(kyf,jy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(kzf,jy) = derv2(kzf,jy) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(kxf,jz) = derv2(kxf,jz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(kyf,jz) = derv2(kyf,jz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(kzf,jz) = derv2(kzf,jz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lklocal) then
    derv2(jxf,kx) = derv2(jxf,kx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(jyf,kx) = derv2(jyf,kx) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(jzf,kx) = derv2(jzf,kx) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(jxf,ky) = derv2(jxf,ky) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(jyf,ky) = derv2(jyf,ky) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(jzf,ky) = derv2(jzf,ky) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(jxf,kz) = derv2(jxf,kz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(jyf,kz) = derv2(jyf,kz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(jzf,kz) = derv2(jzf,kz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  K-L
!     
  if (lklocal) then
    derv2(lxf,kx) = derv2(lxf,kx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(lyf,kx) = derv2(lyf,kx) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(lzf,kx) = derv2(lzf,kx) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(lxf,ky) = derv2(lxf,ky) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(lyf,ky) = derv2(lyf,ky) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(lzf,ky) = derv2(lzf,ky) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(lxf,kz) = derv2(lxf,kz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(lyf,kz) = derv2(lyf,kz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(lzf,kz) = derv2(lzf,kz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lllocal) then
    derv2(kxf,lx) = derv2(kxf,lx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(kyf,lx) = derv2(kyf,lx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(kzf,lx) = derv2(kzf,lx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(kxf,ly) = derv2(kxf,ly) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(kyf,ly) = derv2(kyf,ly) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(kzf,ly) = derv2(kzf,ly) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(kxf,lz) = derv2(kxf,lz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(kyf,lz) = derv2(kyf,lz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(kzf,lz) = derv2(kzf,lz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2sgnosdm')
#endif
!
  return
  end
!
  subroutine add_drv2_1gdm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                           d1trm,d2trm,dr2ijds,d2r2ijdx2,d2r2ijds2,d2r2ijdsdx,ladd2sderv2)
!
!  Adds the second derivatives for a single distance
!  Generalised version that replaces assumption of delta function for 1/2 (d2r/dadb)
!  and uses products of x/y/z
!  NB: Strain products must be handled elsewhere!
!  Distributed memory parallel version.
!
!  On entry : 
!
!  ix(f)             = position of x coordinate for i in second derivative matrix (f for full; otherwise local)
!  iy(f)             = position of y coordinate for i in second derivative matrix (f for full; otherwise local)
!  iz(f)             = position of z coordinate for i in second derivative matrix (f for full; otherwise local)
!  jx(f)             = position of x coordinate for j in second derivative matrix (f for full; otherwise local)
!  jy(f)             = position of y coordinate for j in second derivative matrix (f for full; otherwise local)
!  jz(f)             = position of z coordinate for j in second derivative matrix (f for full; otherwise local)
!  lilocal           = if .true. then i is local to this processor
!  ljlocal           = if .true. then j is local to this processor
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d1trm             = first derivative of energy w.r.t. rij divided by rij
!  d2trm             = derivative of d1trm w.r.t. rij divided by rij
!  ladd2sderv2       = if .true. then add terms to sderv2
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_1g
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
  use current,        only : nstrains, nstrptr
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: ladd2sderv2
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: dr2ijds(6)
  real(dp),    intent(in)                          :: d2r2ijdx2(6)
  real(dp),    intent(in)                          :: d2r2ijds2(6,6)
  real(dp),    intent(in)                          :: d2r2ijdsdx(6,3)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_1gdm')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) - d1trm*d2r2ijdx2(1) - d2trm*xij*xij
    derv2(jyf,ix) = derv2(jyf,ix) - d1trm*d2r2ijdx2(6) - d2trm*yij*xij
    derv2(jzf,ix) = derv2(jzf,ix) - d1trm*d2r2ijdx2(5) - d2trm*zij*xij
    derv2(jxf,iy) = derv2(jxf,iy) - d1trm*d2r2ijdx2(6) - d2trm*xij*yij
    derv2(jyf,iy) = derv2(jyf,iy) - d1trm*d2r2ijdx2(2) - d2trm*yij*yij
    derv2(jzf,iy) = derv2(jzf,iy) - d1trm*d2r2ijdx2(4) - d2trm*zij*yij
    derv2(jxf,iz) = derv2(jxf,iz) - d1trm*d2r2ijdx2(5) - d2trm*xij*zij
    derv2(jyf,iz) = derv2(jyf,iz) - d1trm*d2r2ijdx2(4) - d2trm*yij*zij
    derv2(jzf,iz) = derv2(jzf,iz) - d1trm*d2r2ijdx2(3) - d2trm*zij*zij
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) - d1trm*d2r2ijdx2(1) - d2trm*xij*xij
    derv2(iyf,jx) = derv2(iyf,jx) - d1trm*d2r2ijdx2(6) - d2trm*yij*xij
    derv2(izf,jx) = derv2(izf,jx) - d1trm*d2r2ijdx2(5) - d2trm*zij*xij
    derv2(ixf,jy) = derv2(ixf,jy) - d1trm*d2r2ijdx2(6) - d2trm*xij*yij
    derv2(iyf,jy) = derv2(iyf,jy) - d1trm*d2r2ijdx2(2) - d2trm*yij*yij
    derv2(izf,jy) = derv2(izf,jy) - d1trm*d2r2ijdx2(4) - d2trm*zij*yij
    derv2(ixf,jz) = derv2(ixf,jz) - d1trm*d2r2ijdx2(5) - d2trm*xij*zij
    derv2(iyf,jz) = derv2(iyf,jz) - d1trm*d2r2ijdx2(4) - d2trm*yij*zij
    derv2(izf,jz) = derv2(izf,jz) - d1trm*d2r2ijdx2(3) - d2trm*zij*zij
  endif
!
  if (lstr) then
!
!  Mixed derivatives
!
    if (lilocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ijdsdx(ks,1) - xij*d2trm*dr2ijds(ks)
        derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ijdsdx(ks,2) - yij*d2trm*dr2ijds(ks)
        derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ijdsdx(ks,3) - zij*d2trm*dr2ijds(ks)
      enddo
    endif
!
    if (ljlocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(jx,kl) = derv3(jx,kl) + d1trm*d2r2ijdsdx(ks,1) + xij*d2trm*dr2ijds(ks)
        derv3(jy,kl) = derv3(jy,kl) + d1trm*d2r2ijdsdx(ks,2) + yij*d2trm*dr2ijds(ks)
        derv3(jz,kl) = derv3(jz,kl) + d1trm*d2r2ijdsdx(ks,3) + zij*d2trm*dr2ijds(ks)
      enddo
    endif
!
!  Strain-strain
!
    if (ladd2sderv2) then
      do kk = 1,nstrains
        ks = nstrptr(kk)
        do kl = 1,nstrains
          kt = nstrptr(kl)
          sderv2(kl,kk) = sderv2(kl,kk) + d1trm*d2r2ijds2(kt,ks)
        enddo
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_1gdm')
#endif
!
  return
  end
!
  subroutine add_drv2_2sgdm(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d1trm,dr2ik,dr2jl,dr2ikds,dr2jlds, &
                            d2trm,d2r2ikjl,d2r2ikjlds2,d2r2dsdik,d2r2dsdjl,ladd2sderv2)
!
!  Adds the second derivatives for a single distance that depends on two components, i-k and j-l
!  Generalised version where matrices with the terms are passed in
!  Distributed memory parallel version.
!
!  On entry : 
!
!  ix(f)             = position of x coordinate for i in second derivative matrix (with f is in full matrix; without is local)
!  iy(f)             = position of y coordinate for i in second derivative matrix (with f is in full matrix; without is local)
!  iz(f)             = position of z coordinate for i in second derivative matrix (with f is in full matrix; without is local)
!  jx(f)             = position of x coordinate for j in second derivative matrix (with f is in full matrix; without is local)
!  jy(f)             = position of y coordinate for j in second derivative matrix (with f is in full matrix; without is local)
!  jz(f)             = position of z coordinate for j in second derivative matrix (with f is in full matrix; without is local)
!  kx(f)             = position of x coordinate for k in second derivative matrix (with f is in full matrix; without is local)
!  ky(f)             = position of y coordinate for k in second derivative matrix (with f is in full matrix; without is local)
!  kz(f)             = position of z coordinate for k in second derivative matrix (with f is in full matrix; without is local)
!  lx(f)             = position of x coordinate for l in second derivative matrix (with f is in full matrix; without is local)
!  ly(f)             = position of y coordinate for l in second derivative matrix (with f is in full matrix; without is local)
!  lz(f)             = position of z coordinate for l in second derivative matrix (with f is in full matrix; without is local)
!  lilocal           = if .true. then i is local to this processor
!  ljlocal           = if .true. then j is local to this processor
!  lklocal           = if .true. then k is local to this processor
!  lllocal           = if .true. then l is local to this processor
!  ladd2sderv2       = if .true. then add terms to sderv2
!  d1trm             = (1/r dE/dr)
!  d2trm             = (1/r d/dr (1/r dE/dr))
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_2sg
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
  use current,        only : nstrains, nstrptr
  use derivatives
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  logical,     intent(in)                          :: ladd2sderv2
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: dr2ik(3)
  real(dp),    intent(in)                          :: d2r2ikjl(3,3)
  real(dp),    intent(in)                          :: d2r2ikjlds2(6,6)
  real(dp),    intent(in)                          :: d2r2dsdik(6,3)
  real(dp),    intent(in)                          :: d2r2dsdjl(6,3)
  real(dp),    intent(in)                          :: dr2jl(3)
  real(dp),    intent(in)                          :: dr2ikds(6)
  real(dp),    intent(in)                          :: dr2jlds(6)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_2sgdm')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(jyf,ix) = derv2(jyf,ix) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(jzf,ix) = derv2(jzf,ix) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(jxf,iy) = derv2(jxf,iy) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(jyf,iy) = derv2(jyf,iy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(jzf,iy) = derv2(jzf,iy) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(jxf,iz) = derv2(jxf,iz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(jyf,iz) = derv2(jyf,iz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(jzf,iz) = derv2(jzf,iz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(iyf,jx) = derv2(iyf,jx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(izf,jx) = derv2(izf,jx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ixf,jy) = derv2(ixf,jy) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(iyf,jy) = derv2(iyf,jy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(izf,jy) = derv2(izf,jy) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(ixf,jz) = derv2(ixf,jz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(iyf,jz) = derv2(iyf,jz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(izf,jz) = derv2(izf,jz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  I-L
!
  if (lilocal) then
    derv2(lxf,ix) = derv2(lxf,ix) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(lyf,ix) = derv2(lyf,ix) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(lzf,ix) = derv2(lzf,ix) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(lxf,iy) = derv2(lxf,iy) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(lyf,iy) = derv2(lyf,iy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(lzf,iy) = derv2(lzf,iy) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(lxf,iz) = derv2(lxf,iz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(lyf,iz) = derv2(lyf,iz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(lzf,iz) = derv2(lzf,iz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lllocal) then
    derv2(ixf,lx) = derv2(ixf,lx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(iyf,lx) = derv2(iyf,lx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(izf,lx) = derv2(izf,lx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ixf,ly) = derv2(ixf,ly) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(iyf,ly) = derv2(iyf,ly) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(izf,ly) = derv2(izf,ly) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(ixf,lz) = derv2(ixf,lz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(iyf,lz) = derv2(iyf,lz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(izf,lz) = derv2(izf,lz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  J-K
!
  if (ljlocal) then
    derv2(kxf,jx) = derv2(kxf,jx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(kyf,jx) = derv2(kyf,jx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(kzf,jx) = derv2(kzf,jx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(kxf,jy) = derv2(kxf,jy) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(kyf,jy) = derv2(kyf,jy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(kzf,jy) = derv2(kzf,jy) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(kxf,jz) = derv2(kxf,jz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(kyf,jz) = derv2(kyf,jz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(kzf,jz) = derv2(kzf,jz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lklocal) then
    derv2(jxf,kx) = derv2(jxf,kx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(jyf,kx) = derv2(jyf,kx) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(jzf,kx) = derv2(jzf,kx) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(jxf,ky) = derv2(jxf,ky) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(jyf,ky) = derv2(jyf,ky) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(jzf,ky) = derv2(jzf,ky) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(jxf,kz) = derv2(jxf,kz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(jyf,kz) = derv2(jyf,kz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(jzf,kz) = derv2(jzf,kz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  K-L
!     
  if (lklocal) then
    derv2(lxf,kx) = derv2(lxf,kx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(lyf,kx) = derv2(lyf,kx) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(lzf,kx) = derv2(lzf,kx) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(lxf,ky) = derv2(lxf,ky) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(lyf,ky) = derv2(lyf,ky) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(lzf,ky) = derv2(lzf,ky) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(lxf,kz) = derv2(lxf,kz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(lyf,kz) = derv2(lyf,kz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(lzf,kz) = derv2(lzf,kz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lllocal) then
    derv2(kxf,lx) = derv2(kxf,lx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(kyf,lx) = derv2(kyf,lx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(kzf,lx) = derv2(kzf,lx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(kxf,ly) = derv2(kxf,ly) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(kyf,ly) = derv2(kyf,ly) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(kzf,ly) = derv2(kzf,ly) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(kxf,lz) = derv2(kxf,lz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(kyf,lz) = derv2(kyf,lz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(kzf,lz) = derv2(kzf,lz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  Strain terms
!
  if (lstr) then
    if (lilocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(ix,kl) = derv3(ix,kl) - d2trm*dr2jlds(ks)*dr2ik(1) - d1trm*d2r2dsdik(ks,1)
        derv3(iy,kl) = derv3(iy,kl) - d2trm*dr2jlds(ks)*dr2ik(2) - d1trm*d2r2dsdik(ks,2)
        derv3(iz,kl) = derv3(iz,kl) - d2trm*dr2jlds(ks)*dr2ik(3) - d1trm*d2r2dsdik(ks,3)
      enddo
    endif
!
    if (lklocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(kx,kl) = derv3(kx,kl) + d2trm*dr2jlds(ks)*dr2ik(1) + d1trm*d2r2dsdik(ks,1)
        derv3(ky,kl) = derv3(ky,kl) + d2trm*dr2jlds(ks)*dr2ik(2) + d1trm*d2r2dsdik(ks,2)
        derv3(kz,kl) = derv3(kz,kl) + d2trm*dr2jlds(ks)*dr2ik(3) + d1trm*d2r2dsdik(ks,3)
      enddo
    endif
!
    if (ljlocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(jx,kl) = derv3(jx,kl) - d2trm*dr2ikds(ks)*dr2jl(1) - d1trm*d2r2dsdjl(ks,1)
        derv3(jy,kl) = derv3(jy,kl) - d2trm*dr2ikds(ks)*dr2jl(2) - d1trm*d2r2dsdjl(ks,2)
        derv3(jz,kl) = derv3(jz,kl) - d2trm*dr2ikds(ks)*dr2jl(3) - d1trm*d2r2dsdjl(ks,3)
      enddo
    endif
!
    if (lllocal) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(lx,kl) = derv3(lx,kl) + d2trm*dr2ikds(ks)*dr2jl(1) + d1trm*d2r2dsdjl(ks,1)
        derv3(ly,kl) = derv3(ly,kl) + d2trm*dr2ikds(ks)*dr2jl(2) + d1trm*d2r2dsdjl(ks,2)
        derv3(lz,kl) = derv3(lz,kl) + d2trm*dr2ikds(ks)*dr2jl(3) + d1trm*d2r2dsdjl(ks,3)
      enddo
    endif
!
    if (ladd2sderv2) then
      do kk = 1,nstrains
        ks = nstrptr(kk)
        do kl = 1,nstrains
          kt = nstrptr(kl)
          sderv2(kl,kk) = sderv2(kl,kk) + d2trm*(dr2ikds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2ikds(ks)) + &
                                          d1trm*d2r2ikjlds2(kt,ks)
        enddo
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_2sgdm')
#endif
!
  return
  end
!
  subroutine add_drv2_1p(ix,iy,iz,jx,jy,jz,xij,yij,zij,d1trm,d2trm,d2r2ijdx2)
!
!  Adds the second derivatives for a single distance: Phonon version
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d1trm             = first derivative of energy w.r.t. rij divided by rij
!  d2trm             = derivative of d1trm w.r.t. rij divided by rij
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/21 Created from add_drv2_1
!   5/23 Exclusion of duplicate on-diagonal blocks added
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
  use datatypes
  use derivatives
  use m_strain,       only : real1strterm
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: d2r2ijdx2(6)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_1p')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (jx.gt.ix) then
    derv2(jx,ix) = derv2(jx,ix) - d2trm*d2r2ijdx2(1)
    derv2(jy,ix) = derv2(jy,ix) - d2trm*d2r2ijdx2(6)
    derv2(jz,ix) = derv2(jz,ix) - d2trm*d2r2ijdx2(5)
    derv2(jx,iy) = derv2(jx,iy) - d2trm*d2r2ijdx2(6)
    derv2(jy,iy) = derv2(jy,iy) - d2trm*d2r2ijdx2(2)
    derv2(jz,iy) = derv2(jz,iy) - d2trm*d2r2ijdx2(4)
    derv2(jx,iz) = derv2(jx,iz) - d2trm*d2r2ijdx2(5)
    derv2(jy,iz) = derv2(jy,iz) - d2trm*d2r2ijdx2(4)
    derv2(jz,iz) = derv2(jz,iz) - d2trm*d2r2ijdx2(3)
    derv2(jx,ix) = derv2(jx,ix) - d1trm
    derv2(jy,iy) = derv2(jy,iy) - d1trm
    derv2(jz,iz) = derv2(jz,iz) - d1trm
  elseif (jx.lt.ix) then
    derv2(ix,jx) = derv2(ix,jx) - d2trm*d2r2ijdx2(1)
    derv2(iy,jx) = derv2(iy,jx) - d2trm*d2r2ijdx2(6)
    derv2(iz,jx) = derv2(iz,jx) - d2trm*d2r2ijdx2(5)
    derv2(ix,jy) = derv2(ix,jy) - d2trm*d2r2ijdx2(6)
    derv2(iy,jy) = derv2(iy,jy) - d2trm*d2r2ijdx2(2)
    derv2(iz,jy) = derv2(iz,jy) - d2trm*d2r2ijdx2(4)
    derv2(ix,jz) = derv2(ix,jz) - d2trm*d2r2ijdx2(5)
    derv2(iy,jz) = derv2(iy,jz) - d2trm*d2r2ijdx2(4)
    derv2(iz,jz) = derv2(iz,jz) - d2trm*d2r2ijdx2(3)
    derv2(ix,jx) = derv2(ix,jx) - d1trm
    derv2(iy,jy) = derv2(iy,jy) - d1trm
    derv2(iz,jz) = derv2(iz,jz) - d1trm
  endif
#ifdef TRACE
  call trace_out('add_drv2_1p')
#endif
!
  return
  end
!
  subroutine add_drv2_2p(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,xik,yik,zik,xjl,yjl,zjl,d2trm)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l : Phonon version
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  xik               = x component of vector from i to k
!  yik               = y component of vector from i to k
!  zik               = z component of vector from i to k
!  xjl               = x component of vector from i to k
!  yjl               = y component of vector from i to k
!  zjl               = z component of vector from i to k
!  d2trm             = second derivative of energy w.r.t. rik and rjl with division by distances
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/21 Created from add_drv2_2
!   5/23 Diagonal blocks explicitly excluded
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
  use datatypes
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xik
  real(dp),    intent(in)                          :: yik
  real(dp),    intent(in)                          :: zik
  real(dp),    intent(in)                          :: xjl
  real(dp),    intent(in)                          :: yjl
  real(dp),    intent(in)                          :: zjl
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2p')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.gt.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d2trm*xjl*xik
    derv2(jy,ix) = derv2(jy,ix) + d2trm*yjl*xik
    derv2(jz,ix) = derv2(jz,ix) + d2trm*zjl*xik
    derv2(jx,iy) = derv2(jx,iy) + d2trm*xjl*yik
    derv2(jy,iy) = derv2(jy,iy) + d2trm*yjl*yik
    derv2(jz,iy) = derv2(jz,iy) + d2trm*zjl*yik
    derv2(jx,iz) = derv2(jx,iz) + d2trm*xjl*zik
    derv2(jy,iz) = derv2(jy,iz) + d2trm*yjl*zik
    derv2(jz,iz) = derv2(jz,iz) + d2trm*zjl*zik
  elseif (jx.lt.ix) then
    derv2(ix,jx) = derv2(ix,jx) + d2trm*xjl*xik
    derv2(iy,jx) = derv2(iy,jx) + d2trm*xjl*yik
    derv2(iz,jx) = derv2(iz,jx) + d2trm*xjl*zik
    derv2(ix,jy) = derv2(ix,jy) + d2trm*yjl*xik
    derv2(iy,jy) = derv2(iy,jy) + d2trm*yjl*yik
    derv2(iz,jy) = derv2(iz,jy) + d2trm*yjl*zik
    derv2(ix,jz) = derv2(ix,jz) + d2trm*zjl*xik
    derv2(iy,jz) = derv2(iy,jz) + d2trm*zjl*yik
    derv2(iz,jz) = derv2(iz,jz) + d2trm*zjl*zik
  endif
!
!  I-L
!
  if (lx.gt.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d2trm*xjl*xik
    derv2(ly,ix) = derv2(ly,ix) - d2trm*yjl*xik
    derv2(lz,ix) = derv2(lz,ix) - d2trm*zjl*xik
    derv2(lx,iy) = derv2(lx,iy) - d2trm*xjl*yik
    derv2(ly,iy) = derv2(ly,iy) - d2trm*yjl*yik
    derv2(lz,iy) = derv2(lz,iy) - d2trm*zjl*yik
    derv2(lx,iz) = derv2(lx,iz) - d2trm*xjl*zik
    derv2(ly,iz) = derv2(ly,iz) - d2trm*yjl*zik
    derv2(lz,iz) = derv2(lz,iz) - d2trm*zjl*zik
  elseif (lx.lt.ix) then
    derv2(ix,lx) = derv2(ix,lx) - d2trm*xjl*xik
    derv2(iy,lx) = derv2(iy,lx) - d2trm*xjl*yik
    derv2(iz,lx) = derv2(iz,lx) - d2trm*xjl*zik
    derv2(ix,ly) = derv2(ix,ly) - d2trm*yjl*xik
    derv2(iy,ly) = derv2(iy,ly) - d2trm*yjl*yik
    derv2(iz,ly) = derv2(iz,ly) - d2trm*yjl*zik
    derv2(ix,lz) = derv2(ix,lz) - d2trm*zjl*xik
    derv2(iy,lz) = derv2(iy,lz) - d2trm*zjl*yik
    derv2(iz,lz) = derv2(iz,lz) - d2trm*zjl*zik
  endif
!
!  J-K
!
  if (kx.gt.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d2trm*xjl*xik
    derv2(ky,jx) = derv2(ky,jx) - d2trm*xjl*yik
    derv2(kz,jx) = derv2(kz,jx) - d2trm*xjl*zik
    derv2(kx,jy) = derv2(kx,jy) - d2trm*yjl*xik
    derv2(ky,jy) = derv2(ky,jy) - d2trm*yjl*yik
    derv2(kz,jy) = derv2(kz,jy) - d2trm*yjl*zik
    derv2(kx,jz) = derv2(kx,jz) - d2trm*zjl*xik
    derv2(ky,jz) = derv2(ky,jz) - d2trm*zjl*yik
    derv2(kz,jz) = derv2(kz,jz) - d2trm*zjl*zik
  elseif (kx.lt.jx) then
    derv2(jx,kx) = derv2(jx,kx) - d2trm*xjl*xik
    derv2(jy,kx) = derv2(jy,kx) - d2trm*yjl*xik
    derv2(jz,kx) = derv2(jz,kx) - d2trm*zjl*xik
    derv2(jx,ky) = derv2(jx,ky) - d2trm*xjl*yik
    derv2(jy,ky) = derv2(jy,ky) - d2trm*yjl*yik
    derv2(jz,ky) = derv2(jz,ky) - d2trm*zjl*yik
    derv2(jx,kz) = derv2(jx,kz) - d2trm*xjl*zik
    derv2(jy,kz) = derv2(jy,kz) - d2trm*yjl*zik
    derv2(jz,kz) = derv2(jz,kz) - d2trm*zjl*zik
  endif
!
!  K-L
!     
  if (lx.gt.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d2trm*xik*xjl
    derv2(ly,kx) = derv2(ly,kx) + d2trm*xik*yjl
    derv2(lz,kx) = derv2(lz,kx) + d2trm*xik*zjl
    derv2(lx,ky) = derv2(lx,ky) + d2trm*yik*xjl
    derv2(ly,ky) = derv2(ly,ky) + d2trm*yik*yjl
    derv2(lz,ky) = derv2(lz,ky) + d2trm*yik*zjl
    derv2(lx,kz) = derv2(lx,kz) + d2trm*zik*xjl
    derv2(ly,kz) = derv2(ly,kz) + d2trm*zik*yjl
    derv2(lz,kz) = derv2(lz,kz) + d2trm*zik*zjl
  elseif (lx.lt.kx) then
    derv2(kx,lx) = derv2(kx,lx) + d2trm*xik*xjl
    derv2(ky,lx) = derv2(ky,lx) + d2trm*yik*xjl
    derv2(kz,lx) = derv2(kz,lx) + d2trm*zik*xjl
    derv2(kx,ly) = derv2(kx,ly) + d2trm*xik*yjl
    derv2(ky,ly) = derv2(ky,ly) + d2trm*yik*yjl
    derv2(kz,ly) = derv2(kz,ly) + d2trm*zik*yjl
    derv2(kx,lz) = derv2(kx,lz) + d2trm*xik*zjl
    derv2(ky,lz) = derv2(ky,lz) + d2trm*yik*zjl
    derv2(kz,lz) = derv2(kz,lz) + d2trm*zik*zjl
  endif
#ifdef TRACE
  call trace_out('add_drv2_2p')
#endif
!
  return
  end
!
  subroutine add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                          d1trm,d2trm,d2r2ijdx2)
!
!  Adds the second derivatives for a single distance: Phonon version
!
!  On entry : 
!
!  lilocal           = is i local to this node?
!  ljlocal           = is j local to this node?
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d1trm             = first derivative of energy w.r.t. rij divided by rij
!  d2trm             = derivative of d1trm w.r.t. rij divided by rij
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_1p
!   5/23 Exclusion of duplicate on-diagonal blocks added
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
  use datatypes
  use derivatives
  use m_strain,       only : real1strterm
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: d2r2ijdx2(6)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_1pd')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (ixf.ne.jxf) then
    if (lilocal) then
      derv2(jxf,ix) = derv2(jxf,ix) - d2trm*d2r2ijdx2(1)
      derv2(jyf,ix) = derv2(jyf,ix) - d2trm*d2r2ijdx2(6)
      derv2(jzf,ix) = derv2(jzf,ix) - d2trm*d2r2ijdx2(5)
      derv2(jxf,iy) = derv2(jxf,iy) - d2trm*d2r2ijdx2(6)
      derv2(jyf,iy) = derv2(jyf,iy) - d2trm*d2r2ijdx2(2)
      derv2(jzf,iy) = derv2(jzf,iy) - d2trm*d2r2ijdx2(4)
      derv2(jxf,iz) = derv2(jxf,iz) - d2trm*d2r2ijdx2(5)
      derv2(jyf,iz) = derv2(jyf,iz) - d2trm*d2r2ijdx2(4)
      derv2(jzf,iz) = derv2(jzf,iz) - d2trm*d2r2ijdx2(3)
      derv2(jxf,ix) = derv2(jxf,ix) - d1trm
      derv2(jyf,iy) = derv2(jyf,iy) - d1trm
      derv2(jzf,iz) = derv2(jzf,iz) - d1trm
    endif
!
!  Exclude duplication of on-diagonal blocks
!
    if (ljlocal) then
      derv2(ixf,jx) = derv2(ixf,jx) - d2trm*d2r2ijdx2(1)
      derv2(iyf,jx) = derv2(iyf,jx) - d2trm*d2r2ijdx2(6)
      derv2(izf,jx) = derv2(izf,jx) - d2trm*d2r2ijdx2(5)
      derv2(ixf,jy) = derv2(ixf,jy) - d2trm*d2r2ijdx2(6)
      derv2(iyf,jy) = derv2(iyf,jy) - d2trm*d2r2ijdx2(2)
      derv2(izf,jy) = derv2(izf,jy) - d2trm*d2r2ijdx2(4)
      derv2(ixf,jz) = derv2(ixf,jz) - d2trm*d2r2ijdx2(5)
      derv2(iyf,jz) = derv2(iyf,jz) - d2trm*d2r2ijdx2(4)
      derv2(izf,jz) = derv2(izf,jz) - d2trm*d2r2ijdx2(3)
      derv2(ixf,jx) = derv2(ixf,jx) - d1trm
      derv2(iyf,jy) = derv2(iyf,jy) - d1trm
      derv2(izf,jz) = derv2(izf,jz) - d1trm
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_1pd')
#endif
!
  return
  end
!
  subroutine add_drv2_2pd(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xjl,yjl,zjl,d2trm)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l : Phonon version
!
!  On entry : 
!
!  lilocal           = is i local to this node?
!  ljlocal           = is j local to this node?
!  lklocal           = is k local to this node?
!  lllocal           = is l local to this node?
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  xik               = x component of vector from i to k
!  yik               = y component of vector from i to k
!  zik               = z component of vector from i to k
!  xjl               = x component of vector from i to k
!  yjl               = y component of vector from i to k
!  zjl               = z component of vector from i to k
!  d2trm             = second derivative of energy w.r.t. rik and rjl with division by distances
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_2p
!   5/23 Diagonal blocks explicitly excluded
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
  use datatypes
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xik
  real(dp),    intent(in)                          :: yik
  real(dp),    intent(in)                          :: zik
  real(dp),    intent(in)                          :: xjl
  real(dp),    intent(in)                          :: yjl
  real(dp),    intent(in)                          :: zjl
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2pd')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (ixf.ne.jxf) then
    if (lilocal) then
      derv2(jxf,ix) = derv2(jxf,ix) + d2trm*xjl*xik
      derv2(jyf,ix) = derv2(jyf,ix) + d2trm*yjl*xik
      derv2(jzf,ix) = derv2(jzf,ix) + d2trm*zjl*xik
      derv2(jxf,iy) = derv2(jxf,iy) + d2trm*xjl*yik
      derv2(jyf,iy) = derv2(jyf,iy) + d2trm*yjl*yik
      derv2(jzf,iy) = derv2(jzf,iy) + d2trm*zjl*yik
      derv2(jxf,iz) = derv2(jxf,iz) + d2trm*xjl*zik
      derv2(jyf,iz) = derv2(jyf,iz) + d2trm*yjl*zik
      derv2(jzf,iz) = derv2(jzf,iz) + d2trm*zjl*zik
    endif
!
!  Exclude duplication of on-diagonal blocks
!
    if (ljlocal) then
      derv2(ixf,jx) = derv2(ixf,jx) + d2trm*xjl*xik
      derv2(iyf,jx) = derv2(iyf,jx) + d2trm*xjl*yik
      derv2(izf,jx) = derv2(izf,jx) + d2trm*xjl*zik
      derv2(ixf,jy) = derv2(ixf,jy) + d2trm*yjl*xik
      derv2(iyf,jy) = derv2(iyf,jy) + d2trm*yjl*yik
      derv2(izf,jy) = derv2(izf,jy) + d2trm*yjl*zik
      derv2(ixf,jz) = derv2(ixf,jz) + d2trm*zjl*xik
      derv2(iyf,jz) = derv2(iyf,jz) + d2trm*zjl*yik
      derv2(izf,jz) = derv2(izf,jz) + d2trm*zjl*zik
    endif
  endif
!
!  I-L
!
  if (ixf.ne.lxf) then
    if (lilocal) then
      derv2(lxf,ix) = derv2(lxf,ix) - d2trm*xjl*xik
      derv2(lyf,ix) = derv2(lyf,ix) - d2trm*yjl*xik
      derv2(lzf,ix) = derv2(lzf,ix) - d2trm*zjl*xik
      derv2(lxf,iy) = derv2(lxf,iy) - d2trm*xjl*yik
      derv2(lyf,iy) = derv2(lyf,iy) - d2trm*yjl*yik
      derv2(lzf,iy) = derv2(lzf,iy) - d2trm*zjl*yik
      derv2(lxf,iz) = derv2(lxf,iz) - d2trm*xjl*zik
      derv2(lyf,iz) = derv2(lyf,iz) - d2trm*yjl*zik
      derv2(lzf,iz) = derv2(lzf,iz) - d2trm*zjl*zik
    endif
!
!  Exclude duplication of on-diagonal blocks
!
    if (lllocal) then
      derv2(ixf,lx) = derv2(ixf,lx) - d2trm*xjl*xik
      derv2(iyf,lx) = derv2(iyf,lx) - d2trm*xjl*yik
      derv2(izf,lx) = derv2(izf,lx) - d2trm*xjl*zik
      derv2(ixf,ly) = derv2(ixf,ly) - d2trm*yjl*xik
      derv2(iyf,ly) = derv2(iyf,ly) - d2trm*yjl*yik
      derv2(izf,ly) = derv2(izf,ly) - d2trm*yjl*zik
      derv2(ixf,lz) = derv2(ixf,lz) - d2trm*zjl*xik
      derv2(iyf,lz) = derv2(iyf,lz) - d2trm*zjl*yik
      derv2(izf,lz) = derv2(izf,lz) - d2trm*zjl*zik
    endif
  endif
!
!  J-K
!
  if (jxf.ne.kxf) then
    if (ljlocal) then
      derv2(kxf,jx) = derv2(kxf,jx) - d2trm*xjl*xik
      derv2(kyf,jx) = derv2(kyf,jx) - d2trm*xjl*yik
      derv2(kzf,jx) = derv2(kzf,jx) - d2trm*xjl*zik
      derv2(kxf,jy) = derv2(kxf,jy) - d2trm*yjl*xik
      derv2(kyf,jy) = derv2(kyf,jy) - d2trm*yjl*yik
      derv2(kzf,jy) = derv2(kzf,jy) - d2trm*yjl*zik
      derv2(kxf,jz) = derv2(kxf,jz) - d2trm*zjl*xik
      derv2(kyf,jz) = derv2(kyf,jz) - d2trm*zjl*yik
      derv2(kzf,jz) = derv2(kzf,jz) - d2trm*zjl*zik
    endif
!
!  Exclude duplication of on-diagonal blocks
!
    if (lklocal) then
      derv2(jxf,kx) = derv2(jxf,kx) - d2trm*xjl*xik
      derv2(jyf,kx) = derv2(jyf,kx) - d2trm*yjl*xik
      derv2(jzf,kx) = derv2(jzf,kx) - d2trm*zjl*xik
      derv2(jxf,ky) = derv2(jxf,ky) - d2trm*xjl*yik
      derv2(jyf,ky) = derv2(jyf,ky) - d2trm*yjl*yik
      derv2(jzf,ky) = derv2(jzf,ky) - d2trm*zjl*yik
      derv2(jxf,kz) = derv2(jxf,kz) - d2trm*xjl*zik
      derv2(jyf,kz) = derv2(jyf,kz) - d2trm*yjl*zik
      derv2(jzf,kz) = derv2(jzf,kz) - d2trm*zjl*zik
    endif
  endif
!
!  K-L
!     
  if (kxf.ne.lxf) then
    if (lklocal) then
      derv2(lxf,kx) = derv2(lxf,kx) + d2trm*xik*xjl
      derv2(lyf,kx) = derv2(lyf,kx) + d2trm*xik*yjl
      derv2(lzf,kx) = derv2(lzf,kx) + d2trm*xik*zjl
      derv2(lxf,ky) = derv2(lxf,ky) + d2trm*yik*xjl
      derv2(lyf,ky) = derv2(lyf,ky) + d2trm*yik*yjl
      derv2(lzf,ky) = derv2(lzf,ky) + d2trm*yik*zjl
      derv2(lxf,kz) = derv2(lxf,kz) + d2trm*zik*xjl
      derv2(lyf,kz) = derv2(lyf,kz) + d2trm*zik*yjl
      derv2(lzf,kz) = derv2(lzf,kz) + d2trm*zik*zjl
    endif
!
!  Exclude duplication of on-diagonal blocks
!
    if (lllocal) then
      derv2(kxf,lx) = derv2(kxf,lx) + d2trm*xik*xjl
      derv2(kyf,lx) = derv2(kyf,lx) + d2trm*yik*xjl
      derv2(kzf,lx) = derv2(kzf,lx) + d2trm*zik*xjl
      derv2(kxf,ly) = derv2(kxf,ly) + d2trm*xik*yjl
      derv2(kyf,ly) = derv2(kyf,ly) + d2trm*yik*yjl
      derv2(kzf,ly) = derv2(kzf,ly) + d2trm*zik*yjl
      derv2(kxf,lz) = derv2(kxf,lz) + d2trm*xik*zjl
      derv2(kyf,lz) = derv2(kyf,lz) + d2trm*yik*zjl
      derv2(kzf,lz) = derv2(kzf,lz) + d2trm*zik*zjl
    endif
  endif
#ifdef TRACE
  call trace_out('add_drv2_2pd')
#endif
!
  return
  end
!
  subroutine add_drv2_1bp(ix,iy,iz,jx,jy,jz,d1trm,d2trm,d1c,d2c)
!
!  Adds the second derivatives for a single distance : Phonon version
!  Version that uses pre-computed block
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  d1trm             = multiplier for second derivatives
!  d2trm             = multiplier for first derivative products
!  d1c               = first derivatives - Cartesian
!  d2c               = second derivatives - Cartesian/Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/21 Created from add_drv2_1b
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1c(3)
  real(dp),    intent(in)                          :: d2c(3,3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_1bp')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) - d1trm*d2c(1,1) - d2trm*d1c(1)*d1c(1)
    derv2(jy,ix) = derv2(jy,ix) - d1trm*d2c(2,1) - d2trm*d1c(2)*d1c(1)
    derv2(jz,ix) = derv2(jz,ix) - d1trm*d2c(3,1) - d2trm*d1c(3)*d1c(1)
    derv2(jx,iy) = derv2(jx,iy) - d1trm*d2c(1,2) - d2trm*d1c(1)*d1c(2)
    derv2(jy,iy) = derv2(jy,iy) - d1trm*d2c(2,2) - d2trm*d1c(2)*d1c(2)
    derv2(jz,iy) = derv2(jz,iy) - d1trm*d2c(3,2) - d2trm*d1c(3)*d1c(2)
    derv2(jx,iz) = derv2(jx,iz) - d1trm*d2c(1,3) - d2trm*d1c(1)*d1c(3)
    derv2(jy,iz) = derv2(jy,iz) - d1trm*d2c(2,3) - d2trm*d1c(2)*d1c(3)
    derv2(jz,iz) = derv2(jz,iz) - d1trm*d2c(3,3) - d2trm*d1c(3)*d1c(3)
  else
    derv2(ix,jx) = derv2(ix,jx) - d1trm*d2c(1,1) - d2trm*d1c(1)*d1c(1)
    derv2(iy,jx) = derv2(iy,jx) - d1trm*d2c(2,1) - d2trm*d1c(2)*d1c(1)
    derv2(iz,jx) = derv2(iz,jx) - d1trm*d2c(3,1) - d2trm*d1c(3)*d1c(1)
    derv2(ix,jy) = derv2(ix,jy) - d1trm*d2c(1,2) - d2trm*d1c(1)*d1c(2)
    derv2(iy,jy) = derv2(iy,jy) - d1trm*d2c(2,2) - d2trm*d1c(2)*d1c(2)
    derv2(iz,jy) = derv2(iz,jy) - d1trm*d2c(3,2) - d2trm*d1c(3)*d1c(2)
    derv2(ix,jz) = derv2(ix,jz) - d1trm*d2c(1,3) - d2trm*d1c(1)*d1c(3)
    derv2(iy,jz) = derv2(iy,jz) - d1trm*d2c(2,3) - d2trm*d1c(2)*d1c(3)
    derv2(iz,jz) = derv2(iz,jz) - d1trm*d2c(3,3) - d2trm*d1c(3)*d1c(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_1bp')
#endif
!
  return
  end
!
  subroutine add_drv2_2bp(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,d1trm,d2trm,d1ikc,d1jlc,d2c)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l : Phonon version
!  Version that uses pre-computed block
!  NB: Strain derivatives are assumed to be compacted to nstrains only
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d1trm             = multiplier for second derivatives
!  d2trm             = multiplier for first derivative products
!  d1ikc             = first derivatives for i-k - Cartesian
!  d1jlc             = first derivatives for j-l - Cartesian
!  d2c               = second derivatives - Cartesian/Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/21 Created from add_drv2_2b
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1ikc(3)
  real(dp),    intent(in)                          :: d1jlc(3)
  real(dp),    intent(in)                          :: d2c(3,3)
#ifdef TRACE
  call trace_in('add_drv2_2bp')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(jy,ix) = derv2(jy,ix) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(jz,ix) = derv2(jz,ix) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(jx,iy) = derv2(jx,iy) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(jy,iy) = derv2(jy,iy) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(jz,iy) = derv2(jz,iy) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(jx,iz) = derv2(jx,iz) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(jy,iz) = derv2(jy,iz) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(jz,iz) = derv2(jz,iz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(ix,jx) = derv2(ix,jx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(iy,jx) = derv2(iy,jx) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(iz,jx) = derv2(iz,jx) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ix,jy) = derv2(ix,jy) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(iy,jy) = derv2(iy,jy) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(iz,jy) = derv2(iz,jy) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(ix,jz) = derv2(ix,jz) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(iy,jz) = derv2(iy,jz) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(iz,jz) = derv2(iz,jz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  I-L
!
  if (lx.ge.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(ly,ix) = derv2(ly,ix) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(lz,ix) = derv2(lz,ix) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(lx,iy) = derv2(lx,iy) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(ly,iy) = derv2(ly,iy) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(lz,iy) = derv2(lz,iy) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(lx,iz) = derv2(lx,iz) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ly,iz) = derv2(ly,iz) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(lz,iz) = derv2(lz,iz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(ix,lx) = derv2(ix,lx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(iy,lx) = derv2(iy,lx) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(iz,lx) = derv2(iz,lx) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ix,ly) = derv2(ix,ly) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(iy,ly) = derv2(iy,ly) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(iz,ly) = derv2(iz,ly) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(ix,lz) = derv2(ix,lz) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(iy,lz) = derv2(iy,lz) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(iz,lz) = derv2(iz,lz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  J-K
!
  if (kx.ge.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(ky,jx) = derv2(ky,jx) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(kz,jx) = derv2(kz,jx) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(kx,jy) = derv2(kx,jy) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(ky,jy) = derv2(ky,jy) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(kz,jy) = derv2(kz,jy) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(kx,jz) = derv2(kx,jz) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(ky,jz) = derv2(ky,jz) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(kz,jz) = derv2(kz,jz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(jx,kx) = derv2(jx,kx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(jy,kx) = derv2(jy,kx) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(jz,kx) = derv2(jz,kx) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(jx,ky) = derv2(jx,ky) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(jy,ky) = derv2(jy,ky) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(jz,ky) = derv2(jz,ky) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(jx,kz) = derv2(jx,kz) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(jy,kz) = derv2(jy,kz) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(jz,kz) = derv2(jz,kz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  K-L
!     
  if (lx.ge.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(ly,kx) = derv2(ly,kx) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(lz,kx) = derv2(lz,kx) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(lx,ky) = derv2(lx,ky) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(ly,ky) = derv2(ly,ky) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(lz,ky) = derv2(lz,ky) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(lx,kz) = derv2(lx,kz) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ly,kz) = derv2(ly,kz) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(lz,kz) = derv2(lz,kz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(kx,lx) = derv2(kx,lx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(ky,lx) = derv2(ky,lx) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(kz,lx) = derv2(kz,lx) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(kx,ly) = derv2(kx,ly) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(ky,ly) = derv2(ky,ly) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(kz,ly) = derv2(kz,ly) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(kx,lz) = derv2(kx,lz) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(ky,lz) = derv2(ky,lz) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(kz,lz) = derv2(kz,lz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2bp')
#endif
!
  return
  end
!
  subroutine add_drv2_1bxp(ix,iy,iz,jx,jy,jz,d2trm,d1c1,d1c2)
!
!  Adds the second derivatives for a single distance based only on products
!  of first derivatives. Phonon version.
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  d2trm             = multiplier for first derivative products
!  d1c1              = first derivatives - Cartesian for first term
!  d1c2              = first derivatives - Cartesian for second term
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/21 Created from add_drv2_1bx
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1c1(3)
  real(dp),    intent(in)                          :: d1c2(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_1bxp')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) - d2trm*d1c1(1)*d1c2(1) - d2trm*d1c2(1)*d1c1(1)
    derv2(jy,ix) = derv2(jy,ix) - d2trm*d1c1(2)*d1c2(1) - d2trm*d1c2(2)*d1c1(1)
    derv2(jz,ix) = derv2(jz,ix) - d2trm*d1c1(3)*d1c2(1) - d2trm*d1c2(3)*d1c1(1)
    derv2(jx,iy) = derv2(jx,iy) - d2trm*d1c1(1)*d1c2(2) - d2trm*d1c2(1)*d1c1(2)
    derv2(jy,iy) = derv2(jy,iy) - d2trm*d1c1(2)*d1c2(2) - d2trm*d1c2(2)*d1c1(2)
    derv2(jz,iy) = derv2(jz,iy) - d2trm*d1c1(3)*d1c2(2) - d2trm*d1c2(3)*d1c1(2)
    derv2(jx,iz) = derv2(jx,iz) - d2trm*d1c1(1)*d1c2(3) - d2trm*d1c2(1)*d1c1(3)
    derv2(jy,iz) = derv2(jy,iz) - d2trm*d1c1(2)*d1c2(3) - d2trm*d1c2(2)*d1c1(3)
    derv2(jz,iz) = derv2(jz,iz) - d2trm*d1c1(3)*d1c2(3) - d2trm*d1c2(3)*d1c1(3)
  else
    derv2(ix,jx) = derv2(ix,jx) - d2trm*d1c1(1)*d1c2(1) - d2trm*d1c2(1)*d1c1(1)
    derv2(iy,jx) = derv2(iy,jx) - d2trm*d1c1(2)*d1c2(1) - d2trm*d1c2(2)*d1c1(1)
    derv2(iz,jx) = derv2(iz,jx) - d2trm*d1c1(3)*d1c2(1) - d2trm*d1c2(3)*d1c1(1)
    derv2(ix,jy) = derv2(ix,jy) - d2trm*d1c1(1)*d1c2(2) - d2trm*d1c2(1)*d1c1(2)
    derv2(iy,jy) = derv2(iy,jy) - d2trm*d1c1(2)*d1c2(2) - d2trm*d1c2(2)*d1c1(2)
    derv2(iz,jy) = derv2(iz,jy) - d2trm*d1c1(3)*d1c2(2) - d2trm*d1c2(3)*d1c1(2)
    derv2(ix,jz) = derv2(ix,jz) - d2trm*d1c1(1)*d1c2(3) - d2trm*d1c2(1)*d1c1(3)
    derv2(iy,jz) = derv2(iy,jz) - d2trm*d1c1(2)*d1c2(3) - d2trm*d1c2(2)*d1c1(3)
    derv2(iz,jz) = derv2(iz,jz) - d2trm*d1c1(3)*d1c2(3) - d2trm*d1c2(3)*d1c1(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_1bxp')
#endif
!
  return
  end
!
  subroutine add_drv2_2bxp(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,d2trm,d1ikc,d1jlc)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l: Phonon version
!  Version that only adds the product of first derivatives
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d2trm             = multiplier for first derivative products
!  d1ikc             = first derivatives for i-k - Cartesian
!  d1jlc             = first derivatives for j-l - Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/21 Created from add_drv2_2bx
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1ikc(3)
  real(dp),    intent(in)                          :: d1jlc(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2bxp')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(jy,ix) = derv2(jy,ix) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(jz,ix) = derv2(jz,ix) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(jx,iy) = derv2(jx,iy) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(jy,iy) = derv2(jy,iy) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(jz,iy) = derv2(jz,iy) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(jx,iz) = derv2(jx,iz) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(jy,iz) = derv2(jy,iz) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(jz,iz) = derv2(jz,iz) + d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(ix,jx) = derv2(ix,jx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(iy,jx) = derv2(iy,jx) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(iz,jx) = derv2(iz,jx) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ix,jy) = derv2(ix,jy) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(iy,jy) = derv2(iy,jy) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(iz,jy) = derv2(iz,jy) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(ix,jz) = derv2(ix,jz) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(iy,jz) = derv2(iy,jz) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(iz,jz) = derv2(iz,jz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  I-L
!
  if (lx.ge.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(ly,ix) = derv2(ly,ix) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(lz,ix) = derv2(lz,ix) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(lx,iy) = derv2(lx,iy) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(ly,iy) = derv2(ly,iy) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(lz,iy) = derv2(lz,iy) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(lx,iz) = derv2(lx,iz) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ly,iz) = derv2(ly,iz) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(lz,iz) = derv2(lz,iz) - d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(ix,lx) = derv2(ix,lx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(iy,lx) = derv2(iy,lx) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(iz,lx) = derv2(iz,lx) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ix,ly) = derv2(ix,ly) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(iy,ly) = derv2(iy,ly) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(iz,ly) = derv2(iz,ly) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(ix,lz) = derv2(ix,lz) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(iy,lz) = derv2(iy,lz) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(iz,lz) = derv2(iz,lz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  J-K
!
  if (kx.ge.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(ky,jx) = derv2(ky,jx) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(kz,jx) = derv2(kz,jx) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(kx,jy) = derv2(kx,jy) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(ky,jy) = derv2(ky,jy) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(kz,jy) = derv2(kz,jy) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(kx,jz) = derv2(kx,jz) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(ky,jz) = derv2(ky,jz) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(kz,jz) = derv2(kz,jz) - d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(jx,kx) = derv2(jx,kx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(jy,kx) = derv2(jy,kx) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(jz,kx) = derv2(jz,kx) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(jx,ky) = derv2(jx,ky) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(jy,ky) = derv2(jy,ky) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(jz,ky) = derv2(jz,ky) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(jx,kz) = derv2(jx,kz) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(jy,kz) = derv2(jy,kz) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(jz,kz) = derv2(jz,kz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  K-L
!     
  if (lx.ge.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(ly,kx) = derv2(ly,kx) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(lz,kx) = derv2(lz,kx) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(lx,ky) = derv2(lx,ky) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(ly,ky) = derv2(ly,ky) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(lz,ky) = derv2(lz,ky) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(lx,kz) = derv2(lx,kz) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ly,kz) = derv2(ly,kz) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(lz,kz) = derv2(lz,kz) + d2trm*d1jlc(3)*d1ikc(3)
  else
    derv2(kx,lx) = derv2(kx,lx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(ky,lx) = derv2(ky,lx) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(kz,lx) = derv2(kz,lx) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(kx,ly) = derv2(kx,ly) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(ky,ly) = derv2(ky,ly) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(kz,ly) = derv2(kz,ly) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(kx,lz) = derv2(kx,lz) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(ky,lz) = derv2(ky,lz) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(kz,lz) = derv2(kz,lz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2bxp')
#endif
!
  return
  end
!
  subroutine add_drv2_1bpd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,d1trm,d2trm,d1c,d2c)
!
!  Adds the second derivatives for a single distance : Phonon version
!  Version that uses pre-computed block
!  Distributed memory parallel version.
!
!  On entry : 
!
!  lilocal           = is i local to this node?
!  ljlocal           = is j local to this node?
!  ix(f)             = position of x coordinate for i in second derivative matrix
!  iy(f)             = position of y coordinate for i in second derivative matrix
!  iz(f)             = position of z coordinate for i in second derivative matrix
!  jx(f)             = position of x coordinate for j in second derivative matrix
!  jy(f)             = position of y coordinate for j in second derivative matrix
!  jz(f)             = position of z coordinate for j in second derivative matrix
!  d1trm             = multiplier for second derivatives
!  d2trm             = multiplier for first derivative products
!  d1c               = first derivatives - Cartesian
!  d2c               = second derivatives - Cartesian/Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_1bp
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1c(3)
  real(dp),    intent(in)                          :: d2c(3,3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_1bpd')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) - d1trm*d2c(1,1) - d2trm*d1c(1)*d1c(1)
    derv2(jyf,ix) = derv2(jyf,ix) - d1trm*d2c(2,1) - d2trm*d1c(2)*d1c(1)
    derv2(jzf,ix) = derv2(jzf,ix) - d1trm*d2c(3,1) - d2trm*d1c(3)*d1c(1)
    derv2(jxf,iy) = derv2(jxf,iy) - d1trm*d2c(1,2) - d2trm*d1c(1)*d1c(2)
    derv2(jyf,iy) = derv2(jyf,iy) - d1trm*d2c(2,2) - d2trm*d1c(2)*d1c(2)
    derv2(jzf,iy) = derv2(jzf,iy) - d1trm*d2c(3,2) - d2trm*d1c(3)*d1c(2)
    derv2(jxf,iz) = derv2(jxf,iz) - d1trm*d2c(1,3) - d2trm*d1c(1)*d1c(3)
    derv2(jyf,iz) = derv2(jyf,iz) - d1trm*d2c(2,3) - d2trm*d1c(2)*d1c(3)
    derv2(jzf,iz) = derv2(jzf,iz) - d1trm*d2c(3,3) - d2trm*d1c(3)*d1c(3)
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) - d1trm*d2c(1,1) - d2trm*d1c(1)*d1c(1)
    derv2(iyf,jx) = derv2(iyf,jx) - d1trm*d2c(2,1) - d2trm*d1c(2)*d1c(1)
    derv2(izf,jx) = derv2(izf,jx) - d1trm*d2c(3,1) - d2trm*d1c(3)*d1c(1)
    derv2(ixf,jy) = derv2(ixf,jy) - d1trm*d2c(1,2) - d2trm*d1c(1)*d1c(2)
    derv2(iyf,jy) = derv2(iyf,jy) - d1trm*d2c(2,2) - d2trm*d1c(2)*d1c(2)
    derv2(izf,jy) = derv2(izf,jy) - d1trm*d2c(3,2) - d2trm*d1c(3)*d1c(2)
    derv2(ixf,jz) = derv2(ixf,jz) - d1trm*d2c(1,3) - d2trm*d1c(1)*d1c(3)
    derv2(iyf,jz) = derv2(iyf,jz) - d1trm*d2c(2,3) - d2trm*d1c(2)*d1c(3)
    derv2(izf,jz) = derv2(izf,jz) - d1trm*d2c(3,3) - d2trm*d1c(3)*d1c(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_1bpd')
#endif
!
  return
  end
!
  subroutine add_drv2_2bpd(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                           jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d1trm,d2trm,d1ikc,d1jlc,d2c)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l : Phonon version
!  Version that uses pre-computed block
!  Distributed memory parallel version.
!
!  On entry : 
!
!  lilocal           = is i local to this node?
!  ljlocal           = is j local to this node?
!  lklocal           = is k local to this node?
!  lllocal           = is l local to this node?
!  ix(f)             = position of x coordinate for i in second derivative matrix
!  iy(f)             = position of y coordinate for i in second derivative matrix
!  iz(f)             = position of z coordinate for i in second derivative matrix
!  jx(f)             = position of x coordinate for j in second derivative matrix
!  jy(f)             = position of y coordinate for j in second derivative matrix
!  jz(f)             = position of z coordinate for j in second derivative matrix
!  kx(f)             = position of x coordinate for k in second derivative matrix
!  ky(f)             = position of y coordinate for k in second derivative matrix
!  kz(f)             = position of z coordinate for k in second derivative matrix
!  lx(f)             = position of x coordinate for l in second derivative matrix
!  ly(f)             = position of y coordinate for l in second derivative matrix
!  lz(f)             = position of z coordinate for l in second derivative matrix
!  d1trm             = multiplier for second derivatives
!  d2trm             = multiplier for first derivative products
!  d1ikc             = first derivatives for i-k - Cartesian
!  d1jlc             = first derivatives for j-l - Cartesian
!  d2c               = second derivatives - Cartesian/Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_2bp
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1ikc(3)
  real(dp),    intent(in)                          :: d1jlc(3)
  real(dp),    intent(in)                          :: d2c(3,3)
#ifdef TRACE
  call trace_in('add_drv2_2bpd')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(jyf,ix) = derv2(jyf,ix) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(jzf,ix) = derv2(jzf,ix) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(jxf,iy) = derv2(jxf,iy) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(jyf,iy) = derv2(jyf,iy) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(jzf,iy) = derv2(jzf,iy) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(jxf,iz) = derv2(jxf,iz) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(jyf,iz) = derv2(jyf,iz) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(jzf,iz) = derv2(jzf,iz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(iyf,jx) = derv2(iyf,jx) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(izf,jx) = derv2(izf,jx) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ixf,jy) = derv2(ixf,jy) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(iyf,jy) = derv2(iyf,jy) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(izf,jy) = derv2(izf,jy) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(ixf,jz) = derv2(ixf,jz) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(iyf,jz) = derv2(iyf,jz) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(izf,jz) = derv2(izf,jz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  I-L
!
  if (lilocal) then
    derv2(lxf,ix) = derv2(lxf,ix) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(lyf,ix) = derv2(lyf,ix) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(lzf,ix) = derv2(lzf,ix) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(lxf,iy) = derv2(lxf,iy) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(lyf,iy) = derv2(lyf,iy) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(lzf,iy) = derv2(lzf,iy) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(lxf,iz) = derv2(lxf,iz) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(lyf,iz) = derv2(lyf,iz) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(lzf,iz) = derv2(lzf,iz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
  if (lllocal) then
    derv2(ixf,lx) = derv2(ixf,lx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(iyf,lx) = derv2(iyf,lx) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(izf,lx) = derv2(izf,lx) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ixf,ly) = derv2(ixf,ly) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(iyf,ly) = derv2(iyf,ly) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(izf,ly) = derv2(izf,ly) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(ixf,lz) = derv2(ixf,lz) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(iyf,lz) = derv2(iyf,lz) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(izf,lz) = derv2(izf,lz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  J-K
!
  if (ljlocal) then
    derv2(kxf,jx) = derv2(kxf,jx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(kyf,jx) = derv2(kyf,jx) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(kzf,jx) = derv2(kzf,jx) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(kxf,jy) = derv2(kxf,jy) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(kyf,jy) = derv2(kyf,jy) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(kzf,jy) = derv2(kzf,jy) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(kxf,jz) = derv2(kxf,jz) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(kyf,jz) = derv2(kyf,jz) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(kzf,jz) = derv2(kzf,jz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
  if (lklocal) then
    derv2(jxf,kx) = derv2(jxf,kx) - d1trm*d2c(1,1) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(jyf,kx) = derv2(jyf,kx) - d1trm*d2c(2,1) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(jzf,kx) = derv2(jzf,kx) - d1trm*d2c(3,1) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(jxf,ky) = derv2(jxf,ky) - d1trm*d2c(1,2) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(jyf,ky) = derv2(jyf,ky) - d1trm*d2c(2,2) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(jzf,ky) = derv2(jzf,ky) - d1trm*d2c(3,2) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(jxf,kz) = derv2(jxf,kz) - d1trm*d2c(1,3) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(jyf,kz) = derv2(jyf,kz) - d1trm*d2c(2,3) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(jzf,kz) = derv2(jzf,kz) - d1trm*d2c(3,3) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  K-L
!     
  if (lklocal) then
    derv2(lxf,kx) = derv2(lxf,kx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(lyf,kx) = derv2(lyf,kx) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(lzf,kx) = derv2(lzf,kx) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(lxf,ky) = derv2(lxf,ky) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(lyf,ky) = derv2(lyf,ky) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(lzf,ky) = derv2(lzf,ky) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(lxf,kz) = derv2(lxf,kz) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(lyf,kz) = derv2(lyf,kz) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(lzf,kz) = derv2(lzf,kz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
  if (lllocal) then
    derv2(kxf,lx) = derv2(kxf,lx) + d1trm*d2c(1,1) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(kyf,lx) = derv2(kyf,lx) + d1trm*d2c(1,2) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(kzf,lx) = derv2(kzf,lx) + d1trm*d2c(1,3) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(kxf,ly) = derv2(kxf,ly) + d1trm*d2c(2,1) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(kyf,ly) = derv2(kyf,ly) + d1trm*d2c(2,2) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(kzf,ly) = derv2(kzf,ly) + d1trm*d2c(2,3) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(kxf,lz) = derv2(kxf,lz) + d1trm*d2c(3,1) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(kyf,lz) = derv2(kyf,lz) + d1trm*d2c(3,2) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(kzf,lz) = derv2(kzf,lz) + d1trm*d2c(3,3) + d2trm*d1jlc(3)*d1ikc(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2bpd')
#endif
!
  return
  end
!
  subroutine add_drv2_1bxpd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,d2trm,d1c1,d1c2)
!
!  Adds the second derivatives for a single distance based only on products
!  of first derivatives. Phonon version.
!  Distributed memory parallel version.
!
!  On entry : 
!
!  lilocal           = is i local to this node?
!  ljlocal           = is j local to this node?
!  ix(f)             = position of x coordinate for i in second derivative matrix
!  iy(f)             = position of y coordinate for i in second derivative matrix
!  iz(f)             = position of z coordinate for i in second derivative matrix
!  jx(f)             = position of x coordinate for j in second derivative matrix
!  jy(f)             = position of y coordinate for j in second derivative matrix
!  jz(f)             = position of z coordinate for j in second derivative matrix
!  d2trm             = multiplier for first derivative products
!  d1c1              = first derivatives - Cartesian for first term
!  d1c2              = first derivatives - Cartesian for second term
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_1bxp
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1c1(3)
  real(dp),    intent(in)                          :: d1c2(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_1bxpd')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) - d2trm*d1c1(1)*d1c2(1) - d2trm*d1c2(1)*d1c1(1)
    derv2(jyf,ix) = derv2(jyf,ix) - d2trm*d1c1(2)*d1c2(1) - d2trm*d1c2(2)*d1c1(1)
    derv2(jzf,ix) = derv2(jzf,ix) - d2trm*d1c1(3)*d1c2(1) - d2trm*d1c2(3)*d1c1(1)
    derv2(jxf,iy) = derv2(jxf,iy) - d2trm*d1c1(1)*d1c2(2) - d2trm*d1c2(1)*d1c1(2)
    derv2(jyf,iy) = derv2(jyf,iy) - d2trm*d1c1(2)*d1c2(2) - d2trm*d1c2(2)*d1c1(2)
    derv2(jzf,iy) = derv2(jzf,iy) - d2trm*d1c1(3)*d1c2(2) - d2trm*d1c2(3)*d1c1(2)
    derv2(jxf,iz) = derv2(jxf,iz) - d2trm*d1c1(1)*d1c2(3) - d2trm*d1c2(1)*d1c1(3)
    derv2(jyf,iz) = derv2(jyf,iz) - d2trm*d1c1(2)*d1c2(3) - d2trm*d1c2(2)*d1c1(3)
    derv2(jzf,iz) = derv2(jzf,iz) - d2trm*d1c1(3)*d1c2(3) - d2trm*d1c2(3)*d1c1(3)
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) - d2trm*d1c1(1)*d1c2(1) - d2trm*d1c2(1)*d1c1(1)
    derv2(iyf,jx) = derv2(iyf,jx) - d2trm*d1c1(2)*d1c2(1) - d2trm*d1c2(2)*d1c1(1)
    derv2(izf,jx) = derv2(izf,jx) - d2trm*d1c1(3)*d1c2(1) - d2trm*d1c2(3)*d1c1(1)
    derv2(ixf,jy) = derv2(ixf,jy) - d2trm*d1c1(1)*d1c2(2) - d2trm*d1c2(1)*d1c1(2)
    derv2(iyf,jy) = derv2(iyf,jy) - d2trm*d1c1(2)*d1c2(2) - d2trm*d1c2(2)*d1c1(2)
    derv2(izf,jy) = derv2(izf,jy) - d2trm*d1c1(3)*d1c2(2) - d2trm*d1c2(3)*d1c1(2)
    derv2(ixf,jz) = derv2(ixf,jz) - d2trm*d1c1(1)*d1c2(3) - d2trm*d1c2(1)*d1c1(3)
    derv2(iyf,jz) = derv2(iyf,jz) - d2trm*d1c1(2)*d1c2(3) - d2trm*d1c2(2)*d1c1(3)
    derv2(izf,jz) = derv2(izf,jz) - d2trm*d1c1(3)*d1c2(3) - d2trm*d1c2(3)*d1c1(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_1bxpd')
#endif
!
  return
  end
!
  subroutine add_drv2_2bxpd(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d2trm,d1ikc,d1jlc)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l: Phonon version
!  Version that only adds the product of first derivatives
!  Distributed memory parallel version.
!
!  On entry : 
!
!  lilocal           = is i local to this node?
!  ljlocal           = is j local to this node?
!  lklocal           = is k local to this node?
!  lllocal           = is l local to this node?
!  ix(f)             = position of x coordinate for i in second derivative matrix
!  iy(f)             = position of y coordinate for i in second derivative matrix
!  iz(f)             = position of z coordinate for i in second derivative matrix
!  jx(f)             = position of x coordinate for j in second derivative matrix
!  jy(f)             = position of y coordinate for j in second derivative matrix
!  jz(f)             = position of z coordinate for j in second derivative matrix
!  kx(f)             = position of x coordinate for k in second derivative matrix
!  ky(f)             = position of y coordinate for k in second derivative matrix
!  kz(f)             = position of z coordinate for k in second derivative matrix
!  lx(f)             = position of x coordinate for l in second derivative matrix
!  ly(f)             = position of y coordinate for l in second derivative matrix
!  lz(f)             = position of z coordinate for l in second derivative matrix
!  d2trm             = multiplier for first derivative products
!  d1ikc             = first derivatives for i-k - Cartesian
!  d1jlc             = first derivatives for j-l - Cartesian
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_2bxp
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: d1ikc(3)
  real(dp),    intent(in)                          :: d1jlc(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2bxpd')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(jyf,ix) = derv2(jyf,ix) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(jzf,ix) = derv2(jzf,ix) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(jxf,iy) = derv2(jxf,iy) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(jyf,iy) = derv2(jyf,iy) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(jzf,iy) = derv2(jzf,iy) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(jxf,iz) = derv2(jxf,iz) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(jyf,iz) = derv2(jyf,iz) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(jzf,iz) = derv2(jzf,iz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(iyf,jx) = derv2(iyf,jx) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(izf,jx) = derv2(izf,jx) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(ixf,jy) = derv2(ixf,jy) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(iyf,jy) = derv2(iyf,jy) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(izf,jy) = derv2(izf,jy) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(ixf,jz) = derv2(ixf,jz) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(iyf,jz) = derv2(iyf,jz) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(izf,jz) = derv2(izf,jz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  I-L
!
  if (lilocal) then
    derv2(lxf,ix) = derv2(lxf,ix) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(lyf,ix) = derv2(lyf,ix) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(lzf,ix) = derv2(lzf,ix) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(lxf,iy) = derv2(lxf,iy) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(lyf,iy) = derv2(lyf,iy) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(lzf,iy) = derv2(lzf,iy) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(lxf,iz) = derv2(lxf,iz) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(lyf,iz) = derv2(lyf,iz) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(lzf,iz) = derv2(lzf,iz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
  if (lllocal) then
    derv2(ixf,lx) = derv2(ixf,lx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(iyf,lx) = derv2(iyf,lx) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(izf,lx) = derv2(izf,lx) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(ixf,ly) = derv2(ixf,ly) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(iyf,ly) = derv2(iyf,ly) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(izf,ly) = derv2(izf,ly) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(ixf,lz) = derv2(ixf,lz) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(iyf,lz) = derv2(iyf,lz) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(izf,lz) = derv2(izf,lz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  J-K
!
  if (ljlocal) then
    derv2(kxf,jx) = derv2(kxf,jx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(kyf,jx) = derv2(kyf,jx) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(kzf,jx) = derv2(kzf,jx) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(kxf,jy) = derv2(kxf,jy) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(kyf,jy) = derv2(kyf,jy) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(kzf,jy) = derv2(kzf,jy) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(kxf,jz) = derv2(kxf,jz) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(kyf,jz) = derv2(kyf,jz) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(kzf,jz) = derv2(kzf,jz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
  if (lklocal) then
    derv2(jxf,kx) = derv2(jxf,kx) - d2trm*d1jlc(1)*d1ikc(1)
    derv2(jyf,kx) = derv2(jyf,kx) - d2trm*d1jlc(2)*d1ikc(1)
    derv2(jzf,kx) = derv2(jzf,kx) - d2trm*d1jlc(3)*d1ikc(1)
    derv2(jxf,ky) = derv2(jxf,ky) - d2trm*d1jlc(1)*d1ikc(2)
    derv2(jyf,ky) = derv2(jyf,ky) - d2trm*d1jlc(2)*d1ikc(2)
    derv2(jzf,ky) = derv2(jzf,ky) - d2trm*d1jlc(3)*d1ikc(2)
    derv2(jxf,kz) = derv2(jxf,kz) - d2trm*d1jlc(1)*d1ikc(3)
    derv2(jyf,kz) = derv2(jyf,kz) - d2trm*d1jlc(2)*d1ikc(3)
    derv2(jzf,kz) = derv2(jzf,kz) - d2trm*d1jlc(3)*d1ikc(3)
  endif
!
!  K-L
!     
  if (lklocal) then
    derv2(lxf,kx) = derv2(lxf,kx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(lyf,kx) = derv2(lyf,kx) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(lzf,kx) = derv2(lzf,kx) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(lxf,ky) = derv2(lxf,ky) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(lyf,ky) = derv2(lyf,ky) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(lzf,ky) = derv2(lzf,ky) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(lxf,kz) = derv2(lxf,kz) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(lyf,kz) = derv2(lyf,kz) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(lzf,kz) = derv2(lzf,kz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
!
  if (lllocal) then
    derv2(kxf,lx) = derv2(kxf,lx) + d2trm*d1jlc(1)*d1ikc(1)
    derv2(kyf,lx) = derv2(kyf,lx) + d2trm*d1jlc(1)*d1ikc(2)
    derv2(kzf,lx) = derv2(kzf,lx) + d2trm*d1jlc(1)*d1ikc(3)
    derv2(kxf,ly) = derv2(kxf,ly) + d2trm*d1jlc(2)*d1ikc(1)
    derv2(kyf,ly) = derv2(kyf,ly) + d2trm*d1jlc(2)*d1ikc(2)
    derv2(kzf,ly) = derv2(kzf,ly) + d2trm*d1jlc(2)*d1ikc(3)
    derv2(kxf,lz) = derv2(kxf,lz) + d2trm*d1jlc(3)*d1ikc(1)
    derv2(kyf,lz) = derv2(kyf,lz) + d2trm*d1jlc(3)*d1ikc(2)
    derv2(kzf,lz) = derv2(kzf,lz) + d2trm*d1jlc(3)*d1ikc(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2bxpd')
#endif
!
  return
  end
!
  subroutine add_drv2_1gp(ix,iy,iz,jx,jy,jz,xij,yij,zij,d1trm,d2trm,d2r2ijdx2)
!
!  Adds the second derivatives for a single distance: Phonon version
!  Generalised version that replaces assumption of delta function for 1/2 (d2r/dadb)
!  and uses products of x/y/z
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d1trm             = first derivative of energy w.r.t. rij divided by rij
!  d2trm             = derivative of d1trm w.r.t. rij divided by rij
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/21 Created from add_drv2_1g
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: d2r2ijdx2(6)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_1gp')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) - d1trm*d2r2ijdx2(1) - d2trm*xij*xij
    derv2(jy,ix) = derv2(jy,ix) - d1trm*d2r2ijdx2(6) - d2trm*yij*xij
    derv2(jz,ix) = derv2(jz,ix) - d1trm*d2r2ijdx2(5) - d2trm*zij*xij
    derv2(jx,iy) = derv2(jx,iy) - d1trm*d2r2ijdx2(6) - d2trm*xij*yij
    derv2(jy,iy) = derv2(jy,iy) - d1trm*d2r2ijdx2(2) - d2trm*yij*yij
    derv2(jz,iy) = derv2(jz,iy) - d1trm*d2r2ijdx2(4) - d2trm*zij*yij
    derv2(jx,iz) = derv2(jx,iz) - d1trm*d2r2ijdx2(5) - d2trm*xij*zij
    derv2(jy,iz) = derv2(jy,iz) - d1trm*d2r2ijdx2(4) - d2trm*yij*zij
    derv2(jz,iz) = derv2(jz,iz) - d1trm*d2r2ijdx2(3) - d2trm*zij*zij
  else
    derv2(ix,jx) = derv2(ix,jx) - d1trm*d2r2ijdx2(1) - d2trm*xij*xij
    derv2(iy,jx) = derv2(iy,jx) - d1trm*d2r2ijdx2(6) - d2trm*yij*xij
    derv2(iz,jx) = derv2(iz,jx) - d1trm*d2r2ijdx2(5) - d2trm*zij*xij
    derv2(ix,jy) = derv2(ix,jy) - d1trm*d2r2ijdx2(6) - d2trm*xij*yij
    derv2(iy,jy) = derv2(iy,jy) - d1trm*d2r2ijdx2(2) - d2trm*yij*yij
    derv2(iz,jy) = derv2(iz,jy) - d1trm*d2r2ijdx2(4) - d2trm*zij*yij
    derv2(ix,jz) = derv2(ix,jz) - d1trm*d2r2ijdx2(5) - d2trm*xij*zij
    derv2(iy,jz) = derv2(iy,jz) - d1trm*d2r2ijdx2(4) - d2trm*yij*zij
    derv2(iz,jz) = derv2(iz,jz) - d1trm*d2r2ijdx2(3) - d2trm*zij*zij
  endif
#ifdef TRACE
  call trace_out('add_drv2_1gp')
#endif
!
  return
  end
!
  subroutine add_drv2_2sgp(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,d1trm,dr2ik,dr2jl,d2trm,d2r2ikjl)
!
!  Adds the second derivatives for a single distance that depends on two components, i-k and j-l
!  Generalised version where matrices with the terms are passed in. Phonon version
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d1trm             = (1/r dE/dr)
!  d2trm             = (1/r d/dr (1/r dE/dr))
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/21 Created from add_drv2_2sg
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: dr2ik(3)
  real(dp),    intent(in)                          :: d2r2ikjl(3,3)
  real(dp),    intent(in)                          :: dr2jl(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2sgp')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(jy,ix) = derv2(jy,ix) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(jz,ix) = derv2(jz,ix) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(jx,iy) = derv2(jx,iy) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(jy,iy) = derv2(jy,iy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(jz,iy) = derv2(jz,iy) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(jx,iz) = derv2(jx,iz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(jy,iz) = derv2(jy,iz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(jz,iz) = derv2(jz,iz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(ix,jx) = derv2(ix,jx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(iy,jx) = derv2(iy,jx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(iz,jx) = derv2(iz,jx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ix,jy) = derv2(ix,jy) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(iy,jy) = derv2(iy,jy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(iz,jy) = derv2(iz,jy) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(ix,jz) = derv2(ix,jz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(iy,jz) = derv2(iy,jz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(iz,jz) = derv2(iz,jz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  I-L
!
  if (lx.ge.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(ly,ix) = derv2(ly,ix) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(lz,ix) = derv2(lz,ix) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(lx,iy) = derv2(lx,iy) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(ly,iy) = derv2(ly,iy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(lz,iy) = derv2(lz,iy) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(lx,iz) = derv2(lx,iz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ly,iz) = derv2(ly,iz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(lz,iz) = derv2(lz,iz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(ix,lx) = derv2(ix,lx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(iy,lx) = derv2(iy,lx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(iz,lx) = derv2(iz,lx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ix,ly) = derv2(ix,ly) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(iy,ly) = derv2(iy,ly) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(iz,ly) = derv2(iz,ly) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(ix,lz) = derv2(ix,lz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(iy,lz) = derv2(iy,lz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(iz,lz) = derv2(iz,lz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  J-K
!
  if (kx.ge.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(ky,jx) = derv2(ky,jx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(kz,jx) = derv2(kz,jx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(kx,jy) = derv2(kx,jy) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(ky,jy) = derv2(ky,jy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(kz,jy) = derv2(kz,jy) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(kx,jz) = derv2(kx,jz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(ky,jz) = derv2(ky,jz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(kz,jz) = derv2(kz,jz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(jx,kx) = derv2(jx,kx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(jy,kx) = derv2(jy,kx) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(jz,kx) = derv2(jz,kx) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(jx,ky) = derv2(jx,ky) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(jy,ky) = derv2(jy,ky) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(jz,ky) = derv2(jz,ky) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(jx,kz) = derv2(jx,kz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(jy,kz) = derv2(jy,kz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(jz,kz) = derv2(jz,kz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  K-L
!     
  if (lx.ge.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(ly,kx) = derv2(ly,kx) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(lz,kx) = derv2(lz,kx) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(lx,ky) = derv2(lx,ky) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(ly,ky) = derv2(ly,ky) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(lz,ky) = derv2(lz,ky) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(lx,kz) = derv2(lx,kz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ly,kz) = derv2(ly,kz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(lz,kz) = derv2(lz,kz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(kx,lx) = derv2(kx,lx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(ky,lx) = derv2(ky,lx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(kz,lx) = derv2(kz,lx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(kx,ly) = derv2(kx,ly) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(ky,ly) = derv2(ky,ly) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(kz,ly) = derv2(kz,ly) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(kx,lz) = derv2(kx,lz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(ky,lz) = derv2(ky,lz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(kz,lz) = derv2(kz,lz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2sgp')
#endif
!
  return
  end
!
  subroutine add_drv2_2sgnosp(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,d1trm,dr2ik,dr2jl,d2trm,d2r2ikjl)
!
!  Adds the second derivatives for a single distance that depends on two components, i-k and j-l
!  Generalised version where matrices with the terms are passed in. Phonon version
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d1trm             = (1/r dE/dr)
!  d2trm             = (1/r d/dr (1/r dE/dr))
!
!  On exit :
!
!  Second derivatives have been added to
!
!   3/21 Created from add_drv2_2sgnos
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: dr2ik(3)
  real(dp),    intent(in)                          :: d2r2ikjl(3,3)
  real(dp),    intent(in)                          :: dr2jl(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2sgnosp')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.ge.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(jy,ix) = derv2(jy,ix) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(jz,ix) = derv2(jz,ix) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(jx,iy) = derv2(jx,iy) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(jy,iy) = derv2(jy,iy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(jz,iy) = derv2(jz,iy) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(jx,iz) = derv2(jx,iz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(jy,iz) = derv2(jy,iz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(jz,iz) = derv2(jz,iz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(ix,jx) = derv2(ix,jx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(iy,jx) = derv2(iy,jx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(iz,jx) = derv2(iz,jx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ix,jy) = derv2(ix,jy) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(iy,jy) = derv2(iy,jy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(iz,jy) = derv2(iz,jy) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(ix,jz) = derv2(ix,jz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(iy,jz) = derv2(iy,jz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(iz,jz) = derv2(iz,jz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  I-L
!
  if (lx.ge.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(ly,ix) = derv2(ly,ix) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(lz,ix) = derv2(lz,ix) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(lx,iy) = derv2(lx,iy) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(ly,iy) = derv2(ly,iy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(lz,iy) = derv2(lz,iy) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(lx,iz) = derv2(lx,iz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ly,iz) = derv2(ly,iz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(lz,iz) = derv2(lz,iz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(ix,lx) = derv2(ix,lx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(iy,lx) = derv2(iy,lx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(iz,lx) = derv2(iz,lx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ix,ly) = derv2(ix,ly) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(iy,ly) = derv2(iy,ly) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(iz,ly) = derv2(iz,ly) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(ix,lz) = derv2(ix,lz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(iy,lz) = derv2(iy,lz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(iz,lz) = derv2(iz,lz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  J-K
!
  if (kx.ge.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(ky,jx) = derv2(ky,jx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(kz,jx) = derv2(kz,jx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(kx,jy) = derv2(kx,jy) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(ky,jy) = derv2(ky,jy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(kz,jy) = derv2(kz,jy) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(kx,jz) = derv2(kx,jz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(ky,jz) = derv2(ky,jz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(kz,jz) = derv2(kz,jz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(jx,kx) = derv2(jx,kx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(jy,kx) = derv2(jy,kx) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(jz,kx) = derv2(jz,kx) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(jx,ky) = derv2(jx,ky) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(jy,ky) = derv2(jy,ky) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(jz,ky) = derv2(jz,ky) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(jx,kz) = derv2(jx,kz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(jy,kz) = derv2(jy,kz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(jz,kz) = derv2(jz,kz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  K-L
!     
  if (lx.ge.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(ly,kx) = derv2(ly,kx) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(lz,kx) = derv2(lz,kx) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(lx,ky) = derv2(lx,ky) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(ly,ky) = derv2(ly,ky) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(lz,ky) = derv2(lz,ky) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(lx,kz) = derv2(lx,kz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ly,kz) = derv2(ly,kz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(lz,kz) = derv2(lz,kz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  else
    derv2(kx,lx) = derv2(kx,lx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(ky,lx) = derv2(ky,lx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(kz,lx) = derv2(kz,lx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(kx,ly) = derv2(kx,ly) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(ky,ly) = derv2(ky,ly) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(kz,ly) = derv2(kz,ly) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(kx,lz) = derv2(kx,lz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(ky,lz) = derv2(ky,lz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(kz,lz) = derv2(kz,lz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2sgnosp')
#endif
!
  return
  end
!
  subroutine add_drv2_1gpd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                           xij,yij,zij,d1trm,d2trm,d2r2ijdx2)
!
!  Adds the second derivatives for a single distance: Phonon version
!  Generalised version that replaces assumption of delta function for 1/2 (d2r/dadb)
!  and uses products of x/y/z
!  Distrubuted memory parallel version.
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d1trm             = first derivative of energy w.r.t. rij divided by rij
!  d2trm             = derivative of d1trm w.r.t. rij divided by rij
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_1gp
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: d2r2ijdx2(6)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_1gpd')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) - d1trm*d2r2ijdx2(1) - d2trm*xij*xij
    derv2(jyf,ix) = derv2(jyf,ix) - d1trm*d2r2ijdx2(6) - d2trm*yij*xij
    derv2(jzf,ix) = derv2(jzf,ix) - d1trm*d2r2ijdx2(5) - d2trm*zij*xij
    derv2(jxf,iy) = derv2(jxf,iy) - d1trm*d2r2ijdx2(6) - d2trm*xij*yij
    derv2(jyf,iy) = derv2(jyf,iy) - d1trm*d2r2ijdx2(2) - d2trm*yij*yij
    derv2(jzf,iy) = derv2(jzf,iy) - d1trm*d2r2ijdx2(4) - d2trm*zij*yij
    derv2(jxf,iz) = derv2(jxf,iz) - d1trm*d2r2ijdx2(5) - d2trm*xij*zij
    derv2(jyf,iz) = derv2(jyf,iz) - d1trm*d2r2ijdx2(4) - d2trm*yij*zij
    derv2(jzf,iz) = derv2(jzf,iz) - d1trm*d2r2ijdx2(3) - d2trm*zij*zij
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) - d1trm*d2r2ijdx2(1) - d2trm*xij*xij
    derv2(iyf,jx) = derv2(iyf,jx) - d1trm*d2r2ijdx2(6) - d2trm*yij*xij
    derv2(izf,jx) = derv2(izf,jx) - d1trm*d2r2ijdx2(5) - d2trm*zij*xij
    derv2(ixf,jy) = derv2(ixf,jy) - d1trm*d2r2ijdx2(6) - d2trm*xij*yij
    derv2(iyf,jy) = derv2(iyf,jy) - d1trm*d2r2ijdx2(2) - d2trm*yij*yij
    derv2(izf,jy) = derv2(izf,jy) - d1trm*d2r2ijdx2(4) - d2trm*zij*yij
    derv2(ixf,jz) = derv2(ixf,jz) - d1trm*d2r2ijdx2(5) - d2trm*xij*zij
    derv2(iyf,jz) = derv2(iyf,jz) - d1trm*d2r2ijdx2(4) - d2trm*yij*zij
    derv2(izf,jz) = derv2(izf,jz) - d1trm*d2r2ijdx2(3) - d2trm*zij*zij
  endif
#ifdef TRACE
  call trace_out('add_drv2_1gpd')
#endif
!
  return
  end
!
  subroutine add_drv2_2sgpd(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d1trm,dr2ik,dr2jl,d2trm,d2r2ikjl)
!
!  Adds the second derivatives for a single distance that depends on two components, i-k and j-l
!  Generalised version where matrices with the terms are passed in. Phonon version
!  Distrubuted memory parallel version.
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d1trm             = (1/r dE/dr)
!  d2trm             = (1/r d/dr (1/r dE/dr))
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_2sgp
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: dr2ik(3)
  real(dp),    intent(in)                          :: d2r2ikjl(3,3)
  real(dp),    intent(in)                          :: dr2jl(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2sgpd')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(jyf,ix) = derv2(jyf,ix) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(jzf,ix) = derv2(jzf,ix) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(jxf,iy) = derv2(jxf,iy) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(jyf,iy) = derv2(jyf,iy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(jzf,iy) = derv2(jzf,iy) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(jxf,iz) = derv2(jxf,iz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(jyf,iz) = derv2(jyf,iz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(jzf,iz) = derv2(jzf,iz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(iyf,jx) = derv2(iyf,jx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(izf,jx) = derv2(izf,jx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ixf,jy) = derv2(ixf,jy) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(iyf,jy) = derv2(iyf,jy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(izf,jy) = derv2(izf,jy) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(ixf,jz) = derv2(ixf,jz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(iyf,jz) = derv2(iyf,jz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(izf,jz) = derv2(izf,jz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  I-L
!
  if (lilocal) then
    derv2(lxf,ix) = derv2(lxf,ix) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(lyf,ix) = derv2(lyf,ix) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(lzf,ix) = derv2(lzf,ix) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(lxf,iy) = derv2(lxf,iy) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(lyf,iy) = derv2(lyf,iy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(lzf,iy) = derv2(lzf,iy) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(lxf,iz) = derv2(lxf,iz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(lyf,iz) = derv2(lyf,iz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(lzf,iz) = derv2(lzf,iz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lllocal) then
    derv2(ixf,lx) = derv2(ixf,lx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(iyf,lx) = derv2(iyf,lx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(izf,lx) = derv2(izf,lx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ixf,ly) = derv2(ixf,ly) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(iyf,ly) = derv2(iyf,ly) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(izf,ly) = derv2(izf,ly) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(ixf,lz) = derv2(ixf,lz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(iyf,lz) = derv2(iyf,lz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(izf,lz) = derv2(izf,lz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  J-K
!
  if (ljlocal) then
    derv2(kxf,jx) = derv2(kxf,jx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(kyf,jx) = derv2(kyf,jx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(kzf,jx) = derv2(kzf,jx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(kxf,jy) = derv2(kxf,jy) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(kyf,jy) = derv2(kyf,jy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(kzf,jy) = derv2(kzf,jy) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(kxf,jz) = derv2(kxf,jz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(kyf,jz) = derv2(kyf,jz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(kzf,jz) = derv2(kzf,jz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lklocal) then
    derv2(jxf,kx) = derv2(jxf,kx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(jyf,kx) = derv2(jyf,kx) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(jzf,kx) = derv2(jzf,kx) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(jxf,ky) = derv2(jxf,ky) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(jyf,ky) = derv2(jyf,ky) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(jzf,ky) = derv2(jzf,ky) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(jxf,kz) = derv2(jxf,kz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(jyf,kz) = derv2(jyf,kz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(jzf,kz) = derv2(jzf,kz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  K-L
!     
  if (lklocal) then
    derv2(lxf,kx) = derv2(lxf,kx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(lyf,kx) = derv2(lyf,kx) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(lzf,kx) = derv2(lzf,kx) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(lxf,ky) = derv2(lxf,ky) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(lyf,ky) = derv2(lyf,ky) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(lzf,ky) = derv2(lzf,ky) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(lxf,kz) = derv2(lxf,kz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(lyf,kz) = derv2(lyf,kz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(lzf,kz) = derv2(lzf,kz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lllocal) then
    derv2(kxf,lx) = derv2(kxf,lx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(kyf,lx) = derv2(kyf,lx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(kzf,lx) = derv2(kzf,lx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(kxf,ly) = derv2(kxf,ly) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(kyf,ly) = derv2(kyf,ly) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(kzf,ly) = derv2(kzf,ly) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(kxf,lz) = derv2(kxf,lz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(kyf,lz) = derv2(kyf,lz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(kzf,lz) = derv2(kzf,lz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2sgpd')
#endif
!
  return
  end
!
  subroutine add_drv2_2sgnospd(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                               jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,d1trm,dr2ik,dr2jl,d2trm,d2r2ikjl)
!
!  Adds the second derivatives for a single distance that depends on two components, i-k and j-l
!  Generalised version where matrices with the terms are passed in. Phonon version
!  Distrubuted memory parallel version.
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  d1trm             = (1/r dE/dr)
!  d2trm             = (1/r d/dr (1/r dE/dr))
!
!  On exit :
!
!  Second derivatives have been added to
!
!   4/22 Created from add_drv2_2sgnosp
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
  use derivatives
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: kxf
  integer(i4), intent(in)                          :: kyf
  integer(i4), intent(in)                          :: kzf
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  integer(i4), intent(in)                          :: lxf
  integer(i4), intent(in)                          :: lyf
  integer(i4), intent(in)                          :: lzf
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  logical,     intent(in)                          :: lklocal
  logical,     intent(in)                          :: lllocal
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: dr2ik(3)
  real(dp),    intent(in)                          :: d2r2ikjl(3,3)
  real(dp),    intent(in)                          :: dr2jl(3)
!
!  Local variables
!
#ifdef TRACE
  call trace_in('add_drv2_2sgnospd')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (lilocal) then
    derv2(jxf,ix) = derv2(jxf,ix) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(jyf,ix) = derv2(jyf,ix) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(jzf,ix) = derv2(jzf,ix) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(jxf,iy) = derv2(jxf,iy) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(jyf,iy) = derv2(jyf,iy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(jzf,iy) = derv2(jzf,iy) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(jxf,iz) = derv2(jxf,iz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(jyf,iz) = derv2(jyf,iz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(jzf,iz) = derv2(jzf,iz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (ljlocal) then
    derv2(ixf,jx) = derv2(ixf,jx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(iyf,jx) = derv2(iyf,jx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(izf,jx) = derv2(izf,jx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(ixf,jy) = derv2(ixf,jy) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(iyf,jy) = derv2(iyf,jy) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(izf,jy) = derv2(izf,jy) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(ixf,jz) = derv2(ixf,jz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(iyf,jz) = derv2(iyf,jz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(izf,jz) = derv2(izf,jz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  I-L
!
  if (lilocal) then
    derv2(lxf,ix) = derv2(lxf,ix) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(lyf,ix) = derv2(lyf,ix) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(lzf,ix) = derv2(lzf,ix) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(lxf,iy) = derv2(lxf,iy) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(lyf,iy) = derv2(lyf,iy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(lzf,iy) = derv2(lzf,iy) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(lxf,iz) = derv2(lxf,iz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(lyf,iz) = derv2(lyf,iz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(lzf,iz) = derv2(lzf,iz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lllocal) then
    derv2(ixf,lx) = derv2(ixf,lx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(iyf,lx) = derv2(iyf,lx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(izf,lx) = derv2(izf,lx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(ixf,ly) = derv2(ixf,ly) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(iyf,ly) = derv2(iyf,ly) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(izf,ly) = derv2(izf,ly) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(ixf,lz) = derv2(ixf,lz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(iyf,lz) = derv2(iyf,lz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(izf,lz) = derv2(izf,lz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  J-K
!
  if (ljlocal) then
    derv2(kxf,jx) = derv2(kxf,jx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(kyf,jx) = derv2(kyf,jx) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(kzf,jx) = derv2(kzf,jx) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(kxf,jy) = derv2(kxf,jy) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(kyf,jy) = derv2(kyf,jy) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(kzf,jy) = derv2(kzf,jy) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(kxf,jz) = derv2(kxf,jz) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(kyf,jz) = derv2(kyf,jz) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(kzf,jz) = derv2(kzf,jz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lklocal) then
    derv2(jxf,kx) = derv2(jxf,kx) - d1trm*d2r2ikjl(1,1) - d2trm*dr2ik(1)*dr2jl(1)
    derv2(jyf,kx) = derv2(jyf,kx) - d1trm*d2r2ikjl(2,1) - d2trm*dr2ik(1)*dr2jl(2)
    derv2(jzf,kx) = derv2(jzf,kx) - d1trm*d2r2ikjl(3,1) - d2trm*dr2ik(1)*dr2jl(3)
    derv2(jxf,ky) = derv2(jxf,ky) - d1trm*d2r2ikjl(1,2) - d2trm*dr2ik(2)*dr2jl(1)
    derv2(jyf,ky) = derv2(jyf,ky) - d1trm*d2r2ikjl(2,2) - d2trm*dr2ik(2)*dr2jl(2)
    derv2(jzf,ky) = derv2(jzf,ky) - d1trm*d2r2ikjl(3,2) - d2trm*dr2ik(2)*dr2jl(3)
    derv2(jxf,kz) = derv2(jxf,kz) - d1trm*d2r2ikjl(1,3) - d2trm*dr2ik(3)*dr2jl(1)
    derv2(jyf,kz) = derv2(jyf,kz) - d1trm*d2r2ikjl(2,3) - d2trm*dr2ik(3)*dr2jl(2)
    derv2(jzf,kz) = derv2(jzf,kz) - d1trm*d2r2ikjl(3,3) - d2trm*dr2ik(3)*dr2jl(3)
  endif
!
!  K-L
!     
  if (lklocal) then
    derv2(lxf,kx) = derv2(lxf,kx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(lyf,kx) = derv2(lyf,kx) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(lzf,kx) = derv2(lzf,kx) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(lxf,ky) = derv2(lxf,ky) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(lyf,ky) = derv2(lyf,ky) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(lzf,ky) = derv2(lzf,ky) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(lxf,kz) = derv2(lxf,kz) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(lyf,kz) = derv2(lyf,kz) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(lzf,kz) = derv2(lzf,kz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
!
  if (lllocal) then
    derv2(kxf,lx) = derv2(kxf,lx) + d1trm*d2r2ikjl(1,1) + d2trm*dr2ik(1)*dr2jl(1)
    derv2(kyf,lx) = derv2(kyf,lx) + d1trm*d2r2ikjl(1,2) + d2trm*dr2ik(2)*dr2jl(1)
    derv2(kzf,lx) = derv2(kzf,lx) + d1trm*d2r2ikjl(1,3) + d2trm*dr2ik(3)*dr2jl(1)
    derv2(kxf,ly) = derv2(kxf,ly) + d1trm*d2r2ikjl(2,1) + d2trm*dr2ik(1)*dr2jl(2)
    derv2(kyf,ly) = derv2(kyf,ly) + d1trm*d2r2ikjl(2,2) + d2trm*dr2ik(2)*dr2jl(2)
    derv2(kzf,ly) = derv2(kzf,ly) + d1trm*d2r2ikjl(2,3) + d2trm*dr2ik(3)*dr2jl(2)
    derv2(kxf,lz) = derv2(kxf,lz) + d1trm*d2r2ikjl(3,1) + d2trm*dr2ik(1)*dr2jl(3)
    derv2(kyf,lz) = derv2(kyf,lz) + d1trm*d2r2ikjl(3,2) + d2trm*dr2ik(2)*dr2jl(3)
    derv2(kzf,lz) = derv2(kzf,lz) + d1trm*d2r2ikjl(3,3) + d2trm*dr2ik(3)*dr2jl(3)
  endif
#ifdef TRACE
  call trace_out('add_drv2_2sgnospd')
#endif
!
  return
  end
!
  subroutine add_drv2_1ut(ix,iy,iz,jx,jy,jz,xij,yij,zij,d1trm,d2trm,dr2ijds,d2r2ijdx2, &
                          d2r2ijds2,d2r2ijdsdx)
!
!  Adds the second derivatives for a single distance to upper triangle
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  xij               = x component of vector from i to j
!  yij               = y component of vector from i to j
!  zij               = z component of vector from i to j
!  d1trm             = first derivative of energy w.r.t. rij divided by rij
!  d2trm             = derivative of d1trm w.r.t. rij divided by rij
!
!  On exit :
!
!  Second derivatives have been added to
!
!  11/20 Created
!   9/21 Now adds to upper triangular
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
!  Julian Gale, CIC, Curtin University, September 2021
!
  use datatypes
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  real(dp),    intent(in)                          :: d1trm
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: dr2ijds(6)
  real(dp),    intent(in)                          :: d2r2ijdx2(6)
  real(dp),    intent(in)                          :: d2r2ijds2(6,6)
  real(dp),    intent(in)                          :: d2r2ijdsdx(6,3)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_1ut')
#endif
!----------------------------------
!  Second derivatives w.r.t. rij  |
!----------------------------------
  if (ix.ge.jx) then
    derv2(jx,ix) = derv2(jx,ix) - d2trm*d2r2ijdx2(1)
    derv2(jy,ix) = derv2(jy,ix) - d2trm*d2r2ijdx2(6)
    derv2(jz,ix) = derv2(jz,ix) - d2trm*d2r2ijdx2(5)
    derv2(jx,iy) = derv2(jx,iy) - d2trm*d2r2ijdx2(6)
    derv2(jy,iy) = derv2(jy,iy) - d2trm*d2r2ijdx2(2)
    derv2(jz,iy) = derv2(jz,iy) - d2trm*d2r2ijdx2(4)
    derv2(jx,iz) = derv2(jx,iz) - d2trm*d2r2ijdx2(5)
    derv2(jy,iz) = derv2(jy,iz) - d2trm*d2r2ijdx2(4)
    derv2(jz,iz) = derv2(jz,iz) - d2trm*d2r2ijdx2(3)
    derv2(jx,ix) = derv2(jx,ix) - d1trm
    derv2(jy,iy) = derv2(jy,iy) - d1trm
    derv2(jz,iz) = derv2(jz,iz) - d1trm
  else
    derv2(ix,jx) = derv2(ix,jx) - d2trm*d2r2ijdx2(1)
    derv2(iy,jx) = derv2(iy,jx) - d2trm*d2r2ijdx2(6)
    derv2(iz,jx) = derv2(iz,jx) - d2trm*d2r2ijdx2(5)
    derv2(ix,jy) = derv2(ix,jy) - d2trm*d2r2ijdx2(6)
    derv2(iy,jy) = derv2(iy,jy) - d2trm*d2r2ijdx2(2)
    derv2(iz,jy) = derv2(iz,jy) - d2trm*d2r2ijdx2(4)
    derv2(ix,jz) = derv2(ix,jz) - d2trm*d2r2ijdx2(5)
    derv2(iy,jz) = derv2(iy,jz) - d2trm*d2r2ijdx2(4)
    derv2(iz,jz) = derv2(iz,jz) - d2trm*d2r2ijdx2(3)
    derv2(ix,jx) = derv2(ix,jx) - d1trm
    derv2(iy,jy) = derv2(iy,jy) - d1trm
    derv2(iz,jz) = derv2(iz,jz) - d1trm
  endif
!
  if (lstr) then
!
!  Mixed derivatives
!
    if (ix.ne.jx) then
      do kl = 1,nstrains
        ks = nstrptr(kl)
        derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2ijdsdx(ks,1) - xij*d2trm*dr2ijds(ks)
        derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2ijdsdx(ks,2) - yij*d2trm*dr2ijds(ks)
        derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2ijdsdx(ks,3) - zij*d2trm*dr2ijds(ks)
        derv3(jx,kl) = derv3(jx,kl) + d1trm*d2r2ijdsdx(ks,1) + xij*d2trm*dr2ijds(ks)
        derv3(jy,kl) = derv3(jy,kl) + d1trm*d2r2ijdsdx(ks,2) + yij*d2trm*dr2ijds(ks)
        derv3(jz,kl) = derv3(jz,kl) + d1trm*d2r2ijdsdx(ks,3) + zij*d2trm*dr2ijds(ks)
      enddo
    endif
!
!  Strain-strain
!
    do kk = 1,nstrains
      ks = nstrptr(kk)
      do kl = 1,nstrains
        kt = nstrptr(kl)
        sderv2(kl,kk) = sderv2(kl,kk) + d2trm*dr2ijds(kt)*dr2ijds(ks) + d1trm*d2r2ijds2(kt,ks)
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('add_drv2_1ut')
#endif
!
  return
  end
!
  subroutine add_drv2_2ut(ix,iy,iz,kx,ky,kz,jx,jy,jz,lx,ly,lz,xik,yik,zik,xjl,yjl,zjl,d2trm,dr2ikds,dr2jlds)
!
!  Adds the second derivatives for a pair of distances, i-k and j-l to upper triangle
!
!  On entry : 
!
!  ix                = position of x coordinate for i in second derivative matrix
!  iy                = position of y coordinate for i in second derivative matrix
!  iz                = position of z coordinate for i in second derivative matrix
!  jx                = position of x coordinate for j in second derivative matrix
!  jy                = position of y coordinate for j in second derivative matrix
!  jz                = position of z coordinate for j in second derivative matrix
!  kx                = position of x coordinate for k in second derivative matrix
!  ky                = position of y coordinate for k in second derivative matrix
!  kz                = position of z coordinate for k in second derivative matrix
!  lx                = position of x coordinate for l in second derivative matrix
!  ly                = position of y coordinate for l in second derivative matrix
!  lz                = position of z coordinate for l in second derivative matrix
!  xik               = x component of vector from i to k
!  yik               = y component of vector from i to k
!  zik               = z component of vector from i to k
!  xjl               = x component of vector from i to k
!  yjl               = y component of vector from i to k
!  zjl               = z component of vector from i to k
!  d2trm             = second derivative of energy w.r.t. rik and rjl with division by distances
!
!  On exit :
!
!  Second derivatives have been added to
!
!  11/20 Created
!   9/21 Now adds to upper triangular
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
!  Julian Gale, CIC, Curtin University, September 2021
!
  use datatypes
  use current,        only : nstrains, nstrptr
  use derivatives
  use m_strain,       only : real1strterm
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: kx
  integer(i4), intent(in)                          :: ky
  integer(i4), intent(in)                          :: kz
  integer(i4), intent(in)                          :: lx
  integer(i4), intent(in)                          :: ly
  integer(i4), intent(in)                          :: lz
  real(dp),    intent(in)                          :: d2trm
  real(dp),    intent(in)                          :: xik
  real(dp),    intent(in)                          :: yik
  real(dp),    intent(in)                          :: zik
  real(dp),    intent(in)                          :: xjl
  real(dp),    intent(in)                          :: yjl
  real(dp),    intent(in)                          :: zjl
  real(dp),    intent(in)                          :: dr2ikds(6)
  real(dp),    intent(in)                          :: dr2jlds(6)
!
!  Local variables
!
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
#ifdef TRACE
  call trace_in('add_drv2_2ut')
#endif
!-----------------------------------------------------
!  Add second derivatives for 2 different distances  |
!-----------------------------------------------------
!
!  I-J
!
  if (jx.lt.ix) then
    derv2(jx,ix) = derv2(jx,ix) + d2trm*xjl*xik
    derv2(jy,ix) = derv2(jy,ix) + d2trm*yjl*xik
    derv2(jz,ix) = derv2(jz,ix) + d2trm*zjl*xik
    derv2(jx,iy) = derv2(jx,iy) + d2trm*xjl*yik
    derv2(jy,iy) = derv2(jy,iy) + d2trm*yjl*yik
    derv2(jz,iy) = derv2(jz,iy) + d2trm*zjl*yik
    derv2(jx,iz) = derv2(jx,iz) + d2trm*xjl*zik
    derv2(jy,iz) = derv2(jy,iz) + d2trm*yjl*zik
    derv2(jz,iz) = derv2(jz,iz) + d2trm*zjl*zik
  else
    derv2(ix,jx) = derv2(ix,jx) + d2trm*xjl*xik
    derv2(iy,jx) = derv2(iy,jx) + d2trm*xjl*yik
    derv2(iz,jx) = derv2(iz,jx) + d2trm*xjl*zik
    derv2(ix,jy) = derv2(ix,jy) + d2trm*yjl*xik
    derv2(iy,jy) = derv2(iy,jy) + d2trm*yjl*yik
    derv2(iz,jy) = derv2(iz,jy) + d2trm*yjl*zik
    derv2(ix,jz) = derv2(ix,jz) + d2trm*zjl*xik
    derv2(iy,jz) = derv2(iy,jz) + d2trm*zjl*yik
    derv2(iz,jz) = derv2(iz,jz) + d2trm*zjl*zik
  endif
!
!  I-L
!
  if (lx.lt.ix) then
    derv2(lx,ix) = derv2(lx,ix) - d2trm*xjl*xik
    derv2(ly,ix) = derv2(ly,ix) - d2trm*yjl*xik
    derv2(lz,ix) = derv2(lz,ix) - d2trm*zjl*xik
    derv2(lx,iy) = derv2(lx,iy) - d2trm*xjl*yik
    derv2(ly,iy) = derv2(ly,iy) - d2trm*yjl*yik
    derv2(lz,iy) = derv2(lz,iy) - d2trm*zjl*yik
    derv2(lx,iz) = derv2(lx,iz) - d2trm*xjl*zik
    derv2(ly,iz) = derv2(ly,iz) - d2trm*yjl*zik
    derv2(lz,iz) = derv2(lz,iz) - d2trm*zjl*zik
  else
    derv2(ix,lx) = derv2(ix,lx) - d2trm*xjl*xik
    derv2(iy,lx) = derv2(iy,lx) - d2trm*xjl*yik
    derv2(iz,lx) = derv2(iz,lx) - d2trm*xjl*zik
    derv2(ix,ly) = derv2(ix,ly) - d2trm*yjl*xik
    derv2(iy,ly) = derv2(iy,ly) - d2trm*yjl*yik
    derv2(iz,ly) = derv2(iz,ly) - d2trm*yjl*zik
    derv2(ix,lz) = derv2(ix,lz) - d2trm*zjl*xik
    derv2(iy,lz) = derv2(iy,lz) - d2trm*zjl*yik
    derv2(iz,lz) = derv2(iz,lz) - d2trm*zjl*zik
  endif
!
!  J-K
!
  if (kx.lt.jx) then
    derv2(kx,jx) = derv2(kx,jx) - d2trm*xjl*xik
    derv2(ky,jx) = derv2(ky,jx) - d2trm*xjl*yik
    derv2(kz,jx) = derv2(kz,jx) - d2trm*xjl*zik
    derv2(kx,jy) = derv2(kx,jy) - d2trm*yjl*xik
    derv2(ky,jy) = derv2(ky,jy) - d2trm*yjl*yik
    derv2(kz,jy) = derv2(kz,jy) - d2trm*yjl*zik
    derv2(kx,jz) = derv2(kx,jz) - d2trm*zjl*xik
    derv2(ky,jz) = derv2(ky,jz) - d2trm*zjl*yik
    derv2(kz,jz) = derv2(kz,jz) - d2trm*zjl*zik
  else
    derv2(jx,kx) = derv2(jx,kx) - d2trm*xjl*xik
    derv2(jy,kx) = derv2(jy,kx) - d2trm*yjl*xik
    derv2(jz,kx) = derv2(jz,kx) - d2trm*zjl*xik
    derv2(jx,ky) = derv2(jx,ky) - d2trm*xjl*yik
    derv2(jy,ky) = derv2(jy,ky) - d2trm*yjl*yik
    derv2(jz,ky) = derv2(jz,ky) - d2trm*zjl*yik
    derv2(jx,kz) = derv2(jx,kz) - d2trm*xjl*zik
    derv2(jy,kz) = derv2(jy,kz) - d2trm*yjl*zik
    derv2(jz,kz) = derv2(jz,kz) - d2trm*zjl*zik
  endif
!
!  K-L
!     
  if (lx.lt.kx) then
    derv2(lx,kx) = derv2(lx,kx) + d2trm*xik*xjl
    derv2(ly,kx) = derv2(ly,kx) + d2trm*xik*yjl
    derv2(lz,kx) = derv2(lz,kx) + d2trm*xik*zjl
    derv2(lx,ky) = derv2(lx,ky) + d2trm*yik*xjl
    derv2(ly,ky) = derv2(ly,ky) + d2trm*yik*yjl
    derv2(lz,ky) = derv2(lz,ky) + d2trm*yik*zjl
    derv2(lx,kz) = derv2(lx,kz) + d2trm*zik*xjl
    derv2(ly,kz) = derv2(ly,kz) + d2trm*zik*yjl
    derv2(lz,kz) = derv2(lz,kz) + d2trm*zik*zjl
  else
    derv2(kx,lx) = derv2(kx,lx) + d2trm*xik*xjl
    derv2(ky,lx) = derv2(ky,lx) + d2trm*yik*xjl
    derv2(kz,lx) = derv2(kz,lx) + d2trm*zik*xjl
    derv2(kx,ly) = derv2(kx,ly) + d2trm*xik*yjl
    derv2(ky,ly) = derv2(ky,ly) + d2trm*yik*yjl
    derv2(kz,ly) = derv2(kz,ly) + d2trm*zik*yjl
    derv2(kx,lz) = derv2(kx,lz) + d2trm*xik*zjl
    derv2(ky,lz) = derv2(ky,lz) + d2trm*yik*zjl
    derv2(kz,lz) = derv2(kz,lz) + d2trm*zik*zjl
  endif
!
!  Strain terms
!
  if (lstr) then
    do kl = 1,nstrains
      ks = nstrptr(kl)
      derv3(ix,kl) = derv3(ix,kl) - xik*d2trm*dr2jlds(ks)
      derv3(iy,kl) = derv3(iy,kl) - yik*d2trm*dr2jlds(ks)
      derv3(iz,kl) = derv3(iz,kl) - zik*d2trm*dr2jlds(ks)
      derv3(kx,kl) = derv3(kx,kl) + xik*d2trm*dr2jlds(ks)
      derv3(ky,kl) = derv3(ky,kl) + yik*d2trm*dr2jlds(ks)
      derv3(kz,kl) = derv3(kz,kl) + zik*d2trm*dr2jlds(ks)
    enddo
!
    do kl = 1,nstrains
      ks = nstrptr(kl)
      derv3(jx,kl) = derv3(jx,kl) - xjl*d2trm*dr2ikds(ks)
      derv3(jy,kl) = derv3(jy,kl) - yjl*d2trm*dr2ikds(ks)
      derv3(jz,kl) = derv3(jz,kl) - zjl*d2trm*dr2ikds(ks)
      derv3(lx,kl) = derv3(lx,kl) + xjl*d2trm*dr2ikds(ks)
      derv3(ly,kl) = derv3(ly,kl) + yjl*d2trm*dr2ikds(ks)
      derv3(lz,kl) = derv3(lz,kl) + zjl*d2trm*dr2ikds(ks)
    enddo
!
    do kk = 1,nstrains
      ks = nstrptr(kk)
      do kl = 1,nstrains
        kt = nstrptr(kl)
        sderv2(kl,kk) = sderv2(kl,kk) + d2trm*(dr2ikds(kt)*dr2jlds(ks) + dr2jlds(kt)*dr2ikds(ks))
      enddo
    enddo
  endif
#ifdef TRACE
  call trace_out('add_drv2_2ut')
#endif
!
  return
  end
