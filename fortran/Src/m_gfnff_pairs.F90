module m_gfnff_pairs
! 
!  This module contains the infrastructure to compute the list of pairwise
!  interactions between images of two atoms for GFNFF.
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
!  Julian Gale, Curtin University, October 2020
!
  use datatypes
!
  implicit none

  type t_pairs
!
!  General pair list
!
    integer(i4)                              :: maxpair = 0         ! Maximum number of pairs
    integer(i4)                              :: npair               ! Number of pairs
    real(dp),    dimension(:),       pointer :: r2pair => null()    ! Distance squared for pair
    real(dp),    dimension(:),       pointer :: xpair => null()     ! x component of vector for pair
    real(dp),    dimension(:),       pointer :: ypair => null()     ! y component of vector for pair
    real(dp),    dimension(:),       pointer :: zpair => null()     ! z component of vector for pair
  end type t_pairs

  type t_triads
!
!  General triad list
!
    integer(i4)                              :: maxtriad = 0        ! Maximum number of triads
    integer(i4)                              :: ntriad              ! Number of triads
    real(dp),    dimension(:),       pointer :: r2ikpair => null()  ! Distance squared for i-k pair
    real(dp),    dimension(:),       pointer :: r2jkpair => null()  ! Distance squared for j-k pair
    real(dp),    dimension(:),       pointer :: xikpair => null()   ! x component of vector for i-k pair
    real(dp),    dimension(:),       pointer :: yikpair => null()   ! y component of vector for i-k pair
    real(dp),    dimension(:),       pointer :: zikpair => null()   ! z component of vector for i-k pair
    real(dp),    dimension(:),       pointer :: xjkpair => null()   ! x component of vector for j-k pair
    real(dp),    dimension(:),       pointer :: yjkpair => null()   ! y component of vector for j-k pair
    real(dp),    dimension(:),       pointer :: zjkpair => null()   ! z component of vector for j-k pair
  end type t_triads

  type t_paircell
!
!  Cell pair list
!
    integer(i4)                               :: maxpaircell = 0     ! Maximum number of cell vectors for pairs
    integer(i4)                               :: npaircell           ! Number of cell pair vectors
    integer(i4)                               :: ncentrecell         ! Number of cell image for the centre cell
    real(dp),    dimension(:),       pointer  :: xpaircell => null() ! x component of cell pair list
    real(dp),    dimension(:),       pointer  :: ypaircell => null() ! y component of cell pair list
    real(dp),    dimension(:),       pointer  :: zpaircell => null() ! z component of cell pair list
  end type t_paircell
!
!  Pair list data
!
  type(t_paircell),                      save :: cnhb_paircell       ! Cell for pair list for coordination number in hydrogen bonding
  type(t_pairs),                         save :: cnhb_pairs          ! Pair list for coordination number in hydrogen bonding
  type(t_paircell),                      save :: hb1_paircell        ! Cell for pair list for hydrogen bonding with hbthr1
  type(t_pairs),                         save :: hb1_pairs           ! Pair list for hydrogen bonding with hbthr1
  type(t_pairs),                         save :: hb2_pairs           ! Pair list for hydrogen bonding with hbthr2
  type(t_paircell),                      save :: hb2_paircell        ! Cell for pair list for hydrogen bonding with hbthr2
  type(t_triads),                        save :: hb2_triads          ! Triad list for hydrogen bonding with hbthr2
  type(t_paircell),                      save :: lr_paircell         ! Cell for pair list for long-range interactions
  type(t_pairs),                         save :: lr_pairs            ! Pair list for long-range interactions

CONTAINS

  subroutine changemaxpair(pairs)
!
!  Changes the size of arrays that hold the pair list for GFNFF
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
!  Julian Gale, CIC, Curtin University, October 2020
!
  use reallocate
  implicit none
!
!  Passed variables
!
  type(t_pairs), intent(inout) :: pairs
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(pairs%r2pair,pairs%maxpair,ierror)
  if (ierror.ne.0) call outofmemory('changemaxpair','r2pair')
  call realloc(pairs%xpair,pairs%maxpair,ierror)
  if (ierror.ne.0) call outofmemory('changemaxpair','xpair')
  call realloc(pairs%ypair,pairs%maxpair,ierror)
  if (ierror.ne.0) call outofmemory('changemaxpair','ypair')
  call realloc(pairs%zpair,pairs%maxpair,ierror)
  if (ierror.ne.0) call outofmemory('changemaxpair','zpair')
!
  end subroutine changemaxpair
!
  subroutine changemaxtriad(triads)
!
!  Changes the size of arrays that hold the triad list for GFNFF
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
!  Julian Gale, CIC, Curtin University, October 2020
!
  use reallocate
  implicit none
!
!  Passed variables
!
  type(t_triads), intent(inout) :: triads
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(triads%r2ikpair,triads%maxtriad,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtriad','r2ikpair')
  call realloc(triads%r2jkpair,triads%maxtriad,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtriad','r2jkpair')
  call realloc(triads%xikpair,triads%maxtriad,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtriad','xikpair')
  call realloc(triads%yikpair,triads%maxtriad,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtriad','yikpair')
  call realloc(triads%zikpair,triads%maxtriad,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtriad','zikpair')
  call realloc(triads%xjkpair,triads%maxtriad,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtriad','xjkpair')
  call realloc(triads%yjkpair,triads%maxtriad,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtriad','yjkpair')
  call realloc(triads%zjkpair,triads%maxtriad,ierror)
  if (ierror.ne.0) call outofmemory('changemaxtriad','zjkpair')
!
  end subroutine changemaxtriad
!
  subroutine changemaxpaircell(paircell)
!
!  Changes the size of arrays that hold the cell pair list for GFNFF
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
!  Julian Gale, CIC, Curtin University, October 2020
!
  use reallocate
  implicit none
!
!  Passed variables
!
  type(t_paircell), intent(inout) :: paircell
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(paircell%xpaircell,paircell%maxpaircell,ierror)
  if (ierror.ne.0) call outofmemory('changemaxpaircell','xpaircell')
  call realloc(paircell%ypaircell,paircell%maxpaircell,ierror)
  if (ierror.ne.0) call outofmemory('changemaxpaircell','ypaircell')
  call realloc(paircell%zpaircell,paircell%maxpaircell,ierror)
  if (ierror.ne.0) call outofmemory('changemaxpaircell','zpaircell')
!
  end subroutine changemaxpaircell
!
  subroutine gfnff_setpaircell(cut2,paircell)
!
!  Store linear array of lattice vectors for required cell images
!
!  10/20 Created from rlist2
!   2/21 Number of cells incremented to ensure that all distances
!        are found when cell angles are small.
!   3/22 Angle checks now wrapped for dimensionality
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
  use current
  implicit none
!
!  Passed variables
!
  real(dp),         intent(in)    :: cut2      ! Cutoff squared for interactions
  type(t_paircell), intent(inout) :: paircell  ! Type to hold paircell information
!
!  Local variables
!
  integer(i4)                     :: na
  integer(i4)                     :: nb
  integer(i4)                     :: nc
  integer(i4)                     :: ii
  integer(i4)                     :: jj
  integer(i4)                     :: kk
  logical                         :: lminimal
  real(dp)                        :: cut
  real(dp)                        :: xcdi
  real(dp)                        :: ycdi
  real(dp)                        :: zcdi
  real(dp)                        :: xcdj
  real(dp)                        :: ycdj
  real(dp)                        :: zcdj
  real(dp)                        :: xcrd
  real(dp)                        :: ycrd
  real(dp)                        :: zcrd
!
!  Check cell type
!
  lminimal =  .true.
  if (ndim.eq.3) then
    if (alpha.lt.98.0_dp) lminimal =  .false.
    if (beta.lt.98.0_dp)  lminimal =  .false.
    if (gamma.lt.98.0_dp) lminimal =  .false.
  elseif (ndim.eq.2) then
    if (alpha.lt.98.0_dp) lminimal =  .false.
  endif
!
  if (ndim.eq.3) then
!*************
!  3-D case  *
!*************
!
!  Find cell image dimensions
!
    cut = sqrt(cut2)
    if (lminimal) then
      na = (cut/a) + 1
      nb = (cut/b) + 1
      nc = (cut/c) + 1
    else
      na = nint(cut/a) + 2
      nb = nint(cut/b) + 2
      nc = nint(cut/c) + 2
    endif
!
!  Check memory
!
    paircell%npaircell = (2*na + 1)*(2*nb + 1)*(2*nc + 1)
    if (paircell%npaircell.gt.paircell%maxpaircell) then
      paircell%maxpaircell = paircell%npaircell 
      call changemaxpaircell(paircell)
    endif
!
    xcdi = - dble(na+1)*r1x
    ycdi = - dble(na+1)*r1y
    zcdi = - dble(na+1)*r1z
    paircell%npaircell = 0
!
!  Loop over unit cells
!
    do ii = -na,na
      xcdi = xcdi + r1x
      ycdi = ycdi + r1y
      zcdi = zcdi + r1z
      xcdj = xcdi - dble(nb+1)*r2x
      ycdj = ycdi - dble(nb+1)*r2y
      zcdj = zcdi - dble(nb+1)*r2z
      do jj = -nb,nb
        xcdj = xcdj + r2x
        ycdj = ycdj + r2y
        zcdj = zcdj + r2z
        xcrd = xcdj - dble(nc+1)*r3x
        ycrd = ycdj - dble(nc+1)*r3y
        zcrd = zcdj - dble(nc+1)*r3z
        do kk = -nc,nc
          paircell%npaircell = paircell%npaircell + 1
          xcrd = xcrd + r3x
          ycrd = ycrd + r3y
          zcrd = zcrd + r3z
          paircell%xpaircell(paircell%npaircell) = xcrd
          paircell%ypaircell(paircell%npaircell) = ycrd
          paircell%zpaircell(paircell%npaircell) = zcrd
          if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) paircell%ncentrecell = paircell%npaircell
        enddo
      enddo
    enddo
  elseif (ndim.eq.2) then
!*************
!  2-D case  *
!*************
!
!  Find cell image dimensions
!
    cut = sqrt(cut2)
    if (lminimal) then
      na = (cut/a) + 1
      nb = (cut/b) + 1
    else
      na = nint(cut/a) + 2
      nb = nint(cut/b) + 2
    endif
!
!  Check memory
!
    paircell%npaircell = (2*na + 1)*(2*nb + 1)
    if (paircell%npaircell.gt.paircell%maxpaircell) then
      paircell%maxpaircell = paircell%npaircell
      call changemaxpaircell(paircell)
    endif
!
    xcdi = - dble(na+1)*r1x
    ycdi = - dble(na+1)*r1y
    paircell%npaircell = 0
!
!  Loop over unit cells
!
    do ii = -na,na
      xcdi = xcdi + r1x
      ycdi = ycdi + r1y
      xcrd = xcdi - dble(nb+1)*r2x
      ycrd = ycdi - dble(nb+1)*r2y
      do jj = -nb,nb
        xcrd = xcrd + r2x
        ycrd = ycrd + r2y
        paircell%npaircell = paircell%npaircell + 1
        paircell%xpaircell(paircell%npaircell) = xcrd
        paircell%ypaircell(paircell%npaircell) = ycrd
        paircell%zpaircell(paircell%npaircell) = 0.0_dp
        if (ii.eq.0.and.jj.eq.0) paircell%ncentrecell = paircell%npaircell
      enddo
    enddo
  elseif (ndim.eq.1) then
!*************
!  1-D case  *
!*************
!
!  Find cell image dimensions
!
    cut = sqrt(cut2)
    if (lminimal) then
      na = (cut/a) + 1
    else
      na = nint(cut/a) + 1
    endif
!
!  Check memory
!
    paircell%npaircell = (2*na + 1)
    if (paircell%npaircell.gt.paircell%maxpaircell) then
      paircell%maxpaircell = paircell%npaircell
      call changemaxpaircell(paircell)
    endif
    xcrd = - dble(na+1)*r1x
    paircell%npaircell = 0
!
!  Loop over unit cells
!
    do ii = -na,na
      xcrd = xcrd + r1x
      paircell%npaircell = paircell%npaircell + 1
      paircell%xpaircell(paircell%npaircell) = xcrd
      paircell%ypaircell(paircell%npaircell) = 0.0_dp
      paircell%zpaircell(paircell%npaircell) = 0.0_dp
    enddo
    paircell%ncentrecell = na + 1
  else
!*************
!  0-D case  *
!*************
!
!  Check memory
!
    paircell%npaircell = 1
    if (paircell%npaircell.gt.paircell%maxpaircell) then
      paircell%maxpaircell = paircell%npaircell
      call changemaxpaircell(paircell)
    endif
    paircell%ncentrecell = 1
    paircell%xpaircell(1) = 0.0_dp
    paircell%ypaircell(1) = 0.0_dp
    paircell%zpaircell(1) = 0.0_dp
  endif
!
  return
  end
!
  subroutine gfnff_getpairs(i,j,cut2,paircell,pairs)
!
!  Computes the pair list for GFNFF 
!
!  On entry : 
!
!  cut2            = cutoff squared
!  paircell        = cell information for search over images
!
!  On exit :
!
!  pairs           = type with pair information
!
!  10/20 Created
!   3/21 Intent of pairs corrected to inout
!   6/21 Use of xclat/yclat/zclat changed to xalat/yalat/zalat for benefit of MD
!   7/21 In non-MD case xclat/yclat/zclat are now used to allow for symmetry
!   7/21 Coordinates now use xclat/yclat/zclat for all cases
!   7/21 Coordinates now use inbox arrays
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
!  Julian Gale, CIC, Curtin University, July 2021
!
  use datatypes
  use current,        only : ndim
  use spatial,        only : xinbox, yinbox, zinbox
  implicit none
!
!  Passed variables
!
  integer(i4),                 intent(in)        :: i
  integer(i4),                 intent(in)        :: j
  real(dp),                    intent(in)        :: cut2
  type(t_paircell),            intent(in)        :: paircell
  type(t_pairs),               intent(inout)     :: pairs
!
!  Local variables
!
  integer(i4)                                    :: ii
  real(dp)                                       :: r2
  real(dp)                                       :: xji
  real(dp)                                       :: yji
  real(dp)                                       :: zji
  real(dp)                                       :: xji0
  real(dp)                                       :: yji0
  real(dp)                                       :: zji0
!
!  Set number of pairs to zero
!
  pairs%npair = 0
!******************************************
!  Calculate pair lists for atoms i and j *
!******************************************
!
!  Set centre cell coordinate differences
!
  xji0 = xinbox(j) - xinbox(i)
  yji0 = yinbox(j) - yinbox(i)
  zji0 = zinbox(j) - zinbox(i)
!********
!  PBC  *
!********
  if (ndim.gt.0) then
!
!  Loop over unit cells
!
    do ii = 1,paircell%npaircell
!
!  Exclude self term
!
      if (i.ne.j.or.ii.ne.paircell%ncentrecell) then
        xji = xji0 + paircell%xpaircell(ii)
        yji = yji0 + paircell%ypaircell(ii)
        zji = zji0 + paircell%zpaircell(ii)
        r2 = xji*xji + yji*yji + zji*zji
        if (r2.lt.cut2) then
          if (pairs%npair.ge.pairs%maxpair) then
            pairs%maxpair = pairs%maxpair + 4
            call changemaxpair(pairs)
          endif
          pairs%npair = pairs%npair + 1
          pairs%r2pair(pairs%npair) = r2
          pairs%xpair(pairs%npair) = xji
          pairs%ypair(pairs%npair) = yji
          pairs%zpair(pairs%npair) = zji
        endif
      endif
    enddo
  else
!********
!  0-D  *
!********
    r2 = xji0*xji0 + yji0*yji0 + zji0*zji0
!
!  Check that pair is within cutoff and not a self-term
!
    if (r2.lt.cut2.and.i.ne.j) then
      pairs%npair = pairs%npair + 1
      if (pairs%npair.ge.pairs%maxpair) then
        pairs%maxpair = pairs%maxpair + 1
        call changemaxpair(pairs)
      endif
      pairs%r2pair(pairs%npair) = r2
      pairs%xpair(pairs%npair) = xji0
      pairs%ypair(pairs%npair) = yji0
      pairs%zpair(pairs%npair) = zji0
    endif
  endif
!
  end subroutine gfnff_getpairs

  subroutine gfnff_gettriads(ndim,i,j,k,r2ij,xij,yij,zij,cut2,paircell,triads)
!
!  Computes the a triad list for GFNFF for 3 atoms where the 
!  combined sum of the distances squared is less than cut2.
!
!  On entry : 
!
!  ndim            = number of dimensions
!  i, j            = atoms whose position is fixed for triad
!  k               = atom whose image is being searched for
!  r2ij            = distance squared from i to j
!  xij             = x component of vector from i to j
!  yij             = y component of vector from i to j
!  zij             = z component of vector from i to j
!  cut2            = cutoff squared
!  paircell        = cell information for search over images
!
!  On exit :
!
!  triads           = type with triad information
!
!  10/20 Created
!   3/21 Intent of triads corrected to inout
!   6/21 Use of xclat/yclat/zclat changed to xalat/yalat/zalat for benefit of MD
!   7/21 In non-MD case xclat/yclat/zclat are now used to allow for symmetry
!   7/21 Coordinates now use xclat/yclat/zclat for all cases
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
!  Julian Gale, CIC, Curtin University, July 2021
!
  use datatypes
  use spatial,        only : xinbox, yinbox, zinbox
  implicit none
!
!  Passed variables
!
  integer(i4),                 intent(in)        :: ndim
  integer(i4),                 intent(in)        :: i
  integer(i4),                 intent(in)        :: j
  integer(i4),                 intent(in)        :: k
  real(dp),                    intent(in)        :: cut2
  real(dp),                    intent(in)        :: r2ij
  real(dp),                    intent(in)        :: xij
  real(dp),                    intent(in)        :: yij
  real(dp),                    intent(in)        :: zij
  type(t_paircell),            intent(in)        :: paircell
  type(t_triads),              intent(inout)     :: triads
!
!  Local variables
!
  integer(i4)                                    :: ii
  real(dp)                                       :: r2ik
  real(dp)                                       :: r2jk
  real(dp)                                       :: xki
  real(dp)                                       :: yki
  real(dp)                                       :: zki
  real(dp)                                       :: xki0
  real(dp)                                       :: yki0
  real(dp)                                       :: zki0
  real(dp)                                       :: xkj
  real(dp)                                       :: ykj
  real(dp)                                       :: zkj
  real(dp)                                       :: xkj0
  real(dp)                                       :: ykj0
  real(dp)                                       :: zkj0
!
!  Set number of triads to zero
!
  triads%ntriad = 0
!
!  Check that r2ij doesn't exceed the cutoff already
!
  if (r2ij.gt.cut2) return
!**********************************************
!  Calculate triad lists for atoms i, j and k *
!**********************************************
!
!  Set centre cell coordinate differences
!
  xki0 = xinbox(k) - xinbox(i)
  yki0 = yinbox(k) - yinbox(i)
  zki0 = zinbox(k) - zinbox(i)
!
  xkj0 = xki0 - xij
  ykj0 = yki0 - yij
  zkj0 = zki0 - zij
!********
!  PBC  *
!********
  if (ndim.gt.0) then
!
!  Loop over unit cells
!
    do ii = 1,paircell%npaircell
!
!  Exclude self terms
!
      if (i.eq.k.and.ii.eq.paircell%ncentrecell) cycle
      if (j.eq.k.and.ii.eq.paircell%ncentrecell) cycle
!
      xki = xki0 + paircell%xpaircell(ii)
      yki = yki0 + paircell%ypaircell(ii)
      zki = zki0 + paircell%zpaircell(ii)
      r2ik = xki*xki + yki*yki + zki*zki
!
!  Check that we are still within the cutoff after first 2 distances
!
      if (r2ij+r2ik.gt.cut2) cycle
!
      xkj = xkj0 + paircell%xpaircell(ii)
      ykj = ykj0 + paircell%ypaircell(ii)
      zkj = zkj0 + paircell%zpaircell(ii)
      r2jk = xkj*xkj + ykj*ykj + zkj*zkj
!
!  Check that we are still within the cutoff with all distances
!
      if (r2ij+r2ik+r2jk.gt.cut2) cycle
!
      if (triads%ntriad.ge.triads%maxtriad) then
        triads%maxtriad = triads%maxtriad + 4
        call changemaxtriad(triads)
      endif
      triads%ntriad = triads%ntriad + 1
!
      triads%r2ikpair(triads%ntriad) = r2ik
      triads%xikpair(triads%ntriad) = xki
      triads%yikpair(triads%ntriad) = yki
      triads%zikpair(triads%ntriad) = zki
!
      triads%r2jkpair(triads%ntriad) = r2jk
      triads%xjkpair(triads%ntriad) = xkj
      triads%yjkpair(triads%ntriad) = ykj
      triads%zjkpair(triads%ntriad) = zkj
    enddo
  else
!********
!  0-D  *
!********
    r2ik = xki0*xki0 + yki0*yki0 + zki0*zki0
!
!  Check that triad is within cutoff and does not include self-terms
!
    if (r2ij+r2ik.lt.cut2.and.i.ne.k.and.j.ne.k) then
      r2jk = xkj0*xkj0 + ykj0*ykj0 + zkj0*zkj0
      if (r2ij+r2ik+r2jk.lt.cut2) then
        if (triads%ntriad.ge.triads%maxtriad) then
          triads%maxtriad = triads%maxtriad + 1
          call changemaxtriad(triads)
        endif
        triads%ntriad = triads%ntriad + 1
!
        triads%r2ikpair(triads%ntriad) = r2ik
        triads%xikpair(triads%ntriad) = xki0
        triads%yikpair(triads%ntriad) = yki0
        triads%zikpair(triads%ntriad) = zki0
!
        triads%r2jkpair(triads%ntriad) = r2jk
        triads%xjkpair(triads%ntriad) = xkj0
        triads%yjkpair(triads%ntriad) = ykj0
        triads%zjkpair(triads%ntriad) = zkj0
      endif
    endif
  endif
!
  end subroutine gfnff_gettriads

end module m_gfnff_pairs
