module m_cells
! 
!  This module contains the infrastructure to compute cell images within a cutoff
!
!  Julian Gale, Curtin University, September 2021
!
  use datatypes
!
  implicit none

  type t_cells
!
!  Cell pair list
!
    integer(i4)                               :: maxcell = 0     ! Maximum number of cell vectors for pairs
    integer(i4)                               :: ncell           ! Number of cell pair vectors
    integer(i4)                               :: ncentrecell     ! Number of cell image for the centre cell
    real(dp),    dimension(:),       pointer  :: xcell => null() ! x component of cell pair list
    real(dp),    dimension(:),       pointer  :: ycell => null() ! y component of cell pair list
    real(dp),    dimension(:),       pointer  :: zcell => null() ! z component of cell pair list
  end type t_cells

CONTAINS

  subroutine changemaxcell(cell)
!
!  Changes the size of arrays that hold the cell list
!
!  Julian Gale, CIC, Curtin University, September 2021
!
  use reallocate
  implicit none
!
!  Passed variables
!
  type(t_cells), intent(inout) :: cell
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(cell%xcell,cell%maxcell,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcell','xcell')
  call realloc(cell%ycell,cell%maxcell,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcell','ycell')
  call realloc(cell%zcell,cell%maxcell,ierror)
  if (ierror.ne.0) call outofmemory('changemaxcell','zcell')
!
  end subroutine changemaxcell
!
  subroutine setcells(cut2,cell)
!
!  Store linear array of lattice vectors for required cell images
!
!   9/21 Created 
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
  use current
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),         intent(in)    :: cut2      ! Cutoff squared for interactions
  type(t_cells),    intent(inout) :: cell      ! Type to hold cell information
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
#ifdef TRACE
  call trace_in('setcells')
#endif
!
!  Check cell type
!
  lminimal =  .true.
  if (alpha.lt.98.0_dp) lminimal =  .false.
  if (beta.lt.98.0_dp)  lminimal =  .false.
  if (gamma.lt.98.0_dp) lminimal =  .false.
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
    cell%ncell = (2*na + 1)*(2*nb + 1)*(2*nc + 1)
    if (cell%ncell.gt.cell%maxcell) then
      cell%maxcell = cell%ncell 
      call changemaxcell(cell)
    endif
!
    xcdi = - dble(na+1)*r1x
    ycdi = - dble(na+1)*r1y
    zcdi = - dble(na+1)*r1z
    cell%ncell = 0
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
          cell%ncell = cell%ncell + 1
          xcrd = xcrd + r3x
          ycrd = ycrd + r3y
          zcrd = zcrd + r3z
          cell%xcell(cell%ncell) = xcrd
          cell%ycell(cell%ncell) = ycrd
          cell%zcell(cell%ncell) = zcrd
          if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cell%ncentrecell = cell%ncell
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
    cell%ncell = (2*na + 1)*(2*nb + 1)
    if (cell%ncell.gt.cell%maxcell) then
      cell%maxcell = cell%ncell
      call changemaxcell(cell)
    endif
!
    xcdi = - dble(na+1)*r1x
    ycdi = - dble(na+1)*r1y
    cell%ncell = 0
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
        cell%ncell = cell%ncell + 1
        cell%xcell(cell%ncell) = xcrd
        cell%ycell(cell%ncell) = ycrd
        cell%zcell(cell%ncell) = 0.0_dp
        if (ii.eq.0.and.jj.eq.0) cell%ncentrecell = cell%ncell
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
    cell%ncell = (2*na + 1)
    if (cell%ncell.gt.cell%maxcell) then
      cell%maxcell = cell%ncell
      call changemaxcell(cell)
    endif
    xcrd = - dble(na+1)*r1x
    cell%ncell = 0
!
!  Loop over unit cells
!
    do ii = -na,na
      xcrd = xcrd + r1x
      cell%ncell = cell%ncell + 1
      cell%xcell(cell%ncell) = xcrd
      cell%ycell(cell%ncell) = 0.0_dp
      cell%zcell(cell%ncell) = 0.0_dp
    enddo
    cell%ncentrecell = na + 1
  else
!*************
!  0-D case  *
!*************
!
!  Check memory
!
    cell%ncell = 1
    if (cell%ncell.gt.cell%maxcell) then
      cell%maxcell = cell%ncell
      call changemaxcell(cell)
    endif
    cell%ncentrecell = 1
    cell%xcell(1) = 0.0_dp
    cell%ycell(1) = 0.0_dp
    cell%zcell(1) = 0.0_dp
  endif
#ifdef TRACE
  call trace_out('setcells')
#endif
!
  return
  end

end module m_cells
