  subroutine getinbox
!
!  Set coordinates to be in the box in xinbox/yinbox/zinbox
!
!   7/21 Created
!   1/22 Corrected for 1 and 2-D cases
!   3/22 Further correction for 1- and 2-D cases
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
  use current,        only : numat, ndim
  use current,        only : xfrac, yfrac, zfrac, rv
  use current,        only : xclat, yclat, zclat
  use spatial,        only : xinbox
  use spatial,        only : yinbox
  use spatial,        only : zinbox
  use iochannels
  use gulp_gfnff
  use mdlogic,        only : lmd
  use neighbours
  use spatial
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
  real(dp)                                       :: x
  real(dp)                                       :: y
  real(dp)                                       :: z
#ifdef TRACE
  call trace_in('getinbox')
#endif
!
!  If this is a spatial run then coordinates will be set elsewhere
!
  if (lspatialok) return
!
!  If MD or non-periodic then the atoms are already in the box in xclat/yclat/zclat
!
  if (lmd.or.ndim.eq.0) then
    xinbox(1:numat) = xclat(1:numat)
    yinbox(1:numat) = yclat(1:numat)
    zinbox(1:numat) = zclat(1:numat)
    return
  endif
!******************************
!  Find image of atom in box  *
!******************************
  if (ndim.eq.3) then
    do i = 1,numat
      x = xfrac(i)
      y = yfrac(i)
      z = zfrac(i)
!
      x = mod(x+100.0_dp,1.0_dp)
      y = mod(y+100.0_dp,1.0_dp)
      z = mod(z+100.0_dp,1.0_dp)
!
      xinbox(i) = rv(1,1)*x + rv(1,2)*y + rv(1,3)*z
      yinbox(i) = rv(2,1)*x + rv(2,2)*y + rv(2,3)*z
      zinbox(i) = rv(3,1)*x + rv(3,2)*y + rv(3,3)*z
    enddo
  elseif (ndim.eq.2) then
    do i = 1,numat
      x = xfrac(i)
      y = yfrac(i)
!
      x = mod(x+100.0_dp,1.0_dp)
      y = mod(y+100.0_dp,1.0_dp)
!
      xinbox(i) = rv(1,1)*x + rv(1,2)*y
      yinbox(i) = rv(2,1)*x + rv(2,2)*y
      zinbox(i) = zclat(i)
    enddo
  elseif (ndim.eq.1) then
    do i = 1,numat
      x = xfrac(i)
!
      x = mod(x+100.0_dp,1.0_dp)
!
      xinbox(i) = rv(1,1)*x
      yinbox(i) = yclat(i)
      zinbox(i) = zclat(i)
    enddo
  endif
#ifdef TRACE
  call trace_out('getinbox')
#endif
!
  end subroutine getinbox
