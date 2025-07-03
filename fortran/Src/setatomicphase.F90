  subroutine setatomicphase(xk,yk,zk)
!
!  Calculates the atomic phase factors for later use
!
!   8/22 Created
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
!  Julian Gale, CIC, Curtin University, August 2022
!
  use g_constants
  use control
  use current
  use datatypes
  use frequencies
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
  real(dp),    intent(in)                      :: xk
  real(dp),    intent(in)                      :: yk
  real(dp),    intent(in)                      :: zk
!
!  Local variables
!
  integer(i4)                                  :: i
  real(dp)                                     :: cosk
  real(dp)                                     :: sink
  real(dp)                                     :: kvf(3,3)
  real(dp)                                     :: xkv
  real(dp)                                     :: ykv
  real(dp)                                     :: zkv
#ifdef TRACE
  call trace_in('setatomicphase')
#endif
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
!********************
!  Loop over atoms  *
!********************
  do i = 1,numat
!
!  Compute phase factor
!
    cosk = - (xkv*xclat(i) + ykv*yclat(i) + zkv*zclat(i))
    sink = sin(cosk)
    cosk = cos(cosk)
    aphase(i) = dcmplx(cosk,sink)
  enddo
#ifdef TRACE
  call trace_out('setatomicphase')
#endif
!
  return
  end
