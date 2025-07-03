  subroutine changemaxlambda
!
!  Alters the size of the arrays associated with maxlambda
!
!  11/21 Created
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
!  Julian Gale, CIC, Curtin University, November 2021
!
  use m_ti
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxlambda = 0
!
  call realloc(dUdlambda,maxlambda,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat','dUdlambda')
!
!  Initialise new parts of data arrays
!
  if (maxlambda.gt.oldmaxlambda) then
    do i = oldmaxlambda+1,maxlambda
      dUdlambda(i) = 0.0_dp
    enddo
  endif
!
!  Save current value of maxlambda for next call
!
  oldmaxlambda = maxlambda
!
  return
  end
