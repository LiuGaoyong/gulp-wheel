  subroutine initmaxobsdefaults(i)
!
!  Initialises the arrays associated with maxobs
!
!   9/10 Created from changemax routine
!  11/22 nobstyle added
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
!  Julian Gale, CIC, Curtin University, November 2022
!
  use configurations,   only : maxcfg
  use observables
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: i
!
!  Initialise new parts of data arrays
!
  if (i.ge.1.and.i.le.maxobs) then
    nobstyle(i) = 0
    weight(i) = 0.0_dp
    fparameter(i) = 0.0_dp
    freaction(1:maxcfg,i) = 0.0_dp
  endif
!
  return
  end
