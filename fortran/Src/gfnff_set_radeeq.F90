  subroutine gfnff_set_radeeq(numat,radeeq)
!
!  Initialises radeeq value 
!
!  10/21 Created
!  12/21 Changed to use gfnff_eeq_alp
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
!  Julian Gale, CIC, Curtin University, December 2021
!
  use datatypes
  use gulp_gfnff,   only : gfnff_eeq_alp
  implicit none
!
!  Passed variables
!
  integer(i4),       intent(in)  :: numat
  real(dp),          intent(out) :: radeeq
!
!  Local variables
!
  integer(i4)                    :: i
  real(dp)                       :: alpmax
!
!  Set GFNFF cutoff 
!
  alpmax = 1.0d-12
  do i = 1,numat
    alpmax = max(alpmax,sqrt(gfnff_eeq_alp(i)))
  enddo
  radeeq = 4.0_dp*alpmax
!
  return
  end
