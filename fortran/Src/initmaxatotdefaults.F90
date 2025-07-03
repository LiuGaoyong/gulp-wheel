  subroutine initmaxatotdefaults(i)
!
!  Initialises the arrays associated with maxatot
!
!   9/10 Created from changemaxatot
!   7/15 External potential added
!   3/18 Sign option added to translate
!  11/21 Modifications for TI added
!  11/21 Einstein positions changed
!   7/23 lbsmatcfg added
!   7/23 nregionnocfg added
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
  use configurations
  use moldyn,        only : lfix
  use scan,          only : ltranat, ltranatminus
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: i
!
!  Initialise new parts of data arrays
!
  if (i.ge.1.and.i.le.maxatot) then
    lbsmatcfg(i) = .false.
    leinsteinat(i) = .false.
    ltibeinsteinat(i) = .false.
    ltifeinsteinat(i) = .false.
    lfix(i) = .false.
    lopfi(3*i-2:3*i) = .false.
    ltdforcecfg(1:3,i) = .false.
    lqmatom(i) = .false.
    lsliceatom(i) = .false.
    ltranat(i) = .false.
    ltranatminus(i) = .false.
    nregionnocfg(i) = 1
    extpotcfg(i) = 0.0_dp
    forcecfg(1:3,i) = 0.0_dp
    tdforcecfg(1:3,1:3,i) = 0.0_dp
    keinsteinat(i) = 0.0_dp
    xceinsteinat(i) = 0.0_dp
    yceinsteinat(i) = 0.0_dp
    zceinsteinat(i) = 0.0_dp
    xfeinsteinat(i) = 0.0_dp
    yfeinsteinat(i) = 0.0_dp
    zfeinsteinat(i) = 0.0_dp
  endif
!
  return
  end
