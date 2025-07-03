  subroutine changemaxd2q
!
!  Alters the size of the arrays associated with maxd2q or maxd2qu
!
!   1/22 dqds moved here from changemaxatloc
!   2/22 d2edqdq changed to have maxd2qu as right hand dimension
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
!  Julian Gale, CIC, Curtin University, February 2022
!
  use derivatives
  use current,       only : maxat
  use reallocate
  implicit none
!
  integer(i4) :: ierror
!
  call realloc(dqdxyz,maxd2q,maxd2qu,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd2q','dqdxyz')
  call realloc(dqds,6_i4,maxd2qu,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd2q','dqds')
  call realloc(d2edqc,maxd2,maxd2qu,ierror)
  if (ierror.gt.0) call outofmemory('changemaxd2q','d2edqc')
  call realloc(d2edqdq,maxat,maxd2qu,ierror)
  if (ierror.gt.0) call outofmemory('changemaxd2q','d2edqdq')
!
  return
  end
