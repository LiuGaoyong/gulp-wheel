  subroutine changemaxvalbond
!
!  Alters the size of the arrays associated with maxvalbond
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
  use bondvalence
  use reallocate
  implicit none
!
  integer(i4)       :: ierror
!
  call realloc(nVBspecB1,maxvalbond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbond','nVBspecB1')
  call realloc(nVBspecB2,maxvalbond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbond','nVBspecB2')
  call realloc(nVBtypeB1,maxvalbond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbond','nVBtypeB1')
  call realloc(nVBtypeB2,maxvalbond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbond','nVBtypeB2')
  call realloc(nvalbondtype,maxvalbond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbond','nvalbondtype')
  call realloc(VBparB,maxVBparB,maxvalbond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbond','VBparB')
  call realloc(VBwgtB,maxvalbond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbond','VBwgtB')
  call realloc(rVBmax,maxvalbond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbond','rVBmax')
  call realloc(rVBmin,maxvalbond,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbond','rVBmin')
!
  return
  end
