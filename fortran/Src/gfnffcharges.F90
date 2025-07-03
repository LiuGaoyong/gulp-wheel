  subroutine gfnffcharges
!
!  Outputs the charges for GFNFF which are already computed
!
!   6/23 Created from reaxffcharges.f90
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
!  Julian Gale, CIC, Curtin University, June 2023
!
  use datatypes
  use current
  use iochannels
  use parallel 
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
#ifdef TRACE
  call trace_in('gfnffcharges')
#endif
!*******************
!  Output results  *
!*******************
  if (ioproc) then
    write(ioout,'(//,''  Final charges from pGFNFF :'',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nasym
      write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
#ifdef TRACE
  call trace_out('gfnffcharges')
#endif
!
  return
  end
