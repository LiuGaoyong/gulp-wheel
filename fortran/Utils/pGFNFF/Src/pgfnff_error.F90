  subroutine pgfnff_error(errorstring,iline)
!
!  Outputs error message in standard form
!
!  Julian Gale, CIC, Curtin University, July 2006
!
  use m_pgfnff_types
  use m_io
  implicit none
!
!  Passed variables
!
  character(len=*), intent(in) :: errorstring
  integer(i4),      intent(in) :: iline
!******************
!  Output header  *
!******************
  write(ioout,'(/)')
  write(ioout,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
  write(ioout,'(''!! ERROR : '',a)') trim(errorstring)
  if (iline.gt.0.and.iline.lt.1000000) then
    write(ioout,'(''!!       : Error is apparently on line '',i6)') iline
  endif
  write(ioout,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'',/)')
!
  return
  end
