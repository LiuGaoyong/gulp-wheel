  subroutine pgfnff_warn(warningstring,iline)
!
!  Outputs warning message in standard form
!
!  Julian Gale, CIC, Curtin University, August 2011
!
  use m_pgfnff_types
  use m_io
  implicit none
!
!  Passed variables
!
  character(len=*), intent(in) :: warningstring
  integer(i4),      intent(in) :: iline
!******************
!  Output header  *
!******************
  write(ioout,'(/)')
  write(ioout,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
  write(ioout,'(''!! WARNING : '',a)') trim(warningstring)
  if (iline.gt.0.and.iline.lt.1000000) then
    write(ioout,'(''!!         : Warning is apparently on line '',i6)') iline
  endif
  write(ioout,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'',/)')
!
  return
  end
