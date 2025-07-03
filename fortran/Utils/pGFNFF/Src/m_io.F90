!**************
!  I/O Module *
!**************
!
!  Defines the I/O channels
!
!  Julian Gale, CIC, Curtin University, August 2013
!
  module m_io
    use m_pgfnff_types
    integer(i4),                    save :: ioin  = 5_i4  ! Input channel
    integer(i4),                    save :: ioout = 6_i4  ! Output channel
  end module m_io
