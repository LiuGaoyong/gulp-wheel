  module m_pgfnff_types
!
!  i4 and i2 define the size of integer*4 and integer*2 in the code
!
    use, intrinsic :: iso_c_binding
!
!  Integers
!
    integer, parameter :: i4  = selected_int_kind(5)
    integer, parameter :: i2  = selected_int_kind(3)
!
!  Floats
!
    integer, parameter :: dp  = kind(1.0d0)
    integer, parameter :: dpc = kind((1.0d0,1.0d0))

  end module m_pgfnff_types
