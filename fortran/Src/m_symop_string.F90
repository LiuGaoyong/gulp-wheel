module m_symop_string
contains
  subroutine as_fraction(real_number, numerator, denominator, denominator_limit)
!
!  Convert a real number (input) to a good approximation as a rational number
!  expressed as numerator / denominator
!
!  Peter Spackman, Curtin University, January 2021
!
  use configurations
  implicit none
!
!  Passed variables
!
  real(dp),    intent(in)  :: real_number
  integer(i4), intent(in)  :: denominator_limit
  integer(i4), intent(out) :: numerator, denominator
!
!  Local variables
!
  integer(i4), dimension(3) :: h, k
  integer(i4)               :: a, d, i, n, x
  logical                   :: neg
  real(dp)                  :: f
!
  h = [0, 1, 0]
  k = [1, 0, 0]
  x = 1
  d = 1
  n = 1
  neg = .false.
!
  f = real_number
!
  if (denominator_limit .le. 1) then
    numerator = 0
    denominator = 1
    return
  endif
!
  if (real_number .lt. 0) then
    neg = .true.
    f = - real_number
  endif
!
  do while (f .ne. floor(f)) 
    n = n * 2
    f = f * 2
  enddo
  d = f
!
  do i = 1, 64
    if (n .ne. 0) then
      a = d / n
    else
      a = 0
    endif

    if ((i .gt. 1) .and. (a .eq. 0)) then
      exit
    endif
!
    x = d
    d = n
    n = modulo(x, n)
    x = a
!
    if ((k(2) * a + k(1)) .ge. denominator_limit) then
      x = (denominator_limit - k(1)) / k(2)
      if (((x * 2) .ge. a) .or. (k(2) .ge. denominator_limit)) then
        h(3) = x * h(2) + h(1)
        h(1) = h(2)
        h(2) = h(3)
        k(3) = x * k(2) + k(1)
        k(1) = k(2)
        k(2) = k(3)
      endif
      exit
    endif
!
    h(3) = x * h(2) + h(1)
    h(1) = h(2)
    h(2) = h(3)
    k(3) = x * k(2) + k(1)
    k(1) = k(2)
    k(2) = k(3)
  enddo

  denominator = k(2)
  numerator = h(2)
  if (neg) then
    numerator = - numerator
  endif

  end subroutine
!
  function operation_as_xyz_string(rot, trans) result(str)
!
!  Convert a symmetry operation expressed as a rotation matrix + translation vector
!  to an x,y,z formatted string for output into CIF or other crystallography
!  formats
!
!  Peter Spackman, Curtin University, January 2021
!
  use configurations
  implicit none
!
!  Passed variables
!
  real(dp), dimension(3,3), intent(in) :: rot
  real(dp), dimension(3),   intent(in) :: trans
  character(len=32)                    :: str
!
!  Local variables
!
  character(len=3)                     :: tmp_str
  character, dimension(3)              :: symbols
  integer(i4)                          :: i, j, denominator, numerator
  real(dp)                             :: t, r
!
  symbols = ['x', 'y', 'z']
!
  str = ""
  do i = 1,3
    t = trans(i)
    if (t .ne. 0.0) then
      call as_fraction(t, numerator, denominator, 12)
      write(tmp_str, '(i0,a,i0)') numerator, '/', denominator
      str = trim(str)//trim(tmp_str)
    endif
    do j = 1, 3
      r = rot(i, j)
      if (r .gt. 0.0) then
        str = trim(str)//"+"
        str = trim(str)//symbols(j)
      elseif (r .lt. 0.0) then
        str = trim(str)//"-"
        str = trim(str)//symbols(j)
      endif
    enddo
    if (i .lt. 3) str = trim(str)// ","
  end do
  end function
end module m_symop_string
