module m_precondition
!
!  Preconditioning data for iterative charge solution
!
!   7/21 maxprecon and maxpreconu separated
!
  use datatypes
!
  implicit none
!
  logical,                              save :: luseprecon = .false.      ! If true then use full inverse for preconditioning
  logical,                              save :: lpreconsaved = .false.    ! Indicates whether the preconditioning matrix has been saved yet for use
  integer(i4),                          save :: nprecon = 0               ! Dimension of current precondition matrix
  integer(i4),                          save :: npreconu = 0              ! Upper dimension of current precondition matrix
  integer(i4),                          save :: maxprecon = 0             ! Maximum dimension of precondition matrix
  integer(i4),                          save :: maxpreconu = 0            ! Maximum upper dimension of precondition matrix
  real(dp),    dimension(:,:), pointer, save :: precon_matrix => null()   ! Preconditioning matrix

CONTAINS

  subroutine changemaxprecon
!
!  Changes the size of arrays that hold the preconditioning matrix
!
!  Julian Gale, CIC, Curtin University, June 2021
!
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(precon_matrix,maxprecon,maxpreconu,ierror)
  if (ierror.ne.0) call outofmemory('changemaxprecon','precon_matrix')
!
  end subroutine changemaxprecon

end module m_precondition
