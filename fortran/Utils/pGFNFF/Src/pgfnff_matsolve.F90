  subroutine pgfnff_matsolve(numat,nfrag,ntot,maxA,A,x,q,qfrag)
!
!  10/21 Created
!
!  Julian Gale, Curtin University
!
  use m_pgfnff_types
  use m_io
!
  implicit none
  integer(i4),       intent(in)    :: numat                ! Number of atoms
  integer(i4),       intent(in)    :: nfrag                ! Number of fragments
  integer(i4),       intent(in)    :: ntot                 ! Total size of problem 
  integer(i4),       intent(in)    :: maxA                 ! Left-hand dimension of A
  real(dp),          intent(inout) :: A(maxA,ntot)         ! Matrix
  real(dp),          intent(out)   :: x(ntot)              ! RHS
  real(dp),          intent(out)   :: q(numat)             ! Charges
  real(dp),          intent(in)    :: qfrag(numat)         ! Fragment charges
!
!  Local variables
!
  integer(i4)                    :: ilaenv
  integer(i4)                    :: info
  integer(i4), allocatable, save :: ipiv(:)
  integer(i4)                    :: nb
  integer(i4)                    :: lwrk
  integer(i4)                    :: ierror
  real(dp),    allocatable, save :: wrk(:)
!*****************
!  Matrix solve  *
!*****************
!
!  Set work space for solve
!
  nb = ilaenv(1_i4,'DSYTRF','U',ntot,-1_i4,-1_i4,-1_i4)
  lwrk = nb*ntot
  allocate(wrk(lwrk),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : wrk '')')
    stop
  endif
  allocate(ipiv(ntot),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ipiv '')')
    stop
  endif
!
!  Solve for charges using lapack
!
  call dsytrf('U',ntot,A,maxA,ipiv,wrk,lwrk,info)
!
!  Check info for success
!
  if (info.ne.0) then
    call pgfnff_error('solve for charges failed in pgfnff_matsolve',0_i4)
    stop
  endif
  call dsytrs('U',ntot,1_i4,A,maxA,ipiv,x,ntot,info)
!
!  Check info for success
!
  if (info.ne.0) then
    call pgfnff_error('solve for charges failed in pgfnff_matsolve',0_i4)
    stop
  endif
!
  q(1:numat) = x(1:numat)
  if (numat.eq.1) q(1) = qfrag(1)
!
  deallocate(wrk,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : wrk '')')
    stop
  endif
  deallocate(ipiv,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ipiv '')')
    stop
  endif

  end subroutine pgfnff_matsolve

