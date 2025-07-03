!*********************************************************************
!  Serial or parallel matrix solve of serial or block cyclic matrix  *
!*********************************************************************
  subroutine matrix_solve(lsymmetric,n,ldm,nblock,matrix,x,nproc0,ifail)
!
!  Solves the system of equations Ax = b using Scalapack/Blacs
!  for a matrix in parallel or Lapack in serial.
!
!  NB: In parallel x only needs to be defined on nproc0 on input, but will
!      be returned on all processes
!
!   6/21 Created
!   7/21 Parallel version modified
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
  use iochannels
  use parallel
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                      :: lsymmetric    ! Is this a symmetric matrix?
  integer(i4), intent(in)                      :: n             ! Size of matrix
  integer(i4), intent(in)                      :: ldm           ! First dimension of matrix
  integer(i4), intent(in)                      :: nblock        ! Block size if in parallel
  real(dp),    intent(inout)                   :: matrix(ldm,*) ! Symmetric matrix to be solved
  real(dp),    intent(inout)                   :: x(*)          ! On input = RHS: On output = solution 
  integer(i4), intent(in)                      :: nproc0        ! Processor that has the first element of matrix
  integer(i4), intent(out)                     :: ifail         ! Flag indicating status of call
!
!  Local variables in Scalapack/Blacs integer precision
!
#ifdef MPI
  integer                                      :: ifails
  integer                                      :: idesm(9)
  integer                                      :: idesx(9)
  integer                                      :: ldms
  integer                                      :: lwrk
  integer                                      :: nb
  integer                                      :: np
  integer                                      :: np0
  integer                                      :: ns
  integer                                      :: nrhs
  integer                                      :: one
#else
  integer(i4)                                  :: lwrk
#endif
  integer,     dimension(:), allocatable       :: ipivot
!
!  Local variables in GULP precision
!
  integer(i4)                                  :: ilaenv
  integer(i4)                                  :: status
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: t1p
  real(dp)                                     :: t2p
  real(dp),    dimension(:), allocatable       :: wrk
#ifdef TRACE
  call trace_in('matrix_solve')
#endif
!
  t1p = g_cpu_time()
!
!  Initialise status flag
!
  ifail = 0
!
  if (nprocs.gt.1) then
#ifdef MPI
!***************************************
!  Parallel inversion using scalapack  *
!***************************************
    allocate(ipivot(4*n),stat=status)
    if (status/=0) call outofmemory('matrix_solve','ipivot')
!
!  Set local block size
!
    nb = nblock
    np = nprocs
    np0 = nproc0
    ifails = 0
    ldms = ldm
    ns = n
    nrhs = 1
    one = 1
!
!  Set up Blacs descriptor for the matrix
!
    call descinit( idesm, ns, ns, nb, nb, 0, np0, iBlacsContext, ldms, ifails )
!
    ifail = ifails
!
!  Check for failure
!
    if (ifail.ne.0) then
      call outerror('Blacs has failed to initialise a descriptor',0_i4)
      call stopnow('matrix_solve')
    endif
!
!  Set up Blacs descriptor for the solution
!
    call descinit( idesx, ns, one, nb, nb, 0, np0, iBlacsContext, ldms, ifails )
!
    ifail = ifails
!
!  Check for failure
!
    if (ifail.ne.0) then
      call outerror('Blacs has failed to initialise a descriptor',0_i4)
      call stopnow('matrix_solve')
    endif
!
!  Parallel solve for matrix
!
    call pdgesv(ns,nrhs,matrix,one,one,idesm,ipivot,x,one,one,idesx,ifails)  
    ifail = ifails
    if (ifail.ne.0) then
      call outerror('matrix solve failed in matrix_solve',0_i4)
      call stopnow('matrix_solve')
    endif
!
!  Global broadcast of solution
!
    call sendall(x,n,nproc0,"matrix_solve","x")
!
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('matrix_solve','ipivot')
#else
    call outerror('Parallel matrix solve called without MPI',0_i4)
    call stopnow('matrix_solve')
#endif
  else
!**********************************
!  Serial inversion using lapack  *
!**********************************
!
!  Allocate workspace for solve
!
    if (lsymmetric) then
      lwrk = n*ilaenv(1_i4,'DSYTRF','U',n,-1_i4,-1_i4,-1_i4)
    else
      lwrk = n*ilaenv(1_i4,'DGETRF','U',n,-1_i4,-1_i4,-1_i4)
    endif
!
    allocate(ipivot(n),stat=status)
    if (status/=0) call outofmemory('matrix_solve','ipivot')
    allocate(wrk(lwrk),stat=status)
    if (status/=0) call outofmemory('matrix_solve','wrk')
!     
!  Call lapack to solve for symmetric matrix
!     
    if (lsymmetric) then
      call dsysv('U',n,1_i4,matrix,ldm,ipivot,x,n,wrk,lwrk,ifail)
    else
      call dgesv('U',n,1_i4,matrix,ldm,ipivot,x,n,wrk,lwrk,ifail)
    endif
    if (ifail.ne.0) then
      call outerror('matrix solve failed in matrix_solve',0_i4)
      call stopnow('matrix_solve')
    endif
!
!  Free workspace
!
    deallocate(wrk,stat=status)
    if (status/=0) call deallocate_error('matrix_solve','wrk')
    deallocate(ipivot,stat=status)
    if (status/=0) call deallocate_error('matrix_solve','ipivot')
  endif
!
  t2p = g_cpu_time()
!
#ifdef TRACE
  call trace_out('matrix_solve')
#endif
  return
  end
