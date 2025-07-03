  subroutine dbcgsas(n,nloc,nlocptr,A,maxA,b,x,tol,itmax,iter,err)
!
!  Dense version of bcg solve with special format for charge on SAS
!
!   6/21 Created from bcgsolve
!
  use datatypes
  use iochannels
  use control,       only : keyword
  use parallel,      only : ioproc
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif

  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: n                  ! Total dimension of problem (excluding constraint)
  integer(i4), intent(in)    :: nloc               ! Number of columns of A local to this node
  integer(i4), intent(in)    :: nlocptr(nloc)      ! Pointer from local columns to global number
  integer(i4), intent(in)    :: itmax              ! Maximum number of iterations
  integer(i4), intent(out)   :: iter               ! Number of iterations taken
  integer(i4), intent(in)    :: maxA               ! LHS dimension of listA and A
  real(dp),    intent(in)    :: A(maxA,nloc)       ! Values of non-zero elements of column of A
  real(dp),    intent(in)    :: tol                ! Tolerance
  real(dp),    intent(in)    :: b(n+1)             ! Left-hand vector
  real(dp),    intent(inout) :: x(n+1)             ! Solution vector : guess on input, actual on output
  real(dp),    intent(out)   :: err                ! Error flag
!
!  Local variables
!
  integer(i4)                :: j
  integer(i4)                :: nm
  logical                    :: converged
  real(dp)                   :: ak
  real(dp)                   :: akden
  real(dp)                   :: bk
  real(dp)                   :: bkden
  real(dp)                   :: bknum
  real(dp)                   :: bnrm
  real(dp)                   :: err1(1)
  real(dp)                   :: ddot
  real(dp),    allocatable   :: p(:)
  real(dp),    allocatable   :: pp(:)
  real(dp),    allocatable   :: r(:)
  real(dp),    allocatable   :: rr(:)
  real(dp),    allocatable   :: z(:)
  real(dp),    allocatable   :: zz(:)
#ifdef TRACE
  call trace_in('dbcgsas')
#endif
!
  nm = n + 1
!
!  Allocate workspace
!
  allocate(p(nm))
  allocate(pp(nm))
  allocate(r(nm))
  allocate(rr(nm))
  allocate(z(nm))
  allocate(zz(nm))

  call denseAxVsas(n,nloc,nlocptr,A,maxA,x,r,.false.)
!
  do j = 1,nm
    r(j)  = b(j) - r(j)
    rr(j) = r(j)
  enddo
!
  bnrm = ddot(nm,b,1_i4,b,1_i4)
  call denseAdiagpreconsas(n,nloc,nlocptr,A,maxA,r,z)
!
!  Main loop
!
  iter = 0
  converged = .false.
  do while (iter.le.itmax.and..not.converged)
    iter = iter + 1
    call denseAdiagpreconsas(n,nloc,nlocptr,A,maxA,rr,zz)
    bknum = ddot(nm,z,1_i4,rr,1_i4)
    if (iter.eq.1) then
      do j = 1,nm
        p(j) = z(j)
        pp(j) = zz(j)
      enddo
    else
      bk = bknum/bkden
      do j = 1,nm
        p(j) = bk*p(j) + z(j)
        pp(j) = bk*pp(j) + zz(j)
      enddo
    endif
    bkden = bknum
    call denseAxVsas(n,nloc,nlocptr,A,maxA,p,z,.false.)
    akden = ddot(nm,z,1_i4,pp,1_i4)
    ak = bknum/akden
    call denseAxVsas(n,nloc,nlocptr,A,maxA,pp,zz,.true.)
!
    call daxpy(nm,ak,p,1_i4,x,1_i4)
    call daxpy(nm,-ak,z,1_i4,r,1_i4)
    call daxpy(nm,-ak,zz,1_i4,rr,1_i4)
!
    call denseAdiagpreconsas(n,nloc,nlocptr,A,maxA,r,z)
!
    err = ddot(nm,r,1_i4,r,1_i4)/bnrm
    err = sqrt(err)
!
!  Send error value from a single node to ensure consistent behaviour
!
    err1(1) = err
    call sendall(err1,1_i4,0_i4,"dbcgsas","err")
    err = err1(1)
!
    if (ioproc.and.index(keyword,'verb').ne.0) then
      write(ioout,'('' BCG iteration = '',i4,'' Error = '',f16.12)') iter,err
    endif
    converged = (err.lt.tol)
  enddo
!
!  Free workspace
!
  deallocate(zz)
  deallocate(z)
  deallocate(rr)
  deallocate(r)
  deallocate(pp)
  deallocate(p)
#ifdef TRACE
  call trace_out('dbcgsas')
#endif

  end subroutine dbcgsas

  subroutine denseAdiagpreconsas(n,nloc,nlocptr,A,maxA,Vin,Vout)
!
!  Perform the precondition of a vector by the diagonal of a dense matrix
!
!  Code is specific to the purpose here in that the last row and column 
!  are assumed to be 1, except the diagonal element which is zero. 
!  This avoids having to explicitly store the constraint terms.
!
  use datatypes
  use parallel,      only : nprocs
  use times,         only : tsum
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
  integer(i4),  intent(in)    :: n
  integer(i4),  intent(in)    :: nloc
  integer(i4),  intent(in)    :: nlocptr(nloc)
  integer(i4),  intent(in)    :: maxA
  real(dp),     intent(in)    :: A(maxA,nloc)
  real(dp),     intent(in)    :: Vin(n+1)
  real(dp),     intent(out)   :: Vout(n+1)
!
  integer(i4)                 :: i
  integer(i4)                 :: io
  integer(i4)                 :: nm
  real(dp),     allocatable   :: Vloc(:)
  real(dp)                    :: g_cpu_time
  real(dp)                    :: tsuml
#ifdef TRACE
  call trace_in('denseAdiagpreconsas')
#endif
!
  nm = n + 1
!
!  This code is specific to the form of the matrix being used here!!
!
  Vout(1:nm) = 0.0_dp
  do i = 1,nloc
    io = nlocptr(i)
    if (io.le.n) then
      Vout(io) = Vin(io)/A(io,i)
    else
      Vout(io) = Vin(io)
    endif
  enddo
  if (nprocs.gt.1) then
!
!  Globalization of Vout
!
    tsuml = g_cpu_time()
    allocate(Vloc(nm))
    call sumall(Vout,Vloc,nm,"denseAdiagpreconsas","Vout") 
    Vout(1:nm) = Vloc(1:nm)
    deallocate(Vloc)
    tsum = tsum + g_cpu_time() - tsuml
  endif
#ifdef TRACE
  call trace_out('denseAdiagpreconsas')
#endif

  end subroutine denseAdiagpreconsas

  subroutine denseAxVsas(n,nloc,nlocptr,A,maxA,Vin,Vout,ltranspose)
!
!  Perform the multiplication of a dense matrix by a vector
!
!  Code is specific to the purpose here in that the last row and column 
!  are assumed to be 1, except the diagonal element which is zero. 
!  This avoids having to explicitly store the constraint terms. 
!
  use datatypes
  use parallel,      only : nprocs, procid
  use times,         only : tsum
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
  integer(i4),  intent(in)    :: n
  integer(i4),  intent(in)    :: nloc
  integer(i4),  intent(in)    :: nlocptr(nloc)
  integer(i4),  intent(in)    :: maxA
  logical,      intent(in)    :: ltranspose
  real(dp),     intent(in)    :: A(maxA,*)
  real(dp),     intent(in)    :: Vin(n+1)
  real(dp),     intent(out)   :: Vout(n+1)

  integer(i4)                 :: i
  integer(i4)                 :: il
  integer(i4)                 :: j
  real(dp),     allocatable   :: Vloc(:)
!
  real(dp)                    :: g_cpu_time
  real(dp)                    :: ddot
  real(dp)                    :: tsuml
#ifdef TRACE
  call trace_in('denseAxVsas')
#endif
!
  if (ltranspose) then
    if (nprocs.eq.1) then
      do i = 1,n
        Vout(i) = ddot(n,A(i,1),maxA,Vin,1_i4)
      enddo
    else
      Vout(1:n) = 0.0_dp
      do il = 1,nloc
        i = nlocptr(il)
        do j = 1,n
          Vout(j) = Vout(j) + A(j,il)*Vin(i)
        enddo
      enddo
      if (procid.eq.0) then
        do j = 1,n
          Vout(j) = Vout(j) + A(j,nloc+1)*Vin(n)
        enddo
      endif
    endif
  else
    if (nprocs.eq.1) then
      do i = 1,n
        Vout(i) = ddot(n,A(1,i),1_i4,Vin,1_i4)
      enddo
    else
      Vout(1:n) = 0.0_dp
      do il = 1,nloc
        i = nlocptr(il)
        Vout(i) = ddot(n,A(1,il),1_i4,Vin,1_i4)
      enddo
      if (procid.eq.0) then
        Vout(n) = ddot(n,A(1,nloc+1),1_i4,Vin,1_i4)
      endif
    endif
  endif
  if (nprocs.gt.1) then
    allocate(Vloc(n))
!
!  Globalization of Vout
!
    tsuml = g_cpu_time()
    call sumall(Vout,Vloc,n,"denseAxV","Vout") 
    Vout(1:n) = Vloc(1:n)
    tsum = tsum + g_cpu_time() - tsuml
    deallocate(Vloc)
  endif
#ifdef TRACE
  call trace_out('denseAxVsas')
#endif

  end subroutine denseAxVsas
!
  subroutine dcgsolve(n,nloc,nlocptr,A,maxA,b,x,tol,itmax,iter,err)
!
!  Dense version of cgsolve
!
!   4/21 m argument added
!
  use datatypes
  use iochannels
  use control,       only : keyword
  use parallel,      only : ioproc
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif

  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: n                  ! Total dimension of problem
  integer(i4), intent(in)    :: nloc               ! Number of columns of A local to this node
  integer(i4), intent(in)    :: nlocptr(nloc)      ! Pointer from local columns to global number
  integer(i4), intent(in)    :: itmax              ! Maximum number of iterations
  integer(i4), intent(out)   :: iter               ! Number of iterations taken
  integer(i4), intent(in)    :: maxA               ! LHS dimension of listA and A
  real(dp),    intent(in)    :: A(maxA,nloc)       ! Values of non-zero elements of column of A
  real(dp),    intent(in)    :: tol                ! Tolerance
  real(dp),    intent(in)    :: b(n)               ! Left-hand vector
  real(dp),    intent(inout) :: x(n)               ! Solution vector : guess on input, actual on output
  real(dp),    intent(out)   :: err                ! Error flag
!
!  Local variables
!
  integer(i4)                :: j
  logical                    :: converged
  real(dp)                   :: ak
  real(dp)                   :: akden
  real(dp)                   :: bk
  real(dp)                   :: bkden
  real(dp)                   :: bknum
  real(dp)                   :: bnrm
  real(dp)                   :: err1(1)
  real(dp)                   :: ddot
  real(dp),    allocatable   :: p(:)
  real(dp),    allocatable   :: r(:)
  real(dp),    allocatable   :: z(:)
#ifdef TRACE
  call trace_in('dcgsolve')
#endif
!
!  Allocate workspace
!
  allocate(p(n))
  allocate(r(n))
  allocate(z(n))
!
  call denseAxV(n,nloc,nlocptr,A,maxA,x,r,.false.)
!
  do j = 1,n
    r(j)  = b(j) - r(j)
  enddo
!
  bnrm = ddot(n,b,1_i4,b,1_i4)
  call denseAdiagprecon(n,nloc,nlocptr,A,maxA,r,z)
!
!  Main loop
!
  iter = 0
  converged = .false.
  do while (iter.le.itmax.and..not.converged)
    iter = iter + 1
    bknum = ddot(n,z,1_i4,r,1_i4)
    if (iter.eq.1) then
      do j = 1,n
        p(j) = z(j)
      enddo
    else
      bk = bknum/bkden
      do j = 1,n
        p(j) = bk*p(j) + z(j)
      enddo
    endif
    bkden = bknum
    call denseAxV(n,nloc,nlocptr,A,maxA,p,z,.false.)
    akden = ddot(n,z,1_i4,p,1_i4)
    ak = bknum/akden
!
    call daxpy(n,ak,p,1_i4,x,1_i4)
    call daxpy(n,-ak,z,1_i4,r,1_i4)
!
    call denseAdiagprecon(n,nloc,nlocptr,A,maxA,r,z)
!
    err = ddot(n,r,1_i4,r,1_i4)/bnrm
    err = sqrt(err)
!
!  Send error value from a single node to ensure consistent behaviour
!
    err1(1) = err
    call sendall(err1,1_i4,0_i4,"dcgsolve","err")
    err = err1(1)
!
    if (ioproc.and.index(keyword,'verb').ne.0) then
      write(ioout,'('' CG iteration = '',i4,'' Error = '',f16.12)') iter,err
    endif
    converged = (err.lt.tol)
  enddo
!
!  Free workspace
!
  deallocate(z)
  deallocate(r)
  deallocate(p)
#ifdef TRACE
  call trace_out('dcgsolve')
#endif

  end subroutine dcgsolve

  subroutine dbcgsolve(n,nloc,nlocptr,A,maxA,b,x,tol,itmax,iter,err)
!
!  Dense version of bcgsolve
!
!   4/21 m argument added
!
  use datatypes
  use iochannels
  use control,       only : keyword
  use parallel,      only : ioproc
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif

  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: n                  ! Total dimension of problem 
  integer(i4), intent(in)    :: nloc               ! Number of columns of A local to this node
  integer(i4), intent(in)    :: nlocptr(nloc)      ! Pointer from local columns to global number
  integer(i4), intent(in)    :: itmax              ! Maximum number of iterations
  integer(i4), intent(out)   :: iter               ! Number of iterations taken
  integer(i4), intent(in)    :: maxA               ! LHS dimension of listA and A
  real(dp),    intent(in)    :: A(maxA,nloc)       ! Values of non-zero elements of column of A
  real(dp),    intent(in)    :: tol                ! Tolerance
  real(dp),    intent(in)    :: b(n)               ! Left-hand vector
  real(dp),    intent(inout) :: x(n)               ! Solution vector : guess on input, actual on output
  real(dp),    intent(out)   :: err                ! Error flag
!
!  Local variables
!
  integer(i4)                :: j
  logical                    :: converged
  real(dp)                   :: ak
  real(dp)                   :: akden
  real(dp)                   :: bk
  real(dp)                   :: bkden
  real(dp)                   :: bknum
  real(dp)                   :: bnrm
  real(dp)                   :: err1(1)
  real(dp)                   :: ddot
  real(dp),    allocatable   :: p(:)
  real(dp),    allocatable   :: pp(:)
  real(dp),    allocatable   :: r(:)
  real(dp),    allocatable   :: rr(:)
  real(dp),    allocatable   :: z(:)
  real(dp),    allocatable   :: zz(:)
#ifdef TRACE
  call trace_in('dbcgsolve')
#endif
!
!  Allocate workspace
!
  allocate(p(n))
  allocate(pp(n))
  allocate(r(n))
  allocate(rr(n))
  allocate(z(n))
  allocate(zz(n))

  call denseAxV(n,nloc,nlocptr,A,maxA,x,r,.false.)
!
  do j = 1,n
    r(j)  = b(j) - r(j)
    rr(j) = r(j)
  enddo
!
  bnrm = ddot(n,b,1_i4,b,1_i4)
  call denseAdiagprecon(n,nloc,nlocptr,A,maxA,r,z)
!
!  Main loop
!
  iter = 0
  converged = .false.
  do while (iter.le.itmax.and..not.converged)
    iter = iter + 1
    call denseAdiagprecon(n,nloc,nlocptr,A,maxA,rr,zz)
    bknum = ddot(n,z,1_i4,rr,1_i4)
    if (iter.eq.1) then
      do j = 1,n
        p(j) = z(j)
        pp(j) = zz(j)
      enddo
    else
      bk = bknum/bkden
      do j = 1,n
        p(j) = bk*p(j) + z(j)
        pp(j) = bk*pp(j) + zz(j)
      enddo
    endif
    bkden = bknum
    call denseAxV(n,nloc,nlocptr,A,maxA,p,z,.false.)
    akden = ddot(n,z,1_i4,pp,1_i4)
    ak = bknum/akden
    call denseAxV(n,nloc,nlocptr,A,maxA,pp,zz,.true.)
!
    call daxpy(n,ak,p,1_i4,x,1_i4)
    call daxpy(n,-ak,z,1_i4,r,1_i4)
    call daxpy(n,-ak,zz,1_i4,rr,1_i4)
!
    call denseAdiagprecon(n,nloc,nlocptr,A,maxA,r,z)
!
    err = ddot(n,r,1_i4,r,1_i4)/bnrm
    err = sqrt(err)
!
!  Send error value from a single node to ensure consistent behaviour
!
    err1(1) = err
    call sendall(err1,1_i4,0_i4,"dbcgsolve","err")
    err = err1(1)
!
    if (ioproc.and.index(keyword,'verb').ne.0) then
      write(ioout,'('' BCG iteration = '',i4,'' Error = '',f16.12)') iter,err
    endif
    converged = (err.lt.tol)
  enddo
!
!  Free workspace
!
  deallocate(zz)
  deallocate(z)
  deallocate(rr)
  deallocate(r)
  deallocate(pp)
  deallocate(p)
#ifdef TRACE
  call trace_out('dbcgsolve')
#endif

  end subroutine dbcgsolve
!
  subroutine denseAdiagprecon(n,nloc,nlocptr,A,maxA,Vin,Vout)
!
!  Perform the precondition of a vector by the diagonal of a dense matrix
!
  use datatypes
  use m_precondition
  use parallel,      only : nprocs
  use times,         only : tsum
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
  integer(i4),  intent(in)    :: n
  integer(i4),  intent(in)    :: nloc
  integer(i4),  intent(in)    :: nlocptr(nloc)
  integer(i4),  intent(in)    :: maxA
  real(dp),     intent(in)    :: A(maxA,nloc)
  real(dp),     intent(in)    :: Vin(n)
  real(dp),     intent(out)   :: Vout(n)
!
  integer(i4)                 :: i
  integer(i4)                 :: io
  real(dp),     allocatable   :: Vloc(:)
  real(dp)                    :: g_cpu_time
  real(dp)                    :: tsuml
#ifdef TRACE
  call trace_in('denseAdiagprecon')
#endif
!
  Vout(1:n) = 0.0_dp
!
  if (luseprecon.and.lpreconsaved) then
    call denseAxV(n,nloc,nlocptr,precon_matrix,maxprecon,Vin,Vout,.false.)
  else
    do i = 1,nloc
      io = nlocptr(i)
      if (abs(A(io,i)).gt.1.0d-6) then
        Vout(io) = Vin(io)/A(io,i)
      else
        Vout(io) = Vin(io)
      endif
    enddo
!
    if (nprocs.gt.1) then
!
!  Globalization of Vout
!
      tsuml = g_cpu_time()
      allocate(Vloc(n))
      call sumall(Vout,Vloc,n,"denseAdiagprecon","Vout") 
      Vout(1:n) = Vloc(1:n)
      deallocate(Vloc)
      tsum = tsum + g_cpu_time() - tsuml
    endif
  endif
#ifdef TRACE
  call trace_out('denseAdiagprecon')
#endif

  end subroutine denseAdiagprecon

  subroutine denseAxV(n,nloc,nlocptr,A,maxA,Vin,Vout,ltranspose)
!
!  Perform the multiplication of a dense matrix by a vector
!
  use datatypes
  use parallel,      only : nprocs
  use times,         only : tsum
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
  integer(i4),  intent(in)    :: n
  integer(i4),  intent(in)    :: nloc
  integer(i4),  intent(in)    :: nlocptr(nloc)
  integer(i4),  intent(in)    :: maxA
  logical,      intent(in)    :: ltranspose
  real(dp),     intent(in)    :: A(maxA,*)
  real(dp),     intent(in)    :: Vin(n)
  real(dp),     intent(out)   :: Vout(n)

  integer(i4)                 :: i
  integer(i4)                 :: il
  integer(i4)                 :: j
  real(dp),     allocatable   :: Vloc(:)
!
  real(dp)                    :: g_cpu_time
  real(dp)                    :: ddot
  real(dp)                    :: tsuml
#ifdef TRACE
  call trace_in('denseAxV')
#endif
!
  if (ltranspose) then
    if (nprocs.eq.1) then
      do i = 1,n
        Vout(i) = ddot(n,A(i,1),maxA,Vin,1_i4)
      enddo
    else
      Vout(1:n) = 0.0_dp
      do il = 1,nloc
        i = nlocptr(il)
        do j = 1,n
          Vout(j) = Vout(j) + A(j,il)*Vin(i)
        enddo
      enddo
    endif
  else
    if (nprocs.eq.1) then
      do i = 1,n
        Vout(i) = ddot(n,A(1,i),1_i4,Vin,1_i4)
      enddo
    else
      Vout(1:n) = 0.0_dp
      do il = 1,nloc
        i = nlocptr(il)
        Vout(i) = ddot(n,A(1,il),1_i4,Vin,1_i4)
      enddo
    endif
  endif
  if (nprocs.gt.1) then
    allocate(Vloc(n))
!
!  Globalization of Vout
!
    tsuml = g_cpu_time()
    call sumall(Vout,Vloc,n,"denseAxV","Vout") 
    Vout(1:n) = Vloc(1:n)
    tsum = tsum + g_cpu_time() - tsuml
    deallocate(Vloc)
  endif
#ifdef TRACE
  call trace_out('denseAxV')
#endif

  end subroutine denseAxV
