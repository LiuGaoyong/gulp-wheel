  subroutine pdiaggd(mcv,mcvloc,maxd2,derv2,maxeigr,eigr,freq,fscale,lvectors,lprint,ifail)
!
!  Calculate eigenvalues and eigenvectors of the dynamical matrix
!  and subsequent calculation of properties of eigenvectors.
!
!  Gamma point only version with distributed memory
!
!  11/16 Created from pdiaggd
!  12/16 maxeigr added as a separate argument to maxd2
!   2/17 Corrected so that idesce is not set up when eigenvectors are not needed
!        otherwise there will be an error due to the dimensions of the array
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!   2/18 Trace added
!   7/20 Option to use 2D block cyclic form added
!   7/20 Option to use the ELPA library added
!
!  On entry :
!
!  mcv      = number of modes (global)
!  mcvloc   = number of modes (local)
!  maxd2    = left-hand dimension of derv2
!  maxeigr  = left-hand dimension of eigr
!  derv2    = mass-weighted dynamical matrix
!  fscale   = scale factor for frequencies to convert to wavenumbers
!  lvectors = if .true. then calculate eigenvectors
!  lprint   = if .true. print warnings
!
!  On exit :
!
!  eigr     = eigenvectors of dynamical matrix (if lvectors is true)
!  freq     = frequencies of vibration in wavenumbers
!  ifail    = flag indicating success or failure
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, July 2020
!
  use g_constants
  use current
  use element
#ifdef MPI
  use maths,         only : lblock2D, ldivide_and_conquer
#endif
#ifdef ELPA
  use maths,         only : lelpa
#endif
  use parallel
  use species
  use times
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(out)                    :: ifail
  integer(i4),  intent(in)                     :: maxd2
  integer(i4),  intent(in)                     :: maxeigr
  integer(i4),  intent(in)                     :: mcv
  integer(i4),  intent(in)                     :: mcvloc
  logical,      intent(in)                     :: lprint
  logical,      intent(in)                     :: lvectors
  real(dp),     intent(in)                     :: derv2(maxd2,*)
  real(dp),     intent(out)                    :: eigr(maxeigr,*)
  real(dp),     intent(in)                     :: fscale
  real(dp),     intent(out)                    :: freq(mcv)
#ifdef MPI
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: status
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: root
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp),    dimension(:,:), allocatable     :: wrk2
  real(dp),    dimension(:),   allocatable     :: wrk
!
!  Blacs / Scalapack integers
!
  integer                                      :: iBlacsContext2D
  integer                                      :: iBlacsContextSys
  integer                                      :: idescd(9)
  integer                                      :: idescd2D(9)
  integer                                      :: idesce(9)
  integer                                      :: ifails
  integer                                      :: ld
  integer                                      :: liwrk
  integer                                      :: lwrk
  integer                                      :: mcol2D
  integer                                      :: mrow2D
  integer                                      :: n
  integer                                      :: nb
  integer                                      :: nloc
  integer                                      :: npcol2D
  integer                                      :: npcol2Dloc
  integer                                      :: nprow2D
  integer                                      :: nprow2Dloc
  integer                                      :: nsize
  integer                                      :: numroc
  integer,     dimension(:),   allocatable     :: iwrk
  real*8,      dimension(:,:), allocatable     :: d2D
  real*8,      dimension(:,:), allocatable     :: e2D
#ifdef TRACE
  call trace_in('pdiaggd')
#endif
!
  t1 = g_cpu_time()
!
  n = mcv
  nloc = mcvloc
  nb = 3*nblocksize
!
#ifdef ELPA
  if (lelpa) then
    call eigensolve_elpa(mcv,3_i4*nblocksize,derv2,maxd2,eigr,maxeigr,freq,lvectors)
  elseif (lblock2D) then
#else
  if (lblock2D) then
#endif
!-------------
!  Block 2D  |
!-------------
!
!  Get system context
!
    call blacs_get(0, 0, iBlacsContextSys)
!
!  Set up Blacs descriptors for arrays in current form
!
    ld = maxd2
    call descinit(idescd, n, n, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - idescd ',0_i4)
      call stopnow('pdiaggd')
    endif
!
    if (lvectors.or.ldivide_and_conquer) then
      ld = maxeigr
      call descinit(idesce, n, n, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
      if (ifails.ne.0) then
        call outerror('initialisation in descinit failed - idesce ',0_i4)
        call stopnow('pdiaggd')
      endif
    endif
!
!  Try to set as close to a square grid as possible
!
    do npcol2D = nint(sqrt(dble(nprocs))),2,-1
      if(mod(nprocs,npcol2D) == 0 ) exit
    enddo
    nprow2D = nprocs/npcol2D
!
!  Initialise new blacs grid
!
    iBlacsContext2D = iBlacsContextSys
    call blacs_gridinit(iBlacsContext2D, 'C', nprow2D, npcol2D)
    call blacs_gridinfo(iBlacsContext2D, nprow2D, npcol2D, nprow2Dloc, npcol2Dloc)
!
!  Find local size in 2D
!
    mrow2D = numroc(n, nb, nprow2Dloc, 0, nprow2D)
    mcol2D = numroc(n, nb, npcol2Dloc, 0, npcol2D)
!
!  Set up Blacs descriptors for arrays in 2D form
!
    ld = mrow2D
    call descinit(idescd2D, n, n, nb, nb, 0, 0, iBlacsContext2D, ld, ifails)
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - idescd2D ',0_i4)
      call stopnow('pdiaggd')
    endif
!
!  Allocate arrays for 2D storage of D and E
!
    allocate(d2D(mrow2D,mcol2D),stat=status)
    if (status/=0) call outofmemory('pdiaggd','d2D')
    allocate(e2D(mrow2D,mcol2D),stat=status)
    if (status/=0) call outofmemory('pdiaggd','e2D')
!
!  Copy D to 2D form
!
    call pdgemr2d(n, n, derv2, 1, 1, idescd, d2D, 1, 1, idescd2D, iBlacsContext)
!
!  Initial allocation of workspace until memory requirements are found
!
    allocate(wrk(1),stat=status)
    if (status/=0) call outofmemory('pdiaggd','wrk')
    if (ldivide_and_conquer) then
      allocate(iwrk(1),stat=status)
      if (status/=0) call outofmemory('pdiaggd','iwrk')
    endif
!
!  Call eigensolver to find memory requirements
!
    liwrk = -1
    lwrk = -1
    if (ldivide_and_conquer) then
      call pdsyevd( 'V', 'U', n, d2D, 1, 1, idescd2D, freq, e2D, 1, 1, idescd2D, wrk, lwrk, iwrk, liwrk, ifails )
    else
      if (lvectors) then
        call pdsyev( 'V', 'U', n, d2D, 1, 1, idescd2D, freq, e2D, 1, 1, idescd2D, wrk, lwrk, ifails )
      else
        call pdsyev( 'N', 'U', n, d2D, 1, 1, idescd2D, freq, e2D, 1, 1, idescd2D, wrk, lwrk, ifails )
      endif
    endif
    lwrk = nint(wrk(1))
!
!  Re-size wrk
!
    if (ldivide_and_conquer) then
      liwrk = iwrk(1)
      deallocate(iwrk,stat=status)
      if (status/=0) call deallocate_error('pdiaggd','iwrk')
    endif
    deallocate(wrk,stat=status)
    if (status/=0) call deallocate_error('pdiaggd','wrk')
!
    allocate(wrk(lwrk),stat=status)
    if (status/=0) call outofmemory('pdiaggd','wrk')
    if (ldivide_and_conquer) then
      allocate(iwrk(liwrk),stat=status)
      if (status/=0) call outofmemory('pdiaggd','iwrk')
    endif
!
!  Call eigensolver for real
!
    if (ldivide_and_conquer) then
      call pdsyevd( 'V', 'U', n, d2D, 1, 1, idescd2D, freq, e2D, 1, 1, idescd2D, wrk, lwrk, iwrk, liwrk, ifails )
    else
      if (lvectors) then
        call pdsyev( 'V', 'U', n, d2D, 1, 1, idescd2D, freq, e2D, 1, 1, idescd2D, wrk, lwrk, ifails )
      else
        call pdsyev( 'N', 'U', n, d2D, 1, 1, idescd2D, freq, e2D, 1, 1, idescd2D, wrk, lwrk, ifails )
      endif
    endif
!
!  Copy e2D back to 1D form
!
    if (lvectors.or.ldivide_and_conquer) then
      call pdgemr2d(n, n, e2D, 1, 1, idescd2D, eigr, 1, 1, idesce, iBlacsContext)
    endif
!
!  Deallocate arrays for 2D storage of D and E
!
    deallocate(e2D,stat=status)
    if (status/=0) call deallocate_error('pdiaggd','e2D')
    deallocate(d2D,stat=status)
    if (status/=0) call deallocate_error('pdiaggd','d2D')
!
!  Shutdown blacs context
!
    call blacs_gridexit(iBlacsContext2D)
  else
!-------------
!  Block 1D  |
!-------------
!
!  Set up Blacs descriptors for arrays
!
    ld = maxd2
    call descinit(idescd, n, n, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - idescd ',0_i4)
      call stopnow('pdiaggd')
    endif
!
!  Only need eigr BLACS descriptor if vectors are to be computed
!
    if (lvectors) then
      ld = maxeigr
      call descinit(idesce, n, n, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
      if (ifails.ne.0) then
        call outerror('initialisation in descinit failed - idesce ',0_i4)
        call stopnow('pdiaggd')
      endif
    endif
!
!  Find size of work space array
!
    if (lvectors) then
      nsize = max((nb*(nb-1))/2,(n+nloc)*nb) + nb*nb
      lwrk = 5*n + n*max(1,nloc) + max(2*n,nsize) + 1
    else
      lwrk = 5*n + max(nb*(n+1),3*nb) + 1
    endif
!
!  Allocate work space
!
    allocate(wrk(lwrk),stat=status)
    if (status/=0) call outofmemory('pdiaggd','wrk')
    allocate(wrk2(maxd2,mcvloc),stat=status)
    if (status/=0) call outofmemory('pdiaggd','wrk2')
!
!  Make copy of derv2 to avoid overwriting
!
    do i = 1,mcvloc
      do j = 1,mcv
        wrk2(j,i) = derv2(j,i)
      enddo
    enddo
!
!  Call Scalapack eigensolver
!
    if (lvectors) then
      call pdsyev( 'V', 'U', n, wrk2, 1, 1, idescd, freq, eigr, 1, 1, idesce, wrk, lwrk, ifails )
    else
      call pdsyev( 'N', 'U', n, wrk2, 1, 1, idescd, freq, wrk2, 1, 1, idescd, wrk, lwrk, ifails )
    endif
    ifail = ifails
!
!  Free workspace
!
    deallocate(wrk2,stat=status)
    if (status/=0) call deallocate_error('pdiaggd','wrk2')
    deallocate(wrk,stat=status)
    if (status/=0) call deallocate_error('pdiaggd','wrk')
  endif
!
  t2 = g_cpu_time()
  tdiag = tdiag + t2 - t1
!
!  Trap diagonalisation failure
!
  if (ifail.gt.0) then
    call outerror('diagonalisation of dynamical matrix failed',0_i4)
    call stopnow('pdiaggd')
  endif
!
!  Convert frequency units - imaginary freqs denoted by negative no.
!
  do i = 1,mcv
    root = freq(i)
    if (root.ge.0.0_dp) then
      freq(i) = sqrt(root)*fscale
    else
      root = abs(root)
      freq(i) = - sqrt(root)*fscale
    endif
  enddo
#else
! 
!  Dummy assignment of values to avoid compiler warning
! 
  ifail = 0
  eigr(1:maxeigr,1) = 0.0_dp
  freq(1:mcvloc) = 0.0_dp
!
  call outerror('pdiaggd called when not compiled with MPI',0_i4)
  call stopnow('pdiaggd')
#endif
#ifdef TRACE
  call trace_out('pdiaggd')
#endif
!
  return
  end
