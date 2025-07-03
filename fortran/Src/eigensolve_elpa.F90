  subroutine eigensolve_elpa(ng,ngb,d_matrix,maxld,e_matrix,maxle,eig,lvectors)
!
!  Wrapper to interface GULP to the ELPA eigensolver
!  Involves mapping 1D distribution to 2D and back again
!
!  ng       = global size of problem
!  ngb      = blocksize for problem
!  d_matrix = square symmetric matrix with initial matrix
!  maxld    = left-hand dimension of d_matrix
!  e_matrix = square symmetric matrix with eigenvectors
!  maxle    = left-hand dimension of e_matrix
!  eig      = eigenvalue array
!  lvectors = indicates whether eigenvectors are to be computed
!
!   7/20 Created
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use iochannels
  use parallel
#ifdef ELPA
  use maths,         only : lelpa, elpatype
  use elpa
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),      intent(in)                 :: ng
  integer(i4),      intent(in)                 :: ngb
  integer(i4),      intent(in)                 :: maxld
  integer(i4),      intent(in)                 :: maxle
  logical,          intent(in)                 :: lvectors
  real(dp),         intent(in)                 :: d_matrix(maxld,ng)
  real(dp),         intent(out)                :: e_matrix(maxle,ng)
  real(dp),         intent(out)                :: eig(ng)
#ifdef ELPA
!
!  Local variables
!
  integer(i4)                                  :: status
!
!  Local variables in Scalapack/Blacs/MPI integer precision
!
  integer                                      :: iBlacsContext2D
  integer                                      :: iBlacsContextSys
  integer                                      :: ifails
  integer                                      :: idescd(9)
  integer                                      :: idescd2D(9)
  integer                                      :: idesce(9)
  integer                                      :: ld
  integer                                      :: mcol2D
  integer                                      :: mrow2D
  integer                                      :: npcol2D
  integer                                      :: npcol2Dloc
  integer                                      :: nprow2D
  integer                                      :: nprow2Dloc
  integer                                      :: nb
  integer                                      :: nv
  integer                                      :: numroc
  real*8,      dimension(:,:), allocatable     :: d2D
  real*8,      dimension(:,:), allocatable     :: e2D
!
!  ELPA specific variables
!
  class(elpa_t),                       pointer :: t_elpa
  integer                                      :: elpa_status
  logical                                      :: lfailed
!
  if (nprocs.eq.1) then
    call outerror('eigensolve with parallel ELPA called in serial',0_i4)
    call stopnow('eigensolve_elpa')
  endif
!
  if (lelpa) then
!--------------------------------------------------------
!  Pre-ELPA : Convert 1D data decomposition to 2D form  |
!--------------------------------------------------------
    nb = ngb
    nv = ng
!
!  Get system context
!
    call blacs_get(0, 0, iBlacsContextSys)
!
!  Set up Blacs descriptors for arrays in current form
!
    ld = maxld
    call descinit(idescd, nv, nv, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - idescd ',0_i4)
      call stopnow('eigensolve_elpa')
    endif
!
    if (lvectors) then
      ld = maxle
      call descinit(idesce, nv, nv, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
      if (ifails.ne.0) then
        call outerror('initialisation in descinit failed - idesce ',0_i4)
        call stopnow('eigensolve_elpa')
      endif
    endif
!
!  Try to set as close to a square grid as possible
!
    do npcol2D = nint(sqrt(dble(nprocs))),2,-1
      if (mod(nprocs,npcol2D) == 0 ) exit
    enddo
    nprow2D = nprocs/npcol2D
!
!  Initialise 2D processor grid
!
    iBlacsContext2D = iBlacsContextSys
    call BLACS_Gridinit(iBlacsContext2D,'C',nprow2D,npcol2D)
    call BLACS_Gridinfo(iBlacsContext2D,nprow2D,npcol2D,nprow2Dloc,npcol2Dloc)
!
!  Determine sizes of local matrices for 2D grid
!
    mrow2D = numroc(nv, nb, nprow2Dloc, 0, nprow2D)
    mcol2D = numroc(nv, nb, npcol2Dloc, 0, npcol2D)
!
!  Allocate arrays for 2D storage of D and E
!
    allocate(d2D(mrow2D,mcol2D),stat=status)
    if (status/=0) call outofmemory('eigensolve_elpa','d2D')
    if (lvectors) then
      allocate(e2D(mrow2D,mcol2D),stat=status)
      if (status/=0) call outofmemory('eigensolve_elpa','e2D')
    endif
!
!  Set up Blacs descriptors for arrays in 2D form
!
    ld = mrow2D
    call descinit(idescd2D, nv, nv, nb, nb, 0, 0, iBlacsContext2D, ld, ifails)
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - idescd2D ',0_i4)
      call stopnow('eigensolve_elpa')
    endif
!
!  Copy D to 2D form
!
    call pdgemr2d(nv, nv, d_matrix, 1, 1, idescd, d2D, 1, 1, idescd2D, iBlacsContext)
!-------------------------------------
!  ELPA Solve of eigenvalue problem  |
!-------------------------------------
!
!  Initialise ELPA
!
    if (elpa_init(20200417).ne.elpa_ok) then
      call outerror('Initialisation of ELPA has failed - API version not supported',0_i4)
      call stopnow('eigensolve_elpa')
    endif
!
!  Allocate ELPA data type
!
    t_elpa => elpa_allocate(elpa_status)
    if (elpa_status.ne.elpa_ok) then
      call outerror('Initialisation of ELPA has failed - data type not allocated',0_i4)
      call stopnow('eigensolve_elpa')
    endif
!
!  Set information regarding matrix and MPI distribution
!
    lfailed = .false.
    call t_elpa%set("mpi_comm_parent",MPI_Comm_GULP,elpa_status)      ! MPI communicator launching ELPA
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("na",nv,elpa_status)                              ! Size of global matrix
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("nblk",nb,elpa_status)                            ! Blocksize
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    if (lvectors) then
      call t_elpa%set("nev",nv,elpa_status)                           ! Number of eigenvectors required
    else
      call t_elpa%set("nev",0,elpa_status)                            ! Number of eigenvectors required
    endif
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("local_nrows",mrow2D,elpa_status)                 ! Number of local rows
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("local_ncols",mcol2D,elpa_status)                 ! Number of local columns
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("process_row",nprow2Dloc,elpa_status)             ! Local processor row
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("process_col",npcol2Dloc,elpa_status)             ! Local processor column
    if (elpa_status.ne.elpa_ok) lfailed = .true.
!
!  Check whether an error has occured in setting the above quantities
!
    if (lfailed) then
      call outerror('Initialisation of ELPA has failed - could not set matrix or MPI info',0_i4)
      call stopnow('eigensolve_elpa')
    endif
!
!  Setup the ELPA object
!
    elpa_status = t_elpa%setup()
    if (elpa_status.ne.elpa_ok) then
      call outerror('Initialisation of ELPA has failed - ELPA object setup has failed',0_i4)
      call stopnow('eigensolve_elpa')
    endif
!
!  Set ELPA solver to be used
!
    if (elpatype.eq.1) then
      call t_elpa%set("solver",elpa_solver_1stage,elpa_status)
    else
      call t_elpa%set("solver",elpa_solver_2stage,elpa_status)
    endif
    if (elpa_status.ne.elpa_ok) then
      call outerror('Initialisation of ELPA has failed - ELPA solver type not set',0_i4)
      call stopnow('eigensolve_elpa')
    endif
!
!  Call ELPA solver
!
    if (lvectors) then
      call t_elpa%eigenvectors(d2D, eig, e2D, elpa_status)
    else
      call t_elpa%eigenvalues(d2D, eig, elpa_status)
    endif
    if (elpa_status.ne.elpa_ok) then
      call outerror('Solve by ELPA has failed',0_i4)
      call stopnow('eigensolve_elpa')
    endif
!
!  Deallocate ELPA memory
!
    call elpa_deallocate(t_elpa,elpa_status)
    if (elpa_status.ne.elpa_ok) then
      call outerror('Deallocation of ELPA has failed',0_i4)
      call stopnow('eigensolve_elpa')
    endif
!
!  Shutdown ELPA
!
    call elpa_uninit()
!--------------------------------------------------------------
!  Post-ELPA : Convert 2D data decomposition back to 1D form  |
!--------------------------------------------------------------
    if (lvectors) then
!
!  Copy e2D back to 1D form
!
      call pdgemr2d(nv, nv, e2D, 1, 1, idescd2D, e_matrix, 1, 1, idesce, iBlacsContext)
!
!  Deallocate arrays for 2D storage of D and E
!
      deallocate(e2D,stat=status)
      if (status/=0) call deallocate_error('eigensolve_elpa','e2D')
    endif
    deallocate(d2D,stat=status)
    if (status/=0) call deallocate_error('eigensolve_elpa','d2D')
!
!  Shutdown blacs context
!
    call blacs_gridexit(iBlacsContext2D)
  endif
#else
  call outerror('eigensolve with parallel ELPA called but GULP was not compiled with -e',0_i4)
  call stopnow('eigensolve_elpa')
#endif
  return
  end
!
  subroutine c_eigensolve_elpa(ng,ngb,d_matrix,maxld,e_matrix,maxle,eig,lvectors)
!
!  Wrapper to interface GULP to the ELPA eigensolver
!  Involves mapping 1D distribution to 2D and back again
!  Complex version.
!
!  ng       = global size of problem
!  ngb      = blocksize for problem
!  d_matrix = square symmetric matrix with initial matrix
!  maxld    = left-hand dimension of d_matrix
!  e_matrix = square symmetric matrix with eigenvectors
!  maxle    = left-hand dimension of e_matrix
!  eig      = eigenvalue array
!  lvectors = indicates whether eigenvectors are to be computed
!
!   7/20 Created
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
!  Julian Gale, CIC, Curtin University, July 2020
!
  use iochannels
#ifdef ELPA
  use maths,         only : lelpa, elpatype
#endif
  use parallel
#ifdef ELPA
  use elpa
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),      intent(in)                 :: ng
  integer(i4),      intent(in)                 :: ngb
  integer(i4),      intent(in)                 :: maxld
  integer(i4),      intent(in)                 :: maxle
  logical,          intent(in)                 :: lvectors
  complex(dpc),     intent(in)                 :: d_matrix(maxld,ng)
  complex(dpc),     intent(out)                :: e_matrix(maxle,ng)
  real(dp),         intent(out)                :: eig(ng)
#ifdef ELPA
!
!  Local variables
!
  integer(i4)                                  :: status
!
!  Local variables in Scalapack/Blacs/MPI integer precision
!
  integer                                      :: iBlacsContext2D
  integer                                      :: iBlacsContextSys
  integer                                      :: ifails
  integer                                      :: idescd(9)
  integer                                      :: idescd2D(9)
  integer                                      :: idesce(9)
  integer                                      :: ld
  integer                                      :: mcol2D
  integer                                      :: mrow2D
  integer                                      :: npcol2D
  integer                                      :: npcol2Dloc
  integer                                      :: nprow2D
  integer                                      :: nprow2Dloc
  integer                                      :: nb
  integer                                      :: nv
  integer                                      :: numroc
  complex*16,  dimension(:,:), allocatable     :: d2D
  complex*16,  dimension(:,:), allocatable     :: e2D
!
!  ELPA specific variables
!
  class(elpa_t),                       pointer :: t_elpa
  integer                                      :: elpa_status
  logical                                      :: lfailed
!
  if (nprocs.eq.1) then
    call outerror('eigensolve with parallel ELPA called in serial',0_i4)
    call stopnow('c_eigensolve_elpa')
  endif
!
  if (lelpa) then
!--------------------------------------------------------
!  Pre-ELPA : Convert 1D data decomposition to 2D form  |
!--------------------------------------------------------
    nb = ngb
    nv = ng
!
!  Get system context
!
    call blacs_get(0, 0, iBlacsContextSys)
!
!  Set up Blacs descriptors for arrays in current form
!
    ld = maxld
    call descinit(idescd, nv, nv, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - idescd ',0_i4)
      call stopnow('c_eigensolve_elpa')
    endif
!
    if (lvectors) then
      ld = maxle
      call descinit(idesce, nv, nv, nb, nb, 0, 0, iBlacsContext, ld, ifails)
!
      if (ifails.ne.0) then
        call outerror('initialisation in descinit failed - idesce ',0_i4)
        call stopnow('c_eigensolve_elpa')
      endif
    endif
!
!  Try to set as close to a square grid as possible
!
    do npcol2D = nint(sqrt(dble(nprocs))),2,-1
      if (mod(nprocs,npcol2D) == 0 ) exit
    enddo
    nprow2D = nprocs/npcol2D
!
!  Initialise 2D processor grid
!
    iBlacsContext2D = iBlacsContextSys
    call BLACS_Gridinit(iBlacsContext2D,'C',nprow2D,npcol2D)
    call BLACS_Gridinfo(iBlacsContext2D,nprow2D,npcol2D,nprow2Dloc,npcol2Dloc)
!
!  Determine sizes of local matrices for 2D grid
!
    mrow2D = numroc(nv, nb, nprow2Dloc, 0, nprow2D)
    mcol2D = numroc(nv, nb, npcol2Dloc, 0, npcol2D)
!
!  Allocate arrays for 2D storage of D and E
!
    allocate(d2D(mrow2D,mcol2D),stat=status)
    if (status/=0) call outofmemory('c_eigensolve_elpa','d2D')
    if (lvectors) then
      allocate(e2D(mrow2D,mcol2D),stat=status)
      if (status/=0) call outofmemory('c_eigensolve_elpa','e2D')
    endif
!
!  Set up Blacs descriptors for arrays in 2D form
!
    ld = mrow2D
    call descinit(idescd2D, nv, nv, nb, nb, 0, 0, iBlacsContext2D, ld, ifails)
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed - idescd2D ',0_i4)
      call stopnow('c_eigensolve_elpa')
    endif
!
!  Copy D to 2D form
!
    call pzgemr2d(nv, nv, d_matrix, 1, 1, idescd, d2D, 1, 1, idescd2D, iBlacsContext)
!-------------------------------------
!  ELPA Solve of eigenvalue problem  |
!-------------------------------------
!
!  Initialise ELPA
!
    if (elpa_init(20200417).ne.elpa_ok) then
      call outerror('Initialisation of ELPA has failed - API version not supported',0_i4)
      call stopnow('c_eigensolve_elpa')
    endif
!
!  Allocate ELPA data type
!
    t_elpa => elpa_allocate(elpa_status)
    if (elpa_status.ne.elpa_ok) then
      call outerror('Initialisation of ELPA has failed - data type not allocated',0_i4)
      call stopnow('c_eigensolve_elpa')
    endif
!
!  Set information regarding matrix and MPI distribution
!
    lfailed = .false.
    call t_elpa%set("mpi_comm_parent",MPI_Comm_GULP,elpa_status)      ! MPI communicator launching ELPA
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("na",nv,elpa_status)                              ! Size of global matrix
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("nblk",nb,elpa_status)                            ! Blocksize
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    if (lvectors) then
      call t_elpa%set("nev",nv,elpa_status)                           ! Number of eigenvectors required
    else
      call t_elpa%set("nev",0,elpa_status)                            ! Number of eigenvectors required
    endif
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("local_nrows",mrow2D,elpa_status)                 ! Number of local rows
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("local_ncols",mcol2D,elpa_status)                 ! Number of local columns
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("process_row",nprow2Dloc,elpa_status)             ! Local processor row
    if (elpa_status.ne.elpa_ok) lfailed = .true.
    call t_elpa%set("process_col",npcol2Dloc,elpa_status)             ! Local processor column
    if (elpa_status.ne.elpa_ok) lfailed = .true.
!
!  Check whether an error has occured in setting the above quantities
!
    if (lfailed) then
      call outerror('Initialisation of ELPA has failed - could not set matrix or MPI info',0_i4)
      call stopnow('c_eigensolve_elpa')
    endif
!
!  Setup the ELPA object
!
    elpa_status = t_elpa%setup()
    if (elpa_status.ne.elpa_ok) then
      call outerror('Initialisation of ELPA has failed - ELPA object setup has failed',0_i4)
      call stopnow('c_eigensolve_elpa')
    endif
!
!  Set ELPA solver to be used
!
    if (elpatype.eq.1) then
      call t_elpa%set("solver",elpa_solver_1stage,elpa_status)
    else
      call t_elpa%set("solver",elpa_solver_2stage,elpa_status)
    endif
    if (elpa_status.ne.elpa_ok) then
      call outerror('Initialisation of ELPA has failed - ELPA solver type not set',0_i4)
      call stopnow('c_eigensolve_elpa')
    endif
!
!  Call ELPA solver
!
    if (lvectors) then
      call t_elpa%eigenvectors(d2D, eig, e2D, elpa_status)
    else
      call t_elpa%eigenvalues(d2D, eig, elpa_status)
    endif
    if (elpa_status.ne.elpa_ok) then
      call outerror('Solve by ELPA has failed',0_i4)
      call stopnow('c_eigensolve_elpa')
    endif
!
!  Deallocate ELPA memory
!
    call elpa_deallocate(t_elpa,elpa_status)
    if (elpa_status.ne.elpa_ok) then
      call outerror('Deallocation of ELPA has failed',0_i4)
      call stopnow('c_eigensolve_elpa')
    endif
!
!  Shutdown ELPA
!
    call elpa_uninit()
!--------------------------------------------------------------
!  Post-ELPA : Convert 2D data decomposition back to 1D form  |
!--------------------------------------------------------------
    if (lvectors) then
!
!  Copy e2D back to 1D form
!
      call pzgemr2d(nv, nv, e2D, 1, 1, idescd2D, e_matrix, 1, 1, idesce, iBlacsContext)
!
!  Deallocate arrays for 2D storage of D and E
!
      deallocate(e2D,stat=status)
      if (status/=0) call deallocate_error('c_eigensolve_elpa','e2D')
    endif
    deallocate(d2D,stat=status)
    if (status/=0) call deallocate_error('c_eigensolve_elpa','d2D')
!
!  Shutdown blacs context
!
    call blacs_gridexit(iBlacsContext2D)
  endif
#else
  call outerror('eigensolve with parallel ELPA called but GULP was not compiled with -e',0_i4)
  call stopnow('c_eigensolve_elpa')
#endif
  return
  end
