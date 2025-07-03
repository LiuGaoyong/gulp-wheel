  subroutine pdiagg(mcv,maxd2,derv2,eigr,freq,fscale,lvectors,lprint,ifail)
!
!  Calculate eigenvalues and eigenvectors of the dynamical matrix
!  and subsequent calculation of properties of eigenvectors.
!
!  Gamma point only version!!
!
!   3/02 Created from peigeng
!  10/04 Eispack call replaced by lapack
!   4/05 Modified so that eigr is used as workspace for dsyev so the
!        the dynamical matrix is preversed in derv2 for lvectors =
!        .true. case
!   6/09 Changes added for PDF output and lprint added to arguments
!   2/18 Trace added
!   7/20 Divide and conquer added and normalisation removed as its not needed
!   9/20 Wrapping of ELPA 2stage calls with ifdef added
!
!  On entry :
!
!  mcv      = number of modes
!  maxd2    = left-hand dimension of derv2 and eigr
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
!  Julian Gale, CIC, Curtin University, September 2020
!
  use g_constants
  use current
  use element
  use maths,        only : ldivide_and_conquer, leig_mrrr, eigtol
#ifdef ELPA
  use maths,        only : leig_2stage
#endif
  use parallel
  use species
  use times
#ifdef TRACE
  use trace,        only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(out)                    :: ifail
  integer(i4),  intent(in)                     :: maxd2
  integer(i4),  intent(in)                     :: mcv
  logical,      intent(in)                     :: lprint
  logical,      intent(in)                     :: lvectors
  real(dp),     intent(in)                     :: derv2(maxd2,*)
  real(dp),     intent(out)                    :: eigr(maxd2,*)
  real(dp),     intent(in)                     :: fscale
  real(dp),     intent(out)                    :: freq(*)
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: il
  integer(i4)                                  :: iu
  integer(i4)                                  :: ilaenv
  integer(i4),    dimension(:),   allocatable  :: isup
  integer(i4),    dimension(:),   allocatable  :: iw2
  integer(i4)                                  :: nb
  integer(i4)                                  :: nf
  integer(i4)                                  :: nsize1
  integer(i4)                                  :: nsize2
  integer(i4)                                  :: status
  real(dp),       dimension(:,:), allocatable  :: d2loc
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: root
  real(dp)                                     :: t1
  real(dp)                                     :: t2
  real(dp)                                     :: vl
  real(dp)                                     :: vu
  real(dp),       dimension(:),   allocatable  :: w1
#ifdef TRACE
  call trace_in('pdiagg')
#endif
!
!  Call eigen solver
!
  t1 = g_cpu_time()
  if (lvectors) then
!
!  Copy dynamical matrix to eigenvector array to avoid overwriting
!
    eigr(1:mcv,1:mcv) = derv2(1:mcv,1:mcv)
  endif
  if (ldivide_and_conquer) then
!-----------------------
!  Divide and conquer  |
!-----------------------
!
!  Initial dummy allocation of workspace
!
    allocate(w1(1),stat=status)
    if (status/=0) call outofmemory('pdiagg','w1')
    allocate(iw2(1),stat=status)
    if (status/=0) call outofmemory('pdiagg','iw2')
!
!  Set parameters for eigensolver by calling to obtain memory
!
    nsize1 = -1
    nsize2 = -1
    if (lvectors) then
#ifdef ELPA
      if (leig_2stage) then
        call dsyevd_2stage('V','U',mcv,eigr,maxd2,freq,w1,nsize1,iw2,nsize2,ifail)
      else
#endif
        call dsyevd('V','U',mcv,eigr,maxd2,freq,w1,nsize1,iw2,nsize2,ifail)
#ifdef ELPA
      endif
#endif
    else
#ifdef ELPA
      if (leig_2stage) then
        call dsyevd_2stage('N','U',mcv,derv2,maxd2,freq,w1,nsize1,iw2,nsize2,ifail)
      else
#endif
        call dsyevd('N','U',mcv,derv2,maxd2,freq,w1,nsize1,iw2,nsize2,ifail)
#ifdef ELPA
      endif
#endif
    endif
!
    nsize1 = nint(w1(1))
    nsize2 = iw2(1)
!
!  Deallocate workspace
!
    deallocate(iw2,stat=status)
    if (status/=0) call deallocate_error('pdiagg','iw2')
    deallocate(w1,stat=status)
    if (status/=0) call deallocate_error('pdiagg','w1')
!
!  Real allocation of workspace
!
    allocate(w1(nsize1),stat=status)
    if (status/=0) call outofmemory('pdiagg','w1')
    allocate(iw2(nsize2),stat=status)
    if (status/=0) call outofmemory('pdiagg','iw2')
!
!  Call eigensolver
!
    if (lvectors) then
#ifdef ELPA
      if (leig_2stage) then
        call dsyevd_2stage('V','U',mcv,eigr,maxd2,freq,w1,nsize1,iw2,nsize2,ifail)
      else
#endif
        call dsyevd('V','U',mcv,eigr,maxd2,freq,w1,nsize1,iw2,nsize2,ifail)
#ifdef ELPA
      endif
#endif
    else
#ifdef ELPA
      if (leig_2stage) then
        call dsyevd_2stage('N','U',mcv,derv2,maxd2,freq,w1,nsize1,iw2,nsize2,ifail)
      else
#endif
        call dsyevd('N','U',mcv,derv2,maxd2,freq,w1,nsize1,iw2,nsize2,ifail)
#ifdef ELPA
      endif
#endif
    endif
!
    deallocate(iw2,stat=status)
    if (status/=0) call deallocate_error('pdiagg','iw2')
    deallocate(w1,stat=status)
    if (status/=0) call deallocate_error('pdiagg','w1')
  elseif (leig_mrrr) then
!--------
!  MRRR |
!--------
!
!  Allocation of workspace for copy of matrix if needed
!
    allocate(d2loc(mcv,mcv),stat=status)
    if (status/=0) call outofmemory('pdiagg','d2loc')
!
    d2loc(1:mcv,1:mcv) = derv2(1:mcv,1:mcv)
!
!  Initial dummy allocation of workspace
!
    allocate(isup(2*mcv),stat=status)
    if (status/=0) call outofmemory('pdiagg','isup')
    allocate(w1(1),stat=status)
    if (status/=0) call outofmemory('pdiagg','w1')
    allocate(iw2(1),stat=status)
    if (status/=0) call outofmemory('pdiagg','iw2')
!
!  Set parameters for eigensolver by calling to obtain memory
!
    nsize1 = -1
    nsize2 = -1
    if (lvectors) then
#ifdef ELPA
      if (leig_2stage) then
        call dsyevr_2stage('V','A','U',mcv,d2loc,mcv,vl,vu,il,iu,eigtol,nf,freq,eigr,maxd2,isup,w1,nsize1,iw2,nsize2,ifail)
      else
#endif
        call dsyevr('V','A','U',mcv,d2loc,mcv,vl,vu,il,iu,eigtol,nf,freq,eigr,maxd2,isup,w1,nsize1,iw2,nsize2,ifail)
#ifdef ELPA
      endif
#endif
    else
#ifdef ELPA
      if (leig_2stage) then
        call dsyevr_2stage('N','A','U',mcv,d2loc,mcv,vl,vu,il,iu,eigtol,nf,freq,eigr,maxd2,isup,w1,nsize1,iw2,nsize2,ifail)
      else
#endif
        call dsyevr('N','A','U',mcv,d2loc,mcv,vl,vu,il,iu,eigtol,nf,freq,eigr,maxd2,isup,w1,nsize1,iw2,nsize2,ifail)
#ifdef ELPA
      endif
#endif
    endif
!
    nsize1 = nint(w1(1))
    nsize2 = iw2(1)
!
!  Deallocate workspace
!
    deallocate(iw2,stat=status)
    if (status/=0) call deallocate_error('pdiagg','iw2')
    deallocate(w1,stat=status)
    if (status/=0) call deallocate_error('pdiagg','w1')
!
!  Real allocation of workspace
!
    allocate(w1(nsize1),stat=status)
    if (status/=0) call outofmemory('pdiagg','w1')
    allocate(iw2(nsize2),stat=status)
    if (status/=0) call outofmemory('pdiagg','iw2')
!
!  Call eigensolver
!
    if (lvectors) then
#ifdef ELPA
      if (leig_2stage) then
        call dsyevr_2stage('V','A','U',mcv,d2loc,mcv,vl,vu,il,iu,eigtol,nf,freq,eigr,maxd2,isup,w1,nsize1,iw2,nsize2,ifail)
      else
#endif
        call dsyevr('V','A','U',mcv,d2loc,mcv,vl,vu,il,iu,eigtol,nf,freq,eigr,maxd2,isup,w1,nsize1,iw2,nsize2,ifail)
#ifdef ELPA
      endif
#endif
    else
#ifdef ELPA
      if (leig_2stage) then
        call dsyevr_2stage('N','A','U',mcv,d2loc,mcv,vl,vu,il,iu,eigtol,nf,freq,eigr,maxd2,isup,w1,nsize1,iw2,nsize2,ifail)
      else
#endif
        call dsyevr('N','A','U',mcv,d2loc,mcv,vl,vu,il,iu,eigtol,nf,freq,eigr,maxd2,isup,w1,nsize1,iw2,nsize2,ifail)
#ifdef ELPA
      endif
#endif
    endif
!
    deallocate(iw2,stat=status)
    if (status/=0) call deallocate_error('pdiagg','iw2')
    deallocate(w1,stat=status)
    if (status/=0) call deallocate_error('pdiagg','w1')
    deallocate(isup,stat=status)
    if (status/=0) call deallocate_error('pdiagg','isup')
!
    deallocate(d2loc,stat=status)
    if (status/=0) call deallocate_error('pdiagg','d2loc')
#ifdef ELPA
  elseif (leig_2stage) then
!------------------------
!  Standard with 2stage |
!------------------------
!
!  Initial dummy allocation of workspace
!
    allocate(w1(1),stat=status)
    if (status/=0) call outofmemory('pdiagg','w1')
!
!  Set parameters for eigensolver by calling to obtain memory
!
    nsize1 = -1
    if (lvectors) then
      call dsyev_2stage('V','U',mcv,eigr,maxd2,freq,w1,nsize1,ifail)
    else
      call dsyev_2stage('N','U',mcv,derv2,maxd2,freq,w1,nsize1,ifail)
    endif
!
    nsize1 = nint(w1(1))
!
!  Deallocate workspace
!
    deallocate(w1,stat=status)
    if (status/=0) call deallocate_error('pdiagg','w1')
!
!  Real allocation of workspace
!
    allocate(w1(nsize1),stat=status)
    if (status/=0) call outofmemory('pdiagg','w1')
!
!  Call eigensolver
!
    if (lvectors) then
      call dsyev_2stage('V','U',mcv,eigr,maxd2,freq,w1,nsize1,ifail)
    else
      call dsyev_2stage('N','U',mcv,derv2,maxd2,freq,w1,nsize1,ifail)
    endif
!
    deallocate(w1,stat=status)
    if (status/=0) call deallocate_error('pdiagg','w1')
#endif
  else
!-------------
!  Standard  |
!-------------
!
!  Set parameters for eigensolver
!
    nb = ilaenv(1_i4,'DSYTRD','U',mcv,-1_i4,-1_i4,-1_i4)
    nsize1 = max(1,(nb+2)*mcv)
    allocate(w1(nsize1+10),stat=status)
    if (status/=0) call outofmemory('pdiagg','w1')
!
!  Call eigensolver
!
    if (lvectors) then
#ifdef ELPA
      if (leig_2stage) then
        call dsyev_2stage('V','U',mcv,eigr,maxd2,freq,w1,nsize1,ifail)
      else
#endif
        call dsyev('V','U',mcv,eigr,maxd2,freq,w1,nsize1,ifail)
#ifdef ELPA
      endif
#endif
    else
#ifdef ELPA
      if (leig_2stage) then
        call dsyev_2stage('N','U',mcv,derv2,maxd2,freq,w1,nsize1,ifail)
      else
#endif
        call dsyev('N','U',mcv,derv2,maxd2,freq,w1,nsize1,ifail)
#ifdef ELPA
      endif
#endif
    endif
    deallocate(w1,stat=status)
    if (status/=0) call deallocate_error('pdiagg','w1')
  endif
  t2 = g_cpu_time()
  tdiag = tdiag + t2 - t1
!
!  Trap diagonalisation failure
!
  if (ifail.gt.0) then
    call outerror('diagonalisation of dynamical matrix failed',0_i4)
    call stopnow('pdiagg')
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
#ifdef TRACE
  call trace_out('pdiagg')
#endif
!
  return
  end
