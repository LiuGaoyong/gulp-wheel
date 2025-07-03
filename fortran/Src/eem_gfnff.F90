  subroutine eem_gfnff(ecoul,cn,dcn,d2cn,lgrad1,lgrad2,lphoncall)
!
!  Subroutine for calculating the charges for the GFNFF force field
!  => Allows for fragments and coordination number dependence
!
!  Periodic boundary condition version
!  Now uses only the asymmetric unit
!
!   9/20 Created from eem
!  10/20 Derivatives of coordination number added
!  11/20 Second derivatives started
!   2/21 Second derivatives of self energy added w.r.t. coordination number
!   2/21 Cutoff no longer passed to some coordination number routines
!   3/21 Phonon call argument added
!   4/21 Call to dbcgsolve has an extra argument added
!   5/21 Exclusion of atoms from neem list based on region numbers removed
!   6/21 Improved handling of qf initialisation added
!   6/21 Use of direct solve added
!   6/21 Use of inverse matrix as preconditioner added
!   7/21 Non-zero fragment charges now implemented 
!   7/21 Symmetry handled by reducing full set of charges at the end
!   7/21 Call to matrix inversion routine updated
!  10/21 Modules rearranged
!   1/22 Calls to gfnff_drv2_dcn_dq modified
!   7/23 nregionnocfg added
!   7/23 Update of qlcfg restricted to non-supercell calls
!   7/23 extpot added
!   8/23 Electric fields added
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
!  Copyright Curtin University 2023
!
!  Julian Gale, CIC, Curtin University, August 2023
!
  use control
  use configurations
  use current
  use derivatives
  use eemdata
  use element
  use energies,        only : eregion2region, siteenergy
  use field,           only : lfieldcfg, ntdfieldcfg
  use gulp_gfnff,      only : gfnff_eeq_alp, gfnff_eeq_chi, gfnff_eeq_cnf, gfnff_eeq_gam
  use gulp_gfnff,      only : nfrag, qfrag, nfraglist
  use iochannels
  use m_precondition
  use parallel
  use partial
  use symmetry
  use times
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),                     intent(out)    :: ecoul
  real(dp),                     intent(in)     :: cn(numat)
  real(dp),                     intent(in)     :: dcn(numat)
  real(dp),                     intent(in)     :: d2cn(*)
  logical,                      intent(in)     :: lgrad1
  logical,                      intent(in)     :: lgrad2
  logical,                      intent(in)     :: lphoncall
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ia
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: iter
  integer(i4)                                  :: j
  integer(i4)                                  :: k
  integer(i4)                                  :: n
  integer(i4)                                  :: ninfrag
  integer(i4)                                  :: nloc
  integer(i4), dimension(:), allocatable       :: nlocptr
  integer(i4)                                  :: nmax
  integer(i4)                                  :: nmaxu
  integer(i4)                                  :: nr
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: status
  integer(i4),                           save  :: ncf_last = 0
  logical                                      :: lneedinverse
  logical                                      :: lqiter
  logical                                      :: lprint
  real(dp)                                     :: cni
  real(dp)                                     :: g_cpu_time
  real(dp),    dimension(:,:), allocatable     :: emat
  real(dp)                                     :: err
  real(dp)                                     :: degselfdcn
  real(dp)                                     :: d2egselfdcn2
  real(dp)                                     :: d2egselfdcndq
  real(dp)                                     :: egself
  real(dp)                                     :: egself_before
  real(dp)                                     :: qi
  real(dp)                                     :: reqv
  real(dp)                                     :: sr_cni
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp),                  parameter         :: tsqrt2pi = 0.797884560802866_dp
  real(dp),    dimension(:), allocatable       :: vfield
  real(dp),    dimension(:), allocatable       :: z
#ifdef TRACE
  call trace_in('eem_gfnff')
#endif
!
!  If this is a parallel run then call distributed memory version
!
  if (nprocs.gt.1) then
    call eemd_gfnff(ecoul,cn,dcn,d2cn,lgrad1,lgrad2,lphoncall)
#ifdef TRACE
    call trace_out('eem_gfnff')
#endif
    return
  endif
!
  time1 = g_cpu_time()
!
  ecoul = 0.0_dp
!
!  Set flag for iterative charges - can't be used for all algorithms
!
  lqiter = literativeQ
!
!  Set flag as to whether we need to compute the inverse of the EEM matrix
!
  lneedinverse = (lgrad2.or.ldcharge)
!
!  Preconditioner flags
!
  if (luseprecon.and..not.lpreconsaved) lneedinverse = .true.
!
  if (lneedinverse) lqiter = .false.
!
  neem = 0
  neemrptr(1:numat) = 0
!
!  Check elements 
!
  do i = 1,numat
    neem = neem + 1
    neemptr(neem) = i
    neemrptr(i) = neem
  enddo
!
!  Check the memory for the linear arrays
!
  if (numat+nfrag.gt.maxat) then
    maxat = numat + nfrag
    call changemaxat
  endif
!
!  Allocate the memory for the main matrix
!
  nmaxu = numat + nfrag
  nmax = numat + nfrag
!
  allocate(emat(nmax,nmaxu),stat=status)
  if (status/=0) call outofmemory('eem_gfnff','emat')
!
  if (luseprecon) then
    if (nmax.gt.maxprecon) then
      maxprecon = nmax
      maxpreconu = nmax
      call changemaxprecon
    endif
  endif
!
!  Set the pointer to where the electronegativity should be as well
!
  do i = 1,nfrag
    neemptr(neem+i) = numat + i
  enddo
!
!  Allocate local memory that depends on neem
!
  allocate(z(max(neem+nfrag,numat)),stat=status)
  if (status/=0) call outofmemory('eem_gfnff','z')
!
!  If iterative then set up pointers to local elements
!
  if (lqiter) then
    nloc = neem + nfrag
    allocate(nlocptr(nloc),stat=status)
    if (status/=0) call outofmemory('eem_gfnff','nlocptr')
    do i = 1,nloc
      nlocptr(i) = i
    enddo
  endif
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    allocate(vfield(numat),stat=status)
    if (status/=0) call outofmemory('eem_gfnff','vfield')
    vfield(1:numat) = 0.0_dp
    call electricfieldpotl(vfield)
  endif
!
!  Calculate reciprocal lattice vectors
!
  if (lewald.and.ndim.gt.1) call kindex
!
!  Zero right hand vector
!
  z(1:neem+nfrag) = 0.0_dp
!
!  Ensure emat is full initialised
!
  emat(1:nmax,1:nmaxu) = 0.0_dp
!*******************************************
!  Generate electrostatic potential terms  *
!*******************************************
  if (.not.lnoqeem) then
    call genpot(emat,nmax,z,1_i4)
  endif
!********************************
!  Form matrix of coefficients  *
!********************************
  do ii = 1,neem
    i = neemptr(ii)
    emat(ii,ii) = emat(ii,ii) + (tsqrt2pi/sqrt(gfnff_eeq_alp(i)) + gfnff_eeq_gam(i))*angstoev
  enddo
!
!  Add external potential
!
  do i = 1,neem
    ii = neemptr(i)
    z(i) = z(i) - extpot(nrelf2a(ii))
  enddo
!
  do ii = 1,neem
    i = neemptr(ii)
    z(ii) = z(ii) + gfnff_eeq_chi(i) + gfnff_eeq_cnf(i)*sqrt(cn(i))
  enddo
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    do i = 1,neem
      ii = neemptr(i)
      z(i) = z(i) + vfield(ii)   
    enddo
  endif
!
!  Handle fragments
!
  do i = 1,nfrag
    z(neem+i) = qfrag(i)
    do j = 1,neem
      ii = neemptr(j)
      if (nfraglist(ii).eq.i) then
        emat(j,neem+i) = occuf(ii)
        emat(neem+i,j) = 1.0_dp
      endif
    enddo
    emat(neem+i,neem+i) = 0.0_dp
  enddo
!
!  Debugging output
!
  if (index(keyword,'debu').ne.0.and.ioproc) then
    write(ioout,'(/,''  EEM(GFNFF) Matrix (eV) :'',/)')
    do i = 1,neem
      write(ioout,'(10(1x,f9.5))')(emat(j,i),j=1,neem),(emat(neem+k,i),k=1,nfrag),z(i)
    enddo
    do i = 1,nfrag
      write(ioout,'(10(1x,f9.5))')(emat(j,neem+i),j=1,neem+nfrag),z(neem+i)
    enddo
  endif
  if (lqiter) then
!***************************
!  Iterative charge solve  *
!***************************
!
!  Set initial guess for solution to fragment electronegativity if configuration has changed
!
    if (ncf.ne.ncf_last) then
      do i = 1,nfrag
        ninfrag = 0
        qf(neem+i) = 0.0_dp
        do j = 1,neem
          ii = neemptr(j)
          if (nfraglist(ii).eq.i) then
            ninfrag = ninfrag + 1
            qf(neem+i) = qf(neem+i) + z(j)
          endif
        enddo
        qf(neem+i) = qf(neem+i)/dble(ninfrag)
      enddo
      ncf_last = ncf
    endif
!
!  Solve using iterative route
!
    call dcgsolve(neem+nfrag,nloc,nlocptr,emat,nmax,z,qf,qitertol,nqitermax,iter,err)
!
    if (ioproc) then
      if (index(keyword,'verb').ne.0) then
        write(ioout,'('' Number of iterations / error in dbcgsolve = '',i4,1x,f16.14)') iter,err
      endif
    endif
!
!  Was iterative solution successful?
!
    if (iter.ge.nqitermax) then
      call outerror('iterative charge solution failed in eem',0_i4)
      call stopnow('eem_gfnff')
    endif
  elseif (.not.lneedinverse) then
!*****************
!  Direct solve  *
!*****************
    ifail = 0
    n = neem + nfrag
    qf(1:n) = z(1:n)
    call matrix_solve(.true.,n,nmax,nblocksize,emat,qf,0_i4,ifail)
  else
!******************
!  Invert matrix  *
!******************
    ifail = 0
    n = neem + nfrag
!************************
!  Symmetric inversion  *
!************************
    call matrix_inversion_library(n,1_i4,nmax,nblocksize,emat,0_i4,.true.,ifail)
!
!  Was inversion successful?
!
    if (ifail.ne.0) then
      call outerror('matrix inversion failed in GFNFF_EEM',0_i4)
      call stopnow('eem_gfnff')
    endif
!  
!  Multiply inverse matrix and chi matrix to get charges
!
    do i = 1,neem + nfrag
      ii = neemptr(i)
      qf(ii) = 0.0_dp
      do j = 1,neem + nfrag
        qf(ii) = qf(ii) + z(j)*emat(j,i)
      enddo
    enddo
!
!  Save inverse matrix for use in preconditioning
!
    if (luseprecon) then
      lpreconsaved = .true.
      precon_matrix(1:n,1:n) = emat(1:n,1:n)
    endif
  endif
!
!  Transfer charges to qa
!
  do i = 1,nasym
    nr = nrela2f(i)
    qa(i) = qf(nr)
  enddo
!
!  Store charges in configurational array (unless this is a temporary supercell call)
!
  if (nascfg(ncf).eq.nasym) then
    do i = 1,nasym
      qlcfg(nsft+i) = qa(i)
    enddo
  endif
!*********************************
!  Calculate charge derivatives  *
!*********************************
  if (lgrad2.or.ldcharge) then
    lprint = (index(keyword,'dcha').ne.0) 
    call dcharge_gfnff(lprint,emat,nmax,.true.,cn,dcn)
  endif
!**************************
!  Calculate self energy  *
!**************************
  egself = 0.0_dp
  do ii = 1,neem
    i = neemptr(ii)
    ia = nrelf2a(i)
    qi = qf(i)
    reqv = occuf(i)
    cni = cn(i)
    sr_cni = sqrt(cni)
    egself_before = egself
!
    egself = egself - reqv*qi*(gfnff_eeq_chi(i) + gfnff_eeq_cnf(i)*sr_cni) + &
      0.5_dp*reqv*qi*qi*(gfnff_eeq_gam(i) + tsqrt2pi/sqrt(gfnff_eeq_alp(i)))*angstoev
!
!  Add external potential for site
!
    egself = egself + qi*reqv*extpot(ia)
!
    nregioni = nregionno(ia)
    eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + egself - egself_before
!
    siteenergy(i) = siteenergy(i) + egself - egself_before
!
!  Derivatives of self energy due to coordination number
!
    if (lgrad1.and.sr_cni.gt.1.0d-12) then
      degselfdcn = - 0.5_dp*reqv*qi*gfnff_eeq_cnf(i)/sr_cni
      if (lgrad2) then
        d2egselfdcn2 = 0.25_dp*reqv*qi*gfnff_eeq_cnf(i)/sr_cni**3
        if (lphoncall) then
          call gfnff_drv2_dcn_selfp(i,degselfdcn,d2egselfdcn2,dcn,d2cn)
        else
          call gfnff_drv2_dcn_self(i,degselfdcn,d2egselfdcn2,dcn,d2cn,.true.)
        endif
!
!  Charge-CN products
!
        d2egselfdcndq = - 0.5_dp*reqv*gfnff_eeq_cnf(i)/sr_cni
        if (lphoncall) then
          call gfnff_drv2_dcn_dqp(i,i,d2egselfdcndq,dcn)
        else
          call gfnff_drv2_dcn_dq(i,i,d2egselfdcndq,dcn)
        endif
      else
        call gfnff_drv_dcn(i,degselfdcn,cn,dcn)
      endif
    endif
  enddo
!*******************
!  Output results  *
!*******************
  if ((index(keyword,'verb').ne.0).and.ioproc) then
    write(ioout,'(//,''  Final charges from GFNFF :'',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nasym
      write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Self energy       = '',f16.6,'' eV'')') egself
    write(ioout,'(''--------------------------------------------------------------------------------'')')
  endif
!
!  Add up energy components
!
  ecoul = ecoul + egself
!
!  Free local memory 
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    deallocate(vfield,stat=status)
    if (status/=0) call deallocate_error('eemd_gfnff','vfield')
  endif
  if (lqiter) then
    deallocate(nlocptr,stat=status)
    if (status/=0) call deallocate_error('eem_gfnff','nlocptr')
  endif
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('eem_gfnff','z')
  deallocate(emat,stat=status)
  if (status/=0) call deallocate_error('eem_gfnff','emat')
!
!  Timing
!
  time2 = g_cpu_time()
  teem = teem + time2 - time1
#ifdef TRACE
  call trace_out('eem_gfnff')
#endif
!
  return
  end
