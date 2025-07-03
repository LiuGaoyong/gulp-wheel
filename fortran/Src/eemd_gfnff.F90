  subroutine eemd_gfnff(ecoul,cn,dcn,d2cn,lgrad1,lgrad2,lphoncall)
! 
!  Subroutine for calculating the charges for the GFNFF force field
!  => Allows for fragments and coordination number dependence
!
!  Periodic boundary condition version
!
!  Distributed memory version for parallel case.
!
!   4/21 Created from eemd and gfnff_eem
!   4/21 Call to dbcgsolve has an extra argument added
!   5/21 Exclusion of atoms from neem list based on region numbers removed
!   6/21 Call to dbcgsolve changed to dcgsolve
!   7/21 Symmetry handled by reducing full set of charges at the end
!   7/21 Parallel matrix inversion and solve added
!   8/21 Call to sumall changed to sumone
!  10/21 Modules rearranged
!   1/22 Call to dcharged_gfnff added
!   1/22 Call to gfnff_drv2_dcn_dq modified
!   1/22 lgrad2 now passed to dcharge
!   2/22 Matrix inversion solve added in parallel
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
  use energies
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
#ifdef MPI
  include 'mpif.h'
#endif
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
  integer(i4)                                  :: ie
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ifrag
  integer(i4)                                  :: ii
  integer(i4)                                  :: il
  integer(i4)                                  :: iloc
  integer(i4)                                  :: iopr
  integer(i4)                                  :: iter
  integer(i4)                                  :: j
  integer(i4)                                  :: n
  integer(i4),                           save  :: ncf_last = 0
  integer(i4)                                  :: ninfrag
  integer(i4)                                  :: nloc
  integer(i4)                                  :: nloceem
  integer(i4), dimension(:), allocatable       :: nlocptr
  integer(i4)                                  :: nmax
  integer(i4)                                  :: nr
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: status
#ifdef MPI
  integer                                      :: MPIerror
  integer                                      :: ntag
  integer                                      :: nnode
  integer                                      :: ntmp
  integer                                      :: Request
  integer(i4), dimension(:), allocatable       :: StatMPI       ! Array for status from MPI
#endif
  logical                                      :: lneedinverse
  logical                                      :: lprint
  logical                                      :: lqiter
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: cni
  real(dp)                                     :: degselfdcn
  real(dp)                                     :: d2egselfdcn2
  real(dp)                                     :: d2egselfdcndq
  real(dp)                                     :: egself
  real(dp)                                     :: egselftot
  real(dp)                                     :: egself_before
  real(dp)                                     :: err
  real(dp),    dimension(:,:), allocatable     :: emat
  real(dp)                                     :: qi
  real(dp)                                     :: reqv
  real(dp)                                     :: sr_cni
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp),                  parameter         :: tsqrt2pi = 0.797884560802866_dp
#ifdef MPI
  real(dp),    dimension(:), allocatable       :: tmp
#endif
  real(dp),    dimension(:), allocatable       :: vfield
  real(dp),    dimension(:), allocatable       :: z
#ifdef TRACE
  call trace_in('eemd_gfnff')
#endif
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
!  Allocate local memory
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
  nmax = numat + nfrag
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
  if (status/=0) call outofmemory('eemd_gfnff','z')
  allocate(nlocptr(neem+nfrag),stat=status)
  if (status/=0) call outofmemory('eemd_gfnff','nlocptr')
!
!  Set up pointers to local elements for parallelisation
!
  nloc = 0
  nloceem = 0
  do i = procid+1,neem+nfrag,nprocs
    nloc = nloc + 1
    nlocptr(nloc) = i
    if (i.le.neem) nloceem = nloceem + 1
  enddo
!
!  Allocate the memory for the potential array
!
  allocate(emat(nmax,nloc),stat=status)
  if (status/=0) call outofmemory('eemd_gfnff','emat')
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    allocate(vfield(numat),stat=status)
    if (status/=0) call outofmemory('eemd_gfnff','vfield')
    vfield(1:numat) = 0.0_dp
    call electricfieldpotl(vfield)
  endif
!
  if (luseprecon) then
    if (nmax.gt.maxprecon.or.nloc.gt.maxpreconu) then
      maxprecon = max(nmax,maxprecon)
      maxpreconu = max(nloc,maxpreconu)
      call changemaxprecon
    endif
  endif
!
!  Calculate reciprocal lattice vectors
!
  if (lewald.and.ndim.gt.1) call kindex
!
!  Zero right hand vector
!
  z(1:neem+nfrag) = 0.0_dp
!*******************************************
!  Generate electrostatic potential terms  *
!*******************************************
!
!  Ensure emat is full initialised
!
  emat(1:nmax,1:nloc) = 0.0_dp
  if (.not.lnoqeem) then
    call genpotd(nloceem,nlocptr,emat,nmax,z)
  endif
!********************************
!  Form matrix of coefficients  *
!********************************
  if (nloc.gt.0) then
    do il = 1,nloceem
      i = nlocptr(il)
      ii = neemptr(i)
      emat(i,il) = emat(i,il) + (tsqrt2pi/sqrt(gfnff_eeq_alp(ii)) + gfnff_eeq_gam(ii))*angstoev
    enddo
!
!  Handle fragments
!
    do il = 1,nloc
      i = nlocptr(il)
      if (i.gt.neem) then
        ifrag = i - neem 
        do j = 1,neem
          ii = neemptr(j)
          if (nfraglist(ii).eq.ifrag) then
            emat(j,il) = occuf(ii)
          endif
        enddo
        emat(neem+ifrag,il) = 0.0_dp
      else
        do j = 1,nfrag
          ii = neemptr(i)
          if (nfraglist(ii).eq.j) then
            emat(neem+j,il) = 1.0_dp
          endif
        enddo
      endif
    enddo
  endif
!
!  Fragment charges
!
  do i = 1,nfrag
    z(neem+i) = qfrag(i)
  enddo
!
!  Add external potential
!
  do i = 1,neem
    ii = neemptr(i)
    z(i) = z(i) - extpot(nrelf2a(ii))
  enddo
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
!  Debugging output
!
  if (index(keyword,'debu').ne.0) then
    if (ioproc) then
      write(ioout,'(/,''  EEM(GFNFF) Matrix :'',/)')
    endif
#ifdef MPI
    call mpbarrier
    if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
      ntmp = neem + nfrag
      ntag = 1
      allocate(tmp(ntmp),stat=status)
      if (status/=0) call outofmemory('eemd_gfnff','tmp')
      allocate(StatMPI(MPI_Status_Size),stat=status)
      if (status/=0) call outofmemory('eemd_gfnff','StatMPI')
!
      iopr = 0
      il = 0
      do i = 1,neem+nfrag
        if (iopr.ne.0_i4) then
!
!  Post receive
!
          if (ioproc) then
            nnode = iopr
            call MPI_IRecv(tmp,ntmp,MPI_double_precision,nnode, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
!
!  Pass data to ioproc for writing
!
          if (procid.eq.iopr) then
            il = il + 1
            tmp(1:nmax) = emat(1:nmax,il)
!
!  Post send
!
            call MPI_ISend(tmp,ntmp,MPI_double_precision,0, &
                           ntag,MPI_Comm_World,Request,MPIerror)
          endif
          if (ioproc.or.procid.eq.iopr) then
            call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          endif
          if (ioproc) then
!
!  Write on I/O node
!
            write(ioout,'(i5,10(1x,f9.5))') i,(tmp(j),j=1,neem+nfrag),z(i)
          endif
        elseif (ioproc) then
          il = il + 1
          write(ioout,'(i5,10(1x,f9.5))') i,(emat(j,il),j=1,neem+nfrag),z(i)
        endif
        iopr = iopr + 1
        iopr = mod(iopr,nprocs)
        call mpbarrier
      enddo
!
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('eemd_gfnff','StatMPI')
      deallocate(tmp,stat=status)
      if (status/=0) call deallocate_error('eemd_gfnff','tmp')
    else
#endif
      iopr = 0
      il = 0
      do i = 1,neem+nfrag
        if (procid.eq.iopr) then
          il = il + 1
          write(ioout,'(i5,10(1x,f9.5))') i,(emat(j,il),j=1,neem+nfrag),z(i)
        endif
        iopr = iopr + 1
        iopr = mod(iopr,nprocs)
      enddo
#ifdef MPI
    endif
#endif
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
      call outerror('iterative charge solution failed in eemd_gfnff',0_i4)
      call stopnow('eemd_gfnff')
    endif
  elseif (.not.lneedinverse) then
!*****************
!  Direct solve  *
!*****************
    ifail = 0
    n = neem + nfrag
    qf(1:n) = z(1:n)
    qf(1:neem+nfrag) = z(1:neem+nfrag)
    call matrix_solve(.true.,n,nmax,1_i4,emat,qf,0_i4,ifail)
  else
!******************
!  Invert matrix  *
!******************
    ifail = 0
    n = neem + nfrag
!************************
!  Symmetric inversion  *
!************************
    call matrix_inversion_library(n,1_i4,nmax,1_i4,emat,0_i4,.true.,ifail)
!
!  Was inversion successful?
!
    if (ifail.ne.0) then
      call outerror('matrix inversion failed in GFNFF EEM',0_i4)
      call stopnow('eemd_gfnff')
    endif
!
!  Multiply inverse matrix and chi matrix to get charges
!
    qf(1:neem+nfrag) = 0.0_dp
    do i = 1,nloc
      ii = nlocptr(i)
      qf(ii) = 0.0_dp
      do j = 1,n
        qf(ii) = qf(ii) + z(j)*emat(j,i)
      enddo
    enddo
    if (nprocs.gt.1) then
!
!  Global sum of charges
!
      call sumall(qf,qa,n,"eemd_gfnff","qf")
      qf(1:neem+nfrag) = qa(1:neem+nfrag)
    endif
!
!  Save inverse matrix for use in preconditioning
!
    if (luseprecon) then
      lpreconsaved = .true.
      precon_matrix(1:n,1:nloc) = emat(1:n,1:nloc)
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
    call dcharged_gfnff(lprint,emat,nmax,.true.,lgrad2,cn,dcn)
  endif
!**************************
!  Calculate self energy  *
!**************************
  egself = 0.0_dp
  do ie = 1,neem
    i = neemptr(ie)
    iloc = atom2local(i)
    if (iloc.ne.0) then
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
            call gfnff_drv2_dcn_selfpd_p1(iloc,i,degselfdcn,d2egselfdcn2,dcn,d2cn)
          else
            call gfnff_drv2_dcn_selfd_p1(iloc,i,degselfdcn,d2egselfdcn2,dcn,d2cn,.true.)
          endif
!
!  Charge-CN products
!
          d2egselfdcndq = - 0.5_dp*reqv*gfnff_eeq_cnf(i)/sr_cni
          if (lphoncall) then
            call gfnff_drv2_dcn_dqpd_p1(iloc,i,d2egselfdcndq,dcn)
          else
            call gfnff_drv2_dcn_dqd_p1(iloc,i,d2egselfdcndq,dcn)
          endif
        else
          call gfnff_drv_dcn(i,degselfdcn,cn,dcn)
        endif
      endif
    endif
  enddo
  if (lgrad2) then
!
!  Second pass loop for distributed memory second derivatives
!
    do ii = 1,neem
      i = neemptr(ii)
      qi = qf(i)
      reqv = occuf(i)
      cni = cn(i)
      sr_cni = sqrt(cni)
!
!  Derivatives of self energy due to coordination number
!
      if (sr_cni.gt.1.0d-12) then
        degselfdcn = - 0.5_dp*reqv*qi*gfnff_eeq_cnf(i)/sr_cni
        d2egselfdcn2 = 0.25_dp*reqv*qi*gfnff_eeq_cnf(i)/sr_cni**3
        if (lphoncall) then
          call gfnff_drv2_dcn_selfpd_p2(i,degselfdcn,d2egselfdcn2,dcn,d2cn)
        else
          call gfnff_drv2_dcn_selfd_p2(i,degselfdcn,d2egselfdcn2,dcn,d2cn,.true.)
        endif
!
!  Charge-CN products
!
        d2egselfdcndq = - 0.5_dp*reqv*gfnff_eeq_cnf(i)/sr_cni
        if (lphoncall) then
          call gfnff_drv2_dcn_dqpd_p2(i,d2egselfdcndq,dcn)
        else
          call gfnff_drv2_dcn_dqd_p2(i,d2egselfdcndq,dcn)
        endif
      endif
    enddo
  endif
!*******************
!  Output results  *
!*******************
!
!  Sum self energies
!
  if ((index(keyword,'debu').ne.0)) then
    call sumone(egself,egselftot,"eemd_gfnff","egself")
    if (ioproc) then
      write(ioout,'(//,''  Final charges from GFNFF :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Self energy       = '',f16.6,'' eV'')') egselftot
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
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
  deallocate(emat,stat=status)
  if (status/=0) call deallocate_error('eemd_gfnff','emat')
  deallocate(nlocptr,stat=status)
  if (status/=0) call deallocate_error('eemd_gfnff','nlocptr')
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('eemd_gfnff','z')
!
!  Timing
!
  time2 = g_cpu_time()
  teem = teem + time2 - time1
#ifdef TRACE
  call trace_out('eemd_gfnff')
#endif
!
  return
  end
