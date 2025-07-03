module m_gulp_interface

  implicit none
  logical   :: linitialised = .false.   ! Flag to ensure that initialisation is called before energy
  private
  public    :: init_gulp     ! Initialisation routine
  public    :: gulp_energy   ! Routine for energy and optionally gradients
#ifdef MPI
  public    :: init_gulp_mpi ! Initialisation routine for MPI in GULP
#endif

contains

#ifdef MPI
  subroutine init_gulp_mpi(MPI_comm_in)
!
!  Performs setup tasks for GULP in parallel
!
!  12/23 Created from initcomms.F90
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
!  Julian Gale, CIC, Curtin University, December 2023
!
  use datatypes
  use parallel
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: MPI_comm_in     ! MPI communicator that GULP is being called from
!
!  Local variables
!
  integer                 :: ierr
  integer                 :: lprocid
  integer                 :: lnprocs
!
!  Assign GULP communicator based on that passed in
!
  MPI_comm_GULP = MPI_comm_in
!
!  Find number of processors
!
  call MPI_comm_rank(MPI_comm_GULP,lprocid,ierr)
  call MPI_comm_size(MPI_comm_GULP,lnprocs,ierr)
!
  procid  = lprocid
  nprocs  = lnprocs
!
  end subroutine init_gulp_mpi
#endif

  subroutine init_gulp(lgulpoutput,dimensions,natoms,cell,atomicsymbols,xyz,charges,keyword_in,library_in)
!
!  Performs setup tasks for GULP when called from another code
!
!  11/23 Created from gulpsetup.F90
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
!  Julian Gale, CIC, Curtin University, November 2023
!
  use configurations
  use control
  use current
  use element,         only : maxele, atmass, rvdw
  use general
  use iochannels,      only : iotmp, ioout
  use parallel
  use shifts,          only : nshcfg
  use species
#ifdef TRACE
  use trace,           only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  character(len=*), intent(in)  :: keyword_in            ! Keyword string to control GULP
  character(len=*), intent(in)  :: library_in            ! Name of library file containing potential parameters for GULP
  integer(i4),      intent(in)  :: dimensions            ! Number of periodic dimensions (0 -> 3)
  integer(i4),      intent(in)  :: natoms                ! Number of atoms 
  logical,          intent(in)  :: lgulpoutput           ! If true then allow GULP to print output
  real(dp),         intent(in)  :: cell(3,dimensions)    ! Cell vectors (if dimensions > 0)      (in Angstroms)
  real(dp),         intent(in)  :: xyz(3,natoms)         ! Atomic coordinates (Cartesian)  (in Angstroms)
  real(dp),         intent(in)  :: charges(natoms)       ! Atomic charges (in a.u.)
  character(len=*), intent(in)  :: atomicsymbols(natoms) ! Atomic symbols (that GULP will allow; element symbol, plus optionally a number)
!
  character(len=40)             :: hostname
  character(len=maxwordlength)  :: word
  integer(i4)                   :: i
  integer(i4)                   :: iline
  integer(i4)                   :: j
  integer(i4)                   :: hlength
  integer(i4)                   :: status
  logical                       :: lfound
  logical                       :: librarypresent
  logical                       :: lopened
  logical                       :: lsymbol
  real(dp)                      :: area
  real(dp)                      :: vol
  real(dp)                      :: volume
#ifdef MPI
  include 'mpif.h'
#endif
#ifdef TRACE
  call trace_in('init_gulp')
#endif
#ifdef MPI
  if (lgulpoutput) then
    ioproc = (procid.eq.0)
  endif
#else
  if (lgulpoutput) then
    ioproc = .true.
  endif
  nprocs = 1
#endif
!*************************
!  Nullify all pointers  *
!*************************
  call nullpointer
!*************************
!  Initialise memory     *
!*************************
  call initmemory
!*************************
!  Initialise variables  *
!*************************
  call initial
!******************
!  Output header  *
!******************
  version = 6.3_dp
  if (lgulpoutput) then
    call gulp_banner
  endif
!**************************
!  Get local information  *
!**************************
  call local
!****************************
!  Set element information  *
!****************************
  call setele(lopened)
!************************************************************
!  Update array sizes to accommodate structure if required  *
!************************************************************
  if (natoms.gt.maxat) then
    maxat = natoms
    call changemaxat
  endif
  if (natoms.gt.maxatot) then
    maxatot = natoms
    call changemaxatot
  endif
!***********************************************************
!  Create temporary input file that contains library file  *
!***********************************************************
  librarypresent = (library_in.ne.' ')
  call channels(iotmp,.true.)
  if (librarypresent) then
    write(iotmp,'(a)') 'library '//trim(library_in)
  else
    write(iotmp,'(a)') 'start'
  endif
  rewind(iotmp)
!*****************
!  Set keywords  *
!*****************
  keyword = keyword_in
  call setkeyword(.true.)
!******************
!  Set structure  *
!******************
  ncf = 1
  ncfg = 1
  ndimen(ncfg) = dimensions
  lvecin(ncfg) = .true.
  if (dimensions.gt.0) then
    rvcfg(1:3,1:dimensions,1) = cell(1:3,1:dimensions)
!
!  Check volume/area/length is sensible
!
    if (dimensions.eq.3) then
      vol = volume(rvcfg(1,1,1))
      if (vol.lt.0.5_dp) then
        call outerror('Unit cell volume is too small - check cell parameters',0_i4)
        call stopnow('init_gulp')
      endif
    elseif (dimensions.eq.2) then
      vol = area(rvcfg(1,1,1))
      if (vol.lt.0.5_dp) then
        call outerror('Surface cell area is too small - check cell parameters',0_i4)
        call stopnow('init_gulp')
      endif
    elseif (dimensions.eq.1) then
      vol = rvcfg(1,1,1)
      if (vol.lt.0.5_dp) then
        call outerror('Polymer cell length is too small - check cell parameters',0_i4)
        call stopnow('init_gulp')
      endif
    endif
  endif
  nshcfg(ncfg) = 1
!
  nasym = natoms
  numat = natoms
  nascfg(ncfg) = natoms
!
!  Convert coordinates from Cartesian to fractional if periodic
!
  if (dimensions.gt.0) then
    do i = 1,natoms
      call cart2frac(ndimen(ncfg),xyz(1,i),xyz(2,i),xyz(3,i),cell,xcfg(i),ycfg(i),zcfg(i),icosx(i),icosy(i),icosz(i))
    enddo
  else
    do i = 1,natoms
      xcfg(i) = xyz(1,i)
      ycfg(i) = xyz(2,i)
      zcfg(i) = xyz(3,i)
    enddo
  endif
!
!  Check symbols
!
  do i = 1,natoms
    word = atomicsymbols(i)
    call stolc(word,maxwordlength)
    call worsy(word,lsymbol,.false.)
    if (.not.lsymbol) then
      call outerror('Invalid atomic symbol specified for atom ',i)
      call stopnow('init_gulp')
    endif
    word = atomicsymbols(i)
    call ltont(word,natcfg(i),ntypcfg(i))
  enddo
!****************
!  Set charges  *
!****************
  qlcfg(1:natoms) = charges(1:natoms)
!************************************************************************
!  Set other atomic properties not handled by the interface at present  *
!************************************************************************
  occucfg(1:natoms) = 1.0_dp
  radcfg(1:natoms) = 0.0_dp
!**************************
!  Set number of species  *
!**************************
  do i = 1,natoms
    lfound = .false.
    j = 0
    do while (.not.lfound.and.j.lt.nspec)
      j = j + 1
      lfound = (natcfg(i).eq.natspec(j).and.(ntypcfg(i).eq.ntypspec(j).or.ntypspec(j).eq.0))
    enddo
    if (lfound) then
!
!  Existing species
!
      nspecptrcfg(i) = j
    else
!
!  New species - add to list
!
      nspec = nspec + 1
      if (nspec.gt.maxspec) then
        maxspec = nspec + 20
        call changemaxspec
      endif
      natspec(nspec) = natcfg(i)
      ntypspec(nspec) = ntypcfg(i)
      qlspec(nspec) = qlcfg(i)
      radspec(nspec) = 0.0_dp
      lbrspec(nspec) = .false.
      lqinspec(nspec) = .false.
      if (natspec(nspec).le.maxele) then
        massspec(nspec) = atmass(natspec(nspec))
      else
        massspec(nspec) = 0.0_dp
      endif
      if (natspec(nspec).le.maxele) then
        vdwspec(nspec) = rvdw(natspec(nspec))
      else
        vdwspec(nspec) = rvdw(natspec(nspec)-maxele)
      endif
!
!  Loop over previous species to check whether charge should be
!  set by any of these as qlcfg may not contain anything
!
      do j = 1,nspec - 1
        if (natspec(nspec).eq.natspec(j).and.ntypspec(j).eq.0) then
          qlspec(nspec) = qlspec(j)
        endif
      enddo
!
!  Set pointer from atom to species
!
      nspecptrcfg(i) = nspec
    endif
  enddo
!**************************
!  Process library input  *
!**************************
  iline = 0
  call inword(iline)
!
!  Update charges in case they were input by species in the library
!
  do i = 1,numat
    if (lqinspec(nspecptr(i))) then
      qlcfg(i) = qlspec(nspecptr(i))
    endif
  enddo
!*******************************************************
!  Set keywords again in case any were in the library  *
!*******************************************************
  call setkeyword(.true.)
  close(iotmp,status='delete')
!****************************************
!  Set charge equilibration parameters  *
!****************************************
  call seteem
!***************************
!  Output keyword details  *
!***************************
  if (lgulpoutput.and.ioproc) then
    call outkey
  endif
!**********************************************************************************
!  Set flag for stress output to force strain derivative calculation if required  *
!**********************************************************************************
  if (ndim.gt.0) lstressout = .true.
!*********************************************************
!  Re-initialise memory as keywords may change settings  *
!*********************************************************
  call initmemory
!**************
!  Site name  *
!**************
  if (lgulpoutput.and.ioproc) then
    if (site.ne.' ') then
      write(ioout,'(''* '',a76,'' *'')') site(1:76)
      write(ioout,'(''********************************************************************************'')')
    endif
    call datetime(1_i4)
    write(ioout,'(''  Number of CPUs = '',i5,/)') nprocs
    hostname = ' '
    call get_environment_variable('HOSTNAME',hostname,hlength,status)
    if (status.eq.0) then
      write(ioout,'(''  Host name      = '',a40,/)') hostname
    elseif (status.eq.1) then
      call get_environment_variable('HOST',hostname,hlength,status)
      if (status.eq.0) then
        write(ioout,'(''  Host name      = '',a40,/)') hostname
      endif
    endif
  endif
!*************************
!  GFNFF initialisation  *
!*************************
  if (lgfnff) then
#ifndef NOGFNFF
    call gulp_gfnff_init
    call pgfnff_init
#endif
  endif
!**************************
!  One off initial set up *
!**************************
  call setcfg
!*******************
!  Species output  *
!*******************
  if (lgulpoutput.and.ioproc) then
    call outspec
  endif
!*****************************
!  Electronegativity output  *
!*****************************
  if (lgulpoutput.and.ioproc) then
    call outeem
  endif
!************************
!  Output general info  *
!************************
  if (lgulpoutput.and.ioproc) call outgen
!*******************************
!  Output polarisability info  *
!*******************************
  if (lgulpoutput.and.ioproc) call outpolar
!*********************
!  Check potentials  *
!*********************
  call checkpot
!**********************
!  Output potentials  *
!**********************
  if (lgulpoutput.and.ioproc) then
    call outpot
  endif
!*****************************
!  Check options for method  *
!*****************************
  call methodok
#ifdef TRACE
  call trace_out('init_gulp')
#endif
!
  linitialised = .true.
!
  return
  end subroutine init_gulp
!
  subroutine gulp_energy(cell,xyz,charges,enrgy,gradients,strainderivatives,lgrad1,ierror)
!
!  Returns the energy and optionally the first derivatives from GULP
!
!  11/23 Created
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
!  Julian Gale, CIC, Curtin University, November 2023
!
  use control
  use current
  use derivatives
  use general
  use parallel
#ifdef TRACE
  use trace,       only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(out)   :: ierror                ! Error flag
  logical,     intent(in)    :: lgrad1                ! If true then compute gradients
  real(dp),    intent(in)    :: cell(3,*)             ! Cell vectors for periodic systems
  real(dp),    intent(in)    :: xyz(3,*)              ! Cartesian coordinates in Angstroms
  real(dp),    intent(out)   :: charges(*)            ! Charges from energy calculation
  real(dp),    intent(out)   :: enrgy                 ! Energy (in eV)
  real(dp),    intent(out)   :: gradients(3,*)        ! First derivatives of energy (in eV/Ang)
  real(dp),    intent(out)   :: strainderivatives(*)  ! Strain derivatives (in eV) for periodic systems according to Voight convention
!
!  Local variables
!
  integer(i4)                :: i
#ifdef TRACE
  call trace_in('gulp_energy')
#endif
!
  ierror = 0
!
!  Check whether initialisation has been performed
!
  if (.not.linitialised) then
    call outerror('GULP energy called before initialisation!',0_i4)
    ierror = 1
    enrgy = 1.0d8
    return
  endif
!*********************
!  Update structure  *
!*********************
  if (ndim.gt.0) then
    rv(1:3,1:ndim) = cell(1:3,1:ndim)
  endif
!  
!  Convert coordinates from Cartesian to fractional if periodic
!
  do i = 1,numat
    call cart2frac(ndim,xyz(1,i),xyz(2,i),xyz(3,i),cell,xfrac(i),yfrac(i),zfrac(i),icosx(i),icosy(i),icosz(i))
    xafrac(i) = xfrac(i)
    yafrac(i) = yfrac(i)
    zafrac(i) = zfrac(i)
  enddo
!
  if (ndim.gt.0) then
!
!  Make sure all fractional coords are between 0 and 1
!
    if (lmodco) then
      if (ndim.eq.3) then
        do i = 1,numat
          xafrac(i) = mod(xafrac(i)+1000.0_dp,1.0_dp)
          yafrac(i) = mod(yafrac(i)+1000.0_dp,1.0_dp)
          zafrac(i) = mod(zafrac(i)+1000.0_dp,1.0_dp)
        enddo
      elseif (ndim.eq.2) then
        do i = 1,numat
          xafrac(i) = mod(xafrac(i)+1000.0_dp,1.0_dp)
          yafrac(i) = mod(yafrac(i)+1000.0_dp,1.0_dp)
        enddo
      elseif (ndim.eq.1) then
        do i = 1,numat
          xafrac(i) = mod(xafrac(i)+1000.0_dp,1.0_dp)
        enddo
      endif
    endif
  endif
!
  call x0tostr
!***********
!  Energy  *
!***********
!
!  Check core-shell distances
!
  call cutscheck
!
!  Evaluate function and first derivatives
!
  call energy(enrgy,lgrad1,.false.)
!
!  Complete strain derivatives
!
  if (lgrad1.and.ndim.gt.0) call strfin(.false.)
!
!  First derivatives
!
  if (lgrad1) then
    gradients(1,1:numat) = xdrv(1:numat)
    gradients(2,1:numat) = ydrv(1:numat)
    gradients(3,1:numat) = zdrv(1:numat)
    if (ndim.gt.0) then
      strainderivatives(1:nstrains) = strderv(1:nstrains)
    endif
  endif
!
!  Return charges in case they have been changed
!
  charges(1:numat) = qf(1:numat)
#ifdef TRACE
  call trace_out('gulp_energy')
#endif
!
  return
  end

end module m_gulp_interface
