  subroutine outkpt
!
!  Outputs kpoints for a structure
!
!  nkpt = total number of k points across all structures
!  xkpt = fractional x component of k point
!  ykpt = fractional y component of k point
!  zkpt = fractional z component of k point
!  wkpt = weight of each k point
!
!   2/02 Normalisation of K points corrected for multiple configurations
!   5/15 Output of kpt file added
!  11/15 Output of weights changed to give more decimal places
!   6/17 Module files renamed to gulp_files
!   7/19 Trap added for polarisation and k points
!   3/21 Trap added for non-gamma phonons for GFNFF
!   7/21 Output of k points stopped for molecules
!  10/21 lgfnff moved to control
!  10/21 lnongamma initialised for molecules
!   4/23 Gamma point check added for ReaxFF
!   7/23 Non-gamma point phonons enabled for new fc algorithm
!   8/23 Now to be called by all processors rather than just ioproc
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
!  Julian Gale, CIC, Curtin University, July 2023
!
  use configurations, only : ncfg
  use control
  use current
  use gulp_files,     only : lkpt
  use iochannels
  use ksample
  use parallel,       only : ioproc
  use polarise,       only : lpolar
  implicit none
!
!  Local variables
!
  integer(i4)                               :: i
  integer(i4)                               :: nk
  integer(i4)                               :: status
  logical                                   :: lfout
  logical                                   :: lpout
  logical                                   :: lnongamma
  real(dp), dimension(:), allocatable       :: sumpercfg
!
  lfout = .true.
  lpout = (ioproc.and.index(keyword,'nokp').eq.0)
  nk = 0
  allocate(sumpercfg(ncfg),stat=status)
  if (status/=0) call outofmemory('outkpt','sumpercfg')
  sumpercfg(1:ncfg) = 0.0_dp
  do i = 1,nkpt
    sumpercfg(nkptcfg(i)) = sumpercfg(nkptcfg(i)) + wkpt(i)
  enddo
  do i = 1,ncfg
    if (sumpercfg(i).gt.0.0_dp) then
      sumpercfg(i) = 1.0_dp/sumpercfg(i)
    else
      sumpercfg(i) = 1.0_dp
    endif
  enddo
!
!  Normalise weights
!
  do i = 1,nkpt
    wkpt(i) = wkpt(i)*sumpercfg(nkptcfg(i))
  enddo
  deallocate(sumpercfg,stat=status)
  if (status/=0) call deallocate_error('outkpt','sumpercfg')
!
!  Initialise flag for non-gamma point calculation
!
  lnongamma = .false.
!
!  Loop over k points and output those relating to this configuration
!
  if (ndim.gt.0) then
    lnongamma = .false.
    do i = 1,nkpt
      if (nkptcfg(i).eq.ncf) then
        if (lfout.and.lpout) then
          lfout = .false.
          write(ioout,'(''  Brillouin zone sampling points :'',/)')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  Point number          x          y          z            Weight'')')
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        endif
        nk = nk + 1
        if (lpout) then
          write(ioout,'(4x,i5,10x,3(f9.6,2x),2x,f11.8)') nk,xkpt(i),ykpt(i),zkpt(i),wkpt(i)
        endif
        if (abs(xkpt(i)).gt.1.0d-4) lnongamma = .true.
        if (abs(ykpt(i)).gt.1.0d-4) lnongamma = .true.
        if (abs(zkpt(i)).gt.1.0d-4) lnongamma = .true.
      endif
    enddo
    if (.not.lfout) then
      write(ioout,'(''--------------------------------------------------------------------------------'',/)')
    endif
  endif
!
!  If EEM/QEq/GFNFF then non-gamma point phonons are not allowed yet unless using the new FC algorithm
!
  lkptgammaonlycfg(ncf) = (.not.lnongamma)
  if (leem.or.lgfnff.or.lreaxFF) then
    if (.not.lkptgammaonlycfg(ncf)) lusefcalg = .true.
  endif
  if (.not.lusefcalg) then
    if (leem.and.lnongamma) then
      call outerror('Only gamma point phonons allowed for EEM/QEq',0_i4)
      call stopnow('outkpt')
    endif
    if (lgfnff.and.lnongamma) then
      call outerror('Only gamma point phonons allowed for GFNFF',0_i4)
      call stopnow('outkpt')
    endif
    if (lreaxFF.and.lnongamma) then
      call outerror('Only gamma point phonons allowed for ReaxFF',0_i4)
      call stopnow('outkpt')
    endif
  endif
!
!  If polarisation is being used then non-gamma point phonons are not allowed yet unless using the new FC algorithm
!
  if (lpolar.and.lnongamma) then
    call outerror('Only gamma point phonons allowed for point ion polarisation',0_i4)
    call stopnow('outkpt')
  endif
!
!  Output k points to a file if required
!
  if (lkpt) call outkptfile
!
  return
  end
