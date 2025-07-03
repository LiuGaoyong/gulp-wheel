program main
!
!  This is a simple program that illustrates the calling of GULP as an energy/gradient engine in serial
!
!  Julian Gale, Curtin University, December 2023
!
  use datatypes
  use m_gulp_interface

  implicit none

  character(len=4), allocatable :: atomicsymbols(:)
  character(len=80)             :: keywords
  character(len=80)             :: libraryfile
  integer(i4)                   :: ierror
  integer(i4)                   :: natoms
  integer(i4)                   :: ndim
  integer(i4)                   :: i
  integer(i4)                   :: j
  logical                       :: lgradients
  logical                       :: lgulpoutput
  real(dp)                      :: energy
  real(dp)                      :: cell(3,3)
  real(dp)                      :: strainderivatives(6)
  real(dp),         allocatable :: charges(:)
  real(dp),         allocatable :: gradients(:,:)
  real(dp),         allocatable :: xyz(:,:)
!
!  Control parameters
!
  lgulpoutput = .false.         ! If .true. then GULP will output calculation details
  lgradients  = .true.          ! If .true. then GULP will compute gradients as well as the energy
  keywords = ' '                ! Insert any keywords to pass to GULP here (or in the library file)
  libraryfile = 'nacl.lib  '    ! Insert the name of the library file with force field info here
!
!  Set system details
!
  natoms = 8                    ! Number of atoms in the system
  ndim = 3                      ! Dimensionality of system; 0 => cluster, 1 => polymer, 2 => slab, 3 => solid
!
  allocate(xyz(3,natoms))         ! Array for Cartesian coordinates of the atoms (Angstrom)
  allocate(gradients(3,natoms))   ! Array for Cartesian gradients acting on the atoms (eV/Angstrom)
  allocate(charges(natoms))       ! Array to pass the charges of the atoms (a.u.)
  allocate(atomicsymbols(natoms)) ! Atomic symbols in GULP format (element symbol + optionally a number)
!
!  Set atom types
!
  atomicsymbols(1:4) = 'Na '
  atomicsymbols(5:8) = 'Cl '
!
!  Set unit cell vectors
!
  cell(1:3,1:3) = 0.0_dp          ! Cell vectors (depending on dimensionality) (in Angstrom)
  cell(1,1) = 5.64_dp
  cell(2,2) = 5.64_dp
  cell(3,3) = 5.64_dp
!
!  Set atomic coordinates
!
  xyz(1,1) = 0.00_dp
  xyz(2,1) = 0.00_dp
  xyz(3,1) = 0.00_dp
!
  xyz(1,2) = 0.00_dp
  xyz(2,2) = 2.82_dp
  xyz(3,2) = 2.82_dp
!
  xyz(1,3) = 2.82_dp
  xyz(2,3) = 0.00_dp
  xyz(3,3) = 2.82_dp
!
  xyz(1,4) = 2.82_dp
  xyz(2,4) = 2.82_dp
  xyz(3,4) = 0.00_dp
!
  xyz(1,5) = 2.82_dp
  xyz(2,5) = 2.82_dp
  xyz(3,5) = 2.82_dp
!
  xyz(1,6) = 2.82_dp
  xyz(2,6) = 0.00_dp
  xyz(3,6) = 0.00_dp
!
  xyz(1,7) = 0.00_dp
  xyz(2,7) = 2.82_dp
  xyz(3,7) = 0.00_dp
!
  xyz(1,8) = 0.00_dp
  xyz(2,8) = 0.00_dp
  xyz(3,8) = 2.82_dp
!
!  Set charges - note here they are set via the library file
!
  charges(1:8) = 0.0_dp
!
!  Call to initialise GULP for this structure
!
  call init_gulp(lgulpoutput,ndim,natoms,cell,atomicsymbols,xyz,charges,keywords,libraryfile)
!
  write(6,'(/)')
  write(6,'(''################################################################################'')')
  write(6,'(''#  GULP called as a subroutine : Configuration 1 '')')
  write(6,'(''################################################################################'')')
!
!  Call GULP to get the energy and gradients
!
  call gulp_energy(cell,xyz,charges,energy,gradients,strainderivatives,lgradients,ierror)
!
  write(6,'(/,''  Cell: '')') 
  do i = 1,ndim
    write(6,'(i4,3(1x,f12.6))') i,(cell(j,i),j=1,3)
  enddo
  write(6,'(/,''  Coordinates: '')') 
  do i = 1,natoms
    write(6,'(i4,3(1x,f12.6))') i,(xyz(j,i),j=1,3)
  enddo
  write(6,'(/,''  Energy = '',f12.6,'' eV'')') energy
  write(6,'(/,''  Charges : '')') 
  do i = 1,natoms
    write(6,'(i4,1x,f12.6)') i,charges(i)
  enddo
  if (lgradients) then
    write(6,'(/,''  Gradients : '')') 
    do i = 1,natoms
      write(6,'(i4,3(1x,f12.6))') i,(gradients(j,i),j=1,3)
    enddo
    write(6,'(/,''  Strain gradients : eV '')') 
    write(6,'(4x,6(1x,f12.6))') (strainderivatives(j),j=1,6)
  endif
!
  write(6,'(/)')
  write(6,'(''################################################################################'')')
  write(6,'(''#  GULP called as a subroutine : Configuration 2 '')')
  write(6,'(''################################################################################'')')
!
!  Change some coordinates for example.....
!
  xyz(1,2) = 0.02_dp
  xyz(2,2) = 2.80_dp
  xyz(3,2) = 2.83_dp
!
!  ...and perhaps the cell parameters....
!
  cell(2,2) = 5.66_dp
!
!  Recompute the energy....
!
  call gulp_energy(cell,xyz,charges,energy,gradients,strainderivatives,lgradients,ierror)
!
  write(6,'(/,''  Cell: '')') 
  do i = 1,ndim
    write(6,'(i4,3(1x,f12.6))') i,(cell(j,i),j=1,3)
  enddo
  write(6,'(/,''  Coordinates: '')') 
  do i = 1,natoms
    write(6,'(i4,3(1x,f12.6))') i,(xyz(j,i),j=1,3)
  enddo
  write(6,'(/,''  Energy = '',f12.6,'' eV'')') energy
  write(6,'(/,''  Charges : '')') 
  do i = 1,natoms
    write(6,'(i4,1x,f12.6)') i,charges(i)
  enddo
  if (lgradients) then
    write(6,'(/,''  Gradients : eV/Ang '')') 
    do i = 1,natoms
      write(6,'(i4,3(1x,f12.6))') i,(gradients(j,i),j=1,3)
    enddo
    write(6,'(/,''  Strain gradients : eV '')') 
    write(6,'(4x,6(1x,f12.6))') (strainderivatives(j),j=1,6)
  endif
  write(6,'(/)')

end program main
