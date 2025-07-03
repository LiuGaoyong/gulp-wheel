  subroutine outjdftx(iout)
!
!  Write out a file with the structure in the form suitable
!  for the JDFTx code.
!
!   2/22 Created 
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
!  Copyright Curtin University 2022
!
!  Julian Gale, CIC, Curtin University, February 2022
!
  use configurations
  use current
  use element
  use g_constants,       only : autoangs
  use gulp_files
  use general
  use shells
  use species
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                      :: iout
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ii
  real(dp)                                     :: xcent
  real(dp)                                     :: ycent
  real(dp)                                     :: zcent
  real(dp)                                     :: xmax
  real(dp)                                     :: ymax
  real(dp)                                     :: zmax
  real(dp)                                     :: xmin
  real(dp)                                     :: ymin
  real(dp)                                     :: zmin
!
!  If file name has been given then open file
!
  if (jdftxfile(1:1).ne.' ') then
    open(iout,file=jdftxfile,status='unknown')
  endif
!
  if (ndim.eq.3) then
!*******************
!  Bulk structure  *
!*******************
!
!  Lattice
!
    write(iout,'(''lattice \'')')
    write(iout,'(2x,3f12.6,'' \'')') (rv(i,1)/autoangs,i=1,3)
    write(iout,'(2x,3f12.6,'' \'')') (rv(i,2)/autoangs,i=1,3)
    write(iout,'(2x,3f12.6)') (rv(i,3)/autoangs,i=1,3)
!
!  Coordinates of cores
!
    write(iout,'(/,''coords-type Cartesian'')')
    do ii = 1,ncore
      i = ncoptr(ii)
      write(iout,'(''ion '',a2,3(1x,f12.6),'' 1'')') atsym(nat(i)),xclat(i)/autoangs,yclat(i)/autoangs, &
        zclat(i)/autoangs
    enddo
    write(iout,'(/)')
  elseif (ndim.eq.2) then
!****************
!  Slab output  *
!****************
!
!  Find centre of system
!
    zcent = 0.0_dp
    do ii = 1,ncore
      i = ncoptr(ii)
      zcent = zcent + zclat(i)
    enddo
    zcent = zcent/dble(ncore)
!
!  Find maximum extent in surface normal direction
!
    zmax = 0.0_dp
    zmin = 0.0_dp
    do ii = 1,ncore
      i = ncoptr(ii)
      zmax = max(zmax,zclat(i))
      zmin = min(zmin,zclat(i))
    enddo
!
!  Dummy lattice - size + 10 Ang
!
    write(iout,'(''lattice \'')')
    write(iout,'(2x,3f12.6,'' \'')') rv(1,1)/autoangs,rv(2,1)/autoangs,0.0_dp
    write(iout,'(2x,3f12.6,'' \'')') rv(1,2)/autoangs,rv(2,2)/autoangs,0.0_dp
    write(iout,'(2x,3f12.6)') 0.0_dp,0.0_dp,(zmax-zmin+10.0_dp)/autoangs
!
!  Coordinates of cores
!
    write(iout,'(/,''coords-type Cartesian'')')
    do ii = 1,ncore
      i = ncoptr(ii)
      write(iout,'(''ion '',a2,3(1x,f12.6),'' 1'')') atsym(nat(i)),xclat(i)/autoangs,&
        yclat(i)/autoangs,(zclat(i)-zcent)/autoangs
    enddo
    write(iout,'(/)')
    write(iout,'(''coulomb-interaction slab 001'')')
    write(iout,'(''coulomb-truncation-embed 0.0 0.0 0.0'')')
    write(iout,'(/)')
  elseif (ndim.eq.1) then
!*******************
!  Polymer output  *
!*******************
!
!  Find centre of system
!
    ycent = 0.0_dp
    zcent = 0.0_dp
    do ii = 1,ncore
      i = ncoptr(ii)
      ycent = ycent + yclat(i)
      zcent = zcent + zclat(i)
    enddo
    ycent = ycent/dble(ncore)
    zcent = zcent/dble(ncore)
!
!  Find maximum extent in each direction
!
    ymax = 0.0_dp
    zmax = 0.0_dp
    ymin = 0.0_dp
    zmin = 0.0_dp
    do ii = 1,ncore
      i = ncoptr(ii)
      ymax = max(ymax,yclat(i))
      zmax = max(zmax,zclat(i))
      ymin = min(ymin,yclat(i))
      zmin = min(zmin,zclat(i))
    enddo
!
!  Dummy lattice - size + 10 Ang
!
    write(iout,'(''lattice \'')')
    write(iout,'(2x,3f12.6,'' \'')') rv(1,1)/autoangs,0.0_dp,0.0_dp
    write(iout,'(2x,3f12.6,'' \'')') 0.0_dp,(ymax-ymin+10.0_dp)/autoangs,0.0_dp
    write(iout,'(2x,3f12.6)') 0.0_dp,0.0_dp,(zmax-zmin+10.0_dp)/autoangs
!
!  Coordinates of cores
!
    write(iout,'(/,''coords-type Cartesian'')')
    do ii = 1,ncore
      i = ncoptr(ii)
      write(iout,'(''ion '',a2,3(1x,f12.6),'' 1'')') atsym(nat(i)),xclat(i)/autoangs,&
        (yclat(i)-ycent)/autoangs,(zclat(i)-zcent)/autoangs
    enddo
    write(iout,'(/)')
    write(iout,'(''coulomb-interaction wire 100'')')
    write(iout,'(''coulomb-truncation-embed 0.0 0.0 0.0'')')
    write(iout,'(/)')
  elseif (ndim.eq.0) then
!*******************
!  Cluster output  *
!*******************
!
!  Find centre of system
!
    xcent = 0.0_dp
    ycent = 0.0_dp
    zcent = 0.0_dp
    do ii = 1,ncore
      i = ncoptr(ii)
      xcent = xcent + xclat(i)
      ycent = ycent + yclat(i)
      zcent = zcent + zclat(i)
    enddo
    xcent = xcent/dble(ncore)
    ycent = ycent/dble(ncore)
    zcent = zcent/dble(ncore)
!
!  Find maximum extent in each direction
!
    xmax = 0.0_dp
    ymax = 0.0_dp
    zmax = 0.0_dp
    xmin = 0.0_dp
    ymin = 0.0_dp
    zmin = 0.0_dp
    do ii = 1,ncore
      i = ncoptr(ii)
      xmax = max(xmax,xclat(i))
      ymax = max(ymax,yclat(i))
      zmax = max(zmax,zclat(i))
      xmin = min(xmin,xclat(i))
      ymin = min(ymin,yclat(i))
      zmin = min(zmin,zclat(i))
    enddo
!
!  Dummy lattice - size + 10 Ang
!
    write(iout,'(''lattice \'')')
    write(iout,'(2x,3f12.6,'' \'')') (xmax-xmin+10.0_dp)/autoangs,0.0_dp,0.0_dp
    write(iout,'(2x,3f12.6,'' \'')') 0.0_dp,(ymax-ymin+10.0_dp)/autoangs,0.0_dp
    write(iout,'(2x,3f12.6)') 0.0_dp,0.0_dp,(zmax-zmin+10.0_dp)/autoangs
!
!  Coordinates of cores
!
    write(iout,'(/,''coords-type Cartesian'')')
    do ii = 1,ncore
      i = ncoptr(ii)
      write(iout,'(''ion '',a2,3(1x,f12.6),'' 1'')') atsym(nat(i)),(xclat(i)-xcent)/autoangs,& 
        (yclat(i)-ycent)/autoangs,(zclat(i)-zcent)/autoangs
    enddo
    write(iout,'(/)')
    write(iout,'(''coulomb-interaction isolated'')')
    write(iout,'(''coulomb-truncation-embed 0.0 0.0 0.0'')')
    write(iout,'(/)')
  endif
!
  return
  end
