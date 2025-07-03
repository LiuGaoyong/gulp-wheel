program madelung
!
!  Compute Madelung constant for dielectrics in the presence of an anisotropic dielectric constant
!  Uses the equation given by S.T. Murphy and N.D.M. Hine, PRB, 87, 094111 (2013)
!
!  Conditions of use:
!
!  Madelung is available free of charge to academic institutions
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
!  Julian Gale, Curtin University, April 2021
!
  implicit none

  character(len=132) :: line
  integer*4          :: i
  integer*4          :: ifail
  integer*4          :: j
  integer*4          :: k
  integer*4          :: nreal
  integer*4          :: nrecip
  integer*4          :: ndim
  logical            :: converged
  real*8             :: angstoev
  real*8             :: dielectric
  real*8             :: kv(3,3)
  real*8             :: rv(3,3)
  real*8             :: eps(3,3)
  real*8             :: epsinv(3,3)
  real*8             :: erfc
  real*8             :: eta
  real*8             :: deteps
  real*8             :: ebgd
  real*8             :: elast
  real*8             :: ereal
  real*8             :: erecip
  real*8             :: eself
  real*8             :: etot
  real*8             :: g
  real*8             :: g2
  real*8             :: g1x
  real*8             :: g1y
  real*8             :: g1z
  real*8             :: g2x
  real*8             :: g2y
  real*8             :: g2z
  real*8             :: g3x
  real*8             :: g3y
  real*8             :: g3z
  real*8             :: pi
  real*8             :: q
  real*8             :: q2
  real*8             :: r
  real*8             :: r1x
  real*8             :: r1y
  real*8             :: r1z
  real*8             :: r2x
  real*8             :: r2y
  real*8             :: r2z
  real*8             :: r3x
  real*8             :: r3y
  real*8             :: r3z
  real*8             :: r2
  real*8             :: seta
  real*8             :: volume
  real*8             :: x
  real*8             :: y
  real*8             :: z
  real*8             :: tmp(3)
  real*8             :: wrk(6)
!
!  Set constants
!
  ndim = 3
  pi = 4.0d0*atan(1.0d0)
!********************
!  Read input data  *
!********************
!
!  Read cell
!
  read(5,'(a)') line
  read(5,*) rv(1,1),rv(2,1),rv(3,1)
  read(5,*) rv(1,2),rv(2,2),rv(3,2)
  read(5,*) rv(1,3),rv(2,3),rv(3,3)
!
!  Read dielectric constant
!
  read(5,'(a)') line
  read(5,*) eps(1,1),eps(2,1),eps(3,1)
  read(5,*) eps(1,2),eps(2,2),eps(3,2)
  read(5,*) eps(1,3),eps(2,3),eps(3,3)
!********************
!  Output           *
!********************
  write(6,'(/,''  Madeling calculation with anisotropic dielectric constant : '')')
  write(6,'(/,''  Lattice vectors (Angstroms): '',/)')
  write(6,'(2x,3f12.6)') rv(1,1),rv(2,1),rv(3,1)
  write(6,'(2x,3f12.6)') rv(1,2),rv(2,2),rv(3,2)
  write(6,'(2x,3f12.6)') rv(1,3),rv(2,3),rv(3,3)
  write(6,'(/,''  Dielectric tensor: '',/)')
  write(6,'(2x,3f12.6)') eps(1,1),eps(2,1),eps(3,1)
  write(6,'(2x,3f12.6)') eps(1,2),eps(2,2),eps(3,2)
  write(6,'(2x,3f12.6)') eps(1,3),eps(2,3),eps(3,3)
!
!  Charge
!
  read(5,*) q
  write(6,'(/,''  Charge = '',f12.6,'' a.u.'')') q
  q2 = q*q
!
!  Set initial values
!
  nreal = 4
  nrecip = 4
  eta = 0.1d0
!
  seta = sqrt(eta)
  write(6,'(/,''  Eta    = '',f12.6,'' Angstroms**-2 '')') eta
!
!  Compute inverse dielectric constant matrix
!
  epsinv(1:3,1:3) = eps(1:3,1:3)
  call matrix_inversion(epsinv,ndim,ndim,wrk,ifail)
!
!  Compute determinant of dielectric constant tensor
!
  deteps = eps(1,1)*(eps(2,2)*eps(3,3) - eps(3,2)*eps(2,3)) + &
           eps(1,2)*(eps(2,3)*eps(3,1) - eps(3,3)*eps(2,1)) + &
           eps(1,3)*(eps(2,1)*eps(3,2) - eps(3,1)*eps(2,2))
!************************************
!  Compute self energy              *
!************************************
  eself = - 2.0d0*seta/(sqrt(pi)*deteps**(1.0d0/2.0d0))
!************************************
!  Compute real space energy        *
!************************************
  elast = 0.0d0
  converged = .false.
  do while (.not.converged) 
!
!  Loop over lattice images
!
    ereal = 0.0d0
    do i = -nreal,nreal
      r1x = dble(i)*rv(1,1)
      r1y = dble(i)*rv(2,1)
      r1z = dble(i)*rv(3,1)
      do j = -nreal,nreal
        r2x = dble(j)*rv(1,2)
        r2y = dble(j)*rv(2,2)
        r2z = dble(j)*rv(3,2)
        do k = -nreal,nreal
          r3x = dble(k)*rv(1,3)
          r3y = dble(k)*rv(2,3)
          r3z = dble(k)*rv(3,3)
          x = r1x + r2x + r3x
          y = r1y + r2y + r3y
          z = r1z + r2z + r3z
          tmp(1) = epsinv(1,1)*x + epsinv(1,2)*y + epsinv(1,3)*z
          tmp(2) = epsinv(2,1)*x + epsinv(2,2)*y + epsinv(2,3)*z
          tmp(3) = epsinv(3,1)*x + epsinv(3,2)*y + epsinv(3,3)*z
          r2 = tmp(1)*x + tmp(2)*y + tmp(3)*z
          if (r2.gt.0.0001d0) then
            r = sqrt(r2)
            ereal = ereal + erfc(seta*r)/r
          endif
        enddo
      enddo
    enddo
    converged = (abs(ereal-elast).lt.1.0d-6) 
    elast = ereal
    nreal = nreal + 1
  enddo
! DEBUG
!  write(6,'(/,''  Nreal   = '',i5)') nreal
  ereal = ereal/sqrt(deteps)
!************************************
!  Compute reciprocal space energy  *
!************************************
  erecip = 0.0d0
!
!  Determine reciprocal lattice vectors
!
  kv(1:3,1:3) = rv(1:3,1:3)
  call matrix_inversion(kv,ndim,ndim,wrk,ifail)
  kv = 2.0d0*pi*kv
  write(6,'(/,''  Reciprocal lattice vectors (1/Angstroms): '',/)')
  write(6,'(2x,3f12.6)') kv(1,1),kv(2,1),kv(3,1)
  write(6,'(2x,3f12.6)') kv(1,2),kv(2,2),kv(3,2)
  write(6,'(2x,3f12.6)') kv(1,3),kv(2,3),kv(3,3)
!
!  Determine volume
!
  volume = rv(1,1)*(rv(2,2)*rv(3,3) - rv(3,2)*rv(2,3)) + &
           rv(1,2)*(rv(2,3)*rv(3,1) - rv(3,3)*rv(2,1)) + &
           rv(1,3)*(rv(2,1)*rv(3,2) - rv(3,1)*rv(2,2))
  write(6,'(/,''  Volume = '',f12.3,'' Angstroms**3 '',/)') volume
!
  elast = 0.0d0
  converged = .false.
  do while (.not.converged) 
!
!  Loop over reciprocal space images
!
    erecip = 0.0d0
    do i = -nrecip,nrecip
      g1x = dble(i)*kv(1,1)
      g1y = dble(i)*kv(2,1)
      g1z = dble(i)*kv(3,1)
      do j = -nrecip,nrecip
        g2x = dble(j)*kv(1,2)
        g2y = dble(j)*kv(2,2)
        g2z = dble(j)*kv(3,2)
        do k = -nrecip,nrecip
          g3x = dble(k)*kv(1,3)
          g3y = dble(k)*kv(2,3)
          g3z = dble(k)*kv(3,3)
          x = g1x + g2x + g3x
          y = g1y + g2y + g3y
          z = g1z + g2z + g3z
          tmp(1) = eps(1,1)*x + eps(1,2)*y + eps(1,3)*z
          tmp(2) = eps(2,1)*x + eps(2,2)*y + eps(2,3)*z
          tmp(3) = eps(3,1)*x + eps(3,2)*y + eps(3,3)*z
          g2 = tmp(1)*x + tmp(2)*y + tmp(3)*z
          if (g2.gt.0.0001d0) then
            erecip = erecip + exp(-g2/(4.0d0*eta))/g2
          endif
        enddo
      enddo
    enddo
    converged = (abs(erecip-elast).lt.1.0d-6) 
    elast = erecip
    nrecip = nrecip + 1
  enddo
! DEBUG
!  write(6,'(/,''  Nrecip  = '',i5)') nrecip
!
!  Multiply by constants
!
  erecip = 4.0d0*pi*erecip/volume
!************************************
!  Compute background energy        *
!************************************
  ebgd = - pi/(eta*volume)
!
!  Convert units and multiply by half
!
  angstoev = 14.3997584d0
  ebgd  = ebgd*angstoev*0.5d0*q2
  eself = eself*angstoev*0.5d0*q2
  ereal = ereal*angstoev*0.5d0*q2
  erecip = erecip*angstoev*0.5d0*q2
!
  etot = ereal + eself + ebgd + erecip
!
  write(6,'('' Energy contributions to interaction : '',/)')
  write(6,'('' Energy - background = '',f12.6)') ebgd
  write(6,'('' Energy - self       = '',f12.6)') eself
  write(6,'('' Energy - real       = '',f12.6)') ereal
  write(6,'('' Energy - recip      = '',f12.6)') erecip
  write(6,'('' Energy - TOTAL      = '',f12.6)') etot
  write(6,'(/,'' Energy correction   = '',f12.6,'' eV'',/)') -etot

end
