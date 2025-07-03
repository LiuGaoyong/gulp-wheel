  subroutine outcfg(etot,lgrad1)
!
!  Subroutine for generating CFG file readable by Ovito and AtomEye.
!
!   10/21 Created
!   10/21 Shell handling added
!    9/23 H0 transposed for correct output
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
!  Roman Groger, Czech Acad. Sci.
!  Julian Gale, CIC, Curtin University, September 2023
!
  use current
  use derivatives,     only : xdrv, ydrv, zdrv
  use element
  use energies,        only : siteenergy
  use gulp_files
  use parallel
  use shells
  use symmetry
  implicit none
!
!  Passed variables
!
  real(dp), intent(in)    :: etot
  logical,  intent(in)    :: lgrad1
!
!  Local variables
!
  integer(i4)             :: i
  integer(i4)             :: iaux
  integer(i4)             :: ii
  integer(i4)             :: naux
  integer(i4)             :: ifail
  integer(i4),       save :: iout = 10
  logical                 :: lcfg_floc
  logical                 :: lcfg_fxloc
  logical                 :: lcfg_fyloc
  logical                 :: lcfg_fzloc
  real(dp)                :: oldat
  real(dp)                :: rvi(3,3)
  real(dp)                :: s(3)
  real(dp)                :: wrk(6)
!
!  If not I/O proc then return
!
  if (.not.ioproc) return
!
!  Set local variants of flags depending on whether lgrad1 is true or not
!
  if (lgrad1) then
    lcfg_floc  = lcfg_f
    lcfg_fxloc = lcfg_fx
    lcfg_fyloc = lcfg_fy
    lcfg_fzloc = lcfg_fz
  else
    lcfg_floc  = .false.
    lcfg_fxloc = .false.
    lcfg_fyloc = .false.
    lcfg_fzloc = .false.
  endif
!
!  If name has been given then open file
!
  open(iout,file=cfgfile,status='replace')
! 
!  Total energy as a comment
!
  write(iout,'("# Energy = ",e20.12," eV")') etot
!   
!  Number of particles and lattice parameter
!
  write(iout,'("Number of particles = ",i10)') ncore
  write(iout,'("A = 1.0 Angstrom (basic length-scale)")')
! 
!  Cell vectors
!
  write(iout,'("H0(1,1) = ",e20.12, " A")') rv(1,1)
  write(iout,'("H0(1,2) = ",e20.12, " A")') rv(2,1)
  write(iout,'("H0(1,3) = ",e20.12, " A")') rv(3,1)
  write(iout,'("H0(2,1) = ",e20.12, " A")') rv(1,2)
  write(iout,'("H0(2,2) = ",e20.12, " A")') rv(2,2)
  write(iout,'("H0(2,3) = ",e20.12, " A")') rv(3,2)
  write(iout,'("H0(3,1) = ",e20.12, " A")') rv(1,3)
  write(iout,'("H0(3,2) = ",e20.12, " A")') rv(2,3)
  write(iout,'("H0(3,3) = ",e20.12, " A")') rv(3,3)
! 
!  No velocities
!
  write(iout,'(".NO_VELOCITY.")')
!
!  Auxiliary fields
!
  naux = count((/ lcfg_e, lcfg_fxloc, lcfg_fyloc, lcfg_fzloc, lcfg_floc, lcfg_q, lcfg_dq /))
  write(iout,'("entry_count = ",i2)') 3+naux
  iaux = 0
  if (lcfg_e) then
    write(iout,'("auxiliary[",i2,"] = energy")') iaux
    iaux = iaux + 1
  endif
  if (lcfg_fxloc) then
    write(iout,'("auxiliary[",i2,"] = fx")') iaux
    iaux = iaux + 1
  endif
  if (lcfg_fyloc) then
    write(iout,'("auxiliary[",i2,"] = fy")') iaux
    iaux = iaux + 1
  endif
  if (lcfg_fzloc) then
    write(iout,'("auxiliary[",i2,"] = fz")') iaux
    iaux = iaux + 1
  endif
  if (lcfg_floc) then
    write(iout,'("auxiliary[",i2,"] = f")') iaux
    iaux = iaux + 1
  endif
  if (lcfg_q) then
    write(iout,'("auxiliary[",i2,"] = q")') iaux
    iaux = iaux + 1
  endif
  if (lcfg_dq) then
    write(iout,'("auxiliary[",i2,"] = dq")') iaux
    iaux = iaux + 1
  endif
! 
!  Properties of atoms
!
  oldat = -1
  rvi = rv
  call matrix_inversion(rvi, 3_i4, 3_i4, wrk, ifail)
  do i = 1,ncore
    ii = ncoptr(i)
    if (nat(ii).ne.oldat) then
      write(iout,"(f6.2,/,a)") atmass(nat(ii)), atsym(nat(ii))
      oldat = nat(ii)
    endif
!
!  Mandatory fields (position of atom)
!
    s = matmul(rvi,(/ xclat(ii), yclat(ii), zclat(ii) /))
    write(iout,'(3(1x,f15.12),$)') s
!
!  Auxiliary fields
!
    if (lcfg_e) write(iout,'(1x,f16.12,$)') siteenergy(ii)
    if (lcfg_fxloc) write(iout,'(1x,f16.12,$)') -xdrv(ii)
    if (lcfg_fyloc) write(iout,'(1x,f16.12,$)') -ydrv(ii)
    if (lcfg_fzloc) write(iout,'(1x,f16.12,$)') -zdrv(ii)
    if (lcfg_floc) write(iout,'(1x,f16.12,$)') sqrt(xdrv(ii)**2 + ydrv(ii)**2 + zdrv(ii)**2)
    if (lcfg_q) write(iout,'(1x,f12.8,$)') qf(ii)
    if (lcfg_dq) write(iout,'(1x,f12.8,$)') qf(ii)-q0(nat(ii))
    write(iout,*)
  enddo

  close(iout)

end subroutine outcfg
