  subroutine readreaxfflib
!
!  Reads a ffield file for the original reaxff code
!
!   4/23 Created from utility routine
!   5/23 ReaxFF torsion wildcard handling added
!   8/23 Correction to species handling for torsions
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
  use datatypes
  use configurations
  use control
  use element,          only : atsym
  use g_constants,      only : kcaltoev
  use reaxFFdata
  use species
  implicit none
!
!  Local variables
!
  character(len=132)              :: line
  character(len=2)                :: symbol
  integer(i4)                     :: i
  integer(i4)                     :: ind
  integer(i4)                     :: ind1
  integer(i4)                     :: ind2
  integer(i4)                     :: irxin
  integer(i4)                     :: j
  integer(i4)                     :: n
  integer(i4)                     :: n1
  integer(i4)                     :: n2
  integer(i4)                     :: n3
  integer(i4)                     :: n4
  integer(i4)                     :: nx
  integer(i4)                     :: nval3
  integer(i4)                     :: nangle
  integer(i4)                     :: nbond
  integer(i4)                     :: nhb
  integer(i4)                     :: nod
  integer(i4)                     :: npar
  integer(i4)                     :: ninspec(4)
  integer(i4)                     :: nrxspec
  integer(i4), allocatable        :: nrx2specptr(:)
  integer(i4), allocatable        :: nspec2rxptr(:)
  integer(i4)                     :: ntors
  logical                         :: lfound
  real(dp),    allocatable        :: reaxff_par(:)
  real(dp)                        :: dummy1
  real(dp)                        :: dummy2
  real(dp)                        :: dummy3
  real(dp)                        :: dummy4
  real(dp)                        :: reaxff1_angle(2)
  real(dp)                        :: reaxff_chi
  real(dp)                        :: reaxff_mu
  real(dp)                        :: reaxff_gamma
  real(dp)                        :: reaxff1_under
  real(dp)                        :: reaxff1_lonepair(2)
  real(dp)                        :: reaxff1_morse(4)
  real(dp)                        :: reaxff1_over(4)
  real(dp)                        :: reaxff1_radii(3)
  real(dp)                        :: reaxff1_valence(4)
  real(dp)                        :: reaxff2_bo(9)
  real(dp)                        :: reaxff2_bond(5)
  real(dp)                        :: reaxff2_morse(6)
  real(dp)                        :: reaxff2_over
  real(dp)                        :: reaxff3_angle(6)
  real(dp)                        :: reaxff3_hbond(4)
  real(dp)                        :: reaxff3_conj
  real(dp)                        :: reaxff3_penalty
  real(dp)                        :: reaxff4_torsion(5)
  real(dp)                        :: mass
  real(dp)                        :: units
!
!  Check file has been specifed
!
  if (.not.lreaxff_ffield) return
  if (reaxff_ffield(1:1).eq.' ') then
    call outerror('reaxff library file not specified',0_i4)
    call stopnow('readreaxfflib')
  endif
!
!  Set units as ReaxFF uses kcal
!
  units = kcaltoev
!
!  Open ReaxFF library file
!
  irxin = 91
  open(irxin,file=reaxff_ffield,status='old',form='formatted',err=100)
!
!  Read ReaxFF file
!
!  Dummy line
  read(irxin,'(a)',err=200,end=200) line
!***************
!  Parameters  *
!***************
  read(irxin,*) npar
  allocate(reaxff_par(max(npar,39_i4)))
  reaxff_par = 0.0_dp
!
!  Check parameter numbers
!
  if (npar.ne.39) then
    call outwarning('number of reaxff parameters is not 39',0_i4)
  endif
!
!  Read parameters
!
  do i = 1,npar
    read(irxin,*) reaxff_par(i)
  enddo
!
!  Assign to GULP parameters
!
  reaxFFcutoffVDW = reaxff_par(13)
  reaxFFcutoff = reaxff_par(13)
  reaxFFtol = 0.01_dp*reaxff_par(30)
  reaxFFlam(1) = reaxff_par(1)
  reaxFFlam(2) = reaxff_par(2)
!
  reaxFFlam(6)  = reaxff_par(33)
  reaxFFlam(31) = reaxff_par(32)
  reaxFFlam(7)  = reaxff_par(7)
  reaxFFlam(8)  = reaxff_par(9)
  reaxFFlam(9)  = reaxff_par(10)
!
  reaxFFlam(15) = reaxff_par(15)
  reaxFFlam(16) = reaxff_par(34)
  reaxFFlam(17) = reaxff_par(17)
  reaxFFlam(18) = reaxff_par(18)
!
  reaxFFlam(20) = reaxff_par(20)
  reaxFFlam(21) = reaxff_par(21)
  reaxFFlam(22) = reaxff_par(22)
!
  reaxFFlam(24) = reaxff_par(24)
  reaxFFlam(25) = reaxff_par(25)
  reaxFFlam(26) = reaxff_par(26)
  reaxFFlam(27) = reaxff_par(28)
!
  reaxFFlam(28) = reaxff_par(29)
!
  reaxFFlam(29) = reaxff_par(16)
!************************
!  Species information  *
!************************
!
!  Read number of species in library
!
  read(irxin,*) nrxspec
!
  allocate(nrx2specptr(nrxspec))
  allocate(nspec2rxptr(nspec))
!
!  3 x dummy line
!
  read(irxin,'(a)') line
  read(irxin,'(a)') line
  read(irxin,'(a)') line
!
!  Assign species 0 to be X
!
  nx = 0
!
  nspec2rxptr(1:nspec) = 0
!
!  Loop over species
!
  do n = 1,nrxspec
    nrx2specptr(n) = 0
!
!  First line
!
    read(irxin,*) symbol,reaxff1_radii(1),reaxff1_valence(1),mass,reaxff1_morse(3),reaxff1_morse(2), &
      reaxff_gamma,reaxff1_radii(2),reaxff1_valence(3)
!
!  Find whether this species is present in the system
!
    lfound = .false.
    j = 0
    do while (.not.lfound.and.j.lt.nspec) 
      j = j + 1
      lfound = (atsym(natspec(j)).eq.symbol)
    enddo
    if (lfound) then
      nspec2rxptr(j) = n
    endif
!
!  Compute terms dependent on this line
!
    reaxff1_lonepair(1) = 0.5_dp*(reaxff1_valence(3) - reaxff1_valence(1))
!
!  Second line
!
    read(irxin,*) reaxff1_morse(1),reaxff1_morse(4),reaxff1_valence(4),reaxff1_under,dummy1,reaxff_chi,reaxff_mu,dummy2
!
!  Third line
!
    read(irxin,*) reaxff1_radii(3),reaxff1_lonepair(2),dummy1,reaxff1_over(2),reaxff1_over(1),reaxff1_over(3),dummy2,dummy3
!
!  Fourth line
!
    read(irxin,*) reaxff1_over(4),reaxff1_angle(1),dummy1,reaxff1_valence(2),reaxff1_angle(2),dummy2,dummy3,dummy4
!
!  Set pointer to X species if present
!
    if (symbol.eq.'X ') nx = n
!
    if (lfound.and.n.ne.nx) then
      nreaxFFspec = nreaxFFspec + 1
      nrx2specptr(n) = nreaxFFspec
      if (nreaxFFspec.gt.maxreaxFFspec) then
        maxreaxFFspec = nreaxFFspec
        call changemaxreaxFFspec
      endif
      natreaxFFspec(nreaxFFspec) = natspec(j)
      ntypreaxFFspec(nreaxFFspec) = 0
      if (natreaxFFspec(nreaxFFspec).le.10) lreaxFFunder(nreaxFFspec) = .true.
      symbolreaxFFspec(nreaxFFspec) = symbol
      reaxFFr(1,nreaxFFspec) = reaxff1_radii(1)
      reaxFFr(2,nreaxFFspec) = reaxff1_radii(2)
      reaxFFr(3,nreaxFFspec) = reaxff1_radii(3)
      reaxFFval(1,nreaxFFspec) = reaxff1_valence(1)
      reaxFFval(2,nreaxFFspec) = reaxff1_valence(2)
      reaxFFval(3,nreaxFFspec) = reaxff1_valence(3)
      reaxFFval(4,nreaxFFspec) = reaxff1_valence(4)
      reaxFFpboc(1,nreaxFFspec) = reaxff1_over(1)
      reaxFFpboc(2,nreaxFFspec) = reaxff1_over(2)
      reaxFFpboc(3,nreaxFFspec) = reaxff1_over(3)
      reaxFFoc1(nreaxFFspec) = reaxff1_over(4)
      reaxFFuc1(nreaxFFspec) = reaxff1_under*units
      reaxFFlp(1,nreaxFFspec) = reaxff1_lonepair(1)
      reaxFFlp(2,nreaxFFspec) = reaxff1_lonepair(2)*units
      reaxFFlp(3,nreaxFFspec) = 0.0_dp   ! Parameter below is not currently used
      reaxFFval1(1,nreaxFFspec) = reaxff1_angle(1)
      reaxFFval1(2,nreaxFFspec) = reaxff1_angle(2)
      reaxFFalpha(nreaxFFspec) = reaxff1_morse(1)
      reaxFFeps(nreaxFFspec) = reaxff1_morse(2)*units
      reaxFFrvdw(nreaxFFspec) = reaxff1_morse(3)
      reaxFFgammaw(nreaxFFspec) = reaxff1_morse(4)
      reaxFFchi(nreaxFFspec) = reaxff_chi               ! NB: This term is already in eV
      reaxFFgamma(nreaxFFspec) = reaxff_gamma
      reaxFFmu(nreaxFFspec) = reaxff_mu                 ! NB: This term is already in eV
    endif
  enddo
!
!  Check if any species are missing
!
  do n = 1,nspec
    if (nspec2rxptr(n).eq.0) then
      call outerror('species is missing from reaxff library file',0_i4)
      call stopnow('readreaxfflib')
    endif
  enddo
!
!  Read number of bonds
!
  read(irxin,*) nbond
!
!  Dummy line
!
  read(irxin,'(a)') line
!
!  Read bond parameters
!
  do n = 1,nbond
!
!  Read first line
!
    read(irxin,*) ninspec(1),ninspec(2),reaxff2_bond(1),reaxff2_bond(2),reaxff2_bond(3), &
      reaxff2_bond(4),reaxff2_bo(5),reaxff2_bo(7),reaxff2_bo(6),reaxff2_over
!
!  Read second line
!
    read(irxin,*) reaxff2_bond(5),reaxff2_bo(3),reaxff2_bo(4),dummy1,reaxff2_bo(1),reaxff2_bo(2), &
      reaxff2_bo(8),reaxff2_bo(9)
!
!  If reaxff2_bo(3) = 1 needs to be set to 0 for GULP since this is a dummy value
!
    if (abs(reaxff2_bo(3)-1.0_dp).lt.1.0d-12) reaxff2_bo(3) = 0.0_dp
!
    n1 = nrx2specptr(ninspec(1))
    n2 = nrx2specptr(ninspec(2))
    if (n1.gt.0.and.n2.gt.0) then
      if (n1.ge.n2) then
        ind = n1*(n1-1)/2 + n2
      else
        ind = n2*(n2-1)/2 + n1
      endif
      lreaxFFpboOK(ind) = .true.
!
      reaxFFpbo(1,ind) = reaxff2_bo(1)
      reaxFFpbo(2,ind) = reaxff2_bo(2)
      reaxFFpbo(3,ind) = reaxff2_bo(3)
      reaxFFpbo(4,ind) = reaxff2_bo(4)
      reaxFFpbo(5,ind) = reaxff2_bo(5)
      reaxFFpbo(6,ind) = reaxff2_bo(6)
!
      lreaxFFbocorrect(1,ind) = (reaxff2_bo(8).gt.0.001_dp)
      lreaxFFbocorrect(2,ind) = (reaxff2_bo(7).gt.0.001_dp)
!
      reaxFFDe(1,ind)  = reaxff2_bond(1)*units
      reaxFFDe(2,ind)  = reaxff2_bond(2)*units
      reaxFFDe(3,ind)  = reaxff2_bond(3)*units
      reaxFFpbe(1,ind) = reaxff2_bond(4)
      reaxFFpbe(2,ind) = reaxff2_bond(5)
!
      reaxFFoc2(ind)  = reaxff2_over           ! NB: This term is already in eV
!
      reaxFFpen2(1,ind) = reaxff2_bo(9)*units
      reaxFFpen2(2,ind) = reaxff_par(14)
      reaxFFpen2(3,ind) = 1.0_dp
    endif
  enddo
!
!  Read number of off-diagonal terms
!
  read(irxin,*) nod
!
!  Read off-diagonal parameters
!
  do n = 1,nod
    read(irxin,*) ninspec(1),ninspec(2),reaxff2_morse(1),reaxff2_morse(3),reaxff2_morse(2), &
      reaxff2_morse(4),reaxff2_morse(5),reaxff2_morse(6)
!
    n1 = nrx2specptr(ninspec(1))
    n2 = nrx2specptr(ninspec(2))
    if (n1.gt.0.and.n2.gt.0) then
      if (n1.ge.n2) then
        ind = n1*(n1-1)/2 + n2
      else
        ind = n2*(n2-1)/2 + n1
      endif
      lreaxFFmorseinput(ind) = .true.
      reaxFFmorse(1,ind)  = reaxff2_morse(1)*units
      reaxFFmorse(2,ind)  = reaxff2_morse(2)
      reaxFFmorse(3,ind)  = reaxff2_morse(3)
      reaxFFmorse(4,ind)  = reaxff2_morse(4)
      reaxFFmorse(5,ind)  = reaxff2_morse(5)
      reaxFFmorse(6,ind)  = reaxff2_morse(6)
    endif
  enddo
!
!  Read number of angle terms
!
  read(irxin,*) nangle
!
!  Read in angle and penalty parameters
!
  do n = 1,nangle
    read(irxin,*) ninspec(2),ninspec(1),ninspec(3),reaxff3_angle(1),reaxff3_angle(2), &
      reaxff3_angle(3),reaxff3_conj,reaxff3_angle(5),reaxff3_penalty,reaxff3_angle(4)
    n1 = nrx2specptr(ninspec(1))
    n2 = nrx2specptr(ninspec(2))
    n3 = nrx2specptr(ninspec(3))
    if (n1.gt.0.and.n2.gt.0.and.n3.gt.0) then
      if (n2.ge.n3) then
        ind = n2*(n2-1)/2 + n3
      else
        ind = n3*(n3-1)/2 + n2
      endif
      nreaxFFval3(ind,n1) = nreaxFFval3(ind,n1) + 1
      if (nreaxFFval3(ind,n1).gt.maxreaxFFval3) then
        maxreaxFFval3 = nreaxFFval3(ind,n1)
        call changemaxreaxFFval3
      endif
      nval3 = nreaxFFval3(ind,n1)
!
      reaxFFval3(1,nval3,ind,n1) = reaxff3_angle(1)
      reaxFFval3(2,nval3,ind,n1) = reaxff3_angle(2)*units
      reaxFFval3(3,nval3,ind,n1) = reaxff3_angle(3)
      reaxFFval3(4,nval3,ind,n1) = reaxff3_angle(4)
      reaxFFval3(5,nval3,ind,n1) = reaxff3_angle(5)
      reaxFFval3(6,nval3,ind,n1) = reaxff3_angle(6)
!
      reaxFFpen3(ind,n1) = reaxff3_penalty*units
!
      reaxFFconj3(1,ind,n1) = reaxff3_conj*units
      reaxFFconj3(2,ind,n1) = reaxff_par(3)
      reaxFFconj3(3,ind,n1) = reaxff_par(39)
      reaxFFconj3(4,ind,n1) = reaxff_par(31)
    endif
  enddo
!
!  Read number of torsion terms
!
  read(irxin,*) ntors
!
!  Read torsion parameters
!
  do n = 1,ntors
    read(irxin,*) (ninspec(i),i=1,4),reaxff4_torsion(1),reaxff4_torsion(2),reaxff4_torsion(3), &
      reaxff4_torsion(4),reaxff4_torsion(5),dummy1,dummy2
    if (ninspec(1).eq.0.or.ninspec(4).eq.0) then
!
!  Wildcard torsion 
!
      n2 = nrx2specptr(ninspec(2))
      n3 = nrx2specptr(ninspec(3))
      if (n2.gt.0.and.n3.gt.0) then
        ind2 = 1
        if (n2.ge.n3) then
          ind1 = n2*(n2-1)/2 + n3
        else
          ind1 = n3*(n3-1)/2 + n2
        endif
        reaxFFtor4(1,ind2,ind1) = reaxff4_torsion(1)*units
        reaxFFtor4(2,ind2,ind1) = reaxff4_torsion(2)*units
        reaxFFtor4(3,ind2,ind1) = reaxff4_torsion(3)*units
        reaxFFtor4(4,ind2,ind1) = reaxff4_torsion(4)
        reaxFFtor4(5,ind2,ind1) = reaxff4_torsion(5)*units
      endif
    else
!
!  Specific torsion
!
      n1 = nrx2specptr(ninspec(1))
      n2 = nrx2specptr(ninspec(2))
      n3 = nrx2specptr(ninspec(3))
      n4 = nrx2specptr(ninspec(4))
      if (n1.gt.0.and.n2.gt.0.and.n3.gt.0.and.n4.gt.0) then
        if (n2.ge.n3) then
          ind1 = n2*(n2-1)/2 + n3
        else
          ind1 = n3*(n3-1)/2 + n2
        endif
        if (n1.ge.n4) then
          ind2 = n1*(n1-1)/2 + n4
        else
          ind2 = n4*(n4-1)/2 + n1
        endif
        reaxFFtor4(1,ind2,ind1) = reaxff4_torsion(1)*units
        reaxFFtor4(2,ind2,ind1) = reaxff4_torsion(2)*units
        reaxFFtor4(3,ind2,ind1) = reaxff4_torsion(3)*units
        reaxFFtor4(4,ind2,ind1) = reaxff4_torsion(4)
        reaxFFtor4(5,ind2,ind1) = reaxff4_torsion(5)*units
      endif
    endif
  enddo
!
!  Read number of hydrogen bond terms
!
  read(irxin,*) nhb
!
!  Read hydrogen bond parameters
!
  do n = 1,nhb
    read(irxin,*) ninspec(2),ninspec(1),ninspec(3),(reaxff3_hbond(i),i=1,4)
    n1 = nrx2specptr(ninspec(1))
    n2 = nrx2specptr(ninspec(2))
    n3 = nrx2specptr(ninspec(3))
    if (n1.gt.0.and.n2.gt.0.and.n3.gt.0) then
      reaxFFhb3(1,n3,n2,n1) = reaxff3_hbond(1)
      reaxFFhb3(2,n3,n2,n1) = reaxff3_hbond(2)*units
      reaxFFhb3(3,n3,n2,n1) = reaxff3_hbond(3)
      reaxFFhb3(4,n3,n2,n1) = reaxff3_hbond(4)
    endif
  enddo
!
!  Deallocate memory
!
  deallocate(nspec2rxptr)
  deallocate(nrx2specptr)
  deallocate(reaxff_par)
  return
!
!  Error traps
!
100 continue
  call outerror('reaxff library file not found',0_i4)
  call stopnow('readreaxfflib')
200 continue
  call outerror('error reading reaxff library file',0_i4)
  call stopnow('readreaxfflib')

end subroutine readreaxfflib
