  subroutine settranptr(ntran2try)
!
!  Sets up pointer to atoms that can be transformed
!
!  10/21 Created from setswapptr
!   4/22 Wildcard flag passed to lmatch
!   6/23 Referencing of nspecptr corrected to use asymmetric unit
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
!  Julian Gale, CIC, Curtin University, June 2023
!
  use current
  use montecarlo
  use molecule,     only : natmol
  use reallocate
  use species
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in) :: ntran2try
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: ns
  integer(i4)              :: status
  logical,            save :: lfirstime = .true.
  logical                  :: lfound1
  logical                  :: lfound2
  logical                  :: lmatch
!
!  Initialise memory on first call
!
  if (lfirstime) then
    maxtranable = numat
    allocate(nptrtranable(maxtranable,maxmctrans),stat=status)
    if (status/=0) call outofmemory('settranptr','nptrtranable')
    lfirstime = .false.
  endif
!
  ntranable(ntran2try) = 0
!
!  Find species that match the current transformation
!
  lfound1 = .false.
  ns = 0
  do while (.not.lfound1.and.ns.lt.nspec) 
    ns = ns + 1
    lfound1 = lmatch(natspec(ns),ntypspec(ns),nmctrannat(1,ntran2try),nmctrantype(1,ntran2try),.false.)
  enddo
  if (lfound1) then
    nmctranspec(1,ntran2try) = ns
  else
    nmctranspec(1,ntran2try) = 0
  endif
  lfound2 = .false.
  ns = 0
  do while (.not.lfound2.and.ns.lt.nspec) 
    ns = ns + 1
    lfound2 = lmatch(natspec(ns),ntypspec(ns),nmctrannat(2,ntran2try),nmctrantype(2,ntran2try),.false.)
  enddo
  if (lfound2) then
    nmctranspec(2,ntran2try) = ns
  else
    nmctranspec(2,ntran2try) = 0
  endif
  if (lfound1.and.lfound2) then
!
!  Find atoms that match one of the species
!
    do i = 1,numat
      if (natmol(i).eq.0) then
        if (nspecptr(nrelf2a(i)).eq.nmctranspec(1,ntran2try)) then
          ntranable(ntran2try) = ntranable(ntran2try) + 1
          nptrtranable(ntranable(ntran2try),ntran2try) = i
        elseif (nspecptr(nrelf2a(i)).eq.nmctranspec(2,ntran2try)) then
          ntranable(ntran2try) = ntranable(ntran2try) + 1
          nptrtranable(ntranable(ntran2try),ntran2try) = i
        endif
      endif
    enddo
  endif
!
  return
  end
