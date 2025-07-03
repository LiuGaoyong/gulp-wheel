  subroutine mctran(mode,ntran2try)
!
!  MC routine for transforming atoms. 
!
!  mode = if mode = 1, choose atoms to transform
!         if mode = 2, then create new trial transformation
!         if mode = 3, then undo previous transformation
!
!  10/21 Created from mcswap
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use current
  use general
  use genetic,       only : iseed
  use montecarlo
  use parallel
  use reallocate
  use species
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                  :: mode
  integer(i4), intent(in)                  :: ntran2try    ! Which if the nmctrans to use
!
!  Local variables
!
  integer(i4),                        save :: ntran1 = 0
  integer(i4),                        save :: ntran1spec = 0
  integer(i4),                        save :: ntran2spec = 0
!
  integer(i4)                              :: ns
  integer(i4)                              :: nt
  integer(i4)                              :: nspec1
  real(dp)                                 :: randnum
  real(dp)                                 :: GULP_random
#ifdef TRACE
  call trace_in('mctran')
#endif
!
  if (mode.eq.3) then
!**************************************
!  Mode 3 : Undo last transformation  *
!**************************************
!
!  Check that ntran1 has been set, otherwise there is nothing to undo
!
    if (ntran1.eq.0) return
!
!  Swap back species properties
!
    nt = nptrtranable(ntran1,ntran2try)
    ns = nmctranspec(ntran1spec,ntran2try)
    iatn(nt) = natspec(ns)
    nat(nt) = natspec(ns)
    natype(nt) = ntypspec(ns)
    nftype(nt) = ntypspec(ns)
    qf(nt) = qlspec(ns)
    rada(nt) = radspec(ns)
    radf(nt) = radspec(ns)
    nspecptr(nt) = ns
  elseif (mode.eq.1) then
!********************************
!  Mode 1 : New transformation  *
!********************************
!
!  Choose atom to transform
!
    randnum = GULP_random(iseed,1_i4)
    ntran1 = ntranable(ntran2try)*randnum + 1_i4
    if (ntran1.gt.ntranable(ntran2try)) ntran1 = ntranable(ntran2try)
!
!  Set trial atom pointer
!
    ntrialatom = 1
    nt = nptrtranable(ntran1,ntran2try)
    nptrtrialatom(ntrialatom) = nt
  elseif (mode.eq.2) then
!**********************************
!  Mode 2 : Apply transformation  *
!**********************************
!
!  Apply transformation to configuration array
!
    nt = nptrtranable(ntran1,ntran2try)
    nspec1 = nspecptr(nt)
    if (nmctranspec(1,ntran2try).eq.nspec1) then
      ntran1spec = 1
      ntran2spec = 2
    elseif (nmctranspec(2,ntran2try).eq.nspec1) then
      ntran1spec = 2
      ntran2spec = 1
    else
      call outerror('error has occured in Monte Carlo transformation',0_i4)
      call stopnow('mctran')
    endif
!
!  Swap species properties
!
    ns = nmctranspec(ntran2spec,ntran2try)
    iatn(nt) = natspec(ns)
    nat(nt) = natspec(ns)
    natype(nt) = ntypspec(ns)
    nftype(nt) = ntypspec(ns)
    qf(nt) = qlspec(ns)
    rada(nt) = radspec(ns)
    radf(nt) = radspec(ns)
    nspecptr(nt) = ns
  endif
#ifdef TRACE
  call trace_out('mctran')
#endif
!
  return
  end
