  subroutine specword(nru,word,lwordok,iline,line,l55,l1000)
!
!  Processes input for species
!
!  nru = input fortran channel
!
!   6/01 Initialisation of line added for benefit of some compilers
!  10/02 Handling of missing charge case added
!  11/04 linspec handling modified to be species-wise
!  10/05 Modified to ensure that original case of library symbols is
!        preserved
!  11/05 Logical added to indicate whether charge was input or not
!   5/06 Setting of individual species mass added
!  12/08 Module input renamed to gulpinput
!   4/13 Bug in breathing shell handling fixed
!   1/19 Worsy changed to use word length parameter
!   8/19 Trap for shell accessing rion
!  11/21 Modifications for TI added
!  12/21 ltianyqspec added
!   6/23 AA-CLP force field parameters added
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
  use control,       only : laaclp
  use element
  use gulpinput
  use gulp_lengths
  use library
  use m_aaclp
  use parallel
  use species
  implicit none
!
!  Passed variables
!
  character(len=maxwordlength) :: word
  character(len=maxwordlength) :: wordoriginal
  character(len=maxlinelength) :: line
  integer(i4)                  :: iline
  integer(i4)                  :: nru
  logical                      :: l55
  logical                      :: l1000
  logical                      :: lwordok
!
!  Local variables
!
  integer(i4)                  :: i
  integer(i4)                  :: inat
  integer(i4)                  :: ind2
  integer(i4)                  :: itype
  integer(i4)                  :: ns
  logical                      :: lbreathing
  logical                      :: lfound
  logical                      :: lsymbol
  logical                      :: ltib
  logical                      :: ltif
!
!  Initialise local variables
!
  if (index(word,'spec').eq.1) goto 100
  if (index(word,'aaclp_s').eq.1) goto 110
  return
!********************************************************
!  Species info - at the moment charge and radius only  *
!********************************************************
100 continue
!
!  Check for TI flags
!
  ltib = .false.
  ltif = .false.
  if (nword.gt.1) then
    i = 2
    do while (i.le.nword)
      if (index(words(i),'flam').eq.1) then
        ltif = .true.
      elseif (index(words(i),'blam').eq.1) then
        ltib = .true.
      endif
      i = i + 1
    enddo
  endif
!
!  Sanity check lambda options
!
  if (ltib.and.ltif) then
    call outerror('Lambda cannot be forward and backward for the same species',iline)
    call stopnow('specword')
  endif
!
!  If either ltib or ltif are true then set ltianyqspec
!
  if (ltif.or.ltib) ltianyqspec = .true.
105 line = '  '
  read(nru,'(a)',end=108) line
  iline = iline + 1
  call linepro(nru,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 105
!
!  Check for old fashion specification of number of atoms
!
  if (nword.eq.0.and.nfloat.eq.1) goto 105
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.false.)
    if (.not.lsymbol) then
      call stolc(word,maxwordlength)
      if (index(word,'end').eq.1) then
        lwordok = .true.
        return
      else
        l55 = .true.
        goto 108
      endif
    endif
  endif
!
!  Find atomic number and type number
!
  lbreathing = .false.
  if (nword.eq.0) then
    inat = int(floats(1))
    if (inat.gt.100) inat = inat - 100 + maxele
    itype = 0
  elseif (nword.eq.1) then
    call ltont(words(1),inat,itype)
  elseif (nword.ge.2) then
    call ltont(words(1),inat,itype)
    wordoriginal = words(2)
    call stolc(words(2),maxword)
    if (index(words(2),'bco').eq.1.or.index(words(2),'bsh').eq.1) then
      lbreathing = .true.
      ind2 = 2
    else
      ind2 = 1
    endif
    if (index(words(2),'she').eq.ind2) inat = inat + maxele
  endif
!
!  Check whether species has already been included in list
!
  lfound = .false.
  ns = 0
  do while (ns.lt.nspec.and..not.lfound) 
    ns = ns + 1
    lfound = (inat.eq.natspec(ns).and.itype.eq.ntypspec(ns))
  enddo
  if (.not.lfound) then
!
!  If this is a new species increment counter and check memory
!
    nspec = nspec + 1
    if (nspec.gt.maxspec) then
      maxspec = nspec + 20
      call changemaxspec
    endif
    ns = nspec
  endif
  linspec(ns) = .true.
  natspec(ns) = inat
  ntypspec(ns) = itype
  ltibqspec(ns) = ltib
  ltifqspec(ns) = ltif
  if (lbreathing) lbrspec(ns) = .true.
  if (nword.eq.0) then
    qlspec(ns) = floats(2)
    if (nfloat.gt.2) then
      radspec(ns) = floats(3)
    else
      radspec(ns) = 0.0_dp
    endif
    lqinspec(ns) = .true.
  elseif (nword.eq.1) then
    if (nfloat.ge.2) then
      qlspec(ns) = floats(1)
      radspec(ns) = floats(2)
      lqinspec(ns) = .true.
    elseif (nfloat.eq.1) then
      qlspec(ns) = floats(1)
      radspec(ns) = 0.0_dp
      lqinspec(ns) = .true.
    else
      qlspec(ns) = 0.0_dp
      radspec(ns) = 0.0_dp
      lqinspec(ns) = .false.
    endif
  elseif (nword.eq.2) then
    if (nfloat.ge.2) then
      qlspec(ns) = floats(1)
      radspec(ns) = floats(2)
      lqinspec(ns) = .true.
    else
      if (nfloat.eq.1) then
        qlspec(ns) = floats(1)
        lqinspec(ns) = .true.
      else
        qlspec(ns) = 0.0_dp
        lqinspec(ns) = .false.
      endif
      if (lbrspec(ns)) then
        if (inat.gt.maxele) then
          radspec(ns) = rion(inat-maxele)
        else
          radspec(ns) = rion(inat)
        endif
      else
        radspec(ns) = 0.0_dp
      endif
    endif
    if (index(words(2),'cor').eq.0.and.index(words(2),'she').eq.0) then
!
!  Must be library symbol
!
      libspec(ns) = wordoriginal(1:16)
    endif
  elseif (nword.ge.3) then
    if (nfloat.ge.2) then
      qlspec(ns) = floats(1)
      radspec(ns) = floats(2)
      lqinspec(ns) = .true.
    else
      if (nfloat.eq.1) then
        qlspec(ns) = floats(1)
        lqinspec(ns) = .true.
      else
        qlspec(ns) = 0.0_dp
        lqinspec(ns) = .false.
      endif
      if (lbrspec(ns)) then
        if (inat.gt.maxele) then
          radspec(ns) = rion(inat-maxele)
        else
          radspec(ns) = rion(inat)
        endif
      else
        radspec(ns) = 0.0_dp
      endif
    endif
    libspec(ns) = words(3)(1:16)
  endif
!
!  Set species dependent mass
!
  if (.not.lmassinspec(ns)) then
    if (inat.le.maxele) then
      massspec(ns) = atmass(inat)
    else
      massspec(ns) = 0.0_dp
    endif
  endif
!
  goto 105
108 if (.not.l55) l1000 = .true.
  lwordok = .true.
  return
!***********************************************
!  AA-CLP force field information for species  *
!***********************************************
110 continue
  laaclp = .true.
  if (nfloat.gt.0) then
!
!  AA-CLP coefficients
!
    if (nfloat.ge.4) then
      aaclp_qscale = abs(floats(1))
      aaclp_dscale = abs(floats(2))
      aaclp_pscale = abs(floats(3))
      aaclp_rscale = abs(floats(4))
    elseif (nfloat.eq.3) then
      aaclp_qscale = abs(floats(1))
      aaclp_dscale = abs(floats(2))
      aaclp_pscale = abs(floats(3))
    elseif (nfloat.eq.2) then
      aaclp_qscale = abs(floats(1))
      aaclp_dscale = abs(floats(2))
    elseif (nfloat.eq.1) then
      aaclp_qscale = abs(floats(1))
    endif
  endif
115 line = '  '
  read(nru,'(a)',end=118) line
  iline = iline + 1
  call linepro(nru,line,iline)
!
!  Check whether input is symbol or option
!
  if ((nword+nfloat).eq.0) goto 115
!
!  Check for old fashion specification of number of atoms
!
  if (nword.eq.0.and.nfloat.eq.1) goto 115
  if (nword.gt.0) then
    word = words(1)
    call worsy(word,lsymbol,.false.)
    if (.not.lsymbol) then
      call stolc(word,maxwordlength)
      if (index(word,'end').eq.1) then
        lwordok = .true.
        return
      else
        l55 = .true.
        goto 118
      endif
    endif
  endif
!
!  Find atomic number and type number
!
  if (nword.eq.0) then
    inat = int(floats(1))
    if (inat.gt.100) inat = inat - 100 + maxele
    itype = 0
  elseif (nword.eq.1) then
    call ltont(words(1),inat,itype)
  elseif (nword.ge.2) then
    call ltont(words(1),inat,itype)
    wordoriginal = words(2)
    call stolc(words(2),maxword)
    if (index(words(2),'she').eq.1) inat = inat + maxele
  endif
!
!  Check whether species has already been included in list
!
  lfound = .false.
  ns = 0
  do while (ns.lt.nspec.and..not.lfound) 
    ns = ns + 1
    lfound = (inat.eq.natspec(ns).and.itype.eq.ntypspec(ns))
  enddo
  if (.not.lfound) then
!
!  If this is a new species increment counter and check memory
!
    nspec = nspec + 1
    if (nspec.gt.maxspec) then
      maxspec = nspec + 20
      call changemaxspec
    endif
    ns = nspec
  endif
  linspec(ns) = .true.
  natspec(ns) = inat
  ntypspec(ns) = itype
  ltibqspec(ns) = .false.
  ltifqspec(ns) = .false.
  lbrspec(ns) = .false.
  qlspec(ns) = 0.0_dp
  radspec(ns) = 0.0_dp
!
  if (nfloat.ge.4) then
    aaclp_pol_spec(ns) = abs(floats(1))
    aaclp_ip_spec(ns) = abs(floats(2))
    aaclp_b_spec(ns) = abs(floats(3))
    aaclp_h_spec(ns) = floats(4)
  else
    call outerror('Insufficient parameters given for AA-CLP species',iline)
    call stopnow('genword')
  endif
!
  goto 115
118 if (.not.l55) l1000 = .true.
  lwordok = .true.
  return
!
  end
