  subroutine getpotsymbol1(iline,llibrary,nvar1,itype1,sym1,nword1,nbeg,lvalidpot)
!
!  Gets species symbols for one-body input
!
!   7/06 Created from getpotsymbol2
!   7/06 Error in non-library call corrected
!   8/06 ltype0 flag added to okspec call
!   9/06 Literal symbols now returned 
!  11/06 Words are saved for literal symbols before corruption by okspec call
!  12/08 Module input renamed to gulpinput
!   3/09 Breathing shells now handled when setting atomic number
!   6/09 Module name changed from three to m_three
!   2/18 Trace added
!   1/19 maxwordlength changes added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, January 2019
!
  use datatypes  
  use element,       only : maxele
  use gulpinput,     only : nfloat, nword, floats, words
  use gulp_lengths
  use species,       only : natspec, ntypspec, nspec
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  character(len=5), intent(out) :: sym1
  integer(i4),      intent(in)  :: iline
  integer(i4),      intent(out) :: itype1
  integer(i4),      intent(out) :: nbeg
  integer(i4),      intent(out) :: nvar1
  integer(i4),      intent(in)  :: nword1
  logical,          intent(in)  :: llibrary
  logical,          intent(out) :: lvalidpot
!
!  Local variables
!
  character(len=maxwordlength)  :: wordsave
  integer(i4)                   :: ilp1
  logical                       :: lok1
  logical                       :: ltype0
#ifdef TRACE
  call trace_in('getpotsymbol1')
#endif
!
  lvalidpot = .true.
  nbeg = 0
  sym1 = ' '
  if (nword.gt.0) then
!
!  Symbols used in input
!
    wordsave = words(nword1+1)
    if (nword.eq.nword1+1) then
      if (llibrary) then
        call okspec(lok1,words(nword1+1),ilp1,.true.,ltype0)
        if (.not.lok1) then
          lvalidpot = .false.
        endif
        if (ilp1.gt.0.and.ilp1.le.nspec) then
          nvar1 = natspec(ilp1)
          if (ltype0) then
            itype1 = 0
          else
            itype1 = ntypspec(ilp1)
          endif
        elseif (ilp1.eq.-1) then
          nvar1 = maxele
          itype1 = 0
        endif
      else
        call ltont(words(nword1+1),nvar1,itype1)
      endif
    elseif (nword.eq.nword1+2) then
      if (llibrary) then
        call okspec(lok1,words(nword1+1),ilp1,.true.,ltype0)
        if (.not.lok1) then
          lvalidpot = .false.
        endif
        if (ilp1.gt.0.and.ilp1.le.nspec) then
          nvar1 = natspec(ilp1)
          if (ltype0) then
            itype1 = 0
          else
            itype1 = ntypspec(ilp1)
          endif
        else 
          nvar1 = maxele
          itype1 = 0
        endif
      else
        call ltont(words(nword1+1),nvar1,itype1)
      endif
      if ((index(words(nword1+2),'s').eq.1).or.(index(words(nword1+2),'S').eq.1)) then
        nvar1 = nvar1 + maxele
      elseif ((index(words(nword1+2),'bs').eq.1).or.(index(words(nword1+2),'BS').eq.1)) then
        nvar1 = nvar1 + maxele
      elseif ((index(words(nword1+2),'bS').eq.1).or.(index(words(nword1+2),'Bs').eq.1)) then
        nvar1 = nvar1 + maxele
      endif
    else
      call outerror('Incorrect species input for one-body term',iline)
      call stopnow('getpotsymbol1')
    endif
    nbeg = 0
    sym1 = wordsave(1:5)
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    itype1 = 0
    nbeg = 1
    nfloat = nfloat - 1
    call label(nvar1,itype1,sym1)
  endif
#ifdef TRACE
  call trace_out('getpotsymbol1')
#endif
!
  return
  end
!
  subroutine getpotsymbol2(iline,llibrary,nvar1,itype1,sym1,nvar2,itype2,sym2,nbeg,lvalidpot)
!
!  Gets species symbols for two-body potential input
!
!  10/04 Created from potword21
!  10/05 Correction to handling of atom type added
!   8/06 ltype0 flag added to okspec call
!   9/06 Literal symbols now returned 
!  11/06 Words are saved for literal symbols before corruption by okspec call
!  12/08 Module input renamed to gulpinput
!   3/09 Breathing shells now handled when setting atomic number
!   5/11 Incorrect reference to ilp1 instead of ilp2 corrected
!   2/18 Trace added
!   1/19 maxwordlength changes added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, January 2019
!
  use datatypes  
  use element,       only : maxele
  use gulpinput,     only : nfloat, nword, floats, words
  use gulp_lengths
  use species,       only : natspec, ntypspec, nspec
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  character(len=5), intent(out) :: sym1
  character(len=5), intent(out) :: sym2
  integer(i4),      intent(in)  :: iline
  integer(i4),      intent(out) :: itype1
  integer(i4),      intent(out) :: itype2
  integer(i4),      intent(out) :: nbeg
  integer(i4),      intent(out) :: nvar1
  integer(i4),      intent(out) :: nvar2
  logical,          intent(in)  :: llibrary
  logical,          intent(out) :: lvalidpot
!
!  Local variables
!
  character(len=maxwordlength)  :: word
  character(len=maxwordlength)  :: wordsave1
  character(len=maxwordlength)  :: wordsave2
  integer(i4)                   :: ilp1
  integer(i4)                   :: ilp2
  logical                       :: lok1
  logical                       :: lok2
  logical                       :: ltype01
  logical                       :: ltype02
#ifdef TRACE
  call trace_in('getpotsymbol2')
#endif
!
  lvalidpot = .true.
  nbeg = 0
  sym1 = ' '
  sym2 = ' '
  if (nword.gt.0) then
!
!  Symbols used in input
!
    if (nword.eq.2) then
      wordsave1 = words(1)
      wordsave2 = words(2)
      if (llibrary) then
        call okspec(lok1,words(1),ilp1,.true.,ltype01)
        call okspec(lok2,words(2),ilp2,.true.,ltype02)
        if (.not.lok1.or..not.lok2) then
          lvalidpot = .false.
        endif
        if (ilp1.gt.0.and.ilp1.le.nspec) then
          nvar1 = natspec(ilp1)
          if (ltype01) then
            itype1 = 0
          else
            itype1 = ntypspec(ilp1)
          endif
        elseif (ilp1.eq.-1) then
          nvar1 = maxele
          itype1 = 0
        endif
        if (ilp2.gt.0.and.ilp2.le.nspec) then
          nvar2 = natspec(ilp2)
          if (ltype02) then
            itype2 = 0
          else
            itype2 = ntypspec(ilp2)
          endif
        elseif (ilp2.eq.-1) then
          nvar2 = maxele
          itype2 = 0
        endif
      else
        call ltont(words(1),nvar1,itype1)
        call ltont(words(2),nvar2,itype2)
      endif
      sym1 = wordsave1(1:5)
      sym2 = wordsave2(1:5)
    elseif (nword.eq.3) then
      wordsave1 = words(1)
      if (llibrary) then
        call okspec(lok1,words(1),ilp1,.true.,ltype01)
        if (.not.lok1) then
          lvalidpot = .false.
        endif
        if (ilp1.gt.0.and.ilp1.le.nspec) then
          nvar1 = natspec(ilp1)
          if (ltype01) then
            itype1 = 0
          else
            itype1 = ntypspec(ilp1)
          endif
        elseif (ilp1.eq.-1) then
          nvar1 = maxele
          itype1 = 0
        endif
      else
        call ltont(words(1),nvar1,itype1)
      endif
      sym1 = wordsave1(1:5)
      word = words(2)
      call stolc(word,maxwordlength)
      if (index(word,'cor').eq.1) then
        wordsave2 = words(3)
        if (llibrary) then
          call okspec(lok2,words(3),ilp2,.true.,ltype02)
          if (.not.lok2) then
            lvalidpot = .false.
          endif
          if (ilp2.gt.0.and.ilp2.le.nspec) then
            nvar2 = natspec(ilp2)
            if (ltype02) then
              itype2 = 0
            else
              itype2 = ntypspec(ilp2)
            endif
          else
            nvar2 = maxele
            itype2 = 0
          endif
        else
          call ltont(words(3),nvar2,itype2)
        endif
        sym2 = wordsave2(1:5)
      elseif (index(word,'she').eq.1) then
        nvar1 = nvar1 + maxele
        wordsave2 = words(3)
        if (llibrary) then
          call okspec(lok2,words(3),ilp2,.true.,ltype02)
          if (.not.lok2) then
            lvalidpot = .false.
          endif
          if (ilp2.gt.0.and.ilp2.le.nspec) then
            nvar2 = natspec(ilp2)
            if (ltype02) then
              itype2 = 0
            else
              itype2 = ntypspec(ilp2)
            endif
          else
            nvar2 = maxele
            itype2 = 0
          endif
        else
          call ltont(words(3),nvar2,itype2)
        endif
        sym2 = wordsave2(1:5)
      else
        wordsave2 = words(2)
        if (llibrary) then
          call okspec(lok2,words(2),ilp2,.true.,ltype02)
          if (.not.lok2) then
            lvalidpot = .false.
          endif
          if (ilp2.gt.0.and.ilp2.le.nspec) then
            nvar2 = natspec(ilp2)
            if (ltype02) then
              itype2 = 0
            else
              itype2 = ntypspec(ilp2)
            endif
          else
            nvar2 = maxele
            itype2 = 0
          endif
        else
          call ltont(words(2),nvar2,itype2)
        endif
        word = words(3)
        call stolc(word,maxwordlength)
        if (index(word,'she').eq.1) then
          nvar2 = nvar2 + maxele
        endif
        sym2 = wordsave2(1:5)
      endif
    elseif (nword.eq.4) then
      wordsave1 = words(1)
      wordsave2 = words(3)
      if (llibrary) then
        call okspec(lok1,words(1),ilp1,.true.,ltype01)
        call okspec(lok2,words(3),ilp2,.true.,ltype02)
        if (.not.lok1.or..not.lok2) then
          lvalidpot = .false.
        endif
        if (ilp1.gt.0.and.ilp1.le.nspec) then
          nvar1 = natspec(ilp1)
          if (ltype01) then
            itype1 = 0
          else
            itype1 = ntypspec(ilp1)
          endif
        else 
          nvar1 = maxele
          itype1 = 0
        endif
        if (ilp2.gt.0.and.ilp2.le.nspec) then
          nvar2 = natspec(ilp2)
          if (ltype02) then
            itype2 = 0
          else
            itype2 = ntypspec(ilp2)
          endif
        else
          nvar2 = maxele
          itype2 = 0
        endif
      else
        call ltont(words(1),nvar1,itype1)
        call ltont(words(3),nvar2,itype2)
      endif
      if ((index(words(2),'s').eq.1).or.(index(words(2),'S').eq.1)) then
        nvar1 = nvar1 + maxele
      elseif ((index(words(2),'bs').eq.1).or.(index(words(2),'BS').eq.1)) then
        nvar1 = nvar1 + maxele
      elseif ((index(words(2),'bS').eq.1).or.(index(words(2),'Bs').eq.1)) then
        nvar1 = nvar1 + maxele
      endif
      if ((index(words(4),'s').eq.1).or.(index(words(4),'S').eq.1)) then
        nvar2 = nvar2 + maxele
      elseif ((index(words(4),'bs').eq.1).or.(index(words(4),'BS').eq.1)) then
        nvar2 = nvar2 + maxele
      elseif ((index(words(4),'bS').eq.1).or.(index(words(4),'Bs').eq.1)) then
        nvar2 = nvar2 + maxele
      endif
      sym1 = wordsave1(1:5)
      sym2 = wordsave2(1:5)
    else
      call outerror('Incorrect species input for two-body potential',iline)
      call stopnow('getpotsymbol2')
    endif
    nbeg = 0
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    nvar2 = int(floats(2))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    if (nvar2.gt.100) nvar2 = nvar2 - 100 + maxele
    itype1 = 0
    itype2 = 0
    nbeg = 2
    nfloat = nfloat - 2
    call label(nvar1,itype1,sym1)
    call label(nvar2,itype2,sym2)
  endif
#ifdef TRACE
  call trace_out('getpotsymbol2')
#endif
!
  return
  end
!
  subroutine getpotsymbols2(iline,llibrary,maxpair,npair,nvars,itypes,syms,nbeg,lvalidpot)
!
!  Gets species symbols for two-body potential input.
!  NB: Library version that can return multiple species pairs
!
!   3/22 Created from getpotsymbol2
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
!  Julian Gale, CIC, Curtin University, March 2022
!
  use datatypes  
  use element,       only : maxele
  use gulpinput,     only : nfloat, nword, floats, words
  use gulp_lengths
  use species,       only : natspec, ntypspec, nspec
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),       intent(in)  :: maxpair
  character(len=5),  intent(out) :: syms(2,maxpair)
  integer(i4),       intent(in)  :: iline
  integer(i4),       intent(out) :: itypes(2,maxpair)
  integer(i4),       intent(out) :: nbeg
  integer(i4),       intent(out) :: npair
  integer(i4),       intent(out) :: nvars(2,maxpair)
  logical,           intent(in)  :: llibrary
  logical,           intent(out) :: lvalidpot
!
!  Local variables
!
  character(len=5)               :: lab1
  character(len=5)               :: lab2
  character(len=maxwordlength)   :: word
  integer(i4)                    :: ics
  integer(i4),              save :: maxilp = 20
  integer(i4)                    :: n
  integer(i4)                    :: n1
  integer(i4)                    :: n2
  integer(i4)                    :: nilp1
  integer(i4)                    :: nilp2
  integer(i4), allocatable, save :: ilps1(:)
  integer(i4), allocatable, save :: ilps2(:)
  integer(i4), allocatable, save :: itypes1(:)
  integer(i4), allocatable, save :: itypes2(:)
  integer(i4), allocatable, save :: nvars1(:)
  integer(i4), allocatable, save :: nvars2(:)
  logical                        :: lok1
  logical                        :: lok2
  logical                        :: lsame
  logical,     allocatable, save :: ltypes1(:)
  logical,     allocatable, save :: ltypes2(:)
#ifdef TRACE
  call trace_in('getpotsymbols2')
#endif
!
!  This routine is currently only designed for use with libraries and so check that this is the case
!
  if (.not.llibrary) then
    call outerror('getpotsymbols2 should only be called when reading a library',iline)
    call stopnow('getpotsymbols2')
  endif
!
  lvalidpot = .true.
  nbeg = 0
  npair = 0
!
!  Allocate local workspace
!
  allocate(ilps1(maxilp))
  allocate(ilps2(maxilp))
  allocate(itypes1(maxilp))
  allocate(itypes2(maxilp))
  allocate(ltypes1(maxilp))
  allocate(ltypes2(maxilp))
  allocate(nvars1(maxilp))
  allocate(nvars2(maxilp))
!
  if (nword.gt.0) then
!
!  Symbols used in input
!
    if (nword.ge.2) then
      call okspecs(lok1,words(1),maxilp,nilp1,ilps1,.true.,ltypes1)
!
!  Check dimensions of species arrays
!
      if (nilp1.gt.maxilp) then
        deallocate(nvars2)
        deallocate(nvars1)
        deallocate(ltypes2)
        deallocate(ltypes1)
        deallocate(itypes2)
        deallocate(itypes1)
        deallocate(ilps2)
        deallocate(ilps1)
!
        maxilp = nilp1
!
        allocate(ilps1(maxilp))
        allocate(ilps2(maxilp))
        allocate(itypes1(maxilp))
        allocate(itypes2(maxilp))
        allocate(ltypes1(maxilp))
        allocate(ltypes2(maxilp))
        allocate(nvars1(maxilp))
        allocate(nvars2(maxilp))
!
!  Repeat species calls with correct array sizes
!
        call okspecs(lok1,words(1),maxilp,nilp1,ilps1,.true.,ltypes1)
      endif
!
      if (.not.lok1) then
        lvalidpot = .false.
        goto 100
      endif
!
      if (nword.eq.2) then
        lsame = (words(1).eq.words(2))
        call okspecs(lok2,words(2),maxilp,nilp2,ilps2,.true.,ltypes2)
      elseif (nword.eq.3) then
!
!  Which is the other symbol?
!
        word = words(2)
        call stolc(word,maxwordlength)
        if (index(word,'co').eq.1.or.index(word,'bc').eq.1.or.index(word,'bc').eq.1.or.index(word,'bs').eq.1) then
          ics = 2
        else
          ics = 3
        endif
!
        lsame = (words(1).eq.words(5-ics))
        call okspecs(lok2,words(5-ics),maxilp,nilp2,ilps2,.true.,ltypes2)
      else
        lsame = (words(1).eq.words(3))
        call okspecs(lok2,words(3),maxilp,nilp2,ilps2,.true.,ltypes2)
      endif
      if (.not.lok2) then
        lvalidpot = .false.
        goto 100
      endif
!
!  Check dimensions of species arrays
!
      if (nilp2.gt.maxilp) then
        deallocate(nvars2)
        deallocate(nvars1)
        deallocate(ltypes2)
        deallocate(ltypes1)
        deallocate(itypes2)
        deallocate(itypes1)
        deallocate(ilps2)
        deallocate(ilps1)
!
        maxilp = nilp2
!
        allocate(ilps1(maxilp))
        allocate(ilps2(maxilp))
        allocate(itypes1(maxilp))
        allocate(itypes2(maxilp))
        allocate(ltypes1(maxilp))
        allocate(ltypes2(maxilp))
        allocate(nvars1(maxilp))
        allocate(nvars2(maxilp))
!
!  Repeat species calls with correct array sizes
!
        call okspecs(lok1,words(1),maxilp,nilp1,ilps1,.true.,ltypes1)
        if (nword.eq.2) then
          lsame = (words(1).eq.words(2))
          call okspecs(lok2,words(2),maxilp,nilp2,ilps2,.true.,ltypes2)
        elseif (nword.eq.3) then
!
!  Which is the other symbol?
!
          word = words(2)
          call stolc(word,maxwordlength)
          if (index(word,'co').eq.1.or.index(word,'bc').eq.1.or.index(word,'bc').eq.1.or.index(word,'bs').eq.1) then
            ics = 2
          else
            ics = 3
          endif
!
          lsame = (words(1).eq.words(5-ics))
          call okspecs(lok2,words(5-ics),maxilp,nilp2,ilps2,.true.,ltypes2)
        else
          lsame = (words(1).eq.words(3))
          call okspecs(lok2,words(3),maxilp,nilp2,ilps2,.true.,ltypes2)
        endif
      endif
!
      do n = 1,nilp1
        if (ilps1(n).gt.0.and.ilps1(n).le.nspec) then
          nvars1(n) = natspec(ilps1(n))
          if (ltypes1(n)) then
            itypes1(n) = 0
          else
            itypes1(n) = ntypspec(ilps1(n))
          endif
        elseif (ilps1(n).eq.-1) then
          nvars1(n) = maxele
          itypes1(n) = 0
        endif
      enddo
      do n = 1,nilp2
        if (ilps2(n).gt.0.and.ilps2(n).le.nspec) then
          nvars2(n) = natspec(ilps2(n))
          if (ltypes2(n)) then
            itypes2(n) = 0
          else
            itypes2(n) = ntypspec(ilps2(n))
          endif
        elseif (ilps2(n).eq.-1) then
          nvars2(n) = maxele
          itypes2(n) = 0
        endif
      enddo
      if (lsame) then
!
!  Case where both species are the same
!
        do n1 = 1,nilp1
          call label(nvars1(n1),itypes1(n1),lab1)
          do n2 = 1,n1
            npair = npair + 1
            if (npair.le.maxpair) then
              call label(nvars1(n2),itypes1(n2),lab2)
              syms(1,npair) = lab1
              syms(2,npair) = lab2
              nvars(1,npair) = nvars1(n1)
              nvars(2,npair) = nvars2(n2)
              itypes(1,npair) = itypes1(n1)
              itypes(2,npair) = itypes2(n2)
            endif
          enddo
        enddo
      else
!
!  Case where both species are different
!
        do n1 = 1,nilp1
          call label(nvars1(n1),itypes1(n1),lab1)
          do n2 = 1,nilp2
            npair = npair + 1
            if (npair.le.maxpair) then
              call label(nvars2(n2),itypes2(n2),lab2)
              syms(1,npair) = lab1
              syms(2,npair) = lab2
              nvars(1,npair) = nvars1(n1)
              nvars(2,npair) = nvars2(n2)
              itypes(1,npair) = itypes1(n1)
              itypes(2,npair) = itypes2(n2)
            endif
          enddo
        enddo
      endif
!
!  Handle shell types
!
      if (npair.le.maxpair) then
        if (nword.ge.4) then
          do n1 = 1,npair
            if ((index(words(2),'s').eq.1).or.(index(words(2),'S').eq.1)) then
              nvars(1,n1) = nvars(1,n1) + maxele
            elseif ((index(words(2),'bs').eq.1).or.(index(words(2),'BS').eq.1)) then
              nvars(1,n1) = nvars(1,n1) + maxele
            elseif ((index(words(2),'bS').eq.1).or.(index(words(2),'Bs').eq.1)) then
              nvars(1,n1) = nvars(1,n1) + maxele
            endif
            if ((index(words(4),'s').eq.1).or.(index(words(4),'S').eq.1)) then
              nvars(2,n1) = nvars(2,n1) + maxele
            elseif ((index(words(4),'bs').eq.1).or.(index(words(4),'BS').eq.1)) then
              nvars(2,n1) = nvars(2,n1) + maxele
            elseif ((index(words(4),'bS').eq.1).or.(index(words(4),'Bs').eq.1)) then
              nvars(2,n1) = nvars(2,n1) + maxele
            endif
          enddo
        elseif (nword.eq.3.and.ics.eq.2) then
          do n1 = 1,npair
            if ((index(words(2),'s').eq.1).or.(index(words(2),'S').eq.1)) then
              nvars(1,n1) = nvars(1,n1) + maxele
            elseif ((index(words(2),'bs').eq.1).or.(index(words(2),'BS').eq.1)) then
              nvars(1,n1) = nvars(1,n1) + maxele
            elseif ((index(words(2),'bS').eq.1).or.(index(words(2),'Bs').eq.1)) then
              nvars(1,n1) = nvars(1,n1) + maxele
            endif
          enddo
        elseif (nword.eq.3.and.ics.eq.3) then
          do n1 = 1,npair
            if ((index(words(3),'s').eq.1).or.(index(words(3),'S').eq.1)) then
              nvars(2,n1) = nvars(2,n1) + maxele
            elseif ((index(words(3),'bs').eq.1).or.(index(words(3),'BS').eq.1)) then
              nvars(2,n1) = nvars(2,n1) + maxele
            elseif ((index(words(3),'bS').eq.1).or.(index(words(3),'Bs').eq.1)) then
              nvars(2,n1) = nvars(2,n1) + maxele
            endif
          enddo
        endif
      endif
    endif
    nbeg = 0
  else
!
!  Numeric input
!
    npair = 1
    nvars(1,1) = int(floats(1))
    nvars(2,1) = int(floats(2))
    if (nvars(1,1).gt.100) nvars(1,1) = nvars(1,1) - 100 + maxele
    if (nvars(2,1).gt.100) nvars(2,1) = nvars(2,1) - 100 + maxele
    itypes(1,1) = 0
    itypes(2,1) = 0
    nbeg = 2
    nfloat = nfloat - 2
    call label(nvars(1,1),itypes(1,1),lab1)
    call label(nvars(2,1),itypes(2,1),lab2)
    syms(1,1) = lab1
    syms(2,1) = lab2
  endif
!
!  Return point to ensure memory is freed
!
100 continue
!
!  Free local workspace
!
  deallocate(nvars2)
  deallocate(nvars1)
  deallocate(ltypes2)
  deallocate(ltypes1)
  deallocate(itypes2)
  deallocate(itypes1)
  deallocate(ilps2)
  deallocate(ilps1)
#ifdef TRACE
  call trace_out('getpotsymbols2')
#endif
!
  return
  end
!
  subroutine getpotsymbol3(iline,llibrary,nvar1,itype1,sym1,nvar2,itype2,sym2,nvar3,itype3,sym3,nbeg,lvalidpot)
!
!  Gets species symbols for three-body potential input
!
!  10/04 Created from potword3
!  10/05 Handling of atom types modified
!   8/06 ltype0 flag added to okspec call
!   9/06 Literal symbols now returned 
!  11/06 Words are saved for literal symbols before corruption by okspec call
!  12/08 Module input renamed to gulpinput
!   3/10 Modified so that all three species types are return arguments
!   2/18 Trace added
!   1/19 maxwordlength changes added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, January 2019
!
  use datatypes  
  use element,       only : maxele
  use gulpinput,     only : nfloat, nword, floats, words
  use gulp_lengths
  use species,       only : natspec, ntypspec, nspec
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  character(len=5), intent(out) :: sym1
  character(len=5), intent(out) :: sym2
  character(len=5), intent(out) :: sym3
  integer(i4),      intent(in)  :: iline
  integer(i4),      intent(out) :: itype1
  integer(i4),      intent(out) :: itype2
  integer(i4),      intent(out) :: itype3
  integer(i4),      intent(out) :: nbeg
  integer(i4),      intent(out) :: nvar1
  integer(i4),      intent(out) :: nvar2
  integer(i4),      intent(out) :: nvar3
  logical,          intent(in)  :: llibrary
  logical,          intent(out) :: lvalidpot
!
!  Local variables
!
  character(len=maxwordlength)  :: word
  character(len=maxwordlength)  :: wordsave
  integer(i4)                   :: i
  integer(i4)                   :: ilp1
  integer(i4)                   :: ityp(3)
  integer(i4)                   :: nptr
  integer(i4)                   :: nvr(3)
  logical                       :: lok1
  logical                       :: ltype0
#ifdef TRACE
  call trace_in('getpotsymbol3')
#endif
!
  lvalidpot = .true.
  nbeg = 0
!
  if (nword.gt.0) then
    if (nword.ge.3) then
      nptr = 0
      do i = 1,nword
        wordsave = words(i)
        word = words(i)
        call stolc(word,maxwordlength)
        if (index(word,'she').eq.1.and.nptr.gt.0) then
          nvr(nptr) = nvr(nptr) + maxele
        elseif (index(word,'cor').eq.0) then
          nptr = nptr + 1
          if (llibrary) then
            call okspec(lok1,words(i),ilp1,.true.,ltype0)
            if (.not.lok1) then
              lvalidpot = .false.
            endif
            if (ilp1.gt.0.and.ilp1.le.nspec) then
              nvr(nptr) = natspec(ilp1)
              if (ltype0) then
                ityp(nptr) = 0
              else
                ityp(nptr) = ntypspec(ilp1)
              endif
            else 
              nvr(nptr) = maxele
              ityp(nptr) = 0
            endif
          else
            call ltont(words(i),nvr(nptr),ityp(nptr))
          endif
          if (nptr.eq.1) then
            sym1 = wordsave(1:5)
          elseif (nptr.eq.2) then
            sym2 = wordsave(1:5)
          elseif (nptr.eq.3) then
            sym3 = wordsave(1:5)
          endif
        endif
      enddo
      if (nptr.ne.3) then
        call outerror('Incorrect potential species input for three-body potential',iline)
        call stopnow('getpotsymbol3')
      endif
    else
      call outerror('Incorrect potential species input for three-body potential',iline)
      call stopnow('getpotsymbol3')
    endif
    nbeg = 0
    nvar1 = nvr(1)
    nvar2 = nvr(2)
    nvar3 = nvr(3)
    itype1 = ityp(1)
    itype2 = ityp(2)
    itype3 = ityp(3)
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    nvar2 = int(floats(2))
    nvar3 = int(floats(3))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    if (nvar2.gt.100) nvar2 = nvar2 - 100 + maxele
    if (nvar3.gt.100) nvar3 = nvar3 - 100 + maxele
    itype1 = 0
    itype2 = 0
    itype3 = 0
    nbeg = 3
    nfloat = nfloat - 3
    call label(nvar1,itype1,sym1)
    call label(nvar2,itype2,sym2)
    call label(nvar3,itype3,sym3)
  endif
#ifdef TRACE
  call trace_out('getpotsymbol3')
#endif
!
  return
  end
!
  subroutine getpotsymbols3(iline,llibrary,maxtrio,ntrio,nvars,itypes,syms,nbeg,lvalidpot)
!
!  Gets species symbols for three-body potential input
!  NB: Library version that can return multiple species trios
!
!   3/22 Created from getpotsymbol3
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
!  Julian Gale, CIC, Curtin University, March 2022
!
  use datatypes  
  use element,       only : maxele
  use gulpinput,     only : nfloat, nword, floats, words
  use gulp_lengths
  use species,       only : natspec, ntypspec, nspec
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),       intent(in)  :: maxtrio
  character(len=5),  intent(out) :: syms(3,maxtrio)
  integer(i4),       intent(in)  :: iline
  integer(i4),       intent(out) :: itypes(3,maxtrio)
  integer(i4),       intent(out) :: ntrio
  integer(i4),       intent(out) :: nbeg
  integer(i4),       intent(out) :: nvars(3,maxtrio)
  logical,           intent(in)  :: llibrary
  logical,           intent(out) :: lvalidpot
!
!  Local variables
!
  character(len=5)               :: lab1
  character(len=5)               :: lab2
  character(len=5)               :: lab3
  character(len=maxwordlength)   :: word
  character(len=maxwordlength)   :: wordsave
  character(len=maxwordlength)   :: wordsave2
  character(len=maxwordlength)   :: wordsave3
  integer(i4)                    :: i
  integer(i4),              save :: maxilp = 20
  integer(i4)                    :: n1
  integer(i4)                    :: n2
  integer(i4)                    :: n3
  integer(i4)                    :: nilp(3)
  integer(i4), allocatable, save :: ilp(:,:)
  integer(i4), allocatable, save :: ityp(:,:)
  integer(i4), allocatable, save :: nvr(:,:)
  integer(i4)                    :: nptr
  logical                        :: lok1
  logical                        :: lsame23
  logical,     allocatable, save :: ltyp(:,:)
#ifdef TRACE
  call trace_in('getpotsymbols3')
#endif
!
!  This routine is currently only designed for use with libraries and so check that this is the case
!
  if (.not.llibrary) then
    call outerror('getpotsymbols3 should only be called when reading a library',iline)
    call stopnow('getpotsymbols3')
  endif
!
  lvalidpot = .true.
  nbeg = 0
  ntrio = 0
!
!  Allocate local workspace
!
  allocate(ilp(maxilp,3))
  allocate(ityp(maxilp,3))
  allocate(ltyp(maxilp,3))
  allocate(nvr(maxilp,3))
!
  if (nword.gt.0) then
    if (nword.ge.3) then
      nptr = 0
      do i = 1,nword
        wordsave = words(i)
        word = words(i)
        call stolc(word,maxwordlength)
        if (index(word,'she').eq.1.and.nptr.gt.0) then
          nvr(1:nilp(nptr),nptr) = nvr(1:nilp(nptr),nptr) + maxele
        elseif (index(word,'cor').eq.0) then
          nptr = nptr + 1
          if (nptr.eq.2) then
            wordsave2 = words(i)
          elseif (nptr.eq.3) then
            wordsave3 = words(i)
            lsame23 = (wordsave2.eq.wordsave3)
          endif
          call okspecs(lok1,words(i),maxilp,nilp(nptr),ilp(1,nptr),.true.,ltyp(1,nptr))
!
!  Check dimensions of species arrays
!
          if (nilp(nptr).gt.maxilp) then
            deallocate(nvr)
            deallocate(ltyp)
            deallocate(ityp)
            deallocate(ilp)
!
            maxilp = nilp(nptr)
!
            allocate(ilp(maxilp,3))
            allocate(ityp(maxilp,3))
            allocate(ltyp(maxilp,3))
            allocate(nvr(maxilp,3))
!
!  Repeat species calls with correct array sizes
!
            call okspecs(lok1,words(i),maxilp,nilp(nptr),ilp(1,nptr),.true.,ltyp(1,nptr))
          endif
!
          if (.not.lok1) then
            lvalidpot = .false.
          endif
!
          if (nilp(nptr).gt.0) then
            do n1 = 1,nilp(nptr)
              if (ilp(n1,nptr).gt.0.and.ilp(n1,nptr).le.nspec) then
                nvr(n1,nptr) = natspec(ilp(n1,nptr))
                if (ltyp(n1,nptr)) then
                  ityp(n1,nptr) = 0
                else
                  ityp(n1,nptr) = ntypspec(ilp(n1,nptr))
                endif
              else 
                nvr(n1,nptr) = maxele
                ityp(n1,nptr) = 0
              endif
            enddo
          endif
        endif
      enddo
      if (nptr.ne.3) then
        call outerror('Incorrect potential species input for three-body potential',iline)
        call stopnow('getpotsymbols3')
      endif
    else
      call outerror('Incorrect potential species input for three-body potential',iline)
      call stopnow('getpotsymbols3')
    endif
    nbeg = 0
!
!  Generate trios
!
    do n1 = 1,nilp(1)
      call label(nvr(n1,1),ityp(n1,1),lab1)
      do n2 = 1,nilp(2)
        call label(nvr(n2,2),ityp(n2,2),lab2)
        if (lsame23) then
          do n3 = 1,n2
            ntrio = ntrio + 1
            if (ntrio.le.maxtrio) then
              call label(nvr(n3,3),ityp(n3,3),lab3)
              syms(1,ntrio) = lab1
              syms(2,ntrio) = lab2
              syms(3,ntrio) = lab3
              nvars(1,ntrio) = nvr(n1,1)
              nvars(2,ntrio) = nvr(n2,2)
              nvars(3,ntrio) = nvr(n3,3)
              itypes(1,ntrio) = ityp(n1,1)
              itypes(2,ntrio) = ityp(n2,2)
              itypes(3,ntrio) = ityp(n3,3)
            endif
          enddo
        else
          do n3 = 1,nilp(3)
            ntrio = ntrio + 1
            if (ntrio.le.maxtrio) then
              call label(nvr(n3,3),ityp(n3,3),lab3)
              syms(1,ntrio) = lab1
              syms(2,ntrio) = lab2
              syms(3,ntrio) = lab3
              nvars(1,ntrio) = nvr(n1,1)
              nvars(2,ntrio) = nvr(n2,2)
              nvars(3,ntrio) = nvr(n3,3)
              itypes(1,ntrio) = ityp(n1,1)
              itypes(2,ntrio) = ityp(n2,2)
              itypes(3,ntrio) = ityp(n3,3)
            endif
          enddo
        endif
      enddo
    enddo
  else
!
!  Numeric input
!
    ntrio = 1
    nvars(1,1) = int(floats(1))
    nvars(2,1) = int(floats(2))
    nvars(3,1) = int(floats(3))
    if (nvars(1,1).gt.100) nvars(1,1) = nvars(1,1) - 100 + maxele
    if (nvars(2,1).gt.100) nvars(2,1) = nvars(2,1) - 100 + maxele
    if (nvars(3,1).gt.100) nvars(3,1) = nvars(3,1) - 100 + maxele
    itypes(1,1) = 0
    itypes(2,1) = 0
    itypes(3,1) = 0
    nbeg = 3
    nfloat = nfloat - 3
    call label(nvars(1,1),itypes(1,1),syms(1,1))
    call label(nvars(2,1),itypes(2,1),syms(2,1))
    call label(nvars(3,1),itypes(3,1),syms(3,1))
  endif
!
!  Free local workspace
!
  deallocate(nvr)
  deallocate(ltyp)
  deallocate(ityp)
  deallocate(ilp)
#ifdef TRACE
  call trace_out('getpotsymbols3')
#endif
!
  return
  end
!
  subroutine getpotsymbol4(iline,llibrary,nvar1,itype1,sym1,nvar2,itype2,sym2, &
                           nvar3,itype3,sym3,nvar4,itype4,sym4,nbeg,lvalidpot)
!
!  Gets species symbols for four-body potential input
!
!  10/04 Created from potword4
!  10/05 Correction to handling of atom type added
!   8/06 ltype0 flag added to okspec call
!   9/06 Literal symbols now returned
!  11/06 Words are saved for literal symbols before corruption by okspec call
!  12/08 Module input renamed to gulpinput
!   2/18 Trace added
!   1/19 maxwordlength changes added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, January 2019
!
  use datatypes  
  use element,       only : maxele
  use gulpinput,     only : nfloat, nword, floats, words
  use gulp_lengths
  use species,       only : natspec, ntypspec, nspec
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  character(len=5), intent(out) :: sym1
  character(len=5), intent(out) :: sym2
  character(len=5), intent(out) :: sym3
  character(len=5), intent(out) :: sym4
  integer(i4),      intent(in)  :: iline
  integer(i4),      intent(out) :: itype1
  integer(i4),      intent(out) :: itype2
  integer(i4),      intent(out) :: itype3
  integer(i4),      intent(out) :: itype4
  integer(i4),      intent(out) :: nbeg
  integer(i4),      intent(out) :: nvar1
  integer(i4),      intent(out) :: nvar2
  integer(i4),      intent(out) :: nvar3
  integer(i4),      intent(out) :: nvar4
  logical,          intent(in)  :: llibrary
  logical,          intent(out) :: lvalidpot
!
!  Local variables
!
  character(len=maxwordlength)  :: word
  character(len=maxwordlength)  :: wordsave
  integer(i4)                   :: i
  integer(i4)                   :: ilp1
  integer(i4)                   :: ityp(4)
  integer(i4)                   :: nptr
  integer(i4)                   :: nvr(4)
  logical                       :: lok1
  logical                       :: ltype0
#ifdef TRACE
  call trace_in('getpotsymbol4')
#endif
!
  lvalidpot = .true.
  nbeg = 0
!
  if (nword.gt.0) then
    if (nword.ge.4) then
      nptr = 0
      do i = 1,nword
        wordsave = words(i)
        word = words(i)
        call stolc(word,maxwordlength)
        if (index(word,'she').eq.1.and.nptr.gt.0) then
          nvr(nptr) = nvr(nptr) + maxele
        elseif (index(word,'cor').eq.0) then
          nptr = nptr + 1
          if (llibrary) then
            call okspec(lok1,words(i),ilp1,.true.,ltype0)
            if (.not.lok1) then
              lvalidpot = .false.
            endif
            if (ilp1.gt.0.and.ilp1.le.nspec) then
              nvr(nptr) = natspec(ilp1)
              if (ltype0) then
                ityp(nptr) = 0
              else
                ityp(nptr) = ntypspec(ilp1)
              endif
            else 
              nvr(nptr) = maxele
              ityp(nptr) = 0
            endif
          else
            call ltont(words(i),nvr(nptr),ityp(nptr))
          endif
          if (nptr.eq.1) then
            sym1 = wordsave(1:5)
          elseif (nptr.eq.2) then
            sym2 = wordsave(1:5)
          elseif (nptr.eq.3) then
            sym3 = wordsave(1:5)
          elseif (nptr.eq.4) then
            sym4 = wordsave(1:5)
          endif
        endif
      enddo
      if (nptr.ne.4) then
        call outerror('Incorrect potential species input for four-body potential',iline)
        call stopnow('getpotsymbol4')
      endif
    else
      call outerror('Incorrect potential species input for four-body potential',iline)
      call stopnow('getpotsymbol4')
    endif
    nbeg = 0
    nvar1 = nvr(1)
    nvar2 = nvr(2)
    nvar3 = nvr(3)
    nvar4 = nvr(4)
    itype1 = ityp(1)
    itype2 = ityp(2)
    itype3 = ityp(3)
    itype4 = ityp(4)
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    nvar2 = int(floats(2))
    nvar3 = int(floats(3))
    nvar4 = int(floats(4))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    if (nvar2.gt.100) nvar2 = nvar2 - 100 + maxele
    if (nvar3.gt.100) nvar3 = nvar3 - 100 + maxele
    if (nvar4.gt.100) nvar4 = nvar4 - 100 + maxele
    itype1 = 0
    itype2 = 0
    itype3 = 0
    itype4 = 0
    nbeg = 4
    nfloat = nfloat - 4
    call label(nvar1,itype1,sym1)
    call label(nvar2,itype2,sym2)
    call label(nvar3,itype3,sym3)
    call label(nvar4,itype4,sym4)
  endif
#ifdef TRACE
  call trace_out('getpotsymbol4')
#endif
!
  return
  end
!
  subroutine getpotsymbols4(iline,llibrary,maxqrt,nqrt,nvars,itypes,syms,nbeg,lvalidpot)
!
!  Gets species symbols for four-body potential input
!  NB: Library version that can return multiple species quartets
!
!   3/22 Created from getpotsymbols3
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
!  Julian Gale, CIC, Curtin University, March 2022
!
  use datatypes  
  use element,       only : maxele
  use gulpinput,     only : nfloat, nword, floats, words
  use gulp_lengths
  use species,       only : natspec, ntypspec, nspec
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),       intent(in)  :: maxqrt
  character(len=5),  intent(out) :: syms(4,maxqrt)
  integer(i4),       intent(in)  :: iline
  integer(i4),       intent(out) :: itypes(4,maxqrt)
  integer(i4),       intent(out) :: nqrt
  integer(i4),       intent(out) :: nbeg
  integer(i4),       intent(out) :: nvars(4,maxqrt)
  logical,           intent(in)  :: llibrary
  logical,           intent(out) :: lvalidpot
!
!  Local variables
!
  character(len=5)               :: lab1
  character(len=5)               :: lab2
  character(len=5)               :: lab3
  character(len=5)               :: lab4
  character(len=maxwordlength)   :: word
  character(len=maxwordlength)   :: wordsave
  character(len=maxwordlength)   :: wordsave1
  character(len=maxwordlength)   :: wordsave2
  character(len=maxwordlength)   :: wordsave3
  character(len=maxwordlength)   :: wordsave4
  integer(i4)                    :: i
  integer(i4),              save :: maxilp = 20
  integer(i4)                    :: maxn4
  integer(i4)                    :: n1
  integer(i4)                    :: n2
  integer(i4)                    :: n3
  integer(i4)                    :: n4
  integer(i4)                    :: nilp(4)
  integer(i4), allocatable, save :: ilp(:,:)
  integer(i4), allocatable, save :: ityp(:,:)
  integer(i4), allocatable, save :: nvr(:,:)
  integer(i4)                    :: nptr
  logical                        :: lok1
  logical                        :: lsame14
  logical                        :: lsame23
  logical,     allocatable, save :: ltyp(:,:)
#ifdef TRACE
  call trace_in('getpotsymbols4')
#endif
!
!  This routine is currently only designed for use with libraries and so check that this is the case
!
  if (.not.llibrary) then
    call outerror('getpotsymbols4 should only be called when reading a library',iline)
    call stopnow('getpotsymbols4')
  endif
!
  lvalidpot = .true.
  nbeg = 0
  nqrt = 0
!
!  Allocate local workspace
!
  allocate(ilp(maxilp,4))
  allocate(ityp(maxilp,4))
  allocate(ltyp(maxilp,4))
  allocate(nvr(maxilp,4))
!
  if (nword.gt.0) then
    if (nword.ge.4) then
      nptr = 0
      do i = 1,nword
        wordsave = words(i)
        word = words(i)
        call stolc(word,maxwordlength)
        if (index(word,'she').eq.1.and.nptr.gt.0) then
          nvr(1:nilp(nptr),nptr) = nvr(1:nilp(nptr),nptr) + maxele
        elseif (index(word,'cor').eq.0) then
          nptr = nptr + 1
          if (nptr.eq.1) then
            wordsave1 = words(i)
          elseif (nptr.eq.2) then
            wordsave2 = words(i)
          elseif (nptr.eq.3) then
            wordsave3 = words(i)
            lsame23 = (wordsave2.eq.wordsave3)
          elseif (nptr.eq.4) then
            wordsave4 = words(i)
            lsame14 = (wordsave1.eq.wordsave4)
          endif
          call okspecs(lok1,words(i),maxilp,nilp(nptr),ilp(1,nptr),.true.,ltyp(1,nptr))
!
!
!  Check dimensions of species arrays
!
          if (nilp(nptr).gt.maxilp) then
            deallocate(nvr)
            deallocate(ltyp)
            deallocate(ityp)
            deallocate(ilp)
!
            maxilp = nilp(nptr)
!
            allocate(ilp(maxilp,4))
            allocate(ityp(maxilp,4))
            allocate(ltyp(maxilp,4))
            allocate(nvr(maxilp,4))
!
!  Repeat species calls with correct array sizes
!
            call okspecs(lok1,words(i),maxilp,nilp(nptr),ilp(1,nptr),.true.,ltyp(1,nptr))
          endif
!
          if (.not.lok1) then
            lvalidpot = .false.
          endif
!
          if (nilp(nptr).gt.0) then
            do n1 = 1,nilp(nptr)
              if (ilp(n1,nptr).gt.0.and.ilp(n1,nptr).le.nspec) then
                nvr(n1,nptr) = natspec(ilp(n1,nptr))
                if (ltyp(n1,nptr)) then
                  ityp(n1,nptr) = 0
                else
                  ityp(n1,nptr) = ntypspec(ilp(n1,nptr))
                endif
              else 
                nvr(n1,nptr) = maxele
                ityp(n1,nptr) = 0
              endif
            enddo
          endif
        endif
      enddo
      if (nptr.ne.4) then
        call outerror('Incorrect potential species input for four-body potential',iline)
        call stopnow('getpotsymbols4')
      endif
    else
      call outerror('Incorrect potential species input for four-body potential',iline)
      call stopnow('getpotsymbols4')
    endif
    nbeg = 0
!
!  Generate quartets
!
    do n2 = 1,nilp(2)
      call label(nvr(n2,2),ityp(n2,2),lab2)
      if (lsame23) then
        do n3 = 1,n2
          call label(nvr(n3,3),ityp(n3,3),lab3)
          do n1 = 1,nilp(1)
            call label(nvr(n1,1),ityp(n1,1),lab1)
            if (n3.eq.n2.and.lsame14) then
              maxn4 = n1
            else
              maxn4 = nilp(4)
            endif
            do n4 = 1,maxn4
              nqrt = nqrt + 1
              if (nqrt.le.maxqrt) then
                call label(nvr(n4,4),ityp(n4,4),lab4)
                syms(1,nqrt) = lab1
                syms(2,nqrt) = lab2
                syms(3,nqrt) = lab3
                syms(4,nqrt) = lab4
                nvars(1,nqrt) = nvr(n1,1)
                nvars(2,nqrt) = nvr(n2,2)
                nvars(3,nqrt) = nvr(n3,3)
                nvars(4,nqrt) = nvr(n4,4)
                itypes(1,nqrt) = ityp(n1,1)
                itypes(2,nqrt) = ityp(n2,2)
                itypes(3,nqrt) = ityp(n3,3)
                itypes(4,nqrt) = ityp(n4,4)
              endif
            enddo
          enddo
        enddo
      else
        do n3 = 1,nilp(3)
          call label(nvr(n3,3),ityp(n3,3),lab3)
          do n1 = 1,nilp(1)
            call label(nvr(n1,1),ityp(n1,1),lab1)
            do n4 = 1,nilp(4)
              nqrt = nqrt + 1
              if (nqrt.le.maxqrt) then
                call label(nvr(n4,4),ityp(n4,4),lab4)
                syms(1,nqrt) = lab1
                syms(2,nqrt) = lab2
                syms(3,nqrt) = lab3
                syms(4,nqrt) = lab4
                nvars(1,nqrt) = nvr(n1,1)
                nvars(2,nqrt) = nvr(n2,2)
                nvars(3,nqrt) = nvr(n3,3)
                nvars(4,nqrt) = nvr(n4,4)
                itypes(1,nqrt) = ityp(n1,1)
                itypes(2,nqrt) = ityp(n2,2)
                itypes(3,nqrt) = ityp(n3,3)
                itypes(4,nqrt) = ityp(n4,4)
              endif
            enddo
          enddo
        enddo
      endif
    enddo
  else
!
!  Numeric input
!
    nqrt = 1
    nvars(1,1) = int(floats(1))
    nvars(2,1) = int(floats(2))
    nvars(3,1) = int(floats(3))
    nvars(4,1) = int(floats(4))
    if (nvars(1,1).gt.100) nvars(1,1) = nvars(1,1) - 100 + maxele
    if (nvars(2,1).gt.100) nvars(2,1) = nvars(2,1) - 100 + maxele
    if (nvars(3,1).gt.100) nvars(3,1) = nvars(3,1) - 100 + maxele
    if (nvars(4,1).gt.100) nvars(4,1) = nvars(4,1) - 100 + maxele
    itypes(1,1) = 0
    itypes(2,1) = 0
    itypes(3,1) = 0
    itypes(4,1) = 0
    nbeg = 4
    nfloat = nfloat - 4
    call label(nvars(1,1),itypes(1,1),syms(1,1))
    call label(nvars(2,1),itypes(2,1),syms(2,1))
    call label(nvars(3,1),itypes(3,1),syms(3,1))
    call label(nvars(4,1),itypes(4,1),syms(4,1))
  endif
!
!  Free local workspace
!
  deallocate(nvr)
  deallocate(ltyp)
  deallocate(ityp)
  deallocate(ilp)
#ifdef TRACE
  call trace_out('getpotsymbols4')
#endif
!
  return
  end
!
  subroutine getpotsymbol6(iline,llibrary,nvar1,itype1,sym1,nvar2,itype2,sym2,nvar3,itype3,sym3, &
                           nvar4,itype4,sym4,nvar5,itype5,sym5,nvar6,itype6,sym6,nbeg,lvalidpot)
!
!  Gets species symbols for four-body potential input
!
!  11/04 Created from getpotsymbol4
!  10/05 Correction to handling of atom types added
!   8/06 ltype0 flag added to okspec call
!   9/06 Literal symbols now returned
!  11/06 Words are saved for literal symbols before corruption by okspec call
!  12/08 Module input renamed to gulpinput
!   3/10 ltype0 case corrected to set type number to zero
!   2/18 Trace added
!   1/19 maxwordlength changes added
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, January 2019
!
  use datatypes  
  use element,       only : maxele
  use gulpinput,     only : nfloat, nword, floats, words
  use gulp_lengths
  use species,       only : natspec, ntypspec, nspec
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  character(len=5), intent(out) :: sym1
  character(len=5), intent(out) :: sym2
  character(len=5), intent(out) :: sym3
  character(len=5), intent(out) :: sym4
  character(len=5), intent(out) :: sym5
  character(len=5), intent(out) :: sym6
  integer(i4),      intent(in)  :: iline
  integer(i4),      intent(out) :: itype1
  integer(i4),      intent(out) :: itype2
  integer(i4),      intent(out) :: itype3
  integer(i4),      intent(out) :: itype4
  integer(i4),      intent(out) :: itype5
  integer(i4),      intent(out) :: itype6
  integer(i4),      intent(out) :: nbeg
  integer(i4),      intent(out) :: nvar1
  integer(i4),      intent(out) :: nvar2
  integer(i4),      intent(out) :: nvar3
  integer(i4),      intent(out) :: nvar4
  integer(i4),      intent(out) :: nvar5
  integer(i4),      intent(out) :: nvar6
  logical,          intent(in)  :: llibrary
  logical,          intent(out) :: lvalidpot
!
!  Local variables
!
  character(len=maxwordlength)  :: word
  character(len=maxwordlength)  :: wordsave
  integer(i4)                   :: i
  integer(i4)                   :: ilp1
  integer(i4)                   :: ityp(6)
  integer(i4)                   :: nptr
  integer(i4)                   :: nvr(6)
  logical                       :: lok1
  logical                       :: ltype0
#ifdef TRACE
  call trace_in('getpotsymbol6')
#endif
!
  lvalidpot = .true.
  nbeg = 0
!
  if (nword.gt.0) then
    if (nword.ge.6) then
      nptr = 0
      do i = 1,nword
        wordsave = words(i)
        word = words(i)
        call stolc(word,maxwordlength)
        if (index(word,'she').eq.1.and.nptr.gt.0) then
          nvr(nptr) = nvr(nptr) + maxele
        elseif (index(word,'cor').eq.0) then
          nptr = nptr + 1
          if (llibrary) then
            call okspec(lok1,words(i),ilp1,.true.,ltype0)
            if (.not.lok1) then
              lvalidpot = .false.
            endif
            if (ilp1.gt.0.and.ilp1.le.nspec) then
              nvr(nptr) = natspec(ilp1)
              if (ltype0) then
                ityp(nptr) = 0
              else
                ityp(nptr) = ntypspec(ilp1)
              endif
            else 
              nvr(nptr) = maxele
              ityp(nptr) = 0
            endif
          else
            call ltont(words(i),nvr(nptr),ityp(nptr))
          endif
          if (nptr.eq.1) then
            sym1 = wordsave(1:5)
          elseif (nptr.eq.2) then
            sym2 = wordsave(1:5)
          elseif (nptr.eq.3) then
            sym3 = wordsave(1:5)
          elseif (nptr.eq.4) then
            sym4 = wordsave(1:5)
          elseif (nptr.eq.5) then
            sym5 = wordsave(1:5)
          elseif (nptr.eq.6) then
            sym6 = wordsave(1:5)
          endif
        endif
      enddo
      if (nptr.ne.6) then
        call outerror('Incorrect potential species input for six-body potential',iline)
        call stopnow('getpotsymbol6')
      endif
    else
      call outerror('Incorrect potential species input for six-body potential',iline)
      call stopnow('getpotsymbol6')
    endif
    nbeg = 0
    nvar1 = nvr(1)
    nvar2 = nvr(2)
    nvar3 = nvr(3)
    nvar4 = nvr(4)
    nvar5 = nvr(5)
    nvar6 = nvr(6)
    itype1 = ityp(1)
    itype2 = ityp(2)
    itype3 = ityp(3)
    itype4 = ityp(4)
    itype5 = ityp(5)
    itype6 = ityp(6)
  else
!
!  Numeric input
!
    nvar1 = int(floats(1))
    nvar2 = int(floats(2))
    nvar3 = int(floats(3))
    nvar4 = int(floats(4))
    nvar5 = int(floats(5))
    nvar6 = int(floats(6))
    if (nvar1.gt.100) nvar1 = nvar1 - 100 + maxele
    if (nvar2.gt.100) nvar2 = nvar2 - 100 + maxele
    if (nvar3.gt.100) nvar3 = nvar3 - 100 + maxele
    if (nvar4.gt.100) nvar4 = nvar4 - 100 + maxele
    if (nvar5.gt.100) nvar5 = nvar5 - 100 + maxele
    if (nvar6.gt.100) nvar6 = nvar6 - 100 + maxele
    itype1 = 0
    itype2 = 0
    itype3 = 0
    itype4 = 0
    itype5 = 0
    itype6 = 0
    nbeg = 6
    nfloat = nfloat - 6
    call label(nvar1,itype1,sym1)
    call label(nvar2,itype2,sym2)
    call label(nvar3,itype3,sym3)
    call label(nvar4,itype4,sym4)
    call label(nvar5,itype5,sym5)
    call label(nvar6,itype6,sym6)
  endif
#ifdef TRACE
  call trace_out('getpotsymbol6')
#endif
!
  return
  end
