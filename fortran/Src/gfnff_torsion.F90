  subroutine gfnff_torsion(n4ty,rkfor,rn,phi0,r12,r13,r14,r23,r24,r34,efor,e1d,e2d,lgrad1,lgrad2)
!
!  Calculates four-body energy, first, and second derivatives with respect to the 
!  six interatomic distances that make the four-body term.
!
!  GFNFF version that supports the following cases:
!
!  rn = 1, 2, 3, 6
!  rn*phi0 = multiple of pi
!
!   7/21 New version created that handles the special cases for GFNFF to avoid 1/sin(phi) problem
!   7/21 n = 3/6 cases corrected
!
!  nptr4   = pointer to potential number
!  n4ty    = pointer to type of four-body potential
!  r12     = distance between atoms 1 and 2
!  r13     = distance between atoms 1 and 3
!  r14     = distance between atoms 1 and 4
!  r23     = distance between atoms 2 and 3
!  r24     = distance between atoms 2 and 4
!  r34     = distance between atoms 3 and 4
!  efor    = contribution to four-body energy
!  e1d     = array of first derivative terms
!  e2d     = array of second derivative terms
!  rkfor   = first parameter associated with potential type
!  phi0    = second parameter associated with potential type
!  rn      = third parameter associated with potential type
!  lgrad1  = if .true. calculate the first derivatives
!  lgrad2  = if .true. calculate the second derivatives
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
!  Julian Gale, CIC, Curtin University, July 2021
!
  use g_constants
  use iochannels
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: n4ty
  logical,     intent(in)    :: lgrad1
  logical,     intent(in)    :: lgrad2
  real(dp),    intent(out)   :: e1d(6)
  real(dp),    intent(out)   :: e2d(21)
  real(dp),    intent(out)   :: efor
  real(dp),    intent(in)    :: phi0  
  real(dp),    intent(in)    :: r12
  real(dp),    intent(in)    :: r13
  real(dp),    intent(in)    :: r14
  real(dp),    intent(in)    :: r23
  real(dp),    intent(in)    :: r24
  real(dp),    intent(in)    :: r34
  real(dp),    intent(in)    :: rkfor
  real(dp),    intent(in)    :: rn
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: ii
  integer(i4)                :: j
  integer(i4)                :: k
  logical                    :: lgrad1loc
  logical                    :: lgrad2loc
  real(dp)                   :: cos1
  real(dp)                   :: cos2
  real(dp)                   :: cos3
  real(dp)                   :: cosp2
  real(dp)                   :: cosp4
  real(dp)                   :: cosp6
  real(dp)                   :: cosphi
  real(dp)                   :: cos11d(6),cos21d(6),cos31d(6),cosp1d(6)
  real(dp)                   :: cos12d(21),cos22d(21),cos32d(21),cosp2d(21)
  real(dp)                   :: cos13d(56),cos33d(56),cosp3d(56)
  real(dp)                   :: e1
  real(dp)                   :: e2
  real(dp)                   :: phase_shift
  real(dp)                   :: r122
  real(dp)                   :: r132
  real(dp)                   :: r142
  real(dp)                   :: r232
  real(dp)                   :: r242
  real(dp)                   :: r342
  real(dp)                   :: rr12
  real(dp)                   :: rr122
  real(dp)                   :: rr124
  real(dp)                   :: rr13
  real(dp)                   :: rr132
  real(dp)                   :: rr134
  real(dp)                   :: rr14
  real(dp)                   :: rr142
  real(dp)                   :: rr144
  real(dp)                   :: rr23
  real(dp)                   :: rr232
  real(dp)                   :: rr234
  real(dp)                   :: rr24
  real(dp)                   :: rr242
  real(dp)                   :: rr244
  real(dp)                   :: rr34
  real(dp)                   :: rr342
  real(dp)                   :: rr344
  real(dp)                   :: rsin1
  real(dp)                   :: rsin3
  real(dp)                   :: rtan1
  real(dp)                   :: rtan3
  real(dp)                   :: sin11d(6),sin21d(6),sin31d(6)
  real(dp)                   :: sin12d(21),sin22d(21),sin32d(21)
  real(dp)                   :: sin13d(56),sin33d(56)
  real(dp)                   :: sin1
  real(dp)                   :: sin3
#ifdef TRACE
  call trace_in('gfnff_torsion')
#endif
!
!  Zero terms
!
  efor = 0.0_dp
!
!  Find sign based on phase shift
!
  phase_shift = cos(pi - rn*phi0)
  if (abs(phase_shift).lt.0.9999_dp) then
    call outerror('gfnff_torsion called with n*phi0 not being a multiple of pi',0_i4)
    call stopnow('gfnff_torsion')
  endif
!
  if (lgrad1) then
    do i = 1,6
      e1d(i) = 0.0_dp
      cos11d(i) = 0.0_dp
      cos21d(i) = 0.0_dp
      cos31d(i) = 0.0_dp
      cosp1d(i) = 0.0_dp
      sin11d(i) = 0.0_dp
      sin21d(i) = 0.0_dp
      sin31d(i) = 0.0_dp
    enddo
    if (lgrad2) then
      do i = 1,21
        e2d(i) = 0.0_dp
        cos12d(i) = 0.0_dp
        cos22d(i) = 0.0_dp
        cos32d(i) = 0.0_dp
        cosp2d(i) = 0.0_dp
        sin12d(i) = 0.0_dp
        sin22d(i) = 0.0_dp
        sin32d(i) = 0.0_dp
      enddo
    endif
  endif
!
!  Set up local constants
!
  r122 = r12*r12
  r132 = r13*r13
  r142 = r14*r14
  r232 = r23*r23
  r242 = r24*r24
  r342 = r34*r34
  rr12 = 1.0_dp/r12
  rr13 = 1.0_dp/r13
  rr14 = 1.0_dp/r14
  rr23 = 1.0_dp/r23
  rr24 = 1.0_dp/r24
  rr34 = 1.0_dp/r34
  rr122 = rr12*rr12
  rr132 = rr13*rr13
  rr142 = rr14*rr14
  rr232 = rr23*rr23
  rr242 = rr24*rr24
  rr342 = rr34*rr34
  if (lgrad2) then
    rr124 = rr122*rr122
    rr134 = rr132*rr132
    rr144 = rr142*rr142
    rr234 = rr232*rr232
    rr244 = rr242*rr242
    rr344 = rr342*rr342
  endif
!
!  Set local derivative flags
!
  lgrad1loc = lgrad1
  lgrad2loc = lgrad2
!***************************************
!  Set up potential independent terms  *
!***************************************
!$$$$$$$$$$$$$$$$$
!  Cosine terms  $
!$$$$$$$$$$$$$$$$$
!
!  Cosine theta 1 = 1-2-3, 2 = 2-3-4 and 3 = 3-2-4
!
  cos1 = 0.5_dp*(r232+r122-r132)/(r12*r23)
  cos2 = 0.5_dp*(r232+r342-r242)/(r23*r34)
  cos3 = 0.5_dp*(r232+r242-r342)/(r23*r24)
!
!  Check for angles which are 0 or 180 degrees with potentials
!  which cannot cope with these. For those that can the four-
!  body contribution must go to zero when the angle approaches
!  this limit - hence we can just return having set all the
!  derivatives and energy to zero
!
  if (abs(cos1).ge.0.99999999_dp) then
    lgrad1loc = .false.
    lgrad2loc = .false.
!
! Energy can be unstable when torsion is undefined so return
!
    efor = 0.0_dp
#ifdef TRACE
    call trace_out('gfnff_torsion')
#endif
    return
  endif
  if (abs(cos2).ge.0.99999999_dp) then
    lgrad1loc = .false.
    lgrad2loc = .false.
!
! Energy can be unstable when torsion is undefined so return
!
    efor = 0.0_dp
#ifdef TRACE
    call trace_out('gfnff_torsion')
#endif
    return
  endif
  if (lgrad1loc) then
!
!  First
!
!  1 = r12
!  2 = r13
!  3 = r14
!  4 = r23
!  5 = r24
!  6 = r34
!
    cos11d(1) = rr12*rr23 - cos1*rr122
    cos11d(2) = - rr12*rr23
    cos11d(4) = rr12*rr23 - cos1*rr232
!
    cos21d(4) = rr23*rr34 - cos2*rr232
    cos21d(5) = - rr23*rr34
    cos21d(6) = rr23*rr34 - cos2*rr342
!
    cos31d(4) = rr23*rr24 - cos3*rr232
    cos31d(5) = rr23*rr24 - cos3*rr242
    cos31d(6) = - rr23*rr24
!
    if (lgrad2loc) then
!
!  Second
!
!  1 = 11  7 = 22 12 = 33 16 = 44 19 = 55 21 = 66
!  2 = 21  8 = 32 13 = 43 17 = 54 20 = 65
!  3 = 31  9 = 42 14 = 53 18 = 64
!  4 = 41 10 = 52 15 = 63
!  5 = 51 11 = 62
!  6 = 61 
!
      cos12d(1) = - 2.0_dp*rr122*rr12*rr23 + 3.0_dp*cos1*rr124
      cos12d(2) = rr122*rr12*rr23
      cos12d(4) = rr12*rr23*(cos1*rr12*rr23 - rr122 - rr232)
      cos12d(9) = rr232*rr23*rr12
      cos12d(16) = - 2.0_dp*rr232*rr23*rr12 + 3.0_dp*cos1*rr234
!
      cos22d(16) = - 2.0_dp*rr232*rr23*rr34 + 3.0_dp*cos2*rr234
      cos22d(17) = rr232*rr23*rr34
      cos22d(18) = rr23*rr34*(cos2*rr23*rr34 - rr232 - rr342)
      cos22d(20) = rr342*rr34*rr23
      cos22d(21) = - 2.0_dp*rr342*rr34*rr23 + 3.0_dp*cos2*rr344
!
      cos32d(16) = - 2.0_dp*rr232*rr23*rr24 + 3.0_dp*cos3*rr234
      cos32d(17) = rr23*rr24*(cos3*rr23*rr24 - rr232 - rr242)
      cos32d(18) = rr232*rr23*rr24
      cos32d(19) = - 2.0_dp*rr242*rr24*rr23 + 3.0_dp*cos3*rr244
      cos32d(20) = rr242*rr24*rr23
    endif
  endif
!$$$$$$$$$$$$$$$
!  Sine terms  $
!$$$$$$$$$$$$$$$
  sin1 = sqrt(1.0_dp - cos1*cos1)
  sin3 = sqrt(1.0_dp - cos3*cos3)
  rsin1 = 1.0_dp/sin1
  rsin3 = 1.0_dp/sin3
!
  if (lgrad1loc) then
!
!  First derivatives
!
    rtan1 = cos1*rsin1
    rtan3 = cos3*rsin3
    do i = 1,6
      sin11d(i) = - rtan1*cos11d(i)
      sin31d(i) = - rtan3*cos31d(i)
    enddo
    if (lgrad2loc) then
!
!  Second derivatives
!
      ii = 0
      do i = 1,6
        do j = i,6
          ii = ii + 1
          sin12d(ii) = - rtan1*cos12d(ii) - rsin1*cos11d(i)*cos11d(j)*(1.0_dp + rtan1*rtan1)
          sin32d(ii) = - rtan3*cos32d(ii) - rsin3*cos31d(i)*cos31d(j)*(1.0_dp + rtan3*rtan3)
        enddo
      enddo
    endif
  endif
!$$$$$$$$$$$$$$
!  Phi terms  $
!$$$$$$$$$$$$$$
  call getphi(r12,r14,r24,1_i4,3_i4,5_i4,cos1,cos3,sin1,sin3,cos11d,cos31d,sin11d,sin31d, &
              cos12d,cos32d,sin12d,sin32d,cos13d,cos33d,sin13d,sin33d,cosphi,cosp1d, &
              cosp2d,cosp3d,lgrad1loc,lgrad2loc,.false.)
!******************************
!  Potential dependent terms  *
!******************************
  if (n4ty.eq.1) then
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Standard torsional potential  $
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if (rn.eq.1) then
!---------
! n = 1  |
!---------
      efor = rkfor*(1.0_dp + phase_shift*cosphi)
      if (lgrad1loc) then
        e1 = rkfor*phase_shift
        do k = 1,6
          e1d(k) = e1*cosp1d(k)
        enddo
        if (lgrad2loc) then
          do k = 1,21
            e2d(k) = e1*cosp2d(k)
          enddo
        endif
      endif
    elseif (rn.eq.2) then
!---------
! n = 2  |
!---------
      efor = rkfor*(1.0_dp + phase_shift*(2.0_dp*cosphi**2 - 1.0_dp))
      if (lgrad1loc) then
        e1 = 4.0_dp*rkfor*phase_shift*cosphi
        do k = 1,6
          e1d(k) = e1*cosp1d(k)
        enddo
        if (lgrad2loc) then
          e2 = 4.0_dp*rkfor*phase_shift
          ii = 0
          do i = 1,6
            do j = i,6
              ii = ii + 1
              e2d(ii) = e1*cosp2d(ii) + e2*cosp1d(i)*cosp1d(j)
            enddo
          enddo
        endif
      endif
    elseif (rn.eq.3) then
!---------
! n = 3  |
!---------
      cosp2 = cosphi**2
      efor = rkfor*(1.0_dp + phase_shift*(4.0_dp*cosp2 - 3.0_dp)*cosphi)
      if (lgrad1loc) then
        e1 = rkfor*phase_shift*(12.0_dp*cosp2 - 3.0_dp)
        do k = 1,6
          e1d(k) = e1*cosp1d(k)
        enddo
        if (lgrad2loc) then
          e2 = 24.0_dp*rkfor*phase_shift*cosphi
          ii = 0
          do i = 1,6
            do j = i,6
              ii = ii + 1
              e2d(ii) = e1*cosp2d(ii) + e2*cosp1d(i)*cosp1d(j)
            enddo
          enddo
        endif
      endif
    elseif (rn.eq.6) then
!---------
! n = 6  |
!---------
      cosp2 = cosphi**2
      cosp4 = cosp2**2
      cosp6 = cosp4*cosp2
      efor = rkfor*(1.0_dp + phase_shift*(32.0_dp*cosp6 - 48.0_dp*cosp4 + 18.0_dp*cosp2 - 1.0_dp))
      if (lgrad1loc) then
        e1 = rkfor*phase_shift*(192.0_dp*cosp4 - 192.0_dp*cosp2 + 36.0_dp)*cosphi
        do k = 1,6
          e1d(k) = e1*cosp1d(k)
        enddo
        if (lgrad2loc) then
          e2 = rkfor*phase_shift*(960.0_dp*cosp4 - 576.0_dp*cosp2 + 36.0_dp)
          ii = 0
          do i = 1,6
            do j = i,6
              ii = ii + 1
              e2d(ii) = e1*cosp2d(ii) + e2*cosp1d(i)*cosp1d(j)
            enddo
          enddo
        endif
      endif
    else
      call outerror('gfnff_torsion called with invalid value of n',0_i4)
      call stopnow('gfnff_torsion')
    endif
  endif
#ifdef TRACE
  call trace_out('gfnff_torsion')
#endif
!
  return
  end
