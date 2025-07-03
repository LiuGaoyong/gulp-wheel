!***********************
!  Gamma - GFNFF form  *
!***********************
  subroutine gfnff_gamma(zab,r,gam,dgam,d2gam,lgrad1,lgrad2)
!
!  Subroutine calculates the GFNFF Coulomb term erf(zeta*r)/r
! 
!   9/20 Created
!   2/21 Correction to second derivatives
!   3/21 Tapering added
!
!  On input :
!
!  zab    =  gamma exponent for A-B pair
!  r      =  distance between A and B
!  lgrad1 = if true then compute first derivatives
!  lgrad2 = if true then compute second derivatives
!
!  On exit :
!
!  gam        =  erf(zab*r)/r 
!  dgam       =  if lgrad1 then this is the first derivative w.r.t. r
!  d2gam      =  if lgrad2 then this is the second derivative w.r.t. r
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
!  Julian Gale, CIC, Curtin University, March 2021
!
  use datatypes
  use element,       only : rqeq, rqeqtaper
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp), intent(in)     :: r
  real(dp), intent(in)     :: zab
  real(dp), intent(out)    :: dgam
  real(dp), intent(out)    :: d2gam
  real(dp), intent(out)    :: gam
  logical,  intent(in)     :: lgrad1
  logical,  intent(in)     :: lgrad2
!
!  Local variables
!
  real(dp)                 :: etrm
  real(dp)                 :: expzr
  real(dp)                 :: g_derf
  real(dp)                 :: gamrr
  real(dp)                 :: dgamrr
  real(dp)                 :: d2gamrr
  real(dp)                 :: rr
  real(dp)                 :: tpfn
  real(dp)                 :: dtpfn
  real(dp)                 :: d2tpfn
  real(dp)                 :: d3tpfn
  real(dp)                 :: tworootpi
#ifdef TRACE
  call trace_in('gfnff_gamma')
#endif
!
  gam = 0.0_dp
  dgam = 0.0_dp
  d2gam = 0.0_dp
!
!  Exclude self term
!
  if (r.lt.1.0d-10) then
#ifdef TRACE
    call trace_out('gfnff_gamma')
#endif
    return
  endif
!
  rr = 1.0_dp/r
  etrm = g_derf(zab*r)
  gam = etrm*rr
  if (lgrad1) then
    tworootpi = 1.0_dp/sqrt(atan(1.0_dp))
    expzr = tworootpi*exp(-(zab*r)**2)
    dgam = expzr*zab*rr - gam*rr
    if (lgrad2) then
      d2gam = - 2.0_dp*zab*expzr*(rr*rr + zab*zab) + 2.0_dp*gam*rr*rr
    endif
  endif
!
!  Taper function
!
  if (r.gt.rqeqtaper) then
    call mdftaper(r,rqeqtaper,rqeq,tpfn,dtpfn,d2tpfn,d3tpfn,lgrad1,lgrad2,.false.)
!
!  Compute pure Coulomb values
!
    gamrr = rr
    dgamrr = - rr**2
    d2gamrr = 2.0_dp*rr**3
!
    d2gam = d2gam*tpfn + d2gamrr*(1.0_dp - tpfn) + 2.0_dp*(dgam - dgamrr)*dtpfn + (gam - gamrr)*d2tpfn
    dgam = dgam*tpfn + dgamrr*(1.0_dp - tpfn) + (gam - gamrr)*dtpfn
    gam = gam*tpfn + gamrr*(1.0_dp - tpfn)
  endif
#ifdef TRACE
  call trace_out('gfnff_gamma')
#endif
!
  return
  end subroutine gfnff_gamma
