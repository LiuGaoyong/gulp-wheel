  subroutine einsteintrial(eeinstein,ntrialatom,nptrtrialatom)
!
!  Subroutine for calculating the energy and derivatives due to the Einstein model
!  due to a trial subset of atoms
!
!   1/08 Created 
!   2/18 Trace added
!   3/21 Generalised to all dimensionalities
!  12/21 Einstein positions changed
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
!  Julian Gale, CIC, Curtin University, December 2021
!
  use current
  use configurations
  use parallel
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: ntrialatom
  integer(i4), intent(in)    :: nptrtrialatom(ntrialatom)
  real(dp), intent(out)      :: eeinstein
!
!  Local variables
!
  integer(i4)                :: i
  integer(i4)                :: nt
  real(dp)                   :: ki
  real(dp)                   :: xf
  real(dp)                   :: yf
  real(dp)                   :: zf
  real(dp)                   :: xi
  real(dp)                   :: yi
  real(dp)                   :: zi
!
!  Initialise energy
!
  eeinstein = 0.0_dp
!
!  If there are no Einstein model atoms return
!
  if (.not.leinstein) return
#ifdef TRACE
  call trace_in('einsteintrial')
#endif
!
!  Loop over asymmetric unit performing integral
!
  do nt = 1+procid,ntrialatom,nprocs
    i = nptrtrialatom(nt)
    if (leinsteinat(nsft+i)) then
!
!  Find nearest image in fraction coordinates - must use image of symmetry
!  related atom to avoid problems where lattice points are only defined for
!  the asymmetric unit.
!
      if (ndim.eq.3) then
        xf = xafrac(i) - xfeinsteinat(nsft + i)
        yf = yafrac(i) - yfeinsteinat(nsft + i)
        zf = zafrac(i) - zfeinsteinat(nsft + i)
        xf = mod(xf+10.0_dp,1.0_dp)
        yf = mod(yf+10.0_dp,1.0_dp)
        zf = mod(zf+10.0_dp,1.0_dp)
        if (xf.gt.0.5_dp) xf = xf - 1.0_dp
        if (yf.gt.0.5_dp) yf = yf - 1.0_dp
        if (zf.gt.0.5_dp) zf = zf - 1.0_dp
!
!  Convert coordinate difference into Cartesian frame
!
        xi = xf*rv(1,1) + yf*rv(1,2) + zf*rv(1,3)
        yi = xf*rv(2,1) + yf*rv(2,2) + zf*rv(2,3)
        zi = xf*rv(3,1) + yf*rv(3,2) + zf*rv(3,3)
      elseif (ndim.eq.2) then
        xf = xafrac(i) - xfeinsteinat(nsft + i)
        yf = yafrac(i) - yfeinsteinat(nsft + i)
        zi = zafrac(i) - zfeinsteinat(nsft + i)
        xf = mod(xf+10.0_dp,1.0_dp)
        yf = mod(yf+10.0_dp,1.0_dp)
        if (xf.gt.0.5_dp) xf = xf - 1.0_dp
        if (yf.gt.0.5_dp) yf = yf - 1.0_dp
!
!  Convert coordinate difference into Cartesian frame
!
        xi = xf*rv(1,1) + yf*rv(1,2)
        yi = xf*rv(2,1) + yf*rv(2,2)
      elseif (ndim.eq.1) then
        xf = xafrac(i) - xfeinsteinat(nsft + i)
        yi = yafrac(i) - yfeinsteinat(nsft + i)
        zi = zafrac(i) - zfeinsteinat(nsft + i)
        xf = mod(xf+10.0_dp,1.0_dp)
        if (xf.gt.0.5_dp) xf = xf - 1.0_dp
!
!  Convert coordinate difference into Cartesian frame
!
        xi = xf*rv(1,1)
      else
        xi = xafrac(i) - xfeinsteinat(nsft + i)
        yi = yafrac(i) - yfeinsteinat(nsft + i)
        zi = zafrac(i) - zfeinsteinat(nsft + i)
      endif
!
!  Evaluate energy and derivatives
!
      ki = keinsteinat(nsft + i)
      eeinstein = eeinstein + 0.5_dp*ki*(xi*xi + yi*yi + zi*zi)
    endif
  enddo
#ifdef TRACE
  call trace_out('einsteintrial')
#endif
!
  return
  end
