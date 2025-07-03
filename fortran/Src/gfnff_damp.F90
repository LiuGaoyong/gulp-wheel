  subroutine gfnff_damp(r2,r2cut,damp,ddampdr,d2dampdr2,lgrad1,lgrad2)
!
!  Damping function for GFNFF
!
  use datatypes
  implicit none
!
!  Passed variables
!
  logical,     intent(in)  :: lgrad1
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(in)  :: r2
  real(dp),    intent(in)  :: r2cut     ! Cutoff factor
  real(dp),    intent(out) :: damp      ! Damping factor
  real(dp),    intent(out) :: ddampdr   ! First derivative of damping factor (1/r(d(damp)/dr)
  real(dp),    intent(out) :: d2dampdr2 ! Second derivative of damping factor (1/r d/dr(1/r(d(damp)/dr))
!
!  Local variables
!
  real(dp)                 :: dr2ratiodr
  real(dp)                 :: d2r2ratiodr2
  real(dp)                 :: r2ratio
!
  r2ratio = r2*r2/r2cut
  damp = 1.0_dp/(1.0_dp + r2ratio)
  if (lgrad1) then
    dr2ratiodr = 4.0_dp*r2ratio/r2
    ddampdr = - dr2ratiodr*damp**2
    if (lgrad2) then
      d2r2ratiodr2 = 8.0_dp/r2cut
      d2dampdr2 = damp*damp*(2.0_dp*dr2ratiodr*dr2ratiodr*damp - d2r2ratiodr2)
    endif
  endif

  end subroutine gfnff_damp
