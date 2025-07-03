  subroutine gfnff_radnoq(nati,natj,cni,cnj,radij,dradijdcni,dradijdcnj,lgrad1)
!
!  Computes the coordination number dependent radii for a pair of atoms i and j for GFNFF
!
!  On entry : 
!
!  nati            = atomic number for i
!  natj            = atomic number for j
!  cni             = coordination number for i
!  cnj             = coordination number for j
!  radij           = shift in radii for i and j (in Ang)
!  lgrad1          = if true then compute first derivatives
!
!  On exit :
!
!  radij           = sum of radii for i and j (in Ang)
!  dradijdcni      = derivative of radij w.r.t. cni
!  dradijdcnj      = derivative of radij w.r.t. cnj
!
!   8/20 Created
!  10/20 First derivatives added
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use datatypes
  use gulp_gfnff
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: nati
  integer(i4), intent(in)                          :: natj
  logical,     intent(in)                          :: lgrad1
  real(dp),    intent(in)                          :: cni
  real(dp),    intent(in)                          :: cnj
  real(dp),    intent(inout)                       :: radij
  real(dp),    intent(out)                         :: dradijdcni
  real(dp),    intent(out)                         :: dradijdcnj
!
!  Local variables
!
  real(dp)                                         :: den
  real(dp)                                         :: f12
  real(dp)                                         :: f1
  real(dp)                                         :: f2
  real(dp)                                         :: ri
  real(dp)                                         :: rj
!
!  Compute standard radii with coordination corrections
!
  ri = gfnff_rad_cn(1,nati) + gfnff_rad_cn(2,nati)*cni
  rj = gfnff_rad_cn(1,natj) + gfnff_rad_cn(2,natj)*cnj
!
  den = abs(gfnff_rad_cn(3,nati)-gfnff_rad_cn(3,natj))
  f1 = gfnff_rad_cn(4,nati) + gfnff_rad_cn(4,natj)
  f2 = gfnff_rad_cn(5,nati) + gfnff_rad_cn(5,natj)
  f12 = 1.0_dp - f1*den - f2*den**2
  radij = (radij + ri + rj)*f12
!
  if (lgrad1) then
    dradijdcni = gfnff_rad_cn(2,nati)*f12
    dradijdcnj = gfnff_rad_cn(2,natj)*f12
  endif
!
  return
  end
