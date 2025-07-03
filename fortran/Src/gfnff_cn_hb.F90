  subroutine gfnff_cn_hb(cn)
!
!  Computes the coordination number for GFNFF for hydrogen bonds
!
!  On exit :
!
!  cn              = coordination number
!
!  10/20 Created
!  10/20 Handling of images for PBC when A=B corrected
!   2/21 Cutoff added based on error function going to zero
!
!  Julian Gale, CIC, Curtin University, February 2021
!
  use datatypes
  use current
  use gulp_gfnff
  use m_gfnff_nbr
  use m_gfnff_pairs
  implicit none
!
!  Passed variables
!
  real(dp),                         intent(out)    :: cn(numat)
!
!  Local variables
!
  integer(i4)                                      :: hbA
  integer(i4)                                      :: hbB
  integer(i4)                                      :: hbH
  integer(i4)                                      :: natB
  integer(i4)                                      :: natH
  integer(i4)                                      :: nh
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: np
  logical                                          :: lABok
  real(dp)                                         :: cncut2
  real(dp)                                         :: ctrm
  real(dp),                          parameter     :: rcov_scale = 1.78_dp
  real(dp)                                         :: r0
  real(dp)                                         :: rhb
  real(dp)                                         :: r2ab
  real(dp)                                         :: xab
  real(dp)                                         :: yab
  real(dp)                                         :: zab
!
!  Initialise lattice vectors for this cutoff
!
  call gfnff_setpaircell(gfnff_cnhbthr,cnhb_paircell)
!
!  Initialise coordination number
!
  cn(1:numat)  = 0.0_dp
!
  do ni = 1,nbond_hb_nr
    hbA = nbond_hb_AH(1,ni)
    hbH = nbond_hb_AH(2,ni)
    nh  = nbond_hb_AH(3,ni)
    natH = nat(hbH)
    do nj = 1,nbond_hb_Bn(ni)
      hbB = nbond_hb_B(nj,ni)
      natB = nat(hbB)
      r0 = rcov_scale*(gfnff_rcov(natH) + gfnff_rcov(natB))
!
!  Compute maximum radius for non-zero interactions
!
      cncut2 = min(r0*r0*gfnff_cnerfcut_hb,gfnff_cnhbthr)
!
!  Find interactions for this pair
!
      call gfnff_getpairs(hbH,hbB,cncut2,cnhb_paircell,cnhb_pairs)
!
!  Loop over valid interactions 
!
      do np = 1,cnhb_pairs%npair
!
!  Check that if B is an image of A then A=B case is excluded for the same image
!
        if (hbA.eq.hbB) then
          xab = cnhb_pairs%xpair(np) + xbnbr(nh,hbA)
          yab = cnhb_pairs%ypair(np) + ybnbr(nh,hbA)
          zab = cnhb_pairs%zpair(np) + zbnbr(nh,hbA)
          r2ab = xab**2 + yab**2 + zab**2
          lABok = (r2ab.gt.1.0d-2)
        else
          lABok = .true.
        endif
!
        if (lABok) then
          rhb = sqrt(cnhb_pairs%r2pair(np))
          ctrm = 0.5_dp*(1.0_dp + erf(-gfnff_kn_hb*(rhb-r0)/r0))
          cn(hbH) = cn(hbH) + ctrm
          cn(hbB) = cn(hbB) + ctrm
        endif
      enddo
    enddo
  enddo
!
  return
  end
