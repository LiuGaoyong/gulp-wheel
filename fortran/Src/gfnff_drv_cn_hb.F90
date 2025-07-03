!**********************************************************************************************************
!  Routines that add the HB coordination number derivatives from GFNFF to the appropriate arrays in GULP  *
!**********************************************************************************************************
  subroutine gfnff_drv_dcn_hb(i,dEdcni)
!
!  Computes the derivatives of the coordination number for GFNFF for hydrogen bonds for atom i
!
!  On entry : 
!
!  i               = atom whose coordination number is to contribute to the derivatives
!  dEdcni          = first derivative of energy w.r.t. coordination number of i
!
!  On exit :
!
!  Derivatives have been computed and added to the appropriate arrays
!
!  10/20 Created
!  10/20 Handling of images for PBC when A=B corrected
!  11/20 Separate routine created for first derivatives only
!   2/21 Cutoff added based on error function going to zero
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
!  Julian Gale, CIC, Curtin University, February 2021
!
  use datatypes
  use current
  use derivatives
  use g_constants,    only : pi
  use gulp_gfnff
  use iochannels
  use m_gfnff_nbr
  use m_gfnff_pairs
  use m_strain,       only : real1strterm
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: dEdcni
!
!  Local variables
!
  integer(i4)                                      :: hbA
  integer(i4)                                      :: hbB
  integer(i4)                                      :: hbH
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: natB
  integer(i4)                                      :: natH
  integer(i4)                                      :: nh
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: np
  logical                                          :: lABok
  real(dp)                                         :: cncut2
  real(dp)                                         :: dctrmdrhb
  real(dp)                                         :: dEdrhb
  real(dp)                                         :: dr2ds(6)
  real(dp)                                         :: d2r2dx2(6)
  real(dp)                                         :: d2r2ds2(6,6)
  real(dp)                                         :: d2r2dsdx(6,3)
  real(dp)                                         :: dtrm
  real(dp)                                         :: dr
  real(dp),                          parameter     :: rcov_scale = 1.78_dp
  real(dp)                                         :: r0
  real(dp)                                         :: rhb
  real(dp)                                         :: r2ab
  real(dp)                                         :: xab
  real(dp)                                         :: yab
  real(dp)                                         :: zab
#ifdef TRACE
  call trace_in('gfnff_drv_dcn_hb')
#endif
!
  do ni = 1,nbond_hb_nr
    hbH = nbond_hb_AH(2,ni)
! DEBUG - look for faster approach here!
    if (hbH.ne.i) cycle
    hbA = nbond_hb_AH(1,ni)
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
!
          dr = (rhb - r0)/r0
          dtrm = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr)**2)/(sqrt(pi)*r0)
          dctrmdrhb = - dtrm/rhb
!
          dEdrhb = dEdcni*dctrmdrhb
!
!  First derivatives
!
          xdrv(hbH) = xdrv(hbH) - dEdrhb*cnhb_pairs%xpair(np)
          ydrv(hbH) = ydrv(hbH) - dEdrhb*cnhb_pairs%ypair(np)
          zdrv(hbH) = zdrv(hbH) - dEdrhb*cnhb_pairs%zpair(np)
!
          xdrv(hbB) = xdrv(hbB) + dEdrhb*cnhb_pairs%xpair(np)
          ydrv(hbB) = ydrv(hbB) + dEdrhb*cnhb_pairs%ypair(np)
          zdrv(hbB) = zdrv(hbB) + dEdrhb*cnhb_pairs%zpair(np)
!
          if (lstr) then
            call real1strterm(ndim,cnhb_pairs%xpair(np),cnhb_pairs%ypair(np),cnhb_pairs%zpair(np),0.0_dp,0.0_dp,0.0_dp, &
                              dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          endif
          if (lstr) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + dEdrhb*dr2ds(ks)
            enddo
          endif
        endif
      enddo
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv_dcn_hb')
#endif
!
  return
  end
!
  subroutine gfnff_drv_dcnd_hb(ofctij,i,dEdcni)
!
!  Computes the derivatives of the coordination number for GFNFF for hydrogen bonds for atom i
!  Version for part of distributed memory parallel call
!
!  On entry : 
!
!  ofctij          = correction factor depending on i=j in calling routine
!  i               = atom whose coordination number is to contribute to the derivatives
!  dEdcni          = first derivative of energy w.r.t. coordination number of i
!
!  On exit :
!
!  Derivatives have been computed and added to the appropriate arrays
!
!   3/22 Created from gfnff_drv_dcn_hb
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
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
  use current
  use derivatives
  use g_constants,    only : pi
  use gulp_gfnff
  use iochannels
  use m_gfnff_nbr
  use m_gfnff_pairs
  use m_strain,       only : real1strterm
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: dEdcni
  real(dp),    intent(in)                          :: ofctij
!
!  Local variables
!
  integer(i4)                                      :: hbA
  integer(i4)                                      :: hbB
  integer(i4)                                      :: hbH
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: natB
  integer(i4)                                      :: natH
  integer(i4)                                      :: nh
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: np
  logical                                          :: lABok
  real(dp)                                         :: cncut2
  real(dp)                                         :: dctrmdrhb
  real(dp)                                         :: dEdrhb
  real(dp)                                         :: dr2ds(6)
  real(dp)                                         :: d2r2dx2(6)
  real(dp)                                         :: d2r2ds2(6,6)
  real(dp)                                         :: d2r2dsdx(6,3)
  real(dp)                                         :: dtrm
  real(dp)                                         :: dr
  real(dp),                          parameter     :: rcov_scale = 1.78_dp
  real(dp)                                         :: r0
  real(dp)                                         :: rhb
  real(dp)                                         :: r2ab
  real(dp)                                         :: xab
  real(dp)                                         :: yab
  real(dp)                                         :: zab
#ifdef TRACE
  call trace_in('gfnff_drv_dcnd_hb')
#endif
!
  do ni = 1,nbond_hb_nr
    hbH = nbond_hb_AH(2,ni)
! DEBUG - look for faster approach here!
    if (hbH.ne.i) cycle
    hbA = nbond_hb_AH(1,ni)
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
!
          dr = (rhb - r0)/r0
          dtrm = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr)**2)/(sqrt(pi)*r0)
          dctrmdrhb = - dtrm/rhb
!
          dEdrhb = 0.5_dp*dEdcni*dctrmdrhb*ofctij
!
!  First derivatives
!
          xdrv(hbH) = xdrv(hbH) - dEdrhb*cnhb_pairs%xpair(np)
          ydrv(hbH) = ydrv(hbH) - dEdrhb*cnhb_pairs%ypair(np)
          zdrv(hbH) = zdrv(hbH) - dEdrhb*cnhb_pairs%zpair(np)
!
          xdrv(hbB) = xdrv(hbB) + dEdrhb*cnhb_pairs%xpair(np)
          ydrv(hbB) = ydrv(hbB) + dEdrhb*cnhb_pairs%ypair(np)
          zdrv(hbB) = zdrv(hbB) + dEdrhb*cnhb_pairs%zpair(np)
!
          if (lstr) then
            call real1strterm(ndim,cnhb_pairs%xpair(np),cnhb_pairs%ypair(np),cnhb_pairs%zpair(np),0.0_dp,0.0_dp,0.0_dp, &
                              dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
          endif
          if (lstr) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + dEdrhb*dr2ds(ks)
            enddo
          endif
        endif
      enddo
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv_dcnd_hb')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_hb(hbH,i,j,xij,yij,zij,dEdcnh,d2Edcnh2,d2Edcnhdr,d2Edcnhdcni,d2Edcnhdcnj,dcn,lgrad2)
!
!  Computes the derivatives of the coordination number for GFNFF for hydrogen bonds for atom i-j pair
!
!  On entry : 
!
!  hbH             = hydrogen atom whose coordination number is to contribute to the derivatives
!  i               = atom whose coordination number is to contribute to the derivatives via rij
!  j               = atom whose coordination number is to contribute to the derivatives via rij
!  xij             = x component of vector from i to j
!  yij             = y component of vector from i to j
!  zij             = z component of vector from i to j
!  dEdcnh          = first derivative of energy w.r.t. HB coordination number
!  d2Edcnh2        = second derivative of energy w.r.t. HB coordination number
!  d2Edcnhdr       = second derivative of energy w.r.t. HB coordination number and rij
!  d2Edcnhdcni     = second derivative of energy w.r.t. coordination numbers of hbH and i
!  d2Edcnhdcnj     = second derivative of energy w.r.t. coordination numbers of hbH and j
!  dcn             = derivative of log coordination numbers w.r.t. coordination numbers
!  lgrad2          = if .true. then compute the second derivatives
!
!  On exit :
!
!  Derivatives have been computed and added to the appropriate arrays
!
!  11/20 Created from gfnff_drv_dcn_hb and gfnff_drv2_dcn
!   2/21 Correction to one set of second derivatives
!   2/21 Cutoff added based on error function going to zero
!   2/21 Second derivative term corrected for d2ctrmdrhb2
!   2/21 cut no longer passed in for coordination number term
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
!  Julian Gale, CIC, Curtin University, February 2021
!
  use datatypes
  use current
  use derivatives
  use g_constants,    only : pi
  use gulp_gfnff
  use iochannels
  use m_gfnff_nbr
  use m_gfnff_pairs
  use m_strain,       only : real1strterm
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: hbH
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: dEdcnh
  real(dp),    intent(in)                          :: d2Edcnh2
  real(dp),    intent(in)                          :: d2Edcnhdr
  real(dp),    intent(in)                          :: d2Edcnhdcni
  real(dp),    intent(in)                          :: d2Edcnhdcnj
  real(dp),    intent(in)                          :: dcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: hbA
  integer(i4)                                      :: hbB
  integer(i4)                                      :: hbBx
  integer(i4)                                      :: hbBy
  integer(i4)                                      :: hbBz
  integer(i4)                                      :: hbHx
  integer(i4)                                      :: hbHy
  integer(i4)                                      :: hbHz
  integer(i4)                                      :: indh
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: natB
  integer(i4)                                      :: natH
  integer(i4)                                      :: nh
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nk
  integer(i4)                                      :: np
  integer(i4)                                      :: np2
  logical                                          :: lABok
  logical                                          :: lABok2
  real(dp)                                         :: cncut2
  real(dp)                                         :: dctrmdrhb
  real(dp)                                         :: dctrmdrhb2
  real(dp)                                         :: d2ctrmdrhb2
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: dEdrhb
  real(dp)                                         :: d2Edrhb2
  real(dp)                                         :: dr2ds(6)
  real(dp)                                         :: dr2ds2(6)
  real(dp)                                         :: d2r2dx2(6)
  real(dp)                                         :: d2r2ds2(6,6)
  real(dp)                                         :: d2r2dsdx(6,3)
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
  real(dp)                                         :: dtrm
  real(dp)                                         :: dtrm1a
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: dr
  real(dp)                                         :: dr2
  real(dp)                                         :: drik
  real(dp)                                         :: drjk
  real(dp)                                         :: ofctij
  real(dp),                          parameter     :: rcov_scale = 1.78_dp
  real(dp)                                         :: r0
  real(dp)                                         :: r0ik
  real(dp)                                         :: r0jk
  real(dp)                                         :: rik
  real(dp)                                         :: rjk
  real(dp)                                         :: rhb
  real(dp)                                         :: rhb2
  real(dp)                                         :: r2ab
  real(dp)                                         :: r2ab2
  real(dp)                                         :: xab
  real(dp)                                         :: yab
  real(dp)                                         :: zab
  real(dp)                                         :: xab2
  real(dp)                                         :: yab2
  real(dp)                                         :: zab2
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_hb')
#endif
!
!  Set up rij related terms if needed
!
  if (lgrad2.and.lstr) then
    call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2,d2r2ijdsdx,d2r2ijds2,.false.)
  endif
!
!  Set up for second derivatives
!
  if (lgrad2) then
    indi = 3*(i-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
!
    indj = 3*(j-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
!
    indh = 3*(hbH-1)
    hbHx = indh + 1
    hbHy = indh + 2
    hbHz = indh + 3
!
    dlogcnidcni = dcn(i)
    dlogcnjdcnj = dcn(j)
  endif
!
!  Set factor for double counting
!
  if (i.eq.j) then
    ofctij = 1.0_dp
  else
    ofctij = 0.5_dp
  endif
!
  do ni = 1,nbond_hb_nr
! DEBUG - look for faster approach here!
    if (hbH.ne.nbond_hb_AH(2,ni)) cycle
    hbA = nbond_hb_AH(1,ni)
    nh  = nbond_hb_AH(3,ni)
    natH = nat(hbH)
    do nj = 1,nbond_hb_Bn(ni)
      hbB = nbond_hb_B(nj,ni)
      if (lgrad2) then
        indh = 3*(hbB-1)
        hbBx = indh + 1
        hbBy = indh + 2
        hbBz = indh + 3
      endif
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
!
          dr = (rhb - r0)/r0
          dtrm = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr)**2)/(sqrt(pi)*r0)
          dctrmdrhb = - dtrm/rhb
!
          dEdrhb = dEdcnh*dctrmdrhb
!----------------------
!  First derivatives  |
!----------------------
          xdrv(hbH) = xdrv(hbH) - dEdrhb*cnhb_pairs%xpair(np)*ofctij
          ydrv(hbH) = ydrv(hbH) - dEdrhb*cnhb_pairs%ypair(np)*ofctij
          zdrv(hbH) = zdrv(hbH) - dEdrhb*cnhb_pairs%zpair(np)*ofctij
!
          xdrv(hbB) = xdrv(hbB) + dEdrhb*cnhb_pairs%xpair(np)*ofctij
          ydrv(hbB) = ydrv(hbB) + dEdrhb*cnhb_pairs%ypair(np)*ofctij
          zdrv(hbB) = zdrv(hbB) + dEdrhb*cnhb_pairs%zpair(np)*ofctij
!
          if (lstr.or.lgrad2) then
            call real1strterm(ndim,cnhb_pairs%xpair(np),cnhb_pairs%ypair(np),cnhb_pairs%zpair(np),0.0_dp,0.0_dp,0.0_dp, &
                              dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,lgrad2)
          endif
          if (lstr) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + dEdrhb*dr2ds(ks)*ofctij
            enddo
          endif
!-----------------------
!  Second derivatives  |
!-----------------------
          if (lgrad2) then
!
!  HB coordination number only
!
            d2ctrmdrhb2 = dtrm*(1.0_dp + 2.0_dp*rhb*dr*(gfnff_kn_hb**2)/r0)/rhb**3
            d2Edrhb2 = d2Edcnh2*dctrmdrhb*dctrmdrhb + dEdcnh*d2ctrmdrhb2 ! (d2E/dr2)
!
            if (hbB.gt.hbH) then
              derv2(hbBx,hbHx) = derv2(hbBx,hbHx) - d2Edrhb2*d2r2dx2(1)
              derv2(hbBy,hbHx) = derv2(hbBy,hbHx) - d2Edrhb2*d2r2dx2(6)
              derv2(hbBz,hbHx) = derv2(hbBz,hbHx) - d2Edrhb2*d2r2dx2(5)
              derv2(hbBx,hbHy) = derv2(hbBx,hbHy) - d2Edrhb2*d2r2dx2(6)
              derv2(hbBy,hbHy) = derv2(hbBy,hbHy) - d2Edrhb2*d2r2dx2(2)
              derv2(hbBz,hbHy) = derv2(hbBz,hbHy) - d2Edrhb2*d2r2dx2(4)
              derv2(hbBx,hbHz) = derv2(hbBx,hbHz) - d2Edrhb2*d2r2dx2(5)
              derv2(hbBy,hbHz) = derv2(hbBy,hbHz) - d2Edrhb2*d2r2dx2(4)
              derv2(hbBz,hbHz) = derv2(hbBz,hbHz) - d2Edrhb2*d2r2dx2(3)
              derv2(hbBx,hbHx) = derv2(hbBx,hbHx) - dEdrhb
              derv2(hbBy,hbHy) = derv2(hbBy,hbHy) - dEdrhb
              derv2(hbBz,hbHz) = derv2(hbBz,hbHz) - dEdrhb
            else
              derv2(hbHx,hbBx) = derv2(hbHx,hbBx) - d2Edrhb2*d2r2dx2(1)
              derv2(hbHy,hbBx) = derv2(hbHy,hbBx) - d2Edrhb2*d2r2dx2(6)
              derv2(hbHz,hbBx) = derv2(hbHz,hbBx) - d2Edrhb2*d2r2dx2(5)
              derv2(hbHx,hbBy) = derv2(hbHx,hbBy) - d2Edrhb2*d2r2dx2(6)
              derv2(hbHy,hbBy) = derv2(hbHy,hbBy) - d2Edrhb2*d2r2dx2(2)
              derv2(hbHz,hbBy) = derv2(hbHz,hbBy) - d2Edrhb2*d2r2dx2(4)
              derv2(hbHx,hbBz) = derv2(hbHx,hbBz) - d2Edrhb2*d2r2dx2(5)
              derv2(hbHy,hbBz) = derv2(hbHy,hbBz) - d2Edrhb2*d2r2dx2(4)
              derv2(hbHz,hbBz) = derv2(hbHz,hbBz) - d2Edrhb2*d2r2dx2(3)
              derv2(hbHx,hbBx) = derv2(hbHx,hbBx) - dEdrhb
              derv2(hbHy,hbBy) = derv2(hbHy,hbBy) - dEdrhb
              derv2(hbHz,hbBz) = derv2(hbHz,hbBz) - dEdrhb
            endif
!
            if (lstr) then
!
!  Mixed derivatives
!
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(hbHx,kl) = derv3(hbHx,kl) - dEdrhb*d2r2dsdx(ks,1) - cnhb_pairs%xpair(np)*d2Edrhb2*dr2ds(ks)
                derv3(hbHy,kl) = derv3(hbHy,kl) - dEdrhb*d2r2dsdx(ks,2) - cnhb_pairs%ypair(np)*d2Edrhb2*dr2ds(ks)
                derv3(hbHz,kl) = derv3(hbHz,kl) - dEdrhb*d2r2dsdx(ks,3) - cnhb_pairs%zpair(np)*d2Edrhb2*dr2ds(ks)
                derv3(hbBx,kl) = derv3(hbBx,kl) + dEdrhb*d2r2dsdx(ks,1) + cnhb_pairs%xpair(np)*d2Edrhb2*dr2ds(ks)
                derv3(hbBy,kl) = derv3(hbBy,kl) + dEdrhb*d2r2dsdx(ks,2) + cnhb_pairs%ypair(np)*d2Edrhb2*dr2ds(ks)
                derv3(hbBz,kl) = derv3(hbBz,kl) + dEdrhb*d2r2dsdx(ks,3) + cnhb_pairs%zpair(np)*d2Edrhb2*dr2ds(ks)
              enddo
!
!  Strain-strain
!
              do kk = 1,nstrains
                ks = nstrptr(kk)
                do kl = 1,nstrains
                  kt = nstrptr(kl)
                  sderv2(kl,kk) = sderv2(kl,kk) + d2Edrhb2*dr2ds(kt)*dr2ds(ks) + dEdrhb*d2r2ds2(kt,ks)
                enddo
              enddo
            endif
!
            if (np.gt.1) then
!
!  Loop over pairs of valid interactions
!
              do np2 = 1,np-1
!
!  Check that if B is an image of A then A=B case is excluded for the same image
!
                if (hbA.eq.hbB) then
                  xab2 = cnhb_pairs%xpair(np2) + xbnbr(nh,hbA)
                  yab2 = cnhb_pairs%ypair(np2) + ybnbr(nh,hbA)
                  zab2 = cnhb_pairs%zpair(np2) + zbnbr(nh,hbA)
                  r2ab2 = xab2**2 + yab2**2 + zab2**2
                  lABok2 = (r2ab2.gt.1.0d-2)
                else
                  lABok2 = .true.
                endif
!
                if (lABok2) then
                  rhb2 = sqrt(cnhb_pairs%r2pair(np2))
!
                  dr2 = (rhb2 - r0)/r0
                  dtrm1a = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr2)**2)/(sqrt(pi)*r0)
                  dctrmdrhb2 = - dtrm1a/rhb2
                  d2Edrhb2 = d2Edcnh2*dctrmdrhb*dctrmdrhb2
!
                  if (hbB.gt.hbH) then
                    derv2(hbBx,hbHx) = derv2(hbBx,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBy,hbHx) = derv2(hbBy,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBz,hbHx) = derv2(hbBz,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBx,hbHy) = derv2(hbBx,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBy,hbHy) = derv2(hbBy,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBz,hbHy) = derv2(hbBz,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBx,hbHz) = derv2(hbBx,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBy,hbHz) = derv2(hbBy,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBz,hbHz) = derv2(hbBz,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  else
                    derv2(hbHx,hbBx) = derv2(hbHx,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHy,hbBx) = derv2(hbHy,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHz,hbBx) = derv2(hbHz,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHx,hbBy) = derv2(hbHx,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHy,hbBy) = derv2(hbHy,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHz,hbBy) = derv2(hbHz,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHx,hbBz) = derv2(hbHx,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHy,hbBz) = derv2(hbHy,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHz,hbBz) = derv2(hbHz,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  endif
!
                  if (lstr) then
                    call real1strterm(ndim,cnhb_pairs%xpair(np2),cnhb_pairs%ypair(np2),cnhb_pairs%zpair(np2), &
                                      0.0_dp,0.0_dp,0.0_dp,dr2ds2,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
!  Mixed derivatives
!
                    do kl = 1,nstrains
                      ks = nstrptr(kl)
                      derv3(hbHx,kl) = derv3(hbHx,kl) - cnhb_pairs%xpair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                      - cnhb_pairs%xpair(np2)*d2Edrhb2*dr2ds(ks)
                      derv3(hbHy,kl) = derv3(hbHy,kl) - cnhb_pairs%ypair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                      - cnhb_pairs%ypair(np2)*d2Edrhb2*dr2ds(ks)
                      derv3(hbHz,kl) = derv3(hbHz,kl) - cnhb_pairs%zpair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                      - cnhb_pairs%zpair(np2)*d2Edrhb2*dr2ds(ks)
                      derv3(hbBx,kl) = derv3(hbBx,kl) + cnhb_pairs%xpair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                      + cnhb_pairs%xpair(np2)*d2Edrhb2*dr2ds(ks)
                      derv3(hbBy,kl) = derv3(hbBy,kl) + cnhb_pairs%ypair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                      + cnhb_pairs%ypair(np2)*d2Edrhb2*dr2ds(ks)
                      derv3(hbBz,kl) = derv3(hbBz,kl) + cnhb_pairs%zpair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                      + cnhb_pairs%zpair(np2)*d2Edrhb2*dr2ds(ks)
                    enddo
!
!  Strain-strain
!
                    do kk = 1,nstrains
                      ks = nstrptr(kk)
                      do kl = 1,nstrains
                        kt = nstrptr(kl)
                        sderv2(kl,kk) = sderv2(kl,kk) + d2Edrhb2*dr2ds2(kt)*dr2ds(ks)
                      enddo
                    enddo
                  endif
                endif
              enddo
            endif
!
!  HB coordination number and rij
!
            d2trm = d2Edcnhdr*dctrmdrhb
!
            if (hbH.ge.i) then
              derv2(hbHx,ix) = derv2(hbHx,ix) + d2trm*cnhb_pairs%xpair(np)*xij
              derv2(hbHy,ix) = derv2(hbHy,ix) + d2trm*cnhb_pairs%ypair(np)*xij
              derv2(hbHz,ix) = derv2(hbHz,ix) + d2trm*cnhb_pairs%zpair(np)*xij
              derv2(hbHx,iy) = derv2(hbHx,iy) + d2trm*cnhb_pairs%xpair(np)*yij
              derv2(hbHy,iy) = derv2(hbHy,iy) + d2trm*cnhb_pairs%ypair(np)*yij
              derv2(hbHz,iy) = derv2(hbHz,iy) + d2trm*cnhb_pairs%zpair(np)*yij
              derv2(hbHx,iz) = derv2(hbHx,iz) + d2trm*cnhb_pairs%xpair(np)*zij
              derv2(hbHy,iz) = derv2(hbHy,iz) + d2trm*cnhb_pairs%ypair(np)*zij
              derv2(hbHz,iz) = derv2(hbHz,iz) + d2trm*cnhb_pairs%zpair(np)*zij
            else
              derv2(ix,hbHx) = derv2(ix,hbHx) + d2trm*cnhb_pairs%xpair(np)*xij
              derv2(iy,hbHx) = derv2(iy,hbHx) + d2trm*cnhb_pairs%xpair(np)*yij
              derv2(iz,hbHx) = derv2(iz,hbHx) + d2trm*cnhb_pairs%xpair(np)*zij
              derv2(ix,hbHy) = derv2(ix,hbHy) + d2trm*cnhb_pairs%ypair(np)*xij
              derv2(iy,hbHy) = derv2(iy,hbHy) + d2trm*cnhb_pairs%ypair(np)*yij
              derv2(iz,hbHy) = derv2(iz,hbHy) + d2trm*cnhb_pairs%ypair(np)*zij
              derv2(ix,hbHz) = derv2(ix,hbHz) + d2trm*cnhb_pairs%zpair(np)*xij
              derv2(iy,hbHz) = derv2(iy,hbHz) + d2trm*cnhb_pairs%zpair(np)*yij
              derv2(iz,hbHz) = derv2(iz,hbHz) + d2trm*cnhb_pairs%zpair(np)*zij
            endif
!
            if (hbH.ge.j) then
              derv2(hbHx,jx) = derv2(hbHx,jx) - d2trm*cnhb_pairs%xpair(np)*xij
              derv2(hbHy,jx) = derv2(hbHy,jx) - d2trm*cnhb_pairs%ypair(np)*xij
              derv2(hbHz,jx) = derv2(hbHz,jx) - d2trm*cnhb_pairs%zpair(np)*xij
              derv2(hbHx,jy) = derv2(hbHx,jy) - d2trm*cnhb_pairs%xpair(np)*yij
              derv2(hbHy,jy) = derv2(hbHy,jy) - d2trm*cnhb_pairs%ypair(np)*yij
              derv2(hbHz,jy) = derv2(hbHz,jy) - d2trm*cnhb_pairs%zpair(np)*yij
              derv2(hbHx,jz) = derv2(hbHx,jz) - d2trm*cnhb_pairs%xpair(np)*zij
              derv2(hbHy,jz) = derv2(hbHy,jz) - d2trm*cnhb_pairs%ypair(np)*zij
              derv2(hbHz,jz) = derv2(hbHz,jz) - d2trm*cnhb_pairs%zpair(np)*zij
            else
              derv2(jx,hbHx) = derv2(jx,hbHx) - d2trm*cnhb_pairs%xpair(np)*xij
              derv2(jy,hbHx) = derv2(jy,hbHx) - d2trm*cnhb_pairs%xpair(np)*yij
              derv2(jz,hbHx) = derv2(jz,hbHx) - d2trm*cnhb_pairs%xpair(np)*zij
              derv2(jx,hbHy) = derv2(jx,hbHy) - d2trm*cnhb_pairs%ypair(np)*xij
              derv2(jy,hbHy) = derv2(jy,hbHy) - d2trm*cnhb_pairs%ypair(np)*yij
              derv2(jz,hbHy) = derv2(jz,hbHy) - d2trm*cnhb_pairs%ypair(np)*zij
              derv2(jx,hbHz) = derv2(jx,hbHz) - d2trm*cnhb_pairs%zpair(np)*xij
              derv2(jy,hbHz) = derv2(jy,hbHz) - d2trm*cnhb_pairs%zpair(np)*yij
              derv2(jz,hbHz) = derv2(jz,hbHz) - d2trm*cnhb_pairs%zpair(np)*zij
            endif
!
            if (hbB.ge.i) then
              derv2(hbBx,ix) = derv2(hbBx,ix) - d2trm*cnhb_pairs%xpair(np)*xij
              derv2(hbBy,ix) = derv2(hbBy,ix) - d2trm*cnhb_pairs%ypair(np)*xij
              derv2(hbBz,ix) = derv2(hbBz,ix) - d2trm*cnhb_pairs%zpair(np)*xij
              derv2(hbBx,iy) = derv2(hbBx,iy) - d2trm*cnhb_pairs%xpair(np)*yij
              derv2(hbBy,iy) = derv2(hbBy,iy) - d2trm*cnhb_pairs%ypair(np)*yij
              derv2(hbBz,iy) = derv2(hbBz,iy) - d2trm*cnhb_pairs%zpair(np)*yij
              derv2(hbBx,iz) = derv2(hbBx,iz) - d2trm*cnhb_pairs%xpair(np)*zij
              derv2(hbBy,iz) = derv2(hbBy,iz) - d2trm*cnhb_pairs%ypair(np)*zij
              derv2(hbBz,iz) = derv2(hbBz,iz) - d2trm*cnhb_pairs%zpair(np)*zij
            else
              derv2(ix,hbBx) = derv2(ix,hbBx) - d2trm*cnhb_pairs%xpair(np)*xij
              derv2(iy,hbBx) = derv2(iy,hbBx) - d2trm*cnhb_pairs%xpair(np)*yij
              derv2(iz,hbBx) = derv2(iz,hbBx) - d2trm*cnhb_pairs%xpair(np)*zij
              derv2(ix,hbBy) = derv2(ix,hbBy) - d2trm*cnhb_pairs%ypair(np)*xij
              derv2(iy,hbBy) = derv2(iy,hbBy) - d2trm*cnhb_pairs%ypair(np)*yij
              derv2(iz,hbBy) = derv2(iz,hbBy) - d2trm*cnhb_pairs%ypair(np)*zij
              derv2(ix,hbBz) = derv2(ix,hbBz) - d2trm*cnhb_pairs%zpair(np)*xij
              derv2(iy,hbBz) = derv2(iy,hbBz) - d2trm*cnhb_pairs%zpair(np)*yij
              derv2(iz,hbBz) = derv2(iz,hbBz) - d2trm*cnhb_pairs%zpair(np)*zij
            endif
!
            if (hbB.ge.j) then
              derv2(hbBx,jx) = derv2(hbBx,jx) + d2trm*cnhb_pairs%xpair(np)*xij
              derv2(hbBy,jx) = derv2(hbBy,jx) + d2trm*cnhb_pairs%ypair(np)*xij
              derv2(hbBz,jx) = derv2(hbBz,jx) + d2trm*cnhb_pairs%zpair(np)*xij
              derv2(hbBx,jy) = derv2(hbBx,jy) + d2trm*cnhb_pairs%xpair(np)*yij
              derv2(hbBy,jy) = derv2(hbBy,jy) + d2trm*cnhb_pairs%ypair(np)*yij
              derv2(hbBz,jy) = derv2(hbBz,jy) + d2trm*cnhb_pairs%zpair(np)*yij
              derv2(hbBx,jz) = derv2(hbBx,jz) + d2trm*cnhb_pairs%xpair(np)*zij
              derv2(hbBy,jz) = derv2(hbBy,jz) + d2trm*cnhb_pairs%ypair(np)*zij
              derv2(hbBz,jz) = derv2(hbBz,jz) + d2trm*cnhb_pairs%zpair(np)*zij
            else
              derv2(jx,hbBx) = derv2(jx,hbBx) + d2trm*cnhb_pairs%xpair(np)*xij
              derv2(jy,hbBx) = derv2(jy,hbBx) + d2trm*cnhb_pairs%xpair(np)*yij
              derv2(jz,hbBx) = derv2(jz,hbBx) + d2trm*cnhb_pairs%xpair(np)*zij
              derv2(jx,hbBy) = derv2(jx,hbBy) + d2trm*cnhb_pairs%ypair(np)*xij
              derv2(jy,hbBy) = derv2(jy,hbBy) + d2trm*cnhb_pairs%ypair(np)*yij
              derv2(jz,hbBy) = derv2(jz,hbBy) + d2trm*cnhb_pairs%ypair(np)*zij
              derv2(jx,hbBz) = derv2(jx,hbBz) + d2trm*cnhb_pairs%zpair(np)*xij
              derv2(jy,hbBz) = derv2(jy,hbBz) + d2trm*cnhb_pairs%zpair(np)*yij
              derv2(jz,hbBz) = derv2(jz,hbBz) + d2trm*cnhb_pairs%zpair(np)*zij
            endif
!
!  Strain terms
!
            if (lstr) then
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(ix,kl) = derv3(ix,kl) - xij*d2trm*dr2ds(ks)
                derv3(iy,kl) = derv3(iy,kl) - yij*d2trm*dr2ds(ks)
                derv3(iz,kl) = derv3(iz,kl) - zij*d2trm*dr2ds(ks)
                derv3(jx,kl) = derv3(jx,kl) + xij*d2trm*dr2ds(ks)
                derv3(jy,kl) = derv3(jy,kl) + yij*d2trm*dr2ds(ks)
                derv3(jz,kl) = derv3(jz,kl) + zij*d2trm*dr2ds(ks)
              enddo
!
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(hbHx,kl) = derv3(hbHx,kl) - cnhb_pairs%xpair(np)*d2trm*dr2ijds(ks)
                derv3(hbHy,kl) = derv3(hbHy,kl) - cnhb_pairs%ypair(np)*d2trm*dr2ijds(ks)
                derv3(hbHz,kl) = derv3(hbHz,kl) - cnhb_pairs%zpair(np)*d2trm*dr2ijds(ks)
                derv3(hbBx,kl) = derv3(hbBx,kl) + cnhb_pairs%xpair(np)*d2trm*dr2ijds(ks)
                derv3(hbBy,kl) = derv3(hbBy,kl) + cnhb_pairs%ypair(np)*d2trm*dr2ijds(ks)
                derv3(hbBz,kl) = derv3(hbBz,kl) + cnhb_pairs%zpair(np)*d2trm*dr2ijds(ks)
              enddo
!
              do kk = 1,nstrains
                ks = nstrptr(kk)
                do kl = 1,nstrains
                  kt = nstrptr(kl)
                  sderv2(kl,kk) = sderv2(kl,kk) + d2trm*(dr2ds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2ds(ks))
                enddo
              enddo
            endif
!
!  HB coordination number and coordination number of i
!
            do n = 1,nnbr_cn(i)
              nk = nbrno_cn(n,i)
              rik = rnbr(nk,i)
              k = nbrno(nk,i)
              r0ik = gfnff_rcov(nat(i)) + gfnff_rcov(nat(k))
              drik = (rik - r0ik)/r0ik
              dcnidrik = gfnff_kn_cn*exp(-(gfnff_kn_cn*drik)**2)/(sqrt(pi)*r0ik*rik)
!
              indk = 3*(k-1)
              kx = indk + 1
              ky = indk + 2
              kz = indk + 3
!
              d1trm = d2Edcnhdcni*dcnidrik*dctrmdrhb*dlogcnidcni
!
              if (lstr) then
                call real1strterm(ndim,xnbr(nk,i),ynbr(nk,i),znbr(nk,i),0.0_dp,0.0_dp,0.0_dp, &
                                  dr2ikds,d2r2ikdx2,d2r2ikdsdx,d2r2ikds2,.false.)
              endif
!
              if (hbH.ge.i) then
                derv2(hbHx,ix) = derv2(hbHx,ix) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
                derv2(hbHy,ix) = derv2(hbHy,ix) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
                derv2(hbHz,ix) = derv2(hbHz,ix) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
                derv2(hbHx,iy) = derv2(hbHx,iy) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
                derv2(hbHy,iy) = derv2(hbHy,iy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
                derv2(hbHz,iy) = derv2(hbHz,iy) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
                derv2(hbHx,iz) = derv2(hbHx,iz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
                derv2(hbHy,iz) = derv2(hbHy,iz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
                derv2(hbHz,iz) = derv2(hbHz,iz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
              else
                derv2(ix,hbHx) = derv2(ix,hbHx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
                derv2(iy,hbHx) = derv2(iy,hbHx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
                derv2(iz,hbHx) = derv2(iz,hbHx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
                derv2(ix,hbHy) = derv2(ix,hbHy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
                derv2(iy,hbHy) = derv2(iy,hbHy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
                derv2(iz,hbHy) = derv2(iz,hbHy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
                derv2(ix,hbHz) = derv2(ix,hbHz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
                derv2(iy,hbHz) = derv2(iy,hbHz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
                derv2(iz,hbHz) = derv2(iz,hbHz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
              endif
!
              if (hbH.ge.k) then
                derv2(hbHx,kx) = derv2(hbHx,kx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
                derv2(hbHy,kx) = derv2(hbHy,kx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
                derv2(hbHz,kx) = derv2(hbHz,kx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
                derv2(hbHx,ky) = derv2(hbHx,ky) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
                derv2(hbHy,ky) = derv2(hbHy,ky) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
                derv2(hbHz,ky) = derv2(hbHz,ky) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
                derv2(hbHx,kz) = derv2(hbHx,kz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
                derv2(hbHy,kz) = derv2(hbHy,kz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
                derv2(hbHz,kz) = derv2(hbHz,kz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
              else
                derv2(kx,hbHx) = derv2(kx,hbHx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
                derv2(ky,hbHx) = derv2(ky,hbHx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
                derv2(kz,hbHx) = derv2(kz,hbHx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
                derv2(kx,hbHy) = derv2(kx,hbHy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
                derv2(ky,hbHy) = derv2(ky,hbHy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
                derv2(kz,hbHy) = derv2(kz,hbHy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
                derv2(kx,hbHz) = derv2(kx,hbHz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
                derv2(ky,hbHz) = derv2(ky,hbHz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
                derv2(kz,hbHz) = derv2(kz,hbHz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
              endif
!
              if (hbB.ge.i) then
                derv2(hbBx,ix) = derv2(hbBx,ix) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
                derv2(hbBy,ix) = derv2(hbBy,ix) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
                derv2(hbBz,ix) = derv2(hbBz,ix) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
                derv2(hbBx,iy) = derv2(hbBx,iy) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
                derv2(hbBy,iy) = derv2(hbBy,iy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
                derv2(hbBz,iy) = derv2(hbBz,iy) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
                derv2(hbBx,iz) = derv2(hbBx,iz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
                derv2(hbBy,iz) = derv2(hbBy,iz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
                derv2(hbBz,iz) = derv2(hbBz,iz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
              else
                derv2(ix,hbBx) = derv2(ix,hbBx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
                derv2(iy,hbBx) = derv2(iy,hbBx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
                derv2(iz,hbBx) = derv2(iz,hbBx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
                derv2(ix,hbBy) = derv2(ix,hbBy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
                derv2(iy,hbBy) = derv2(iy,hbBy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
                derv2(iz,hbBy) = derv2(iz,hbBy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
                derv2(ix,hbBz) = derv2(ix,hbBz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
                derv2(iy,hbBz) = derv2(iy,hbBz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
                derv2(iz,hbBz) = derv2(iz,hbBz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
              endif
!
              if (hbB.ge.k) then
                derv2(hbBx,kx) = derv2(hbBx,kx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
                derv2(hbBy,kx) = derv2(hbBy,kx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
                derv2(hbBz,kx) = derv2(hbBz,kx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
                derv2(hbBx,ky) = derv2(hbBx,ky) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
                derv2(hbBy,ky) = derv2(hbBy,ky) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
                derv2(hbBz,ky) = derv2(hbBz,ky) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
                derv2(hbBx,kz) = derv2(hbBx,kz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
                derv2(hbBy,kz) = derv2(hbBy,kz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
                derv2(hbBz,kz) = derv2(hbBz,kz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
              else
                derv2(kx,hbBx) = derv2(kx,hbBx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
                derv2(ky,hbBx) = derv2(ky,hbBx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
                derv2(kz,hbBx) = derv2(kz,hbBx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
                derv2(kx,hbBy) = derv2(kx,hbBy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
                derv2(ky,hbBy) = derv2(ky,hbBy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
                derv2(kz,hbBy) = derv2(kz,hbBy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
                derv2(kx,hbBz) = derv2(kx,hbBz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
                derv2(ky,hbBz) = derv2(ky,hbBz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
                derv2(kz,hbBz) = derv2(kz,hbBz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
              endif
!
!  Strain terms
!
              if (lstr) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(ix,kl) = derv3(ix,kl) - xnbr(nk,i)*d1trm*dr2ds(ks)
                  derv3(iy,kl) = derv3(iy,kl) - ynbr(nk,i)*d1trm*dr2ds(ks)
                  derv3(iz,kl) = derv3(iz,kl) - znbr(nk,i)*d1trm*dr2ds(ks)
                  derv3(kx,kl) = derv3(kx,kl) + xnbr(nk,i)*d1trm*dr2ds(ks)
                  derv3(ky,kl) = derv3(ky,kl) + ynbr(nk,i)*d1trm*dr2ds(ks)
                  derv3(kz,kl) = derv3(kz,kl) + znbr(nk,i)*d1trm*dr2ds(ks)
                enddo
!
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(hbHx,kl) = derv3(hbHx,kl) - cnhb_pairs%xpair(np)*d1trm*dr2ikds(ks)
                  derv3(hbHy,kl) = derv3(hbHy,kl) - cnhb_pairs%ypair(np)*d1trm*dr2ikds(ks)
                  derv3(hbHz,kl) = derv3(hbHz,kl) - cnhb_pairs%zpair(np)*d1trm*dr2ikds(ks)
                  derv3(hbBx,kl) = derv3(hbBx,kl) + cnhb_pairs%xpair(np)*d1trm*dr2ikds(ks)
                  derv3(hbBy,kl) = derv3(hbBy,kl) + cnhb_pairs%ypair(np)*d1trm*dr2ikds(ks)
                  derv3(hbBz,kl) = derv3(hbBz,kl) + cnhb_pairs%zpair(np)*d1trm*dr2ikds(ks)
                enddo
!
                do kk = 1,nstrains
                  ks = nstrptr(kk)
                  do kl = 1,nstrains
                    kt = nstrptr(kl)
                    sderv2(kl,kk) = sderv2(kl,kk) + d1trm*(dr2ds(kt)*dr2ikds(ks) + dr2ikds(kt)*dr2ds(ks))
                  enddo
                enddo
              endif
            enddo
!
!  HB coordination number and coordination number of j
!
            do n = 1,nnbr_cn(j)
              nk = nbrno_cn(n,j)
              rjk = rnbr(nk,j)
              k = nbrno(nk,j)
              r0jk = gfnff_rcov(nat(j)) + gfnff_rcov(nat(k))
              drjk = (rjk - r0jk)/r0jk
              dcnjdrjk = gfnff_kn_cn*exp(-(gfnff_kn_cn*drjk)**2)/(sqrt(pi)*r0jk*rjk)
!
              indk = 3*(k-1)
              kx = indk + 1
              ky = indk + 2
              kz = indk + 3
!
              d1trm = d2Edcnhdcnj*dcnjdrjk*dctrmdrhb*dlogcnjdcnj
!
              if (lstr) then
                call real1strterm(ndim,xnbr(nk,j),ynbr(nk,j),znbr(nk,j),0.0_dp,0.0_dp,0.0_dp, &
                                  dr2jkds,d2r2jkdx2,d2r2jkdsdx,d2r2jkds2,.false.)
              endif
!
              if (hbH.ge.j) then
                derv2(hbHx,jx) = derv2(hbHx,jx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
                derv2(hbHy,jx) = derv2(hbHy,jx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
                derv2(hbHz,jx) = derv2(hbHz,jx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
                derv2(hbHx,jy) = derv2(hbHx,jy) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
                derv2(hbHy,jy) = derv2(hbHy,jy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
                derv2(hbHz,jy) = derv2(hbHz,jy) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
                derv2(hbHx,jz) = derv2(hbHx,jz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
                derv2(hbHy,jz) = derv2(hbHy,jz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
                derv2(hbHz,jz) = derv2(hbHz,jz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
              else
                derv2(jx,hbHx) = derv2(jx,hbHx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
                derv2(jy,hbHx) = derv2(jy,hbHx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
                derv2(jz,hbHx) = derv2(jz,hbHx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
                derv2(jx,hbHy) = derv2(jx,hbHy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
                derv2(jy,hbHy) = derv2(jy,hbHy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
                derv2(jz,hbHy) = derv2(jz,hbHy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
                derv2(jx,hbHz) = derv2(jx,hbHz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
                derv2(jy,hbHz) = derv2(jy,hbHz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
                derv2(jz,hbHz) = derv2(jz,hbHz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
              endif
!
              if (hbH.ge.k) then
                derv2(hbHx,kx) = derv2(hbHx,kx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
                derv2(hbHy,kx) = derv2(hbHy,kx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
                derv2(hbHz,kx) = derv2(hbHz,kx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
                derv2(hbHx,ky) = derv2(hbHx,ky) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
                derv2(hbHy,ky) = derv2(hbHy,ky) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
                derv2(hbHz,ky) = derv2(hbHz,ky) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
                derv2(hbHx,kz) = derv2(hbHx,kz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
                derv2(hbHy,kz) = derv2(hbHy,kz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
                derv2(hbHz,kz) = derv2(hbHz,kz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
              else
                derv2(kx,hbHx) = derv2(kx,hbHx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
                derv2(ky,hbHx) = derv2(ky,hbHx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
                derv2(kz,hbHx) = derv2(kz,hbHx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
                derv2(kx,hbHy) = derv2(kx,hbHy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
                derv2(ky,hbHy) = derv2(ky,hbHy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
                derv2(kz,hbHy) = derv2(kz,hbHy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
                derv2(kx,hbHz) = derv2(kx,hbHz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
                derv2(ky,hbHz) = derv2(ky,hbHz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
                derv2(kz,hbHz) = derv2(kz,hbHz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
              endif
!
              if (hbB.ge.j) then
                derv2(hbBx,jx) = derv2(hbBx,jx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
                derv2(hbBy,jx) = derv2(hbBy,jx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
                derv2(hbBz,jx) = derv2(hbBz,jx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
                derv2(hbBx,jy) = derv2(hbBx,jy) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
                derv2(hbBy,jy) = derv2(hbBy,jy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
                derv2(hbBz,jy) = derv2(hbBz,jy) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
                derv2(hbBx,jz) = derv2(hbBx,jz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
                derv2(hbBy,jz) = derv2(hbBy,jz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
                derv2(hbBz,jz) = derv2(hbBz,jz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
              else
                derv2(jx,hbBx) = derv2(jx,hbBx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
                derv2(jy,hbBx) = derv2(jy,hbBx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
                derv2(jz,hbBx) = derv2(jz,hbBx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
                derv2(jx,hbBy) = derv2(jx,hbBy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
                derv2(jy,hbBy) = derv2(jy,hbBy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
                derv2(jz,hbBy) = derv2(jz,hbBy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
                derv2(jx,hbBz) = derv2(jx,hbBz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
                derv2(jy,hbBz) = derv2(jy,hbBz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
                derv2(jz,hbBz) = derv2(jz,hbBz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
              endif
!
              if (hbB.ge.k) then
                derv2(hbBx,kx) = derv2(hbBx,kx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
                derv2(hbBy,kx) = derv2(hbBy,kx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
                derv2(hbBz,kx) = derv2(hbBz,kx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
                derv2(hbBx,ky) = derv2(hbBx,ky) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
                derv2(hbBy,ky) = derv2(hbBy,ky) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
                derv2(hbBz,ky) = derv2(hbBz,ky) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
                derv2(hbBx,kz) = derv2(hbBx,kz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
                derv2(hbBy,kz) = derv2(hbBy,kz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
                derv2(hbBz,kz) = derv2(hbBz,kz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
              else
                derv2(kx,hbBx) = derv2(kx,hbBx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
                derv2(ky,hbBx) = derv2(ky,hbBx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
                derv2(kz,hbBx) = derv2(kz,hbBx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
                derv2(kx,hbBy) = derv2(kx,hbBy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
                derv2(ky,hbBy) = derv2(ky,hbBy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
                derv2(kz,hbBy) = derv2(kz,hbBy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
                derv2(kx,hbBz) = derv2(kx,hbBz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
                derv2(ky,hbBz) = derv2(ky,hbBz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
                derv2(kz,hbBz) = derv2(kz,hbBz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
              endif
!
!  Strain terms
!
              if (lstr) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(jx,kl) = derv3(jx,kl) - xnbr(nk,j)*d1trm*dr2ds(ks)
                  derv3(jy,kl) = derv3(jy,kl) - ynbr(nk,j)*d1trm*dr2ds(ks)
                  derv3(jz,kl) = derv3(jz,kl) - znbr(nk,j)*d1trm*dr2ds(ks)
                  derv3(kx,kl) = derv3(kx,kl) + xnbr(nk,j)*d1trm*dr2ds(ks)
                  derv3(ky,kl) = derv3(ky,kl) + ynbr(nk,j)*d1trm*dr2ds(ks)
                  derv3(kz,kl) = derv3(kz,kl) + znbr(nk,j)*d1trm*dr2ds(ks)
                enddo
!
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(hbHx,kl) = derv3(hbHx,kl) - cnhb_pairs%xpair(np)*d1trm*dr2jkds(ks)
                  derv3(hbHy,kl) = derv3(hbHy,kl) - cnhb_pairs%ypair(np)*d1trm*dr2jkds(ks)
                  derv3(hbHz,kl) = derv3(hbHz,kl) - cnhb_pairs%zpair(np)*d1trm*dr2jkds(ks)
                  derv3(hbBx,kl) = derv3(hbBx,kl) + cnhb_pairs%xpair(np)*d1trm*dr2jkds(ks)
                  derv3(hbBy,kl) = derv3(hbBy,kl) + cnhb_pairs%ypair(np)*d1trm*dr2jkds(ks)
                  derv3(hbBz,kl) = derv3(hbBz,kl) + cnhb_pairs%zpair(np)*d1trm*dr2jkds(ks)
                enddo
!
                do kk = 1,nstrains
                  ks = nstrptr(kk)
                  do kl = 1,nstrains
                    kt = nstrptr(kl)
                    sderv2(kl,kk) = sderv2(kl,kk) + d1trm*(dr2ds(kt)*dr2jkds(ks) + dr2jkds(kt)*dr2ds(ks))
                  enddo
                enddo
              endif
            enddo
          endif
        endif
      enddo
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_hb')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnd_hb(hbH,i,j,xij,yij,zij,dEdcnh,d2Edcnh2,d2Edcnhdr,d2Edcnhdcni,d2Edcnhdcnj,dcn)
!
!  Computes the second derivatives of the coordination number for GFNFF for hydrogen bonds for atom i-j pair
!  Distributed memory parallel version
!
!  On entry : 
!
!  hbH             = hydrogen atom whose coordination number is to contribute to the derivatives
!  i               = atom whose coordination number is to contribute to the derivatives via rij
!  j               = atom whose coordination number is to contribute to the derivatives via rij
!  xij             = x component of vector from i to j
!  yij             = y component of vector from i to j
!  zij             = z component of vector from i to j
!  dEdcnh          = first derivative of energy w.r.t. HB coordination number
!  d2Edcnh2        = second derivative of energy w.r.t. HB coordination number
!  d2Edcnhdr       = second derivative of energy w.r.t. HB coordination number and rij
!  d2Edcnhdcni     = second derivative of energy w.r.t. coordination numbers of hbH and i
!  d2Edcnhdcnj     = second derivative of energy w.r.t. coordination numbers of hbH and j
!  dcn             = derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  Derivatives have been computed and added to the appropriate arrays
!
!   3/22 Created from gfnff_drv2_dcn_hb
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
  use current
  use derivatives
  use g_constants,    only : pi
  use gulp_gfnff
  use iochannels
  use m_gfnff_nbr
  use m_gfnff_pairs
  use m_strain,       only : real1strterm
  use parallel
  use spatial
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: hbH
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: dEdcnh
  real(dp),    intent(in)                          :: d2Edcnh2
  real(dp),    intent(in)                          :: d2Edcnhdr
  real(dp),    intent(in)                          :: d2Edcnhdcni
  real(dp),    intent(in)                          :: d2Edcnhdcnj
  real(dp),    intent(in)                          :: dcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: hbA
  integer(i4)                                      :: hbB
  integer(i4)                                      :: hbBloc
  integer(i4)                                      :: hbBx
  integer(i4)                                      :: hbBy
  integer(i4)                                      :: hbBz
  integer(i4)                                      :: hbBxf
  integer(i4)                                      :: hbByf
  integer(i4)                                      :: hbBzf
  integer(i4)                                      :: hbHloc
  integer(i4)                                      :: hbHx
  integer(i4)                                      :: hbHy
  integer(i4)                                      :: hbHz
  integer(i4)                                      :: hbHxf
  integer(i4)                                      :: hbHyf
  integer(i4)                                      :: hbHzf
  integer(i4)                                      :: iloc
  integer(i4)                                      :: indh
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: ixf
  integer(i4)                                      :: iyf
  integer(i4)                                      :: izf
  integer(i4)                                      :: jloc
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: jxf
  integer(i4)                                      :: jyf
  integer(i4)                                      :: jzf
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: n
  integer(i4)                                      :: natB
  integer(i4)                                      :: natH
  integer(i4)                                      :: nh
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nk
  integer(i4)                                      :: np
  integer(i4)                                      :: np2
  logical                                          :: lABok
  logical                                          :: lABok2
  logical                                          :: lilocal
  logical                                          :: ljlocal
  logical                                          :: lklocal
  logical                                          :: lhbBlocal
  logical                                          :: lhbHlocal
  real(dp)                                         :: cncut2
  real(dp)                                         :: dctrmdrhb
  real(dp)                                         :: dctrmdrhb2
  real(dp)                                         :: d2ctrmdrhb2
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: dEdrhb
  real(dp)                                         :: d2Edrhb2
  real(dp)                                         :: dr2ds(6)
  real(dp)                                         :: dr2ds2(6)
  real(dp)                                         :: d2r2dx2(6)
  real(dp)                                         :: d2r2ds2(6,6)
  real(dp)                                         :: d2r2dsdx(6,3)
  real(dp)                                         :: dr2ijds(6)
  real(dp)                                         :: d2r2ijdx2(6)
  real(dp)                                         :: d2r2ijds2(6,6)
  real(dp)                                         :: d2r2ijdsdx(6,3)
  real(dp)                                         :: dr2ikds(6)
  real(dp)                                         :: d2r2ikdx2(6)
  real(dp)                                         :: d2r2ikds2(6,6)
  real(dp)                                         :: d2r2ikdsdx(6,3)
  real(dp)                                         :: dr2jkds(6)
  real(dp)                                         :: d2r2jkdx2(6)
  real(dp)                                         :: d2r2jkds2(6,6)
  real(dp)                                         :: d2r2jkdsdx(6,3)
  real(dp)                                         :: dtrm
  real(dp)                                         :: dtrm1a
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: dr
  real(dp)                                         :: dr2
  real(dp)                                         :: drik
  real(dp)                                         :: drjk
  real(dp),                          parameter     :: rcov_scale = 1.78_dp
  real(dp)                                         :: ofctij
  real(dp)                                         :: r0
  real(dp)                                         :: r0ik
  real(dp)                                         :: r0jk
  real(dp)                                         :: rik
  real(dp)                                         :: rjk
  real(dp)                                         :: rhb
  real(dp)                                         :: rhb2
  real(dp)                                         :: r2ab
  real(dp)                                         :: r2ab2
  real(dp)                                         :: xab
  real(dp)                                         :: yab
  real(dp)                                         :: zab
  real(dp)                                         :: xab2
  real(dp)                                         :: yab2
  real(dp)                                         :: zab2
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnd_hb')
#endif
!
!  Set up rij related terms if needed
!
  if (lstr) then
    call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ijds,d2r2ijdx2,d2r2ijdsdx,d2r2ijds2,.false.)
  endif
!
!  Is i or j or hbH local
!
  iloc = atom2local(i)
  jloc = atom2local(j)
  hbHloc = atom2local(hbH)
  lilocal = (iloc.ne.0)
  ljlocal = (jloc.ne.0)
  lhbHlocal = (hbHloc.ne.0)
!
!  Set up for second derivatives
!
  if (lilocal) then
    indi = 3*(iloc-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
  endif
!
  indi = 3*(i-1)
  ixf = indi + 1
  iyf = indi + 2
  izf = indi + 3
!
  if (ljlocal) then
    indj = 3*(jloc-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
  endif
!
  indj = 3*(j-1)
  jxf = indj + 1
  jyf = indj + 2
  jzf = indj + 3
!
  if (lhbHlocal) then
    indh = 3*(hbHloc-1)
    hbHx = indh + 1
    hbHy = indh + 2
    hbHz = indh + 3
  endif
!
  indh = 3*(hbH-1)
  hbHxf = indh + 1
  hbHyf = indh + 2
  hbHzf = indh + 3
!
  dlogcnidcni = dcn(i)
  dlogcnjdcnj = dcn(j)
!
!  Set factor for double counting
!
  if (i.eq.j) then
    ofctij = 1.0_dp
  else
    ofctij = 0.5_dp
  endif
!
  do ni = 1,nbond_hb_nr
! DEBUG - look for faster approach here!
    if (hbH.ne.nbond_hb_AH(2,ni)) cycle
    hbA = nbond_hb_AH(1,ni)
    nh  = nbond_hb_AH(3,ni)
    natH = nat(hbH)
    do nj = 1,nbond_hb_Bn(ni)
      hbB = nbond_hb_B(nj,ni)
      hbBloc = atom2local(hbB)
      lhbBlocal = (hbBloc.ne.0)
!
      if (lhbBlocal) then
        indh = 3*(hbBloc-1)
        hbBx = indh + 1
        hbBy = indh + 2
        hbBz = indh + 3
      endif
!
      indh = 3*(hbB-1)
      hbBxf = indh + 1
      hbByf = indh + 2
      hbBzf = indh + 3
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
!
          dr = (rhb - r0)/r0
          dtrm = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr)**2)/(sqrt(pi)*r0)
          dctrmdrhb = - dtrm/rhb
!
          dEdrhb = dEdcnh*dctrmdrhb
!----------------------
!  First derivatives  |
!----------------------
          if (lilocal) then
            xdrv(hbH) = xdrv(hbH) - dEdrhb*cnhb_pairs%xpair(np)*ofctij
            ydrv(hbH) = ydrv(hbH) - dEdrhb*cnhb_pairs%ypair(np)*ofctij
            zdrv(hbH) = zdrv(hbH) - dEdrhb*cnhb_pairs%zpair(np)*ofctij
!
            xdrv(hbB) = xdrv(hbB) + dEdrhb*cnhb_pairs%xpair(np)*ofctij
            ydrv(hbB) = ydrv(hbB) + dEdrhb*cnhb_pairs%ypair(np)*ofctij
            zdrv(hbB) = zdrv(hbB) + dEdrhb*cnhb_pairs%zpair(np)*ofctij
          endif
!
          call real1strterm(ndim,cnhb_pairs%xpair(np),cnhb_pairs%ypair(np),cnhb_pairs%zpair(np),0.0_dp,0.0_dp,0.0_dp, &
                            dr2ds,d2r2dx2,d2r2dsdx,d2r2ds2,.true.)
          if (lstr) then
            if (lilocal) then
              do kl = 1,nstrains
                ks = nstrptr(kl)
                rstrd(kl) = rstrd(kl) + dEdrhb*dr2ds(ks)*ofctij
              enddo
            endif
          endif
!
!  Check whether any atoms are local for the next section
!
          if (lhbBlocal.or.lhbHlocal) then
!-----------------------
!  Second derivatives  |
!-----------------------
!
!  HB coordination number only
!
            d2ctrmdrhb2 = dtrm*(1.0_dp + 2.0_dp*rhb*dr*(gfnff_kn_hb**2)/r0)/rhb**3
            d2Edrhb2 = d2Edcnh2*dctrmdrhb*dctrmdrhb + dEdcnh*d2ctrmdrhb2 ! (d2E/dr2)
!
            if (lhbHlocal) then
              derv2(hbBxf,hbHx) = derv2(hbBxf,hbHx) - d2Edrhb2*d2r2dx2(1)
              derv2(hbByf,hbHx) = derv2(hbByf,hbHx) - d2Edrhb2*d2r2dx2(6)
              derv2(hbBzf,hbHx) = derv2(hbBzf,hbHx) - d2Edrhb2*d2r2dx2(5)
              derv2(hbBxf,hbHy) = derv2(hbBxf,hbHy) - d2Edrhb2*d2r2dx2(6)
              derv2(hbByf,hbHy) = derv2(hbByf,hbHy) - d2Edrhb2*d2r2dx2(2)
              derv2(hbBzf,hbHy) = derv2(hbBzf,hbHy) - d2Edrhb2*d2r2dx2(4)
              derv2(hbBxf,hbHz) = derv2(hbBxf,hbHz) - d2Edrhb2*d2r2dx2(5)
              derv2(hbByf,hbHz) = derv2(hbByf,hbHz) - d2Edrhb2*d2r2dx2(4)
              derv2(hbBzf,hbHz) = derv2(hbBzf,hbHz) - d2Edrhb2*d2r2dx2(3)
              derv2(hbBxf,hbHx) = derv2(hbBxf,hbHx) - dEdrhb
              derv2(hbByf,hbHy) = derv2(hbByf,hbHy) - dEdrhb
              derv2(hbBzf,hbHz) = derv2(hbBzf,hbHz) - dEdrhb
            endif
!
            if (lhbBlocal) then
              derv2(hbHxf,hbBx) = derv2(hbHxf,hbBx) - d2Edrhb2*d2r2dx2(1)
              derv2(hbHyf,hbBx) = derv2(hbHyf,hbBx) - d2Edrhb2*d2r2dx2(6)
              derv2(hbHzf,hbBx) = derv2(hbHzf,hbBx) - d2Edrhb2*d2r2dx2(5)
              derv2(hbHxf,hbBy) = derv2(hbHxf,hbBy) - d2Edrhb2*d2r2dx2(6)
              derv2(hbHyf,hbBy) = derv2(hbHyf,hbBy) - d2Edrhb2*d2r2dx2(2)
              derv2(hbHzf,hbBy) = derv2(hbHzf,hbBy) - d2Edrhb2*d2r2dx2(4)
              derv2(hbHxf,hbBz) = derv2(hbHxf,hbBz) - d2Edrhb2*d2r2dx2(5)
              derv2(hbHyf,hbBz) = derv2(hbHyf,hbBz) - d2Edrhb2*d2r2dx2(4)
              derv2(hbHzf,hbBz) = derv2(hbHzf,hbBz) - d2Edrhb2*d2r2dx2(3)
              derv2(hbHxf,hbBx) = derv2(hbHxf,hbBx) - dEdrhb
              derv2(hbHyf,hbBy) = derv2(hbHyf,hbBy) - dEdrhb
              derv2(hbHzf,hbBz) = derv2(hbHzf,hbBz) - dEdrhb
            endif
!
            if (lstr) then
!
!  Mixed derivatives
!
              if (lhbBlocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(hbBx,kl) = derv3(hbBx,kl) + dEdrhb*d2r2dsdx(ks,1) + cnhb_pairs%xpair(np)*d2Edrhb2*dr2ds(ks)
                  derv3(hbBy,kl) = derv3(hbBy,kl) + dEdrhb*d2r2dsdx(ks,2) + cnhb_pairs%ypair(np)*d2Edrhb2*dr2ds(ks)
                  derv3(hbBz,kl) = derv3(hbBz,kl) + dEdrhb*d2r2dsdx(ks,3) + cnhb_pairs%zpair(np)*d2Edrhb2*dr2ds(ks)
                enddo
              endif
              if (lhbHlocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(hbHx,kl) = derv3(hbHx,kl) - dEdrhb*d2r2dsdx(ks,1) - cnhb_pairs%xpair(np)*d2Edrhb2*dr2ds(ks)
                  derv3(hbHy,kl) = derv3(hbHy,kl) - dEdrhb*d2r2dsdx(ks,2) - cnhb_pairs%ypair(np)*d2Edrhb2*dr2ds(ks)
                  derv3(hbHz,kl) = derv3(hbHz,kl) - dEdrhb*d2r2dsdx(ks,3) - cnhb_pairs%zpair(np)*d2Edrhb2*dr2ds(ks)
                enddo
              endif
!
!  Strain-strain
!
              if (lhbHlocal) then
                do kk = 1,nstrains
                  ks = nstrptr(kk)
                  do kl = 1,nstrains
                    kt = nstrptr(kl)
                    sderv2(kl,kk) = sderv2(kl,kk) + d2Edrhb2*dr2ds(kt)*dr2ds(ks) + dEdrhb*d2r2ds2(kt,ks)
                  enddo
                enddo
              endif
            endif
!
            if (np.gt.1) then
!
!  Loop over pairs of valid interactions
!
              do np2 = 1,np-1
!
!  Check that if B is an image of A then A=B case is excluded for the same image
!
                if (hbA.eq.hbB) then
                  xab2 = cnhb_pairs%xpair(np2) + xbnbr(nh,hbA)
                  yab2 = cnhb_pairs%ypair(np2) + ybnbr(nh,hbA)
                  zab2 = cnhb_pairs%zpair(np2) + zbnbr(nh,hbA)
                  r2ab2 = xab2**2 + yab2**2 + zab2**2
                  lABok2 = (r2ab2.gt.1.0d-2)
                else
                  lABok2 = .true.
                endif
!
                if (lABok2) then
                  rhb2 = sqrt(cnhb_pairs%r2pair(np2))
!
                  dr2 = (rhb2 - r0)/r0
                  dtrm1a = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr2)**2)/(sqrt(pi)*r0)
                  dctrmdrhb2 = - dtrm1a/rhb2
                  d2Edrhb2 = d2Edcnh2*dctrmdrhb*dctrmdrhb2
!
                  if (lhbHlocal) then
                    derv2(hbBxf,hbHx) = derv2(hbBxf,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbByf,hbHx) = derv2(hbByf,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBzf,hbHx) = derv2(hbBzf,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBxf,hbHy) = derv2(hbBxf,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbByf,hbHy) = derv2(hbByf,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBzf,hbHy) = derv2(hbBzf,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBxf,hbHz) = derv2(hbBxf,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbByf,hbHz) = derv2(hbByf,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBzf,hbHz) = derv2(hbBzf,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  endif
!
                  if (lhbBlocal) then
                    derv2(hbHxf,hbBx) = derv2(hbHxf,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHyf,hbBx) = derv2(hbHyf,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHzf,hbBx) = derv2(hbHzf,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHxf,hbBy) = derv2(hbHxf,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHyf,hbBy) = derv2(hbHyf,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHzf,hbBy) = derv2(hbHzf,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHxf,hbBz) = derv2(hbHxf,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHyf,hbBz) = derv2(hbHyf,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHzf,hbBz) = derv2(hbHzf,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  endif
!
                  if (lstr) then
                    call real1strterm(ndim,cnhb_pairs%xpair(np2),cnhb_pairs%ypair(np2),cnhb_pairs%zpair(np2), &
                                      0.0_dp,0.0_dp,0.0_dp,dr2ds2,d2r2dx2,d2r2dsdx,d2r2ds2,.false.)
!
!  Mixed derivatives
!
                    if (lhbBlocal) then
                      do kl = 1,nstrains
                        ks = nstrptr(kl)
                        derv3(hbBx,kl) = derv3(hbBx,kl) + cnhb_pairs%xpair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                        + cnhb_pairs%xpair(np2)*d2Edrhb2*dr2ds(ks)
                        derv3(hbBy,kl) = derv3(hbBy,kl) + cnhb_pairs%ypair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                        + cnhb_pairs%ypair(np2)*d2Edrhb2*dr2ds(ks)
                        derv3(hbBz,kl) = derv3(hbBz,kl) + cnhb_pairs%zpair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                        + cnhb_pairs%zpair(np2)*d2Edrhb2*dr2ds(ks)
                      enddo
                    endif
                    if (lhbHlocal) then
                      do kl = 1,nstrains
                        ks = nstrptr(kl)
                        derv3(hbHx,kl) = derv3(hbHx,kl) - cnhb_pairs%xpair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                        - cnhb_pairs%xpair(np2)*d2Edrhb2*dr2ds(ks)
                        derv3(hbHy,kl) = derv3(hbHy,kl) - cnhb_pairs%ypair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                        - cnhb_pairs%ypair(np2)*d2Edrhb2*dr2ds(ks)
                        derv3(hbHz,kl) = derv3(hbHz,kl) - cnhb_pairs%zpair(np)*d2Edrhb2*dr2ds2(ks)  &
                                                        - cnhb_pairs%zpair(np2)*d2Edrhb2*dr2ds(ks)
                      enddo
                    endif
!
!  Strain-strain
!
                    if (lhbHlocal) then
                      do kk = 1,nstrains
                        ks = nstrptr(kk)
                        do kl = 1,nstrains
                          kt = nstrptr(kl)
                          sderv2(kl,kk) = sderv2(kl,kk) + d2Edrhb2*dr2ds2(kt)*dr2ds(ks)
                        enddo
                      enddo
                    endif
                  endif
                endif
              enddo
            endif
          endif
!
!  HB coordination number and rij
!
          d2trm = d2Edcnhdr*dctrmdrhb
!
          if (lilocal) then
            derv2(hbHxf,ix) = derv2(hbHxf,ix) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbHyf,ix) = derv2(hbHyf,ix) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbHzf,ix) = derv2(hbHzf,ix) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbHxf,iy) = derv2(hbHxf,iy) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbHyf,iy) = derv2(hbHyf,iy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbHzf,iy) = derv2(hbHzf,iy) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbHxf,iz) = derv2(hbHxf,iz) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbHyf,iz) = derv2(hbHyf,iz) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbHzf,iz) = derv2(hbHzf,iz) + d2trm*cnhb_pairs%zpair(np)*zij
!
            derv2(hbBxf,ix) = derv2(hbBxf,ix) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbByf,ix) = derv2(hbByf,ix) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbBzf,ix) = derv2(hbBzf,ix) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbBxf,iy) = derv2(hbBxf,iy) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbByf,iy) = derv2(hbByf,iy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbBzf,iy) = derv2(hbBzf,iy) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbBxf,iz) = derv2(hbBxf,iz) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbByf,iz) = derv2(hbByf,iz) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbBzf,iz) = derv2(hbBzf,iz) - d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
          if (lhbHlocal) then
            derv2(ixf,hbHx) = derv2(ixf,hbHx) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(iyf,hbHx) = derv2(iyf,hbHx) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(izf,hbHx) = derv2(izf,hbHx) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(ixf,hbHy) = derv2(ixf,hbHy) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(iyf,hbHy) = derv2(iyf,hbHy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(izf,hbHy) = derv2(izf,hbHy) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(ixf,hbHz) = derv2(ixf,hbHz) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(iyf,hbHz) = derv2(iyf,hbHz) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(izf,hbHz) = derv2(izf,hbHz) + d2trm*cnhb_pairs%zpair(np)*zij
!
            derv2(jxf,hbHx) = derv2(jxf,hbHx) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(jyf,hbHx) = derv2(jyf,hbHx) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(jzf,hbHx) = derv2(jzf,hbHx) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(jxf,hbHy) = derv2(jxf,hbHy) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(jyf,hbHy) = derv2(jyf,hbHy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(jzf,hbHy) = derv2(jzf,hbHy) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(jxf,hbHz) = derv2(jxf,hbHz) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(jyf,hbHz) = derv2(jyf,hbHz) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(jzf,hbHz) = derv2(jzf,hbHz) - d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
          if (lhbBlocal) then
            derv2(ixf,hbBx) = derv2(ixf,hbBx) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(iyf,hbBx) = derv2(iyf,hbBx) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(izf,hbBx) = derv2(izf,hbBx) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(ixf,hbBy) = derv2(ixf,hbBy) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(iyf,hbBy) = derv2(iyf,hbBy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(izf,hbBy) = derv2(izf,hbBy) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(ixf,hbBz) = derv2(ixf,hbBz) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(iyf,hbBz) = derv2(iyf,hbBz) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(izf,hbBz) = derv2(izf,hbBz) - d2trm*cnhb_pairs%zpair(np)*zij
!
            derv2(jxf,hbBx) = derv2(jxf,hbBx) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(jyf,hbBx) = derv2(jyf,hbBx) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(jzf,hbBx) = derv2(jzf,hbBx) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(jxf,hbBy) = derv2(jxf,hbBy) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(jyf,hbBy) = derv2(jyf,hbBy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(jzf,hbBy) = derv2(jzf,hbBy) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(jxf,hbBz) = derv2(jxf,hbBz) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(jyf,hbBz) = derv2(jyf,hbBz) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(jzf,hbBz) = derv2(jzf,hbBz) + d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
          if (ljlocal) then
            derv2(hbHxf,jx) = derv2(hbHxf,jx) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbHyf,jx) = derv2(hbHyf,jx) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbHzf,jx) = derv2(hbHzf,jx) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbHxf,jy) = derv2(hbHxf,jy) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbHyf,jy) = derv2(hbHyf,jy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbHzf,jy) = derv2(hbHzf,jy) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbHxf,jz) = derv2(hbHxf,jz) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbHyf,jz) = derv2(hbHyf,jz) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbHzf,jz) = derv2(hbHzf,jz) - d2trm*cnhb_pairs%zpair(np)*zij
!
            derv2(hbBxf,jx) = derv2(hbBxf,jx) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbByf,jx) = derv2(hbByf,jx) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbBzf,jx) = derv2(hbBzf,jx) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbBxf,jy) = derv2(hbBxf,jy) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbByf,jy) = derv2(hbByf,jy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbBzf,jy) = derv2(hbBzf,jy) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbBxf,jz) = derv2(hbBxf,jz) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbByf,jz) = derv2(hbByf,jz) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbBzf,jz) = derv2(hbBzf,jz) + d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
!  Strain terms
!
          if (lstr) then
            if (lilocal) then
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(ix,kl) = derv3(ix,kl) - xij*d2trm*dr2ds(ks)
                derv3(iy,kl) = derv3(iy,kl) - yij*d2trm*dr2ds(ks)
                derv3(iz,kl) = derv3(iz,kl) - zij*d2trm*dr2ds(ks)
              enddo
            endif
!
            if (ljlocal) then
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(jx,kl) = derv3(jx,kl) + xij*d2trm*dr2ds(ks)
                derv3(jy,kl) = derv3(jy,kl) + yij*d2trm*dr2ds(ks)
                derv3(jz,kl) = derv3(jz,kl) + zij*d2trm*dr2ds(ks)
              enddo
            endif
!
            if (lhbBlocal) then
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(hbBx,kl) = derv3(hbBx,kl) + cnhb_pairs%xpair(np)*d2trm*dr2ijds(ks)
                derv3(hbBy,kl) = derv3(hbBy,kl) + cnhb_pairs%ypair(np)*d2trm*dr2ijds(ks)
                derv3(hbBz,kl) = derv3(hbBz,kl) + cnhb_pairs%zpair(np)*d2trm*dr2ijds(ks)
              enddo
            endif
!
            if (lhbHlocal) then
              do kl = 1,nstrains
                ks = nstrptr(kl)
                derv3(hbHx,kl) = derv3(hbHx,kl) - cnhb_pairs%xpair(np)*d2trm*dr2ijds(ks)
                derv3(hbHy,kl) = derv3(hbHy,kl) - cnhb_pairs%ypair(np)*d2trm*dr2ijds(ks)
                derv3(hbHz,kl) = derv3(hbHz,kl) - cnhb_pairs%zpair(np)*d2trm*dr2ijds(ks)
              enddo
            endif
!
            if (lilocal) then
              do kk = 1,nstrains
                ks = nstrptr(kk)
                do kl = 1,nstrains
                  kt = nstrptr(kl)
                  sderv2(kl,kk) = sderv2(kl,kk) + d2trm*(dr2ds(kt)*dr2ijds(ks) + dr2ijds(kt)*dr2ds(ks))
                enddo
              enddo
            endif
          endif
!
!  HB coordination number and coordination number of i
!
          do n = 1,nnbr_cn(i)
            nk = nbrno_cn(n,i)
            rik = rnbr(nk,i)
            k = nbrno(nk,i)
!
            kloc = atom2local(k)
            lklocal = (kloc.ne.0)
!
!  Check whether any atoms are local
!
            if (.not.lilocal.and..not.lklocal.and..not.lhbBlocal.and..not.lhbHlocal) cycle
!
            r0ik = gfnff_rcov(nat(i)) + gfnff_rcov(nat(k))
            drik = (rik - r0ik)/r0ik
            dcnidrik = gfnff_kn_cn*exp(-(gfnff_kn_cn*drik)**2)/(sqrt(pi)*r0ik*rik)
!
            if (lklocal) then
              indk = 3*(kloc-1)
              kx = indk + 1
              ky = indk + 2
              kz = indk + 3
            endif
!
            indk = 3*(k-1)
            kxf = indk + 1
            kyf = indk + 2
            kzf = indk + 3
!
            d1trm = d2Edcnhdcni*dcnidrik*dctrmdrhb*dlogcnidcni
!
            if (lstr) then
              call real1strterm(ndim,xnbr(nk,i),ynbr(nk,i),znbr(nk,i),0.0_dp,0.0_dp,0.0_dp, &
                                dr2ikds,d2r2ikdx2,d2r2ikdsdx,d2r2ikds2,.false.)
            endif
!
            if (lilocal) then
              derv2(hbHxf,ix) = derv2(hbHxf,ix) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbHyf,ix) = derv2(hbHyf,ix) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbHzf,ix) = derv2(hbHzf,ix) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbHxf,iy) = derv2(hbHxf,iy) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbHyf,iy) = derv2(hbHyf,iy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbHzf,iy) = derv2(hbHzf,iy) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbHxf,iz) = derv2(hbHxf,iz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbHyf,iz) = derv2(hbHyf,iz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbHzf,iz) = derv2(hbHzf,iz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
!
              derv2(hbBxf,ix) = derv2(hbBxf,ix) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbByf,ix) = derv2(hbByf,ix) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbBzf,ix) = derv2(hbBzf,ix) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbBxf,iy) = derv2(hbBxf,iy) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbByf,iy) = derv2(hbByf,iy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbBzf,iy) = derv2(hbBzf,iy) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbBxf,iz) = derv2(hbBxf,iz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbByf,iz) = derv2(hbByf,iz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbBzf,iz) = derv2(hbBzf,iz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
            if (lhbHlocal) then
              derv2(ixf,hbHx) = derv2(ixf,hbHx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(iyf,hbHx) = derv2(iyf,hbHx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(izf,hbHx) = derv2(izf,hbHx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(ixf,hbHy) = derv2(ixf,hbHy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(iyf,hbHy) = derv2(iyf,hbHy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(izf,hbHy) = derv2(izf,hbHy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(ixf,hbHz) = derv2(ixf,hbHz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(iyf,hbHz) = derv2(iyf,hbHz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(izf,hbHz) = derv2(izf,hbHz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
!
              derv2(kxf,hbHx) = derv2(kxf,hbHx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(kyf,hbHx) = derv2(kyf,hbHx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(kzf,hbHx) = derv2(kzf,hbHx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(kxf,hbHy) = derv2(kxf,hbHy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(kyf,hbHy) = derv2(kyf,hbHy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(kzf,hbHy) = derv2(kzf,hbHy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(kxf,hbHz) = derv2(kxf,hbHz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(kyf,hbHz) = derv2(kyf,hbHz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(kzf,hbHz) = derv2(kzf,hbHz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
            if (lhbBlocal) then
              derv2(ixf,hbBx) = derv2(ixf,hbBx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(iyf,hbBx) = derv2(iyf,hbBx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(izf,hbBx) = derv2(izf,hbBx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(ixf,hbBy) = derv2(ixf,hbBy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(iyf,hbBy) = derv2(iyf,hbBy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(izf,hbBy) = derv2(izf,hbBy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(ixf,hbBz) = derv2(ixf,hbBz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(iyf,hbBz) = derv2(iyf,hbBz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(izf,hbBz) = derv2(izf,hbBz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
!
              derv2(kxf,hbBx) = derv2(kxf,hbBx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(kyf,hbBx) = derv2(kyf,hbBx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(kzf,hbBx) = derv2(kzf,hbBx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(kxf,hbBy) = derv2(kxf,hbBy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(kyf,hbBy) = derv2(kyf,hbBy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(kzf,hbBy) = derv2(kzf,hbBy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(kxf,hbBz) = derv2(kxf,hbBz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(kyf,hbBz) = derv2(kyf,hbBz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(kzf,hbBz) = derv2(kzf,hbBz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
            if (lklocal) then
              derv2(hbHxf,kx) = derv2(hbHxf,kx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbHyf,kx) = derv2(hbHyf,kx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbHzf,kx) = derv2(hbHzf,kx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbHxf,ky) = derv2(hbHxf,ky) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbHyf,ky) = derv2(hbHyf,ky) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbHzf,ky) = derv2(hbHzf,ky) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbHxf,kz) = derv2(hbHxf,kz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbHyf,kz) = derv2(hbHyf,kz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbHzf,kz) = derv2(hbHzf,kz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
!
              derv2(hbBxf,kx) = derv2(hbBxf,kx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbByf,kx) = derv2(hbByf,kx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbBzf,kx) = derv2(hbBzf,kx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbBxf,ky) = derv2(hbBxf,ky) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbByf,ky) = derv2(hbByf,ky) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbBzf,ky) = derv2(hbBzf,ky) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbBxf,kz) = derv2(hbBxf,kz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbByf,kz) = derv2(hbByf,kz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbBzf,kz) = derv2(hbBzf,kz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
!  Strain terms
!
            if (lstr) then
              if (lilocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(ix,kl) = derv3(ix,kl) - xnbr(nk,i)*d1trm*dr2ds(ks)
                  derv3(iy,kl) = derv3(iy,kl) - ynbr(nk,i)*d1trm*dr2ds(ks)
                  derv3(iz,kl) = derv3(iz,kl) - znbr(nk,i)*d1trm*dr2ds(ks)
                enddo
              endif
!
              if (lklocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(kx,kl) = derv3(kx,kl) + xnbr(nk,i)*d1trm*dr2ds(ks)
                  derv3(ky,kl) = derv3(ky,kl) + ynbr(nk,i)*d1trm*dr2ds(ks)
                  derv3(kz,kl) = derv3(kz,kl) + znbr(nk,i)*d1trm*dr2ds(ks)
                enddo
              endif
!
              if (lhbBlocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(hbBx,kl) = derv3(hbBx,kl) + cnhb_pairs%xpair(np)*d1trm*dr2ikds(ks)
                  derv3(hbBy,kl) = derv3(hbBy,kl) + cnhb_pairs%ypair(np)*d1trm*dr2ikds(ks)
                  derv3(hbBz,kl) = derv3(hbBz,kl) + cnhb_pairs%zpair(np)*d1trm*dr2ikds(ks)
                enddo
              endif
!
              if (lhbHlocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(hbHx,kl) = derv3(hbHx,kl) - cnhb_pairs%xpair(np)*d1trm*dr2ikds(ks)
                  derv3(hbHy,kl) = derv3(hbHy,kl) - cnhb_pairs%ypair(np)*d1trm*dr2ikds(ks)
                  derv3(hbHz,kl) = derv3(hbHz,kl) - cnhb_pairs%zpair(np)*d1trm*dr2ikds(ks)
                enddo
              endif
!
              if (lilocal) then
                do kk = 1,nstrains
                  ks = nstrptr(kk)
                  do kl = 1,nstrains
                    kt = nstrptr(kl)
                    sderv2(kl,kk) = sderv2(kl,kk) + d1trm*(dr2ds(kt)*dr2ikds(ks) + dr2ikds(kt)*dr2ds(ks))
                  enddo
                enddo
              endif
            endif
          enddo
!
!  HB coordination number and coordination number of j
!
          do n = 1,nnbr_cn(j)
            nk = nbrno_cn(n,j)
            rjk = rnbr(nk,j)
            k = nbrno(nk,j)
!
            kloc = atom2local(k)
            lklocal = (kloc.ne.0)
!
!  Check whether any atoms are local
!
            if (.not.ljlocal.and..not.lklocal.and..not.lhbBlocal.and..not.lhbHlocal) cycle
!
            r0jk = gfnff_rcov(nat(j)) + gfnff_rcov(nat(k))
            drjk = (rjk - r0jk)/r0jk
            dcnjdrjk = gfnff_kn_cn*exp(-(gfnff_kn_cn*drjk)**2)/(sqrt(pi)*r0jk*rjk)
!
            if (lklocal) then
              indk = 3*(kloc-1)
              kx = indk + 1
              ky = indk + 2
              kz = indk + 3
            endif
!
            indk = 3*(k-1)
            kxf = indk + 1
            kyf = indk + 2
            kzf = indk + 3
!
            d1trm = d2Edcnhdcnj*dcnjdrjk*dctrmdrhb*dlogcnjdcnj
!
            if (lstr) then
              call real1strterm(ndim,xnbr(nk,j),ynbr(nk,j),znbr(nk,j),0.0_dp,0.0_dp,0.0_dp, &
                                dr2jkds,d2r2jkdx2,d2r2jkdsdx,d2r2jkds2,.false.)
            endif
!
            if (ljlocal) then
              derv2(hbHxf,jx) = derv2(hbHxf,jx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbHyf,jx) = derv2(hbHyf,jx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbHzf,jx) = derv2(hbHzf,jx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbHxf,jy) = derv2(hbHxf,jy) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbHyf,jy) = derv2(hbHyf,jy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbHzf,jy) = derv2(hbHzf,jy) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbHxf,jz) = derv2(hbHxf,jz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbHyf,jz) = derv2(hbHyf,jz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbHzf,jz) = derv2(hbHzf,jz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
!
              derv2(hbBxf,jx) = derv2(hbBxf,jx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbByf,jx) = derv2(hbByf,jx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbBzf,jx) = derv2(hbBzf,jx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbBxf,jy) = derv2(hbBxf,jy) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbByf,jy) = derv2(hbByf,jy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbBzf,jy) = derv2(hbBzf,jy) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbBxf,jz) = derv2(hbBxf,jz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbByf,jz) = derv2(hbByf,jz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbBzf,jz) = derv2(hbBzf,jz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
            if (lhbHlocal) then
              derv2(jxf,hbHx) = derv2(jxf,hbHx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(jyf,hbHx) = derv2(jyf,hbHx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(jzf,hbHx) = derv2(jzf,hbHx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(jxf,hbHy) = derv2(jxf,hbHy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(jyf,hbHy) = derv2(jyf,hbHy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(jzf,hbHy) = derv2(jzf,hbHy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(jxf,hbHz) = derv2(jxf,hbHz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(jyf,hbHz) = derv2(jyf,hbHz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(jzf,hbHz) = derv2(jzf,hbHz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
!
              derv2(kxf,hbHx) = derv2(kxf,hbHx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(kyf,hbHx) = derv2(kyf,hbHx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(kzf,hbHx) = derv2(kzf,hbHx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(kxf,hbHy) = derv2(kxf,hbHy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(kyf,hbHy) = derv2(kyf,hbHy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(kzf,hbHy) = derv2(kzf,hbHy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(kxf,hbHz) = derv2(kxf,hbHz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(kyf,hbHz) = derv2(kyf,hbHz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(kzf,hbHz) = derv2(kzf,hbHz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
            if (lhbBlocal) then
              derv2(jxf,hbBx) = derv2(jxf,hbBx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(jyf,hbBx) = derv2(jyf,hbBx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(jzf,hbBx) = derv2(jzf,hbBx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(jxf,hbBy) = derv2(jxf,hbBy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(jyf,hbBy) = derv2(jyf,hbBy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(jzf,hbBy) = derv2(jzf,hbBy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(jxf,hbBz) = derv2(jxf,hbBz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(jyf,hbBz) = derv2(jyf,hbBz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(jzf,hbBz) = derv2(jzf,hbBz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
!
              derv2(kxf,hbBx) = derv2(kxf,hbBx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(kyf,hbBx) = derv2(kyf,hbBx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(kzf,hbBx) = derv2(kzf,hbBx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(kxf,hbBy) = derv2(kxf,hbBy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(kyf,hbBy) = derv2(kyf,hbBy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(kzf,hbBy) = derv2(kzf,hbBy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(kxf,hbBz) = derv2(kxf,hbBz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(kyf,hbBz) = derv2(kyf,hbBz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(kzf,hbBz) = derv2(kzf,hbBz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
            if (lklocal) then
              derv2(hbHxf,kx) = derv2(hbHxf,kx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbHyf,kx) = derv2(hbHyf,kx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbHzf,kx) = derv2(hbHzf,kx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbHxf,ky) = derv2(hbHxf,ky) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbHyf,ky) = derv2(hbHyf,ky) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbHzf,ky) = derv2(hbHzf,ky) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbHxf,kz) = derv2(hbHxf,kz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbHyf,kz) = derv2(hbHyf,kz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbHzf,kz) = derv2(hbHzf,kz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
!
              derv2(hbBxf,kx) = derv2(hbBxf,kx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbByf,kx) = derv2(hbByf,kx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbBzf,kx) = derv2(hbBzf,kx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbBxf,ky) = derv2(hbBxf,ky) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbByf,ky) = derv2(hbByf,ky) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbBzf,ky) = derv2(hbBzf,ky) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbBxf,kz) = derv2(hbBxf,kz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbByf,kz) = derv2(hbByf,kz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbBzf,kz) = derv2(hbBzf,kz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
!  Strain terms
!
            if (lstr) then
              if (ljlocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(jx,kl) = derv3(jx,kl) - xnbr(nk,j)*d1trm*dr2ds(ks)
                  derv3(jy,kl) = derv3(jy,kl) - ynbr(nk,j)*d1trm*dr2ds(ks)
                  derv3(jz,kl) = derv3(jz,kl) - znbr(nk,j)*d1trm*dr2ds(ks)
                enddo
              endif
!
              if (lklocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(kx,kl) = derv3(kx,kl) + xnbr(nk,j)*d1trm*dr2ds(ks)
                  derv3(ky,kl) = derv3(ky,kl) + ynbr(nk,j)*d1trm*dr2ds(ks)
                  derv3(kz,kl) = derv3(kz,kl) + znbr(nk,j)*d1trm*dr2ds(ks)
                enddo
              endif
!
              if (lhbBlocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(hbBx,kl) = derv3(hbBx,kl) + cnhb_pairs%xpair(np)*d1trm*dr2jkds(ks)
                  derv3(hbBy,kl) = derv3(hbBy,kl) + cnhb_pairs%ypair(np)*d1trm*dr2jkds(ks)
                  derv3(hbBz,kl) = derv3(hbBz,kl) + cnhb_pairs%zpair(np)*d1trm*dr2jkds(ks)
                enddo
              endif
!
              if (lhbHlocal) then
                do kl = 1,nstrains
                  ks = nstrptr(kl)
                  derv3(hbHx,kl) = derv3(hbHx,kl) - cnhb_pairs%xpair(np)*d1trm*dr2jkds(ks)
                  derv3(hbHy,kl) = derv3(hbHy,kl) - cnhb_pairs%ypair(np)*d1trm*dr2jkds(ks)
                  derv3(hbHz,kl) = derv3(hbHz,kl) - cnhb_pairs%zpair(np)*d1trm*dr2jkds(ks)
                enddo
              endif
!
              if (ljlocal) then
                do kk = 1,nstrains
                  ks = nstrptr(kk)
                  do kl = 1,nstrains
                    kt = nstrptr(kl)
                    sderv2(kl,kk) = sderv2(kl,kk) + d1trm*(dr2ds(kt)*dr2jkds(ks) + dr2jkds(kt)*dr2ds(ks))
                  enddo
                enddo
              endif
            endif
          enddo
        endif
      enddo
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnd_hb')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcn_hbp(hbH,i,j,xij,yij,zij,dEdcnh,d2Edcnh2,d2Edcnhdr,d2Edcnhdcni,d2Edcnhdcnj,dcn)
!
!  Computes the derivatives of the coordination number for GFNFF for hydrogen bonds for atom i-j pair
!  Phonon version.
!
!  On entry : 
!
!  hbH             = hydrogen atom whose coordination number is to contribute to the derivatives
!  i               = atom whose coordination number is to contribute to the derivatives via rij
!  j               = atom whose coordination number is to contribute to the derivatives via rij
!  xij             = x component of vector from i to j
!  yij             = y component of vector from i to j
!  zij             = z component of vector from i to j
!  dEdcnh          = first derivative of energy w.r.t. HB coordination number
!  d2Edcnh2        = second derivative of energy w.r.t. HB coordination number
!  d2Edcnhdr       = second derivative of energy w.r.t. HB coordination number and rij
!  d2Edcnhdcni     = second derivative of energy w.r.t. coordination numbers of hbH and i
!  d2Edcnhdcnj     = second derivative of energy w.r.t. coordination numbers of hbH and j
!  dcn             = derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  Derivatives have been computed and added to the appropriate arrays
!
!   3/21 Created from gfnff_drv2_dcn_hb
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
  use current
  use derivatives
  use g_constants,    only : pi
  use gulp_gfnff
  use iochannels
  use m_gfnff_nbr
  use m_gfnff_pairs
  use spatial
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: hbH
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: dEdcnh
  real(dp),    intent(in)                          :: d2Edcnh2
  real(dp),    intent(in)                          :: d2Edcnhdr
  real(dp),    intent(in)                          :: d2Edcnhdcni
  real(dp),    intent(in)                          :: d2Edcnhdcnj
  real(dp),    intent(in)                          :: dcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: hbA
  integer(i4)                                      :: hbB
  integer(i4)                                      :: hbBx
  integer(i4)                                      :: hbBy
  integer(i4)                                      :: hbBz
  integer(i4)                                      :: hbHx
  integer(i4)                                      :: hbHy
  integer(i4)                                      :: hbHz
  integer(i4)                                      :: indh
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: k
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: n
  integer(i4)                                      :: natB
  integer(i4)                                      :: natH
  integer(i4)                                      :: nh
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nk
  integer(i4)                                      :: np
  integer(i4)                                      :: np2
  logical                                          :: lABok
  logical                                          :: lABok2
  real(dp)                                         :: cncut2
  real(dp)                                         :: dctrmdrhb
  real(dp)                                         :: dctrmdrhb2
  real(dp)                                         :: d2ctrmdrhb2
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: dEdrhb
  real(dp)                                         :: d2Edrhb2
  real(dp)                                         :: d2r2dx2(6)
  real(dp)                                         :: dtrm
  real(dp)                                         :: dtrm1a
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: dr
  real(dp)                                         :: dr2
  real(dp)                                         :: drik
  real(dp)                                         :: drjk
  real(dp),                          parameter     :: rcov_scale = 1.78_dp
  real(dp)                                         :: r0
  real(dp)                                         :: r0ik
  real(dp)                                         :: r0jk
  real(dp)                                         :: rik
  real(dp)                                         :: rjk
  real(dp)                                         :: rhb
  real(dp)                                         :: rhb2
  real(dp)                                         :: r2ab
  real(dp)                                         :: r2ab2
  real(dp)                                         :: xab
  real(dp)                                         :: yab
  real(dp)                                         :: zab
  real(dp)                                         :: xab2
  real(dp)                                         :: yab2
  real(dp)                                         :: zab2
#ifdef TRACE
  call trace_in('gfnff_drv2_dcn_hbp')
#endif
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ix = indi + 1
  iy = indi + 2
  iz = indi + 3
!
  indj = 3*(j-1)
  jx = indj + 1
  jy = indj + 2
  jz = indj + 3
!
  indh = 3*(hbH-1)
  hbHx = indh + 1
  hbHy = indh + 2
  hbHz = indh + 3
!
  dlogcnidcni = dcn(i)
  dlogcnjdcnj = dcn(j)
!
  do ni = 1,nbond_hb_nr
! DEBUG - look for faster approach here!
    if (hbH.ne.nbond_hb_AH(2,ni)) cycle
    hbA = nbond_hb_AH(1,ni)
    nh  = nbond_hb_AH(3,ni)
    natH = nat(hbH)
    do nj = 1,nbond_hb_Bn(ni)
      hbB = nbond_hb_B(nj,ni)
      indh = 3*(hbB-1)
      hbBx = indh + 1
      hbBy = indh + 2
      hbBz = indh + 3
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
!
          dr = (rhb - r0)/r0
          dtrm = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr)**2)/(sqrt(pi)*r0)
          dctrmdrhb = - dtrm/rhb
!
          dEdrhb = dEdcnh*dctrmdrhb
!
          d2r2dx2(1) = cnhb_pairs%xpair(np)*cnhb_pairs%xpair(np)
          d2r2dx2(2) = cnhb_pairs%ypair(np)*cnhb_pairs%ypair(np)
          d2r2dx2(3) = cnhb_pairs%zpair(np)*cnhb_pairs%zpair(np)
          d2r2dx2(4) = cnhb_pairs%ypair(np)*cnhb_pairs%zpair(np)
          d2r2dx2(5) = cnhb_pairs%xpair(np)*cnhb_pairs%zpair(np)
          d2r2dx2(6) = cnhb_pairs%xpair(np)*cnhb_pairs%ypair(np)
!-----------------------
!  Second derivatives  |
!-----------------------
!
!  HB coordination number only
!
          d2ctrmdrhb2 = dtrm*(1.0_dp + 2.0_dp*rhb*dr*(gfnff_kn_hb**2)/r0)/rhb**3
          d2Edrhb2 = d2Edcnh2*dctrmdrhb*dctrmdrhb + dEdcnh*d2ctrmdrhb2 ! (d2E/dr2)
!
          if (hbB.gt.hbH) then
            derv2(hbBx,hbHx) = derv2(hbBx,hbHx) - d2Edrhb2*d2r2dx2(1)
            derv2(hbBy,hbHx) = derv2(hbBy,hbHx) - d2Edrhb2*d2r2dx2(6)
            derv2(hbBz,hbHx) = derv2(hbBz,hbHx) - d2Edrhb2*d2r2dx2(5)
            derv2(hbBx,hbHy) = derv2(hbBx,hbHy) - d2Edrhb2*d2r2dx2(6)
            derv2(hbBy,hbHy) = derv2(hbBy,hbHy) - d2Edrhb2*d2r2dx2(2)
            derv2(hbBz,hbHy) = derv2(hbBz,hbHy) - d2Edrhb2*d2r2dx2(4)
            derv2(hbBx,hbHz) = derv2(hbBx,hbHz) - d2Edrhb2*d2r2dx2(5)
            derv2(hbBy,hbHz) = derv2(hbBy,hbHz) - d2Edrhb2*d2r2dx2(4)
            derv2(hbBz,hbHz) = derv2(hbBz,hbHz) - d2Edrhb2*d2r2dx2(3)
            derv2(hbBx,hbHx) = derv2(hbBx,hbHx) - dEdrhb
            derv2(hbBy,hbHy) = derv2(hbBy,hbHy) - dEdrhb
            derv2(hbBz,hbHz) = derv2(hbBz,hbHz) - dEdrhb
          else
            derv2(hbHx,hbBx) = derv2(hbHx,hbBx) - d2Edrhb2*d2r2dx2(1)
            derv2(hbHy,hbBx) = derv2(hbHy,hbBx) - d2Edrhb2*d2r2dx2(6)
            derv2(hbHz,hbBx) = derv2(hbHz,hbBx) - d2Edrhb2*d2r2dx2(5)
            derv2(hbHx,hbBy) = derv2(hbHx,hbBy) - d2Edrhb2*d2r2dx2(6)
            derv2(hbHy,hbBy) = derv2(hbHy,hbBy) - d2Edrhb2*d2r2dx2(2)
            derv2(hbHz,hbBy) = derv2(hbHz,hbBy) - d2Edrhb2*d2r2dx2(4)
            derv2(hbHx,hbBz) = derv2(hbHx,hbBz) - d2Edrhb2*d2r2dx2(5)
            derv2(hbHy,hbBz) = derv2(hbHy,hbBz) - d2Edrhb2*d2r2dx2(4)
            derv2(hbHz,hbBz) = derv2(hbHz,hbBz) - d2Edrhb2*d2r2dx2(3)
            derv2(hbHx,hbBx) = derv2(hbHx,hbBx) - dEdrhb
            derv2(hbHy,hbBy) = derv2(hbHy,hbBy) - dEdrhb
            derv2(hbHz,hbBz) = derv2(hbHz,hbBz) - dEdrhb
          endif
!
          if (np.gt.1) then
!
!  Loop over pairs of valid interactions
!
            do np2 = 1,np-1
!
!  Check that if B is an image of A then A=B case is excluded for the same image
!
              if (hbA.eq.hbB) then
                xab2 = cnhb_pairs%xpair(np2) + xbnbr(nh,hbA)
                yab2 = cnhb_pairs%ypair(np2) + ybnbr(nh,hbA)
                zab2 = cnhb_pairs%zpair(np2) + zbnbr(nh,hbA)
                r2ab2 = xab2**2 + yab2**2 + zab2**2
                lABok2 = (r2ab2.gt.1.0d-2)
              else
                lABok2 = .true.
              endif
!
              if (lABok2) then
                rhb2 = sqrt(cnhb_pairs%r2pair(np2))
!
                dr2 = (rhb2 - r0)/r0
                dtrm1a = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr2)**2)/(sqrt(pi)*r0)
                dctrmdrhb2 = - dtrm1a/rhb2
                d2Edrhb2 = d2Edcnh2*dctrmdrhb*dctrmdrhb2
!
                if (hbB.gt.hbH) then
                  derv2(hbBx,hbHx) = derv2(hbBx,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbBy,hbHx) = derv2(hbBy,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbBz,hbHx) = derv2(hbBz,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbBx,hbHy) = derv2(hbBx,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbBy,hbHy) = derv2(hbBy,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbBz,hbHy) = derv2(hbBz,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbBx,hbHz) = derv2(hbBx,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbBy,hbHz) = derv2(hbBy,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbBz,hbHz) = derv2(hbBz,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                else
                  derv2(hbHx,hbBx) = derv2(hbHx,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbHy,hbBx) = derv2(hbHy,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbHz,hbBx) = derv2(hbHz,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbHx,hbBy) = derv2(hbHx,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbHy,hbBy) = derv2(hbHy,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbHz,hbBy) = derv2(hbHz,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbHx,hbBz) = derv2(hbHx,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbHy,hbBz) = derv2(hbHy,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  derv2(hbHz,hbBz) = derv2(hbHz,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                endif
!
              endif
            enddo
          endif
!
!  HB coordination number and rij
!
          d2trm = d2Edcnhdr*dctrmdrhb
!
          if (hbH.ge.i) then
            derv2(hbHx,ix) = derv2(hbHx,ix) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbHy,ix) = derv2(hbHy,ix) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbHz,ix) = derv2(hbHz,ix) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbHx,iy) = derv2(hbHx,iy) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbHy,iy) = derv2(hbHy,iy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbHz,iy) = derv2(hbHz,iy) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbHx,iz) = derv2(hbHx,iz) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbHy,iz) = derv2(hbHy,iz) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbHz,iz) = derv2(hbHz,iz) + d2trm*cnhb_pairs%zpair(np)*zij
          else
            derv2(ix,hbHx) = derv2(ix,hbHx) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(iy,hbHx) = derv2(iy,hbHx) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(iz,hbHx) = derv2(iz,hbHx) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(ix,hbHy) = derv2(ix,hbHy) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(iy,hbHy) = derv2(iy,hbHy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(iz,hbHy) = derv2(iz,hbHy) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(ix,hbHz) = derv2(ix,hbHz) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(iy,hbHz) = derv2(iy,hbHz) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(iz,hbHz) = derv2(iz,hbHz) + d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
          if (hbH.ge.j) then
            derv2(hbHx,jx) = derv2(hbHx,jx) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbHy,jx) = derv2(hbHy,jx) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbHz,jx) = derv2(hbHz,jx) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbHx,jy) = derv2(hbHx,jy) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbHy,jy) = derv2(hbHy,jy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbHz,jy) = derv2(hbHz,jy) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbHx,jz) = derv2(hbHx,jz) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbHy,jz) = derv2(hbHy,jz) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbHz,jz) = derv2(hbHz,jz) - d2trm*cnhb_pairs%zpair(np)*zij
          else
            derv2(jx,hbHx) = derv2(jx,hbHx) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(jy,hbHx) = derv2(jy,hbHx) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(jz,hbHx) = derv2(jz,hbHx) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(jx,hbHy) = derv2(jx,hbHy) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(jy,hbHy) = derv2(jy,hbHy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(jz,hbHy) = derv2(jz,hbHy) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(jx,hbHz) = derv2(jx,hbHz) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(jy,hbHz) = derv2(jy,hbHz) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(jz,hbHz) = derv2(jz,hbHz) - d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
          if (hbB.ge.i) then
            derv2(hbBx,ix) = derv2(hbBx,ix) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbBy,ix) = derv2(hbBy,ix) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbBz,ix) = derv2(hbBz,ix) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbBx,iy) = derv2(hbBx,iy) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbBy,iy) = derv2(hbBy,iy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbBz,iy) = derv2(hbBz,iy) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbBx,iz) = derv2(hbBx,iz) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbBy,iz) = derv2(hbBy,iz) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbBz,iz) = derv2(hbBz,iz) - d2trm*cnhb_pairs%zpair(np)*zij
          else
            derv2(ix,hbBx) = derv2(ix,hbBx) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(iy,hbBx) = derv2(iy,hbBx) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(iz,hbBx) = derv2(iz,hbBx) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(ix,hbBy) = derv2(ix,hbBy) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(iy,hbBy) = derv2(iy,hbBy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(iz,hbBy) = derv2(iz,hbBy) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(ix,hbBz) = derv2(ix,hbBz) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(iy,hbBz) = derv2(iy,hbBz) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(iz,hbBz) = derv2(iz,hbBz) - d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
          if (hbB.ge.j) then
            derv2(hbBx,jx) = derv2(hbBx,jx) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbBy,jx) = derv2(hbBy,jx) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbBz,jx) = derv2(hbBz,jx) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbBx,jy) = derv2(hbBx,jy) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbBy,jy) = derv2(hbBy,jy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbBz,jy) = derv2(hbBz,jy) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbBx,jz) = derv2(hbBx,jz) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbBy,jz) = derv2(hbBy,jz) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbBz,jz) = derv2(hbBz,jz) + d2trm*cnhb_pairs%zpair(np)*zij
          else
            derv2(jx,hbBx) = derv2(jx,hbBx) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(jy,hbBx) = derv2(jy,hbBx) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(jz,hbBx) = derv2(jz,hbBx) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(jx,hbBy) = derv2(jx,hbBy) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(jy,hbBy) = derv2(jy,hbBy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(jz,hbBy) = derv2(jz,hbBy) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(jx,hbBz) = derv2(jx,hbBz) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(jy,hbBz) = derv2(jy,hbBz) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(jz,hbBz) = derv2(jz,hbBz) + d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
!  HB coordination number and coordination number of i
!
          do n = 1,nnbr_cn(i)
            nk = nbrno_cn(n,i)
            rik = rnbr(nk,i)
            k = nbrno(nk,i)
            r0ik = gfnff_rcov(nat(i)) + gfnff_rcov(nat(k))
            drik = (rik - r0ik)/r0ik
            dcnidrik = gfnff_kn_cn*exp(-(gfnff_kn_cn*drik)**2)/(sqrt(pi)*r0ik*rik)
!
            indk = 3*(k-1)
            kx = indk + 1
            ky = indk + 2
            kz = indk + 3
!
            d1trm = d2Edcnhdcni*dcnidrik*dctrmdrhb*dlogcnidcni
!
            if (hbH.ge.i) then
              derv2(hbHx,ix) = derv2(hbHx,ix) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbHy,ix) = derv2(hbHy,ix) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbHz,ix) = derv2(hbHz,ix) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbHx,iy) = derv2(hbHx,iy) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbHy,iy) = derv2(hbHy,iy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbHz,iy) = derv2(hbHz,iy) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbHx,iz) = derv2(hbHx,iz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbHy,iz) = derv2(hbHy,iz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbHz,iz) = derv2(hbHz,iz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            else
              derv2(ix,hbHx) = derv2(ix,hbHx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(iy,hbHx) = derv2(iy,hbHx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(iz,hbHx) = derv2(iz,hbHx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(ix,hbHy) = derv2(ix,hbHy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(iy,hbHy) = derv2(iy,hbHy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(iz,hbHy) = derv2(iz,hbHy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(ix,hbHz) = derv2(ix,hbHz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(iy,hbHz) = derv2(iy,hbHz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(iz,hbHz) = derv2(iz,hbHz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
            if (hbH.ge.k) then
              derv2(hbHx,kx) = derv2(hbHx,kx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbHy,kx) = derv2(hbHy,kx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbHz,kx) = derv2(hbHz,kx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbHx,ky) = derv2(hbHx,ky) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbHy,ky) = derv2(hbHy,ky) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbHz,ky) = derv2(hbHz,ky) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbHx,kz) = derv2(hbHx,kz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbHy,kz) = derv2(hbHy,kz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbHz,kz) = derv2(hbHz,kz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            else
              derv2(kx,hbHx) = derv2(kx,hbHx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(ky,hbHx) = derv2(ky,hbHx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(kz,hbHx) = derv2(kz,hbHx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(kx,hbHy) = derv2(kx,hbHy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(ky,hbHy) = derv2(ky,hbHy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(kz,hbHy) = derv2(kz,hbHy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(kx,hbHz) = derv2(kx,hbHz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(ky,hbHz) = derv2(ky,hbHz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(kz,hbHz) = derv2(kz,hbHz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
            if (hbB.ge.i) then
              derv2(hbBx,ix) = derv2(hbBx,ix) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbBy,ix) = derv2(hbBy,ix) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbBz,ix) = derv2(hbBz,ix) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbBx,iy) = derv2(hbBx,iy) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbBy,iy) = derv2(hbBy,iy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbBz,iy) = derv2(hbBz,iy) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbBx,iz) = derv2(hbBx,iz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbBy,iz) = derv2(hbBy,iz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbBz,iz) = derv2(hbBz,iz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            else
              derv2(ix,hbBx) = derv2(ix,hbBx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(iy,hbBx) = derv2(iy,hbBx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(iz,hbBx) = derv2(iz,hbBx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(ix,hbBy) = derv2(ix,hbBy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(iy,hbBy) = derv2(iy,hbBy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(iz,hbBy) = derv2(iz,hbBy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(ix,hbBz) = derv2(ix,hbBz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(iy,hbBz) = derv2(iy,hbBz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(iz,hbBz) = derv2(iz,hbBz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
            if (hbB.ge.k) then
              derv2(hbBx,kx) = derv2(hbBx,kx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbBy,kx) = derv2(hbBy,kx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbBz,kx) = derv2(hbBz,kx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbBx,ky) = derv2(hbBx,ky) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbBy,ky) = derv2(hbBy,ky) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbBz,ky) = derv2(hbBz,ky) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbBx,kz) = derv2(hbBx,kz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbBy,kz) = derv2(hbBy,kz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbBz,kz) = derv2(hbBz,kz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            else
              derv2(kx,hbBx) = derv2(kx,hbBx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(ky,hbBx) = derv2(ky,hbBx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(kz,hbBx) = derv2(kz,hbBx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(kx,hbBy) = derv2(kx,hbBy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(ky,hbBy) = derv2(ky,hbBy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(kz,hbBy) = derv2(kz,hbBy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(kx,hbBz) = derv2(kx,hbBz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(ky,hbBz) = derv2(ky,hbBz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(kz,hbBz) = derv2(kz,hbBz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
          enddo
!
!  HB coordination number and coordination number of j
!
          do n = 1,nnbr_cn(j)
            nk = nbrno_cn(n,j)
            rjk = rnbr(nk,j)
            k = nbrno(nk,j)
            r0jk = gfnff_rcov(nat(j)) + gfnff_rcov(nat(k))
            drjk = (rjk - r0jk)/r0jk
            dcnjdrjk = gfnff_kn_cn*exp(-(gfnff_kn_cn*drjk)**2)/(sqrt(pi)*r0jk*rjk)
!
            indk = 3*(k-1)
            kx = indk + 1
            ky = indk + 2
            kz = indk + 3
!
            d1trm = d2Edcnhdcnj*dcnjdrjk*dctrmdrhb*dlogcnjdcnj
!
            if (hbH.ge.j) then
              derv2(hbHx,jx) = derv2(hbHx,jx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbHy,jx) = derv2(hbHy,jx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbHz,jx) = derv2(hbHz,jx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbHx,jy) = derv2(hbHx,jy) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbHy,jy) = derv2(hbHy,jy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbHz,jy) = derv2(hbHz,jy) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbHx,jz) = derv2(hbHx,jz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbHy,jz) = derv2(hbHy,jz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbHz,jz) = derv2(hbHz,jz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            else
              derv2(jx,hbHx) = derv2(jx,hbHx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(jy,hbHx) = derv2(jy,hbHx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(jz,hbHx) = derv2(jz,hbHx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(jx,hbHy) = derv2(jx,hbHy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(jy,hbHy) = derv2(jy,hbHy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(jz,hbHy) = derv2(jz,hbHy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(jx,hbHz) = derv2(jx,hbHz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(jy,hbHz) = derv2(jy,hbHz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(jz,hbHz) = derv2(jz,hbHz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
            if (hbH.ge.k) then
              derv2(hbHx,kx) = derv2(hbHx,kx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbHy,kx) = derv2(hbHy,kx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbHz,kx) = derv2(hbHz,kx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbHx,ky) = derv2(hbHx,ky) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbHy,ky) = derv2(hbHy,ky) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbHz,ky) = derv2(hbHz,ky) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbHx,kz) = derv2(hbHx,kz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbHy,kz) = derv2(hbHy,kz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbHz,kz) = derv2(hbHz,kz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            else
              derv2(kx,hbHx) = derv2(kx,hbHx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(ky,hbHx) = derv2(ky,hbHx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(kz,hbHx) = derv2(kz,hbHx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(kx,hbHy) = derv2(kx,hbHy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(ky,hbHy) = derv2(ky,hbHy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(kz,hbHy) = derv2(kz,hbHy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(kx,hbHz) = derv2(kx,hbHz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(ky,hbHz) = derv2(ky,hbHz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(kz,hbHz) = derv2(kz,hbHz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
            if (hbB.ge.j) then
              derv2(hbBx,jx) = derv2(hbBx,jx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbBy,jx) = derv2(hbBy,jx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbBz,jx) = derv2(hbBz,jx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbBx,jy) = derv2(hbBx,jy) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbBy,jy) = derv2(hbBy,jy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbBz,jy) = derv2(hbBz,jy) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbBx,jz) = derv2(hbBx,jz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbBy,jz) = derv2(hbBy,jz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbBz,jz) = derv2(hbBz,jz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            else
              derv2(jx,hbBx) = derv2(jx,hbBx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(jy,hbBx) = derv2(jy,hbBx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(jz,hbBx) = derv2(jz,hbBx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(jx,hbBy) = derv2(jx,hbBy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(jy,hbBy) = derv2(jy,hbBy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(jz,hbBy) = derv2(jz,hbBy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(jx,hbBz) = derv2(jx,hbBz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(jy,hbBz) = derv2(jy,hbBz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(jz,hbBz) = derv2(jz,hbBz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
            if (hbB.ge.k) then
              derv2(hbBx,kx) = derv2(hbBx,kx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbBy,kx) = derv2(hbBy,kx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbBz,kx) = derv2(hbBz,kx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbBx,ky) = derv2(hbBx,ky) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbBy,ky) = derv2(hbBy,ky) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbBz,ky) = derv2(hbBz,ky) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbBx,kz) = derv2(hbBx,kz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbBy,kz) = derv2(hbBy,kz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbBz,kz) = derv2(hbBz,kz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            else
              derv2(kx,hbBx) = derv2(kx,hbBx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(ky,hbBx) = derv2(ky,hbBx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(kz,hbBx) = derv2(kz,hbBx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(kx,hbBy) = derv2(kx,hbBy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(ky,hbBy) = derv2(ky,hbBy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(kz,hbBy) = derv2(kz,hbBy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(kx,hbBz) = derv2(kx,hbBz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(ky,hbBz) = derv2(ky,hbBz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(kz,hbBz) = derv2(kz,hbBz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
          enddo
        endif
      enddo
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcn_hbp')
#endif
!
  return
  end
!
  subroutine gfnff_drv2_dcnd_hbp(hbH,i,j,xij,yij,zij,dEdcnh,d2Edcnh2,d2Edcnhdr,d2Edcnhdcni,d2Edcnhdcnj,dcn)
!
!  Computes the derivatives of the coordination number for GFNFF for hydrogen bonds for atom i-j pair
!  Phonon version. Distributed memory parallel.
!
!  On entry : 
!
!  hbH             = hydrogen atom whose coordination number is to contribute to the derivatives
!  i               = atom whose coordination number is to contribute to the derivatives via rij
!  j               = atom whose coordination number is to contribute to the derivatives via rij
!  xij             = x component of vector from i to j
!  yij             = y component of vector from i to j
!  zij             = z component of vector from i to j
!  dEdcnh          = first derivative of energy w.r.t. HB coordination number
!  d2Edcnh2        = second derivative of energy w.r.t. HB coordination number
!  d2Edcnhdr       = second derivative of energy w.r.t. HB coordination number and rij
!  d2Edcnhdcni     = second derivative of energy w.r.t. coordination numbers of hbH and i
!  d2Edcnhdcnj     = second derivative of energy w.r.t. coordination numbers of hbH and j
!  dcn             = derivative of log coordination numbers w.r.t. coordination numbers
!
!  On exit :
!
!  Derivatives have been computed and added to the appropriate arrays
!
!   4/22 Created from gfnff_drv2_dcn_hbp
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
!  Julian Gale, CIC, Curtin University, April 2022
!
  use datatypes
  use current
  use derivatives
  use g_constants,    only : pi
  use gulp_gfnff
  use iochannels
  use m_gfnff_nbr
  use m_gfnff_pairs
  use parallel
  use spatial
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: hbH
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: j
  real(dp),    intent(in)                          :: dEdcnh
  real(dp),    intent(in)                          :: d2Edcnh2
  real(dp),    intent(in)                          :: d2Edcnhdr
  real(dp),    intent(in)                          :: d2Edcnhdcni
  real(dp),    intent(in)                          :: d2Edcnhdcnj
  real(dp),    intent(in)                          :: dcn(*)
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
!
!  Local variables
!
  integer(i4)                                      :: hbA
  integer(i4)                                      :: hbB
  integer(i4)                                      :: hbBloc
  integer(i4)                                      :: hbBx
  integer(i4)                                      :: hbBy
  integer(i4)                                      :: hbBz
  integer(i4)                                      :: hbBxf
  integer(i4)                                      :: hbByf
  integer(i4)                                      :: hbBzf
  integer(i4)                                      :: hbHloc
  integer(i4)                                      :: hbHx
  integer(i4)                                      :: hbHy
  integer(i4)                                      :: hbHz
  integer(i4)                                      :: hbHxf
  integer(i4)                                      :: hbHyf
  integer(i4)                                      :: hbHzf
  integer(i4)                                      :: iloc
  integer(i4)                                      :: indh
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: ixf
  integer(i4)                                      :: iyf
  integer(i4)                                      :: izf
  integer(i4)                                      :: jloc
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: jxf
  integer(i4)                                      :: jyf
  integer(i4)                                      :: jzf
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: n
  integer(i4)                                      :: natB
  integer(i4)                                      :: natH
  integer(i4)                                      :: nh
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nk
  integer(i4)                                      :: np
  integer(i4)                                      :: np2
  logical                                          :: lABok
  logical                                          :: lABok2
  logical                                          :: lilocal
  logical                                          :: ljlocal
  logical                                          :: lklocal
  logical                                          :: lhbBlocal
  logical                                          :: lhbHlocal
  real(dp)                                         :: cncut2
  real(dp)                                         :: dctrmdrhb
  real(dp)                                         :: dctrmdrhb2
  real(dp)                                         :: d2ctrmdrhb2
  real(dp)                                         :: dcnidrik
  real(dp)                                         :: dcnjdrjk
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcnjdcnj
  real(dp)                                         :: dEdrhb
  real(dp)                                         :: d2Edrhb2
  real(dp)                                         :: d2r2dx2(6)
  real(dp)                                         :: dtrm
  real(dp)                                         :: dtrm1a
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: dr
  real(dp)                                         :: dr2
  real(dp)                                         :: drik
  real(dp)                                         :: drjk
  real(dp),                          parameter     :: rcov_scale = 1.78_dp
  real(dp)                                         :: r0
  real(dp)                                         :: r0ik
  real(dp)                                         :: r0jk
  real(dp)                                         :: rik
  real(dp)                                         :: rjk
  real(dp)                                         :: rhb
  real(dp)                                         :: rhb2
  real(dp)                                         :: r2ab
  real(dp)                                         :: r2ab2
  real(dp)                                         :: xab
  real(dp)                                         :: yab
  real(dp)                                         :: zab
  real(dp)                                         :: xab2
  real(dp)                                         :: yab2
  real(dp)                                         :: zab2
#ifdef TRACE
  call trace_in('gfnff_drv2_dcnd_hbp')
#endif
!
!  Is i or j or hbH local
!
  iloc = atom2local(i)
  jloc = atom2local(j)
  hbHloc = atom2local(hbH)
  lilocal = (iloc.ne.0)
  ljlocal = (jloc.ne.0)
  lhbHlocal = (hbHloc.ne.0)
!
!  Set up for second derivatives
!
  indi = 3*(i-1)
  ixf = indi + 1
  iyf = indi + 2
  izf = indi + 3
  if (lilocal) then
    indi = 3*(iloc-1)
    ix = indi + 1
    iy = indi + 2
    iz = indi + 3
  endif
!
  indj = 3*(j-1)
  jxf = indj + 1
  jyf = indj + 2
  jzf = indj + 3
  if (ljlocal) then
    indj = 3*(jloc-1)
    jx = indj + 1
    jy = indj + 2
    jz = indj + 3
  endif
!
  indh = 3*(hbH-1)
  hbHxf = indh + 1
  hbHyf = indh + 2
  hbHzf = indh + 3
  if (lhbHlocal) then
    indh = 3*(hbHloc-1)
    hbHx = indh + 1
    hbHy = indh + 2
    hbHz = indh + 3
  endif
!
  dlogcnidcni = dcn(i)
  dlogcnjdcnj = dcn(j)
!
  do ni = 1,nbond_hb_nr
! DEBUG - look for faster approach here!
    if (hbH.ne.nbond_hb_AH(2,ni)) cycle
    hbA = nbond_hb_AH(1,ni)
    nh  = nbond_hb_AH(3,ni)
    natH = nat(hbH)
    do nj = 1,nbond_hb_Bn(ni)
      hbB = nbond_hb_B(nj,ni)
      hbBloc = atom2local(hbB)
      lhbBlocal = (hbBloc.ne.0)
!
      if (lhbBlocal) then
        indh = 3*(hbBloc-1)
        hbBx = indh + 1
        hbBy = indh + 2
        hbBz = indh + 3
      endif
!
      indh = 3*(hbB-1)
      hbBxf = indh + 1
      hbByf = indh + 2
      hbBzf = indh + 3
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
!
          dr = (rhb - r0)/r0
          dtrm = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr)**2)/(sqrt(pi)*r0)
          dctrmdrhb = - dtrm/rhb
!
          dEdrhb = dEdcnh*dctrmdrhb
!
!  Check whether any atoms are local for the next section
!
          if (lhbBlocal.or.lhbHlocal) then
            d2r2dx2(1) = cnhb_pairs%xpair(np)*cnhb_pairs%xpair(np)
            d2r2dx2(2) = cnhb_pairs%ypair(np)*cnhb_pairs%ypair(np)
            d2r2dx2(3) = cnhb_pairs%zpair(np)*cnhb_pairs%zpair(np)
            d2r2dx2(4) = cnhb_pairs%ypair(np)*cnhb_pairs%zpair(np)
            d2r2dx2(5) = cnhb_pairs%xpair(np)*cnhb_pairs%zpair(np)
            d2r2dx2(6) = cnhb_pairs%xpair(np)*cnhb_pairs%ypair(np)
!-----------------------
!  Second derivatives  |
!-----------------------
!
!  HB coordination number only
!
            d2ctrmdrhb2 = dtrm*(1.0_dp + 2.0_dp*rhb*dr*(gfnff_kn_hb**2)/r0)/rhb**3
            d2Edrhb2 = d2Edcnh2*dctrmdrhb*dctrmdrhb + dEdcnh*d2ctrmdrhb2 ! (d2E/dr2)
!
            if (lhbHlocal) then
              derv2(hbBxf,hbHx) = derv2(hbBxf,hbHx) - d2Edrhb2*d2r2dx2(1)
              derv2(hbByf,hbHx) = derv2(hbByf,hbHx) - d2Edrhb2*d2r2dx2(6)
              derv2(hbBzf,hbHx) = derv2(hbBzf,hbHx) - d2Edrhb2*d2r2dx2(5)
              derv2(hbBxf,hbHy) = derv2(hbBxf,hbHy) - d2Edrhb2*d2r2dx2(6)
              derv2(hbByf,hbHy) = derv2(hbByf,hbHy) - d2Edrhb2*d2r2dx2(2)
              derv2(hbBzf,hbHy) = derv2(hbBzf,hbHy) - d2Edrhb2*d2r2dx2(4)
              derv2(hbBxf,hbHz) = derv2(hbBxf,hbHz) - d2Edrhb2*d2r2dx2(5)
              derv2(hbByf,hbHz) = derv2(hbByf,hbHz) - d2Edrhb2*d2r2dx2(4)
              derv2(hbBzf,hbHz) = derv2(hbBzf,hbHz) - d2Edrhb2*d2r2dx2(3)
              derv2(hbBxf,hbHx) = derv2(hbBxf,hbHx) - dEdrhb
              derv2(hbByf,hbHy) = derv2(hbByf,hbHy) - dEdrhb
              derv2(hbBzf,hbHz) = derv2(hbBzf,hbHz) - dEdrhb
            endif
!
            if (lhbBlocal) then
              derv2(hbHxf,hbBx) = derv2(hbHxf,hbBx) - d2Edrhb2*d2r2dx2(1)
              derv2(hbHyf,hbBx) = derv2(hbHyf,hbBx) - d2Edrhb2*d2r2dx2(6)
              derv2(hbHzf,hbBx) = derv2(hbHzf,hbBx) - d2Edrhb2*d2r2dx2(5)
              derv2(hbHxf,hbBy) = derv2(hbHxf,hbBy) - d2Edrhb2*d2r2dx2(6)
              derv2(hbHyf,hbBy) = derv2(hbHyf,hbBy) - d2Edrhb2*d2r2dx2(2)
              derv2(hbHzf,hbBy) = derv2(hbHzf,hbBy) - d2Edrhb2*d2r2dx2(4)
              derv2(hbHxf,hbBz) = derv2(hbHxf,hbBz) - d2Edrhb2*d2r2dx2(5)
              derv2(hbHyf,hbBz) = derv2(hbHyf,hbBz) - d2Edrhb2*d2r2dx2(4)
              derv2(hbHzf,hbBz) = derv2(hbHzf,hbBz) - d2Edrhb2*d2r2dx2(3)
              derv2(hbHxf,hbBx) = derv2(hbHxf,hbBx) - dEdrhb
              derv2(hbHyf,hbBy) = derv2(hbHyf,hbBy) - dEdrhb
              derv2(hbHzf,hbBz) = derv2(hbHzf,hbBz) - dEdrhb
            endif
!
            if (np.gt.1) then
!
!  Loop over pairs of valid interactions
!
              do np2 = 1,np-1
!
!  Check that if B is an image of A then A=B case is excluded for the same image
!
                if (hbA.eq.hbB) then
                  xab2 = cnhb_pairs%xpair(np2) + xbnbr(nh,hbA)
                  yab2 = cnhb_pairs%ypair(np2) + ybnbr(nh,hbA)
                  zab2 = cnhb_pairs%zpair(np2) + zbnbr(nh,hbA)
                  r2ab2 = xab2**2 + yab2**2 + zab2**2
                  lABok2 = (r2ab2.gt.1.0d-2)
                else
                  lABok2 = .true.
                endif
!
                if (lABok2) then
                  rhb2 = sqrt(cnhb_pairs%r2pair(np2))
!
                  dr2 = (rhb2 - r0)/r0
                  dtrm1a = gfnff_kn_hb*exp(-(gfnff_kn_hb*dr2)**2)/(sqrt(pi)*r0)
                  dctrmdrhb2 = - dtrm1a/rhb2
                  d2Edrhb2 = d2Edcnh2*dctrmdrhb*dctrmdrhb2
!
                  if (lhbHlocal) then
                    derv2(hbBxf,hbHx) = derv2(hbBxf,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbByf,hbHx) = derv2(hbByf,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBzf,hbHx) = derv2(hbBzf,hbHx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBxf,hbHy) = derv2(hbBxf,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbByf,hbHy) = derv2(hbByf,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBzf,hbHy) = derv2(hbBzf,hbHy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBxf,hbHz) = derv2(hbBxf,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbByf,hbHz) = derv2(hbByf,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbBzf,hbHz) = derv2(hbBzf,hbHz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  endif
!
                  if (lhbBlocal) then
                    derv2(hbHxf,hbBx) = derv2(hbHxf,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHyf,hbBx) = derv2(hbHyf,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHzf,hbBx) = derv2(hbHzf,hbBx) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHxf,hbBy) = derv2(hbHxf,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHyf,hbBy) = derv2(hbHyf,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHzf,hbBy) = derv2(hbHzf,hbBy) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHxf,hbBz) = derv2(hbHxf,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHyf,hbBz) = derv2(hbHyf,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                    derv2(hbHzf,hbBz) = derv2(hbHzf,hbBz) - d2Edrhb2*cnhb_pairs%xpair(np2)*cnhb_pairs%xpair(np)
                  endif
                endif
              enddo
            endif
          endif
!
!  HB coordination number and rij
!
          d2trm = d2Edcnhdr*dctrmdrhb
!
          if (lilocal) then
            derv2(hbHxf,ix) = derv2(hbHxf,ix) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbHyf,ix) = derv2(hbHyf,ix) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbHzf,ix) = derv2(hbHzf,ix) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbHxf,iy) = derv2(hbHxf,iy) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbHyf,iy) = derv2(hbHyf,iy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbHzf,iy) = derv2(hbHzf,iy) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbHxf,iz) = derv2(hbHxf,iz) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbHyf,iz) = derv2(hbHyf,iz) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbHzf,iz) = derv2(hbHzf,iz) + d2trm*cnhb_pairs%zpair(np)*zij
!
            derv2(hbBxf,ix) = derv2(hbBxf,ix) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbByf,ix) = derv2(hbByf,ix) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbBzf,ix) = derv2(hbBzf,ix) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbBxf,iy) = derv2(hbBxf,iy) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbByf,iy) = derv2(hbByf,iy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbBzf,iy) = derv2(hbBzf,iy) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbBxf,iz) = derv2(hbBxf,iz) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbByf,iz) = derv2(hbByf,iz) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbBzf,iz) = derv2(hbBzf,iz) - d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
          if (lhbHlocal) then
            derv2(ixf,hbHx) = derv2(ixf,hbHx) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(iyf,hbHx) = derv2(iyf,hbHx) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(izf,hbHx) = derv2(izf,hbHx) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(ixf,hbHy) = derv2(ixf,hbHy) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(iyf,hbHy) = derv2(iyf,hbHy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(izf,hbHy) = derv2(izf,hbHy) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(ixf,hbHz) = derv2(ixf,hbHz) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(iyf,hbHz) = derv2(iyf,hbHz) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(izf,hbHz) = derv2(izf,hbHz) + d2trm*cnhb_pairs%zpair(np)*zij
!
            derv2(jxf,hbHx) = derv2(jxf,hbHx) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(jyf,hbHx) = derv2(jyf,hbHx) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(jzf,hbHx) = derv2(jzf,hbHx) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(jxf,hbHy) = derv2(jxf,hbHy) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(jyf,hbHy) = derv2(jyf,hbHy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(jzf,hbHy) = derv2(jzf,hbHy) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(jxf,hbHz) = derv2(jxf,hbHz) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(jyf,hbHz) = derv2(jyf,hbHz) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(jzf,hbHz) = derv2(jzf,hbHz) - d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
          if (lhbBlocal) then
            derv2(ixf,hbBx) = derv2(ixf,hbBx) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(iyf,hbBx) = derv2(iyf,hbBx) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(izf,hbBx) = derv2(izf,hbBx) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(ixf,hbBy) = derv2(ixf,hbBy) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(iyf,hbBy) = derv2(iyf,hbBy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(izf,hbBy) = derv2(izf,hbBy) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(ixf,hbBz) = derv2(ixf,hbBz) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(iyf,hbBz) = derv2(iyf,hbBz) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(izf,hbBz) = derv2(izf,hbBz) - d2trm*cnhb_pairs%zpair(np)*zij
!
            derv2(jxf,hbBx) = derv2(jxf,hbBx) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(jyf,hbBx) = derv2(jyf,hbBx) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(jzf,hbBx) = derv2(jzf,hbBx) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(jxf,hbBy) = derv2(jxf,hbBy) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(jyf,hbBy) = derv2(jyf,hbBy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(jzf,hbBy) = derv2(jzf,hbBy) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(jxf,hbBz) = derv2(jxf,hbBz) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(jyf,hbBz) = derv2(jyf,hbBz) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(jzf,hbBz) = derv2(jzf,hbBz) + d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
          if (ljlocal) then
            derv2(hbBxf,jx) = derv2(hbBxf,jx) + d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbByf,jx) = derv2(hbByf,jx) + d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbBzf,jx) = derv2(hbBzf,jx) + d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbBxf,jy) = derv2(hbBxf,jy) + d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbByf,jy) = derv2(hbByf,jy) + d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbBzf,jy) = derv2(hbBzf,jy) + d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbBxf,jz) = derv2(hbBxf,jz) + d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbByf,jz) = derv2(hbByf,jz) + d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbBzf,jz) = derv2(hbBzf,jz) + d2trm*cnhb_pairs%zpair(np)*zij
!
            derv2(hbHxf,jx) = derv2(hbHxf,jx) - d2trm*cnhb_pairs%xpair(np)*xij
            derv2(hbHyf,jx) = derv2(hbHyf,jx) - d2trm*cnhb_pairs%ypair(np)*xij
            derv2(hbHzf,jx) = derv2(hbHzf,jx) - d2trm*cnhb_pairs%zpair(np)*xij
            derv2(hbHxf,jy) = derv2(hbHxf,jy) - d2trm*cnhb_pairs%xpair(np)*yij
            derv2(hbHyf,jy) = derv2(hbHyf,jy) - d2trm*cnhb_pairs%ypair(np)*yij
            derv2(hbHzf,jy) = derv2(hbHzf,jy) - d2trm*cnhb_pairs%zpair(np)*yij
            derv2(hbHxf,jz) = derv2(hbHxf,jz) - d2trm*cnhb_pairs%xpair(np)*zij
            derv2(hbHyf,jz) = derv2(hbHyf,jz) - d2trm*cnhb_pairs%ypair(np)*zij
            derv2(hbHzf,jz) = derv2(hbHzf,jz) - d2trm*cnhb_pairs%zpair(np)*zij
          endif
!
!  HB coordination number and coordination number of i
!
          do n = 1,nnbr_cn(i)
            nk = nbrno_cn(n,i)
            rik = rnbr(nk,i)
            k = nbrno(nk,i)
!
            kloc = atom2local(k)
            lklocal = (kloc.ne.0)
!
!  Check whether any atoms are local
!
            if (.not.lilocal.and..not.lklocal.and..not.lhbBlocal.and..not.lhbHlocal) cycle
!
            r0ik = gfnff_rcov(nat(i)) + gfnff_rcov(nat(k))
            drik = (rik - r0ik)/r0ik
            dcnidrik = gfnff_kn_cn*exp(-(gfnff_kn_cn*drik)**2)/(sqrt(pi)*r0ik*rik)
!
            indk = 3*(k-1)
            kxf = indk + 1
            kyf = indk + 2
            kzf = indk + 3
            if (lklocal) then
              indk = 3*(kloc-1)
              kx = indk + 1
              ky = indk + 2
              kz = indk + 3
            endif
!
            d1trm = d2Edcnhdcni*dcnidrik*dctrmdrhb*dlogcnidcni
!
            if (lilocal) then
              derv2(hbHxf,ix) = derv2(hbHxf,ix) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbHyf,ix) = derv2(hbHyf,ix) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbHzf,ix) = derv2(hbHzf,ix) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbHxf,iy) = derv2(hbHxf,iy) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbHyf,iy) = derv2(hbHyf,iy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbHzf,iy) = derv2(hbHzf,iy) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbHxf,iz) = derv2(hbHxf,iz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbHyf,iz) = derv2(hbHyf,iz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbHzf,iz) = derv2(hbHzf,iz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
!
              derv2(hbBxf,ix) = derv2(hbBxf,ix) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbByf,ix) = derv2(hbByf,ix) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbBzf,ix) = derv2(hbBzf,ix) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbBxf,iy) = derv2(hbBxf,iy) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbByf,iy) = derv2(hbByf,iy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbBzf,iy) = derv2(hbBzf,iy) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbBxf,iz) = derv2(hbBxf,iz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbByf,iz) = derv2(hbByf,iz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbBzf,iz) = derv2(hbBzf,iz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
            if (lhbHlocal) then
              derv2(ixf,hbHx) = derv2(ixf,hbHx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(iyf,hbHx) = derv2(iyf,hbHx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(izf,hbHx) = derv2(izf,hbHx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(ixf,hbHy) = derv2(ixf,hbHy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(iyf,hbHy) = derv2(iyf,hbHy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(izf,hbHy) = derv2(izf,hbHy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(ixf,hbHz) = derv2(ixf,hbHz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(iyf,hbHz) = derv2(iyf,hbHz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(izf,hbHz) = derv2(izf,hbHz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
!
              derv2(kxf,hbHx) = derv2(kxf,hbHx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(kyf,hbHx) = derv2(kyf,hbHx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(kzf,hbHx) = derv2(kzf,hbHx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(kxf,hbHy) = derv2(kxf,hbHy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(kyf,hbHy) = derv2(kyf,hbHy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(kzf,hbHy) = derv2(kzf,hbHy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(kxf,hbHz) = derv2(kxf,hbHz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(kyf,hbHz) = derv2(kyf,hbHz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(kzf,hbHz) = derv2(kzf,hbHz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
            if (lhbBlocal) then
              derv2(ixf,hbBx) = derv2(ixf,hbBx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(iyf,hbBx) = derv2(iyf,hbBx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(izf,hbBx) = derv2(izf,hbBx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(ixf,hbBy) = derv2(ixf,hbBy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(iyf,hbBy) = derv2(iyf,hbBy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(izf,hbBy) = derv2(izf,hbBy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(ixf,hbBz) = derv2(ixf,hbBz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(iyf,hbBz) = derv2(iyf,hbBz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(izf,hbBz) = derv2(izf,hbBz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
!
              derv2(kxf,hbBx) = derv2(kxf,hbBx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(kyf,hbBx) = derv2(kyf,hbBx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(kzf,hbBx) = derv2(kzf,hbBx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(kxf,hbBy) = derv2(kxf,hbBy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(kyf,hbBy) = derv2(kyf,hbBy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(kzf,hbBy) = derv2(kzf,hbBy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(kxf,hbBz) = derv2(kxf,hbBz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(kyf,hbBz) = derv2(kyf,hbBz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(kzf,hbBz) = derv2(kzf,hbBz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
!
            if (lklocal) then
              derv2(hbHxf,kx) = derv2(hbHxf,kx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbHyf,kx) = derv2(hbHyf,kx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbHzf,kx) = derv2(hbHzf,kx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbHxf,ky) = derv2(hbHxf,ky) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbHyf,ky) = derv2(hbHyf,ky) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbHzf,ky) = derv2(hbHzf,ky) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbHxf,kz) = derv2(hbHxf,kz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbHyf,kz) = derv2(hbHyf,kz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbHzf,kz) = derv2(hbHzf,kz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
!
              derv2(hbBxf,kx) = derv2(hbBxf,kx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,i)
              derv2(hbByf,kx) = derv2(hbByf,kx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,i)
              derv2(hbBzf,kx) = derv2(hbBzf,kx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,i)
              derv2(hbBxf,ky) = derv2(hbBxf,ky) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,i)
              derv2(hbByf,ky) = derv2(hbByf,ky) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,i)
              derv2(hbBzf,ky) = derv2(hbBzf,ky) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,i)
              derv2(hbBxf,kz) = derv2(hbBxf,kz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,i)
              derv2(hbByf,kz) = derv2(hbByf,kz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,i)
              derv2(hbBzf,kz) = derv2(hbBzf,kz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,i)
            endif
          enddo
!
!  HB coordination number and coordination number of j
!
          do n = 1,nnbr_cn(j)
            nk = nbrno_cn(n,j)
            rjk = rnbr(nk,j)
            k = nbrno(nk,j)
!
            kloc = atom2local(k)
            lklocal = (kloc.ne.0)
!
!  Check whether any atoms are local
!
            if (.not.ljlocal.and..not.lklocal.and..not.lhbBlocal.and..not.lhbHlocal) cycle
!
            r0jk = gfnff_rcov(nat(j)) + gfnff_rcov(nat(k))
            drjk = (rjk - r0jk)/r0jk
            dcnjdrjk = gfnff_kn_cn*exp(-(gfnff_kn_cn*drjk)**2)/(sqrt(pi)*r0jk*rjk)
!
            indk = 3*(k-1)
            kxf = indk + 1
            kyf = indk + 2
            kzf = indk + 3
            if (lklocal) then
              indk = 3*(kloc-1)
              kx = indk + 1
              ky = indk + 2
              kz = indk + 3
            endif
!
            d1trm = d2Edcnhdcnj*dcnjdrjk*dctrmdrhb*dlogcnjdcnj
!
            if (ljlocal) then
              derv2(hbHxf,jx) = derv2(hbHxf,jx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbHyf,jx) = derv2(hbHyf,jx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbHzf,jx) = derv2(hbHzf,jx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbHxf,jy) = derv2(hbHxf,jy) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbHyf,jy) = derv2(hbHyf,jy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbHzf,jy) = derv2(hbHzf,jy) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbHxf,jz) = derv2(hbHxf,jz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbHyf,jz) = derv2(hbHyf,jz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbHzf,jz) = derv2(hbHzf,jz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
!
              derv2(hbBxf,jx) = derv2(hbBxf,jx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbByf,jx) = derv2(hbByf,jx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbBzf,jx) = derv2(hbBzf,jx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbBxf,jy) = derv2(hbBxf,jy) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbByf,jy) = derv2(hbByf,jy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbBzf,jy) = derv2(hbBzf,jy) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbBxf,jz) = derv2(hbBxf,jz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbByf,jz) = derv2(hbByf,jz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbBzf,jz) = derv2(hbBzf,jz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
            if (lhbHlocal) then
              derv2(jxf,hbHx) = derv2(jxf,hbHx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(jyf,hbHx) = derv2(jyf,hbHx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(jzf,hbHx) = derv2(jzf,hbHx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(jxf,hbHy) = derv2(jxf,hbHy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(jyf,hbHy) = derv2(jyf,hbHy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(jzf,hbHy) = derv2(jzf,hbHy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(jxf,hbHz) = derv2(jxf,hbHz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(jyf,hbHz) = derv2(jyf,hbHz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(jzf,hbHz) = derv2(jzf,hbHz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
!
              derv2(kxf,hbHx) = derv2(kxf,hbHx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(kyf,hbHx) = derv2(kyf,hbHx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(kzf,hbHx) = derv2(kzf,hbHx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(kxf,hbHy) = derv2(kxf,hbHy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(kyf,hbHy) = derv2(kyf,hbHy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(kzf,hbHy) = derv2(kzf,hbHy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(kxf,hbHz) = derv2(kxf,hbHz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(kyf,hbHz) = derv2(kyf,hbHz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(kzf,hbHz) = derv2(kzf,hbHz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
            if (lhbBlocal) then
              derv2(jxf,hbBx) = derv2(jxf,hbBx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(jyf,hbBx) = derv2(jyf,hbBx) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(jzf,hbBx) = derv2(jzf,hbBx) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(jxf,hbBy) = derv2(jxf,hbBy) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(jyf,hbBy) = derv2(jyf,hbBy) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(jzf,hbBy) = derv2(jzf,hbBy) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(jxf,hbBz) = derv2(jxf,hbBz) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(jyf,hbBz) = derv2(jyf,hbBz) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(jzf,hbBz) = derv2(jzf,hbBz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
!
              derv2(kxf,hbBx) = derv2(kxf,hbBx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(kyf,hbBx) = derv2(kyf,hbBx) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(kzf,hbBx) = derv2(kzf,hbBx) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(kxf,hbBy) = derv2(kxf,hbBy) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(kyf,hbBy) = derv2(kyf,hbBy) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(kzf,hbBy) = derv2(kzf,hbBy) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(kxf,hbBz) = derv2(kxf,hbBz) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(kyf,hbBz) = derv2(kyf,hbBz) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(kzf,hbBz) = derv2(kzf,hbBz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
!
            if (lklocal) then
              derv2(hbBxf,kx) = derv2(hbBxf,kx) + d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbByf,kx) = derv2(hbByf,kx) + d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbBzf,kx) = derv2(hbBzf,kx) + d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbBxf,ky) = derv2(hbBxf,ky) + d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbByf,ky) = derv2(hbByf,ky) + d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbBzf,ky) = derv2(hbBzf,ky) + d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbBxf,kz) = derv2(hbBxf,kz) + d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbByf,kz) = derv2(hbByf,kz) + d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbBzf,kz) = derv2(hbBzf,kz) + d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            else
              derv2(hbHxf,kx) = derv2(hbHxf,kx) - d1trm*cnhb_pairs%xpair(np)*xnbr(nk,j)
              derv2(hbHyf,kx) = derv2(hbHyf,kx) - d1trm*cnhb_pairs%ypair(np)*xnbr(nk,j)
              derv2(hbHzf,kx) = derv2(hbHzf,kx) - d1trm*cnhb_pairs%zpair(np)*xnbr(nk,j)
              derv2(hbHxf,ky) = derv2(hbHxf,ky) - d1trm*cnhb_pairs%xpair(np)*ynbr(nk,j)
              derv2(hbHyf,ky) = derv2(hbHyf,ky) - d1trm*cnhb_pairs%ypair(np)*ynbr(nk,j)
              derv2(hbHzf,ky) = derv2(hbHzf,ky) - d1trm*cnhb_pairs%zpair(np)*ynbr(nk,j)
              derv2(hbHxf,kz) = derv2(hbHxf,kz) - d1trm*cnhb_pairs%xpair(np)*znbr(nk,j)
              derv2(hbHyf,kz) = derv2(hbHyf,kz) - d1trm*cnhb_pairs%ypair(np)*znbr(nk,j)
              derv2(hbHzf,kz) = derv2(hbHzf,kz) - d1trm*cnhb_pairs%zpair(np)*znbr(nk,j)
            endif
          enddo
        endif
      enddo
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_drv2_dcnd_hbp')
#endif
!
  return
  end
