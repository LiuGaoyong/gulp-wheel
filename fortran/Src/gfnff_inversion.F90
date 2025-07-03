  subroutine gfnff_inversion(ndim,numat,nnbr,maxnbr,rnbr,xnbr,ynbr,znbr,i,n1,n2,n3,phi,phi1dx,phi1ds,phi2dx,phi2ds,phi2dsdx, &
                             lgrad1,lgrad2)
!
!  Compute improper torsion angle
!  NB: GFNFF uses the angle between the cross product of 2 vectors and a third one
!      which is not the conventional torsion angle.
!  The derivatives are returned in Cartesian form for the vectors ij, jk and il
!
!  Uses vectors from the neighbour list to ensure correct images of atoms
!
!  11/20 Sign change for j->i to i->j moved to dvndij (instead of dvndji)
!
!  Julian Gale, Curtin University, November 2021
!
  use datatypes
  use current,      only : nstrains
  use m_strain,     only : real1strterm, cartstrterm
  use symmetry,     only : lstr
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: n1
  integer(i4), intent(in)  :: n2
  integer(i4), intent(in)  :: n3
  integer(i4), intent(in)  :: nnbr(numat)
  integer(i4), intent(in)  :: maxnbr
  logical,     intent(in)  :: lgrad1
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(in)  :: rnbr(maxnbr,numat)
  real(dp),    intent(in)  :: xnbr(maxnbr,numat)
  real(dp),    intent(in)  :: ynbr(maxnbr,numat)
  real(dp),    intent(in)  :: znbr(maxnbr,numat)
  real(dp),    intent(out) :: phi
  real(dp),    intent(out) :: phi1dx(3,3)              ! If lgrad1 is true then this contains the Cartesian derivatives of phi
  real(dp),    intent(out) :: phi1ds(6)                ! If lgrad1 is true then this contains the strain derivatives of phi
  real(dp),    intent(out) :: phi2dx(3,3,6)            ! If lgrad2 is true then this contains the Cartesian second derivatives of phi
  real(dp),    intent(out) :: phi2ds(6,6)              ! If lgrad2 is true then this contains the strain second derivatives of phi
  real(dp),    intent(out) :: phi2dsdx(6,3,3)          ! If lgrad2 is true then this contains the mixed strain-Cartesian second derivatives of phi
!
!  Local variables
!
  integer(i4)              :: ii
  integer(i4)              :: ind
  integer(i4)              :: ix
  integer(i4)              :: jj
  integer(i4)              :: jx
  integer(i4)              :: ks
  integer(i4)              :: kt
  logical                  :: ldistOK
  real(dp)                 :: vil(3)
  real(dp)                 :: vjk(3)
  real(dp)                 :: vji(3)
  real(dp)                 :: vn(3)
  real(dp)                 :: dphidrnv
  real(dp)                 :: d2phidrnv2
  real(dp)                 :: drnvdx(3,3)        ! First derivatives of rnv for 3 vectors
  real(dp)                 :: drnvds(6)          ! First derivatives of rnv w.r.t. strain
  real(dp)                 :: d2rnvdx2(3,3,6)    ! Second derivatives of rnv for 3 vectors
  real(dp)                 :: d2rnvds2(6,6)      ! Second derivatives of rnv for 2 strains
  real(dp)                 :: d2rnvdsdx(6,3,3)   ! Second derivatives of rnv for 1 strain and 1 Cartesian component
  real(dp)                 :: dvndij(3,3)        ! Derivatives of vn components (right index) w.r.t. vji components (left index)
  real(dp)                 :: dvndjk(3,3)        ! Derivatives of vn components (right index) w.r.t. vjk components (left index)
  real(dp)                 :: dvnds(6,3)         ! Derivatives of vn components (right index) w.r.t. strain components (left index)
  real(dp)                 :: d2vnds2(6,6,3)     ! Second derivatives of vn components (right index) w.r.t. two strain components (left indices)
  real(dp)                 :: d2vndsdij(6,3,3)   ! Second derivatives of vn components (right index) w.r.t. strain and Cartesian components (left/middle indices)
  real(dp)                 :: d2vndsdjk(6,3,3)   ! Second derivatives of vn components (right index) w.r.t. strain and Cartesian components (left/middle indices)
  real(dp)                 :: d2vndijdjk(3,3,3)  ! Second derivatives of vn components (right index) w.r.t. vij/vjk components (middle/left index)
  real(dp)                 :: dvilds(6,3)
  real(dp)                 :: d2vildsdx(6,3,3)
  real(dp)                 :: d2vilds2(6,6,3)
  real(dp)                 :: dvjids(6,3)
  real(dp)                 :: d2vjidsdx(6,3,3)
  real(dp)                 :: d2vjids2(6,6,3)
  real(dp)                 :: dvjkds(6,3)
  real(dp)                 :: d2vjkdsdx(6,3,3)
  real(dp)                 :: d2vjkds2(6,6,3)
  real(dp)                 :: dvilvilds(6)
  real(dp)                 :: dvilvnds(6)
  real(dp)                 :: dvnvilds(6)
  real(dp)                 :: dvnvnds(6)
  real(dp)                 :: dvnvildij(3)
  real(dp)                 :: dvnvildjk(3)
  real(dp)                 :: dvnvndij(3)
  real(dp)                 :: dvnvndjk(3)
  real(dp)                 :: ril
  real(dp)                 :: ril2
  real(dp)                 :: rn
  real(dp)                 :: rn2
  real(dp)                 :: rrnil
  real(dp)                 :: rnv
  real(dp)                 :: rnvrril2
  real(dp)                 :: rnvrrn2
!
!  Check validity of n1, n2 and n3
!
  if (n1.gt.nnbr(i).or.n2.gt.nnbr(i).or.n3.gt.nnbr(i)) then
    call outerror('invalid arguments for bonds to i in gfnff_inversion',0_i4)
    call stopnow('gfnff_inversion')
  endif
!
!  j -> i vector
!
  vji(1) = - xnbr(n1,i)
  vji(2) = - ynbr(n1,i)
  vji(3) = - znbr(n1,i)
!
!  j -> k vector
!
  vjk(1) = xnbr(n2,i) - xnbr(n1,i)
  vjk(2) = ynbr(n2,i) - ynbr(n1,i)
  vjk(3) = znbr(n2,i) - znbr(n1,i)
!
!  i -> l vector
!
  vil(1) = xnbr(n3,i)
  vil(2) = ynbr(n3,i)
  vil(3) = znbr(n3,i)
!
!  Compute cross product of rji and rjk
!
  vn(1) = vji(2)*vjk(3) - vji(3)*vjk(2)
  vn(2) = vji(3)*vjk(1) - vji(1)*vjk(3)
  vn(3) = vji(1)*vjk(2) - vji(2)*vjk(1)
!
  rn2 = vn(1)**2 + vn(2)**2 + vn(3)**2
  rn = sqrt(rn2)
!
  ril = rnbr(n3,i)
  ril2 = ril**2
!
  ldistOK = (rn*ril.gt.1.0d-12)
  if (ldistOK) then
    rrnil = 1.0_dp/(rn*ril)
    rnv = (vn(1)*vil(1) + vn(2)*vil(2) + vn(3)*vil(3))*rrnil
  else
    rnv = 0.0_dp
  endif
!
  phi = asin(rnv)
!
  if (lgrad1) then
    phi1dx(1:3,1:3) = 0.0_dp
    phi1ds(1:6) = 0.0_dp
    if (lgrad2) then
      phi2dx(1:3,1:3,1:6) = 0.0_dp
      phi2ds(1:6,1:6) = 0.0_dp
      phi2dsdx(1:6,1:3,1:3) = 0.0_dp
    endif
    if (ldistOK) then
!
!  Compute derivatives of vectors for cross products
!
!  NB: Sign changed due to differentiating w.r.t. vij rather than vji
!
      dvndij(1,1) = 0.0_dp
      dvndij(2,1) = - vjk(3)
      dvndij(3,1) = vjk(2)
!
      dvndij(1,2) = vjk(3)
      dvndij(2,2) = 0.0_dp
      dvndij(3,2) = - vjk(1)
!
      dvndij(1,3) = - vjk(2)
      dvndij(2,3) = vjk(1)
      dvndij(3,3) = 0.0_dp
!
      dvndjk(1,1) = 0.0_dp
      dvndjk(2,1) = - vji(3)
      dvndjk(3,1) = vji(2)
!
      dvndjk(1,2) = vji(3)
      dvndjk(2,2) = 0.0_dp
      dvndjk(3,2) = - vji(1)
!
      dvndjk(1,3) = - vji(2)
      dvndjk(2,3) = vji(1)
      dvndjk(3,3) = 0.0_dp
!
      if (lgrad2) then
        d2vndijdjk(1:3,1:3,1:3) = 0.0_dp
!
        d2vndijdjk(3,2,1) = - 1.0_dp
        d2vndijdjk(2,3,1) = 1.0_dp
!
        d2vndijdjk(3,1,2) = 1.0_dp
        d2vndijdjk(1,3,2) = - 1.0_dp
!
        d2vndijdjk(2,1,3) = - 1.0_dp
        d2vndijdjk(1,2,3) = 1.0_dp
      endif
!
!  Compute derivatives of phi w.r.t. vector components
!
      dphidrnv = 1.0_dp/sqrt(1.0_dp - rnv**2)
      rnvrrn2 = rnv/rn2
!
!  Form terms useful for Cartesian derivatives
!
      do ix = 1,3
        dvnvildij(ix)  = dvndij(ix,1)*vil(1) + dvndij(ix,2)*vil(2) + dvndij(ix,3)*vil(3)
        dvnvildjk(ix)  = dvndjk(ix,1)*vil(1) + dvndjk(ix,2)*vil(2) + dvndjk(ix,3)*vil(3)
        dvnvndij(ix)   = dvndij(ix,1)*vn(1) + dvndij(ix,2)*vn(2) + dvndij(ix,3)*vn(3)
        dvnvndjk(ix)   = dvndjk(ix,1)*vn(1) + dvndjk(ix,2)*vn(2) + dvndjk(ix,3)*vn(3)
      enddo
!
!  Derivatives of rnv:
!
!  i-j
!
      drnvdx(1,1) = dvnvildij(1)*rrnil - rnvrrn2*dvnvndij(1)
      drnvdx(2,1) = dvnvildij(2)*rrnil - rnvrrn2*dvnvndij(2)
      drnvdx(3,1) = dvnvildij(3)*rrnil - rnvrrn2*dvnvndij(3)
!
!  j-k
!
      drnvdx(1,2) = dvnvildjk(1)*rrnil - rnvrrn2*dvnvndjk(1)
      drnvdx(2,2) = dvnvildjk(2)*rrnil - rnvrrn2*dvnvndjk(2)
      drnvdx(3,2) = dvnvildjk(3)*rrnil - rnvrrn2*dvnvndjk(3)
!
!  i-l
!
      rnvrril2 = rnv/ril**2
      drnvdx(1,3) = (rrnil*vn(1) - rnvrril2*vil(1))
      drnvdx(2,3) = (rrnil*vn(2) - rnvrril2*vil(2))
      drnvdx(3,3) = (rrnil*vn(3) - rnvrril2*vil(3))
!
      do ii = 1,3
        phi1dx(1:3,ii) = dphidrnv*drnvdx(1:3,ii)
      enddo
!
!  Second derivatives for Cartesian coordinates
!
      if (lgrad2) then
        d2phidrnv2 = rnv/(1.0_dp - rnv**2)**1.5_dp
!
!  Second derivatives of rnv
!
!  i-j/i-j
!
        do ix = 1,3
          do jx = 1,3
            d2rnvdx2(jx,ix,1) = - rrnil*dvnvildij(ix)*dvnvndij(jx)/rn2 &
                                - rrnil*dvnvildij(jx)*dvnvndij(ix)/rn2 &
                                + 3.0_dp*rnvrrn2*dvnvndij(ix)*dvnvndij(jx)/rn2  &
                                - rnvrrn2*(dvndij(ix,1)*dvndij(jx,1) + dvndij(ix,2)*dvndij(jx,2) + dvndij(ix,3)*dvndij(jx,3))
          enddo
        enddo
!
!  j-k/i-j
!
        do ix = 1,3
          do jx = 1,3
            d2rnvdx2(jx,ix,2) = rrnil*(d2vndijdjk(jx,ix,1)*vil(1) + d2vndijdjk(jx,ix,2)*vil(2) + d2vndijdjk(jx,ix,3)*vil(3)) &
                                - rrnil*dvnvildij(ix)*dvnvndjk(jx)/rn2  &
                                - rrnil*dvnvildjk(jx)*dvnvndij(ix)/rn2  &
                                + 3.0_dp*rnvrrn2*dvnvndij(ix)*dvnvndjk(jx)/rn2 &
                                - rnvrrn2*(dvndij(ix,1)*dvndjk(jx,1) + dvndij(ix,2)*dvndjk(jx,2) + dvndij(ix,3)*dvndjk(jx,3) &
                                + d2vndijdjk(jx,ix,1)*vn(1) + d2vndijdjk(jx,ix,2)*vn(2) + d2vndijdjk(jx,ix,3)*vn(3))
          enddo
        enddo
!
!  i-l/i-j
!
        do ix = 1,3
          do jx = 1,3
            d2rnvdx2(jx,ix,3) = rrnil*dvndij(ix,jx) &
                                - rrnil*dvnvildij(ix)*vil(jx)/ril2   &
                                - rrnil*vn(jx)*dvnvndij(ix)/rn2      &
                                + rnv*dvnvndij(ix)*vil(jx)/(rn2*ril2)
          enddo
        enddo
!
!  j-k/j-k
!
        do ix = 1,3
          do jx = 1,3
            d2rnvdx2(jx,ix,4) = - rrnil*dvnvildjk(ix)*dvnvndjk(jx)/rn2   &
                                - rrnil*dvnvildjk(jx)*dvnvndjk(ix)/rn2   &
                                + 3.0_dp*rnvrrn2*dvnvndjk(ix)*dvnvndjk(jx)/rn2   &
                                - rnvrrn2*(dvndjk(ix,1)*dvndjk(jx,1) + dvndjk(ix,2)*dvndjk(jx,2) + dvndjk(ix,3)*dvndjk(jx,3))
          enddo
        enddo
!
!  i-l/j-k
!
        do ix = 1,3
          do jx = 1,3
            d2rnvdx2(jx,ix,5) = rrnil*dvndjk(ix,jx) &
                                - rrnil*dvnvildjk(ix)*vil(jx)/ril2   &
                                - rrnil*vn(jx)*dvnvndjk(ix)/rn2      &
                                + rnv*dvnvndjk(ix)*vil(jx)/(rn2*ril2)
          enddo
        enddo
!
!  i-l/i-l
!
        do ix = 1,3
          do jx = 1,3
            d2rnvdx2(jx,ix,6) = - rrnil*(vn(ix)*vil(jx) + vn(jx)*vil(ix))/ril2 &
                                + 3.0_dp*rnvrril2*vil(ix)*vil(jx)/ril2
            if (ix.eq.jx) then
              d2rnvdx2(jx,ix,6) = d2rnvdx2(jx,ix,6) - rnvrril2
            endif
          enddo
        enddo
!
!  Second derivatives of phi
!
        ind = 0
        do ii = 1,3
          do jj = ii,3
            ind = ind + 1
            do ix = 1,3
              do jx = 1,3
                phi2dx(jx,ix,ind) = d2phidrnv2*drnvdx(jx,jj)*drnvdx(ix,ii) + dphidrnv*d2rnvdx2(jx,ix,ind)
              enddo
            enddo
          enddo
        enddo
      endif
!
      if (lstr) then
!
!  Strain derivatives of vn components
!
        call cartstrterm(ndim,vji(1),vji(2),vji(3),0.0_dp,0.0_dp,0.0_dp,dvjids,d2vjidsdx,d2vjids2,lgrad2)
        call cartstrterm(ndim,vjk(1),vjk(2),vjk(3),0.0_dp,0.0_dp,0.0_dp,dvjkds,d2vjkdsdx,d2vjkds2,lgrad2)
!
!  Switch sign due to vji -> vij
!
        if (lgrad2) d2vjidsdx = - d2vjidsdx
!
        do ks = 1,nstrains
          dvnds(ks,1) = dvjids(ks,2)*vjk(3) - dvjids(ks,3)*vjk(2) + vji(2)*dvjkds(ks,3) - vji(3)*dvjkds(ks,2)
          dvnds(ks,2) = dvjids(ks,3)*vjk(1) - dvjids(ks,1)*vjk(3) + vji(3)*dvjkds(ks,1) - vji(1)*dvjkds(ks,3)
          dvnds(ks,3) = dvjids(ks,1)*vjk(2) - dvjids(ks,2)*vjk(1) + vji(1)*dvjkds(ks,2) - vji(2)*dvjkds(ks,1)
        enddo
!
!  Strain derivatives of vil components
!
        call cartstrterm(ndim,vil(1),vil(2),vil(3),0.0_dp,0.0_dp,0.0_dp,dvilds,d2vildsdx,d2vilds2,lgrad2)
!
!  Use terms for forming strain derivatives
!
        do ks = 1,nstrains
          dvilvilds(ks) = dvilds(ks,1)*vil(1) + dvilds(ks,2)*vil(2) + dvilds(ks,3)*vil(3)
          dvnvilds(ks)  = dvnds(ks,1)*vil(1) + dvnds(ks,2)*vil(2) + dvnds(ks,3)*vil(3)
          dvilvnds(ks)  = vn(1)*dvilds(ks,1) + vn(2)*dvilds(ks,2) + vn(3)*dvilds(ks,3)
          dvnvnds(ks)   = dvnds(ks,1)*vn(1) + dvnds(ks,2)*vn(2) + dvnds(ks,3)*vn(3)
        enddo
!
!  Strain derivatives of rnv
!
        do ks = 1,nstrains
          drnvds(ks) =  rrnil*dvnvilds(ks) - rnvrrn2*dvnvnds(ks) + rrnil*dvilvnds(ks) - rnvrril2*dvilvilds(ks)
        enddo
!
        do ks = 1,nstrains
          phi1ds(ks) = dphidrnv*drnvds(ks)
        enddo
!
        if (lgrad2) then
!
!  Strain-strain second derivatives
!
          do ks = 1,nstrains
            do kt = 1,nstrains
              d2vnds2(kt,ks,1) = d2vjids2(kt,ks,2)*vjk(3) - d2vjids2(kt,ks,3)*vjk(2) &
                                 + vji(2)*d2vjkds2(kt,ks,3) - vji(3)*d2vjkds2(kt,ks,2) &
                                 + dvjids(kt,2)*dvjkds(ks,3) - dvjids(kt,3)*dvjkds(ks,2) &
                                 + dvjids(ks,2)*dvjkds(kt,3) - dvjids(ks,3)*dvjkds(kt,2)
              d2vnds2(kt,ks,2) = d2vjids2(kt,ks,3)*vjk(1) - d2vjids2(kt,ks,1)*vjk(3) &
                                 + vji(3)*d2vjkds2(kt,ks,1) - vji(1)*d2vjkds2(kt,ks,3) &
                                 + dvjids(kt,3)*dvjkds(ks,1) - dvjids(kt,1)*dvjkds(ks,3) &
                                 + dvjids(ks,3)*dvjkds(kt,1) - dvjids(ks,1)*dvjkds(kt,3)
              d2vnds2(kt,ks,3) = d2vjids2(kt,ks,1)*vjk(2) - d2vjids2(kt,ks,2)*vjk(1) &
                                 + vji(1)*d2vjkds2(kt,ks,2) - vji(2)*d2vjkds2(kt,ks,1) &
                                 + dvjids(kt,1)*dvjkds(ks,2) - dvjids(kt,2)*dvjkds(ks,1) &
                                 + dvjids(ks,1)*dvjkds(kt,2) - dvjids(ks,2)*dvjkds(kt,1)
            enddo
          enddo
!
          do ks = 1,nstrains
            do kt = 1,nstrains
              d2rnvds2(kt,ks) = rrnil*(d2vnds2(kt,ks,1)*vil(1) + d2vnds2(kt,ks,2)*vil(2) + d2vnds2(kt,ks,3)*vil(3) + &
                                       vn(1)*d2vilds2(kt,ks,1) + vn(2)*d2vilds2(kt,ks,2) + vn(3)*d2vilds2(kt,ks,3) + &
                                       dvnds(kt,1)*dvilds(ks,1) + dvnds(kt,2)*dvilds(ks,2) + dvnds(kt,3)*dvilds(ks,3) + &
                                       dvnds(ks,1)*dvilds(kt,1) + dvnds(ks,2)*dvilds(kt,2) + dvnds(ks,3)*dvilds(kt,3)) &
                              - rnvrrn2*(vn(1)*d2vnds2(kt,ks,1) + vn(2)*d2vnds2(kt,ks,2) + vn(3)*d2vnds2(kt,ks,3) + &
                                     dvnds(kt,1)*dvnds(ks,1) + dvnds(kt,2)*dvnds(ks,2) + dvnds(kt,3)*dvnds(ks,3)) &
                              - rnvrril2*(vil(1)*d2vilds2(kt,ks,1) + vil(2)*d2vilds2(kt,ks,2) + vil(3)*d2vilds2(kt,ks,3) + &
                                     dvilds(kt,1)*dvilds(ks,1) + dvilds(kt,2)*dvilds(ks,2) + dvilds(kt,3)*dvilds(ks,3)) &
                              - rrnil*(((dvilvnds(kt) + dvnvilds(kt))*dvnvnds(ks) + &
                                        (dvilvnds(ks) + dvnvilds(ks))*dvnvnds(kt))/rn2 + &
                                       ((dvilvnds(kt) + dvnvilds(kt))*dvilvilds(ks) + &
                                        (dvilvnds(ks) + dvnvilds(ks))*dvilvilds(kt))/ril2) &
                              + rnvrrn2*(dvnvnds(kt)*dvilvilds(ks) + dvnvnds(ks)*dvilvilds(kt))/ril2 &
                              + 3.0_dp*(rnvrrn2*dvnvnds(kt)*dvnvnds(ks)/rn2 + rnvrril2*dvilvilds(kt)*dvilvilds(ks)/ril2)
            enddo
          enddo
!
          do ks = 1,nstrains
            do kt = 1,nstrains
              phi2ds(kt,ks) = dphidrnv*d2rnvds2(kt,ks) + d2phidrnv2*drnvds(kt)*drnvds(ks)
            enddo
          enddo
!
          do ix = 1,3
            do ks = 1,nstrains
              d2vndsdij(ks,ix,1) = d2vjidsdx(ks,ix,2)*vjk(3) - d2vjidsdx(ks,ix,3)*vjk(2)
              d2vndsdij(ks,ix,2) = d2vjidsdx(ks,ix,3)*vjk(1) - d2vjidsdx(ks,ix,1)*vjk(3)
              d2vndsdij(ks,ix,3) = d2vjidsdx(ks,ix,1)*vjk(2) - d2vjidsdx(ks,ix,2)*vjk(1)
              d2vndsdjk(ks,ix,1) = vji(2)*d2vjkdsdx(ks,ix,3) - vji(3)*d2vjkdsdx(ks,ix,2)
              d2vndsdjk(ks,ix,2) = vji(3)*d2vjkdsdx(ks,ix,1) - vji(1)*d2vjkdsdx(ks,ix,3)
              d2vndsdjk(ks,ix,3) = vji(1)*d2vjkdsdx(ks,ix,2) - vji(2)*d2vjkdsdx(ks,ix,1)
            enddo
          enddo
!
!  NB: Sign changed due to differentiating w.r.t. vij rather than vji
!
          do ks = 1,nstrains
            d2vndsdij(ks,2,1) = d2vndsdij(ks,2,1) - dvjkds(ks,3)
            d2vndsdij(ks,3,1) = d2vndsdij(ks,3,1) + dvjkds(ks,2)
            d2vndsdij(ks,1,2) = d2vndsdij(ks,1,2) + dvjkds(ks,3)
            d2vndsdij(ks,3,2) = d2vndsdij(ks,3,2) - dvjkds(ks,1)
            d2vndsdij(ks,1,3) = d2vndsdij(ks,1,3) - dvjkds(ks,2)
            d2vndsdij(ks,2,3) = d2vndsdij(ks,2,3) + dvjkds(ks,1)
!
            d2vndsdjk(ks,2,1) = d2vndsdjk(ks,2,1) - dvjids(ks,3)
            d2vndsdjk(ks,3,1) = d2vndsdjk(ks,3,1) + dvjids(ks,2)
            d2vndsdjk(ks,1,2) = d2vndsdjk(ks,1,2) + dvjids(ks,3)
            d2vndsdjk(ks,3,2) = d2vndsdjk(ks,3,2) - dvjids(ks,1)
            d2vndsdjk(ks,1,3) = d2vndsdjk(ks,1,3) - dvjids(ks,2)
            d2vndsdjk(ks,2,3) = d2vndsdjk(ks,2,3) + dvjids(ks,1)
          enddo
!
          do ix = 1,3
            do ks = 1,nstrains
              d2rnvdsdx(ks,ix,1) = rrnil*(d2vndsdij(ks,ix,1)*vil(1) + d2vndsdij(ks,ix,2)*vil(2) + &
                                          d2vndsdij(ks,ix,3)*vil(3) + dvilds(ks,1)*dvndij(ix,1) + &
                                          dvilds(ks,2)*dvndij(ix,2) + dvilds(ks,3)*dvndij(ix,3)) &
                                 - rrnil*(dvilvnds(ks) + dvnvilds(ks))*dvnvndij(ix)/rn2 &
                                 - rrnil*dvnvnds(ks)*dvnvildij(ix)/rn2 &
                                 + 3.0_dp*rnvrrn2*dvnvnds(ks)*dvnvndij(ix)/rn2 &
                                 - rnvrrn2*(dvndij(ix,1)*dvnds(ks,1) + dvndij(ix,2)*dvnds(ks,2) + &
                                            dvndij(ix,3)*dvnds(ks,3) + vn(1)*d2vndsdij(ks,ix,1) + &
                                            vn(2)*d2vndsdij(ks,ix,2) + vn(3)*d2vndsdij(ks,ix,3)) &
                                 - rrnil*dvnvildij(ix)*dvilvilds(ks)/ril2 &
                                 + rnvrrn2*dvilvilds(ks)*dvnvndij(ix)/ril2
!
              d2rnvdsdx(ks,ix,2) = rrnil*(d2vndsdjk(ks,ix,1)*vil(1) + d2vndsdjk(ks,ix,2)*vil(2) + &
                                          d2vndsdjk(ks,ix,3)*vil(3) + dvilds(ks,1)*dvndjk(ix,1) + &
                                          dvilds(ks,2)*dvndjk(ix,2) + dvilds(ks,3)*dvndjk(ix,3)) &
                                 - rrnil*(dvilvnds(ks) + dvnvilds(ks))*dvnvndjk(ix)/rn2 &
                                 - rrnil*dvnvnds(ks)*dvnvildjk(ix)/rn2 &
                                 + 3.0_dp*rnvrrn2*dvnvnds(ks)*dvnvndjk(ix)/rn2 &
                                 - rnvrrn2*(dvndjk(ix,1)*dvnds(ks,1) + dvndjk(ix,2)*dvnds(ks,2) + &
                                            dvndjk(ix,3)*dvnds(ks,3) + vn(1)*d2vndsdjk(ks,ix,1) + &
                                            vn(2)*d2vndsdjk(ks,ix,2) + vn(3)*d2vndsdjk(ks,ix,3)) &
                                 - rrnil*dvnvildjk(ix)*dvilvilds(ks)/ril2 &
                                 + rnvrrn2*dvilvilds(ks)*dvnvndjk(ix)/ril2
!
              d2rnvdsdx(ks,ix,3) = rrnil*(d2vildsdx(ks,ix,1)*vn(1) + d2vildsdx(ks,ix,2)*vn(2) + &
                                          d2vildsdx(ks,ix,3)*vn(3) + dvnds(ks,ix)) &
                                 - rrnil*(dvilvnds(ks) + dvnvilds(ks))*vil(ix)/ril2 &
                                 - rrnil*vn(ix)*dvnvnds(ks)/rn2 &
                                 + rnvrrn2*dvnvnds(ks)*vil(ix)/ril2 &
                                 - rnvrril2*(vil(1)*d2vildsdx(ks,ix,1) + vil(2)*d2vildsdx(ks,ix,2) + &
                                             vil(3)*d2vildsdx(ks,ix,3) + dvilds(ks,ix)) &
                                 - rrnil*dvilvilds(ks)*vn(ix)/ril2 &
                                 + 3.0_dp*rnvrril2*dvilvilds(ks)*vil(ix)/ril2 
            enddo
          enddo
!
          do ix = 1,3
            do ks = 1,nstrains
              phi2dsdx(ks,ix,1:3) = dphidrnv*d2rnvdsdx(ks,ix,1:3) + d2phidrnv2*drnvds(ks)*drnvdx(ix,1:3)
            enddo
          enddo
        endif
      endif
    endif
  endif

  end subroutine gfnff_inversion

  subroutine gfnff_inversion_phon(numat,nnbr,maxnbr,rnbr,xnbr,ynbr,znbr,i,n1,n2,n3,phi,phi1dx,phi2dx)
!
!  Compute improper torsion angle: Phonon version
!  NB: GFNFF uses the angle between the cross product of 2 vectors and a third one
!      which is not the conventional torsion angle.
!  The derivatives are returned in Cartesian form for the vectors ij, jk and il
!
!  Uses vectors from the neighbour list to ensure correct images of atoms
!
!  11/20 Sign change for j->i to i->j moved to dvndij (instead of dvndji)
!   3/21 Phonon version created from gfnff_inversion
!
!  Julian Gale, Curtin University, November 2021
!
  use datatypes
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: n1
  integer(i4), intent(in)  :: n2
  integer(i4), intent(in)  :: n3
  integer(i4), intent(in)  :: nnbr(numat)
  integer(i4), intent(in)  :: maxnbr
  real(dp),    intent(in)  :: rnbr(maxnbr,numat)
  real(dp),    intent(in)  :: xnbr(maxnbr,numat)
  real(dp),    intent(in)  :: ynbr(maxnbr,numat)
  real(dp),    intent(in)  :: znbr(maxnbr,numat)
  real(dp),    intent(out) :: phi
  real(dp),    intent(out) :: phi1dx(3,3)              ! Contains the Cartesian derivatives of phi
  real(dp),    intent(out) :: phi2dx(3,3,6)            ! Contains the Cartesian second derivatives of phi
!
!  Local variables
!
  integer(i4)              :: ii
  integer(i4)              :: ind
  integer(i4)              :: ix
  integer(i4)              :: jj
  integer(i4)              :: jx
  logical                  :: ldistOK
  real(dp)                 :: vil(3)
  real(dp)                 :: vjk(3)
  real(dp)                 :: vji(3)
  real(dp)                 :: vn(3)
  real(dp)                 :: dphidrnv
  real(dp)                 :: d2phidrnv2
  real(dp)                 :: drnvdx(3,3)        ! First derivatives of rnv for 3 vectors
  real(dp)                 :: d2rnvdx2(3,3,6)    ! Second derivatives of rnv for 3 vectors
  real(dp)                 :: dvndij(3,3)        ! Derivatives of vn components (right index) w.r.t. vji components (left index)
  real(dp)                 :: dvndjk(3,3)        ! Derivatives of vn components (right index) w.r.t. vjk components (left index)
  real(dp)                 :: d2vndijdjk(3,3,3)  ! Second derivatives of vn components (right index) w.r.t. vij/vjk components (middle/left index)
  real(dp)                 :: dvnvildij(3)
  real(dp)                 :: dvnvildjk(3)
  real(dp)                 :: dvnvndij(3)
  real(dp)                 :: dvnvndjk(3)
  real(dp)                 :: ril
  real(dp)                 :: ril2
  real(dp)                 :: rn
  real(dp)                 :: rn2
  real(dp)                 :: rrnil
  real(dp)                 :: rnv
  real(dp)                 :: rnvrril2
  real(dp)                 :: rnvrrn2
!
!  Check validity of n1, n2 and n3
!
  if (n1.gt.nnbr(i).or.n2.gt.nnbr(i).or.n3.gt.nnbr(i)) then
    call outerror('invalid arguments for bonds to i in gfnff_inversion_phon',0_i4)
    call stopnow('gfnff_inversion_phon')
  endif
!
!  j -> i vector
!
  vji(1) = - xnbr(n1,i)
  vji(2) = - ynbr(n1,i)
  vji(3) = - znbr(n1,i)
!
!  j -> k vector
!
  vjk(1) = xnbr(n2,i) - xnbr(n1,i)
  vjk(2) = ynbr(n2,i) - ynbr(n1,i)
  vjk(3) = znbr(n2,i) - znbr(n1,i)
!
!  i -> l vector
!
  vil(1) = xnbr(n3,i)
  vil(2) = ynbr(n3,i)
  vil(3) = znbr(n3,i)
!
!  Compute cross product of rji and rjk
!
  vn(1) = vji(2)*vjk(3) - vji(3)*vjk(2)
  vn(2) = vji(3)*vjk(1) - vji(1)*vjk(3)
  vn(3) = vji(1)*vjk(2) - vji(2)*vjk(1)
!
  rn2 = vn(1)**2 + vn(2)**2 + vn(3)**2
  rn = sqrt(rn2)
!
  ril = rnbr(n3,i)
  ril2 = ril**2
!
  ldistOK = (rn*ril.gt.1.0d-12)
  if (ldistOK) then
    rrnil = 1.0_dp/(rn*ril)
    rnv = (vn(1)*vil(1) + vn(2)*vil(2) + vn(3)*vil(3))*rrnil
  else
    rnv = 0.0_dp
  endif
!
  phi = asin(rnv)
!
  phi1dx(1:3,1:3) = 0.0_dp
  phi2dx(1:3,1:3,1:6) = 0.0_dp
  if (ldistOK) then
!
!  Compute derivatives of vectors for cross products
!
!  NB: Sign changed due to differentiating w.r.t. vij rather than vji
!
    dvndij(1,1) = 0.0_dp
    dvndij(2,1) = - vjk(3)
    dvndij(3,1) = vjk(2)
!
    dvndij(1,2) = vjk(3)
    dvndij(2,2) = 0.0_dp
    dvndij(3,2) = - vjk(1)
!
    dvndij(1,3) = - vjk(2)
    dvndij(2,3) = vjk(1)
    dvndij(3,3) = 0.0_dp
!
    dvndjk(1,1) = 0.0_dp
    dvndjk(2,1) = - vji(3)
    dvndjk(3,1) = vji(2)
!
    dvndjk(1,2) = vji(3)
    dvndjk(2,2) = 0.0_dp
    dvndjk(3,2) = - vji(1)
!
    dvndjk(1,3) = - vji(2)
    dvndjk(2,3) = vji(1)
    dvndjk(3,3) = 0.0_dp
!
    d2vndijdjk(1:3,1:3,1:3) = 0.0_dp
!
    d2vndijdjk(3,2,1) = - 1.0_dp
    d2vndijdjk(2,3,1) = 1.0_dp
!
    d2vndijdjk(3,1,2) = 1.0_dp
    d2vndijdjk(1,3,2) = - 1.0_dp
!
    d2vndijdjk(2,1,3) = - 1.0_dp
    d2vndijdjk(1,2,3) = 1.0_dp
!
!  Compute derivatives of phi w.r.t. vector components
!
    dphidrnv = 1.0_dp/sqrt(1.0_dp - rnv**2)
    rnvrrn2 = rnv/rn2
!
!  Form terms useful for Cartesian derivatives
!
    do ix = 1,3
      dvnvildij(ix)  = dvndij(ix,1)*vil(1) + dvndij(ix,2)*vil(2) + dvndij(ix,3)*vil(3)
      dvnvildjk(ix)  = dvndjk(ix,1)*vil(1) + dvndjk(ix,2)*vil(2) + dvndjk(ix,3)*vil(3)
      dvnvndij(ix)   = dvndij(ix,1)*vn(1) + dvndij(ix,2)*vn(2) + dvndij(ix,3)*vn(3)
      dvnvndjk(ix)   = dvndjk(ix,1)*vn(1) + dvndjk(ix,2)*vn(2) + dvndjk(ix,3)*vn(3)
    enddo
!
!  Derivatives of rnv:
!
!  i-j
!
    drnvdx(1,1) = dvnvildij(1)*rrnil - rnvrrn2*dvnvndij(1)
    drnvdx(2,1) = dvnvildij(2)*rrnil - rnvrrn2*dvnvndij(2)
    drnvdx(3,1) = dvnvildij(3)*rrnil - rnvrrn2*dvnvndij(3)
!
!  j-k
!
    drnvdx(1,2) = dvnvildjk(1)*rrnil - rnvrrn2*dvnvndjk(1)
    drnvdx(2,2) = dvnvildjk(2)*rrnil - rnvrrn2*dvnvndjk(2)
    drnvdx(3,2) = dvnvildjk(3)*rrnil - rnvrrn2*dvnvndjk(3)
!
!  i-l
!
    rnvrril2 = rnv/ril**2
    drnvdx(1,3) = (rrnil*vn(1) - rnvrril2*vil(1))
    drnvdx(2,3) = (rrnil*vn(2) - rnvrril2*vil(2))
    drnvdx(3,3) = (rrnil*vn(3) - rnvrril2*vil(3))
!
    do ii = 1,3
      phi1dx(1:3,ii) = dphidrnv*drnvdx(1:3,ii)
    enddo
!
!  Second derivatives for Cartesian coordinates
!
    d2phidrnv2 = rnv/(1.0_dp - rnv**2)**1.5_dp
!
!  Second derivatives of rnv
!
!  i-j/i-j
!
    do ix = 1,3
      do jx = 1,3
        d2rnvdx2(jx,ix,1) = - rrnil*dvnvildij(ix)*dvnvndij(jx)/rn2 &
                            - rrnil*dvnvildij(jx)*dvnvndij(ix)/rn2 &
                            + 3.0_dp*rnvrrn2*dvnvndij(ix)*dvnvndij(jx)/rn2  &
                            - rnvrrn2*(dvndij(ix,1)*dvndij(jx,1) + dvndij(ix,2)*dvndij(jx,2) + dvndij(ix,3)*dvndij(jx,3))
      enddo
    enddo
!
!  j-k/i-j
!
    do ix = 1,3
      do jx = 1,3
        d2rnvdx2(jx,ix,2) = rrnil*(d2vndijdjk(jx,ix,1)*vil(1) + d2vndijdjk(jx,ix,2)*vil(2) + d2vndijdjk(jx,ix,3)*vil(3)) &
                            - rrnil*dvnvildij(ix)*dvnvndjk(jx)/rn2  &
                            - rrnil*dvnvildjk(jx)*dvnvndij(ix)/rn2  &
                            + 3.0_dp*rnvrrn2*dvnvndij(ix)*dvnvndjk(jx)/rn2 &
                            - rnvrrn2*(dvndij(ix,1)*dvndjk(jx,1) + dvndij(ix,2)*dvndjk(jx,2) + dvndij(ix,3)*dvndjk(jx,3) &
                            + d2vndijdjk(jx,ix,1)*vn(1) + d2vndijdjk(jx,ix,2)*vn(2) + d2vndijdjk(jx,ix,3)*vn(3))
      enddo
    enddo
!
!  i-l/i-j
!
    do ix = 1,3
      do jx = 1,3
        d2rnvdx2(jx,ix,3) = rrnil*dvndij(ix,jx) &
                            - rrnil*dvnvildij(ix)*vil(jx)/ril2   &
                            - rrnil*vn(jx)*dvnvndij(ix)/rn2      &
                            + rnv*dvnvndij(ix)*vil(jx)/(rn2*ril2)
      enddo
    enddo
!
!  j-k/j-k
!
    do ix = 1,3
      do jx = 1,3
        d2rnvdx2(jx,ix,4) = - rrnil*dvnvildjk(ix)*dvnvndjk(jx)/rn2   &
                            - rrnil*dvnvildjk(jx)*dvnvndjk(ix)/rn2   &
                            + 3.0_dp*rnvrrn2*dvnvndjk(ix)*dvnvndjk(jx)/rn2   &
                            - rnvrrn2*(dvndjk(ix,1)*dvndjk(jx,1) + dvndjk(ix,2)*dvndjk(jx,2) + dvndjk(ix,3)*dvndjk(jx,3))
      enddo
    enddo
!
!  i-l/j-k
!
    do ix = 1,3
      do jx = 1,3
        d2rnvdx2(jx,ix,5) = rrnil*dvndjk(ix,jx) &
                            - rrnil*dvnvildjk(ix)*vil(jx)/ril2   &
                            - rrnil*vn(jx)*dvnvndjk(ix)/rn2      &
                            + rnv*dvnvndjk(ix)*vil(jx)/(rn2*ril2)
      enddo
    enddo
!
!  i-l/i-l
!
    do ix = 1,3
      do jx = 1,3
        d2rnvdx2(jx,ix,6) = - rrnil*(vn(ix)*vil(jx) + vn(jx)*vil(ix))/ril2 &
                            + 3.0_dp*rnvrril2*vil(ix)*vil(jx)/ril2
        if (ix.eq.jx) then
          d2rnvdx2(jx,ix,6) = d2rnvdx2(jx,ix,6) - rnvrril2
        endif
      enddo
    enddo
!
!  Second derivatives of phi
!
    ind = 0
    do ii = 1,3
      do jj = ii,3
        ind = ind + 1
        do ix = 1,3
          do jx = 1,3
            phi2dx(jx,ix,ind) = d2phidrnv2*drnvdx(jx,jj)*drnvdx(ix,ii) + dphidrnv*d2rnvdx2(jx,ix,ind)
          enddo
        enddo
      enddo
    enddo
  endif

  end subroutine gfnff_inversion_phon
