!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  HB calculation without lists  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gfnff_hbpd(nnbr,maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xkv,ykv,zkv)
!
!  Compute the hydrogen bond energy without the use of list based methods
!  so that periodic boundary conditions can be handled correctly.
!  Distribributed memory parallel version.
!
!   4/22 Created from gfnff_hbd
!
!  Julian Gale, Curtin University, April 2022
!
  use datatypes
  use current
  use derivatives
  use gulp_gfnff
  use m_gfnff_nbr3
  use m_gfnff_pairs
  use parallel
  implicit none
!
!  Passed variables
!
  integer(i4),                         intent(in)  :: nnbr(numat)
  integer(i4),                         intent(in)  :: maxnbr
  integer(i4),                         intent(in)  :: nbrno(maxnbr,numat)
  real(dp),                            intent(in)  :: rbnbr(maxnbr,numat)
  real(dp),                            intent(in)  :: xbnbr(maxnbr,numat)
  real(dp),                            intent(in)  :: ybnbr(maxnbr,numat)
  real(dp),                            intent(in)  :: zbnbr(maxnbr,numat)
  real(dp),                            intent(in)  :: xkv
  real(dp),                            intent(in)  :: ykv
  real(dp),                            intent(in)  :: zkv
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ii
  integer(i4)                                      :: iloc
  integer(i4)                                      :: ilast
  integer(i4)                                      :: indi
  integer(i4)                                      :: indj
  integer(i4)                                      :: indk
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: ixf
  integer(i4)                                      :: iyf
  integer(i4)                                      :: izf
  integer(i4)                                      :: j
  integer(i4)                                      :: jj
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
  integer(i4)                                      :: nh
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nk
  integer(i4)                                      :: np
  integer(i4)                                      :: np2
  logical                                          :: lbonded3
  logical                                          :: lijbonded
  logical                                          :: likbonded
  logical                                          :: ljkbonded
  logical                                          :: lilocal
  logical                                          :: ljlocal
  logical                                          :: lklocal
  logical                                          :: loneneighC
  logical                                          :: loneneighN
  real(dp)                                         :: acid
  real(dp)                                         :: base
  real(dp)                                         :: caa
  real(dp)                                         :: cbb
  real(dp)                                         :: const
  real(dp)                                         :: dadrik
  real(dp)                                         :: dadrjk
  real(dp)                                         :: dabdrik
  real(dp)                                         :: dabdrjk
  real(dp)                                         :: dbdrik
  real(dp)                                         :: dbdrjk
  real(dp)                                         :: d2adrik2
  real(dp)                                         :: d2adrikdrjk
  real(dp)                                         :: d2adrjk2
  real(dp)                                         :: d2abdrik2
  real(dp)                                         :: d2abdrikdrjk
  real(dp)                                         :: d2abdrjk2
  real(dp)                                         :: d2bdrik2
  real(dp)                                         :: d2bdrikdrjk
  real(dp)                                         :: d2bdrjk2
  real(dp)                                         :: damp
  real(dp)                                         :: dampl
  real(dp)                                         :: ddampl
  real(dp)                                         :: d2dampl
  real(dp)                                         :: dampo
  real(dp)                                         :: ddampodrij
  real(dp)                                         :: ddampodrik
  real(dp)                                         :: ddampodrjk
  real(dp)                                         :: d2dampo(6)
  real(dp)                                         :: damps
  real(dp)                                         :: ddamps
  real(dp)                                         :: d2damps
  real(dp)                                         :: denom
  real(dp)                                         :: deijkdrij
  real(dp)                                         :: deijkdrik
  real(dp)                                         :: deijkdrjk
  real(dp)                                         :: d2trm
  real(dp)                                         :: expo
  real(dp)                                         :: dexpodrij
  real(dp)                                         :: dexpodrik
  real(dp)                                         :: dexpodrjk
  real(dp)                                         :: d2expo(6)
  real(dp)                                         :: eijk
  real(dp)                                         :: fct
  real(dp)                                         :: qhba
  real(dp)                                         :: qhbb
  real(dp)                                         :: qhbx
  real(dp)                                         :: drdampdrij
  real(dp)                                         :: d2rdampdrij2
  real(dp)                                         :: drdampdrjk
  real(dp)                                         :: d2rdampdrjk2
  real(dp)                                         :: d2r2dx23(6,3)
  real(dp)                                         :: rij
  real(dp)                                         :: rik
  real(dp)                                         :: rjk
  real(dp)                                         :: r2ij
  real(dp)                                         :: r2ik
  real(dp)                                         :: r2jk
  real(dp)                                         :: r4ik
  real(dp)                                         :: r4jk
  real(dp)                                         :: radij
  real(dp)                                         :: rsum
  real(dp)                                         :: ratio1
  real(dp)                                         :: dratio1drij
  real(dp)                                         :: dratio1drjk
  real(dp)                                         :: d2ratio1drij2
  real(dp)                                         :: d2ratio1drjk2
  real(dp)                                         :: ratio2
  real(dp)                                         :: ratio3
  real(dp)                                         :: dratio3drij
  real(dp)                                         :: dratio3drjk
  real(dp)                                         :: d2ratio3drij2
  real(dp)                                         :: d2ratio3drjk2
  real(dp)                                         :: rdamp
  real(dp)                                         :: shortcut
  real(dp)                                         :: xij
  real(dp)                                         :: yij
  real(dp)                                         :: zij
  real(dp)                                         :: xik
  real(dp)                                         :: yik
  real(dp)                                         :: zik
  real(dp)                                         :: xjk
  real(dp)                                         :: yjk
  real(dp)                                         :: zjk
!
  nhb1 = 0
  nhb2 = 0
!
  call gfnff_setpaircell(gfnff_hbthr1,hb1_paircell)
  call gfnff_setpaircell(gfnff_hbthr2,hb2_paircell)
!
!  Loop over possible A/B atoms
!
  do ii = 1,n_gfnff_hb_AB
    i = n_gfnff_hb_ABptr(ii)
    iloc = atom2local(i)
    lilocal = (iloc.ne.0)
!
!  A charge scaled term
!
    qhba = gfnff_hb_ABq(i)
!
!  Set up bonded neighbour info for i
!
    call gfnff_get_n3atoms(numat,i,1_i4)
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
    do jj = 1,ii
      j = n_gfnff_hb_ABptr(jj)
!
!  Check on A and B as a pair
!
      if (gfnff_hb_base(i)*gfnff_hb_acid(j).lt.1.d-6.and.gfnff_hb_acid(i)*gfnff_hb_base(j).lt.1.d-6) cycle
!
      jloc = atom2local(j)
      ljlocal = (jloc.ne.0)
!
!  Set parameters
!
      radij = (gfnff_rad(nat(i)) + gfnff_rad(nat(j)))
!
!  B charge scaled term
!
      qhbb = gfnff_hb_ABq(j)
!
!  Correction for i = j
!
      if (i.eq.j) then
        fct = 0.5_dp
      else
        fct = 1.0_dp
      endif
      fct = fct*gfnff_autoangs**3
!
!  Find interactions for this pair
!
      call gfnff_getpairs(i,j,gfnff_hbthr1,hb1_paircell,hb1_pairs)
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
!  Loop over valid interactions (self term is excluded by getpairs)
!
      do np = 1,hb1_pairs%npair
        r2ij = hb1_pairs%r2pair(np)
        rij = sqrt(r2ij)
        xij = hb1_pairs%xpair(np)
        yij = hb1_pairs%ypair(np)
        zij = hb1_pairs%zpair(np)
!
!  Are i and j bonded?
!
        call gfnff_bond_check(j,xij,yij,zij,1_i4,lijbonded)
!
        if (.not.lijbonded) then
!
!  Search for atoms bonded to A and B that are H for HB
!
          do ni = 1,nnbr(i)
            k = nbrno(ni,i)
            nh = n_gfnff_hb_Hrptr(k)
            if (nh.eq.0) cycle
            nhb2 = nhb2 + 1
!!!!!!!!!!!!!!!!!!!!
!  Compute energy  !
!!!!!!!!!!!!!!!!!!!!
!
!  Carbonyl case R-C=O...H_A
!
            if (nnbr(j).eq.1) then
              loneneighC = (nat(nbrno(1,j)).eq.6)
              loneneighN = (nat(nbrno(1,j)).eq.7)
            else
              loneneighC = .false.
              loneneighN = .false.
            endif
            if (nat(j).eq.8.and.loneneighC) then
              const = gfnff_hb_acid(i)*qhba*gfnff_hb_base(j)*qhbb*gfnff_hb_scale_coh
              call gfnff_hb2_eg3pd(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                   maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const,radij,xkv,ykv,zkv)
!
!  Nitro case R-N=O...H_A
!
            elseif (nat(j).eq.8.and.loneneighN) then
              const = gfnff_hb_acid(i)*qhba*gfnff_hb_base(j)*qhbb*gfnff_hb_scale_coh
              call gfnff_hb2_eg3pd(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                   maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const,radij,xkv,ykv,zkv)
!
!  N hetero aromatic
!
            elseif (nat(j).eq.7.and.nnbr(j).eq.2) then
              const = gfnff_hb_acid(i)*qhba*gfnff_hb_base(j)*qhbb*gfnff_hb_scale_gen
              call gfnff_hb2_eg2pd_rnr(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                       maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const,radij,xkv,ykv,zkv)
            else
!
!  Default
!
              const = gfnff_hb_acid(i)*qhba*gfnff_hb_base(j)*qhbb*gfnff_hb_scale_gen
              call gfnff_hb2_eg2pd(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                   maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const,radij,xkv,ykv,zkv)
            endif
          enddo
! DEBUG - check that H isn't bonded to image of i!!
          if (j.ne.i) then
            do nj = 1,nnbr(j)
              k = nbrno(nj,j)
              nh = n_gfnff_hb_Hrptr(k)
              if (nh.eq.0) cycle
              nhb2 = nhb2 + 1
!!!!!!!!!!!!!!!!!!!!
!  Compute energy  !
!!!!!!!!!!!!!!!!!!!!
!
!  Carbonyl case R-C=O...H_A
!
              if (nnbr(i).eq.1) then
                loneneighC = (nat(nbrno(1,i)).eq.6)
                loneneighN = (nat(nbrno(1,i)).eq.7)
              else
                loneneighC = .false.
                loneneighN = .false.
              endif
              if (nat(i).eq.8.and.loneneighC) then
                const = gfnff_hb_base(i)*qhba*gfnff_hb_acid(j)*qhbb*gfnff_hb_scale_coh
                call gfnff_hb2_eg3pd(ljlocal,lilocal,j,i,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf,nj,nnbr, &
                                     maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,-xij,-yij,-zij,rij,const,radij,xkv,ykv,zkv)
!
!  Nitro case R-N=O...H_A
!
              elseif (nat(i).eq.8.and.loneneighN) then
                const = gfnff_hb_base(i)*qhba*gfnff_hb_acid(j)*qhbb*gfnff_hb_scale_coh
                call gfnff_hb2_eg3pd(ljlocal,lilocal,j,i,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf,nj,nnbr, &
                                     maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,-xij,-yij,-zij,rij,const,radij,xkv,ykv,zkv)
!
!  N hetero aromatic
!
              elseif (nat(i).eq.7.and.nnbr(i).eq.2) then
                const = gfnff_hb_base(i)*qhba*gfnff_hb_acid(j)*qhbb*gfnff_hb_scale_gen
                call gfnff_hb2_eg2pd_rnr(ljlocal,lilocal,j,i,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf,nj,nnbr, &
                                         maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,-xij,-yij,-zij,rij,const,radij,xkv,ykv,zkv)
              else
!
!  Default
!
                const = gfnff_hb_base(i)*qhba*gfnff_hb_acid(j)*qhbb*gfnff_hb_scale_gen
                call gfnff_hb2_eg2pd(ljlocal,lilocal,j,i,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf,nj,nnbr, &
                                     maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,-xij,-yij,-zij,rij,const,radij,xkv,ykv,zkv)
              endif
            enddo
          endif
        endif
!
!  Now loop over H atoms excluding those found to have a bond to A or B
!
        do nh = 1,n_gfnff_hb_H
          k = n_gfnff_hb_Hptr(nh)
          kloc = atom2local(k)
          lklocal = (kloc.ne.0)
!
!  Are any of the atoms local to this node?
!
          if (.not.lilocal.and..not.ljlocal.and..not.lklocal) cycle
!
!  Find images
!
          call gfnff_gettriads(ndim,i,j,k,r2ij,xij,yij,zij,gfnff_hbthr2,hb2_paircell,hb2_triads)
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
!  Loop over valid interactions (self terms are excluded by gettriads)
!
          do np2 = 1,hb2_triads%ntriad
            xik = hb2_triads%xikpair(np2)
            yik = hb2_triads%yikpair(np2)
            zik = hb2_triads%zikpair(np2)
!
!  Are i and k bonded?
!
            call gfnff_bond_check(k,xik,yik,zik,1_i4,likbonded)
            if (likbonded) cycle
!
            r2jk = hb2_triads%r2jkpair(np2)
            rjk = sqrt(r2jk)
!
!  Are j and k bonded?
!
            ljkbonded = .false.
            do nj = 1,nnbr(j)
              if (k.ne.nbrno(nj,j)) cycle
!
!  Atom number matches, but does the distance?
!
              if (abs(rjk-rbnbr(nj,j)).lt.1.0d-2) then
                ljkbonded = .true.
                exit
              endif
            enddo
            if (ljkbonded) cycle
!
            r2ik = hb2_triads%r2ikpair(np2)
            rik = sqrt(r2ik)
!
            xjk = hb2_triads%xjkpair(np2)
            yjk = hb2_triads%yjkpair(np2)
            zjk = hb2_triads%zjkpair(np2)
!
            nhb1 = nhb1 + 1
!
            rsum = rik + rjk + 1.0d-12
!
!  Out-of-line damping
!
            expo = (gfnff_hb_a_cut/radij)*(rsum/rij - 1.0_dp)
            if (expo.gt.15.0_dp) cycle
            ratio2 = exp(expo)
            dampo = 2.0_dp/(1.0_dp + ratio2)
!
!  Long-range damping
!
            ratio1 = (r2ij/gfnff_hb_long_cut)**gfnff_hb_alp
            dampl = 1.0_dp/(1.0_dp + ratio1)
!
!  Short-range damping
!
            shortcut = gfnff_hb_short_cut*radij
            ratio3 = (shortcut/r2ij)**gfnff_hb_alp
            damps = 1.0_dp/(1.0_dp + ratio3)
!
            damp = damps*dampl
            rdamp = damp/r2ij/rij
!
!  Donor-acceptor term
!
            r4ik = r2ik*r2ik
            r4jk = r2jk*r2jk
            denom = 1.0_dp/(r4ik + r4jk)
!
            caa = qhba*gfnff_hb_base(i)
            cbb = qhbb*gfnff_hb_base(j)
!
            base = (caa*r4ik + cbb*r4jk)*denom
            acid = (gfnff_hb_acid(j)*r4ik + gfnff_hb_acid(i)*r4jk)*denom
!
!  Energy
!
            const = fct*gfnff_hb_ABq(k)
            eijk = - const*acid*base*rdamp*dampo
!
!----------------------
!  First derivatives  |
!----------------------
            dexpodrij = - (gfnff_hb_a_cut/radij)*ratio2*rsum/rij**3
            dexpodrik = (gfnff_hb_a_cut/radij)*ratio2/(rij*rik)
            dexpodrjk = (gfnff_hb_a_cut/radij)*ratio2/(rij*rjk)
!
            ddampodrij = - 0.5_dp*dampo*dampo*dexpodrij
            ddampodrik = - 0.5_dp*dampo*dampo*dexpodrik
            ddampodrjk = - 0.5_dp*dampo*dampo*dexpodrjk
!
!  Damping derivatives
!
            dratio1drij = 2.0_dp*gfnff_hb_alp*ratio1/r2ij
            dratio3drij = - 2.0_dp*gfnff_hb_alp*ratio3/r2ij
            ddampl = - dampl*dampl*dratio1drij
            ddamps = - damps*damps*dratio3drij
            drdampdrij = (ddamps*dampl + damps*ddampl)/r2ij/rij - 3.0_dp*rdamp/r2ij
!
!  Derivatives of acid*base
!
            dadrik = 4.0_dp*denom*(gfnff_hb_acid(j) - acid)*r2ik
            dbdrik = 4.0_dp*denom*(caa - base)*r2ik
            dadrjk = 4.0_dp*denom*(gfnff_hb_acid(i) - acid)*r2jk
            dbdrjk = 4.0_dp*denom*(cbb - base)*r2jk
            dabdrik = dadrik*base + acid*dbdrik
            dabdrjk = dadrjk*base + acid*dbdrjk
!
            deijkdrij = - const*acid*base*ddampodrij*rdamp &         ! (1/r(dE/dr))
                        - const*acid*base*dampo*drdampdrij
            deijkdrik = - const*acid*base*ddampodrik*rdamp &         ! (1/r(dE/dr))
                        - const*dabdrik*rdamp*dampo
            deijkdrjk = - const*acid*base*ddampodrjk*rdamp &         ! (1/r(dE/dr))
                        - const*dabdrjk*rdamp*dampo
!
            d2r2dx23(1,1) = xij*xij
            d2r2dx23(2,1) = yij*yij
            d2r2dx23(3,1) = zij*zij
            d2r2dx23(4,1) = yij*zij
            d2r2dx23(5,1) = xij*zij
            d2r2dx23(6,1) = xij*yij
!
            d2r2dx23(1,2) = xik*xik
            d2r2dx23(2,2) = yik*yik
            d2r2dx23(3,2) = zik*zik
            d2r2dx23(4,2) = yik*zik
            d2r2dx23(5,2) = xik*zik
            d2r2dx23(6,2) = xik*yik
!
            d2r2dx23(1,3) = xjk*xjk
            d2r2dx23(2,3) = yjk*yjk
            d2r2dx23(3,3) = zjk*zjk
            d2r2dx23(4,3) = yjk*zjk
            d2r2dx23(5,3) = xjk*zjk
            d2r2dx23(6,3) = xjk*yjk
!-----------------------
!  Second derivatives  |
!-----------------------
            d2expo(1) = - dexpodrij*(gfnff_hb_a_cut/radij)*rsum/rij**3 - 3.0_dp*dexpodrij/r2ij
            d2expo(2) = dexpodrij*(gfnff_hb_a_cut/radij)/(rij*rik) - dexpodrik/r2ij
            d2expo(3) = dexpodrij*(gfnff_hb_a_cut/radij)/(rij*rjk) - dexpodrjk/r2ij
            d2expo(4) = dexpodrik*(gfnff_hb_a_cut/radij)/(rij*rik) - dexpodrik/r2ik
            d2expo(5) = dexpodrik*(gfnff_hb_a_cut/radij)/(rij*rjk)
            d2expo(6) = dexpodrjk*(gfnff_hb_a_cut/radij)/(rij*rjk) - dexpodrjk/r2jk
!
            d2dampo(1) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrij - d2expo(1))
            d2dampo(2) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrik - d2expo(2))
            d2dampo(3) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrjk - d2expo(3))
            d2dampo(4) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrik - d2expo(4))
            d2dampo(5) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrjk - d2expo(5))
            d2dampo(6) = 0.5_dp*dampo*dampo*(dampo*dexpodrjk*dexpodrjk - d2expo(6))
!
            d2ratio1drij2 = 4.0_dp*gfnff_hb_alp*ratio1*(gfnff_hb_alp - 1.0_dp)/r2ij**2
            d2ratio3drij2 = 4.0_dp*gfnff_hb_alp*ratio3*(gfnff_hb_alp + 1.0_dp)/r2ij**2
            d2dampl = - dampl*dampl*d2ratio1drij2 + 2.0_dp*dampl*dampl*dampl*dratio1drij**2
            d2damps = - damps*damps*d2ratio3drij2 + 2.0_dp*damps*damps*damps*dratio3drij**2
            d2rdampdrij2 = (d2damps*dampl + 2.0_dp*ddamps*ddampl + damps*d2dampl)/r2ij/rij + &
                           15.0_dp*rdamp/r2ij**2 - 6.0_dp*(ddamps*dampl + damps*ddampl)/r2ij/r2ij/rij
!
            d2adrik2 = 8.0_dp*denom*(gfnff_hb_acid(j) - acid)*(1.0_dp - 2.0_dp*denom*r2ik*r2ik) - 4.0_dp*denom*dadrik*r2ik
            d2bdrik2 = 8.0_dp*denom*(caa - base)*(1.0_dp - 2.0_dp*denom*r2ik*r2ik) - 4.0_dp*denom*dbdrik*r2ik
            d2adrjk2 = 8.0_dp*denom*(gfnff_hb_acid(i) - acid)*(1.0_dp - 2.0_dp*denom*r2jk*r2jk) - 4.0_dp*denom*dadrjk*r2jk
            d2bdrjk2 = 8.0_dp*denom*(cbb - base)*(1.0_dp - 2.0_dp*denom*r2jk*r2jk) - 4.0_dp*denom*dbdrjk*r2jk
            d2adrikdrjk = - 16.0_dp*denom*denom*(gfnff_hb_acid(j) - acid)*r2ik*r2jk - 4.0_dp*denom*dadrjk*r2ik
            d2bdrikdrjk = - 16.0_dp*denom*denom*(caa - base)*r2ik*r2jk - 4.0_dp*denom*dbdrjk*r2ik
!
            d2abdrik2 = d2adrik2*base + 2.0_dp*dadrik*dbdrik + acid*d2bdrik2
            d2abdrjk2 = d2adrjk2*base + 2.0_dp*dadrjk*dbdrjk + acid*d2bdrjk2
            d2abdrikdrjk = d2adrikdrjk*base + dadrik*dbdrjk + dadrjk*dbdrik + acid*d2bdrikdrjk
!
!  i-j / i-j
!
            d2trm = - const*acid*base*(dampo*d2rdampdrij2 + 2.0_dp*drdampdrij*ddampodrij + rdamp*d2dampo(1))
            call add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                              deijkdrij,d2trm,d2r2dx23(1,1))
!
!  i-j / i-k
!
            d2trm = - const*(dabdrik*(drdampdrij*dampo + rdamp*ddampodrij) + &
                             acid*base*(drdampdrij*ddampodrik + rdamp*d2dampo(2)))
            call add_drv2_2pd(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                              ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm)
!
!  i-j / j-k
!
            d2trm = - const*(dabdrjk*(drdampdrij*dampo + rdamp*ddampodrij) + &
                             acid*base*(drdampdrij*ddampodrjk + rdamp*d2dampo(3)))
            call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                              jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm)
!
!  i-k / i-k
!
            d2trm = - const*rdamp*(d2abdrik2*dampo + 2.0_dp*dabdrik*ddampodrik + acid*base*d2dampo(4))
            call add_drv2_1pd(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik, &
                              deijkdrik,d2trm,d2r2dx23(1,2))
!
!  i-k / j-k
!
            d2trm = - const*rdamp*(d2abdrikdrjk*dampo + dabdrik*ddampodrjk + dabdrjk*ddampodrik + acid*base*d2dampo(5))
            call add_drv2_2pd(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                              jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm)
!
!  j-k / j-k
!
            d2trm = - const*rdamp*(d2abdrjk2*dampo + 2.0_dp*dabdrjk*ddampodrjk + acid*base*d2dampo(6))
            call add_drv2_1pd(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk, &
                              deijkdrjk,d2trm,d2r2dx23(1,3))
!
          enddo ! End loop over valid images of this hydrogen
        enddo   ! End loop over valid hydrogens
!
      enddo     ! End loop over images
    enddo       ! End loop over j
  enddo         ! End loop over i
!
!  Halogen bonds
!
  nhb3 = 0
  ilast = 0
  do ii = 1,n_gfnff_xb_AB
    i  = n_gfnff_xb_ABptr(1,ii)
    iloc = atom2local(i)
    lilocal = (iloc.ne.0)
!
    j  = n_gfnff_xb_ABptr(2,ii)
    nk = n_gfnff_xb_ABptr(3,ii)   ! X
    k  = nbrno(nk,i)
!
    jloc = atom2local(j)
    kloc = atom2local(k)
    ljlocal = (jloc.ne.0)
    lklocal = (kloc.ne.0)
!
!  Are any of the atoms local to this node?
!
    if (.not.lilocal.and..not.ljlocal.and..not.lklocal) cycle
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
!  Set up 3 bond neighbour info
!
    if (i.ne.ilast) then
      ilast = i
      call gfnff_get_n3atoms(numat,i,3_i4)
    endif
!
    rik = rbnbr(nk,i)
    xik = xbnbr(nk,i)
    yik = ybnbr(nk,i)
    zik = zbnbr(nk,i)
!
!  Set parameters
!
    radij = (gfnff_rad(nat(i)) + gfnff_rad(nat(j)))
!
!  Find interactions for this pair
!
    call gfnff_getpairs(i,j,gfnff_hbthr2,hb2_paircell,hb2_pairs)
!
    if (hb2_pairs%npair.gt.0) then
!
!  Charge-scaling terms
!
      qhbx = gfnff_xb_ABq(1,k)
      qhbb = gfnff_xb_ABq(2,j)
    endif
!
    r2ik = rik*rik
!
!  Loop over valid interactions (self term is excluded by getpairs)
!
    do np = 1,hb2_pairs%npair
      r2ij = hb2_pairs%r2pair(np)
      xij = hb2_pairs%xpair(np)
      yij = hb2_pairs%ypair(np)
      zij = hb2_pairs%zpair(np)
!
!  Check whether j is allowed based on being more than 3 bonds away from i
!
      call gfnff_bond_check(j,xij,yij,zij,3_i4,lbonded3)
      if (lbonded3) cycle
!
      rij = sqrt(r2ij)
!
      xjk = xik - xij
      yjk = yik - yij
      zjk = zik - zij
      r2jk = xjk**2 + yjk**2 + zjk**2
      rjk = sqrt(r2jk)
!
      nhb3 = nhb3 + 1
!
!  Out-of-line damping
!
      expo = gfnff_xb_a_cut*((rik + rjk)/rij - 1.0_dp)
      if (expo.gt.15.0_dp) cycle
      ratio2 = exp(expo)
      dampo = 2.0_dp/(1.0_dp + ratio2)
!
!  Long-range damping
!
      ratio1 = (r2jk/gfnff_xb_long_cut)**gfnff_hb_alp
      dampl = 1.0_dp/(1.0_dp + ratio1)
!
!  Short-range damping
!
      shortcut = gfnff_xb_short_cut*radij
      ratio3 = (shortcut/r2jk)**gfnff_hb_alp
      damps = 1.0_dp/(1.0_dp + ratio3)
!
      damp = damps*dampl
      rdamp = damp/r2jk/rjk
!
      const = qhbb*qhbx*gfnff_xb_scale(nat(k))
!
!  Energy
!
      eijk = - rdamp*dampo*const
!----------------------
!  First derivatives  |
!----------------------
      dexpodrij = - gfnff_xb_a_cut*ratio2*(rik + rjk)/rij**3
      dexpodrik = gfnff_xb_a_cut*ratio2/(rij*rik)
      dexpodrjk = gfnff_xb_a_cut*ratio2/(rij*rjk)
!
      ddampodrij = - 0.5_dp*dampo*dampo*dexpodrij
      ddampodrik = - 0.5_dp*dampo*dampo*dexpodrik
      ddampodrjk = - 0.5_dp*dampo*dampo*dexpodrjk
!
!  Damping derivatives
!
      dratio1drjk = 2.0_dp*gfnff_hb_alp*ratio1/r2jk
      dratio3drjk = - 2.0_dp*gfnff_hb_alp*ratio3/r2jk
      ddampl = - dampl*dampl*dratio1drjk
      ddamps = - damps*damps*dratio3drjk
      drdampdrjk = (ddamps*dampl + damps*ddampl)/r2jk/rjk - 3.0_dp*rdamp/r2jk
!
      deijkdrij = - rdamp*ddampodrij*const                           ! (1/r(dE/dr))
      deijkdrik = - rdamp*ddampodrik*const                           ! (1/r(dE/dr))
      deijkdrjk = - drdampdrjk*dampo*const - rdamp*ddampodrjk*const  ! (1/r(dE/dr))
!-----------------------
!  Second derivatives  |
!-----------------------
      d2expo(1) = - dexpodrij*gfnff_xb_a_cut*(rik + rjk)/rij**3 - 3.0_dp*dexpodrij/r2ij
      d2expo(2) = dexpodrij*gfnff_xb_a_cut/(rij*rik) - dexpodrik/r2ij
      d2expo(3) = dexpodrij*gfnff_xb_a_cut/(rij*rjk) - dexpodrjk/r2ij
      d2expo(4) = dexpodrik*gfnff_xb_a_cut/(rij*rik) - dexpodrik/r2ik
      d2expo(5) = dexpodrik*gfnff_xb_a_cut/(rij*rjk)
      d2expo(6) = dexpodrjk*gfnff_xb_a_cut/(rij*rjk) - dexpodrjk/r2jk
!
      d2dampo(1) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrij - d2expo(1))
      d2dampo(2) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrik - d2expo(2))
      d2dampo(3) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrjk - d2expo(3))
      d2dampo(4) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrik - d2expo(4))
      d2dampo(5) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrjk - d2expo(5))
      d2dampo(6) = 0.5_dp*dampo*dampo*(dampo*dexpodrjk*dexpodrjk - d2expo(6))
!
      d2ratio1drjk2 = 4.0_dp*gfnff_hb_alp*ratio1*(gfnff_hb_alp - 1.0_dp)/r2jk**2
      d2ratio3drjk2 = 4.0_dp*gfnff_hb_alp*ratio3*(gfnff_hb_alp + 1.0_dp)/r2jk**2
      d2dampl = - dampl*dampl*d2ratio1drjk2 + 2.0_dp*dampl*dampl*dampl*dratio1drjk**2
      d2damps = - damps*damps*d2ratio3drjk2 + 2.0_dp*damps*damps*damps*dratio3drjk**2
      d2rdampdrjk2 = (d2damps*dampl + 2.0_dp*ddamps*ddampl + damps*d2dampl)/r2jk/rjk + &
                     15.0_dp*rdamp/r2jk**2 - 6.0_dp*(ddamps*dampl + damps*ddampl)/r2jk/r2jk/rjk
!
!  i-j / i-j
!
      d2trm = - const*rdamp*d2dampo(1)
      call add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                        d2trm,d2r2dx23(1,1))
!
!  i-j / i-k
!
      d2trm = - const*rdamp*d2dampo(2)
      call add_drv2_2pd(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm)
!
!  i-j / j-k
!
      d2trm = - const*(drdampdrjk*ddampodrij + rdamp*d2dampo(3))
      call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm)
!
!  i-k / i-k
!
      d2trm = - const*rdamp*d2dampo(4)
      call add_drv2_1pd(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,deijkdrik, &
                        d2trm,d2r2dx23(1,2))
!
!  i-k / j-k
!
      d2trm = - const*(drdampdrjk*ddampodrik + rdamp*d2dampo(5))
      call add_drv2_2pd(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm)
!
!  j-k / j-k
!
      d2trm = - const*(dampo*d2rdampdrjk2 + 2.0_dp*drdampdrjk*ddampodrjk + rdamp*d2dampo(6))
      call add_drv2_1pd(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,deijkdrjk, &
                        d2trm,d2r2dx23(1,3))
    enddo
  enddo
!
  end subroutine gfnff_hbpd
!
  subroutine gfnff_hb2_eg2pd(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                             maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const_in,radij,xkv,ykv,zkv)
!
!  Calculates the hydrogen bonding energy and derivatives for the GFNFF force fields for nhb2
!  Distributed memory parallel phonon version.
!
!  On entry : 
!
!  lilocal         = is i local to this node?
!  ljlocal         = is j local to this node?
!  i               = first atom of AB to which k is bonded
!  j               = second atom of AB
!  ix,  iy,  iz    = second derivative indices for i if local
!  ixf, iyf, izf   = second derivative indices for i (global)
!  jx,  jy,  jz    = second derivative indices for j if local
!  jxf, jyf, jzf   = second derivative indices for j (global)
!  ni              = pointer to k in bonding lists of i
!  rij             = distance from i to j (A to B)
!  xij             = x component of vector from i to j (A to B)
!  yij             = y component of vector from i to j (A to B)
!  zij             = z component of vector from i to j (A to B)
!  const_in        = energy factor
!
!   4/22 Created from gfnff_hb2_eg2d
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
  use iochannels
  use gulp_gfnff
  use neighbours
  use parallel
  use spatialbo
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: ni
  integer(i4), intent(in)                          :: maxnbr
  integer(i4), intent(in)                          :: nnbr(maxat)
  integer(i4), intent(in)                          :: nbrno(maxnbr,maxat)
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: const_in
  real(dp),    intent(in)                          :: radij
  real(dp),    intent(in)                          :: rij
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: xkv
  real(dp),    intent(in)                          :: ykv
  real(dp),    intent(in)                          :: zkv
  real(dp),    intent(in)                          :: rbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: xbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: ybnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: zbnbr(maxnbr,maxat)
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: indm
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: l
  integer(i4)                                      :: lloc
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: lxf
  integer(i4)                                      :: lyf
  integer(i4)                                      :: lzf
  integer(i4)                                      :: m
  integer(i4)                                      :: mloc
  integer(i4)                                      :: mx
  integer(i4)                                      :: my
  integer(i4)                                      :: mz
  integer(i4)                                      :: mxf
  integer(i4)                                      :: myf
  integer(i4)                                      :: mzf
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  integer(i4)                                      :: nk
  integer(i4)                                      :: status
  logical                                          :: lklocal
  logical                                          :: lllocal
  logical                                          :: lmlocal
  real(dp)                                         :: const
  real(dp)                                         :: damp
  real(dp)                                         :: dampl
  real(dp)                                         :: ddampl
  real(dp)                                         :: d2dampl
  real(dp)                                         :: dampo
  real(dp)                                         :: ddampodrij
  real(dp)                                         :: ddampodrik
  real(dp)                                         :: ddampodrjk
  real(dp)                                         :: d2dampo(6)
  real(dp)                                         :: damps
  real(dp)                                         :: ddamps
  real(dp)                                         :: d2damps
  real(dp)                                         :: d2trm
  real(dp)                                         :: expo
  real(dp)                                         :: dexpodrij
  real(dp)                                         :: dexpodrik
  real(dp)                                         :: dexpodrjk
  real(dp)                                         :: d2expo(6)
  real(dp)                                         :: expo_nb
  real(dp)                                         :: dexpo_nbdrij
  real(dp)                                         :: dexpo_nbdril
  real(dp)                                         :: dexpo_nbdrjl
  real(dp)                                         :: d2expo_nb(6)
  real(dp)                                         :: deijkdrij
  real(dp)                                         :: deijkdrik
  real(dp)                                         :: deijkdril
  real(dp)                                         :: deijkdrjk
  real(dp)                                         :: deijkdrjl
  real(dp)                                         :: hbnbcutloc
  real(dp)                                         :: dampo_nbloc
  real(dp)                                         :: dampo_nb_tot
  real(dp)                                         :: dampo_nb_totm1
  real(dp)                                         :: dampo_nb_totm2
  real(dp),    dimension(:),     allocatable, save :: dampo_nb
  real(dp),    dimension(:,:),   allocatable, save :: ddampo_nbdr
  real(dp),    dimension(:,:),   allocatable, save :: d2dampo_nbdr2
  real(dp)                                         :: p_ab
  real(dp)                                         :: p_bh
  real(dp)                                         :: qhdampo
  real(dp)                                         :: r2ij
  real(dp)                                         :: rik
  real(dp)                                         :: r2ik
  real(dp)                                         :: ril
  real(dp)                                         :: r2il
  real(dp)                                         :: rsum
  real(dp)                                         :: rsum2
  real(dp)                                         :: rjk
  real(dp)                                         :: r2jk
  real(dp)                                         :: rjl
  real(dp)                                         :: r2jl
  real(dp)                                         :: ratio1
  real(dp)                                         :: ratio2
  real(dp)                                         :: ratio2_nb
  real(dp)                                         :: ratio3
  real(dp)                                         :: dratio1drij
  real(dp)                                         :: dratio3drij
  real(dp)                                         :: d2ratio1drij2
  real(dp)                                         :: d2ratio3drij2
  real(dp)                                         :: rdamp
  real(dp)                                         :: drdampdrij
  real(dp)                                         :: drdampdrjk
  real(dp)                                         :: d2rdampdrij2
  real(dp)                                         :: d2rdampdrijdrjk
  real(dp)                                         :: d2rdampdrjk2
  real(dp)                                         :: d2r2dx23(6,3)
  real(dp)                                         :: d2r2dx23l(6,3)
  real(dp)                                         :: rijdamp
  real(dp)                                         :: rjkdamp
  real(dp)                                         :: shortcut
  real(dp)                                         :: xik
  real(dp)                                         :: yik
  real(dp)                                         :: zik
  real(dp)                                         :: xil
  real(dp)                                         :: yil
  real(dp)                                         :: zil
  real(dp)                                         :: xim
  real(dp)                                         :: yim
  real(dp)                                         :: zim
  real(dp)                                         :: xjk
  real(dp)                                         :: yjk
  real(dp)                                         :: zjk
  real(dp)                                         :: xjl
  real(dp)                                         :: yjl
  real(dp)                                         :: zjl
  real(dp)                                         :: xjm
  real(dp)                                         :: yjm
  real(dp)                                         :: zjm
#ifdef TRACE
  call trace_in('gfnff_hb2_eg2pd')
#endif
!
  allocate(dampo_nb(maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg2pd','dampo_nb')
  allocate(ddampo_nbdr(3,maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg2pd','ddampo_nbdr')
  allocate(d2dampo_nbdr2(6,maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg2pd','d2dampo_nbdr2')
!
  k = nbrno(ni,i)   ! H
  kloc = atom2local(k)
  lklocal = (kloc.ne.0)
!
  p_bh = (1.0_dp + gfnff_hbabmix)*gfnff_autoangs**3
  p_ab = - gfnff_hbabmix*gfnff_autoangs**3
!
  r2ij = rij**2
!
!  A-H distance (i-k)
!
  xik = xbnbr(ni,i)
  yik = ybnbr(ni,i)
  zik = zbnbr(ni,i)
  rik = rbnbr(ni,i)
  r2ik = rik**2
!
!  B-H distance (j-k)
!
  xjk = xik - xij
  yjk = yik - yij
  zjk = zik - zij
  r2jk = xjk**2 + yjk**2 + zjk**2
  rjk = sqrt(r2jk)
!
  rsum = rik + rjk + 1.d-12
!
!  Out-of-line damping : A-H...B
!
  expo = (gfnff_hb_a_cut/radij)*(rsum/rij - 1.0_dp)
  if (expo.gt.15.0_dp) goto 1000
  ratio2 = exp(expo)
  dampo = 2.0_dp/(1.0_dp + ratio2)
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
!  Out-of-line damping : A...neighbour(B)-B
!
  if (nat(j).eq.7.and.nnbr(j).eq.1) then
    hbnbcutloc = 2.0_dp*gfnff_autoangs
  else
    hbnbcutloc = gfnff_hb_nb_cut
  endif
!
!  Loop over neighbours of B
!
  dampo_nb_tot = 1.0_dp
  do nj = 1,nnbr(j)
!
!  Compute A - neighbour of B distance
!
    xil = xij + xbnbr(nj,j)
    yil = yij + ybnbr(nj,j)
    zil = zij + zbnbr(nj,j)
    r2il = xil**2 + yil**2 + zil**2
    ril = sqrt(r2il)
    rjl = rbnbr(nj,j)
    r2jl = rjl**2
!
    rsum2 = ril + rjl + 1.0d-12
    expo_nb = (hbnbcutloc/radij)*(rsum2/rij - 1.0_dp)
    ratio2_nb = exp(-expo_nb)
    dampo_nbloc = (2.0_dp/(1.0_dp + ratio2_nb))
    dampo_nb(nj) = dampo_nbloc - 1.0_dp
    dampo_nb_tot = dampo_nb_tot*dampo_nb(nj)
!
    dexpo_nbdrij = (hbnbcutloc/radij)*ratio2_nb*rsum2/rij**3
    dexpo_nbdril = - (hbnbcutloc/radij)*ratio2_nb/(rij*ril)
    dexpo_nbdrjl = - (hbnbcutloc/radij)*ratio2_nb/(rij*rjl)
!
    ddampo_nbdr(1,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrij
    ddampo_nbdr(2,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdril
    ddampo_nbdr(3,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrjl
!
    d2expo_nb(1) = dexpo_nbdrij*(hbnbcutloc/radij)*rsum2/rij**3 - 3.0_dp*dexpo_nbdrij/r2ij
    d2expo_nb(2) = - dexpo_nbdrij*(hbnbcutloc/radij)/(rij*ril) - dexpo_nbdril/r2ij
    d2expo_nb(3) = - dexpo_nbdrij*(hbnbcutloc/radij)/(rij*rjl) - dexpo_nbdrjl/r2ij
    d2expo_nb(4) = - dexpo_nbdril*(hbnbcutloc/radij)/(rij*ril) - dexpo_nbdril/r2il
    d2expo_nb(5) = - dexpo_nbdril*(hbnbcutloc/radij)/(rij*rjl)
    d2expo_nb(6) = - dexpo_nbdrjl*(hbnbcutloc/radij)/(rij*rjl) - dexpo_nbdrjl/r2jl
!
    d2dampo_nbdr2(1,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrij*dexpo_nbdrij - d2expo_nb(1))
    d2dampo_nbdr2(2,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrij*dexpo_nbdril - d2expo_nb(2))
    d2dampo_nbdr2(3,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrij*dexpo_nbdrjl - d2expo_nb(3))
    d2dampo_nbdr2(4,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdril*dexpo_nbdril - d2expo_nb(4))
    d2dampo_nbdr2(5,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdril*dexpo_nbdrjl - d2expo_nb(5))
    d2dampo_nbdr2(6,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrjl*dexpo_nbdrjl - d2expo_nb(6))
  enddo
!
!  Long-range damping
!
  ratio1 = (r2ij/gfnff_hb_long_cut)**gfnff_hb_alp
  dampl = 1.0_dp/(1.0_dp + ratio1)
!
!  Short-range damping
!
  shortcut = gfnff_hb_short_cut*radij
  ratio3 = (shortcut/r2ij)**gfnff_hb_alp
  damps = 1.0_dp/(1.0_dp + ratio3)
!
  damp = damps*dampl
  rijdamp = damp*(p_ab/r2ij/rij)
  rjkdamp = damp*(p_bh/r2jk/rjk)
  rdamp   = rjkdamp + rijdamp
!
  const = const_in*gfnff_hb_ABq(k)
!
  qhdampo = dampo*dampo_nb_tot
!----------------------
!  First derivatives  |
!----------------------
  dexpodrij = - (gfnff_hb_a_cut/radij)*ratio2*rsum/rij**3
  dexpodrik = (gfnff_hb_a_cut/radij)*ratio2/(rij*rik)
  dexpodrjk = (gfnff_hb_a_cut/radij)*ratio2/(rij*rjk)
!
  ddampodrij = - 0.5_dp*dampo*dampo*dexpodrij
  ddampodrik = - 0.5_dp*dampo*dampo*dexpodrik
  ddampodrjk = - 0.5_dp*dampo*dampo*dexpodrjk
!
!  Damping derivatives
!
  dratio1drij = 2.0_dp*gfnff_hb_alp*ratio1/r2ij
  dratio3drij = - 2.0_dp*gfnff_hb_alp*ratio3/r2ij
  ddampl = - dampl*dampl*dratio1drij
  ddamps = - damps*damps*dratio3drij
  drdampdrij = (ddamps*dampl + damps*ddampl)*(p_ab/r2ij/rij + p_bh/r2jk/rjk) - 3.0_dp*rijdamp/r2ij
  drdampdrjk = - 3.0_dp*rjkdamp/r2jk
!
  deijkdrij = - const*rdamp*dampo_nb_tot*ddampodrij - const*qhdampo*drdampdrij
  deijkdrik = - const*rdamp*dampo_nb_tot*ddampodrik
  deijkdrjk = - const*rdamp*dampo_nb_tot*ddampodrjk - const*qhdampo*drdampdrjk
!
  d2r2dx23(1,1) = xij*xij
  d2r2dx23(2,1) = yij*yij
  d2r2dx23(3,1) = zij*zij
  d2r2dx23(4,1) = yij*zij
  d2r2dx23(5,1) = xij*zij
  d2r2dx23(6,1) = xij*yij
!
  d2r2dx23(1,2) = xik*xik
  d2r2dx23(2,2) = yik*yik
  d2r2dx23(3,2) = zik*zik
  d2r2dx23(4,2) = yik*zik
  d2r2dx23(5,2) = xik*zik
  d2r2dx23(6,2) = xik*yik
!
  d2r2dx23(1,3) = xjk*xjk
  d2r2dx23(2,3) = yjk*yjk
  d2r2dx23(3,3) = zjk*zjk
  d2r2dx23(4,3) = yjk*zjk
  d2r2dx23(5,3) = xjk*zjk
  d2r2dx23(6,3) = xjk*yjk
!-----------------------
!  Second derivatives  |
!-----------------------
  d2expo(1) = - dexpodrij*(gfnff_hb_a_cut/radij)*rsum/rij**3 - 3.0_dp*dexpodrij/r2ij
  d2expo(2) = dexpodrij*(gfnff_hb_a_cut/radij)/(rij*rik) - dexpodrik/r2ij
  d2expo(3) = dexpodrij*(gfnff_hb_a_cut/radij)/(rij*rjk) - dexpodrjk/r2ij
  d2expo(4) = dexpodrik*(gfnff_hb_a_cut/radij)/(rij*rik) - dexpodrik/r2ik
  d2expo(5) = dexpodrik*(gfnff_hb_a_cut/radij)/(rij*rjk)
  d2expo(6) = dexpodrjk*(gfnff_hb_a_cut/radij)/(rij*rjk) - dexpodrjk/r2jk
!
  d2dampo(1) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrij - d2expo(1))
  d2dampo(2) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrik - d2expo(2))
  d2dampo(3) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrjk - d2expo(3))
  d2dampo(4) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrik - d2expo(4))
  d2dampo(5) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrjk - d2expo(5))
  d2dampo(6) = 0.5_dp*dampo*dampo*(dampo*dexpodrjk*dexpodrjk - d2expo(6))
!
  d2ratio1drij2 = 4.0_dp*gfnff_hb_alp*ratio1*(gfnff_hb_alp - 1.0_dp)/r2ij**2
  d2ratio3drij2 = 4.0_dp*gfnff_hb_alp*ratio3*(gfnff_hb_alp + 1.0_dp)/r2ij**2
  d2dampl = - dampl*dampl*d2ratio1drij2 + 2.0_dp*dampl*dampl*dampl*dratio1drij**2
  d2damps = - damps*damps*d2ratio3drij2 + 2.0_dp*damps*damps*damps*dratio3drij**2
  d2rdampdrij2 = (d2damps*dampl + 2.0_dp*ddamps*ddampl + damps*d2dampl)*(p_ab/r2ij/rij + p_bh/r2jk/rjk) + &
                 15.0_dp*rijdamp/r2ij**2 - 6.0_dp*(ddamps*dampl + damps*ddampl)*p_ab/r2ij/r2ij/rij
  d2rdampdrijdrjk = - 3.0_dp*(ddamps*dampl + damps*ddampl)*p_bh/r2jk/r2jk/rjk
  d2rdampdrjk2 = 15.0_dp*rjkdamp/r2jk**2
!
!  i-j / i-j
!
  d2trm = - const*dampo_nb_tot*(dampo*d2rdampdrij2 + 2.0_dp*drdampdrij*ddampodrij + rdamp*d2dampo(1))
  call add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                    d2trm,d2r2dx23(1,1))
!
!  i-j / i-k
!
  d2trm = - const*dampo_nb_tot*(drdampdrij*ddampodrik + rdamp*d2dampo(2))
  call add_drv2_2pd(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm)
!
!  i-j / j-k
!
  d2trm = - const*dampo_nb_tot*(dampo*d2rdampdrijdrjk + drdampdrij*ddampodrjk + drdampdrjk*ddampodrij + &
            rdamp*d2dampo(3))
  call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm)
!
!  i-k / i-k
!
  d2trm = - const*dampo_nb_tot*rdamp*d2dampo(4)
  call add_drv2_1pd(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,deijkdrik, &
                    d2trm,d2r2dx23(1,2))
!
!  i-k / j-k
!
  d2trm = - const*dampo_nb_tot*(drdampdrjk*ddampodrik + rdamp*d2dampo(5))
  call add_drv2_2pd(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm)
!
!  j-k / j-k
!
  d2trm = - const*dampo_nb_tot*(dampo*d2rdampdrjk2 + 2.0_dp*drdampdrjk*ddampodrjk + rdamp*d2dampo(6))
  call add_drv2_1pd(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,deijkdrjk, &
                    d2trm,d2r2dx23(1,3))
!
!  Compute A - neighbour of B distance contributions to derivatives
!  
  do nj = 1,nnbr(j)
    l = nbrno(nj,j)
    lloc = atom2local(l)
    lllocal = (lloc.ne.0) 
!
    xjl = xbnbr(nj,j)
    yjl = ybnbr(nj,j)
    zjl = zbnbr(nj,j)
!
    xil = xij + xjl
    yil = yij + yjl
    zil = zij + zjl
!
!  Compute product of factors except the current one
!
    dampo_nb_totm1 = 1.0_dp
    do nj2 = 1,nnbr(j)
      if (nj2.ne.nj) then
        dampo_nb_totm1 = dampo_nb_totm1*dampo_nb(nj2)
      endif
    enddo
!
    deijkdrij = - const*rdamp*dampo_nb_totm1*dampo*ddampo_nbdr(1,nj)   ! (1/r(dE/dr))
    deijkdril = - const*rdamp*dampo_nb_totm1*dampo*ddampo_nbdr(2,nj)   ! (1/r(dE/dr))
    deijkdrjl = - const*rdamp*dampo_nb_totm1*dampo*ddampo_nbdr(3,nj)   ! (1/r(dE/dr))
!
    d2r2dx23l(1,1) = xij*xij
    d2r2dx23l(2,1) = yij*yij
    d2r2dx23l(3,1) = zij*zij
    d2r2dx23l(4,1) = yij*zij
    d2r2dx23l(5,1) = xij*zij
    d2r2dx23l(6,1) = xij*yij
!
    d2r2dx23l(1,2) = xil*xil
    d2r2dx23l(2,2) = yil*yil
    d2r2dx23l(3,2) = zil*zil
    d2r2dx23l(4,2) = yil*zil
    d2r2dx23l(5,2) = xil*zil
    d2r2dx23l(6,2) = xil*yil
!
    d2r2dx23l(1,3) = xjl*xjl
    d2r2dx23l(2,3) = yjl*yjl
    d2r2dx23l(3,3) = zjl*zjl
    d2r2dx23l(4,3) = yjl*zjl
    d2r2dx23l(5,3) = xjl*zjl
    d2r2dx23l(6,3) = xjl*yjl
!
    indl = 3*(l-1)
    lxf = indl + 1
    lyf = indl + 2
    lzf = indl + 3
    if (lllocal) then
      indl = 3*(lloc-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
    endif
!---------------------------------------------------------------------------
!  Contribution from second derivatives of dampo_nb for a single distance  |
!---------------------------------------------------------------------------
!
!  i-j / i-j
!
    d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(1,nj)
    call add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                      d2trm,d2r2dx23l(1,1))
!
!  i-j / i-l
!
    d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(2,nj)
    call add_drv2_2pd(lilocal,ljlocal,lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xil,yil,zil,d2trm)
!
!  i-j / j-l
!
    d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(3,nj)
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm)
!
!  i-l / i-l
!
    d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(4,nj)
    call add_drv2_1pd(lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xil,yil,zil,deijkdril, &
                      d2trm,d2r2dx23l(1,2))
!
!  i-l / j-l
!
    d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(5,nj)
    call add_drv2_2pd(lilocal,lllocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xil,yil,zil,xjl,yjl,zjl,d2trm)
!
!  j-l / j-l
!
    d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(6,nj)
    call add_drv2_1pd(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,deijkdrjl, &
                      d2trm,d2r2dx23l(1,3))
!------------------------------------------------
!  Contribution from dampo_nb and 1 other term  |
!------------------------------------------------
!
!  i-j / i-j
!
    d2trm = - const*dampo_nb_totm1*ddampo_nbdr(1,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
    call add_drv2_2pd(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm)
!
!  i-j / i-l
!
    d2trm = - const*dampo_nb_totm1*ddampo_nbdr(2,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
    call add_drv2_2pd(lilocal,ljlocal,lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xil,yil,zil,d2trm)
!
!  i-j / j-l
!
    d2trm = - const*dampo_nb_totm1*ddampo_nbdr(3,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm)
!
!  i-k / i-j
!
    d2trm = - const*dampo_nb_totm1*ddampo_nbdr(1,nj)*rdamp*ddampodrik
    call add_drv2_2pd(lilocal,lklocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xik,yik,zik,xij,yij,zij,d2trm)
!
!  i-k / i-l
!
    d2trm = - const*dampo_nb_totm1*ddampo_nbdr(2,nj)*rdamp*ddampodrik
    call add_drv2_2pd(lilocal,lklocal,lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xil,yil,zil,d2trm)
!
!  i-k / j-l
!
    d2trm = - const*dampo_nb_totm1*ddampo_nbdr(3,nj)*rdamp*ddampodrik
    call add_drv2_2pd(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xjl,yjl,zjl,d2trm)
!
!  j-k / i-j
!
    d2trm = - const*dampo_nb_totm1*ddampo_nbdr(1,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
    call add_drv2_2pd(ljlocal,lklocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjk,yjk,zjk,xij,yij,zij,d2trm)
!
!  j-k / i-l
!
    d2trm = - const*dampo_nb_totm1*ddampo_nbdr(2,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
    call add_drv2_2pd(ljlocal,lklocal,lilocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xil,yil,zil,d2trm)
!
!  j-k / j-l
!
    d2trm = - const*dampo_nb_totm1*ddampo_nbdr(3,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm)
!
!  Loop over second neighbour of B for products of first derivatives
!
    do nk = 1,nj-1
      m = nbrno(nk,j)
      mloc = atom2local(m)
      lmlocal = (mloc.ne.0) 
!
      if (.not.lilocal.and..not.ljlocal.and..not.lllocal.and..not.lmlocal) cycle
!
      indm = 3*(m-1)
      mxf = indm + 1
      myf = indm + 2
      mzf = indm + 3
      if (lmlocal) then
        indm = 3*(mloc-1)
        mx = indm + 1
        my = indm + 2
        mz = indm + 3
      endif
!
      xjm = xbnbr(nk,j)
      yjm = ybnbr(nk,j)
      zjm = zbnbr(nk,j)
!
      xim = xij + xjm
      yim = yij + yjm
      zim = zij + zjm
!
!  Compute product of factors except the current one
!
      dampo_nb_totm2 = 1.0_dp
      do nj2 = 1,nnbr(j)
        if (nj2.ne.nj.and.nj2.ne.nk) then
          dampo_nb_totm2 = dampo_nb_totm2*dampo_nb(nj2)
        endif
      enddo
!
!  i-j / i-j
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(1,nk)
      call add_drv2_2pd(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm)
!
!  i-j / i-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(2,nk)
      call add_drv2_2pd(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm)
!
!  i-j / j-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(3,nk)
      call add_drv2_2pd(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm)
!
!  i-l / i-j
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(1,nk)
      call add_drv2_2pd(lilocal,lllocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xil,yil,zil,xij,yij,zij,d2trm)
!
!  i-l / i-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(2,nk)
      call add_drv2_2pd(lilocal,lllocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                        ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xil,yil,zil,xim,yim,zim,d2trm)
!
!  i-l / j-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(3,nk)
      call add_drv2_2pd(lilocal,lllocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xil,yil,zil,xjm,yjm,zjm,d2trm)
!
!  j-l / i-j
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(1,nk)
      call add_drv2_2pd(ljlocal,lllocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjl,yjl,zjl,xij,yij,zij,d2trm)
!
!  j-l / i-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(2,nk)
      call add_drv2_2pd(ljlocal,lllocal,lilocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xim,yim,zim,d2trm)
!
!  j-l / j-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(3,nk)
      call add_drv2_2pd(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xjm,yjm,zjm,d2trm)
    enddo
  enddo
!
!  Exit point to ensure memory is deallocated
!
1000 continue
!
  deallocate(d2dampo_nbdr2,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg2pd','d2dampo_nbdr2')
  deallocate(ddampo_nbdr,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg2pd','ddampo_nbdr')
  deallocate(dampo_nb,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg2pd','dampo_nb')
!
#ifdef TRACE
  call trace_out('gfnff_hb2_eg2pd')
#endif
!
  return
  end
!
  subroutine gfnff_hb2_eg2pd_rnr(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                 maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const_in,radij,xkv,ykv,zkv)
!
!  Calculates the hydrogen bonding energy and derivatives for the GFNFF force fields for nhb2
!  Lone pair version. Distributed memory parallel phonon version.
!
!  On entry : 
!
!  lilocal         = is i local to this node?
!  ljlocal         = is j local to this node?
!  i               = first atom of AB to which k is bonded
!  j               = second atom of AB
!  ix,  iy,  iz    = second derivative indices for i if local
!  ixf, iyf, izf   = second derivative indices for i (global)
!  jx,  jy,  jz    = second derivative indices for j if local
!  jxf, jyf, jzf   = second derivative indices for j (global)
!  ni              = pointer to k in bonding lists of i
!  rij             = distance from i to j (A to B)
!  xij             = x component of vector from i to j (A to B)
!  yij             = y component of vector from i to j (A to B)
!  zij             = z component of vector from i to j (A to B)
!  const_in        = energy factor
!
!   4/22 Created from gfnff_hb2_eg2d_rnr
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
  use gulp_gfnff
  use iochannels
  use neighbours
  use parallel
  use spatialbo
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: ni
  integer(i4), intent(in)                          :: maxnbr
  integer(i4), intent(in)                          :: nnbr(maxat)
  integer(i4), intent(in)                          :: nbrno(maxnbr,maxat)
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: const_in
  real(dp),    intent(in)                          :: radij
  real(dp),    intent(in)                          :: rij
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: xkv
  real(dp),    intent(in)                          :: ykv
  real(dp),    intent(in)                          :: zkv
  real(dp),    intent(in)                          :: rbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: xbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: ybnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: zbnbr(maxnbr,maxat)
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: indm
  integer(i4)                                      :: indn
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: l
  integer(i4)                                      :: lloc
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: lxf
  integer(i4)                                      :: lyf
  integer(i4)                                      :: lzf
  integer(i4)                                      :: m
  integer(i4)                                      :: mloc
  integer(i4)                                      :: mx
  integer(i4)                                      :: my
  integer(i4)                                      :: mz
  integer(i4)                                      :: mxf
  integer(i4)                                      :: myf
  integer(i4)                                      :: mzf
  integer(i4)                                      :: n
  integer(i4)                                      :: nloc
  integer(i4)                                      :: nx
  integer(i4)                                      :: ny
  integer(i4)                                      :: nz
  integer(i4)                                      :: nxf
  integer(i4)                                      :: nyf
  integer(i4)                                      :: nzf
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  integer(i4)                                      :: nk
  integer(i4)                                      :: status
  logical                                          :: lijmlocal
  logical                                          :: lklocal
  logical                                          :: lllocal
  logical                                          :: lmlocal
  logical                                          :: lnlocal
  logical                                          :: llpok
  real(dp)                                         :: const
  real(dp)                                         :: d2trm
  real(dp)                                         :: damp
  real(dp)                                         :: dampl
  real(dp)                                         :: ddampl
  real(dp)                                         :: d2dampl
  real(dp)                                         :: damps
  real(dp)                                         :: ddamps
  real(dp)                                         :: d2damps
  real(dp)                                         :: dampo
  real(dp)                                         :: ddampodrij
  real(dp)                                         :: ddampodrik
  real(dp)                                         :: ddampodrjk
  real(dp)                                         :: d2dampo(6)
  real(dp)                                         :: dampo_lp
  real(dp)                                         :: ddampo_lpdralp
  real(dp)                                         :: ddampo_lpdrij
  real(dp)                                         :: d2dampo_lp(3)
  real(dp)                                         :: dampo_nbloc
  real(dp)                                         :: dampo_nb_tot
  real(dp)                                         :: dampo_nb_totm1
  real(dp)                                         :: dampo_nb_totm2
  real(dp),    dimension(:),     allocatable, save :: dampo_nb
  real(dp),    dimension(:,:),   allocatable, save :: ddampo_nbdr
  real(dp),    dimension(:,:),   allocatable, save :: d2dampo_nbdr2
  real(dp)                                         :: dr2alp(3)
  real(dp)                                         :: dr2alpij(3)
  real(dp)                                         :: dxblp(3)
  real(dp)                                         :: d2xblp(3,3)
  real(dp)                                         :: dyblp(3)
  real(dp)                                         :: d2yblp(3,3)
  real(dp)                                         :: dzblp(3)
  real(dp)                                         :: d2zblp(3,3)
  real(dp)                                         :: d2r2dx23(6,3)
  real(dp)                                         :: d2r2dx23m(6,3)
  real(dp)                                         :: d2r2dx2ra(6)
  real(dp)                                         :: expo
  real(dp)                                         :: dexpodrij
  real(dp)                                         :: dexpodrik
  real(dp)                                         :: dexpodrjk
  real(dp)                                         :: d2expo(6)
  real(dp)                                         :: expo_lp
  real(dp)                                         :: dexpo_lpdrij
  real(dp)                                         :: dexpo_lpdralp
  real(dp)                                         :: d2expo_lp(6)
  real(dp)                                         :: expo_nb
  real(dp)                                         :: dexpo_nbdrij
  real(dp)                                         :: dexpo_nbdrim
  real(dp)                                         :: dexpo_nbdrjm
  real(dp)                                         :: d2expo_nb(6)
  real(dp)                                         :: deijkdrij
  real(dp)                                         :: deijkdrik
  real(dp)                                         :: deijkdrim
  real(dp)                                         :: deijkdrjk
  real(dp)                                         :: deijkdrjm
  real(dp)                                         :: dtrmalp
  real(dp)                                         :: d2trmalp2
  real(dp)                                         :: hbnbcutloc
  real(dp)                                         :: p_ab
  real(dp)                                         :: p_bh
  real(dp)                                         :: qhdampo
  real(dp)                                         :: r2ij
  real(dp)                                         :: rik
  real(dp)                                         :: r2ik
  real(dp)                                         :: rim
  real(dp)                                         :: r2im
  real(dp)                                         :: ralp
  real(dp)                                         :: ralp2
  real(dp)                                         :: rblp
  real(dp)                                         :: rblp2
  real(dp)                                         :: rblpsum
  real(dp)                                         :: rblpsum2
  real(dp)                                         :: rsum
  real(dp)                                         :: rsum2
  real(dp)                                         :: rlp_dist
  real(dp)                                         :: rjk
  real(dp)                                         :: r2jk
  real(dp)                                         :: rjm
  real(dp)                                         :: r2jm
  real(dp)                                         :: ralpprblp
  real(dp)                                         :: ratio1
  real(dp)                                         :: dratio1drij
  real(dp)                                         :: d2ratio1drij2
  real(dp)                                         :: ratio2
  real(dp)                                         :: ratio2_lp
  real(dp)                                         :: ratio2_nb
  real(dp)                                         :: ratio3
  real(dp)                                         :: dratio3drij
  real(dp)                                         :: d2ratio3drij2
  real(dp)                                         :: rdamp
  real(dp)                                         :: drdampdrij
  real(dp)                                         :: d2rdampdrij2
  real(dp)                                         :: d2rdampdrijdrjk
  real(dp)                                         :: drdampdrjk
  real(dp)                                         :: d2rdampdrjk2
  real(dp)                                         :: rhblpcut
  real(dp)                                         :: rijdamp
  real(dp)                                         :: rjkdamp
  real(dp)                                         :: shortcut
  real(dp)                                         :: xalp
  real(dp)                                         :: yalp
  real(dp)                                         :: zalp
  real(dp)                                         :: xblp
  real(dp)                                         :: yblp
  real(dp)                                         :: zblp
  real(dp)                                         :: xblpsum
  real(dp)                                         :: yblpsum
  real(dp)                                         :: zblpsum
  real(dp)                                         :: xik
  real(dp)                                         :: yik
  real(dp)                                         :: zik
  real(dp)                                         :: xim
  real(dp)                                         :: yim
  real(dp)                                         :: zim
  real(dp)                                         :: xin
  real(dp)                                         :: yin
  real(dp)                                         :: zin
  real(dp)                                         :: xjk
  real(dp)                                         :: yjk
  real(dp)                                         :: zjk
  real(dp)                                         :: xjl
  real(dp)                                         :: yjl
  real(dp)                                         :: zjl
  real(dp)                                         :: xjm
  real(dp)                                         :: yjm
  real(dp)                                         :: zjm
  real(dp)                                         :: xjn
  real(dp)                                         :: yjn
  real(dp)                                         :: zjn
  real(dp)                                         :: d2r2dx2mix(3,3)
  real(dp)                                         :: d2r2blpdx2(6)
#ifdef TRACE
  call trace_in('gfnff_hb2_eg2pd_rnr')
#endif
!
  allocate(dampo_nb(maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg2pd_rnr','dampo_nb')
  allocate(ddampo_nbdr(3,maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg2pd_rnr','ddampo_nbdr')
  allocate(d2dampo_nbdr2(6,maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg2pd_rnr','d2dampo_nbdr2')
!
  k = nbrno(ni,i)   ! H
  kloc = atom2local(k)
  lklocal = (kloc.ne.0)
!
  p_bh = (1.0_dp + gfnff_hbabmix)*gfnff_autoangs**3
  p_ab = - gfnff_hbabmix*gfnff_autoangs**3
!
  r2ij = rij**2
!
  rlp_dist = (0.50_dp - 0.018_dp*gfnff_repulsion_z(nat(j)))*gfnff_autoangs
  rhblpcut = 56.0_dp*gfnff_autoangs
!
!  Find lone pair position using negative sum of vectors for neighbours to B
!
  xblpsum = 0.0_dp
  yblpsum = 0.0_dp
  zblpsum = 0.0_dp
  do nj = 1,nnbr(j)
    xblpsum = xblpsum + xbnbr(nj,j)
    yblpsum = yblpsum + ybnbr(nj,j)
    zblpsum = zblpsum + zbnbr(nj,j)
  enddo
  rblpsum2 = xblpsum**2 + yblpsum**2 + zblpsum**2
  rblpsum = sqrt(rblpsum2)
  if (rblpsum.gt.1.0d-10) then
    llpok = .true.
    xblp = - rlp_dist*xblpsum/rblpsum
    yblp = - rlp_dist*yblpsum/rblpsum
    zblp = - rlp_dist*zblpsum/rblpsum
  else
    llpok = .false.
    xblp = 0.0_dp
    yblp = 0.0_dp
    zblp = 0.0_dp
  endif
!
!  A-H distance (i-k)
!
  xik = xbnbr(ni,i)
  yik = ybnbr(ni,i)
  zik = zbnbr(ni,i)
  rik = rbnbr(ni,i)
  r2ik = rik**2
!
!  B-H distance (j-k)
!
  xjk = xik - xij
  yjk = yik - yij
  zjk = zik - zij
  r2jk = xjk**2 + yjk**2 + zjk**2
  rjk = sqrt(r2jk)
!
  rsum = rik + rjk + 1.d-12
!
!  Out-of-line damping : A-H...B
!
  expo = (gfnff_hb_a_cut/radij)*(rsum/rij - 1.0_dp)
  if (expo.gt.15.0_dp) goto 1000
  ratio2 = exp(expo)
  dampo = 2.0_dp/(1.0_dp + ratio2)
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
!  Compute distances from A and B to LP
!
  xalp = xblp + xij
  yalp = yblp + yij
  zalp = zblp + zij
  ralp2 = xalp**2 + yalp**2 + zalp**2
  ralp = sqrt(ralp2)
!
  rblp2 = xblp**2 + yblp**2 + zblp**2
  rblp  = sqrt(rblp2)
!
!  Out-of-line damping : A...LP-B
!
  ralpprblp = ralp + rblp + 1.0d-12
  expo_lp = (rhblpcut/radij)*(ralpprblp/rij - 1.0_dp)
  ratio2_lp = exp(expo_lp)
  dampo_lp = 2.0_dp/(1.0_dp + ratio2_lp)
!
!  Loop over neighbours of B
!
  if (nat(j).eq.7.and.nnbr(j).eq.1) then
    hbnbcutloc = 2.0_dp*gfnff_autoangs
  else
    hbnbcutloc = gfnff_hb_nb_cut
  endif
  dampo_nb_tot = 1.0_dp
  if (llpok) then
    do nj = 1,nnbr(j)
!
!  Compute A - neighbour of B distance
!
      xim = xij + xbnbr(nj,j)
      yim = yij + ybnbr(nj,j)
      zim = zij + zbnbr(nj,j)
      r2im = xim**2 + yim**2 + zim**2
      rim = sqrt(r2im)
      rjm = rbnbr(nj,j)
      r2jm = rjm**2
!
      rsum2 = rim + rjm + 1.0d-12
      expo_nb = (hbnbcutloc/radij)*(rsum2/rij - 1.0_dp)
      ratio2_nb = exp(-expo_nb)
      dampo_nbloc = (2.0_dp/(1.0_dp + ratio2_nb))
      dampo_nb(nj) = dampo_nbloc - 1.0_dp
      dampo_nb_tot = dampo_nb_tot*dampo_nb(nj)
!
      dexpo_nbdrij = (hbnbcutloc/radij)*ratio2_nb*rsum2/rij**3
      dexpo_nbdrim = - (hbnbcutloc/radij)*ratio2_nb/(rij*rim)
      dexpo_nbdrjm = - (hbnbcutloc/radij)*ratio2_nb/(rij*rjm)
!
      ddampo_nbdr(1,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrij
      ddampo_nbdr(2,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrim
      ddampo_nbdr(3,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrjm
!
      d2expo_nb(1) = dexpo_nbdrij*(hbnbcutloc/radij)*rsum2/rij**3 - 3.0_dp*dexpo_nbdrij/r2ij
      d2expo_nb(2) = - dexpo_nbdrij*(hbnbcutloc/radij)/(rij*rim) - dexpo_nbdrim/r2ij
      d2expo_nb(3) = - dexpo_nbdrij*(hbnbcutloc/radij)/(rij*rjm) - dexpo_nbdrjm/r2ij
      d2expo_nb(4) = - dexpo_nbdrim*(hbnbcutloc/radij)/(rij*rim) - dexpo_nbdrim/r2im
      d2expo_nb(5) = - dexpo_nbdrim*(hbnbcutloc/radij)/(rij*rjm)
      d2expo_nb(6) = - dexpo_nbdrjm*(hbnbcutloc/radij)/(rij*rjm) - dexpo_nbdrjm/r2jm
!
      d2dampo_nbdr2(1,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrij*dexpo_nbdrij - d2expo_nb(1))
      d2dampo_nbdr2(2,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrij*dexpo_nbdrim - d2expo_nb(2))
      d2dampo_nbdr2(3,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrij*dexpo_nbdrjm - d2expo_nb(3))
      d2dampo_nbdr2(4,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrim*dexpo_nbdrim - d2expo_nb(4))
      d2dampo_nbdr2(5,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrim*dexpo_nbdrjm - d2expo_nb(5))
      d2dampo_nbdr2(6,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrjm*dexpo_nbdrjm - d2expo_nb(6))
    enddo
  endif
!
!  Long-range damping
!
  ratio1 = (r2ij/gfnff_hb_long_cut)**gfnff_hb_alp
  dampl = 1.0_dp/(1.0_dp + ratio1)
!
!  Short-range damping
!
  shortcut = gfnff_hb_short_cut*radij
  ratio3 = (shortcut/r2ij)**gfnff_hb_alp
  damps = 1.0_dp/(1.0_dp + ratio3)
!
  damp = damps*dampl
  rijdamp = damp*(p_ab/r2ij/rij)
  rjkdamp = damp*(p_bh/r2jk/rjk)
  rdamp   = rjkdamp + rijdamp
!
  qhdampo = dampo*dampo_nb_tot*dampo_lp
!
  const = const_in*gfnff_hb_ABq(k)
!----------------------
!  First derivatives  |
!----------------------
  dexpodrij = - (gfnff_hb_a_cut/radij)*ratio2*rsum/rij**3
  dexpodrik = (gfnff_hb_a_cut/radij)*ratio2/(rij*rik)
  dexpodrjk = (gfnff_hb_a_cut/radij)*ratio2/(rij*rjk)
!
  ddampodrij = - 0.5_dp*dampo*dampo*dexpodrij
  ddampodrik = - 0.5_dp*dampo*dampo*dexpodrik
  ddampodrjk = - 0.5_dp*dampo*dampo*dexpodrjk
!
  dexpo_lpdrij = - (rhblpcut/radij)*ratio2_lp*ralpprblp/rij**3
  dexpo_lpdralp = (rhblpcut/radij)*ratio2_lp/(rij*ralp)
  ddampo_lpdrij = - 0.5_dp*dampo_lp*dampo_lp*dexpo_lpdrij
  ddampo_lpdralp = - 0.5_dp*dampo_lp*dampo_lp*dexpo_lpdralp
!---------------------------------------
!  Derivatives of damping and qhdampo  |
!---------------------------------------
  dratio1drij = 2.0_dp*gfnff_hb_alp*ratio1/r2ij
  dratio3drij = - 2.0_dp*gfnff_hb_alp*ratio3/r2ij
  ddampl = - dampl*dampl*dratio1drij
  ddamps = - damps*damps*dratio3drij
  drdampdrij = (ddamps*dampl + damps*ddampl)*(p_ab/r2ij/rij + p_bh/r2jk/rjk) - 3.0_dp*rijdamp/r2ij
  drdampdrjk = - 3.0_dp*rjkdamp/r2jk
!
  deijkdrij = - const*rdamp*dampo_nb_tot*(dampo_lp*ddampodrij + dampo*ddampo_lpdrij) - const*qhdampo*drdampdrij
  deijkdrik = - const*rdamp*dampo_nb_tot*dampo_lp*ddampodrik
  deijkdrjk = - const*rdamp*dampo_nb_tot*dampo_lp*ddampodrjk - const*qhdampo*drdampdrjk
!
  d2r2dx23(1,1) = xij*xij
  d2r2dx23(2,1) = yij*yij
  d2r2dx23(3,1) = zij*zij
  d2r2dx23(4,1) = yij*zij
  d2r2dx23(5,1) = xij*zij
  d2r2dx23(6,1) = xij*yij
!
  d2r2dx23(1,2) = xik*xik
  d2r2dx23(2,2) = yik*yik
  d2r2dx23(3,2) = zik*zik
  d2r2dx23(4,2) = yik*zik
  d2r2dx23(5,2) = xik*zik
  d2r2dx23(6,2) = xik*yik
!
  d2r2dx23(1,3) = xjk*xjk
  d2r2dx23(2,3) = yjk*yjk
  d2r2dx23(3,3) = zjk*zjk
  d2r2dx23(4,3) = yjk*zjk
  d2r2dx23(5,3) = xjk*zjk
  d2r2dx23(6,3) = xjk*yjk
!
!  Derivatives of ralp : rij component
!
  dtrmalp = - const*rdamp*dampo_nb_tot*dampo*ddampo_lpdralp
!-----------------------
!  Second derivatives  |
!-----------------------
!
!  Set up terms for second derivatives
!
  d2expo_lp(1) = - dexpo_lpdrij*(rhblpcut/radij)*ralpprblp/rij**3 - 3.0_dp*dexpo_lpdrij/r2ij
  d2expo_lp(2) = dexpo_lpdrij*(rhblpcut/radij)/(rij*ralp) - dexpo_lpdralp/r2ij
  d2expo_lp(3) = dexpo_lpdralp*(rhblpcut/radij)/(rij*ralp) - dexpo_lpdralp/ralp2
!
  d2dampo_lp(1) = 0.5_dp*dampo_lp*dampo_lp*(dampo_lp*dexpo_lpdrij*dexpo_lpdrij - d2expo_lp(1))
  d2dampo_lp(2) = 0.5_dp*dampo_lp*dampo_lp*(dampo_lp*dexpo_lpdrij*dexpo_lpdralp - d2expo_lp(2))
  d2dampo_lp(3) = 0.5_dp*dampo_lp*dampo_lp*(dampo_lp*dexpo_lpdralp*dexpo_lpdralp - d2expo_lp(3))
!
  d2trmalp2 = - const*rdamp*dampo_nb_tot*dampo*d2dampo_lp(3)
!
!  i-j / i-j for ralp
!
  d2r2dx2ra(1) = xalp*xalp
  d2r2dx2ra(2) = yalp*yalp
  d2r2dx2ra(3) = zalp*zalp
  d2r2dx2ra(4) = yalp*zalp
  d2r2dx2ra(5) = xalp*zalp
  d2r2dx2ra(6) = xalp*yalp
!
  call add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xalp,yalp,zalp, &
                    dtrmalp,d2trmalp2,d2r2dx2ra)
!
!  Mixed terms
!
!  i-j / i-j
!
  d2trm = - const*dampo_nb_tot*(dampo*rdamp*d2dampo_lp(2) + dampo*drdampdrij*ddampo_lpdralp + &
            rdamp*ddampodrij*ddampo_lpdralp)
  call add_drv2_2pd(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xalp,yalp,zalp,d2trm)
!
!  Set up terms for second derivatives of non-lp term
!
  d2expo(1) = - dexpodrij*(gfnff_hb_a_cut/radij)*rsum/rij**3 - 3.0_dp*dexpodrij/r2ij
  d2expo(2) = dexpodrij*(gfnff_hb_a_cut/radij)/(rij*rik) - dexpodrik/r2ij
  d2expo(3) = dexpodrij*(gfnff_hb_a_cut/radij)/(rij*rjk) - dexpodrjk/r2ij
  d2expo(4) = dexpodrik*(gfnff_hb_a_cut/radij)/(rij*rik) - dexpodrik/r2ik
  d2expo(5) = dexpodrik*(gfnff_hb_a_cut/radij)/(rij*rjk)
  d2expo(6) = dexpodrjk*(gfnff_hb_a_cut/radij)/(rij*rjk) - dexpodrjk/r2jk
!
  d2dampo(1) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrij - d2expo(1))
  d2dampo(2) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrik - d2expo(2))
  d2dampo(3) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrjk - d2expo(3))
  d2dampo(4) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrik - d2expo(4))
  d2dampo(5) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrjk - d2expo(5))
  d2dampo(6) = 0.5_dp*dampo*dampo*(dampo*dexpodrjk*dexpodrjk - d2expo(6))
!
  d2ratio1drij2 = 4.0_dp*gfnff_hb_alp*ratio1*(gfnff_hb_alp - 1.0_dp)/r2ij**2
  d2ratio3drij2 = 4.0_dp*gfnff_hb_alp*ratio3*(gfnff_hb_alp + 1.0_dp)/r2ij**2
  d2dampl = - dampl*dampl*d2ratio1drij2 + 2.0_dp*dampl*dampl*dampl*dratio1drij**2
  d2damps = - damps*damps*d2ratio3drij2 + 2.0_dp*damps*damps*damps*dratio3drij**2
  d2rdampdrij2 = (d2damps*dampl + 2.0_dp*ddamps*ddampl + damps*d2dampl)*(p_ab/r2ij/rij + p_bh/r2jk/rjk) + &
                 15.0_dp*rijdamp/r2ij**2 - 6.0_dp*(ddamps*dampl + damps*ddampl)*p_ab/r2ij/r2ij/rij
  d2rdampdrijdrjk = - 3.0_dp*(ddamps*dampl + damps*ddampl)*p_bh/r2jk/r2jk/rjk
  d2rdampdrjk2 = 15.0_dp*rjkdamp/r2jk**2
!
!  i-j / i-j
!
  d2trm = - const*dampo_nb_tot*(dampo_lp*(dampo*d2rdampdrij2 + 2.0_dp*drdampdrij*ddampodrij + rdamp*d2dampo(1)) + &
            dampo*rdamp*d2dampo_lp(1) + 2.0_dp*dampo*drdampdrij*ddampo_lpdrij + &
            2.0_dp*rdamp*ddampodrij*ddampo_lpdrij)
  call add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                    deijkdrij,d2trm,d2r2dx23(1,1))
!
!  i-j / i-k
!
  d2trm = - const*dampo_nb_tot*(dampo_lp*(drdampdrij*ddampodrik + rdamp*d2dampo(2)) + &
            rdamp*ddampodrik*ddampo_lpdrij)
  call add_drv2_2pd(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm)
!
!  i-j (ralp)/ i-k 
!
  d2trm = - const*dampo_nb_tot*rdamp*ddampodrik*ddampo_lpdralp
  call add_drv2_2pd(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xalp,yalp,zalp,xik,yik,zik,d2trm)

!  i-j / j-k
!
  d2trm = - const*dampo_nb_tot*(dampo_lp*(dampo*d2rdampdrijdrjk + drdampdrij*ddampodrjk + drdampdrjk*ddampodrij + &
            rdamp*d2dampo(3)) + (rdamp*ddampodrjk + dampo*drdampdrjk)*ddampo_lpdrij)
  call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm)
!
!  i-j (ralp)/ j-k 
!
  d2trm = - const*dampo_nb_tot*(rdamp*ddampodrjk + dampo*drdampdrjk)*ddampo_lpdralp
  call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xalp,yalp,zalp,xjk,yjk,zjk,d2trm)

!  i-k / i-k
!
  d2trm = - const*dampo_nb_tot*dampo_lp*rdamp*d2dampo(4)
  call add_drv2_1pd(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik, &
                    deijkdrik,d2trm,d2r2dx23(1,2))
!
!  i-k / j-k
!
  d2trm = - const*dampo_nb_tot*dampo_lp*(drdampdrjk*ddampodrik + rdamp*d2dampo(5))
  call add_drv2_2pd(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm)
!
!  j-k / j-k
!
  d2trm = - const*dampo_nb_tot*dampo_lp*(dampo*d2rdampdrjk2 + 2.0_dp*drdampdrjk*ddampodrjk + rdamp*d2dampo(6))
  call add_drv2_1pd(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk, &
                    deijkdrjk,d2trm,d2r2dx23(1,3))
!
!  Compute derivatives of rblp components w.r.t. neighbour distance
!
  dxblp(1) = xblpsum*(rlp_dist/rblpsum**3)*xblpsum
  dxblp(2) = xblpsum*(rlp_dist/rblpsum**3)*yblpsum
  dxblp(3) = xblpsum*(rlp_dist/rblpsum**3)*zblpsum
!
  dyblp(1) = yblpsum*(rlp_dist/rblpsum**3)*xblpsum
  dyblp(2) = yblpsum*(rlp_dist/rblpsum**3)*yblpsum
  dyblp(3) = yblpsum*(rlp_dist/rblpsum**3)*zblpsum
!
  dzblp(1) = zblpsum*(rlp_dist/rblpsum**3)*xblpsum
  dzblp(2) = zblpsum*(rlp_dist/rblpsum**3)*yblpsum
  dzblp(3) = zblpsum*(rlp_dist/rblpsum**3)*zblpsum
!
  dxblp(1) = dxblp(1) - rlp_dist/rblpsum
  dyblp(2) = dyblp(2) - rlp_dist/rblpsum
  dzblp(3) = dzblp(3) - rlp_dist/rblpsum
!
  dr2alp(1) = xalp*dxblp(1) + yalp*dyblp(1) + zalp*dzblp(1)
  dr2alp(2) = xalp*dxblp(2) + yalp*dyblp(2) + zalp*dzblp(2)
  dr2alp(3) = xalp*dxblp(3) + yalp*dyblp(3) + zalp*dzblp(3)
!
  d2xblp(1,1) = - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*xblpsum*xblpsum + 3.0_dp*(rlp_dist/rblpsum**3)*xblpsum
  d2xblp(2,1) = - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*xblpsum*xblpsum + (rlp_dist/rblpsum**3)*yblpsum
  d2xblp(3,1) = - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*xblpsum*xblpsum + (rlp_dist/rblpsum**3)*zblpsum
  d2xblp(1,2) = - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*yblpsum*xblpsum + (rlp_dist/rblpsum**3)*yblpsum
  d2xblp(2,2) = - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*yblpsum*xblpsum + (rlp_dist/rblpsum**3)*xblpsum
  d2xblp(3,2) = - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*yblpsum*xblpsum
  d2xblp(1,3) = - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*zblpsum*xblpsum + (rlp_dist/rblpsum**3)*zblpsum
  d2xblp(2,3) = - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*zblpsum*xblpsum
  d2xblp(3,3) = - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*zblpsum*xblpsum + (rlp_dist/rblpsum**3)*xblpsum
!
  d2yblp(1,1) = - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*xblpsum*yblpsum + (rlp_dist/rblpsum**3)*yblpsum
  d2yblp(2,1) = - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*xblpsum*yblpsum + (rlp_dist/rblpsum**3)*xblpsum
  d2yblp(3,1) = - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*xblpsum*yblpsum 
  d2yblp(1,2) = - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*yblpsum*yblpsum + (rlp_dist/rblpsum**3)*xblpsum
  d2yblp(2,2) = - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*yblpsum*yblpsum + 3.0_dp*(rlp_dist/rblpsum**3)*yblpsum
  d2yblp(3,2) = - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*yblpsum*yblpsum + (rlp_dist/rblpsum**3)*zblpsum
  d2yblp(1,3) = - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*zblpsum*yblpsum 
  d2yblp(2,3) = - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*zblpsum*yblpsum + (rlp_dist/rblpsum**3)*zblpsum
  d2yblp(3,3) = - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*zblpsum*yblpsum + (rlp_dist/rblpsum**3)*yblpsum
!
  d2zblp(1,1) = - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*xblpsum*zblpsum + (rlp_dist/rblpsum**3)*zblpsum
  d2zblp(2,1) = - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*xblpsum*zblpsum
  d2zblp(3,1) = - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*xblpsum*zblpsum + (rlp_dist/rblpsum**3)*xblpsum
  d2zblp(1,2) = - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*yblpsum*zblpsum
  d2zblp(2,2) = - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*yblpsum*zblpsum + (rlp_dist/rblpsum**3)*zblpsum
  d2zblp(3,2) = - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*yblpsum*zblpsum + (rlp_dist/rblpsum**3)*yblpsum
  d2zblp(1,3) = - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*zblpsum*zblpsum + (rlp_dist/rblpsum**3)*xblpsum
  d2zblp(2,3) = - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*zblpsum*zblpsum + (rlp_dist/rblpsum**3)*yblpsum
  d2zblp(3,3) = - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*zblpsum*zblpsum + 3.0_dp*(rlp_dist/rblpsum**3)*zblpsum
!
  d2r2blpdx2(1) = xalp*d2xblp(1,1) + dxblp(1)*dxblp(1) + yalp*d2yblp(1,1) + dyblp(1)*dyblp(1) + &
                  zalp*d2zblp(1,1) + dzblp(1)*dzblp(1)
  d2r2blpdx2(2) = xalp*d2xblp(2,2) + dxblp(2)*dxblp(2) + yalp*d2yblp(2,2) + dyblp(2)*dyblp(2) + &
                  zalp*d2zblp(2,2) + dzblp(2)*dzblp(2)
  d2r2blpdx2(3) = xalp*d2xblp(3,3) + dxblp(3)*dxblp(3) + yalp*d2yblp(3,3) + dyblp(3)*dyblp(3) + &
                  zalp*d2zblp(3,3) + dzblp(3)*dzblp(3)
  d2r2blpdx2(6) = xalp*d2xblp(2,1) + dxblp(2)*dxblp(1) + yalp*d2yblp(2,1) + dyblp(2)*dyblp(1) + &
                  zalp*d2zblp(2,1) + dzblp(2)*dzblp(1)
  d2r2blpdx2(5) = xalp*d2xblp(3,1) + dxblp(3)*dxblp(1) + yalp*d2yblp(3,1) + dyblp(3)*dyblp(1) + &
                  zalp*d2zblp(3,1) + dzblp(3)*dzblp(1)
  d2r2blpdx2(4) = xalp*d2xblp(3,2) + dxblp(3)*dxblp(2) + yalp*d2yblp(3,2) + dyblp(3)*dyblp(2) + &
                  zalp*d2zblp(3,2) + dzblp(3)*dzblp(2)
!
!  Derivatives of dampo_lp with respect to ralp and rblp
!
  do nj = 1,nnbr(j)
    l = nbrno(nj,j)
    lloc = atom2local(l)
    lllocal = (lloc.ne.0)
!
    indl = 3*(l-1)
    lxf = indl + 1
    lyf = indl + 2
    lzf = indl + 3
    if (lllocal) then
      indl = 3*(lloc-1)
      lx = indl + 1
      ly = indl + 2
      lz = indl + 3
    endif
    xjl = xbnbr(nj,j)
    yjl = ybnbr(nj,j)
    zjl = zbnbr(nj,j)
!
!  j-l/j-l for rblp: NB use of different routine to allow for explicit definition of 1/2 dr2/dadb
!
    call add_drv2_1gpd(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,dr2alp(1),dr2alp(2),dr2alp(3), &
                       dtrmalp,d2trmalp2,d2r2blpdx2)
!
!  Set up terms for mixed ij/jl contributions
!
    dr2alpij(1) = xalp
    dr2alpij(2) = yalp
    dr2alpij(3) = zalp
!
    d2r2dx2mix(1,1) = dxblp(1)
    d2r2dx2mix(2,1) = dxblp(2)
    d2r2dx2mix(3,1) = dxblp(3)
    d2r2dx2mix(1,2) = dyblp(1)
    d2r2dx2mix(2,2) = dyblp(2)
    d2r2dx2mix(3,2) = dyblp(3)
    d2r2dx2mix(1,3) = dzblp(1)
    d2r2dx2mix(2,3) = dzblp(2)
    d2r2dx2mix(3,3) = dzblp(3)
!
!  i-j/j-l for rblp mixed terms: NB use of different routine due to special case
!
    call add_drv2_2sgpd(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,dtrmalp,dr2alpij,dr2alp,d2trmalp2,d2r2dx2mix)
!
!  Mixed terms
!
!  i-j / j-l (ralp)
!
    d2trm = - const*dampo_nb_tot*(dampo*rdamp*d2dampo_lp(2) + dampo*drdampdrij*ddampo_lpdralp + &
              rdamp*ddampodrij*ddampo_lpdralp)
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,dr2alp(1),dr2alp(2),dr2alp(3), &
                      d2trm)
!
!  j-l (ralp)/ i-k
!
    d2trm = - const*dampo_nb_tot*rdamp*ddampodrik*ddampo_lpdralp
    call add_drv2_2pd(ljlocal,lllocal,lilocal,lklocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,dr2alp(1),dr2alp(2),dr2alp(3),xik,yik,zik,d2trm)
!
!  j-l (ralp)/ j-k
!
    d2trm = - const*dampo_nb_tot*(rdamp*ddampodrjk + dampo*drdampdrjk)*ddampo_lpdralp
    call add_drv2_2pd(ljlocal,lllocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,dr2alp(1),dr2alp(2),dr2alp(3),xjk,yjk,zjk,d2trm)
!
    do nk = 1,nj-1
      m = nbrno(nk,j)
      mloc = atom2local(m)
      lmlocal = (mloc.ne.0)
!
      indm = 3*(m-1)
      mxf = indm + 1
      myf = indm + 2
      mzf = indm + 3
      if (lmlocal) then
        indm = 3*(mloc-1)
        mx = indm + 1
        my = indm + 2
        mz = indm + 3
      endif
!
!  j-l/j-m for rblp mixed Cartesian terms: NB use of different routine due to special case
!
      d2r2dx2mix(1,1) = d2r2blpdx2(1)
      d2r2dx2mix(2,1) = d2r2blpdx2(6)
      d2r2dx2mix(3,1) = d2r2blpdx2(5)
      d2r2dx2mix(1,2) = d2r2blpdx2(6)
      d2r2dx2mix(2,2) = d2r2blpdx2(2)
      d2r2dx2mix(3,2) = d2r2blpdx2(4)
      d2r2dx2mix(1,3) = d2r2blpdx2(5)
      d2r2dx2mix(2,3) = d2r2blpdx2(4)
      d2r2dx2mix(3,3) = d2r2blpdx2(3)
!
      call add_drv2_2sgnospd(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                             jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,dtrmalp,dr2alp,dr2alp, &
                             d2trmalp2,d2r2dx2mix)
    enddo
  enddo
!
  if (llpok) then
!
!  Compute A - neighbour of B distance contributions to derivatives
!
    do nj = 1,nnbr(j)
      m = nbrno(nj,j)
      mloc = atom2local(m)
      lmlocal = (mloc.ne.0)
!
      xjm = xbnbr(nj,j)
      yjm = ybnbr(nj,j)
      zjm = zbnbr(nj,j)
!
      xim = xij + xjm
      yim = yij + yjm
      zim = zij + zjm
!
!  Compute product of factors except the current one
!
      dampo_nb_totm1 = 1.0_dp
      do nj2 = 1,nnbr(j)
        if (nj2.ne.nj) then
          dampo_nb_totm1 = dampo_nb_totm1*dampo_nb(nj2)
        endif
      enddo
!
      deijkdrij = - const*rdamp*dampo_nb_totm1*dampo*dampo_lp*ddampo_nbdr(1,nj)   ! (1/r(dE/dr))
      deijkdrim = - const*rdamp*dampo_nb_totm1*dampo*dampo_lp*ddampo_nbdr(2,nj)   ! (1/r(dE/dr))
      deijkdrjm = - const*rdamp*dampo_nb_totm1*dampo*dampo_lp*ddampo_nbdr(3,nj)   ! (1/r(dE/dr))
!
      d2r2dx23m(1,1) = xij*xij
      d2r2dx23m(2,1) = yij*yij
      d2r2dx23m(3,1) = zij*zij
      d2r2dx23m(4,1) = yij*zij
      d2r2dx23m(5,1) = xij*zij
      d2r2dx23m(6,1) = xij*yij
!
      d2r2dx23m(1,2) = xim*xim
      d2r2dx23m(2,2) = yim*yim
      d2r2dx23m(3,2) = zim*zim
      d2r2dx23m(4,2) = yim*zim
      d2r2dx23m(5,2) = xim*zim
      d2r2dx23m(6,2) = xim*yim
!
      d2r2dx23m(1,3) = xjm*xjm
      d2r2dx23m(2,3) = yjm*yjm
      d2r2dx23m(3,3) = zjm*zjm
      d2r2dx23m(4,3) = yjm*zjm
      d2r2dx23m(5,3) = xjm*zjm
      d2r2dx23m(6,3) = xjm*yjm
!
      indm = 3*(m-1)
      mxf = indm + 1
      myf = indm + 2
      mzf = indm + 3
      if (lmlocal) then
        indm = 3*(mloc-1)
        mx = indm + 1
        my = indm + 2
        mz = indm + 3
      endif
!---------------------------------------------------------------------------
!  Contribution from second derivatives of dampo_nb for a single distance  |
!---------------------------------------------------------------------------
!
!  i-j / i-j
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(1,nj)
      call add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                        deijkdrij,d2trm,d2r2dx23m(1,1))
!
!  i-j / i-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(2,nj)
      call add_drv2_2pd(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm)
!
!  i-j / j-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(3,nj)
      call add_drv2_2pd(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm)
!
!  i-m / i-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(4,nj)
      call add_drv2_1pd(lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xim,yim,zim, &
                        deijkdrim,d2trm,d2r2dx23m(1,2))
!
!  i-m / j-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(5,nj)
      call add_drv2_2pd(lilocal,lmlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xim,yim,zim,xjm,yjm,zjm,d2trm)
!
!  j-m / j-m
!
      d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(6,nj)
      call add_drv2_1pd(ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjm,yjm,zjm, &
                        deijkdrjm,d2trm,d2r2dx23m(1,3))
!------------------------------------------------
!  Contribution from dampo_nb and 1 other term  |
!------------------------------------------------
!
!  i-j / i-j
!
      d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(1,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
      call add_drv2_2pd(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm)
!
!  i-j / i-m
!
      d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(2,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
      call add_drv2_2pd(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm)
!
!  i-j / j-m
!
      d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(3,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
      call add_drv2_2pd(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm)
!
!  i-k / i-j
!
      d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(1,nj)*rdamp*ddampodrik
      call add_drv2_2pd(lilocal,lklocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xik,yik,zik,xij,yij,zij,d2trm)
!
!  i-k / i-m
!
      d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(2,nj)*rdamp*ddampodrik
      call add_drv2_2pd(lilocal,lklocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xim,yim,zim,d2trm)
!
!  i-k / j-m
!
      d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(3,nj)*rdamp*ddampodrik
      call add_drv2_2pd(lilocal,lklocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xjm,yjm,zjm,d2trm)
!
!  j-k / i-j
!
      d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(1,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
      call add_drv2_2pd(ljlocal,lklocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjk,yjk,zjk,xij,yij,zij,d2trm)
!
!  j-k / i-m
!
      d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(2,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
      call add_drv2_2pd(ljlocal,lklocal,lilocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xim,yim,zim,d2trm)
!
!  j-k / j-m
!
      d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(3,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
      call add_drv2_2pd(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,d2trm)
!--------------------------------------------------------
!  Mixed derivatives between dampo_lp and dampo_nb_tot  |
!--------------------------------------------------------
!
!  i-j / i-j
!
      d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(1,nj)*ddampo_lpdrij
      call add_drv2_2pd(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm)
!
!  i-m / i-j
!
      d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(2,nj)*ddampo_lpdrij
      call add_drv2_2pd(lilocal,lmlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xim,yim,zim,xij,yij,zij,d2trm)
!
!  j-m / i-j
!
      d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(3,nj)*ddampo_lpdrij
      call add_drv2_2pd(ljlocal,lmlocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjm,yjm,zjm,xij,yij,zij,d2trm)
!
!  i-j / i-j (ralp)
!
      d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(1,nj)*ddampo_lpdralp
      call add_drv2_2pd(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xalp,yalp,zalp,d2trm)
!
!  i-m / i-j (ralp)
!
      d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(2,nj)*ddampo_lpdralp
      call add_drv2_2pd(lilocal,lmlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xim,yim,zim,xalp,yalp,zalp,d2trm)
!
!  j-m / i-j (ralp)
!
      d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(3,nj)*ddampo_lpdralp
      call add_drv2_2pd(ljlocal,lmlocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjm,yjm,zjm,xalp,yalp,zalp,d2trm)
!
      lijmlocal = (lilocal.or.ljlocal.or.lmlocal)
!
!  Loop over dampo_lp j-l contributions
!
      do nk = 1,nnbr(j)
        l = nbrno(nk,j)
        lloc = atom2local(l)
        lllocal = (lloc.ne.0)
!
        if (.not.lijmlocal.and..not.lllocal) cycle
!
        indl = 3*(l-1)
        lxf = indl + 1
        lyf = indl + 2
        lzf = indl + 3
        if (lllocal) then
          indl = 3*(lloc-1)
          lx = indl + 1
          ly = indl + 2
          lz = indl + 3
        endif
!
        xjl = xbnbr(nk,j)
        yjl = ybnbr(nk,j)
        zjl = zbnbr(nk,j)
!
!  j-l (ralp)/ i-j
!
        d2trm = - const*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(1,nj)*ddampo_lpdralp
        call add_drv2_2pd(ljlocal,lllocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,dr2alp(1),dr2alp(2),dr2alp(3), &
                          xij,yij,zij,d2trm)
!
!  j-l (ralp)/ i-m
!
        d2trm = - const*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(2,nj)*ddampo_lpdralp
        call add_drv2_2pd(ljlocal,lllocal,lilocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,dr2alp(1),dr2alp(2),dr2alp(3), &
                          xim,yim,zim,d2trm)
!
!  j-l (ralp)/ j-m
!
        d2trm = - const*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(3,nj)*ddampo_lpdralp
        call add_drv2_2pd(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,dr2alp(1),dr2alp(2),dr2alp(3), &
                          xjm,yjm,zjm,d2trm)
      enddo
!----------------------------------------------------------------------
!  Loop over second neighbour of B for products of first derivatives  |
!----------------------------------------------------------------------
      do nk = 1,nj-1
        n = nbrno(nk,j)
        nloc = atom2local(n)
        lnlocal = (nloc.ne.0)
!
        if (.not.lijmlocal.and..not.lnlocal) cycle
!
        indn = 3*(n-1)
        nxf = indn + 1
        nyf = indn + 2
        nzf = indn + 3
        if (lnlocal) then
          indn = 3*(nloc-1)
          nx = indn + 1
          ny = indn + 2
          nz = indn + 3
        endif
!
        xjn = xbnbr(nk,j)
        yjn = ybnbr(nk,j)
        zjn = zbnbr(nk,j)
!
        xin = xij + xjn
        yin = yij + yjn
        zin = zij + zjn
!
!  Compute product of factors except the current one
!
        dampo_nb_totm2 = 1.0_dp
        do nj2 = 1,nnbr(j)
          if (nj2.ne.nj.and.nj2.ne.nk) then
            dampo_nb_totm2 = dampo_nb_totm2*dampo_nb(nj2)
          endif
        enddo
!
!  i-j / i-j
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(1,nj)*ddampo_nbdr(1,nk)
        call add_drv2_2pd(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm)
!
!  i-j / i-n
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(1,nj)*ddampo_nbdr(2,nk)
        call add_drv2_2pd(lilocal,ljlocal,lilocal,lnlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xij,yij,zij,xin,yin,zin,d2trm)
!
!  i-j / j-n
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(1,nj)*ddampo_nbdr(3,nk)
        call add_drv2_2pd(lilocal,ljlocal,ljlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xij,yij,zij,xjn,yjn,zjn,d2trm)
!
!  i-m / i-j
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(2,nj)*ddampo_nbdr(1,nk)
        call add_drv2_2pd(lilocal,lmlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xim,yim,zim,xij,yij,zij,d2trm)
!
!  i-m / i-n
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(2,nj)*ddampo_nbdr(2,nk)
        call add_drv2_2pd(lilocal,lmlocal,lilocal,lnlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                          ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xim,yim,zim,xin,yin,zin,d2trm)
!
!  i-m / j-n
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(2,nj)*ddampo_nbdr(3,nk)
        call add_drv2_2pd(lilocal,lmlocal,ljlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                          jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xim,yim,zim,xjn,yjn,zjn,d2trm)
!
!  j-m / i-j
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(3,nj)*ddampo_nbdr(1,nk)
        call add_drv2_2pd(ljlocal,lmlocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjm,yjm,zjm,xij,yij,zij,d2trm)
!
!  j-m / i-n
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(3,nj)*ddampo_nbdr(2,nk)
        call add_drv2_2pd(ljlocal,lmlocal,lilocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                          ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xin,yin,zin,d2trm)
!
!  j-m / j-n
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(3,nj)*ddampo_nbdr(3,nk)
        call add_drv2_2pd(ljlocal,lmlocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                          jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xjn,yjn,zjn,d2trm)
      enddo
    enddo
  endif
!
!  Exit point to ensure memory is deallocated
!
  1000 continue
!
  deallocate(d2dampo_nbdr2,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg2pd_rnr','d2dampo_nbdr2')
  deallocate(ddampo_nbdr,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg2pd_rnr','ddampo_nbdr')
  deallocate(dampo_nb,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg2pd_rnr','dampo_nb')
!
#ifdef TRACE
  call trace_out('gfnff_hb2_eg2pd_rnr')
#endif
!
  return
  end
!
  subroutine gfnff_hb2_eg3pd(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                             maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const_in,radij,xkv,ykv,zkv)
!
!  Calculates the hydrogen bonding energy and derivatives for the GFNFF force fields for nhb2
!  Special case version for carbonyls/nitryls with two in plane lone pairs.
!  Distributed memory parallel phonon version.
!
!  On entry : 
!
!  lilocal         = is i local to this node?
!  ljlocal         = is j local to this node?
!  i               = first atom of AB to which k is bonded
!  j               = second atom of AB
!  ix,  iy,  iz    = second derivative indices for i if local
!  ixf, iyf, izf   = second derivative indices for i (global)
!  jx,  jy,  jz    = second derivative indices for j if local
!  jxf, jyf, jzf   = second derivative indices for j (global)
!  ni              = pointer to k in bonding lists of i
!  rij             = distance from i to j (A to B)
!  xij             = x component of vector from i to j (A to B)
!  yij             = y component of vector from i to j (A to B)
!  zij             = z component of vector from i to j (A to B)
!  const_in        = energy factor
!
!   4/22 Created from gfnff_hb2_eg3d
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
  use neighbours
  use numbers,        only : third
  use parallel
  use spatialbo
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  integer(i4), intent(in)                          :: ix
  integer(i4), intent(in)                          :: iy
  integer(i4), intent(in)                          :: iz
  integer(i4), intent(in)                          :: ixf
  integer(i4), intent(in)                          :: iyf
  integer(i4), intent(in)                          :: izf
  integer(i4), intent(in)                          :: j
  integer(i4), intent(in)                          :: jx
  integer(i4), intent(in)                          :: jy
  integer(i4), intent(in)                          :: jz
  integer(i4), intent(in)                          :: jxf
  integer(i4), intent(in)                          :: jyf
  integer(i4), intent(in)                          :: jzf
  integer(i4), intent(in)                          :: ni
  integer(i4), intent(in)                          :: maxnbr
  integer(i4), intent(in)                          :: nnbr(maxat)
  integer(i4), intent(in)                          :: nbrno(maxnbr,maxat)
  logical,     intent(in)                          :: lilocal
  logical,     intent(in)                          :: ljlocal
  real(dp),    intent(in)                          :: const_in
  real(dp),    intent(in)                          :: radij
  real(dp),    intent(in)                          :: rij
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: xkv
  real(dp),    intent(in)                          :: ykv
  real(dp),    intent(in)                          :: zkv
  real(dp),    intent(in)                          :: rbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: xbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: ybnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: zbnbr(maxnbr,maxat)
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: indm
  integer(i4)                                      :: indn
  integer(i4)                                      :: k
  integer(i4)                                      :: kloc
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: l
  integer(i4)                                      :: lloc
  integer(i4)                                      :: lx
  integer(i4)                                      :: ly
  integer(i4)                                      :: lz
  integer(i4)                                      :: lxf
  integer(i4)                                      :: lyf
  integer(i4)                                      :: lzf
  integer(i4)                                      :: m
  integer(i4)                                      :: mloc
  integer(i4)                                      :: mx
  integer(i4)                                      :: my
  integer(i4)                                      :: mz
  integer(i4)                                      :: mxf
  integer(i4)                                      :: myf
  integer(i4)                                      :: mzf
  integer(i4)                                      :: n
  integer(i4)                                      :: nloc
  integer(i4)                                      :: nx
  integer(i4)                                      :: ny
  integer(i4)                                      :: nz
  integer(i4)                                      :: nxf
  integer(i4)                                      :: nyf
  integer(i4)                                      :: nzf
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  integer(i4)                                      :: nk
  integer(i4)                                      :: nl
  integer(i4)                                      :: nm
  integer(i4)                                      :: nt
  integer(i4)                                      :: nt2
  integer(i4)                                      :: nt3
  integer(i4)                                      :: ntorsion
  integer(i4)                                      :: status
  logical                                          :: lijklmlocal
  logical                                          :: lijmlocal
  logical                                          :: lklocal
  logical                                          :: lllocal
  logical                                          :: lmlocal
  logical                                          :: lnlocal
  real(dp)                                         :: const
  real(dp)                                         :: cosa
  real(dp)                                         :: costheta0
  real(dp)                                         :: dcosdrjk
  real(dp)                                         :: dcosdrjl
  real(dp)                                         :: dcosdrkl
  real(dp)                                         :: d2cos(6)
  real(dp)                                         :: d2trm
  real(dp)                                         :: damp
  real(dp)                                         :: dampl
  real(dp)                                         :: ddampl
  real(dp)                                         :: d2dampl
  real(dp)                                         :: damps
  real(dp)                                         :: ddamps
  real(dp)                                         :: d2damps
  real(dp)                                         :: dampo
  real(dp)                                         :: ddampodrij
  real(dp)                                         :: ddampodrik
  real(dp)                                         :: ddampodrjk
  real(dp)                                         :: d2dampo(6)
  real(dp)                                         :: dampo_nbloc
  real(dp)                                         :: dampo_nb_tot
  real(dp)                                         :: dampo_nb_totm1
  real(dp)                                         :: dampo_nb_totm2
  real(dp),    dimension(:),     allocatable, save :: dampo_nb
  real(dp),    dimension(:,:),   allocatable, save :: ddampo_nbdr
  real(dp),    dimension(:,:),   allocatable, save :: d2dampo_nbdr2
  real(dp)                                         :: dfangledcosa
  real(dp)                                         :: d2fangledcosa2
  real(dp)                                         :: d2r2dx23(6,3)
  real(dp)                                         :: d2r2dx23a(6,3)
  real(dp)                                         :: d2r2dx26t(6,6)
  real(dp)                                         :: d2r2dx23m(6,3)
  real(dp)                                         :: d2r2dx23n(6,3)
  real(dp)                                         :: expo
  real(dp)                                         :: dexpodrij
  real(dp)                                         :: dexpodrik
  real(dp)                                         :: dexpodrjk
  real(dp)                                         :: d2expo(6)
  real(dp)                                         :: expo_nb
  real(dp)                                         :: dexpo_nbdrij
  real(dp)                                         :: dexpo_nbdrim
  real(dp)                                         :: dexpo_nbdrjm
  real(dp)                                         :: d2expo_nb(6)
  real(dp)                                         :: eijk
  real(dp)                                         :: deijkdr
  real(dp)                                         :: deijkdrij
  real(dp)                                         :: deijkdrik
  real(dp)                                         :: deijkdrim
  real(dp)                                         :: deijkdrjk
  real(dp)                                         :: deijkdrjm
  real(dp)                                         :: e1d(6)
  real(dp)                                         :: e2d(21)
  real(dp)                                         :: fangle
  real(dp)                                         :: dfangledrjk
  real(dp)                                         :: dfangledrjl
  real(dp)                                         :: dfangledrkl
  real(dp)                                         :: fconst
  real(dp),    dimension(:),     allocatable, save :: ft
  real(dp),    dimension(:,:),   allocatable, save :: dft
  real(dp),    dimension(:,:),   allocatable, save :: d2ft
  real(dp)                                         :: ftloc
  real(dp)                                         :: ftors
  real(dp)                                         :: ftorsm1
  real(dp)                                         :: ftorsm2
  real(dp)                                         :: hbnbcutloc
  real(dp)                                         :: ktheta
  real(dp)                                         :: p_ab
  real(dp)                                         :: p_bh
  real(dp)                                         :: phi0
  real(dp)                                         :: qhdampo
  real(dp)                                         :: rn
  real(dp)                                         :: r2ij
  real(dp)                                         :: rik
  real(dp)                                         :: r2ik
  real(dp)                                         :: rim
  real(dp)                                         :: r2im
  real(dp)                                         :: rsum
  real(dp)                                         :: rsum2
  real(dp)                                         :: rjk
  real(dp)                                         :: r2jk
  real(dp)                                         :: rrjk
  real(dp)                                         :: rjl
  real(dp)                                         :: r2jl
  real(dp)                                         :: rrjl
  real(dp)                                         :: rjm
  real(dp)                                         :: r2jm
  real(dp)                                         :: rkl
  real(dp)                                         :: r2kl
  real(dp)                                         :: rkm
  real(dp)                                         :: r2km
  real(dp)                                         :: rlm
  real(dp)                                         :: ratio1
  real(dp)                                         :: dratio1drij
  real(dp)                                         :: d2ratio1drij2
  real(dp)                                         :: ratio2
  real(dp)                                         :: ratio2_nb
  real(dp)                                         :: ratio3
  real(dp)                                         :: dratio3drij
  real(dp)                                         :: d2ratio3drij2
  real(dp)                                         :: rdamp
  real(dp)                                         :: rkfor
  real(dp)                                         :: drdampdrij
  real(dp)                                         :: drdampdrjk
  real(dp)                                         :: d2rdampdrij2
  real(dp)                                         :: d2rdampdrijdrjk
  real(dp)                                         :: d2rdampdrjk2
  real(dp)                                         :: rijdamp
  real(dp)                                         :: rjkdamp
  real(dp)                                         :: shortcut
  real(dp)                                         :: theta0
  real(dp)                                         :: xik
  real(dp)                                         :: yik
  real(dp)                                         :: zik
  real(dp)                                         :: xim
  real(dp)                                         :: yim
  real(dp)                                         :: zim
  real(dp)                                         :: xin
  real(dp)                                         :: yin
  real(dp)                                         :: zin
  real(dp)                                         :: xjk
  real(dp)                                         :: yjk
  real(dp)                                         :: zjk
  real(dp)                                         :: xjl
  real(dp)                                         :: yjl
  real(dp)                                         :: zjl
  real(dp)                                         :: xjm
  real(dp)                                         :: yjm
  real(dp)                                         :: zjm
  real(dp)                                         :: xjn
  real(dp)                                         :: yjn
  real(dp)                                         :: zjn
  real(dp)                                         :: xkl
  real(dp)                                         :: ykl
  real(dp)                                         :: zkl
  real(dp)                                         :: xkm
  real(dp)                                         :: ykm
  real(dp)                                         :: zkm
  real(dp)                                         :: xkn
  real(dp)                                         :: ykn
  real(dp)                                         :: zkn
  real(dp)                                         :: xlm
  real(dp)                                         :: ylm
  real(dp)                                         :: zlm
  real(dp)                                         :: xln
  real(dp)                                         :: yln
  real(dp)                                         :: zln
#ifdef TRACE
  call trace_in('gfnff_hb2_eg3pd')
#endif
!
  allocate(dampo_nb(maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg3pd','dampo_nb')
  allocate(ddampo_nbdr(3,maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg3pd','ddampo_nbdr')
  allocate(d2dampo_nbdr2(6,maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg3pd','d2dampo_nbdr2')
!
  k = nbrno(ni,i)   ! H
  kloc = atom2local(k)
  lklocal = (kloc.ne.0)
!
  p_bh = (1.0_dp + gfnff_hbabmix)*gfnff_autoangs**3
  p_ab = - gfnff_hbabmix*gfnff_autoangs**3
!
  r2ij = rij**2
!
!  A-H distance (i-k)
!
  xik = xbnbr(ni,i)
  yik = ybnbr(ni,i)
  zik = zbnbr(ni,i)
  rik = rbnbr(ni,i)
  r2ik = rik**2
!
!  B-H distance (j-k)
!
  xjk = xik - xij
  yjk = yik - yij
  zjk = zik - zij
  r2jk = xjk**2 + yjk**2 + zjk**2
  rjk = sqrt(r2jk)
!
  rsum = rik + rjk + 1.d-12
!
!  Out-of-line damping : A-H...B
!
  expo = (gfnff_hb_a_cut/radij)*(rsum/rij - 1.0_dp)
  if (expo.gt.15.0_dp) goto 1000 
  ratio2 = exp(expo)
  dampo = 2.0_dp/(1.0_dp + ratio2)
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
!  Loop over neighbours of B
!
  if (nat(j).eq.7.and.nnbr(j).eq.1) then
    hbnbcutloc = 2.0_dp*gfnff_autoangs
  else
    hbnbcutloc = gfnff_hb_nb_cut
  endif
!
  dampo_nb_tot = 1.0_dp
  do nj = 1,nnbr(j)
!
!  Compute A - neighbour of B distance
!
    xim = xij + xbnbr(nj,j)
    yim = yij + ybnbr(nj,j)
    zim = zij + zbnbr(nj,j)
    r2im = xim**2 + yim**2 + zim**2
    rim = sqrt(r2im)
    rjm = rbnbr(nj,j)
    r2jm = rjm**2
!
    rsum2 = rim + rjm + 1.0d-12
    expo_nb = (hbnbcutloc/radij)*(rsum2/rij - 1.0_dp)
    ratio2_nb = exp(-expo_nb)
    dampo_nbloc = (2.0_dp/(1.0_dp + ratio2_nb))
    dampo_nb(nj) = dampo_nbloc - 1.0_dp
    dampo_nb_tot = dampo_nb_tot*dampo_nb(nj)
!
    dexpo_nbdrij = (hbnbcutloc/radij)*ratio2_nb*rsum2/rij**3
    dexpo_nbdrim = - (hbnbcutloc/radij)*ratio2_nb/(rij*rim)
    dexpo_nbdrjm = - (hbnbcutloc/radij)*ratio2_nb/(rij*rjm)
!
    ddampo_nbdr(1,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrij
    ddampo_nbdr(2,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrim
    ddampo_nbdr(3,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrjm
!
    d2expo_nb(1) = dexpo_nbdrij*(hbnbcutloc/radij)*rsum2/rij**3 - 3.0_dp*dexpo_nbdrij/r2ij
    d2expo_nb(2) = - dexpo_nbdrij*(hbnbcutloc/radij)/(rij*rim) - dexpo_nbdrim/r2ij
    d2expo_nb(3) = - dexpo_nbdrij*(hbnbcutloc/radij)/(rij*rjm) - dexpo_nbdrjm/r2ij
    d2expo_nb(4) = - dexpo_nbdrim*(hbnbcutloc/radij)/(rij*rim) - dexpo_nbdrim/r2im
    d2expo_nb(5) = - dexpo_nbdrim*(hbnbcutloc/radij)/(rij*rjm)
    d2expo_nb(6) = - dexpo_nbdrjm*(hbnbcutloc/radij)/(rij*rjm) - dexpo_nbdrjm/r2jm
!
    d2dampo_nbdr2(1,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrij*dexpo_nbdrij - d2expo_nb(1))
    d2dampo_nbdr2(2,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrij*dexpo_nbdrim - d2expo_nb(2))
    d2dampo_nbdr2(3,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrij*dexpo_nbdrjm - d2expo_nb(3))
    d2dampo_nbdr2(4,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrim*dexpo_nbdrim - d2expo_nb(4))
    d2dampo_nbdr2(5,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrim*dexpo_nbdrjm - d2expo_nb(5))
    d2dampo_nbdr2(6,nj) = 0.5_dp*dampo_nbloc*dampo_nbloc*(dampo_nbloc*dexpo_nbdrjm*dexpo_nbdrjm - d2expo_nb(6))
  enddo
!
!  Long-range damping
!
  ratio1 = (r2ij/gfnff_hb_long_cut)**gfnff_hb_alp
  dampl = 1.0_dp/(1.0_dp + ratio1)
!
!  Short-range damping
!
  shortcut = gfnff_hb_short_cut*radij
  ratio3 = (shortcut/r2ij)**gfnff_hb_alp
  damps = 1.0_dp/(1.0_dp + ratio3)
!
  damp = damps*dampl
  rijdamp = damp*(p_ab/r2ij/rij)
  rjkdamp = damp*(p_bh/r2jk/rjk)
  rdamp   = rjkdamp + rijdamp
!
  qhdampo = dampo*dampo_nb_tot
!
  const = const_in*gfnff_hb_ABq(k)
!--------------------
!  Torsional factor |
!--------------------
!
!  Set number of torsions:
!    - since j is a O= then the number of neighbours, l (C/N) must be one
!    - number of torsions is therefore the number of neighbours of l - 1
!
  l = nbrno(1,j)
  lloc = atom2local(l)
  lllocal = (lloc.ne.0)
!
  indl = 3*(l-1)
  lxf = indl + 1
  lyf = indl + 2
  lzf = indl + 3
  if (lllocal) then
    indl = 3*(lloc-1)
    lx = indl + 1
    ly = indl + 2
    lz = indl + 3
  endif
!
!  Distances for l
!
  xjl = xbnbr(1,j)
  yjl = ybnbr(1,j)
  zjl = zbnbr(1,j)
  rjl = rbnbr(1,j)
  r2jl = rjl**2
!
  xkl = xjl - xjk
  ykl = yjl - yjk
  zkl = zjl - zjk
  r2kl = xkl**2 + ykl**2 + zkl**2
  rkl = sqrt(r2kl)
!
!  Set parameters
!
  rn = 2.0_dp
  phi0 = 0.5_dp*pi
  rkfor = 0.5_dp*(1.0_dp - gfnff_tors_hb)
!
!  Set up storage for torsions
!
  ntorsion = nnbr(l) - 1
  allocate(ft(ntorsion),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg3pd','ft')
  allocate(dft(6,ntorsion),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg3pd','dft')
  allocate(d2ft(21,ntorsion),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg3pd','d2ft')
!
!  Compute torsional factor
!
  ftors = 1.0_dp
  ntorsion = 0
  do nl = 1,nnbr(l)
!
!  Exclude case where neighbour is j
!
    m = nbrno(nl,l)
    if (m.eq.j) cycle
    ntorsion = ntorsion + 1
!
!  Compute distances for m
!
    xkm = xbnbr(nl,l) + xkl
    ykm = ybnbr(nl,l) + ykl
    zkm = zbnbr(nl,l) + zkl
    r2km = xkm**2 + ykm**2 + zkm**2
    rkm = sqrt(r2km)
!
    rlm = rbnbr(nl,l)
!
    xjm = xbnbr(nl,l) + xbnbr(1,j)
    yjm = ybnbr(nl,l) + ybnbr(1,j)
    zjm = zbnbr(nl,l) + zbnbr(1,j)
    r2jm = xjm**2 + yjm**2 + zjm**2
    rjm = sqrt(r2jm)
!
    call gfnff_torsion(1_i4,rkfor,rn,phi0,rjk,rkl,rkm,rjl,rjm,rlm,ftloc,e1d,e2d,.true.,.true.)
    ftloc = ftloc + gfnff_tors_hb
!
    ft(ntorsion) = ftloc
    ftors = ftors*ftloc
    dft(1:6,ntorsion) = e1d(1:6)
    d2ft(1:21,ntorsion) = e2d(1:21)
  enddo
!-----------------
!  Angle factor  |
!-----------------
!
!  Set parameters
!
  theta0 = 2.0_dp*pi*third
  costheta0 = cos(theta0)
  ktheta = (1.0_dp - gfnff_bend_hb)/(1.0_dp - costheta0)**2
!
!  Compute angle k-j-l
!
  cosa = 0.5_dp*(r2jk + r2jl - r2kl)/(rjk*rjl)
  cosa = dble(min(1.0_dp,max(-1.0_dp,cosa)))
!
  fangle = 1.0_dp - ktheta*(cosa - costheta0)**2
!-----------
!  Energy  |
!-----------
  fconst = const*fangle*ftors
  eijk =  - fconst*rdamp*qhdampo
!----------------------
!  First derivatives  |
!----------------------
  dexpodrij = - (gfnff_hb_a_cut/radij)*ratio2*rsum/rij**3
  dexpodrik = (gfnff_hb_a_cut/radij)*ratio2/(rij*rik)
  dexpodrjk = (gfnff_hb_a_cut/radij)*ratio2/(rij*rjk)
!
  ddampodrij = - 0.5_dp*dampo*dampo*dexpodrij
  ddampodrik = - 0.5_dp*dampo*dampo*dexpodrik
  ddampodrjk = - 0.5_dp*dampo*dampo*dexpodrjk
!--------------------------------------
!  Derivatives of damping and qhdampo  |
!--------------------------------------
  dratio1drij = 2.0_dp*gfnff_hb_alp*ratio1/r2ij
  dratio3drij = - 2.0_dp*gfnff_hb_alp*ratio3/r2ij
  ddampl = - dampl*dampl*dratio1drij
  ddamps = - damps*damps*dratio3drij
  drdampdrij = (ddamps*dampl + damps*ddampl)*(p_ab/r2ij/rij + p_bh/r2jk/rjk) - 3.0_dp*rijdamp/r2ij
  drdampdrjk = - 3.0_dp*rjkdamp/r2jk
!
  deijkdrij = - fconst*rdamp*dampo_nb_tot*ddampodrij - fconst*qhdampo*drdampdrij
  deijkdrik = - fconst*rdamp*dampo_nb_tot*ddampodrik
  deijkdrjk = - fconst*rdamp*dampo_nb_tot*ddampodrjk - fconst*qhdampo*drdampdrjk
!
  d2r2dx23(1,1) = xij*xij
  d2r2dx23(2,1) = yij*yij
  d2r2dx23(3,1) = zij*zij
  d2r2dx23(4,1) = yij*zij
  d2r2dx23(5,1) = xij*zij
  d2r2dx23(6,1) = xij*yij
!
  d2r2dx23(1,2) = xik*xik
  d2r2dx23(2,2) = yik*yik
  d2r2dx23(3,2) = zik*zik
  d2r2dx23(4,2) = yik*zik
  d2r2dx23(5,2) = xik*zik
  d2r2dx23(6,2) = xik*yik
!
  d2r2dx23(1,3) = xjk*xjk
  d2r2dx23(2,3) = yjk*yjk
  d2r2dx23(3,3) = zjk*zjk
  d2r2dx23(4,3) = yjk*zjk
  d2r2dx23(5,3) = xjk*zjk
  d2r2dx23(6,3) = xjk*yjk
!-----------------------
!  Second derivatives  |
!-----------------------
  d2expo(1) = - dexpodrij*(gfnff_hb_a_cut/radij)*rsum/rij**3 - 3.0_dp*dexpodrij/r2ij
  d2expo(2) = dexpodrij*(gfnff_hb_a_cut/radij)/(rij*rik) - dexpodrik/r2ij
  d2expo(3) = dexpodrij*(gfnff_hb_a_cut/radij)/(rij*rjk) - dexpodrjk/r2ij
  d2expo(4) = dexpodrik*(gfnff_hb_a_cut/radij)/(rij*rik) - dexpodrik/r2ik
  d2expo(5) = dexpodrik*(gfnff_hb_a_cut/radij)/(rij*rjk)
  d2expo(6) = dexpodrjk*(gfnff_hb_a_cut/radij)/(rij*rjk) - dexpodrjk/r2jk
!
  d2dampo(1) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrij - d2expo(1))
  d2dampo(2) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrik - d2expo(2))
  d2dampo(3) = 0.5_dp*dampo*dampo*(dampo*dexpodrij*dexpodrjk - d2expo(3))
  d2dampo(4) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrik - d2expo(4))
  d2dampo(5) = 0.5_dp*dampo*dampo*(dampo*dexpodrik*dexpodrjk - d2expo(5))
  d2dampo(6) = 0.5_dp*dampo*dampo*(dampo*dexpodrjk*dexpodrjk - d2expo(6))
!
  d2ratio1drij2 = 4.0_dp*gfnff_hb_alp*ratio1*(gfnff_hb_alp - 1.0_dp)/r2ij**2
  d2ratio3drij2 = 4.0_dp*gfnff_hb_alp*ratio3*(gfnff_hb_alp + 1.0_dp)/r2ij**2
  d2dampl = - dampl*dampl*d2ratio1drij2 + 2.0_dp*dampl*dampl*dampl*dratio1drij**2
  d2damps = - damps*damps*d2ratio3drij2 + 2.0_dp*damps*damps*damps*dratio3drij**2
  d2rdampdrij2 = (d2damps*dampl + 2.0_dp*ddamps*ddampl + damps*d2dampl)*(p_ab/r2ij/rij + p_bh/r2jk/rjk) + &
                 15.0_dp*rijdamp/r2ij**2 - 6.0_dp*(ddamps*dampl + damps*ddampl)*p_ab/r2ij/r2ij/rij
  d2rdampdrijdrjk = - 3.0_dp*(ddamps*dampl + damps*ddampl)*p_bh/r2jk/r2jk/rjk
  d2rdampdrjk2 = 15.0_dp*rjkdamp/r2jk**2
!
!  i-j / i-j
!
  d2trm = - fconst*dampo_nb_tot*(dampo*d2rdampdrij2 + 2.0_dp*drdampdrij*ddampodrij + rdamp*d2dampo(1))
  call add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                    d2trm,d2r2dx23(1,1))
!
!  i-j / i-k
!
  d2trm = - fconst*dampo_nb_tot*(drdampdrij*ddampodrik + rdamp*d2dampo(2))
  call add_drv2_2pd(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm)
!
!  i-j / j-k
!
  d2trm = - fconst*dampo_nb_tot*(dampo*d2rdampdrijdrjk + drdampdrij*ddampodrjk + drdampdrjk*ddampodrij + &
            rdamp*d2dampo(3))
  call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm)
!
!  i-k / i-k
!
  d2trm = - fconst*dampo_nb_tot*rdamp*d2dampo(4)
  call add_drv2_1pd(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,deijkdrik, &
                    d2trm,d2r2dx23(1,2))
!
!  i-k / j-k
!
  d2trm = - fconst*dampo_nb_tot*(drdampdrjk*ddampodrik + rdamp*d2dampo(5))
  call add_drv2_2pd(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm)
!
!  j-k / j-k
!
  d2trm = - fconst*dampo_nb_tot*(dampo*d2rdampdrjk2 + 2.0_dp*drdampdrjk*ddampodrjk + rdamp*d2dampo(6))
  call add_drv2_1pd(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,deijkdrjk, &
                    d2trm,d2r2dx23(1,3))
!--------------------------
!  Derivatives of fangle  |
!--------------------------
  dfangledcosa = - 2.0_dp*ktheta*(cosa - costheta0)
!
  rrjk = 1.0_dp/rjk
  rrjl = 1.0_dp/rjl
!
  dcosdrjk = rrjk*rrjl - cosa*rrjk**2
  dcosdrjl = rrjk*rrjl - cosa*rrjl**2
  dcosdrkl = - rrjk*rrjl
!
!  j-k
!
  dfangledrjk = - const*ftors*rdamp*qhdampo*dfangledcosa*dcosdrjk      ! (1/r(dEdr))
!
!  j-l
!
  dfangledrjl = - const*ftors*rdamp*qhdampo*dfangledcosa*dcosdrjl      ! (1/r(dEdr))
!
!  k-l
!
  dfangledrkl = - const*ftors*rdamp*qhdampo*dfangledcosa*dcosdrkl      ! (1/r(dEdr))
!
  d2r2dx23a(1,1) = xjk*xjk
  d2r2dx23a(2,1) = yjk*yjk
  d2r2dx23a(3,1) = zjk*zjk
  d2r2dx23a(4,1) = yjk*zjk
  d2r2dx23a(5,1) = xjk*zjk
  d2r2dx23a(6,1) = xjk*yjk
!
  d2r2dx23a(1,2) = xjl*xjl
  d2r2dx23a(2,2) = yjl*yjl
  d2r2dx23a(3,2) = zjl*zjl
  d2r2dx23a(4,2) = yjl*zjl
  d2r2dx23a(5,2) = xjl*zjl
  d2r2dx23a(6,2) = xjl*yjl
!
  d2r2dx23a(1,3) = xkl*xkl
  d2r2dx23a(2,3) = ykl*ykl
  d2r2dx23a(3,3) = zkl*zkl
  d2r2dx23a(4,3) = ykl*zkl
  d2r2dx23a(5,3) = xkl*zkl
  d2r2dx23a(6,3) = xkl*ykl
!
!  Second derivatives of angle terms
!
  d2fangledcosa2 = - 2.0_dp*ktheta
!
  d2cos(1) = - rrjk*rrjk*(dcosdrjk + rrjk*rrjl - 2.0_dp*cosa*rrjk**2)
  d2cos(2) = - rrjk*rrjl**3 - rrjk*rrjk*dcosdrjl
  d2cos(3) = rrjl*rrjk**3
  d2cos(4) = - rrjl*rrjl*(dcosdrjl + rrjk*rrjl - 2.0_dp*cosa*rrjl**2)
  d2cos(5) = rrjk*rrjl**3
  d2cos(6) = 0.0_dp
!
!  j-k / j-k
!
  d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(1) + d2fangledcosa2*dcosdrjk*dcosdrjk)
  call add_drv2_1pd(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,dfangledrjk, &
                    d2trm,d2r2dx23a(1,1))
!
!  j-k / j-l
!
  d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(2) + d2fangledcosa2*dcosdrjk*dcosdrjl)
  call add_drv2_2pd(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                    jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm)
!
!  j-k / k-l
!
  d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(3) + d2fangledcosa2*dcosdrjk*dcosdrkl)
  call add_drv2_2pd(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                    kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm)
!
!  j-l / j-l
!
  d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(4) + d2fangledcosa2*dcosdrjl*dcosdrjl)
  call add_drv2_1pd(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,dfangledrjl, &
                    d2trm,d2r2dx23a(1,2))
!
!  j-l / k-l
!
  d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(5) + d2fangledcosa2*dcosdrjl*dcosdrkl)
  call add_drv2_2pd(ljlocal,lllocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                    kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xkl,ykl,zkl,d2trm)
!
!  k-l / k-l
!
  d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(6) + d2fangledcosa2*dcosdrkl*dcosdrkl)
  call add_drv2_1pd(lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,dfangledrkl, &
                    d2trm,d2r2dx23a(1,3))
!--------------------------------------------------------
!  Second derivatives of fangle mixed with other terms  |
!--------------------------------------------------------
!
!  i-j / j-k
!
  d2trm = - const*ftors*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*dfangledcosa*dcosdrjk
  call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm)
!
!  i-j / j-l
!
  d2trm = - const*ftors*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*dfangledcosa*dcosdrjl
  call add_drv2_2pd(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm)
!
!  i-j / k-l
!
  d2trm = - const*ftors*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*dfangledcosa*dcosdrkl
  call add_drv2_2pd(lilocal,ljlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                    kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xkl,ykl,zkl,d2trm)
!
!  i-k / j-k
!
  d2trm = - const*ftors*dampo_nb_tot*rdamp*ddampodrik*dfangledcosa*dcosdrjk
  call add_drv2_2pd(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm)
!
!  i-k / j-l
!
  d2trm = - const*ftors*dampo_nb_tot*rdamp*ddampodrik*dfangledcosa*dcosdrjl
  call add_drv2_2pd(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                    jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xjl,yjl,zjl,d2trm)
!
!  i-k / k-l
!
  d2trm = - const*ftors*dampo_nb_tot*rdamp*ddampodrik*dfangledcosa*dcosdrkl
  call add_drv2_2pd(lilocal,lklocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                    kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xkl,ykl,zkl,d2trm)
!
!  j-k / j-k
!
  d2trm = - const*ftors*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*dfangledcosa*dcosdrjk
  call add_drv2_2pd(ljlocal,lklocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                    jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,xjk,yjk,zjk,d2trm)
!
!  j-k / j-l
!
  d2trm = - const*ftors*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*dfangledcosa*dcosdrjl
  call add_drv2_2pd(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                    jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm)
!
!  j-k / k-l
!
  d2trm = - const*ftors*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*dfangledcosa*dcosdrkl
  call add_drv2_2pd(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                    kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm)
!
!  Compute A - neighbour of B distance contributions to derivatives
!
  do nj = 1,nnbr(j)
    m = nbrno(nj,j)
    mloc = atom2local(m)
    lmlocal = (mloc.ne.0)
!
    xjm = xbnbr(nj,j)
    yjm = ybnbr(nj,j)
    zjm = zbnbr(nj,j)
!
    xim = xij + xjm
    yim = yij + yjm
    zim = zij + zjm
!
!  Compute product of factors except the current one
!
    dampo_nb_totm1 = 1.0_dp
    do nj2 = 1,nnbr(j)
      if (nj2.ne.nj) then
        dampo_nb_totm1 = dampo_nb_totm1*dampo_nb(nj2)
      endif
    enddo
!
    deijkdrij = - fconst*rdamp*dampo_nb_totm1*dampo*ddampo_nbdr(1,nj)   ! (1/r(dE/dr))
    deijkdrim = - fconst*rdamp*dampo_nb_totm1*dampo*ddampo_nbdr(2,nj)   ! (1/r(dE/dr))
    deijkdrjm = - fconst*rdamp*dampo_nb_totm1*dampo*ddampo_nbdr(3,nj)   ! (1/r(dE/dr))
!
    d2r2dx23m(1,1) = xij*xij
    d2r2dx23m(2,1) = yij*yij
    d2r2dx23m(3,1) = zij*zij
    d2r2dx23m(4,1) = yij*zij
    d2r2dx23m(5,1) = xij*zij
    d2r2dx23m(6,1) = xij*yij
!
    d2r2dx23m(1,2) = xim*xim
    d2r2dx23m(2,2) = yim*yim
    d2r2dx23m(3,2) = zim*zim
    d2r2dx23m(4,2) = yim*zim
    d2r2dx23m(5,2) = xim*zim
    d2r2dx23m(6,2) = xim*yim
!
    d2r2dx23m(1,3) = xjm*xjm
    d2r2dx23m(2,3) = yjm*yjm
    d2r2dx23m(3,3) = zjm*zjm
    d2r2dx23m(4,3) = yjm*zjm
    d2r2dx23m(5,3) = xjm*zjm
    d2r2dx23m(6,3) = xjm*yjm
!
    indm = 3*(m-1)
    mxf = indm + 1
    myf = indm + 2
    mzf = indm + 3
    if (lmlocal) then
      indm = 3*(mloc-1)
      mx = indm + 1
      my = indm + 2
      mz = indm + 3
    endif
!
    lijmlocal = (lilocal.or.ljlocal.or.lmlocal)
!---------------------------------------------------------------------------
!  Contribution from second derivatives of dampo_nb for a single distance  |
!---------------------------------------------------------------------------
!
!  i-j / i-j
!
    d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(1,nj)
    call add_drv2_1pd(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                      d2trm,d2r2dx23m(1,1))
!
!  i-j / i-m
!
    d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(2,nj)
    call add_drv2_2pd(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm)
!
!  i-j / j-m
!
    d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(3,nj)
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm)
!
!  i-m / i-m
!
    d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(4,nj)
    call add_drv2_1pd(lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xim,yim,zim,deijkdrim, &
                      d2trm,d2r2dx23m(1,2))
!
!  i-m / j-m
!
    d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(5,nj)
    call add_drv2_2pd(lilocal,lmlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xim,yim,zim,xjm,yjm,zjm,d2trm)
!
!  j-m / j-m
!
    d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(6,nj)
    call add_drv2_1pd(ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjm,yjm,zjm,deijkdrjm, &
                      d2trm,d2r2dx23m(1,3))
!------------------------------------------------
!  Contribution from dampo_nb and 1 other term  |
!------------------------------------------------
!
!  i-j / i-j
!
    d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(1,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
    call add_drv2_2pd(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm)
!
!  i-j / i-m
!
    d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(2,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
    call add_drv2_2pd(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm)
!
!  i-j / j-m
!
    d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(3,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm)
!
!  i-k / i-j
!
    d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(1,nj)*rdamp*ddampodrik
    call add_drv2_2pd(lilocal,lklocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xik,yik,zik,xij,yij,zij,d2trm)
!
!  i-k / i-m
!
    d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(2,nj)*rdamp*ddampodrik
    call add_drv2_2pd(lilocal,lklocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xim,yim,zim,d2trm)
!
!  i-k / j-m
!
    d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(3,nj)*rdamp*ddampodrik
    call add_drv2_2pd(lilocal,lklocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xjm,yjm,zjm,d2trm)
!
!  j-k / i-j
!
    d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(1,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
    call add_drv2_2pd(ljlocal,lklocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjk,yjk,zjk,xij,yij,zij,d2trm)
!
!  j-k / i-m
!
    d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(2,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
    call add_drv2_2pd(ljlocal,lklocal,lilocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xim,yim,zim,d2trm)
!
!  j-k / j-m
!
    d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(3,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,d2trm)
!------------------------------------------
!  Contribution from dampo_nb and fangle  |
!------------------------------------------
!
!  i-j / j-k
!
    d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(1,nj)*dfangledcosa*dcosdrjk
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm)
!
!  i-j / j-l
!
    d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(1,nj)*dfangledcosa*dcosdrjl
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm)
!
!  i-j / k-l
!
    d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(1,nj)*dfangledcosa*dcosdrkl
    call add_drv2_2pd(lilocal,ljlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xkl,ykl,zkl,d2trm)
!
!  i-m / j-k
!
    d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(2,nj)*dfangledcosa*dcosdrjk
    call add_drv2_2pd(lilocal,lmlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xim,yim,zim,xjk,yjk,zjk,d2trm)
!
!  i-m / j-l
!
    d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(2,nj)*dfangledcosa*dcosdrjl
    call add_drv2_2pd(lilocal,lmlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xim,yim,zim,xjl,yjl,zjl,d2trm)
!
!  i-m / k-l
!
    d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(2,nj)*dfangledcosa*dcosdrkl
    call add_drv2_2pd(lilocal,lmlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xim,yim,zim,xkl,ykl,zkl,d2trm)
!
!  j-m / j-k
!
    d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(3,nj)*dfangledcosa*dcosdrjk
    call add_drv2_2pd(ljlocal,lmlocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjm,yjm,zjm,xjk,yjk,zjk,d2trm)
!
!  j-m / j-l
!
    d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(3,nj)*dfangledcosa*dcosdrjl
    call add_drv2_2pd(ljlocal,lmlocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjm,yjm,zjm,xjl,yjl,zjl,d2trm)
!
!  j-m / k-l
!
    d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(3,nj)*dfangledcosa*dcosdrkl
    call add_drv2_2pd(ljlocal,lmlocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjm,yjm,zjm,xkl,ykl,zkl,d2trm)
!
!  Loop over second neighbour of B for products of first derivatives
!
    do nk = 1,nj-1
      n = nbrno(nk,j)
      nloc = atom2local(n)
      lnlocal = (nloc.ne.0)
!
      if (.not.lijmlocal.and..not.lnlocal) cycle
!
      indn = 3*(n-1)
      nxf = indn + 1
      nyf = indn + 2
      nzf = indn + 3
      if (lnlocal) then
        indn = 3*(nloc-1)
        nx = indn + 1
        ny = indn + 2
        nz = indn + 3
      endif
!
      xjn = xbnbr(nk,j)
      yjn = ybnbr(nk,j)
      zjn = zbnbr(nk,j)
!
      xin = xij + xjn
      yin = yij + yjn
      zin = zij + zjn
!
!  Compute product of factors except the current one
!
      dampo_nb_totm2 = 1.0_dp
      do nj2 = 1,nnbr(j)
        if (nj2.ne.nj.and.nj2.ne.nk) then
          dampo_nb_totm2 = dampo_nb_totm2*dampo_nb(nj2)
        endif
      enddo
!
!  i-j / i-j
!
      d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(1,nk)
      call add_drv2_2pd(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm)
!
!  i-j / i-n
!
      d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(2,nk)
      call add_drv2_2pd(lilocal,ljlocal,lilocal,lnlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xij,yij,zij,xin,yin,zin,d2trm)
!
!  i-j / j-n
!
      d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(3,nk)
      call add_drv2_2pd(lilocal,ljlocal,ljlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xij,yij,zij,xjn,yjn,zjn,d2trm)
!
!  i-m / i-j
!
      d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(1,nk)
      call add_drv2_2pd(lilocal,lmlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xim,yim,zim,xij,yij,zij,d2trm)
!
!  i-m / i-n
!
      d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(2,nk)
      call add_drv2_2pd(lilocal,lmlocal,lilocal,lnlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                        ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xim,yim,zim,xin,yin,zin,d2trm)
!
!  i-m / j-n
!
      d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(3,nk)
      call add_drv2_2pd(lilocal,lmlocal,ljlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xim,yim,zim,xjn,yjn,zjn,d2trm)
!
!  j-m / i-j
!
      d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(1,nk)
      call add_drv2_2pd(ljlocal,lmlocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjm,yjm,zjm,xij,yij,zij,d2trm)
!
!  j-m / i-n
!
      d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(2,nk)
      call add_drv2_2pd(ljlocal,lmlocal,lilocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xin,yin,zin,d2trm)
!
!  j-m / j-n
!
      d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(3,nk)
      call add_drv2_2pd(ljlocal,lmlocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xjn,yjn,zjn,d2trm)
    enddo
  enddo
!---------------------------
!  Derivatives of torsion  |
!---------------------------
  nt = 0
  do nl = 1,nnbr(l)
!
!  Exclude case where neighbour is j
!
    m = nbrno(nl,l)
    if (m.eq.j) cycle
    nt = nt + 1
!
    mloc = atom2local(m)
    lmlocal = (mloc.ne.0)
!
!  Compute distances for m
!
    xkm = xbnbr(nl,l) + xkl
    ykm = ybnbr(nl,l) + ykl
    zkm = zbnbr(nl,l) + zkl
!
    xjm = xbnbr(nl,l) + xbnbr(1,j)
    yjm = ybnbr(nl,l) + ybnbr(1,j)
    zjm = zbnbr(nl,l) + zbnbr(1,j)
!
    xlm = xbnbr(nl,l)
    ylm = ybnbr(nl,l)
    zlm = zbnbr(nl,l)
!
!  Compute product of factors except the current one
!
    ftorsm1 = 1.0_dp
    do nt2 = 1,ntorsion
      if (nt2.ne.nt) then
        ftorsm1 = ftorsm1*ft(nt2)
      endif
    enddo
    deijkdr = - const*fangle*rdamp*qhdampo*ftorsm1
!
    e1d(1:6) = dft(1:6,nt)*deijkdr
!
    d2r2dx26t(1,1) = xjk*xjk
    d2r2dx26t(2,1) = yjk*yjk
    d2r2dx26t(3,1) = zjk*zjk
    d2r2dx26t(4,1) = yjk*zjk
    d2r2dx26t(5,1) = xjk*zjk
    d2r2dx26t(6,1) = xjk*yjk
!
    d2r2dx26t(1,2) = xkl*xkl
    d2r2dx26t(2,2) = ykl*ykl
    d2r2dx26t(3,2) = zkl*zkl
    d2r2dx26t(4,2) = ykl*zkl
    d2r2dx26t(5,2) = xkl*zkl
    d2r2dx26t(6,2) = xkl*ykl
!
    d2r2dx26t(1,3) = xkm*xkm
    d2r2dx26t(2,3) = ykm*ykm
    d2r2dx26t(3,3) = zkm*zkm
    d2r2dx26t(4,3) = ykm*zkm
    d2r2dx26t(5,3) = xkm*zkm
    d2r2dx26t(6,3) = xkm*ykm
!
    d2r2dx26t(1,4) = xjl*xjl
    d2r2dx26t(2,4) = yjl*yjl
    d2r2dx26t(3,4) = zjl*zjl
    d2r2dx26t(4,4) = yjl*zjl
    d2r2dx26t(5,4) = xjl*zjl
    d2r2dx26t(6,4) = xjl*yjl
!
    d2r2dx26t(1,5) = xjm*xjm
    d2r2dx26t(2,5) = yjm*yjm
    d2r2dx26t(3,5) = zjm*zjm
    d2r2dx26t(4,5) = yjm*zjm
    d2r2dx26t(5,5) = xjm*zjm
    d2r2dx26t(6,5) = xjm*yjm
!
    d2r2dx26t(1,6) = xlm*xlm
    d2r2dx26t(2,6) = ylm*ylm
    d2r2dx26t(3,6) = zlm*zlm
    d2r2dx26t(4,6) = ylm*zlm
    d2r2dx26t(5,6) = xlm*zlm
    d2r2dx26t(6,6) = xlm*ylm
!-----------------------------------
!  Second derivatives of torsions  |
!-----------------------------------
    e2d(1:21) = d2ft(1:21,nt)*deijkdr
    indm = 3*(m-1)
    mxf = indm + 1
    myf = indm + 2
    mzf = indm + 3
    if (lmlocal) then
      indm = 3*(mloc-1)
      mx = indm + 1
      my = indm + 2
      mz = indm + 3
    endif
!
!  j-k / j-k
!
    call add_drv2_1pd(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,e1d(1), &
                      e2d(1),d2r2dx26t(1,1))
!
!  j-k / k-l
!
    call add_drv2_2pd(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,e2d(2))
!
!  j-k / k-m
!
    call add_drv2_2pd(ljlocal,lklocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xkm,ykm,zkm,e2d(3))
!
!  j-k / j-l
!
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,e2d(4))
!
!  j-k / j-m
!
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,e2d(5))
!
!  j-k / l-m
!
    call add_drv2_2pd(ljlocal,lklocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xlm,ylm,zlm,e2d(6))
!
!  k-l / k-l
!
    call add_drv2_1pd(lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,e1d(2), &
                      e2d(7),d2r2dx26t(1,2))
!
!  k-l / k-m
!
    call add_drv2_2pd(lklocal,lllocal,lklocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xkm,ykm,zkm,e2d(8))
!
!  k-l / j-l
!
    call add_drv2_2pd(lklocal,lllocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xjl,yjl,zjl,e2d(9))
!
!  k-l / j-m
!
    call add_drv2_2pd(lklocal,lllocal,ljlocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xjm,yjm,zjm,e2d(10))
!
!  k-l / l-m
!
    call add_drv2_2pd(lklocal,lllocal,lllocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xlm,ylm,zlm,e2d(11))
!
!  k-m / k-m
!
    call add_drv2_1pd(lklocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xkm,ykm,zkm,e1d(3), &
                      e2d(12),d2r2dx26t(1,3))
!
!  k-m / j-l
!
    call add_drv2_2pd(lklocal,lmlocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkm,ykm,zkm,xjl,yjl,zjl,e2d(13))
!
!  k-m / j-m
!
    call add_drv2_2pd(lklocal,lmlocal,ljlocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xkm,ykm,zkm,xjm,yjm,zjm,e2d(14))
!     
!  k-m / l-m
! 
    call add_drv2_2pd(lklocal,lmlocal,lllocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xkm,ykm,zkm,xlm,ylm,zlm,e2d(15))
!     
!  j-l / j-l
! 
    call add_drv2_1pd(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,e1d(4), &
                      e2d(16),d2r2dx26t(1,4))
!
!  j-l / j-m
!
    call add_drv2_2pd(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xjm,yjm,zjm,e2d(17))
!
!  j-l / l-m
!
    call add_drv2_2pd(ljlocal,lllocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xlm,ylm,zlm,e2d(18))
!
!  j-m / j-m
!
    call add_drv2_1pd(ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjm,yjm,zjm,e1d(5), &
                      e2d(19),d2r2dx26t(1,5))
!
!  j-m / l-m
!
    call add_drv2_2pd(ljlocal,lmlocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjm,yjm,zjm,xlm,ylm,zlm,e2d(20))
!
!  l-m / l-m
!
    call add_drv2_1pd(lllocal,lmlocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xlm,ylm,zlm,e1d(6), &
                      e2d(21),d2r2dx26t(1,6))
!-------------------------------------------------------
!  Second derivatives of ftors mixed with other terms  |
!-------------------------------------------------------
!
!  i-j / j-k
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(1,nt)
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm)
!
!  i-j / k-l
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(2,nt)
    call add_drv2_2pd(lilocal,ljlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xkl,ykl,zkl,d2trm)
!
!  i-j / k-m
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(3,nt)
    call add_drv2_2pd(lilocal,ljlocal,lklocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xkm,ykm,zkm,d2trm)
!
!  i-j / j-l
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(4,nt)
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm)
!
!  i-j / j-m
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(5,nt)
    call add_drv2_2pd(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm)
!
!  i-j / l-m
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(6,nt)
    call add_drv2_2pd(lilocal,ljlocal,lllocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xlm,ylm,zlm,d2trm)
!
!  i-k / j-k
!
    d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(1,nt)
    call add_drv2_2pd(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm)
!
!  i-k / k-l
!
    d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(2,nt)
    call add_drv2_2pd(lilocal,lklocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xkl,ykl,zkl,d2trm)
!
!  i-k / k-m
!
    d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(3,nt)
    call add_drv2_2pd(lilocal,lklocal,lklocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xkm,ykm,zkm,d2trm)
!
!  i-k / j-l
!
    d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(4,nt)
    call add_drv2_2pd(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xjl,yjl,zjl,d2trm)
!
!  i-k / j-m
!
    d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(5,nt)
    call add_drv2_2pd(lilocal,lklocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xjm,yjm,zjm,d2trm)
!
!  i-k / l-m
!
    d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(6,nt)
    call add_drv2_2pd(lilocal,lklocal,lllocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xlm,ylm,zlm,d2trm)
!
!  j-k / j-k
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(1,nt)
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,xjk,yjk,zjk,d2trm)
!
!  j-k / k-l
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(2,nt)
    call add_drv2_2pd(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm)
!
!  j-k / k-m
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(3,nt)
    call add_drv2_2pd(ljlocal,lklocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xkm,ykm,zkm,d2trm)
!
!  j-k / j-l
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(4,nt)
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm)
!
!  j-k / j-m
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(5,nt)
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,d2trm)
!
!  j-k / l-m
!
    d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(6,nt)
    call add_drv2_2pd(ljlocal,lklocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xlm,ylm,zlm,d2trm)
!--------------------------------------------------
!  Second derivatives of ftors mixed with fangle  |
!--------------------------------------------------
!
!  j-k / j-k
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(1,nt)
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,xjk,yjk,zjk,d2trm)
!
!  j-k / k-l
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(2,nt)
    call add_drv2_2pd(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm)
!
!  j-k / k-m
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(3,nt)
    call add_drv2_2pd(ljlocal,lklocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xkm,ykm,zkm,d2trm)
!
!  j-k / j-l
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(4,nt)
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm)
!
!  j-k / j-m
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(5,nt)
    call add_drv2_2pd(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,d2trm)
!
!  j-k / l-m
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(6,nt)
    call add_drv2_2pd(ljlocal,lklocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xlm,ylm,zlm,d2trm)
!
!  j-l / j-k
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(1,nt)
    call add_drv2_2pd(ljlocal,lllocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjl,yjl,zjl,xjk,yjk,zjk,d2trm)
!
!  j-l / k-l
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(2,nt)
    call add_drv2_2pd(ljlocal,lllocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xkl,ykl,zkl,d2trm)
!
!  j-l / k-m 
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(3,nt)
    call add_drv2_2pd(ljlocal,lllocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xkm,ykm,zkm,d2trm)
!
!  j-l / j-l
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(4,nt)
    call add_drv2_2pd(ljlocal,lllocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xjl,yjl,zjl,d2trm)
!
!  j-l / j-m
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(5,nt)
    call add_drv2_2pd(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xjm,yjm,zjm,d2trm)
!
!  j-l / l-m
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(6,nt)
    call add_drv2_2pd(ljlocal,lllocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xlm,ylm,zlm,d2trm)
!
!  k-l / j-k
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(1,nt)
    call add_drv2_2pd(lklocal,lllocal,ljlocal,lklocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xkl,ykl,zkl,xjk,yjk,zjk,d2trm)
!
!  k-l / k-l
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(2,nt)
    call add_drv2_2pd(lklocal,lllocal,lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xkl,ykl,zkl,d2trm)
!
!  k-l / k-m
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(3,nt)
    call add_drv2_2pd(lklocal,lllocal,lklocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xkm,ykm,zkm,d2trm)
!
!  k-l / j-l
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(4,nt)
    call add_drv2_2pd(lklocal,lllocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xjl,yjl,zjl,d2trm)
!
!  k-l / j-m
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(5,nt)
    call add_drv2_2pd(lklocal,lllocal,ljlocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xjm,yjm,zjm,d2trm)
!
!  k-l / l-m
!
    d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(6,nt)
    call add_drv2_2pd(lklocal,lllocal,lllocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                      lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xlm,ylm,zlm,d2trm)
!
    lijklmlocal = (lilocal.or.ljlocal.or.lklocal.or.lllocal.or.lmlocal)
!
!  Loop over second neighbour of B for products of first derivatives
!
    nt2 = 0
    do nm = 1,nl-1
!
!  Exclude case where neighbour is j
!
      n = nbrno(nm,l)
      if (n.eq.j) cycle
      nt2 = nt2 + 1
!
      nloc = atom2local(n)
      lnlocal = (nloc.ne.0)
!
      if (.not.lijklmlocal.and..not.lnlocal) cycle
!
      indn = 3*(n-1)
      nxf = indn + 1
      nyf = indn + 2
      nzf = indn + 3
      if (lnlocal) then
        indn = 3*(nloc-1)
        nx = indn + 1
        ny = indn + 2
        nz = indn + 3
      endif
!
!  Compute product of factors except the current two
!
      ftorsm2 = 1.0_dp
      do nt3 = 1,ntorsion
        if (nt3.ne.nt.and.nt3.ne.nt2) then
          ftorsm2 = ftorsm2*ft(nt3)
        endif
      enddo
!
!  Compute distances for n
!
      xkn = xbnbr(nm,l) + xkl
      ykn = ybnbr(nm,l) + ykl
      zkn = zbnbr(nm,l) + zkl
!
      xjn = xbnbr(nm,l) + xbnbr(1,j)
      yjn = ybnbr(nm,l) + ybnbr(1,j)
      zjn = zbnbr(nm,l) + zbnbr(1,j)
!
      xln = xbnbr(nm,l)
      yln = ybnbr(nm,l)
      zln = zbnbr(nm,l)
!
      d2r2dx26t(1,1) = xkn*xkn
      d2r2dx26t(2,1) = ykn*ykn
      d2r2dx26t(3,1) = zkn*zkn
      d2r2dx26t(4,1) = ykn*zkn
      d2r2dx26t(5,1) = xkn*zkn
      d2r2dx26t(6,1) = xkn*ykn
!
      d2r2dx26t(1,2) = xjn*xjn
      d2r2dx26t(2,2) = yjn*yjn
      d2r2dx26t(3,2) = zjn*zjn
      d2r2dx26t(4,2) = yjn*zjn
      d2r2dx26t(5,2) = xjn*zjn
      d2r2dx26t(6,2) = xjn*yjn
!
      d2r2dx26t(1,3) = xln*xln
      d2r2dx26t(2,3) = yln*yln
      d2r2dx26t(3,3) = zln*zln
      d2r2dx26t(4,3) = yln*zln
      d2r2dx26t(5,3) = xln*zln
      d2r2dx26t(6,3) = xln*yln
!
!  j-k / j-k
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(1,nt2)
      call add_drv2_2pd(ljlocal,lklocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,xjk,yjk,zjk,d2trm)
!
!  j-k / k-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(2,nt2)
      call add_drv2_2pd(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm)
!
!  j-k / k-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(3,nt2)
      call add_drv2_2pd(ljlocal,lklocal,lklocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xjk,yjk,zjk,xkn,ykn,zkn,d2trm)
!
!  j-k / j-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(4,nt2)
      call add_drv2_2pd(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm)
!
!  j-k / j-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(5,nt2)
      call add_drv2_2pd(ljlocal,lklocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjk,yjk,zjk,xjn,yjn,zjn,d2trm)
!
!  j-k / l-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(6,nt2)
      call add_drv2_2pd(ljlocal,lklocal,lllocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xjk,yjk,zjk,xln,yln,zln,d2trm)
!
!  k-l / j-k
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(1,nt2)
      call add_drv2_2pd(lklocal,lllocal,ljlocal,lklocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xkl,ykl,zkl,xjk,yjk,zjk,d2trm)
!
!  k-l / k-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(2,nt2)
      call add_drv2_2pd(lklocal,lllocal,lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xkl,ykl,zkl,d2trm)
!
!  k-l / k-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(3,nt2)
      call add_drv2_2pd(lklocal,lllocal,lklocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                        kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xkn,ykn,zkn,d2trm)
!
!  k-l / j-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(4,nt2)
      call add_drv2_2pd(lklocal,lllocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xjl,yjl,zjl,d2trm)
!
!  k-l / j-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(5,nt2)
      call add_drv2_2pd(lklocal,lllocal,ljlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                        jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xjn,yjn,zjn,d2trm)
!
!  k-l / l-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(6,nt2)
      call add_drv2_2pd(lklocal,lllocal,lllocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                        lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xln,yln,zln,d2trm)
!
!  k-m / j-k
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(1,nt2)
      call add_drv2_2pd(lklocal,lmlocal,ljlocal,lklocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xkm,ykm,zkm,xjk,yjk,zjk,d2trm)
!
!  k-m / k-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(2,nt2)
      call add_drv2_2pd(lklocal,lmlocal,lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkm,ykm,zkm,xkl,ykl,zkl,d2trm)
!
!  k-m / k-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(3,nt2)
      call add_drv2_2pd(lklocal,lmlocal,lklocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                        kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xkm,ykm,zkm,xkn,ykn,zkn,d2trm)
!
!  k-m / j-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(4,nt2)
      call add_drv2_2pd(lklocal,lmlocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkm,ykm,zkm,xjl,yjl,zjl,d2trm)
!
!  k-m / j-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(5,nt2)
      call add_drv2_2pd(lklocal,lmlocal,ljlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xkm,ykm,zkm,xjn,yjn,zjn,d2trm)
!
!  k-m / l-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(6,nt2)
      call add_drv2_2pd(lklocal,lmlocal,lllocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                        lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xkm,ykm,zkm,xln,yln,zln,d2trm)
!
!  j-l / j-k
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(1,nt2)
      call add_drv2_2pd(ljlocal,lllocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjl,yjl,zjl,xjk,yjk,zjk,d2trm)
!
!  j-l / k-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(2,nt2)
      call add_drv2_2pd(ljlocal,lllocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xkl,ykl,zkl,d2trm)
!
!  j-l / k-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(3,nt2)
      call add_drv2_2pd(ljlocal,lllocal,lklocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xjl,yjl,zjl,xkn,ykn,zkn,d2trm)
!
!  j-l / j-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(4,nt2)
      call add_drv2_2pd(ljlocal,lllocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xjl,yjl,zjl,d2trm)
!
!  j-l / j-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(5,nt2)
      call add_drv2_2pd(ljlocal,lllocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjl,yjl,zjl,xjn,yjn,zjn,d2trm)
!
!  j-l / l-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(6,nt2)
      call add_drv2_2pd(ljlocal,lllocal,lllocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xjl,yjl,zjl,xln,yln,zln,d2trm)
!
!  j-m / j-k
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(1,nt2)
      call add_drv2_2pd(ljlocal,lmlocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjm,yjm,zjm,xjk,yjk,zjk,d2trm)
!
!  j-m / k-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(2,nt2)
      call add_drv2_2pd(ljlocal,lmlocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjm,yjm,zjm,xkl,ykl,zkl,d2trm)
!
!  j-m / k-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(3,nt2)
      call add_drv2_2pd(ljlocal,lmlocal,lklocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xkn,ykn,zkn,d2trm)
!
!  j-m / j-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(4,nt2)
      call add_drv2_2pd(ljlocal,lmlocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjm,yjm,zjm,xjl,yjl,zjl,d2trm)
!
!  j-m / j-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(5,nt2)
      call add_drv2_2pd(ljlocal,lmlocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xjn,yjn,zjn,d2trm)
!
!  j-m / l-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(6,nt2)
      call add_drv2_2pd(ljlocal,lmlocal,lllocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                        lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xln,yln,zln,d2trm)
!
!  l-m / j-k
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(1,nt2)
      call add_drv2_2pd(lllocal,lmlocal,ljlocal,lklocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xlm,ylm,zlm,xjk,yjk,zjk,d2trm)
!
!  l-m / k-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(2,nt2)
      call add_drv2_2pd(lllocal,lmlocal,lklocal,lllocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xlm,ylm,zlm,xkl,ykl,zkl,d2trm)
!
!  l-m / k-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(3,nt2)
      call add_drv2_2pd(lllocal,lmlocal,lklocal,lnlocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                        kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xlm,ylm,zlm,xkn,ykn,zkn,d2trm)
!
!  l-m / j-l
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(4,nt2)
      call add_drv2_2pd(lllocal,lmlocal,ljlocal,lllocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xlm,ylm,zlm,xjl,yjl,zjl,d2trm)
!
!  l-m / j-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(5,nt2)
      call add_drv2_2pd(lllocal,lmlocal,ljlocal,lnlocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                        jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xlm,ylm,zlm,xjn,yjn,zjn,d2trm)
!
!  l-m / l-n
!
      d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(6,nt2)
      call add_drv2_2pd(lllocal,lmlocal,lllocal,lnlocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                        lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xlm,ylm,zlm,xln,yln,zln,d2trm)
    enddo
!------------------------------------------------------------
!  Second derivatives from product of torsion and dampo_nb  |
!------------------------------------------------------------
    do nj = 1,nnbr(j)
      n = nbrno(nj,j)
      nloc = atom2local(n)
      lnlocal = (nloc.ne.0)
!
      if (.not.lijklmlocal.and..not.lnlocal) cycle
!
      indn = 3*(n-1)
      nxf = indn + 1
      nyf = indn + 2
      nzf = indn + 3
      if (lnlocal) then
        indn = 3*(nloc-1)
        nx = indn + 1
        ny = indn + 2
        nz = indn + 3
      endif
!
      xjn = xbnbr(nj,j)
      yjn = ybnbr(nj,j)
      zjn = zbnbr(nj,j)
!
      xin = xij + xjn
      yin = yij + yjn
      zin = zij + zjn
!
!  Compute product of factors except the current one
!
      dampo_nb_totm1 = 1.0_dp
      do nj2 = 1,nnbr(j)
        if (nj2.ne.nj) then
          dampo_nb_totm1 = dampo_nb_totm1*dampo_nb(nj2)
        endif
      enddo
!
      d2r2dx23n(1,1) = xij*xij
      d2r2dx23n(2,1) = yij*yij
      d2r2dx23n(3,1) = zij*zij
      d2r2dx23n(4,1) = yij*zij
      d2r2dx23n(5,1) = xij*zij
      d2r2dx23n(6,1) = xij*yij
!
      d2r2dx23n(1,2) = xin*xin
      d2r2dx23n(2,2) = yin*yin
      d2r2dx23n(3,2) = zin*zin
      d2r2dx23n(4,2) = yin*zin
      d2r2dx23n(5,2) = xin*zin
      d2r2dx23n(6,2) = xin*yin
!
      d2r2dx23n(1,3) = xjn*xjn
      d2r2dx23n(2,3) = yjn*yjn
      d2r2dx23n(3,3) = zjn*zjn
      d2r2dx23n(4,3) = yjn*zjn
      d2r2dx23n(5,3) = xjn*zjn
      d2r2dx23n(6,3) = xjn*yjn
!
!  i-j / j-k
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(1,nt)
      call add_drv2_2pd(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm)
!
!  i-j / k-l
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(2,nt)
      call add_drv2_2pd(lilocal,ljlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xkl,ykl,zkl,d2trm)
!
!  i-j / k-m
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(3,nt)
      call add_drv2_2pd(lilocal,ljlocal,lklocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xkm,ykm,zkm,d2trm)
!
!  i-j / j-l
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(4,nt)
      call add_drv2_2pd(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm)
!
!  i-j / j-m
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(5,nt)
      call add_drv2_2pd(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm)
!
!  i-j / l-m
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(6,nt)
      call add_drv2_2pd(lilocal,ljlocal,lllocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xlm,ylm,zlm,d2trm)
!
!  i-n / j-k
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(1,nt)
      call add_drv2_2pd(lilocal,lnlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xin,yin,zin,xjk,yjk,zjk,d2trm)
!
!  i-n / k-l
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(2,nt)
      call add_drv2_2pd(lilocal,lnlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xin,yin,zin,xkl,ykl,zkl,d2trm)
!
!  i-n / k-m
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(3,nt)
      call add_drv2_2pd(lilocal,lnlocal,lklocal,lmlocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                        kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xin,yin,zin,xkm,ykm,zkm,d2trm)
!
!  i-n / j-l
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(4,nt)
      call add_drv2_2pd(lilocal,lnlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xin,yin,zin,xjl,yjl,zjl,d2trm)
!
!  i-n / j-m
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(5,nt)
      call add_drv2_2pd(lilocal,lnlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xin,yin,zin,xjm,yjm,zjm,d2trm)
!
!  i-n / l-m
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(6,nt)
      call add_drv2_2pd(lilocal,lnlocal,lllocal,lmlocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                        lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xin,yin,zin,xlm,ylm,zlm,d2trm)
!
!  j-n / j-k
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(1,nt)
      call add_drv2_2pd(ljlocal,lnlocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjn,yjn,zjn,xjk,yjk,zjk,d2trm)
!
!  j-n / k-l
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(2,nt)
      call add_drv2_2pd(ljlocal,lnlocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjn,yjn,zjn,xkl,ykl,zkl,d2trm)
!       
!  j-n / k-m
!   
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(3,nt)
      call add_drv2_2pd(ljlocal,lnlocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                        kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjn,yjn,zjn,xkm,ykm,zkm,d2trm)
!
!  j-n / j-l
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(4,nt)
      call add_drv2_2pd(ljlocal,lnlocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjn,yjn,zjn,xjl,yjl,zjl,d2trm)
!
!  j-n / j-m
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(5,nt)
      call add_drv2_2pd(ljlocal,lnlocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                        jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjn,yjn,zjn,xjm,yjm,zjm,d2trm)
!
!  j-n / l-m
!
      d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(6,nt)
      call add_drv2_2pd(ljlocal,lnlocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                        lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjn,yjn,zjn,xlm,ylm,zlm,d2trm)
    enddo
  enddo
!
  deallocate(d2ft,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg3pd','d2ft')
  deallocate(dft,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg3pd','dft')
  deallocate(ft,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg3pd','ft')
!
!  Exit point to ensure memory is deallocated apart from ft
!
  1000 continue
!
  deallocate(d2dampo_nbdr2,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg3pd','d2dampo_nbdr2')
  deallocate(ddampo_nbdr,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg3pd','ddampo_nbdr')
  deallocate(dampo_nb,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg3pd','dampo_nb')
!
#ifdef TRACE
  call trace_out('gfnff_hb2_eg3pd')
#endif
!
  return
  end
