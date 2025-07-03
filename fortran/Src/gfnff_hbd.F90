!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  HB calculation without lists  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gfnff_ehbd(ehb,exb,nnbr,maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,lgrad1,lgrad2)
!
!  Compute the hydrogen bond energy without the use of list based methods
!  so that periodic boundary conditions can be handled correctly.
!  Distribributed memory parallel version.
!
!   3/22 Created from gfnff_ehb
!   7/23 nregionnocfg added
!
!  Julian Gale, Curtin University, July 2023
!
  use datatypes
  use configurations, only : lsliceatom, nregions
  use control,        only : lseok
  use current
  use derivatives
  use energies,       only : eattach, esregion12, esregion2
  use energies,       only : eregion2region, siteenergy
  use gulp_gfnff
  use m_gfnff_nbr3
  use m_gfnff_pairs
  use m_strain,       only : real1strterm
  use numbers,        only : third
  use parallel
  use symmetry,       only : lstr
  implicit none
!
!  Passed variables
!
  integer(i4),                         intent(in)  :: nnbr(numat)
  integer(i4),                         intent(in)  :: maxnbr
  integer(i4),                         intent(in)  :: nbrno(maxnbr,numat)
  logical,                             intent(in)  :: lgrad1
  logical,                             intent(in)  :: lgrad2
  real(dp),                            intent(in)  :: rbnbr(maxnbr,numat)
  real(dp),                            intent(in)  :: xbnbr(maxnbr,numat)
  real(dp),                            intent(in)  :: ybnbr(maxnbr,numat)
  real(dp),                            intent(in)  :: zbnbr(maxnbr,numat)
  real(dp),                            intent(out) :: ehb
  real(dp),                            intent(out) :: exb
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
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: ks
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
  integer(i4)                                      :: nregioni           ! Region number for i
  integer(i4)                                      :: nregionj           ! Region number for j
  integer(i4)                                      :: nregionk           ! Region number for k
  logical                                          :: lattach
  logical                                          :: lbonded3
  logical                                          :: lijbonded
  logical                                          :: likbonded
  logical                                          :: ljkbonded
  logical                                          :: lilocal
  logical                                          :: ljlocal
  logical                                          :: lklocal
  logical                                          :: loneneighC
  logical                                          :: loneneighN
  logical                                          :: lreg12
  logical                                          :: lreg2trio
  logical                                          :: lslicei
  logical                                          :: lslicej
  logical                                          :: lslicek
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
  real(dp)                                         :: dr2ds3(6,3)
  real(dp)                                         :: d2r2dx23(6,3)
  real(dp)                                         :: d2r2ds23(6,6,3)
  real(dp)                                         :: d2r2dsdx3(6,3,3)
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
!  Initialise the energies
!
  ehb = 0.0_dp
  exb = 0.0_dp
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
!  Only do energy and first derivatives if i is local
!
    if (.not.lilocal.and..not.lgrad2) cycle
!
    nregioni = nregionno(nrelf2a(i))
    lslicei = lsliceatom(nsft+nrelf2a(i))
!
!  A charge scaled term
!
    qhba = gfnff_hb_ABq(i)
!
!  Set up bonded neighbour info for i
!
    call gfnff_get_n3atoms(numat,i,1_i4)
!
    if (lgrad2) then
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
      nregionj = nregionno(nrelf2a(j))
      lslicej = lsliceatom(nsft+nrelf2a(j))
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
      if (lgrad2) then
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
            nregionk = nregionno(nrelf2a(k))
            lslicek = lsliceatom(nsft+nrelf2a(k))
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
              call gfnff_hb2_eg3d(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                  maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const,radij,ehb, &
                                  lgrad1,lgrad2)
!
!  Nitro case R-N=O...H_A
!
            elseif (nat(j).eq.8.and.loneneighN) then
              const = gfnff_hb_acid(i)*qhba*gfnff_hb_base(j)*qhbb*gfnff_hb_scale_coh
              call gfnff_hb2_eg3d(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                  maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const,radij,ehb, &
                                  lgrad1,lgrad2)
!
!  N hetero aromatic
!
            elseif (nat(j).eq.7.and.nnbr(j).eq.2) then
              const = gfnff_hb_acid(i)*qhba*gfnff_hb_base(j)*qhbb*gfnff_hb_scale_gen
              call gfnff_hb2_eg2d_rnr(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                      maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const,radij,ehb, &
                                      lgrad1,lgrad2)
            else
!
!  Default
!
              const = gfnff_hb_acid(i)*qhba*gfnff_hb_base(j)*qhbb*gfnff_hb_scale_gen
              call gfnff_hb2_eg2d(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                  maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const,radij,ehb, &
                                  lgrad1,lgrad2)
            endif
          enddo
! DEBUG - check that H isn't bonded to image of i!!
          if (j.ne.i) then
            do nj = 1,nnbr(j)
              k = nbrno(nj,j)
              nh = n_gfnff_hb_Hrptr(k)
              if (nh.eq.0) cycle
              nregionk = nregionno(nrelf2a(k))
              lslicek = lsliceatom(nsft+nrelf2a(k))
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
                call gfnff_hb2_eg3d(ljlocal,lilocal,j,i,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf,nj,nnbr, &
                                    maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,-xij,-yij,-zij,rij,const,radij,ehb, &
                                    lgrad1,lgrad2)
!
!  Nitro case R-N=O...H_A
!
              elseif (nat(i).eq.8.and.loneneighN) then
                const = gfnff_hb_base(i)*qhba*gfnff_hb_acid(j)*qhbb*gfnff_hb_scale_coh
                call gfnff_hb2_eg3d(ljlocal,lilocal,j,i,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf,nj,nnbr, &
                                    maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,-xij,-yij,-zij,rij,const,radij,ehb, &
                                    lgrad1,lgrad2)
!
!  N hetero aromatic
!
              elseif (nat(i).eq.7.and.nnbr(i).eq.2) then
                const = gfnff_hb_base(i)*qhba*gfnff_hb_acid(j)*qhbb*gfnff_hb_scale_gen
                call gfnff_hb2_eg2d_rnr(ljlocal,lilocal,j,i,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf,nj,nnbr, &
                                        maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,-xij,-yij,-zij,rij,const,radij,ehb, &
                                        lgrad1,lgrad2)
              else
!
!  Default
!
                const = gfnff_hb_base(i)*qhba*gfnff_hb_acid(j)*qhbb*gfnff_hb_scale_gen
                call gfnff_hb2_eg2d(ljlocal,lilocal,j,i,jx,jy,jz,jxf,jyf,jzf,ix,iy,iz,ixf,iyf,izf,nj,nnbr, &
                                    maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,-xij,-yij,-zij,rij,const,radij,ehb, &
                                    lgrad1,lgrad2)
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
          nregionk = nregionno(nrelf2a(k))
          lslicek = lsliceatom(nsft+nrelf2a(k))
!
!  Setup region quantities
!
          lreg12    = .false.
          lreg2trio = .false.
          if (lseok.and.nregions(ncf).gt.1) then
            lreg2trio = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1)
            if (.not.lreg2trio) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1)
          endif
          lattach = .true.
          if (lslicei.and.lslicej.and.lslicek) lattach = .false.
          if (.not.lslicei.and..not.lslicej.and..not.lslicek) lattach = .false.
!
!  Find images
!
          call gfnff_gettriads(ndim,i,j,k,r2ij,xij,yij,zij,gfnff_hbthr2,hb2_paircell,hb2_triads)
!
          if (lgrad2) then
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
            if (lgrad1) then
              xjk = hb2_triads%xjkpair(np2)
              yjk = hb2_triads%yjkpair(np2)
              zjk = hb2_triads%zjkpair(np2)
            endif
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
!  Handle regions
!
            if (lilocal) then
              if (lreg2trio) then
                esregion2 = esregion2 + eijk
              elseif (lreg12) then
                esregion12 = esregion12 + eijk
              else
                ehb = ehb + eijk
              endif
!
!  Region - region energy
!
              eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*eijk
              eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + third*eijk
              eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + third*eijk
!
!  Attachment energy
!
              if (lattach) eattach = eattach + eijk
!
!  Site energies
!
              siteenergy(i) = siteenergy(i) + third*eijk
              siteenergy(j) = siteenergy(j) + third*eijk
              siteenergy(k) = siteenergy(k) + third*eijk
            endif
!
            if (lgrad1) then
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
              if (lilocal) then
                xdrv(i) = xdrv(i) - deijkdrij*xij - deijkdrik*xik
                ydrv(i) = ydrv(i) - deijkdrij*yij - deijkdrik*yik
                zdrv(i) = zdrv(i) - deijkdrij*zij - deijkdrik*zik
!
                xdrv(j) = xdrv(j) + deijkdrij*xij - deijkdrjk*xjk
                ydrv(j) = ydrv(j) + deijkdrij*yij - deijkdrjk*yjk
                zdrv(j) = zdrv(j) + deijkdrij*zij - deijkdrjk*zjk
!
                xdrv(k) = xdrv(k) + deijkdrik*xik + deijkdrjk*xjk
                ydrv(k) = ydrv(k) + deijkdrik*yik + deijkdrjk*yjk
                zdrv(k) = zdrv(k) + deijkdrik*zik + deijkdrjk*zjk
              endif
!
              if (lstr.or.lgrad2) then
                call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,1),d2r2dx23(1,1),d2r2dsdx3(1,1,1), &
                                  d2r2ds23(1,1,1),lgrad2)
                call real1strterm(ndim,xik,yik,zik,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,2),d2r2dx23(1,2),d2r2dsdx3(1,1,2), &
                                  d2r2ds23(1,1,2),lgrad2)
                call real1strterm(ndim,xjk,yjk,zjk,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,3),d2r2dx23(1,3),d2r2dsdx3(1,1,3), &
                                  d2r2ds23(1,1,3),lgrad2)
              endif
              if (lstr) then
                if (lilocal) then
                  do kl = 1,nstrains
                    ks = nstrptr(kl)
                    rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds3(ks,1) + deijkdrik*dr2ds3(ks,2) + deijkdrjk*dr2ds3(ks,3)
                  enddo
                endif
              endif
              if (lgrad2) then
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
                call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                                  deijkdrij,d2trm,dr2ds3(1,1),d2r2dx23(1,1),d2r2ds23(1,1,1),d2r2dsdx3(1,1,1),lilocal)
!
!  i-j / i-k
!
                d2trm = - const*(dabdrik*(drdampdrij*dampo + rdamp*ddampodrij) + &
                                 acid*base*(drdampdrij*ddampodrik + rdamp*d2dampo(2)))
                call add_drv2_2dm(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                                  ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm, &
                                  dr2ds3(1,1),dr2ds3(1,2),lilocal)
!
!  i-j / j-k
!
                d2trm = - const*(dabdrjk*(drdampdrij*dampo + rdamp*ddampodrij) + &
                                 acid*base*(drdampdrij*ddampodrjk + rdamp*d2dampo(3)))
                call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                                  jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm, &
                                  dr2ds3(1,1),dr2ds3(1,3),lilocal)
!
!  i-k / i-k
!
                d2trm = - const*rdamp*(d2abdrik2*dampo + 2.0_dp*dabdrik*ddampodrik + acid*base*d2dampo(4))
                call add_drv2_1dm(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik, &
                                  deijkdrik,d2trm,dr2ds3(1,2),d2r2dx23(1,2),d2r2ds23(1,1,2),d2r2dsdx3(1,1,2),lilocal)
!
!  i-k / j-k
!
                d2trm = - const*rdamp*(d2abdrikdrjk*dampo + dabdrik*ddampodrjk + dabdrjk*ddampodrik + acid*base*d2dampo(5))
                call add_drv2_2dm(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                                  jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm, &
                                  dr2ds3(1,2),dr2ds3(1,3),lilocal)
!
!  j-k / j-k
!
                d2trm = - const*rdamp*(d2abdrjk2*dampo + 2.0_dp*dabdrjk*ddampodrjk + acid*base*d2dampo(6))
                call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk, &
                                  deijkdrjk,d2trm,dr2ds3(1,3),d2r2dx23(1,3),d2r2ds23(1,1,3),d2r2dsdx3(1,1,3),lilocal)
              endif
            endif
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
!  Only do energy and first derivatives if i is local
!
    if (.not.lilocal.and..not.lgrad2) cycle
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
    nregioni = nregionno(nrelf2a(i))
    nregionj = nregionno(nrelf2a(j))
    nregionk = nregionno(nrelf2a(k))
    lslicei = lsliceatom(nsft+nrelf2a(i))
    lslicej = lsliceatom(nsft+nrelf2a(j))
    lslicek = lsliceatom(nsft+nrelf2a(k))
!
    if (lgrad2) then
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
    endif
!
!  Setup region quantities
!
    lreg12    = .false.
    lreg2trio = .false.
    if (lseok.and.nregions(ncf).gt.1) then
      lreg2trio = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1)
      if (.not.lreg2trio) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1)
    endif
    lattach = .true.
    if (lslicei.and.lslicej.and.lslicek) lattach = .false.
    if (.not.lslicei.and..not.lslicej.and..not.lslicek) lattach = .false.
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
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xik,yik,zik,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,2),d2r2dx23(1,2),d2r2dsdx3(1,1,2), &
                        d2r2ds23(1,1,2),lgrad2)
    endif
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
      if (lilocal) then
!
!  Handle regions
!
        if (lreg2trio) then
          esregion2 = esregion2 + eijk
        elseif (lreg12) then
          esregion12 = esregion12 + eijk
        else
          exb = exb + eijk
        endif
!
!  Region - region energy
!
        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*eijk
        eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + third*eijk
        eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + third*eijk
!
!  Attachment energy
!
        if (lattach) eattach = eattach + eijk
!
!  Site energies
!
        siteenergy(i) = siteenergy(i) + third*eijk
        siteenergy(j) = siteenergy(j) + third*eijk
        siteenergy(k) = siteenergy(k) + third*eijk
      endif
!
      if (lgrad1) then
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
!
        if (lilocal) then
          xdrv(i) = xdrv(i) - deijkdrij*xij - deijkdrik*xik
          ydrv(i) = ydrv(i) - deijkdrij*yij - deijkdrik*yik
          zdrv(i) = zdrv(i) - deijkdrij*zij - deijkdrik*zik
!
          xdrv(j) = xdrv(j) + deijkdrij*xij - deijkdrjk*xjk
          ydrv(j) = ydrv(j) + deijkdrij*yij - deijkdrjk*yjk
          zdrv(j) = zdrv(j) + deijkdrij*zij - deijkdrjk*zjk
!
          xdrv(k) = xdrv(k) + deijkdrik*xik + deijkdrjk*xjk
          ydrv(k) = ydrv(k) + deijkdrik*yik + deijkdrjk*yjk
          zdrv(k) = zdrv(k) + deijkdrik*zik + deijkdrjk*zjk
        endif
!
        if (lstr.or.lgrad2) then
          call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,1),d2r2dx23(1,1),d2r2dsdx3(1,1,1), &
                            d2r2ds23(1,1,1),lgrad2)
          call real1strterm(ndim,xjk,yjk,zjk,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,3),d2r2dx23(1,3),d2r2dsdx3(1,1,3), &
                            d2r2ds23(1,1,3),lgrad2)
        endif
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds3(ks,1) + deijkdrik*dr2ds3(ks,2) + deijkdrjk*dr2ds3(ks,3)
            enddo
          endif
        endif
        if (lgrad2) then
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
          call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                            d2trm,dr2ds3(1,1),d2r2dx23(1,1), &
                            d2r2ds23(1,1,1),d2r2dsdx3(1,1,1),lilocal)
!
!  i-j / i-k
!
          d2trm = - const*rdamp*d2dampo(2)
          call add_drv2_2dm(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm, &
                            dr2ds3(1,1),dr2ds3(1,2),lilocal)
!
!  i-j / j-k
!
          d2trm = - const*(drdampdrjk*ddampodrij + rdamp*d2dampo(3))
          call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm, &
                            dr2ds3(1,1),dr2ds3(1,3),lilocal)
!
!  i-k / i-k
!
          d2trm = - const*rdamp*d2dampo(4)
          call add_drv2_1dm(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,deijkdrik, &
                            d2trm,dr2ds3(1,2),d2r2dx23(1,2), &
                            d2r2ds23(1,1,2),d2r2dsdx3(1,1,2),lilocal)
!
!  i-k / j-k
!
          d2trm = - const*(drdampdrjk*ddampodrik + rdamp*d2dampo(5))
          call add_drv2_2dm(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm, &
                            dr2ds3(1,2),dr2ds3(1,3),lilocal)
!
!  j-k / j-k
!
          d2trm = - const*(dampo*d2rdampdrjk2 + 2.0_dp*drdampdrjk*ddampodrjk + rdamp*d2dampo(6))
          call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,deijkdrjk, &
                            d2trm,dr2ds3(1,3),d2r2dx23(1,3),d2r2ds23(1,1,3),d2r2dsdx3(1,1,3),lilocal)
        endif
      endif
    enddo
  enddo
!
  end subroutine gfnff_ehbd
!
  subroutine gfnff_hb2_eg2d(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                            maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const_in,radij,ehb, &
                            lgrad1,lgrad2)
!
!  Calculates the hydrogen bonding energy and derivatives for the GFNFF force fields for nhb2
!  Distributed memory parallel version.
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
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!
!  On exit :
!
!  ehb             = the value of the energy contribution
!
!   3/22 Created from gfnff_hb2_eg2
!   7/23 nregionnocfg added
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
!  Copyright Curtin University 2023
!
!  Julian Gale, CIC, Curtin University, July 2023
!
  use datatypes
  use configurations, only : nregions, lsliceatom
  use control,        only : lseok
  use current
  use derivatives
  use energies,       only : eattach, esregion12, esregion2
  use energies,       only : eregion2region, siteenergy
  use iochannels
  use gulp_gfnff
  use m_strain,       only : real1strterm
  use neighbours
  use numbers,        only : third
  use parallel
  use spatialbo
  use symmetry,       only : lstr
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
  real(dp),    intent(inout)                       :: ehb
  real(dp),    intent(in)                          :: const_in
  real(dp),    intent(in)                          :: radij
  real(dp),    intent(in)                          :: rij
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: rbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: xbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: ybnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: zbnbr(maxnbr,maxat)
  logical,     intent(in)                          :: lgrad1
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: indm
  integer(i4)                                      :: k
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: ks
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
  integer(i4)                                      :: nregioni           ! Region number for i
  integer(i4)                                      :: nregionj           ! Region number for j
  integer(i4)                                      :: nregionk           ! Region number for k
  integer(i4)                                      :: status
  logical                                          :: lattach
  logical                                          :: lklocal
  logical                                          :: lllocal
  logical                                          :: lmlocal
  logical                                          :: lreg12
  logical                                          :: lreg2trio
  logical                                          :: lslicei
  logical                                          :: lslicej
  logical                                          :: lslicek
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
  real(dp)                                         :: eijk
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
  real(dp)                                         :: dr2ds3(6,3)
  real(dp)                                         :: dr2ds3l(6,3)
  real(dp)                                         :: dr2ds3m(6,3)
  real(dp)                                         :: d2r2dx23(6,3)
  real(dp)                                         :: d2r2ds23(6,6,3)
  real(dp)                                         :: d2r2dsdx3(6,3,3)
  real(dp)                                         :: d2r2dx23l(6,3)
  real(dp)                                         :: d2r2ds23l(6,6,3)
  real(dp)                                         :: d2r2dsdx3l(6,3,3)
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
  call trace_in('gfnff_hb2_eg2d')
#endif
!
  allocate(dampo_nb(maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg2d','dampo_nb')
  if (lgrad1) then
    allocate(ddampo_nbdr(3,maxnbr),stat=status)
    if (status/=0) call outofmemory('gfnff_hb2_eg2d','ddampo_nbdr')
    if (lgrad2) then
      allocate(d2dampo_nbdr2(6,maxnbr),stat=status)
      if (status/=0) call outofmemory('gfnff_hb2_eg2d','d2dampo_nbdr2')
    endif
  endif
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
  nregioni = nregionno(nrelf2a(i))
  nregionj = nregionno(nrelf2a(j))
  nregionk = nregionno(nrelf2a(k))
  lslicei = lsliceatom(nsft+nrelf2a(i))
  lslicej = lsliceatom(nsft+nrelf2a(j))
  lslicek = lsliceatom(nsft+nrelf2a(k))
!
  if (lgrad2) then
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
  endif
!
!  Set region 2 trio flag
!
  lreg12    = .false.
  lreg2trio = .false.
  if (lseok.and.nregions(ncf).gt.1) then
    lreg2trio = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1)
    if (.not.lreg2trio) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1)
  endif
  lattach = .true.
  if (lslicei.and.lslicej.and.lslicek) lattach = .false.
  if (.not.lslicei.and..not.lslicej.and..not.lslicek) lattach = .false.
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
    if (lgrad1) then
      dexpo_nbdrij = (hbnbcutloc/radij)*ratio2_nb*rsum2/rij**3
      dexpo_nbdril = - (hbnbcutloc/radij)*ratio2_nb/(rij*ril)
      dexpo_nbdrjl = - (hbnbcutloc/radij)*ratio2_nb/(rij*rjl)
!
      ddampo_nbdr(1,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrij
      ddampo_nbdr(2,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdril
      ddampo_nbdr(3,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrjl
!
      if (lgrad2) then
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
      endif
    endif
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
!
!  Energy
!
  eijk =  - const*rdamp*qhdampo
  if (lilocal) then
!
!  Handle regions
!
    if (lreg2trio) then
      esregion2 = esregion2 + eijk
    elseif (lreg12) then
      esregion12 = esregion12 + eijk
    else
      ehb = ehb + eijk
    endif 
!
!  Region - region energy
!
    eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*eijk
    eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + third*eijk
    eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + third*eijk
!
!  Attachment energy
!
    if (lattach) eattach = eattach + eijk
!
!  Site energies
!
    siteenergy(i) = siteenergy(i) + third*eijk
    siteenergy(j) = siteenergy(j) + third*eijk
    siteenergy(k) = siteenergy(k) + third*eijk
  endif
!
  if (lgrad1) then
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
    if (lilocal) then
      xdrv(i) = xdrv(i) - deijkdrij*xij - deijkdrik*xik
      ydrv(i) = ydrv(i) - deijkdrij*yij - deijkdrik*yik
      zdrv(i) = zdrv(i) - deijkdrij*zij - deijkdrik*zik
!
      xdrv(j) = xdrv(j) + deijkdrij*xij - deijkdrjk*xjk
      ydrv(j) = ydrv(j) + deijkdrij*yij - deijkdrjk*yjk
      zdrv(j) = zdrv(j) + deijkdrij*zij - deijkdrjk*zjk
!
      xdrv(k) = xdrv(k) + deijkdrik*xik + deijkdrjk*xjk
      ydrv(k) = ydrv(k) + deijkdrik*yik + deijkdrjk*yjk
      zdrv(k) = zdrv(k) + deijkdrik*zik + deijkdrjk*zjk
    endif
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,1),d2r2dx23(1,1), &
                        d2r2dsdx3(1,1,1),d2r2ds23(1,1,1),lgrad2)
      call real1strterm(ndim,xik,yik,zik,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,2),d2r2dx23(1,2), &
                        d2r2dsdx3(1,1,2),d2r2ds23(1,1,2),lgrad2)
      call real1strterm(ndim,xjk,yjk,zjk,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,3),d2r2dx23(1,3), &
                        d2r2dsdx3(1,1,3),d2r2ds23(1,1,3),lgrad2)
    endif
    if (lstr) then
      if (lilocal) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds3(ks,1) + deijkdrik*dr2ds3(ks,2) + deijkdrjk*dr2ds3(ks,3)
        enddo
      endif
    endif
    if (lgrad2) then
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
      call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                        d2trm,dr2ds3(1,1),d2r2dx23(1,1),d2r2ds23(1,1,1),d2r2dsdx3(1,1,1),lilocal)
!
!  i-j / i-k
!
      d2trm = - const*dampo_nb_tot*(drdampdrij*ddampodrik + rdamp*d2dampo(2))
      call add_drv2_2dm(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm, &
                        dr2ds3(1,1),dr2ds3(1,2),lilocal)
!
!  i-j / j-k
!
      d2trm = - const*dampo_nb_tot*(dampo*d2rdampdrijdrjk + drdampdrij*ddampodrjk + drdampdrjk*ddampodrij + &
                rdamp*d2dampo(3))
      call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm, &
                        dr2ds3(1,1),dr2ds3(1,3),lilocal)
!
!  i-k / i-k
!
      d2trm = - const*dampo_nb_tot*rdamp*d2dampo(4)
      call add_drv2_1dm(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,deijkdrik, &
                        d2trm,dr2ds3(1,2),d2r2dx23(1,2),d2r2ds23(1,1,2),d2r2dsdx3(1,1,2),lilocal)
!
!  i-k / j-k
!
      d2trm = - const*dampo_nb_tot*(drdampdrjk*ddampodrik + rdamp*d2dampo(5))
      call add_drv2_2dm(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm, &
                        dr2ds3(1,2),dr2ds3(1,3),lilocal)
!
!  j-k / j-k
!
      d2trm = - const*dampo_nb_tot*(dampo*d2rdampdrjk2 + 2.0_dp*drdampdrjk*ddampodrjk + rdamp*d2dampo(6))
      call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,deijkdrjk, &
                        d2trm,dr2ds3(1,3),d2r2dx23(1,3),d2r2ds23(1,1,3),d2r2dsdx3(1,1,3),lilocal)
    endif
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
      if (lilocal) then
        xdrv(i) = xdrv(i) - deijkdrij*xij - deijkdril*xil
        ydrv(i) = ydrv(i) - deijkdrij*yij - deijkdril*yil
        zdrv(i) = zdrv(i) - deijkdrij*zij - deijkdril*zil
!
        xdrv(j) = xdrv(j) + deijkdrij*xij - deijkdrjl*xjl
        ydrv(j) = ydrv(j) + deijkdrij*yij - deijkdrjl*yjl
        zdrv(j) = zdrv(j) + deijkdrij*zij - deijkdrjl*zjl
!
        xdrv(l) = xdrv(l) + deijkdril*xil + deijkdrjl*xjl
        ydrv(l) = ydrv(l) + deijkdril*yil + deijkdrjl*yjl
        zdrv(l) = zdrv(l) + deijkdril*zil + deijkdrjl*zjl
      endif
!
      if (lstr.or.lgrad2) then
        call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3l(1,1),d2r2dx23l(1,1), &
                          d2r2dsdx3l(1,1,1),d2r2ds23l(1,1,1),lgrad2)
        call real1strterm(ndim,xil,yil,zil,0.0_dp,0.0_dp,0.0_dp,dr2ds3l(1,2),d2r2dx23l(1,2), &
                          d2r2dsdx3l(1,1,2),d2r2ds23l(1,1,2),lgrad2)
        call real1strterm(ndim,xjl,yjl,zjl,0.0_dp,0.0_dp,0.0_dp,dr2ds3l(1,3),d2r2dx23l(1,3), &
                          d2r2dsdx3l(1,1,3),d2r2ds23l(1,1,3),lgrad2)
      endif
      if (lstr) then
        if (lilocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds3l(ks,1) + deijkdril*dr2ds3l(ks,2) + deijkdrjl*dr2ds3l(ks,3)
          enddo
        endif
      endif
!
      if (lgrad2) then
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
        call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                          d2trm,dr2ds3l(1,1),d2r2dx23l(1,1),d2r2ds23l(1,1,1),d2r2dsdx3l(1,1,1),lilocal)
!
!  i-j / i-l
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(2,nj)
        call add_drv2_2dm(lilocal,ljlocal,lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xil,yil,zil,d2trm, &
                          dr2ds3l(1,1),dr2ds3l(1,2),lilocal)
!
!  i-j / j-l
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(3,nj)
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm, &
                          dr2ds3l(1,1),dr2ds3l(1,3),lilocal)
!
!  i-l / i-l
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(4,nj)
        call add_drv2_1dm(lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xil,yil,zil,deijkdril, &
                          d2trm,dr2ds3l(1,2),d2r2dx23l(1,2),d2r2ds23l(1,1,2),d2r2dsdx3l(1,1,2),lilocal)
!
!  i-l / j-l
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(5,nj)
        call add_drv2_2dm(lilocal,lllocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xil,yil,zil,xjl,yjl,zjl,d2trm, &
                          dr2ds3l(1,2),dr2ds3l(1,3),lilocal)
!
!  j-l / j-l
!
        d2trm = - const*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(6,nj)
        call add_drv2_1dm(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,deijkdrjl, &
                          d2trm,dr2ds3l(1,3),d2r2dx23l(1,3),d2r2ds23l(1,1,3),d2r2dsdx3l(1,1,3),lilocal)
!------------------------------------------------
!  Contribution from dampo_nb and 1 other term  |
!------------------------------------------------
!
!  i-j / i-j
!
        d2trm = - const*dampo_nb_totm1*ddampo_nbdr(1,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
        call add_drv2_2dm(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm, &
                          dr2ds3(1,1),dr2ds3l(1,1),lilocal)
!
!  i-j / i-l
!
        d2trm = - const*dampo_nb_totm1*ddampo_nbdr(2,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
        call add_drv2_2dm(lilocal,ljlocal,lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xil,yil,zil,d2trm, &
                          dr2ds3(1,1),dr2ds3l(1,2),lilocal)
!
!  i-j / j-l
!
        d2trm = - const*dampo_nb_totm1*ddampo_nbdr(3,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm, &
                          dr2ds3(1,1),dr2ds3l(1,3),lilocal)
!
!  i-k / i-j
!
        d2trm = - const*dampo_nb_totm1*ddampo_nbdr(1,nj)*rdamp*ddampodrik
        call add_drv2_2dm(lilocal,lklocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xik,yik,zik,xij,yij,zij,d2trm, &
                          dr2ds3(1,2),dr2ds3l(1,1),lilocal)
!
!  i-k / i-l
!
        d2trm = - const*dampo_nb_totm1*ddampo_nbdr(2,nj)*rdamp*ddampodrik
        call add_drv2_2dm(lilocal,lklocal,lilocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xil,yil,zil,d2trm, &
                          dr2ds3(1,2),dr2ds3l(1,2),lilocal)
!
!  i-k / j-l
!
        d2trm = - const*dampo_nb_totm1*ddampo_nbdr(3,nj)*rdamp*ddampodrik
        call add_drv2_2dm(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xjl,yjl,zjl,d2trm, &
                          dr2ds3(1,2),dr2ds3l(1,3),lilocal)
!
!  j-k / i-j
!
        d2trm = - const*dampo_nb_totm1*ddampo_nbdr(1,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
        call add_drv2_2dm(ljlocal,lklocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjk,yjk,zjk,xij,yij,zij,d2trm, &
                          dr2ds3(1,3),dr2ds3l(1,1),lilocal)
!
!  j-k / i-l
!
        d2trm = - const*dampo_nb_totm1*ddampo_nbdr(2,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
        call add_drv2_2dm(ljlocal,lklocal,lilocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xil,yil,zil,d2trm, &
                          dr2ds3(1,3),dr2ds3l(1,2),lilocal)
!
!  j-k / j-l
!
        d2trm = - const*dampo_nb_totm1*ddampo_nbdr(3,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm, &
                          dr2ds3(1,3),dr2ds3l(1,3),lilocal)
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
          if (lstr) then
            call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3m(1,1),d2r2dx23(1,1), &
                              d2r2dsdx3(1,1,1),d2r2ds23(1,1,1),.false.)
            call real1strterm(ndim,xim,yim,zim,0.0_dp,0.0_dp,0.0_dp,dr2ds3m(1,2),d2r2dx23(1,2), &
                              d2r2dsdx3(1,1,2),d2r2ds23(1,1,2),.false.)
            call real1strterm(ndim,xjm,yjm,zjm,0.0_dp,0.0_dp,0.0_dp,dr2ds3m(1,3),d2r2dx23(1,3), &
                              d2r2dsdx3(1,1,3),d2r2ds23(1,1,3),.false.)
          endif
!
!  i-j / i-j
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(1,nk)
          call add_drv2_2dm(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm, &
                            dr2ds3l(1,1),dr2ds3m(1,1),lilocal)
!
!  i-j / i-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(2,nk)
          call add_drv2_2dm(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm, &
                            dr2ds3l(1,1),dr2ds3m(1,2),lilocal)
!
!  i-j / j-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(3,nk)
          call add_drv2_2dm(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm, &
                            dr2ds3l(1,1),dr2ds3m(1,3),lilocal)
!
!  i-l / i-j
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(1,nk)
          call add_drv2_2dm(lilocal,lllocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xil,yil,zil,xij,yij,zij,d2trm, &
                            dr2ds3l(1,2),dr2ds3m(1,1),lilocal)
!
!  i-l / i-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(2,nk)
          call add_drv2_2dm(lilocal,lllocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                            ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xil,yil,zil,xim,yim,zim,d2trm, &
                            dr2ds3l(1,2),dr2ds3m(1,2),lilocal)
!
!  i-l / j-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(3,nk)
          call add_drv2_2dm(lilocal,lllocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xil,yil,zil,xjm,yjm,zjm,d2trm, &
                            dr2ds3l(1,2),dr2ds3m(1,3),lilocal)
!
!  j-l / i-j
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(1,nk)
          call add_drv2_2dm(ljlocal,lllocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjl,yjl,zjl,xij,yij,zij,d2trm, &
                            dr2ds3l(1,3),dr2ds3m(1,1),lilocal)
!
!  j-l / i-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(2,nk)
          call add_drv2_2dm(ljlocal,lllocal,lilocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xim,yim,zim,d2trm, &
                            dr2ds3l(1,3),dr2ds3m(1,2),lilocal)
!
!  j-l / j-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(3,nk)
          call add_drv2_2dm(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xjm,yjm,zjm,d2trm, &
                            dr2ds3l(1,3),dr2ds3m(1,3),lilocal)
        enddo
      endif
    enddo
  endif
!
!  Exit point to ensure memory is deallocated
!
  1000 continue
!
  if (lgrad1) then
    if (lgrad2) then
      deallocate(d2dampo_nbdr2,stat=status)
      if (status/=0) call deallocate_error('gfnff_hb2_eg2d','d2dampo_nbdr2')
    endif
    deallocate(ddampo_nbdr,stat=status)
    if (status/=0) call deallocate_error('gfnff_hb2_eg2d','ddampo_nbdr')
  endif
  deallocate(dampo_nb,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg2d','dampo_nb')
!
#ifdef TRACE
  call trace_out('gfnff_hb2_eg2d')
#endif
!
  return
  end
!
  subroutine gfnff_hb2_eg2d_rnr(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                                maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const_in,radij,ehb, &
                                lgrad1,lgrad2)
!
!  Calculates the hydrogen bonding energy and derivatives for the GFNFF force fields for nhb2
!  Lone pair version. Distributed memory parallel version.
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
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!
!  On exit :
!
!  ehb             = the value of the energy contribution
!
!   3/22 Created from gfnff_hb2_eg2_rnr
!   7/23 nregionnocfg added
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
!  Copyright Curtin University 2023
!
!  Julian Gale, CIC, Curtin University, July 2023
!
  use datatypes
  use configurations, only : nregions, lsliceatom
  use control,        only : lseok
  use current
  use derivatives
  use energies,       only : eattach, esregion12, esregion2
  use energies,       only : eregion2region, siteenergy
  use gulp_gfnff
  use iochannels
  use m_strain,       only : real1strterm, cartstrterm
  use neighbours
  use numbers,        only : third
  use parallel
  use spatialbo
  use symmetry,       only : lstr
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
  real(dp),    intent(inout)                       :: ehb
  real(dp),    intent(in)                          :: const_in
  real(dp),    intent(in)                          :: radij
  real(dp),    intent(in)                          :: rij
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: rbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: xbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: ybnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: zbnbr(maxnbr,maxat)
  logical,     intent(in)                          :: lgrad1
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: indm
  integer(i4)                                      :: indn
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: km
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
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
  integer(i4)                                      :: nregioni           ! Region number for i
  integer(i4)                                      :: nregionj           ! Region number for j
  integer(i4)                                      :: nregionk           ! Region number for k
  integer(i4)                                      :: status
  logical                                          :: lattach
  logical                                          :: lijmlocal
  logical                                          :: lklocal
  logical                                          :: lllocal
  logical                                          :: lmlocal
  logical                                          :: lnlocal
  logical                                          :: llpok
  logical                                          :: lreg12
  logical                                          :: lreg2trio
  logical                                          :: lslicei
  logical                                          :: lslicej
  logical                                          :: lslicek
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
  real(dp)                                         :: dr2ds3(6,3)
  real(dp)                                         :: d2r2dx23(6,3)
  real(dp)                                         :: d2r2ds23(6,6,3)
  real(dp)                                         :: d2r2dsdx3(6,3,3)
  real(dp)                                         :: dr2ds3m(6,3)
  real(dp)                                         :: d2r2dx23m(6,3)
  real(dp)                                         :: d2r2ds23m(6,6,3)
  real(dp)                                         :: d2r2dsdx3m(6,3,3)
  real(dp)                                         :: dr2ds3n(6,3)
  real(dp)                                         :: d2r2dx23n(6,3)
  real(dp)                                         :: d2r2ds23n(6,6,3)
  real(dp)                                         :: d2r2dsdx3n(6,3,3)
  real(dp)                                         :: dr2dsra(6)
  real(dp)                                         :: d2r2dx2ra(6)
  real(dp)                                         :: d2r2ds2ra(6,6)
  real(dp)                                         :: d2r2dsdxra(6,3)
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
  real(dp)                                         :: eijk
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
  real(dp)                                         :: trms
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
  real(dp)                                         :: d2r2ds2mix(6,6)
  real(dp)                                         :: d2r2dsdij(6,3)
  real(dp)                                         :: d2r2dsdjl(6,3)
  real(dp)                                         :: dblpds(6,3)
  real(dp)                                         :: dblpsumds(6,3)
  real(dp)                                         :: d2blpdsdx(6,3,3)
  real(dp)                                         :: d2blpds2(6,6,3)
  real(dp)                                         :: d2r2blpds2(6,6)
  real(dp)                                         :: d2r2blpdx2(6)
  real(dp)                                         :: d2r2blpdsdx(6,3)
  real(dp)                                         :: dr2blpds(6)
  real(dp)                                         :: dr2blpsumds(6)
  real(dp)                                         :: dtrmsumds(6)
  real(dp)                                         :: dxyzijds(6,3)
  real(dp)                                         :: dxyzds(6,3)
  real(dp)                                         :: dxyzsumds(6,3)
  real(dp)                                         :: d2xyzdsdx(6,3,3)
  real(dp)                                         :: d2xyzds2(6,6,3)
#ifdef TRACE
  call trace_in('gfnff_hb2_eg2d_rnr')
#endif
!
  allocate(dampo_nb(maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg2d_rnr','dampo_nb')
  if (lgrad1) then
    allocate(ddampo_nbdr(3,maxnbr),stat=status)
    if (status/=0) call outofmemory('gfnff_hb2_eg2d_rnr','ddampo_nbdr')
    if (lgrad2) then
      allocate(d2dampo_nbdr2(6,maxnbr),stat=status)
      if (status/=0) call outofmemory('gfnff_hb2_eg2d_rnr','d2dampo_nbdr2')
    endif
  endif
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
  nregioni = nregionno(nrelf2a(i))
  nregionj = nregionno(nrelf2a(j))
  nregionk = nregionno(nrelf2a(k))
  lslicei = lsliceatom(nsft+nrelf2a(i))
  lslicej = lsliceatom(nsft+nrelf2a(j))
  lslicek = lsliceatom(nsft+nrelf2a(k))
!
  if (lgrad2) then
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
  endif
!
!  Set region 2 trio flag
!
  lreg12    = .false. 
  lreg2trio = .false. 
  if (lseok.and.nregions(ncf).gt.1) then
    lreg2trio = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1)
    if (.not.lreg2trio) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1)
  endif
  lattach = .true.
  if (lslicei.and.lslicej.and.lslicek) lattach = .false. 
  if (.not.lslicei.and..not.lslicej.and..not.lslicek) lattach = .false. 
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
      if (lgrad1) then
        dexpo_nbdrij = (hbnbcutloc/radij)*ratio2_nb*rsum2/rij**3
        dexpo_nbdrim = - (hbnbcutloc/radij)*ratio2_nb/(rij*rim)
        dexpo_nbdrjm = - (hbnbcutloc/radij)*ratio2_nb/(rij*rjm)
!
        ddampo_nbdr(1,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrij
        ddampo_nbdr(2,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrim
        ddampo_nbdr(3,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrjm
!
        if (lgrad2) then
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
        endif
      endif
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
!
!  Energy
!
  eijk = - const*rdamp*qhdampo
  if (lilocal) then
!
!  Handle regions
!
    if (lreg2trio) then
      esregion2 = esregion2 + eijk
    elseif (lreg12) then
      esregion12 = esregion12 + eijk
    else
      ehb = ehb + eijk
    endif
!
!  Region - region energy
!
    eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*eijk
    eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + third*eijk
    eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + third*eijk
!
!  Attachment energy
!
    if (lattach) eattach = eattach + eijk
!
!  Site energies
!
    siteenergy(i) = siteenergy(i) + third*eijk
    siteenergy(j) = siteenergy(j) + third*eijk
    siteenergy(k) = siteenergy(k) + third*eijk
  endif
!
  if (lgrad1) then
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
    if (lilocal) then
      xdrv(i) = xdrv(i) - deijkdrij*xij - deijkdrik*xik
      ydrv(i) = ydrv(i) - deijkdrij*yij - deijkdrik*yik
      zdrv(i) = zdrv(i) - deijkdrij*zij - deijkdrik*zik
!
      xdrv(j) = xdrv(j) + deijkdrij*xij - deijkdrjk*xjk
      ydrv(j) = ydrv(j) + deijkdrij*yij - deijkdrjk*yjk
      zdrv(j) = zdrv(j) + deijkdrij*zij - deijkdrjk*zjk
!
      xdrv(k) = xdrv(k) + deijkdrik*xik + deijkdrjk*xjk
      ydrv(k) = ydrv(k) + deijkdrik*yik + deijkdrjk*yjk
      zdrv(k) = zdrv(k) + deijkdrik*zik + deijkdrjk*zjk
    endif
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,1),d2r2dx23(1,1),d2r2dsdx3(1,1,1), &
                        d2r2ds23(1,1,1),lgrad2)
      call real1strterm(ndim,xik,yik,zik,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,2),d2r2dx23(1,2),d2r2dsdx3(1,1,2), &
                        d2r2ds23(1,1,2),lgrad2)
      call real1strterm(ndim,xjk,yjk,zjk,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,3),d2r2dx23(1,3),d2r2dsdx3(1,1,3), &
                        d2r2ds23(1,1,3),lgrad2)
    endif
    if (lstr) then
      if (lilocal) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds3(ks,1) + deijkdrik*dr2ds3(ks,2) + deijkdrjk*dr2ds3(ks,3)
        enddo
      endif
    endif
!
!  Derivatives of ralp : rij component
!
    dtrmalp = - const*rdamp*dampo_nb_tot*dampo*ddampo_lpdralp
!
    if (lilocal) then
      xdrv(i) = xdrv(i) - dtrmalp*xalp
      ydrv(i) = ydrv(i) - dtrmalp*yalp
      zdrv(i) = zdrv(i) - dtrmalp*zalp
!
      xdrv(j) = xdrv(j) + dtrmalp*xalp
      ydrv(j) = ydrv(j) + dtrmalp*yalp
      zdrv(j) = zdrv(j) + dtrmalp*zalp
    endif
!
    if (lstr) then
      call cartstrterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dxyzijds,d2xyzdsdx,d2xyzds2,lgrad2)
!
      do kl = 1,nstrains
        ks = nstrptr(kl)
        dr2dsra(kl) = xalp*dxyzijds(ks,1) + yalp*dxyzijds(ks,2) + zalp*dxyzijds(ks,3)
      enddo
!
      if (lilocal) then
        do kl = 1,nstrains
          rstrd(kl) = rstrd(kl) + dtrmalp*dr2dsra(kl)
        enddo
      endif
    endif
    if (lgrad2) then
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
      if (lstr) then
        do kl = 1,nstrains
          d2r2dsdxra(kl,1) = xalp*d2xyzdsdx(kl,1,1) + yalp*d2xyzdsdx(kl,1,2) + zalp*d2xyzdsdx(kl,1,3) + dxyzijds(kl,1)
          d2r2dsdxra(kl,2) = xalp*d2xyzdsdx(kl,2,1) + yalp*d2xyzdsdx(kl,2,2) + zalp*d2xyzdsdx(kl,2,3) + dxyzijds(kl,2)
          d2r2dsdxra(kl,3) = xalp*d2xyzdsdx(kl,3,1) + yalp*d2xyzdsdx(kl,3,2) + zalp*d2xyzdsdx(kl,3,3) + dxyzijds(kl,3)
!
          do km = 1,nstrains
            d2r2ds2ra(km,kl) = xalp*d2xyzds2(km,kl,1) + dxyzijds(km,1)*dxyzijds(kl,1) + &
                               yalp*d2xyzds2(km,kl,2) + dxyzijds(km,2)*dxyzijds(kl,2) + &
                               zalp*d2xyzds2(km,kl,3) + dxyzijds(km,3)*dxyzijds(kl,3)
          enddo
        enddo
      endif
!
      call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xalp,yalp,zalp, &
                        dtrmalp,d2trmalp2,dr2dsra,d2r2dx2ra,d2r2ds2ra,d2r2dsdxra,lilocal)
!
!  Mixed terms
!
!  i-j / i-j
!
      d2trm = - const*dampo_nb_tot*(dampo*rdamp*d2dampo_lp(2) + dampo*drdampdrij*ddampo_lpdralp + &
                rdamp*ddampodrij*ddampo_lpdralp)
      call add_drv2_2dm(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xalp,yalp,zalp,d2trm, &
                        dr2ds3(1,1),dr2dsra,lilocal)
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
      call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                        deijkdrij,d2trm,dr2ds3(1,1),d2r2dx23(1,1),d2r2ds23(1,1,1),d2r2dsdx3(1,1,1),lilocal)
!
!  i-j / i-k
!
      d2trm = - const*dampo_nb_tot*(dampo_lp*(drdampdrij*ddampodrik + rdamp*d2dampo(2)) + &
                rdamp*ddampodrik*ddampo_lpdrij)
      call add_drv2_2dm(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm, &
                        dr2ds3(1,1),dr2ds3(1,2),lilocal)
!
!  i-j (ralp)/ i-k 
!
      d2trm = - const*dampo_nb_tot*rdamp*ddampodrik*ddampo_lpdralp
      call add_drv2_2dm(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xalp,yalp,zalp,xik,yik,zik,d2trm, &
                        dr2dsra,dr2ds3(1,2),lilocal)
!
!  i-j / j-k
!
      d2trm = - const*dampo_nb_tot*(dampo_lp*(dampo*d2rdampdrijdrjk + drdampdrij*ddampodrjk + drdampdrjk*ddampodrij + &
                rdamp*d2dampo(3)) + (rdamp*ddampodrjk + dampo*drdampdrjk)*ddampo_lpdrij)
      call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm, &
                        dr2ds3(1,1),dr2ds3(1,3),lilocal)
!
!  i-j (ralp)/ j-k 
!
      d2trm = - const*dampo_nb_tot*(rdamp*ddampodrjk + dampo*drdampdrjk)*ddampo_lpdralp
      call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xalp,yalp,zalp,xjk,yjk,zjk,d2trm, &
                        dr2dsra,dr2ds3(1,3),lilocal)
!
!  i-k / i-k
!
      d2trm = - const*dampo_nb_tot*dampo_lp*rdamp*d2dampo(4)
      call add_drv2_1dm(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik, &
                        deijkdrik,d2trm,dr2ds3(1,2),d2r2dx23(1,2),d2r2ds23(1,1,2),d2r2dsdx3(1,1,2),lilocal)
!
!  i-k / j-k
!
      d2trm = - const*dampo_nb_tot*dampo_lp*(drdampdrjk*ddampodrik + rdamp*d2dampo(5))
      call add_drv2_2dm(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm, &
                        dr2ds3(1,2),dr2ds3(1,3),lilocal)
!
!  j-k / j-k
!
      d2trm = - const*dampo_nb_tot*dampo_lp*(dampo*d2rdampdrjk2 + 2.0_dp*drdampdrjk*ddampodrjk + rdamp*d2dampo(6))
      call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk, &
                        deijkdrjk,d2trm,dr2ds3(1,3),d2r2dx23(1,3),d2r2ds23(1,1,3),d2r2dsdx3(1,1,3),lilocal)
    endif
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
    if (lgrad2) then
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
    endif
!
!  First pass for sum of strain derivatives for blpsum
!
    if (lstr.and.lgrad2) then
      dblpsumds(1:nstrains,1:3) = 0.0_dp
      dr2blpsumds(1:nstrains) = 0.0_dp
      dxyzsumds(1:nstrains,1:3) = 0.0_dp
      dtrmsumds(1:nstrains) = 0.0_dp
!
      do nj = 1,nnbr(j)
        l = nbrno(nj,j)
        xjl = xbnbr(nj,j)
        yjl = ybnbr(nj,j)
        zjl = zbnbr(nj,j)
!
        call cartstrterm(ndim,xjl,yjl,zjl,0.0_dp,0.0_dp,0.0_dp,dxyzds,d2xyzdsdx,d2xyzds2,lgrad2)
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          dblpds(kl,1:3) = - (rlp_dist/rblpsum)*dxyzds(kl,1:3)
          trms = xblpsum*dxyzds(kl,1) + yblpsum*dxyzds(kl,2) + zblpsum*dxyzds(kl,3)
          dtrmsumds(kl) = dtrmsumds(kl) + trms
          dblpds(kl,1) = dblpds(kl,1) + (rlp_dist/rblpsum**3)*xblpsum*trms
          dblpds(kl,2) = dblpds(kl,2) + (rlp_dist/rblpsum**3)*yblpsum*trms
          dblpds(kl,3) = dblpds(kl,3) + (rlp_dist/rblpsum**3)*zblpsum*trms
!
          dr2blpsumds(kl) = dr2blpsumds(kl) + xalp*dblpds(kl,1) + yalp*dblpds(kl,2) + zalp*dblpds(kl,3)
          dblpsumds(kl,1:3) = dblpsumds(kl,1:3) + dblpds(kl,1:3)
          dxyzsumds(kl,1:3) = dxyzsumds(kl,1:3) + dxyzds(kl,1:3)
        enddo
      enddo
!
!  Overall strain-strain terms
!
      if (lilocal) then
        do kk = 1,nstrains
          ks = nstrptr(kk)
          do kl = 1,nstrains
            kt = nstrptr(kl)
            sderv2(kl,kk) = sderv2(kl,kk) + d2trmalp2*dr2blpsumds(kt)*dr2blpsumds(ks) &
                                          + dtrmalp*(dblpsumds(kt,1)*dblpsumds(ks,1) + &
                                                     dblpsumds(kt,2)*dblpsumds(ks,2) + &
                                                     dblpsumds(kt,3)*dblpsumds(ks,3))
!
            sderv2(kl,kk) = sderv2(kl,kk) + dtrmalp*(rlp_dist/rblpsum**3)*( &
                                            xalp*dxyzsumds(kl,1)*dtrmsumds(kk) + &
                                            yalp*dxyzsumds(kl,2)*dtrmsumds(kk) + &
                                            zalp*dxyzsumds(kl,3)*dtrmsumds(kk) + &
                                            xalp*dxyzsumds(kk,1)*dtrmsumds(kl) + &
                                            yalp*dxyzsumds(kk,2)*dtrmsumds(kl) + &
                                            zalp*dxyzsumds(kk,3)*dtrmsumds(kl) + &
                                            xalp*xblpsum*(dxyzsumds(kk,1)*dxyzsumds(kl,1) + &
                                                          dxyzsumds(kk,2)*dxyzsumds(kl,2) + &
                                                          dxyzsumds(kk,3)*dxyzsumds(kl,3)) + &
                                            yalp*yblpsum*(dxyzsumds(kk,1)*dxyzsumds(kl,1) + &
                                                          dxyzsumds(kk,2)*dxyzsumds(kl,2) + &
                                                          dxyzsumds(kk,3)*dxyzsumds(kl,3)) + &
                                            zalp*zblpsum*(dxyzsumds(kk,1)*dxyzsumds(kl,1) + &
                                                          dxyzsumds(kk,2)*dxyzsumds(kl,2) + &
                                                          dxyzsumds(kk,3)*dxyzsumds(kl,3)))
!
            sderv2(kl,kk) = sderv2(kl,kk) - 3.0_dp*dtrmalp*(rlp_dist/rblpsum**5)*  &
                                            dtrmsumds(kk)*dtrmsumds(kl)*( &
                                            xalp*xblpsum + yalp*yblpsum + zalp*zblpsum)
          enddo
        enddo
      endif
    endif
!
!  Derivatives of dampo_lp with respect to ralp and rblp
!
    do nj = 1,nnbr(j)
      l = nbrno(nj,j)
      lloc = atom2local(l)
      lllocal = (lloc.ne.0)
!
      if (lgrad2) then
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
      endif
      xjl = xbnbr(nj,j)
      yjl = ybnbr(nj,j)
      zjl = zbnbr(nj,j)
!
      if (lilocal) then
        xdrv(j) = xdrv(j) - dtrmalp*dr2alp(1)
        ydrv(j) = ydrv(j) - dtrmalp*dr2alp(2)
        zdrv(j) = zdrv(j) - dtrmalp*dr2alp(3)
!
        xdrv(l) = xdrv(l) + dtrmalp*dr2alp(1)
        ydrv(l) = ydrv(l) + dtrmalp*dr2alp(2)
        zdrv(l) = zdrv(l) + dtrmalp*dr2alp(3)
      endif
!
      if (lstr) then
        call cartstrterm(ndim,xjl,yjl,zjl,0.0_dp,0.0_dp,0.0_dp,dxyzds,d2xyzdsdx,d2xyzds2,lgrad2)
!
        do kl = 1,nstrains
          ks = nstrptr(kl)
          dblpds(kl,1:3) = - (rlp_dist/rblpsum)*dxyzds(kl,1:3)
          trms = xblpsum*dxyzds(kl,1) + yblpsum*dxyzds(kl,2) + zblpsum*dxyzds(kl,3)
          dblpds(kl,1) = dblpds(kl,1) + (rlp_dist/rblpsum**3)*xblpsum*trms
          dblpds(kl,2) = dblpds(kl,2) + (rlp_dist/rblpsum**3)*yblpsum*trms
          dblpds(kl,3) = dblpds(kl,3) + (rlp_dist/rblpsum**3)*zblpsum*trms
          dr2blpds(kl) = xalp*dblpds(kl,1) + yalp*dblpds(kl,2) + zalp*dblpds(kl,3)
        enddo
        if (lilocal) then
          do kl = 1,nstrains
            rstrd(kl) = rstrd(kl) + dtrmalp*dr2blpds(kl)
          enddo
        endif
      endif
!
      if (lgrad2) then
        if (lstr) then
          do kl = 1,nstrains
            trms = dtrmsumds(kl)
            d2blpdsdx(kl,1,1) = - (rlp_dist/rblpsum)*d2xyzdsdx(kl,1,1) + (rlp_dist/rblpsum**3)*( &
                                dxyzsumds(kl,1)*xblpsum + xblpsum*(dxyzsumds(kl,1) + &
                                xblpsum*d2xyzdsdx(kl,1,1) + yblpsum*d2xyzdsdx(kl,1,2) + &
                                zblpsum*d2xyzdsdx(kl,1,3)) + trms) &
                                - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*xblpsum*trms
            d2blpdsdx(kl,2,1) = - (rlp_dist/rblpsum)*d2xyzdsdx(kl,2,1) + (rlp_dist/rblpsum**3)*( &
                                dxyzsumds(kl,1)*yblpsum + xblpsum*(dxyzsumds(kl,2) + &
                                xblpsum*d2xyzdsdx(kl,2,1) + yblpsum*d2xyzdsdx(kl,2,2) + &
                                zblpsum*d2xyzdsdx(kl,2,3))) &
                                - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*yblpsum*trms
            d2blpdsdx(kl,3,1) = - (rlp_dist/rblpsum)*d2xyzdsdx(kl,3,1) + (rlp_dist/rblpsum**3)*( &
                                dxyzsumds(kl,1)*zblpsum + xblpsum*(dxyzsumds(kl,3) + &
                                xblpsum*d2xyzdsdx(kl,3,1) + yblpsum*d2xyzdsdx(kl,3,2) + &
                                zblpsum*d2xyzdsdx(kl,3,3))) &
                                - 3.0_dp*(rlp_dist/rblpsum**5)*xblpsum*zblpsum*trms
!
            d2blpdsdx(kl,1,2) = - (rlp_dist/rblpsum)*d2xyzdsdx(kl,1,2) + (rlp_dist/rblpsum**3)*( &
                                dxyzsumds(kl,2)*xblpsum + yblpsum*(dxyzsumds(kl,1) + &
                                xblpsum*d2xyzdsdx(kl,1,1) + yblpsum*d2xyzdsdx(kl,1,2) + &
                                zblpsum*d2xyzdsdx(kl,1,3))) &
                                - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*xblpsum*trms
            d2blpdsdx(kl,2,2) = - (rlp_dist/rblpsum)*d2xyzdsdx(kl,2,2) + (rlp_dist/rblpsum**3)*( &
                                dxyzsumds(kl,2)*yblpsum + yblpsum*(dxyzsumds(kl,2) + &
                                xblpsum*d2xyzdsdx(kl,2,1) + yblpsum*d2xyzdsdx(kl,2,2) + &
                                zblpsum*d2xyzdsdx(kl,2,3)) + trms) &
                                - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*yblpsum*trms
            d2blpdsdx(kl,3,2) = - (rlp_dist/rblpsum)*d2xyzdsdx(kl,3,2) + (rlp_dist/rblpsum**3)*( &
                                dxyzsumds(kl,2)*zblpsum + yblpsum*(dxyzsumds(kl,3) + &
                                xblpsum*d2xyzdsdx(kl,3,1) + yblpsum*d2xyzdsdx(kl,3,2) + &
                                zblpsum*d2xyzdsdx(kl,3,3))) &
                                - 3.0_dp*(rlp_dist/rblpsum**5)*yblpsum*zblpsum*trms
!
            d2blpdsdx(kl,1,3) = - (rlp_dist/rblpsum)*d2xyzdsdx(kl,1,3) + (rlp_dist/rblpsum**3)*( &
                                dxyzsumds(kl,3)*xblpsum + zblpsum*(dxyzsumds(kl,1) + &
                                xblpsum*d2xyzdsdx(kl,1,1) + yblpsum*d2xyzdsdx(kl,1,2) + &
                                zblpsum*d2xyzdsdx(kl,1,3))) &
                                - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*xblpsum*trms
            d2blpdsdx(kl,2,3) = - (rlp_dist/rblpsum)*d2xyzdsdx(kl,2,3) + (rlp_dist/rblpsum**3)*( &
                                dxyzsumds(kl,3)*yblpsum + zblpsum*(dxyzsumds(kl,2) + &
                                xblpsum*d2xyzdsdx(kl,2,1) + yblpsum*d2xyzdsdx(kl,2,2) + &
                                zblpsum*d2xyzdsdx(kl,2,3))) &
                                - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*yblpsum*trms
            d2blpdsdx(kl,3,3) = - (rlp_dist/rblpsum)*d2xyzdsdx(kl,3,3) + (rlp_dist/rblpsum**3)*( &
                                dxyzsumds(kl,3)*zblpsum + zblpsum*(dxyzsumds(kl,3) + &
                                xblpsum*d2xyzdsdx(kl,3,1) + yblpsum*d2xyzdsdx(kl,3,2) + &
                                zblpsum*d2xyzdsdx(kl,3,3)) + trms) &
                                - 3.0_dp*(rlp_dist/rblpsum**5)*zblpsum*zblpsum*trms
!
            d2r2blpdsdx(kl,1) = xalp*d2blpdsdx(kl,1,1) + dxblp(1)*dblpsumds(kl,1) + &
                                yalp*d2blpdsdx(kl,1,2) + dyblp(1)*dblpsumds(kl,2) + &
                                zalp*d2blpdsdx(kl,1,3) + dzblp(1)*dblpsumds(kl,3)
            d2r2blpdsdx(kl,2) = xalp*d2blpdsdx(kl,2,1) + dxblp(2)*dblpsumds(kl,1) + &
                                yalp*d2blpdsdx(kl,2,2) + dyblp(2)*dblpsumds(kl,2) + &
                                zalp*d2blpdsdx(kl,2,3) + dzblp(2)*dblpsumds(kl,3)
            d2r2blpdsdx(kl,3) = xalp*d2blpdsdx(kl,3,1) + dxblp(3)*dblpsumds(kl,1) + &
                                yalp*d2blpdsdx(kl,3,2) + dyblp(3)*dblpsumds(kl,2) + &
                                zalp*d2blpdsdx(kl,3,3) + dzblp(3)*dblpsumds(kl,3)
!
            do km = 1,nstrains
              d2blpds2(km,kl,1) = - (rlp_dist/rblpsum)*d2xyzds2(km,kl,1) + (rlp_dist/rblpsum**3)*( &
                                    xblpsum*(xblpsum*d2xyzds2(km,kl,1) + yblpsum*d2xyzds2(km,kl,2) + &
                                             zblpsum*d2xyzds2(km,kl,3)))
              d2blpds2(km,kl,2) = - (rlp_dist/rblpsum)*d2xyzds2(km,kl,2) + (rlp_dist/rblpsum**3)*( &
                                    yblpsum*(xblpsum*d2xyzds2(km,kl,1) + yblpsum*d2xyzds2(km,kl,2) + &
                                             zblpsum*d2xyzds2(km,kl,3)))
              d2blpds2(km,kl,3) = - (rlp_dist/rblpsum)*d2xyzds2(km,kl,3) + (rlp_dist/rblpsum**3)*( &
                                    zblpsum*(xblpsum*d2xyzds2(km,kl,1) + yblpsum*d2xyzds2(km,kl,2) + &
                                             zblpsum*d2xyzds2(km,kl,3)))
              d2r2blpds2(km,kl) = xalp*d2blpds2(km,kl,1) + yalp*d2blpds2(km,kl,2) + zalp*d2blpds2(km,kl,3)
            enddo
          enddo
        else
          d2r2blpds2 = 0.0_dp
          d2r2blpdsdx = 0.0_dp
        endif
!
!  j-l/j-l for rblp: NB use of different routine to allow for explicit definition of 1/2 dr2/dadb
!
        call add_drv2_1gdm(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,dr2alp(1),dr2alp(2),dr2alp(3), &
                           dtrmalp,d2trmalp2,dr2blpsumds,d2r2blpdx2,d2r2blpds2,d2r2blpdsdx,lilocal)
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
        if (lstr) then
          do kl = 1,nstrains
            d2r2dsdij(kl,1) = dblpds(kl,1)
            d2r2dsdij(kl,2) = dblpds(kl,2)
            d2r2dsdij(kl,3) = dblpds(kl,3)
!
            d2r2dsdjl(kl,1) = dxblp(1)*dxyzijds(kl,1) + &
                              dyblp(1)*dxyzijds(kl,2) + &
                              dzblp(1)*dxyzijds(kl,3)
            d2r2dsdjl(kl,2) = dxblp(2)*dxyzijds(kl,1) + &
                              dyblp(2)*dxyzijds(kl,2) + &
                              dzblp(2)*dxyzijds(kl,3)
            d2r2dsdjl(kl,3) = dxblp(3)*dxyzijds(kl,1) + &
                              dyblp(3)*dxyzijds(kl,2) + &
                              dzblp(3)*dxyzijds(kl,3)
!
            do km = 1,nstrains
              d2r2ds2mix(km,kl) = dxyzijds(km,1)*dblpds(kl,1) + dxyzijds(kl,1)*dblpds(km,1) + &
                                  dxyzijds(km,2)*dblpds(kl,2) + dxyzijds(kl,2)*dblpds(km,2) + &
                                  dxyzijds(km,3)*dblpds(kl,3) + dxyzijds(kl,3)*dblpds(km,3)
            enddo
          enddo
        else
          d2r2dsdij = 0.0d0
          d2r2dsdjl = 0.0d0
          d2r2ds2mix = 0.0d0
        endif
!
!  i-j/j-l for rblp mixed terms: NB use of different routine due to special case
!
        call add_drv2_2sgdm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,dtrmalp,dr2alpij,dr2alp,dr2dsra, &
                            dr2blpds,d2trmalp2,d2r2dx2mix,d2r2ds2mix,d2r2dsdij,d2r2dsdjl,lilocal)
!
!  Mixed terms
!
!  i-j / j-l (ralp)
!
        d2trm = - const*dampo_nb_tot*(dampo*rdamp*d2dampo_lp(2) + dampo*drdampdrij*ddampo_lpdralp + &
                  rdamp*ddampodrij*ddampo_lpdralp)
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,dr2alp(1),dr2alp(2),dr2alp(3), &
                          d2trm,dr2ds3(1,1),dr2blpds,lilocal)
!
!  j-l (ralp)/ i-k
!
        d2trm = - const*dampo_nb_tot*rdamp*ddampodrik*ddampo_lpdralp
        call add_drv2_2dm(ljlocal,lllocal,lilocal,lklocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,dr2alp(1),dr2alp(2),dr2alp(3),xik,yik,zik,d2trm, &
                          dr2blpds,dr2ds3(1,2),lilocal)
!
!  j-l (ralp)/ j-k
!
        d2trm = - const*dampo_nb_tot*(rdamp*ddampodrjk + dampo*drdampdrjk)*ddampo_lpdralp
        call add_drv2_2dm(ljlocal,lllocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,dr2alp(1),dr2alp(2),dr2alp(3),xjk,yjk,zjk,d2trm, &
                          dr2blpds,dr2ds3(1,3),lilocal)
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
          call add_drv2_2sgnosdm(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                                 jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,dtrmalp,dr2alp,dr2alp, &
                                 d2trmalp2,d2r2dx2mix)
        enddo
      endif
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
        if (lilocal) then
          xdrv(i) = xdrv(i) - deijkdrij*xij - deijkdrim*xim
          ydrv(i) = ydrv(i) - deijkdrij*yij - deijkdrim*yim
          zdrv(i) = zdrv(i) - deijkdrij*zij - deijkdrim*zim
!
          xdrv(j) = xdrv(j) + deijkdrij*xij - deijkdrjm*xjm
          ydrv(j) = ydrv(j) + deijkdrij*yij - deijkdrjm*yjm
          zdrv(j) = zdrv(j) + deijkdrij*zij - deijkdrjm*zjm
!
          xdrv(m) = xdrv(m) + deijkdrim*xim + deijkdrjm*xjm
          ydrv(m) = ydrv(m) + deijkdrim*yim + deijkdrjm*yjm
          zdrv(m) = zdrv(m) + deijkdrim*zim + deijkdrjm*zjm
        endif
!
        if (lstr.or.lgrad2) then
          call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3m(1,1),d2r2dx23m(1,1),d2r2dsdx3m(1,1,1), &
                            d2r2ds23m(1,1,1),lgrad2)
          call real1strterm(ndim,xim,yim,zim,0.0_dp,0.0_dp,0.0_dp,dr2ds3m(1,2),d2r2dx23m(1,2),d2r2dsdx3m(1,1,2), &
                            d2r2ds23m(1,1,2),lgrad2)
          call real1strterm(ndim,xjm,yjm,zjm,0.0_dp,0.0_dp,0.0_dp,dr2ds3m(1,3),d2r2dx23m(1,3),d2r2dsdx3m(1,1,3), &
                            d2r2ds23m(1,1,3),lgrad2)
        endif
        if (lstr) then
          if (lilocal) then
            do kl = 1,nstrains
              ks = nstrptr(kl)
              rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds3m(ks,1) + deijkdrim*dr2ds3m(ks,2) + deijkdrjm*dr2ds3m(ks,3)
            enddo
          endif
        endif
        if (lgrad2) then
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
          call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij, &
                            deijkdrij,d2trm,dr2ds3m(1,1),d2r2dx23m(1,1), &
                            d2r2ds23m(1,1,1),d2r2dsdx3m(1,1,1),lilocal)
!
!  i-j / i-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(2,nj)
          call add_drv2_2dm(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm, &
                            dr2ds3m(1,1),dr2ds3m(1,2),lilocal)
!
!  i-j / j-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(3,nj)
          call add_drv2_2dm(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm, &
                            dr2ds3m(1,1),dr2ds3m(1,3),lilocal)
!
!  i-m / i-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(4,nj)
          call add_drv2_1dm(lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xim,yim,zim, &
                            deijkdrim,d2trm,dr2ds3m(1,2),d2r2dx23m(1,2), &
                            d2r2ds23m(1,1,2),d2r2dsdx3m(1,1,2),lilocal)
!
!  i-m / j-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(5,nj)
          call add_drv2_2dm(lilocal,lmlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xim,yim,zim,xjm,yjm,zjm,d2trm, &
                            dr2ds3m(1,2),dr2ds3m(1,3),lilocal)
!
!  j-m / j-m
!
          d2trm = - const*rdamp*dampo*dampo_nb_totm1*dampo_lp*d2dampo_nbdr2(6,nj)
          call add_drv2_1dm(ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjm,yjm,zjm, &
                            deijkdrjm,d2trm,dr2ds3m(1,3),d2r2dx23m(1,3), &
                            d2r2ds23m(1,1,3),d2r2dsdx3m(1,1,3),lilocal)
!------------------------------------------------
!  Contribution from dampo_nb and 1 other term  |
!------------------------------------------------
!
!  i-j / i-j
!
          d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(1,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
          call add_drv2_2dm(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm, &
                            dr2ds3(1,1),dr2ds3m(1,1),lilocal)
!
!  i-j / i-m
!
          d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(2,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
          call add_drv2_2dm(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm, &
                            dr2ds3(1,1),dr2ds3m(1,2),lilocal)
!
!  i-j / j-m
!
          d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(3,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
          call add_drv2_2dm(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm, &
                            dr2ds3(1,1),dr2ds3m(1,3),lilocal)
!
!  i-k / i-j
!
          d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(1,nj)*rdamp*ddampodrik
          call add_drv2_2dm(lilocal,lklocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xik,yik,zik,xij,yij,zij,d2trm, &
                            dr2ds3(1,2),dr2ds3m(1,1),lilocal)
!
!  i-k / i-m
!
          d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(2,nj)*rdamp*ddampodrik
          call add_drv2_2dm(lilocal,lklocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                            ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xim,yim,zim,d2trm, &
                            dr2ds3(1,2),dr2ds3m(1,2),lilocal)
!
!  i-k / j-m
!
          d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(3,nj)*rdamp*ddampodrik
          call add_drv2_2dm(lilocal,lklocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xjm,yjm,zjm,d2trm, &
                            dr2ds3(1,2),dr2ds3m(1,3),lilocal)
!
!  j-k / i-j
!
          d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(1,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
          call add_drv2_2dm(ljlocal,lklocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjk,yjk,zjk,xij,yij,zij,d2trm, &
                            dr2ds3(1,3),dr2ds3m(1,1),lilocal)
!
!  j-k / i-m
!
          d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(2,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
          call add_drv2_2dm(ljlocal,lklocal,lilocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xim,yim,zim,d2trm, &
                            dr2ds3(1,3),dr2ds3m(1,2),lilocal)
!
!  j-k / j-m
!
          d2trm = - const*dampo_nb_totm1*dampo_lp*ddampo_nbdr(3,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
          call add_drv2_2dm(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,d2trm, &
                            dr2ds3(1,3),dr2ds3m(1,3),lilocal)
!--------------------------------------------------------
!  Mixed derivatives between dampo_lp and dampo_nb_tot  |
!--------------------------------------------------------
!
!  i-j / i-j
!
          d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(1,nj)*ddampo_lpdrij
          call add_drv2_2dm(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm, &
                            dr2ds3m(1,1),dr2ds3(1,1),lilocal)
!
!  i-m / i-j
!
          d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(2,nj)*ddampo_lpdrij
          call add_drv2_2dm(lilocal,lmlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xim,yim,zim,xij,yij,zij,d2trm, &
                            dr2ds3m(1,2),dr2ds3(1,1),lilocal)
!
!  j-m / i-j
!
          d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(3,nj)*ddampo_lpdrij
          call add_drv2_2dm(ljlocal,lmlocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjm,yjm,zjm,xij,yij,zij,d2trm, &
                            dr2ds3m(1,3),dr2ds3(1,1),lilocal)
!
!  i-j / i-j (ralp)
!
          d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(1,nj)*ddampo_lpdralp
          call add_drv2_2dm(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xalp,yalp,zalp,d2trm, &
                            dr2ds3m(1,1),dr2dsra,lilocal)
!
!  i-m / i-j (ralp)
!
          d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(2,nj)*ddampo_lpdralp
          call add_drv2_2dm(lilocal,lmlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xim,yim,zim,xalp,yalp,zalp,d2trm, &
                            dr2ds3m(1,2),dr2dsra,lilocal)
!
!  j-m / i-j (ralp)
!
          d2trm = - const*dampo_nb_totm1*dampo*rdamp*ddampo_nbdr(3,nj)*ddampo_lpdralp
          call add_drv2_2dm(ljlocal,lmlocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjm,yjm,zjm,xalp,yalp,zalp,d2trm, &
                            dr2ds3m(1,3),dr2dsra,lilocal)
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
            if (lstr) then
              call cartstrterm(ndim,xjl,yjl,zjl,0.0_dp,0.0_dp,0.0_dp,dxyzds,d2xyzdsdx,d2xyzds2,.false.)
!
              do kl = 1,nstrains
                ks = nstrptr(kl)
                dblpds(kl,1:3) = - (rlp_dist/rblpsum)*dxyzds(kl,1:3)
                trms = xblpsum*dxyzds(kl,1) + yblpsum*dxyzds(kl,2) + zblpsum*dxyzds(kl,3)
                dblpds(kl,1) = dblpds(kl,1) + (rlp_dist/rblpsum**3)*xblpsum*trms
                dblpds(kl,2) = dblpds(kl,2) + (rlp_dist/rblpsum**3)*yblpsum*trms
                dblpds(kl,3) = dblpds(kl,3) + (rlp_dist/rblpsum**3)*zblpsum*trms
!
                dr2blpds(kl) = xalp*dblpds(kl,1) + yalp*dblpds(kl,2) + zalp*dblpds(kl,3)
              enddo
            endif
!
!  j-l (ralp)/ i-j
!
            d2trm = - const*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(1,nj)*ddampo_lpdralp
            call add_drv2_2dm(ljlocal,lllocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                              ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,dr2alp(1),dr2alp(2),dr2alp(3), &
                              xij,yij,zij,d2trm,dr2blpds,dr2ds3m(1,1),lilocal)
!
!  j-l (ralp)/ i-m
!
            d2trm = - const*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(2,nj)*ddampo_lpdralp
            call add_drv2_2dm(ljlocal,lllocal,lilocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                              ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,dr2alp(1),dr2alp(2),dr2alp(3), &
                              xim,yim,zim,d2trm,dr2blpds,dr2ds3m(1,2),lilocal)
!
!  j-l (ralp)/ j-m
!
            d2trm = - const*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(3,nj)*ddampo_lpdralp
            call add_drv2_2dm(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                              jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,dr2alp(1),dr2alp(2),dr2alp(3), &
                              xjm,yjm,zjm,d2trm,dr2blpds,dr2ds3m(1,3),lilocal)
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
            if (lstr) then
              call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3n(1,1),d2r2dx23n(1,1), &
                                d2r2dsdx3n(1,1,1),d2r2ds23n(1,1,1),.false.)
              call real1strterm(ndim,xin,yin,zin,0.0_dp,0.0_dp,0.0_dp,dr2ds3n(1,2),d2r2dx23n(1,2), &
                                d2r2dsdx3n(1,1,2),d2r2ds23n(1,1,2),.false.)
              call real1strterm(ndim,xjn,yjn,zjn,0.0_dp,0.0_dp,0.0_dp,dr2ds3n(1,3),d2r2dx23n(1,3), &
                                d2r2dsdx3n(1,1,3),d2r2ds23n(1,1,3),.false.)
            endif
!
!  i-j / i-j
!
            d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(1,nj)*ddampo_nbdr(1,nk)
            call add_drv2_2dm(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                              ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm, &
                              dr2ds3m(1,1),dr2ds3n(1,1),lilocal)
!
!  i-j / i-n
!
            d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(1,nj)*ddampo_nbdr(2,nk)
            call add_drv2_2dm(lilocal,ljlocal,lilocal,lnlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                              ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xij,yij,zij,xin,yin,zin,d2trm, &
                              dr2ds3m(1,1),dr2ds3n(1,2),lilocal)
!
!  i-j / j-n
!
            d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(1,nj)*ddampo_nbdr(3,nk)
            call add_drv2_2dm(lilocal,ljlocal,ljlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                              jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xij,yij,zij,xjn,yjn,zjn,d2trm, &
                              dr2ds3m(1,1),dr2ds3n(1,3),lilocal)
!
!  i-m / i-j
!
            d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(2,nj)*ddampo_nbdr(1,nk)
            call add_drv2_2dm(lilocal,lmlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                              ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xim,yim,zim,xij,yij,zij,d2trm, &
                              dr2ds3m(1,2),dr2ds3n(1,1),lilocal)
!
!  i-m / i-n
!
            d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(2,nj)*ddampo_nbdr(2,nk)
            call add_drv2_2dm(lilocal,lmlocal,lilocal,lnlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                              ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xim,yim,zim,xin,yin,zin,d2trm, &
                              dr2ds3m(1,2),dr2ds3n(1,2),lilocal)
!
!  i-m / j-n
!
            d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(2,nj)*ddampo_nbdr(3,nk)
            call add_drv2_2dm(lilocal,lmlocal,ljlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                              jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xim,yim,zim,xjn,yjn,zjn,d2trm, &
                              dr2ds3m(1,2),dr2ds3n(1,3),lilocal)
!
!  j-m / i-j
!
            d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(3,nj)*ddampo_nbdr(1,nk)
            call add_drv2_2dm(ljlocal,lmlocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                              ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjm,yjm,zjm,xij,yij,zij,d2trm, &
                              dr2ds3m(1,3),dr2ds3n(1,1),lilocal)
!
!  j-m / i-n
!
            d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(3,nj)*ddampo_nbdr(2,nk)
            call add_drv2_2dm(ljlocal,lmlocal,lilocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                              ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xin,yin,zin,d2trm, &
                              dr2ds3m(1,3),dr2ds3n(1,2),lilocal)
!
!  j-m / j-n
!
            d2trm = - const*rdamp*dampo*dampo_nb_totm2*dampo_lp*ddampo_nbdr(3,nj)*ddampo_nbdr(3,nk)
            call add_drv2_2dm(ljlocal,lmlocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                              jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xjn,yjn,zjn,d2trm, &
                              dr2ds3m(1,3),dr2ds3n(1,3),lilocal)
          enddo
        endif
      enddo
    endif
  endif
!
!  Exit point to ensure memory is deallocated
!
  1000 continue
!
  if (lgrad1) then
    if (lgrad2) then
      deallocate(d2dampo_nbdr2,stat=status)
      if (status/=0) call deallocate_error('gfnff_hb2_eg2d_rnr','d2dampo_nbdr2')
    endif
    deallocate(ddampo_nbdr,stat=status)
    if (status/=0) call deallocate_error('gfnff_hb2_eg2d_rnr','ddampo_nbdr')
  endif
  deallocate(dampo_nb,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg2d_rnr','dampo_nb')
!
#ifdef TRACE
  call trace_out('gfnff_hb2_eg2d_rnr')
#endif
!
  return
  end
!
  subroutine gfnff_hb2_eg3d(lilocal,ljlocal,i,j,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,ni,nnbr, &
                            maxnbr,nbrno,rbnbr,xbnbr,ybnbr,zbnbr,xij,yij,zij,rij,const_in,radij,ehb, &
                            lgrad1,lgrad2)
!
!  Calculates the hydrogen bonding energy and derivatives for the GFNFF force fields for nhb2
!  Special case version for carbonyls/nitryls with two in plane lone pairs.
!  Distributed memory parallel version.
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
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!
!  On exit :
!
!  ehb             = the value of the energy contribution
!
!   3/22 Created from gfnff_hb2_eg3
!   7/23 nregionnocfg added
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
!  Copyright Curtin University 2023
!
!  Julian Gale, CIC, Curtin University, July 2023
!
  use datatypes
  use configurations, only : nregions, lsliceatom
  use control,        only : lseok
  use current
  use derivatives
  use energies,       only : eattach, esregion12, esregion2
  use energies,       only : eregion2region, siteenergy
  use g_constants,    only : pi
  use gulp_gfnff
  use iochannels
  use m_strain,       only : real1strterm, realstrterms
  use neighbours
  use numbers,        only : third
  use parallel
  use spatialbo
  use symmetry,       only : lstr
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
  real(dp),    intent(inout)                       :: ehb
  real(dp),    intent(in)                          :: const_in
  real(dp),    intent(in)                          :: radij
  real(dp),    intent(in)                          :: rij
  real(dp),    intent(in)                          :: xij
  real(dp),    intent(in)                          :: yij
  real(dp),    intent(in)                          :: zij
  real(dp),    intent(in)                          :: rbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: xbnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: ybnbr(maxnbr,maxat)
  real(dp),    intent(in)                          :: zbnbr(maxnbr,maxat)
  logical,     intent(in)                          :: lgrad1
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: indk
  integer(i4)                                      :: indl
  integer(i4)                                      :: indm
  integer(i4)                                      :: indn
  integer(i4)                                      :: k
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: ks
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
  integer(i4)                                      :: nregioni           ! Region number for i
  integer(i4)                                      :: nregionj           ! Region number for j
  integer(i4)                                      :: nregionk           ! Region number for k
  integer(i4)                                      :: nt
  integer(i4)                                      :: nt2
  integer(i4)                                      :: nt3
  integer(i4)                                      :: ntorsion
  integer(i4)                                      :: status
  logical                                          :: lattach
  logical                                          :: lijklmlocal
  logical                                          :: lijmlocal
  logical                                          :: lklocal
  logical                                          :: lllocal
  logical                                          :: lmlocal
  logical                                          :: lnlocal
  logical                                          :: lreg12
  logical                                          :: lreg2trio
  logical                                          :: lslicei
  logical                                          :: lslicej
  logical                                          :: lslicek
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
  real(dp)                                         :: dr2ds3(6,6)
  real(dp)                                         :: d2r2dx23(6,3)
  real(dp)                                         :: d2r2ds23(6,6,3)
  real(dp)                                         :: d2r2dsdx3(6,3,3)
  real(dp)                                         :: dr2ds3a(6,3)
  real(dp)                                         :: d2r2dx23a(6,3)
  real(dp)                                         :: d2r2ds23a(6,6,3)
  real(dp)                                         :: d2r2dsdx3a(6,3,3)
  real(dp)                                         :: dr2ds6t(6,6)
  real(dp)                                         :: d2r2dx26t(6,6)
  real(dp)                                         :: d2r2ds26t(6,6,6)
  real(dp)                                         :: d2r2dsdx6t(6,3,6)
  real(dp)                                         :: dr2ds6t2(6,3)
  real(dp)                                         :: dr2ds3m(6,3)
  real(dp)                                         :: d2r2dx23m(6,3)
  real(dp)                                         :: d2r2ds23m(6,6,3)
  real(dp)                                         :: d2r2dsdx3m(6,3,3)
  real(dp)                                         :: dr2ds3n(6,3)
  real(dp)                                         :: d2r2dx23n(6,3)
  real(dp)                                         :: d2r2ds23n(6,6,3)
  real(dp)                                         :: d2r2dsdx3n(6,3,3)
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
  call trace_in('gfnff_hb2_eg3d')
#endif
!
  allocate(dampo_nb(maxnbr),stat=status)
  if (status/=0) call outofmemory('gfnff_hb2_eg3d','dampo_nb')
  if (lgrad1) then
    allocate(ddampo_nbdr(3,maxnbr),stat=status)
    if (status/=0) call outofmemory('gfnff_hb2_eg3d','ddampo_nbdr')
    if (lgrad2) then
      allocate(d2dampo_nbdr2(6,maxnbr),stat=status)
      if (status/=0) call outofmemory('gfnff_hb2_eg3d','d2dampo_nbdr2')
    endif
  endif
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
  nregioni = nregionno(nrelf2a(i))
  nregionj = nregionno(nrelf2a(j))
  nregionk = nregionno(nrelf2a(k))
  lslicei = lsliceatom(nsft+nrelf2a(i))
  lslicej = lsliceatom(nsft+nrelf2a(j))
  lslicek = lsliceatom(nsft+nrelf2a(k))
!
  if (lgrad2) then
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
  endif
!
!  Set region 2 trio flag
!
  lreg12    = .false. 
  lreg2trio = .false. 
  if (lseok.and.nregions(ncf).gt.1) then
    lreg2trio = (nregioni.gt.1.and.nregionj.gt.1.and.nregionk.gt.1)
    if (.not.lreg2trio) lreg12 = (nregioni.gt.1.or.nregionj.gt.1.or.nregionk.gt.1)
  endif
  lattach = .true.
  if (lslicei.and.lslicej.and.lslicek) lattach = .false. 
  if (.not.lslicei.and..not.lslicej.and..not.lslicek) lattach = .false. 
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
    if (lgrad1) then
      dexpo_nbdrij = (hbnbcutloc/radij)*ratio2_nb*rsum2/rij**3
      dexpo_nbdrim = - (hbnbcutloc/radij)*ratio2_nb/(rij*rim)
      dexpo_nbdrjm = - (hbnbcutloc/radij)*ratio2_nb/(rij*rjm)
!
      ddampo_nbdr(1,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrij
      ddampo_nbdr(2,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrim
      ddampo_nbdr(3,nj) = - 0.5_dp*dampo_nbloc*dampo_nbloc*dexpo_nbdrjm
!
      if (lgrad2) then
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
      endif
    endif
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
  if (lgrad2) then
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
  if (status/=0) call outofmemory('gfnff_hb2_eg3d','ft')
  if (lgrad1) then
    allocate(dft(6,ntorsion),stat=status)
    if (status/=0) call outofmemory('gfnff_hb2_eg3d','dft')
    if (lgrad2) then
      allocate(d2ft(21,ntorsion),stat=status)
      if (status/=0) call outofmemory('gfnff_hb2_eg3d','d2ft')
    endif
  endif
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
    call gfnff_torsion(1_i4,rkfor,rn,phi0,rjk,rkl,rkm,rjl,rjm,rlm,ftloc,e1d,e2d,lgrad1,lgrad2)
    ftloc = ftloc + gfnff_tors_hb
!
    ft(ntorsion) = ftloc
    ftors = ftors*ftloc
    if (lgrad1) then
      dft(1:6,ntorsion) = e1d(1:6)
      if (lgrad2) then
        d2ft(1:21,ntorsion) = e2d(1:21)
      endif
    endif
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
  if (lilocal) then
!
!  Handle regions
!
    if (lreg2trio) then
      esregion2 = esregion2 + eijk
    elseif (lreg12) then
      esregion12 = esregion12 + eijk
    else
      ehb = ehb + eijk
    endif
!
!  Region - region energy
!
    eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + third*eijk
    eregion2region(nregionk,nregioni) = eregion2region(nregionk,nregioni) + third*eijk
    eregion2region(nregionk,nregionj) = eregion2region(nregionk,nregionj) + third*eijk
!
!  Attachment energy
!
    if (lattach) eattach = eattach + eijk
!
!  Site energies
!
    siteenergy(i) = siteenergy(i) + third*eijk
    siteenergy(j) = siteenergy(j) + third*eijk
    siteenergy(k) = siteenergy(k) + third*eijk
  endif
!
  if (lgrad1) then
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
    if (lilocal) then
      xdrv(i) = xdrv(i) - deijkdrij*xij - deijkdrik*xik
      ydrv(i) = ydrv(i) - deijkdrij*yij - deijkdrik*yik
      zdrv(i) = zdrv(i) - deijkdrij*zij - deijkdrik*zik
!
      xdrv(j) = xdrv(j) + deijkdrij*xij - deijkdrjk*xjk
      ydrv(j) = ydrv(j) + deijkdrij*yij - deijkdrjk*yjk
      zdrv(j) = zdrv(j) + deijkdrij*zij - deijkdrjk*zjk
!
      xdrv(k) = xdrv(k) + deijkdrik*xik + deijkdrjk*xjk
      ydrv(k) = ydrv(k) + deijkdrik*yik + deijkdrjk*yjk
      zdrv(k) = zdrv(k) + deijkdrik*zik + deijkdrjk*zjk
    endif
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,1),d2r2dx23(1,1),d2r2dsdx3(1,1,1), &
                        d2r2ds23(1,1,1),lgrad2)
      call real1strterm(ndim,xik,yik,zik,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,2),d2r2dx23(1,2),d2r2dsdx3(1,1,2), &
                        d2r2ds23(1,1,2),lgrad2)
      call real1strterm(ndim,xjk,yjk,zjk,0.0_dp,0.0_dp,0.0_dp,dr2ds3(1,3),d2r2dx23(1,3),d2r2dsdx3(1,1,3), &
                        d2r2ds23(1,1,3),lgrad2)
    endif
    if (lstr) then
      if (lilocal) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds3(ks,1) + deijkdrik*dr2ds3(ks,2) + deijkdrjk*dr2ds3(ks,3)
        enddo
      endif
    endif
!
    if (lgrad2) then
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
      call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                        d2trm,dr2ds3(1,1),d2r2dx23(1,1),d2r2ds23(1,1,1),d2r2dsdx3(1,1,1),lilocal)
!
!  i-j / i-k
!
      d2trm = - fconst*dampo_nb_tot*(drdampdrij*ddampodrik + rdamp*d2dampo(2))
      call add_drv2_2dm(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xik,yik,zik,d2trm, &
                        dr2ds3(1,1),dr2ds3(1,2),lilocal)
!
!  i-j / j-k
!
      d2trm = - fconst*dampo_nb_tot*(dampo*d2rdampdrijdrjk + drdampdrij*ddampodrjk + drdampdrjk*ddampodrij + &
                rdamp*d2dampo(3))
      call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm, &
                        dr2ds3(1,1),dr2ds3(1,3),lilocal)
!
!  i-k / i-k
!
      d2trm = - fconst*dampo_nb_tot*rdamp*d2dampo(4)
      call add_drv2_1dm(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,deijkdrik, &
                        d2trm,dr2ds3(1,2),d2r2dx23(1,2),d2r2ds23(1,1,2),d2r2dsdx3(1,1,2),lilocal)
!
!  i-k / j-k
!
      d2trm = - fconst*dampo_nb_tot*(drdampdrjk*ddampodrik + rdamp*d2dampo(5))
      call add_drv2_2dm(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm, &
                        dr2ds3(1,2),dr2ds3(1,3),lilocal)
!
!  j-k / j-k
!
      d2trm = - fconst*dampo_nb_tot*(dampo*d2rdampdrjk2 + 2.0_dp*drdampdrjk*ddampodrjk + rdamp*d2dampo(6))
      call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,deijkdrjk, &
                        d2trm,dr2ds3(1,3),d2r2dx23(1,3),d2r2ds23(1,1,3),d2r2dsdx3(1,1,3),lilocal)
    endif
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
    if (lilocal) then
      xdrv(j) = xdrv(j) - dfangledrjk*xjk
      ydrv(j) = ydrv(j) - dfangledrjk*yjk
      zdrv(j) = zdrv(j) - dfangledrjk*zjk
      xdrv(k) = xdrv(k) + dfangledrjk*xjk
      ydrv(k) = ydrv(k) + dfangledrjk*yjk
      zdrv(k) = zdrv(k) + dfangledrjk*zjk
    endif
!
!  j-l
!
    dfangledrjl = - const*ftors*rdamp*qhdampo*dfangledcosa*dcosdrjl      ! (1/r(dEdr))
    if (lilocal) then
      xdrv(j) = xdrv(j) - dfangledrjl*xjl
      ydrv(j) = ydrv(j) - dfangledrjl*yjl
      zdrv(j) = zdrv(j) - dfangledrjl*zjl
      xdrv(l) = xdrv(l) + dfangledrjl*xjl
      ydrv(l) = ydrv(l) + dfangledrjl*yjl
      zdrv(l) = zdrv(l) + dfangledrjl*zjl
    endif
!
!  k-l
!
    dfangledrkl = - const*ftors*rdamp*qhdampo*dfangledcosa*dcosdrkl      ! (1/r(dEdr))
    if (lilocal) then
      xdrv(k) = xdrv(k) - dfangledrkl*xkl
      ydrv(k) = ydrv(k) - dfangledrkl*ykl
      zdrv(k) = zdrv(k) - dfangledrkl*zkl
      xdrv(l) = xdrv(l) + dfangledrkl*xkl
      ydrv(l) = ydrv(l) + dfangledrkl*ykl
      zdrv(l) = zdrv(l) + dfangledrkl*zkl
    endif
!
    if (lstr.or.lgrad2) then
      call real1strterm(ndim,xjk,yjk,zjk,0.0_dp,0.0_dp,0.0_dp,dr2ds3a(1,1),d2r2dx23a(1,1), &
                        d2r2dsdx3a(1,1,1),d2r2ds23a(1,1,1),lgrad2)
      call real1strterm(ndim,xjl,yjl,zjl,0.0_dp,0.0_dp,0.0_dp,dr2ds3a(1,2),d2r2dx23a(1,2), &
                        d2r2dsdx3a(1,1,2),d2r2ds23a(1,1,2),lgrad2)
      call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds3a(1,3),d2r2dx23a(1,3), &
                        d2r2dsdx3a(1,1,3),d2r2ds23a(1,1,3),lgrad2)
    endif
    if (lstr) then
      if (lilocal) then
        do kl = 1,nstrains
          ks = nstrptr(kl)
          rstrd(kl) = rstrd(kl) + dfangledrjk*dr2ds3a(ks,1) + dfangledrjl*dr2ds3a(ks,2) + dfangledrkl*dr2ds3a(ks,3)
        enddo
      endif
    endif
!
    if (lgrad2) then
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
      call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,dfangledrjk, &
                        d2trm,dr2ds3a(1,1),d2r2dx23a(1,1),d2r2ds23a(1,1,1),d2r2dsdx3a(1,1,1),lilocal)
!
!  j-k / j-l
!
      d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(2) + d2fangledcosa2*dcosdrjk*dcosdrjl)
      call add_drv2_2dm(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm, &
                        dr2ds3a(1,1),dr2ds3a(1,2),lilocal)
!
!  j-k / k-l
!
      d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(3) + d2fangledcosa2*dcosdrjk*dcosdrkl)
      call add_drv2_2dm(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm, &
                        dr2ds3a(1,1),dr2ds3a(1,3),lilocal)
!
!  j-l / j-l
!
      d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(4) + d2fangledcosa2*dcosdrjl*dcosdrjl)
      call add_drv2_1dm(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,dfangledrjl, &
                        d2trm,dr2ds3a(1,2),d2r2dx23a(1,2),d2r2ds23a(1,1,2),d2r2dsdx3a(1,1,2),lilocal)
!
!  j-l / k-l
!
      d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(5) + d2fangledcosa2*dcosdrjl*dcosdrkl)
      call add_drv2_2dm(ljlocal,lllocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xkl,ykl,zkl,d2trm, &
                        dr2ds3a(1,2),dr2ds3a(1,3),lilocal)
!
!  k-l / k-l
!
      d2trm = - const*ftors*rdamp*qhdampo*(dfangledcosa*d2cos(6) + d2fangledcosa2*dcosdrkl*dcosdrkl)
      call add_drv2_1dm(lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,dfangledrkl, &
                        d2trm,dr2ds3a(1,3),d2r2dx23a(1,3),d2r2ds23a(1,1,3),d2r2dsdx3a(1,1,3),lilocal)
!--------------------------------------------------------
!  Second derivatives of fangle mixed with other terms  |
!--------------------------------------------------------
!
!  i-j / j-k
!
      d2trm = - const*ftors*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*dfangledcosa*dcosdrjk
      call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm, &
                        dr2ds3(1,1),dr2ds3a(1,1),lilocal)
!
!  i-j / j-l
!
      d2trm = - const*ftors*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*dfangledcosa*dcosdrjl
      call add_drv2_2dm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm, &
                        dr2ds3(1,1),dr2ds3a(1,2),lilocal)
!
!  i-j / k-l
!
      d2trm = - const*ftors*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*dfangledcosa*dcosdrkl
      call add_drv2_2dm(lilocal,ljlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xkl,ykl,zkl,d2trm, &
                        dr2ds3(1,1),dr2ds3a(1,3),lilocal)
!
!  i-k / j-k
!
      d2trm = - const*ftors*dampo_nb_tot*rdamp*ddampodrik*dfangledcosa*dcosdrjk
      call add_drv2_2dm(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm, &
                        dr2ds3(1,2),dr2ds3a(1,1),lilocal)
!
!  i-k / j-l
!
      d2trm = - const*ftors*dampo_nb_tot*rdamp*ddampodrik*dfangledcosa*dcosdrjl
      call add_drv2_2dm(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xjl,yjl,zjl,d2trm, &
                        dr2ds3(1,2),dr2ds3a(1,2),lilocal)
!
!  i-k / k-l
!
      d2trm = - const*ftors*dampo_nb_tot*rdamp*ddampodrik*dfangledcosa*dcosdrkl
      call add_drv2_2dm(lilocal,lklocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xkl,ykl,zkl,d2trm, &
                        dr2ds3(1,2),dr2ds3a(1,3),lilocal)
!
!  j-k / j-k
!
      d2trm = - const*ftors*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*dfangledcosa*dcosdrjk
      call add_drv2_2dm(ljlocal,lklocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,xjk,yjk,zjk,d2trm, &
                        dr2ds3(1,3),dr2ds3a(1,1),lilocal)
!
!  j-k / j-l
!
      d2trm = - const*ftors*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*dfangledcosa*dcosdrjl
      call add_drv2_2dm(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm, &
                        dr2ds3(1,3),dr2ds3a(1,2),lilocal)
!
!  j-k / k-l
!
      d2trm = - const*ftors*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*dfangledcosa*dcosdrkl
      call add_drv2_2dm(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                        kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm, &
                        dr2ds3(1,3),dr2ds3a(1,3),lilocal)
    endif
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
      if (lilocal) then
        xdrv(i) = xdrv(i) - deijkdrij*xij - deijkdrim*xim
        ydrv(i) = ydrv(i) - deijkdrij*yij - deijkdrim*yim
        zdrv(i) = zdrv(i) - deijkdrij*zij - deijkdrim*zim
!
        xdrv(j) = xdrv(j) + deijkdrij*xij - deijkdrjm*xjm
        ydrv(j) = ydrv(j) + deijkdrij*yij - deijkdrjm*yjm
        zdrv(j) = zdrv(j) + deijkdrij*zij - deijkdrjm*zjm
!
        xdrv(m) = xdrv(m) + deijkdrim*xim + deijkdrjm*xjm
        ydrv(m) = ydrv(m) + deijkdrim*yim + deijkdrjm*yjm
        zdrv(m) = zdrv(m) + deijkdrim*zim + deijkdrjm*zjm
      endif
!
      if (lstr.or.lgrad2) then
        call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3m(1,1),d2r2dx23m(1,1),d2r2dsdx3m(1,1,1), &
                          d2r2ds23m(1,1,1),lgrad2)
        call real1strterm(ndim,xim,yim,zim,0.0_dp,0.0_dp,0.0_dp,dr2ds3m(1,2),d2r2dx23m(1,2),d2r2dsdx3m(1,1,2), &
                          d2r2ds23m(1,1,2),lgrad2)
        call real1strterm(ndim,xjm,yjm,zjm,0.0_dp,0.0_dp,0.0_dp,dr2ds3m(1,3),d2r2dx23m(1,3),d2r2dsdx3m(1,1,3), &
                          d2r2ds23m(1,1,3),lgrad2)
      endif
      if (lstr) then
        if (lilocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + deijkdrij*dr2ds3m(ks,1) + deijkdrim*dr2ds3m(ks,2) + deijkdrjm*dr2ds3m(ks,3)
          enddo
        endif
      endif
!
      if (lgrad2) then
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
        call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,deijkdrij, &
                          d2trm,dr2ds3m(1,1),d2r2dx23m(1,1),d2r2ds23m(1,1,1),d2r2dsdx3m(1,1,1),lilocal)
!
!  i-j / i-m
!
        d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(2,nj)
        call add_drv2_2dm(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm, &
                          dr2ds3m(1,1),dr2ds3m(1,2),lilocal)
!
!  i-j / j-m
!
        d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(3,nj)
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm, &
                          dr2ds3m(1,1),dr2ds3m(1,3),lilocal)
!
!  i-m / i-m
!
        d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(4,nj)
        call add_drv2_1dm(lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xim,yim,zim,deijkdrim, &
                          d2trm,dr2ds3m(1,2),d2r2dx23m(1,2),d2r2ds23m(1,1,2),d2r2dsdx3m(1,1,2),lilocal)
!
!  i-m / j-m
!
        d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(5,nj)
        call add_drv2_2dm(lilocal,lmlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xim,yim,zim,xjm,yjm,zjm,d2trm, &
                          dr2ds3m(1,2),dr2ds3m(1,3),lilocal)
!
!  j-m / j-m
!
        d2trm = - fconst*rdamp*dampo*dampo_nb_totm1*d2dampo_nbdr2(6,nj)
        call add_drv2_1dm(ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjm,yjm,zjm,deijkdrjm, &
                          d2trm,dr2ds3m(1,3),d2r2dx23m(1,3),d2r2ds23m(1,1,3),d2r2dsdx3m(1,1,3),lilocal)
!------------------------------------------------
!  Contribution from dampo_nb and 1 other term  |
!------------------------------------------------
!
!  i-j / i-j
!
        d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(1,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
        call add_drv2_2dm(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm, &
                          dr2ds3(1,1),dr2ds3m(1,1),lilocal)
!
!  i-j / i-m
!
        d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(2,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
        call add_drv2_2dm(lilocal,ljlocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xim,yim,zim,d2trm, &
                          dr2ds3(1,1),dr2ds3m(1,2),lilocal)
!
!  i-j / j-m
!
        d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(3,nj)*(rdamp*ddampodrij + drdampdrij*dampo)
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm, &
                          dr2ds3(1,1),dr2ds3m(1,3),lilocal)
!
!  i-k / i-j
!
        d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(1,nj)*rdamp*ddampodrik
        call add_drv2_2dm(lilocal,lklocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xik,yik,zik,xij,yij,zij,d2trm, &
                          dr2ds3(1,2),dr2ds3m(1,1),lilocal)
!
!  i-k / i-m
!
        d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(2,nj)*rdamp*ddampodrik
        call add_drv2_2dm(lilocal,lklocal,lilocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xim,yim,zim,d2trm, &
                          dr2ds3(1,2),dr2ds3m(1,2),lilocal)
!
!  i-k / j-m
!
        d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(3,nj)*rdamp*ddampodrik
        call add_drv2_2dm(lilocal,lklocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xjm,yjm,zjm,d2trm, &
                          dr2ds3(1,2),dr2ds3m(1,3),lilocal)
!
!  j-k / i-j
!
        d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(1,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
        call add_drv2_2dm(ljlocal,lklocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjk,yjk,zjk,xij,yij,zij,d2trm, &
                          dr2ds3(1,3),dr2ds3m(1,1),lilocal)
!
!  j-k / i-m
!
        d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(2,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
        call add_drv2_2dm(ljlocal,lklocal,lilocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xim,yim,zim,d2trm, &
                          dr2ds3(1,3),dr2ds3m(1,2),lilocal)
!
!  j-k / j-m
!
        d2trm = - fconst*dampo_nb_totm1*ddampo_nbdr(3,nj)*(rdamp*ddampodrjk + drdampdrjk*dampo)
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,d2trm, &
                          dr2ds3(1,3),dr2ds3m(1,3),lilocal)
!------------------------------------------
!  Contribution from dampo_nb and fangle  |
!------------------------------------------
!
!  i-j / j-k
!
        d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(1,nj)*dfangledcosa*dcosdrjk
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm, &
                          dr2ds3m(1,1),dr2ds3a(1,1),lilocal)
!
!  i-j / j-l
!
        d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(1,nj)*dfangledcosa*dcosdrjl
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm, &
                          dr2ds3m(1,1),dr2ds3a(1,2),lilocal)
!
!  i-j / k-l
!
        d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(1,nj)*dfangledcosa*dcosdrkl
        call add_drv2_2dm(lilocal,ljlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xkl,ykl,zkl,d2trm, &
                          dr2ds3m(1,1),dr2ds3a(1,3),lilocal)
!
!  i-m / j-k
!
        d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(2,nj)*dfangledcosa*dcosdrjk
        call add_drv2_2dm(lilocal,lmlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xim,yim,zim,xjk,yjk,zjk,d2trm, &
                          dr2ds3m(1,2),dr2ds3a(1,1),lilocal)
!
!  i-m / j-l
!
        d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(2,nj)*dfangledcosa*dcosdrjl
        call add_drv2_2dm(lilocal,lmlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xim,yim,zim,xjl,yjl,zjl,d2trm, &
                          dr2ds3m(1,2),dr2ds3a(1,2),lilocal)
!
!  i-m / k-l
!
        d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(2,nj)*dfangledcosa*dcosdrkl
        call add_drv2_2dm(lilocal,lmlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xim,yim,zim,xkl,ykl,zkl,d2trm, &
                          dr2ds3m(1,2),dr2ds3a(1,3),lilocal)
!
!  j-m / j-k
!
        d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(3,nj)*dfangledcosa*dcosdrjk
        call add_drv2_2dm(ljlocal,lmlocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjm,yjm,zjm,xjk,yjk,zjk,d2trm, &
                          dr2ds3m(1,3),dr2ds3a(1,1),lilocal)
!
!  j-m / j-l
!
        d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(3,nj)*dfangledcosa*dcosdrjl
        call add_drv2_2dm(ljlocal,lmlocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjm,yjm,zjm,xjl,yjl,zjl,d2trm, &
                          dr2ds3m(1,3),dr2ds3a(1,2),lilocal)
!
!  j-m / k-l
!
        d2trm = - const*ftors*dampo_nb_totm1*rdamp*dampo*ddampo_nbdr(3,nj)*dfangledcosa*dcosdrkl
        call add_drv2_2dm(ljlocal,lmlocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjm,yjm,zjm,xkl,ykl,zkl,d2trm, &
                          dr2ds3m(1,3),dr2ds3a(1,3),lilocal)
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
          if (lstr) then
            call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3n(1,1),d2r2dx23n(1,1), &
                              d2r2dsdx3n(1,1,1),d2r2ds23n(1,1,1),.false.)
            call real1strterm(ndim,xin,yin,zin,0.0_dp,0.0_dp,0.0_dp,dr2ds3n(1,2),d2r2dx23n(1,2), &
                              d2r2dsdx3n(1,1,2),d2r2ds23n(1,1,2),.false.)
            call real1strterm(ndim,xjn,yjn,zjn,0.0_dp,0.0_dp,0.0_dp,dr2ds3n(1,3),d2r2dx23n(1,3), &
                              d2r2dsdx3n(1,1,3),d2r2ds23n(1,1,3),.false.)
          endif
!
!  i-j / i-j
!
          d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(1,nk)
          call add_drv2_2dm(lilocal,ljlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xij,yij,zij,xij,yij,zij,d2trm, &
                            dr2ds3m(1,1),dr2ds3n(1,1),lilocal)
!
!  i-j / i-n
!
          d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(2,nk)
          call add_drv2_2dm(lilocal,ljlocal,lilocal,lnlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xij,yij,zij,xin,yin,zin,d2trm, &
                            dr2ds3m(1,1),dr2ds3n(1,2),lilocal)
!
!  i-j / j-n
!
          d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(1,nj)*ddampo_nbdr(3,nk)
          call add_drv2_2dm(lilocal,ljlocal,ljlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xij,yij,zij,xjn,yjn,zjn,d2trm, &
                            dr2ds3m(1,1),dr2ds3n(1,3),lilocal)
!
!  i-m / i-j
!
          d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(1,nk)
          call add_drv2_2dm(lilocal,lmlocal,lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xim,yim,zim,xij,yij,zij,d2trm, &
                            dr2ds3m(1,2),dr2ds3n(1,1),lilocal)
!
!  i-m / i-n
!
          d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(2,nk)
          call add_drv2_2dm(lilocal,lmlocal,lilocal,lnlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                            ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xim,yim,zim,xin,yin,zin,d2trm, &
                            dr2ds3m(1,2),dr2ds3n(1,2),lilocal)
!
!  i-m / j-n
!
          d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(2,nj)*ddampo_nbdr(3,nk)
          call add_drv2_2dm(lilocal,lmlocal,ljlocal,lnlocal,ix,iy,iz,ixf,iyf,izf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xim,yim,zim,xjn,yjn,zjn,d2trm, &
                            dr2ds3m(1,2),dr2ds3n(1,3),lilocal)
!
!  j-m / i-j
!
          d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(1,nk)
          call add_drv2_2dm(ljlocal,lmlocal,lilocal,ljlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xjm,yjm,zjm,xij,yij,zij,d2trm, &
                            dr2ds3m(1,3),dr2ds3n(1,1),lilocal)
!
!  j-m / i-n
!
          d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(2,nk)
          call add_drv2_2dm(ljlocal,lmlocal,lilocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xin,yin,zin,d2trm, &
                            dr2ds3m(1,3),dr2ds3n(1,2),lilocal)
!
!  j-m / j-n
!
          d2trm = - fconst*rdamp*dampo*dampo_nb_totm2*ddampo_nbdr(3,nj)*ddampo_nbdr(3,nk)
          call add_drv2_2dm(ljlocal,lmlocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xjn,yjn,zjn,d2trm, &
                            dr2ds3m(1,3),dr2ds3n(1,3),lilocal)
        enddo
      endif
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
      if (lilocal) then
        xdrv(j) = xdrv(j) - xjk*e1d(1) - xjl*e1d(4) - xjm*e1d(5)
        ydrv(j) = ydrv(j) - yjk*e1d(1) - yjl*e1d(4) - yjm*e1d(5)
        zdrv(j) = zdrv(j) - zjk*e1d(1) - zjl*e1d(4) - zjm*e1d(5)
!
        xdrv(k) = xdrv(k) - xkl*e1d(2) + xjk*e1d(1) - xkm*e1d(3)
        ydrv(k) = ydrv(k) - ykl*e1d(2) + yjk*e1d(1) - ykm*e1d(3)
        zdrv(k) = zdrv(k) - zkl*e1d(2) + zjk*e1d(1) - zkm*e1d(3)
!
        xdrv(l) = xdrv(l) + xkl*e1d(2) - xlm*e1d(6) + xjl*e1d(4)
        ydrv(l) = ydrv(l) + ykl*e1d(2) - ylm*e1d(6) + yjl*e1d(4)
        zdrv(l) = zdrv(l) + zkl*e1d(2) - zlm*e1d(6) + zjl*e1d(4)
!
        xdrv(m) = xdrv(m) + xlm*e1d(6) + xkm*e1d(3) + xjm*e1d(5)
        ydrv(m) = ydrv(m) + ylm*e1d(6) + ykm*e1d(3) + yjm*e1d(5)
        zdrv(m) = zdrv(m) + zlm*e1d(6) + zkm*e1d(3) + zjm*e1d(5)
      endif
!
      if (lstr.or.lgrad2) then
        call real1strterm(ndim,xjk,yjk,zjk,0.0_dp,0.0_dp,0.0_dp,dr2ds6t(1,1),d2r2dx26t(1,1),d2r2dsdx6t(1,1,1), &
                          d2r2ds26t(1,1,1),lgrad2)
        call real1strterm(ndim,xkl,ykl,zkl,0.0_dp,0.0_dp,0.0_dp,dr2ds6t(1,2),d2r2dx26t(1,2),d2r2dsdx6t(1,1,2), &
                          d2r2ds26t(1,1,2),lgrad2)
        call real1strterm(ndim,xkm,ykm,zkm,0.0_dp,0.0_dp,0.0_dp,dr2ds6t(1,3),d2r2dx26t(1,3),d2r2dsdx6t(1,1,3), &
                          d2r2ds26t(1,1,3),lgrad2)
        call real1strterm(ndim,xjl,yjl,zjl,0.0_dp,0.0_dp,0.0_dp,dr2ds6t(1,4),d2r2dx26t(1,4),d2r2dsdx6t(1,1,4), &
                          d2r2ds26t(1,1,4),lgrad2)
        call real1strterm(ndim,xjm,yjm,zjm,0.0_dp,0.0_dp,0.0_dp,dr2ds6t(1,5),d2r2dx26t(1,5),d2r2dsdx6t(1,1,5), &
                          d2r2ds26t(1,1,5),lgrad2)
        call real1strterm(ndim,xlm,ylm,zlm,0.0_dp,0.0_dp,0.0_dp,dr2ds6t(1,6),d2r2dx26t(1,6),d2r2dsdx6t(1,1,6), &
                          d2r2ds26t(1,1,6),lgrad2)
      endif
      if (lstr) then
        if (lilocal) then
          do kl = 1,nstrains
            ks = nstrptr(kl)
            rstrd(kl) = rstrd(kl) + e1d(1)*dr2ds6t(ks,1)
            rstrd(kl) = rstrd(kl) + e1d(2)*dr2ds6t(ks,2)
            rstrd(kl) = rstrd(kl) + e1d(3)*dr2ds6t(ks,3)
            rstrd(kl) = rstrd(kl) + e1d(4)*dr2ds6t(ks,4)
            rstrd(kl) = rstrd(kl) + e1d(5)*dr2ds6t(ks,5)
            rstrd(kl) = rstrd(kl) + e1d(6)*dr2ds6t(ks,6)
          enddo
        endif
      endif
      if (lgrad2) then
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
        call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,e1d(1), &
                          e2d(1),dr2ds6t(1,1),d2r2dx26t(1,1),d2r2ds26t(1,1,1),d2r2dsdx6t(1,1,1),lilocal)
!
!  j-k / k-l
!
        call add_drv2_2dm(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,e2d(2), &
                          dr2ds6t(1,1),dr2ds6t(1,2),lilocal)
!
!  j-k / k-m
!
        call add_drv2_2dm(ljlocal,lklocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xkm,ykm,zkm,e2d(3), &
                          dr2ds6t(1,1),dr2ds6t(1,3),lilocal)
!
!  j-k / j-l
!
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,e2d(4), &
                          dr2ds6t(1,1),dr2ds6t(1,4),lilocal)
!
!  j-k / j-m
!
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,e2d(5), &
                          dr2ds6t(1,1),dr2ds6t(1,5),lilocal)
!
!  j-k / l-m
!
        call add_drv2_2dm(ljlocal,lklocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xlm,ylm,zlm,e2d(6), &
                          dr2ds6t(1,1),dr2ds6t(1,6),lilocal)
!
!  k-l / k-l
!
        call add_drv2_1dm(lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,e1d(2), &
                          e2d(7),dr2ds6t(1,2),d2r2dx26t(1,2),d2r2ds26t(1,1,2),d2r2dsdx6t(1,1,2),lilocal)
!
!  k-l / k-m
!
        call add_drv2_2dm(lklocal,lllocal,lklocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xkm,ykm,zkm,e2d(8), &
                          dr2ds6t(1,2),dr2ds6t(1,3),lilocal)
!
!  k-l / j-l
!
        call add_drv2_2dm(lklocal,lllocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xjl,yjl,zjl,e2d(9), &
                          dr2ds6t(1,2),dr2ds6t(1,4),lilocal)
!
!  k-l / j-m
!
        call add_drv2_2dm(lklocal,lllocal,ljlocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xjm,yjm,zjm,e2d(10), &
                          dr2ds6t(1,2),dr2ds6t(1,5),lilocal)
!
!  k-l / l-m
!
        call add_drv2_2dm(lklocal,lllocal,lllocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xlm,ylm,zlm,e2d(11), &
                          dr2ds6t(1,2),dr2ds6t(1,6),lilocal)
!
!  k-m / k-m
!
        call add_drv2_1dm(lklocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xkm,ykm,zkm,e1d(3), &
                          e2d(12),dr2ds6t(1,3),d2r2dx26t(1,3),d2r2ds26t(1,1,3),d2r2dsdx6t(1,1,3),lilocal)
!
!  k-m / j-l
!
        call add_drv2_2dm(lklocal,lmlocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkm,ykm,zkm,xjl,yjl,zjl,e2d(13), &
                          dr2ds6t(1,3),dr2ds6t(1,4),lilocal)
!
!  k-m / j-m
!
        call add_drv2_2dm(lklocal,lmlocal,ljlocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xkm,ykm,zkm,xjm,yjm,zjm,e2d(14), &
                          dr2ds6t(1,3),dr2ds6t(1,5),lilocal)
!     
!  k-m / l-m
! 
        call add_drv2_2dm(lklocal,lmlocal,lllocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xkm,ykm,zkm,xlm,ylm,zlm,e2d(15), &
                          dr2ds6t(1,3),dr2ds6t(1,6),lilocal)
!     
!  j-l / j-l
! 
        call add_drv2_1dm(ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,e1d(4), &
                          e2d(16),dr2ds6t(1,4),d2r2dx26t(1,4),d2r2ds26t(1,1,4),d2r2dsdx6t(1,1,4),lilocal)
!
!  j-l / j-m
!
        call add_drv2_2dm(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xjm,yjm,zjm,e2d(17), &
                          dr2ds6t(1,4),dr2ds6t(1,5),lilocal)
!
!  j-l / l-m
!
        call add_drv2_2dm(ljlocal,lllocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xlm,ylm,zlm,e2d(18), &
                          dr2ds6t(1,4),dr2ds6t(1,6),lilocal)
!
!  j-m / j-m
!
        call add_drv2_1dm(ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjm,yjm,zjm,e1d(5), &
                          e2d(19),dr2ds6t(1,5),d2r2dx26t(1,5),d2r2ds26t(1,1,5),d2r2dsdx6t(1,1,5),lilocal)
!
!  j-m / l-m
!
        call add_drv2_2dm(ljlocal,lmlocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjm,yjm,zjm,xlm,ylm,zlm,e2d(20), &
                          dr2ds6t(1,5),dr2ds6t(1,6),lilocal)
!
!  l-m / l-m
!
        call add_drv2_1dm(lllocal,lmlocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xlm,ylm,zlm,e1d(6), &
                          e2d(21),dr2ds6t(1,6),d2r2dx26t(1,6),d2r2ds26t(1,1,6),d2r2dsdx6t(1,1,6),lilocal)
!-------------------------------------------------------
!  Second derivatives of ftors mixed with other terms  |
!-------------------------------------------------------
!
!  i-j / j-k
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(1,nt)
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm, &
                          dr2ds3(1,1),dr2ds6t(1,1),lilocal)
!
!  i-j / k-l
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(2,nt)
        call add_drv2_2dm(lilocal,ljlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xkl,ykl,zkl,d2trm, &
                          dr2ds3(1,1),dr2ds6t(1,2),lilocal)
!
!  i-j / k-m
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(3,nt)
        call add_drv2_2dm(lilocal,ljlocal,lklocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xkm,ykm,zkm,d2trm, &
                          dr2ds3(1,1),dr2ds6t(1,3),lilocal)
!
!  i-j / j-l
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(4,nt)
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm, &
                          dr2ds3(1,1),dr2ds6t(1,4),lilocal)
!
!  i-j / j-m
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(5,nt)
        call add_drv2_2dm(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm, &
                          dr2ds3(1,1),dr2ds6t(1,5),lilocal)
!
!  i-j / l-m
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrij*dampo + rdamp*ddampodrij)*ftorsm1*dft(6,nt)
        call add_drv2_2dm(lilocal,ljlocal,lllocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xlm,ylm,zlm,d2trm, &
                          dr2ds3(1,1),dr2ds6t(1,6),lilocal)
!
!  i-k / j-k
!
        d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(1,nt)
        call add_drv2_2dm(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xik,yik,zik,xjk,yjk,zjk,d2trm, &
                          dr2ds3(1,2),dr2ds6t(1,1),lilocal)
!
!  i-k / k-l
!
        d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(2,nt)
        call add_drv2_2dm(lilocal,lklocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xkl,ykl,zkl,d2trm, &
                          dr2ds3(1,2),dr2ds6t(1,2),lilocal)
!
!  i-k / k-m
!
        d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(3,nt)
        call add_drv2_2dm(lilocal,lklocal,lklocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xkm,ykm,zkm,d2trm, &
                          dr2ds3(1,2),dr2ds6t(1,3),lilocal)
!
!  i-k / j-l
!
        d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(4,nt)
        call add_drv2_2dm(lilocal,lklocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xik,yik,zik,xjl,yjl,zjl,d2trm, &
                          dr2ds3(1,2),dr2ds6t(1,4),lilocal)
!
!  i-k / j-m
!
        d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(5,nt)
        call add_drv2_2dm(lilocal,lklocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xjm,yjm,zjm,d2trm, &
                          dr2ds3(1,2),dr2ds6t(1,5),lilocal)
!
!  i-k / l-m
!
        d2trm = - const*fangle*dampo_nb_tot*rdamp*ddampodrik*ftorsm1*dft(6,nt)
        call add_drv2_2dm(lilocal,lklocal,lllocal,lmlocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xik,yik,zik,xlm,ylm,zlm,d2trm, &
                          dr2ds3(1,2),dr2ds6t(1,6),lilocal)
!
!  j-k / j-k
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(1,nt)
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,xjk,yjk,zjk,d2trm, &
                          dr2ds3(1,3),dr2ds6t(1,1),lilocal)
!
!  j-k / k-l
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(2,nt)
        call add_drv2_2dm(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm, &
                          dr2ds3(1,3),dr2ds6t(1,2),lilocal)
!
!  j-k / k-m
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(3,nt)
        call add_drv2_2dm(ljlocal,lklocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xkm,ykm,zkm,d2trm, &
                          dr2ds3(1,3),dr2ds6t(1,3),lilocal)
!
!  j-k / j-l
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(4,nt)
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm, &
                          dr2ds3(1,3),dr2ds6t(1,4),lilocal)
!
!  j-k / j-m
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(5,nt)
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,d2trm, &
                          dr2ds3(1,3),dr2ds6t(1,5),lilocal)
!
!  j-k / l-m
!
        d2trm = - const*fangle*dampo_nb_tot*(drdampdrjk*dampo + rdamp*ddampodrjk)*ftorsm1*dft(6,nt)
        call add_drv2_2dm(ljlocal,lklocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xlm,ylm,zlm,d2trm, &
                          dr2ds3(1,3),dr2ds6t(1,6),lilocal)
!--------------------------------------------------
!  Second derivatives of ftors mixed with fangle  |
!--------------------------------------------------
!
!  j-k / j-k
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(1,nt)
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,xjk,yjk,zjk,d2trm, &
                          dr2ds3a(1,1),dr2ds6t(1,1),lilocal)
!
!  j-k / k-l
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(2,nt)
        call add_drv2_2dm(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm, &
                          dr2ds3a(1,1),dr2ds6t(1,2),lilocal)
!
!  j-k / k-m
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(3,nt)
        call add_drv2_2dm(ljlocal,lklocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xkm,ykm,zkm,d2trm, &
                          dr2ds3a(1,1),dr2ds6t(1,3),lilocal)
!
!  j-k / j-l
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(4,nt)
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm, &
                          dr2ds3a(1,1),dr2ds6t(1,4),lilocal)
!
!  j-k / j-m
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(5,nt)
        call add_drv2_2dm(ljlocal,lklocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xjm,yjm,zjm,d2trm, &
                          dr2ds3a(1,1),dr2ds6t(1,5),lilocal)
!
!  j-k / l-m
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjk*ftorsm1*dft(6,nt)
        call add_drv2_2dm(ljlocal,lklocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjk,yjk,zjk,xlm,ylm,zlm,d2trm, &
                          dr2ds3a(1,1),dr2ds6t(1,6),lilocal)
!
!  j-l / j-k
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(1,nt)
        call add_drv2_2dm(ljlocal,lllocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjl,yjl,zjl,xjk,yjk,zjk,d2trm, &
                          dr2ds3a(1,2),dr2ds6t(1,1),lilocal)
!
!  j-l / k-l
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(2,nt)
        call add_drv2_2dm(ljlocal,lllocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xkl,ykl,zkl,d2trm, &
                          dr2ds3a(1,2),dr2ds6t(1,2),lilocal)
!
!  j-l / k-m 
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(3,nt)
        call add_drv2_2dm(ljlocal,lllocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xkm,ykm,zkm,d2trm, &
                          dr2ds3a(1,2),dr2ds6t(1,3),lilocal)
!
!  j-l / j-l
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(4,nt)
        call add_drv2_2dm(ljlocal,lllocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xjl,yjl,zjl,d2trm, &
                          dr2ds3a(1,2),dr2ds6t(1,4),lilocal)
!
!  j-l / j-m
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(5,nt)
        call add_drv2_2dm(ljlocal,lllocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xjm,yjm,zjm,d2trm, &
                          dr2ds3a(1,2),dr2ds6t(1,5),lilocal)
!
!  j-l / l-m
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrjl*ftorsm1*dft(6,nt)
        call add_drv2_2dm(ljlocal,lllocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjl,yjl,zjl,xlm,ylm,zlm,d2trm, &
                          dr2ds3a(1,2),dr2ds6t(1,6),lilocal)
!
!  k-l / j-k
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(1,nt)
        call add_drv2_2dm(lklocal,lllocal,ljlocal,lklocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xkl,ykl,zkl,xjk,yjk,zjk,d2trm, &
                          dr2ds3a(1,3),dr2ds6t(1,1),lilocal)
!
!  k-l / k-l
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(2,nt)
        call add_drv2_2dm(lklocal,lllocal,lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xkl,ykl,zkl,d2trm, &
                          dr2ds3a(1,3),dr2ds6t(1,2),lilocal)
!
!  k-l / k-m
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(3,nt)
        call add_drv2_2dm(lklocal,lllocal,lklocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xkm,ykm,zkm,d2trm, &
                          dr2ds3a(1,3),dr2ds6t(1,3),lilocal)
!
!  k-l / j-l
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(4,nt)
        call add_drv2_2dm(lklocal,lllocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xjl,yjl,zjl,d2trm, &
                          dr2ds3a(1,3),dr2ds6t(1,4),lilocal)
!
!  k-l / j-m
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(5,nt)
        call add_drv2_2dm(lklocal,lllocal,ljlocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xjm,yjm,zjm,d2trm, &
                          dr2ds3a(1,3),dr2ds6t(1,5),lilocal)
!
!  k-l / l-m
!
        d2trm = - const*rdamp*qhdampo*dfangledcosa*dcosdrkl*ftorsm1*dft(6,nt)
        call add_drv2_2dm(lklocal,lllocal,lllocal,lmlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                          lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xkl,ykl,zkl,xlm,ylm,zlm,d2trm, &
                          dr2ds3a(1,3),dr2ds6t(1,6),lilocal)
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
          call real1strterm(ndim,xkn,ykn,zkn,0.0_dp,0.0_dp,0.0_dp,dr2ds6t2(1,1),d2r2dx26t(1,1),d2r2dsdx6t(1,1,1), &
                            d2r2ds26t(1,1,1),.false.)
          call real1strterm(ndim,xjn,yjn,zjn,0.0_dp,0.0_dp,0.0_dp,dr2ds6t2(1,2),d2r2dx26t(1,2),d2r2dsdx6t(1,1,2), &
                            d2r2ds26t(1,1,2),.false.)
          call real1strterm(ndim,xln,yln,zln,0.0_dp,0.0_dp,0.0_dp,dr2ds6t2(1,3),d2r2dx26t(1,3),d2r2dsdx6t(1,1,3), &
                            d2r2ds26t(1,1,3),.false.)
!
!  j-k / j-k
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(1,nt2)
          call add_drv2_2dm(ljlocal,lklocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjk,yjk,zjk,xjk,yjk,zjk,d2trm, &
                            dr2ds6t(1,1),dr2ds6t(1,1),lilocal)
!
!  j-k / k-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(2,nt2)
          call add_drv2_2dm(ljlocal,lklocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xkl,ykl,zkl,d2trm, &
                            dr2ds6t(1,1),dr2ds6t(1,2),lilocal)
!
!  j-k / k-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(3,nt2)
          call add_drv2_2dm(ljlocal,lklocal,lklocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xjk,yjk,zjk,xkn,ykn,zkn,d2trm, &
                            dr2ds6t(1,1),dr2ds6t2(1,1),lilocal)
!
!  j-k / j-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(4,nt2)
          call add_drv2_2dm(ljlocal,lklocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjk,yjk,zjk,xjl,yjl,zjl,d2trm, &
                            dr2ds6t(1,1),dr2ds6t(1,4),lilocal)
!
!  j-k / j-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(5,nt2)
          call add_drv2_2dm(ljlocal,lklocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjk,yjk,zjk,xjn,yjn,zjn,d2trm, &
                            dr2ds6t(1,1),dr2ds6t2(1,2),lilocal)
!
!  j-k / l-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(1,nt)*dft(6,nt2)
          call add_drv2_2dm(ljlocal,lklocal,lllocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf, &
                            lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xjk,yjk,zjk,xln,yln,zln,d2trm, &
                            dr2ds6t(1,1),dr2ds6t2(1,3),lilocal)
!
!  k-l / j-k
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(1,nt2)
          call add_drv2_2dm(lklocal,lllocal,ljlocal,lklocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xkl,ykl,zkl,xjk,yjk,zjk,d2trm, &
                            dr2ds6t(1,2),dr2ds6t(1,1),lilocal)
!
!  k-l / k-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(2,nt2)
          call add_drv2_2dm(lklocal,lllocal,lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                            kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xkl,ykl,zkl,d2trm, &
                            dr2ds6t(1,2),dr2ds6t(1,2),lilocal)
!
!  k-l / k-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(3,nt2)
          call add_drv2_2dm(lklocal,lllocal,lklocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                            kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xkn,ykn,zkn,d2trm, &
                            dr2ds6t(1,2),dr2ds6t2(1,1),lilocal)
!
!  k-l / j-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(4,nt2)
          call add_drv2_2dm(lklocal,lllocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkl,ykl,zkl,xjl,yjl,zjl,d2trm, &
                            dr2ds6t(1,2),dr2ds6t(1,4),lilocal)
!
!  k-l / j-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(5,nt2)
          call add_drv2_2dm(lklocal,lllocal,ljlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xjn,yjn,zjn,d2trm, &
                            dr2ds6t(1,2),dr2ds6t2(1,2),lilocal)
!
!  k-l / l-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(2,nt)*dft(6,nt2)
          call add_drv2_2dm(lklocal,lllocal,lllocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf, &
                            lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xkl,ykl,zkl,xln,yln,zln,d2trm, &
                            dr2ds6t(1,2),dr2ds6t2(1,3),lilocal)
!
!  k-m / j-k
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(1,nt2)
          call add_drv2_2dm(lklocal,lmlocal,ljlocal,lklocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xkm,ykm,zkm,xjk,yjk,zjk,d2trm, &
                            dr2ds6t(1,3),dr2ds6t(1,1),lilocal)
!
!  k-m / k-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(2,nt2)
          call add_drv2_2dm(lklocal,lmlocal,lklocal,lllocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                            kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xkm,ykm,zkm,xkl,ykl,zkl,d2trm, &
                            dr2ds6t(1,3),dr2ds6t(1,2),lilocal)
!
!  k-m / k-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(3,nt2)
          call add_drv2_2dm(lklocal,lmlocal,lklocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                            kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xkm,ykm,zkm,xkn,ykn,zkn,d2trm, &
                            dr2ds6t(1,3),dr2ds6t2(1,1),lilocal)
!
!  k-m / j-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(4,nt2)
          call add_drv2_2dm(lklocal,lmlocal,ljlocal,lllocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xkm,ykm,zkm,xjl,yjl,zjl,d2trm, &
                            dr2ds6t(1,3),dr2ds6t(1,4),lilocal)
!
!  k-m / j-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(5,nt2)
          call add_drv2_2dm(lklocal,lmlocal,ljlocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xkm,ykm,zkm,xjn,yjn,zjn,d2trm, &
                            dr2ds6t(1,3),dr2ds6t2(1,2),lilocal)
!
!  k-m / l-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(3,nt)*dft(6,nt2)
          call add_drv2_2dm(lklocal,lmlocal,lllocal,lnlocal,kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf, &
                            lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xkm,ykm,zkm,xln,yln,zln,d2trm, &
                            dr2ds6t(1,3),dr2ds6t2(1,3),lilocal)
!
!  j-l / j-k
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(1,nt2)
          call add_drv2_2dm(ljlocal,lllocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjl,yjl,zjl,xjk,yjk,zjk,d2trm, &
                            dr2ds6t(1,4),dr2ds6t(1,1),lilocal)
!
!  j-l / k-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(2,nt2)
          call add_drv2_2dm(ljlocal,lllocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xkl,ykl,zkl,d2trm, &
                            dr2ds6t(1,4),dr2ds6t(1,2),lilocal)
!
!  j-l / k-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(3,nt2)
          call add_drv2_2dm(ljlocal,lllocal,lklocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xjl,yjl,zjl,xkn,ykn,zkn,d2trm, &
                            dr2ds6t(1,4),dr2ds6t2(1,1),lilocal)
!
!  j-l / j-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(4,nt2)
          call add_drv2_2dm(ljlocal,lllocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjl,yjl,zjl,xjl,yjl,zjl,d2trm, &
                            dr2ds6t(1,4),dr2ds6t(1,4),lilocal)
!
!  j-l / j-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(5,nt2)
          call add_drv2_2dm(ljlocal,lllocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjl,yjl,zjl,xjn,yjn,zjn,d2trm, &
                            dr2ds6t(1,4),dr2ds6t2(1,2),lilocal)
!
!  j-l / l-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(4,nt)*dft(6,nt2)
          call add_drv2_2dm(ljlocal,lllocal,lllocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf, &
                            lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xjl,yjl,zjl,xln,yln,zln,d2trm, &
                            dr2ds6t(1,4),dr2ds6t2(1,3),lilocal)
!
!  j-m / j-k
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(1,nt2)
          call add_drv2_2dm(ljlocal,lmlocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjm,yjm,zjm,xjk,yjk,zjk,d2trm, &
                            dr2ds6t(1,5),dr2ds6t(1,1),lilocal)
!
!  j-m / k-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(2,nt2)
          call add_drv2_2dm(ljlocal,lmlocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjm,yjm,zjm,xkl,ykl,zkl,d2trm, &
                            dr2ds6t(1,5),dr2ds6t(1,2),lilocal)
!
!  j-m / k-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(3,nt2)
          call add_drv2_2dm(ljlocal,lmlocal,lklocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xkn,ykn,zkn,d2trm, &
                            dr2ds6t(1,5),dr2ds6t2(1,1),lilocal)
!
!  j-m / j-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(4,nt2)
          call add_drv2_2dm(ljlocal,lmlocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjm,yjm,zjm,xjl,yjl,zjl,d2trm, &
                            dr2ds6t(1,5),dr2ds6t(1,4),lilocal)
!
!  j-m / j-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(5,nt2)
          call add_drv2_2dm(ljlocal,lmlocal,ljlocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xjn,yjn,zjn,d2trm, &
                            dr2ds6t(1,5),dr2ds6t2(1,2),lilocal)
!
!  j-m / l-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(5,nt)*dft(6,nt2)
          call add_drv2_2dm(ljlocal,lmlocal,lllocal,lnlocal,jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf, &
                            lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xjm,yjm,zjm,xln,yln,zln,d2trm, &
                            dr2ds6t(1,5),dr2ds6t2(1,3),lilocal)
!
!  l-m / j-k
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(1,nt2)
          call add_drv2_2dm(lllocal,lmlocal,ljlocal,lklocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xlm,ylm,zlm,xjk,yjk,zjk,d2trm, &
                            dr2ds6t(1,6),dr2ds6t(1,1),lilocal)
!
!  l-m / k-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(2,nt2)
          call add_drv2_2dm(lllocal,lmlocal,lklocal,lllocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                            kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xlm,ylm,zlm,xkl,ykl,zkl,d2trm, &
                            dr2ds6t(1,6),dr2ds6t(1,2),lilocal)
!
!  l-m / k-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(3,nt2)
          call add_drv2_2dm(lllocal,lmlocal,lklocal,lnlocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                            kx,ky,kz,kxf,kyf,kzf,nx,ny,nz,nxf,nyf,nzf,xlm,ylm,zlm,xkn,ykn,zkn,d2trm, &
                            dr2ds6t(1,6),dr2ds6t2(1,1),lilocal)
!
!  l-m / j-l
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(4,nt2)
          call add_drv2_2dm(lllocal,lmlocal,ljlocal,lllocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xlm,ylm,zlm,xjl,yjl,zjl,d2trm, &
                            dr2ds6t(1,6),dr2ds6t(1,4),lilocal)
!
!  l-m / j-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(5,nt2)
          call add_drv2_2dm(lllocal,lmlocal,ljlocal,lnlocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                            jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf,xlm,ylm,zlm,xjn,yjn,zjn,d2trm, &
                            dr2ds6t(1,6),dr2ds6t2(1,2),lilocal)
!
!  l-m / l-n
!
          d2trm = - const*fangle*rdamp*qhdampo*ftorsm2*dft(6,nt)*dft(6,nt2)
          call add_drv2_2dm(lllocal,lmlocal,lllocal,lnlocal,lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf, &
                            lx,ly,lz,lxf,lyf,lzf,nx,ny,nz,nxf,nyf,nzf,xlm,ylm,zlm,xln,yln,zln,d2trm, &
                            dr2ds6t(1,6),dr2ds6t2(1,3),lilocal)
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
          call real1strterm(ndim,xij,yij,zij,0.0_dp,0.0_dp,0.0_dp,dr2ds3n(1,1),d2r2dx23n(1,1),d2r2dsdx3n(1,1,1), &
                            d2r2ds23n(1,1,1),.false.)
          call real1strterm(ndim,xin,yin,zin,0.0_dp,0.0_dp,0.0_dp,dr2ds3n(1,2),d2r2dx23n(1,2),d2r2dsdx3n(1,1,2), &
                            d2r2ds23n(1,1,2),.false.)
          call real1strterm(ndim,xjn,yjn,zjn,0.0_dp,0.0_dp,0.0_dp,dr2ds3n(1,3),d2r2dx23n(1,3),d2r2dsdx3n(1,1,3), &
                            d2r2ds23n(1,1,3),.false.)
!
!  i-j / j-k
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(1,nt)
          call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xij,yij,zij,xjk,yjk,zjk,d2trm, &
                            dr2ds3n(1,1),dr2ds6t(1,1),lilocal)
!
!  i-j / k-l
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(2,nt)
          call add_drv2_2dm(lilocal,ljlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xkl,ykl,zkl,d2trm, &
                            dr2ds3n(1,1),dr2ds6t(1,2),lilocal)
!
!  i-j / k-m
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(3,nt)
          call add_drv2_2dm(lilocal,ljlocal,lklocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xkm,ykm,zkm,d2trm, &
                            dr2ds3n(1,1),dr2ds6t(1,3),lilocal)
!
!  i-j / j-l
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(4,nt)
          call add_drv2_2dm(lilocal,ljlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xij,yij,zij,xjl,yjl,zjl,d2trm, &
                            dr2ds3n(1,1),dr2ds6t(1,4),lilocal)
!
!  i-j / j-m
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(5,nt)
          call add_drv2_2dm(lilocal,ljlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xjm,yjm,zjm,d2trm, &
                            dr2ds3n(1,1),dr2ds6t(1,5),lilocal)
!
!  i-j / l-m
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(1,nj)*dft(6,nt)
          call add_drv2_2dm(lilocal,ljlocal,lllocal,lmlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                            lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xij,yij,zij,xlm,ylm,zlm,d2trm, &
                            dr2ds3n(1,1),dr2ds6t(1,6),lilocal)
!
!  i-n / j-k
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(1,nt)
          call add_drv2_2dm(lilocal,lnlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xin,yin,zin,xjk,yjk,zjk,d2trm, &
                            dr2ds3n(1,2),dr2ds6t(1,1),lilocal)
!
!  i-n / k-l
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(2,nt)
          call add_drv2_2dm(lilocal,lnlocal,lklocal,lllocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                            kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xin,yin,zin,xkl,ykl,zkl,d2trm, &
                            dr2ds3n(1,2),dr2ds6t(1,2),lilocal)
!
!  i-n / k-m
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(3,nt)
          call add_drv2_2dm(lilocal,lnlocal,lklocal,lmlocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                            kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xin,yin,zin,xkm,ykm,zkm,d2trm, &
                            dr2ds3n(1,2),dr2ds6t(1,3),lilocal)
!
!  i-n / j-l
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(4,nt)
          call add_drv2_2dm(lilocal,lnlocal,ljlocal,lllocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xin,yin,zin,xjl,yjl,zjl,d2trm, &
                            dr2ds3n(1,2),dr2ds6t(1,4),lilocal)
!
!  i-n / j-m
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(5,nt)
          call add_drv2_2dm(lilocal,lnlocal,ljlocal,lmlocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xin,yin,zin,xjm,yjm,zjm,d2trm, &
                            dr2ds3n(1,2),dr2ds6t(1,5),lilocal)
!
!  i-n / l-m
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(2,nj)*dft(6,nt)
          call add_drv2_2dm(lilocal,lnlocal,lllocal,lmlocal,ix,iy,iz,ixf,iyf,izf,nx,ny,nz,nxf,nyf,nzf, &
                            lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xin,yin,zin,xlm,ylm,zlm,d2trm, &
                            dr2ds3n(1,2),dr2ds6t(1,6),lilocal)
!
!  j-n / j-k
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(1,nt)
          call add_drv2_2dm(ljlocal,lnlocal,ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                            jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xjn,yjn,zjn,xjk,yjk,zjk,d2trm, &
                            dr2ds3n(1,3),dr2ds6t(1,1),lilocal)
!
!  j-n / k-l
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(2,nt)
          call add_drv2_2dm(ljlocal,lnlocal,lklocal,lllocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                            kx,ky,kz,kxf,kyf,kzf,lx,ly,lz,lxf,lyf,lzf,xjn,yjn,zjn,xkl,ykl,zkl,d2trm, &
                            dr2ds3n(1,3),dr2ds6t(1,2),lilocal)
!       
!  j-n / k-m
!   
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(3,nt)
          call add_drv2_2dm(ljlocal,lnlocal,lklocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                            kx,ky,kz,kxf,kyf,kzf,mx,my,mz,mxf,myf,mzf,xjn,yjn,zjn,xkm,ykm,zkm,d2trm, &
                            dr2ds3n(1,3),dr2ds6t(1,3),lilocal)
!
!  j-n / j-l
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(4,nt)
          call add_drv2_2dm(ljlocal,lnlocal,ljlocal,lllocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                            jx,jy,jz,jxf,jyf,jzf,lx,ly,lz,lxf,lyf,lzf,xjn,yjn,zjn,xjl,yjl,zjl,d2trm, &
                            dr2ds3n(1,3),dr2ds6t(1,4),lilocal)
!
!  j-n / j-m
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(5,nt)
          call add_drv2_2dm(ljlocal,lnlocal,ljlocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                            jx,jy,jz,jxf,jyf,jzf,mx,my,mz,mxf,myf,mzf,xjn,yjn,zjn,xjm,yjm,zjm,d2trm, &
                            dr2ds3n(1,3),dr2ds6t(1,5),lilocal)
!
!  j-n / l-m
!
          d2trm = - const*fangle*rdamp*dampo*dampo_nb_totm1*ftorsm1*ddampo_nbdr(3,nj)*dft(6,nt)
          call add_drv2_2dm(ljlocal,lnlocal,lllocal,lmlocal,jx,jy,jz,jxf,jyf,jzf,nx,ny,nz,nxf,nyf,nzf, &
                            lx,ly,lz,lxf,lyf,lzf,mx,my,mz,mxf,myf,mzf,xjn,yjn,zjn,xlm,ylm,zlm,d2trm, &
                            dr2ds3n(1,3),dr2ds6t(1,6),lilocal)
        enddo
      endif
    enddo
  endif
!
  if (lgrad1) then
    if (lgrad2) then
      deallocate(d2ft,stat=status)
      if (status/=0) call deallocate_error('gfnff_hb2_eg3d','d2ft')
    endif
    deallocate(dft,stat=status)
    if (status/=0) call deallocate_error('gfnff_hb2_eg3d','dft')
  endif
  deallocate(ft,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg3d','ft')
!
!  Exit point to ensure memory is deallocated apart from ft
!
  1000 continue
!
  if (lgrad1) then
    if (lgrad2) then
      deallocate(d2dampo_nbdr2,stat=status)
      if (status/=0) call deallocate_error('gfnff_hb2_eg3d','d2dampo_nbdr2')
    endif
    deallocate(ddampo_nbdr,stat=status)
    if (status/=0) call deallocate_error('gfnff_hb2_eg3d','ddampo_nbdr')
  endif
  deallocate(dampo_nb,stat=status)
  if (status/=0) call deallocate_error('gfnff_hb2_eg3d','dampo_nb')
!
#ifdef TRACE
  call trace_out('gfnff_hb2_eg3d')
#endif
!
  return
  end
