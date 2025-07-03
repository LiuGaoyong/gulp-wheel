  subroutine setreaxffQ(nboatom,nboatomRptr,nbos,nbosptr,qreaxFF,eself,loutputQ,ldoQdrv)
!
!  Subroutine for performing electronegativity equalisation calcns
!  in order to determine reaxFF charges.
!
!   1/08 Created based on eem & reaxff
!   4/08 Qtot now properly initialised
!   4/08 qreaxFF returned to full atom array at the end
!   6/08 Shell structure correction added along with iterative solution
!   6/08 Option to fix charges in ReaxFF added
!  10/09 qr12 term added
!  11/09 qr12 term removed
!  12/09 Solution for charges modified to use lapack solution rather 
!        than matrix inversion.
!   4/10 Setting of nfi corrected
!   9/11 Use of ni/nj for nat of i and j replaced by nati/natj for
!        consistency with setreaxffQiter.
!  10/11 Site energy terms added
!  12/11 nbosptr added to arguments of setreaxff
!  12/11 Dimension of nbos corrected to numat from nboatom
!  12/11 lreaxFFqfix and reaxFFgamma referenced by species
!   8/12 Modified to trap segmentation fault for the case where there
!        are no atoms in a structure!
!   8/13 Declaration of qreaxFF changed to inout
!   7/15 External potential added
!   8/15 Argument added to flag whether charges should be output
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   9/19 ReaxFF charge parameters now given by species
!   3/20 Location of angstoev changed to current
!  12/20 Atom number format increased to allow for larger values
!  12/22 Calculation of charge derivatives added
!   3/23 Correction to dqds
!   4/23 Sign of strain term corrected
!   6/23 Electrostatic cutoff squared now passed to rsearch subroutines
!   7/23 extpot added
!   8/23 Electric fields added
!   8/23 Transfer of charges to regular arrays added
!   9/23 dqdVxyz added
!  10/23 Electric field charge derivatives corrected for symmetry-adapted case
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
!  Julian Gale, CIC, Curtin University, October 2023
!
  use control
  use current
  use derivatives
  use element
  use energies,         only : siteenergy
  use field,            only : lfieldcfg, ntdfieldcfg
  use iochannels
  use m_strain,         only : twostrterms
  use numbers,          only : third
  use parallel
  use realvectors
  use reaxFFdata,       only : reaxFFcutoffQ, reaxFFqfix, reaxFFchi, reaxFFmu, reaxFFgamma
  use reaxFFdata,       only : reaxFFqdamp, nreaxFFqiter, reaxFFqconverged, reaxFFqconverged1
  use reaxFFdata,       only : reaxFFshell, lreaxFFqfix
  use times
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                      :: nboatom                ! Number of bond order atoms
  integer(i4), intent(in)                      :: nboatomRptr(nboatom)   ! Pointer to bond order atoms
  integer(i4), intent(in)                      :: nbos(numat)            ! Number of ReaxFF species for atoms
  integer(i4), intent(in)                      :: nbosptr(numat)         ! Pointer from atom to species
  logical,     intent(in)                      :: loutputQ               ! If true then output charges 
  logical,     intent(in)                      :: ldoQdrv                ! If true then compute charge derivatives
  real(dp),    intent(inout)                   :: qreaxFF(numat+1)       ! At end, contains array of charges
  real(dp),    intent(out)                     :: eself                  ! At end, contains self energy
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ifail
  integer(i4)                                  :: igroup
  integer(i4)                                  :: ii
  integer(i4)                                  :: imin
  integer(i4)                                  :: ind
  integer(i4)                                  :: indi
  integer(i4)                                  :: indj
  integer(i4), dimension(:), allocatable       :: ipivot
  integer(i4)                                  :: iresid
  integer(i4)                                  :: is
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: jmax
  integer(i4)                                  :: k
  integer(i4)                                  :: n
  integer(i4)                                  :: nd
  integer(i4)                                  :: nfree
  integer(i4), dimension(:), allocatable       :: nfreeptr
  integer(i4)                                  :: nfi
  integer(i4)                                  :: nfj
  integer(i4)                                  :: niter
  integer(i4)                                  :: nitermax
  integer(i4)                                  :: nmolonly
  integer(i4)                                  :: nor
  integer(i4)                                  :: nr
  integer(i4)                                  :: ns
  integer(i4)                                  :: nspeci
  integer(i4)                                  :: nspecj
  integer(i4)                                  :: status
  logical                                      :: lfixQi
  logical                                      :: lfixQj
  logical                                      :: literative
  logical                                      :: lmixed
  logical                                      :: lself
  logical                                      :: lstrain
  real(dp),    dimension(:), allocatable       :: dpacked
  real(dp),    dimension(:), allocatable       :: dpackedsave
  real(dp)                                     :: cut2
  real(dp)                                     :: damp
  real(dp),    dimension(:,:), allocatable     :: dqs
  real(dp),    dimension(:),   allocatable     :: dqx
  real(dp),    dimension(:),   allocatable     :: dqy
  real(dp),    dimension(:),   allocatable     :: dqz
  real(dp)                                     :: esite
  real(dp)                                     :: dgamdrij
  real(dp)                                     :: fieldx
  real(dp)                                     :: fieldy
  real(dp)                                     :: fieldz
  real(dp)                                     :: gam
  real(dp)                                     :: gam1
  real(dp)                                     :: gammai
  real(dp)                                     :: gammaj
  real(dp)                                     :: gammaij
  real(dp)                                     :: gamthird
  real(dp)                                     :: qdiff
  real(dp)                                     :: qdiff1
  real(dp)                                     :: qi
  real(dp)                                     :: qj
  real(dp)                                     :: qsum
  real(dp)                                     :: qtot
  real(dp),    dimension(:), allocatable       :: qreaxFFfree
  real(dp),    dimension(:), allocatable       :: qreaxFFsave
  real(dp)                                     :: rij
  real(dp)                                     :: rij2
  real(dp)                                     :: tp
  real(dp)                                     :: dtpdr
  real(dp)                                     :: d2tpdr2
  real(dp)                                     :: xji
  real(dp)                                     :: yji
  real(dp)                                     :: zji
  real(dp),    dimension(:), allocatable       :: vfield
  real(dp),    dimension(:), allocatable       :: z
!
!  Set local strain flag
!
  lstrain = (ndim.gt.0)
!
  if (ldoQdrv) then
!
!  Check the memory for the charge derivatives if required
!
    if (numat.gt.maxd2qu) then
      maxd2qu = numat
      call changemaxd2q
    endif
    if (3*numat.gt.maxd2q) then
      maxd2q = 3*numat
      call changemaxd2q
    endif
!
!  Zero dq/dxyz and dq/ds
!
    do i = 1,numat
      do j = 1,3*numat
        dqdxyz(j,i) = 0.0_dp
      enddo
      dqdVxyz(1:3,i) = 0.0_dp
    enddo
    if (ndim.gt.0) then
      do i = 1,numat
        do j = 1,nstrains
          dqds(j,i) = 0.0_dp
        enddo
      enddo
    endif
  endif
!
!  Allocate pointer array for atoms to act on
!
  allocate(nfreeptr(nboatom),stat=status)
  if (status/=0) call outofmemory('setreaxffQ','nfreeptr')
!
!  Find out whether an iterative solution is needed due to shell structure or whether there are any fixed charges
!
  literative = .false.
  qsum = 0.0_dp
  nfree = 0
  do i = 1,nboatom
    ii = nboatomRptr(i)
    if (ii.gt.0) then
      nspeci = nbosptr(ii)
      if (abs(reaxFFshell(1,nspeci)).gt.1.0d-12) then
        literative = .true.
      endif
      if (lreaxFFqfix(nspeci)) then
        qsum = qsum - reaxFFqfix(nspeci)*occuf(ii)
      else
        nfree = nfree + 1
        nfreeptr(nfree) = i
      endif
    endif
  enddo
!
!  Allocate local memory for charge derivatives
!
  if (ldoQdrv) then
    allocate(dqx(nfree+1),stat=status)
    if (status/=0) call outofmemory('setreaxffQ','dqx')
    allocate(dqy(nfree+1),stat=status)
    if (status/=0) call outofmemory('setreaxffQ','dqy')
    allocate(dqz(nfree+1),stat=status)
    if (status/=0) call outofmemory('setreaxffQ','dqz')
    if (lstrain) then
      allocate(dqs(nfree+1,nstrains),stat=status)
      if (status/=0) call outofmemory('setreaxffQ','dqs')
    endif
  endif
!
!  Set electric field if required
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    allocate(vfield(numat),stat=status)
    if (status/=0) call outofmemory('setreaxffQ','vfield') 
    vfield(1:numat) = 0.0_dp
    call electricfieldpotl(vfield)
    call electricfieldparts(fieldx,fieldy,fieldz)
  endif
!
!  Decide on total EEM fragment charge - for cluster
!  there is no constraint so use sum of initial charges.
!  For periodic system, charge must be equal to the
!  negative sum of the non-EEM ion charges.
!
  qtot = totalcharge + qsum
!
!  Allocate local memory that depends on nboatom
!
  n = nfree + 1
  allocate(z(n),stat=status)
  if (status/=0) call outofmemory('setreaxffQ','z')
  allocate(qreaxFFfree(nfree+1),stat=status)
  if (status/=0) call outofmemory('setreaxffQ','qreaxFFfree')
  allocate(qreaxFFsave(nboatom+1),stat=status)
  if (status/=0) call outofmemory('setreaxffQ','qreaxFFsave')
  allocate(dpacked(n*(n+1)/2),stat=status)
  if (status/=0) call outofmemory('setreaxffQ','dpacked')
  if (literative) then
    allocate(dpackedsave(n*(n+1)/2),stat=status)
    if (status/=0) call outofmemory('setreaxffQ','dpackedsave')
  endif
!
!  Form right hand vector
!
  z(1:n) = 0.0_dp
  do i = 1,nfree
    ii = nboatomRptr(nfreeptr(i))
    nspeci = nbosptr(ii)
    z(i) = z(i) - reaxFFchi(nspeci) - extpot(nrelf2a(ii))
  enddo
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    do i = 1,nfree
      ii = nboatomRptr(nfreeptr(i))
      z(i) = z(i) + vfield(ii)   ! NB Sign is opposite to extpot and chi
    enddo
  endif
!
  z(nfree+1) = qtot
!*******************************************
!  Generate electrostatic potential terms  *
!*******************************************
!
!  Initialise i-j matrix elements
!
  dpacked(1:n*(n+1)/2) = 0.0_dp
!
!  Set cutoffs
!
  cut2 = reaxFFcutoffQ**2
!
!  Set lower bound for i loop
!
  if (ndim.gt.0) then
    imin = 1
    ind = 0
    nfi = 0
  else
    imin = 2
    if (nboatom.gt.0) then
      if (nboatomRptr(1).gt.0) then
        if (lreaxFFqfix(nbosptr(nboatomRptr(1)))) then
          ind = 0
        else
          ind = 1
        endif
      else
        ind = 0
      endif
!
!  Set nfi allowing for whether skipped first atom is fixed or not
!
      lfixQi = lreaxFFqfix(nbosptr(nboatomRptr(1)))
      if (lfixQi) then
        nfi = 0
      else
        nfi = 1
      endif
    else
      ind = 0
      nfi = 1
    endif
  endif
!
!  Loop over first atom
!
  do ii = imin,nboatom
    i = nboatomRptr(ii)
    if (i.gt.0) then
      nspeci = nbosptr(i)
      lfixQi = lreaxFFqfix(nspeci)
      if (.not.lfixQi) nfi = nfi + 1
!
!  Does i have a reaxFF species?
!
      if (nbos(i).gt.0.or.lfixQi) then
        gammai = reaxFFgamma(nspeci)
!
!  Set upper bound for j loop
!
        if (ndim.gt.0) then
          jmax = ii
        else
          jmax = ii - 1
        endif
!
!  Loop over second atom
!
        nfj = 0
        jloop: do jj = 1,jmax
          j = nboatomRptr(jj)
          if (j.gt.0) then
            nspecj = nbosptr(j)
            lfixQj = lreaxFFqfix(nspecj)
            if (.not.lfixQj) nfj = nfj + 1
!
!  Does j have a reaxFF species?
!
            if (nbos(j).gt.0.or.lfixQj) then
!
!  Skip loop for pairs of fixed atoms
!
              if (lfixQi.and.lfixQj) then
                cycle jloop
              elseif (.not.lfixQi.and..not.lfixQj) then
!
!  Only increment matrix element for pairs of free atoms
!
                ind = ind + 1
                lmixed = .false.
              else
                lmixed = .true.
              endif
!
              gammaj = reaxFFgamma(nspecj)
!
!  Compute basic interatomic vector
!
              xji = xclat(j) - xclat(i)
              yji = yclat(j) - yclat(i)
              zji = zclat(j) - zclat(i)
!
!  Find valid vectors
!
              nor = 0
              if (ndim.eq.3) then
                call rsearch3D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,cut2)
              elseif (ndim.eq.2) then
                call rsearch2D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,cut2)
              elseif (ndim.eq.1) then
                call rsearch1D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,cut2)
              elseif (ndim.eq.0) then
                rij2 = xji*xji + yji*yji + zji*zji
                if (rij2.lt.cut2.and.rij2.gt.1.0d-12) then
                  nor = nor + 1
                endif
              endif
!
!  If no distances then cycle
!
              if (nor.eq.0) cycle jloop
!
!  Compute Coulomb shielding parameters 
!
              gammaij = sqrt(gammai*gammaj)
              gammaij = 1.0_dp/(gammaij**3)
!
!  Loop over valid distances and calculate contributions
!
              if (ndim.gt.0) then
                do nd = 1,nor
                  rij2 = dist(nd)
                  rij = sqrt(rij2)
!
!  Compute taper function
!
                  call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.false.,.false.)
!                     
!  Pure real space tapered lreaxFF case
!                     
                  gam = tp/(rij2*rij + gammaij)**third
                  if (lmixed) then
                    if (lfixQi) then
                      z(nfj) = z(nfj) - gam*reaxFFqfix(nspeci)*angstoev
                    else
                      z(nfi) = z(nfi) - gam*reaxFFqfix(nspecj)*angstoev
                    endif
                  else
                    dpacked(ind) = dpacked(ind) + gam
                  endif
                enddo
              else
                rij = sqrt(rij2)
!
!  Compute taper function
!
                call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.false.,.false.)
!                     
!  Pure real space tapered lreaxFF case
!                     
                gam = tp/(rij2*rij + gammaij)**third
                if (lmixed) then
                  if (lfixQi) then
                    z(nfj) = z(nfj) - gam*reaxFFqfix(nspeci)*angstoev
                  else
                    z(nfi) = z(nfi) - gam*reaxFFqfix(nspecj)*angstoev
                  endif
                else
                  dpacked(ind) = dpacked(ind) + gam
                endif
              endif
            endif
          endif
!
!  End of loop over j
!
        enddo jloop
!
!  For 0-D case we need to increment ind by 1 to allow for diagonal element that is skipped
!
        if (ndim.eq.0.and..not.lfixQi) then
          ind = ind + 1
        endif
      endif
    endif
!
!  End of loop over i
!
  enddo
!
!  Scale matrix elements by conversion factor
!
  n = nfree + 1
  do i = 1,n*(n+1)/2
    dpacked(i) = dpacked(i)*angstoev
  enddo
!     
!  Allocate workspace for inversion
!     
  allocate(ipivot(n),stat=status)
  if (status/=0) call outofmemory('setreaxffQ','ipivot')
!
!  If this an iterative run then save the original dpacked matrix to avoid recomputing it
!
  if (literative) then
    do i = 1,n*(n+1)/2
      dpackedsave(i) = dpacked(i)
    enddo
    nitermax = nreaxFFqiter
  else
    nitermax = 1
  endif
!****************************
!  Start of iterative loop  *
!****************************
  do niter = 1,nitermax
!
!  If this is not the first time then copy saved dpacked back
!
    if (niter.gt.1) then
      n = nfree + 1
      do i = 1,n*(n+1)/2
        dpacked(i) = dpackedsave(i)
      enddo
    endif
!
!  Save qreaxFF for later comparison
!
    qreaxFFsave(1:nboatom) = qreaxFF(1:nboatom) 
!********************************
!  Form matrix of coefficients  *
!********************************
    do i = 1,nfree
      ii = nboatomRptr(nfreeptr(i))
      ind = i*(i+1)/2
      nspeci = nbosptr(ii)
      dpacked(ind) = dpacked(ind) + 2.0_dp*reaxFFmu(nspeci)*occuf(ii)
!
!  Shell structure option
!
      if (abs(reaxFFshell(1,nspeci)).gt.1.0d-12) then
        if (qreaxFF(ii).gt.reaxFFshell(2,nspeci)) then
          dpacked(ind) = dpacked(ind) + 4.0_dp*reaxFFshell(1,nspeci)*occuf(ii)*(qreaxFF(ii) - reaxFFshell(2,nspeci))**3
        endif
      endif
    enddo
    ind = nfree*(nfree + 1)/2 
    do i = 1,nfree
      dpacked(ind+i) = 1.0_dp
    enddo
!
!  Debugging output
!
    if (index(keyword,'debu').ne.0.and.ioproc) then
      write(ioout,'(/,''  ReaxFF Charge Matrix :'',/)')
      ind = 0
      do i = 1,nfree + 1
        write(ioout,'(1x,f9.5,4x,9(1x,f9.5))') z(i),(dpacked(j),j=ind+1,ind+i)
        ind = ind + i
      enddo
    endif
    if (nfree.gt.0) then
!******************
!  Invert matrix  *
!******************
      ifail = 0
!     
!  Factorise matrix
!     
      n = nfree + 1
      call dsptrf('U',n,dpacked,ipivot,ifail)
      if (ifail.eq.0) then
!
!  Copy z to qreaxFFfree
!
        qreaxFFfree(1:n) = z(1:n)
!     
!  Solve for new charges
!     
        call dsptrs('U',n,1_i4,dpacked,ipivot,qreaxFFfree,n,ifail)
      endif
!
!  Was inversion successful?
!
      if (ifail.ne.0) then
        call outerror('charge solution failed in setreaxffQ',0_i4)
        call stopnow('setreaxffQ')
      endif
    endif
!
!  Expand back to full array
!
    nfi = 0
    do i = 1,nboatom
      ii = nboatomRptr(i)
      if (ii.gt.0) then
        nspeci = nbosptr(ii)
        if (lreaxFFqfix(nspeci)) then
          qreaxFF(i) = reaxFFqfix(nspeci)
        else
          nfi = nfi + 1
          qreaxFF(i) = qreaxFFfree(nfi)
        endif
      endif
    enddo
!
    if (literative) then
!
!  Check for convergence
!
      qdiff  = 0.0_dp
      qdiff1 = 0.0_dp
      do i = 1,nboatom
        qdiff  = qdiff + abs(qreaxFF(i) - qreaxFFsave(i))
        qdiff1 = max(qdiff1,abs(qreaxFF(i) - qreaxFFsave(i)))
      enddo
      qdiff = qdiff/dble(nboatom)
      if (index(keyword,'debu').ne.0.and.ioproc) then
        write(ioout,'(''  ** Cycle : '',i4,'' Qdiffs : '',2f10.8)') niter,qdiff,qdiff1
      endif
      if (qdiff.lt.reaxFFqconverged.and.qdiff1.lt.reaxFFqconverged1) exit 
!
!  If not converged then damp charges for next iteration except for first iteration
!
      if (niter.gt.1) then
        damp = reaxFFqdamp
        do i = 1,nboatom
          qreaxFF(i) = damp*qreaxFF(i) + (1.0_dp - damp)*qreaxFFsave(i)
        enddo
      endif
      if (index(keyword,'debu').ne.0.and.ioproc) then
        write(ioout,'(//,''  Charges for ReaxFF during iteration :'',/)')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
        write(ioout,'(''--------------------------------------------------------------------------------'')')
        do i = 1,nboatom
          ii = nboatomRptr(i)
          if (ii.gt.0) then
            write(ioout,'(6x,i4,18x,i2,16x,f10.7)') i,nat(ii),qreaxFF(i)
          endif
        enddo
        write(ioout,'(''--------------------------------------------------------------------------------'',/)')
      endif
    endif
!**************************
!  End of iterative loop  *
!**************************
  enddo
!***********************
!  Charge derivatives  *
!***********************
  if (ldoQdrv) then
!
!  Build matrix derivatives
!
!  Loop over first atom
!
    nfi = 0
    do ii = 1,nboatom
      i = nboatomRptr(ii)
      indi = 3*(i-1)
      if (i.gt.0) then
        nspeci = nbosptr(i)
        lfixQi = lreaxFFqfix(nspeci)
        qi = qreaxFF(i)
        if (.not.lfixQi) nfi = nfi + 1
!
!  Does i have a reaxFF species?
!
        if (nbos(i).gt.0.or.lfixQi) then
          gammai = reaxFFgamma(nspeci)
!
!  Initialise charge derivatives
!
          dqx(1:n) = 0.0_dp
          dqy(1:n) = 0.0_dp
          dqz(1:n) = 0.0_dp
          if (lstrain) then
            dqs(1:n,1:nstrains) = 0.0_dp
          endif
!
!  Set upper bound for j loop
!
          if (ndim.gt.0) then
            jmax = ii
          else
            jmax = ii - 1
          endif
!
!  Loop over second atom
!
          nfj = 0
          jqloop: do jj = 1,nboatom
            j = nboatomRptr(jj)
            if (j.gt.0) then
              nspecj = nbosptr(j)
              lfixQj = lreaxFFqfix(nspecj)
              qj = qreaxFF(j)
              if (.not.lfixQj) nfj = nfj + 1
!
!  Does j have a reaxFF species?
!
              if (nbos(j).gt.0.or.lfixQj) then
!
!  Skip loop for pairs of fixed atoms
!
                if (lfixQi.and.lfixQj) then
                  cycle jqloop
                elseif (.not.lfixQi.and..not.lfixQj) then
!
!  Only increment matrix element for pairs of free atoms
!
                  ind = ind + 1
                  lmixed = .false.
                else
                  lmixed = .true.
                endif
!
                gammaj = reaxFFgamma(nspecj)
!
!  Compute basic interatomic vector
!
                xji = xclat(j) - xclat(i)
                yji = yclat(j) - yclat(i)
                zji = zclat(j) - zclat(i)
!
!  Find valid vectors
!
                nor = 0
                if (ndim.eq.3) then
                  call rsearch3D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,cut2)
                elseif (ndim.eq.2) then
                  call rsearch2D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,cut2)
                elseif (ndim.eq.1) then
                  call rsearch1D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,cut2)
                elseif (ndim.eq.0) then
                  rij2 = xji*xji + yji*yji + zji*zji
                  if (rij2.lt.cut2.and.rij2.gt.1.0d-12) then
                    nor = nor + 1
                  endif
                endif
!
!  If no distances then cycle
!
                if (nor.eq.0) cycle jqloop
!
!  Generate strain products
!
                if (lstrain) then
                  call twostrterms(ndim,maxdis,nor,xtmp,ytmp,ztmp,0.0_dp,0.0_dp,0.0_dp,dr2ds,rpd,d2r2dsdx,d2r2ds2,.false.)
                endif
!
!  Compute Coulomb shielding parameters 
!
                gammaij = sqrt(gammai*gammaj)
                gammaij = 1.0_dp/(gammaij**3)
!
!  Loop over valid distances and calculate contributions
!
                if (ndim.gt.0) then
                  do nd = 1,nor
                    rij2 = dist(nd)
                    rij = sqrt(rij2)
                    xji = xtmp(nd)
                    yji = ytmp(nd)
                    zji = ztmp(nd)
!
!  Compute taper function
!
                    call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.true.,.false.)
!                     
!  Pure real space tapered lreaxFF case
!                     
                    gam1 = 1.0_dp/(rij2*rij + gammaij)
                    gamthird = gam1**third
                    dgamdrij = gamthird*(dtpdr - tp*rij*gam1)
!
                    if (lmixed) then
                      if (lfixQi) then
                        dqx(nfj) = dqx(nfj) - dgamdrij*xji*qi
                        dqy(nfj) = dqy(nfj) - dgamdrij*yji*qi
                        dqz(nfj) = dqz(nfj) - dgamdrij*zji*qi
                      else
                        dqx(nfi) = dqx(nfi) - dgamdrij*xji*qj
                        dqy(nfi) = dqy(nfi) - dgamdrij*yji*qj
                        dqz(nfi) = dqz(nfi) - dgamdrij*zji*qj
                      endif
                    else
                      dqx(nfi) = dqx(nfi) + dgamdrij*xji*qj
                      dqy(nfi) = dqy(nfi) + dgamdrij*yji*qj
                      dqz(nfi) = dqz(nfi) + dgamdrij*zji*qj
                      dqx(nfj) = dqx(nfj) + dgamdrij*xji*qi
                      dqy(nfj) = dqy(nfj) + dgamdrij*yji*qi
                      dqz(nfj) = dqz(nfj) + dgamdrij*zji*qi
                      if (lstrain) then
                        do is = 1,nstrains
                          ns = nstrptr(is)
                          dqs(nfi,is) = dqs(nfi,is) + qj*dgamdrij*dr2ds(nd,ns)
                        enddo
                      endif
                    endif
                  enddo
                else
                  rij = sqrt(rij2)
!
!  Compute taper function
!
                  call p7reaxFFqtaper(rij,reaxFFcutoffQ,tp,dtpdr,d2tpdr2,.true.,.false.)
!                     
!  Pure real space tapered lreaxFF case
!                     
                  gam1 = 1.0_dp/(rij2*rij + gammaij)
                  gamthird = gam1**third
                  dgamdrij = gamthird*(dtpdr - tp*rij*gam1)
!
                  if (lmixed) then
                    if (lfixQi) then
                      dqx(nfj) = dqx(nfj) - dgamdrij*xji*qi
                      dqy(nfj) = dqy(nfj) - dgamdrij*yji*qi
                      dqz(nfj) = dqz(nfj) - dgamdrij*zji*qi
                    else
                      dqx(nfi) = dqx(nfi) - dgamdrij*xji*qj
                      dqy(nfi) = dqy(nfi) - dgamdrij*yji*qj
                      dqz(nfi) = dqz(nfi) - dgamdrij*zji*qj
                    endif
                  else
                    dqx(nfi) = dqx(nfi) + dgamdrij*xji*qj
                    dqy(nfi) = dqy(nfi) + dgamdrij*yji*qj
                    dqz(nfi) = dqz(nfi) + dgamdrij*zji*qj
                    dqx(nfj) = dqx(nfj) + dgamdrij*xji*qi
                    dqy(nfj) = dqy(nfj) + dgamdrij*yji*qi
                    dqz(nfj) = dqz(nfj) + dgamdrij*zji*qi
                    if (lstrain) then
                      do is = 1,nstrains
                        ns = nstrptr(is)
                        dqs(nfi,is) = dqs(nfi,is) + qj*dgamdrij*dr2ds(nd,ns)
                      enddo
                    endif
                  endif
                endif
              endif
            endif
!
!  End of loop over j
!
          enddo jqloop
!
!  For 0-D case we need to increment ind by 1 to allow for diagonal element that is skipped
!
          if (ndim.eq.0.and..not.lfixQi) then
            ind = ind + 1
          endif
!
!  Solve for charge derivatives and convert units
!
!  X
!
          qreaxFFfree(1:n) = dqx(1:n)*angstoev
!
!  Solve for X charge derivatives
!
          call dsptrs('U',n,1_i4,dpacked,ipivot,qreaxFFfree,n,ifail)
          do j = 1,nfree
            dqdxyz(indi+1,j) = dqdxyz(indi+1,j) + qreaxFFfree(j)
          enddo
!
!  Y
!
          qreaxFFfree(1:n) = dqy(1:n)*angstoev
!
!  Solve for Y charge derivatives
!
          call dsptrs('U',n,1_i4,dpacked,ipivot,qreaxFFfree,n,ifail)
          do j = 1,nfree
            dqdxyz(indi+2,j) = dqdxyz(indi+2,j) + qreaxFFfree(j)
          enddo
!
!  Z
!
          qreaxFFfree(1:n) = dqz(1:n)*angstoev
!
!  Solve for Z charge derivatives
!
          call dsptrs('U',n,1_i4,dpacked,ipivot,qreaxFFfree,n,ifail)
          do j = 1,nfree
            dqdxyz(indi+3,j) = dqdxyz(indi+3,j) + qreaxFFfree(j)
          enddo
!
!  Strains
!
          if (lstrain) then
            do ns = 1,nstrains
              qreaxFFfree(1:n) = dqs(1:n,ns)*angstoev
!
!  Solve for strain charge derivatives
!
!  NB: Sign for terms for strain needs to be negative
!
              call dsptrs('U',n,1_i4,dpacked,ipivot,qreaxFFfree,n,ifail)
              do j = 1,nfree
                dqds(ns,j) = dqds(ns,j) - qreaxFFfree(j)
              enddo
            enddo
          endif
!
!  Solve for charge-electric field derivatives
!
          qreaxFFfree(1:n) = 0.0_dp
          qreaxFFfree(nfi) = 1.0_dp
          call dsptrs('U',n,1_i4,dpacked,ipivot,qreaxFFfree,n,ifail)
!
          if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
!
!  Electric field
!
            do jj = 1,nfree
              j = nboatomRptr(nfreeptr(jj))
              indj = 3*(j - 1)
              dqdxyz(indj+1,i) = dqdxyz(indj+1,i) + fieldx*qreaxFFfree(jj)
              dqdxyz(indj+2,i) = dqdxyz(indj+2,i) + fieldy*qreaxFFfree(jj)
              dqdxyz(indj+3,i) = dqdxyz(indj+3,i) + fieldz*qreaxFFfree(jj)
            enddo
          endif
!   
!  Compute electric field derivatives of charges
! 
          do jj = 1,nfree
            j = nboatomRptr(nfreeptr(jj))
            dqdVxyz(1,i) = dqdVxyz(1,i) + xclat(j)*qreaxFFfree(jj)
            dqdVxyz(2,i) = dqdVxyz(2,i) + yclat(j)*qreaxFFfree(jj)
            dqdVxyz(3,i) = dqdVxyz(3,i) + zclat(j)*qreaxFFfree(jj)
          enddo
        endif
      endif
!
!  End of loop over i
!
    enddo
  endif
!     
!  Free workspace for inversion
!     
  deallocate(ipivot,stat=status)  
  if (status/=0) call deallocate_error('setreaxffQ','ipivot')
!**************************
!  Calculate self energy  *
!**************************
  eself = 0.0_dp
  do i = 1,nfree
    ii = nfreeptr(i)
    qi = qreaxFF(ii)
    ii = nboatomRptr(ii)
    nspeci = nbosptr(ii)
    esite = qi*occuf(ii)*(reaxFFchi(nspeci)+extpot(nrelf2a(ii))+qi*reaxFFmu(nspeci))
    eself = eself + esite
    siteenergy(ii) = siteenergy(ii) + esite
!
!  Shell structure correction
!
    if (abs(reaxFFshell(1,nspeci)).gt.1.0d-12) then
      if (qi.gt.reaxFFshell(2,nspeci)) then
        esite = occuf(ii)*reaxFFshell(1,nspeci)*(qi - reaxFFshell(2,nspeci))**4
        eself = eself + esite
        siteenergy(ii) = siteenergy(ii) + esite
      endif
    endif
  enddo
!*******************
!  Output results  *
!*******************
  if ((loutputQ.or.index(keyword,'debu').ne.0).and.ioproc) then
!
!  Charges
!
    write(ioout,'(//,''  Final charges from ReaxFF :'',/)')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    do i = 1,nboatom
      ii = nboatomRptr(i)
      if (ii.gt.0) then
        write(ioout,'(4x,i6,18x,i2,16x,f10.7)') i,nat(ii),qreaxFF(i)
      endif
    enddo
    write(ioout,'(''--------------------------------------------------------------------------------'',/)')
  endif
  if ((loutputQ.or.index(keyword,'debu').ne.0).and.ioproc.and.ldoQdrv) then
!
!  Charge derivatives
!
    if (ndim.eq.0) then
      write(ioout,'(/,''  First derivatives of ReaxFF charge distribution : (Angstroms**-1)'',/)')
    else
      write(ioout,'(/,''  First derivatives of ReaxFF charge distribution : '',/)')
      write(ioout,'(''  Strain :'',/)')
      write(ioout,'(''  Atom   '',6(i10))') (j,j=1,nstrains)
      do i = 1,numat
        write(ioout,'(i6,4x,6f10.6)') i,(dqds(j,i),j=1,nstrains)
      enddo
      write(ioout,'(/,''  Coordinate (Angstroms**-1) :'',/)')
    endif
    igroup = numat/6
    iresid = numat - igroup*6
    indi = 0
    if (igroup.gt.0) then
      do i = 1,igroup
        write(ioout,'(''  Atom   '',6(i10))') (indi+j,j=1,6)
        indj = 0
        do j = 1,numat
          write(ioout,'(i6,'' x '',4x,6f10.6)') j,(dqdxyz(indj+1,indi+k),k=1,6)
          write(ioout,'(i6,'' y '',4x,6f10.6)') j,(dqdxyz(indj+2,indi+k),k=1,6)
          write(ioout,'(i6,'' z '',4x,6f10.6)') j,(dqdxyz(indj+3,indi+k),k=1,6)
          indj = indj + 3
        enddo
        indi = indi + 6
        write(ioout,'(/)')
      enddo
    endif
    if (iresid.gt.0) then
      write(ioout,'(''  Atom   '',6(i10))') (indi+j,j=1,iresid)
      indj = 0
      do j = 1,numat
        write(ioout,'(i6,'' x '',4x,6f10.6)') j,(dqdxyz(indj+1,indi+k),k=1,iresid)
        write(ioout,'(i6,'' y '',4x,6f10.6)') j,(dqdxyz(indj+2,indi+k),k=1,iresid)
        write(ioout,'(i6,'' z '',4x,6f10.6)') j,(dqdxyz(indj+3,indi+k),k=1,iresid)
        indj = indj + 3
      enddo
    endif
    write(ioout,'(/)')
    write(ioout,'(''  First derivatives of charge distribution for electric field : (Angstroms/V)'',/)')
    write(ioout,'(''  Atom         X             Y            Z'')')
    do i = 1,numat
      write(ioout,'(i6,3(4x,f10.6))') i,(dqdVxyz(k,i),k=1,3)
    enddo
    write(ioout,'(/)')
  endif
! 
!  Transfer charges to regular arrays
!   
  qf(1:numat) = 0.0_dp
  do ii = 1,nboatom
    i = nboatomRptr(ii)
    qf(i) = qreaxFF(ii)
  enddo
  qa(1:nasym) = 0.0_dp
  do i = 1,numat
    nr = nrelf2a(i)
    qa(nr) = qf(i)
  enddo
!
!  Free local memory 
!
  if (literative) then
    deallocate(dpackedsave,stat=status)
    if (status/=0) call deallocate_error('setreaxffQ','dpackedsave')
  endif
  deallocate(dpacked,stat=status)
  if (status/=0) call deallocate_error('setreaxffQ','dpacked')
  deallocate(qreaxFFsave,stat=status)
  if (status/=0) call deallocate_error('setreaxffQ','qreaxFFsave')
  deallocate(qreaxFFfree,stat=status)
  if (status/=0) call deallocate_error('setreaxffQ','qreaxFFfree')
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('setreaxffQ','z')
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    deallocate(vfield,stat=status)
    if (status/=0) call deallocate_error('setreaxffQ','vfield')
  endif
  deallocate(nfreeptr,stat=status)
  if (status/=0) call deallocate_error('setreaxffQ','nfreeptr')
!
  if (ldoQdrv) then
    if (lstrain) then
      deallocate(dqs,stat=status)
      if (status/=0) call deallocate_error('setreaxffQ','dqs')
    endif
    deallocate(dqz,stat=status)
    if (status/=0) call deallocate_error('setreaxffQ','dqz')
    deallocate(dqy,stat=status)
    if (status/=0) call deallocate_error('setreaxffQ','dqy')
    deallocate(dqx,stat=status)
    if (status/=0) call deallocate_error('setreaxffQ','dqx')
  endif
!
  return
  end
