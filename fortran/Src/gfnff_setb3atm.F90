  subroutine gfnff_setb3atm(ndim,numat,nnbr,maxnbr,nbrno,xnbr,ynbr,znbr)
!
!  This routine replaces the algorithm in the original GFNFF that determines the number of covalent
!  bonds between atoms and the list of 3-body atoms for ATM dispersion. A different
!  algorithm is needed due to the need to handle multiple images of the same atom with PBC
!
!  10/20 Created
!   5/21 Modified to accelerate
!
!  Julian Gale, Curtin University, May 2021
!
  use datatypes
  use m_gfnff_nbr3
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: nnbr(numat)
  integer(i4), intent(in)  :: maxnbr
  integer(i4), intent(in)  :: nbrno(maxnbr,numat)
  real(dp),    intent(in)  :: xnbr(maxnbr,numat)
  real(dp),    intent(in)  :: ynbr(maxnbr,numat)
  real(dp),    intent(in)  :: znbr(maxnbr,numat)
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: inb
  integer(i4)              :: j
  integer(i4)              :: jnb
  integer(i4)              :: k
  integer(i4)              :: knb
  integer(i4)              :: l
  integer(i4)              :: lnb
  integer(i4)              :: nb
  integer(i4)              :: nb3min
  integer(i4)              :: ni
  integer(i4)              :: ninb
  integer(i4)              :: nj
  integer(i4)              :: njnb
  integer(i4)              :: nk
  integer(i4)              :: nknb
  logical                  :: liisk
  logical                  :: lik1bond
  logical                  :: lil2bond
  logical                  :: lil3bond
  logical                  :: ljisl
  real(dp)                 :: rik(3)
  real(dp)                 :: rik2
  real(dp)                 :: ril(3)
  real(dp)                 :: ril2
  real(dp)                 :: rjl(3)
  real(dp)                 :: rjl2
!
!  Initialise arrays to store 1-4 (3 bond) connection info
!
  call changemaxb3atm
  nb3atm = 0
  nb3min = 0
!******************************
!  Set third neighbour shell  *
!******************************
  do i = 1,numat
!
!  Set up bonded neighbour info for i
!
    call gfnff_get_n3atoms(numat,i,2_i4)
!
    do ni = 1,nnbr(i)
      j = nbrno(ni,i)
!
!  Loop over neighbours of j
!
      do nj = 1,nnbr(j)
        k = nbrno(nj,j)
!
!  Check on k being i
!
        if (k.eq.i) then
          rik(1) = xnbr(ni,i) + xnbr(nj,j)
          rik(2) = ynbr(ni,i) + ynbr(nj,j)
          rik(3) = znbr(ni,i) + znbr(nj,j)
          rik2 = abs(rik(1)) + abs(rik(2)) + abs(rik(3))
          liisk = (rik2.lt.1.0d-2)
          if (liisk) cycle
        endif
        lik1bond = .false.    ! Flag for case where i is bonded to k
!
!  Loop over neighbours of k
!
        do nk = 1,nnbr(k)
          l = nbrno(nk,k)
          if (l.eq.j) then
!
!  Check on l being j
!
            rjl(1) = xnbr(nj,j) + xnbr(nk,k)
            rjl(2) = ynbr(nj,j) + ynbr(nk,k)
            rjl(3) = znbr(nj,j) + znbr(nk,k)
            rjl2 = abs(rjl(1)) + abs(rjl(2)) + abs(rjl(3))
            ljisl = (rjl2.lt.1.0d-2)
            if (ljisl) cycle
          endif
!
!  Reduce search by excluding repeats from the opposite direction
!
          if (l.lt.i) cycle
!
          if (l.eq.i) then
!
!  Check on l being i => i is bonded to k
!
            ril(1) = xnbr(ni,i) + xnbr(nj,j) + xnbr(nk,k)
            ril(2) = ynbr(ni,i) + ynbr(nj,j) + ynbr(nk,k)
            ril(3) = znbr(ni,i) + znbr(nj,j) + znbr(nk,k)
            ril2 = abs(ril(1)) + abs(ril(2)) + abs(ril(3))
            lik1bond = (ril2.lt.1.0d-2)
            if (lik1bond) cycle
          endif
          if (.not.lik1bond) then
            ril(1) = xnbr(ni,i) + xnbr(nj,j) + xnbr(nk,k)
            ril(2) = ynbr(ni,i) + ynbr(nj,j) + ynbr(nk,k)
            ril(3) = znbr(ni,i) + znbr(nj,j) + znbr(nk,k)
!
!  Final checks on whether l is bonded to i or j
!
            call gfnff_bond_check(l,ril(1),ril(2),ril(3),2_i4,lil2bond)
            if (lil2bond) cycle
!
!  Check on 3 bonds via another route
!
            lil3bond = .false.
            do nb = nb3min+1,nb3atm
              inb  = nb3list(1,nb)
              ninb = nb3list(2,nb)
              njnb = nb3list(3,nb)
              nknb = nb3list(4,nb)
              jnb = nbrno(ninb,inb)
              knb = nbrno(njnb,jnb)
              lnb = nbrno(nknb,knb)
!
!  Compute i-l vector and test against one being added
!
              if (i.ne.l) then
                if (inb.eq.i.and.lnb.eq.l) then
                  ril2 = abs(rilb3(1,nb) - ril(1)) + abs(rilb3(2,nb) - ril(2)) + abs(rilb3(3,nb) - ril(3))
                  lil3bond = (ril2.lt.1.0d-2)
                  if (lil3bond) exit
                endif
              else
                if (inb.eq.i.and.lnb.eq.l) then
                  ril2 = abs(rilb3(1,nb) - ril(1)) + abs(rilb3(2,nb) - ril(2)) + abs(rilb3(3,nb) - ril(3))
                  lil3bond = (ril2.lt.1.0d-2)
                  if (lil3bond) exit
                  ril2 = abs(rilb3(1,nb) + ril(1)) + abs(rilb3(2,nb) + ril(2)) + abs(rilb3(3,nb) + ril(3))
                  lil3bond = (ril2.lt.1.0d-2)
                  if (lil3bond) exit
                endif
              endif
            enddo
            if (.not.lil3bond) then
!
!  Valid 1-4 pair
!
              nb3atm = nb3atm + 1
              if (nb3atm.gt.maxb3atm) then
                maxb3atm = nb3atm + 6
                call changemaxb3atm
              endif
              nb3list(1,nb3atm) = i
              nb3list(2,nb3atm) = ni
              nb3list(3,nb3atm) = nj
              nb3list(4,nb3atm) = nk
              rilb3(1,nb3atm) = ril(1)
              rilb3(2,nb3atm) = ril(2)
              rilb3(3,nb3atm) = ril(3)
            endif
          endif
        enddo
      enddo
    enddo
    nb3min = nb3atm
  enddo

  end subroutine gfnff_setb3atm
