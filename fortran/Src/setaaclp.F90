  subroutine setaaclp
!
!  Setup tasks for the AA-CLP force field of Gavezzotti
!
!   6/23 Created from setc6
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
!  Julian Gale, CIC, Curtin University, June 2023
!
  use control
  use current
  use g_constants, only : kjmtoev
  use m_aaclp
  use m_eht,       only : eht_mol, ehtd
  use molecule,    only : nmol, nmollist, nmolatom, nmolptr
  use parallel,    only : nprocs, procid
  use species
  use two
  implicit none
!
!  Local variables
!
  integer(i4)                                     :: i
  integer(i4)                                     :: ii
  integer(i4)                                     :: nati
  integer(i4)                                     :: natj
  integer(i4)                                     :: natk
  integer(i4)                                     :: nb
  integer(i4)                                     :: nm
  integer(i4)                                     :: noxygen
  integer(i4)                                     :: noxygenptr(4)
  integer(i4)                                     :: ntypi
  integer(i4)                                     :: ntypj
  integer(i4)                                     :: ntypk
  integer(i4)                                     :: nh
  integer(i4)                                     :: nh1
  integer(i4)                                     :: nh2
  integer(i4)                                     :: nr
  integer(i4)                                     :: nr1
  integer(i4)                                     :: nr2
  integer(i4)                                     :: status
  real(dp)                                        :: deltaq
  real(dp)                                        :: eeht
  real(dp)                                        :: qh
  real(dp)                                        :: qo
  real(dp)                                        :: qr
  real(dp),     dimension(:),   allocatable, save :: qsum
!************************************************
!  Generate Mulliken charges for each molecule  *
!************************************************
  if (nprocs.gt.1) then
    allocate(qsum(numat),stat=status)
    if (status/=0) call outofmemory('setaaclp','qsum')
!
!  Distribute molecules across nodes rather than Hamiltonian
!
    qsum(1:numat) = 0.0_dp
    do nm = 1+procid,nmol,nprocs
      call eht_mol(nm,eeht,.false.)
!
!  Globalise charges
!
      do ii = 1,nmolatom(nm)
        i = nmollist(nmolptr(nm)+ii)
        qsum(i) = qf(i)
      enddo
    enddo
!
!  Globalise charges
!
    call sumall(qsum,qf,numat,"setaaclp","qf")
!
!  Set values in qa
!
    do i = 1,nasym
      qa(i) = qf(nrela2f(i))
    enddo
!
    deallocate(qsum,stat=status)
    if (status/=0) call deallocate_error('setaaclp','qsum')
  else
    do nm = 1,nmol
      call eht_mol(nm,eeht,.false.)
    enddo
  endif
!**************************************
!  Special case charge modifications  *
!**************************************
  do i = 1,numat
    nati = nat(i)
    ntypi = nftype(i)
    if (nati.eq.8) then
!
!  Oxygen special cases
!
      if (ntypi.eq.5) then
!
!  ROH groups
!
        if (nbonds(i).ne.2) then
          call outerror('Oxygen type O5 in AA-CLP has the wrong number of bonds',0_i4)
          call stopnow('setaaclp')
        endif
        natj = nat(nbonded(1,i))
        natk = nat(nbonded(2,i))
        if (natj.eq.1.and.natk.ne.1) then
          nh = nbonded(1,i)
          nr = nbonded(2,i)
        elseif (natj.ne.1.and.natk.eq.1) then
          nr = nbonded(1,i)
          nh = nbonded(2,i)
        endif
!
        qo = qf(i)
        qh = qf(nh)
        qr = qf(nr)
!
        qf(i) = 1.2_dp*qo
        qf(nh) = 1.1_dp*qh
        qf(nr) = qr - 0.1_dp*qh - 0.2_dp*qo
      elseif (ntypi.eq.2) then
!
!  HOH groups (water)
!
        if (nbonds(i).ne.2) then
          call outerror('Oxygen type O2 in AA-CLP has the wrong number of bonds',0_i4)
          call stopnow('setaaclp')
        endif
        natj = nat(nbonded(1,i))
        natk = nat(nbonded(2,i))
        if (natj.ne.1.or.natk.ne.1) then
          call outerror('Oxygen type O2 in AA-CLP should be bonded to 2 hydrogens',0_i4)
          call stopnow('setaaclp')
        endif
        nh1 = nbonded(1,i)
        nh2 = nbonded(2,i)
!
        qf(i) = -1.54_dp
        qf(nh1) = 0.77_dp
        qf(nh2) = 0.77_dp
      endif
    elseif (nati.eq.7) then
!
!  Nitrogen special cases
!
      if (ntypi.eq.4) then
!
!  Nitro groups
!
        if (nbonds(i).le.2) then
          call outerror('Nitrogen type N5 in AA-CLP has the wrong number of bonds',0_i4)
          call stopnow('setaaclp')
        endif
!
!  Find oxygens within group
!
        noxygen = 0
        do nb = 1,nbonds(i)
          natj = nat(nbonded(nb,i))
          ntypj = nftype(nbonded(nb,i))
          if (natj.eq.8.and.ntypj.eq.6) then
            noxygen = noxygen + 1
            noxygenptr(noxygen) = nb
          endif
        enddo
!
        if (noxygen.lt.1) then
          call outerror('Nitrogen type N5 in AA-CLP has no bonds to oxygen',0_i4)
          call stopnow('setaaclp')
        endif
!
        deltaq = -0.20_dp*qf(i)
        qf(i) = qf(i) + deltaq
        deltaq = deltaq/dble(noxygen)
        do nb = 1,noxygen
          qf(noxygenptr(nb)) = qf(noxygenptr(nb)) - deltaq
        enddo
      endif
    elseif (nati.eq.16) then
!
!  Sulphur special cases
!
      if (ntypi.eq.1) then
!
!  RSR' groups
!
        if (nbonds(i).ne.2) then
          call outerror('Sulphur type S1 in AA-CLP has the wrong number of bonds',0_i4)
          call stopnow('setaaclp')
        endif
        nr1 = nbonded(1,i)
        nr2 = nbonded(2,i)
        natj = nat(nr1)
        ntypj = nftype(nr1)
        natk = nat(nr2)
        ntypk = nftype(nr2)
!
        if (ntypj.eq.3.and.ntypk.eq.3) then
          deltaq = 0.1_dp*qf(i)
        elseif (ntypj.eq.3.and.ntypk.eq.4) then
          deltaq = 0.4_dp*qf(i)
        elseif (ntypj.eq.4.and.ntypk.eq.3) then
          deltaq = 0.4_dp*qf(i)
        elseif (ntypj.eq.4.and.ntypk.eq.4) then
          deltaq = 1.5_dp*qf(i)
        endif
!
        qf(i) = qf(i) + deltaq
        qf(nr1) = qf(nr1) - 0.5_dp*deltaq
        qf(nr2) = qf(nr2) - 0.5_dp*deltaq
      elseif (ntypi.eq.1) then
!
!  Sulfone/sulphoxide groups
!
        if (nbonds(i).le.2) then
          call outerror('Sulphur type S3 in AA-CLP has the wrong number of bonds',0_i4)
          call stopnow('setaaclp')
        endif
!
!  Find oxygens within group
!
        noxygen = 0
        do nb = 1,nbonds(i)
          natj = nat(nbonded(nb,i))
          ntypj = nftype(nbonded(nb,i))
          if (natj.eq.8.and.ntypj.eq.7) then
            noxygen = noxygen + 1
            noxygenptr(noxygen) = nb
          endif
        enddo
!
        if (noxygen.lt.1) then
          call outerror('Sulphur type S3 in AA-CLP has no bonds to oxygen',0_i4)
          call stopnow('setaaclp')
        endif
!
        deltaq = -0.33_dp*qf(i)
        qf(i) = qf(i) + deltaq
        deltaq = deltaq/dble(noxygen)
        do nb = 1,noxygen
          qf(noxygenptr(nb)) = qf(noxygenptr(nb)) - deltaq
        enddo
      endif
    elseif (nati.eq.9) then
!
!  Fluorine
!
      if (ntypi.eq.1) then
        if (nbonds(i).eq.1) then
!
!  Find carbon bonds
!
          noxygen = 0
          do nb = 1,nbonds(i)
            natj = nat(nbonded(nb,i))
            ntypj = nftype(nbonded(nb,i))
            if (natj.eq.6) then
              noxygen = noxygen + 1
              noxygenptr(noxygen) = nb
            endif
          enddo
!
!  Only do case of a single C-F bond
!
          if (noxygen.eq.1) then
            deltaq = -0.25_dp*qf(i)
            qf(i) = qf(i) + deltaq
            qf(noxygenptr(1)) = qf(noxygenptr(1)) - deltaq
          endif
        endif
      endif
    endif
  enddo
!
!  Update asymmetric unit charges
!
  do i = 1,nasym
    qa(i) = qf(nrela2f(i))
  enddo
!***********************************
!  Scale charges by AA-CLP factor  *
!***********************************
  do i = 1,nasym
    qa(i) = aaclp_qscale*qa(i)
  enddo
  do i = 1,numat
    qf(i) = aaclp_qscale*qf(i)
  enddo
!*************************************************************
!  Set up scale factors while handling units (kJ/mol -> eV)  *
!*************************************************************
  aaclp_dscale_ev = aaclp_dscale*kjmtoev
  aaclp_pscale_ev = aaclp_pscale*kjmtoev
  aaclp_rscale_ev = aaclp_rscale*kjmtoev
!
  return
  end
