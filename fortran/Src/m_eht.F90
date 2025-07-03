  module m_eht
!
!  Module for Extended Huckel Theory (EHT) calculation of charges
!
    use datatypes
    use control,     only : ldebug
    use g_constants
    use iochannels,  only : ioout, lioproconly
    use parallel
!
!  Data for EHT
!
    integer(i4), parameter :: maxele_eht = 103
    integer(i4),      save :: natorb(maxele_eht)                   ! Number of orbitals per atom
    integer(i4),      save :: norbital_l(9)                        ! Orbital l value
    integer(i4),      save :: norbital_m(9)                        ! Orbital m value
    integer(i4),      save :: npqn(maxele_eht)                     ! Principal quantum number 
    integer(i4),      save :: nvalenceelectron(maxele_eht)         ! Number of valence electrons for each element
    real(dp),         save :: vsip_s(maxele_eht)                   ! Valence State Ionisation Potentials for s orbitals
    real(dp),         save :: vsip_sin(maxele_eht)                 ! Valence State Ionisation Potentials for s orbitals (initial values)
    real(dp),         save :: vsip_p(maxele_eht)                   ! Valence State Ionisation Potentials for p orbitals
    real(dp),         save :: vsip_pin(maxele_eht)                 ! Valence State Ionisation Potentials for p orbitals (initial values)
    real(dp),         save :: zeta_s(maxele_eht)                   ! Slater exponents for s orbitals
    real(dp),         save :: zeta_sin(maxele_eht)                 ! Slater exponents for s orbitals (initial values)
    real(dp),         save :: zeta_p(maxele_eht)                   ! Slater exponents for p orbitals
    real(dp),         save :: zeta_pin(maxele_eht)                 ! Slater exponents for p orbitals (initial values)
    real(dp),         save :: scale_eht = 1.75_dp                  ! Scale factor for EHT Hamiltonian
    real(dp),         save :: qscale_eht = 1.0_dp                  ! Scale factor for EHT charges
    real(dp),         save :: rmax_overlap = 10.0_dp               ! Maximum distance for calculating overlap integrals
    real(dp),         save :: fct(30)                              ! Factorials
!
!  Huckel k point sampling
!
    integer(i4),                        save :: maxnehtk = 1           ! Maximum number of k points for EHT
    integer(i4),                        save :: nehtk                  ! Number of k points for EHT
    logical,                            save :: lgamma_eht = .false.   ! If true, use gamma point even for periodic cases
    real(dp),                           save :: kspace_eht = 0.10_dp   ! Default k point spacing for EHT
    real(dp),    dimension(:), pointer, save :: wehtk        => null() ! Weight of k point for EHT
    real(dp),    dimension(:), pointer, save :: xehtk        => null() ! X component of k point for EHT
    real(dp),    dimension(:), pointer, save :: yehtk        => null() ! Y component of k point for EHT
    real(dp),    dimension(:), pointer, save :: zehtk        => null() ! Z component of k point for EHT
!
    integer(i4), private   :: i
!
!  Element information for orbitals
!
    data (nvalenceelectron(i),i=1,maxele_eht)/1,2,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2, &
                3,4,5,6,7,8,9,10,11,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12,3, &
                4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,4,5,6,7,8, &
                9,10,11,12,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/ 
    data norbital_l/0,1,1,1,2,2,2,2,2/
    data norbital_m/0,1,-1,0,0,1,-1,2,-2/
    data npqn/2*1,8*2,8*3,18*4,18*5,32*6,17*0/
    data natorb/2*1,0,7*4,1,7*4,4,4,9*9,7*4,2*4,9*9,7*4,2*2,14*8,9*9,7*4,17*0/
!
!  Valence state ionisation potentials for Extended Huckel Theory
!
    data (vsip_sin(i),i=1,maxele_eht)/-13.60_dp,-23.40_dp,-5.40_dp,-10.00_dp,-15.20_dp, &
              -21.40_dp,-26.00_dp,-32.30_dp,-40.00_dp,-43.20_dp, -5.10_dp, -9.00_dp, &
              -12.30_dp,-17.30_dp,-18.60_dp,-20.00_dp,-26.30_dp,  0.00_dp, -4.34_dp, &
               -7.00_dp, -8.87_dp, -8.97_dp, -8.81_dp, -8.66_dp, -9.75_dp, -9.10_dp, &
               -9.21_dp, -9.17_dp,-11.40_dp,-12.41_dp,-14.58_dp,-16.00_dp,-16.22_dp, &
              -20.50_dp,-22.07_dp,  0.00_dp, -4.18_dp, -6.62_dp,  0.00_dp, -9.87_dp, &
              -10.10_dp, -8.34_dp,-10.07_dp,-10.40_dp, -8.09_dp, -7.32_dp,  0.00_dp, &
              -11.80_dp,-12.60_dp,-16.16_dp,-18.80_dp,-20.80_dp,-18.00_dp,  0.00_dp, &
               -3.88_dp,  0.00_dp, -7.67_dp, -4.97_dp,  0.00_dp,  0.00_dp,  0.00_dp, &
               -4.86_dp,  0.00_dp, -5.44_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp, &
                0.00_dp, -5.35_dp, -6.05_dp, -4.84_dp,-10.10_dp, -8.26_dp, -9.36_dp, &
               -8.17_dp,-11.36_dp, -9.077_dp,-10.92_dp,-13.68_dp,-11.60_dp,-15.70_dp,&
              -15.19_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp, -5.19_dp, &
               -5.39_dp, -5.41_dp, -5.50_dp, -5.60_dp,  0.00_dp,  0.00_dp,  0.00_dp, &
                0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp, -6.74_dp/
    data (vsip_pin(i),i=1,maxele_eht)/ -0.00_dp, -0.00_dp, -3.50_dp, -6.00_dp, -8.50_dp, &
              -11.40_dp,-13.40_dp,-14.80_dp,-18.10_dp,-20.00_dp, -3.00_dp, -4.50_dp, &
               -6.50_dp, -9.20_dp,-14.00_dp,-11.00_dp,-14.20_dp,  0.00_dp, -2.73_dp, &
               -4.00_dp, -2.75_dp, -5.44_dp, -5.52_dp, -5.24_dp, -5.89_dp, -5.32_dp, &
               -5.29_dp, -6.27_dp, -6.06_dp, -6.53_dp, -6.75_dp, -9.00_dp,-12.16_dp, &
              -14.40_dp,-13.10_dp,  0.00_dp, -2.60_dp, -3.92_dp,  0.00_dp, -6.76_dp, &
               -6.86_dp, -5.24_dp, -5.40_dp, -6.87_dp, -4.57_dp, -3.75_dp,  0.00_dp, &
               -8.20_dp, -6.19_dp, -8.32_dp,-11.70_dp,-14.80_dp,-12.70_dp,  0.00_dp, &
               -2.49_dp,  0.00_dp, -5.01_dp, -4.97_dp,  0.00_dp,  0.00_dp,  0.00_dp, &
               -4.86_dp,  0.00_dp, -5.44_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp, &
                0.00_dp, -5.35_dp, -6.05_dp, -4.84_dp, -6.86_dp, -5.17_dp, -5.96_dp, &
               -4.81_dp, -4.50_dp, -5.475_dp, -5.55_dp, -8.47_dp, -5.80_dp, -8.00_dp,&
               -7.79_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp, -5.19_dp, &
               -5.39_dp, -5.41_dp, -5.50_dp, -5.60_dp,  0.00_dp,  0.00_dp,  0.00_dp, &
                0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp,  0.00_dp, -6.74_dp/
!
!  Exponents for Slater functions
!
    data (zeta_sin(i),i=1,maxele_eht)/ 1.300_dp, 1.688_dp, 0.650_dp, 0.975_dp, 1.300_dp, &
                1.625_dp, 1.950_dp, 2.275_dp, 2.425_dp, 2.879_dp, 0.733_dp, 1.100_dp, &
                1.167_dp, 1.383_dp, 1.750_dp, 2.122_dp, 2.183_dp, 0.000_dp, 0.874_dp, &
                1.200_dp, 1.300_dp, 1.075_dp, 1.300_dp, 1.700_dp, 0.970_dp, 1.900_dp, &
                2.000_dp, 1.825_dp, 2.200_dp, 2.010_dp, 1.770_dp, 2.160_dp, 2.230_dp, &
                2.440_dp, 2.588_dp, 0.000_dp, 0.997_dp, 1.214_dp, 0.000_dp, 1.817_dp, &
                1.890_dp, 1.960_dp, 2.018_dp, 2.080_dp, 2.135_dp, 2.190_dp, 0.000_dp, &
                0.000_dp, 1.903_dp, 2.120_dp, 2.323_dp, 2.510_dp, 2.679_dp, 0.000_dp, &
                1.060_dp, 0.000_dp, 2.140_dp, 1.799_dp, 0.000_dp, 0.000_dp, 0.000_dp, &
                1.400_dp, 0.000_dp, 1.369_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, &
                0.000_dp, 1.540_dp, 1.666_dp, 1.527_dp, 2.280_dp, 2.341_dp, 2.398_dp, &
                2.452_dp, 2.500_dp, 2.554_dp, 2.602_dp, 2.649_dp, 2.300_dp, 2.350_dp, &
                2.560_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, &
                1.834_dp, 0.000_dp, 1.914_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, &
                0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp/ 
    data (zeta_pin(i),i=1,maxele_eht)/ 0.000_dp, 0.000_dp, 0.650_dp, 0.975_dp, 1.300_dp, &
                1.625_dp, 1.950_dp, 2.275_dp, 2.425_dp, 2.879_dp, 0.733_dp, 1.100_dp, &
                1.167_dp, 1.383_dp, 1.300_dp, 1.827_dp, 1.733_dp, 0.000_dp, 0.874_dp, &
                1.200_dp, 1.300_dp, 0.675_dp, 1.300_dp, 1.700_dp, 0.970_dp, 1.900_dp, &
                2.000_dp, 1.125_dp, 2.200_dp, 1.700_dp, 1.550_dp, 1.850_dp, 1.890_dp, &
                2.070_dp, 2.131_dp, 0.000_dp, 0.997_dp, 1.214_dp, 0.000_dp, 1.776_dp, &
                1.850_dp, 1.900_dp, 1.984_dp, 2.040_dp, 2.100_dp, 2.152_dp, 0.000_dp, &
                0.000_dp, 1.677_dp, 1.820_dp, 1.999_dp, 2.160_dp, 2.322_dp, 0.000_dp, &
                1.060_dp, 0.000_dp, 2.080_dp, 1.799_dp, 0.000_dp, 0.000_dp, 0.000_dp, &
                1.400_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, &
                0.000_dp, 1.540_dp, 1.666_dp, 0.000_dp, 2.241_dp, 2.309_dp, 2.372_dp, &
                2.429_dp, 2.200_dp, 2.554_dp, 2.584_dp, 2.631_dp, 1.600_dp, 2.060_dp, &
                2.072_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, &
                1.834_dp, 0.000_dp, 1.914_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, &
                0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp, 0.000_dp/ 
!
!  Declare public and private attributes to control visibility external to module
!
    public :: eht, ehtd, eht_mol
    private :: ss, coeff, coeffs, aintgs, bintgs, binm, overlap, getrotationmatrix
    private :: changemaxnehtk

  contains

    subroutine eht(eeht,lmain)
!
!  Subroutine for Extended Huckel Theory calculation of charges
!
!   5/23 Created
!   6/23 Electrostatic cutoff squared now passed to rsearch subroutines
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
    use current
    use element
    use maths,         only : ldivide_and_conquer
    use realvectors
    use times
#ifdef TRACE
    use trace,         only : trace_in, trace_out
#endif
    implicit none
!
!  Passed variables
!
    logical,     intent(in)                      :: lmain   ! If true then this is a main call for output
    real(dp),    intent(inout)                   :: eeht
!
!  Local variables
!
    integer(i4)                                  :: i
    integer(i4)                                  :: ifail
    integer(i4)                                  :: ilaenv
    integer(i4),  dimension(:),   allocatable    :: iw2
    integer(i4)                                  :: j
    integer(i4)                                  :: k
    integer(i4)                                  :: m
    integer(i4)                                  :: mm
    integer(i4)                                  :: n
    integer(i4)                                  :: nb
    integer(i4)                                  :: nk
    integer(i4)                                  :: nn
    integer(i4)                                  :: nati
    integer(i4)                                  :: natj
    integer(i4)                                  :: nelectron
    integer(i4)                                  :: nmolonly
    integer(i4)                                  :: no
    integer(i4)                                  :: nocc
    integer(i4)                                  :: nor
    integer(i4)                                  :: norbi
    integer(i4)                                  :: norbj
    integer(i4)                                  :: norbsofari
    integer(i4)                                  :: norbsofari_next
    integer(i4)                                  :: norbsofarj
    integer(i4)                                  :: norbsofarj_next
    integer(i4)                                  :: norbtot
    integer(i4)                                  :: nsize1
    integer(i4)                                  :: nsize2
    integer(i4)                                  :: nsize3
    integer(i4)                                  :: status
    logical,                                save :: lfirst = .true.
    logical                                      :: lgamma        ! If true then this is a gamma point calculation or non-periodic
    logical                                      :: lself
    complex(dpc)                                 :: cHij
    complex(dpc)                                 :: cSij
    complex(dpc), dimension(:),   allocatable    :: cw1
    complex(dpc), dimension(:,:), allocatable    :: cHmat
    complex(dpc), dimension(:,:), allocatable    :: cSmat
    real(dp)                                     :: g_cpu_time
    real(dp)                                     :: cut2
    real(dp),     dimension(:,:), allocatable    :: Hmat
    real(dp),     dimension(:,:), allocatable    :: Smat
    real(dp),     dimension(:,:), allocatable    :: eig
    real(dp)                                     :: Hij
    real(dp)                                     :: oci
    real(dp)                                     :: ocj
    real(dp)                                     :: occ
    real(dp)                                     :: ofct
    complex(dpc)                                 :: phase
    real(dp)                                     :: qtrm
    real(dp)                                     :: r2
    real(dp)                                     :: rkr
    real(dp)                                     :: Spair(9,9)    ! Overlap matrix between a pair of atoms
    real(dp)                                     :: time1
    real(dp)                                     :: time2
    real(dp)                                     :: vsipi(4)
    real(dp)                                     :: vsipj(4)
    real(dp),     dimension(:),   allocatable    :: w1
    real(dp)                                     :: xal
    real(dp)                                     :: yal
    real(dp)                                     :: zal
    real(dp)                                     :: xcrd
    real(dp)                                     :: ycrd
    real(dp)                                     :: zcrd
    real(dp)                                     :: xk
    real(dp)                                     :: yk
    real(dp)                                     :: zk
    real(dp)                                     :: xkc
    real(dp)                                     :: ykc
    real(dp)                                     :: zkc
#ifdef TRACE
    call trace_in('eht')
#endif
!
    time1 = g_cpu_time()
!
!  On first call initialise EHT k point array memory
!
    if (lfirst) then
      lfirst = .false.
      call changemaxnehtk
    endif
!
!  Set radius squared for overlap integrals
!
    cut2 = rmax_overlap**2
!   
!  Set up factorials - stored as g_factorial(n) in location(n+1)
!   
    fct(1) = 1.0_dp
    do i = 1,29
      fct(i+1) = fct(i)*dble(i)
    enddo
!
!  Find total number of orbitals
!
    norbtot = 0
    nelectron = 0
    do i = 1,numat
      nati = nat(i)
      norbi = natorb(nati)
      norbtot = norbtot + norbi
      nelectron = nelectron + nvalenceelectron(nati)
!
!  Check that no elements have more than s and p orbitals as d/f are not currently supported
!
      if (norbi.gt.4) then
        call outerror('Elements with d and/or f orbitals are not currently supported for EHT',0_i4)
        call stopnow('eht')
      elseif (norbi.eq.0) then
        call outerror('Elements without orbitals are not currently supported for EHT',0_i4)
        call stopnow('eht')
      endif
    enddo
!           
!  Ensure we have k vectors
!        
    if (ndim.eq.3) then
      call kvector3D
    elseif (ndim.eq.2) then
      call kvector2D
    elseif (ndim.eq.1) then
      call kvector1D
    endif
!
    if (ndim.ne.0) then
      if (lgamma_eht) then
!
!  Force gamma point
!
        nehtk = 1
        wehtk(1) = 1.0_dp
        xehtk(1) = 0.0_dp
        yehtk(1) = 0.0_dp
        zehtk(1) = 0.0_dp
      else
!
!  Set k points
!
        call setehtkpt(ndim,kv)
      endif
    endif
!
!  Is this a gamma point calculation?
!
    if (ndim.eq.0) then
      lgamma = .true.
    else
      if (nehtk.eq.1) then
        lgamma = (abs(xehtk(1))+abs(yehtk(1))+abs(zehtk(1)).lt.1.0d-6)
      else
        lgamma = .false.
      endif
    endif
    if (lgamma) then
!---------------------------------------------------------------------------------------------------
!  Gamma point only                                                                                |
!---------------------------------------------------------------------------------------------------
!
!  Allocate Hamiltonian and overlap matrices
!
      allocate(Hmat(norbtot,norbtot),stat=status)
      if (status/=0) call outofmemory('eht','Hmat')
      allocate(Smat(norbtot,norbtot),stat=status)
      if (status/=0) call outofmemory('eht','Smat')
      allocate(eig(norbtot,1),stat=status)
      if (status/=0) call outofmemory('eht','eig')
!
!  Initialise Hmat
!
      Hmat(1:norbtot,1:norbtot) = 0.0_dp
      Smat(1:norbtot,1:norbtot) = 0.0_dp
!**************************************
!  Build Extended Huckel Hamiltonian  *
!**************************************
!
!  Loop over first atom
!
      norbsofari_next = 0
      do i = 1,numat
        xal = xclat(i)
        yal = yclat(i)
        zal = zclat(i)
        nati = nat(i)
        norbi = natorb(nati)
        norbsofari = norbsofari_next
        norbsofari_next = norbsofari + norbi
        oci = occuf(i)
!
        vsipi(1) = vsip_s(nati)
        if (norbi.gt.1) then
          vsipi(2) = vsip_p(nati)
          vsipi(3) = vsip_p(nati)
          vsipi(4) = vsip_p(nati)
        endif
!
!  Start of second atom loop
!
        norbsofarj_next = 0
        jloop: do j = 1,i
          natj = nat(j)
          norbj = natorb(natj)
          norbsofarj = norbsofarj_next
          norbsofarj_next = norbsofarj + norbj
          ocj = occuf(j)
!
          xcrd = xclat(j) - xal
          ycrd = yclat(j) - yal
          zcrd = zclat(j) - zal
!***********************
!  Find valid vectors  *
!***********************
          nor = 0
          if (ndim.eq.3) then
            call rsearch3D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
          elseif (ndim.eq.2) then
            call rsearch2D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
          elseif (ndim.eq.1) then
            call rsearch1D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
          elseif (i.ne.j) then   ! 0-D so exclude self term
            r2 = xcrd**2 + ycrd**2 + zcrd**2
            if (r2.lt.cut2) then
              nor = 1
              dist(1) = r2
              xtmp(1) = xcrd
              ytmp(1) = ycrd
              ztmp(1) = zcrd
            endif
          endif
          if (i.eq.j) then
!***************
!  Self terms  *
!***************
            ofct = oci*oci
            mm = norbsofari
            do m = 1,norbi
              mm = mm + 1
              Hij = ofct*vsipi(m)
              Hmat(mm,mm) = Hmat(mm,mm) + Hij
              Smat(mm,mm) = Smat(mm,mm) + 1.0_dp
            enddo
          endif
!
          if (nor.eq.0) cycle jloop
!
          ofct = oci*ocj
!
          vsipj(1) = vsip_s(natj)
          if (norbj.gt.1) then
            vsipj(2) = vsip_p(natj)
            vsipj(3) = vsip_p(natj)
            vsipj(4) = vsip_p(natj)
          endif
!******************************
!  Loop over valid distances  *
!******************************
          do k = 1,nor
!
!  Sqrt distances
!
            dist(k) = sqrt(dist(k))
!*********************************************************
!  Calculate overlap matrix contribution from this pair  *
!*********************************************************
            call overlap(nati,natj,dist(k),xtmp(k),ytmp(k),ztmp(k),Spair)
!************************************
!  Add terms to Hamiltonian matrix  *
!************************************
            mm = norbsofari
            do m = 1,norbi
              mm = mm + 1
              nn = norbsofarj
              do n = 1,norbj
                nn = nn + 1
!
                Hij = 0.5_dp*ofct*scale_eht*Spair(m,n)*(vsipi(m) + vsipj(n))
!
                Hmat(nn,mm) = Hmat(nn,mm) + Hij
                Hmat(mm,nn) = Hmat(mm,nn) + Hij
!
                Smat(nn,mm) = Smat(nn,mm) + Spair(m,n)
                Smat(mm,nn) = Smat(mm,nn) + Spair(m,n)
              enddo
            enddo
          enddo
        enddo jloop
      enddo
!
!  Output matrices if desired
!
      if (ldebug.and.ioproc) then
        write(ioout,'(/,''  Extended Huckel Theory Hamiltonian matrix (eV) :'',/)')
        do i = 1,norbtot
          write(ioout,'(12f11.6)')(Hmat(j,i),j=1,norbtot)
        enddo
        write(ioout,'(/,''  Extended Huckel Theory overlap matrix :'',/)')
        do i = 1,norbtot
          write(ioout,'(12f11.6)')(Smat(j,i),j=1,norbtot)
        enddo
      endif
!*****************************************************
!  Find eigenvalues and eigenvectors of Hamiltonian  *
!*****************************************************
      if (ldivide_and_conquer) then
!
!  Initial dummy allocation of workspace
!
        allocate(w1(10),stat=status)
        if (status/=0) call outofmemory('eht','w1')
        allocate(iw2(10),stat=status)
        if (status/=0) call outofmemory('eht','iw2')
!
!  Set parameters for eigensolver by calling to obtain memory
!
        nsize1 = -1
        nsize2 = -1
!
        call dsygvd(1_i4,'V','U',norbtot,Hmat,norbtot,Smat,norbtot,eig(1,1),w1,nsize1,iw2,nsize2,ifail)
!
        nsize1 = nint(w1(1))
        nsize2 = iw2(1)
!
!  Deallocate workspace
!
        deallocate(iw2,stat=status)
        if (status/=0) call deallocate_error('eht','iw2')
        deallocate(w1,stat=status)
        if (status/=0) call deallocate_error('eht','w1')
!
!  Real allocation of workspace
!
        allocate(w1(nsize1),stat=status)
        if (status/=0) call outofmemory('eht','w1')
        allocate(iw2(nsize2),stat=status)
        if (status/=0) call outofmemory('eht','iw2')
!
!  Call eigensolver
!
        call dsygvd(1_i4,'V','U',norbtot,Hmat,norbtot,Smat,norbtot,eig(1,1),w1,nsize1,iw2,nsize2,ifail)
!
        deallocate(iw2,stat=status)
        if (status/=0) call deallocate_error('eht','iw2')
        deallocate(w1,stat=status)
        if (status/=0) call deallocate_error('eht','w1')
      else
!
!  Set parameters for eigensolver
!
        nb = ilaenv(1_i4,'DSYTRD','U',norbtot,-1_i4,-1_i4,-1_i4)
        nsize1 = max(1_i4,(nb+2_i4)*norbtot)
!
        allocate(w1(nsize1+10_i4),stat=status)
        if (status/=0) call outofmemory('eht','w1')
!
!  Call eigensolver
!
        call dsygv(1_i4,'V','U',norbtot,Hmat,norbtot,Smat,norbtot,eig(1,1),w1,nsize1,ifail)
!
        deallocate(w1,stat=status)
        if (status/=0) call deallocate_error('eht','w1')
      endif
!
!  Eigenvalues
!
      if (lmain.and.ioproc) then
        write(ioout,'(/,''  Extended Huckel Theory Orbital Eigenvalues (eV): '',/)')
        do i = 1,norbtot
          write(ioout,'(2x,i5,2x,f12.6)') i,eig(i,1)
        enddo
      endif
!
!  Compute total energy
!
      nocc = nelectron/2
      eeht = 0.0_dp
      do no = 1,nocc
        eeht = eeht + 2.0_dp*eig(no,1)
      enddo
      if (2*nocc.ne.nelectron) then
!
!  Odd electron case
!
        eeht = eeht + eig(nocc+1,1)
      endif
!*********************************
!  Compute Mulliken populations  *
!*********************************
      norbsofari_next = 0
      do i = 1,numat
        nati = nat(i)
        norbi = natorb(nati)
        norbsofari = norbsofari_next
        norbsofari_next = norbsofari + norbi
!
        qf(i) = 0.0_dp
!
        do no = 1,nocc
          mm = norbsofari
          do m = 1,norbi
            mm = mm + 1
!
!  Orbital population - NB Hmat now contains the eigenvectors
!
            qf(i) = qf(i) + 2.0_dp*Hmat(mm,no)**2
!
!  Overlap contribution - compute in 2 stages as only half of Smat remains after eigensolve
!
            do j = 1,mm-1
              qf(i) = qf(i) + 2.0_dp*Hmat(mm,no)*Hmat(j,no)*Smat(mm,j)
            enddo
            do j = mm+1,norbtot
              qf(i) = qf(i) + 2.0_dp*Hmat(mm,no)*Hmat(j,no)*Smat(j,mm)
            enddo
          enddo
        enddo
        if (2*nocc.ne.nelectron) then
          mm = norbsofari
          do m = 1,norbi
            mm = mm + 1
            qf(i) = qf(i) + Hmat(mm,nocc+1)**2
          enddo
          do j = 1,mm-1
            qf(i) = qf(i) + Hmat(mm,nocc+1)*Hmat(j,nocc+1)*Smat(mm,j)
          enddo
          do j = mm+1,norbtot
            qf(i) = qf(i) + Hmat(mm,nocc+1)*Hmat(j,nocc+1)*Smat(j,mm)
          enddo
        endif
      enddo
!
!  Compute total energy
!
      eeht = 0.0_dp
      do i = 1,nocc
        eeht = eeht + 2.0_dp*eig(i,1)
      enddo
    else
!---------------------------------------------------------------------------------------------------
!  Periodic systems with k points                                                                  |
!---------------------------------------------------------------------------------------------------
!
!  Assume that periodic systems will be closed shell for now
!
!  NB: This algorithm assumes that the system is a wide gap insulator, such that only the lowest
!      nelectron/2 bands will be occupied
!
      nocc = nelectron/2
!
!  Check just in case
!
      if (2*nocc.ne.nelectron) then
        call outerror('Periodic open shell systems are not currently supported for EHT',0_i4)
        call stopnow('eht')
      endif
!
!  Allocate Hamiltonian and overlap matrices
!
      allocate(cHmat(norbtot,norbtot),stat=status)
      if (status/=0) call outofmemory('eht','cHmat')
      allocate(cSmat(norbtot,norbtot),stat=status)
      if (status/=0) call outofmemory('eht','cSmat')
      allocate(eig(norbtot,nehtk),stat=status)
      if (status/=0) call outofmemory('eht','eig')
!
!  Initialise charges
!
      qf(1:numat) = 0.0_dp
      eeht = 0.0_dp
!************************
!  First over k points  *
!************************
      do nk = 1,nehtk
!
!  Generate Cartesian k vectors
!
        xk = xehtk(nk)
        yk = yehtk(nk)
        zk = zehtk(nk)
        if (ndim.eq.3) then
          xkc = xk*kv(1,1) + yk*kv(1,2) + zk*kv(1,3)
          ykc = xk*kv(2,1) + yk*kv(2,2) + zk*kv(2,3)
          zkc = xk*kv(3,1) + yk*kv(3,2) + zk*kv(3,3)
        elseif (ndim.eq.2) then
          xkc = xk*kv(1,1) + yk*kv(1,2)
          ykc = xk*kv(2,1) + yk*kv(2,2)
          zkc = 0.0_dp
        elseif (ndim.eq.1) then
          xkc = xk*kv(1,1)
          ykc = 0.0_dp
          zkc = 0.0_dp
        endif
!
!  Initialise Hmat
!
        cHmat(1:norbtot,1:norbtot) = 0.0_dpc
        cSmat(1:norbtot,1:norbtot) = 0.0_dpc
!**************************************
!  Build Extended Huckel Hamiltonian  *
!**************************************
!
!  Loop over first atom
!
        norbsofari_next = 0
        do i = 1,numat
          xal = xclat(i)
          yal = yclat(i)
          zal = zclat(i)
          nati = nat(i)
          norbi = natorb(nati)
          norbsofari = norbsofari_next
          norbsofari_next = norbsofari + norbi
          oci = occuf(i)
!
          vsipi(1) = vsip_s(nati)
          if (norbi.gt.1) then
            vsipi(2) = vsip_p(nati)
            vsipi(3) = vsip_p(nati)
            vsipi(4) = vsip_p(nati)
          endif
!
!  Start of second atom loop
!
          norbsofarj_next = 0
          jloopk1: do j = 1,i
            natj = nat(j)
            norbj = natorb(natj)
            norbsofarj = norbsofarj_next
            norbsofarj_next = norbsofarj + norbj
            ocj = occuf(j)
!
            xcrd = xclat(j) - xal
            ycrd = yclat(j) - yal
            zcrd = zclat(j) - zal
!***********************
!  Find valid vectors  *
!***********************
            nor = 0
            if (ndim.eq.3) then
              call rsearch3D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (ndim.eq.2) then
              call rsearch2D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (ndim.eq.1) then
              call rsearch1D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            endif
            if (i.eq.j) then
!***************
!  Self terms  *
!***************
              ofct = oci*oci
              mm = norbsofari
              do m = 1,norbi
                mm = mm + 1
                Hij = ofct*vsipi(m)
                cHmat(mm,mm) = cHmat(mm,mm) + dcmplx(Hij,0.0_dp)
                cSmat(mm,mm) = cSmat(mm,mm) + dcmplx(1.0_dp,0.0_dp)
              enddo
            endif
!
            if (nor.eq.0) cycle jloopk1
!
            ofct = oci*ocj
!
            vsipj(1) = vsip_s(natj)
            if (norbj.gt.1) then
              vsipj(2) = vsip_p(natj)
              vsipj(3) = vsip_p(natj)
              vsipj(4) = vsip_p(natj)
            endif
!******************************
!  Loop over valid distances  *
!******************************
            do k = 1,nor
!
!  Check for self term
!
              if (i.eq.j.and.dist(k).lt.1.0d-2) cycle
!
!  Sqrt distances
!
              dist(k) = sqrt(dist(k))
!*********************************************************
!  Calculate overlap matrix contribution from this pair  *
!*********************************************************
              call overlap(nati,natj,dist(k),xtmp(k),ytmp(k),ztmp(k),Spair)
!
!  Compute phase factor
!
              rkr = xkc*xtmp(k) + ykc*ytmp(k) + zkc*ztmp(k)
              phase = dcmplx(cos(rkr),-sin(rkr))
!************************************
!  Add terms to Hamiltonian matrix  *
!************************************
              mm = norbsofari
              do m = 1,norbi
                mm = mm + 1
                nn = norbsofarj
                do n = 1,norbj
                  nn = nn + 1
!
                  Hij = 0.5_dp*ofct*scale_eht*Spair(m,n)*(vsipi(m) + vsipj(n))
                  cHij = dcmplx(Hij,0.0_dp)*phase
                  cSij = dcmplx(Spair(m,n),0.0_dp)*phase
!
                  cHmat(nn,mm) = cHmat(nn,mm) + cHij
                  cHmat(mm,nn) = cHmat(mm,nn) + dconjg(cHij)
!
                  cSmat(nn,mm) = cSmat(nn,mm) + cSij
                  cSmat(mm,nn) = cSmat(mm,nn) + dconjg(cSij)
                enddo
              enddo
            enddo
          enddo jloopk1
        enddo
!
!  Output matrices if desired
!
        if (lmain.and.ioproc) then
          write(ioout,'(/,''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  EHT K point : '',i4,3(1x,f8.5),2x,'' weight = '',f8.5)') nk,xehtk(nk),yehtk(nk),zehtk(nk),wehtk(nk)
          write(ioout,'(''--------------------------------------------------------------------------------'',/)')
        endif
        if (ldebug.and.ioproc) then
          write(ioout,'(''  Extended Huckel Theory Hamiltonian matrix (eV) :'',/)')
          write(ioout,'(''  Real: '',/)')
          do i = 1,norbtot
            write(ioout,'(12f11.6)')(dble(cHmat(j,i)),j=1,norbtot)
          enddo
          write(ioout,'(/,''  Imaginary: '',/)')
          do i = 1,norbtot
            write(ioout,'(12f11.6)')(dimag(cHmat(j,i)),j=1,norbtot)
          enddo
          write(ioout,'(/,''  Extended Huckel Theory overlap matrix :'',/)')
          write(ioout,'(''  Real: '',/)')
          do i = 1,norbtot
            write(ioout,'(12f11.6)')(dble(cSmat(j,i)),j=1,norbtot)
          enddo
          write(ioout,'(/,''  Imaginary: '',/)')
          do i = 1,norbtot
            write(ioout,'(12f11.6)')(dimag(cSmat(j,i)),j=1,norbtot)
          enddo
        endif
!*****************************************************
!  Find eigenvalues and eigenvectors of Hamiltonian  *
!*****************************************************
        if (ldivide_and_conquer) then
!
!  Initial dummy allocation of workspace
!
          allocate(w1(10),stat=status)
          if (status/=0) call outofmemory('eht','w1')
          allocate(iw2(10),stat=status)
          if (status/=0) call outofmemory('eht','iw2')
          allocate(cw1(10),stat=status)
          if (status/=0) call outofmemory('eht','cw1')
!
!  Set parameters for eigensolver by calling to obtain memory
!
          nsize1 = -1
          nsize2 = -1
          nsize3 = -1
!
          call zhegvd(1_i4,'V','U',norbtot,cHmat,norbtot,cSmat,norbtot,eig(1,nk),cw1,nsize1,w1,nsize2,iw2,nsize3,ifail)
!
          nsize1 = nint(dble(cw1(1)))
          nsize2 = nint(w1(1))
          nsize3 = iw2(1)
!
!  Deallocate workspace
!
          deallocate(cw1,stat=status)
          if (status/=0) call deallocate_error('eht','cw1')
          deallocate(iw2,stat=status)
          if (status/=0) call deallocate_error('eht','iw2')
          deallocate(w1,stat=status)
          if (status/=0) call deallocate_error('eht','w1')
!
!  Real allocation of workspace
!
          allocate(cw1(nsize1),stat=status)
          if (status/=0) call outofmemory('eht','cw1')
          allocate(w1(nsize2),stat=status)
          if (status/=0) call outofmemory('eht','w1')
          allocate(iw2(nsize3),stat=status)
          if (status/=0) call outofmemory('eht','iw2')
!
!  Call eigensolver
!
          call zhegvd(1_i4,'V','U',norbtot,cHmat,norbtot,cSmat,norbtot,eig(1,nk),cw1,nsize1,w1,nsize2,iw2,nsize3,ifail)
!
          deallocate(iw2,stat=status)
          if (status/=0) call deallocate_error('eht','iw2')
          deallocate(w1,stat=status)
          if (status/=0) call deallocate_error('eht','w1')
          deallocate(cw1,stat=status)
          if (status/=0) call deallocate_error('eht','cw1')
        else
!
!  Set parameters for eigensolver
!
          nsize1 = 3*norbtot - 2
          nb = ilaenv(1_i4,'ZHETRD','U',norbtot,-1_i4,-1_i4,-1_i4)
          nsize2 = max(1_i4,(nb+1_i4)*norbtot)
!
          allocate(w1(nsize1),stat=status)
          if (status/=0) call outofmemory('eht','w1')
          allocate(cw1(nsize2),stat=status)
          if (status/=0) call outofmemory('eht','cw1')
!
!  Call eigensolver
!
          call zhegv(1_i4,'V','U',norbtot,cHmat,norbtot,cSmat,norbtot,eig(1,nk),cw1,nsize2,w1,ifail)
!
          deallocate(cw1,stat=status)
          if (status/=0) call deallocate_error('eht','cw1')
          deallocate(w1,stat=status)
          if (status/=0) call deallocate_error('eht','w1')
        endif
!
!  Eigenvalues
!
        if (lmain.and.ioproc) then
          write(ioout,'(/,''  Extended Huckel Theory Orbital Eigenvalues (eV): '',/)')
          do i = 1,norbtot
            write(ioout,'(2x,i5,2x,f12.6)') i,eig(i,nk)
          enddo
        endif
!*********************************
!  Compute Mulliken populations  *
!*********************************
!
!  Loop over first atom
!
        norbsofari_next = 0
        do i = 1,numat
          xal = xclat(i)
          yal = yclat(i)
          zal = zclat(i)
          nati = nat(i)
          norbi = natorb(nati)
          norbsofari = norbsofari_next
          norbsofari_next = norbsofari + norbi
          oci = occuf(i)
!
!  Start of second atom loop
!
          norbsofarj_next = 0
          jloopk2: do j = 1,i
            natj = nat(j)
            norbj = natorb(natj)
            norbsofarj = norbsofarj_next
            norbsofarj_next = norbsofarj + norbj
            ocj = occuf(j)
!
            xcrd = xclat(j) - xal
            ycrd = yclat(j) - yal
            zcrd = zclat(j) - zal
!***********************
!  Find valid vectors  *
!***********************
            nor = 0
            if (ndim.eq.3) then
              call rsearch3D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (ndim.eq.2) then
              call rsearch2D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (ndim.eq.1) then
              call rsearch1D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            endif
!
            ofct = oci*ocj
            occ = 2.0_dp*wehtk(nk)*ofct
!
            if (i.eq.j) then
!***************
!  Self terms  *
!***************
              do no = 1,nocc
                mm = norbsofari
                do m = 1,norbi
                  mm = mm + 1
                  cHij = cHmat(mm,no)
                  qf(i) = qf(i) + occ*dble(dconjg(cHij)*cHij)
                enddo
              enddo
            endif
!
            if (nor.eq.0) cycle jloopk2
!******************************
!  Loop over valid distances  *
!******************************
            do k = 1,nor
!
!  Check for self term
!
              if (i.eq.j.and.dist(k).lt.1.0d-2) cycle
!
!  Sqrt distances
!
              dist(k) = sqrt(dist(k))
!*********************************************************
!  Calculate overlap matrix contribution from this pair  *
!*********************************************************
              call overlap(nati,natj,dist(k),xtmp(k),ytmp(k),ztmp(k),Spair)
!
!  Compute phase factor
!
              rkr = xkc*xtmp(k) + ykc*ytmp(k) + zkc*ztmp(k)
              phase = dcmplx(cos(rkr),-sin(rkr))
!*********************************
!  Add contributions to charges  *
!*********************************
              mm = norbsofari
              do m = 1,norbi
                mm = mm + 1
                nn = norbsofarj
                do n = 1,norbj
                  nn = nn + 1
                  do no = 1,nocc
!
!  Divide overlap density between pair of atoms
!
                    qtrm = occ*dble(dconjg(cHmat(nn,no))*cHmat(mm,no)*phase)*Spair(m,n)
                    qf(i) = qf(i) + qtrm
                    qf(j) = qf(j) + qtrm
                  enddo
                enddo
              enddo
            enddo
          enddo jloopk2
        enddo
!
!  Add contribution to total energy from this k point
!
        do i = 1,nocc
          eeht = eeht + 2.0_dp*wehtk(nk)*eig(i,nk)
        enddo
      enddo
    endif
!
!  Convert orbital densities to charges
!
    do i = 1,numat
      qf(i) = nvalenceelectron(nat(i)) - qf(i)
    enddo
!
!  Update asymmetric unit charges
!
    do i = 1,nasym
     qa(i) = qf(nrela2f(i))
    enddo
!*******************
!  Output results  *
!*******************
    if ((lmain.or.ldebug).and.ioproc) then
      write(ioout,'(//,''  Final Mulliken charges from Extended Huckel Theory :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Extended Huckel Theory Energy = '',f12.6,'' eV'')') eeht
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
!
!  Free local memory
!
    if (lgamma) then
      deallocate(eig,stat=status)
      if (status/=0) call deallocate_error('eht','eig')
      deallocate(Smat,stat=status)
      if (status/=0) call deallocate_error('eht','Smat')
      deallocate(Hmat,stat=status)
      if (status/=0) call deallocate_error('eht','Hmat')
    else
      deallocate(eig,stat=status)
      if (status/=0) call deallocate_error('eht','eig')
      deallocate(cSmat,stat=status)
      if (status/=0) call deallocate_error('eht','cSmat')
      deallocate(cHmat,stat=status)
      if (status/=0) call deallocate_error('eht','cHmat')
    endif
!
!  Timing
!
    time2 = g_cpu_time()
    teht = teht + time2 - time1
#ifdef TRACE
    call trace_out('eht')
#endif
!
    return
    end
!
    subroutine ehtd(eeht,lmain)
!
!  Subroutine for Extended Huckel Theory calculation of charges
!  Parallel version with distributed memory.
!
!  NB: The algorithm changes the order of atoms to place all heavy atoms first
!      followed by hydrogens. This ensures that a blocksize of 4 will align
!      with the boundaries between the orbitals of different atoms.
!
!   6/23 Created from eht
!   6/23 Electrostatic cutoff squared now passed to rsearch subroutines
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
    use current
    use element
    use realvectors
    use times
#ifdef TRACE
    use trace,         only : trace_in, trace_out
#endif
    implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
    logical,     intent(in)                         :: lmain   ! If true then this is a main call for output
    real(dp),    intent(inout)                      :: eeht
#ifdef MPI
!
!  Local variables
!
    character(len=4)                                :: cinfo
    integer(i4)                                     :: i
    integer(i4)                                     :: ii
    integer(i4)                                     :: ilaenv
    integer(i4)                                     :: iloc
    integer(i4),  dimension(:),   allocatable, save :: iw2
    integer(i4)                                     :: j
    integer(i4)                                     :: jj
    integer(i4)                                     :: k
    integer(i4)                                     :: m
    integer(i4)                                     :: mm
    integer(i4)                                     :: n
    integer(i4)                                     :: nk
    integer(i4)                                     :: nn
    integer(i4)                                     :: nati
    integer(i4)                                     :: natj
    integer(i4)                                     :: natomloc
    integer(i4)                                     :: natomtot
    integer(i4),  dimension(:),   allocatable, save :: natomn2o      ! Pointer from new order of atoms to old (global)
    integer(i4),  dimension(:),   allocatable, save :: natoml2gold   ! Pointer to local atoms in original order
    integer(i4),  dimension(:),   allocatable, save :: natoml2gnew   ! Pointer to local atoms in new order
    integer(i4)                                     :: nelectron
    integer(i4)                                     :: nheavy
    integer(i4)                                     :: nhyd
    integer(i4)                                     :: nhydloc
    integer(i4)                                     :: nmolonly
    integer(i4)                                     :: no
    integer(i4)                                     :: nocc
    integer(i4)                                     :: noccloc
    integer(i4)                                     :: nor
    integer(i4)                                     :: norb
    integer(i4)                                     :: norbi
    integer(i4)                                     :: norbj
    integer(i4)                                     :: norbloc
    integer(i4),  dimension(:),   allocatable, save :: norbg2l
    integer(i4),  dimension(:),   allocatable, save :: norbl2g
    integer(i4),  dimension(:),   allocatable, save :: norb2node
    integer(i4)                                     :: norbsofari
    integer(i4)                                     :: norbsofari_next
    integer(i4)                                     :: norbsofarj
    integer(i4)                                     :: norbsofarj_next
    integer(i4)                                     :: norbtot
    integer(i4)                                     :: nproc2add
    integer(i4)                                     :: nsize1
    integer(i4)                                     :: nsize2
    integer(i4)                                     :: nsize3
    integer(i4)                                     :: status
    logical,                                   save :: lfirst = .true.
    logical                                         :: lgamma        ! If true then this is a gamma point calculation or non-periodic
    logical                                         :: lself
    complex(dpc)                                    :: cHij
    complex(dpc)                                    :: cSij
    complex(dpc), dimension(:),   allocatable, save :: cw1
    complex(dpc), dimension(:,:), allocatable, save :: cHmat
    complex(dpc), dimension(:,:), allocatable, save :: cSmat
    complex(dpc), dimension(:,:), allocatable, save :: cZmat
    real(dp)                                        :: g_cpu_time
    real(dp)                                        :: cut2
    real(dp),     dimension(:,:), allocatable, save :: Hmat
    real(dp),     dimension(:,:), allocatable, save :: Smat
    real(dp),     dimension(:,:), allocatable, save :: Zmat
    real(dp),     dimension(:,:), allocatable, save :: eig
    real(dp)                                        :: Hij
    real(dp)                                        :: oci
    real(dp)                                        :: ocj
    real(dp)                                        :: occ
    real(dp)                                        :: ofct
    complex(dpc)                                    :: phase
    real(dp)                                        :: qtrm
    real(dp)                                        :: r2
    real(dp)                                        :: rkr
    real(dp)                                        :: Spair(9,9)    ! Overlap matrix between a pair of atoms
    real(dp)                                        :: time1
    real(dp)                                        :: time2
    real(dp)                                        :: vsipi(4)
    real(dp)                                        :: vsipj(4)
    real(dp),     dimension(:),   allocatable, save :: w1
    real(dp)                                        :: xal
    real(dp)                                        :: yal
    real(dp)                                        :: zal
    real(dp)                                        :: xcrd
    real(dp)                                        :: ycrd
    real(dp)                                        :: zcrd
    real(dp)                                        :: xk
    real(dp)                                        :: yk
    real(dp)                                        :: zk
    real(dp)                                        :: xkc
    real(dp)                                        :: ykc
    real(dp)                                        :: zkc
!
!  Blacs / Scalapack integers / reals
!
    integer,      dimension(:),   allocatable, save :: icluster
    integer                                      :: idesc(9)
    integer                                      :: ifail
    integer,      dimension(:),   allocatable, save :: ifails
    integer                                      :: il
    integer                                      :: info
    integer                                      :: iu
    integer                                      :: ld
    integer                                      :: liwrk
    integer                                      :: lrwrk
    integer                                      :: lwrk
    integer                                      :: MPIerror
    integer                                      :: na
    integer                                      :: nb
    integer                                      :: nefound
    integer                                      :: nvfound
    integer                                      :: nloc
    integer                                      :: ntag
    integer                                      :: nnode
    integer                                      :: ntmp
    integer                                      :: Request
    integer(i4),  dimension(:),   allocatable, save :: StatMPI       ! Array for status from MPI
    real*8,       dimension(:),   allocatable, save :: dtmp
    real*8                                       :: etol
    real*8,       dimension(:),   allocatable, save :: gap
    real*8                                       :: orfac
    real*8                                       :: pdlamch
    real*8                                       :: vl
    real*8                                       :: vu
#ifdef TRACE
    call trace_in('ehtd')
#endif
!
    time1 = g_cpu_time()
!
!  On first call initialise EHT k point array memory
!
    if (lfirst) then
      lfirst = .false.
      call changemaxnehtk
    endif
!
!  Set radius squared for overlap integrals
!
    cut2 = rmax_overlap**2
!   
!  Set up factorials - stored as g_factorial(n) in location(n+1)
!   
    fct(1) = 1.0_dp
    do i = 1,29
      fct(i+1) = fct(i)*dble(i)
    enddo
!
!  Find total number of orbitals
!
    norbtot = 0
    nelectron = 0
    nheavy = 0
    nhyd = 0
    do i = 1,numat
      nati = nat(i)
      norbi = natorb(nati)
      norbtot = norbtot + norbi
      nelectron = nelectron + nvalenceelectron(nati)
!
!  Check that no elements have more than s and p orbitals as d/f are not currently supported
!
      if (norbi.gt.4) then
        call outerror('Elements with d and/or f orbitals are not currently supported for EHT',0_i4)
        call stopnow('ehtd')
      elseif (norbi.eq.0) then
        call outerror('Elements without orbitals are not currently supported for EHT',0_i4)
        call stopnow('ehtd')
      endif
!
      if (norbi.eq.1) then
        nhyd = nhyd + 1
      else
        nheavy = nheavy + 1
      endif
    enddo
!           
!  Ensure we have k vectors
!        
    if (ndim.eq.3) then
      call kvector3D
    elseif (ndim.eq.2) then
      call kvector2D
    elseif (ndim.eq.1) then
      call kvector1D
    endif
!
    allocate(natomn2o(numat),stat=status)
    if (status/=0) call outofmemory('ehtd','natomn2o')
    allocate(natoml2gnew(numat),stat=status)
    if (status/=0) call outofmemory('ehtd','natoml2gnew')
    allocate(natoml2gold(numat),stat=status)
    if (status/=0) call outofmemory('ehtd','natoml2gold')
    allocate(norbl2g(norbtot),stat=status)
    if (status/=0) call outofmemory('ehtd','norbl2g')
    allocate(norbg2l(norbtot),stat=status)
    if (status/=0) call outofmemory('ehtd','norbg2l')
    allocate(norb2node(norbtot),stat=status)
    if (status/=0) call outofmemory('ehtd','norb2node')
!
!  Find parallel distribution of orbitals
!
!  NB: Heavy atoms with s and p first, followed by hydrogens
!      Use a blocksize of 4 to ensure that boundaries align as much as possible
!
    natomloc = 0
    norbloc = 0
    nproc2add = 0
!
!  First pass: Assign heavy atoms
!
    norb = 0
    natomtot = 0
    do i = 1,numat
      nati = nat(i)
      norbi = natorb(nati)
      if (norbi.eq.4) then
        natomtot = natomtot + 1
        natomn2o(natomtot) = i
        if (procid.eq.nproc2add) then
          natomloc = natomloc + 1
          do j = 1,4
            norbl2g(norbloc+j) = norb + j
            norbg2l(norb+j) = norbloc + j
            norb2node(norb+j) = nproc2add
          enddo
          norbloc = norbloc + 4
          natoml2gold(natomloc) = i
          natoml2gnew(natomloc) = natomtot
        else
          do j = 1,4
            norb2node(norb+j) = nproc2add
            norbg2l(norb+j) = 0
          enddo
        endif
!
!  Increment processor pointer (with wrapping)
!
        nproc2add = nproc2add + 1
        if (nproc2add.ge.nprocs) nproc2add = 0
!
        norb = norb + norbi
      endif
    enddo
!
!  Second pass: Assign hydrogen atoms
!
    nhydloc = 0
    do i = 1,numat
      nati = nat(i)
      norbi = natorb(nati)
      if (norbi.eq.1) then
        natomtot = natomtot + 1
        natomn2o(natomtot) = i
        if (procid.eq.nproc2add) then
          natomloc = natomloc + 1
          norbloc = norbloc + 1
          norbl2g(norbloc) = norb + 1
          norbg2l(norb+1) = norbloc
          natoml2gold(natomloc) = i
          natoml2gnew(natomloc) = natomtot
        else
          norbg2l(norb+1) = 0
        endif
        norb = norb + 1
        norb2node(norb) = nproc2add
        nhydloc = nhydloc + 1
!
!  If number of hydrogens has reached 4 then increment processor (with wrapping) and reset counter
!
        if (nhydloc.eq.4) then
          nhydloc = 0
          nproc2add = nproc2add + 1
          if (nproc2add.ge.nprocs) nproc2add = 0
        endif
      endif
    enddo
!
!  Find out how many occupied orbitals are local for a blocksize of 4
!
    nocc = nelectron/2  ! NB: Need to handle odd numbers of electrons
    noccloc = 0
    nhydloc = 0   ! Here this is being used as a counter rather than for the number of hydrogen atoms
    nproc2add = 0
    do no = 1,nocc
      if (procid.eq.nproc2add) then
        noccloc = noccloc + 1
      endif
      nhydloc = nhydloc + 1
!
!  If number of local occupied orbitals has been incremented by 4 then move to next processor (with wrapping) and reset counter
!
      if (nhydloc.eq.4) then
        nhydloc = 0
        nproc2add = nproc2add + 1
        if (nproc2add.ge.nprocs) nproc2add = 0
      endif
    enddo
!  
!  Set up Blacs descriptors for arrays
!   
    na = norbtot
    nloc = norbloc 
    nb = 4
    ld = norbtot
    call descinit(idesc, na, na, nb, nb, 0, 0, iBlacsContext, ld, ifail)
!     
    if (ifail.ne.0) then
      call outerror('initialisation in descinit failed - idesc ',0_i4)
      call stopnow('ehtd')
    endif
!   
!  Set tolerance for eigenvalues
!   
    etol = 2.0_dp*pdlamch(idesc(2),'U')
    orfac = 1.0d-3
!
!  Allocate arrays for parallel eigensolvers
!
    allocate(icluster(2*nprocs),stat=status)
    if (status/=0) call outofmemory('ehtd','icluster')
    allocate(ifails(norbtot),stat=status)
    if (status/=0) call outofmemory('ehtd','ifails')
    allocate(gap(nprocs),stat=status)
    if (status/=0) call outofmemory('ehtd','gap')
!
    if (ndim.ne.0) then
      if (lgamma_eht) then
!
!  Force gamma point
!
        nehtk = 1
        wehtk(1) = 1.0_dp
        xehtk(1) = 0.0_dp
        yehtk(1) = 0.0_dp
        zehtk(1) = 0.0_dp
      else
!
!  Set k points
!
        call setehtkpt(ndim,kv)
      endif
    endif
!
!  Is this a gamma point calculation?
!
    if (ndim.eq.0) then
      lgamma = .true.
    else
      if (nehtk.eq.1) then
        lgamma = (abs(xehtk(1))+abs(yehtk(1))+abs(zehtk(1)).lt.1.0d-6)
      else
        lgamma = .false.
      endif
    endif
!
!  Initialise charges and energy
!
    qf(1:numat) = 0.0_dp
    eeht = 0.0_dp
!
    if (lgamma) then
!---------------------------------------------------------------------------------------------------
!  Gamma point only                                                                                |
!---------------------------------------------------------------------------------------------------
!
!  Allocate Hamiltonian and overlap matrices
!
      allocate(Hmat(norbtot,max(1_i4,norbloc)),stat=status)
      if (status/=0) call outofmemory('ehtd','Hmat')
      allocate(Smat(norbtot,max(1_i4,norbloc)),stat=status)
      if (status/=0) call outofmemory('ehtd','Smat')
      allocate(Zmat(norbtot,max(1_i4,norbloc)),stat=status)
      if (status/=0) call outofmemory('ehtd','Zmat')
      allocate(eig(norbtot,1),stat=status)
      if (status/=0) call outofmemory('ehtd','eig')
!
!  Initialise Hmat
!
      if (norbloc.gt.0) then
        Hmat(1:norbtot,1:norbloc) = 0.0_dp
        Smat(1:norbtot,1:norbloc) = 0.0_dp
!**************************************
!  Build Extended Huckel Hamiltonian  *
!**************************************
!
!  Loop over first atom
!
        norbsofari_next = 0
        do ii = 1,natomloc
          i = natoml2gold(ii)
          xal = xclat(i)
          yal = yclat(i)
          zal = zclat(i)
          nati = nat(i)
          norbi = natorb(nati)
          norbsofari = norbsofari_next
          norbsofari_next = norbsofari + norbi
          oci = occuf(i)
!
          vsipi(1) = vsip_s(nati)
          if (norbi.gt.1) then
            vsipi(2) = vsip_p(nati)
            vsipi(3) = vsip_p(nati)
            vsipi(4) = vsip_p(nati)
          endif
!
!  Start of second atom loop
!
          norbsofarj_next = 0
          jloop: do jj = 1,numat
            j = natomn2o(jj)
            natj = nat(j)
            norbj = natorb(natj)
            norbsofarj = norbsofarj_next
            norbsofarj_next = norbsofarj + norbj
            ocj = occuf(j)
!
            xcrd = xclat(j) - xal
            ycrd = yclat(j) - yal
            zcrd = zclat(j) - zal
!***********************
!  Find valid vectors  *
!***********************
            nor = 0
            if (ndim.eq.3) then
              call rsearch3D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (ndim.eq.2) then
              call rsearch2D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (ndim.eq.1) then
              call rsearch1D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (i.ne.j) then   ! 0-D so exclude self term
              r2 = xcrd**2 + ycrd**2 + zcrd**2
              if (r2.lt.cut2) then
                nor = 1
                dist(1) = r2
                xtmp(1) = xcrd
                ytmp(1) = ycrd
                ztmp(1) = zcrd
              endif
            endif
            if (i.eq.j) then
!***************
!  Self terms  *
!***************
              ofct = oci*oci
              mm = norbsofari
              nn = norbsofarj
              do m = 1,norbi
                mm = mm + 1
                nn = nn + 1
                Hij = ofct*vsipi(m)
                Hmat(nn,mm) = Hmat(nn,mm) + Hij
                Smat(nn,mm) = Smat(nn,mm) + 1.0_dp
              enddo
            endif
!
            if (nor.eq.0) cycle jloop
!
            ofct = oci*ocj
!
            vsipj(1) = vsip_s(natj)
            if (norbj.gt.1) then
              vsipj(2) = vsip_p(natj)
              vsipj(3) = vsip_p(natj)
              vsipj(4) = vsip_p(natj)
            endif
!******************************
!  Loop over valid distances  *
!******************************
            do k = 1,nor
!
!  Sqrt distances
!
              dist(k) = sqrt(dist(k))
!*********************************************************
!  Calculate overlap matrix contribution from this pair  *
!*********************************************************
              call overlap(nati,natj,dist(k),xtmp(k),ytmp(k),ztmp(k),Spair)
!************************************
!  Add terms to Hamiltonian matrix  *
!************************************
              mm = norbsofari
              do m = 1,norbi
                mm = mm + 1
                nn = norbsofarj
                do n = 1,norbj
                  nn = nn + 1
!
                  Hij = 0.5_dp*ofct*scale_eht*Spair(m,n)*(vsipi(m) + vsipj(n))
!
                  Hmat(nn,mm) = Hmat(nn,mm) + Hij
!
                  Smat(nn,mm) = Smat(nn,mm) + Spair(m,n)
                enddo
              enddo
            enddo
          enddo jloop
        enddo
      endif
!
!  Output matrices if desired
!
      if (ldebug) then
        if (ioproc) then
          write(ioout,'(/,''  Extended Huckel Theory Hamiltonian matrix (eV) :'',/)')
        endif
!
        if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
          ntmp = norbtot
          ntag = 1
          allocate(dtmp(norbtot),stat=status)
          if (status/=0) call outofmemory('ehtd','dtmp')
          allocate(StatMPI(MPI_Status_Size),stat=status)
          if (status/=0) call outofmemory('ehtd','StatMPI')
        endif
        call mpbarrier
        do i = 1,norbtot
          iloc = norbg2l(i)
          if (lioproconly.and.norb2node(i).ne.0_i4) then
!
!  Post receive
!
            if (ioproc) then
              nnode = norb2node(i)
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_GULP,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (iloc.gt.0) then
              dtmp(1:norbtot) = Hmat(1:norbtot,iloc)
!
!  Post send
!
              call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_GULP,Request,MPIerror)
            endif
            if (ioproc.or.iloc.gt.0) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
            if (ioproc) then
!
!  Write on I/O node
!
              write(ioout,'(12f11.6)')(dtmp(j),j=1,norbtot)
            endif
          else
            if (iloc.gt.0) then
              write(ioout,'(12f11.6)')(Hmat(j,iloc),j=1,norbtot)
            endif
          endif
          call mpbarrier
        enddo
!
        if (ioproc) then
          write(ioout,'(/,''  Extended Huckel Theory overlap matrix :'',/)')
        endif
        call mpbarrier
        do i = 1,norbtot
          iloc = norbg2l(i)
          if (lioproconly.and.norb2node(i).ne.0_i4) then
!
!  Post receive
!  
            if (ioproc) then
              nnode = norb2node(i)
              call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                             ntag,MPI_Comm_GULP,Request,MPIerror)
            endif
!
!  Pass data to ioproc for writing
!
            if (iloc.gt.0) then
              dtmp(1:norbtot) = Smat(1:norbtot,iloc)
!
!  Post send
!                  
              call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                             ntag,MPI_Comm_GULP,Request,MPIerror)
            endif
            if (ioproc.or.iloc.gt.0) then
              call MPI_WaitAll(1,Request,StatMPI,MPIerror)
            endif
            if (ioproc) then
!
!  Write on I/O node
!
              write(ioout,'(12f11.6)')(dtmp(j),j=1,norbtot)
            endif
          else
            if (iloc.gt.0) then
              write(ioout,'(12f11.6)')(Smat(j,iloc),j=1,norbtot)
            endif
          endif
          call mpbarrier
        enddo
!
        if (lioproconly) then
          deallocate(StatMPI,stat=status)
          if (status/=0) call deallocate_error('ehtd','StatMPI')
          deallocate(dtmp,stat=status)
          if (status/=0) call deallocate_error('ehtd','dtmp')
        endif
      endif
!*****************************************************
!  Find eigenvalues and eigenvectors of Hamiltonian  *
!*****************************************************
!
!  Initial dummy allocation of workspace
!
      allocate(w1(10),stat=status)
      if (status/=0) call outofmemory('ehtd','w1')
      allocate(iw2(10),stat=status)
      if (status/=0) call outofmemory('ehtd','iw2')
!
!  Set parameters for eigensolver by calling to obtain memory
!
      na = norbtot
      nsize1 = -1
      nsize2 = -1
!
      call pdsygvx(1,'V','A','U',na,Hmat,1,1,idesc,Smat,1,1,idesc,vl,vu,il,iu,etol,nefound,nvfound, &
                   eig(1,1),orfac,Zmat,1,1,idesc,w1,nsize1,iw2,nsize2,ifails,icluster,gap,info)
      if (info.ne.0) then
        write(cinfo,'(i4)') info
        call outerror('Error calling pdsygvx in EHTD: Info = '//cinfo,0_i4)
        call stopnow('ehtd')
      endif
!
      nsize1 = nint(w1(1))
      nsize2 = iw2(1)
!
!  Deallocate workspace
!
      deallocate(iw2,stat=status)
      if (status/=0) call deallocate_error('ehtd','iw2')
      deallocate(w1,stat=status)
      if (status/=0) call deallocate_error('ehtd','w1')
!
!  Real allocation of workspace
!
      allocate(w1(nsize1),stat=status)
      if (status/=0) call outofmemory('ehtd','w1')
      allocate(iw2(nsize2),stat=status)
      if (status/=0) call outofmemory('ehtd','iw2')
!
!  Call eigensolver
!
      call pdsygvx(1,'V','A','U',na,Hmat,1,1,idesc,Smat,1,1,idesc,vl,vu,il,iu,etol,nefound,nvfound, &
                   eig(1,1),orfac,Zmat,1,1,idesc,w1,nsize1,iw2,nsize2,ifails,icluster,gap,info)
!
      if (info.ne.0) then
        write(cinfo,'(i4)') info
        call outerror('Error calling pdsygvx in EHTD: Info = '//cinfo,0_i4)
        call stopnow('ehtd')
      endif
!
      deallocate(iw2,stat=status)
      if (status/=0) call deallocate_error('ehtd','iw2')
      deallocate(w1,stat=status)
      if (status/=0) call deallocate_error('ehtd','w1')
!
!  Eigenvalues
!
      if (lmain.and.ioproc) then
        write(ioout,'(/,''  Extended Huckel Theory Orbital Eigenvalues (eV): '',/)')
        do i = 1,norbtot
          write(ioout,'(2x,i5,2x,f12.6)') i,eig(i,1)
        enddo
      endif
!
!  Compute total energy
!
      do no = 1,nocc
        eeht = eeht + 2.0_dp*eig(no,1)
      enddo
      if (2*nocc.ne.nelectron) then
!
!  Odd electron case
!
        eeht = eeht + eig(nocc+1,1)
      endif
!********************************************************************************************************
!  Compute Mulliken populations while minimising communication by rebuilding overlap matrix on the fly  *
!********************************************************************************************************
!
!  Loop over first atom
!
      if (noccloc.gt.0) then
        norbsofari_next = 0
        do ii = 1,numat
          i = natomn2o(ii)
          xal = xclat(i)
          yal = yclat(i)
          zal = zclat(i)
          nati = nat(i)
          norbi = natorb(nati)
          norbsofari = norbsofari_next
          norbsofari_next = norbsofari + norbi
          oci = occuf(i)
!
!  Start of second atom loop
!
          norbsofarj_next = 0
          jloop2: do jj = 1,ii
            j = natomn2o(jj)
            natj = nat(j)
            norbj = natorb(natj)
            norbsofarj = norbsofarj_next
            norbsofarj_next = norbsofarj + norbj
            ocj = occuf(j)
!
            xcrd = xclat(j) - xal
            ycrd = yclat(j) - yal
            zcrd = zclat(j) - zal
!***********************
!  Find valid vectors  *
!***********************
            nor = 0
            if (ndim.eq.3) then
              call rsearch3D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (ndim.eq.2) then
              call rsearch2D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (ndim.eq.1) then
              call rsearch1D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
            elseif (i.ne.j) then   ! 0-D so exclude self term
              r2 = xcrd**2 + ycrd**2 + zcrd**2
              if (r2.lt.cut2) then
                nor = 1
                dist(1) = r2
                xtmp(1) = xcrd
                ytmp(1) = ycrd
                ztmp(1) = zcrd
              endif
            endif
!
            ofct = oci*ocj
            occ = 2.0_dp*ofct
!
            if (i.eq.j) then
!***************
!  Self terms  *
!***************
              do no = 1,noccloc
                mm = norbsofari
                do m = 1,norbi
                  mm = mm + 1
                  qtrm = occ*Zmat(mm,no)**2
                  qf(i) = qf(i) + qtrm
                enddo
              enddo
            endif
!
            if (nor.eq.0) cycle jloop2
!******************************
!  Loop over valid distances  *
!******************************
            do k = 1,nor
!
!  Check for self term
!
              if (i.eq.j.and.dist(k).lt.1.0d-2) cycle
!
!  Sqrt distances
!
              dist(k) = sqrt(dist(k))
!*********************************************************
!  Calculate overlap matrix contribution from this pair  *
!*********************************************************
              call overlap(nati,natj,dist(k),xtmp(k),ytmp(k),ztmp(k),Spair)
!*********************************
!  Add contributions to charges  *
!*********************************
              mm = norbsofari
              do m = 1,norbi
                mm = mm + 1
                nn = norbsofarj
                do n = 1,norbj
                  nn = nn + 1
                  do no = 1,noccloc
!
!  Divide overlap density between pair of atoms
!
                    qtrm = occ*Zmat(nn,no)*Zmat(mm,no)*Spair(m,n)
                    qf(i) = qf(i) + qtrm
                    qf(j) = qf(j) + qtrm
                  enddo
                enddo
              enddo
            enddo
          enddo jloop2
        enddo
      endif
    else
!---------------------------------------------------------------------------------------------------
!  Periodic systems with k points                                                                  |
!---------------------------------------------------------------------------------------------------
!
!  Assume that periodic systems will be closed shell for now
!
!  NB: This algorithm assumes that the system is a wide gap insulator, such that only the lowest
!      nelectron/2 bands will be occupied
!
!  Check just in case
!
      if (2*nocc.ne.nelectron) then
        call outerror('Periodic open shell systems are not currently supported for EHT',0_i4)
        call stopnow('ehtd')
      endif
!
!  Allocate Hamiltonian and overlap matrices
!
      allocate(cHmat(norbtot,max(1_i4,norbloc)),stat=status)
      if (status/=0) call outofmemory('ehtd','cHmat')
      allocate(cSmat(norbtot,max(1_i4,norbloc)),stat=status)
      if (status/=0) call outofmemory('ehtd','cSmat')
      allocate(cZmat(norbtot,max(1_i4,norbloc)),stat=status)
      if (status/=0) call outofmemory('ehtd','cZmat')
      allocate(eig(norbtot,nehtk),stat=status)
      if (status/=0) call outofmemory('ehtd','eig')
!************************
!  First over k points  *
!************************
      do nk = 1,nehtk
!
!  Generate Cartesian k vectors
!
        xk = xehtk(nk)
        yk = yehtk(nk)
        zk = zehtk(nk)
        if (ndim.eq.3) then
          xkc = xk*kv(1,1) + yk*kv(1,2) + zk*kv(1,3)
          ykc = xk*kv(2,1) + yk*kv(2,2) + zk*kv(2,3)
          zkc = xk*kv(3,1) + yk*kv(3,2) + zk*kv(3,3)
        elseif (ndim.eq.2) then
          xkc = xk*kv(1,1) + yk*kv(1,2)
          ykc = xk*kv(2,1) + yk*kv(2,2)
          zkc = 0.0_dp
        elseif (ndim.eq.1) then
          xkc = xk*kv(1,1)
          ykc = 0.0_dp
          zkc = 0.0_dp
        endif
        if (norbloc.gt.0) then
!
!  Initialise Hmat
!
          cHmat(1:norbtot,1:norbloc) = 0.0_dpc
          cSmat(1:norbtot,1:norbloc) = 0.0_dpc
!**************************************
!  Build Extended Huckel Hamiltonian  *
!**************************************
!
!  Loop over first atom
!
          norbsofari_next = 0
          do ii = 1,natomloc
            i = natoml2gold(ii)
            xal = xclat(i)
            yal = yclat(i)
            zal = zclat(i)
            nati = nat(i)
            norbi = natorb(nati)
            norbsofari = norbsofari_next
            norbsofari_next = norbsofari + norbi
            oci = occuf(i)
!
            vsipi(1) = vsip_s(nati)
            if (norbi.gt.1) then
              vsipi(2) = vsip_p(nati)
              vsipi(3) = vsip_p(nati)
              vsipi(4) = vsip_p(nati)
            endif
!
!  Start of second atom loop
!
            norbsofarj_next = 0
            jloopk1: do jj = 1,numat
              j = natomn2o(jj)
              natj = nat(j)
              norbj = natorb(natj)
              norbsofarj = norbsofarj_next
              norbsofarj_next = norbsofarj + norbj
              ocj = occuf(j)
!
              xcrd = xclat(j) - xal
              ycrd = yclat(j) - yal
              zcrd = zclat(j) - zal
!***********************
!  Find valid vectors  *
!***********************
              nor = 0
              if (ndim.eq.3) then
                call rsearch3D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
              elseif (ndim.eq.2) then
                call rsearch2D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
              elseif (ndim.eq.1) then
                call rsearch1D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
              endif
              if (i.eq.j) then
!***************
!  Self terms  *
!***************
                ofct = oci*oci
                mm = norbsofari
                nn = norbsofarj
                do m = 1,norbi
                  mm = mm + 1
                  nn = nn + 1
                  Hij = ofct*vsipi(m)
                  cHmat(nn,mm) = cHmat(nn,mm) + dcmplx(Hij,0.0_dp)
                  cSmat(nn,mm) = cSmat(nn,mm) + dcmplx(1.0_dp,0.0_dp)
                enddo
              endif
!
              if (nor.eq.0) cycle jloopk1
!
              ofct = oci*ocj
!
              vsipj(1) = vsip_s(natj)
              if (norbj.gt.1) then
                vsipj(2) = vsip_p(natj)
                vsipj(3) = vsip_p(natj)
                vsipj(4) = vsip_p(natj)
              endif
!******************************
!  Loop over valid distances  *
!******************************
              do k = 1,nor
!
!  Check for self term
!
                if (i.eq.j.and.dist(k).lt.1.0d-2) cycle
!
!  Sqrt distances
!
                dist(k) = sqrt(dist(k))
!*********************************************************
!  Calculate overlap matrix contribution from this pair  *
!*********************************************************
                call overlap(nati,natj,dist(k),xtmp(k),ytmp(k),ztmp(k),Spair)
!
!  Compute phase factor
!
                rkr = xkc*xtmp(k) + ykc*ytmp(k) + zkc*ztmp(k)
                phase = dcmplx(cos(rkr),-sin(rkr))
!************************************
!  Add terms to Hamiltonian matrix  *
!************************************
                mm = norbsofari
                do m = 1,norbi
                  mm = mm + 1
                  nn = norbsofarj
                  do n = 1,norbj
                    nn = nn + 1
!
                    Hij = 0.5_dp*ofct*scale_eht*Spair(m,n)*(vsipi(m) + vsipj(n))
                    cHij = dcmplx(Hij,0.0_dp)*phase
                    cSij = dcmplx(Spair(m,n),0.0_dp)*phase
!
                    cHmat(nn,mm) = cHmat(nn,mm) + cHij
!
                    cSmat(nn,mm) = cSmat(nn,mm) + cSij
                  enddo
                enddo
              enddo
            enddo jloopk1
          enddo
        endif
!
!  Output matrices if desired
!
        if (lmain.and.ioproc) then
          write(ioout,'(/,''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  EHT K point : '',i4,3(1x,f8.5),2x,'' weight = '',f8.5)') nk,xehtk(nk),yehtk(nk),zehtk(nk),wehtk(nk)
          write(ioout,'(''--------------------------------------------------------------------------------'',/)')
        endif
        if (ldebug) then
          if (ioproc) then
            write(ioout,'(''  Extended Huckel Theory Hamiltonian matrix (eV) :'',/)')
            write(ioout,'(''  Real: '',/)')
          endif
          if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
            ntmp = norbtot
            ntag = 1
            allocate(dtmp(norbtot),stat=status)
            if (status/=0) call outofmemory('ehtd','dtmp')
            allocate(StatMPI(MPI_Status_Size),stat=status)
            if (status/=0) call outofmemory('ehtd','StatMPI')
          endif
          call mpbarrier
          do i = 1,norbtot
            iloc = norbg2l(i)
            if (lioproconly.and.norb2node(i).ne.0_i4) then
!
!  Post receive
!
              if (ioproc) then
                nnode = norb2node(i)
                call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_GULP,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (iloc.gt.0) then
                dtmp(1:norbtot) = dble(cHmat(1:norbtot,iloc))
!
!  Post send
!
                call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                               ntag,MPI_Comm_GULP,Request,MPIerror)
              endif
              if (ioproc.or.iloc.gt.0) then
                call MPI_WaitAll(1,Request,StatMPI,MPIerror)
              endif
              if (ioproc) then
!
!  Write on I/O node
!
                write(ioout,'(12f11.6)')(dtmp(j),j=1,norbtot)
              endif
            else
              if (iloc.gt.0) then
                write(ioout,'(12f11.6)')(dble(cHmat(j,iloc)),j=1,norbtot)
              endif
            endif
            call mpbarrier
          enddo
          if (ioproc) then
            write(ioout,'(/,''  Imaginary: '',/)')
          endif
          do i = 1,norbtot
            iloc = norbg2l(i)
            if (lioproconly.and.norb2node(i).ne.0_i4) then
!
!  Post receive
!
              if (ioproc) then
                nnode = norb2node(i)
                call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_GULP,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (iloc.gt.0) then
                dtmp(1:norbtot) = dimag(cHmat(1:norbtot,iloc))
!
!  Post send
!
                call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                               ntag,MPI_Comm_GULP,Request,MPIerror)
              endif
              if (ioproc.or.iloc.gt.0) then
                call MPI_WaitAll(1,Request,StatMPI,MPIerror)
              endif
              if (ioproc) then
!
!  Write on I/O node
!
                write(ioout,'(12f11.6)')(dtmp(j),j=1,norbtot)
              endif
            else
              if (iloc.gt.0) then
                write(ioout,'(12f11.6)')(dimag(cHmat(j,iloc)),j=1,norbtot)
              endif
            endif
            call mpbarrier
          enddo
          if (ioproc) then
            write(ioout,'(/,''  Extended Huckel Theory overlap matrix :'',/)')
            write(ioout,'(''  Real: '',/)')
          endif
          call mpbarrier
          do i = 1,norbtot
            iloc = norbg2l(i)
            if (lioproconly.and.norb2node(i).ne.0_i4) then
!             
!  Post receive
!  
              if (ioproc) then
                nnode = norb2node(i)
                call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_GULP,Request,MPIerror)
              endif
!               
!  Pass data to ioproc for writing
!           
              if (iloc.gt.0) then
                dtmp(1:norbtot) = dble(cSmat(1:norbtot,iloc))
!         
!  Post send
!           
                call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                               ntag,MPI_Comm_GULP,Request,MPIerror)
              endif
              if (ioproc.or.iloc.gt.0) then
                call MPI_WaitAll(1,Request,StatMPI,MPIerror)
              endif
              if (ioproc) then
!
!  Write on I/O node
!
                write(ioout,'(12f11.6)')(dtmp(j),j=1,norbtot)
              endif
            else
              if (iloc.gt.0) then
                write(ioout,'(12f11.6)')(dble(cSmat(j,iloc)),j=1,norbtot)
              endif
            endif
            call mpbarrier
          enddo
          if (ioproc) then
            write(ioout,'(/,''  Imaginary: '',/)')
          endif
          do i = 1,norbtot
            iloc = norbg2l(i)
            if (lioproconly.and.norb2node(i).ne.0_i4) then
!
!  Post receive
!
              if (ioproc) then
                nnode = norb2node(i)  
                call MPI_IRecv(dtmp,ntmp,MPI_double_precision,nnode, &
                               ntag,MPI_Comm_GULP,Request,MPIerror)
              endif
!
!  Pass data to ioproc for writing
!
              if (iloc.gt.0) then
                dtmp(1:norbtot) = dimag(cSmat(1:norbtot,iloc))
!
!  Post send
!
                call MPI_ISend(dtmp,ntmp,MPI_double_precision,0, &
                               ntag,MPI_Comm_GULP,Request,MPIerror)
              endif
              if (ioproc.or.iloc.gt.0) then
                call MPI_WaitAll(1,Request,StatMPI,MPIerror)
              endif
              if (ioproc) then
!
!  Write on I/O node
!
                write(ioout,'(12f11.6)')(dtmp(j),j=1,norbtot)
              endif
            else
              if (iloc.gt.0) then
                write(ioout,'(12f11.6)')(dimag(cSmat(j,iloc)),j=1,norbtot)
              endif
            endif
            call mpbarrier
          enddo
!
          if (lioproconly) then
            deallocate(StatMPI,stat=status)
            if (status/=0) call deallocate_error('ehtd','StatMPI')
            deallocate(dtmp,stat=status)
            if (status/=0) call deallocate_error('ehtd','dtmp')
          endif
        endif
!*****************************************************
!  Find eigenvalues and eigenvectors of Hamiltonian  *
!*****************************************************
        if (nk.eq.1) then
!
!  Initial dummy allocation of workspace for call to find sizes
!  NB: Size must be greater than 1 since other elements can be referenced. 10 is used here for safety.
!
          allocate(w1(10),stat=status)
          if (status/=0) call outofmemory('ehtd','w1')
          allocate(iw2(10),stat=status)
          if (status/=0) call outofmemory('ehtd','iw2')
          allocate(cw1(10),stat=status)
          if (status/=0) call outofmemory('ehtd','cw1')
!
!  Set parameters for eigensolver by calling to obtain memory
!
          na = norbtot
          nsize1 = -1
          nsize2 = -1
          nsize3 = -1
!
          call pzhegvx(1,'V','A','U',na,cHmat,1,1,idesc,cSmat,1,1,idesc,vl,vu,il,iu,etol,nefound,nvfound, &
                       eig(1,nk),orfac,cZmat,1,1,idesc,cw1,nsize1,w1,nsize2,iw2,nsize3,ifails,icluster,gap,info)
!
          if (info.ne.0) then
            write(cinfo,'(i4)') info
            call outerror('Error calling pzhegvx in EHTD: Info = '//cinfo,0_i4)
            call stopnow('ehtd')
          endif
!
          nsize1 = nint(dble(cw1(1)))
          nsize2 = nint(w1(1)) + 8_i4*norbtot   ! NB: The extra memory is added to handle clusters of eigenvalues
          nsize3 = iw2(1)
!
!  Deallocate workspace
!
          deallocate(cw1,stat=status)
          if (status/=0) call deallocate_error('ehtd','cw1')
          deallocate(iw2,stat=status)
          if (status/=0) call deallocate_error('ehtd','iw2')
          deallocate(w1,stat=status)
          if (status/=0) call deallocate_error('ehtd','w1')
        endif
!
!  Real allocation of workspace
!
        allocate(cw1(nsize1),stat=status)
        if (status/=0) call outofmemory('ehtd','cw1')
        allocate(w1(nsize2),stat=status)
        if (status/=0) call outofmemory('ehtd','w1')
        allocate(iw2(nsize3),stat=status)
        if (status/=0) call outofmemory('ehtd','iw2')
!
!  Call eigensolver
!
        call pzhegvx(1,'V','A','U',na,cHmat,1,1,idesc,cSmat,1,1,idesc,vl,vu,il,iu,etol,nefound,nvfound, &
                     eig(1,nk),orfac,cZmat,1,1,idesc,cw1,nsize1,w1,nsize2,iw2,nsize3,ifails,icluster,gap,info)
!
        if (info.ne.0) then
          write(cinfo,'(i4)') info
          call outerror('Error calling pzhegvx in EHTD: Info = '//cinfo,0_i4)
          call stopnow('ehtd')
        endif
!
        deallocate(iw2,stat=status)
        if (status/=0) call deallocate_error('ehtd','iw2')
        deallocate(w1,stat=status)
        if (status/=0) call deallocate_error('ehtd','w1')
        deallocate(cw1,stat=status)
        if (status/=0) call deallocate_error('ehtd','cw1')
!
!  Checks
!
        if (nefound.lt.nocc) then
          call outerror('Insufficient eigenvalues converged for EHT',0_i4)
          call stopnow('eht')
        endif
        if (nvfound.lt.nocc) then
          call outerror('Insufficient eigenvectors converged for EHT',0_i4)
          call stopnow('eht')
        endif
!
!  Eigenvalues
!
        if (lmain.and.ioproc) then
          write(ioout,'(/,''  Extended Huckel Theory Orbital Eigenvalues (eV): '',/)')
          do i = 1,norbtot
            write(ioout,'(2x,i5,2x,f12.6)') i,eig(i,nk)
          enddo
        endif
!********************************************************************************************************
!  Compute Mulliken populations while minimising communication by rebuilding overlap matrix on the fly  *
!********************************************************************************************************
        if (noccloc.gt.0) then
!
!  Loop over first atom
!
          norbsofari_next = 0
          do ii = 1,numat
            i = natomn2o(ii)
            xal = xclat(i)
            yal = yclat(i)
            zal = zclat(i)
            nati = nat(i)
            norbi = natorb(nati)
            norbsofari = norbsofari_next
            norbsofari_next = norbsofari + norbi
            oci = occuf(i)
!
!  Start of second atom loop
!
            norbsofarj_next = 0
            jloopk2: do jj = 1,ii
              j = natomn2o(jj)
              natj = nat(j)
              norbj = natorb(natj)
              norbsofarj = norbsofarj_next
              norbsofarj_next = norbsofarj + norbj
              ocj = occuf(j)
!
              xcrd = xclat(j) - xal
              ycrd = yclat(j) - yal
              zcrd = zclat(j) - zal
!***********************
!  Find valid vectors  *
!***********************
              nor = 0
              if (ndim.eq.3) then
                call rsearch3D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
              elseif (ndim.eq.2) then
                call rsearch2D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
              elseif (ndim.eq.1) then
                call rsearch1D(xcrd,ycrd,zcrd,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
              endif
!
              ofct = oci*ocj
              occ = 2.0_dp*wehtk(nk)*ofct
!
              if (i.eq.j) then
!***************
!  Self terms  *
!***************
                do no = 1,noccloc
                  mm = norbsofari
                  do m = 1,norbi
                    mm = mm + 1
                    cHij = cZmat(mm,no)
                    qf(i) = qf(i) + occ*dble(dconjg(cHij)*cHij)
                  enddo
                enddo
              endif
!
              if (nor.eq.0) cycle jloopk2
!******************************
!  Loop over valid distances  *
!******************************
              do k = 1,nor
!
!  Check for self term
!
                if (i.eq.j.and.dist(k).lt.1.0d-2) cycle
!
!  Sqrt distances
!
                dist(k) = sqrt(dist(k))
!*********************************************************
!  Calculate overlap matrix contribution from this pair  *
!*********************************************************
                call overlap(nati,natj,dist(k),xtmp(k),ytmp(k),ztmp(k),Spair)
!
!  Compute phase factor
!
                rkr = xkc*xtmp(k) + ykc*ytmp(k) + zkc*ztmp(k)
                phase = dcmplx(cos(rkr),-sin(rkr))
!*********************************
!  Add contributions to charges  *
!*********************************
                mm = norbsofari
                do m = 1,norbi
                  mm = mm + 1
                  nn = norbsofarj
                  do n = 1,norbj
                    nn = nn + 1
                    do no = 1,noccloc
!
!  Divide overlap density between pair of atoms
!
                      qtrm = occ*dble(dconjg(cZmat(nn,no))*cZmat(mm,no)*phase)*Spair(m,n)
                      qf(i) = qf(i) + qtrm
                      qf(j) = qf(j) + qtrm
                    enddo
                  enddo
                enddo
              enddo
            enddo jloopk2
          enddo
        endif
!
!  Add contribution to total energy from this k point
!
        do i = 1,nocc
          eeht = eeht + 2.0_dp*wehtk(nk)*eig(i,nk)
        enddo
      enddo
    endif
!
!  Global sum of charges while swapping atom order back from the local orbital-based one to the standard one
!
    allocate(dtmp(norbtot),stat=status)
    if (status/=0) call outofmemory('ehtd','dtmp')
!
    dtmp(1:numat) = qf(1:numat)
    call sumall(dtmp,qf,numat,'ehtd','qf')
!
    deallocate(dtmp,stat=status)
    if (status/=0) call deallocate_error('ehtd','dtmp')
!
!  Convert orbital densities to charges
!
    do i = 1,numat
      qf(i) = nvalenceelectron(nat(i)) - qf(i)
    enddo
!
!  Update asymmetric unit charges
!
    do i = 1,nasym
     qa(i) = qf(nrela2f(i))
    enddo
!*******************
!  Output results  *
!*******************
    if ((lmain.or.ldebug).and.ioproc) then
      write(ioout,'(//,''  Final Mulliken charges from Extended Huckel Theory :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Extended Huckel Theory Energy = '',f12.6,'' eV'')') eeht
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
!
!  Free local memory
!
    if (lgamma) then
      deallocate(eig,stat=status)
      if (status/=0) call deallocate_error('ehtd','eig')
      deallocate(Zmat,stat=status)
      if (status/=0) call deallocate_error('ehtd','Zmat')
      deallocate(Smat,stat=status)
      if (status/=0) call deallocate_error('ehtd','Smat')
      deallocate(Hmat,stat=status)
      if (status/=0) call deallocate_error('ehtd','Hmat')
    else
      deallocate(eig,stat=status)
      if (status/=0) call deallocate_error('ehtd','eig')
      deallocate(cZmat,stat=status)
      if (status/=0) call deallocate_error('ehtd','cZmat')
      deallocate(cSmat,stat=status)
      if (status/=0) call deallocate_error('ehtd','cSmat')
      deallocate(cHmat,stat=status)
      if (status/=0) call deallocate_error('ehtd','cHmat')
    endif
!
    deallocate(gap,stat=status)
    if (status/=0) call deallocate_error('ehtd','gap')
    deallocate(ifails,stat=status)
    if (status/=0) call deallocate_error('ehtd','ifails')
    deallocate(icluster,stat=status)
    if (status/=0) call deallocate_error('ehtd','icluster')
    deallocate(norb2node,stat=status)
    if (status/=0) call deallocate_error('ehtd','norb2node')
    deallocate(norbg2l,stat=status)
    if (status/=0) call deallocate_error('ehtd','norbg2l')
    deallocate(norbl2g,stat=status)
    if (status/=0) call deallocate_error('ehtd','norbl2g')
    deallocate(natoml2gold,stat=status)
    if (status/=0) call deallocate_error('ehtd','natoml2gold')
    deallocate(natoml2gnew,stat=status)
    if (status/=0) call deallocate_error('ehtd','natoml2gnew')
    deallocate(natomn2o,stat=status)
    if (status/=0) call deallocate_error('ehtd','natomn2o')
!
!  Timing
!
    time2 = g_cpu_time()
    teht = teht + time2 - time1
#ifdef TRACE
    call trace_out('ehtd')
#endif
#else           
    call outerror('ehtd called when not compiled with MPI',0_i4)
    call stopnow('ehtd')
#endif 
!
    return
    end
!
    subroutine eht_mol(nmoltodo,eeht,lmain)
!
!  Subroutine for Extended Huckel Theory calculation of charges of a molecule sub-system
!
!   6/23 Created from eht_mol
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
    use current
    use element
    use molecule
    use times
#ifdef TRACE
    use trace,         only : trace_in, trace_out
#endif
    implicit none
!
!  Passed variables
!
    integer(i4), intent(in)                      :: nmoltodo  ! Number of molecule whose charges are to be found
    logical,     intent(in)                      :: lmain     ! If true then this is a main call for output
    real(dp),    intent(inout)                   :: eeht
!
!  Local variables
!
    integer(i4)                                  :: i
    integer(i4)                                  :: ii
    integer(i4)                                  :: ifail
    integer(i4)                                  :: ilaenv
    integer(i4)                                  :: indi
    integer(i4)                                  :: indj
    integer(i4)                                  :: ix
    integer(i4)                                  :: iy
    integer(i4)                                  :: iz
    integer(i4)                                  :: j
    integer(i4)                                  :: jj
    integer(i4)                                  :: jx
    integer(i4)                                  :: jy
    integer(i4)                                  :: jz
    integer(i4)                                  :: m
    integer(i4)                                  :: mm
    integer(i4)                                  :: n
    integer(i4)                                  :: nb
    integer(i4)                                  :: nn
    integer(i4)                                  :: nati
    integer(i4)                                  :: natj
    integer(i4)                                  :: nelectron
    integer(i4)                                  :: no
    integer(i4)                                  :: nocc
    integer(i4)                                  :: norbi
    integer(i4)                                  :: norbj
    integer(i4)                                  :: norbsofari
    integer(i4)                                  :: norbsofari_next
    integer(i4)                                  :: norbsofarj
    integer(i4)                                  :: norbsofarj_next
    integer(i4)                                  :: norbtot
    integer(i4)                                  :: nsize1
    integer(i4)                                  :: status
    real(dp)                                     :: g_cpu_time
    real(dp),     dimension(:,:), allocatable    :: Hmat
    real(dp),     dimension(:,:), allocatable    :: Smat
    real(dp),     dimension(:,:), allocatable    :: eig
    real(dp)                                     :: Hij
    real(dp)                                     :: oci
    real(dp)                                     :: ocj
    real(dp)                                     :: ofct
    real(dp)                                     :: r2
    real(dp)                                     :: Spair(9,9)    ! Overlap matrix between a pair of atoms
    real(dp)                                     :: time1
    real(dp)                                     :: time2
    real(dp)                                     :: vsipi(4)
    real(dp)                                     :: vsipj(4)
    real(dp),     dimension(:),   allocatable    :: w1
    real(dp)                                     :: xal
    real(dp)                                     :: yal
    real(dp)                                     :: zal
    real(dp)                                     :: xcrd
    real(dp)                                     :: ycrd
    real(dp)                                     :: zcrd
#ifdef TRACE
    call trace_in('eht_mol')
#endif
    time1 = g_cpu_time()
!
!  Check molecule number
!
    if (nmoltodo.gt.nmol.or.nmoltodo.lt.1) then
      call outerror('Molecule number is out of range in EHT_mol',0_i4)
      call stopnow('eht_mol')
    endif
!   
!  Set up factorials - stored as g_factorial(n) in location(n+1)
!   
    fct(1) = 1.0_dp
    do i = 1,29
      fct(i+1) = fct(i)*dble(i)
    enddo
!
!  Find total number of orbitals
!
    norbtot = 0
    nelectron = 0
    do ii = 1,nmolatom(nmoltodo)
      i = nmollist(nmolptr(nmoltodo)+ii)
      nati = nat(i)
      norbi = natorb(nati)
      norbtot = norbtot + norbi
      nelectron = nelectron + nvalenceelectron(nati)
!
!  Check that no elements have more than s and p orbitals as d/f are not currently supported
!
      if (norbi.gt.4) then
        call outerror('Elements with d and/or f orbitals are not currently supported for EHT',0_i4)
        call stopnow('eht_mol')
      elseif (norbi.eq.0) then
        call outerror('Elements without orbitals are not currently supported for EHT',0_i4)
        call stopnow('eht_mol')
      endif
    enddo
!
!  Allocate Hamiltonian and overlap matrices
!
    allocate(Hmat(norbtot,norbtot),stat=status)
    if (status/=0) call outofmemory('eht_mol','Hmat')
    allocate(Smat(norbtot,norbtot),stat=status)
    if (status/=0) call outofmemory('eht_mol','Smat')
    allocate(eig(norbtot,1),stat=status)
    if (status/=0) call outofmemory('eht_mol','eig')
!
!  Initialise Hmat
!
    Hmat(1:norbtot,1:norbtot) = 0.0_dp
    Smat(1:norbtot,1:norbtot) = 0.0_dp
!**************************************
!  Build Extended Huckel Hamiltonian  *
!**************************************
!
!  Loop over first atom
!
    norbsofari_next = 0
    do ii = 1,nmolatom(nmoltodo)
      i = nmollist(nmolptr(nmoltodo)+ii)
      xal = xclat(i)
      yal = yclat(i)
      zal = zclat(i)
!
!  Correct coordinates of i for cell wrapping
!
      indi = nmolind(i)
      call mindtoijk(indi,ix,iy,iz)
      if (ix.ne.0) then
        xal = xal + ix*rv(1,1)
        yal = yal + ix*rv(2,1)
        zal = zal + ix*rv(3,1)
      endif
      if (iy.ne.0) then
        xal = xal + iy*rv(1,2)
        yal = yal + iy*rv(2,2)
        zal = zal + iy*rv(3,2)
      endif
      if (iz.ne.0) then
        xal = xal + iz*rv(1,3)
        yal = yal + iz*rv(2,3)
        zal = zal + iz*rv(3,3)
      endif
!
      nati = nat(i)
      norbi = natorb(nati)
      norbsofari = norbsofari_next
      norbsofari_next = norbsofari + norbi
      oci = occuf(i)
!
      vsipi(1) = vsip_s(nati)
      if (norbi.gt.1) then
        vsipi(2) = vsip_p(nati)
        vsipi(3) = vsip_p(nati)
        vsipi(4) = vsip_p(nati)
      endif
!
!  Start of second atom loop
!
      norbsofarj_next = 0
      jloop: do jj = 1,ii
        j = nmollist(nmolptr(nmoltodo)+jj)
        natj = nat(j)
        norbj = natorb(natj)
        norbsofarj = norbsofarj_next
        norbsofarj_next = norbsofarj + norbj
        ocj = occuf(j)
!
        xcrd = xclat(j) - xal
        ycrd = yclat(j) - yal
        zcrd = zclat(j) - zal
!
!  Correct coordinates of j for cell wrapping
!
        indj = nmolind(j)
        call mindtoijk(indj,jx,jy,jz)
        if (jx.ne.0) then
          xcrd = xcrd + jx*rv(1,1)
          ycrd = ycrd + jx*rv(2,1)
          zcrd = zcrd + jx*rv(3,1)
        endif
        if (jy.ne.0) then
          xcrd = xcrd + jy*rv(1,2)
          ycrd = ycrd + jy*rv(2,2)
          zcrd = zcrd + jy*rv(3,2)
        endif
        if (jz.ne.0) then
          xcrd = xcrd + jz*rv(1,3)
          ycrd = ycrd + jz*rv(2,3)
          zcrd = zcrd + jz*rv(3,3)
        endif
!
        if (i.eq.j) then
!***************
!  Self terms  *
!***************
          ofct = oci*oci
          mm = norbsofari
          do m = 1,norbi
            mm = mm + 1
            Hij = ofct*vsipi(m)
            Hmat(mm,mm) = Hmat(mm,mm) + Hij
            Smat(mm,mm) = Smat(mm,mm) + 1.0_dp
          enddo
          cycle jloop
        else
          r2 = xcrd**2 + ycrd**2 + zcrd**2
        endif
!
        ofct = oci*ocj
!
        vsipj(1) = vsip_s(natj)
        if (norbj.gt.1) then
          vsipj(2) = vsip_p(natj)
          vsipj(3) = vsip_p(natj)
          vsipj(4) = vsip_p(natj)
        endif
!
!  Sqrt distance
!
        r2 = sqrt(r2)
!*********************************************************
!  Calculate overlap matrix contribution from this pair  *
!*********************************************************
        call overlap(nati,natj,r2,xcrd,ycrd,zcrd,Spair)
!************************************
!  Add terms to Hamiltonian matrix  *
!************************************
        mm = norbsofari
        do m = 1,norbi
          mm = mm + 1
          nn = norbsofarj
          do n = 1,norbj
            nn = nn + 1
!
            Hij = 0.5_dp*ofct*scale_eht*Spair(m,n)*(vsipi(m) + vsipj(n))
!
            Hmat(nn,mm) = Hmat(nn,mm) + Hij
            Hmat(mm,nn) = Hmat(mm,nn) + Hij
!
            Smat(nn,mm) = Smat(nn,mm) + Spair(m,n)
            Smat(mm,nn) = Smat(mm,nn) + Spair(m,n)
          enddo
        enddo
      enddo jloop
    enddo
!
!  Output matrices if desired
!
    if (ldebug.and.ioproc) then
      write(ioout,'(/,''  Extended Huckel Theory Hamiltonian matrix (eV) for molecule '',i4,'' :'',/)') nmoltodo
      do i = 1,norbtot
        write(ioout,'(12f11.6)')(Hmat(j,i),j=1,norbtot)
      enddo
      write(ioout,'(/,''  Extended Huckel Theory overlap matrix for molecule '',i4,'' :'',/)') nmoltodo
      do i = 1,norbtot
        write(ioout,'(12f11.6)')(Smat(j,i),j=1,norbtot)
      enddo
    endif
!*****************************************************
!  Find eigenvalues and eigenvectors of Hamiltonian  *
!*****************************************************
!
!  Set parameters for eigensolver
!
    nb = ilaenv(1_i4,'DSYTRD','U',norbtot,-1_i4,-1_i4,-1_i4)
    nsize1 = max(1_i4,(nb+2_i4)*norbtot)
!
    allocate(w1(nsize1+10_i4),stat=status)
    if (status/=0) call outofmemory('eht_mol','w1')
!
!  Call eigensolver
!
    call dsygv(1_i4,'V','U',norbtot,Hmat,norbtot,Smat,norbtot,eig(1,1),w1,nsize1,ifail)
!
    deallocate(w1,stat=status)
    if (status/=0) call deallocate_error('eht_mol','w1')
!
!  Eigenvalues
!
    if (lmain.and.ioproc) then
      write(ioout,'(/,''  Extended Huckel Theory Orbital Eigenvalues (eV): '',/)')
      do i = 1,norbtot
        write(ioout,'(2x,i5,2x,f12.6)') i,eig(i,1)
      enddo
    endif
!
!  Compute total energy
!
    nocc = nelectron/2
    eeht = 0.0_dp
    do no = 1,nocc
      eeht = eeht + 2.0_dp*eig(no,1)
    enddo
    if (2*nocc.ne.nelectron) then
!
!  Odd electron case
!
      eeht = eeht + eig(nocc+1,1)
    endif
!*********************************
!  Compute Mulliken populations  *
!*********************************
    norbsofari_next = 0
    do ii = 1,nmolatom(nmoltodo)
      i = nmollist(nmolptr(nmoltodo)+ii)
      nati = nat(i)
      norbi = natorb(nati)
      norbsofari = norbsofari_next
      norbsofari_next = norbsofari + norbi
!
      qf(i) = 0.0_dp
!
      do no = 1,nocc
        mm = norbsofari
        do m = 1,norbi
          mm = mm + 1
!
!  Orbital population - NB Hmat now contains the eigenvectors
!
          qf(i) = qf(i) + 2.0_dp*Hmat(mm,no)**2
!
!  Overlap contribution - compute in 2 stages as only half of Smat remains after eigensolve
!
          do j = 1,mm-1
            qf(i) = qf(i) + 2.0_dp*Hmat(mm,no)*Hmat(j,no)*Smat(mm,j)
          enddo
          do j = mm+1,norbtot
            qf(i) = qf(i) + 2.0_dp*Hmat(mm,no)*Hmat(j,no)*Smat(j,mm)
          enddo
        enddo
      enddo
      if (2*nocc.ne.nelectron) then
        mm = norbsofari
        do m = 1,norbi
          mm = mm + 1
          qf(i) = qf(i) + Hmat(mm,nocc+1)**2
        enddo
        do j = 1,mm-1
          qf(i) = qf(i) + Hmat(mm,nocc+1)*Hmat(j,nocc+1)*Smat(mm,j)
        enddo
        do j = mm+1,norbtot
          qf(i) = qf(i) + Hmat(mm,nocc+1)*Hmat(j,nocc+1)*Smat(j,mm)
        enddo
      endif
    enddo
!
!  Compute total energy
!
    eeht = 0.0_dp
    do i = 1,nocc
      eeht = eeht + 2.0_dp*eig(i,1)
    enddo
!
!  Convert orbital densities to charges
!
    do ii = 1,nmolatom(nmoltodo)
      i = nmollist(nmolptr(nmoltodo)+ii)
      qf(i) = nvalenceelectron(nat(i)) - qf(i)
    enddo
!
!  Update asymmetric unit charges
!
    do ii = 1,nmolatom(nmoltodo)
      i = nmollist(nmolptr(nmoltodo)+ii)
      qa(nrelf2a(i)) = qf(i)
    enddo
!*******************
!  Output results  *
!*******************
    if ((lmain.or.ldebug).and.ioproc) then
      write(ioout,'(//,''  Final Mulliken charges from Extended Huckel Theory for molecule '',i4,'' :'',/)') nmoltodo
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do ii = 1,nmolatom(nmoltodo)
        i = nmollist(nmolptr(nmoltodo)+ii)
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,nat(i),qf(i)
      enddo
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
!
!  Free local memory
!
    deallocate(eig,stat=status)
    if (status/=0) call deallocate_error('eht_mol','eig')
    deallocate(Smat,stat=status)
    if (status/=0) call deallocate_error('eht_mol','Smat')
    deallocate(Hmat,stat=status)
    if (status/=0) call deallocate_error('eht_mol','Hmat')
!
!  Timing
!
    time2 = g_cpu_time()
    teht = teht + time2 - time1
#ifdef TRACE
    call trace_out('eht_mol')
#endif
!
    return
    end
!
    subroutine setehtkpt(ndim,kv)
!
!  Sets up k points according to shrinking factors for EHT solve
!
!  nehtk = total number of k points 
!  xehtk = fractional x component of k point
!  yehtk = fractional y component of k point
!  zehtk = fractional z component of k point
!
!  Julian Gale, CIC, Curtin University, May 2023
!
    implicit none
!
!  Passed variables
!
    integer(i4),   intent(in)   :: ndim      ! Dimensionality of system
    real(dp),      intent(in)   :: kv(3,3)   ! Reciprocal lattice vectors
!
!  Local variables
!
    integer(i4)                 :: i
    integer(i4)                 :: j
    integer(i4)                 :: k
    integer(i4)                 :: nx
    integer(i4)                 :: nx2
    integer(i4)                 :: ny
    integer(i4)                 :: nz
    real(dp)                    :: delx
    real(dp)                    :: dely
    real(dp)                    :: delz
    real(dp)                    :: ka
    real(dp)                    :: kb
    real(dp)                    :: kc
    real(dp)                    :: pi
    real(dp)                    :: rks
    real(dp)                    :: rx
    real(dp)                    :: ry
    real(dp)                    :: rz
    real(dp)                    :: sumwehtk
    real(dp)                    :: wx
    real(dp)                    :: xadd
    real(dp)                    :: yadd
    real(dp)                    :: zadd
!
!  Initial call to ensure arrays are dimensioned as per initial value of maxnkh
!
    call changemaxnehtk
!
    pi = 4.0_dp*atan(1.0_dp)
!
!  Compute reciprocal lattice vector lengths
!
    ka = sqrt(kv(1,1)**2 + kv(2,1)**2 + kv(3,1)**2)
    kb = sqrt(kv(1,2)**2 + kv(2,2)**2 + kv(3,2)**2)
    kc = sqrt(kv(1,3)**2 + kv(2,3)**2 + kv(3,3)**2)
!
!  Multiply ks by 2*pi to match the reciprocal lattice vectors
!
    rks = 0.5_dp/(kspace_eht*pi)
!
    nx = ka*rks + 1
    ny = kb*rks + 1
    nz = kc*rks + 1
!
    if (ndim.eq.2) then
      nz = 1
    elseif (ndim.eq.1) then
      ny = 1
      nz = 1
    endif
!
    delx = 1.0_dp/float(nx)
    dely = 1.0_dp/float(ny)
    delz = 1.0_dp/float(nz)
!
    xadd = 0.5_dp*delx
    yadd = 0.5_dp*dely
    zadd = 0.5_dp*delz
!
    if (ndim.le.2) then
      delz = 0.0_dp
      zadd = 0.0_dp
    endif
    if (ndim.eq.1) then
      dely = 0.0_dp
      yadd = 0.0_dp
    endif
!
!  Apply time reversal symmetry in x
!
    nx2 = nx/2
    if (2*nx2.ne.nx) nx2 = nx2 + 1
!*****************************************************
!  Generate k points across the full Brillouin zone  *
!*****************************************************
!
!  Check there is enough space to insert new K points
!
    if (nx*ny*nz.gt.maxnehtk) then
      maxnehtk = nx*ny*nz
      call changemaxnehtk
    endif
    nehtk = 0
    rx = - delx + xadd
    sumwehtk = 0.0_dp
    do i = 1,nx2
      rx = rx + delx
      if (abs(rx).lt.1.0d-8.or.abs(rx-0.5_dp).lt.1.0d-8) then
        wx = 0.5_dp
      else
        wx = 1.0_dp
      endif
      ry = - dely + yadd
      do j = 1,ny
        ry = ry + dely
        rz = - delz + zadd
        do k = 1,nz
          rz = rz + delz
          nehtk = nehtk + 1
!
!  Add in new point
!
          xehtk(nehtk) = rx
          yehtk(nehtk) = ry
          zehtk(nehtk) = rz
          wehtk(nehtk) = wx
          sumwehtk = sumwehtk + wx
        enddo
      enddo
    enddo
!
!  Normalise weights
!
    do i = 1,nehtk
      wehtk(i) = wehtk(i)/sumwehtk
    enddo
!
    return
    end subroutine setehtkpt
!
    subroutine changemaxnehtk
!
!  Changes the size of arrays that hold k points for EHT solve
!
!  Julian Gale, CIC, Curtin University, May 2023
!
    use reallocate
    implicit none
!
!  Local variables
!
    integer(i4)              :: ierror
!
    call realloc(wehtk,maxnehtk,ierror)
    if (ierror.ne.0) call outofmemory('changemaxnehtk','wehtk')
    call realloc(xehtk,maxnehtk,ierror)
    if (ierror.ne.0) call outofmemory('changemaxnehtk','xehtk')
    call realloc(yehtk,maxnehtk,ierror)
    if (ierror.ne.0) call outofmemory('changemaxnehtk','yehtk')
    call realloc(zehtk,maxnehtk,ierror)
    if (ierror.ne.0) call outofmemory('changemaxnehtk','zehtk')
!
    end subroutine changemaxnehtk
!
    subroutine overlap(ni,nj,r,rx,ry,rz,di)
!
!  Calculates the overlap integrals according to CNDO type approach.
!
    implicit none
!
!  Passed variables
!
    integer(i4), intent(in)    :: ni     ! Atomic number of atom i
    integer(i4), intent(in)    :: nj     ! Atomic number of atom j
    real(dp),    intent(in)    :: r
    real(dp),    intent(in)    :: rx
    real(dp),    intent(in)    :: ry
    real(dp),    intent(in)    :: rz
    real(dp),    intent(out)   :: di(9,9)
!
!  Local variables
!
    integer(i4)                :: i
    integer(i4)                :: ii
    integer(i4)                :: j
    integer(i4)                :: jj
    integer(i4)                :: kk
    integer(i4)                :: maxl
    integer(i4)                :: n11
    integer(i4)                :: n22
    integer(i4)                :: norbi
    integer(i4)                :: norbj
    real(dp)                   :: t(9,9)       ! Transformation for rotation of frame
    real(dp)                   :: td(9,9,3)    ! Cartesian derivatives of transformation
    real(dp)                   :: tds(9,9,6)   ! Strain derivatives of transformation
    real(dp)                   :: temp(9,9)
    real(dp)                   :: rau
    real(dp)                   :: rr
    real(dp)                   :: rra
    real(dp)                   :: rrx
    real(dp)                   :: rry
    real(dp)                   :: rrz
    real(dp)                   :: s
    real(dp)                   :: z1r
    real(dp)                   :: z2r
    real(dp)                   :: zetai(3)
    real(dp)                   :: zetaj(3)
!
!  Initialise local variables
!
    norbi = natorb(ni)
    norbj = natorb(nj)
    maxl = max0(norbital_l(norbi),norbital_l(norbj))
!
    zetai(1) = zeta_s(ni)
    zetai(2) = zeta_p(ni)
    zetai(3) = 0.0_dp      ! Zeta for d set to zero as d-block is currently not fully handled
    zetaj(1) = zeta_s(nj)
    zetaj(2) = zeta_p(nj)
    zetaj(3) = 0.0_dp      ! Zeta for d set to zero as d-block is currently not fully handled
!
!  Units of r need to be a.u.
!
    rr = 1.0_dp/r
    rrx = rx*rr
    rry = ry*rr
    rrz = rz*rr
    rau = r/0.529167_dp
    rra = rr*0.529167_dp
!
!  Calculate rotation matrices
!
    call getrotationmatrix(t,maxl,rau,rrx,rry,rrz,td,tds,.false.,.false.)
!
!  Loop through pairs of basis functions
!
    do i = 1,norbi
      ii = norbital_l(i) + 1
      z1r = zetai(ii)*rau
      n11 = 2*npqn(ni) + 1
      do j = 1,norbj
        jj = norbital_l(j) + 1
        if (norbital_m(i).ne.norbital_m(j)) then
          di(i,j) = 0.0_dp
        elseif (norbital_m(i).lt.0) then
          di(i,j) = di(i-1,j-1)
        else
          z2r = zetaj(jj)*rau
          n22 = 2*npqn(nj) + 1
          call ss(s,npqn(ni),norbital_l(i),norbital_m(i),npqn(nj),norbital_l(j),z1r,z2r,rra)
          di(i,j) = dsqrt((z1r**n11*z2r**n22)/(fct(n11)*fct(n22)))*(-1.0_dp)**(norbital_l(j)+norbital_m(j))*s
        endif
      enddo
    enddo
!
!  Rotate from diatomic to molecular coordinates
!
    do i = 1,norbi
      do j = 1,norbj
        temp(i,j) = 0.0_dp
        do kk = 1,norbj
          temp(i,j) = temp(i,j) + t(j,kk)*di(i,kk)
        enddo
      enddo
    enddo
!
    do i = 1,norbi
      do j = 1,norbj
        di(i,j) = 0.0_dp
        do kk = 1,norbi
          di(i,j) = di(i,j) + t(i,kk)*temp(kk,j)
        enddo
      enddo
    enddo
!
    return
    end subroutine overlap
!
    subroutine getrotationmatrix(t,maxl,r,x,y,z,td,tds,lgrad1,lstr)
!
!   Find the rotational transformation matrix and optionally its derivatives
!
!   This routine has now been modified to include calculation of the derivative
!   of the transformation matrix in order to allow correct calculation of
!   the analytical gradients.
!
!   Strain derivatives of transformation matrix added for s and p fns
!   NB: d-orbital terms need checking for case when x and y = 0
!
    implicit none
!
!  Passed variables
!
    integer(i4), intent(in)    :: maxl       ! Maximum angular momentum
    logical,     intent(in)    :: lgrad1
    logical,     intent(in)    :: lstr
    real(dp),    intent(out)   :: t(9,9)     ! Transformation matrix
    real(dp),    intent(out)   :: td(9,9,3)  ! Cartesian derivatives of transformation matrix
    real(dp),    intent(out)   :: tds(9,9,6) ! Strain derivatives of transformation matrix
    real(dp),    intent(in)    :: r          ! Distance in a.u.
    real(dp),    intent(in)    :: x
    real(dp),    intent(in)    :: y
    real(dp),    intent(in)    :: z
!
!  Local variables
!
    integer(i4)                :: i
    integer(i4)                :: j
    integer(i4)                :: k
    integer(i4)                :: m
    real(dp)                   :: cosp
    real(dp)                   :: cos2p
    real(dp)                   :: cost
    real(dp)                   :: cos2t
    real(dp)                   :: d
    real(dp)                   :: dd
    real(dp)                   :: ddd
    real(dp)                   :: dsinp(3)
    real(dp)                   :: dsin2p(3)
    real(dp)                   :: dsinps(6)
    real(dp)                   :: dsint(3)
    real(dp)                   :: dsin2t(3)
    real(dp)                   :: dsints(6)
    real(dp)                   :: dcosp(3)
    real(dp)                   :: dcos2p(3)
    real(dp)                   :: dcosps(6)
    real(dp)                   :: dcost(3)
    real(dp)                   :: dcos2t(3)
    real(dp)                   :: dcosts(6)
    real(dp)                   :: rr
    real(dp)                   :: rsint
    real(dp)                   :: rsint3
    real(dp)                   :: sinp
    real(dp)                   :: sin2p
    real(dp)                   :: sint
    real(dp)                   :: sin2t
    real(dp)                   :: sqrt3
!
    cost = z
    if (abs(cost).gt.1.0_dp) then
      cost = sign(1.0_dp,cost)
      sint = 0.0_dp
    else
      sint = dsqrt(1.0_dp-cost**2)
    endif
    rr = 1.0_dp/r
    if (sint.gt.1.0d-6) then
      cosp = x/sint
      sinp = y/sint
    else
      cosp = 1.0_dp
      sinp = 0.0_dp
    endif
    do i = 1,9
      do j = 1,9
        t(j,i) = 0.0_dp
      enddo
    enddo
    t(1,1) = 1.0_dp
!
!  For gradient calculations calculate angular derivatives in cartesian coordinates.
!
    if (lgrad1) then
      do i = 1,3
        do j = 1,9
          do k = 1,9
            td(k,j,i) = 0.0_dp
          enddo
        enddo
      enddo
      ddd = (x**2 + y**2)
      dd = sqrt(ddd)*r
      d = dd*ddd
      dcost(1) = - x*z*rr
      dcost(2) = - y*z*rr
!
!  Cartesian derivatives of sines and cosines
!
      if (abs(sint).le.1.0d-6) then
        dcosp(1) = 0.0_dp
        dcosp(2) = 0.0_dp
        dsinp(2) = 0.0_dp
        dsint(1) = 0.0_dp
        dsint(2) = 0.0_dp
      else
        dcosp(1) = y**2/d
        dcosp(2) = -x*y/d
        dsint(1) = -dcost(1)*z*r/dd
        dsint(2) = -dcost(2)*z*r/dd
        dsinp(2) = x**2/d
      endif
      dcosp(3) = 0.0_dp
      dsinp(1) = dcosp(2)
      dsinp(3) = 0.0_dp
      dcost(3) = ddd*rr
      dsint(3) = -dd*z*rr*rr
      if (lstr) then
        do i = 1,6
          do j = 1,4
            do k = 1,4
              tds(k,j,i) = 0.0_dp
            enddo
          enddo
        enddo
!
!  Strain derivatives of sines and cosines
!
        dcosts(1) = - z*x*x
        dcosts(2) = - z*y*y
        dcosts(3) = - z*z*z+z
        dcosts(4) = - z*y*z + 0.5_dp*y
        dcosts(5) = - z*x*z + 0.5_dp*x
        dcosts(6) = - z*x*y
        if (abs(sint).le.1.0d-6) then
          do i = 1,6
            dsints(i) = 0.0_dp
            dcosps(i) = 0.0_dp
            dsinps(i) = 0.0_dp
          enddo
        else
          rsint = 1.0_dp/sint
          dsints(1) = - sint*x*x+rsint*x*x
          dsints(2) = - sint*y*y+rsint*y*y
          dsints(3) = - sint*z*z
          dsints(4) = - sint*y*z+0.5d0*rsint*y*z
          dsints(5) = - sint*x*z+0.5d0*rsint*x*z
          dsints(6) = - sint*x*y+rsint*x*y
          rsint3 = rsint*rsint*rsint
          dsinps(1) = - y*rsint3*x*x
          dsinps(2) = - y*rsint3*y*y+rsint*y
          dsinps(3) = 0.0_dp
          dsinps(4) = - y*rsint3*y*z+rsint*z
          dsinps(5) = - y*rsint3*x*z
          dsinps(6) = - 2.0_dp*y*rsint3*x*y + rsint*x
          dcosps(1) = - x*rsint3*x*x + rsint*x
          dcosps(2) = - x*rsint3*y*y
          dcosps(3) = 0.0_dp
          dcosps(4) = - x*rsint3*y*z
          dcosps(5) = - x*rsint3*x*z + rsint*z
          dcosps(6) = - 2.0_dp*x*rsint3*x*y + rsint*y
          dsinps(4) = 0.5_dp*dsinps(4)
          dsinps(5) = 0.5_dp*dsinps(5)
          dsinps(6) = 0.5_dp*dsinps(6)
          dcosps(4) = 0.5_dp*dcosps(4)
          dcosps(5) = 0.5_dp*dcosps(5)
          dcosps(6) = 0.5_dp*dcosps(6)
        endif
      endif
    endif
!
    if (maxl.gt.1) then
      cos2t = cost**2 - sint**2
      sin2t = 2.0_dp*sint*cost
      cos2p = cosp**2 - sinp**2
      sin2p = 2.0_dp*sinp*cosp

      if (lgrad1) then
        do m = 1,3
          dsin2p(m) = 2.0_dp*(sinp*dcosp(m) + cosp*dsinp(m))
          dcos2p(m) = 4.0_dp*cosp*dcosp(m)
          dsin2t(m) = 2.0_dp*(sint*dcost(m) + cost*dsint(m))
          dcos2t(m) = 4.0_dp*cost*dcost(m)
        enddo
      endif
!
!  Transformation matrix elements for d functions
!
      sqrt3 = dsqrt(3.0_dp)
!      t(5,5) = (3.0d0*cost**2-1.0d0)/2.0d0
!      t(5,6) = -sqrt3*sin2t/2.0d0
!      t(5,8) = sqrt3*sint**2/2.0d0
!      t(6,5) = sqrt3*sin2t*cosp/2.0d0
!      t(6,6) = cos2t*cosp
!      t(6,7) = -cost*sinp
!      t(6,8) =-t(6,5)/sqrt3
!      t(6,9) = sint*sinp
!      t(7,5) = sqrt3*sin2t*sinp/2.0d0
!      t(7,6) = cos2t*sinp
!      t(7,7) = cost*cosp
!      t(7,8) = -t(7,5)/sqrt3
!      t(7,9) = -sint*cosp
!      t(8,5) = sqrt3*sint**2*cos2p/2.0d0
!      t(8,6) = sin2t*cos2p/2.0d0
!      t(8,7) = -sint*sin2p
!      t(8,8) = (1.0d0+cost**2)*cos2p/2.0d0
!      t(8,9) = -cost*sin2p
!      t(9,5) = sqrt3*sint**2*sin2p/2.0d0
!      t(9,6) = sin2t*sin2p/2.0d0
!      t(9,7) = sint*cos2p
!      t(9,8) = (1.0d0+cost**2)*sin2p/2.0d0
!      t(9,9) = cost*cos2p

!      if (lgrad1) then
!        do m = 1,3
!          td(5,5,m)=3.0d0*cost*dcost(m)
!          td(5,6,m)=-sqrt3*dsin2t(m)/2.0d0
!          td(5,8,m)=sqrt3*sint*dsint(m)
!          td(6,5,m)=sqrt3*(sin2t*dcosp(m)+cosp*dsin2t(m))/2.0d0
!          td(6,6,m)=cos2t*dcosp(m)+cosp*dcos2t(m)
!          td(6,7,m)=-(cost*dsinp(m)+sinp*dcost(m))
!          td(6,8,m)=-td(6,5,m)/sqrt3
!          td(6,9,m)=sint*dsinp(m)+sinp*dsint(m)
!          td(7,5,m)=sqrt3*(sin2t*dsinp(m)+sinp*dsin2t(m))/2.0d0
!          td(7,6,m)=cos2t*dsinp(m)+sinp*dcos2t(m)
!          td(7,7,m)=cost*dcosp(m)+cosp*dcost(m)
!          td(7,8,m)=-td(7,5,m)/sqrt3
!          td(7,9,m)=-(sint*dcosp(m)+cosp*dsint(m))
!          td(8,5,m)=sqrt3*(sint**2*dcos2p(m)/2.0d0+ cos2p*sint*dsint(m))
!          td(8,6,m)=(sin2t*dcos2p(m)+cos2p*dsin2t(m))/2.0d0
!          td(8,7,m)=-(sint*dsin2p(m)+sin2p*dsint(m))
!          td(8,8,m)=(1.0d0+cost**2)*dcos2p(m)/2.0d0+ cost*cos2p*dcost(m)
!          td(8,9,m)=-(cost*dsin2p(m)+sin2p*dcost(m))
!          td(9,5,m)=sqrt3*(sint**2*dsin2p(m)/2.0d0+ sin2p*sint*dsint(m))
!          td(9,6,m)=(sin2t*dsin2p(m)+sin2p*dsin2t(m))/2.0d0
!          td(9,7,m)=sint*dcos2p(m)+cos2p*dsint(m)
!          td(9,8,m)=(1.0d0+cost**2)*dsin2p(m)/2.0d0+ sin2p*cost*dcost(m)
!          td(9,9,m)=cost*dcos2p(m)+cos2p*dcost(m)
!        enddo
!      endif
    endif
    if (maxl.gt.0) then
!
!  Transformation matrix elements for p functions
!
      t(2,2) = cost*cosp
      t(2,3) = - sinp
      t(2,4) = sint*cosp
      t(3,2) = cost*sinp
      t(3,3) = cosp
      t(3,4) = sint*sinp
      t(4,2) = - sint
      t(4,4) = cost

      if (lgrad1) then
!
!  Cartesian derivatives
!
        do m = 1,3
          td(2,2,m) = cost*dcosp(m) + cosp*dcost(m)
          td(2,3,m) = - dsinp(m)
          td(2,4,m) = sint*dcosp(m) + cosp*dsint(m)
          td(3,2,m) = cost*dsinp(m) + sinp*dcost(m)
          td(3,3,m) = dcosp(m)
          td(3,4,m) = sint*dsinp(m) + sinp*dsint(m)
          td(4,2,m) = - dsint(m)
          td(4,4,m) = dcost(m)
        enddo
!
!  Handle terms which become inaccurate as products when x/y become small
!
        if (sint.le.1.0d-6) then
          td(2,4,1) = (1.0_dp - x*x)*rr
          td(3,4,2) = (1.0_dp - y*y)*rr
          td(4,2,1) = - (1.0_dp - x*x)*rr
          td(4,3,2) = (1.0_dp - y*y)*rr
        endif
!
!  Strain derivatives
!
        do m = 1,6
          tds(2,2,m) = dcosts(m)*cosp + cost*dcosps(m)
          tds(2,3,m) = - dsinps(m)
          tds(2,4,m) = dsints(m)*cosp + sint*dcosps(m)
          tds(3,2,m) = dcosts(m)*sinp + cost*dsinps(m)
          tds(3,3,m) = dcosps(m)
          tds(3,4,m) = dsints(m)*sinp + sint*dsinps(m)
          tds(4,2,m) = - dsints(m)
          tds(4,4,m) = dcosts(m)
        enddo
      endif
    endif
!
    return
    end subroutine getrotationmatrix
!
    subroutine ss(s,nn1,ll1,mm,nn2,ll2,alpha,beta,rr)
!
!  Calculates reduced overlap integral
!
!  rr = 1/r
!
    implicit none
!     
!  Passed variables
!       
    integer(i4), intent(in)    :: ll1        ! Angular momentum quantum number for orbital 1
    integer(i4), intent(in)    :: ll2        ! Angular momentum quantum number for orbital 2
    integer(i4), intent(in)    :: mm         ! Angular momentum projection quantum number
    integer(i4), intent(in)    :: nn1        ! Principal quantum number for orbital 1
    integer(i4), intent(in)    :: nn2        ! Principal quantum number for orbital 2
    real(dp),    intent(out)   :: s          ! Overlap integral
    real(dp),    intent(in)    :: rr         ! Distance in a.u.
    real(dp),    intent(in)    :: alpha      ! Zeta times distance for orbital 1
    real(dp),    intent(in)    :: beta       ! Zeta times distance for orbital 2
!     
!  Local variables
!     
    integer(i4)                :: i
    integer(i4)                :: ii
    integer(i4)                :: j
    integer(i4)                :: k
    integer(i4)                :: l1
    integer(i4)                :: l2
    integer(i4)                :: llim
    integer(i4)                :: m
    integer(i4)                :: n1
    integer(i4)                :: n2
    integer(i4)                :: nni1
    integer(i4)                :: ul2
    integer(i4)                :: ulim
    real(dp)                   :: a(21)
    real(dp)                   :: b(21)
    real(dp)                   :: ai1
    real(dp)                   :: ai11
    real(dp)                   :: bi
    real(dp)                   :: coff
    real(dp)                   :: d
    real(dp)                   :: p
    real(dp)                   :: pt
    real(dp)                   :: x
    real(dp)                   :: zeta1
    real(dp)                   :: zeta2
!
    n1 = nn1
    l1 = ll1
    m  = mm
    n2 = nn2
    l2 = ll2
    p  = 0.5_dp*(alpha + beta)
    pt = 0.5_dp*(alpha - beta)
    x = 0.0_dp
    zeta1 = alpha*rr
    zeta2 = beta*rr
    m = iabs(m)
!
!  Reverse quantum numbers if necessary
!
    if ((l2.lt.l1).or.((l2.eq.l1).and.(n2.lt.n1))) then
      k  = n1
      n1 = n2
      n2 = k
      k  = l1
      l1 = l2
      l2 = k
      pt = - pt
    endif
!
!  Trap for enormously long distances which would cause aintgs or bintgs to crash with an overflow
!
    if (p.gt.86.0_dp.or.pt.gt.86.0_dp) then
      s = 0.0_dp
      return
    endif
!
    k = mod((n1+n2-l1-l2),2_i4)
!
!  Find a and b integrals
!
    call aintgs(p,n1+n2+1,a)
    call bintgs(pt,n1+n2+1,b)
    llim = 1
    ulim = n1 + n2 + 1
    if ((l1.gt.0).or.(l2.gt.0)) then
!
!  Begin section used for overlaps involving non-s functions
!
      do i = llim,ulim
        ul2 = ulim/2 + mod(ulim,2_i4)
        ai1 = a(i)
        ai11 = a(i+1)
        do j = llim,ul2
          ii = 2*j + mod(k+i+1_i4,2_i4) - 1
          coff = coeff(n1,n2,l1,l2,m,i,ii)
          bi = b(ii)
          x = x + coff*ai1*bi
        enddo
      enddo
      d = (fct(m+2)/8.0_dp)**2*dsqrt(real(2*l1+1)*fct(l1-m+1)*real(2*l2+1)*fct(l2-m+1)/(4.0_dp*fct(l1+m+1)*fct(l2+m+1)))
      s = x*d
    else
!
!  Begin section used for overlap integrals involving s functions
!
      do i = llim,ulim
        nni1 = n1 + n2 - i + 2
        coff = coeffs(n1,n2,i-1)
        x = x + coff*a(i)*b(nni1)
      enddo
      s = 0.5_dp*x
    endif
    return
    end
!
    subroutine aintgs(x,k,a)
!
!  Aintgs forms the "a" integrals for the overlap calculation
!
    implicit none
!       
!  Passed variables
!       
    integer(i4), intent(in)    :: k
    real(dp),    intent(out)   :: a(*)
    real(dp),    intent(in)    :: x
!       
!  Local variables
!     
    integer(i4)                :: i
    real(dp)                   :: c
    real(dp)                   :: rx
!
    c = exp(-x)
    rx = 1.0_dp/x
    a(1) = c*rx
!
    do i = 1,k
      a(i+1) = (a(i)*dble(i) + c)*rx
    enddo
!
    return
    end subroutine aintgs
!
    subroutine bintgs(x,k,b)
!
!  Bintgs forms the "b" integrals for the overlap calculation
!
    implicit none
!
!  Passed variables
!
    integer(i4), intent(in)    :: k
    real(dp),    intent(out)   :: b(*)
    real(dp),    intent(in)    :: x
!   
!  Local variables
!   
    integer(i4)                :: i
    integer(i4)                :: io
    integer(i4)                :: last
    integer(i4)                :: m
    logical                    :: leqn1
    real(dp)                   :: absx
    real(dp)                   :: expmx
    real(dp)                   :: expx
    real(dp)                   :: rx
    real(dp)                   :: xf
    real(dp)                   :: y
!
    io = 0
    absx = abs(x)
    if (absx.gt.3.0_dp) then
      expx = exp(x)
      expmx = 1.0_dp/expx
      rx = 1.0_dp/x
      b(1) = (expx - expmx)*rx
      do i = 1,k
        b(i+1) = (dble(i)*b(i)+(-1.0_dp)**i*expx - expmx)*rx
      enddo
    elseif (absx.le.1.0d-6) then
      do i = io,k
        b(i+1) = dble(2*mod(i+1_i4,2_i4))/(dble(i)+1.0_dp)
      enddo
    else
      leqn1 = .true.
      if (absx.le.0.5_dp) then
        last = 6
      elseif (absx.le.1.0_dp) then
        if (k.le.5) then
          leqn1 = .false.
        else
          last = 7
        endif
      elseif (absx.le.2.0_dp) then
        if (k.le.7) then
          leqn1 = .false.
        else
          last = 12
        endif
      else
        if (k.le.10) then
          leqn1 = .false.
        else
          last = 15
        endif
      endif
      if (leqn1) then
        do i = io,k
          y = 0.0_dp
          do m = io,last
            xf = 1.0_dp
!
!  NB: Factorial index is m+1 here
!
            if (m.ne.0) xf = fct(m+1)
            y = y + (-x)**m*(2*mod(m+i+1_i4,2_i4))/(xf*(m+i+1))
          enddo
          b(i+1) = y
        enddo
      else
        expx = exp(x)
        expmx = 1.0_dp/expx
        rx = 1.0_dp/x
        b(1) = (expx-expmx)*rx
        do i = 1,k
          b(i+1) = (dble(i)*b(i)+(-1.0_dp)**i*expx-expmx)*rx
        enddo
      endif
    endif
!
    return
    end subroutine bintgs
!
    real(dp) function coeff(na,nb,la,lb,m,kl,ll)
    implicit none
!
!  Passed variables
!
    integer(i4), intent(in)    :: la
    integer(i4), intent(in)    :: lb
    integer(i4), intent(in)    :: ll
    integer(i4), intent(in)    :: kl
    integer(i4), intent(in)    :: m
    integer(i4), intent(in)    :: na
    integer(i4), intent(in)    :: nb
!
!  Local variables
!
    integer(i4)                :: a
    integer(i4)                :: au
    integer(i4)                :: b
    integer(i4)                :: ba
    integer(i4)                :: be
    integer(i4)                :: bl
    integer(i4)                :: bma
    integer(i4)                :: bmi
    integer(i4)                :: c
    integer(i4)                :: ca
    integer(i4)                :: cd
    integer(i4)                :: ce
    integer(i4)                :: cl
    integer(i4)                :: cma
    integer(i4)                :: d
    integer(i4)                :: dma
    integer(i4)                :: i
    integer(i4)                :: ia
    integer(i4)                :: iaa
    integer(i4)                :: ie
    integer(i4)                :: ijcd
    integer(i4)                :: il
    integer(i4)                :: ima
    integer(i4)                :: j
    integer(i4)                :: ja
    integer(i4)                :: jcd
    integer(i4)                :: je
    integer(i4)                :: jl
    integer(i4)                :: jma
    integer(i4)                :: k
    integer(i4)                :: kkll
    integer(i4)                :: l
    integer(i4)                :: l2
    integer(i4)                :: lam
    integer(i4)                :: lbm
    integer(i4)                :: namu
    integer(i4)                :: nbmv
    integer(i4)                :: nnll
    integer(i4)                :: p
    integer(i4)                :: u
    integer(i4)                :: ul
    integer(i4)                :: v
    integer(i4)                :: vl
    real(dp),            save  :: co(3,3,3)
    real(dp)                   :: f
!
    data co/8.0_dp,0.0_dp,-4.0_dp,0.0_dp,4.0_dp,0.0_dp,0.0_dp,0.0_dp,4.0_dp,0.0_dp,8.0_dp,0.0_dp,0.0_dp,0.0_dp, &
            12.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,12.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/
!
    coeff = 0.0_dp
    k = kl - 1
    l = ll - 1
    kkll = k + l
    nnll = na + nb - la - lb
    if (kkll.lt.nnll) return
    lam = la - m + 1
    lbm = lb - m + 1
    do ul = 1,lam
      vlloop: do vl = 1,lbm
        f = 0.0_dp
        u = ul - 1
        v = vl - 1
        p = na + nb - u - v
        iaa = mod(l+p-k,2_i4)
        if (iaa.ne.0) cycle vlloop
        namu = na - m - u
        nbmv = nb - m - v
        ima = min0(l,u)
        jma = min0(l,v)
        cma = min0(l,namu)
        dma = min0(l,nbmv)
        l2 = int(dble(l)/2.0_dp)
        bma = min0(l2,m)
        bmi = int(dble(l-ima-jma-cma-dma)/2.0_dp)
        bmi = max0(0_i4,bmi)
        ba = bmi + 1
        be = bma + 1
        if (be.lt.ba) cycle vlloop
        blloop: do bl = ba,be
          b = bl - 1
          ijcd = l - 2*b
          ie = min0(u,ijcd) + 1
          ia = ijcd - jma - cma - dma
          ia = max0(0_i4,ia) + 1
          if (ie.lt.ia) cycle blloop
          illoop: do il = ia,ie
            i = il - 1
            jcd = ijcd - i
            je = min0(v,jcd) + 1
            ja = jcd - cma - dma
            ja = max0(0_i4,ja) + 1
            if (je.lt.ja) cycle illoop
            jlloop: do jl=ja,je
              j = jl - 1
              cd = jcd - j
              ce = min0(namu,cd) + 1
              ca = cd - dma
              ca = max0(0_i4,ca) + 1
              if (ce.lt.ca) cycle jlloop
              clloop: do cl = ca,ce
                c = cl - 1
                d = cd - c
                au = mod(p-k+i+j-c-d,2_i4)
                if (au.ne.0) cycle clloop
                a = (p-k+i+j-c-d)/2
                if (a.lt.0) cycle clloop
                if (a.gt.m) cycle clloop
                f = f + binm(m,a)*binm(m,b)*binm(u,i)*binm(v,j)*binm(namu,c)*binm(nbmv,d)*(-1.0_dp)**(a+b+j+d)
              enddo clloop
            enddo jlloop
          enddo illoop
        enddo blloop
        coeff = coeff + co(la+1,m+1,u+1)*co(lb+1,m+1,v+1)*f
      enddo vlloop
    enddo
!
    return
    end function coeff
!
    real(dp) function coeffs(na,nb,k)
    implicit none
!   
!  Passed variables
!   
    integer(i4), intent(in)    :: na
    integer(i4), intent(in)    :: nb
    integer(i4), intent(in)    :: k
!   
!  Local variables 
!   
    integer(i4)                :: i
    integer(i4)                :: ia
    integer(i4)                :: ie
    integer(i4)                :: il
    integer(i4)                :: j
    integer(i4)                :: je
    integer(i4)                :: l
!
    coeffs = 0.0_dp
    l = na + nb - k
    ie = min0(l,na) + 1
    je = min0(l,nb)
    ia = l - je + 1
    do il = ia,ie
      i = il - 1
      j = l - i
      coeffs = coeffs + binm(na,i)*binm(nb,j)*(-1.0_dp)**j
    enddo
!
    return
    end function coeffs
!
    real(dp) function binm(n,i)
    implicit none
!
!  Passed variables
!
    integer(i4), intent(in)    :: n
    integer(i4), intent(in)    :: i
!
    binm = fct(n+1)/(fct(n-i+1)*fct(i+1))
!
    return
    end function binm

  end module m_eht
