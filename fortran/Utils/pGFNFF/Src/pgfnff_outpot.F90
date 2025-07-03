  subroutine pgfnff_outpot(numat,maxnbr,nnbr_bond,nbrno_bond)
!
!  Outputs GFNFF interatomic potential information 
!  NB: Should only be called if GFNFF is to be used and this is the I/O processor
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  use m_pgfnff_disp
  use m_pgfnff_nbr_lib
  use m_pgfnff_topo
  implicit none
!
!  Passed variables
!
  integer(i4),                  intent(in)     :: numat                    ! Number of atoms
  integer(i4),                  intent(in)     :: maxnbr                   ! Maximum number of neighbours
  integer(i4),                  intent(in)     :: nnbr_bond(numat)         ! Number of neighbours
  integer(i4),                  intent(in)     :: nbrno_bond(maxnbr,numat) ! Pointer to neighbours
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: ni
  real(dp)                                     :: r0
  real(dp)                                     :: v2
  real(dp)                                     :: v3
!
!***************************************
!  Output pGFNFF potential parameters  *
!***************************************
  write(ioout,'(/,''  pGFNFF potential parameters: '',/)')
!
!  Bond energy
!
  write(ioout,'(''  Bond energy: '',/)')
  write(ioout,'(''----------------------------------------------------------------------'')')
  write(ioout,'(''  I     J            r0 (initial)      Exponent        Pre-exponential'')')
  write(ioout,'(''                        (Ang)           (Ang^-2)            (eV) '')')
  write(ioout,'(''----------------------------------------------------------------------'')')
  do i = 1,numat
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
      r0 = vbnbr(1,ni,i)
      v2 = vbnbr(2,ni,i)
      v3 = vbnbr(3,ni,i)
      write(ioout,'(i4,1x,i4,5x,f8.5,1x,f12.6,1x,f12.6)')  i,j,r0,v2,v3
    enddo
  enddo
!
  return
  end
