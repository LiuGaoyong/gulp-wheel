  subroutine gfnff_cn(numat,cut2,logcn,dlogcndcn,d2logcndcn2,lgrad1,lgrad2)
!
!  Computes the coordination number for GFNFF
!
!  On entry : 
!
!  cut2            = cutoff squared
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!
!  On exit :
!
!  logcn           = log coordination numbers
!  dlogcndcn       = log coordination number derivatives w.r.t. coordination number if lgrad1 is true
!  d2logcndcn2     = log coordination number second derivatives w.r.t. coordination number if lgrad2 is true
!
!   8/20 Created
!  10/20 Strain added
!  10/20 Derivative with respect to cn now return rather than Cartesian/strain versions
!  11/20 Second derivatives added
!   2/21 Checking for non-zero terms added with storage of result in nnbr_cn / nbrno_cn
!   2/21 Derivative terms for coordination number computed here to avoid repetition
!
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, February 2021
!
  use datatypes
  use current,       only : nat
  use gulp_gfnff
  use m_gfnff_nbr
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: numat
  logical,     intent(in)                          :: lgrad1
  logical,     intent(in)                          :: lgrad2
  real(dp),    intent(in)                          :: cut2
  real(dp),    intent(out)                         :: logcn(*)
  real(dp),    intent(out)                         :: dlogcndcn(*)
  real(dp),    intent(out)                         :: d2logcndcn2(*)
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: j
  integer(i4)                                      :: ni
  real(dp)                                         :: cni
  real(dp)                                         :: cut
  real(dp)                                         :: dr
  real(dp)                                         :: erfCN
  real(dp)                                         :: exptrm
  real(dp)                                         :: logcnmax
  real(dp)                                         :: pi
  real(dp)                                         :: r0
  real(dp)                                         :: rij
!
  pi = 4.0_dp*atan(1.0_dp)
!
  cut = sqrt(cut2)
!
!  Check whether neighbour list allows for this cutoff
!
  if (cut.gt.cutnbr) call gfnff_getnbr(cut2,.false.)
!
!  Initialise number of non-zero terms in neighbour list
!
  nnbr_cn(1:numat) = 0
!************************************************************
!  Loop over pairs of atoms to compute coordination number  *
!************************************************************
  logcnmax = log(1.0_dp + exp(gfnff_cnmax))
  do i = 1,numat
    cni = 0.0_dp
!
!  Compute total coordination number for i, cni
!
    do ni = 1,nnbr(i)
      rij = rnbr(ni,i)
      if (rij.le.cut) then
        j = nbrno(ni,i)
        r0 = gfnff_rcov(nat(i)) + gfnff_rcov(nat(j))
        if (rij.lt.gfnff_cnerfcut_cn*r0) then
          dr = (rij - r0)/r0
          erfCN = 0.5_dp*(1.0_dp + erf(gfnff_kn_cn*dr))
          cni = cni + erfCN
!
!  Store pointer to non-zero term
!
          nnbr_cn(i) = nnbr_cn(i) + 1
          nbrno_cn(nnbr_cn(i),i) = ni
!
!  Store derivative terms if needed
!
          if (lgrad1) then
            d1cndr_cn(nnbr_cn(i),i) = gfnff_kn_cn*exp(-(gfnff_kn_cn*dr)**2)/(sqrt(pi)*r0*rij)
            if (lgrad2) then
              d2cndr_cn(nnbr_cn(i),i) = - d1cndr_cn(nnbr_cn(i),i)*(1.0_dp/rij + 2.0_dp*gfnff_kn_cn*gfnff_kn_cn*dr/r0)/rij
            endif
          endif
        endif
      endif
    enddo
!
!  Now create a function of this coordination number
!
    exptrm = exp(gfnff_cnmax - cni)
    logcn(i) = logcnmax - log(1.0_dp + exptrm)
    if (lgrad1) then
      dlogcndcn(i) = exptrm/(1.0_dp + exptrm)
      if (lgrad2) then
        d2logcndcn2(i) = dlogcndcn(i)*(dlogcndcn(i) - 1.0_dp)
      endif
    endif
  enddo
!
  return
  end
!
  subroutine gfnff_dcn(i,cut,logcn,sumdlogcn)
!
!  Computes the coordination number for GFNFF and the sum of derivatives
!
!  On entry : 
!
!  i               = atom number who coordination number is to be computed
!  cut2            = cutoff squared
!
!  On exit :
!
!  logcn           = log coordination number for i
!  sumdlogcn       = sum of log coordination number derivatives for i
!
!   8/20 Created
!
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, August 2020
!
  use datatypes
  use current,       only : nat
  use gulp_gfnff
  use m_gfnff_nbr
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i
  real(dp),    intent(in)                          :: cut
  real(dp),    intent(out)                         :: logcn
  real(dp),    intent(out)                         :: sumdlogcn
!
!  Local variables
!
  integer(i4)                                      :: j
  integer(i4)                                      :: ni
  real(dp)                                         :: cni
  real(dp)                                         :: dcnidrij
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcni(3)
  real(dp)                                         :: dtrm
  real(dp)                                         :: dr
  real(dp)                                         :: erfCN
  real(dp)                                         :: exptrm
  real(dp)                                         :: logcnmax
  real(dp)                                         :: pi
  real(dp)                                         :: r0
  real(dp)                                         :: rij
  real(dp)                                         :: sumdlogcni
!
  pi = 4.0_dp*atan(1.0_dp)
!************************************************************
!  Loop over pairs of atoms to compute coordination number  *
!************************************************************
  logcnmax = log(1.0_dp + exp(gfnff_cnmax))
  cni = 0.0_dp
!
!  Compute total coordination number for i, cni
!
  do ni = 1,nnbr(i)
    rij = rnbr(ni,i)
    if (rij.le.cut) then
      j = nbrno(ni,i)
      r0 = gfnff_rcov(nat(i)) + gfnff_rcov(nat(j))
      dr = (rij - r0)/r0
      erfCN = 0.5_dp*(1.0_dp + erf(gfnff_kn_cn*dr))
      cni = cni + erfCN
    endif
  enddo
!
!  Now create a function of this coordination number
!
  exptrm = exp(gfnff_cnmax - cni)
  logcn = logcnmax - log(1.0_dp + exptrm)
!
  sumdlogcn = 0.0_dp
  sumdlogcni = 0.0_dp
  dlogcnidcni = exptrm/(1.0_dp + exptrm)
!
!  Compute derivatives of total coordination number for i
!
  dlogcni(1:3) = 0.0_dp
  do ni = 1,nnbr(i)
    rij = rnbr(ni,i)
    if (rij.le.cut) then
      j = nbrno(ni,i)
      r0 = gfnff_rcov(nat(i)) + gfnff_rcov(nat(j))
      dr = (rij - r0)/r0
      dcnidrij = gfnff_kn_cn*exp(-(gfnff_kn_cn*dr)**2)/(sqrt(pi)*r0)
      dtrm = dlogcnidcni*dcnidrij
!
!  i-j contribution
!
      sumdlogcn  = sumdlogcn  + abs(dtrm)
!
!  Self term
!
      dcnidrij = dcnidrij/rij
      dlogcni(1) = dlogcni(1) - dlogcnidcni*dcnidrij*xnbr(ni,i)
      dlogcni(2) = dlogcni(2) - dlogcnidcni*dcnidrij*ynbr(ni,i)
      dlogcni(3) = dlogcni(3) - dlogcnidcni*dcnidrij*znbr(ni,i)
    endif
  enddo
  sumdlogcn  = sumdlogcn  + sqrt(dlogcni(1)**2+dlogcni(2)**2+dlogcni(3)**2)
!
  return
  end
