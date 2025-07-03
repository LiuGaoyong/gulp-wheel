  subroutine reaxFF_bomds(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                          d1BOp,d1BOp_pi,d1BOp_pipi,BOij_s,BOij_pi,BOij_pipi,d1BOi_s,d1BOi_pi,d1BOi_pipi, &
                          nonzero,nonzeroi,nonzeroptr,lgrad1)
!
!  Calculates the bond order terms and first derivatives
!  Sparse derivative version.
!
!  On entry : 
!
!  lgrad1              = if .true. calculate the first derivatives
!
!  On exit :
!
!  nonzero             = number of nonzero derivatives
!  nonzeroi            = number of nonzero derivatives for neighbours of i only
!  nonzeroptr(nonzero) = pointer to real index of nonzero terms
!  BOij_s              = corrected sigma bond order
!  BOij_pi             = corrected pi bond order
!  BOij_pipi           = corrected pi-pi bond order
!  + corresponding derivatives as required
!
!   4/23 Created from reaxff_bos
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
!  Julian Gale, CIC, Curtin University, April 2023
!
  use datatypes
  use current
  use iochannels
  use neighbours
  use reaxFFdata
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                        :: i
  integer(i4), intent(in)                        :: j
  integer(i4), intent(in)                        :: nspeci
  integer(i4), intent(in)                        :: nspecj
  integer(i4), intent(in)                        :: ni
  integer(i4), intent(in)                        :: nj
  integer(i4), intent(in)                        :: mneigh
  integer(i4), intent(in)                        :: nneigh(*)
  integer(i4), intent(out)                       :: nonzero
  integer(i4), intent(out)                       :: nonzeroi
  integer(i4), intent(out)                       :: nonzeroptr(*)
  real(dp),    intent(in)                        :: deltapi
  real(dp),    intent(in)                        :: deltapj
  real(dp),    intent(in)                        :: BOp(maxneigh,*)
  real(dp),    intent(in)                        :: BOp_pi(maxneigh,*)
  real(dp),    intent(in)                        :: BOp_pipi(maxneigh,*)
  real(dp),    intent(in)                        :: d1BOp(maxneigh,*)
  real(dp),    intent(in)                        :: d1BOp_pi(maxneigh,*)
  real(dp),    intent(in)                        :: d1BOp_pipi(maxneigh,*)
  real(dp),    intent(out)                       :: BOij_s
  real(dp),    intent(out)                       :: BOij_pi
  real(dp),    intent(out)                       :: BOij_pipi
  real(dp),    intent(out)                       :: d1BOi_s(*)
  real(dp),    intent(out)                       :: d1BOi_pi(*)
  real(dp),    intent(out)                       :: d1BOi_pipi(*)
  logical,     intent(in)                        :: lgrad1
!
!  Local variables
!
  integer(i4)                                    :: ind
  integer(i4)                                    :: k
  integer(i4)                                    :: nk
  integer(i4)                                    :: nneighij
  integer(i4)                                    :: nneighij2
  integer(i4)                                    :: status
  logical                                        :: lf1squared
  real(dp)                                       :: BOij
  real(dp)                                       :: BOpij
  real(dp)                                       :: BOpij_s
  real(dp)                                       :: BOpij_pi
  real(dp),    dimension(:),   allocatable, save :: df1dr
  real(dp),    dimension(:),   allocatable, save :: df45dr
  real(dp)                                       :: dBOpijdr
  real(dp)                                       :: df1ddpi
  real(dp)                                       :: df1ddpj
  real(dp)                                       :: df4ddpi
  real(dp)                                       :: df4dbo
  real(dp)                                       :: df5ddpj
  real(dp)                                       :: df5dbo
  real(dp)                                       :: d2f1ddpi2
  real(dp)                                       :: d2f1ddpiddpj
  real(dp)                                       :: d2f1ddpj2
  real(dp)                                       :: d2f4ddpi2
  real(dp)                                       :: d2f4ddpidbo
  real(dp)                                       :: d2f4dbo2
  real(dp)                                       :: d2f5ddpj2
  real(dp)                                       :: d2f5ddpjdbo
  real(dp)                                       :: d2f5dbo2
  real(dp)                                       :: f1
  real(dp)                                       :: f4
  real(dp)                                       :: f5
  real(dp)                                       :: fs
  real(dp)                                       :: fpi
  real(dp)                                       :: fpipi
#ifdef TRACE
  call trace_in('reaxFF_bomds')
#endif
!
!  Set flag for f1squared or not
!
!  Original version didn't have square, but new one does
!
  lf1squared = .true.
!
!  Sizes for arrays
!
  nneighij = nneigh(i) + nneigh(j)
  nneighij2 = nneighij*(nneighij+1)/2
!
!  Branch according to which corrections need to be applied to the bond order
!
  if (nspeci.gt.nspecj) then
    ind = nspeci*(nspeci - 1)/2 + nspecj
  else
    ind = nspecj*(nspecj - 1)/2 + nspeci
  endif
  if (lreaxFFbocorrect(1,ind).and.lreaxFFbocorrect(2,ind)) then
!***********************************
!  Both corrections to be applied  *
!***********************************
!
!  Allocate local memory
!
    if (lgrad1) then
      allocate(df1dr(nneighij),stat=status)
      if (status/=0) call outofmemory('reaxFF_bo','df1dr')
      allocate(df45dr(nneighij),stat=status)
      if (status/=0) call outofmemory('reaxFF_bo','df45dr')
    endif
!
!  Calculate Bij and Bji - loop over all other neighbours
!
    if (lgrad1) then
      d1BOi_s(1:nneighij) = 0.0_dp
      d1BOi_pi(1:nneighij) = 0.0_dp
      d1BOi_pipi(1:nneighij) = 0.0_dp
      df1dr(1:nneighij) = 0.0_dp
      df45dr(1:nneighij) = 0.0_dp
    endif
!
!  Set local variable for bond order prime of currrent bond
!
    BOpij_s  = BOp(ni,i)
    BOpij_pi = BOp_pi(ni,i) + BOp_pipi(ni,i)
    BOpij = BOpij_s + BOpij_pi
!
!  Compute functions of delta'
!
    call reaxFF_f1(nspeci,nspecj,deltapi,deltapj,f1,df1ddpi,df1ddpj,d2f1ddpi2,d2f1ddpiddpj,d2f1ddpj2,lgrad1,.false.)
    call reaxFF_f45(nspeci,nspecj,BOpij,deltapi,f4,df4ddpi,df4dbo,d2f4ddpi2,d2f4ddpidbo,d2f4dbo2,lgrad1,.false.)
    call reaxFF_f45(nspecj,nspeci,BOpij,deltapj,f5,df5ddpj,df5dbo,d2f5ddpj2,d2f5ddpjdbo,d2f5dbo2,lgrad1,.false.)
!
!  Calculate corrected bond orders for i-j
!
    BOij = BOpij*f1*f4*f5
    if (lf1squared) then
      BOij_pi   = BOp_pi(ni,i)*f1*f1*f4*f5
      BOij_pipi = BOp_pipi(ni,i)*f1*f1*f4*f5
    else
      BOij_pi   = BOp_pi(ni,i)*f1*f4*f5
      BOij_pipi = BOp_pipi(ni,i)*f1*f4*f5
    endif
    BOij_s = BOij - BOij_pi - BOij_pipi
!
!  Compute derivatives of combined bond order term
!
    if (lgrad1) then
      if (lf1squared) then
        BOpij_s = BOpij_s + (1.0_dp - f1)*BOpij_pi
      endif
!
!  First derivatives of f1
!
      nonzero = 0
      do k = 1,nneigh(i)
        nonzero = nonzero + 1
        nonzeroptr(nonzero) = k
        dBOpijdr = (d1BOp(k,i) + d1BOp_pi(k,i) + d1BOp_pipi(k,i))
        df1dr(nonzero) = df1dr(nonzero) + df1ddpi*dBOpijdr
      enddo
      nonzeroi = nonzero
      do k = 1,nneigh(j)
        nonzero = nonzero + 1
        nonzeroptr(nonzero) = k
        dBOpijdr = (d1BOp(k,j) + d1BOp_pi(k,j) + d1BOp_pipi(k,j))
        df1dr(nonzero) = df1dr(nonzero) + df1ddpj*dBOpijdr
      enddo
!
!  First derivatives of f4 and f5
!
      dBOpijdr = (d1BOp(ni,i) + d1BOp_pi(ni,i) + d1BOp_pipi(ni,i))
      df45dr(ni) = df45dr(ni) + f5*df4dbo*dBOpijdr
      df45dr(nneigh(i)+nj) = df45dr(nneigh(i)+nj) + f4*df5dbo*dBOpijdr
!
      do k = 1,nneigh(i)
        dBOpijdr = (d1BOp(k,i) + d1BOp_pi(k,i) + d1BOp_pipi(k,i))
        df45dr(k) = df45dr(k) + f5*df4ddpi*dBOpijdr
      enddo
!
      do k = 1,nneigh(j)
        dBOpijdr = (d1BOp(k,j) + d1BOp_pi(k,j) + d1BOp_pipi(k,j))
        df45dr(nneigh(i)+k) = df45dr(nneigh(i)+k) + f4*df5ddpj*dBOpijdr
      enddo
!
!  Derivatives of BO' 
!
      if (lf1squared) then
        d1BOi_s(ni)    = d1BOi_s(ni)    + (d1BOp(ni,i) + (1.0_dp - f1)*(d1BOp_pi(ni,i) + d1BOp_pipi(ni,i)))*f1*f4*f5
        d1BOi_pi(ni)   = d1BOi_pi(ni)   + d1BOp_pi(ni,i)*f1*f1*f4*f5
        d1BOi_pipi(ni) = d1BOi_pipi(ni) + d1BOp_pipi(ni,i)*f1*f1*f4*f5
!
!  Derivatives of f1, f4 & f5 - note that derivatives of BOp are the same as those of deltap with respect to the same bond
!
        fs = f4*f5*(BOpij_s - f1*BOpij_pi)
        fpi = 2.0_dp*f4*f5*f1*BOp_pi(ni,i)
        fpipi = 2.0_dp*f4*f5*f1*BOp_pipi(ni,i)
!
        do nk = 1,nonzero
          d1BOi_s(nk)    = d1BOi_s(nk)    + fs*df1dr(nk)
          d1BOi_pi(nk)   = d1BOi_pi(nk)   + fpi*df1dr(nk)
          d1BOi_pipi(nk) = d1BOi_pipi(nk) + fpipi*df1dr(nk)
!
          d1BOi_s(nk)    = d1BOi_s(nk)    + BOpij_s*f1*df45dr(nk)
          d1BOi_pi(nk)   = d1BOi_pi(nk)   + BOp_pi(ni,i)*f1*f1*df45dr(nk)
          d1BOi_pipi(nk) = d1BOi_pipi(nk) + BOp_pipi(ni,i)*f1*f1*df45dr(nk)
        enddo
      else
!
!  Old ReaxFF scheme - no f1 squared
!
        d1BOi_s(ni)    = d1BOi_s(ni)    + d1BOp(ni,i)*f1*f4*f5
        d1BOi_pi(ni)   = d1BOi_pi(ni)   + d1BOp_pi(ni,i)*f1*f4*f5
        d1BOi_pipi(ni) = d1BOi_pipi(ni) + d1BOp_pipi(ni,i)*f1*f4*f5
!
!  Derivatives of f1, f4 & f5 - note that derivatives of BOp are the same as those of deltap with respect to the same bond
!
        do nk = 1,nonzero
          d1BOi_s(nk)    = d1BOi_s(nk)    + BOp(ni,i)*f4*f5*df1dr(nk)
          d1BOi_pi(nk)   = d1BOi_pi(nk)   + BOp_pi(ni,i)*f4*f5*df1dr(nk)
          d1BOi_pipi(nk) = d1BOi_pipi(nk) + BOp_pipi(ni,i)*f4*f5*df1dr(nk)
!
          d1BOi_s(nk)    = d1BOi_s(nk)    + BOp(ni,i)*f1*df45dr(nk)
          d1BOi_pi(nk)   = d1BOi_pi(nk)   + BOp_pi(ni,i)*f1*df45dr(nk)
          d1BOi_pipi(nk) = d1BOi_pipi(nk) + BOp_pipi(ni,i)*f1*df45dr(nk)
        enddo
      endif
    endif
!
!  Free local memory
!
    if (lgrad1) then
      deallocate(df45dr,stat=status)
      if (status/=0) call deallocate_error('reaxFF_bo','df45dr')
      deallocate(df1dr,stat=status)
      if (status/=0) call deallocate_error('reaxFF_bo','df1dr')
    endif
  elseif (lreaxFFbocorrect(1,ind).and..not.lreaxFFbocorrect(2,ind)) then
!***************************************************
!  Only overcoordination correction to be applied  *
!***************************************************
!
!  Allocate local memory
!
    if (lgrad1) then
      allocate(df1dr(nneighij),stat=status)
      if (status/=0) call outofmemory('reaxFF_bo','df1dr')
    endif
!
!  Calculate Bij and Bji - loop over all other neighbours
!
    if (lgrad1) then
      d1BOi_s(1:nneighij) = 0.0_dp
      d1BOi_pi(1:nneighij) = 0.0_dp
      d1BOi_pipi(1:nneighij) = 0.0_dp
      df1dr(1:nneighij) = 0.0_dp
    endif
!
!  Set local variable for bond order prime of currrent bond
!
    BOpij_s  = BOp(ni,i)
    BOpij_pi = BOp_pi(ni,i) + BOp_pipi(ni,i)
    BOpij = BOpij_s + BOpij_pi
!
!  Compute functions of delta'
!
    call reaxFF_f1(nspeci,nspecj,deltapi,deltapj,f1,df1ddpi,df1ddpj,d2f1ddpi2,d2f1ddpiddpj,d2f1ddpj2,lgrad1,.false.)
!
!  Calculate corrected bond orders for i-j
!
    BOij = BOpij*f1
    if (lf1squared) then
      BOij_pi   = BOp_pi(ni,i)*f1*f1
      BOij_pipi = BOp_pipi(ni,i)*f1*f1
    else
      BOij_pi   = BOp_pi(ni,i)*f1
      BOij_pipi = BOp_pipi(ni,i)*f1
    endif
    BOij_s = BOij - BOij_pi - BOij_pipi
!
!  Compute derivatives of combined bond order term
!
    if (lgrad1) then
      if (lf1squared) then
        BOpij_s = BOpij_s + (1.0_dp - f1)*BOpij_pi
      endif
!
!  First derivatives of f1 
!
      nonzero = 0
      do k = 1,nneigh(i)
        dBOpijdr = (d1BOp(k,i) + d1BOp_pi(k,i) + d1BOp_pipi(k,i))
        nonzero = nonzero + 1
        nonzeroptr(nonzero) = k
        df1dr(nonzero) = df1dr(nonzero) + df1ddpi*dBOpijdr
      enddo
      nonzeroi = nonzero
      do k = 1,nneigh(j)
        dBOpijdr = (d1BOp(k,j) + d1BOp_pi(k,j) + d1BOp_pipi(k,j))
        nonzero = nonzero + 1
        nonzeroptr(nonzero) = k
        df1dr(nonzero) = df1dr(nonzero) + df1ddpj*dBOpijdr
      enddo
!
!  Derivatives of BO' 
!
      if (lf1squared) then
        d1BOi_s(ni)    = d1BOi_s(ni)    + (d1BOp(ni,i) + (1.0_dp - f1)*(d1BOp_pi(ni,i) + d1BOp_pipi(ni,i)))*f1
        d1BOi_pi(ni)   = d1BOi_pi(ni)   + d1BOp_pi(ni,i)*f1*f1
        d1BOi_pipi(ni) = d1BOi_pipi(ni) + d1BOp_pipi(ni,i)*f1*f1
!
!  Derivatives of f1 - note that derivatives of BOp are the same as those of deltap with respect to the same bond
!
        do nk = 1,nonzero
          d1BOi_s(nk)    = d1BOi_s(nk)    + (BOpij_s - f1*BOpij_pi)*df1dr(nk)
          d1BOi_pi(nk)   = d1BOi_pi(nk)   + 2.0_dp*f1*BOp_pi(ni,i)*df1dr(nk)
          d1BOi_pipi(nk) = d1BOi_pipi(nk) + 2.0_dp*f1*BOp_pipi(ni,i)*df1dr(nk)
        enddo
      else
!
!  Old ReaxFF scheme - no f1 squared
!
        d1BOi_s(ni)    = d1BOi_s(ni)    + d1BOp(ni,i)*f1
        d1BOi_pi(ni)   = d1BOi_pi(ni)   + d1BOp_pi(ni,i)*f1
        d1BOi_pipi(ni) = d1BOi_pipi(ni) + d1BOp_pipi(ni,i)*f1
!
!  Derivatives of f1 - note that derivatives of BOp are the same as those of deltap with respect to the same bond
!
        do nk = 1,nonzero
          d1BOi_s(nk)    = d1BOi_s(nk)    + BOp(ni,i)*df1dr(nk)
          d1BOi_pi(nk)   = d1BOi_pi(nk)   + BOp_pi(ni,i)*df1dr(nk)
          d1BOi_pipi(nk) = d1BOi_pipi(nk) + BOp_pipi(ni,i)*df1dr(nk)
        enddo
      endif
    endif
!
!  Free local memory
!
    if (lgrad1) then
      deallocate(df1dr,stat=status)
      if (status/=0) call deallocate_error('reaxFF_bo','df1dr')
    endif
  elseif (.not.lreaxFFbocorrect(1,ind).and.lreaxFFbocorrect(2,ind)) then
!**************************************
!  Only 1-3 correction to be applied  *
!**************************************
!
!  Allocate local memory
!
    if (lgrad1) then
      allocate(df45dr(nneighij),stat=status)
      if (status/=0) call outofmemory('reaxFF_bo','df45dr')
    endif
!
!  Calculate Bij and Bji - loop over all other neighbours
!
    if (lgrad1) then
      d1BOi_s(1:nneighij) = 0.0_dp
      d1BOi_pi(1:nneighij) = 0.0_dp
      d1BOi_pipi(1:nneighij) = 0.0_dp
      df45dr(1:nneighij) = 0.0_dp
    endif
!
!  Set local variable for bond order prime of currrent bond
!
    BOpij_s  = BOp(ni,i)
    BOpij_pi = BOp_pi(ni,i) + BOp_pipi(ni,i)
    BOpij = BOpij_s + BOpij_pi
!
!  Compute functions of delta'
!
    call reaxFF_f45(nspeci,nspecj,BOpij,deltapi,f4,df4ddpi,df4dbo,d2f4ddpi2,d2f4ddpidbo,d2f4dbo2,lgrad1,.false.)
    call reaxFF_f45(nspecj,nspeci,BOpij,deltapj,f5,df5ddpj,df5dbo,d2f5ddpj2,d2f5ddpjdbo,d2f5dbo2,lgrad1,.false.)
!
!  Calculate corrected bond orders for i-j
!
    BOij = BOpij*f4*f5
    BOij_pi   = BOp_pi(ni,i)*f4*f5
    BOij_pipi = BOp_pipi(ni,i)*f4*f5
    BOij_s = BOij - BOij_pi - BOij_pipi
!
!  Compute derivatives of combined bond order term
!
    if (lgrad1) then
!
!  First derivatives of f4 and f5
!
      dBOpijdr = (d1BOp(ni,i) + d1BOp_pi(ni,i) + d1BOp_pipi(ni,i))
      df45dr(ni) = df45dr(ni) + f5*df4dbo*dBOpijdr
      df45dr(nneigh(i)+nj) = df45dr(nneigh(i)+nj) + f4*df5dbo*dBOpijdr
!
      nonzero = 0
      do k = 1,nneigh(i)
        dBOpijdr = (d1BOp(k,i) + d1BOp_pi(k,i) + d1BOp_pipi(k,i))
        nonzero = nonzero + 1
        nonzeroptr(nonzero) = k
        df45dr(nonzero) = df45dr(nonzero) + f5*df4ddpi*dBOpijdr
      enddo
      nonzeroi = nonzero
!
      do k = 1,nneigh(j)
        dBOpijdr = (d1BOp(k,j) + d1BOp_pi(k,j) + d1BOp_pipi(k,j))
        nonzero = nonzero + 1
        nonzeroptr(nonzero) = k
        df45dr(nonzero) = df45dr(nonzero) + f4*df5ddpj*dBOpijdr
      enddo
!
!  Derivatives of BO' 
!
      d1BOi_s(ni)    = d1BOi_s(ni)    + d1BOp(ni,i)*f4*f5
      d1BOi_pi(ni)   = d1BOi_pi(ni)   + d1BOp_pi(ni,i)*f4*f5
      d1BOi_pipi(ni) = d1BOi_pipi(ni) + d1BOp_pipi(ni,i)*f4*f5
!
!  Derivatives of f4 & f5 - note that derivatives of BOp are the same as those of deltap with respect to the same bond
!
      do nk = 1,nonzero
        d1BOi_s(nk)    = d1BOi_s(nk)    + BOp(ni,i)*df45dr(nk)
        d1BOi_pi(nk)   = d1BOi_pi(nk)   + BOp_pi(ni,i)*df45dr(nk)
        d1BOi_pipi(nk) = d1BOi_pipi(nk) + BOp_pipi(ni,i)*df45dr(nk)
      enddo
    endif
!
!  Free local memory
!
    if (lgrad1) then
      deallocate(df45dr,stat=status)
      if (status/=0) call deallocate_error('reaxFF_bo','df45dr')
    endif
  else
!*********************************
!  No corrections to be applied  *
!*********************************
    nonzero = 1
    nonzeroptr(1) = ni
    nonzeroi = nonzero
!  
!  Set local variable for bond order prime of currrent bond
!  
    BOpij_s  = BOp(ni,i)
    BOpij_pi = BOp_pi(ni,i) + BOp_pipi(ni,i)
    BOpij = BOpij_s + BOpij_pi
!  
!  Calculate corrected bond orders for i-j 
!  
    BOij      = BOpij 
    BOij_pi   = BOp_pi(ni,i) 
    BOij_pipi = BOp_pipi(ni,i)
    BOij_s = BOij - BOij_pi - BOij_pipi
!
!  Compute derivatives of combined bond order term
!
    if (lgrad1) then
!
!  Derivatives of BO' 
!
      d1BOi_s(1)    = d1BOp(ni,i)
      d1BOi_pi(1)   = d1BOp_pi(ni,i)
      d1BOi_pipi(1) = d1BOp_pipi(ni,i)
    endif
  endif
#ifdef TRACE
  call trace_out('reaxFF_bomds')
#endif
!
  return
  end
