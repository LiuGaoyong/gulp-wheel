module m_gfnff_c6
  use datatypes
!
  implicit none
!
!  Variables for dispersion model with D4
!
  integer(i4), dimension(:),       pointer, save :: d4_nref => null()
  real(dp),    dimension(:,:),     pointer, save :: d4_cn => null()
  real(dp),    dimension(:,:,:,:), pointer, save :: d4_c6 => null()
  real(dp),    dimension(:),       pointer, save :: d4_c9 => null()
  real(dp),    dimension(:),       pointer, save :: d4_r0 => null()
  real(dp),    dimension(:),       pointer, save :: d4_sqrtZr4r2 => null()
  real(dp),    dimension(:),       pointer, save :: d4_zeta_c6 => null()
!
  integer(i4),                              save :: maxc6ele = 1
  integer(i4),                              save :: maxc6ref = 7

CONTAINS

  subroutine gfnff_cnc6(numat,nat,cn,gw,dgwdcn,d2gwdcn2,lgrad1,lgrad2)
!
!  Set up factors that modify the C6 coefficents based on coordination number
!  for use in GFN-FF
!
  use datatypes
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: nat(numat)
  logical,     intent(in)  :: lgrad1
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(in)  :: cn(numat)
  real(dp),    intent(out) :: gw(maxc6ref,numat)
  real(dp),    intent(out) :: dgwdcn(maxc6ref,numat)     ! First derivative of gw w.r.t. cn
  real(dp),    intent(out) :: d2gwdcn2(maxc6ref,numat)   ! Second derivative of gw w.r.t. cn
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: iref
  integer(i4)              :: nati
  logical                  :: lgwkok
  real(dp)                 :: d1norm
  real(dp)                 :: d2norm
  real(dp)                 :: expw
  real(dp)                 :: d1expw
  real(dp)                 :: d2expw
  real(dp)                 :: gwk
  real(dp)                 :: gwloc
  real(dp)                 :: dgwloc
  real(dp)                 :: norm
  real(dp)                 :: wf
!
!  Initialise gw
!
  gw(1:maxc6ref,1:numat) = 0.0_dp
  if (lgrad1) then
    dgwdcn(1:maxc6ref,1:numat) = 0.0_dp
    if (lgrad2) then
      d2gwdcn2(1:maxc6ref,1:numat) = 0.0_dp
    endif
  endif
!
  wf = 4.0_dp
!
  do i = 1,numat
    nati = nat(i)
    norm = 0.0_dp
    d1norm = 0.0_dp
    d2norm = 0.0_dp
    if (lgrad2) then
      do iref = 1,d4_nref(nati)
        gwloc = exp(-wf*(cn(i)-d4_cn(iref,nati))**2)
        norm = norm + gwloc
        dgwloc = 2.0_dp*wf*(d4_cn(iref,nati) - cn(i))*gwloc
        d1norm = d1norm + dgwloc
        d2norm = d2norm + 2.0_dp*wf*((d4_cn(iref,nati) - cn(i))*dgwloc - gwloc)
      enddo
    elseif (lgrad1) then
      do iref = 1,d4_nref(nati)
        gwloc = exp(-wf*(cn(i)-d4_cn(iref,nati))**2)
        norm = norm + gwloc
        dgwloc = 2.0_dp*wf*(d4_cn(iref,nati) - cn(i))*gwloc
        d1norm = d1norm + dgwloc
      enddo
    else
      do iref = 1,d4_nref(nati)
        gwloc = exp(-wf*(cn(i)-d4_cn(iref,nati))**2)
        norm = norm + gwloc
      enddo
    endif
    norm = 1.0_dp/norm
    do iref = 1,d4_nref(nati)
      expw = exp(-wf*(cn(i)-d4_cn(iref,nati))**2)
!
      gwk = expw*norm
      lgwkok = .true.
      if (gwk /= gwk) then
        lgwkok = .false.
        if (maxval(d4_cn(1:d4_nref(nati),nati)).eq.d4_cn(iref,nati)) then
          gwk = 1.0_dp
        else
          gwk = 0.0_dp
        endif
      endif
      gw(iref,i) = gwk
!
      if (lgwkok) then
        if (lgrad1) then
          d1expw = 2.0_dp*wf*(d4_cn(iref,nati) - cn(i))*expw
          dgwdcn(iref,i) = d1expw*norm - expw*d1norm*norm**2
          if (lgrad2) then
            d2expw = 2.0_dp*wf*((d4_cn(iref,nati) - cn(i))*d1expw - expw)
            d2gwdcn2(iref,i) = d2expw*norm - 2.0_dp*d1expw*d1norm*norm**2 + &
                               2.0_dp*expw*d1norm*d1norm*norm**3 - expw*d2norm*norm**2
          endif
        endif
      endif
    enddo
  enddo

  end subroutine gfnff_cnc6

  subroutine gfnff_setc6(numat,nati,natj,i,j,gw,dgwdcn,d2gwdcn2,c6,dc6dcni,dc6dcnj, &
                         d2c6dcni2,d2c6dcnj2,d2c6dcnidcnj,lgrad1,lgrad2)
!
!  Set C6 coefficient for a pair of atoms
!
  use datatypes
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: j
  integer(i4), intent(in)  :: nati
  integer(i4), intent(in)  :: natj
  logical,     intent(in)  :: lgrad1
  logical,     intent(in)  :: lgrad2
  real(dp),    intent(in)  :: gw(maxc6ref,numat)
  real(dp),    intent(in)  :: dgwdcn(maxc6ref,*)
  real(dp),    intent(in)  :: d2gwdcn2(maxc6ref,*)
  real(dp),    intent(out) :: c6
  real(dp),    intent(out) :: dc6dcni
  real(dp),    intent(out) :: dc6dcnj
  real(dp),    intent(out) :: d2c6dcni2
  real(dp),    intent(out) :: d2c6dcnidcnj
  real(dp),    intent(out) :: d2c6dcnj2
!
!  Local variables
!
  integer(i4)              :: ii
  integer(i4)              :: jj
  real(dp)                 :: refc6
!
  c6 = 0.0_dp
  dc6dcni = 0.0_dp
  dc6dcnj = 0.0_dp
  d2c6dcni2 = 0.0_dp
  d2c6dcnidcnj = 0.0_dp
  d2c6dcnj2 = 0.0_dp
!
  if (lgrad2) then
    do ii = 1,d4_nref(nati)
      do jj = 1,d4_nref(natj)
        refc6 = d4_c6(ii,jj,nati,natj)
        c6 = c6 + gw(ii,i)*gw(jj,j)*refc6
        dc6dcni = dc6dcni + dgwdcn(ii,i)*gw(jj,j)*refc6
        dc6dcnj = dc6dcnj + dgwdcn(jj,j)*gw(ii,i)*refc6
        d2c6dcni2 = d2c6dcni2 + d2gwdcn2(ii,i)*gw(jj,j)*refc6
        d2c6dcnj2 = d2c6dcnj2 + d2gwdcn2(jj,j)*gw(ii,i)*refc6
        d2c6dcnidcnj = d2c6dcnidcnj + dgwdcn(ii,i)*dgwdcn(jj,j)*refc6
      enddo
    enddo
  elseif (lgrad1) then
    do ii = 1,d4_nref(nati)
      do jj = 1,d4_nref(natj)
        refc6 = d4_c6(ii,jj,nati,natj)
        c6 = c6 + gw(ii,i)*gw(jj,j)*refc6
        dc6dcni = dc6dcni + dgwdcn(ii,i)*gw(jj,j)*refc6
        dc6dcnj = dc6dcnj + dgwdcn(jj,j)*gw(ii,i)*refc6
      enddo
    enddo
  else
    do ii = 1,d4_nref(nati)
      do jj = 1,d4_nref(natj)
        refc6 = d4_c6(ii,jj,nati,natj)
        c6 = c6 + gw(ii,i)*gw(jj,j)*refc6
      enddo
    enddo
  endif

  end subroutine gfnff_setc6

  subroutine changemaxd4c6
!
!  Alters the size of the arrays associated with maxc6ref
!
!  12/21 Created
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
!  Copyright Curtin University 2021
!
!  Julian Gale, CIC, Curtin University, December 2021
!
  use reallocate
  implicit none
!
  integer(i4) :: ierror
!
  call realloc(d4_nref,maxc6ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd4c6','d4_nref')
  call realloc(d4_cn,maxc6ref,maxc6ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd4c6','d4_cn')
  call realloc(d4_c6,maxc6ref,maxc6ref,maxc6ele,maxc6ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd4c6','d4_c6')
  call realloc(d4_r0,maxc6ele*(maxc6ele+1_i4)/2_i4,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd4c6','d4_r0')
  call realloc(d4_sqrtZr4r2,maxc6ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxd4c6','d4_sqrtZr4r2')
!
  return
  end

end module m_gfnff_c6
