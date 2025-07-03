  subroutine backgroundmd(ebgd,emad,lgrad1)
!
!  Calculates neutralising background contribution to energy
!  for non-charge neutral solids. MD version
!
!  12/21 Created from background
!   6/23 Referencing of nspecptr corrected to use asymmetric unit
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
  use g_constants
  use control,      only : lDoElectrostatics, lwolf, lmadelung, latomicstress
  use current
  use derivatives
  use kspace
  use m_strain,     only : strainddetds, straindet
  use m_ti
  use mdlogic
  use parallel
  use qmedata,      only : maxloop
  use species
#ifdef TRACE
  use trace,        only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,  intent(in)     :: lgrad1
  real(dp), intent(out)    :: ebgd
  real(dp), intent(out)    :: emad
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: nspeci
  real(dp)                 :: asstress
  real(dp)                 :: demad
  real(dp)                 :: dsumqdl
  real(dp)                 :: erealq
  real(dp)                 :: erecipq
  real(dp)                 :: esum
  real(dp)                 :: qi
  real(dp)                 :: reta
  real(dp)                 :: rvol
  real(dp)                 :: sf_i
  real(dp)                 :: sumq
  real(dp)                 :: vol
  real(dp)                 :: volume
#ifdef TRACE
  call trace_in('backgroundmd')
#endif
!
!  Check that electrostatics are being used
!
  emad = 0.0_dp
  if (lDoElectrostatics.and..not.lwolf) then
    if (ndim.eq.3) then
      vol = volume(rv)
      rvol = 1.0_dp/vol
      reta = 1.0_dp/eta
      if (lti) then
!
!  If this TI then recompute the total charge to allow for lambda
!
        sumq = 0.0_dp
        dsumqdl = 0.0_dp
        do i = 1,numat
          qi = qf(i)*occuf(i)
          nspeci = nspecptr(nrelf2a(i))
          sf_i = 1.0_dp
          if (ltifqspec(nspeci)) then
            sf_i = lambda
            dsumqdl = dsumqdl + qi
          elseif (ltibqspec(nspeci)) then
            sf_i = (1.0_dp - lambda)
            dsumqdl = dsumqdl - qi
          endif
          qi = qi*sf_i
          sumq = sumq + qi
        enddo
        dsumqdl = dsumqdl*rtistep
        dUdlambda(nlambdanow) = dUdlambda(nlambdanow) - pi*angstoev*sumq*rvol*reta*dsumqdl
      else
        sumq = totalcharge
      endif
      ebgd = - 0.5_dp*pi*angstoev*(sumq**2)*rvol*reta
      if (lmadelung) then
!
!  Check the system is cubic
!
        if (a.ne.b.or.b.ne.c.or.alpha.ne.90.0_dp.or.beta.ne.90.0_dp.or.gamma.ne.90_dp) then
          call outerror('Madelung correction can only be applied to cubic systems',0_i4)
          call stopnow('background')
        endif
!
        if (lti) then
          dUdlambda(nlambdanow) = dUdlambda(nlambdanow) + 2.837297_dp*angstoev*sumq*dsumqdl/a
        endif
!
        emad = 0.5_dp*2.837297_dp*angstoev*(sumq**2)/a
        demad = - emad/a
      else
        demad = 0.0_dp
      endif
      if (lgrad1) then
        esum = ebgd + demad
        if (latomicstress) then
          asstress = esum/dble(numat)
          do i = 1,numat
            atomicstress(1,i) = atomicstress(1,i) - asstress
            atomicstress(2,i) = atomicstress(2,i) - asstress
            atomicstress(3,i) = atomicstress(3,i) - asstress
          enddo
        endif
        if (lfinitestrain) then
!
!  Volume first derivatives
!
          do i = 1,6
            strderv(i) = strderv(i) - strainddetds(i)*straindet*esum
          enddo
        else
          strderv(1) = strderv(1) - esum
          strderv(2) = strderv(2) - esum
          strderv(3) = strderv(3) - esum
        endif
      endif
    elseif (ndim.eq.2) then
      erecipq = 0.0_dp
      erealq = 0.0_dp
      call recip2Dq(erecipq,sumq,lgrad1,.false.)
      call real2Dq(erealq,erecipq,sumq,lgrad1,.false.)
      ebgd = erecipq + erealq
    elseif (ndim.eq.1) then
      call setmaxcell1D(maxloop(1))
      call real1Dq(ebgd,sumq,lgrad1,.false.)
    endif
  else
    ebgd = 0.0_dp
  endif
#ifdef TRACE
  call trace_out('backgroundmd')
#endif
!
  return
  end
