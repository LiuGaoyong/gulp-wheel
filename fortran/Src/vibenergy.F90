  subroutine vibenergy(mcv,nllkpt,wkpt,freq,maxfd,rnokpt,lprint,fc)
!
!  Calculates the vibrational energetics from the phonon modes and outputs them
!
!   9/13 Created from phonon
!  11/13 Trap for inverting rkt added when T < 1.0d-6
!  12/16 Modified to handle cluster case
!   7/17 lprinloc moved to only wrap output and not calculation
!   8/17 linear now set for cases other than 0-D
!  10/17 fhenergy moved to energies module
!   2/18 Trace added
!   3/19 Multiple temperature ramps added
!   8/19 Correction to handling of zero point energy in equipartition free energy
!   3/21 Handling of Einstein model added
!   7/21 fhenergy_eq added
!  10/21 Terms rearranged and equipartition output changed to entropy
!   4/22 Output changed
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
!  Copyright Curtin University 2022
!
!  Julian Gale, CIC, Curtin University, April 2022
!
  use configurations
  use g_constants
  use control
  use current
  use element
  use energies,          only : fhenergy, fhenergy_eq
  use general
  use iochannels
  use parallel
  use properties,        only : cv, entropy
#ifdef TRACE
  use trace,             only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)                          :: mcv                 ! Number of modes
  integer(i4),  intent(in)                          :: nllkpt              ! Number of kpoints for this configuration
  integer(i4),  intent(in)                          :: maxfd               ! Maximum first dimension of the frequency array
  logical,      intent(in)                          :: lprint              ! If true then output results
  real(dp),     intent(in)                          :: fc                  ! Internal energy
  real(dp),     intent(in)                          :: wkpt(nllkpt)        ! Weights for k points
  real(dp),     intent(in)                          :: freq(maxfd,nllkpt)  ! Vibrational frequencies for the k points
  real(dp),     intent(in)                          :: rnokpt              ! Inverse sum of k point weights
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: mcvmax
  integer(i4)                                       :: mcvmin
  integer(i4)                                       :: nimag
  integer(i4)                                       :: nk
  integer(i4)                                       :: nt
  integer(i4)                                       :: nt0
  integer(i4)                                       :: ntr
  integer(i4)                                       :: status
  logical                                           :: linear
  logical                                           :: lnozero
  logical                                           :: lprinloc
  real(dp)                                          :: cmfact
  real(dp)                                          :: cv2
  real(dp)                                          :: ent2
  real(dp)                                          :: entropy_eq
  real(dp)                                          :: factor
  real(dp)                                          :: freqmin
  real(dp)                                          :: kinenergy
  real(dp)                                          :: logx
  real(dp)                                          :: rk
  real(dp)                                          :: rkt
  real(dp),     dimension(:),     allocatable       :: hw_over_kT
  real(dp)                                          :: tem
  real(dp)                                          :: trm1
  real(dp)                                          :: trmcv
  real(dp)                                          :: trmen
  real(dp)                                          :: trmen_eq
  real(dp)                                          :: trmfe
  real(dp)                                          :: trmfe_eq
  real(dp)                                          :: trmzp
  real(dp)                                          :: uvib_eq
  real(dp),     dimension(:),     allocatable       :: w1
  real(dp),     dimension(:),     allocatable       :: w2
  real(dp),     dimension(:),     allocatable       :: w3
  real(dp)                                          :: wk
  real(dp)                                          :: zpe
#ifdef TRACE
  call trace_in('vibenergy')
#endif
!
!  Set logicals
!
  lprinloc = (lprint.and.ioproc)
  lnozero  = (index(keyword,'noze').ne.0)
!
!  Allocate local pointer arrays
!
  allocate(hw_over_kT(mcv),stat=status)
  if (status/=0) call outofmemory('vibenergy','hw_over_kT')
  allocate(w1(mcv),stat=status)
  if (status/=0) call outofmemory('vibenergy','w1')
  allocate(w2(mcv),stat=status)
  if (status/=0) call outofmemory('vibenergy','w2')
  allocate(w3(mcv),stat=status)
  if (status/=0) call outofmemory('vibenergy','w3')
!
!  Set bounds on vibrational modes to consider
!
  if (ndim.eq.0) then
!
!  Cluster
!
    call moltype(linear)
    if (leinstein) then
      mcvmin = 1
    else
      if (linear) then
        mcvmin = 6
      else
        mcvmin = 7
      endif
    endif
  else
!
!  Periodic system
!
    linear = .false.
    mcvmin = 1
    mcvmax = mcv
  endif
!
!  Check for imaginary modes and exclude them
!
  nimag = 0
  do i = 1,mcv
    if (freq(i,1).lt.-0.5_dp) nimag = nimag + 1
  enddo
  if (linear) then
    mcvmin = mcvmin + max(0,nimag-2)
  else
    mcvmin = mcvmin + max(0,nimag-3)
  endif
  mcvmax = mcv
  if (minmode.ne.1) mcvmin = minmode
  if (maxmode.ne.0) mcvmax = maxmode
!*************************************
!  Output phonon related properties  *
!*************************************
  if (ntemperatureramp.gt.0) then
    do ntr = 1,ntemperatureramp
      tem = temperaturestart(ntr)
!
!  For first ramp need to start from initial temperature
!  For subsequent ramps this would be the same as the end of the previous ramp
!
      if (ntr.eq.1) then
        nt0 = 0
      else
        nt0 = 1
      endif
      do nt = nt0,ntemperaturestep(ntr)
!
!  Set temperature 
!
        tem = temperaturestart(ntr) + dble(nt)*temperaturestep(ntr)
!***************************************
!  Evaluate phonon related properties  *
!***************************************
!
!  Initialise thermodynamic properties
!
        zpe = 0.0_dp
        entropy = 0.0_dp
        entropy_eq = 0.0_dp
        fhenergy = 0.0_dp
        fhenergy_eq = 0.0_dp
        kinenergy = 0.0_dp
        uvib_eq = 0.0_dp
        cv = 0.0_dp
!
!  Set kT and k in eV
!
        rkt = boltz*tem/evtoj
        rk = boltz/evtoj
!
!  Set conversion factor for frequencies to unitless ratio to kT
!
        if (tem.gt.1.0d-6) then
          cmfact = planck*speedl/(boltz*tem)
        else
          cmfact = 0.0_dp
        endif
!
!  Loop over K points
!
        do nk = 1,nllkpt
          wk = wkpt(nk)*rnokpt
          if (tem.gt.1.0d-6) then
!
!  Scale frequencies to hw/kT
!
            do i = 1,mcv
              hw_over_kT(i) = cmfact*freq(i,nk)
            enddo
!
!  Store exp(x) in w1, exp(x)-1 in w2 and 1/(exp(x)-1) in w3
!
            do i = 1,mcv
              if (hw_over_kT(i).lt.12.0_dp) then
                w1(i) = exp(hw_over_kT(i))
                w2(i) = w1(i) - 1.0_dp
                if (abs(w2(i)).gt.0.0_dp) w3(i) = 1.0_dp/w2(i)
              else
                w3(i) = exp(-hw_over_kT(i))
              endif
            enddo
!
!  Zero point energy
!
            trmzp = 0.0_dp
            do i = mcvmin,mcvmax
              if (hw_over_kT(i).gt.cmfact) trmzp = trmzp + hw_over_kT(i)
            enddo
            trmzp = 0.5_dp*wk*rkt*trmzp
!
!  Kinetic energy
!
            do i = mcvmin,mcvmax
              kinenergy = kinenergy + 0.5_dp*wk*rkt*hw_over_kT(i)*(0.5_dp + w3(i))
            enddo
!
!  Entropy and free energy
!
            trmfe = 0.0_dp
            trmen = 0.0_dp
            freqmin = cmfact
            do i = mcvmin,mcvmax
              trm1 = hw_over_kT(i)
              if (trm1.gt.freqmin) then
                trm1 = 1.0_dp - exp(-trm1)
                logx = log(trm1)
                trmfe = trmfe + logx
                trmen = trmen + hw_over_kT(i)*w3(i) - logx
              endif
            enddo
!
            if (lnozero) then
              fhenergy = fhenergy + wk*rkt*trmfe 
            else
              fhenergy = fhenergy + wk*rkt*trmfe + trmzp
              zpe = zpe + trmzp
            endif
            entropy = entropy + wk*rk*trmen
!
!  Heat capacity - constant volume
!
            trmcv = 0.0_dp
            do i = mcvmin,mcvmax
              if (hw_over_kT(i).gt.freqmin) then
                if (hw_over_kT(i).lt.12.0_dp) then
                  trmcv = trmcv + hw_over_kT(i)*hw_over_kT(i)*w1(i)*w3(i)*w3(i)
                else
                  trmcv = trmcv + hw_over_kT(i)*hw_over_kT(i)*w3(i)
                endif
              endif
            enddo
            cv = cv + wk*rk*trmcv
!
!  Equipartition quantities
!
            trmfe_eq = 0.0_dp
            trmen_eq = 0.0_dp
            do i = mcvmin,mcvmax
              trm1 = hw_over_kT(i)
              if (trm1.gt.freqmin) then
                logx = log(trm1)
                trmfe_eq = trmfe_eq + logx
                trmen_eq = trmen_eq + 1.0_dp - logx
              endif
              uvib_eq = uvib_eq + wk*rkt
            enddo
            if (lnozero) then
              fhenergy_eq = fhenergy_eq + wk*rkt*trmfe_eq
            else
              fhenergy_eq = fhenergy_eq + wk*rkt*trmfe_eq + trmzp
            endif
            entropy_eq = entropy_eq + wk*rk*trmen_eq
          elseif (.not.lnozero) then
!
!  Zero point energy
!
            factor = 0.5_dp*wk*planck*speedl/evtoj
            trmzp = 0.0_dp
            do i = mcvmin,mcvmax
              if (freq(i,nk).gt.1.0_dp) trmzp = trmzp + freq(i,nk)
            enddo
            trmzp = factor*trmzp
            zpe = zpe + trmzp
          endif
!
!  End of loop over K points
!
        enddo
!
        if (lprinloc) then
          if (ndim.eq.0) then
            write(ioout,'(''  Vibrational properties (for cluster):  Temperature  =  '',f10.3,'' K'')') tem
          else
            write(ioout,'(''  Phonon properties (per mole of unit cells): Temperature = '',f10.3,'' K'')') tem
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'')')
          write(ioout,'(''  Zero point energy             = '',f15.6,'' eV'')') zpe
          if (tem.gt.1.0d-06) then
            ent2 = entropy*evtoj*avogadro
            cv2 = cv*evtoj*avogadro
            write(ioout,'(''  Entropy (TS)                  = '',f15.6,'' eV'')') entropy*tem
            write(ioout,'(''          (S)                   = '',f15.6,'' J/(mol.K)'')') ent2
            write(ioout,'(''  Helmholtz free-energy         = '',f15.6,'' eV'')') fhenergy + fc
            write(ioout,'(''                                = '',f15.6,'' kJmol-1'')') (fhenergy+fc)*evtoj*avogadro*0.001_dp
            write(ioout,'(''  Free energy  (equipartition)  = '',f15.6,'' eV'')') fc + fhenergy_eq
            write(ioout,'(''  Entropy (TS) (equipartition)  = '',f15.6,'' eV'')') entropy_eq*tem
            write(ioout,'(''  Uvib         (equipartition)  = '',f15.6,'' eV'')') uvib_eq
            write(ioout,'(''  Mean kinetic energy           = '',f15.6,'' eV'')') kinenergy
            write(ioout,'(''  Heat capacity - const volume  = '',f15.6,'' eV/K'')') cv
            write(ioout,'(''                                = '',f15.6,'' J/(mol.K)'')') cv2
          endif
          write(ioout,'(''--------------------------------------------------------------------------------'')')
        endif
      enddo
    enddo
  else
!
!  Single temperature
!
    tem = temperature
!***************************************
!  Evaluate phonon related properties  *
!***************************************
!
!  Initialise thermodynamic properties
!
    zpe = 0.0_dp
    entropy = 0.0_dp
    entropy_eq = 0.0_dp
    fhenergy = 0.0_dp
    fhenergy_eq = 0.0_dp
    kinenergy = 0.0_dp
    uvib_eq = 0.0_dp
    cv = 0.0_dp
!
!  Set kT and k in eV
!
    rkt = boltz*tem/evtoj
    rk = boltz/evtoj
!
!  Set conversion factor for frequencies to unitless ratio to kT
!
    if (tem.gt.1.0d-6) then
      cmfact = planck*speedl/(boltz*tem)
    else
      cmfact = 0.0_dp
    endif
!
!  Loop over K points
!
    do nk = 1,nllkpt
      wk = wkpt(nk)*rnokpt
      if (tem.gt.1.0d-6) then
!
!  Scale frequencies to hw/kT
!
        do i = 1,mcv
          hw_over_kT(i) = cmfact*freq(i,nk)
        enddo
!
!  Store exp(x) in w1, exp(x)-1 in w2 and 1/(exp(x)-1) in w3
!
        do i = 1,mcv
          if (hw_over_kT(i).lt.12.0_dp) then
            w1(i) = exp(hw_over_kT(i))
            w2(i) = w1(i) - 1.0_dp
            if (abs(w2(i)).gt.0.0_dp) w3(i) = 1.0_dp/w2(i)
          else
            w3(i) = exp(-hw_over_kT(i))
          endif
        enddo
!
!  Zero point energy
!
        trmzp = 0.0_dp
        do i = mcvmin,mcvmax
          if (hw_over_kT(i).gt.cmfact) trmzp = trmzp + hw_over_kT(i)
        enddo
        trmzp = 0.5_dp*wk*rkt*trmzp
!
!  Kinetic energy
!
        do i = mcvmin,mcvmax
          kinenergy = kinenergy + 0.5_dp*wk*rkt*hw_over_kT(i)*(0.5_dp + w3(i))
        enddo
!
!  Entropy and free energy
!
        trmfe = 0.0_dp
        trmen = 0.0_dp
        freqmin = cmfact
        do i = mcvmin,mcvmax
          trm1 = hw_over_kT(i)
          if (trm1.gt.freqmin) then
            trm1 = 1.0_dp - exp(-trm1)
            logx = log(trm1)
            trmfe = trmfe + logx
            trmen = trmen + hw_over_kT(i)*w3(i) - logx
          endif
        enddo
!
        if (lnozero) then
          fhenergy = fhenergy + wk*rkt*trmfe
        else
          fhenergy = fhenergy + wk*rkt*trmfe + trmzp
          zpe = zpe + trmzp
        endif
        entropy = entropy + wk*rk*trmen
!
!  Equipartition free energy
!
        trmfe_eq = 0.0_dp
        trmen_eq = 0.0_dp
        do i = mcvmin,mcvmax
          trm1 = hw_over_kT(i)
          if (trm1.gt.freqmin) then
            logx = log(trm1)
            trmfe_eq = trmfe_eq + logx
            trmen_eq = trmen_eq + 1.0_dp - logx
          endif
          uvib_eq = uvib_eq + wk*rkt
        enddo
        if (lnozero) then
          fhenergy_eq = fhenergy_eq + wk*rkt*trmfe_eq
        else
          fhenergy_eq = fhenergy_eq + wk*rkt*trmfe_eq + trmzp
        endif
        entropy_eq = entropy_eq + wk*rk*trmen_eq
!
!  Heat capacity - constant volume
!
        trmcv = 0.0_dp
        do i = mcvmin,mcvmax
          if (hw_over_kT(i).gt.freqmin) then
            if (hw_over_kT(i).lt.12.0_dp) then
              trmcv = trmcv + hw_over_kT(i)*hw_over_kT(i)*w1(i)*w3(i)*w3(i)
            else
              trmcv = trmcv + hw_over_kT(i)*hw_over_kT(i)*w3(i)
            endif
          endif
        enddo
        cv = cv + wk*rk*trmcv
      elseif (.not.lnozero) then
!
!  Zero point energy
!
        factor = 0.5_dp*wk*planck*speedl/evtoj
        trmzp = 0.0_dp
        do i = mcvmin,mcvmax
          if (freq(i,nk).gt.1.0_dp) trmzp = trmzp + freq(i,nk)
        enddo
        trmzp = factor*trmzp
        zpe = zpe + trmzp
      endif
!
!  End of loop over K points
!
    enddo
!
    if (lprinloc) then
      if (ndim.eq.0) then
        write(ioout,'(''  Vibrational properties (for cluster):  Temperature  =  '',f10.3,'' K'')') tem
      else
        write(ioout,'(''  Phonon properties (per mole of unit cells): Temperature = '',f10.3,'' K'')') tem
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''  Zero point energy             = '',f15.6,'' eV'')') zpe
      if (tem.gt.1.0d-06) then
        ent2 = entropy*evtoj*avogadro
        cv2 = cv*evtoj*avogadro
        write(ioout,'(''  Entropy (TS)                  = '',f15.6,'' eV'')') entropy*tem
        write(ioout,'(''          (S)                   = '',f15.6,'' J/(mol.K)'')') ent2
        write(ioout,'(''  Helmholtz free-energy         = '',f15.6,'' eV'')') fhenergy + fc
        write(ioout,'(''                                = '',f15.6,'' kJmol-1'')') (fhenergy+fc)*evtoj*avogadro*0.001_dp
        write(ioout,'(''  Free energy  (equipartition)  = '',f15.6,'' eV'')') fc + fhenergy_eq
        write(ioout,'(''  Entropy (TS) (equipartition)  = '',f15.6,'' eV'')') entropy_eq*tem
        write(ioout,'(''  Uvib         (equipartition)  = '',f15.6,'' eV'')') uvib_eq
        write(ioout,'(''  Mean kinetic energy           = '',f15.6,'' eV'')') kinenergy
        write(ioout,'(''  Heat capacity - const volume  = '',f15.6,'' eV/K'')') cv
        write(ioout,'(''                                = '',f15.6,'' J/(mol.K)'')') cv2
      endif
      write(ioout,'(''--------------------------------------------------------------------------------'')')
    endif
  endif
!
!  Free local memory
!
  deallocate(w3,stat=status)
  if (status/=0) call deallocate_error('vibenergy','w3')
  deallocate(w2,stat=status)
  if (status/=0) call deallocate_error('vibenergy','w2')
  deallocate(w1,stat=status)
  if (status/=0) call deallocate_error('vibenergy','w1')
  deallocate(hw_over_kT,stat=status)
  if (status/=0) call deallocate_error('vibenergy','hw_over_kT')
#ifdef TRACE
  call trace_out('vibenergy')
#endif
!
  return
  end
