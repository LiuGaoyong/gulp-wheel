  subroutine gfnff_outpot
!
!  Outputs pGFNFF parameter information 
!  NB: Should only be called if GFNFF is to be used and this is the I/O processor
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use current
  use gulp_gfnff
  use iochannels
  implicit none
!
!*****************************
!  Output pGFNFF information  *
!*****************************
  write(ioout,'(/,''  pGFNFF forcefield to be used'',/)')
  if (lgfnff_xtbtopo) then
    write(ioout,'(''  pGFNFF parameter charges will use original XTB topology '')') 
  else
    write(ioout,'(''  pGFNFF parameter charges will use new PBC-consistent topology '')') 
    write(ioout,'(''  pGFNFF maximum number of shells for topology - loop 0   = '',i4)') maxtoposhell1
    write(ioout,'(''  pGFNFF maximum number of shells for topology - loop > 0 = '',i4)') maxtoposhell
  endif
  if (lgfnff_fragment_bond) then
    write(ioout,'(''  pGFNFF fragments to be set based on bonds '')') 
  else
    write(ioout,'(''  pGFNFF fragments to be set based on original neighbour list'')') 
  endif
  if (lgfnff_newtopo) then
    write(ioout,'(''  pGFNFF topological large distance = 10**12 '')') 
  else
    write(ioout,'(''  pGFNFF topological large distance = 13 '')') 
  endif
  write(ioout,'(''  pGFNFF topological distances scaled by = '',f8.5)') gfnff_rtopo_scale
  write(ioout,'(''  pGFNFF topological distance  cutoff    = '',f8.5,'' a.u.'')') tdist_thr
  if (lgfnff_topowolf) then
    write(ioout,'(''  pGFNFF topological charges to be computed with a Wolf sum with eta = '',f8.5)') gfnff_wolf_eta
  endif
  write(ioout,'(''  pGFNFF minimum allowed charge in topology set = '',f12.4)') gfnff_q_trap
  write(ioout,'(''  pGFNFF dispersion-coordination number tolerance = '',f10.4)') - log10(gfnff_cnc6tol)
  if (lgfnff_highcn_trap) then
    write(ioout,'(''  pGFNFF high coordination numbers will be trapped'')')
  else
    write(ioout,'(''  pGFNFF high coordination numbers will NOT be trapped'')')
  endif
  if (gfnff_pi_change.eq.0) then
    write(ioout,'(''  pGFNFF pi biradicals modified based on charge of atoms'')')
  elseif (gfnff_pi_change.gt.0) then
    write(ioout,'(''  pGFNFF pi biradicals modified by increase'')')
  else
    write(ioout,'(''  pGFNFF pi biradicals modified by decrease'')')
  endif
  write(ioout,'(''  pGFNFF accuracy       = '',f12.6)') gfnff_accuracy
  write(ioout,'(''  pGFNFF accuracy_disp  = '',f12.6)') gfnff_accuracy_disp
  write(ioout,'(''  pGFNFF accuracy_rep   = '',f12.6)') gfnff_accuracy_rep
  write(ioout,'(''  pGFNFF accuracy_cn    = '',f12.6)') gfnff_accuracy_cn
  write(ioout,'(''  pGFNFF accuracy_hb1   = '',f12.6)') gfnff_accuracy_hb1
  write(ioout,'(''  pGFNFF accuracy_hb2   = '',f12.6)') gfnff_accuracy_hb2
  write(ioout,'(''  pGFNFF taper          = '',f12.6)') gfnff_taper
  write(ioout,'(''  pGFNFF temperature 1  = '',f12.6,'' K'')') gfnff_pi_temp1
  write(ioout,'(''  pGFNFF temperature 2  = '',f12.6,'' K'')') gfnff_pi_temp2
  if (ndim.gt.0) then
    write(ioout,'(''  pGFNFF kspace         = '',f12.6,'' 1/Ang'')') gfnff_ks
  endif
  if (atm_alpha1.gt.0.0_dp) then
    write(ioout,'(''  pGFNFF ATM dispersion =   damped '')') 
  else
    write(ioout,'(''  pGFNFF ATM dispersion =   undamped '')') 
  endif
!
  return
  end
