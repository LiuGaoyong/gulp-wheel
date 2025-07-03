  subroutine gulp_gfnff_init
!
!  Initialises various parameters for GFNFF passed from GULP
!
!  10/21 Created
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
!  Julian Gale, CIC, Curtin University, October 2021
!
  use datatypes
  use gulp_gfnff
  use m_gfnff_c6
  implicit none
!
!  Local variables
!
  integer(i4)          :: ierror
#ifndef NOGFNFF
!
!  Pass parameters from GULP to GFNFF
!
  call pgfnff_init_param_i('pi_change',gfnff_pi_change,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: pi_change',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('pi_temp1',gfnff_pi_temp1,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: pi_temp1',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('pi_temp2',gfnff_pi_temp2,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: pi_temp2',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('ks',gfnff_ks,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: ks',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_i('maxtoposhell',maxtoposhell,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: maxtoposhell',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_i('maxtoposhell1',maxtoposhell1,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: maxtoposhell1',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_l('newtopo',lgfnff_newtopo,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: newtopo',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_l('xtbtopo',lgfnff_xtbtopo,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: xtbtopo',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_l('fragment_bond',lgfnff_fragment_bond,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: fragment_bond',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('tdist_thr',tdist_thr,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: tdist_thr',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('rtopo',gfnff_rtopo_scale,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: rtopo',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('atm_alpha1',atm_alpha1,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: atm_alpha1',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('q_trap',gfnff_q_trap,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: q_trap',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_l('highcn_trap',lgfnff_highcn_trap,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: highcn_trap',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('max_accuracy',max_gfnff_accuracy,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: max_accuracy',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('max_acc_disp',max_gfnff_acc_disp,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: max_acc_disp',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('max_acc_rep',max_gfnff_acc_rep,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: max_acc_rep',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('max_acc_cn',max_gfnff_acc_cn,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: max_acc_cn',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('max_acc_hb1',max_gfnff_acc_hb1,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: max_acc_hb1',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('max_acc_hb2',max_gfnff_acc_hb2,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: max_acc_hb2',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('accuracy_overall',gfnff_accuracy,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: accuracy_overall',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('accuracy_disp',gfnff_accuracy_disp,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: accuracy_disp',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('accuracy_rep',gfnff_accuracy_rep,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: accuracy_rep',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('accuracy_cn',gfnff_accuracy_cn,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: accuracy_cn',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('accuracy_hb1',gfnff_accuracy_hb1,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: accuracy_hb1',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('accuracy_hb2',gfnff_accuracy_hb2,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: accuracy_hb2',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('cnc6tol',gfnff_cnc6tol,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: cnc6tol',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('taper',gfnff_taper,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: taper',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_r('wolf_eta',gfnff_wolf_eta,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: wolf_eta',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
!
  call pgfnff_init_param_l('topowolf',lgfnff_topowolf,ierror)
  if (ierror.lt.0) then
    call outerror('GFNFF parameter not found when passed from GULP: topowolf',0_i4)
    call stopnow('gulp_gfnff_init')
  endif
#else
  call outerror('cannot initialise GFNFF parameters as pGFNFF library is not linked',0_i4)
  call stopnow('gulp_gfnff_init')
#endif
!
!  Set up memory to pass to pGFNFF library 
!
  call changemaxd4c6
!
  return
  end
