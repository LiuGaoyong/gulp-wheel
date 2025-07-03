  subroutine potreal(xji,yji,zji,maxextra,cuts2,qme)
!
!  Routine calculates the long range sum potential between two points
!  For 3-D this invokes the Ewald sum, while for 2-D it uses the Parry sum. 
!  The reciprocal space terms are assumed to have been computed elsewhere.
!
!   7/21 Created from qmatrixelement
!   7/21 Maxextra added
!
!  On entry:
!
!    xji      = difference in X coordinates of two points
!    yji      = difference in Y coordinates of two points
!    zji      = difference in Z coordinates of two points
!    maxextra = pairwise-specific extra looping due to i-j vector
!    cuts2    = core-shell cutoff, if applicable
!
!  On exit:
!
!    qme    = Coulomb potential between points (in Angs**-1)
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
!  Copyright Curtin Univerisity 2021
!
!  Julian Gale, CIC, Curtin University, July 2021
!
  use control
  use current
  use general,       only : cutw, etaw, selfwolf, selfwolfenergy
  use kspace
  use qmedata
  use shells
  use symmetry,      only : lra
#ifdef TRACE
  use trace,         only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: maxextra(3)
  real(dp),    intent(in)    :: xji
  real(dp),    intent(in)    :: yji
  real(dp),    intent(in)    :: zji
  real(dp),    intent(in)    :: cuts2
  real(dp),    intent(out)   :: qme
!
!  Local variables
!
  integer(i4)                :: ii
  integer(i4)                :: iim
  integer(i4)                :: iimn
  integer(i4)                :: iimx
  integer(i4)                :: jj
  integer(i4)                :: jjm
  integer(i4)                :: jjmn
  integer(i4)                :: jjmx
  integer(i4)                :: kk
  integer(i4)                :: kkm
  integer(i4)                :: kkmn
  integer(i4)                :: kkmx
  integer(i4)                :: ml1
  integer(i4)                :: ml2
  integer(i4)                :: ml3
  real(dp)                   :: dqme
  real(dp)                   :: d2qme
  real(dp)                   :: etaloc
  real(dp)                   :: r2
  real(dp)                   :: rl
  real(dp)                   :: rrl
  real(dp)                   :: rr2
  real(dp)                   :: rxi, ryi, rzi
  real(dp)                   :: rxj, ryj, rzj
  real(dp)                   :: rxk, ryk, rzk
  real(dp)                   :: setaloc
  real(dp)                   :: trm
  real(dp)                   :: xcd, ycd, zcd
  real(dp)                   :: xij
  real(dp)                   :: yji2
  real(dp)                   :: zji2
!
!  Functions
!
  real(dp)                   :: g_derfc
#ifdef TRACE
  call trace_in('potreal')
#endif
!
!  Zero matrix element / derivatives
!
  qme = 0.0_dp
  if (lwolf) then
    etaloc = etaw*etaw
    setaloc = etaw
  elseif (lewald) then
    etaloc = eta
    setaloc = seta
  endif
!*************************
!  Real space summation  *
!*************************
  if (.not.lnoreal) then
    if (ndim.gt.1.or.lwolf) then
      if (lminimage) then
!****************************
!  Minimum image algorithm  *
!****************************
!
!  Set up looping extents
!
        if (ndim.eq.3) then
          iimn = -1
          iimx =  1
          jjmn = -1
          jjmx =  1
          kkmn = -1
          kkmx =  1
        elseif (ndim.eq.2) then
          iimn = -1
          iimx =  1
          jjmn = -1
          jjmx =  1
          kkmn =  0
          kkmx =  0
        elseif (ndim.eq.1) then
          iimn = -1
          iimx =  1
          jjmn =  0
          jjmx =  0
          kkmn =  0
          kkmx =  0
        endif
!
!  Find minimum image
!
        if (lra) then
          rr2 = 1.0d10
          rxi = xji + (iimn-1)*r1x
          do ii = iimn,iimx
            rxi = rxi + r1x
            ryi = yji + (jjmn-1)*r2y
            do jj = jjmn,jjmx
              ryi = ryi + r2y
              rzi = zji + (kkmn-1)*r3z
              do kk = kkmn,kkmx
                rzi = rzi + r3z
                r2 = rxi*rxi + ryi*ryi + rzi*rzi
                if (r2.lt.rr2) then
                  rr2 = r2
                  rxk = rxi
                  ryk = ryi
                  rzk = rzi
                  iim = ii
                  jjm = jj
                  kkm = kk
                endif
              enddo
            enddo
          enddo
        else
          rr2 = 1.0d10
          rxi = xji + (iimn-1)*r1x
          ryi = yji + (iimn-1)*r1y
          rzi = zji + (iimn-1)*r1z
          do ii = iimn,iimx
            rxi = rxi + r1x
            ryi = ryi + r1y
            rzi = rzi + r1z
            rxj = rxi + (jjmn-1)*r2x
            ryj = ryi + (jjmn-1)*r2y
            rzj = rzi + (jjmn-1)*r2z
            do jj = jjmn,jjmx
              rxj = rxj + r2x
              ryj = ryj + r2y
              rzj = rzj + r2z
              xcd = rxj + (kkmn-1)*r3x
              ycd = ryj + (kkmn-1)*r3y
              zcd = rzj + (kkmn-1)*r3z
              do kk = kkmn,kkmx
                xcd = xcd + r3x
                ycd = ycd + r3y
                zcd = zcd + r3z
                r2 = xcd*xcd + ycd*ycd + zcd*zcd
                if (r2.lt.rr2) then
                  rr2 = r2
                  rxk = xcd
                  ryk = ycd
                  rzk = zcd
                  iim = ii
                  jjm = jj
                  kkm = kk
                endif
              enddo
            enddo
          enddo
        endif
!
!  Exclude distances outside maximum cutoff
!
        if (rr2.le.rmax2) then
!
!  Trap self term
!
          if (rr2.lt.1.0d-15) then
            if (lwolf) then
              qme = qme - selfwolf - tweatpi
            else
              qme = qme - tweatpi
            endif
          else
            rl = sqrt(rr2)
            rrl = 1.0_dp/rl
            if (rr2.lt.cuts2) then
!
!  Core-shell interaction
!
              qme = qme - rrl
            endif
            trm = g_derfc(setaloc*rl)*rrl
            qme = qme + trm
            if (lwolf) then
              qme = qme - selfwolf + selfwolfenergy*(rl - cutw)
            endif
          endif
        endif
      else
!******************************
!  Periodic 2-D and 3-D case  *
!******************************
        ml1 = maxloop(1) + maxextra(1)
        ml2 = maxloop(2) + maxextra(2)
        ml3 = maxloop(3) + maxextra(3)
!
!  Loop over cell vectors
!
        rxi = xji - (ml1+1)*r1x
        ryi = yji - (ml1+1)*r1y
        rzi = zji - (ml1+1)*r1z
        do ii = -ml1,ml1
          rxi = rxi + r1x
          ryi = ryi + r1y
          rzi = rzi + r1z
          rxj = rxi - (ml2+1)*r2x
          ryj = ryi - (ml2+1)*r2y
          rzj = rzi - (ml2+1)*r2z
          do jj = -ml2,ml2
            rxj = rxj + r2x
            ryj = ryj + r2y
            rzj = rzj + r2z
            rxk = rxj - (ml3+1)*r3x
            ryk = ryj - (ml3+1)*r3y
            rzk = rzj - (ml3+1)*r3z
            do kk = -ml3,ml3
              rxk = rxk + r3x
              ryk = ryk + r3y
              rzk = rzk + r3z
!
!  Calculate distance squared
!
              rr2 = rxk*rxk + ryk*ryk + rzk*rzk
!
!  Exclude distances outside maximum cutoff
!
              if (rr2.le.rmax2) then
!
!  Trap self term
!
                if (rr2.lt.1.0d-15) then
                  if (lwolf) then
                    qme = qme - selfwolf - tweatpi
                  else
                    qme = qme - tweatpi
                  endif
                else
                  rl = sqrt(rr2)
                  rrl = 1.0_dp/rl
                  if (rr2.lt.cuts2) then
!
!  Core-shell interaction
!
                    qme = qme - rrl
                  endif
                  trm = g_derfc(setaloc*rl)*rrl
                  qme = qme + trm
                  if (lwolf) then
                    qme = qme - selfwolf + selfwolfenergy*(rl - cutw)
                  endif
                endif
              endif
!
!  End of loops over lattice vectors
!
            enddo
          enddo
        enddo
      endif
    else
!*************************
!  Cluster and 1-D case  *
!*************************
      if (ndim.eq.1) then
        ml1 = maxloop(1) + maxextra(1)
        call qmatrix1D(xji,yji,zji,.false.,.false.,qme,dqme,d2qme)
      else
        ml1 = 0
        r1x = 0.0_dp
      endif
      xij = xji - (ml1+1)*r1x
      yji2 = yji*yji
      zji2 = zji*zji
      do ii = -ml1,ml1
        xij = xij + r1x
        rr2 = xij*xij + yji2 + zji2
!
!  Exclude core-shell interaction
!
        if (rr2.ge.cuts2.and.rr2.ge.1.0d-10) then
          rl = sqrt(rr2)
          rrl = 1.0_dp/rl
          qme = qme + rrl
        endif
      enddo
    endif
  endif
#ifdef TRACE
  call trace_out('potreal')
#endif
!
  return
  end
