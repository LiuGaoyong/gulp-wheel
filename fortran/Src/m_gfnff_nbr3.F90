module m_gfnff_nbr3
! 
!  This module contains the information on the neighbours of an atom out to 3 bonds
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
!  Copyright Curtin University 2020
!
!  Julian Gale, Curtin University, October 2020
!
  use datatypes
!
  implicit none
!
  integer(i4),                              save :: maxn3atom = 36      ! Maximum number of different atoms in shells
  integer(i4),                              save :: maxn3nbr  = 1       ! Maximum number of images of an atom within the shells
  integer(i4),                              save :: n3atom              ! Number of different atoms in shells
  integer(i4),                              save :: nshell_set = 0      ! Current number of shells for which the list is set up
!
!  Arrays for list
!
  integer(i4), dimension(:),       pointer, save :: n3atomptr => null()    ! Pointer from n3atom to real atom index
  integer(i4), dimension(:),       pointer, save :: n3atomrptr => null()   ! Pointer from real atom index to n3atom
  integer(i4), dimension(:),       pointer, save :: n3nbr => null()        ! Number of images of each n3atom
  integer(i4), dimension(:,:),     pointer, save :: n3nbrshell => null()   ! Index of shell for each image (1, 2 or 3)
  real(dp),    dimension(:,:),     pointer, save :: x3nbr => null()        ! x component of vector to images of each n3atom
  real(dp),    dimension(:,:),     pointer, save :: y3nbr => null()        ! y component of vector to images of each n3atom
  real(dp),    dimension(:,:),     pointer, save :: z3nbr => null()        ! z component of vector to images of each n3atom
!
!  Bonding list
!
  integer(i4),                              save :: maxb3atm = 1           ! Maximum number of 1-4 atoms
  integer(i4),                              save :: nb3atm                 ! Number of 1-4 bonded atom pairs
  integer(i4),                              save :: nbatm                  ! Number of 1-4 bonded atoms x bonds for 1 x bonds for 4
  integer(i4), dimension(:,:),     pointer, save :: nb3list => null()      ! Bonding info that connects 1-4 atoms
  real(dp),    dimension(:,:),     pointer, save :: rilb3 => null()        ! Ril vector for 1-4 atoms

CONTAINS

  subroutine changemaxn3nbr
!
!  Changes the size of arrays that hold the data for neighbours of an atom out to 3 bonds
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)              :: ierror
!
  call realloc(n3nbrshell,maxn3nbr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3nbr','n3nbrshell')
  call realloc(x3nbr,maxn3nbr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3nbr','x3nbr')
  call realloc(y3nbr,maxn3nbr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3nbr','y3nbr')
  call realloc(z3nbr,maxn3nbr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3nbr','z3nbr')
!
  end subroutine changemaxn3nbr

  subroutine changemaxn3atom
!
!  Changes the size of arrays that hold the data for neighbours of an atom out to 3 bonds
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)              :: ierror
!
  call realloc(n3nbr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3atom','n3nbr')
  call realloc(n3atomptr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3atom','n3atomptr')
  call realloc(n3nbrshell,maxn3nbr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3atom','n3nbrshell')
  call realloc(x3nbr,maxn3nbr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3atom','x3nbr')
  call realloc(y3nbr,maxn3nbr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3atom','y3nbr')
  call realloc(z3nbr,maxn3nbr,maxn3atom,ierror)
  if (ierror.ne.0) call outofmemory('changemaxn3atom','z3nbr')
!
  end subroutine changemaxn3atom

  subroutine changemaxb3atm
!
!  Changes the size of arrays that hold the information for 1-4 connected atoms (3 bonds)
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(nb3list,4_i4,maxb3atm,ierror)
  if (ierror.ne.0) call outofmemory('changemaxb3atm','nb3list')
  call realloc(rilb3,3_i4,maxb3atm,ierror)
  if (ierror.ne.0) call outofmemory('changemaxb3atm','rilb3')
!
  end subroutine changemaxb3atm

  subroutine gfnff_get_n3atoms(numat,i,nshell)
!
!  Computes the information on the neighbours of an atom up to the shell specified by nshell
!  NB: nshell must be less than or equal to 3
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
!  Julian Gale, CIC, Curtin University, March 2021
!
  use datatypes
  use m_gfnff_nbr
  implicit none
!
!  Passed variables
!
  integer(i4),                       intent(in)  :: numat
  integer(i4),                       intent(in)  :: i
  integer(i4),                       intent(in)  :: nshell
!
!  Local variables
!
  integer(i4)                                    :: j
  integer(i4)                                    :: j3
  integer(i4)                                    :: k
  integer(i4)                                    :: k3
  integer(i4)                                    :: kt
  integer(i4)                                    :: l
  integer(i4)                                    :: l3
  integer(i4)                                    :: lt
  integer(i4)                                    :: ni
  integer(i4)                                    :: nj
  integer(i4)                                    :: nk
  logical                                        :: lfound
  real(dp)                                       :: diff
  real(dp)                                       :: xik
  real(dp)                                       :: yik
  real(dp)                                       :: zik
  real(dp)                                       :: xil
  real(dp)                                       :: yil
  real(dp)                                       :: zil
!
!  Check nshell argument
!
  if (nshell.lt.1.or.nshell.gt.3) then
    call outerror('gfnff_get_n3atoms called with an invalid value of nshell',0_i4)
    call stopnow('gfnff_get_n3atoms')
  endif
!
!  Initialise values
!
  n3atom = 0
  n3atomrptr(1:numat) = 0
!
!  Store nshell value for future checking
!
  nshell_set = nshell
!**********************************
!  Add first shell of neighbours  *
!**********************************
  do ni = 1,nnbr_bond(i)
    j = nbrno_bond(ni,i)
    j3 = n3atomrptr(j)
    if (j3.eq.0) then
!
!  New atom in shells
!
      n3atom = n3atom + 1
      if (n3atom.gt.maxn3atom) then
        maxn3atom = n3atom + 6
        call changemaxn3atom
      endif
      n3atomptr(n3atom) = j
      n3atomrptr(j) = n3atom
      n3nbr(n3atom) = 1
      n3nbrshell(1,n3atom) = 1
      x3nbr(1,n3atom) = xbnbr(ni,i)
      y3nbr(1,n3atom) = ybnbr(ni,i)
      z3nbr(1,n3atom) = zbnbr(ni,i)
    else
!
!  Add image for this atom
!
      n3nbr(j3) = n3nbr(j3) + 1
      if (n3nbr(j3).gt.maxn3nbr) then
        maxn3nbr = n3nbr(j3) + 1
        call changemaxn3nbr
      endif
      n3nbrshell(n3nbr(j3),n3atom) = 1
      x3nbr(n3nbr(j3),j3) = xbnbr(ni,i)
      y3nbr(n3nbr(j3),j3) = ybnbr(ni,i)
      z3nbr(n3nbr(j3),j3) = zbnbr(ni,i)
    endif
  enddo
!***********************************
!  Add second shell of neighbours  *
!***********************************
  if (nshell.ge.2) then
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
!*****************************************
!  Loop over second shell of neighbours  *
!*****************************************
      do nj = 1,nnbr_bond(j)
        k = nbrno_bond(nj,j)
        k3 = n3atomrptr(k)
!
!  Set vector from i
!
        xik = xbnbr(ni,i) + xbnbr(nj,j)
        yik = ybnbr(ni,i) + ybnbr(nj,j)
        zik = zbnbr(ni,i) + zbnbr(nj,j)
!
!  Exclude i = k for the same image
!
        if (i.eq.k) then
          if ((abs(xik) + abs(yik) + abs(zik)).lt.1.0d-2) cycle
        endif
!
        if (k3.eq.0) then
!
!  New atom in shells
!
          n3atom = n3atom + 1
          if (n3atom.gt.maxn3atom) then
            maxn3atom = n3atom + 6
            call changemaxn3atom
          endif
          n3atomptr(n3atom) = k
          n3atomrptr(k) = n3atom
          n3nbr(n3atom) = 1
          n3nbrshell(1,n3atom) = 2
          x3nbr(1,n3atom) = xik
          y3nbr(1,n3atom) = yik
          z3nbr(1,n3atom) = zik
        else
!
!  Check whether this image is a duplicate
!
          lfound = .false.
          do kt = 1,n3nbr(k3)
            diff = abs(x3nbr(kt,k3) - xik) + abs(y3nbr(kt,k3) - yik) + abs(z3nbr(kt,k3) - zik)
            lfound = (diff.lt.1.0d-2)
            if (lfound) exit
          enddo
          if (lfound) cycle
!
!  Add image for this atom
!
          n3nbr(k3) = n3nbr(k3) + 1
          if (n3nbr(k3).gt.maxn3nbr) then
            maxn3nbr = n3nbr(k3) + 1
            call changemaxn3nbr
          endif
          n3nbrshell(n3nbr(k3),n3atom) = 2
          x3nbr(n3nbr(k3),k3) = xik
          y3nbr(n3nbr(k3),k3) = yik
          z3nbr(n3nbr(k3),k3) = zik
        endif
      enddo
    enddo
  endif
!**********************************
!  Add third shell of neighbours  *
!**********************************
  if (nshell.ge.3) then
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
!*****************************************
!  Loop over second shell of neighbours  *
!*****************************************
      do nj = 1,nnbr_bond(j)
        k = nbrno_bond(nj,j)
        xik = xbnbr(ni,i) + xbnbr(nj,j)
        yik = ybnbr(ni,i) + ybnbr(nj,j)
        zik = zbnbr(ni,i) + zbnbr(nj,j)
!****************************************
!  Loop over third shell of neighbours  *
!****************************************
        do nk = 1,nnbr_bond(k)
          l = nbrno_bond(nk,k)
          l3 = n3atomrptr(l)
!
!  Set vector from i
!
          xil = xik + xbnbr(nk,k)
          yil = yik + ybnbr(nk,k)
          zil = zik + zbnbr(nk,k)
!
!  Exclude i = l for the same image
!
          if (i.eq.l) then
            if ((abs(xil) + abs(yil) + abs(zil)).lt.1.0d-2) cycle
          endif
!
          if (l3.eq.0) then
!
!  New atom in shells 
!
            n3atom = n3atom + 1
            if (n3atom.gt.maxn3atom) then
              maxn3atom = n3atom + 6
              call changemaxn3atom
            endif
            n3atomptr(n3atom) = l
            n3atomrptr(l) = n3atom
            n3nbr(n3atom) = 1
            n3nbrshell(1,n3atom) = 3
            x3nbr(1,n3atom) = xil
            y3nbr(1,n3atom) = yil
            z3nbr(1,n3atom) = zil
          else
!
!  Check whether this image is a duplicate
!
            lfound = .false.
            do lt = 1,n3nbr(l3)
              diff = abs(x3nbr(lt,l3) - xil) + abs(y3nbr(lt,l3) - yil) + abs(z3nbr(lt,l3) - zil)
              lfound = (diff.lt.1.0d-2) 
              if (lfound) exit
            enddo
            if (lfound) cycle
!
!  Add image for this atom
!
            n3nbr(l3) = n3nbr(l3) + 1
            if (n3nbr(l3).gt.maxn3nbr) then
              maxn3nbr = n3nbr(l3) + 1
              call changemaxn3nbr
            endif
            n3nbrshell(n3nbr(l3),n3atom) = 3
            x3nbr(n3nbr(l3),l3) = xil
            y3nbr(n3nbr(l3),l3) = yil
            z3nbr(n3nbr(l3),l3) = zil
          endif
        enddo  ! End of loop over third shell
      enddo    ! End of loop over second shell
    enddo      ! End of loop over first shell
  endif
!
  end subroutine gfnff_get_n3atoms
!
  subroutine gfnff_bond_check(j,xij,yij,zij,nshell,lbonded)
!
!  Checks whether an atom j is within the bonding lists up to 
!  the nshell number of bonds
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
!  Copyright Curtin University 2020
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use datatypes
  use m_gfnff_nbr
  implicit none
!
!  Passed variables
!
  integer(i4),                       intent(in)  :: j
  integer(i4),                       intent(in)  :: nshell
  real(dp),                          intent(in)  :: xij
  real(dp),                          intent(in)  :: yij
  real(dp),                          intent(in)  :: zij
  logical,                           intent(out) :: lbonded
!
!  Local variables
!
  integer(i4)                                    :: j3
  integer(i4)                                    :: ni
  logical                                        :: lfound
  real(dp)                                       :: diff
!
!  Check that nshell is within values for which the list has been set up
!
  if (nshell.gt.nshell_set) then
    call outerror('gfnff_bond_check called with nshell greater than the maximum set',0_i4)
    call stopnow('gfnff_bond_check')
  endif
!
!  Initialise values
!
  lbonded = .false.
  j3 = n3atomrptr(j)
!
!  If j3 is 0 then j is not any shell
!
  if (j3.eq.0) return
!
!  Find out whether the precise vector is in the lists for j
!
  lfound = .false.
  do ni = 1,n3nbr(j3)
    diff = abs(x3nbr(ni,j3) - xij) + abs(y3nbr(ni,j3) - yij) + abs(z3nbr(ni,j3) - zij)
    lfound = (diff.lt.1.0d-2) 
    if (lfound) exit
  enddo
!
!  If vector has been found then is it within the request shell?
!
  if (lfound) lbonded = (n3nbrshell(ni,j3).le.nshell)
!
  end subroutine gfnff_bond_check

  subroutine gfnff_bond_shell(j,xij,yij,zij,nshell)
!
!  Checks whether an atom j is within the bonding lists and if so in which shell
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use datatypes
  use m_gfnff_nbr
  implicit none
!
!  Passed variables
!
  integer(i4),                       intent(in)  :: j
  integer(i4),                       intent(out) :: nshell
  real(dp),                          intent(in)  :: xij
  real(dp),                          intent(in)  :: yij
  real(dp),                          intent(in)  :: zij
!
!  Local variables
!
  integer(i4)                                    :: j3
  integer(i4)                                    :: ni
  logical                                        :: lfound
  real(dp)                                       :: diff
!
!  Initialise values
!
  nshell = 0
  j3 = n3atomrptr(j)
!
!  If j3 is 0 then j is not any shell
!
  if (j3.eq.0) return
!
!  Find out whether the precise vector is in the lists for j
!
  lfound = .false.
  do ni = 1,n3nbr(j3)
    diff = abs(x3nbr(ni,j3) - xij) + abs(y3nbr(ni,j3) - yij) + abs(z3nbr(ni,j3) - zij)
    lfound = (diff.lt.1.0d-2) 
    if (lfound) exit
  enddo
!
!  If vector has been found then set the shell for return
!
  if (lfound) nshell = n3nbrshell(ni,j3)
!
  end subroutine gfnff_bond_shell

end module m_gfnff_nbr3
