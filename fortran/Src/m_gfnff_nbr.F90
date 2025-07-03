module m_gfnff_nbr
! 
!  This module contains the variables that hold the configuration specific
!  neighbour list information for GFNFF in GULP
!
!   8/20 Created
!   6/21 Use of xclat/yclat/zclat changed to xalat/yalat/zalat for benefit of MD
!   7/21 In non-MD case xclat/yclat/zclat are now used to allow for symmetry
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
!  Julian Gale, Curtin University, July 2021
!
  use datatypes
!
  implicit none
!
  integer(i4),                              save :: maxnbr = 6          ! Maximum number of neighbours
  real(dp),                                 save :: cutnbr = 0.0_dp     ! Cutoff used to generate list
!
!  General neighbour list
!
  integer(i4), dimension(:),       pointer, save :: nnbr => null()      ! Number of neighbours
  integer(i4), dimension(:,:),     pointer, save :: nbrno => null()     ! Pointer to atom that is neighbour
  integer(i4), dimension(:,:),     pointer, save :: ncnbr => null()     ! Pointer to cell of neighbour 
  logical,     dimension(:),       pointer, save :: latmdone => null()  ! Flag that atom has been done
  real(dp),    dimension(:,:),     pointer, save :: rnbr => null()      ! Distance to neighbour
  real(dp),    dimension(:,:),     pointer, save :: xnbr => null()      ! x component of distance to neighbour
  real(dp),    dimension(:,:),     pointer, save :: ynbr => null()      ! y component of distance to neighbour
  real(dp),    dimension(:,:),     pointer, save :: znbr => null()      ! z component of distance to neighbour
!
!  Coordination number pointers to data within nnbr / nbrno and derivatives
!
  integer(i4), dimension(:),       pointer, save :: nnbr_cn   => null()      ! Number of neighbours with non-zero coordination number
  integer(i4), dimension(:,:),     pointer, save :: nbrno_cn  => null()      ! Pointer to neighbour that has non-zero coordination number
  real(dp),    dimension(:,:),     pointer, save :: d1cndr_cn => null()      ! First derivatives w.r.t. distance for this coordination number term
  real(dp),    dimension(:,:),     pointer, save :: d2cndr_cn => null()      ! Second derivatives w.r.t. distance for this coordination number term
!
!  Bonding lists
!
  integer(i4), dimension(:),       pointer, save :: nnbr_bond => null()      ! Number of bonded neighbours
  integer(i4), dimension(:,:),     pointer, save :: nbrno_bond => null()     ! Pointer to atom that is a bonded neighbour
  integer(i4), dimension(:,:),     pointer, save :: nbrno_hb   => null()     ! Hydrogen bonding pointer for bond (nr_hb)
  integer(i4), dimension(:,:),     pointer, save :: ncnbr_bond => null()     ! Pointer to cell of bonded neighbour 
  real(dp),    dimension(:,:,:),   pointer, save :: par_gfnff_bond => null() ! Parameters for the interaction with a bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: rbnbr => null()          ! Distance to bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: xbnbr => null()          ! x component of distance to bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: ybnbr => null()          ! y component of distance to bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: zbnbr => null()          ! z component of distance to bonded neighbour

CONTAINS

  subroutine changemaxnbr
!
!  Changes the size of arrays that hold the neighbour list for GFNFF
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
!  Julian Gale, CIC, Curtin University, August 2020
!
  use current,     only : maxat
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(nnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','nnbr')
  call realloc(nbrno,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','nbrno')
  call realloc(ncnbr,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','ncnbr')
  call realloc(latmdone,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','latmdone')
  call realloc(rnbr,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','rnbr')
  call realloc(xnbr,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','xnbr')
  call realloc(ynbr,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','ynbr')
  call realloc(znbr,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','znbr')
!
  call realloc(nnbr_cn,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','nnbr_cn')
  call realloc(nbrno_cn,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','nbrno_cn')
  call realloc(d1cndr_cn,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','d1cndr_cn')
  call realloc(d2cndr_cn,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','d2cndr_cn')
!
  call realloc(nnbr_bond,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','nnbr_bond')
  call realloc(nbrno_bond,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','nbrno_bond')
  call realloc(ncnbr_bond,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','ncnbr_bond')
  call realloc(nbrno_hb,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','nbrno_hb')
  call realloc(par_gfnff_bond,3_i4,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','par_gfnff_bond')
  call realloc(rbnbr,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','rbnbr')
  call realloc(xbnbr,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','xbnbr')
  call realloc(ybnbr,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','ybnbr')
  call realloc(zbnbr,maxnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxnbr','zbnbr')
!
  end subroutine changemaxnbr

  subroutine gfnff_getnbr(cut2,lsetup)
!
!  Computes the neighbour list for GFNFF 
!
!  On entry : 
!
!  cut2            = cutoff squared
!  lsetup          = if .true. then this is a setup call
!
!  On exit :
!
!  Neighbour list is set
!
!   8/20 Created
!   2/21 Switch to using rlist2 variables instead of rlist
!   3/21 Trap for atoms being on top of each other added
!   6/21 Use of xclat/yclat/zclat changed to xalat/yalat/zalat for benefit of MD
!   7/21 Except for setup coordinates now use inbox arrays
!   8/21 Checking of maxnbr moved
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
!  Julian Gale, CIC, Curtin University, August 2021
!
  use datatypes
  use control,        only : keyword
  use current,        only : numat, nat, iimax2, iimid2
  use current,        only : xalat, yalat, zalat
  use current,        only : xclat, yclat, zclat
  use current,        only : xvec2cell, yvec2cell, zvec2cell
  use spatial,        only : lbuffercell
  use spatial,        only : lspatialok
  use spatial,        only : ncellsearch
  use spatial,        only : ncellpernode
  use spatial,        only : ncellnodeptr
  use spatial,        only : nspcell
  use spatial,        only : nspcellat
  use spatial,        only : nspcellatptr
  use spatial,        only : nspcellat1ptr
  use spatial,        only : nspcellatptrcell
  use spatial,        only : xinbox
  use spatial,        only : yinbox
  use spatial,        only : zinbox
  use iochannels
  use mdlogic,        only : lmd
  use neighbours
  use parallel,       only : ioproc
  use spatial
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,     intent(in)                        :: lsetup
  real(dp),    intent(in)                        :: cut2
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ic
  integer(i4)                                    :: ii
  integer(i4)                                    :: imx
  integer(i4)                                    :: imy
  integer(i4)                                    :: imz
  integer(i4)                                    :: ind1
  integer(i4)                                    :: ind2
  integer(i4)                                    :: ix
  integer(i4)                                    :: ixyz
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: j
  integer(i4)                                    :: jc
  integer(i4)                                    :: jj
  integer(i4)                                    :: maxxy
  integer(i4)                                    :: maxx
  integer(i4)                                    :: n1
  integer(i4)                                    :: n1j
  integer(i4)                                    :: ni
  integer(i4)                                    :: nj
  integer(i4)                                    :: nn
  integer(i4)                                    :: nsplower(3)
  integer(i4)                                    :: nspupper(3)
  logical                                        :: lzeroR
  real(dp)                                       :: r2
  real(dp)                                       :: rij
  real(dp),                                 save :: rzerotol = 1.0d-12
  real(dp)                                       :: xi
  real(dp)                                       :: yi
  real(dp)                                       :: zi
  real(dp)                                       :: xji
  real(dp)                                       :: yji
  real(dp)                                       :: zji
  real(dp)                                       :: xji0
  real(dp)                                       :: yji0
  real(dp)                                       :: zji0
#ifdef TRACE
  call trace_in('gfnff_getnbr')
#endif
!
  cutnbr = sqrt(cut2)
!
!  Check that arrays are initialised to the correct size
!
  call changemaxnbr
!
!  Set number of neighbours to zero
!
  nnbr(1:numat) = 0
  latmdone(1:numat) = .false.
  lzeroR = .false.
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
  if (lspatialok) then
    maxxy = nspcell(1)*nspcell(2)
    maxx  = nspcell(1)
!
!  Loop over all spatial cells looking for non-buffer cells
!
    do ixyz = 1,ncellpernode
      if (.not.lbuffercell(ixyz)) then
        ind1 = ncellnodeptr(ixyz)
        ind2 = ind1 - 1
        iz = ind2/maxxy
        ind2 = ind2 - maxxy*iz
        iy = ind2/maxx
        ix = ind2 - maxx*iy + 1
        iy = iy + 1
        iz = iz + 1
!
!  Set cell search bounds
!
        nspupper(1) = min(ix+ncellsearch(1),nspcell(1))
        nspupper(2) = min(iy+ncellsearch(2),nspcell(2))
        nspupper(3) = min(iz+ncellsearch(3),nspcell(3))
        nsplower(1) = max(ix-ncellsearch(1),1)
        nsplower(2) = max(iy-ncellsearch(2),1)
        nsplower(3) = max(iz-ncellsearch(3),1)
!
!  Get number of atoms in this cell
!
        ni = nspcellat(ind1)
        n1 = nspcellat1ptr(ind1)
!
!  Loop over atoms in the cell finding neighbours
!
        do ii = 1,ni
          i = nspcellatptr(n1+ii)
          ic = nspcellatptrcell(n1+ii)
          latmdone(i) = .true.
!
!  Set coordinates
!
          xi = xinbox(i) + xvec2cell(ic)
          yi = yinbox(i) + yvec2cell(ic)
          zi = zinbox(i) + zvec2cell(ic)
!
!  Loop over neighbouring cells
!
          do imz = nsplower(3),nspupper(3)
            do imy = nsplower(2),nspupper(2)
              do imx = nsplower(1),nspupper(1)
                ind2 = (imz-1)*maxxy + (imy-1)*maxx + imx
!
!  Loop over atoms within neighbouring cells
!
                nj = nspcellat(ind2)
                n1j = nspcellat1ptr(ind2)
                do jj = 1,nj
                  j = nspcellatptr(n1j+jj)
!
!  Exclude self term
!   
                  if (i.ne.j.or.ind1.ne.ind2) then
                    jc = nspcellatptrcell(n1j+jj)
!
!  Set centre cell coordinate differences
! 
                    xji = xvec2cell(jc) + xinbox(j) - xi
                    yji = yvec2cell(jc) + yinbox(j) - yi
                    zji = zvec2cell(jc) + zinbox(j) - zi
! 
                    r2 = xji*xji + yji*yji + zji*zji
                    if (r2.lt.cut2) then
                      rij = sqrt(r2)
                      nnbr(i) = nnbr(i) + 1
                      if (nnbr(i).ge.maxnbr) then
                        maxnbr = maxnbr + 2
                        call changemaxnbr
                      endif
                      nbrno(nnbr(i),i) = j
                      ncnbr(nnbr(i),i) = jc
                      if (rij.lt.rzerotol) lzeroR = .true.
                      rnbr(nnbr(i),i) = rij
                      xnbr(nnbr(i),i) = xji
                      ynbr(nnbr(i),i) = yji
                      znbr(nnbr(i),i) = zji
                    endif
                  endif
!
                enddo
              enddo
            enddo
          enddo
!
        enddo
!
      endif
    enddo
  else
    do i = 1,numat
      latmdone(i) = .true.
!
!  Loop over atoms
!
      do j = 1,numat
!
!  Set centre cell coordinate differences
!
        if (lsetup) then
          if (lmd) then
            xji0 = xalat(j) - xalat(i)
            yji0 = yalat(j) - yalat(i)
            zji0 = zalat(j) - zalat(i)
          else
            xji0 = xclat(j) - xclat(i)
            yji0 = yclat(j) - yclat(i)
            zji0 = zclat(j) - zclat(i)
          endif
        else
          xji0 = xinbox(j) - xinbox(i)
          yji0 = yinbox(j) - yinbox(i)
          zji0 = zinbox(j) - zinbox(i)
        endif
!
!  Loop over unit cells
!
        do ii = 1,iimax2
!
!  Exclude self term
!
          if (i.ne.j.or.ii.ne.iimid2) then
            xji = xji0 + xvec2cell(ii)
            yji = yji0 + yvec2cell(ii)
            zji = zji0 + zvec2cell(ii)
            r2 = xji*xji + yji*yji + zji*zji
            if (r2.lt.cut2) then
              rij = sqrt(r2)
              if (rij.lt.rzerotol) lzeroR = .true.
              nnbr(i) = nnbr(i) + 1
              if (nnbr(i).ge.maxnbr) then
                maxnbr = maxnbr + 2
                call changemaxnbr
              endif
              nbrno(nnbr(i),i) = j
              ncnbr(nnbr(i),i) = ii
              rnbr(nnbr(i),i) = rij
              xnbr(nnbr(i),i) = xji
              ynbr(nnbr(i),i) = yji
              znbr(nnbr(i),i) = zji
            endif
          endif
        enddo
      enddo
    enddo
  endif
!*******************
!  Debug printing  *
!*******************
  if (index(keyword,'debu').ne.0) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
      write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    endif
    call mpbarrier
    if (ioproc) then
      do i = 1,numat
        write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nnbr(i)
      enddo
    endif
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  Neighbours of atoms :'',/)')
    endif
    call mpbarrier
    if (ioproc) then
      do i = 1,numat
        write(ioout,'(i4,8(1x,i4))') i,(nbrno(nn,i),nn=1,nnbr(i))
      enddo
    endif
    call mpbarrier
  endif
!
!  Trap zero distances
!
  if (lzeroR) then
    call outerror('atoms are too close together',0_i4)
    call stopnow('gfnff_getnbr')
  endif
#ifdef TRACE
  call trace_out('gfnff_getnbr')
#endif
!
  end subroutine gfnff_getnbr

  subroutine gfnff_update_nbr_bond
!
!  Updates the distances and vectors in the bonding neighbour list for GFNFF 
!
!   9/20 Created
!   2/21 Switch to using rlist2 variables instead of rlist
!   3/21 Spatial option added
!   6/21 Use of xclat/yclat/zclat changed to xalat/yalat/zalat for benefit of MD
!   7/21 In non-MD case xclat/yclat/zclat are now used to allow for symmetry
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
!  Julian Gale, CIC, Curtin University, July 2021
!
  use datatypes
  use current,        only : numat
  use current,        only : xalat, yalat, zalat
  use current,        only : xclat, yclat, zclat
  use current,        only : xvec2cell, yvec2cell, zvec2cell
  use iochannels
  use mdlogic,        only : lmd
  use neighbours
  use spatial
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ii
  integer(i4)                                    :: j
  integer(i4)                                    :: ni
  real(dp)                                       :: rij
  real(dp)                                       :: r2
  real(dp)                                       :: xji
  real(dp)                                       :: yji
  real(dp)                                       :: zji
  real(dp)                                       :: xji0
  real(dp)                                       :: yji0
  real(dp)                                       :: zji0
#ifdef TRACE
  call trace_in('gfnff_update_nbr_bond')
#endif
!*********************************************
!  Update bonding neighbour lists for atoms  *
!*********************************************
!
!  Loop over atoms
!
  do i = 1,numat
!
!  Loop over neighbours
!
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
!
!  Set centre cell coordinate differences
!
      if (lspatialok) then
        xji0 = xinbox(j) - xinbox(i)
        yji0 = yinbox(j) - yinbox(i)
        zji0 = zinbox(j) - zinbox(i)
      else
        if (lmd) then
          xji0 = xalat(j) - xalat(i)
          yji0 = yalat(j) - yalat(i)
          zji0 = zalat(j) - zalat(i)
        else
          xji0 = xclat(j) - xclat(i)
          yji0 = yclat(j) - yclat(i)
          zji0 = zclat(j) - zclat(i)
        endif
      endif
!
!  Add unit cell contribution
!
      ii = ncnbr_bond(ni,i)
      xji = xji0 + xvec2cell(ii)
      yji = yji0 + yvec2cell(ii)
      zji = zji0 + zvec2cell(ii)
!
      r2 = xji*xji + yji*yji + zji*zji
      rij = sqrt(r2)
!
      rbnbr(ni,i) = rij
      xbnbr(ni,i) = xji
      ybnbr(ni,i) = yji
      zbnbr(ni,i) = zji
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_update_nbr_bond')
#endif
!
  end subroutine gfnff_update_nbr_bond

  subroutine gfnff_update_nbr
!
!  Updates the distances and vectors in the neighbour list for GFNFF 
!
!   9/20 Created
!   2/21 Switch to using rlist2 variables instead of rlist
!   3/21 Spatial option added
!   6/21 Use of xclat/yclat/zclat changed to xalat/yalat/zalat for benefit of MD
!   7/21 In non-MD case xclat/yclat/zclat are now used to allow for symmetry
!   7/21 Minimum image added
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
!  Julian Gale, CIC, Curtin University, July 2021
!
  use datatypes
  use current,        only : numat
  use current,        only : xalat, yalat, zalat
  use current,        only : xclat, yclat, zclat
  use current,        only : xvec2cell, yvec2cell, zvec2cell
  use iochannels
  use mdlogic,        only : lmd
  use neighbours
  use spatial
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ii
  integer(i4)                                    :: j
  integer(i4)                                    :: ni
  real(dp)                                       :: rij
  real(dp)                                       :: r2
  real(dp)                                       :: xji
  real(dp)                                       :: yji
  real(dp)                                       :: zji
  real(dp)                                       :: xji0
  real(dp)                                       :: yji0
  real(dp)                                       :: zji0
#ifdef TRACE
  call trace_in('gfnff_update_nbr')
#endif
!*************************************
!  Update neighbour lists for atoms  *
!*************************************
!
!  Loop over atoms
!
  do i = 1,numat
!
!  Loop over neighbours
!
    do ni = 1,nnbr(i)
      j = nbrno(ni,i)
!
!  Set centre cell coordinate differences
!
      if (lspatialok) then
        xji0 = xinbox(j) - xinbox(i)
        yji0 = yinbox(j) - yinbox(i)
        zji0 = zinbox(j) - zinbox(i)
      else
        if (lmd) then
          xji0 = xalat(j) - xalat(i)
          yji0 = yalat(j) - yalat(i)
          zji0 = zalat(j) - zalat(i)
        else
          xji0 = xclat(j) - xclat(i)
          yji0 = yclat(j) - yclat(i)
          zji0 = zclat(j) - zclat(i)
        endif
      endif
!
!  Add unit cell contribution
!
      ii = ncnbr(ni,i)
      xji = xji0 + xvec2cell(ii)
      yji = yji0 + yvec2cell(ii)
      zji = zji0 + zvec2cell(ii)
      r2 = xji*xji + yji*yji + zji*zji
!
      rij = sqrt(r2)
!
      rnbr(ni,i) = rij
      xnbr(ni,i) = xji
      ynbr(ni,i) = yji
      znbr(ni,i) = zji
    enddo
  enddo
#ifdef TRACE
  call trace_out('gfnff_update_nbr')
#endif
!
  end subroutine gfnff_update_nbr

end module m_gfnff_nbr
