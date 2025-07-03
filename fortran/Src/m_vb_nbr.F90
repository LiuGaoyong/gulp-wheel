module m_vb_nbr
! 
!  This module contains the variables that hold the configuration specific
!  neighbour list information for valence bond potentials in GULP
!
!   8/21 Created
!
!  Julian Gale, Curtin University, August 2021
!
  use datatypes
  use m_cells
!
  implicit none
!
  integer(i4),                              save :: maxvbnbr = 6             ! Maximum number of neighbours
  integer(i4),                              save :: maxvalbondatm = 1        ! Maximum number of valbond potentials per atom
!
!  General neighbour list
!
  integer(i4), dimension(:),       pointer, save :: nvbnbr => null()         ! Number of neighbours
  integer(i4), dimension(:,:),     pointer, save :: nvbnbrno => null()       ! Pointer to atom that is neighbour
  integer(i4), dimension(:,:),     pointer, save :: nvbcnbr => null()        ! Pointer to cell of neighbour 
  logical,     dimension(:),       pointer, save :: lvbatmdone => null()     ! Flag that atom has been done
  real(dp),    dimension(:,:),     pointer, save :: rvbnbr => null()         ! Distance to neighbour
  real(dp),    dimension(:,:),     pointer, save :: xvbnbr => null()         ! x component of distance to neighbour
  real(dp),    dimension(:,:),     pointer, save :: yvbnbr => null()         ! y component of distance to neighbour
  real(dp),    dimension(:,:),     pointer, save :: zvbnbr => null()         ! z component of distance to neighbour
!
!  Arrays to helps find potentials for atoms and neighbours
!
  integer(i4), dimension(:),       pointer, save :: nvalbondatm => null()    ! Pointer to atom that is neighbour
  integer(i4), dimension(:,:),     pointer, save :: nvalbondatmptr => null() ! Pointer to atom that is neighbour
  real(dp),    dimension(:),       pointer, save :: rvbcut2 => null()        ! Atom based cutoff squared
!
  type(t_cells),                            save :: vb_cells                 ! Cell image list for VB

CONTAINS

  subroutine changemaxvbnbr
!
!  Changes the size of arrays that hold the neighbour list for valence bond potentials
!
!  Julian Gale, CIC, Curtin University, August 2021
!
  use current,      only : maxat
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(nvbnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvbnbr','nvbnbr')
  call realloc(nvbnbrno,maxvbnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvbnbr','nvbnbrno')
  call realloc(nvbcnbr,maxvbnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvbnbr','nvbcnbr')
  call realloc(lvbatmdone,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvbnbr','lvbatmdone')
  call realloc(rvbcut2,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvbnbr','rvbcut2')
  call realloc(rvbnbr,maxvbnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvbnbr','rvbnbr')
  call realloc(xvbnbr,maxvbnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvbnbr','xvbnbr')
  call realloc(yvbnbr,maxvbnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvbnbr','yvbnbr')
  call realloc(zvbnbr,maxvbnbr,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvbnbr','zvbnbr')
!
  end subroutine changemaxvbnbr

  subroutine changemaxvalbondatm
!
!  Changes the size of arrays that hold the pointer for atoms to valence bond potentials
!
!  Julian Gale, CIC, Curtin University, August 2021
!
  use current,      only : maxat
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(nvalbondatm,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbondatm','nvalbondatm')
  call realloc(nvalbondatmptr,maxvalbondatm,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxvalbondatm','nvalbondatmptr')
!
  end subroutine changemaxvalbondatm

  subroutine vb_getcut
!
!  Computes the atom-based overall cutoffs for valence bond potentials
!
!  On exit :
!
!  rvbcut2 contains the atom-based cutoffs squared
!
!   8/21 Created
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
  use bondvalence
  use current
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: n
  real(dp)                                       :: cut2
#ifdef TRACE
  call trace_in('vb_getcut')
#endif
!
!  Check that arrays are initialised to the correct size
!
  call changemaxvbnbr
  call changemaxvalbondatm
!
!  Initialise cutoffs to zero
!
  rvbcut2(1:numat) = 0.0_dp
!
!  Initialise pointers to potentials for atoms
!
  rvbcut2(1:numat) = 0.0_dp
  nvalbondatm(1:numat) = 0
!
  do i = 1,numat
    do n = 1,nvalbond
      if (nVBspecB1(n).eq.nat(i)) then
        if (nVBtypeB1(n).eq.nftype(i).or.nVBtypeB1(n).eq.0) then
          nvalbondatm(i) = nvalbondatm(i) + 1
          if (nvalbondatm(i).gt.maxvalbondatm) then
            maxvalbondatm = nvalbondatm(i) + 1
            call changemaxvalbondatm
          endif
          nvalbondatmptr(nvalbondatm(i),i) = n
          cut2 = rVBmax(n)**2
          rvbcut2(i) = max(rvbcut2(i),cut2)
        endif
      endif
    enddo
  enddo
#ifdef TRACE
  call trace_out('vb_getcut')
#endif
!
  end subroutine vb_getcut

  subroutine vb_getnbr
!
!  Computes the neighbour list for valence bond potentials
!
!  On exit :
!
!  Neighbour list is set
!
!   8/21 Created
!   9/21 Cell images now from dedicated list
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
  use bondvalence
  use control,        only : keyword
  use current,        only : numat, nat, nftype
  use current,        only : xvec2cell, yvec2cell, zvec2cell
  use element,        only : maxele
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
  use neighbours
  use parallel,       only : ioproc
  use spatial
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
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
  integer(i4)                                    :: m
  integer(i4)                                    :: maxxy
  integer(i4)                                    :: maxx
  integer(i4)                                    :: n
  integer(i4)                                    :: n1
  integer(i4)                                    :: n1j
  integer(i4)                                    :: ni
  integer(i4)                                    :: nj
  integer(i4)                                    :: nn
  integer(i4)                                    :: nsplower(3)
  integer(i4)                                    :: nspupper(3)
  logical                                        :: lallowed
  logical                                        :: lzeroR
  real(dp)                                       :: cut2
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
  call trace_in('vb_getnbr')
#endif
!
!  Check that arrays are initialised to the correct size
!
  call changemaxvbnbr
!
!  Set cells
!
  cut2 = 0.0_dp
  do i = 1,numat
    cut2 = max(rvbcut2(i),cut2)
  enddo
  call setcells(cut2,vb_cells)
!
!  Set number of neighbours to zero
!
  nvbnbr(1:numat) = 0
  lvbatmdone(1:numat) = .false.
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
          lvbatmdone(i) = .true.
          cut2 = rvbcut2(i)
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
!  Check atom is allowed for a valence bond potential
!
                    lallowed = .false.
                    m = 0
                    do while (.not.lallowed.and.m.lt.nvalbondatm(i))
                      m = m + 1
                      n = nvalbondatmptr(m,i)
                      lallowed = (nVBspecB2(n).eq.maxele.or.(nVBspecB2(n).eq.nat(j).and. &
                                  (nVBtypeB2(n).eq.nftype(j).or.nVBtypeB2(n).eq.0)))
                    enddo
                    if (lallowed) then
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
                        nvbnbr(i) = nvbnbr(i) + 1
                        if (nvbnbr(i).ge.maxvbnbr) then
                          maxvbnbr = maxvbnbr + 2
                          call changemaxvbnbr
                        endif
                        nvbnbrno(nvbnbr(i),i) = j
                        nvbcnbr(nvbnbr(i),i) = jc
                        if (rij.lt.rzerotol) lzeroR = .true.
                        rvbnbr(nvbnbr(i),i) = rij
                        xvbnbr(nvbnbr(i),i) = xji
                        yvbnbr(nvbnbr(i),i) = yji
                        zvbnbr(nvbnbr(i),i) = zji
                      endif
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
      lvbatmdone(i) = .true.
      cut2 = rvbcut2(i)
!
!  Loop over atoms
!
      do j = 1,numat
!
!  Check atom is allowed for a valence bond potential
!
        lallowed = .false.
        m = 0
        do while (.not.lallowed.and.m.lt.nvalbondatm(i))
          m = m + 1
          n = nvalbondatmptr(m,i)
          lallowed = (nVBspecB2(n).eq.maxele.or.(nVBspecB2(n).eq.nat(j).and.(nVBtypeB2(n).eq.nftype(j).or.nVBtypeB2(n).eq.0)))
        enddo
        if (lallowed) then
!
!  Set centre cell coordinate differences
!
          xji0 = xinbox(j) - xinbox(i)
          yji0 = yinbox(j) - yinbox(i)
          zji0 = zinbox(j) - zinbox(i)
!
!  Loop over unit cells
!
          do ii = 1,vb_cells%ncell
!
!  Exclude self term
!
            if (i.ne.j.or.ii.ne.vb_cells%ncentrecell) then
              xji = xji0 + vb_cells%xcell(ii)
              yji = yji0 + vb_cells%ycell(ii)
              zji = zji0 + vb_cells%zcell(ii)
              r2 = xji*xji + yji*yji + zji*zji
              if (r2.lt.cut2) then
                rij = sqrt(r2)
                if (rij.lt.rzerotol) lzeroR = .true.
                nvbnbr(i) = nvbnbr(i) + 1
                if (nvbnbr(i).ge.maxvbnbr) then
                  maxvbnbr = maxvbnbr + 2
                  call changemaxvbnbr
                endif
                nvbnbrno(nvbnbr(i),i) = j
                nvbcnbr(nvbnbr(i),i) = ii
                rvbnbr(nvbnbr(i),i) = rij
                xvbnbr(nvbnbr(i),i) = xji
                yvbnbr(nvbnbr(i),i) = yji
                zvbnbr(nvbnbr(i),i) = zji
              endif
            endif
          enddo
        endif
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
        write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nvbnbr(i)
      enddo
    endif
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  Neighbours of atoms :'',/)')
    endif
    call mpbarrier
    if (ioproc) then
      do i = 1,numat
        write(ioout,'(i4,8(1x,i4))') i,(nvbnbrno(nn,i),nn=1,nvbnbr(i))
      enddo
    endif
    call mpbarrier
  endif
!
!  Trap zero distances
!
  if (lzeroR) then
    call outerror('atoms are too close together',0_i4)
    call stopnow('vb_getnbr')
  endif
#ifdef TRACE
  call trace_out('vb_getnbr')
#endif
!
  end subroutine vb_getnbr

end module m_vb_nbr
