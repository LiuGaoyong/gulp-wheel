  subroutine changemaxatgfnff
!
!  Changes the size of arrays that depend on the maximum number of atoms in GFNFF arrays
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use reallocate
  use gulp_gfnff
  use m_gfnff_nbr3,    only : n3atomrptr
  use m_gfnff_c6,      only : d4_zeta_c6, d4_c9
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
  integer(i4)       :: maxat2
!
  call realloc(n3atomrptr,maxat_gfnff,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat_gfnff','n3atomrptr')
  call realloc(n_gfnff_hb_ABptr,maxat_gfnff,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat_gfnff','n_gfnff_hb_ABptr')
  call realloc(gfnff_hb_ABq,maxat_gfnff,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat_gfnff','gfnff_hb_ABq')
  call realloc(gfnff_hb_acid,maxat_gfnff,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat_gfnff','gfnff_hb_acid')
  call realloc(gfnff_hb_base,maxat_gfnff,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat_gfnff','gfnff_hb_base')
  call realloc(gfnff_xb_ABq,2_i4,maxat_gfnff,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat_gfnff','gfnff_xb_ABq')
  call realloc(n_gfnff_hb_Hrptr,maxat_gfnff,ierror)
  if (ierror.gt.0) call outofmemory('changemaxat_gfnff','n_gfnff_hb_Hrptr')
!
  maxat2 = maxat_gfnff*(maxat_gfnff+1)/2
!
  call realloc(gfnff_repulsion_p,maxat2,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat_gfnff','gfnff_repulsion_p')
  call realloc(d4_c9,maxat_gfnff,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat_gfnff','d4_c9')
  call realloc(d4_zeta_c6,maxat_gfnff,ierror)
  if (ierror.ne.0) call outofmemory('changemaxat_gfnff','d4_zeta_c6')
!
  end subroutine changemaxatgfnff

  subroutine changemaxathbH
!
!  Changes the size of arrays that hold information for hydrogen bonds - nathbH
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use gulp_gfnff
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(n_gfnff_hb_Hptr,maxathbH,ierror)
  if (ierror.ne.0) call outofmemory('changemaxathbH','n_gfnff_hb_Hptr')
!
  end subroutine changemaxathbH
!
  subroutine changemaxatxbAB
!
!  Changes the size of arrays that hold information for hydrogen bonds - natxbAB
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use gulp_gfnff
  use reallocate
  implicit none
! 
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(n_gfnff_xb_ABptr,3_i4,maxatxbAB,ierror)
  if (ierror.ne.0) call outofmemory('changemaxatxbAB','n_gfnff_xb_ABptr')
!
  end subroutine changemaxatxbAB
 
  subroutine changemaxgfnffangles
!
!  Changes the size of arrays that hold information for the angle bends
!
!  Julian Gale, CIC, Curtin University, September 2020
!
  use gulp_gfnff
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(n_gfnff_angleptr,3_i4,max_gfnff_angles,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnffangles','n_gfnff_angleptr')
  call realloc(par_gfnff_angle,2_i4,max_gfnff_angles,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnffangles','par_gfnff_angle')
!
  end subroutine changemaxgfnffangles

  subroutine changemaxgfnfftorsions
!
!  Changes the size of arrays that hold information for the torsions
!
!  Julian Gale, CIC, Curtin University, September 2020
!
  use gulp_gfnff
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(n_gfnff_torsionptr,5_i4,max_gfnff_torsions,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnfftorsions','n_gfnff_torsionptr')
  call realloc(par_gfnff_torsion,2_i4,max_gfnff_torsions,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnfftorsions','par_gfnff_torsion')
!
  end subroutine changemaxgfnfftorsions

  subroutine changemaxbond_hb
!
!  Changes the size of arrays that hold the list of trios for hydrogen bonds 3
!
!  Julian Gale, CIC, Curtin University, September 2020
!
  use gulp_gfnff
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(nbond_hb_AH,3_i4,maxbond_hb_nr,ierror)
  if (ierror.ne.0) call outofmemory('changemaxbond_hb','nbond_hb_AH')
  call realloc(nbond_hb_Bn,maxbond_hb_nr,ierror)
  if (ierror.ne.0) call outofmemory('changemaxbond_hb','nbond_hb_Bn')
  call realloc(nbond_hb_B,maxbond_hb_Bn,maxbond_hb_nr,ierror)
  if (ierror.ne.0) call outofmemory('changemaxbond_hb','nbond_hb_B')
!
  end subroutine changemaxbond_hb
!
  subroutine changemaxgfnffele
!
!  Changes the size of arrays that depend on the maximum number of elements in GFNFF arrays
!
!  Julian Gale, CIC, Curtin University, December 2021
!
  use gulp_gfnff
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror
!
  call realloc(gfnff_rad,max_gfnff_ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnffele','gfnff_rad')
  call realloc(gfnff_rad_cn,5_i4,max_gfnff_ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnffele','gfnff_rad_cn')
  call realloc(gfnff_rcov,max_gfnff_ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnffele','gfnff_rcov')
  call realloc(gfnff_repulsion_a,max_gfnff_ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnffele','gfnff_repulsion_a')
  call realloc(gfnff_repulsion_z,max_gfnff_ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnffele','gfnff_repulsion_z')
  call realloc(gfnff_xb_scale,max_gfnff_ele,ierror)
  if (ierror.ne.0) call outofmemory('changemaxgfnffele','gfnff_xb_scale')
!
  end subroutine changemaxgfnffele
