  subroutine reaxFFds(ereaxFF,lgrad1,lgrad2)
!
!  Calculates the energy and derivatives for the reaxFF force field up to second derivatives.
!  Uses sparsity to speed up calculation.
!  Parallel version. All nodes compute the bond order data for atoms connected to those that are local.
!  NB: In serial reaxFFs should be called!
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!
!  On exit :
!
!  ereaxFF         = the value of the energy contribution
!
!   8/23 Created from reaxFF.F90 and old version of reaxffds.F90
!   8/23 Correction to addressing nonzeroptr when not doing derivatives
!   9/23 Scale corrected for ebond in parallel
!  10/23 neighnoRptr now set for all atoms 
!  10/23 Bond orders now computed for all atoms to avoid problems where neighbour 
!        may not be in numattodo list for extended case
!  10/23 Handling of d1is allocation corrected
!  10/23 d2evaldthetadBO_neigh, d2etorsdthetaddeltaj_neigh and d2etorsdthetadBOpi_neigh made sparse in first dimension
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
!  Julian Gale, CIC, Curtin University, October 2023
!
  use datatypes
  use configurations, only : nregions, lsliceatom, nregiontype, lregionrigid, QMMMmode
  use control,        only : keyword, lseok, lreaxFFQ, latomicstress, literativeQ
  use g_constants,    only : evtokcal
  use current
  use derivatives,    only : xdrv, ydrv, zdrv, rstrd, xregdrv, yregdrv, zregdrv, atomicstress
  use derivatives,    only : derv2, derv3, sderv2
  use energies,       only : eattach, esregion12, esregion2, siteenergy
  use gulp_files,     only : lqbo
  use iochannels
  use m_qderiv,       only : d2charged, d2chargeselfd
  use m_strain,       only : real1strterm
  use neighbours
  use numbers,        only : third
  use optimisation,   only : lfreeze, lopf
  use parallel
  use reallocate
  use reaxFFdata
  use realvectors,    only : dist, xtmp, ytmp, ztmp
  use realvectors,    only : d0i, d0j, d2i2, d2ij, d2j2
  use realvectors,    only : d1qi => d1i
  use realvectors,    only : d1qj => d1j
  use spatialbo
  use symmetry,       only : lstr
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!  Passed variables
!
  real(dp),    intent(out)                         :: ereaxFF
  logical,     intent(in)                          :: lgrad1
  logical,     intent(in)                          :: lgrad2
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ib
  integer(i4)                                      :: ii
  integer(i4)                                      :: iloc
  integer(i4)                                      :: ind
  integer(i4)                                      :: ind0
  integer(i4)                                      :: ind3
  integer(i4)                                      :: indj
  integer(i4)                                      :: indj2
  integer(i4)                                      :: indnj2
  integer(i4)                                      :: ind1
  integer(i4)                                      :: ind1i
  integer(i4)                                      :: ind1i2
  integer(i4)                                      :: ind2
  integer(i4)                                      :: indboij2jk
  integer(i4)                                      :: indbojk
  integer(i4)                                      :: indbojk2
  integer(i4)                                      :: indbojl
  integer(i4)                                      :: indil
  integer(i4)                                      :: indjk
  integer(i4)                                      :: indjk2
  integer(i4)                                      :: indjl
  integer(i4)                                      :: indkl
  integer(i4)                                      :: indni2
  integer(i4)                                      :: indninj2
  integer(i4)                                      :: inds
  integer(i4)                                      :: indsjk
  integer(i4)                                      :: indcosphi(6)
  integer(i4)                                      :: indcosphis(6)
  integer(i4)                                      :: indsinkij(3)
  integer(i4)                                      :: indsinkijs(3)
  integer(i4)                                      :: indsinijl(3)
  integer(i4)                                      :: indsinijls(3)
  integer(i4), dimension(:,:,:), allocatable, save :: ineigh         ! Cell indices for vector
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: ixf
  integer(i4)                                      :: iyf
  integer(i4)                                      :: izf
  integer(i4)                                      :: j
  integer(i4)                                      :: j2
  integer(i4)                                      :: jj
  integer(i4)                                      :: jloc
  integer(i4)                                      :: jx
  integer(i4)                                      :: jy
  integer(i4)                                      :: jz
  integer(i4)                                      :: jxf
  integer(i4)                                      :: jyf
  integer(i4)                                      :: jzf
  integer(i4)                                      :: k
  integer(i4)                                      :: k2
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: kloc
  integer(i4)                                      :: ks
  integer(i4)                                      :: kt
  integer(i4)                                      :: kx
  integer(i4)                                      :: ky
  integer(i4)                                      :: kz
  integer(i4)                                      :: kxf
  integer(i4)                                      :: kyf
  integer(i4)                                      :: kzf
  integer(i4)                                      :: l
  integer(i4)                                      :: m
  integer(i4)                                      :: maxind
  integer(i4),                                save :: maxsparse = 1000  ! Maximum size of d1is
  integer(i4)                                      :: maxsparse2
  integer(i4)                                      :: maxneigh1
  integer(i4)                                      :: maxneigh1a
  integer(i4)                                      :: maxneigh2
  integer(i4)                                      :: maxneigh2a
  integer(i4)                                      :: maxneigh2e     ! Extended dimension
  integer(i4)                                      :: maxneigh2t     ! Up to torsion dimension
  integer(i4)                                      :: maxnonzero
  integer(i4)                                      :: maxnonzero2
  integer(i4)                                      :: mneigh
  integer(i4)                                      :: n
  integer(i4)                                      :: nati
  integer(i4)                                      :: nboij
  integer(i4)                                      :: nboij2
  integer(i4)                                      :: nboik
  integer(i4)                                      :: nbojk
  integer(i4)                                      :: nboj2k2
  integer(i4), dimension(:),     allocatable, save :: nbos
  integer(i4), dimension(:),     allocatable, save :: nbosptr
  integer(i4)                                      :: nboatom
  integer(i4), dimension(:),     allocatable, save :: nboatomptr
  integer(i4), dimension(:),     allocatable, save :: nboatomRptr
  integer(i4), dimension(:),     allocatable, save :: nbodummyptr
  integer(i4)                                      :: nbotodo
  integer(i4), dimension(:),     allocatable, save :: nbotodoptr
  integer(i4)                                      :: ni
  integer(i4)                                      :: nivalid
  integer(i4)                                      :: nj
  integer(i4)                                      :: nj2
  integer(i4)                                      :: nj2i
  integer(i4)                                      :: nk
  integer(i4)                                      :: nks
  integer(i4)                                      :: nk2
  integer(i4)                                      :: nk2j2
  integer(i4)                                      :: nl
  integer(i4)                                      :: nls
  integer(i4)                                      :: nm
  integer(i4)                                      :: nms
  integer(i4)                                      :: nmolonly
  integer(i4)                                      :: nn
  integer(i4)                                      :: nns
  integer(i4)                                      :: nor
  integer(i4), dimension(:,:),   allocatable, save :: neighno
  integer(i4), dimension(:,:),   allocatable, save :: neighnoRptr
  integer(i4), dimension(:),     allocatable, save :: nfreeatom
  integer(i4), dimension(:),     allocatable, save :: nneigh
  integer(i4)                                      :: nneighi1
  integer(i4)                                      :: nneighi1a
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: nneighi2a
  integer(i4)                                      :: nneighi2e
  integer(i4)                                      :: nneighi2t
  integer(i4)                                      :: nneighi22
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregionj
  integer(i4)                                      :: nregionk
  integer(i4)                                      :: nregiontypi
  integer(i4)                                      :: nregiontypj
  integer(i4)                                      :: ns
  integer(i4)                                      :: nspeci
  integer(i4)                                      :: nspecj
  integer(i4)                                      :: nspecj2
  integer(i4)                                      :: nspeck
  integer(i4)                                      :: nspeck2
  integer(i4)                                      :: nspecl
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: nval3
  integer(i4)                                      :: numattodo
  integer(i4)                                      :: numattodo_main
  integer(i4)                                      :: numattodo_main2
  integer(i4), dimension(:),     allocatable, save :: numattodoptr
  integer(i4), dimension(:),     allocatable, save :: numattodonptr
  integer(i4), dimension(:),     allocatable, save :: numattodoRptr
  integer(i4)                                      :: nonzero
  integer(i4)                                      :: nonzero_ij
  integer(i4)                                      :: nonzero_ij_i
  integer(i4)                                      :: nonzero_ij2
  integer(i4)                                      :: nonzero_ij2_i
  integer(i4)                                      :: nonzero_ik
  integer(i4)                                      :: nonzero_ik_i
  integer(i4)                                      :: nonzero_jk
  integer(i4)                                      :: nonzero_jk_j
  integer(i4)                                      :: nonzero_j2k2
  integer(i4)                                      :: nonzero_j2k2_j2
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr
  integer(i4), dimension(:),     allocatable, save :: nonzerorptr
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_ij
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_ik
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_ij2
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_jk
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_j2k2
  integer(i4)                                      :: oldmaxsparse
  integer(i4)                                      :: oldmaxsparse2
  integer(i4)                                      :: status
#ifdef MPI
  integer(i4)                                      :: numtmp
  integer                                          :: MPIerror
  integer                                          :: nloc
  integer                                          :: ntag
  integer                                          :: ntag2
  integer                                          :: nnode
  integer                                          :: ntmp
  integer                                          :: Request(2)
  integer,     dimension(:),     allocatable, save :: nntmp
  integer(i4), dimension(:,:),   allocatable, save :: StatMPI       ! Array for status from MPI
#endif
  logical                                          :: lAnyHB
  logical,     dimension(:),     allocatable, save :: latomdone
  logical                                          :: lattach
  logical                                          :: lanyBOnonzero
  logical                                          :: lBOijOK
  logical                                          :: lBOikOK
  logical                                          :: lBOnonzero
  logical                                          :: lDoEpen
  logical                                          :: lextended
  logical                                          :: lfixQi
  logical                                          :: lfixQj
  logical                                          :: lfound
  logical                                          :: lilocal
  logical                                          :: ljlocal
  logical                                          :: lklocal
  logical                                          :: lmaxneighok
  logical                                          :: lneeded
  logical                                          :: lQMMMelectro
  logical                                          :: lQMMMok
  logical                                          :: loverflow_expuc3
  logical                                          :: lreactivei
  logical                                          :: lreactivej
  logical                                          :: lreg2one
  logical                                          :: lreg2pair
  logical                                          :: lself
  logical                                          :: lslicei
  logical                                          :: lslicej
  logical                                          :: lStrict
  logical,     dimension(:),     allocatable, save :: lopanyneigh
  real(dp)                                         :: dtrm
  real(dp)                                         :: bo8
  real(dp)                                         :: BOij
  real(dp)                                         :: BOij_s
  real(dp)                                         :: BOij_pi
  real(dp)                                         :: BOij_pipi
  real(dp)                                         :: BOij2_s
  real(dp)                                         :: BOij2_pi
  real(dp)                                         :: BOij2_pipi
  real(dp)                                         :: BOijpbe2
  real(dp)                                         :: BOik
  real(dp)                                         :: BOik_s
  real(dp)                                         :: BOik_pi
  real(dp)                                         :: BOik_pipi
  real(dp)                                         :: BOjk_s
  real(dp)                                         :: BOjk_pi
  real(dp)                                         :: BOjk_pipi
  real(dp)                                         :: BOjk2_s
  real(dp)                                         :: BOjk2_pi
  real(dp)                                         :: BOjk2_pipi
  real(dp)                                         :: BOjl
  real(dp)                                         :: BOpij
#ifdef MPI
  real(dp),    dimension(:),     allocatable, save :: botmp
#endif
  real(dp)                                         :: conjtrm1
  real(dp)                                         :: conjtrm2
  real(dp)                                         :: coshalftheta
  real(dp)                                         :: cosp1d(6)
  real(dp)                                         :: cosp2d(21)
  real(dp)                                         :: cosphi
  real(dp)                                         :: cos2phi
  real(dp)                                         :: cos3phi
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: cut2
  real(dp)                                         :: cutq2
  real(dp)                                         :: cutvdw2
  real(dp)                                         :: d1si(6)
  real(dp)                                         :: d1ix
  real(dp)                                         :: d1iy
  real(dp)                                         :: d1iz
  real(dp)                                         :: d1sj(6)
  real(dp)                                         :: d1jx
  real(dp)                                         :: d1jy
  real(dp)                                         :: d1jz
  real(dp)                                         :: d1trm
  real(dp)                                         :: d2trm
  real(dp)                                         :: d2self
  real(dp)                                         :: dcos2phidcosphi
  real(dp)                                         :: d2cos2phidcosphi2
  real(dp)                                         :: dcos3phidcosphi
  real(dp)                                         :: d2cos3phidcosphi2
  real(dp)                                         :: delta_angle
  real(dp)                                         :: ddeltaeddeltai
  real(dp)                                         :: ddeltalpddeltai
  real(dp)                                         :: ddeltalpddeltaj
  real(dp)                                         :: ddeltalpddeltaj2
  real(dp)                                         :: ddeltalplpddeltai
  real(dp)                                         :: d2deltalplpddeltai2
  real(dp)                                         :: ddeltalpcorrdbopi
  real(dp)                                         :: ddeltalpcorrddeltai
  real(dp)                                         :: ddeltalpcorrddeltaj
  real(dp)                                         :: ddeltalpcorrddeltaj2
  real(dp)                                         :: derive0self
  real(dp)                                         :: derive0selfi
  real(dp)                                         :: derive0selfj
  real(dp)                                         :: dehbdBOij
  real(dp)                                         :: d2ehbdBOij2
  real(dp)                                         :: dehbdrij
  real(dp)                                         :: d2ehbdrij2
  real(dp)                                         :: d2ehbdrijdrik
  real(dp)                                         :: d2ehbdrijdrjk
  real(dp)                                         :: dehbdrik
  real(dp)                                         :: d2ehbdrik2
  real(dp)                                         :: d2ehbdrikdrjk
  real(dp)                                         :: dehbdrjk
  real(dp)                                         :: d2ehbdrjk2
  real(dp)                                         :: debpenddeltai
  real(dp)                                         :: d2ebpenddeltai2
  real(dp)                                         :: deloneddeltai
  real(dp)                                         :: d2eloneddeltai2
  real(dp)                                         :: deloneddeltalplp
  real(dp)                                         :: d2eloneddeltalplp2
  real(dp)                                         :: deoverddeltalpcorr
  real(dp)                                         :: deconjdcosphi
  real(dp)                                         :: d2econjdf2
  real(dp)                                         :: deconjdf12
  real(dp)                                         :: deconjdsinkij
  real(dp)                                         :: deconjdsinijl
  real(dp)                                         :: d2econjdcosphidf12
  real(dp)                                         :: d2econjdsinkijdf12
  real(dp)                                         :: d2econjdsinijldf12
  real(dp)                                         :: detorsconjdf
  real(dp)                                         :: detorsdcosphi
  real(dp)                                         :: d2etorsdcosphi2
  real(dp)                                         :: d2etorsdcosphidf
  real(dp)                                         :: d2etorsdcosphidf10
  real(dp)                                         :: d2etorsdcosphidsinkij
  real(dp)                                         :: d2etorsdcosphidsinijl
  real(dp)                                         :: d2etorsdf2
  real(dp)                                         :: detorsdcos1phi
  real(dp)                                         :: detorsdcos2phi
  real(dp)                                         :: detorsdcos3phi
  real(dp)                                         :: detorsdexpV2
  real(dp)                                         :: detorsdexpV2nof10
  real(dp)                                         :: detorsdexpV2notaper
  real(dp)                                         :: detorsdf10
  real(dp)                                         :: detorsdsinkij
  real(dp)                                         :: detorsdsinijl
  real(dp)                                         :: d2etorsdsinkijdsinijl
  real(dp)                                         :: d2etorsdsinkijdf
  real(dp)                                         :: d2etorsdsinkijdf10
  real(dp)                                         :: d2etorsdsinijldf
  real(dp)                                         :: d2etorsdsinijldf10
  real(dp)                                         :: deunderddeltalpcorr
  real(dp)                                         :: deunderdsumdeltabopi
  real(dp)                                         :: d2eunderdsumdeltabopi2
  real(dp)                                         :: d2eunderddeltalpcorrdsumdeltabopi
  real(dp)                                         :: d2deltalpcorrddeltaiddeltaj
  real(dp)                                         :: d2deltalpcorrddeltaj2
  real(dp)                                         :: d2deltalpcorrddeltajddeltaj2
  real(dp)                                         :: d2deltalpcorrddeltaidbopi
  real(dp)                                         :: d2deltalpcorrddeltaidsumdeltabopi
  real(dp)                                         :: d2deltalpcorrddeltajdbopi
  real(dp)                                         :: d2eunderddeltalpcorr2
  real(dp)                                         :: devalddeltai
  real(dp)                                         :: d2evalddeltai2
  real(dp)                                         :: d2evalddeltaidsbo
  real(dp)                                         :: d2evaldthetadtheta0
  real(dp)                                         :: d2evaldf7dtheta
  real(dp)                                         :: d2evaldf7dtheta0
  real(dp)                                         :: d2evaldfijdtheta
  real(dp)                                         :: d2evaldfikdtheta
  real(dp)                                         :: d2evaldfijdtheta0
  real(dp)                                         :: d2evaldfikdtheta0
  real(dp)                                         :: d2evaldthetaddeltai
  real(dp)                                         :: devaldsbo
  real(dp)                                         :: d2evaldsbo2
  real(dp)                                         :: devdwdrij
  real(dp)                                         :: d2evdwdrij2
  real(dp)                                         :: diffBOdelta
  real(dp)                                         :: dtheta0dsbo
  real(dp)                                         :: dthetadrij
  real(dp)                                         :: dthetadrik
  real(dp)                                         :: dthetadrjk
  real(dp)                                         :: devaldtheta
  real(dp)                                         :: devaldtheta0
  real(dp)                                         :: d2evaldtheta2
  real(dp)                                         :: d2evaldtheta02
  real(dp)                                         :: dsboddeltai
  real(dp)                                         :: d2sboddeltai2
  real(dp)                                         :: dr2ds(6,3)
  real(dp)                                         :: d2r2dx2(6,3)
  real(dp)                                         :: d2r2ds2(6,6,3)
  real(dp)                                         :: d2r2dsdx(6,3,3)
  real(dp)                                         :: dexpuc1ddeltalpcorr
  real(dp)                                         :: d2expuc1ddeltalpcorr2
  real(dp)                                         :: dutrm1ddeltalpcorr
  real(dp)                                         :: d2utrm1ddeltalpcorr2
  real(dp)                                         :: dutrm2dsumdeltabopi
  real(dp)                                         :: d2utrm2dsumdeltabopi2
  real(dp)                                         :: ebond
  real(dp)                                         :: ebpen
  real(dp)                                         :: ebpenij
  real(dp)                                         :: ecoa
  real(dp)                                         :: ecoatrm
  real(dp)                                         :: econj
  real(dp)                                         :: ecoul
  real(dp)                                         :: ehb
  real(dp)                                         :: elone
  real(dp)                                         :: eover
  real(dp)                                         :: epen
  real(dp)                                         :: epentrm
  real(dp)                                         :: eself
  real(dp)                                         :: esite
  real(dp)                                         :: etors
  real(dp)                                         :: etrm
  real(dp)                                         :: eunder
  real(dp)                                         :: eval
  real(dp)                                         :: evdw
  real(dp)                                         :: expH
  real(dp)                                         :: exptheta
  real(dp)                                         :: fangle
  real(dp)                                         :: dfangledtheta
  real(dp)                                         :: dfangledtheta0
  real(dp)                                         :: d2fangledtheta2
  real(dp)                                         :: d2fangledthetadtheta0
  real(dp)                                         :: d2fangledtheta02
  real(dp),    dimension(:,:),   allocatable, save :: BO
  real(dp),    dimension(:,:),   allocatable, save :: BO_pi
  real(dp),    dimension(:,:),   allocatable, save :: BO_pipi
  real(dp),    dimension(:,:),   allocatable, save :: BOp
  real(dp),    dimension(:,:),   allocatable, save :: BOp_pi
  real(dp),    dimension(:,:),   allocatable, save :: BOp_pipi
  real(dp),    dimension(:,:),   allocatable, save :: d1BOp
  real(dp),    dimension(:,:),   allocatable, save :: d1BOp_pi
  real(dp),    dimension(:,:),   allocatable, save :: d1BOp_pipi
  real(dp),    dimension(:),     allocatable, save :: d2mix
  real(dp),    dimension(:),     allocatable, save :: d1BOi
  real(dp),    dimension(:),     allocatable, save :: d1BOi_s
  real(dp),    dimension(:),     allocatable, save :: d1BOi_pi
  real(dp),    dimension(:),     allocatable, save :: d1BOi_pipi
  real(dp),    dimension(:),     allocatable, save :: d1BOij2
  real(dp),    dimension(:),     allocatable, save :: d1BOij2_s
  real(dp),    dimension(:),     allocatable, save :: d1BOij2_pi
  real(dp),    dimension(:),     allocatable, save :: d1BOij2_pipi
  real(dp),    dimension(:),     allocatable, save :: d1BOik
  real(dp),    dimension(:),     allocatable, save :: d1BOik_s
  real(dp),    dimension(:),     allocatable, save :: d1BOik_pi
  real(dp),    dimension(:),     allocatable, save :: d1BOik_pipi
  real(dp),    dimension(:),     allocatable, save :: d1BOjk
  real(dp),    dimension(:),     allocatable, save :: d1BOjk_s
  real(dp),    dimension(:),     allocatable, save :: d1BOjk_pi
  real(dp),    dimension(:),     allocatable, save :: d1BOjk_pipi
  real(dp),    dimension(:),     allocatable, save :: d1BOjk2
  real(dp),    dimension(:),     allocatable, save :: d1BOjk2_s
  real(dp),    dimension(:),     allocatable, save :: d1BOjk2_pi
  real(dp),    dimension(:),     allocatable, save :: d1BOjk2_pipi
  real(dp),    dimension(:,:),   allocatable, save :: d2BOp
  real(dp),    dimension(:,:),   allocatable, save :: d2BOp_pi
  real(dp),    dimension(:,:),   allocatable, save :: d2BOp_pipi
  real(dp),    dimension(:),     allocatable, save :: d2BOi
  real(dp),    dimension(:),     allocatable, save :: d2BOi_s
  real(dp),    dimension(:),     allocatable, save :: d2BOi_pi
  real(dp),    dimension(:),     allocatable, save :: d2BOi_pipi
  real(dp),    dimension(:),     allocatable, save :: d2BOjk
  real(dp),    dimension(:),     allocatable, save :: d2BOjk_s
  real(dp),    dimension(:),     allocatable, save :: d2BOjk_pi
  real(dp),    dimension(:),     allocatable, save :: d2BOjk_pipi
  real(dp),    dimension(:),     allocatable, save :: delta
  real(dp),    dimension(:),     allocatable, save :: deltap
  real(dp),    dimension(:),     allocatable, save :: deltalplp ! delta lp for lone pairs
  real(dp),    dimension(:),     allocatable, save :: deltalp   ! delta lp for over coordination
  real(dp),    dimension(:),     allocatable, save :: debpendBO_neigh
  real(dp),    dimension(:),     allocatable, save :: d2ebpendBO2_neigh
  real(dp),    dimension(:),     allocatable, save :: d2ebpenddeltaidBO_neigh
  real(dp),    dimension(:),     allocatable, save :: devaldBO_neigh
  real(dp),    dimension(:),     allocatable, save :: d2evaldBO2_neigh
  real(dp),    dimension(:),     allocatable, save :: d2evalddeltaidBO_neigh
  real(dp),    dimension(:),     allocatable, save :: d2evalddeltaiddeltaj_neigh
  real(dp),    dimension(:),     allocatable, save :: d2evalddeltajk_neigh
  real(dp),    dimension(:),     allocatable, save :: detorsdBOpi_neigh
  real(dp),    dimension(:),     allocatable, save :: devalddeltaj_neigh
  real(dp),    dimension(:),     allocatable, save :: devaldsbo_neigh
  real(dp),    dimension(:),     allocatable, save :: d2evaldsbo_neigh
  real(dp),    dimension(:),     allocatable, save :: d2evaldsbo2_neigh
  real(dp),    dimension(:),     allocatable, save :: d2evaldthetaddeltai_neigh
  real(dp),    dimension(:),     allocatable, save :: d2evaldthetadBOpi_neigh
  real(dp),    dimension(:),     allocatable, save :: d2etorsdBOpi2_neigh
  real(dp),    dimension(:),     allocatable, save :: d2etorsddeltadBOpi_neigh
  real(dp),    dimension(:),     allocatable, save :: d2etorsddeltaidBOpi_neigh
  real(dp),    dimension(:),     allocatable, save :: d2etorsddeltajdBOpi_neigh
  real(dp),    dimension(:,:),   allocatable, save :: d2evalddeltajdBO_neigh
  real(dp),    dimension(:,:),   allocatable, save :: d2etorsdBOdBOpi_neigh
  real(dp),    dimension(:),     allocatable, save :: d2etorsdthetaddeltai
  real(dp)                                         :: dBOpijdr
  real(dp)                                         :: d2BOpijdr2
  real(dp)                                         :: deltapi
  real(dp)                                         :: deltapj
  real(dp)                                         :: deltapj2
  real(dp)                                         :: deltapk
  real(dp)                                         :: deltapk2
  real(dp)                                         :: delta_coa
  real(dp)                                         :: delta_e
  real(dp)                                         :: delta_e0
  real(dp)                                         :: deltalpcorr_i
  real(dp)                                         :: deijdBO_s
  real(dp)                                         :: d2eijdBO2_s
  real(dp)                                         :: deijdBO_pi
  real(dp)                                         :: deijdBO_pipi
  real(dp)                                         :: dexpV2dBO_pi
  real(dp)                                         :: d2expV2dBO_pi2
  real(dp)                                         :: d2expV2df11dBO_pi
  real(dp)                                         :: d2expV2df112
  real(dp)                                         :: dexpV2df11
  real(dp)                                         :: df7dBOij
  real(dp)                                         :: df7dBOik
  real(dp)                                         :: d2f7dBOij2
  real(dp)                                         :: d2f7dBOijdBOik
  real(dp)                                         :: d2f7dBOik2
  real(dp)                                         :: df8ddeltai
  real(dp)                                         :: d2f8ddeltai2
  real(dp)                                         :: df9ddeltai
  real(dp)                                         :: d2f9ddeltai2
  real(dp)                                         :: df10dBOij
  real(dp)                                         :: df10dBOik
  real(dp)                                         :: df10dBOjl
  real(dp)                                         :: d2f10dBOij2
  real(dp)                                         :: d2f10dBOijdBOjl
  real(dp)                                         :: d2f10dBOik2
  real(dp)                                         :: d2f10dBOikdBOij
  real(dp)                                         :: d2f10dBOikdBOjl
  real(dp)                                         :: d2f10dBOjl2
  real(dp)                                         :: df11ddeltai
  real(dp)                                         :: df11ddeltaj
  real(dp)                                         :: df11ddeltaij
  real(dp)                                         :: d2f11ddeltai2
  real(dp)                                         :: d2f11ddeltaij
  real(dp)                                         :: d2f11ddeltaij2
  real(dp)                                         :: d2f11ddeltaj2
  real(dp)                                         :: df12dBOij
  real(dp)                                         :: df12dBOik
  real(dp)                                         :: df12dBOjl
  real(dp)                                         :: d2f12dBOik2
  real(dp)                                         :: d2f12dBOikdBOij
  real(dp)                                         :: d2f12dBOikdBOjl
  real(dp)                                         :: d2f12dBOij2
  real(dp)                                         :: d2f12dBOijdBOjl
  real(dp)                                         :: d2f12dBOjl2
  real(dp)                                         :: df13drij
  real(dp)                                         :: d2f13drij2
  real(dp)                                         :: fij
  real(dp)                                         :: dfijdBOij
  real(dp)                                         :: d2fijdBOij2
  real(dp)                                         :: d3fijdBOij3
  real(dp)                                         :: fik
  real(dp)                                         :: dfikdBOik
  real(dp)                                         :: d2fikdBOik2
  real(dp)                                         :: d3fikdBOik3
  real(dp)                                         :: fjl
  real(dp)                                         :: dfjldBOjl
  real(dp)                                         :: d2fjldBOjl2
  real(dp)                                         :: d3fjldBOjl3
  real(dp)                                         :: fHB
  real(dp)                                         :: dfHBdBO
  real(dp)                                         :: d2fHBdBO2
  real(dp)                                         :: d3fHBdBO3
  real(dp)                                         :: frHB
  real(dp)                                         :: dfrHBdr
  real(dp)                                         :: d2frHBdr2
  real(dp)                                         :: d3frHBdr3
  real(dp)                                         :: dsinkijdr(3)
  real(dp)                                         :: dsinijldr(3)
  real(dp)                                         :: d2sinkijdr2(6)
  real(dp)                                         :: d2sinijldr2(6)
  real(dp)                                         :: d2thetadr2(6)
  real(dp)                                         :: d2theta0dsbo2
  real(dp)                                         :: dVDWtrmdf13
  real(dp)                                         :: dVDWtrmdrij
  real(dp)                                         :: d2VDWtrmdrij2
  real(dp)                                         :: d2deltalpcorrddeltai2
  real(dp)                                         :: d2deltalpcorrdbopi2
  real(dp)                                         :: d2deltalpddeltai2
  real(dp)                                         :: d2deltalpddeltaj2
  real(dp)                                         :: deoverdsumover
  real(dp)                                         :: d2eoverddeltalpcorr2
  real(dp)                                         :: d2eoverdsumoverddeltalpcorr
  real(dp)                                         :: detotdbo_tot
  real(dp)                                         :: detotdbopi
  real(dp)                                         :: detotdbopi_tot
  real(dp)                                         :: detotddeltaj
  real(dp)                                         :: detotddeltaj_tot
  real(dp)                                         :: detotddeltalpcorr
  real(dp)                                         :: d2etotdbo2ik
  real(dp)                                         :: d2etotdbo2_tot
  real(dp)                                         :: d2etotdbodbopi_1
  real(dp)                                         :: d2etotdbodbopi_2
  real(dp)                                         :: d2etotdbodbopi_tot
  real(dp)                                         :: d2etotdbodsumdeltabopi
  real(dp)                                         :: d2etotdbopi2
  real(dp)                                         :: d2etotdbopi2_tot
  real(dp)                                         :: d2etotdbopi2jk
  real(dp)                                         :: d2etotddeltaj2
  real(dp)                                         :: d2etotddeltajdbojk
  real(dp)                                         :: d2etotddeltaj2_tot
  real(dp)                                         :: d2etotddeltajdbo_tot
  real(dp)                                         :: d2etotddeltajdbopijk
  real(dp)                                         :: d2etotddeltajddeltaj2
  real(dp)                                         :: d2etotddeltajddeltaj2_tot
  real(dp)                                         :: d2etotddeltalpcorr2
  real(dp)                                         :: d2etotdsumdeltabopi2_tot
  real(dp)                                         :: dsumdeltabopidbopi
  real(dp)                                         :: dsumdeltabopidbopij2
  real(dp)                                         :: dsumdeltabopidbopik
  real(dp)                                         :: dsumdeltabopiddeltaj
  real(dp)                                         :: dsumdeltabopiddeltaj2
  real(dp)                                         :: d2sumdeltabopiddeltaj2
  real(dp)                                         :: d2sumdeltabopiddeltajdbopi
  real(dp)                                         :: eij
  real(dp)                                         :: expbo8
  real(dp)                                         :: expcoa2
  real(dp)                                         :: expcoa3a
  real(dp)                                         :: d1expcoa3adBOij
  real(dp)                                         :: d1expcoa3addeltaj
  real(dp)                                         :: d2expcoa3addeltaj2
  real(dp)                                         :: d2expcoa3addeltajdBOij
  real(dp)                                         :: d2expcoa3adBOij2
  real(dp)                                         :: d2expcoa3abdBOdBO
  real(dp)                                         :: d2expcoa3abddeltadBO
  real(dp)                                         :: d2expcoa3abddeltajddeltak
  real(dp)                                         :: d2expcoa3a4adBOij
  real(dp)                                         :: d2expcoa3a4addeltaj
  real(dp)                                         :: d2expcoa3a4bdBOij
  real(dp)                                         :: d2expcoa3a4bddeltaj
  real(dp)                                         :: expcoa3b
  real(dp)                                         :: d1expcoa3bdBOik
  real(dp)                                         :: d2expcoa3bdBOik2
  real(dp)                                         :: d1expcoa3bddeltak
  real(dp)                                         :: d2expcoa3bddeltak2
  real(dp)                                         :: d2expcoa3bddeltakdBOik
  real(dp)                                         :: d2expcoa3b4adBOik
  real(dp)                                         :: d2expcoa3b4addeltak
  real(dp)                                         :: d2expcoa3b4bdBOik
  real(dp)                                         :: d2expcoa3b4bddeltak
  real(dp)                                         :: expcoa4a
  real(dp)                                         :: d1expcoa4a
  real(dp)                                         :: d2expcoa4a
  real(dp)                                         :: d2expcoa4ab
  real(dp)                                         :: expcoa4b
  real(dp)                                         :: d1expcoa4b
  real(dp)                                         :: d2expcoa4b
  real(dp)                                         :: expcoaprod
  real(dp)                                         :: expcoatrm
  real(dp)                                         :: dexpcoatrmddeltai
  real(dp)                                         :: d2expcoatrmddeltai2
  real(dp)                                         :: expdi
  real(dp)                                         :: exphb1
  real(dp)                                         :: dexphb1dBOij
  real(dp)                                         :: d2exphb1dBOij2
  real(dp)                                         :: exphb2
  real(dp)                                         :: dexphb2drik
  real(dp)                                         :: d2exphb2drik2
  real(dp)                                         :: expij
  real(dp)                                         :: explp
  real(dp)                                         :: explp2
  real(dp)                                         :: expoc
  real(dp)                                         :: exppen1
  real(dp)                                         :: d1exppen1
  real(dp)                                         :: d2exppen1
  real(dp)                                         :: exppen2
  real(dp)                                         :: d1exppen2
  real(dp)                                         :: d2exppen2
  real(dp)                                         :: expuc1
  real(dp)                                         :: expuc2
  real(dp)                                         :: expuc3
  real(dp)                                         :: expV2
  real(dp)                                         :: expVDW1
  real(dp)                                         :: expVDW2
  real(dp)                                         :: f7
  real(dp)                                         :: f8
  real(dp)                                         :: f9
  real(dp)                                         :: f10
  real(dp)                                         :: f11
  real(dp)                                         :: f12
  real(dp)                                         :: f13
  real(dp)                                         :: decouldrij
  real(dp)                                         :: d2ecouldrij2
  real(dp)                                         :: detrmdrij
  real(dp)                                         :: d2etrmdrij2
  real(dp)                                         :: dgamdrij
  real(dp)                                         :: d2gamdrij2
  real(dp)                                         :: dgamnoqdrij
  real(dp)                                         :: gam
  real(dp)                                         :: gamnoq
  real(dp)                                         :: gammai
  real(dp)                                         :: gammaj
  real(dp)                                         :: gammaij
  real(dp)                                         :: half
  real(dp)                                         :: otrm1
  real(dp)                                         :: otrm2
  real(dp)                                         :: dotrm2ddeltalpcorr
  real(dp)                                         :: d2otrm2ddeltalpcorr2
  real(dp)                                         :: otrm3
  real(dp)                                         :: dotrm3ddeltalpcorr
  real(dp)                                         :: d2otrm3ddeltalpcorr2
  real(dp)                                         :: ddeltalpcorrdsumdeltabopi
  real(dp)                                         :: d2deltalpcorrdsumdeltabopi2
  real(dp)                                         :: ppen2
  real(dp)                                         :: pcoa1
  real(dp)                                         :: pcoa2
  real(dp)                                         :: pcoa3
  real(dp)                                         :: pcoa4
  real(dp)                                         :: pcot1
  real(dp)                                         :: phb1
  real(dp)                                         :: phb2
  real(dp)                                         :: phb3
  real(dp)                                         :: ptor1
  real(dp)                                         :: pval1
  real(dp)                                         :: pval2
  real(dp)                                         :: qij
  real(dp)                                         :: r0hb
  real(dp)                                         :: reaxFFrhtollower
  real(dp)                                         :: utrm1
  real(dp)                                         :: utrm2
  real(dp)                                         :: rintde
  real(dp)                                         :: drintdeddeltae
  real(dp)                                         :: d2rintdeddeltae2
  real(dp)                                         :: rij
  real(dp)                                         :: rij2
  real(dp)                                         :: ril
  real(dp)                                         :: ril2
  real(dp)                                         :: rkl
  real(dp)                                         :: rkl2
  real(dp)                                         :: rji
  real(dp)                                         :: rjk
  real(dp)                                         :: rjk2
  real(dp)                                         :: rki
  real(dp)                                         :: rki2
  real(dp)                                         :: rkj
  real(dp)                                         :: rkj2
  real(dp)                                         :: rrij
  real(dp)                                         :: rrik
  real(dp)                                         :: rrjk
  real(dp)                                         :: rnlp
  real(dp)                                         :: rstrdloc(6)
  real(dp)                                         :: sbo
  real(dp)                                         :: sbopi
  real(dp)                                         :: sboprod
  real(dp)                                         :: scale
  real(dp)                                         :: sinkij
  real(dp)                                         :: sinijl
  real(dp)                                         :: sinhalftheta
  real(dp)                                         :: sin4halftheta
  real(dp)                                         :: dsin4halfthetadtheta
  real(dp)                                         :: d2sin4halfthetadtheta2
  real(dp)                                         :: smoothH
  real(dp)                                         :: sumBOcoa_ij
  real(dp)                                         :: sumBOcoa_ik
  real(dp)                                         :: sumdeltabopi
  real(dp)                                         :: sumover
  real(dp)                                         :: dsumoverdbo
  real(dp)                                         :: dsumoverdboik
  real(dp)                                         :: dsumoverdboj2
  real(dp)                                         :: ctwo
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp)                                         :: t3
  real(dp)                                         :: t4
  real(dp)                                         :: theta
  real(dp)                                         :: theta0
  real(dp)                                         :: tp
  real(dp)                                         :: dtpdr
  real(dp)                                         :: d2tpdr2
  real(dp)                                         :: tpQ
  real(dp)                                         :: dtpQdr
  real(dp)                                         :: d2tpQdr2
  real(dp)                                         :: trm1lp
  real(dp)                                         :: trm2lp
  real(dp)                                         :: V1
  real(dp)                                         :: V1trm
  real(dp)                                         :: V2
  real(dp)                                         :: V2trm
  real(dp)                                         :: V3
  real(dp)                                         :: V3trm
  real(dp)                                         :: VDWtrm
  real(dp),    dimension(:,:),   allocatable, save :: rneigh
  real(dp),    dimension(:,:),   allocatable, save :: xneigh
  real(dp),    dimension(:,:),   allocatable, save :: yneigh
  real(dp),    dimension(:,:),   allocatable, save :: zneigh
  real(dp),    dimension(:),     allocatable, save :: sum
  real(dp),    dimension(:),     allocatable, save :: sum0
  real(dp),    dimension(:),     allocatable, save :: sumsij
  real(dp),    dimension(:),     allocatable, save :: dsboproddbo
  real(dp),    dimension(:),     allocatable, save :: d2sboproddbo2
  real(dp),    dimension(:),     allocatable, save :: d2sboddeltaidbo
  real(dp)                                         :: xd
  real(dp)                                         :: yd
  real(dp)                                         :: zd
  real(dp)                                         :: xdiff
  real(dp)                                         :: ydiff
  real(dp)                                         :: zdiff
  real(dp)                                         :: xil
  real(dp)                                         :: yil
  real(dp)                                         :: zil
  real(dp)                                         :: xkl
  real(dp)                                         :: ykl
  real(dp)                                         :: zkl
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xjk
  real(dp)                                         :: yjk
  real(dp)                                         :: zjk
  real(dp)                                         :: xki
  real(dp)                                         :: yki
  real(dp)                                         :: zki
  real(dp)                                         :: xkj
  real(dp)                                         :: ykj
  real(dp)                                         :: zkj
!   
!  Sparse arrays - use pointer since size might change during calculation
!   
  real(dp),    dimension(:),     pointer,     save :: d1is => null()
  real(dp),    dimension(:),     pointer,     save :: d2is => null()
  real(dp),    dimension(:,:),   pointer,     save :: d2evaldthetadBO_neigh => null()
  real(dp),    dimension(:,:),   pointer,     save :: d2etorsdthetaddeltaj_neigh => null()
  real(dp),    dimension(:,:),   pointer,     save :: d2etorsdthetadBOpi_neigh => null()
#ifdef TRACE
  call trace_in('reaxffds')
#endif
!
  t1 = g_cpu_time()
! 
!  Check this is a parallel run
! 
  if (nprocs.eq.1) then
    call outerror('reaxFFds called in serial!',0_i4)
    call stopnow('reaxFFds')
  endif
!
!  If lStrict is true then we try to set the behaviour to be the same as the
!  original ReaxFF code and turn off the modifications implemented for smoothness
!
  lStrict = (index(keyword,'stri').eq.1.or.index(keyword,' stri').ne.0)
!
!  Set lower bounds for HB taper to be 0.9 of upper bound for now
!
  reaxFFrhtollower = 0.9_dp*reaxFFrhtol
!
!  Find out whether any pair requires the use of extended derivatives
!
  lextended = .false.
  maxind = nreaxFFspec*(nreaxFFspec+1)/2
  ind = 0
  do while (.not.lextended.and.ind.lt.maxind)
    ind = ind + 1
    if (lreaxFFbocorrect(1,ind).or.lreaxFFbocorrect(2,ind)) lextended = .true.
  enddo
!
  allocate(nboatomptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','nboatomptr')
  allocate(nboatomRptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','nboatomRptr')
  allocate(nbodummyptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','nbodummyptr')
!
!  Allocate memory that does not depend on maxneigh
!
  allocate(numattodoptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','numattodoptr')
  allocate(numattodonptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','numattodonptr')
  allocate(numattodoRptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','numattodoRptr')
  allocate(nbos(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','nbos')
  allocate(nbosptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','nbosptr')
  allocate(latomdone(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','latomdone')
  allocate(lopanyneigh(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','lopanyneigh')
  allocate(nfreeatom(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','nfreeatom')
  allocate(nneigh(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','nneigh')
  allocate(sumsij(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','sumsij')
  allocate(sum(max(3*numat,14_i4)),stat=status)
  if (status/=0) call outofmemory('reaxFFds','sum')
  allocate(sum0(max(3*numat,14_i4)),stat=status)
  if (status/=0) call outofmemory('reaxFFds','sum0')
!
!  Reinitialisation point should maxneigh be increased
!
100 continue
  lmaxneighok = .true.
  if (allocated(rneigh)) then
    deallocate(d2sboddeltaidbo,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2sboddeltaidbo')
    deallocate(d2sboproddbo2,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2sboproddbo2')
    deallocate(dsboproddbo,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','dsboproddbo')
    deallocate(nonzeroptr_j2k2,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_j2k2')
    deallocate(nonzeroptr_jk,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_jk')
    deallocate(nonzeroptr_ij2,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_ij2')
    deallocate(nonzeroptr_ik,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_ik')
    deallocate(nonzeroptr_ij,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_ij')
    deallocate(nonzeroptr,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','nonzeroptr')
    deallocate(nonzerorptr,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','nonzerorptr')
    deallocate(d1BOjk2_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOjk2_pipi')
    deallocate(d1BOjk2_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOjk2_pi')
    deallocate(d1BOjk2_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOjk2_s')
    deallocate(d1BOjk2,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOjk2')
    deallocate(d1BOjk_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOjk_pipi')
    deallocate(d1BOjk_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOjk_pi')
    deallocate(d1BOjk_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOjk_s')
    deallocate(d1BOjk,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOjk')
    deallocate(d2BOp_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOp_pipi')
    deallocate(d2BOp_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOp_pi')
    deallocate(d2BOp,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOp')
    deallocate(d2BOjk_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOjk_pipi')
    deallocate(d2BOjk_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOjk_pi')
    deallocate(d2BOjk_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOjk_s')
    deallocate(d2BOjk,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOjk')
    deallocate(d2BOi_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOi_pipi')
    deallocate(d2BOi_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOi_pi')
    deallocate(d2BOi_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOi_s')
    deallocate(d2BOi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2BOi')
    deallocate(d2mix,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2mix')
    deallocate(d1BOik_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOik_pipi')
    deallocate(d1BOik_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOik_pi')
    deallocate(d1BOik_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOik_s')
    deallocate(d1BOik,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOik')
    deallocate(d1BOij2_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOij2_pipi')
    deallocate(d1BOij2_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOij2_pi')
    deallocate(d1BOij2_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOij2_s')
    deallocate(d1BOij2,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOij2')
    deallocate(d1BOi_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOi_pipi')
    deallocate(d1BOi_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOi_pi')
    deallocate(d1BOi_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOi_s')
    deallocate(d1BOi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOi')
    deallocate(devaldsbo_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','devaldsbo_neigh')
    deallocate(devalddeltaj_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','devalddeltaj_neigh')
    deallocate(d2etorsdthetaddeltai,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2etorsdthetaddeltai')
    deallocate(d2etorsdBOdBOpi_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2etorsdBOdBOpi_neigh')
    deallocate(d2etorsdBOpi2_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2etorsdBOpi2_neigh')
    deallocate(detorsdBOpi_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','detorsdBOpi_neigh')
    deallocate(d2evaldthetaddeltai_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2evaldthetaddeltai_neigh')
    deallocate(d2evaldthetadBOpi_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2evaldthetadBOpi_neigh')
    deallocate(d2evalddeltajdBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2evalddeltajdBO_neigh')
    deallocate(d2etorsddeltajdBOpi_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2etorsddeltajdBOpi_neigh')
    deallocate(d2etorsddeltaidBOpi_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2etorsddeltaidBOpi_neigh')
    deallocate(d2etorsddeltadBOpi_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2etorsddeltadBOpi_neigh')
    deallocate(d2evaldsbo2_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2evaldsbo2_neigh')
    deallocate(d2evaldsbo_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2evaldsbo_neigh')
    deallocate(d2evalddeltajk_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2evalddeltajk_neigh')
    deallocate(d2evalddeltaiddeltaj_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2evalddeltaiddeltaj_neigh')
    deallocate(d2evalddeltaidBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2evalddeltaidBO_neigh')
    deallocate(d2evaldBO2_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2evaldBO2_neigh')
    deallocate(devaldBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','devaldBO_neigh')
    deallocate(d2ebpenddeltaidBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2ebpendddeltaiBO_neigh')
    deallocate(d2ebpendBO2_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d2ebpendBO2_neigh')
    deallocate(debpendBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','debpendBO_neigh')
    deallocate(d1BOp_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOp_pipi')
    deallocate(d1BOp_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOp_pi')
    deallocate(d1BOp,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','d1BOp')
    deallocate(deltalplp,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','deltalplp')
    deallocate(deltalp,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','deltalp')
    deallocate(deltap,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','deltap')
    deallocate(delta,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','delta')
    deallocate(BOp_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','BOp_pipi')
    deallocate(BOp_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','BOp_pi')
    deallocate(BOp,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','BOp')
    deallocate(BO_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','BO_pipi')
    deallocate(BO_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','BO_pi')
    deallocate(BO,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','BO')
    deallocate(zneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','zneigh')
    deallocate(yneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','yneigh')
    deallocate(xneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','xneigh')
    deallocate(rneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','rneigh')
    deallocate(ineigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','ineigh')
    deallocate(neighnoRptr,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','neighnoRptr')
    deallocate(neighno,stat=status)
    if (status/=0) call deallocate_error('reaxFFds','neighno')
  endif
!
!  Initialise Bond Order energy contributions
!
  ebond  = 0.0_dp
  ebpen  = 0.0_dp
  ecoa   = 0.0_dp
  econj  = 0.0_dp
  ecoul  = 0.0_dp
  elone  = 0.0_dp
  eover  = 0.0_dp
  epen   = 0.0_dp
  eself  = 0.0_dp
  etors  = 0.0_dp
  eunder = 0.0_dp
  eval   = 0.0_dp
  evdw   = 0.0_dp
  ehb    = 0.0_dp
!
!  Set parameter for pairwise storage memory
!
!  maxneigh2  is the size of first derivative arrays except d1i
!  maxneigh2e is the size of first derivative array d1i - this is larger to hold torsional derivatives (and extended if required)
!  maxneigh2t is the size of first derivative array d1i - this is larger to hold torsional derivatives
!  maxneigh1a is the size of devaldBO_neigh
!
  maxneigh1   = maxneigh*(maxneigh + 1)/2
  maxneigh1a  = 2*maxneigh1
  maxneigh2   = maxneigh + maxneigh1
  maxneigh2a  = maxneigh1a*(maxneigh1a + 1)/2
  if (lextended) then
    maxneigh2e  = maxneigh + maxneigh1 + maxneigh*(2*maxneigh + 1)*maxneigh
    maxneigh2t  = maxneigh + maxneigh1 + maxneigh*(maxneigh + 1)*maxneigh
  else
    maxneigh2e  = maxneigh + maxneigh1 + maxneigh*(maxneigh + 1)*maxneigh
    maxneigh2t  = maxneigh + maxneigh1 + maxneigh*(maxneigh + 1)*maxneigh
  endif
! 
!  Set sizes for arrays that are passed to reaxFF_bos 
! 
  maxnonzero  = 2*maxneigh                         
  maxnonzero2 = maxnonzero*(maxnonzero+1)/2
!
!  Allocate local memory
!
  allocate(neighno(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','neighno')
  allocate(neighnoRptr(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','neighnoRptr')
  allocate(ineigh(3_i4,maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','ineigh')
  allocate(rneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','rneigh')
  allocate(xneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','xneigh')
  allocate(yneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','yneigh')
  allocate(zneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','zneigh')
  allocate(BO(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','BO')
  allocate(BO_pi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','BO_pi')
  allocate(BO_pipi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','BO_pipi')
  allocate(BOp(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','BOp')
  allocate(BOp_pi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','BOp_pi')
  allocate(BOp_pipi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','BOp_pipi')
  allocate(delta(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','delta')
  allocate(deltap(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','deltap')
  allocate(deltalp(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','deltalp')
  allocate(deltalplp(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','deltalplp')
  if (lgrad1) then
    allocate(nonzerorptr(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzerorptr')
    allocate(nonzeroptr(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr')
    allocate(nonzeroptr_ij(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_ij')
    allocate(nonzeroptr_ij2(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_ij2')
    allocate(nonzeroptr_ik(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_ik')
    allocate(nonzeroptr_jk(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_jk')
    allocate(nonzeroptr_j2k2(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_j2k2')
    allocate(d1BOp(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOp')
    allocate(d1BOp_pi(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOp_pi')
    allocate(d1BOp_pipi(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOp_pipi')
    allocate(debpendBO_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','debpendBO_neigh')
    allocate(devaldBO_neigh(maxneigh1a),stat=status)
    if (status/=0) call outofmemory('reaxFFds','devaldBO_neigh')
    allocate(detorsdBOpi_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','detorsdBOpi_neigh')
    allocate(devalddeltaj_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','devalddeltaj_neigh')
    allocate(devaldsbo_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','devaldsbo_neigh')
    allocate(d1BOi(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOi')
    allocate(d1BOi_s(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOi_s')
    allocate(d1BOi_pi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOi_pi')
    allocate(d1BOi_pipi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOi_pipi')
    allocate(d1BOij2(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOij2')
    allocate(d1BOij2_s(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOij2_s')
    allocate(d1BOij2_pi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOij2_pi')
    allocate(d1BOij2_pipi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOij2_pipi')
    allocate(d1BOik(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOik')
    allocate(d1BOik_s(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOik_s')
    allocate(d1BOik_pi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOik_pi')
    allocate(d1BOik_pipi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOik_pipi')
    allocate(d1BOjk(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk')
    allocate(d1BOjk_s(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk_s')
    allocate(d1BOjk_pi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk_pi')
    allocate(d1BOjk_pipi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk_pipi')
    allocate(d1BOjk2(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk2')
    allocate(d1BOjk2_s(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk2_s')
    allocate(d1BOjk2_pi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk2_pi')
    allocate(d1BOjk2_pipi(maxnonzero),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk2_pipi')
    allocate(dsboproddbo(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','dsboproddbo')
  else
    allocate(nonzerorptr(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzerorptr')
    allocate(nonzeroptr(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr')
    allocate(nonzeroptr_ij(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_ij')
    allocate(nonzeroptr_ij2(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_ij2')
    allocate(nonzeroptr_ik(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_ik')
    allocate(nonzeroptr_jk(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_jk')
    allocate(nonzeroptr_j2k2(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','nonzeroptr_j2k2')
    allocate(d1BOp(1,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOp')
    allocate(d1BOp_pi(1,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOp_pi')
    allocate(d1BOp_pipi(1,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOp_pipi')
    allocate(debpendBO_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','debpendBO_neigh')
    allocate(devaldBO_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','devaldBO_neigh')
    allocate(detorsdBOpi_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','detorsdBOpi_neigh')
    allocate(devalddeltaj_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','devalddeltaj_neigh')
    allocate(devaldsbo_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','devaldsbo_neigh')
    allocate(d1BOi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOi')
    allocate(d1BOi_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOi_s')
    allocate(d1BOi_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOi_pi')
    allocate(d1BOi_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOi_pipi')
    allocate(d1BOij2(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOij2')
    allocate(d1BOij2_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOij2_s')
    allocate(d1BOij2_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOij2_pi')
    allocate(d1BOij2_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOij2_pipi')
    allocate(d1BOik(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOik')
    allocate(d1BOik_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOik_s')
    allocate(d1BOik_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOik_pi')
    allocate(d1BOik_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOik_pipi')
    allocate(d1BOjk(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk')
    allocate(d1BOjk_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk_s')
    allocate(d1BOjk_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk_pi')
    allocate(d1BOjk_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk_pipi')
    allocate(d1BOjk2(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk2')
    allocate(d1BOjk2_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk2_s')
    allocate(d1BOjk2_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk2_pi')
    allocate(d1BOjk2_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d1BOjk2_pipi')
    allocate(dsboproddbo(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','dsboproddbo')
  endif
  if (lgrad2) then
    allocate(d2mix(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2max')
    allocate(d2BOi(maxnonzero2),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOi')
    allocate(d2BOi_s(maxnonzero2),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOi_s')
    allocate(d2BOi_pi(maxnonzero2),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOi_pi')
    allocate(d2BOi_pipi(maxnonzero2),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOi_pipi')
    allocate(d2BOjk(maxnonzero2),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOjk')
    allocate(d2BOjk_s(maxnonzero2),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOjk_s')
    allocate(d2BOjk_pi(maxnonzero2),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOjk_pi')
    allocate(d2BOjk_pipi(maxnonzero2),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOjk_pipi')
    allocate(d2BOp(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOp')
    allocate(d2BOp_pi(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOp_pi')
    allocate(d2BOp_pipi(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOp_pipi')
    allocate(d2ebpenddeltaidBO_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2ebpenddeltaidBO_neigh')
    allocate(d2ebpendBO2_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2ebpendBO2_neigh')
    allocate(d2etorsdBOpi2_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsdBOpi2_neigh')
    allocate(d2etorsdBOdBOpi_neigh(maxneigh1a,maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsdBOdBOpi_neigh')
    allocate(d2etorsdthetaddeltai(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsdthetaddeltai')
    allocate(d2evaldBO2_neigh(maxneigh2a),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldBO2_neigh')
    allocate(d2evalddeltaidBO_neigh(maxneigh1a),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evalddeltaidBO_neigh')
    allocate(d2evalddeltaiddeltaj_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evalddeltaiddeltaj_neigh')
    allocate(d2evalddeltajk_neigh(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evalddeltajk_neigh')
    allocate(d2evaldsbo_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldsbo_neigh')
    allocate(d2evaldsbo2_neigh(maxneigh1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldsbo2_neigh')
    allocate(d2etorsddeltadBOpi_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsddeltadBOpi_neigh')
    allocate(d2etorsddeltaidBOpi_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsddeltaidBOpi_neigh')
    allocate(d2etorsddeltajdBOpi_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsddeltajdBOpi_neigh')
    allocate(d2evalddeltajdBO_neigh(maxneigh,maxneigh1a),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evalddeltajdBO_neigh')
    allocate(d2evaldthetadBOpi_neigh(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldthetadBOpi_neigh')
    allocate(d2evaldthetaddeltai_neigh(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldthetaddeltai_neigh')
    allocate(d2sboproddbo2(maxneigh1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2sboproddbo2')
    allocate(d2sboddeltaidbo(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2sboddeltaidbo')
  else
    allocate(d2mix(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2mix')
    allocate(d2BOi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOi')
    allocate(d2BOi_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOi_s')
    allocate(d2BOi_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOi_pi')
    allocate(d2BOi_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOi_pipi')
    allocate(d2BOjk(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOjk')
    allocate(d2BOjk_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOjk_s')
    allocate(d2BOjk_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOjk_pi')
    allocate(d2BOjk_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOjk_pipi')
    allocate(d2BOp(1,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOp')
    allocate(d2BOp_pi(1,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOp_pi')
    allocate(d2BOp_pipi(1,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2BOp_pipi')
    allocate(d2ebpenddeltaidBO_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2ebpenddeltaidBO_neigh')
    allocate(d2ebpendBO2_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2ebpendBO2_neigh')
    allocate(d2etorsdBOpi2_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsdBOpi2_neigh')
    allocate(d2etorsdBOdBOpi_neigh(1,1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsdBOdBOpi_neigh')
    allocate(d2etorsdthetaddeltai(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsdthetaddeltai')
    allocate(d2evaldBO2_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldBO2_neigh')
    allocate(d2evalddeltaidBO_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evalddeltaidBO_neigh')
    allocate(d2evalddeltaiddeltaj_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evalddeltaiddeltaj_neigh')
    allocate(d2evalddeltajk_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evalddeltajk_neigh')
    allocate(d2evaldsbo_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldsbo_neigh')
    allocate(d2evaldsbo2_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldsbo2_neigh')
    allocate(d2etorsddeltadBOpi_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsddeltadBOpi_neigh')
    allocate(d2etorsddeltaidBOpi_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsddeltaidBOpi_neigh')
    allocate(d2etorsddeltajdBOpi_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2etorsddeltajdBOpi_neigh')
    allocate(d2evalddeltajdBO_neigh(1,1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evalddeltajdBO_neigh')
    allocate(d2evaldthetadBOpi_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldthetadBOpi_neigh')
    allocate(d2evaldthetaddeltai_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2evaldthetaddeltai_neigh')
    allocate(d2sboproddbo2(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2sboproddbo2')
    allocate(d2sboddeltaidbo(1),stat=status)
    if (status/=0) call outofmemory('reaxFFds','d2sboddeltaidbo')
  endif
!
!  Initial allocation of sparse arrays
!
  if (lgrad2) then
    call realloc(d1is,maxsparse,status)
    if (status.ne.0) call outofmemory('reaxFFds','d1is')
    maxsparse2 = maxsparse*(maxsparse+1)/2
    call realloc(d2is,maxsparse2,status)
    if (status.ne.0) call outofmemory('reaxFFds','d2is')
    call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
    if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
    call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
    if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
    call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
    if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
  elseif (lgrad1) then
    call realloc(d1is,maxsparse,status)
    if (status.ne.0) call outofmemory('reaxFFds','d1is')
  endif
!****************************
!  Find list of free atoms  *
!****************************
  if (lfreeze) then
    ii = 0
    do i = 1,numat
      if (lopf(nrelf2a(i))) then
        ii = ii + 1
        nfreeatom(i) = ii
      else
        nfreeatom(i) = 0
      endif
    enddo
  else
    do i = 1,numat
      nfreeatom(i) = i
    enddo
  endif
!*************************************
!  Find cut-off radii for all atoms  *
!*************************************
  nbosptr(1:numat) = 0
  nboatom = 0
  do i = 1,numat
    nati = nat(i)
    ntypi = nftype(i)
!
!  Set dummy pointer for derivative routines since neighbour list is not compressed yet for ReaxFF
!
    nbodummyptr(i) = i
!
!  Check repulsive terms and build atom pointers to species
!
    nbos(i) = 0
    do j = 1,nreaxFFspec
      if (nati.eq.natreaxFFspec(j).and.(ntypi.eq.ntypreaxFFspec(j).or.ntypreaxFFspec(j).eq.0)) then
        nbos(i) = nbos(i) + 1
        nbosptr(i) = j
      endif
    enddo
!
!  Check number of species for now
!
    if (nbos(i).gt.1) then
      call outerror('Multiple species per atom not yet allowed for in reaxFF',0_i4)
      call stopnow('reaxFFds')
    elseif (nbos(i).eq.1) then
      nboatom = nboatom + 1
      nboatomptr(i) = nboatom
      nboatomRptr(nboatom) = i
    endif
  enddo
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
  call getReaxFFneighbour(maxneigh,nbosptr,nneigh,neighno,rneigh, &
                          xneigh,yneigh,zneigh,ineigh,latomdone,lmaxneighok)
!**********************************************************************
!  If maxneigh has been exceeded, increase value and return to start  *
!**********************************************************************
  if (.not.lmaxneighok) then
    do i = 1,numat
      if (nneigh(i).gt.maxneigh) maxneigh = nneigh(i)
    enddo
    if (index(keyword,'verb').ne.0.and.ioproc) then
      write(ioout,'(/,''  Increasing maxneigh to '',i6)') maxneigh
    endif
    goto 100
  endif
!
!  Set pointer to atoms that are needed on this node
!
  numattodo = 0
  numattodoRptr(1:numat) = 0
!
  do i = 1,numat
!
!  Find if any atom, neighbour or neighbour of neighbour is needed on this node
!
!  First check on whether i is local
!
    lneeded = (atom2local(i).ne.0)
    numattodonptr(i) = atom2node(i)
!
!  Second, if i is not local check neighbours of i
!
    nj = 0
    do while (.not.lneeded.and.nj.lt.nneigh(i))
      nj = nj + 1
      j = neighno(nj,i)
      lneeded = (atom2local(j).ne.0)
    enddo
!
!  Third, if i and neighbours are not local check neighbours of j
!
    nj = 0
    do while (.not.lneeded.and.nj.lt.nneigh(i))
      nj = nj + 1
      j = neighno(nj,i)
      nk = 0
      do while (.not.lneeded.and.nk.lt.nneigh(j))
        nk = nk + 1
        k = neighno(nk,j)
        lneeded = (atom2local(k).ne.0)
      enddo
    enddo
!
!  Atom and shell is needed locally
!
    if (lneeded) then
      numattodo = numattodo + 1
      numattodoptr(numattodo) = i
      numattodoRptr(i) = numattodo
    endif
  enddo
!
!  Now include neighbours of atoms
!
  numattodo_main = numattodo
  do ni = 1,numattodo_main
    i = numattodoptr(ni)
    do nj = 1,nneigh(i)
      j = neighno(nj,i)
      if (numattodoRptr(j).eq.0) then
        numattodo = numattodo + 1
        numattodoptr(numattodo) = j
        numattodoRptr(j) = numattodo
      endif
    enddo
  enddo
!
!  Now include neighbours of neighbours of atoms - needed for torsions
!
  numattodo_main2 = numattodo
  do ni = numattodo_main+1,numattodo_main2
    i = numattodoptr(ni)
    do nj = 1,nneigh(i)
      j = neighno(nj,i)
      if (numattodoRptr(j).eq.0) then
        numattodo = numattodo + 1
        numattodoptr(numattodo) = j
        numattodoRptr(j) = numattodo
      endif
    enddo
  enddo
!*********************************************
!  Calculate numbers of neighbours for each  *
!*********************************************
  do i = 1,numat
!
!  Set initial value for lopanyneigh - this variable indicates whether an atom has 
!  any neighbours for which derivatives are required
!
    if (.not.lfreeze) then
      lopanyneigh(i) = .true.
    else
      lopanyneigh(i) = lopf(nrelf2a(i))
    endif
    do n = 1,nneigh(i)
      j = neighno(n,i)
      rij = rneigh(n,i)
!
!  Check whether atom is free to optimise
!
      if (lopf(nrelf2a(j))) then
        lopanyneigh(i) = .true.
      endif
    enddo
  enddo
  if (index(keyword,'debu').ne.0.and.ioproc) then
    write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
    write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    do i = 1,numat
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneigh(i)
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do i = 1,numat
      write(ioout,'(i4,8(1x,i4))') i,(neighno(nn,i),nn=1,nneigh(i))
    enddo
  endif
!   
!  Check size of sparse arrays
!   
  if (maxnonzero.gt.maxsparse) then
    if (lgrad2) then 
      maxsparse = maxnonzero
      call realloc(d1is,maxsparse,status)
      if (status.ne.0) call outofmemory('reaxFFds','d1is')
      maxsparse2 = maxsparse*(maxsparse+1)/2
      call realloc(d2is,maxsparse2,status)
      if (status.ne.0) call outofmemory('reaxFFds','d2is')
      call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
      if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
      call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
      call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
    elseif (lgrad1) then
      maxsparse = maxnonzero
      call realloc(d1is,maxsparse,status)
      if (status.ne.0) call outofmemory('reaxFFds','d1is')
    endif
  endif
!*************************************************************************
!  Loop over pairs of atoms to compute bond order prime and delta prime  *
!*************************************************************************
  allocate(nbotodoptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFds','nbotodoptr')
!
  if (lextended) then
    nbotodo = numat
    do i = 1,numat
      nbotodoptr(i) = i
    enddo
  else
    nbotodo = numattodo
    do i = 1,numattodo
      nbotodoptr(i) = numattodoptr(i)
    enddo
  endif
!
!  mneigh contains the real value of maxneigh need after removing neighbours
!  whose bond order tolerance is below the allowed threshold
!
  mneigh = 0
  do ii = 1,nbotodo
    i = nbotodoptr(ii)
!
!  Set variables relating to i
!
    if (nbos(i).gt.0) then
      nspeci = nbosptr(i)
!
!  Initialise deltap
!
      deltap(i) = - reaxFFval(1,nspeci)
!
!  Loop over neighbours of i 
!
      nivalid = 0
      do ni = 1,nneigh(i)
        j = neighno(ni,i)
!
!  Does j have a bond order species type?
!
        if (nbos(j).gt.0) then
!
!  Set variables relating to j
!
          nspecj = nbosptr(j)
!
!  Set up i-j quantities
!
          rij = rneigh(ni,i)
!
!  Set index to parameters for pairwise interaction of species
!
          if (nspeci.ge.nspecj) then
            nboij = nspeci*(nspeci - 1)/2 + nspecj
          else
            nboij = nspecj*(nspecj - 1)/2 + nspeci
          endif
          lanyBOnonzero = .false.
          if (rij.lt.reaxFFrmaxpair(nboij)) then
!**********************************
!  Valid bond order contribution  *
!**********************************
!
!  Sigma
!
            call reaxFFbo_sigma(nboij,nspeci,nspecj,rij,BOpij,dBOpijdr,d2BOpijdr2,lgrad1,lgrad2, &
                                .not.lStrict,lBOnonzero)
            if (lBOnonzero) then
              lanyBOnonzero = .true.
              BOp(ni,i) = BOpij
              deltap(i) = deltap(i) + BOpij
              if (lgrad1) then
                rrij = 1.0_dp/rij
                d1BOp(ni,i) = rrij*dBOpijdr
                if (lgrad2) then
                  d2BOp(ni,i) = rrij*rrij*(d2BOpijdr2 - rrij*dBOpijdr)
                endif
              endif
!
!  Pi
!
              call reaxFFbo_pi(nboij,nspeci,nspecj,rij,BOpij,dBOpijdr,d2BOpijdr2,lgrad1,lgrad2, &
                               .not.lStrict)
              BOp_pi(ni,i) = BOpij
              deltap(i) = deltap(i) + BOpij
              if (lgrad1) then
                d1BOp_pi(ni,i) = rrij*dBOpijdr
                if (lgrad2) then
                  d2BOp_pi(ni,i) = rrij*rrij*(d2BOpijdr2 - rrij*dBOpijdr)
                endif
              endif
!
!  Pi_pi
!
              call reaxFFbo_pipi(nboij,nspeci,nspecj,rij,BOpij,dBOpijdr,d2BOpijdr2,lgrad1,lgrad2, &
                                 .not.lStrict)
              BOp_pipi(ni,i) = BOpij
              deltap(i) = deltap(i) + BOpij
              if (lgrad1) then
                d1BOp_pipi(ni,i) = rrij*dBOpijdr
                if (lgrad2) then
                  d2BOp_pipi(ni,i) = rrij*rrij*(d2BOpijdr2 - rrij*dBOpijdr)
                endif
              endif
            else
              BOp(ni,i) = 0.0_dp
              BOp_pi(ni,i) = 0.0_dp
              BOp_pipi(ni,i) = 0.0_dp
              if (lgrad1) then
                d1BOp(ni,i) = 0.0_dp
                d1BOp_pi(ni,i) = 0.0_dp
                d1BOp_pipi(ni,i) = 0.0_dp
                if (lgrad2) then
                  d2BOp(ni,i) = 0.0_dp
                  d2BOp_pi(ni,i) = 0.0_dp
                  d2BOp_pipi(ni,i) = 0.0_dp
                endif
              endif
            endif
          else
            BOp(ni,i) = 0.0_dp
            BOp_pi(ni,i) = 0.0_dp
            BOp_pipi(ni,i) = 0.0_dp
            if (lgrad1) then
              d1BOp(ni,i) = 0.0_dp
              d1BOp_pi(ni,i) = 0.0_dp
              d1BOp_pipi(ni,i) = 0.0_dp
              if (lgrad2) then
                d2BOp(ni,i) = 0.0_dp
                d2BOp_pi(ni,i) = 0.0_dp
                d2BOp_pipi(ni,i) = 0.0_dp
              endif
            endif
          endif
          if (lanyBOnonzero) then
!
!  Valid terms found - move terms to correct location
!
            nivalid = nivalid + 1
            if (nivalid.ne.ni) then
              neighno(nivalid,i) = neighno(ni,i)
              rneigh(nivalid,i) = rneigh(ni,i)
              xneigh(nivalid,i) = xneigh(ni,i)
              yneigh(nivalid,i) = yneigh(ni,i)
              zneigh(nivalid,i) = zneigh(ni,i)
              BOp(nivalid,i) = BOp(ni,i)
              BOp_pi(nivalid,i) = BOp_pi(ni,i)
              BOp_pipi(nivalid,i) = BOp_pipi(ni,i)
              if (lgrad1) then
                d1BOp(nivalid,i) = d1BOp(ni,i)
                d1BOp_pi(nivalid,i) = d1BOp_pi(ni,i)
                d1BOp_pipi(nivalid,i) = d1BOp_pipi(ni,i)
                if (lgrad2) then
                  d2BOp(nivalid,i) = d2BOp(ni,i)
                  d2BOp_pi(nivalid,i) = d2BOp_pi(ni,i)
                  d2BOp_pipi(nivalid,i) = d2BOp_pipi(ni,i)
                endif
              endif
            endif
          endif
        endif
      enddo
!  
!  Now reset number of neighbours to reflect the number that have non-zero bond orders
!           
      nneigh(i) = nivalid
      mneigh = max(mneigh,nivalid)
    endif
!
!  End loop over atoms i
!
  enddo
!
  deallocate(nbotodoptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nbotodoptr')
!
!  Set neighnoRptr
!
  do i = 1,numat
    if (nbos(i).gt.0) then
!
!  Loop over neighbours of i 
!
      nivalid = 0
      do ni = 1,nneigh(i)
        j = neighno(ni,i)
!
!  Does j have a bond order species type?
!
        if (nbos(j).gt.0) then
!
!  Set up i-j quantities
!
          rij = rneigh(ni,i)
          xji = xneigh(ni,i)
          yji = yneigh(ni,i)
          zji = zneigh(ni,i)
!
!  Find i in neighbour list for j
!
          nj = 1
          lfound = .false.
          do while (nj.le.nneigh(j).and..not.lfound)
            if (neighno(nj,j).eq.i) then
              xdiff = xneigh(nj,j) + xji
              ydiff = yneigh(nj,j) + yji
              zdiff = zneigh(nj,j) + zji
              lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
            endif
            if (.not.lfound) nj = nj + 1
          enddo
          if (lfound) then
            neighnoRptr(ni,i) = nj
          else
            call outerror('neighbour lists are inconsistent in reaxFF',0_i4)
            call stopnow('reaxff')
          endif
        endif
      enddo
    endif
!
!  End loop over atoms i
!
  enddo
  if (index(keyword,'debu').ne.0) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  Number of neighbours for atoms after bond order screening:'',/)')
      write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    endif
#ifdef MPI
    if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
      ntmp = 1
      ntag = 1
      allocate(StatMPI(MPI_Status_Size,1),stat=status)
      if (status/=0) call outofmemory('reaxFFds','StatMPI')
!
      do i = 1,numat
        if (procid.eq.numattodonptr(i).and.ioproc) then
!
!  Data is local to I/O node
!
          write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneigh(i)
        elseif (procid.eq.numattodonptr(i).or.ioproc) then
!
!  Post receive
!
          if (ioproc) then
            nnode = numattodonptr(i)
            call MPI_IRecv(nloc,ntmp,MPI_integer,nnode, &
                           ntag,MPI_Comm_GULP,Request,MPIerror)
          else
!
!  Pass data to ioproc for writing
!
            nloc = nneigh(i)
!
!  Post send
!
            call MPI_ISend(nloc,ntmp,MPI_integer,0, &
                           ntag,MPI_Comm_GULP,Request,MPIerror)
          endif
          call MPI_WaitAll(1,Request,StatMPI,MPIerror)
          if (ioproc) then
!
!  Write on I/O node
!
            write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nloc
          endif
        endif
      enddo
!
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('reaxFFds','StatMPI')
    else
#endif
!
      call mpbarrier
      do i = 1,numat
        if (procid.eq.numattodonptr(i)) then
          write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneigh(i)
        endif
        call mpbarrier
      enddo
#ifdef MPI
    endif
#endif
  endif
  if (index(keyword,'verb').ne.0) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Delta and Bond Order prime: '',/)')
    endif
    call mpbarrier
#ifdef MPI
    if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
      ntmp = maxneigh + 1
      ntag = 1
      ntag2 = 2
      allocate(nntmp(ntmp),stat=status)
      if (status/=0) call outofmemory('reaxFFds','nntmp')
      allocate(botmp(2*ntmp),stat=status)
      if (status/=0) call outofmemory('reaxFFds','botmp')
      allocate(StatMPI(MPI_Status_Size,2),stat=status)
      if (status/=0) call outofmemory('reaxFFds','StatMPI')
!
      do i = 1,numat
        if (procid.eq.numattodonptr(i).and.ioproc) then
!
!  Data is local to I/O node
!
          write(ioout,'(i4,f10.6,6(1x,i4,1x,f7.4))') i,deltap(i),(neighno(j,i),(BOp(j,i)+BOp_pi(j,i)+BOp_pipi(j,i)),j=1,nneigh(i))
        elseif (procid.eq.numattodonptr(i).or.ioproc) then
!
!  Post receives
!
          if (ioproc) then
            nnode = numattodonptr(i)
            call MPI_IRecv(nntmp,ntmp,MPI_integer,nnode, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_IRecv(botmp,ntmp,MPI_double_precision,nnode, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          else
!
!  Pass data to ioproc for writing
!
            nntmp(1:ntmp) = 0
            nntmp(1) = nneigh(i)
            do j = 1,nneigh(i)
              nntmp(j+1) = neighno(j,i)
            enddo
            botmp(1:ntmp) = 0
            botmp(1) = deltap(i)
            do j = 1,nneigh(i)
              botmp(1+j) = BOp(j,i) + BOp_pi(j,i) + BOp_pipi(j,i)
            enddo
!
!  Post sends
!
            call MPI_ISend(nntmp,ntmp,MPI_integer,0, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_ISend(botmp,ntmp,MPI_double_precision,0, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          endif
          call MPI_WaitAll(2,Request,StatMPI,MPIerror)
          if (ioproc) then
!
!  Write on I/O node
!
            numtmp = nntmp(1)
            write(ioout,'(i4,f10.6,6(1x,i4,1x,f7.4))') i,botmp(1),(nntmp(1+j),botmp(1+j),j=1,numtmp)
          endif
        endif
      enddo
    else
#endif
      do i = 1,numat
        if (procid.eq.numattodonptr(i)) then
          write(ioout,'(i4,f10.6,6(1x,i4,1x,f7.4))') i,deltap(i),(neighno(j,i),(BOp(j,i)+BOp_pi(j,i)+BOp_pipi(j,i)),j=1,nneigh(i))
        endif
        call mpbarrier
      enddo
#ifdef MPI
    endif
#endif
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Bond Order sigma prime: '',/)')
    endif
    call mpbarrier
#ifdef MPI
    if (lioproconly) then
      do i = 1,numat
        if (procid.eq.numattodonptr(i).and.ioproc) then
!
!  Data is local to I/O node
!
          write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(neighno(j,i),(BOp(j,i)),j=1,nneigh(i))
        elseif (procid.eq.numattodonptr(i).or.ioproc) then
!
!  Post receives
!
          if (ioproc) then
            nnode = numattodonptr(i)
            call MPI_IRecv(nntmp,ntmp,MPI_integer,nnode, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_IRecv(botmp,ntmp,MPI_double_precision,nnode, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          else
!
!  Pass data to ioproc for writing
!
            nntmp(1:ntmp) = 0
            nntmp(1) = nneigh(i)
            do j = 1,nneigh(i)
              nntmp(j+1) = neighno(j,i)
            enddo
            botmp(1:ntmp) = 0
            do j = 1,nneigh(i)
              botmp(1+j) = BOp(j,i)
            enddo
!
!  Post sends
!
            call MPI_ISend(nntmp,ntmp,MPI_integer,0, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_ISend(botmp,ntmp,MPI_double_precision,0, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          endif
          call MPI_WaitAll(2,Request,StatMPI,MPIerror)
          if (ioproc) then
!
!  Write on I/O node
!
            numtmp = nntmp(1)
            write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(nntmp(1+j),botmp(1+j),j=1,numtmp)
          endif
        endif
      enddo
    else
#endif
      do i = 1,numat
        if (procid.eq.numattodonptr(i)) then
          write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(neighno(j,i),(BOp(j,i)),j=1,nneigh(i))
        endif
        call mpbarrier
      enddo
#ifdef MPI
    endif
#endif
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Bond Order pi prime: '',/)')
    endif
    call mpbarrier
#ifdef MPI
    if (lioproconly) then
      do i = 1,numat
        if (procid.eq.numattodonptr(i).and.ioproc) then
!
!  Data is local to I/O node
!
          write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(neighno(j,i),(BOp_pi(j,i)),j=1,nneigh(i))
        elseif (procid.eq.numattodonptr(i).or.ioproc) then
!
!  Post receives
!
          if (ioproc) then
            nnode = numattodonptr(i)
            call MPI_IRecv(nntmp,ntmp,MPI_integer,nnode, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_IRecv(botmp,ntmp,MPI_double_precision,nnode, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          else
!
!  Pass data to ioproc for writing
!
            nntmp(1:ntmp) = 0
            nntmp(1) = nneigh(i)
            do j = 1,nneigh(i)
              nntmp(j+1) = neighno(j,i)
            enddo
            botmp(1:ntmp) = 0
            do j = 1,nneigh(i)
              botmp(1+j) = BOp_pi(j,i)
            enddo
!
!  Post sends
!
            call MPI_ISend(nntmp,ntmp,MPI_integer,0, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_ISend(botmp,ntmp,MPI_double_precision,0, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          endif
          call MPI_WaitAll(2,Request,StatMPI,MPIerror)
          if (ioproc) then
!
!  Write on I/O node
!
            numtmp = nntmp(1)
            write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(nntmp(1+j),botmp(1+j),j=1,numtmp)
          endif
        endif
      enddo
    else
#endif
      do i = 1,numat
        if (procid.eq.numattodonptr(i)) then
          write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(neighno(j,i),(BOp_pi(j,i)),j=1,nneigh(i))
        endif
        call mpbarrier
      enddo
#ifdef MPI
    endif
#endif
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Bond Order pi-pi prime: '',/)')
    endif
    call mpbarrier
#ifdef MPI
    if (lioproconly) then
      do i = 1,numat
        if (procid.eq.numattodonptr(i).and.ioproc) then
!
!  Data is local to I/O node
!
          write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(neighno(j,i),(BOp_pipi(j,i)),j=1,nneigh(i))
        elseif (procid.eq.numattodonptr(i).or.ioproc) then
!
!  Post receives
!
          if (ioproc) then
            nnode = numattodonptr(i)
            call MPI_IRecv(nntmp,ntmp,MPI_integer,nnode, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_IRecv(botmp,ntmp,MPI_double_precision,nnode, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          else
!
!  Pass data to ioproc for writing
!
            nntmp(1:ntmp) = 0
            nntmp(1) = nneigh(i)
            do j = 1,nneigh(i)
              nntmp(j+1) = neighno(j,i)
            enddo
            botmp(1:ntmp) = 0
            do j = 1,nneigh(i)
              botmp(1+j) = BOp_pipi(j,i)
            enddo
!
!  Post sends
!
            call MPI_ISend(nntmp,ntmp,MPI_integer,0, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_ISend(botmp,ntmp,MPI_double_precision,0, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          endif
          call MPI_WaitAll(2,Request,StatMPI,MPIerror)
          if (ioproc) then
!
!  Write on I/O node
!
            numtmp = nntmp(1)
            write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(nntmp(1+j),botmp(1+j),j=1,numtmp)
          endif
        endif
      enddo
!
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('reaxFFds','StatMPI')
      deallocate(botmp,stat=status)
      if (status/=0) call deallocate_error('reaxFFds','botmp')
      deallocate(nntmp,stat=status)
      if (status/=0) call deallocate_error('reaxFFds','nntmp')
    else
#endif
      do i = 1,numat
        if (procid.eq.numattodonptr(i)) then
          write(ioout,'(i4,6(1x,i4,1x,f7.4))') i,(neighno(j,i),(BOp_pipi(j,i)),j=1,nneigh(i))
        endif
        call mpbarrier
      enddo
#ifdef MPI
    endif
#endif
    if (ioproc) then
      write(ioout,'(/)')
    endif
  endif
!*****************************
!  Loop over pairs of atoms  *
!*****************************
  iloop: do ii = 1,nboatom
    i = nboatomptr(ii)
!       
!  If none of the neighhbours are local to this node then cycle
!     
    if (numattodoRptr(i).eq.0) cycle
!     
    lilocal = (atom2local(i).gt.0)
!
!  Set variables relating to i
!
    nspeci = nbosptr(i)
    nregioni = nregionno(nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    lslicei = lsliceatom(nsft + nrelf2a(i))
!
!  Loop over neighbours of i (=> j)
!
    do ni = 1,nneigh(i)
      j = neighno(ni,i)
!
      nregionj = nregionno(nrelf2a(j))
      nregiontypj = nregiontype(nregionj,ncf)
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
      lQMMMok = .true.
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
      endif
!
!  Do we need to do this pair of atoms
!
      if (nbos(j).gt.0) then
!   
!  Initialise sparse derivative storage for neighbours of i
!   
        if (lgrad1) then
          d1is(1:maxnonzero) = 0.0_dp
          if (lgrad2) then
            d2is(1:maxnonzero2) = 0.0_dp
          endif
        endif
!
!  If i = j set scale to be half to correct for double counting
!
        scale = 0.5_dp
!
!  Set variables relating to j
!
        nspecj = nbosptr(j)
        lslicej = lsliceatom(nsft + nrelf2a(j))
!
!  Set up i-j quantities
!
        lreg2one  = .false.
        lreg2pair = .false.
        if (lseok.and.nregions(ncf).gt.1) then
          lreg2pair = (nregioni.eq.nregionj.and.lregionrigid(nregioni,ncf))
          if (.not.lreg2pair) lreg2one = (lregionrigid(nregioni,ncf).neqv.lregionrigid(nregionj,ncf))
        endif
        lattach = (lslicei.and..not.lslicej.or.lslicej.and..not.lslicei)
!
!  Find i in neighbour list for j
!
        nj = neighnoRptr(ni,i)
!
!  Set index to parameters for pairwise interaction of species
!
        if (nspeci.ge.nspecj) then
          nboij = nspeci*(nspeci - 1)/2 + nspecj
        else
          nboij = nspecj*(nspecj - 1)/2 + nspeci
        endif
!
!  Set delta' 
!
        deltapi = deltap(i)
        deltapj = deltap(j)
!***************************************
!  Valid pairwise reaxFF contribution  *
!***************************************
        call reaxFF_bos(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                        d1BOp,d1BOp_pi,d1BOp_pipi,d2BOp,d2BOp_pi,d2BOp_pipi,BOij_s,BOij_pi,BOij_pipi, &
                        d1BOi_s,d1BOi_pi,d1BOi_pipi,d2BOi_s,d2BOi_pi,d2BOi_pipi,nonzero_ij, &
                        nonzero_ij_i,nonzeroptr_ij,lgrad1,lgrad2)
!
!  Convert non-zero indices
!
        if (lgrad1) then
          ind0 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + ni*mneigh
          do nl = nonzero_ij_i+1,nonzero_ij
            nonzeroptr_ij(nl) = ind0 + nonzeroptr_ij(nl)
          enddo
        endif
!
!  Save corrected bond orders
!
        BO(ni,i) = BOij_s + BOij_pi + BOij_pipi
        BO(nj,j) = BOij_s + BOij_pi + BOij_pipi
        BO_pi(ni,i) = BOij_pi
        BO_pi(nj,j) = BOij_pi
        BO_pipi(ni,i) = BOij_pipi
        BO_pipi(nj,j) = BOij_pipi
!
!  If lQMMMok then compute energy terms
!
        if (lreaxff_ebond.and.lQMMMok) then
!
!  Raise BOij to the power of Pbe,2 (in paper says Pbe,1) but I think this is a typo!
!
          BOijpbe2 = BOij_s**(reaxFFpbe(2,nboij))
!
!  Compute exponential factor
!
          expij = exp(reaxFFpbe(1,nboij)*(1.0_dp - BOijpbe2))
!
!  Calculate total i-j potential
!
          eij = - scale*(reaxFFDe(1,nboij)*BOij_s*expij + reaxFFDe(2,nboij)*BOij_pi + reaxFFDe(3,nboij)*BOij_pipi)
!
!  Add to surface energy totals if appropriate
!
          if (lilocal) then
            if (lseok) then
              if (lreg2one) then
                esregion12 = esregion12 + eij
              elseif (lreg2pair) then
                esregion2 = esregion2 + eij
              else
                ebond = ebond + eij
              endif
            else
              ebond = ebond + eij
            endif
            if (lattach) eattach = eattach + eij
!
!  Site energy contributions for ebond
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*eij
            siteenergy(j) = siteenergy(j) + 0.5_dp*eij
          endif
!
!  Derivatives of Bond Order potential energy
!
          if (lgrad1) then
            deijdBO_s    = - scale*reaxFFDe(1,nboij)*expij*(1.0_dp - reaxFFpbe(1,nboij)*reaxFFpbe(2,nboij)*BOijpbe2)
            deijdBO_pi   = - scale*reaxFFDe(2,nboij)
            deijdBO_pipi = - scale*reaxFFDe(3,nboij)
            do nk = 1,nonzero_ij
              d1is(nk) = d1is(nk) + deijdBO_s*d1BOi_s(nk) + deijdBO_pi*d1BOi_pi(nk) + deijdBO_pipi*d1BOi_pipi(nk)
            enddo
!
            if (lgrad2) then
              d2eijdBO2_s = scale*reaxFFDe(1,nboij)*expij*reaxFFpbe(1,nboij)*reaxFFpbe(2,nboij)* &
                (BOij_s**(reaxFFpbe(2,nboij)-1.0_dp))*((1.0_dp - reaxFFpbe(1,nboij)*reaxFFpbe(2,nboij)*BOijpbe2) + &
                reaxFFpbe(2,nboij))
              inds = 0
              do nk = 1,nonzero_ij
                do nl = 1,nk
                  inds = inds + 1
                  d2is(inds) = d2is(inds) + deijdBO_s*d2BOi_s(inds) + deijdBO_pi*d2BOi_pi(inds) + deijdBO_pipi*d2BOi_pipi(inds)
                  d2is(inds) = d2is(inds) + d2eijdBO2_s*d1BOi_s(nk)*d1BOi_s(nl)
                enddo 
              enddo
!
!  Add second derivatives
!
              call d2adddmsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nboatomRptr,d1is,d2is, &
                             nonzero_ij,nonzeroptr_ij,.true.,.false.)
            endif
!
!  Add first derivatives
!
            if (lilocal) then
              call d1addsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nbodummyptr,d1is, &
                           nonzero_ij,nonzeroptr_ij,.true.,.false.,.true.)
            endif
          endif
        endif
!
!  End condition section on i or j being associated with moving atom
!
      endif
    enddo
  enddo iloop
!***************************************************
!  Loop over neighbours of atoms to compute delta  *
!***************************************************
  delta(1:numat) = 0.0_dp
  deltalp(1:numat) = 0.0_dp
  deltalplp(1:numat) = 0.0_dp
  ddeltaeddeltai = 0.5_dp
!
  do ii = 1,numattodo
    i = numattodoptr(ii)
!
!  Set variables relating to i
!
    if (nbos(i).gt.0) then
      nspeci = nbosptr(i)
!
!  Compute delta from refined bond orders
!
      delta(i) = - reaxFFval(1,nspeci)
      do ni = 1,nneigh(i)
        delta(i) = delta(i) + BO(ni,i)
      enddo
!
!  Compute 1/2 x delta_e
!
      delta_e = 0.5_dp*(delta(i) + reaxFFval(1,nspeci) - reaxFFval(3,nspeci))
!
!  Compute deltalp for lone pairs
!
      if (lreaxFFlpsmooth) then
!
!  New smooth form
!
        delta_e0 = nint(delta_e)
        smoothH = 1.0_dp/(1.0_dp + exp(-reaxFFlpsmooth*(delta_e - delta_e0)))
        rintde = delta_e0 - 1.0_dp + smoothH
      else
!
!  Old discontinuous form
!
        rintde = dble(int(delta_e))
      endif
      rnlp = - rintde + exp(-4.0_dp*reaxFFlam(29)*(1.0_dp + delta_e - rintde)**2)
      deltalplp(i) = reaxFFlp(1,nspeci) - rnlp
!
!  Compute deltalp for lone pairs for overcoordination term
!
      if (lreaxFFovsmooth) then
!
!  New smooth form
!
        delta_e0 = nint(delta_e)
        smoothH = 1.0_dp/(1.0_dp + exp(-reaxFFlpsmooth*(delta_e - delta_e0)))
        rintde = delta_e0 - 1.0_dp + smoothH
      else
!
!  Old discontinuous form
!
        rintde = dble(int(delta_e))
      endif
      rnlp = - rintde + exp(-4.0_dp*reaxFFlam(29)*(1.0_dp + delta_e - rintde)**2)
      deltalp(i) = reaxFFlp(1,nspeci) - rnlp
    endif
  enddo
!
  if (index(keyword,'verb').ne.0) then
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Delta and Bond Order: '',/)')
    endif
#ifdef MPI
    if (lioproconly) then
!
!  Allocate temporary workspace for communication
!
      ntmp = maxneigh + 1
      ntag = 1
      ntag2 = 2
      allocate(nntmp(ntmp),stat=status)
      if (status/=0) call outofmemory('reaxFFds','nntmp')
      allocate(botmp(2*ntmp),stat=status)
      if (status/=0) call outofmemory('reaxFFds','botmp')
      allocate(StatMPI(MPI_Status_Size,2),stat=status)
      if (status/=0) call outofmemory('reaxFFds','StatMPI')
!
      do i = 1,numat
        if (procid.eq.numattodonptr(i).and.ioproc) then
!
!  Data is local to I/O node
!
          write(ioout,'(i4,f10.6,6(1x,i4,1x,f7.4))') i,delta(i),(neighno(j,i),BO(j,i),j=1,nneigh(i))
        elseif (procid.eq.numattodonptr(i).or.ioproc) then
!
!  Post receives
!
          if (ioproc) then
            nnode = numattodonptr(i)
            call MPI_IRecv(nntmp,ntmp,MPI_integer,nnode, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_IRecv(botmp,ntmp,MPI_double_precision,nnode, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          else
!
!  Pass data to ioproc for writing
!
            nntmp(1:ntmp) = 0
            nntmp(1) = nneigh(i)
            do j = 1,nneigh(i)
              nntmp(j+1) = neighno(j,i)
            enddo
            botmp(1:ntmp) = 0.0_dp
            botmp(1) = delta(i)
            do j = 1,nneigh(i)
              botmp(1+j) = BO(j,i)
            enddo
!
!  Post sends
!
            call MPI_ISend(nntmp,ntmp,MPI_integer,0, &
                           ntag,MPI_Comm_GULP,Request(1),MPIerror)
            call MPI_ISend(botmp,ntmp,MPI_double_precision,0, &
                           ntag2,MPI_Comm_GULP,Request(2),MPIerror)
          endif
          call MPI_WaitAll(2,Request,StatMPI,MPIerror)
          if (ioproc) then
!
!  Write on I/O node
!
            numtmp = nntmp(1)
            write(ioout,'(i4,f10.6,6(1x,i4,1x,f7.4))') i,botmp(1),(nntmp(1+j),botmp(1+j),j=1,numtmp)
          endif
        endif
      enddo
!
      deallocate(StatMPI,stat=status)
      if (status/=0) call deallocate_error('reaxFFds','StatMPI')
      deallocate(botmp,stat=status)
      if (status/=0) call deallocate_error('reaxFFds','botmp')
      deallocate(nntmp,stat=status)
      if (status/=0) call deallocate_error('reaxFFds','nntmp')
    else
#endif
      call mpbarrier
      do i = 1,numat
        if (procid.eq.numattodonptr(i)) then
          write(ioout,'(i4,f10.6,6(1x,i4,1x,f7.4))') i,delta(i),(neighno(j,i),BO(j,i),j=1,nneigh(i))
        endif
        call mpbarrier
      enddo
#ifdef MPI
    endif
#endif
    call gflush(ioout)
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Delta lone pair: '',/)')
      do i = 1,numat
        write(ioout,'(i4,f10.6)') i,deltalplp(i)
      enddo
      write(ioout,'(/)')
      call gflush(ioout)
    endif
    call mpbarrier
    if (ioproc) then
      write(ioout,'(/,''  ReaxFF: Delta lone pair for overcoordination: '',/)')
      do i = 1,numat
        write(ioout,'(i4,f10.6)') i,deltalp(i)
      enddo
      write(ioout,'(/)')
      call gflush(ioout)
    endif
  endif
!***********************************************************
!  Loop over atoms to compute over and under coordination  *
!***********************************************************
  do ii = 1,nboatom
    i = nboatomptr(ii)
!
!  If none of the neighhbours are local to this node then cycle
!
    if (numattodoRptr(i).eq.0) cycle
!
!  Set variables relating to i
!
    lilocal = (atom2local(i).gt.0)
!
!  Set variables relating to i
!
    if (nbos(i).gt.0) then
      nspeci = nbosptr(i)
      nregioni = nregionno(nrelf2a(i))
      nregiontypi = nregiontype(nregioni,ncf)
      lslicei = lsliceatom(nsft + nrelf2a(i))
!   
!  QM/MM handling : i is a QM atom => exclude
!   
      lQMMMok = .true.                           
      if (QMMMmode(ncf).gt.0) then               
        if (nregiontypi.eq.1) lQMMMok = .false.  
      endif
!
!  Do we need to do this atom based on QMMM scheme?
!
      if (lQMMMok) then
!
!  Setup in case we need derivatives
!
        if (lgrad1) then
!
!  Set total number of distances for neighbours of i
!
          nneighi1   = nneigh(i)*(nneigh(i) + 1)/2 
          nneighi1a  = maxneigh*(nneigh(i) + 1)
          nneighi2   = nneigh(i) + nneighi1
          nneighi2a  = nneighi1a*(nneighi1a + 1)/2
          nneighi22  = nneighi2*(nneighi2 + 1)/2
          nneighi2t  = nneigh(i) + nneighi1 + nneigh(i)*(mneigh+1)*mneigh
          if (lextended) then
            nneighi2e = nneigh(i) + nneighi1 + nneigh(i)*(2*mneigh+1)*mneigh
          else
            nneighi2e = nneighi2t
          endif
!
!  Initialise derivative storage for neighbours of i
!
          d1is(1:maxsparse) = 0.0_dp
          if (lgrad2) then
            d2is(1:maxsparse2) = 0.0_dp
          endif
!
!  Initialise sparse indices
!
          nonzero = 0
          nonzerorptr(1:nneighi2e) = 0
        endif
!
!  Lone pair energy
!
        if (lreaxFF_elone) then
          if (abs(reaxFFlp(2,nspeci)).gt.1.0d-12) then
            explp = exp(-reaxFFlam(30)*deltalplp(i))
            trm1lp = 1.0_dp/(1.0_dp + explp)
            if (lilocal) then
              esite = reaxFFlp(2,nspeci)*deltalplp(i)*trm1lp
              elone = elone + esite
!
!  Site energy contributions for elone
!
              siteenergy(i) = siteenergy(i) + esite
            endif
          endif
        endif
!
!  Bond penalty energy
!
        if (lreaxFF_ebpen.and.lilocal) then
          do ni = 1,nneigh(i)
            j = neighno(ni,i)
            if (nbos(j).gt.0) then
              nspecj = nbosptr(j)
              if (nspeci.ge.nspecj) then
                nboij = nspeci*(nspeci - 1)/2 + nspecj
              else
                nboij = nspecj*(nspecj - 1)/2 + nspeci
              endif
              diffBOdelta = BO(ni,i) - delta(i) - reaxFFpen2(2,nboij)*delta(i)**4 - reaxFFpen2(3,nboij)
              if (diffBOdelta.gt.0.0_dp) then
                ebpenij = reaxFFpen2(1,nboij)*diffBOdelta**2
                ebpen = ebpen + ebpenij
!
!  Site energy contributions for ebpen
!
                siteenergy(i) = siteenergy(i) + 0.5_dp*ebpenij
                siteenergy(j) = siteenergy(j) + 0.5_dp*ebpenij
              endif
            endif
          enddo
        endif
!
!  Compute sum over neighbours for over/under-coordination energy terms
!
        sumdeltabopi = 0.0_dp
        sumover = 0.0_dp
        do ni = 1,nneigh(i)
          j = neighno(ni,i)
          if (nbos(j).gt.0) then
            nspecj = nbosptr(j)
            if (nspeci.ge.nspecj) then
              nboij = nspeci*(nspeci - 1)/2 + nspecj
            else
              nboij = nspecj*(nspecj - 1)/2 + nspeci
            endif
            if (.not.lreaxFFunder(nspeci).or..not.lreaxFFunder(nspecj)) then
              sumdeltabopi = sumdeltabopi + delta(j)*(BO_pi(ni,i) + BO_pipi(ni,i))
            else
              sumdeltabopi = sumdeltabopi + (delta(j) - deltalp(j))*(BO_pi(ni,i) + BO_pipi(ni,i))
            endif
            sumover = sumover + reaxFFoc2(nboij)*reaxFFDe(1,nboij)*BO(ni,i)
          endif
        enddo
!
!  Compute delta_lpcorr_i - expression is different for first row elements
!
        if (.not.lreaxFFunder(nspeci)) then
          expdi = 0.0_dp
          otrm1 = 0.0_dp
          deltalpcorr_i = delta(i)
        else
          expdi = exp(reaxFFlam(31)*sumdeltabopi)
          otrm1 = 1.0_dp/(1.0_dp + reaxFFlam(6)*expdi)
          deltalpcorr_i = delta(i) - deltalp(i)*otrm1
        endif
!
!  Over coordination term
!
        if (lreaxFF_eover) then
          otrm2 = 1.0_dp/(deltalpcorr_i + reaxFFval(1,nspeci) + 1.0d-8)
          expoc = exp(reaxFFoc1(nspeci)*deltalpcorr_i)
          otrm3 = 1.0_dp/(1.0_dp + expoc)
          if (lilocal) then
            esite = sumover*otrm2*deltalpcorr_i*otrm3
            eover = eover + esite
!
!  Site energy contributions for eover
!
            siteenergy(i) = siteenergy(i) + esite
          endif
        else
          expoc = 0.0_dp
          otrm2 = 0.0_dp
          otrm3 = 0.0_dp
        endif
!
!  Under coordination term
!
        if (lreaxFF_eunder) then
          expuc1 = exp(reaxFFlam(7)*deltalpcorr_i)
          expuc2 = exp(-reaxFFoc1(nspeci)*deltalpcorr_i)
          loverflow_expuc3 = (reaxFFlam(9)*sumdeltabopi.gt.80.0_dp)
          if (loverflow_expuc3) then
            expuc3 = 0.0_dp
            utrm2  = 0.0_dp
          else
            expuc3 = exp(reaxFFlam(9)*sumdeltabopi)
            utrm2  = 1.0_dp/(1.0_dp + reaxFFlam(8)*expuc3)
          endif
          utrm1  = 1.0_dp/(1.0_dp + expuc2)
          if (lilocal) then
            esite  = - reaxFFuc1(nspeci)*(1.0_dp - expuc1)*utrm1*utrm2
            eunder = eunder + esite
!
!  Site energy contributions for eunder
!
            siteenergy(i) = siteenergy(i) + esite
          endif
        else
          expuc1 = 0.0_dp
          expuc2 = 0.0_dp
          expuc3 = 0.0_dp
          utrm1 = 0.0_dp
          utrm2 = 0.0_dp
        endif
!**********************************************************
!  Set derivatives of deltalp for use in SBO derivatives  *
!**********************************************************
        if (lgrad1) then
          delta_e = 0.5_dp*(delta(i) + reaxFFval(1,nspeci) - reaxFFval(3,nspeci))
          if (lreaxFFovsmooth) then
!
!  New smooth form
!
            delta_e0 = nint(delta_e)
            expH = exp(-reaxFFlpsmooth*(delta_e - delta_e0))
            smoothH = 1.0_dp/(1.0_dp + expH)
            rintde = delta_e0 - 1.0_dp + smoothH
            drintdeddeltae = reaxFFlpsmooth*expH*smoothH**2
            if (lgrad2) then
              d2rintdeddeltae2 = reaxFFlpsmooth*reaxFFlpsmooth*expH*(smoothH**2)* &
                                 (2.0_dp*expH*smoothH - 1.0_dp)
            endif
          else
!
!  Old discontinuous form
!
            rintde = dble(int(delta_e))
            drintdeddeltae = 0.0_dp
            d2rintdeddeltae2 = 0.0_dp
          endif
!
!  Derivatives of deltalp
!
          trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
          explp2 = exp(-reaxFFlam(29)*trm2lp**2)
          ddeltalpddeltai = 0.5_dp*(drintdeddeltae + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(1.0_dp - drintdeddeltae)) ! 1/2 is to convert from delta_e to delta_i
!
          if (lgrad2) then
            d2deltalpddeltai2 = 0.25_dp*(d2rintdeddeltae2 + reaxFFlam(29)*explp2*( &
                                         8.0_dp*(1.0_dp - drintdeddeltae)**2 - &
                                         4.0_dp*trm2lp*d2rintdeddeltae2 - &
                                         16.0_dp*reaxFFlam(29)*trm2lp*trm2lp*(1.0_dp - drintdeddeltae)**2))
          endif
        endif
!*****************************************************
!  Compute SBO term needed for valence angle energy  *
!*****************************************************
!
!  During pass over neighbours, sbo stores sum of pi bond orders & sboprod is product of exponentiated bond orders
!
        sbo = 0.0_dp
        sbopi = 1.0_dp     ! This value can be used to turn off pi contribution to sbo for debugging
        sboprod = 1.0_dp
        do ni = 1,nneigh(i)
          sbo = sbo + BO_pi(ni,i) + BO_pipi(ni,i)
          bo8 = BO(ni,i)**8
          expbo8 = exp(-bo8)
          sboprod = sboprod*expbo8
        enddo
!
!  Combine sbo contributions
!
        delta_angle = delta(i) + reaxFFval(1,nspeci) - reaxFFval(4,nspeci)
        rnlp = reaxFFlp(1,nspeci) - deltalp(i)
!
        sbo = sbo - (1.0_dp - sboprod)*(delta_angle + reaxFFlam(16)*rnlp)
        sbopi = 1.0_dp
        if (lgrad1) then
          dsboddeltai = - (1.0_dp - sboprod)*(1.0_dp - reaxFFlam(16)*ddeltalpddeltai)
          do nl = 1,nneigh(i)
            dsboproddbo(nl) = - (8.0_dp*BO(nl,i)**7)*sboprod*(delta_angle + reaxFFlam(16)*rnlp)
          enddo
          if (lgrad2) then
            d2sboddeltai2 = (1.0_dp - sboprod)*reaxFFlam(16)*d2deltalpddeltai2
            ind = 0
            do nl = 1,nneigh(i)
              d2sboddeltaidbo(nl) = - (8.0_dp*BO(nl,i)**7)*sboprod*(1.0_dp - reaxFFlam(16)*ddeltalpddeltai)
              do nm = 1,nl
                ind = ind + 1
                d2sboproddbo2(ind) = 64.0_dp*(BO(nl,i)**7)*(BO(nm,i)**7)*sboprod*(delta_angle + reaxFFlam(16)*rnlp)
                if (nl.eq.nm) then
                  d2sboproddbo2(ind) = d2sboproddbo2(ind) - (56.0_dp*BO(nl,i)**6)*sboprod*(delta_angle + reaxFFlam(16)*rnlp)
                endif
              enddo
            enddo
          endif
        endif
!
!  Valence angle and penalty term
!
        if (lgrad1) then
          devaldsbo = 0.0_dp
          devalddeltai = 0.0_dp
          devaldBO_neigh(1:maxneigh1a) = 0.0_dp
          devaldsbo_neigh(1:nneigh(i)) = 0.0_dp
          devalddeltaj_neigh(1:nneigh(i)) = 0.0_dp
          if (lgrad2) then
            d2evaldsbo2 = 0.0_dp
            d2evalddeltai2 = 0.0_dp
            d2evalddeltaidsbo = 0.0_dp
            d2evalddeltaidBO_neigh(1:maxneigh1a) = 0.0_dp
            d2evalddeltaiddeltaj_neigh(1:nneigh(i)) = 0.0_dp
            d2evaldBO2_neigh(1:maxneigh2a) = 0.0_dp
            d2evaldsbo_neigh(1:nneigh(i)) = 0.0_dp
            d2evaldsbo2_neigh(1:nneighi1) = 0.0_dp
            d2evaldthetadBO_neigh(1:maxsparse,1:maxneigh1a) = 0.0_dp
            d2evaldthetadBOpi_neigh(1:nneighi2e) = 0.0_dp
            d2evaldthetaddeltai_neigh(1:nneighi2e) = 0.0_dp
            d2evalddeltajk_neigh(1:nneighi2e) = 0.0_dp
            d2evalddeltajdBO_neigh(1:nneigh(i),1:maxneigh1a) = 0.0_dp
            d2etorsddeltadBOpi_neigh(1:nneigh(i)) = 0.0_dp
            d2etorsddeltaidBOpi_neigh(1:nneigh(i)) = 0.0_dp
            d2etorsddeltajdBOpi_neigh(1:nneigh(i)) = 0.0_dp
          endif
!
!  Set up derivative of delta_lp with respect to delta_e
!
          delta_e = 0.5_dp*(delta(i) + reaxFFval(1,nspeci) - reaxFFval(3,nspeci))
!
          if (lreaxFFlpsmooth) then
!
!  New smooth form
!
            delta_e0 = nint(delta_e)
            expH = exp(-reaxFFlpsmooth*(delta_e - delta_e0))
            smoothH = 1.0_dp/(1.0_dp + expH)
            rintde = delta_e0 - 1.0_dp + smoothH
            drintdeddeltae = reaxFFlpsmooth*expH*smoothH**2
            if (lgrad2) then
              d2rintdeddeltae2 = reaxFFlpsmooth*reaxFFlpsmooth*expH*(smoothH**2)* &
                                 (2.0_dp*expH*smoothH - 1.0_dp)
            endif
          else
!
!  Old discontinuous form
!
            rintde = dble(int(delta_e))
            drintdeddeltae = 0.0_dp
            d2rintdeddeltae2 = 0.0_dp
          endif
!
          trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
          explp2 = exp(-reaxFFlam(29)*trm2lp**2)
          ddeltalplpddeltai = 0.5_dp*(drintdeddeltae + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(1.0_dp - drintdeddeltae)) ! 1/2 is to convert from delta_e to delta_i
!
          if (lgrad2) then
            d2deltalplpddeltai2 = 0.25_dp*(d2rintdeddeltae2 + reaxFFlam(29)*explp2*( &
                                         8.0_dp*(1.0_dp - drintdeddeltae)**2 - &
                                         4.0_dp*trm2lp*d2rintdeddeltae2 - &
                                         16.0_dp*reaxFFlam(29)*trm2lp*trm2lp*(1.0_dp - drintdeddeltae)**2))
          endif
!
          if (lreaxFFovsmooth) then
!
!  New smooth form
!
            delta_e0 = nint(delta_e)
            expH = exp(-reaxFFlpsmooth*(delta_e - delta_e0))
            smoothH = 1.0_dp/(1.0_dp + expH)
            rintde = delta_e0 - 1.0_dp + smoothH
            drintdeddeltae = reaxFFlpsmooth*expH*smoothH**2
            if (lgrad2) then
              d2rintdeddeltae2 = reaxFFlpsmooth*reaxFFlpsmooth*expH*(smoothH**2)* &
                                 (2.0_dp*expH*smoothH - 1.0_dp)
            endif
          else
!
!  Old discontinuous form
!
            rintde = dble(int(delta_e))
            drintdeddeltae = 0.0_dp
            d2rintdeddeltae2 = 0.0_dp
          endif
!
!  Derivatives of deltalp
!
          trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
          explp2 = exp(-reaxFFlam(29)*trm2lp**2)
          ddeltalpddeltai = 0.5_dp*(drintdeddeltae + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(1.0_dp - drintdeddeltae)) ! 1/2 is to convert from delta_e to delta_i
!
          if (lgrad2) then
            d2deltalpddeltai2 = 0.25_dp*(d2rintdeddeltae2 + reaxFFlam(29)*explp2*( &
                                         8.0_dp*(1.0_dp - drintdeddeltae)**2 - &
                                         4.0_dp*trm2lp*d2rintdeddeltae2 - &
                                         16.0_dp*reaxFFlam(29)*trm2lp*trm2lp*(1.0_dp - drintdeddeltae)**2))
          endif
        endif
!
!  Compute delta for 3-body conjugation
!
        delta_coa = delta(i) + reaxFFval(1,nspeci) - reaxFFval(2,nspeci)
!******************
!  Angular terms  *
!******************
!
!  NB: The use of cutoffs is not applied to the hydrogen bond term, as per the original program
!
        do ni = 2,nneigh(i)
          j = neighno(ni,i)
          if (nbos(j).gt.0) then
            nspecj = nbosptr(j)
            rij = rneigh(ni,i)
!
!  Bond order check for i-j
!
            BOij = BO(ni,i) - reaxFFatol
            lBOijOK = (BOij.gt.0.0_dp) 
!
!  Calculate taper function to smooth cutoff for BOij
!
            if (lStrict) then
              fij = 1.0_dp
              dfijdBOij = 0.0_dp
              d2fijdBOij2 = 0.0_dp
            else
              call ataper(.false.,BO(ni,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fij,dfijdBOij,d2fijdBOij2, &
                          d3fijdBOij3,lgrad1,lgrad2,.false.)
            endif
            do nk = 1,ni-1
              k = neighno(nk,i)
              if (nbos(k).gt.0) then
                nspeck = nbosptr(k)
!
!  Bond order check for i-k
!
                BOik = BO(nk,i) - reaxFFatol
                lBOikOK = (BOik.gt.0.0_dp)
!
!  Calculate taper function to smooth cutoff for BOij
!
                if (lStrict) then
                  fik = 1.0_dp
                  dfikdBOik = 0.0_dp
                  d2fikdBOik2 = 0.0_dp
                else
                  call ataper(.false.,BO(nk,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fik,dfikdBOik,d2fikdBOik2, &
                              d3fikdBOik3,lgrad1,lgrad2,.false.)
                endif
!
!  Find triad index
!
                if (nspecj.ge.nspeck) then
                  indsjk = nspecj*(nspecj - 1)/2 + nspeck
                else
                  indsjk = nspeck*(nspeck - 1)/2 + nspecj
                endif
!
!  Compute distance squared between j - k
!
                xjk = xneigh(ni,i) - xneigh(nk,i)
                yjk = yneigh(ni,i) - yneigh(nk,i)
                zjk = zneigh(ni,i) - zneigh(nk,i)
                rjk2 = xjk**2 + yjk**2 + zjk**2
!
!  Compute theta for j - i - k
!
                rjk = sqrt(rjk2)
                call reaxff_theta(rneigh(ni,i),rneigh(nk,i),rjk,theta,dthetadrij,dthetadrik,dthetadrjk, &
                                  d2thetadr2,lgrad1,lgrad2)
                if (lBOijOK.and.lBOikOK.and.(BO(ni,i)*BO(nk,i).gt.reaxFFatol2)) then
!------------------
!  Valence energy |
!------------------
! 
!  Set index for j-k derivatives                   
! 
                  indjk = nneigh(i) + ni*(ni - 1)/2 + nk
!
!  Add indices to sparsity pattern
!
                  if (lgrad1) then
                    if (nonzerorptr(ni).eq.0) then
                      nonzero = nonzero + 1
                      nonzeroptr(nonzero) = ni
                      nonzerorptr(ni) = nonzero
                    endif
                    if (nonzerorptr(nk).eq.0) then
                      nonzero = nonzero + 1
                      nonzeroptr(nonzero) = nk
                      nonzerorptr(nk) = nonzero
                    endif
                    if (nonzerorptr(indjk).eq.0) then
                      nonzero = nonzero + 1
                      nonzeroptr(nonzero) = indjk
                      nonzerorptr(indjk) = nonzero
                    endif
!
!  Check there is room in the sparse arrays and zero any elements added
!
                    if (nonzero.gt.maxsparse) then
                      oldmaxsparse = maxsparse
                      oldmaxsparse2 = maxsparse*(maxsparse+1)/2
                      maxsparse = nonzero + 500
                      call realloc(d1is,maxsparse,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d1is')
                      maxsparse2 = maxsparse*(maxsparse+1)/2
                      call realloc(d2is,maxsparse2,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2is')
                      call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
                      call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
                      call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
!
                      d1is(oldmaxsparse+1:maxsparse) = 0.0_dp
                      d2is(oldmaxsparse2+1:maxsparse2) = 0.0_dp
                      d2evaldthetadBO_neigh(oldmaxsparse+1:maxsparse,1:maxneigh1a) = 0.0_dp
                      d2etorsdthetaddeltaj_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                      d2etorsdthetadBOpi_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                    endif
                  endif
!
!  Loop over three-body valence potentials
!
                  if (lreaxFF_eval) then
                    do nval3 = 1,nreaxFFval3(indsjk,nspeci)
                      pval1 = reaxFFval3(2,nval3,indsjk,nspeci)
                      pval2 = reaxFFval3(3,nval3,indsjk,nspeci)
!
!  Compute f7, f8 and theta0 for this triad
!
                      call reaxFF_f7(nspeci,nspecj,nspeck,nval3,BOij,BOik,f7,df7dBOij,df7dBOik, &
                                     d2f7dBOij2,d2f7dBOijdBOik,d2f7dBOik2,lgrad1,lgrad2)
                      call reaxFF_f8(nspeci,nspecj,nspeck,nval3,delta(i),f8,df8ddeltai,d2f8ddeltai2,lgrad1,lgrad2)
                      call reaxFF_theta0(nspeci,nspecj,nspeck,nval3,sbo,theta0,dtheta0dsbo,d2theta0dsbo2,lgrad1,lgrad2)
!
!  Set angle dependence
!
                      exptheta = exp(-pval2*(theta0 - theta)**2)
                      fangle = (1.0_dp - exptheta)
                      if (lgrad1) then
                        dfangledtheta = - 2.0_dp*pval2*(theta0 - theta)*exptheta
                        dfangledtheta0 = - dfangledtheta
                        if (lgrad2) then
                          d2fangledtheta2 = 2.0_dp*pval2*exptheta*(1.0_dp - 2.0_dp*pval2*(theta0 - theta)**2)
                          d2fangledtheta02 = d2fangledtheta2
                          d2fangledthetadtheta0 = - d2fangledtheta2
                        endif
                      endif
                      if (lilocal) then
!
!  Evaluate valence energy
!
                        eval = eval + fij*fik*f7*f8*pval1*fangle
!
!  Site energy for eval
!
                        esite = fij*fik*f7*f8*pval1*fangle/3.0_dp
                        siteenergy(i) = siteenergy(i) + esite
                        siteenergy(j) = siteenergy(j) + esite
                        siteenergy(k) = siteenergy(k) + esite
                      endif
!
!  Derivatives of valence energy
!
                      if (lgrad1) then
!  
!  Compute derivative terms from f7 and f8
!                 
                        devalddeltai = devalddeltai + pval1*fij*fik*f7*df8ddeltai*fangle
                        devaldBO_neigh(ni) = devaldBO_neigh(ni) + pval1*fij*fik*f8*df7dBOij*fangle
                        devaldBO_neigh(ni) = devaldBO_neigh(ni) + pval1*dfijdBOij*fik*f8*f7*fangle
                        devaldBO_neigh(nk) = devaldBO_neigh(nk) + pval1*fij*fik*f8*df7dBOik*fangle
                        devaldBO_neigh(nk) = devaldBO_neigh(nk) + pval1*fij*dfikdBOik*f8*f7*fangle
!                 
!  Compute derivative terms from (1 - exp(-pval2*(theta0 - theta)**2) - Part 1 derivatives for theta
!                 
                        devaldtheta = pval1*fij*fik*f7*f8*dfangledtheta
                        ind1 = nonzerorptr(ni)
                        ind2 = nonzerorptr(nk)
                        ind3 = nonzerorptr(indjk)
                        d1is(ind1) = d1is(ind1) + devaldtheta*dthetadrij
                        d1is(ind2) = d1is(ind2) + devaldtheta*dthetadrik
                        d1is(ind3) = d1is(ind3) + devaldtheta*dthetadrjk
!
!  Compute derivative terms from (1 - exp(-pval2*(theta0 - theta)**2) - Part 2 derivatives for theta0
!               
                        devaldtheta0 = pval1*fij*fik*f7*f8*dfangledtheta0
                        devaldsbo = devaldsbo + devaldtheta0*dtheta0dsbo
!
                        devalddeltai = devalddeltai + devaldtheta0*dtheta0dsbo*dsboddeltai
                        do nl = 1,nneigh(i)
                          devaldsbo_neigh(nl) = devaldsbo_neigh(nl) + devaldtheta0*dtheta0dsbo*dsboproddbo(nl)
                        enddo
!
                        if (lgrad2) then
!
!  Second derivatives
!
                          d2evalddeltai2 = d2evalddeltai2 + pval1*fij*fik*f7*d2f8ddeltai2*fangle
!
                          ind = ni*(ni+1)/2
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + pval1*fij*fik*f8*d2f7dBOij2*fangle
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*pval1*dfijdBOij*fik*f8*df7dBOij*fangle
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + pval1*d2fijdBOij2*fik*f8*f7*fangle
!
                          ind = nk*(nk+1)/2
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + pval1*fij*fik*f8*d2f7dBOik2*fangle
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*pval1*fij*dfikdBOik*f8*df7dBOik*fangle
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + pval1*fij*d2fikdBOik2*f8*f7*fangle
!
                          ind = ni*(ni-1)/2 + nk
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + pval1*dfijdBOij*dfikdBOik*f8*f7*fangle
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + pval1*dfijdBOij*fik*f8*df7dBOik*fangle
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + pval1*fij*dfikdBOik*f8*df7dBOij*fangle
                          d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + pval1*fij*fik*d2f7dBOijdBOik*f8*fangle
!
                          d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                                                       pval1*fij*fik*df8ddeltai*df7dBOij*fangle
                          d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                                                       pval1*dfijdBOij*fik*df8ddeltai*f7*fangle
                          d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                                                       pval1*fij*fik*df8ddeltai*df7dBOik*fangle
                          d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                                                       pval1*fij*dfikdBOik*df8ddeltai*f7*fangle
!
!  Compute derivative terms from (1 - exp(-pval2*(theta0 - theta)**2) - Part 1 derivatives for theta
!
                          d2evaldtheta2 = pval1*fij*fik*f7*f8*d2fangledtheta2
                          ind = ind1*(ind1+1)/2
                          d2is(ind) = d2is(ind) + devaldtheta*d2thetadr2(1) + d2evaldtheta2*dthetadrij**2
                          if (ind1.gt.ind2) then
                            ind = ind1*(ind1-1)/2 + ind2
                          else
                            ind = ind2*(ind2-1)/2 + ind1
                          endif
                          d2is(ind) = d2is(ind) + devaldtheta*d2thetadr2(2) + d2evaldtheta2*dthetadrij*dthetadrik
                          if (ind1.gt.ind3) then
                            ind = ind1*(ind1-1)/2 + ind3
                          else
                            ind = ind3*(ind3-1)/2 + ind1
                          endif
                          d2is(ind) = d2is(ind) + devaldtheta*d2thetadr2(3) + d2evaldtheta2*dthetadrij*dthetadrjk
                          ind = ind2*(ind2+1)/2
                          d2is(ind) = d2is(ind) + devaldtheta*d2thetadr2(4) + d2evaldtheta2*dthetadrik**2
                          if (ind2.gt.ind3) then
                            ind = ind2*(ind2-1)/2 + ind3
                          else
                            ind = ind3*(ind3-1)/2 + ind2
                          endif
                          d2is(ind) = d2is(ind) + devaldtheta*d2thetadr2(5) + d2evaldtheta2*dthetadrik*dthetadrjk
                          ind = ind3*(ind3+1)/2
                          d2is(ind) = d2is(ind) + devaldtheta*d2thetadr2(6) + d2evaldtheta2*dthetadrjk**2
!
!  Compute derivative terms from (1 - exp(-pval2*(theta0 - theta)**2) - Part 2 derivatives for theta0
!
                          d2evaldtheta02 = pval1*fij*fik*f7*f8*d2fangledtheta02
                          d2evaldsbo2 = d2evaldsbo2 + d2evaldtheta02*dtheta0dsbo**2 + devaldtheta0*d2theta0dsbo2
!
                          d2evalddeltai2 = d2evalddeltai2 + d2evaldtheta02*(dtheta0dsbo*dsboddeltai)**2 &
                                                          + devaldtheta0*d2theta0dsbo2*dsboddeltai**2 &
                                                          + devaldtheta0*dtheta0dsbo*d2sboddeltai2
                          d2evalddeltaidsbo = d2evalddeltaidsbo + d2evaldtheta02*(dtheta0dsbo**2)*dsboddeltai &
                                                                + devaldtheta0*d2theta0dsbo2*dsboddeltai
!
                          ind = 0
                          do nl = 1,nneigh(i)
                            do nm = 1,nl
                              ind = ind + 1
                              d2evaldsbo2_neigh(ind) = d2evaldsbo2_neigh(ind) + devaldtheta0*dtheta0dsbo*d2sboproddbo2(ind) &
                                                                              + (devaldtheta0*d2theta0dsbo2 + &
                                                                                 d2evaldtheta02*dtheta0dsbo**2)* &
                                                                                 dsboproddbo(nl)*dsboproddbo(nm)
                            enddo
                          enddo
!
!  Cross terms for theta0 from sbo
!
                          do nl = 1,nneigh(i)
                            d2evalddeltaidBO_neigh(nl) = d2evalddeltaidBO_neigh(nl) + &
                                                         d2evaldtheta02*(dtheta0dsbo**2)*dsboproddbo(nl)*dsboddeltai + &
                                                         devaldtheta0*d2theta0dsbo2*dsboproddbo(nl)*dsboddeltai + &
                                                         devaldtheta0*dtheta0dsbo*d2sboddeltaidbo(nl)
!
                            d2evaldsbo_neigh(nl) = d2evaldsbo_neigh(nl) + d2evaldtheta02*(dtheta0dsbo**2)*dsboproddbo(nl) &
                                                                        + devaldtheta0*d2theta0dsbo2*dsboproddbo(nl)
                          enddo
!
!  Cross terms between fij, fik and f7 with theta 
!
                          d2evaldfijdtheta = pval1*fik*f7*f8*dfangledtheta
                          d2evaldfikdtheta = pval1*fij*f7*f8*dfangledtheta
                          d2evaldf7dtheta = pval1*fij*fik*f8*dfangledtheta
!
                          d2evaldthetadBO_neigh(ind1,ni) = d2evaldthetadBO_neigh(ind1,ni) + &
                            d2evaldfijdtheta*dthetadrij*dfijdBOij + &
                            d2evaldf7dtheta*dthetadrij*df7dBOij
                          d2evaldthetadBO_neigh(ind2,ni) = d2evaldthetadBO_neigh(ind2,ni) + &
                            d2evaldfijdtheta*dthetadrik*dfijdBOij + &
                            d2evaldf7dtheta*dthetadrik*df7dBOij
                          d2evaldthetadBO_neigh(ind3,ni) = d2evaldthetadBO_neigh(ind3,ni) + &
                            d2evaldfijdtheta*dthetadrjk*dfijdBOij + &
                            d2evaldf7dtheta*dthetadrjk*df7dBOij
                          d2evaldthetadBO_neigh(ind1,nk) = d2evaldthetadBO_neigh(ind1,nk) + &
                            d2evaldfikdtheta*dthetadrij*dfikdBOik + &
                            d2evaldf7dtheta*dthetadrij*df7dBOik
                          d2evaldthetadBO_neigh(ind2,nk) = d2evaldthetadBO_neigh(ind2,nk) + &
                            d2evaldfikdtheta*dthetadrik*dfikdBOik + &
                            d2evaldf7dtheta*dthetadrik*df7dBOik
                          d2evaldthetadBO_neigh(ind3,nk) = d2evaldthetadBO_neigh(ind3,nk) + &
                            d2evaldfikdtheta*dthetadrjk*dfikdBOik + &
                            d2evaldf7dtheta*dthetadrjk*df7dBOik
!
!  Cross terms between fij, fik and f7 with theta0
!
                          d2evaldfijdtheta0 = pval1*fik*f7*f8*dfangledtheta0
                          d2evaldfikdtheta0 = pval1*fij*f7*f8*dfangledtheta0
                          d2evaldf7dtheta0 = pval1*fij*fik*f8*dfangledtheta0
!
                          d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                            d2evaldfijdtheta0*dtheta0dsbo*dsboddeltai*dfijdBOij + &
                            d2evaldf7dtheta0*dtheta0dsbo*dsboddeltai*df7dBOij
                          d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                            d2evaldfikdtheta0*dtheta0dsbo*dsboddeltai*dfikdBOik + &
                            d2evaldf7dtheta0*dtheta0dsbo*dsboddeltai*df7dBOik
!
                          d2evaldsbo_neigh(ni) = d2evaldsbo_neigh(ni) + &
                            d2evaldfijdtheta0*dtheta0dsbo*dfijdBOij + &
                            d2evaldf7dtheta0*dtheta0dsbo*df7dBOij
                          d2evaldsbo_neigh(nk) = d2evaldsbo_neigh(nk) + &
                            d2evaldfikdtheta0*dtheta0dsbo*dfikdBOik + &
                            d2evaldf7dtheta0*dtheta0dsbo*df7dBOik
!
                          do nl = 1,nneigh(i)
                            if (ni.ge.nl) then
                              ind = ni*(ni - 1)/2 + nl
                            else
                              ind = nl*(nl - 1)/2 + ni
                            endif
!
                            if (ni.eq.nl) then
                              d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                                2.0_dp*d2evaldfijdtheta0*dtheta0dsbo*dsboproddbo(nl)*dfijdBOij + &
                                2.0_dp*d2evaldf7dtheta0*dtheta0dsbo*dsboproddbo(nl)*df7dBOij
                            else
                              d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                                d2evaldfijdtheta0*dtheta0dsbo*dsboproddbo(nl)*dfijdBOij + &
                                d2evaldf7dtheta0*dtheta0dsbo*dsboproddbo(nl)*df7dBOij
                            endif
!
                            if (nk.ge.nl) then
                              ind = nk*(nk - 1)/2 + nl
                            else
                              ind = nl*(nl - 1)/2 + nk
                            endif
!
                            if (nk.eq.nl) then
                              d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                                2.0_dp*d2evaldfikdtheta0*dtheta0dsbo*dsboproddbo(nl)*dfikdBOik + &
                                2.0_dp*d2evaldf7dtheta0*dtheta0dsbo*dsboproddbo(nl)*df7dBOik
                            else
                              d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                                d2evaldfikdtheta0*dtheta0dsbo*dsboproddbo(nl)*dfikdBOik + &
                                d2evaldf7dtheta0*dtheta0dsbo*dsboproddbo(nl)*df7dBOik
                            endif
                          enddo
!
!  Cross terms between theta and theta0
!
                          d2evaldthetadtheta0 = pval1*fij*fik*f7*f8*d2fangledthetadtheta0
!
                          d2evaldthetadBOpi_neigh(ni) = d2evaldthetadBOpi_neigh(ni) + &
                            d2evaldthetadtheta0*dtheta0dsbo*dthetadrij
                          d2evaldthetadBOpi_neigh(nk) = d2evaldthetadBOpi_neigh(nk) + &
                            d2evaldthetadtheta0*dtheta0dsbo*dthetadrik
                          d2evaldthetadBOpi_neigh(indjk) = d2evaldthetadBOpi_neigh(indjk) + &
                            d2evaldthetadtheta0*dtheta0dsbo*dthetadrjk
!
                          d2evaldthetaddeltai_neigh(ni) = d2evaldthetaddeltai_neigh(ni) + &
                            d2evaldthetadtheta0*dtheta0dsbo*dthetadrij*dsboddeltai
                          d2evaldthetaddeltai_neigh(nk) = d2evaldthetaddeltai_neigh(nk) + &
                            d2evaldthetadtheta0*dtheta0dsbo*dthetadrik*dsboddeltai
                          d2evaldthetaddeltai_neigh(indjk) = d2evaldthetaddeltai_neigh(indjk) + &
                            d2evaldthetadtheta0*dtheta0dsbo*dthetadrjk*dsboddeltai
!
                          do nl = 1,nneigh(i)
                            d2evaldthetadBO_neigh(ind1,nl) = d2evaldthetadBO_neigh(ind1,nl) + &
                              d2evaldthetadtheta0*dtheta0dsbo*dthetadrij*dsboproddbo(nl)
                            d2evaldthetadBO_neigh(ind2,nl) = d2evaldthetadBO_neigh(ind2,nl) + &
                              d2evaldthetadtheta0*dtheta0dsbo*dthetadrik*dsboproddbo(nl)
                            d2evaldthetadBO_neigh(ind3,nl) = d2evaldthetadBO_neigh(ind3,nl) + &
                              d2evaldthetadtheta0*dtheta0dsbo*dthetadrjk*dsboproddbo(nl)
                          enddo
!
!  Cross terms between f8 and theta/theta0 terms
!
                          d2evalddeltai2 = d2evalddeltai2 + &
                            2.0_dp*pval1*fij*fik*f7*df8ddeltai*dfangledtheta0*dtheta0dsbo*dsboddeltai
                          d2evalddeltaidsbo = d2evalddeltaidsbo + &
                            pval1*fij*fik*f7*df8ddeltai*dfangledtheta0*dtheta0dsbo
                          do nl = 1,nneigh(i)
                            d2evalddeltaidBO_neigh(nl) = d2evalddeltaidBO_neigh(nl) + &
                              pval1*fij*fik*f7*df8ddeltai*dfangledtheta0*dtheta0dsbo*dsboproddbo(nl)
                          enddo
!
                          d2evaldthetaddeltai = pval1*fij*fik*f7*df8ddeltai*dfangledtheta
                          d2evaldthetaddeltai_neigh(ni) = d2evaldthetaddeltai_neigh(ni) + &
                            d2evaldthetaddeltai*dthetadrij
                          d2evaldthetaddeltai_neigh(nk) = d2evaldthetaddeltai_neigh(nk) + &
                            d2evaldthetaddeltai*dthetadrik
                          d2evaldthetaddeltai_neigh(indjk) = d2evaldthetaddeltai_neigh(indjk) + &
                            d2evaldthetaddeltai*dthetadrjk
                        endif
                      endif
!
!  End of loop over three-body valence potentials
!
                    enddo
                  endif
!------------------
!  Penalty energy |
!------------------
                  lDoEpen = (lreaxFF_epen.and.abs(reaxFFpen3(indsjk,nspeci)).gt.1.0d-12)
                  if (lDoEpen) then
!
!  Compute f9 for this triad
!
                    call reaxFF_f9(delta(i),f9,df9ddeltai,d2f9ddeltai2,lgrad1,lgrad2)
!
!  Evaluate penalty energy
!
                    ppen2 = reaxFFlam(20)
                    exppen1 = exp(-ppen2*(BOij - 2.0_dp)**2)
                    exppen2 = exp(-ppen2*(BOik - 2.0_dp)**2)
                    epentrm = reaxFFpen3(indsjk,nspeci)*f9*exppen1*exppen2
                    if (lilocal) then
                      esite = fij*fik*epentrm
                      epen = epen + esite
!
!  Site energy for epen
!
                      esite = esite/3.0_dp
                      siteenergy(i) = siteenergy(i) + esite
                      siteenergy(j) = siteenergy(j) + esite
                      siteenergy(k) = siteenergy(k) + esite
                    endif
!
                    if (lgrad1) then
!  
!  Derivatives of penalty energy with respect to deltai & bond orders -> add to eval terms
!
                      d1exppen1 = - 2.0_dp*ppen2*(BOij - 2.0_dp)*exppen1
                      d1exppen2 = - 2.0_dp*ppen2*(BOik - 2.0_dp)*exppen2
!
                      devalddeltai = devalddeltai + fij*fik*reaxFFpen3(indsjk,nspeci)*df9ddeltai*exppen1*exppen2
!
                      devaldBO_neigh(ni) = devaldBO_neigh(ni) + fij*fik*reaxFFpen3(indsjk,nspeci)*f9*d1exppen1*exppen2
                      devaldBO_neigh(ni) = devaldBO_neigh(ni) + dfijdBOij*fik*epentrm
!
                      devaldBO_neigh(nk) = devaldBO_neigh(nk) + fij*fik*reaxFFpen3(indsjk,nspeci)*f9*exppen1*d1exppen2
                      devaldBO_neigh(nk) = devaldBO_neigh(nk) + dfikdBOik*fij*epentrm
!
                      if (lgrad2) then
                        d2exppen1 = (4.0_dp*ppen2*(BOij - 2.0_dp)**2 - 2.0_dp)*ppen2*exppen1
                        d2exppen2 = (4.0_dp*ppen2*(BOik - 2.0_dp)**2 - 2.0_dp)*ppen2*exppen2
!
                        d2evalddeltai2 = d2evalddeltai2 + fij*fik*reaxFFpen3(indsjk,nspeci)*d2f9ddeltai2*exppen1*exppen2
!
                        d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                          dfijdBOij*fik*reaxFFpen3(indsjk,nspeci)*df9ddeltai*exppen1*exppen2
                        d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                          fij*fik*reaxFFpen3(indsjk,nspeci)*df9ddeltai*d1exppen1*exppen2
                        d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                          fij*dfikdBOik*reaxFFpen3(indsjk,nspeci)*df9ddeltai*exppen1*exppen2
                        d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                          fij*fik*reaxFFpen3(indsjk,nspeci)*df9ddeltai*exppen1*d1exppen2
!
                        ind = ni*(ni+1)/2
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          d2fijdBOij2*fik*epentrm
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          2.0_dp*dfijdBOij*fik*reaxFFpen3(indsjk,nspeci)*f9*d1exppen1*exppen2
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          fij*fik*reaxFFpen3(indsjk,nspeci)*f9*d2exppen1*exppen2
!
                        ind = nk*(nk+1)/2
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          fij*d2fikdBOik2*epentrm
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          2.0_dp*fij*dfikdBOik*reaxFFpen3(indsjk,nspeci)*f9*exppen1*d1exppen2
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          fij*fik*reaxFFpen3(indsjk,nspeci)*f9*exppen1*d2exppen2
!
                        ind = ni*(ni-1)/2 + nk
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          dfijdBOij*dfikdBOik*epentrm
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          dfijdBOij*fik*reaxFFpen3(indsjk,nspeci)*f9*exppen1*d1exppen2
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          fij*dfikdBOik*reaxFFpen3(indsjk,nspeci)*f9*d1exppen1*exppen2
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + &
                          fij*fik*reaxFFpen3(indsjk,nspeci)*f9*d1exppen1*d1exppen2
                      endif
                    endif
                  endif
!-----------------------------
!  3-body conjugation energy |
!-----------------------------
!
!  Set parameters for ecoa
!
                  pcoa1 = reaxFFconj3(1,indsjk,nspeci)
!
!  If pcoa1 is zero then there is no potential to do
!
                  if (lreaxFF_ecoa.and.abs(pcoa1).gt.1.0d-12) then
                    pcoa2 = reaxFFconj3(2,indsjk,nspeci)
                    pcoa3 = reaxFFconj3(3,indsjk,nspeci)
                    pcoa4 = reaxFFconj3(4,indsjk,nspeci)
!
!  Compute sum of bond orders for terminal atoms excluding i-j / i-k contribution
!  => use the fact that delta is already the sum of all bond orders
!
                    sumBOcoa_ij = delta(j) + reaxFFval(1,nspecj) - BOij
                    sumBOcoa_ik = delta(k) + reaxFFval(1,nspeck) - BOik
!
!  Evaluate 3-body conjugation energy
!
                    expcoa2 = exp(pcoa2*delta_coa)
                    expcoatrm = 1.0_dp/(1.0_dp + expcoa2)
!
                    expcoa3a = exp(-pcoa3*sumBOcoa_ij**2)
                    expcoa3b = exp(-pcoa3*sumBOcoa_ik**2)
                    expcoa4a = exp(-pcoa4*(BOij - 1.5_dp)**2)
                    expcoa4b = exp(-pcoa4*(BOik - 1.5_dp)**2)
!
                    expcoaprod = pcoa1*expcoa3a*expcoa3b*expcoa4a*expcoa4b
                    ecoatrm = expcoaprod*expcoatrm
                    if (lilocal) then
                      esite = fij*fik*ecoatrm
                      ecoa = ecoa + esite
!
!  Site energy for ecoa
!
                      esite = esite/3.0_dp
                      siteenergy(i) = siteenergy(i) + esite
                      siteenergy(j) = siteenergy(j) + esite
                      siteenergy(k) = siteenergy(k) + esite
                    endif
!
                    if (lgrad1) then
                      d1expcoa3addeltaj = - 2.0_dp*pcoa3*sumBOcoa_ij*ecoatrm
                      d1expcoa3adBOij = 2.0_dp*pcoa3*sumBOcoa_ij*ecoatrm
                      d1expcoa3bddeltak = - 2.0_dp*pcoa3*sumBOcoa_ik*ecoatrm
                      d1expcoa3bdBOik = 2.0_dp*pcoa3*sumBOcoa_ik*ecoatrm
                      d1expcoa4a = - 2.0_dp*pcoa4*(BOij - 1.5_dp)*ecoatrm
                      d1expcoa4b = - 2.0_dp*pcoa4*(BOik - 1.5_dp)*ecoatrm
!
!  Derivatives of 3-body conjugation energy with respect to deltai & bond orders
!
                      dexpcoatrmddeltai = - pcoa2*expcoa2*expcoatrm   ! Divided by expcoatrm since this is part of ecoatrm
                      devalddeltai = devalddeltai + fij*fik*ecoatrm*dexpcoatrmddeltai
!
                      devaldBO_neigh(ni) = devaldBO_neigh(ni) + fij*fik*d1expcoa4a
                      devaldBO_neigh(ni) = devaldBO_neigh(ni) + dfijdBOij*fik*ecoatrm
                      devaldBO_neigh(ni) = devaldBO_neigh(ni) + fij*fik*d1expcoa3adBOij
!
                      devaldBO_neigh(nk) = devaldBO_neigh(nk) + fij*fik*d1expcoa4b
                      devaldBO_neigh(nk) = devaldBO_neigh(nk) + dfikdBOik*fij*ecoatrm
                      devaldBO_neigh(nk) = devaldBO_neigh(nk) + fij*fik*d1expcoa3bdBOik
!
!  Derivatives of 3-body conjugation energy with respect to deltaj/k
!
                      devalddeltaj_neigh(ni) = devalddeltaj_neigh(ni) + fij*fik*d1expcoa3addeltaj
                      devalddeltaj_neigh(nk) = devalddeltaj_neigh(nk) + fij*fik*d1expcoa3bddeltak
!
                      if (lgrad2) then
!
!  Second derivatives of 3-body conjugation energy with respect to deltai & bond orders
!
                        d2expcoatrmddeltai2 = - dexpcoatrmddeltai*pcoa2*(2.0_dp*expcoa2*expcoatrm - 1.0_dp)
                        d2evalddeltai2 = d2evalddeltai2 + fij*fik*ecoatrm*d2expcoatrmddeltai2
!
                        d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                          dfijdBOij*fik*dexpcoatrmddeltai*ecoatrm
                        d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                          fij*fik*dexpcoatrmddeltai*d1expcoa4a
                        d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                          fij*fik*dexpcoatrmddeltai*d1expcoa3adBOij
                        d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                          fij*dfikdBOik*dexpcoatrmddeltai*ecoatrm
                        d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                          fij*fik*dexpcoatrmddeltai*d1expcoa4b
                        d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                          fij*fik*dexpcoatrmddeltai*d1expcoa3bdBOik
!
                        d2evalddeltaiddeltaj_neigh(ni) = d2evalddeltaiddeltaj_neigh(ni) + &
                          fij*fik*dexpcoatrmddeltai*d1expcoa3addeltaj
                        d2evalddeltaiddeltaj_neigh(nk) = d2evalddeltaiddeltaj_neigh(nk) + &
                          fij*fik*dexpcoatrmddeltai*d1expcoa3bddeltak
!
                        d2expcoa3addeltaj2 = (4.0_dp*pcoa3*sumBOcoa_ij**2 - 2.0_dp)*pcoa3*ecoatrm
                        d2expcoa3addeltajdBOij = - (4.0_dp*pcoa3*sumBOcoa_ij**2 - 2.0_dp)*pcoa3*ecoatrm
                        d2expcoa3adBOij2 = (4.0_dp*pcoa3*sumBOcoa_ij**2 - 2.0_dp)*pcoa3*ecoatrm
                        d2expcoa3bddeltak2 = (4.0_dp*pcoa3*sumBOcoa_ik**2 - 2.0_dp)*pcoa3*ecoatrm
                        d2expcoa3bddeltakdBOik = - (4.0_dp*pcoa3*sumBOcoa_ik**2 - 2.0_dp)*pcoa3*ecoatrm
                        d2expcoa3bdBOik2 = (4.0_dp*pcoa3*sumBOcoa_ik**2 - 2.0_dp)*pcoa3*ecoatrm
                        d2expcoa3abdBOdBO = 4.0_dp*pcoa3*pcoa3*sumBOcoa_ij*sumBOcoa_ik*ecoatrm
                        d2expcoa3abddeltadBO = - 4.0_dp*pcoa3*pcoa3*sumBOcoa_ij*sumBOcoa_ik*ecoatrm
                        d2expcoa3abddeltajddeltak = 4.0_dp*pcoa3*pcoa3*sumBOcoa_ij*sumBOcoa_ik*ecoatrm
!
                        d2expcoa4a = (4.0_dp*pcoa4*(BOij - 1.5_dp)**2 - 2.0_dp)*pcoa4*ecoatrm
                        d2expcoa4b = (4.0_dp*pcoa4*(BOik - 1.5_dp)**2 - 2.0_dp)*pcoa4*ecoatrm
                        d2expcoa4ab = 4.0_dp*pcoa4*pcoa4*(BOij - 1.5_dp)*(BOik - 1.5_dp)*ecoatrm
!
                        d2expcoa3a4addeltaj = 4.0_dp*pcoa3*pcoa4*sumBOcoa_ij*(BOij - 1.5_dp)*ecoatrm
                        d2expcoa3a4adBOij = - 4.0_dp*pcoa3*pcoa4*sumBOcoa_ij*(BOij - 1.5_dp)*ecoatrm
                        d2expcoa3a4bddeltaj = 4.0_dp*pcoa3*pcoa4*sumBOcoa_ij*(BOik - 1.5_dp)*ecoatrm
                        d2expcoa3a4bdBOij = - 4.0_dp*pcoa3*pcoa4*sumBOcoa_ij*(BOik - 1.5_dp)*ecoatrm
                        d2expcoa3b4addeltak = 4.0_dp*pcoa3*pcoa4*sumBOcoa_ik*(BOij - 1.5_dp)*ecoatrm
                        d2expcoa3b4adBOik = - 4.0_dp*pcoa3*pcoa4*sumBOcoa_ik*(BOij - 1.5_dp)*ecoatrm
                        d2expcoa3b4bddeltak = 4.0_dp*pcoa3*pcoa4*sumBOcoa_ik*(BOik - 1.5_dp)*ecoatrm
                        d2expcoa3b4bdBOik = - 4.0_dp*pcoa3*pcoa4*sumBOcoa_ik*(BOik - 1.5_dp)*ecoatrm
!
                        ind = ni*(ni+1)/2
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + d2fijdBOij2*fik*ecoatrm
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*dfijdBOij*fik*d1expcoa4a
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*dfijdBOij*fik*d1expcoa3adBOij
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fij*fik*d2expcoa3a4adBOij
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*d2expcoa4a
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*d2expcoa3adBOij2
!
                        ind = nk*(nk+1)/2
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*d2fikdBOik2*ecoatrm
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fij*dfikdBOik*d1expcoa4b
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fij*dfikdBOik*d1expcoa3bdBOik
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fij*fik*d2expcoa3b4bdBOik
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*d2expcoa4b
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*d2expcoa3bdBOik2
!
                        ind = ni*(ni-1)/2 + nk
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + dfijdBOij*dfikdBOik*ecoatrm
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + dfijdBOij*fik*d1expcoa4b
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + dfijdBOij*fik*d1expcoa3bdBOik
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*d2expcoa3a4bdBOij
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*dfikdBOik*d1expcoa4a
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*dfikdBOik*d1expcoa3adBOij
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*d2expcoa3b4adBOik
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*d2expcoa4ab
                        d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*d2expcoa3abdBOdBO
!
!  Second derivatives of 3-body conjugation energy with respect to deltaj/k
!
                        ind = ni*(ni+1)/2
                        d2evalddeltajk_neigh(ind) = d2evalddeltajk_neigh(ind) + fij*fik*d2expcoa3addeltaj2
                        ind = nk*(nk+1)/2
                        d2evalddeltajk_neigh(ind) = d2evalddeltajk_neigh(ind) + fij*fik*d2expcoa3bddeltak2
                        ind = ni*(ni-1)/2 + nk
                        d2evalddeltajk_neigh(ind) = d2evalddeltajk_neigh(ind) + fij*fik*d2expcoa3abddeltajddeltak
!
                        d2evalddeltajdBO_neigh(ni,ni) = d2evalddeltajdBO_neigh(ni,ni) + fij*fik*d2expcoa3addeltajdBOij
                        d2evalddeltajdBO_neigh(ni,nk) = d2evalddeltajdBO_neigh(ni,nk) + fij*fik*d2expcoa3abddeltadBO
                        d2evalddeltajdBO_neigh(nk,ni) = d2evalddeltajdBO_neigh(nk,ni) + fij*fik*d2expcoa3abddeltadBO
                        d2evalddeltajdBO_neigh(nk,nk) = d2evalddeltajdBO_neigh(nk,nk) + fij*fik*d2expcoa3bddeltakdBOik
!
                        d2evalddeltajdBO_neigh(ni,ni) = d2evalddeltajdBO_neigh(ni,ni) + fij*fik*d2expcoa3a4addeltaj
                        d2evalddeltajdBO_neigh(ni,nk) = d2evalddeltajdBO_neigh(ni,nk) + fij*fik*d2expcoa3a4bddeltaj
                        d2evalddeltajdBO_neigh(nk,ni) = d2evalddeltajdBO_neigh(nk,ni) + fij*fik*d2expcoa3b4addeltak
                        d2evalddeltajdBO_neigh(nk,nk) = d2evalddeltajdBO_neigh(nk,nk) + fij*fik*d2expcoa3b4bddeltak
!
                        d2evalddeltajdBO_neigh(ni,ni) = d2evalddeltajdBO_neigh(ni,ni) + dfijdBOij*fik*d1expcoa3addeltaj
                        d2evalddeltajdBO_neigh(ni,nk) = d2evalddeltajdBO_neigh(ni,nk) + fij*dfikdBOik*d1expcoa3addeltaj
                        d2evalddeltajdBO_neigh(nk,ni) = d2evalddeltajdBO_neigh(nk,ni) + dfijdBOij*fik*d1expcoa3bddeltak
                        d2evalddeltajdBO_neigh(nk,nk) = d2evalddeltajdBO_neigh(nk,nk) + fij*dfikdBOik*d1expcoa3bddeltak
                      endif
                    endif
                  endif
                endif
              endif
            enddo
          endif
        enddo
!********************
!  Torsional terms  *
!********************
!
!  Initialise arrays for summing the derivatives of torsional energy with respect to bond orders of neighbours of i
!
        if (lgrad1) then
          detorsdBOpi_neigh(1:nneigh(i)) = 0.0_dp
          if (lgrad2) then
            d2etorsdBOpi2_neigh(1:nneigh(i)) = 0.0_dp
            d2etorsdBOdBOpi_neigh(1:nneighi1a,1:nneigh(i)) = 0.0_dp
            d2etorsdthetadBOpi_neigh(1:maxsparse,1:nneigh(i)) = 0.0_dp
            d2etorsdthetaddeltai(1:nneighi2e) = 0.0_dp
            d2etorsdthetaddeltaj_neigh(1:maxsparse,1:nneigh(i)) = 0.0_dp
          endif
        endif
!
!  Loop over neighbours of i to find central bonds 
!
        tniloop: do ni = 1,nneigh(i)
          j = neighno(ni,i)
!
!  Impose conditions that j must be a reaxFF atom 
!  NB: Trying to use i > j algorithm which differs from original method
!
          if (j.lt.i) cycle tniloop
!
          if (nbos(j).gt.0) then
!
!  Check bond order for i-j
!
            BOij = BO(ni,i) - reaxFFatol
            if (BOij.gt.0.0_dp) then
              if (lStrict) then                    
                fij = 1.0_dp                       
                dfijdBOij = 0.0_dp                 
                d2fijdBOij2 = 0.0_dp                 
              else                                 
                call ataper(.false.,BO(ni,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fij,dfijdBOij,d2fijdBOij2,d3fijdBOij3, &
                            lgrad1,lgrad2,.false.)
              endif
              nspecj = nbosptr(j)
              if (nspeci.ge.nspecj) then
                ind1 = nspeci*(nspeci - 1)/2 + nspecj
              else
                ind1 = nspecj*(nspecj - 1)/2 + nspeci
              endif
!
!  Loop over neighbours of i to find terminal atom k while excluding k = j
!
              do nk = 1,nneigh(i)
                k = neighno(nk,i)
!  The condition below modified from k.ne.i to handle small cell systems
                if (nbos(k).gt.0.and.nk.ne.ni) then
!
!  Check bond order for i-k
!
                  BOik = BO(nk,i) - reaxFFatol
                  if (BOik.gt.0.0_dp) then
                    if (lStrict) then
                      fik = 1.0_dp
                      dfikdBOik = 0.0_dp
                      d2fikdBOik2 = 0.0_dp
                    else
                      call ataper(.false.,BO(nk,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fik,dfikdBOik,d2fikdBOik2, &
                                  d3fikdBOik3,lgrad1,lgrad2,.false.)
                    endif
                    nspeck = nbosptr(k)
!
!  Loop over neighbours of j to find terminal atom l while excluding i = l and k = l
!
                    do nl = 1,nneigh(j)
                      l = neighno(nl,j)
!
!  Check bond order for j-l
!
                      BOjl = BO(nl,j) - reaxFFatol
                      if (BOjl.gt.0.0_dp.and.BO(ni,i)*BO(nk,i)*BO(nl,j).gt.reaxFFatol3) then
                        if (lStrict) then
                          fjl = 1.0_dp
                          dfjldBOjl = 0.0_dp
                          d2fjldBOjl2 = 0.0_dp                 
                        else
                          call ataper(.false.,BO(nl,j),reaxFFatol,reaxFFtaperscale*reaxFFatol,fjl,dfjldBOjl,d2fjldBOjl2, &
                                      d3fjldBOjl3,lgrad1,lgrad2,.false.)
                        endif
!
                        xkl = xneigh(ni,i) + xneigh(nl,j) - xneigh(nk,i)
                        ykl = yneigh(ni,i) + yneigh(nl,j) - yneigh(nk,i)
                        zkl = zneigh(ni,i) + zneigh(nl,j) - zneigh(nk,i)
                        rkl2 = xkl*xkl + ykl*ykl + zkl*zkl
!  The condition below modified from l.ne.i & l.ne.k to handle small cell systems
                        if (nbos(l).gt.0.and.nl.ne.neighnoRptr(ni,i).and.rkl2.gt.1.0d-8) then
                          nspecl = nbosptr(l)
!
!  We now have a quartet of valid reaxFF species
!
                          if (nspeck.ge.nspecl) then
                            ind2 = nspeck*(nspeck - 1)/2 + nspecl + 1
                          else
                            ind2 = nspecl*(nspecl - 1)/2 + nspeck + 1
                          endif
!
!  Assign local scalars
!
                          V1 = 0.5_dp*reaxFFtor4(1,ind2,ind1)
                          V2 = 0.5_dp*reaxFFtor4(2,ind2,ind1)
                          V3 = 0.5_dp*reaxFFtor4(3,ind2,ind1)
                          ptor1 = reaxFFtor4(4,ind2,ind1)
                          pcot1 = reaxFFtor4(5,ind2,ind1)
!
                          if (V1.eq.0.0_dp.and.V2.eq.0.0_dp.and.V3.eq.0.0_dp) then
!
!  No specific term so try wildcard
!
                            V1 = 0.5_dp*reaxFFtor4(1,1,ind1)
                            V2 = 0.5_dp*reaxFFtor4(2,1,ind1)
                            V3 = 0.5_dp*reaxFFtor4(3,1,ind1)
                            ptor1 = reaxFFtor4(4,1,ind1)
                            pcot1 = reaxFFtor4(5,1,ind1)
                          endif
!
                          if (.not.lreaxFF_etors) then
                            V1 = 0.0_dp
                            V2 = 0.0_dp
                            V3 = 0.0_dp
                          endif
                          if (.not.lreaxFF_econj) then
                            pcot1 = 0.0_dp
                          endif
!
!  Check whether the potentials are non-zero to avoid wasted calculation
!
                          if (abs(V1)+abs(V2)+abs(V3)+abs(pcot1).gt.1.0d-12) then
!
!  Compute unknown distances, j-k and i-l
!
                            xjk =   xneigh(ni,i) - xneigh(nk,i)
                            yjk =   yneigh(ni,i) - yneigh(nk,i)
                            zjk =   zneigh(ni,i) - zneigh(nk,i)
                            rjk2 = xjk**2 + yjk**2 + zjk**2
                            xil = - xneigh(ni,i) - xneigh(nl,j)
                            yil = - yneigh(ni,i) - yneigh(nl,j)
                            zil = - zneigh(ni,i) - zneigh(nl,j)
                            ril2 = xil**2 + yil**2 + zil**2
                            xkl = xil + xneigh(nk,i) 
                            ykl = yil + yneigh(nk,i) 
                            zkl = zil + zneigh(nk,i) 
                            rkl2 = xkl**2 + ykl**2 + zkl**2
!
!  Calculate angles k-i-j and i-j-l
!
                            rjk = sqrt(rjk2)
                            ril = sqrt(ril2)
                            rkl = sqrt(rkl2)
                            call reaxff_sintheta(rneigh(ni,i),rneigh(nk,i),rjk,sinkij,dsinkijdr(1),dsinkijdr(2),dsinkijdr(3), &
                                                 d2sinkijdr2,lgrad1,lgrad2)
                            call reaxff_sintheta(rneigh(ni,i),rneigh(nl,j),ril,sinijl,dsinijldr(1),dsinijldr(2),dsinijldr(3), &
                                                 d2sinijldr2,lgrad1,lgrad2)
!
!  Compute torsion angle
!
                            call reaxff_cosphi(rneigh(nk,i),rjk,rkl,rneigh(ni,i),ril,rneigh(nl,j),cosphi,cosp1d,cosp2d, &
                                               lgrad1,lgrad2)
!
!  Compute multiple angle terms
!
                            cos2phi = 2.0_dp*cosphi*cosphi - 1.0_dp
                            cos3phi = cosphi*(4.0_dp*cosphi*cosphi - 3.0_dp)
!
!  Compute functions f10, f11 and f12
!
                            call reaxFF_f10(BOik,BOij,BOjl,f10,df10dBOik,df10dBOij,df10dBOjl, &
                                            d2f10dBOik2,d2f10dBOikdBOij,d2f10dBOikdBOjl,d2f10dBOij2, &
                                            d2f10dBOijdBOjl,d2f10dBOjl2,lgrad1,lgrad2)
                            call reaxFF_f11(nspeci,nspecj,delta(i),delta(j),f11,df11ddeltaij,d2f11ddeltaij2,lgrad1,lgrad2)
!
                            if (lgrad1) then
                              df11ddeltai = df11ddeltaij
                              df11ddeltaj = df11ddeltaij
                              if (lgrad2) then
                                d2f11ddeltai2 = d2f11ddeltaij2
                                d2f11ddeltaij = d2f11ddeltaij2
                                d2f11ddeltaj2 = d2f11ddeltaij2
                              endif
                            endif
!
                            call reaxFF_f12(BOik,BOij,BOjl,f12,df12dBOik,df12dBOij,df12dBOjl, &
                                            d2f12dBOik2,d2f12dBOikdBOij,d2f12dBOikdBOjl,d2f12dBOij2, &
                                            d2f12dBOijdBOjl,d2f12dBOjl2,lStrict,lgrad1,lgrad2)
!
!  Compute phi related terms
!
                            ctwo = 2.0_dp
                            expV2 = exp(ptor1*(ctwo - BO_pi(ni,i) - f11)**2)
                            V1trm = V1*(1.0_dp + cosphi)
                            V2trm = V2*expV2*(1.0_dp - cos2phi)
                            V3trm = V3*(1.0_dp + cos3phi)
                            conjtrm1 = (cosphi**2 - 1.0_dp)
                            conjtrm2 = (1.0_dp + conjtrm1*sinkij*sinijl)
                            if (lilocal) then
!
!  Compute torsional and conjugation energies
!
                              etors = etors + fij*fik*fjl*f10*sinkij*sinijl*(V1trm + V2trm + V3trm)
                              econj = econj + fij*fik*fjl*f12*pcot1*conjtrm2
!
!  Site energy contributions for etors and econj
!
                              esite = 0.25_dp*fij*fik*fjl*(f10*sinkij*sinijl*(V1trm + V2trm + V3trm) + &
                                                           f12*pcot1*conjtrm2)
                              siteenergy(i) = siteenergy(i) + esite
                              siteenergy(j) = siteenergy(j) + esite
                              siteenergy(k) = siteenergy(k) + esite
                              siteenergy(l) = siteenergy(l) + esite
                            endif
!
                            if (lgrad1) then
                              indbojl = ni*maxneigh + nl   ! 2nd shell index for l in shell of j
!
!  Set indices for cosphi
!
                              if (ni.ge.nk) then
                                indjk = nneigh(i) + ni*(ni - 1)/2 + nk
                              else
                                indjk = nneigh(i) + nk*(nk - 1)/2 + ni
                              endif
                              indkl = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + nk*mneigh + nl
                              indil = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + nl
                              indjl = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + ni*mneigh + nl
!
                              indcosphi(1) = nk
                              indcosphi(2) = indjk
                              indcosphi(3) = indkl
                              indcosphi(4) = ni
                              indcosphi(5) = indil
                              indcosphi(6) = indjl
!
!  Set indices for sines
!
                              indsinkij(1) = ni
                              indsinkij(2) = nk
                              indsinkij(3) = indjk
                              indsinijl(1) = ni
                              indsinijl(2) = indjl
                              indsinijl(3) = indil
!                                                
!  Add indices to sparsity pattern
!                     
                              if (nonzerorptr(ni).eq.0) then                                 
                                nonzero = nonzero + 1
                                nonzeroptr(nonzero) = ni
                                nonzerorptr(ni) = nonzero
                              endif       
                              if (nonzerorptr(nk).eq.0) then
                                nonzero = nonzero + 1
                                nonzeroptr(nonzero) = nk
                                nonzerorptr(nk) = nonzero
                              endif       
                              if (nonzerorptr(indil).eq.0) then 
                                nonzero = nonzero + 1
                                nonzeroptr(nonzero) = indil
                                nonzerorptr(indil) = nonzero       
                              endif
                              if (nonzerorptr(indjk).eq.0) then 
                                nonzero = nonzero + 1
                                nonzeroptr(nonzero) = indjk
                                nonzerorptr(indjk) = nonzero       
                              endif
                              if (nonzerorptr(indjl).eq.0) then 
                                nonzero = nonzero + 1
                                nonzeroptr(nonzero) = indjl
                                nonzerorptr(indjl) = nonzero       
                              endif
                              if (nonzerorptr(indkl).eq.0) then 
                                nonzero = nonzero + 1
                                nonzeroptr(nonzero) = indkl
                                nonzerorptr(indkl) = nonzero       
                              endif
!
!  Check there is room in the sparse arrays
!
                              if (nonzero.gt.maxsparse) then
                                oldmaxsparse = maxsparse
                                oldmaxsparse2 = maxsparse*(maxsparse+1)/2
                                maxsparse = nonzero + 500 
                                call realloc(d1is,maxsparse,status)
                                if (status.ne.0) call outofmemory('reaxFFds','d1is')
                                maxsparse2 = maxsparse*(maxsparse+1)/2
                                call realloc(d2is,maxsparse2,status)
                                if (status.ne.0) call outofmemory('reaxFFds','d2is')
                                call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
                                if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
                                call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
                                if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
                                call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
                                if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
!
                                d1is(oldmaxsparse+1:maxsparse) = 0.0_dp
                                d2is(oldmaxsparse2+1:maxsparse2) = 0.0_dp
                                d2evaldthetadBO_neigh(oldmaxsparse+1:maxsparse,1:maxneigh1a) = 0.0_dp
                                d2etorsdthetaddeltaj_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                                d2etorsdthetadBOpi_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                              endif 
!
!  Set indices for sines and cosines after mapping to sparse form
!
                              do m = 1,6
                                indcosphis(m) = nonzerorptr(indcosphi(m))
                              enddo
!
                              do m = 1,3
                                indsinkijs(m) = nonzerorptr(indsinkij(m))
                                indsinijls(m) = nonzerorptr(indsinijl(m))
                              enddo
!
!  Add on contributions to bond order derivatives for i-j & i-k
!
                              detorsdf10 = fij*fik*fjl*sinkij*sinijl*(V1trm + V2trm + V3trm)
                              devaldBO_neigh(ni) = devaldBO_neigh(ni) + detorsdf10*df10dBOij
                              devaldBO_neigh(nk) = devaldBO_neigh(nk) + detorsdf10*df10dBOik
                              devaldBO_neigh(indbojl) = devaldBO_neigh(indbojl) + detorsdf10*df10dBOjl
!
                              deconjdf12 = fij*fik*fjl*pcot1*conjtrm2
                              devaldBO_neigh(ni) = devaldBO_neigh(ni) + deconjdf12*df12dBOij
                              devaldBO_neigh(nk) = devaldBO_neigh(nk) + deconjdf12*df12dBOik
                              devaldBO_neigh(indbojl) = devaldBO_neigh(indbojl) + deconjdf12*df12dBOjl
!
!  Derivatives of taper functions with respect to bond orders
!
                              detorsconjdf = f10*sinkij*sinijl*(V1trm + V2trm + V3trm) + f12*pcot1*conjtrm2
                              devaldBO_neigh(ni) = devaldBO_neigh(ni) + detorsconjdf*fik*fjl*dfijdBOij
                              devaldBO_neigh(nk) = devaldBO_neigh(nk) + detorsconjdf*fij*fjl*dfikdBOik
                              devaldBO_neigh(indbojl) = devaldBO_neigh(indbojl) + detorsconjdf*fij*fik*dfjldBOjl
!
!  Add on contributions to pi bond order derivatives for i-j from expV2
!
                              detorsdexpV2 = fij*fik*fjl*f10*sinkij*sinijl*V2*(1.0_dp - cos2phi)
                              dexpV2dBO_pi = - 2.0_dp*ptor1*expV2*(ctwo - BO_pi(ni,i) - f11)
                              detorsdBOpi_neigh(ni) = detorsdBOpi_neigh(ni) + detorsdexpV2*dexpV2dBO_pi
!
!  Delta i and j derivatives of f11
!
                              dexpV2df11 = - 2.0_dp*ptor1*expV2*(ctwo - BO_pi(ni,i) - f11)
                              devalddeltai = devalddeltai + detorsdexpV2*dexpV2df11*df11ddeltai
                              devalddeltaj_neigh(ni) = devalddeltaj_neigh(ni) + detorsdexpV2*dexpV2df11*df11ddeltaj
!
!  Derivatives of torsional angles
!
                              dcos2phidcosphi = 4.0_dp*cosphi
                              dcos3phidcosphi = 12.0_dp*cosphi*cosphi - 3.0_dp
!
!  Derivatives of energy with respect to multiple angles
!
                              detorsdcos1phi = f10*sinkij*sinijl*V1
                              detorsdcos2phi = - f10*sinkij*sinijl*V2*expV2
                              detorsdcos3phi = f10*sinkij*sinijl*V3
!
!  Finish torsional angle derivatives + taper derivatives where needed
!
                              detorsdcosphi = (detorsdcos1phi + detorsdcos2phi*dcos2phidcosphi + detorsdcos3phi*dcos3phidcosphi)
                              deconjdcosphi = 2.0_dp*f12*pcot1*cosphi*sinkij*sinijl
!
                              detorsdcosphi = detorsdcosphi + deconjdcosphi
                              detorsdcosphi = fij*fik*fjl*detorsdcosphi
!
                              do m = 1,6
                                d1is(indcosphis(m)) = d1is(indcosphis(m)) + detorsdcosphi*cosp1d(m)
                              enddo
!
!  Derivatives with respect to sine of angles
!
                              detorsdsinkij = f10*sinijl*(V1trm + V2trm + V3trm)
                              detorsdsinijl = f10*sinkij*(V1trm + V2trm + V3trm)
                              deconjdsinkij = f12*pcot1*conjtrm1*sinijl
                              deconjdsinijl = f12*pcot1*conjtrm1*sinkij
                              detorsdsinkij = detorsdsinkij + deconjdsinkij
                              detorsdsinijl = detorsdsinijl + deconjdsinijl
!
                              detorsdsinkij = fij*fik*fjl*detorsdsinkij
                              detorsdsinijl = fij*fik*fjl*detorsdsinijl
!
                              do m = 1,3
                                d1is(indsinkijs(m)) = d1is(indsinkijs(m)) + detorsdsinkij*dsinkijdr(m)
                                d1is(indsinijls(m)) = d1is(indsinijls(m)) + detorsdsinijl*dsinijldr(m)
                              enddo
!
                              if (lgrad2) then
!
!  Second derivatives
!
                                d2etorsdf2 = sinkij*sinijl*(V1trm + V2trm + V3trm)
                                d2econjdf2 = pcot1*conjtrm2
!
                                ind = ni*(ni+1)/2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fik*fjl*d2fijdBOij2*detorsconjdf
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*fjl*d2f10dBOij2*d2etorsdf2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*fjl*d2f12dBOij2*d2econjdf2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fik*fjl*dfijdBOij*df10dBOij*d2etorsdf2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fik*fjl*dfijdBOij*df12dBOij*d2econjdf2
                                ind = nk*(nk+1)/2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fjl*d2fikdBOik2*detorsconjdf
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*fjl*d2f10dBOik2*d2etorsdf2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*fjl*d2f12dBOik2*d2econjdf2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fij*fjl*dfikdBOik*df10dBOik*d2etorsdf2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fij*fjl*dfikdBOik*df12dBOik*d2econjdf2
                                ind = indbojl*(indbojl+1)/2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*d2fjldBOjl2*detorsconjdf
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*fjl*d2f10dBOjl2*d2etorsdf2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*fjl*d2f12dBOjl2*d2econjdf2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fij*fik*dfjldBOjl*df10dBOjl*d2etorsdf2
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + 2.0_dp*fij*fik*dfjldBOjl*df12dBOjl*d2econjdf2
!
!  Second derivatives of expV2 w.r.t. BOpi
!
                                d2expV2dBO_pi2 = ptor1*expV2*(4.0_dp*ptor1*(ctwo - BO_pi(ni,i) - f11)**2 + 2.0_dp)
                                d2expV2df11dBO_pi = ptor1*expV2*(4.0_dp*ptor1*(ctwo - BO_pi(ni,i) - f11)**2 + 2.0_dp)
                                d2expV2df112 = ptor1*expV2*(4.0_dp*ptor1*(ctwo - BO_pi(ni,i) - f11)**2 + 2.0_dp)
                                d2etorsdBOpi2_neigh(ni) = d2etorsdBOpi2_neigh(ni) + detorsdexpV2*d2expV2dBO_pi2
!
!  Mixed expV2 and cosphi second derivatives
!
                                do m = 1,6
                                  d2etorsdthetadBOpi_neigh(indcosphis(m),ni) = d2etorsdthetadBOpi_neigh(indcosphis(m),ni) - &
                                    fij*fik*fjl*f10*sinkij*sinijl*V2*dexpV2dBO_pi*dcos2phidcosphi*cosp1d(m)
                                  d2etorsdthetaddeltai(indcosphi(m)) = d2etorsdthetaddeltai(indcosphi(m)) - &
                                    fij*fik*fjl*f10*sinkij*sinijl*V2*dexpV2df11*df11ddeltai*dcos2phidcosphi*cosp1d(m)
                                  d2etorsdthetaddeltaj_neigh(indcosphis(m),ni) = d2etorsdthetaddeltaj_neigh(indcosphis(m),ni) - &
                                    fij*fik*fjl*f10*sinkij*sinijl*V2*dexpV2df11*df11ddeltaj*dcos2phidcosphi*cosp1d(m)
                                enddo
!
!  Mixed expV2 and sin second derivatives
!
                                do m = 1,3
                                  d2etorsdthetadBOpi_neigh(indsinkijs(m),ni) = d2etorsdthetadBOpi_neigh(indsinkijs(m),ni) + &
                                    fij*fik*fjl*f10*sinijl*V2*dexpV2dBO_pi*(1.0_dp - cos2phi)*dsinkijdr(m)
                                  d2etorsdthetadBOpi_neigh(indsinijls(m),ni) = d2etorsdthetadBOpi_neigh(indsinijls(m),ni) + &
                                    fij*fik*fjl*f10*sinkij*V2*dexpV2dBO_pi*(1.0_dp - cos2phi)*dsinijldr(m)
!
                                  d2etorsdthetaddeltai(indsinkij(m)) = d2etorsdthetaddeltai(indsinkij(m)) + &
                                    fij*fik*fjl*f10*sinijl*V2*dexpV2df11*df11ddeltai*(1.0_dp - cos2phi)*dsinkijdr(m)
                                  d2etorsdthetaddeltai(indsinijl(m)) = d2etorsdthetaddeltai(indsinijl(m)) + &
                                    fij*fik*fjl*f10*sinkij*V2*dexpV2df11*df11ddeltai*(1.0_dp - cos2phi)*dsinijldr(m)
!
                                  d2etorsdthetaddeltaj_neigh(indsinkijs(m),ni) = d2etorsdthetaddeltaj_neigh(indsinkijs(m),ni) + &
                                    fij*fik*fjl*f10*sinijl*V2*dexpV2df11*df11ddeltaj*(1.0_dp - cos2phi)*dsinkijdr(m)
                                  d2etorsdthetaddeltaj_neigh(indsinijls(m),ni) = d2etorsdthetaddeltaj_neigh(indsinijls(m),ni) + &
                                    fij*fik*fjl*f10*sinkij*V2*dexpV2df11*df11ddeltaj*(1.0_dp - cos2phi)*dsinijldr(m)
                                enddo
!
!  Mixed expV2 and BO second derivatives
!
                                d2etorsdBOdBOpi_neigh(ni,ni) = d2etorsdBOdBOpi_neigh(ni,ni) + &
                                  fik*fjl*f10*sinkij*sinijl*V2*dexpV2dBO_pi*(1.0_dp - cos2phi)*dfijdBOij
                                d2etorsdBOdBOpi_neigh(nk,ni) = d2etorsdBOdBOpi_neigh(nk,ni) + &
                                  fij*fjl*f10*sinkij*sinijl*V2*dexpV2dBO_pi*(1.0_dp - cos2phi)*dfikdBOik
                                d2etorsdBOdBOpi_neigh(indbojl,ni) = d2etorsdBOdBOpi_neigh(indbojl,ni) + &
                                  fij*fik*f10*sinkij*sinijl*V2*dexpV2dBO_pi*(1.0_dp - cos2phi)*dfjldBOjl
!
                                d2etorsdBOdBOpi_neigh(ni,ni) = d2etorsdBOdBOpi_neigh(ni,ni) + &
                                  fij*fik*fjl*sinkij*sinijl*V2*dexpV2dBO_pi*(1.0_dp - cos2phi)*df10dBOij
                                d2etorsdBOdBOpi_neigh(nk,ni) = d2etorsdBOdBOpi_neigh(nk,ni) + &
                                  fij*fik*fjl*sinkij*sinijl*V2*dexpV2dBO_pi*(1.0_dp - cos2phi)*df10dBOik
                                d2etorsdBOdBOpi_neigh(indbojl,ni) = d2etorsdBOdBOpi_neigh(indbojl,ni) + &
                                  fij*fik*fjl*sinkij*sinijl*V2*dexpV2dBO_pi*(1.0_dp - cos2phi)*df10dBOjl
!
!  Second derivatives of expV2 w.r.t. f11
!
                                d2evalddeltai2 = d2evalddeltai2 + detorsdexpV2*( &
                                  d2expV2df112*df11ddeltai**2 + dexpV2df11*d2f11ddeltai2)
                                d2evalddeltaiddeltaj_neigh(ni) = d2evalddeltaiddeltaj_neigh(ni) + detorsdexpV2*( &
                                  d2expV2df112*df11ddeltai*df11ddeltaj + dexpV2df11*d2f11ddeltaij)
                                ind = ni*(ni+1)/2
                                d2evalddeltajk_neigh(ind) = d2evalddeltajk_neigh(ind) + detorsdexpV2*( &
                                  d2expV2df112*df11ddeltaj**2 + dexpV2df11*d2f11ddeltaj2)
!
!  Mixed second derivatives of expV2 w.r.t. f11 and BOpi
!
                                d2etorsddeltaidBOpi_neigh(ni) = d2etorsddeltaidBOpi_neigh(ni) + &
                                  detorsdexpV2*d2expV2df11dBO_pi*df11ddeltai
                                d2etorsddeltajdBOpi_neigh(ni) = d2etorsddeltajdBOpi_neigh(ni) + &
                                  detorsdexpV2*d2expV2df11dBO_pi*df11ddeltaj
!
!  Set up torsion terms for second derivatives
!
                                d2cos2phidcosphi2 = 4.0_dp
                                d2cos3phidcosphi2 = 24.0_dp*cosphi
                                d2etorsdcosphi2 = (detorsdcos2phi*d2cos2phidcosphi2 + detorsdcos3phi*d2cos3phidcosphi2) + &
                                                  2.0_dp*f12*pcot1*sinkij*sinijl
                                d2etorsdcosphi2 = fij*fik*fjl*d2etorsdcosphi2
!
!  Mixed second derivatives of expV2 w.r.t. f11 and tapers
!
                                detorsdexpV2notaper = f10*sinkij*sinijl*V2*(1.0_dp - cos2phi)
                                d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                                  fik*fjl*detorsdexpV2notaper*dexpV2df11*df11ddeltai*dfijdBOij
                                d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                                  fij*fjl*detorsdexpV2notaper*dexpV2df11*df11ddeltai*dfikdBOik
                                d2evalddeltaidBO_neigh(indbojl) = d2evalddeltaidBO_neigh(indbojl) + &
                                  fij*fik*detorsdexpV2notaper*dexpV2df11*df11ddeltai*dfjldBOjl
!
                                d2evalddeltajdBO_neigh(ni,ni) = d2evalddeltajdBO_neigh(ni,ni) + &
                                  fik*fjl*detorsdexpV2notaper*dexpV2df11*df11ddeltaj*dfijdBOij
                                d2evalddeltajdBO_neigh(ni,nk) = d2evalddeltajdBO_neigh(ni,nk) + &
                                  fij*fjl*detorsdexpV2notaper*dexpV2df11*df11ddeltaj*dfikdBOik
                                d2evalddeltajdBO_neigh(ni,indbojl) = d2evalddeltajdBO_neigh(ni,indbojl) + &
                                  fij*fik*detorsdexpV2notaper*dexpV2df11*df11ddeltaj*dfjldBOjl
!
!  Mixed second derivatives of expV2 w.r.t. f11 and f10
!
                                detorsdexpV2nof10 = fij*fik*fjl*sinkij*sinijl*V2*(1.0_dp - cos2phi)
                                d2evalddeltaidBO_neigh(ni) = d2evalddeltaidBO_neigh(ni) + &
                                  detorsdexpV2nof10*dexpV2df11*df11ddeltai*df10dBOij
                                d2evalddeltaidBO_neigh(nk) = d2evalddeltaidBO_neigh(nk) + &
                                  detorsdexpV2nof10*dexpV2df11*df11ddeltai*df10dBOik
                                d2evalddeltaidBO_neigh(indbojl) = d2evalddeltaidBO_neigh(indbojl) + &
                                  detorsdexpV2nof10*dexpV2df11*df11ddeltai*df10dBOjl
!
                                d2evalddeltajdBO_neigh(ni,ni) = d2evalddeltajdBO_neigh(ni,ni) + &
                                  detorsdexpV2nof10*dexpV2df11*df11ddeltaj*df10dBOij
                                d2evalddeltajdBO_neigh(ni,nk) = d2evalddeltajdBO_neigh(ni,nk) + &
                                  detorsdexpV2nof10*dexpV2df11*df11ddeltaj*df10dBOik
                                d2evalddeltajdBO_neigh(ni,indbojl) = d2evalddeltajdBO_neigh(ni,indbojl) + &
                                  detorsdexpV2nof10*dexpV2df11*df11ddeltaj*df10dBOjl
!
!  Torsion second derivatives
!
!  1 = nk
!  2 = indjk
!  3 = indkl
!  4 = ni
!  5 = indil
!  6 = indjl
!
                                ind2 = 0
                                do m = 1,6
                                  do n = m,6
                                    ind2 = ind2 + 1
                                    if (indcosphis(m).ge.indcosphis(n)) then
                                      ind = indcosphis(m)*(indcosphis(m)-1)/2 + indcosphis(n)
                                    else
                                      ind = indcosphis(n)*(indcosphis(n)-1)/2 + indcosphis(m)
                                    endif
                                    d2is(ind) = d2is(ind) + detorsdcosphi*cosp2d(ind2) + d2etorsdcosphi2*cosp1d(m)*cosp1d(n)
                                  enddo
                                enddo
!
!  Mixed torsion angle - sine derivatives
!
                                d2etorsdcosphidsinkij = fij*fik*fjl*sinijl*(f10*(V1 - V2*expV2*dcos2phidcosphi + &
                                                        V3*dcos3phidcosphi) + 2.0_dp*f12*pcot1*cosphi)
                                d2etorsdcosphidsinijl = fij*fik*fjl*sinkij*(f10*(V1 - V2*expV2*dcos2phidcosphi + &
                                                        V3*dcos3phidcosphi) + 2.0_dp*f12*pcot1*cosphi)
                                do m = 1,6
                                  do n = 1,3
                                    if (indcosphis(m).eq.indsinkijs(n)) then
                                      ind = indcosphis(m)*(indcosphis(m)+1)/2
                                      d2is(ind) = d2is(ind) + 2.0_dp*d2etorsdcosphidsinkij*dsinkijdr(n)*cosp1d(m)
                                    elseif (indcosphis(m).gt.indsinkijs(n)) then
                                      ind = indcosphis(m)*(indcosphis(m)-1)/2 + indsinkijs(n)
                                      d2is(ind) = d2is(ind) + d2etorsdcosphidsinkij*dsinkijdr(n)*cosp1d(m)
                                    else
                                      ind = indsinkijs(n)*(indsinkijs(n)-1)/2 + indcosphis(m)
                                      d2is(ind) = d2is(ind) + d2etorsdcosphidsinkij*dsinkijdr(n)*cosp1d(m)
                                    endif
                                    if (indcosphis(m).eq.indsinijls(n)) then
                                      ind = indcosphis(m)*(indcosphis(m)+1)/2
                                      d2is(ind) = d2is(ind) + 2.0_dp*d2etorsdcosphidsinijl*dsinijldr(n)*cosp1d(m)
                                    elseif (indcosphis(m).gt.indsinijls(n)) then
                                      ind = indcosphis(m)*(indcosphis(m)-1)/2 + indsinijls(n)
                                      d2is(ind) = d2is(ind) + d2etorsdcosphidsinijl*dsinijldr(n)*cosp1d(m)
                                    else
                                      ind = indsinijls(n)*(indsinijls(n)-1)/2 + indcosphis(m)
                                      d2is(ind) = d2is(ind) + d2etorsdcosphidsinijl*dsinijldr(n)*cosp1d(m)
                                    endif
                                  enddo
                                enddo
!
!  Mixed bond order - bond order second derivatives
!
                                if (ni.ge.nk) then
                                  ind = ni*(ni - 1)/2 + nk
                                else
                                  ind = nk*(nk - 1)/2 + ni
                                endif
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fjl*dfijdBOij*dfikdBOik*detorsconjdf
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fik*fjl*dfijdBOij*(df10dBOik*d2etorsdf2 + &
                                                                                df12dBOik*d2econjdf2)
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fjl*dfikdBOik*(df10dBOij*d2etorsdf2 + &
                                                                                df12dBOij*d2econjdf2)
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*fjl*(d2f10dBOikdBOij*d2etorsdf2 + &
                                                                                d2f12dBOikdBOij*d2econjdf2)
                                ind = indbojl*(indbojl - 1)/2 + ni
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fik*dfjldBOjl*dfijdBOij*detorsconjdf
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*dfjldBOjl*(df10dBOij*d2etorsdf2 + &
                                                                                df12dBOij*d2econjdf2)
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fik*fjl*dfijdBOij*(df10dBOjl*d2etorsdf2 + &
                                                                                df12dBOjl*d2econjdf2)
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*fjl*(d2f10dBOijdBOjl*d2etorsdf2 + &
                                                                                d2f12dBOijdBOjl*d2econjdf2)
                                ind = indbojl*(indbojl - 1)/2 + nk
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*dfjldBOjl*dfikdBOik*detorsconjdf
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*dfjldBOjl*(df10dBOik*d2etorsdf2 + &
                                                                                df12dBOik*d2econjdf2)
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fjl*dfikdBOik*(df10dBOjl*d2etorsdf2 + &
                                                                                df12dBOjl*d2econjdf2)
                                d2evaldBO2_neigh(ind) = d2evaldBO2_neigh(ind) + fij*fik*fjl*(d2f10dBOikdBOjl*d2etorsdf2 + &
                                                                                d2f12dBOikdBOjl*d2econjdf2)
!
!  Second derivatives from sin terms
!
                                d2etorsdsinkijdsinijl = fij*fik*fjl*(f10*(V1trm + V2trm + V3trm) + f12*pcot1*conjtrm1)
!
                                ind2 = 0
                                do m = 1,3
                                  do n = m,3
                                    ind2 = ind2 + 1
                                    if (indsinkijs(m).ge.indsinkijs(n)) then
                                      ind = indsinkijs(m)*(indsinkijs(m)-1)/2 + indsinkijs(n)
                                    else
                                      ind = indsinkijs(n)*(indsinkijs(n)-1)/2 + indsinkijs(m)
                                    endif
                                    d2is(ind) = d2is(ind) + detorsdsinkij*d2sinkijdr2(ind2)
                                    if (indsinijls(m).ge.indsinijls(n)) then
                                      ind = indsinijls(m)*(indsinijls(m)-1)/2 + indsinijls(n)
                                    else
                                      ind = indsinijls(n)*(indsinijls(n)-1)/2 + indsinijls(m)
                                    endif
                                    d2is(ind) = d2is(ind) + detorsdsinijl*d2sinijldr2(ind2)
                                  enddo
                                enddo
!
                                do m = 1,3
                                  do n = 1,3
                                    if (indsinkijs(m).eq.indsinijls(n)) then
                                      ind = indsinkijs(m)*(indsinkijs(m)+1)/2
                                      d2is(ind) = d2is(ind) + 2.0_dp*d2etorsdsinkijdsinijl*dsinkijdr(m)*dsinijldr(n)
                                    else
                                      if (indsinkijs(m).gt.indsinijls(n)) then
                                        ind = indsinkijs(m)*(indsinkijs(m)-1)/2 + indsinijls(n)
                                      else
                                        ind = indsinijls(n)*(indsinijls(n)-1)/2 + indsinkijs(m)
                                      endif
                                      d2is(ind) = d2is(ind) + d2etorsdsinkijdsinijl*dsinkijdr(m)*dsinijldr(n)
                                    endif
                                  enddo
                                enddo
!
!  Mixed sin with fij, fik, fjl, f10 terms
!
                                d2etorsdsinkijdf = (f10*(V1trm + V2trm + V3trm) + f12*pcot1*conjtrm1)*sinijl
                                d2etorsdsinijldf = (f10*(V1trm + V2trm + V3trm) + f12*pcot1*conjtrm1)*sinkij
                                d2etorsdsinkijdf10 = fij*fik*fjl*(V1trm + V2trm + V3trm)*sinijl
                                d2etorsdsinijldf10 = fij*fik*fjl*(V1trm + V2trm + V3trm)*sinkij
                                d2econjdsinkijdf12 = fij*fik*fjl*pcot1*conjtrm1*sinijl
                                d2econjdsinijldf12 = fij*fik*fjl*pcot1*conjtrm1*sinkij
                                do n = 1,3
                                  d2evaldthetadBO_neigh(indsinkijs(n),ni) = d2evaldthetadBO_neigh(indsinkijs(n),ni) + &
                                    d2etorsdsinkijdf*dsinkijdr(n)*dfijdBOij*fik*fjl + &
                                    d2etorsdsinkijdf10*dsinkijdr(n)*df10dBOij + &
                                    d2econjdsinkijdf12*dsinkijdr(n)*df12dBOij
                                  d2evaldthetadBO_neigh(indsinkijs(n),nk) = d2evaldthetadBO_neigh(indsinkijs(n),nk) + &
                                    d2etorsdsinkijdf*dsinkijdr(n)*dfikdBOik*fij*fjl + &
                                    d2etorsdsinkijdf10*dsinkijdr(n)*df10dBOik + &
                                    d2econjdsinkijdf12*dsinkijdr(n)*df12dBOik
                                  d2evaldthetadBO_neigh(indsinkijs(n),indbojl) = d2evaldthetadBO_neigh(indsinkijs(n),indbojl) + &
                                    d2etorsdsinkijdf*dsinkijdr(n)*dfjldBOjl*fij*fik + &
                                    d2etorsdsinkijdf10*dsinkijdr(n)*df10dBOjl + &
                                    d2econjdsinkijdf12*dsinkijdr(n)*df12dBOjl
!
                                  d2evaldthetadBO_neigh(indsinijls(n),ni) = d2evaldthetadBO_neigh(indsinijls(n),ni) + &
                                    d2etorsdsinijldf*dsinijldr(n)*dfijdBOij*fik*fjl + &
                                    d2etorsdsinijldf10*dsinijldr(n)*df10dBOij + &
                                    d2econjdsinijldf12*dsinijldr(n)*df12dBOij
                                  d2evaldthetadBO_neigh(indsinijls(n),nk) = d2evaldthetadBO_neigh(indsinijls(n),nk) + &
                                    d2etorsdsinijldf*dsinijldr(n)*dfikdBOik*fij*fjl + &
                                    d2etorsdsinijldf10*dsinijldr(n)*df10dBOik + &
                                    d2econjdsinijldf12*dsinijldr(n)*df12dBOik
                                  d2evaldthetadBO_neigh(indsinijls(n),indbojl) = d2evaldthetadBO_neigh(indsinijls(n),indbojl) + &
                                    d2etorsdsinijldf*dsinijldr(n)*dfjldBOjl*fij*fik + &
                                    d2etorsdsinijldf10*dsinijldr(n)*df10dBOjl + &
                                    d2econjdsinijldf12*dsinijldr(n)*df12dBOjl
                                enddo
!
!  Mixed cosphi with fij, fik, fjl, f10 terms
!
                                d2etorsdcosphidf = sinkij*sinijl*(f10*(V1 - V2*expV2*dcos2phidcosphi + V3*dcos3phidcosphi) + &
                                                   2.0_dp*f12*pcot1*cosphi)
                                d2etorsdcosphidf10 = fij*fik*fjl*sinkij*sinijl*(V1 - V2*expV2*dcos2phidcosphi + V3*dcos3phidcosphi)
                                d2econjdcosphidf12 = 2.0_dp*fij*fik*fjl*sinkij*sinijl*pcot1*cosphi
                                do n = 1,6
                                  d2evaldthetadBO_neigh(indcosphis(n),ni) = d2evaldthetadBO_neigh(indcosphis(n),ni) + &
                                    d2etorsdcosphidf*cosp1d(n)*dfijdBOij*fik*fjl + &
                                    d2etorsdcosphidf10*cosp1d(n)*df10dBOij + &
                                    d2econjdcosphidf12*cosp1d(n)*df12dBOij
                                  d2evaldthetadBO_neigh(indcosphis(n),nk) = d2evaldthetadBO_neigh(indcosphis(n),nk) + &
                                    d2etorsdcosphidf*cosp1d(n)*dfikdBOik*fij*fjl + &
                                    d2etorsdcosphidf10*cosp1d(n)*df10dBOik + &
                                    d2econjdcosphidf12*cosp1d(n)*df12dBOik
                                  d2evaldthetadBO_neigh(indcosphis(n),indbojl) = d2evaldthetadBO_neigh(indcosphis(n),indbojl) + &
                                    d2etorsdcosphidf*cosp1d(n)*dfjldBOjl*fij*fik + &
                                    d2etorsdcosphidf10*cosp1d(n)*df10dBOjl + &
                                    d2econjdcosphidf12*cosp1d(n)*df12dBOjl
                                enddo
                              endif
                            endif
                          endif
                        endif
!
!  End of check on bond order for j-l
!
                      endif
                    enddo
!
!  End of check on bond order for i-k
!
                  endif
!
                endif
              enddo
!
!  End of check on bond order for i-j
!
            endif
!
          endif
        enddo tniloop
!
!  Now for the derivatives......
!
        if (lgrad1) then
!
!  Set up derivative of deltalp_corr
!
          if (.not.lreaxFFunder(nspeci)) then
            ddeltalpcorrddeltai = 1.0_dp
            if (lgrad2) then
              d2deltalpcorrddeltai2 = 0.0_dp
            endif
          else
            ddeltalpcorrddeltai = 1.0_dp - otrm1*ddeltalpddeltai
            if (lgrad2) then
              d2deltalpcorrddeltai2 = - otrm1*d2deltalpddeltai2
            endif
          endif
!
!  Set up term needed later for derivative of deltalpcorr_i with respect to sumdeltabopi
!
          if (.not.lreaxFFunder(nspeci)) then
            ddeltalpcorrdsumdeltabopi = 0.0_dp
            if (lgrad2) then
              d2deltalpcorrdsumdeltabopi2 = 0.0_dp
              d2deltalpcorrddeltaidsumdeltabopi = 0.0_dp
            endif
          else
            ddeltalpcorrdsumdeltabopi = (otrm1**2)*expdi*reaxFFlam(6)*reaxFFlam(31)*deltalp(i)
            if (lgrad2) then
              d2deltalpcorrdsumdeltabopi2 = (1.0_dp - 2.0_dp*otrm1*expdi*reaxFFlam(6))*ddeltalpcorrdsumdeltabopi* &
                                            reaxFFlam(31)
              d2deltalpcorrddeltaidsumdeltabopi = (otrm1**2)*expdi*reaxFFlam(6)*reaxFFlam(31)*ddeltalpddeltai
            endif
          endif
!---------------------------------------------------------------------------------
!  Set up terms for partial derivatives of eunder                                |
!---------------------------------------------------------------------------------
          dexpuc1ddeltalpcorr = reaxFFlam(7)*expuc1
          dutrm1ddeltalpcorr  = utrm1*utrm1*expuc2*reaxFFoc1(nspeci)
          dutrm2dsumdeltabopi = - utrm2*utrm2*expuc3*reaxFFlam(8)*reaxFFlam(9)
          if (lgrad2) then
            d2expuc1ddeltalpcorr2 = dexpuc1ddeltalpcorr*reaxFFlam(7)
            d2utrm1ddeltalpcorr2  = 2.0_dp*utrm1*utrm1*utrm1*(expuc2*reaxFFoc1(nspeci))**2 - &
                                    utrm1*utrm1*expuc2*reaxFFoc1(nspeci)**2
            d2utrm2dsumdeltabopi2 = - dutrm2dsumdeltabopi*reaxFFlam(9)*(2.0_dp*reaxFFlam(8)*expuc3*utrm2 - 1.0_dp)
          endif
!
!  Compute derivatives of eunder w.r.t. deltalpcorr and sumdeltabopi
!
          deunderddeltalpcorr = - reaxFFuc1(nspeci)*(-dexpuc1ddeltalpcorr*utrm1 + (1.0_dp - expuc1)*dutrm1ddeltalpcorr)*utrm2
          deunderdsumdeltabopi = - reaxFFuc1(nspeci)*(1.0_dp - expuc1)*utrm1*dutrm2dsumdeltabopi
          if (lgrad2) then
            d2eunderddeltalpcorr2 = - reaxFFuc1(nspeci)*(-d2expuc1ddeltalpcorr2*utrm1 + (1.0_dp - expuc1)*d2utrm1ddeltalpcorr2 &
                                                         - 2.0_dp*dexpuc1ddeltalpcorr*dutrm1ddeltalpcorr)*utrm2
            d2eunderddeltalpcorrdsumdeltabopi = - reaxFFuc1(nspeci)*(-dexpuc1ddeltalpcorr*utrm1 + &
                                                  (1.0_dp - expuc1)*dutrm1ddeltalpcorr)*dutrm2dsumdeltabopi
            d2eunderdsumdeltabopi2 = - reaxFFuc1(nspeci)*(1.0_dp - expuc1)*utrm1*d2utrm2dsumdeltabopi2
          endif
!---------------------------------------------------------------------------------
!  Set up terms for partial derivatives of eover                                 |
!---------------------------------------------------------------------------------
          deoverdsumover = otrm2*deltalpcorr_i*otrm3
          dotrm2ddeltalpcorr = - otrm2**2
          dotrm3ddeltalpcorr = - (otrm3**2)*reaxFFoc1(nspeci)*expoc
          deoverddeltalpcorr = sumover*(otrm2*otrm3 + deltalpcorr_i*(dotrm2ddeltalpcorr*otrm3 + otrm2*dotrm3ddeltalpcorr))
!
          if (lgrad2) then
            d2otrm2ddeltalpcorr2 = 2.0_dp*otrm2**3
            d2otrm3ddeltalpcorr2 = dotrm3ddeltalpcorr*reaxFFoc1(nspeci)*(1.0_dp - 2.0_dp*expoc*otrm3)
            d2eoverdsumoverddeltalpcorr = otrm2*otrm3 + deltalpcorr_i*(dotrm2ddeltalpcorr*otrm3 + otrm2*dotrm3ddeltalpcorr)
            d2eoverddeltalpcorr2 = sumover*(deltalpcorr_i*(d2otrm2ddeltalpcorr2*otrm3 + otrm2*d2otrm3ddeltalpcorr2 + &
                                            2.0_dp*dotrm2ddeltalpcorr*dotrm3ddeltalpcorr) + &
                                            2.0_dp*(dotrm2ddeltalpcorr*otrm3 + otrm2*dotrm3ddeltalpcorr))
          endif
!---------------------------------------------------------------------------------
!  Set up terms for partial derivatives of elone                                 |
!---------------------------------------------------------------------------------
          if (lreaxFF_elone.and.abs(reaxFFlp(2,nspeci)).gt.1.0d-12) then
            deloneddeltalplp = reaxFFlp(2,nspeci)*trm1lp*(trm1lp*deltalplp(i)*explp*reaxFFlam(30) + 1.0_dp)
            deloneddeltai = deloneddeltalplp*ddeltalplpddeltai
          else
            deloneddeltalplp = 0.0_dp
            deloneddeltai = 0.0_dp
          endif
          if (lgrad2) then
            if (lreaxFF_elone.and.abs(reaxFFlp(2,nspeci)).gt.1.0d-12) then
              d2eloneddeltalplp2 = reaxFFlp(2,nspeci)*trm1lp*trm1lp*reaxFFlam(30)*explp*( &
                                   2.0_dp - reaxFFlam(30)*deltalplp(i)*(1.0_dp - 2.0_dp*explp*trm1lp))
              d2eloneddeltai2 = d2eloneddeltalplp2*(ddeltalplpddeltai**2) + deloneddeltalplp*d2deltalplpddeltai2
            else
              d2eloneddeltalplp2 = 0.0_dp
              d2eloneddeltai2 = 0.0_dp
            endif
          endif
!---------------------------------------------------------------------------------
!  Set up terms for partial derivatives of ebpen                                 |
!---------------------------------------------------------------------------------
          debpenddeltai = 0.0_dp
          debpendBO_neigh(1:nneigh(i)) = 0.0_dp
          if (lgrad2) then
            d2ebpenddeltai2 = 0.0_dp
            d2ebpendBO2_neigh(1:nneigh(i)) = 0.0_dp
            d2ebpenddeltaidBO_neigh(1:nneigh(i)) = 0.0_dp
          endif
          if (lreaxFF_ebpen) then
            do ni = 1,nneigh(i)
              j = neighno(ni,i)
              if (nbos(j).gt.0) then
                nspecj = nbosptr(j)
                if (nspeci.ge.nspecj) then
                  nboij = nspeci*(nspeci - 1)/2 + nspecj
                else
                  nboij = nspecj*(nspecj - 1)/2 + nspeci
                endif
                diffBOdelta = BO(ni,i) - delta(i) - reaxFFpen2(2,nboij)*delta(i)**4 - reaxFFpen2(3,nboij)
                if (diffBOdelta.gt.0.0_dp) then
                  debpenddeltai = debpenddeltai - 2.0_dp*reaxFFpen2(1,nboij)*diffBOdelta* &
                    (1.0_dp + 4.0_dp*reaxFFpen2(2,nboij)*delta(i)**3)
                  debpendBO_neigh(ni) = 2.0_dp*reaxFFpen2(1,nboij)*diffBOdelta
                  if (lgrad2) then
                    d2ebpenddeltai2 = d2ebpenddeltai2 + 2.0_dp*reaxFFpen2(1,nboij)* &
                      ((1.0_dp + 4.0_dp*reaxFFpen2(2,nboij)*delta(i)**3)**2 - &
                       12.0_dp*diffBOdelta*reaxFFpen2(2,nboij)*delta(i)**2)
                    d2ebpendBO2_neigh(ni) = 2.0_dp*reaxFFpen2(1,nboij)
                    d2ebpenddeltaidBO_neigh(ni) = - 2.0_dp*reaxFFpen2(1,nboij)* &
                       (1.0_dp + 4.0_dp*reaxFFpen2(2,nboij)*delta(i)**3)
                  endif
                endif
              endif
            enddo
          endif
!---------------------------------------------------------------------------------
!  Sum up terms for partial derivatives of energy                                |
!---------------------------------------------------------------------------------
          detotddeltalpcorr = deunderddeltalpcorr + deoverddeltalpcorr
          if (lgrad2) then
            d2etotddeltalpcorr2 = d2eunderddeltalpcorr2 + d2eoverddeltalpcorr2
          endif
!*************************************
!  Loop over neighbours of i (=> j)  *
!*************************************
          do ni = 1,nneigh(i)
            j = neighno(ni,i)
            indj = nneigh(i) + nneigh(i)*(nneigh(i)+1)/2 + (ni-1)*mneigh*(mneigh + 1) + ni*mneigh
!
!  Do we need to do this pair of atoms
!
            if ((lopanyneigh(i).or.lopanyneigh(j)).and.nbos(j).gt.0) then
!
!  Set variables relating to j
!
              nspecj = nbosptr(j)
!
!  Find i in neighbour list for j
!
              nj = neighnoRptr(ni,i)
!
!  Set index to parameters for pairwise interaction of species
!
              if (nspeci.ge.nspecj) then
                nboij = nspeci*(nspeci - 1)/2 + nspecj
              else
                nboij = nspecj*(nspecj - 1)/2 + nspeci
              endif
!
!  Set up derivative of delta_lp for j with respect to delta_e
!
              delta_e = 0.5_dp*(delta(j) + reaxFFval(1,nspecj) - reaxFFval(3,nspecj))
!
              if (lreaxFFovsmooth) then
                delta_e0 = nint(delta_e)
                expH = exp(-reaxFFlpsmooth*(delta_e - delta_e0))
                smoothH = 1.0_dp/(1.0_dp + expH)
                rintde = delta_e0 - 1.0_dp + smoothH
                drintdeddeltae = reaxFFlpsmooth*expH*smoothH**2
                if (lgrad2) then
                  d2rintdeddeltae2 = reaxFFlpsmooth*reaxFFlpsmooth*expH*(smoothH**2)* &
                                     (2.0_dp*expH*smoothH - 1.0_dp)
                endif
              else
                rintde = dble(int(delta_e))
                drintdeddeltae = 0.0_dp
                d2rintdeddeltae2 = 0.0_dp
              endif
!
              trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
              explp2 = exp(-reaxFFlam(29)*trm2lp**2)
              ddeltalpddeltaj = 0.5_dp*(drintdeddeltae + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(1.0_dp - drintdeddeltae)) ! 1/2 is to convert from delta_e to delta_i
!
              if (lgrad2) then
                d2deltalpddeltaj2 = 0.25_dp*(d2rintdeddeltae2 + reaxFFlam(29)*explp2*( &
                                             8.0_dp*(1.0_dp - drintdeddeltae)**2 - &
                                             4.0_dp*trm2lp*d2rintdeddeltae2 - &
                                             16.0_dp*reaxFFlam(29)*trm2lp*trm2lp*(1.0_dp - drintdeddeltae)**2))
              endif
!
!  Set delta' 
!
              deltapi = deltap(i)
              deltapj = deltap(j)
!***************************************
!  Valid pairwise reaxFF contribution  *
!***************************************
              call reaxFF_bos(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                              d1BOp,d1BOp_pi,d1BOp_pipi,d2BOp,d2BOp_pi,d2BOp_pipi,BOij_s,BOij_pi,BOij_pipi, &
                              d1BOi_s,d1BOi_pi,d1BOi_pipi,d2BOi_s,d2BOi_pi,d2BOi_pipi, &
                              nonzero_ij,nonzero_ij_i,nonzeroptr_ij,lgrad1,lgrad2)
!
!  Convert non-zero indices
!
              if (lgrad1) then
                ind0 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + ni*mneigh
                do nl = nonzero_ij_i+1,nonzero_ij
                  nonzeroptr_ij(nl) = ind0 + nonzeroptr_ij(nl)
                enddo
              endif
!
!  Add indices to sparsity pattern if necessary
!
              do nl = 1,nonzero_ij
                l = nonzeroptr_ij(nl)
                if (nonzerorptr(l).eq.0) then
                  nonzero = nonzero + 1
                  nonzeroptr(nonzero) = l
                  nonzerorptr(l) = nonzero
                endif
              enddo
!  
!  Check there is room in the sparse arrays
!             
              if (nonzero.gt.maxsparse) then
                oldmaxsparse = maxsparse
                oldmaxsparse2 = maxsparse*(maxsparse+1)/2
                maxsparse = nonzero + 500
                call realloc(d1is,maxsparse,status)
                if (status.ne.0) call outofmemory('reaxFFds','d1is')
                maxsparse2 = maxsparse*(maxsparse+1)/2
                call realloc(d2is,maxsparse2,status)
                if (status.ne.0) call outofmemory('reaxFFds','d2is')
                call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
                if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
                call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
                if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
                call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
                if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
!
                d1is(oldmaxsparse+1:maxsparse) = 0.0_dp
                d2is(oldmaxsparse2+1:maxsparse2) = 0.0_dp
                d2evaldthetadBO_neigh(oldmaxsparse+1:maxsparse,1:maxneigh1a) = 0.0_dp
                d2etorsdthetaddeltaj_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                d2etorsdthetadBOpi_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
              endif 
!
!  Calculate total bond order derivatives from sum of components
!
              do nk = 1,nonzero_ij
                d1BOi(nk) = d1BOi_s(nk) + d1BOi_pi(nk) + d1BOi_pipi(nk)
              enddo
!
!  Derivatives of S w.r.t. BO pi and deltaj
!
              if (.not.lreaxFFunder(nspeci).or..not.lreaxFFunder(nspecj)) then
                dsumdeltabopidbopi = delta(j)
                dsumdeltabopiddeltaj = BO_pi(ni,i) + BO_pipi(ni,i)
              else
                dsumdeltabopidbopi = delta(j) - deltalp(j)
                dsumdeltabopiddeltaj = (1.0_dp - ddeltalpddeltaj)*(BO_pi(ni,i) + BO_pipi(ni,i))
              endif
              if (lgrad2) then
                if (.not.lreaxFFunder(nspeci).or..not.lreaxFFunder(nspecj)) then
                  d2sumdeltabopiddeltaj2 = 0.0_dp
                  d2sumdeltabopiddeltajdbopi = 1.0_dp
                else
                  d2sumdeltabopiddeltaj2 = - d2deltalpddeltaj2*(BO_pi(ni,i) + BO_pipi(ni,i))
                  d2sumdeltabopiddeltajdbopi = (1.0_dp - ddeltalpddeltaj)
                endif
              endif
!
!  Build derivatives of deltalpcorr w.r.t. BO pi and pi-pi and deltaj, plus mixed terms
!
              ddeltalpcorrdbopi = ddeltalpcorrdsumdeltabopi*dsumdeltabopidbopi
              ddeltalpcorrddeltaj = ddeltalpcorrdsumdeltabopi*dsumdeltabopiddeltaj
              if (lgrad2) then
                d2deltalpcorrdbopi2 = d2deltalpcorrdsumdeltabopi2*dsumdeltabopidbopi**2
                d2deltalpcorrddeltaidbopi = d2deltalpcorrddeltaidsumdeltabopi*dsumdeltabopidbopi
                d2deltalpcorrddeltajdbopi = d2deltalpcorrdsumdeltabopi2*dsumdeltabopiddeltaj*dsumdeltabopidbopi + &
                                            ddeltalpcorrdsumdeltabopi*d2sumdeltabopiddeltajdbopi
                d2deltalpcorrddeltaj2 = ddeltalpcorrdsumdeltabopi*d2sumdeltabopiddeltaj2   ! NB: d2deltalpcorrdsumdeltabopi2 handled via deltaj/deltaj2 term
                d2deltalpcorrddeltaiddeltaj = d2deltalpcorrddeltaidsumdeltabopi*dsumdeltabopiddeltaj 
              endif
!
              dsumoverdbo = reaxFFoc2(nboij)*reaxFFDe(1,nboij)
!
!  Compute energy derivatives
!
              detotdbopi = deunderdsumdeltabopi*dsumdeltabopidbopi                     ! Partial derivative from sumdeltabopi
              detotddeltaj = deunderdsumdeltabopi*dsumdeltabopiddeltaj                 ! Partial derivative from sumdeltabopi
!
              detotdbo_tot = detotddeltalpcorr*ddeltalpcorrddeltai + deoverdsumover*dsumoverdbo + &
                             deloneddeltai + debpenddeltai + debpendBO_neigh(ni)
              detotdbopi_tot = detotddeltalpcorr*ddeltalpcorrdbopi + detotdbopi
              detotddeltaj_tot = detotddeltalpcorr*ddeltalpcorrddeltaj + detotddeltaj
!
              if (lgrad2) then
                d2etotdbopi2 = d2eunderdsumdeltabopi2*dsumdeltabopidbopi**2
                d2etotddeltaj2 = deunderdsumdeltabopi*d2sumdeltabopiddeltaj2   ! NB: d2eunderdsumdeltabopi2 handled via deltaj/deltaj2 term
!
                d2etotdbo2_tot = d2etotddeltalpcorr2*ddeltalpcorrddeltai**2 + &
                                 detotddeltalpcorr*d2deltalpcorrddeltai2 + &
                                 2.0_dp*d2eoverdsumoverddeltalpcorr*dsumoverdbo*ddeltalpcorrddeltai + &
                                 d2eloneddeltai2 + d2ebpenddeltai2 + d2ebpendBO2_neigh(ni) + &
                                 2.0_dp*d2ebpenddeltaidBO_neigh(ni)
                d2etotdbopi2_tot = d2etotddeltalpcorr2*ddeltalpcorrdbopi**2 + &
                                   detotddeltalpcorr*d2deltalpcorrdbopi2 + &
                                   2.0_dp*d2eunderddeltalpcorrdsumdeltabopi*ddeltalpcorrdbopi*dsumdeltabopidbopi + &
                                   d2etotdbopi2
!
                d2etotdbodbopi_tot = d2etotddeltalpcorr2*ddeltalpcorrddeltai*ddeltalpcorrdbopi + &
                                     detotddeltalpcorr*d2deltalpcorrddeltaidbopi + &
                                     d2eunderddeltalpcorrdsumdeltabopi*ddeltalpcorrddeltai*dsumdeltabopidbopi + &
                                     d2eoverdsumoverddeltalpcorr*dsumoverdbo*ddeltalpcorrdbopi
                d2etotddeltaj2_tot = detotddeltalpcorr*d2deltalpcorrddeltaj2 + d2etotddeltaj2
!
                d2etotddeltajdbo_tot = d2etotddeltalpcorr2*ddeltalpcorrddeltai*ddeltalpcorrddeltaj + &
                                       detotddeltalpcorr*d2deltalpcorrddeltaiddeltaj + &
                                       d2eunderddeltalpcorrdsumdeltabopi*dsumdeltabopiddeltaj*ddeltalpcorrddeltai
                d2etotdbodsumdeltabopi = d2etotddeltalpcorr2*ddeltalpcorrddeltai*ddeltalpcorrdsumdeltabopi + &
                                         detotddeltalpcorr*d2deltalpcorrddeltaidsumdeltabopi + &
                                         d2eunderddeltalpcorrdsumdeltabopi*ddeltalpcorrddeltai
!
                d2etotdsumdeltabopi2_tot = d2eunderdsumdeltabopi2 + d2etotddeltalpcorr2*ddeltalpcorrdsumdeltabopi**2 + &
                                           2.0_dp*d2eunderddeltalpcorrdsumdeltabopi*ddeltalpcorrdsumdeltabopi + &
                                           detotddeltalpcorr*d2deltalpcorrdsumdeltabopi2
              endif
!
              do nk = 1,nonzero_ij
                k = nonzeroptr_ij(nk)
                nks = nonzerorptr(k)
!
!  Derivatives w.r.t. BOij (via deltai)
!
                d1is(nks) = d1is(nks) + detotdbo_tot*d1BOi(nk)
!
!  Derivatives w.r.t. BOij from angles and torsions via deltai and BOij
!
                d1is(nks) = d1is(nks) + (devalddeltai + devaldBO_neigh(ni) + devaldsbo_neigh(ni))*d1BOi(nk)
!
!  Derivatives w.r.t. BOij_pi and BOij_pipi
!
                d1is(nks) = d1is(nks) + detotdbopi_tot*(d1BOi_pi(nk) + d1BOi_pipi(nk))
!
!  Derivatives w.r.t. BOij_pi and BOij_pipi from angles
!
                d1is(nks) = d1is(nks) + devaldsbo*(d1BOi_pi(nk) + d1BOi_pipi(nk))*sbopi
!
!  Derivatives of etors with respect to BOpi_ij
!
                d1is(nks) = d1is(nks) + detorsdBOpi_neigh(ni)*d1BOi_pi(nk)
              enddo
!*********************************
!  Compute delta_j contribution  *
!*********************************
              do nj = 1,nneigh(j)
                k = neighno(nj,j)
!
!  Do we need to do this atom?
!
                if (nbos(k).gt.0) then
!
!  Set variables relating to k
!
                  nspeck = nbosptr(k)
!
!  Find j in neighbour list for k
!
                  nk = neighnoRptr(nj,j)
!
!  Set index to parameters for pairwise interaction of species
!
                  if (nspecj.ge.nspeck) then
                    nbojk = nspecj*(nspecj - 1)/2 + nspeck
                  else
                    nbojk = nspeck*(nspeck - 1)/2 + nspecj
                  endif
!
!  Set delta'
!
                  deltapk = deltap(k)
                  call reaxFF_bomds(j,k,nspecj,nspeck,nj,nk,deltapj,deltapk,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                                    d1BOp,d1BOp_pi,d1BOp_pipi,BOjk_s,BOjk_pi,BOjk_pipi,d1BOjk_s,d1BOjk_pi, &
                                    d1BOjk_pipi,nonzero_jk,nonzero_jk_j,nonzeroptr_jk,.true.)
!
!  Convert non-zero indices
!
                  do nl = 1,nonzero_jk_j
                    nonzeroptr_jk(nl) = indj + nonzeroptr_jk(nl)
                  enddo
                  ind1i = nneigh(i) + nneigh(i)*(nneigh(i)+1)/2 + nneigh(i)*mneigh*(mneigh + 1) + &
                         (ni - 1)*mneigh*mneigh + (nj - 1)*mneigh
                  do nl = nonzero_jk_j+1,nonzero_jk
                    nonzeroptr_jk(nl) = ind1i + nonzeroptr_jk(nl)
                  enddo
!                 
!  Add indices to sparsity pattern if necessary
!  
                  do nl = 1,nonzero_jk
                    l = nonzeroptr_jk(nl)
                    if (nonzerorptr(l).eq.0) then
                      nonzero = nonzero + 1 
                      nonzeroptr(nonzero) = l
                      nonzerorptr(l) = nonzero
                    endif
                  enddo
!                   
!  Check there is room in the sparse arrays
!
                  if (nonzero.gt.maxsparse) then
                    oldmaxsparse = maxsparse
                    oldmaxsparse2 = maxsparse*(maxsparse+1)/2
                    maxsparse = nonzero + 500
                    call realloc(d1is,maxsparse,status)
                    if (status.ne.0) call outofmemory('reaxFFds','d1is')
                    maxsparse2 = maxsparse*(maxsparse+1)/2
                    call realloc(d2is,maxsparse2,status)
                    if (status.ne.0) call outofmemory('reaxFFds','d2is')
                    call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
                    if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
                    call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
                    if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
                    call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
                    if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
!
                    d1is(oldmaxsparse+1:maxsparse) = 0.0_dp
                    d2is(oldmaxsparse2+1:maxsparse2) = 0.0_dp
                    d2evaldthetadBO_neigh(oldmaxsparse+1:maxsparse,1:maxneigh1a) = 0.0_dp
                    d2etorsdthetaddeltaj_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                    d2etorsdthetadBOpi_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                  endif
!
!  Calculate total bond order derivatives from sum of components
!
                  do nl = 1,nonzero_jk
                    d1BOjk(nl) = d1BOjk_s(nl) + d1BOjk_pi(nl) + d1BOjk_pipi(nl)
                  enddo
!
                  do nl = 1,nonzero_jk
                    l = nonzeroptr_jk(nl)
                    nls = nonzerorptr(l)
                    d1is(nls) = d1is(nls) + (detotddeltaj_tot + devalddeltaj_neigh(ni))*d1BOjk(nl)
                  enddo
!
!  Derivatives w.r.t. BOjl from torsions 
!
                  indbojk = ni*maxneigh + nj
                  do nl = 1,nonzero_jk
                    l = nonzeroptr_jk(nl)
                    nls = nonzerorptr(l)
                    d1is(nls) = d1is(nls) + devaldBO_neigh(indbojk)*d1BOjk(nl)
                  enddo
                endif
!
!  End of loop over neighbours of j
!
              enddo
!********************************
!  End of delta_j contribution  *
!********************************
              if (lgrad2) then
!
!  Calculate total bond order derivatives from sum of components
!
                do nk = 1,nonzero_ij*(nonzero_ij+1)/2
                  d2BOi(nk) = d2BOi_s(nk) + d2BOi_pi(nk) + d2BOi_pipi(nk)
                enddo
!
                ind3 = ni*(ni+1)/2
                inds = 0
                do nk = 1,nonzero_ij
                  k = nonzeroptr_ij(nk)
                  nks = nonzerorptr(k)
                  do nl = 1,nk
                    l = nonzeroptr_ij(nl)
                    nls = nonzerorptr(l)
                    if (nks.ge.nls) then
                      ind = nks*(nks-1)/2 + nls
                    else
                      ind = nls*(nls-1)/2 + nks
                    endif
                    inds = inds + 1
!
!  Add second derivative contributions that depend on bond order and pi bond orders including angles
!
                    d2is(ind) = d2is(ind) + detotdbo_tot*d2BOi(inds)
                    d2is(ind) = d2is(ind) + detotdbopi_tot*(d2BOi_pi(inds) + d2BOi_pipi(inds))
!
                    d2is(ind) = d2is(ind) + (devalddeltai + devaldBO_neigh(ni) + devaldsbo_neigh(ni))*d2BOi(inds)
                    d2is(ind) = d2is(ind) + devaldsbo*(d2BOi_pi(inds) + d2BOi_pipi(inds))*sbopi
!
                    d2is(ind) = d2is(ind) + d2etotdbo2_tot*d1BOi(nk)*d1BOi(nl)
                    d2is(ind) = d2is(ind) + d2etotdbodbopi_tot*(d1BOi_pi(nk)+d1BOi_pipi(nk))*d1BOi(nl)
                    d2is(ind) = d2is(ind) + d2etotdbodbopi_tot*(d1BOi_pi(nl)+d1BOi_pipi(nl))*d1BOi(nk)
                    d2is(ind) = d2is(ind) + d2etotdbopi2_tot*(d1BOi_pi(nk)+d1BOi_pipi(nk))*(d1BOi_pi(nl)+d1BOi_pipi(nl))
!
                    d2is(ind) = d2is(ind) + d2evaldsbo2*(d1BOi_pi(nk)+d1BOi_pipi(nk))*(d1BOi_pi(nl)+d1BOi_pipi(nl))*sbopi
                    d2is(ind) = d2is(ind) + (d2evalddeltaidsbo + d2evaldsbo_neigh(ni))*(d1BOi_pi(nk)+d1BOi_pipi(nk))*d1BOi(nl)*sbopi
                    d2is(ind) = d2is(ind) + (d2evalddeltaidsbo + d2evaldsbo_neigh(ni))*(d1BOi_pi(nl)+d1BOi_pipi(nl))*d1BOi(nk)*sbopi
!
!  Second derivatives w.r.t. BOij from angles via deltai and BOij
!
                    d2is(ind) = d2is(ind) + (d2evalddeltai2 + d2evaldBO2_neigh(ind3) + d2evaldsbo2_neigh(ind3))*d1BOi(nk)*d1BOi(nl)
                    d2is(ind) = d2is(ind) + 2.0_dp*(d2evalddeltaidBO_neigh(ni))*d1BOi(nk)*d1BOi(nl)
!
!  Add second derivative contributions that depend on pi bond orders from torsions
!
                    d2is(ind) = d2is(ind) + detorsdBOpi_neigh(ni)*d2BOi_pi(inds)
                    d2is(ind) = d2is(ind) + d2etorsdBOpi2_neigh(ni)*d1BOi_pi(nk)*d1BOi_pi(nl)
                  enddo
                enddo
!
                do nk = 1,nonzero_ij
                  k = nonzeroptr_ij(nk)
                  nks = nonzerorptr(k)
                  ind0 = nks*(nks-1)/2
                  do nl = 1,nonzero
                    l = nonzeroptr(nl)
                    if (l.le.nneighi2t) then
                      nls = nonzerorptr(l)
                      if (nks.ge.nls) then
                        ind = ind0 + nls
                      else
                        ind = nls*(nls-1)/2 + nks
                      endif
! 
!  Second derivatives from mixed angle bond order terms
! 
                      d2is(ind) = d2is(ind) + d2evaldthetadBO_neigh(nl,ni)*d1BOi(nk)
                    endif
                  enddo
                  ind = ind0 + nks
                  d2is(ind) = d2is(ind) + d2evaldthetadBO_neigh(nks,ni)*d1BOi(nk)
                enddo
!
!  More second derivatives from mixed angle bond order terms
!
                do nk = 1,nonzero_ij
                  k = nonzeroptr_ij(nk)
                  nks = nonzerorptr(k)
                  ind0 = nks*(nks-1)/2
                  do nl = 1,nonzero
                    l = nonzeroptr(nl)
                    if (l.le.nneighi2) then
                      nls = nonzerorptr(l)
                      if (nks.ge.nls) then
                        ind = ind0 + nls
                      else
                        ind = nls*(nls-1)/2 + nks
                      endif
                      d2is(ind) = d2is(ind) + d2evaldthetadBOpi_neigh(l)*(d1BOi_pi(nk) + d1BOi_pipi(nk))*sbopi
                      d2is(ind) = d2is(ind) + d2evaldthetaddeltai_neigh(l)*d1BOi(nk)
                    endif
                  enddo
                  ind = ind0 + nks
                  d2is(ind) = d2is(ind) + d2evaldthetadBOpi_neigh(k)*(d1BOi_pi(nk) + d1BOi_pipi(nk))*sbopi
                  d2is(ind) = d2is(ind) + d2evaldthetaddeltai_neigh(k)*d1BOi(nk)
                enddo
!
!  Second derivatives from mixed torsion bond order terms
!
                do nk = 1,nonzero_ij
                  k = nonzeroptr_ij(nk)
                  nks = nonzerorptr(k)
                  ind0 = nks*(nks-1)/2
                  do nl = 1,nonzero
                    l = nonzeroptr(nl)
                    nls = nonzerorptr(l)
                    if (nks.ge.nls) then
                      ind = ind0 + nls
                    else
                      ind = nls*(nls-1)/2 + nks
                    endif
                    d2is(ind) = d2is(ind) + d2etorsdthetadBOpi_neigh(nl,ni)*d1BOi_pi(nk)
                    d2is(ind) = d2is(ind) + d2etorsdthetaddeltai(l)*d1BOi(nk)
                  enddo
                  ind = ind0 + nks
                  d2is(ind) = d2is(ind) + d2etorsdthetadBOpi_neigh(nks,ni)*d1BOi_pi(nk)
                  d2is(ind) = d2is(ind) + d2etorsdthetaddeltai(k)*d1BOi(nk)
                enddo
!
                do nk = 1,nonzero_ij
                  k = nonzeroptr_ij(nk)
                  nks = nonzerorptr(k)
                  ind0 = nks*(nks-1)/2
                  do nl = 1,nk
                    l = nonzeroptr_ij(nl)
                    nls = nonzerorptr(l)
                    if (nks.gt.nls) then
                      ind = ind0 + nls
                    else
                      ind = nls*(nls-1)/2 + nks
                    endif
                    ind = ind0 + nls
                    d2is(ind) = d2is(ind) + d2etorsdBOdBOpi_neigh(ni,ni)*d1BOi_pi(nl)*d1BOi(nk)
                    d2is(ind) = d2is(ind) + d2etorsdBOdBOpi_neigh(ni,ni)*d1BOi_pi(nk)*d1BOi(nl)
!
                    d2is(ind) = d2is(ind) + d2etorsddeltaidBOpi_neigh(ni)*d1BOi_pi(nl)*d1BOi(nk)
                    d2is(ind) = d2is(ind) + d2etorsddeltaidBOpi_neigh(ni)*d1BOi_pi(nk)*d1BOi(nl)
                  enddo
                enddo
!--------------------------------------------
!  Loop over second neighbour of i up to j  |
!--------------------------------------------
                do nj = 1,ni-1
                  k = neighno(nj,i)
!
!  Do we need to do this atom?
!
                  if (nbos(k).gt.0) then
!
!  Set variables relating to k
!
                    nspeck = nbosptr(k)
!
!  Find i in neighbour list for k
!
                    nk = neighnoRptr(nj,i)
!
!  Set index to parameters for pairwise interaction of species
!
                    if (nspeci.ge.nspeck) then
                      nboik = nspeci*(nspeci - 1)/2 + nspeck
                    else
                      nboik = nspeck*(nspeck - 1)/2 + nspeci
                    endif
!
!  Set delta'
!
                    deltapk = deltap(k)
!***************************************
!  Valid pairwise reaxFF contribution  *
!***************************************
                    call reaxFF_bomds(i,k,nspeci,nspeck,nj,nk,deltapi,deltapk,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                                      d1BOp,d1BOp_pi,d1BOp_pipi,BOik_s,BOik_pi,BOik_pipi,d1BOik_s,d1BOik_pi, &
                                      d1BOik_pipi,nonzero_ik,nonzero_ik_i,nonzeroptr_ik,.true.)
!
!  Convert non-zero indices
!
                    ind0 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (nj - 1)*(mneigh+1)*mneigh + nj*mneigh
                    do nl = nonzero_ik_i+1,nonzero_ik
                      nonzeroptr_ik(nl) = ind0 + nonzeroptr_ik(nl)
                    enddo
!                 
!  Add indices to sparsity pattern if necessary
!               
                    do nl = 1,nonzero_ik
                      l = nonzeroptr_ik(nl)
                      if (nonzerorptr(l).eq.0) then
                        nonzero = nonzero + 1
                        nonzeroptr(nonzero) = l
                        nonzerorptr(l) = nonzero
                      endif 
                    enddo
!                   
!  Check there is room in the sparse arrays
!  
                    if (nonzero.gt.maxsparse) then
                      oldmaxsparse = maxsparse
                      oldmaxsparse2 = maxsparse*(maxsparse+1)/2
                      maxsparse = nonzero + 500
                      call realloc(d1is,maxsparse,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d1is')
                      maxsparse2 = maxsparse*(maxsparse+1)/2
                      call realloc(d2is,maxsparse2,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2is')
                      call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
                      call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
                      call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
!
                      d1is(oldmaxsparse+1:maxsparse) = 0.0_dp
                      d2is(oldmaxsparse2+1:maxsparse2) = 0.0_dp
                      d2evaldthetadBO_neigh(oldmaxsparse+1:maxsparse,1:maxneigh1a) = 0.0_dp
                      d2etorsdthetaddeltaj_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                      d2etorsdthetadBOpi_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                    endif
!
!  Calculate total bond order derivatives from sum of components
!
                    do nl = 1,nonzero_ik
                      d1BOik(nl) = d1BOik_s(nl) + d1BOik_pi(nl) + d1BOik_pipi(nl)
                    enddo
!
                    if (.not.lreaxFFunder(nspeci).or..not.lreaxFFunder(nspeck)) then
                      dsumdeltabopidbopik = delta(k)
                    else
                      dsumdeltabopidbopik = delta(k) - deltalp(k)
                    endif
!
                    d2etotdbopi2jk = d2etotddeltalpcorr2*ddeltalpcorrdsumdeltabopi**2 + &
                                     detotddeltalpcorr*d2deltalpcorrdsumdeltabopi2 + &
                                     2.0_dp*d2eunderddeltalpcorrdsumdeltabopi*ddeltalpcorrdsumdeltabopi + &
                                     d2eunderdsumdeltabopi2
!
                    d2etotdbopi2jk = d2etotdbopi2jk*dsumdeltabopidbopi*dsumdeltabopidbopik
!
                    dsumoverdboik = reaxFFoc2(nboik)*reaxFFDe(1,nboik)
!
                    d2etotdbodbopi_1 = (d2etotdbodsumdeltabopi+d2eoverdsumoverddeltalpcorr*dsumoverdboik* &
                                        ddeltalpcorrdsumdeltabopi)*dsumdeltabopidbopi
                    d2etotdbodbopi_2 = (d2etotdbodsumdeltabopi+d2eoverdsumoverddeltalpcorr*dsumoverdbo* &
                                        ddeltalpcorrdsumdeltabopi)*dsumdeltabopidbopik
                    d2etotdbo2ik = d2etotddeltalpcorr2*ddeltalpcorrddeltai**2 + &
                                   detotddeltalpcorr*d2deltalpcorrddeltai2 + &
                                   d2eloneddeltai2 + & 
                                   d2eoverdsumoverddeltalpcorr*dsumoverdboik*ddeltalpcorrddeltai + &
                                   d2eoverdsumoverddeltalpcorr*dsumoverdbo*ddeltalpcorrddeltai + &
                                   d2ebpendBO2_neigh(nj) + d2ebpenddeltaidBO_neigh(ni) + &
                                   2.0_dp*d2ebpenddeltaidBO_neigh(nj) + d2ebpenddeltai2
!
!  Add second derivative contributions that depend on bond orders 
!
                    ind3 = ni*(ni-1)/2 + nj
                    do nl = 1,nonzero_ij
                      l = nonzeroptr_ij(nl)
                      nls = nonzerorptr(l)
                      ind0 = nls*(nls-1)/2
                      do nm = 1,nonzero_ik
                        m = nonzeroptr_ik(nm)
                        nms = nonzerorptr(m)
                        if (nls.ge.nms) then
                          ind = ind0 + nms
                        else
                          ind = nms*(nms-1)/2 + nls
                        endif
                        if (nls.eq.nms) then
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etotdbo2ik*d1BOi(nl)*d1BOik(nm)
!                       
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etotdbopi2jk*(d1BOi_pi(nl)+d1BOi_pipi(nl))* &
                                                  (d1BOik_pi(nm)+d1BOik_pipi(nm))
!                       
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etotdbodbopi_1*(d1BOi_pi(nl)+d1BOi_pipi(nl))*d1BOik(nm)
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etotdbodbopi_2*(d1BOik_pi(nm)+d1BOik_pipi(nm))*d1BOi(nl)
!                                             
!  Second derivatives from angle terms
!  
                          d2is(ind) = d2is(ind) + 2.0_dp*(d2evalddeltai2 + d2evaldBO2_neigh(ind3) + &
                                                d2evaldsbo2_neigh(ind3))*d1BOi(nl)*d1BOik(nm)
!                       
                          d2is(ind) = d2is(ind) + 2.0_dp*(d2evalddeltaidBO_neigh(ni))*d1BOi(nl)*d1BOik(nm)
                          d2is(ind) = d2is(ind) + 2.0_dp*(d2evalddeltaidBO_neigh(nj))*d1BOi(nl)*d1BOik(nm)
!
                          d2is(ind) = d2is(ind) + 2.0_dp*d2evaldsbo2*(d1BOi_pi(nl)+d1BOi_pipi(nl))* &
                                                (d1BOik_pi(nm)+d1BOik_pipi(nm))*sbopi
!                       
                          d2is(ind) = d2is(ind) + 2.0_dp*(d2evalddeltaidsbo+d2evaldsbo_neigh(nj))* &
                                                (d1BOi_pi(nl)+d1BOi_pipi(nl))*d1BOik(nm)*sbopi
                          d2is(ind) = d2is(ind) + 2.0_dp*(d2evalddeltaidsbo+d2evaldsbo_neigh(ni))* &
                                                (d1BOik_pi(nm)+d1BOik_pipi(nm))*d1BOi(nl)*sbopi
!
!  Torsion mixed BO and BOpi terms
!
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etorsdBOdBOpi_neigh(nj,ni)*d1BOi_pi(nl)*d1BOik(nm)
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etorsdBOdBOpi_neigh(ni,nj)*d1BOik_pi(nm)*d1BOi(nl)
!
!  Mixed delta i terms with BOij_pi only
!
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etorsddeltaidBOpi_neigh(ni)*d1BOi_pi(nl)*d1BOik(nm)
                        else
                          d2is(ind) = d2is(ind) + d2etotdbo2ik*d1BOi(nl)*d1BOik(nm)
!
                          d2is(ind) = d2is(ind) + d2etotdbopi2jk*(d1BOi_pi(nl)+d1BOi_pipi(nl))*(d1BOik_pi(nm)+d1BOik_pipi(nm))
!
                          d2is(ind) = d2is(ind) + d2etotdbodbopi_1*(d1BOi_pi(nl)+d1BOi_pipi(nl))*d1BOik(nm)
                          d2is(ind) = d2is(ind) + d2etotdbodbopi_2*(d1BOik_pi(nm)+d1BOik_pipi(nm))*d1BOi(nl)
!
!  Second derivatives from angle terms
!
                          d2is(ind) = d2is(ind) + (d2evalddeltai2 + d2evaldBO2_neigh(ind3) + &
                                                d2evaldsbo2_neigh(ind3))*d1BOi(nl)*d1BOik(nm)
!
                          d2is(ind) = d2is(ind) + (d2evalddeltaidBO_neigh(ni))*d1BOi(nl)*d1BOik(nm)
                          d2is(ind) = d2is(ind) + (d2evalddeltaidBO_neigh(nj))*d1BOi(nl)*d1BOik(nm)
!
                          d2is(ind) = d2is(ind) + d2evaldsbo2*(d1BOi_pi(nl)+d1BOi_pipi(nl))*(d1BOik_pi(nm)+d1BOik_pipi(nm))*sbopi
!
                          d2is(ind) = d2is(ind) + (d2evalddeltaidsbo+d2evaldsbo_neigh(nj))* &
                                                (d1BOi_pi(nl)+d1BOi_pipi(nl))*d1BOik(nm)*sbopi
                          d2is(ind) = d2is(ind) + (d2evalddeltaidsbo+d2evaldsbo_neigh(ni))* &
                                                (d1BOik_pi(nm)+d1BOik_pipi(nm))*d1BOi(nl)*sbopi
!
!  Torsion mixed BO and BOpi terms
!
                          d2is(ind) = d2is(ind) + d2etorsdBOdBOpi_neigh(nj,ni)*d1BOi_pi(nl)*d1BOik(nm)
                          d2is(ind) = d2is(ind) + d2etorsdBOdBOpi_neigh(ni,nj)*d1BOik_pi(nm)*d1BOi(nl)
!
!  Mixed delta i terms with BOij_pi only
!
                          d2is(ind) = d2is(ind) + d2etorsddeltaidBOpi_neigh(ni)*d1BOi_pi(nl)*d1BOik(nm)
                          d2is(ind) = d2is(ind) + d2etorsddeltaidBOpi_neigh(nj)*d1BOik_pi(nm)*d1BOi(nl)
                        endif
                      enddo
                    enddo
                  endif
!
!  End of second loop over neighbours of i
!
                enddo
!*******************************************************
!  Compute delta_j contribution to second derivatives  *
!*******************************************************
                indni2 = ni*(ni+1)/2
                do nj = 1,nneigh(j)
                  k = neighno(nj,j)
!
!  Do we need to do this atom?
!
                  if (nbos(k).gt.0) then
!
!  Set variables relating to k
!
                    nspeck = nbosptr(k)
!
!  Find j in neighbour list for k
!
                    nk = neighnoRptr(nj,j)
!
!  Set index to parameters for pairwise interaction of species
!
                    if (nspecj.ge.nspeck) then
                      nbojk = nspecj*(nspecj - 1)/2 + nspeck
                    else
                      nbojk = nspeck*(nspeck - 1)/2 + nspecj
                    endif
!
!  Set delta'
!
                    deltapk = deltap(k)
!
                    call reaxFF_bos(j,k,nspecj,nspeck,nj,nk,deltapj,deltapk,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                                    d1BOp,d1BOp_pi,d1BOp_pipi,d2BOp,d2BOp_pi,d2BOp_pipi,BOjk_s,BOjk_pi,BOjk_pipi, &
                                    d1BOjk_s,d1BOjk_pi,d1BOjk_pipi,d2BOjk_s,d2BOjk_pi,d2BOjk_pipi,nonzero_jk, &
                                    nonzero_jk_j,nonzeroptr_jk,.true.,.true.)
!
!  Convert non-zero indices
!
                    do nl = 1,nonzero_jk_j
                      nonzeroptr_jk(nl) = indj + nonzeroptr_jk(nl)
                    enddo
                    ind1i = nneigh(i) + nneigh(i)*(nneigh(i)+1)/2 + nneigh(i)*mneigh*(mneigh + 1) + &
                            (ni - 1)*mneigh*mneigh + (nj - 1)*mneigh
                    do nl = nonzero_jk_j+1,nonzero_jk
                      nonzeroptr_jk(nl) = ind1i + nonzeroptr_jk(nl)
                    enddo
!                         
!  Add indices to sparsity pattern if necessary 
!                   
                    do nl = 1,nonzero_jk_j
                      l = nonzeroptr_jk(nl)
                      if (nonzerorptr(l).eq.0) then
                        nonzero = nonzero + 1   
                        nonzeroptr(nonzero) = l 
                        nonzerorptr(l) = nonzero
                      endif
                    enddo
!
!  Check there is room in the sparse arrays
!                         
                    if (nonzero.gt.maxsparse) then
                      oldmaxsparse = maxsparse
                      oldmaxsparse2 = maxsparse*(maxsparse+1)/2
                      maxsparse = nonzero + 500
                      call realloc(d1is,maxsparse,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d1is')
                      maxsparse2 = maxsparse*(maxsparse+1)/2
                      call realloc(d2is,maxsparse2,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2is')
                      call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
                      call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
                      call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
                      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
!
                      d1is(oldmaxsparse+1:maxsparse) = 0.0_dp
                      d2is(oldmaxsparse2+1:maxsparse2) = 0.0_dp
                      d2evaldthetadBO_neigh(oldmaxsparse+1:maxsparse,1:maxneigh1a) = 0.0_dp
                      d2etorsdthetaddeltaj_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                      d2etorsdthetadBOpi_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                    endif
!
!  Calculate total bond order derivatives from sum of components
!
                    do nl = 1,nonzero_jk
                      d1BOjk(nl) = d1BOjk_s(nl) + d1BOjk_pi(nl) + d1BOjk_pipi(nl)
                    enddo
                    do nl = 1,nonzero_jk*(nonzero_jk+1)/2
                      d2BOjk(nl) = d2BOjk_s(nl) + d2BOjk_pi(nl) + d2BOjk_pipi(nl)
                    enddo
!
                    inds = 0
                    do nm = 1,nonzero_jk
                      m = nonzeroptr_jk(nm)
                      nms = nonzerorptr(m)
                      ind0 = nms*(nms-1)/2
                      do nn = 1,nm
                        n = nonzeroptr_jk(nn)
                        nns = nonzerorptr(n)
                        if (nms.ge.nns) then
                          ind = ind0 + nns
                        else
                          ind = nns*(nns-1)/2 + nms
                        endif
                        inds = inds + 1
                        d2is(ind) = d2is(ind) + detotddeltaj_tot*d2BOjk(inds)
                        d2is(ind) = d2is(ind) + devalddeltaj_neigh(ni)*d2BOjk(inds)
                      enddo
                    enddo
!
                    d2etotddeltajdbopijk = detotddeltalpcorr*ddeltalpcorrdsumdeltabopi + deunderdsumdeltabopi
                    d2etotddeltajdbopijk = d2etotddeltajdbopijk*d2sumdeltabopiddeltajdbopi
!
                    indbojk = ni*maxneigh + nj
                    indbojk2 = indbojk*(indbojk+1)/2
                    indjk = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + ni*mneigh + nj
                    indjk2 = indjk*(indjk+1)/2
                    indnj2 = nj*(nj+1)/2
!
                    do nm = 1,nonzero_jk
                      m = nonzeroptr_jk(nm)
                      nms = nonzerorptr(m)
                      do nn = 1,nonzero_ij
                        n = nonzeroptr_ij(nn)
                        nns = nonzerorptr(n)
                        if (nms.ge.nns) then
                          ind = nms*(nms-1)/2 + nns
                        else
                          ind = nns*(nns-1)/2 + nms
                        endif
                        if (m.eq.n) then
!
!  Mixed delta j terms with BOij_pi only
!
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etorsddeltajdBOpi_neigh(ni)*d1BOi_pi(nn)*d1BOjk(nm)
!
!  Mixed delta j terms with BOij_pi and BOij_pipi
!
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etotddeltajdbopijk*(d1BOi_pi(nn)+d1BOi_pipi(nn))*d1BOjk(nm)
!
!  Torsion mixed BO and BOpi terms
!
                          d2is(ind) = d2is(ind) + 2.0_dp*d2etorsdBOdBOpi_neigh(indbojk,ni)*d1BOi_pi(nn)*d1BOjk(nm)
                        else
!
!  Mixed delta j terms with BOij_pi only
!
                          d2is(ind) = d2is(ind) + d2etorsddeltajdBOpi_neigh(ni)*d1BOi_pi(nn)*d1BOjk(nm)
!
!  Mixed delta j terms with BOij_pi and BOij_pipi
!
                          d2is(ind) = d2is(ind) + d2etotddeltajdbopijk*(d1BOi_pi(nn)+d1BOi_pipi(nn))*d1BOjk(nm)
!
!  Torsion mixed BO and BOpi terms
!
                          d2is(ind) = d2is(ind) + d2etorsdBOdBOpi_neigh(indbojk,ni)*d1BOi_pi(nn)*d1BOjk(nm)
                        endif
                      enddo
                    enddo
!
                    do nm = 1,nonzero_jk
                      m = nonzeroptr_jk(nm)
                      nms = nonzerorptr(m)
                      ind0 = nms*(nms-1)/2
                      do nn = 1,nonzero
                        n = nonzeroptr(nn)
                        nns = nonzerorptr(n)
                        if (nms.ge.nns) then
                          ind = ind0 + nns
                        else
                          ind = nns*(nns-1)/2 + nms
                        endif
!
!  Mixed delta j terms with angle from torsion
!
                        d2is(ind) = d2is(ind) + d2etorsdthetaddeltaj_neigh(nn,ni)*d1BOjk(nm)
!
!  Torsion second derivatives from mixed angle bond order terms
!
                        d2is(ind) = d2is(ind) + d2evaldthetadBO_neigh(nn,indbojk)*d1BOjk(nm)
                      enddo
!
                      ind = ind0 + nms
!
!  Mixed delta j terms with angle from torsion
!
                      d2is(ind) = d2is(ind) + d2etorsdthetaddeltaj_neigh(nms,ni)*d1BOjk(nm)
!
!  Torsion second derivatives from mixed angle bond order terms
!
                      d2is(ind) = d2is(ind) + d2evaldthetadBO_neigh(nms,indbojk)*d1BOjk(nm)
                    enddo
!
!  Torsion second derivatives
!
                    inds = 0
                    do nm = 1,nonzero_jk
                      m = nonzeroptr_jk(nm)
                      nms = nonzerorptr(m)
                      ind0 = nms*(nms-1)/2
                      do nn = 1,nm
                        n = nonzeroptr_jk(nn)
                        nns = nonzerorptr(n)
                        inds = inds + 1
                        if (nms.ge.nns) then
                          ind = ind0 + nns
                        else
                          ind = nns*(nns-1)/2 + nms
                        endif
                        d2is(ind) = d2is(ind) + devaldBO_neigh(indbojk)*d2BOjk(inds) 
                        d2is(ind) = d2is(ind) + d2evaldBO2_neigh(indbojk2)*d1BOjk(nm)*d1BOjk(nn)
                      enddo
                    enddo
!*********************************************
!  Second loop over neighbours of i (=> j2)  *
!*********************************************
                    do nj2 = 1,nneigh(i)
                      j2 = neighno(nj2,i)
                      indj2 = nneigh(i) + nneigh(i)*(nneigh(i)+1)/2 + (nj2-1)*mneigh*(mneigh + 1) + nj2*mneigh
                      if (ni.ge.nj2) then
                        indninj2 = ni*(ni-1)/2 + nj2
                      else
                        indninj2 = nj2*(nj2-1)/2 + ni
                      endif
!
!  Do we need to do this pair of atoms
!
                      if ((lopanyneigh(i).or.lopanyneigh(j2)).and.nbos(j2).gt.0) then
!
!  Set variables relating to j2
!
                        nspecj2 = nbosptr(j2)
!
!  Find j2 in neighbour list for i
!
                        nj2i = neighnoRptr(nj2,i)
!
!  Set index to parameters for pairwise interaction of species
!
                        if (nspeci.ge.nspecj2) then
                          nboij2 = nspeci*(nspeci - 1)/2 + nspecj2
                        else
                          nboij2 = nspecj2*(nspecj2 - 1)/2 + nspeci
                        endif
!
!  Set delta'
!
                        deltapj2 = deltap(j2)
                        call reaxFF_bomds(i,j2,nspeci,nspecj2,nj2,nj2i,deltapi,deltapj2,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                                          d1BOp,d1BOp_pi,d1BOp_pipi,BOij2_s,BOij2_pi,BOij2_pipi,d1BOij2_s,d1BOij2_pi, &
                                          d1BOij2_pipi,nonzero_ij2,nonzero_ij2_i,nonzeroptr_ij2,.true.)
!  
!  Convert non-zero indices
! 
                        ind0 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (nj2 - 1)*(mneigh+1)*mneigh + nj2*mneigh
                        do nl = nonzero_ij2_i+1,nonzero_ij2
                          nonzeroptr_ij2(nl) = ind0 + nonzeroptr_ij2(nl)
                        enddo
!                        
!  Add indices to sparsity pattern if necessary
!                   
                        do nl = 1,nonzero_ij2 
                          l = nonzeroptr_ij2(nl)
                          if (nonzerorptr(l).eq.0) then
                            nonzero = nonzero + 1   
                            nonzeroptr(nonzero) = l 
                            nonzerorptr(l) = nonzero
                          endif
                        enddo
!                   
!  Check there is room in the sparse arrays
!                         
                        if (nonzero.gt.maxsparse) then
                          oldmaxsparse = maxsparse
                          oldmaxsparse2 = maxsparse*(maxsparse+1)/2
                          maxsparse = nonzero + 500
                          call realloc(d1is,maxsparse,status)
                          if (status.ne.0) call outofmemory('reaxFFds','d1is')
                          maxsparse2 = maxsparse*(maxsparse+1)/2
                          call realloc(d2is,maxsparse2,status)
                          if (status.ne.0) call outofmemory('reaxFFds','d2is')
                          call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
                          if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
                          call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
                          if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
                          call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
                          if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
!
                          d1is(oldmaxsparse+1:maxsparse) = 0.0_dp
                          d2is(oldmaxsparse2+1:maxsparse2) = 0.0_dp
                          d2evaldthetadBO_neigh(oldmaxsparse+1:maxsparse,1:maxneigh1a) = 0.0_dp
                          d2etorsdthetaddeltaj_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                          d2etorsdthetadBOpi_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                        endif
!
!  Calculate total bond order derivatives from sum of components
!
                        do nl = 1,nonzero_ij2
                          d1BOij2(nl) = d1BOij2_s(nl) + d1BOij2_pi(nl) + d1BOij2_pipi(nl)
                        enddo
!
!  Set up derivative of delta_lp for j with respect to delta_e
!
                        delta_e = 0.5_dp*(delta(j2) + reaxFFval(1,nspecj2) - reaxFFval(3,nspecj2))
!
                        if (lreaxFFovsmooth) then
                          delta_e0 = nint(delta_e)
                          expH = exp(-reaxFFlpsmooth*(delta_e - delta_e0))
                          smoothH = 1.0_dp/(1.0_dp + expH)
                          rintde = delta_e0 - 1.0_dp + smoothH
                          drintdeddeltae = reaxFFlpsmooth*expH*smoothH**2
                          if (lgrad2) then
                            d2rintdeddeltae2 = reaxFFlpsmooth*reaxFFlpsmooth*expH*(smoothH**2)* &
                                               (2.0_dp*expH*smoothH - 1.0_dp)
                          endif
                        else
                          rintde = dble(int(delta_e))
                          drintdeddeltae = 0.0_dp
                          d2rintdeddeltae2 = 0.0_dp
                        endif
!
                        trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
                        explp2 = exp(-reaxFFlam(29)*trm2lp**2)
                        ddeltalpddeltaj2 = 0.5_dp*(drintdeddeltae + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(1.0_dp - drintdeddeltae)) ! 1/2 is to convert from delta_e to delta_i
!
!  Derivatives of S w.r.t. deltaj2 and BOij pi and pipi
!
                        if (.not.lreaxFFunder(nspeci).or..not.lreaxFFunder(nspecj2)) then
                          dsumdeltabopiddeltaj2 = BO_pi(nj2,i) + BO_pipi(nj2,i)
                          dsumdeltabopidbopij2 = delta(j2)
                        else                       
                          dsumdeltabopiddeltaj2 = (1.0_dp - ddeltalpddeltaj2)*(BO_pi(nj2,i) + BO_pipi(nj2,i))
                          dsumdeltabopidbopij2 = delta(j2) - deltalp(j2)
                        endif
!
!  Set derivatives that depend on deltaj2
!
                        ddeltalpcorrddeltaj2 = ddeltalpcorrdsumdeltabopi*dsumdeltabopiddeltaj2
                        d2deltalpcorrddeltajddeltaj2 = d2deltalpcorrdsumdeltabopi2*dsumdeltabopiddeltaj*dsumdeltabopiddeltaj2
!
                        d2etotddeltajddeltaj2 = d2eunderdsumdeltabopi2*dsumdeltabopiddeltaj*dsumdeltabopiddeltaj2
                        d2etotddeltajddeltaj2_tot = d2etotddeltalpcorr2*ddeltalpcorrddeltaj*ddeltalpcorrddeltaj2 + &
                                detotddeltalpcorr*d2deltalpcorrddeltajddeltaj2 + &
                                2.0_dp*d2eunderddeltalpcorrdsumdeltabopi*ddeltalpcorrddeltaj*dsumdeltabopiddeltaj2 + &
                                d2etotddeltajddeltaj2
!
                        dsumoverdboj2 = reaxFFoc2(nboij2)*reaxFFDe(1,nboij2)
                        d2etotddeltajdbojk = d2etotddeltajdbo_tot + &
                          d2eoverdsumoverddeltalpcorr*dsumoverdboj2*ddeltalpcorrddeltaj
                        dtrm = d2etotdsumdeltabopi2_tot*dsumdeltabopidbopij2*dsumdeltabopiddeltaj
!
                        indboij2jk = indbojk*(indbojk-1)/2 + nj2
!
                        do nm = 1,nonzero_jk
                          m = nonzeroptr_jk(nm)
                          nms = nonzerorptr(m)
                          do nn = 1,nonzero_ij2
                            n = nonzeroptr_ij2(nn)
                            nns = nonzerorptr(n)
                            if (nms.ge.nns) then
                              ind = nms*(nms-1)/2 + nns
                            else
                              ind = nns*(nns-1)/2 + nms
                            endif
                            if (m.eq.n) then
!
!  Mixed derivatives between delta j and BO pi
!
                              d2is(ind) = d2is(ind) + 2.0_dp*dtrm*(d1BOij2_pi(nn)+d1BOij2_pipi(nn))*d1BOjk(nm)
!
!  Mixed derivatives between delta j and delta i
!
                              d2is(ind) = d2is(ind) + 2.0_dp*d2etotddeltajdbojk*d1BOij2(nn)*d1BOjk(nm)
!
!  Threebody conjugation energy mixed deltaj/BOik terms
!
                              d2is(ind) = d2is(ind) + 2.0_dp*d2evalddeltajdBO_neigh(ni,nj2)*d1BOij2(nn)*d1BOjk(nm)
!
!  Mixed delta i - delta j term
!
                              d2is(ind) = d2is(ind) + 2.0_dp*d2evalddeltaiddeltaj_neigh(ni)*d1BOij2(nn)*d1BOjk(nm)
!
!  Torsion second derivatives from i-j2 and j-k
!
                              d2is(ind) = d2is(ind) + 2.0_dp*d2evaldBO2_neigh(indboij2jk)*d1BOij2(nn)*d1BOjk(nm)
                              d2is(ind) = d2is(ind) + 2.0_dp*d2evalddeltaidBO_neigh(indbojk)*d1BOij2(nn)*d1BOjk(nm)
                            else
!
!  Mixed derivatives between delta j and BO pi
!
                              d2is(ind) = d2is(ind) + dtrm*(d1BOij2_pi(nn)+d1BOij2_pipi(nn))*d1BOjk(nm)
!
!  Mixed derivatives between delta j and delta i
!
                              d2is(ind) = d2is(ind) + d2etotddeltajdbojk*d1BOij2(nn)*d1BOjk(nm)
!
!  Threebody conjugation energy mixed deltaj/BOik terms
!
                              d2is(ind) = d2is(ind) + d2evalddeltajdBO_neigh(ni,nj2)*d1BOij2(nn)*d1BOjk(nm)
!
!  Mixed delta i - delta j term
!
                              d2is(ind) = d2is(ind) + d2evalddeltaiddeltaj_neigh(ni)*d1BOij2(nn)*d1BOjk(nm)
!
!  Torsion second derivatives from i-j2 and j-k
!
                              d2is(ind) = d2is(ind) + d2evaldBO2_neigh(indboij2jk)*d1BOij2(nn)*d1BOjk(nm)
                              d2is(ind) = d2is(ind) + d2evalddeltaidBO_neigh(indbojk)*d1BOij2(nn)*d1BOjk(nm)
                            endif
                          enddo
                        enddo
!**********************************
!  Compute delta_j2 contribution  *
!**********************************
                        do nk2 = 1,nneigh(j2)
                          k2 = neighno(nk2,j2)
!
!  Do we need to do this atom?
!
                          if (nbos(k2).gt.0) then
!
!  Set variables relating to k2
!
                            nspeck2 = nbosptr(k2)
!
!  Find k2 in neighbour list for j2
!
                            nk2j2 = neighnoRptr(nk2,j2)
!
!  Set index to parameters for pairwise interaction of species
!
                            if (nspecj2.ge.nspeck2) then
                              nboj2k2 = nspecj2*(nspecj2 - 1)/2 + nspeck2
                            else
                              nboj2k2 = nspeck2*(nspeck2 - 1)/2 + nspecj2
                            endif
!
!  Set delta'
!
                            deltapk2 = deltap(k2)
                            call reaxFF_bomds(j2,k2,nspecj2,nspeck2,nk2,nk2j2,deltapj2,deltapk2,mneigh,nneigh,BOp,BOp_pi, &
                                              BOp_pipi,d1BOp,d1BOp_pi,d1BOp_pipi,BOjk2_s,BOjk2_pi,BOjk2_pipi,d1BOjk2_s, &
                                              d1BOjk2_pi,d1BOjk2_pipi,nonzero_j2k2,nonzero_j2k2_j2,nonzeroptr_j2k2,.true.)
!  
!  Convert non-zero indices
! 
                            do nl = 1,nonzero_j2k2_j2
                              nonzeroptr_j2k2(nl) = indj2 + nonzeroptr_j2k2(nl)
                            enddo
                            ind1i2 = nneigh(i) + nneigh(i)*(nneigh(i)+1)/2 + nneigh(i)*mneigh*(mneigh + 1) + &
                                    (nj2 - 1)*mneigh*mneigh + (nk2 - 1)*mneigh
                            do nl = nonzero_j2k2_j2+1,nonzero_j2k2
                              nonzeroptr_j2k2(nl) = ind1i2 + nonzeroptr_j2k2(nl)
                            enddo
!                        
!  Add indices to sparsity pattern if necessary
!                   
                            do nl = 1,nonzero_j2k2 
                              l = nonzeroptr_j2k2(nl)
                              if (nonzerorptr(l).eq.0) then
                                nonzero = nonzero + 1   
                                nonzeroptr(nonzero) = l 
                                nonzerorptr(l) = nonzero
                              endif
                            enddo
!                   
!  Check there is room in the sparse arrays
!                         
                            if (nonzero.gt.maxsparse) then
                              oldmaxsparse = maxsparse
                              oldmaxsparse2 = maxsparse*(maxsparse+1)/2
                              maxsparse = nonzero + 500
                              call realloc(d1is,maxsparse,status)
                              if (status.ne.0) call outofmemory('reaxFFds','d1is')
                              maxsparse2 = maxsparse*(maxsparse+1)/2
                              call realloc(d2is,maxsparse2,status)
                              if (status.ne.0) call outofmemory('reaxFFds','d2is')
                              call realloc(d2evaldthetadBO_neigh,maxsparse,maxneigh1a,status)
                              if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
                              call realloc(d2etorsdthetaddeltaj_neigh,maxsparse,maxneigh,status)
                              if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
                              call realloc(d2etorsdthetadBOpi_neigh,maxsparse,maxneigh,status)
                              if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
!
                              d1is(oldmaxsparse+1:maxsparse) = 0.0_dp
                              d2is(oldmaxsparse2+1:maxsparse2) = 0.0_dp
                              d2evaldthetadBO_neigh(oldmaxsparse+1:maxsparse,1:maxneigh1a) = 0.0_dp
                              d2etorsdthetaddeltaj_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                              d2etorsdthetadBOpi_neigh(oldmaxsparse+1:maxsparse,1:maxneigh) = 0.0_dp
                            endif
!
!  Calculate total bond order derivatives from sum of components
!
                            do nl = 1,nonzero_j2k2
                              d1BOjk2(nl) = d1BOjk2_s(nl) + d1BOjk2_pi(nl) + d1BOjk2_pipi(nl)
                            enddo
!
!  Second derivatives from delta j terms
!
                            do nm = 1,nonzero_jk
                              m = nonzeroptr_jk(nm)
                              nms = nonzerorptr(m)
                              do nn = 1,nonzero_j2k2
                                n = nonzeroptr_j2k2(nn)
                                nns = nonzerorptr(n)
                                if (nms.ge.nns) then
                                  ind = nms*(nms-1)/2 + nns
                                else
                                  ind = nns*(nns-1)/2 + nms
                                endif
                                if (nms.eq.nns) then
                                  half = 1.0_dp
                                else
                                  half = 0.5_dp
                                endif
                                d2is(ind) = d2is(ind) + half*d2etotddeltajddeltaj2_tot*d1BOjk(nm)*d1BOjk2(nn)
                                d2is(ind) = d2is(ind) + half*d2evalddeltajk_neigh(indninj2)*d1BOjk(nm)*d1BOjk2(nn)
                              enddo
                            enddo
!
                            if (j.eq.j2) then
!
!  Torsion second derivatives from delta j and bond order of j
!
                              indbojk2 = nj2*maxneigh + nk2
                              do nm = 1,nonzero_jk
                                m = nonzeroptr_jk(nm)
                                nms = nonzerorptr(m)
                                do nn = 1,nonzero_j2k2
                                  n = nonzeroptr_j2k2(nn)
                                  nns = nonzerorptr(n)
                                  if (nms.ge.nns) then
                                    ind = nms*(nms-1)/2 + nns
                                  else
                                    ind = nns*(nns-1)/2 + nms
                                  endif
                                  if (nms.eq.nns) then
                                    d2is(ind) = d2is(ind) + 2.0_dp*d2evalddeltajdBO_neigh(ni,indbojk2)*d1BOjk(nm)*d1BOjk2(nn)
                                    d2is(ind) = d2is(ind) + d2etotddeltaj2_tot*d1BOjk(nm)*d1BOjk2(nn)
                                  else
                                    d2is(ind) = d2is(ind) + d2evalddeltajdBO_neigh(ni,indbojk2)*d1BOjk(nm)*d1BOjk2(nn)
                                    d2is(ind) = d2is(ind) + 0.5_dp*d2etotddeltaj2_tot*d1BOjk(nm)*d1BOjk2(nn)
                                  endif
                                enddo
                              enddo
                            endif
                          endif
!
!  End of loop over neighbours of j2
!
                        enddo
!*********************************
!  End of delta_j2 contribution  *
!*********************************
!
!  End condition section on i or j2 being associated with moving atom
!
                      endif
!
!  End of second loop over neighbours of i
!
                    enddo
                  endif
!
!  End of loop over neighbours of j
!
                enddo
!********************************
!  End of delta_j contribution  *
!********************************
              endif
!
!  End condition section on i or j being associated with moving atom
!
            endif
!
!  End of loop over neighbours of i
!
          enddo
!
!  Add derivatives due to neighbours of i - here we include torsional derivatives
!
          if (lilocal) then
            call d1addsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nbodummyptr,d1is, &
                         nonzero,nonzeroptr,.true.,lextended,.true.)
          endif
          if (lgrad2) then
            call d2adddmsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nboatomRptr, &
                           d1is,d2is,nonzero,nonzeroptr,.true.,lextended)
          endif
!
!  End condition on lgrad1 being true
!
        endif
!
!  End of QMMM condition
!
      endif
!
!  End of condition on i having a ReaxFF species type
!
    endif
  enddo
!**********************************
!  Hydrogen bonding contribution  *
!**********************************
  cut2 = reaxFFrhtol**2
  if (lreaxff_ehb) then
!***************************
!  Standard N^2 algorithm  *
!***************************
!
!  Loop over atoms looking for hydrogen
!
    do ii = 1,nboatom
      i = nboatomptr(ii)
      if (nat(i).eq.1) then
        iloc = atom2local(i)
        lilocal = (iloc.gt.0)
!
!  Set variables relating to i
!
        nspeci = nbosptr(i)
        nregioni = nregionno(nrelf2a(i))
        nregiontypi = nregiontype(nregioni,ncf)
        lslicei = lsliceatom(nsft + nrelf2a(i))
        deltapi = deltap(i)
!
!  QM/MM handling : i is a QM atom => exclude
!
        lQMMMok = .true.
        if (QMMMmode(ncf).gt.0) then
          if (nregiontypi.eq.1) lQMMMok = .false.
        endif
!
!  Do we need to do this atom based on QMMM scheme?
!
        if (lQMMMok) then
!           
!  Setup in case we need derivatives 
!           
!  Initialise derivative storage for neighbours of i
!       
          if (lgrad2) then
            ix = 3*(iloc-1) + 1
            iy = ix + 1
            iz = ix + 2
            ixf = 3*(i-1) + 1
            iyf = ixf + 1
            izf = ixf + 2
          endif
!
!  Loop over atoms that are connected by bond order to i
!
          do ni = 1,nneigh(i)
            j = neighno(ni,i)
            if (nbos(j).gt.0) then
              nspecj = nbosptr(j)
!
              jloc = atom2local(j)
              ljlocal = (jloc.gt.0)
!
!  Is the bond order above the threshold for a hydrogen bond?
!
              if (BO(ni,i).gt.reaxFFhtol) then
!
!  Is there any H-bonding term for this i-j pair?
!
                lAnyHB = .false.
                do nspeck = 1,nreaxFFspec
                  phb1 = reaxFFhb3(2,nspeck,nspecj,nspeci)
                  if (abs(phb1).gt.1.0d-12) lAnyHB = .true.
                enddo
                if (lAnyHB) then
                  xji = xneigh(ni,i)
                  yji = yneigh(ni,i)
                  zji = zneigh(ni,i)
                  rji = rneigh(ni,i)
!                 
!  Calculate taper function to smooth cutoff
!                 
                  if (lStrict) then
                    fHB = 1.0_dp
                    dfHBdBO = 0.0_dp
                    d2fHBdBO2 = 0.0_dp
                  else 
                    call ataper(.false.,BO(ni,i),reaxFFhtol,reaxFFtaperscale*reaxFFhtol,fHB,dfHBdBO,d2fHBdBO2,d3fHBdBO3, &
                                lgrad1,lgrad2,.false.)
                  endif
!
                  if (lgrad1) then
!
!  Set region for j
!
                    nregionj = nregionno(nrelf2a(j))
!
!  Find i in neighbour list for j
!
                    nj = neighnoRptr(ni,i)
!
!  Set delta'
!
                    deltapj = deltap(j)
!
!  Compute bond order derivatives
!
                    call reaxFF_bos(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                                    d1BOp,d1BOp_pi,d1BOp_pipi,d2BOp,d2BOp_pi,d2BOp_pipi,BOij_s,BOij_pi,BOij_pipi, &
                                    d1BOi_s,d1BOi_pi,d1BOi_pipi,d2BOi_s,d2BOi_pi,d2BOi_pipi,nonzero_ij,nonzero_ij_i, &
                                    nonzeroptr_ij,lgrad1,lgrad2)
! 
!  Convert non-zero indices
! 
                    ind0 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + ni*mneigh
                    do nl = nonzero_ij_i+1,nonzero_ij
                      nonzeroptr_ij(nl) = ind0 + nonzeroptr_ij(nl)
                    enddo     
!                 
!  Initialise derivative storage for neighbours of i
!
                    d1is(1:nonzero_ij) = 0.0_dp
                    if (lgrad2) then
                      d2is(1:nonzero_ij*(nonzero_ij+1)/2) = 0.0_dp
                      d2mix(1:nonzero_ij) = 0.0_dp
                    endif
!
!  Calculate total bond order derivatives from sum of components
!
                    do nk = 1,nonzero_ij
                      d1BOi(nk) = d1BOi_s(nk) + d1BOi_pi(nk) + d1BOi_pipi(nk)
                    enddo
                    if (lgrad2) then
                      jx = 3*(jloc-1) + 1
                      jy = jx + 1
                      jz = jx + 2
                      jxf = 3*(j-1) + 1
                      jyf = jxf + 1
                      jzf = jxf + 2
                      do nk = 1,nonzero_ij*(nonzero_ij+1)/2
                        d2BOi(nk) = d2BOi_s(nk) + d2BOi_pi(nk) + d2BOi_pipi(nk)
                      enddo
                    endif
                  endif
!
!  Loop over atoms looking for an acceptor
!
                  kloop: do k = 1,numat
!
!  Does k have a reaxFF species?
!
                    if (nbos(k).gt.0) then
                      nspeck = nbosptr(k)
!
                      kloc = atom2local(k)
                      lklocal = (kloc.gt.0)
!
!  Is this triad valid for a hydrogen bond?
!
                      phb1 = reaxFFhb3(2,nspeck,nspecj,nspeci)
                      if (abs(phb1).lt.1.0d-12) cycle kloop
!
!  Compute basic interatomic vector
!
                      xkj = xclat(k) - xclat(j)
                      ykj = yclat(k) - yclat(j)
                      zkj = zclat(k) - zclat(j)
!
!  Find valid vectors
!
                      nor = 0
                      if (ndim.eq.3) then
                        call rsearch3D(xkj,ykj,zkj,.false.,.false.,j,k,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
                      elseif (ndim.eq.2) then
                        call rsearch2D(xkj,ykj,zkj,.false.,.false.,j,k,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
                      elseif (ndim.eq.1) then
                        call rsearch1D(xkj,ykj,zkj,.false.,.false.,j,k,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,0.0_dp)
                      elseif (ndim.eq.0) then
                        rkj2 = xkj*xkj + ykj*ykj + zkj*zkj
                        if (rkj2.lt.cut2) then
                          nor = nor + 1
                        endif
                      endif
!
!  If no distances then cycle
!
                      if (nor.eq.0) cycle kloop
!
                      if (lgrad2) then
                        kx = 3*(kloc-1) + 1
                        ky = kx + 1
                        kz = kx + 2
                        kxf = 3*(k-1) + 1
                        kyf = kxf + 1
                        kzf = kxf + 2
                      endif
!
!  Set remaining parameters
!
                      r0hb = reaxFFhb3(1,nspeck,nspecj,nspeci)
                      phb2 = reaxFFhb3(3,nspeck,nspecj,nspeci)
                      phb3 = reaxFFhb3(4,nspeck,nspecj,nspeci)
!
!  Set region for k
!
                      nregionk = nregionno(nrelf2a(k))
!
!  Loop over valid distances and calculate contributions
!
                      nloop: do n = 1,nor
                        if (ndim.gt.0) rkj2 = dist(n)
                        if (ndim.gt.0) then
                          xkj = xtmp(n)
                          ykj = ytmp(n)
                          zkj = ztmp(n)
                        endif
!
!  Calculate distances for i - k
!
                        xki = xkj + xji
                        yki = ykj + yji
                        zki = zkj + zji
                        rki2 = xki*xki + yki*yki + zki*zki
!
!  Exclude i = k case and j = k case
!
                        if (i.eq.k.and.rki2.lt.1.0d-12) cycle nloop
                        if (j.eq.k.and.rkj2.lt.1.0d-12) cycle nloop
!
                        rki = sqrt(rki2)
                        rkj = sqrt(rkj2)
!
                        rrij = 1.0_dp/rji
                        rrik = 1.0_dp/rki
                        rrjk = 1.0_dp/rkj
!                         
!  Compute distance taper function 
!
                        if (lStrict) then
                          frHB = 1.0_dp
                          dfrHBdr = 0.0_dp
                          d2frHBdr2 = 0.0_dp
                        else
                          call ataper(.true.,rkj,reaxFFrhtollower,reaxFFrhtol,frHB,dfrHBdr,d2frHBdr2,d3frHBdr3, &
                                      lgrad1,lgrad2,.false.)
!
!  Convert taper derivatives to correct form
!
                          if (lgrad1) then
                            dfrHBdr = rrjk*dfrHBdr
                            if (lgrad2) then
                              d2frHBdr2 = rrjk*rrjk*(d2frHBdr2 - dfrHBdr)
                            endif
                          endif
                        endif
!  
!  Compute theta for j - i - k
!                     
                        call reaxff_theta(rji,rki,rkj,theta,dthetadrij,dthetadrik,dthetadrjk, &
                                          d2thetadr2,lgrad1,lgrad2)
!
!  Compute remaining terms for hydrogen bond energy
!
                        sinhalftheta  = sin(0.5_dp*theta)
                        sin4halftheta = sinhalftheta**4
                        exphb1 = exp(-phb2*BO(ni,i))
                        exphb2 = exp(-phb3*((r0hb*rrik) + (rki/r0hb) - 2.0_dp))
!
                        if (lilocal) then
                          ehb = ehb + fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta
                        endif
!
                        if (lgrad1) then
!**********************
!  First derivatives  *
!**********************
!
!  Hydrogen bond derivatives - general terms
!
                          dexphb1dBOij = - phb2*exphb1
                          dexphb2drik = - phb3*exphb2*(-r0hb*(rrik**3) + rrik/r0hb)
                          coshalftheta = cos(0.5_dp*theta)
                          dsin4halfthetadtheta = (2.0_dp*coshalftheta*sinhalftheta**3)
!
                          if (lstr.or.lgrad2) then
                            call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,1), &
                              d2r2dx2(1,1),d2r2dsdx(1,1,1),d2r2ds2(1,1,1),lgrad2)
                            call real1strterm(ndim,xki,yki,zki,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,2), &
                              d2r2dx2(1,2),d2r2dsdx(1,1,2),d2r2ds2(1,1,2),lgrad2)
                            call real1strterm(ndim,xkj,ykj,zkj,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,3), &
                              d2r2dx2(1,3),d2r2dsdx(1,1,3),d2r2ds2(1,1,3),lgrad2)
                          endif
!
!  Hydrogen bond derivatives with respect to the distance i-j
!
                          dehbdrij = fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrij
                          dehbdrik = fHB*frHB*phb1*(1.0_dp - exphb1)*dexphb2drik*sin4halftheta + &
                                     fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrik
                          dehbdrjk = fHB*dfrHBdr*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta + &
                                     fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrjk
!
                          if (lilocal) then
!
!  Hydrogen bond derivatives with respect to the distance i-j
!
                            xdrv(i) = xdrv(i) - dehbdrij*xji
                            ydrv(i) = ydrv(i) - dehbdrij*yji
                            zdrv(i) = zdrv(i) - dehbdrij*zji
                            xdrv(j) = xdrv(j) + dehbdrij*xji
                            ydrv(j) = ydrv(j) + dehbdrij*yji
                            zdrv(j) = zdrv(j) + dehbdrij*zji
!
                            if (nregioni.ne.nregionj) then
                              xregdrv(nregioni) = xregdrv(nregioni) - dehbdrij*xji
                              yregdrv(nregioni) = yregdrv(nregioni) - dehbdrij*yji
                              zregdrv(nregioni) = zregdrv(nregioni) - dehbdrij*zji
                              xregdrv(nregionj) = xregdrv(nregionj) + dehbdrij*xji
                              yregdrv(nregionj) = yregdrv(nregionj) + dehbdrij*yji
                              zregdrv(nregionj) = zregdrv(nregionj) + dehbdrij*zji
                            endif
!
                            if (lstr) then
                              rstrdloc(1:nstrains) = 0.0_dp
                              do kl = 1,nstrains
                                ns = nstrptr(kl)
                                rstrdloc(kl) = rstrdloc(kl) + dehbdrij*dr2ds(ns,1)
                              enddo
                              do kl = 1,nstrains
                                rstrd(kl) = rstrd(kl) + rstrdloc(kl)
                              enddo
                              if (latomicstress) then
                                do kl = 1,nstrains
                                  atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*rstrdloc(kl)
                                  atomicstress(kl,k) = atomicstress(kl,k) + 0.5_dp*rstrdloc(kl)
                                enddo
                              endif
                            endif
!
!  Hydrogen bond derivatives with respect to the distance i-k
!
                            xdrv(i) = xdrv(i) - dehbdrik*xki
                            ydrv(i) = ydrv(i) - dehbdrik*yki
                            zdrv(i) = zdrv(i) - dehbdrik*zki
                            xdrv(k) = xdrv(k) + dehbdrik*xki
                            ydrv(k) = ydrv(k) + dehbdrik*yki
                            zdrv(k) = zdrv(k) + dehbdrik*zki
!
                            if (nregioni.ne.nregionk) then
                              xregdrv(nregioni) = xregdrv(nregioni) - dehbdrik*xki
                              yregdrv(nregioni) = yregdrv(nregioni) - dehbdrik*yki
                              zregdrv(nregioni) = zregdrv(nregioni) - dehbdrik*zki
                              xregdrv(nregionk) = xregdrv(nregionk) + dehbdrik*xki
                              yregdrv(nregionk) = yregdrv(nregionk) + dehbdrik*yki
                              zregdrv(nregionk) = zregdrv(nregionk) + dehbdrik*zki
                            endif
!
                            if (lstr) then
                              rstrdloc(1:nstrains) = 0.0_dp
                              do kl = 1,nstrains
                                ns = nstrptr(kl)
                                rstrdloc(kl) = rstrdloc(kl) + dehbdrik*dr2ds(ns,2)
                              enddo
                              do kl = 1,nstrains 
                                rstrd(kl) = rstrd(kl) + rstrdloc(kl)
                              enddo
                              if (latomicstress) then 
                                do kl = 1,nstrains
                                  atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*rstrdloc(kl)
                                  atomicstress(kl,k) = atomicstress(kl,k) + 0.5_dp*rstrdloc(kl)
                                enddo
                              endif
                            endif
!                           
!  Hydrogen bond derivatives with respect to the distance j-k
!                           
                            xdrv(j) = xdrv(j) - dehbdrjk*xkj
                            ydrv(j) = ydrv(j) - dehbdrjk*ykj
                            zdrv(j) = zdrv(j) - dehbdrjk*zkj
                            xdrv(k) = xdrv(k) + dehbdrjk*xkj
                            ydrv(k) = ydrv(k) + dehbdrjk*ykj
                            zdrv(k) = zdrv(k) + dehbdrjk*zkj
!
                            if (nregionj.ne.nregionk) then
                              xregdrv(nregionj) = xregdrv(nregionj) - dehbdrjk*xkj
                              yregdrv(nregionj) = yregdrv(nregionj) - dehbdrjk*ykj
                              zregdrv(nregionj) = zregdrv(nregionj) - dehbdrjk*zkj
                              xregdrv(nregionk) = xregdrv(nregionk) + dehbdrjk*xkj
                              yregdrv(nregionk) = yregdrv(nregionk) + dehbdrjk*ykj
                              zregdrv(nregionk) = zregdrv(nregionk) + dehbdrjk*zkj
                            endif
!
                            if (lstr) then
                              rstrdloc(1:nstrains) = 0.0_dp
                              do kl = 1,nstrains
                                ns = nstrptr(kl)
                                rstrdloc(kl) = rstrdloc(kl) + dehbdrjk*dr2ds(ns,3)
                              enddo
                              do kl = 1,nstrains 
                                rstrd(kl) = rstrd(kl) + rstrdloc(kl)
                              enddo
                              if (latomicstress) then 
                                do kl = 1,nstrains
                                  atomicstress(kl,j) = atomicstress(kl,j) + 0.5_dp*rstrdloc(kl)
                                  atomicstress(kl,k) = atomicstress(kl,k) + 0.5_dp*rstrdloc(kl)
                                enddo
                              endif
                            endif
                          endif
!
!  Derivatives of hydrogen bond energy with respect to bond orders -> add to eval terms
!
                          dehbdBOij = phb1*frHB*(dfHBdBO*(1.0_dp - exphb1) - fHB*dexphb1dBOij)*exphb2*sin4halftheta
!
                          do nl = 1,nonzero_ij
                            d1is(nl) = d1is(nl) + dehbdBOij*d1BOi(nl)
                          enddo
!
                          if (lgrad2) then
!**********************
!  Second derivatives *
!**********************
!
!  Hydrogen bond derivatives - general terms
!
                            d2exphb1dBOij2 = phb2*phb2*exphb1
                            d2exphb2drik2 = phb3*phb3*exphb2*(-r0hb*(rrik**3) + rrik/r0hb)**2 - &
                                            phb3*exphb2*(3.0_dp*r0hb*(rrik**5) - (rrik**3)/r0hb)
                            d2sin4halfthetadtheta2 = 3.0_dp*(coshalftheta*sinhalftheta)**2 - sin4halftheta
!
!  Distance-based second derivatives
!
!  i-j / i-j
!
                            d2ehbdrij2 = fHB*frHB*(phb1*(1.0_dp - exphb1)*exphb2)* &
                                         (dsin4halfthetadtheta*d2thetadr2(1) + d2sin4halfthetadtheta2*dthetadrij**2)
                            call add_drv2_1dm(lilocal,ljlocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf,xji,yji,zji,dehbdrij, &
                                              d2ehbdrij2,dr2ds(1,1),d2r2dx2(1,1),d2r2ds2(1,1,1),d2r2dsdx(1,1,1),lilocal)
!
!  i-j / i-k
!
                            d2ehbdrijdrik = fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2* &
                                            (dsin4halfthetadtheta*d2thetadr2(2) + d2sin4halfthetadtheta2*dthetadrij*dthetadrik) + &
                                            fHB*frHB*phb1*(1.0_dp - exphb1)*dexphb2drik*dsin4halfthetadtheta*dthetadrij
                            call add_drv2_2dm(lilocal,ljlocal,lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                                              ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xji,yji,zji,xki,yki,zki,d2ehbdrijdrik, &
                                              dr2ds(1,1),dr2ds(1,2),lilocal)
!
!  i-j / j-k
!
                            d2ehbdrijdrjk = fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2* &
                                            (dsin4halfthetadtheta*d2thetadr2(3) + d2sin4halfthetadtheta2*dthetadrij*dthetadrjk) + &
                                            fHB*dfrHBdr*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrij
                            call add_drv2_2dm(lilocal,ljlocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,jx,jy,jz,jxf,jyf,jzf, &
                                              jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xji,yji,zji,xkj,ykj,zkj,d2ehbdrijdrjk, &
                                              dr2ds(1,1),dr2ds(1,3),lilocal)
!
!  i-k / i-k
!
                            d2ehbdrik2 = fHB*frHB*phb1*(1.0_dp - exphb1)*d2exphb2drik2*sin4halftheta + &
                                         fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2* &
                                         (dsin4halfthetadtheta*d2thetadr2(4) + d2sin4halfthetadtheta2*dthetadrik**2) + &
                                         2.0_dp*fHB*frHB*(phb1*(1.0_dp - exphb1)*dexphb2drik)*dsin4halfthetadtheta*dthetadrik
                            call add_drv2_1dm(lilocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf,xki,yki,zki,dehbdrik, &
                                              d2ehbdrik2,dr2ds(1,2),d2r2dx2(1,2),d2r2ds2(1,1,2),d2r2dsdx(1,1,2),lilocal)
!
!  i-k / j-k
!
                            d2ehbdrikdrjk = fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2* &
                                            (dsin4halfthetadtheta*d2thetadr2(5) + d2sin4halfthetadtheta2*dthetadrik*dthetadrjk) + &
                                            fHB*frHB*phb1*(1.0_dp - exphb1)*dexphb2drik*dsin4halfthetadtheta*dthetadrjk + &
                                            fHB*dfrHBdr*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrik + &
                                            fHB*dfrHBdr*phb1*(1.0_dp - exphb1)*dexphb2drik*sin4halftheta
                            call add_drv2_2dm(lilocal,lklocal,ljlocal,lklocal,ix,iy,iz,ixf,iyf,izf,kx,ky,kz,kxf,kyf,kzf, &
                                              jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xki,yki,zki,xkj,ykj,zkj,d2ehbdrikdrjk, &
                                              dr2ds(1,2),dr2ds(1,3),lilocal)
!
!  j-k / j-k
!
                            d2ehbdrjk2 = fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2* &
                                         (dsin4halfthetadtheta*d2thetadr2(6) + d2sin4halfthetadtheta2*dthetadrjk**2) + &
                                         fHB*d2frHBdr2*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta + &
                                         2.0_dp*fHB*dfrHBdr*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrjk
                            call add_drv2_1dm(ljlocal,lklocal,jx,jy,jz,jxf,jyf,jzf,kx,ky,kz,kxf,kyf,kzf,xkj,ykj,zkj,dehbdrjk, &
                                              d2ehbdrjk2,dr2ds(1,3),d2r2dx2(1,3),d2r2ds2(1,1,3),d2r2dsdx(1,1,3),lilocal)
!
!  Derivatives of hydrogen bond energy with respect to bond orders -> add to eval terms
!
                            d2ehbdBOij2 = phb1*frHB*(d2fHBdBO2*(1.0_dp - exphb1) - fHB*d2exphb1dBOij2 - &
                                                     2.0_dp*dfHBdBO*dexphb1dBOij)*exphb2*sin4halftheta
!
                            inds = 0
                            do nl = 1,nonzero_ij
                              do nm = 1,nl
                                inds = inds + 1
                                d2is(inds) = d2is(inds) + dehbdBOij*d2BOi(inds)
                                d2is(inds) = d2is(inds) + d2ehbdBOij2*d1BOi(nl)*d1BOi(nm)
                              enddo
                            enddo
!
!  Mixed bond order / distance derivatives
!
!  i-j
!
                            d2trm = dfHBdBO*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrij - &
                                    fHB*frHB*phb1*dexphb1dBOij*exphb2*dsin4halfthetadtheta*dthetadrij
                            do nl = 1,nonzero_ij
                              d2mix(nl) = d2trm*d1BOi(nl)
                            enddo
                            call d2addmixdsp(lilocal,lilocal,ljlocal,i,ix,iy,iz,ixf,iyf,izf,maxneigh,mneigh,nneigh,neighno, &
                                             xneigh,yneigh,zneigh,nbodummyptr,d2mix,nonzero_ij,nonzeroptr_ij,ix,iy,iz,ixf,iyf,izf, &
                                             jx,jy,jz,jxf,jyf,jzf,xji,yji,zji,dr2ds(1,1),.true.)
!
!  i-k
!
                            d2trm = dfHBdBO*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrik + &
                                    dfHBdBO*frHB*phb1*(1.0_dp - exphb1)*dexphb2drik*sin4halftheta - &
                                    fHB*frHB*phb1*dexphb1dBOij*exphb2*dsin4halfthetadtheta*dthetadrik - &
                                    fHB*frHB*phb1*dexphb1dBOij*dexphb2drik*sin4halftheta
                            do nl = 1,nonzero_ij
                              d2mix(nl) = d2trm*d1BOi(nl)
                            enddo
                            call d2addmixdsp(lilocal,lilocal,lklocal,i,ix,iy,iz,ixf,iyf,izf,maxneigh,mneigh,nneigh,neighno, &
                                             xneigh,yneigh,zneigh,nbodummyptr,d2mix,nonzero_ij,nonzeroptr_ij,ix,iy,iz,ixf,iyf,izf, &
                                             kx,ky,kz,kxf,kyf,kzf,xki,yki,zki,dr2ds(1,2),.true.)
!
!  j-k
!
                            d2trm = dfHBdBO*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrjk + &
                                    dfHBdBO*dfrHBdr*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta - &
                                    fHB*frHB*phb1*dexphb1dBOij*exphb2*dsin4halfthetadtheta*dthetadrjk - &
                                    fHB*dfrHBdr*phb1*dexphb1dBOij*exphb2*sin4halftheta
                            do nl = 1,nonzero_ij
                              d2mix(nl) = d2trm*d1BOi(nl)
                            enddo
                            call d2addmixdsp(lilocal,ljlocal,lklocal,i,ix,iy,iz,ixf,iyf,izf,maxneigh,mneigh,nneigh,neighno, &
                                             xneigh,yneigh,zneigh,nbodummyptr,d2mix,nonzero_ij,nonzeroptr_ij,jx,jy,jz,jxf,jyf,jzf, &
                                             kx,ky,kz,kxf,kyf,kzf,xkj,ykj,zkj,dr2ds(1,3),.true.)
                          endif
                        endif
                      enddo nloop
!
!  End check on whether k is ReaxFF species
!
                    endif
!
!  End of loop over k
!
                  enddo kloop
!
!  Add derivatives due to i-j pair
!
                  if (lgrad1) then
                    if (lilocal) then
                      call d1addsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nbodummyptr,d1is, &
                                   nonzero_ij,nonzeroptr_ij,.true.,.false.,.true.)
                    endif
                    if (lgrad2) then
                      call d2adddmsp(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nboatomRptr,d1is,d2is, &
                                     nonzero_ij,nonzeroptr_ij,.true.,.false.)
                    endif
                  endif
!
!  End of check as to whether there are any H-bonding terms for this i-j pair
!
                endif
!
!  End of check on bond order threshold
!
              endif
!
!  End of check on whether j is a ReaxFF atom
!
            endif
          enddo
!
!  End of QMMM condition
!
        endif
      endif
    enddo
  endif   ! End of condition on hydrogen bonding energy being included
!***************************************
!  Twobody VDW & Coulomb contribution  *
!***************************************
  t3 = g_cpu_time()
  if (lreaxFF_ecoul.and.lreaxFFQ) then
    if (literativeQ.and..not.lgrad2) then
      call setreaxffQiter(nboatom,nboatomptr,nboatomRptr,nbos,nbosptr,qreaxFF,eself,.false.)
    else
      call setreaxffQ(nboatom,nboatomRptr,nbos,nbosptr,qreaxFF,eself,.false.,lgrad2)
! 
!  This version of the charge solution is not currently parallel and so the self-energy has to be divided by the number of processors
! 
      eself = eself/dble(nprocs)
    endif
  else
    eself = 0.0_dp
    qreaxFF(1:nboatom) = 0.0_dp
  endif
  t4 = g_cpu_time()
!
!  Set cutoffs
!
  cutvdw2 = reaxFFcutoffVDW**2
  cutq2   = reaxFFcutoffQ**2
  cut2    = max(cutvdw2,cutq2)
!
  if (lreaxff_evdw.or.lreaxff_ecoul) then
!***************************
!  Standard N^2 algorithm  *
!***************************
    derive0self = 0.0_dp
    d2self = 0.0_dp
!         
!  Loop over pairs 
!         
    do ii = 1,natomsonnode
      i = node2atom(ii)
      if (lgrad2) then
        ix = 3*(ii-1) + 1
        iy = ix + 1
        iz = ix + 2   
        ixf = 3*(i-1) + 1
        iyf = ixf + 1
        izf = ixf + 2
      endif
!
!  Does i have a reaxFF species or a fixed charge?
!
      nspeci = nbosptr(i)
      lfixQi = lreaxFFqfix(nspeci)
      lreactivei = (nbos(i).gt.0)
!
      if (lreactivei.or.lfixQi) then
        gammai = reaxFFgamma(nspeci)
        ib = nboatomptr(i)
!
!  Set region for i
!
        nregioni = nregionno(nrelf2a(i))
        nregiontypi = nregiontype(nregioni,ncf)
!
!  Loop over second atom
!
        jloop: do jj = 1,nboatom
          j = nboatomRptr(jj)
          if (lgrad2) then
            jx = 3*(j-1) + 1
            jy = jx + 1
            jz = jx + 2
          endif
!
!  Does j have a reaxFF species or a fixed charge?
!
          nspecj = nbosptr(j)
          lfixQj = lreaxFFqfix(nspecj)
          lreactivej = (nbos(j).gt.0)
!
          if (lreactivej.or.lfixQj) then
            gammaj = reaxFFgamma(nspecj)
!
!  Compute basic interatomic vector
!
            xji = xclat(j) - xclat(i)
            yji = yclat(j) - yclat(i)
            zji = zclat(j) - zclat(i)
!
!  Find valid vectors
!
            nor = 0
            if (ndim.eq.3) then
              call rsearch3D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,cutq2)
            elseif (ndim.eq.2) then
              call rsearch2D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,cutq2)
            elseif (ndim.eq.1) then
              call rsearch1D(xji,yji,zji,.false.,.false.,i,j,0_i4,0_i4,0_i4,0_i4,nor,nmolonly,lself,cut2,cutq2)
            elseif (ndim.eq.0) then
              rij2 = xji*xji + yji*yji + zji*zji
              if (rij2.lt.cut2.and.rij2.gt.1.0d-12) then
                nor = nor + 1
                dist(nor) = rij2
                xtmp(nor) = xji
                ytmp(nor) = yji
                ztmp(nor) = zji
              endif
            endif
!
            if (nor.eq.0) then
!
!  If there are no distances just add the self-term before cycling
!
              call d2chargeselfd(ii,i,j,ix,iy,iz,jx,jy,jz,nbosptr)
              cycle jloop
            endif
!
!  Set index for pairwise parameters
!
            if (nspeci.ge.nspecj) then
              ind = nspeci*(nspeci - 1)/2 + nspecj
            else
              ind = nspecj*(nspecj - 1)/2 + nspeci
            endif
!
!  Set region for j
!
            nregionj = nregionno(nrelf2a(j))
            nregiontypj = nregiontype(nregionj,ncf)
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
            lQMMMok = .true.
            if (QMMMmode(ncf).gt.0) then
              if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
            endif
            if (lQMMMok) then
!
!  QM/MM : Set electrostatic embedding flag : If either i or j are QM atoms => exclude electrostatics
!
              lQMMMelectro = (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1))
!
              qij = qreaxFF(ib)*qreaxFF(jj)*angstoev
!
!  Loop over valid distances and calculate contributions
!
              do n = 1,nor
                if (ndim.gt.0) rij2 = dist(n)
                rij = sqrt(rij2)
!
!  Compute taper function(s)
!
                call p7reaxFFvdwtaper(rij,reaxFFcutoffVDW,tp,dtpdr,d2tpdr2,lgrad1,lgrad2)
                if (reaxFFcutoffQ.ne.reaxFFcutoffVDW) then
                  call p7reaxFFqtaper(rij,reaxFFcutoffQ,tpQ,dtpQdr,d2tpQdr2,lgrad1,lgrad2)
                else
                  tpQ      = tp
                  dtpQdr   = dtpdr
                  d2tpQdr2 = d2tpdr2
                endif
                if (lreaxff_evdw.and.lreactivei.and.lreactivej) then
!
!  Compute f13
!
                  call reaxFF_f13(rij,reaxFFgammaVDW(ind),f13,df13drij,d2f13drij2,lgrad1,lgrad2)
!
!  Energy terms
!
                  VDWtrm  = reaxFFalphaVDW(ind)*(1.0_dp - f13/reaxFFr0VDW(ind))
                  expVDW1 = exp(0.5_dp*VDWtrm)
                  expVDW2 = expVDW1*expVDW1
!
!  Add VDW energy
!
                  etrm = reaxFFDeVDW(ind)*(expVDW2 - 2.0_dp*expVDW1)
                  evdw = evdw + 0.5_dp*etrm*tp
!
!  Site energies for evdw
!
                  siteenergy(i) = siteenergy(i) + 0.5_dp*etrm*tp
                endif
                if (lreaxFF_ecoul.and..not.lQMMMelectro) then
!
!  Compute Coulomb shielding parameters
!
                  gammaij  = sqrt(gammai*gammaj)
                  gammaij  = 1.0_dp/(gammaij**3)
!                     
!  Pure real space tapered lreaxFF case
!                     
                  gam = qij/(rij*rij2 + gammaij)**third
                  ecoul = ecoul + 0.5_dp*gam*tpQ
!
!  Site energies for ecoul
!
                  siteenergy(i) = siteenergy(i) + 0.5_dp*gam*tpQ
                endif
!
                if (lgrad1) then
                  if (lreaxff_evdw.and.lreactivei.and.lreactivej) then
!
!  Compute derivatives of VDW energy
!
                    dVDWtrmdf13 = - reaxFFalphaVDW(ind)/reaxFFr0VDW(ind)
                    dVDWtrmdrij = dVDWtrmdf13*df13drij
                    detrmdrij = reaxFFDeVDW(ind)*(expVDW2 - expVDW1)*dVDWtrmdrij
                    devdwdrij = detrmdrij*tp + etrm*dtpdr
                    if (lgrad2) then
                      d2VDWtrmdrij2 = dVDWtrmdf13*d2f13drij2
                      d2etrmdrij2 = reaxFFDeVDW(ind)*((expVDW2 - 0.5_dp*expVDW1)*dVDWtrmdrij**2 + &
                                                            (expVDW2 - expVDW1)*d2VDWtrmdrij2)
                      d2evdwdrij2 = d2etrmdrij2*tp + 2.0_dp*detrmdrij*dtpdr + etrm*d2tpdr2
                    endif
                  else
                    devdwdrij = 0.0_dp
                    d2evdwdrij2 = 0.0_dp
                  endif
                  if (lreaxFF_ecoul.and..not.lQMMMelectro) then
!
!  Compute derivative of Coulomb energy
!           
                    dgamdrij = - (gam*rij)/(rij*rij2 + gammaij)
                    decouldrij = gam*dtpQdr + dgamdrij*tpQ
                    if (lgrad2) then
                      d2gamdrij2 = dgamdrij*(1.0_dp/rij2 - 4.0_dp*rij/(rij*rij2 + gammaij))
                      d2ecouldrij2 = d2gamdrij2*tpQ + 2.0_dp*dgamdrij*dtpQdr + gam*d2tpQdr2
!
!  Mixed charge distance derivatives
!
                      gamnoq = angstoev/(rij*rij2 + gammaij)**third
                      dgamnoqdrij = - (gamnoq*rij)/(rij*rij2 + gammaij)
                      d0i(n) = gamnoq*tpQ*qreaxFF(jj)
                      d0j(n) = gamnoq*tpQ*qreaxFF(ib)
                      d2i2(n) = 0.0_dp
                      if (i.eq.j) then
                        d2ij(n) = 0.5_dp*gamnoq*tpQ
                      else
                        d2ij(n) = gamnoq*tpQ
                      endif
                      d2j2(n) = 0.0_dp
                      d1qi(n) = (dgamnoqdrij*tpQ + gamnoq*dtpQdr)*qreaxFF(jj)
                      d1qj(n) = (dgamnoqdrij*tpQ + gamnoq*dtpQdr)*qreaxFF(ib)
                    endif
                  else
                    decouldrij = 0.0_dp
                    d2ecouldrij2 = 0.0_dp
                    d0i(n) = 0.0_dp
                    d0j(n) = 0.0_dp
                    d2i2(n) = 0.0_dp
                    d2ij(n) = 0.0_dp
                    d2j2(n) = 0.0_dp
                    d1qi(n) = 0.0_dp
                    d1qj(n) = 0.0_dp
                  endif
!
!  Get Cartesian components of vector and derivative
!
                  if (ndim.gt.0) then
                    xji = xtmp(n)
                    yji = ytmp(n)
                    zji = ztmp(n)
                  endif
                  xd = xji*(devdwdrij + decouldrij)
                  yd = yji*(devdwdrij + decouldrij)
                  zd = zji*(devdwdrij + decouldrij)
!
!  Add derivatives to main arrays
!
                  xdrv(i) = xdrv(i) - xd
                  ydrv(i) = ydrv(i) - yd
                  zdrv(i) = zdrv(i) - zd
!
                  if (nregioni.ne.nregionj) then
                    xregdrv(nregioni) = xregdrv(nregioni) - xd
                    yregdrv(nregioni) = yregdrv(nregioni) - yd
                    zregdrv(nregioni) = zregdrv(nregioni) - zd
                  endif
!
                  if (lstr.or.lgrad2) then
                    call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,1), &
                      d2r2dx2(1,1),d2r2dsdx(1,1,1),d2r2ds2(1,1,1),lgrad2)
                  endif
                  if (lstr) then
                    rstrdloc(1:nstrains) = 0.0_dp
                    do kl = 1,nstrains
                      ns = nstrptr(kl)
                      rstrdloc(kl) = rstrdloc(kl) + (devdwdrij + decouldrij)*dr2ds(ns,1)
                    enddo
                    do kl = 1,nstrains 
                      rstrd(kl) = rstrd(kl) + 0.5_dp*rstrdloc(kl)
                    enddo
                    if (latomicstress) then 
                      do kl = 1,nstrains
                        atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*rstrdloc(kl)
                      enddo
                    endif
                  endif
                  if (lgrad2) then
!
!  Second derivatives for VDW and fixed charge part of Coulomb terms
!
                    d1trm = devdwdrij + decouldrij
                    d2trm = d2evdwdrij2 + d2ecouldrij2
!
                    derv2(jx,ix) = derv2(jx,ix) - d2trm*d2r2dx2(1,1)
                    derv2(jy,ix) = derv2(jy,ix) - d2trm*d2r2dx2(6,1)
                    derv2(jz,ix) = derv2(jz,ix) - d2trm*d2r2dx2(5,1)
                    derv2(jx,iy) = derv2(jx,iy) - d2trm*d2r2dx2(6,1)
                    derv2(jy,iy) = derv2(jy,iy) - d2trm*d2r2dx2(2,1)
                    derv2(jz,iy) = derv2(jz,iy) - d2trm*d2r2dx2(4,1)
                    derv2(jx,iz) = derv2(jx,iz) - d2trm*d2r2dx2(5,1)
                    derv2(jy,iz) = derv2(jy,iz) - d2trm*d2r2dx2(4,1)
                    derv2(jz,iz) = derv2(jz,iz) - d2trm*d2r2dx2(3,1)
                    derv2(jx,ix) = derv2(jx,ix) - d1trm
                    derv2(jy,iy) = derv2(jy,iy) - d1trm
                    derv2(jz,iz) = derv2(jz,iz) - d1trm
!
                    if (lstr) then
                      do kl = 1,nstrains
                        ks = nstrptr(kl)
                        derv3(ix,kl) = derv3(ix,kl) - d1trm*d2r2dsdx(ks,1,1) - xji*d2trm*dr2ds(ks,1)
                        derv3(iy,kl) = derv3(iy,kl) - d1trm*d2r2dsdx(ks,2,1) - yji*d2trm*dr2ds(ks,1)
                        derv3(iz,kl) = derv3(iz,kl) - d1trm*d2r2dsdx(ks,3,1) - zji*d2trm*dr2ds(ks,1)
                      enddo
!
                      do kk = 1,nstrains
                        ks = nstrptr(kk)
                        do kl = 1,nstrains
                          kt = nstrptr(kl)
                          sderv2(kl,kk) = sderv2(kl,kk) + 0.5_dp*d2trm*dr2ds(kt,1)*dr2ds(ks,1)
                          sderv2(kl,kk) = sderv2(kl,kk) + 0.5_dp*d1trm*d2r2ds2(kt,ks,1)
                        enddo
                      enddo
                    endif
                  endif
                endif
!
!  End loop over distances
!
              enddo
            endif
          endif
!*****************************************
!  Contribution from charge derivatives  *
!*****************************************
          if (lgrad2.and.lreaxFFQ) then
            derive0selfi = derive0self*qreaxFF(jj)
            derive0selfj = derive0self*qreaxFF(ib)
            d1ix = 0.0_dp
            d1iy = 0.0_dp
            d1iz = 0.0_dp
            d1jx = 0.0_dp
            d1jy = 0.0_dp
            d1jz = 0.0_dp
            do k = 1,nor
              d1ix = d1ix + d1qi(k)*xtmp(k)
              d1iy = d1iy + d1qi(k)*ytmp(k)
              d1iz = d1iz + d1qi(k)*ztmp(k)
              d1jx = d1jx + d1qj(k)*xtmp(k)
              d1jy = d1jy + d1qj(k)*ytmp(k)
              d1jz = d1jz + d1qj(k)*ztmp(k)
            enddo
            if (lstr) then
              do kl = 1,nstrains
                ns = nstrptr(kl)
                d1si(kl) = 0.0_dp
                d1sj(kl) = 0.0_dp
                do k = 1,nor
                  call real1strterm(ndim,xtmp(k),ytmp(k),ztmp(k),0.0_dp,0.0_dp,0.0_dp,dr2ds(1,1), &
                    d2r2dx2(1,1),d2r2dsdx(1,1,1),d2r2ds2(1,1,1),.false.)
                  d1si(kl) = d1si(kl) + d1qi(k)*dr2ds(ns,1)
                  d1sj(kl) = d1sj(kl) + d1qj(k)*dr2ds(ns,1)
                enddo
              enddo
              if (i.eq.j) then
                do kl = 1,nstrains
                  d1si(kl) = 0.5_dp*d1si(kl)
                  d1sj(kl) = 0.5_dp*d1sj(kl)
                enddo
              endif
            endif
            call d2charged(ii,i,j,nor,ix,iy,iz,jx,jy,jz,d0i,d0j,d1ix,d1iy,d1iz, &
                           d1jx,d1jy,d1jz,d1si,d1sj,d2i2,d2ij,d2j2,d2self,derive0selfi, &
                           derive0selfj,.true.,.true.,nbosptr)
          endif
!
!  End of loop over j
!
        enddo jloop
      endif
!
!  End of loop over i
!
    enddo
  endif
!
!  Sum energy contributions
!
  ereaxFF = ebond + ebpen + econj + elone + eover + epen + ecoa + etors + eunder + eval + ehb + evdw + ecoul + eself
  if (index(keyword,'verb').ne.0)  then
!
!  Global sum of energy components for print out purposes
!
    sum0(1)  = ebond
    sum0(2)  = ebpen
    sum0(3)  = elone
    sum0(4)  = eover
    sum0(5)  = eunder
    sum0(6)  = eval
    sum0(7)  = epen
    sum0(8)  = ecoa
    sum0(9)  = etors
    sum0(10) = econj
    sum0(11) = ehb
    sum0(12) = evdw
    sum0(13) = ecoul
    sum0(14) = eself
    call sumall(sum0,sum,14_i4,"reaxFFds","energies")
    ebond  = sum(1)
    ebpen  = sum(2)
    elone  = sum(3)
    eover  = sum(4)
    eunder = sum(5)
    eval   = sum(6)
    epen   = sum(7)
    ecoa   = sum(8)
    etors  = sum(9)
    econj  = sum(10)
    ehb    = sum(11)
    evdw   = sum(12)
    ecoul  = sum(13)
    eself  = sum(14)
  endif
!
  if (index(keyword,'verb').ne.0.and.ioproc)  then
!
!  Print out contributions to energy
!
    write(ioout,'(''  ReaxFF : Energy contributions: '',/)')
    write(ioout,'(''  E(bond)     = '',f16.8,'' eV = '',f16.7,'' kcal'')') ebond,ebond*evtokcal
    write(ioout,'(''  E(bpen)     = '',f16.8,'' eV = '',f16.7,'' kcal'')') ebpen,ebpen*evtokcal
    write(ioout,'(''  E(lonepair) = '',f16.8,'' eV = '',f16.7,'' kcal'')') elone,elone*evtokcal
    write(ioout,'(''  E(over)     = '',f16.8,'' eV = '',f16.7,'' kcal'')') eover,eover*evtokcal
    write(ioout,'(''  E(under)    = '',f16.8,'' eV = '',f16.7,'' kcal'')') eunder,eunder*evtokcal
    write(ioout,'(''  E(val)      = '',f16.8,'' eV = '',f16.7,'' kcal'')') eval,eval*evtokcal
    write(ioout,'(''  E(pen)      = '',f16.8,'' eV = '',f16.7,'' kcal'')') epen,epen*evtokcal
    write(ioout,'(''  E(coa)      = '',f16.8,'' eV = '',f16.7,'' kcal'')') ecoa,ecoa*evtokcal
    write(ioout,'(''  E(tors)     = '',f16.8,'' eV = '',f16.7,'' kcal'')') etors,etors*evtokcal
    write(ioout,'(''  E(conj)     = '',f16.8,'' eV = '',f16.7,'' kcal'')') econj,econj*evtokcal
    write(ioout,'(''  E(hb)       = '',f16.8,'' eV = '',f16.7,'' kcal'')') ehb,ehb*evtokcal
    write(ioout,'(''  E(vdw)      = '',f16.8,'' eV = '',f16.7,'' kcal'')') evdw,evdw*evtokcal
    write(ioout,'(''  E(coulomb)  = '',f16.8,'' eV = '',f16.7,'' kcal'')') ecoul,ecoul*evtokcal
    write(ioout,'(''  E(self)     = '',f16.8,'' eV = '',f16.7,'' kcal'')') eself,eself*evtokcal
  endif
!
!  Output of charge / bond orders
!
  if (lqbo) then
    call outqbo(1_i4,numat,qreaxFF,maxneigh,nneigh,neighno,BO)
  endif
!
!  Free local memory
!
  deallocate(d2sboddeltaidbo,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2sboddeltaidbo')
  deallocate(d2sboproddbo2,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2sboproddbo2')
  deallocate(nonzeroptr_j2k2,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_j2k2')
  deallocate(nonzeroptr_jk,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_jk')
  deallocate(nonzeroptr_ij2,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_ij2')
  deallocate(nonzeroptr_ik,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_ik')
  deallocate(nonzeroptr_ij,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nonzeroptr_ij')
  deallocate(nonzeroptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nonzeroptr')
  deallocate(nonzerorptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nonzerorptr')
  deallocate(d1BOjk2_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOjk2_pipi')
  deallocate(d1BOjk2_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOjk2_pi')
  deallocate(d1BOjk2_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOjk2_s')
  deallocate(d1BOjk2,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOjk2')
  deallocate(d1BOjk_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOjk_pipi')
  deallocate(d1BOjk_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOjk_pi')
  deallocate(d1BOjk_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOjk_s')
  deallocate(d1BOjk,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOjk')
  deallocate(d2BOp_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOp_pipi')
  deallocate(d2BOp_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOp_pi')
  deallocate(d2BOp,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOp')
  deallocate(d2BOjk_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOjk_pipi')
  deallocate(d2BOjk_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOjk_pi')
  deallocate(d2BOjk_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOjk_s')
  deallocate(d2BOjk,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOjk')
  deallocate(d2BOi_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOi_pipi')
  deallocate(d2BOi_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOi_pi')
  deallocate(d2BOi_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOi_s')
  deallocate(d2BOi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2BOi')
  deallocate(d2mix,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2mix')
!
  deallocate(dsboproddbo,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','dsboproddbo')
  deallocate(d1BOik_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOik_pipi')
  deallocate(d1BOik_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOik_pi')
  deallocate(d1BOik_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOik_s')
  deallocate(d1BOik,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOik')
  deallocate(d1BOij2_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOij2_pipi')
  deallocate(d1BOij2_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOij2_pi')
  deallocate(d1BOij2_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOij2_s')
  deallocate(d1BOij2,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOij2')
  deallocate(d1BOi_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOi_pipi')
  deallocate(d1BOi_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOi_pi')
  deallocate(d1BOi_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOi_s')
  deallocate(d1BOi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOi')
  deallocate(devaldsbo_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','devaldsbo_neigh')
  deallocate(devalddeltaj_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','devalddeltaj_neigh')
  deallocate(d2etorsdthetaddeltai,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2etorsdthetaddeltai')
  deallocate(d2etorsdBOdBOpi_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2etorsdBOdBOpi_neigh')
  deallocate(d2etorsdBOpi2_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2etorsdBOpi2_neigh')
  deallocate(detorsdBOpi_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','detorsdBOpi_neigh')
  deallocate(d2evaldthetaddeltai_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2evaldthetaddeltai_neigh')
  deallocate(d2evaldthetadBOpi_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2evaldthetadBOpi_neigh')
  deallocate(d2evalddeltajdBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2evalddeltajdBO_neigh')
  deallocate(d2etorsddeltajdBOpi_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2etorsddeltajdBOpi_neigh')
  deallocate(d2etorsddeltaidBOpi_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2etorsddeltaidBOpi_neigh')
  deallocate(d2etorsddeltadBOpi_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2etorsddeltadBOpi_neigh')
  deallocate(d2evaldsbo2_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2evaldsbo2_neigh')
  deallocate(d2evaldsbo_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2evaldsbo_neigh')
  deallocate(d2evalddeltajk_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2evalddeltajk_neigh')
  deallocate(d2evalddeltaiddeltaj_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2evalddeltaiddeltaj_neigh')
  deallocate(d2evalddeltaidBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2evalddeltaidBO_neigh')
  deallocate(d2evaldBO2_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2evaldBO2_neigh')
  deallocate(devaldBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','devaldBO_neigh')
  deallocate(d2ebpendBO2_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2ebpendBO2_neigh')
  deallocate(d2ebpenddeltaidBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d2ebpenddeltaidBO_neigh')
  deallocate(debpendBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','debpendBO_neigh')
  deallocate(d1BOp_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOp_pipi')
  deallocate(d1BOp_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOp_pi')
  deallocate(d1BOp,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','d1BOp')
  deallocate(deltalplp,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','deltalplp')
  deallocate(deltalp,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','deltalp')
  deallocate(deltap,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','deltap')
  deallocate(delta,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','delta')
  deallocate(BOp_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','BOp_pipi')
  deallocate(BOp_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','BOp_pi')
  deallocate(BOp,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','BOp')
  deallocate(BO_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','BO_pipi')
  deallocate(BO_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','BO_pi')
  deallocate(BO,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','BO')
  deallocate(zneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','zneigh')
  deallocate(yneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','yneigh')
  deallocate(xneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','xneigh')
  deallocate(rneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','rneigh')
  deallocate(ineigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','ineigh')
  deallocate(neighnoRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','neighnoRptr')
  deallocate(neighno,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','neighno')
  deallocate(sum0,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','sum0')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','sum')
  deallocate(sumsij,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','sumsij')
  deallocate(nneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nneigh')
  deallocate(nfreeatom,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nfreeatom')
  deallocate(lopanyneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','lopanyneigh')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','latomdone')
  deallocate(nbosptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nbosptr')
  deallocate(nbos,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nbos')
  deallocate(nbodummyptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nbodummyptr')
  deallocate(nboatomRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nboatomRptr')
  deallocate(nboatomptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','nboatomptr')
  deallocate(numattodoRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','numattodoRptr')
  deallocate(numattodonptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','numattodonptr')
  deallocate(numattodoptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFds','numattodoptr')
!     
  if (lgrad1) then
    call realloc(d1is,0_i4,status)
    if (status.ne.0) call outofmemory('reaxFFds','d1is')
    if (lgrad2) then
      call realloc(d2is,0_i4,status)
      if (status.ne.0) call outofmemory('reaxFFds','d2is')
      call realloc(d2evaldthetadBO_neigh,0_i4,0_i4,status)
      if (status.ne.0) call outofmemory('reaxFFds','d2evaldthetadBO_neigh')
      call realloc(d2etorsdthetaddeltaj_neigh,0_i4,0_i4,status)
      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetaddeltaj_neigh')
      call realloc(d2etorsdthetadBOpi_neigh,0_i4,0_i4,status)
      if (status.ne.0) call outofmemory('reaxFFds','d2etorsdthetadBOpi_neigh')
    endif
  endif
!
  t2 = g_cpu_time()
  treaxFF = treaxFF + t2 - t1 - (t4 - t3)
  teem    = teem + t4 - t3
#ifdef TRACE
  call trace_out('reaxffds')
#endif
!
  return
  end
