  subroutine reaxFFmds(ereaxFF,lgrad1)
!
!  Calculates the energy and derivatives for the reaxFF force field up to first derivatives.
!  Uses sparsity to speed up calculation.
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
!
!  On exit :
!
!  ereaxFF         = the value of the energy contribution
!
!   5/23 Created from reaxFFs.F90
!   6/23 Electrostatic cutoff squared now passed to rsearch subroutines
!   7/23 nregionnocfg added
!   8/23 Overflow trap added for expuc3
!   8/23 Electric fields added
!   9/23 Scale corrected for ebond in parallel
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
!  Julian Gale, CIC, Curtin University, September 2023
!
  use datatypes
  use configurations, only : nregions, lsliceatom, nregiontype, lregionrigid, QMMMmode
  use control,        only : keyword, lseok, lreaxFFQ, latomicstress, literativeQ
  use g_constants,    only : evtokcal
  use current
  use derivatives,    only : xdrv, ydrv, zdrv, rstrd, xregdrv, yregdrv, zregdrv, atomicstress
  use energies,       only : eattach, esregion12, esregion2, siteenergy
  use gulp_files,     only : lqbo
  use iochannels
  use m_strain,       only : real1strterm
  use neighbours
  use numbers,        only : third
  use optimisation,   only : lfreeze, lopf
  use parallel
  use reaxFFdata
  use realvectors,    only : dist, xtmp, ytmp, ztmp
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
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ib
  integer(i4)                                      :: ic
  integer(i4)                                      :: ii
  integer(i4)                                      :: imin
  integer(i4)                                      :: imx
  integer(i4)                                      :: imy
  integer(i4)                                      :: imz
  integer(i4)                                      :: ind
  integer(i4)                                      :: ind0
  integer(i4)                                      :: indj
  integer(i4)                                      :: indn
  integer(i4)                                      :: ind1
  integer(i4)                                      :: ind1i
  integer(i4)                                      :: ind2
  integer(i4)                                      :: indbojk
  integer(i4)                                      :: indbojl
  integer(i4)                                      :: indil
  integer(i4)                                      :: indjk
  integer(i4)                                      :: indjl
  integer(i4)                                      :: indkl
  integer(i4)                                      :: indsjk
  integer(i4)                                      :: indcosphi(6)
  integer(i4)                                      :: indsinkij(3)
  integer(i4)                                      :: indsinijl(3)
  integer(i4), dimension(:,:,:), allocatable, save :: ineigh         ! Cell indices for vector
  integer(i4)                                      :: isx
  integer(i4)                                      :: isy
  integer(i4)                                      :: isz
  integer(i4)                                      :: itmp
  integer(i4)                                      :: ix
  integer(i4)                                      :: ixyz
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jb
  integer(i4)                                      :: jc
  integer(i4)                                      :: jj
  integer(i4)                                      :: jmax
  integer(i4)                                      :: k
  integer(i4)                                      :: kc
  integer(i4)                                      :: kk
  integer(i4)                                      :: kl
  integer(i4)                                      :: l
  integer(i4)                                      :: m
  integer(i4)                                      :: maxind
  integer(i4)                                      :: maxx
  integer(i4)                                      :: maxxy
  integer(i4)                                      :: maxneigh1
  integer(i4)                                      :: maxneigh1a
  integer(i4)                                      :: maxneigh2
  integer(i4)                                      :: maxneigh2a
  integer(i4)                                      :: maxneigh2e     ! Extended dimension
  integer(i4)                                      :: maxneigh2t     ! Up to torsion dimension
  integer(i4)                                      :: mneigh
  integer(i4)                                      :: n
  integer(i4)                                      :: n1i
  integer(i4)                                      :: n1j
  integer(i4)                                      :: n1k
  integer(i4)                                      :: nati
  integer(i4)                                      :: nboij
  integer(i4)                                      :: nbojk
  integer(i4), dimension(:),     allocatable, save :: nbos
  integer(i4), dimension(:),     allocatable, save :: nbosptr
  integer(i4)                                      :: nboatom
  integer(i4), dimension(:),     allocatable, save :: nboatomptr
  integer(i4), dimension(:),     allocatable, save :: nboatomRptr
  integer(i4), dimension(:),     allocatable, save :: nbodummyptr
  integer(i4)                                      :: ncellsearchHB(3)
  integer(i4)                                      :: ni
  integer(i4)                                      :: nivalid
  integer(i4)                                      :: nj
  integer(i4)                                      :: nk
  integer(i4)                                      :: nl
  integer(i4)                                      :: nmin
  integer(i4)                                      :: nmolonly
  integer(i4)                                      :: nn
  integer(i4)                                      :: nn2
  integer(i4)                                      :: nor
  integer(i4)                                      :: nptr
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
  integer(i4)                                      :: nneighi2e2
  integer(i4)                                      :: nneighi2t2
  integer(i4)                                      :: nneighi22
  integer(i4)                                      :: nneighj2
  integer(i4)                                      :: nneighj2t
  integer(i4)                                      :: nneighj2t2
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregionj
  integer(i4)                                      :: nregionk
  integer(i4)                                      :: nregiontypi
  integer(i4)                                      :: nregiontypj
  integer(i4)                                      :: ns
  integer(i4)                                      :: nspeci
  integer(i4)                                      :: nspecj
  integer(i4)                                      :: nspeck
  integer(i4)                                      :: nspecl
  integer(i4)                                      :: nsplower(3)
  integer(i4)                                      :: nspupper(3)
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: numattodo
  integer(i4)                                      :: numattodo_main
  integer(i4)                                      :: numattodo_main2
  integer(i4), dimension(:),     allocatable, save :: numattodoptr
  integer(i4), dimension(:),     allocatable, save :: numattodonptr
  integer(i4), dimension(:),     allocatable, save :: numattodoRptr
  integer(i4)                                      :: nval3
  integer(i4)                                      :: nonzero_ij
  integer(i4)                                      :: nonzero_ij_i
  integer(i4)                                      :: nonzero_jk
  integer(i4)                                      :: nonzero_jk_j
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_ij
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_ik
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_ij2
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_jk
  integer(i4), dimension(:),     allocatable, save :: nonzeroptr_j2k2
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
  logical                                          :: lmaxneighok
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
  real(dp)                                         :: bo8
  real(dp)                                         :: BOij
  real(dp)                                         :: BOij_s
  real(dp)                                         :: BOij_pi
  real(dp)                                         :: BOij_pipi
  real(dp)                                         :: BOijpbe2
  real(dp)                                         :: BOik
  real(dp)                                         :: BOjk_s
  real(dp)                                         :: BOjk_pi
  real(dp)                                         :: BOjk_pipi
  real(dp)                                         :: BOjl
  real(dp)                                         :: BOpij
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
  real(dp)                                         :: d2self
  real(dp)                                         :: dcos2phidcosphi
  real(dp)                                         :: dcos3phidcosphi
  real(dp)                                         :: delta_angle
  real(dp)                                         :: ddeltaeddeltai
  real(dp)                                         :: ddeltalpddeltai
  real(dp)                                         :: ddeltalpddeltaj
  real(dp)                                         :: ddeltalplpddeltai
  real(dp)                                         :: ddeltalpcorrdbopi
  real(dp)                                         :: ddeltalpcorrddeltai
  real(dp)                                         :: ddeltalpcorrddeltaj
  real(dp)                                         :: derive0self
  real(dp)                                         :: dehbdBOij
  real(dp)                                         :: dehbdrij
  real(dp)                                         :: dehbdrik
  real(dp)                                         :: dehbdrjk
  real(dp)                                         :: debpenddeltai
  real(dp)                                         :: deloneddeltai
  real(dp)                                         :: deloneddeltalplp
  real(dp)                                         :: deoverddeltalpcorr
  real(dp)                                         :: deconjdcosphi
  real(dp)                                         :: deconjdf12
  real(dp)                                         :: deconjdsinkij
  real(dp)                                         :: deconjdsinijl
  real(dp)                                         :: detorsconjdf
  real(dp)                                         :: detorsdcosphi
  real(dp)                                         :: detorsdcos1phi
  real(dp)                                         :: detorsdcos2phi
  real(dp)                                         :: detorsdcos3phi
  real(dp)                                         :: detorsdexpV2
  real(dp)                                         :: detorsdf10
  real(dp)                                         :: detorsdsinkij
  real(dp)                                         :: detorsdsinijl
  real(dp)                                         :: deunderddeltalpcorr
  real(dp)                                         :: deunderdsumdeltabopi
  real(dp)                                         :: devalddeltai
  real(dp)                                         :: devaldsbo
  real(dp)                                         :: devdwdrij
  real(dp)                                         :: diffBOdelta
  real(dp)                                         :: dtheta0dsbo
  real(dp)                                         :: dthetadrij
  real(dp)                                         :: dthetadrik
  real(dp)                                         :: dthetadrjk
  real(dp)                                         :: devaldtheta
  real(dp)                                         :: devaldtheta0
  real(dp)                                         :: dsboddeltai
  real(dp)                                         :: dr2ds(6,3)
  real(dp)                                         :: d2r2dx2(6,3)
  real(dp)                                         :: d2r2ds2(6,6,3)
  real(dp)                                         :: d2r2dsdx(6,3,3)
  real(dp)                                         :: dexpuc1ddeltalpcorr
  real(dp)                                         :: dutrm1ddeltalpcorr
  real(dp)                                         :: dutrm2dsumdeltabopi
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
  real(dp),    dimension(:,:),   allocatable, save :: BO
  real(dp),    dimension(:,:),   allocatable, save :: BO_pi
  real(dp),    dimension(:,:),   allocatable, save :: BO_pipi
  real(dp),    dimension(:,:),   allocatable, save :: BOp
  real(dp),    dimension(:,:),   allocatable, save :: BOp_pi
  real(dp),    dimension(:,:),   allocatable, save :: BOp_pipi
  real(dp),    dimension(:,:),   allocatable, save :: d1BOp
  real(dp),    dimension(:,:),   allocatable, save :: d1BOp_pi
  real(dp),    dimension(:,:),   allocatable, save :: d1BOp_pipi
  real(dp),    dimension(:),     allocatable, save :: d1i
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
  real(dp),    dimension(:),     allocatable, save :: delta
  real(dp),    dimension(:),     allocatable, save :: deltap
  real(dp),    dimension(:),     allocatable, save :: deltalplp ! delta lp for lone pairs
  real(dp),    dimension(:),     allocatable, save :: deltalp   ! delta lp for over coordination
  real(dp),    dimension(:),     allocatable, save :: debpendBO_neigh
  real(dp),    dimension(:),     allocatable, save :: devaldBO_neigh
  real(dp),    dimension(:),     allocatable, save :: detorsdBOpi_neigh
  real(dp),    dimension(:),     allocatable, save :: devalddeltaj_neigh
  real(dp),    dimension(:),     allocatable, save :: devaldsbo_neigh
#ifdef MPI
  real(dp),    dimension(:),     allocatable, save :: botmp
#endif
  real(dp)                                         :: dBOpijdr
  real(dp)                                         :: d2BOpijdr2
  real(dp)                                         :: deltapi
  real(dp)                                         :: deltapj
  real(dp)                                         :: deltapk
  real(dp)                                         :: delta_coa
  real(dp)                                         :: delta_e
  real(dp)                                         :: delta_e0
  real(dp)                                         :: deltalpcorr_i
  real(dp)                                         :: deijdBO_s
  real(dp)                                         :: deijdBO_pi
  real(dp)                                         :: deijdBO_pipi
  real(dp)                                         :: dexpV2dBO_pi
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
  real(dp)                                         :: d2f11ddeltaij2
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
  real(dp)                                         :: deoverdsumover
  real(dp)                                         :: detotdbo_tot
  real(dp)                                         :: detotdbopi
  real(dp)                                         :: detotdbopi_tot
  real(dp)                                         :: detotddeltaj
  real(dp)                                         :: detotddeltaj_tot
  real(dp)                                         :: detotddeltalpcorr
  real(dp)                                         :: dsumdeltabopidbopi
  real(dp)                                         :: dsumdeltabopiddeltaj
  real(dp)                                         :: eij
  real(dp)                                         :: expbo8
  real(dp)                                         :: expcoa2
  real(dp)                                         :: expcoa3a
  real(dp)                                         :: d1expcoa3adBOij
  real(dp)                                         :: d1expcoa3addeltaj
  real(dp)                                         :: expcoa3b
  real(dp)                                         :: d1expcoa3bdBOik
  real(dp)                                         :: d1expcoa3bddeltak
  real(dp)                                         :: expcoa4a
  real(dp)                                         :: d1expcoa4a
  real(dp)                                         :: expcoa4b
  real(dp)                                         :: d1expcoa4b
  real(dp)                                         :: expcoaprod
  real(dp)                                         :: expcoatrm
  real(dp)                                         :: dexpcoatrmddeltai
  real(dp)                                         :: expdi
  real(dp)                                         :: exphb1
  real(dp)                                         :: dexphb1dBOij
  real(dp)                                         :: exphb2
  real(dp)                                         :: dexphb2drik
  real(dp)                                         :: expij
  real(dp)                                         :: explp
  real(dp)                                         :: explp2
  real(dp)                                         :: expoc
  real(dp)                                         :: exppen1
  real(dp)                                         :: d1exppen1
  real(dp)                                         :: exppen2
  real(dp)                                         :: d1exppen2
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
  real(dp)                                         :: detrmdrij
  real(dp)                                         :: dgamdrij
  real(dp)                                         :: gam
  real(dp)                                         :: gammai
  real(dp)                                         :: gammaj
  real(dp)                                         :: gammaij
  real(dp)                                         :: otrm1
  real(dp)                                         :: otrm2
  real(dp)                                         :: dotrm2ddeltalpcorr
  real(dp)                                         :: otrm3
  real(dp)                                         :: dotrm3ddeltalpcorr
  real(dp)                                         :: ddeltalpcorrdsumdeltabopi
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
  real(dp)                                         :: rtmp
  real(dp)                                         :: sbo
  real(dp)                                         :: sbopi
  real(dp)                                         :: sboprod
  real(dp)                                         :: scale
  real(dp)                                         :: sinkij
  real(dp)                                         :: sinijl
  real(dp)                                         :: sinhalftheta
  real(dp)                                         :: sin4halftheta
  real(dp)                                         :: dsin4halfthetadtheta
  real(dp)                                         :: smoothH
  real(dp)                                         :: sumBOcoa_ij
  real(dp)                                         :: sumBOcoa_ik
  real(dp)                                         :: sumdeltabopi
  real(dp)                                         :: sumover
  real(dp)                                         :: dsumoverdbo
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
  real(dp)                                         :: xd
  real(dp)                                         :: yd
  real(dp)                                         :: zd
  real(dp)                                         :: xdiff
  real(dp)                                         :: ydiff
  real(dp)                                         :: zdiff
  real(dp)                                         :: xi
  real(dp)                                         :: yi
  real(dp)                                         :: zi
  real(dp)                                         :: xil
  real(dp)                                         :: yil
  real(dp)                                         :: zil
  real(dp)                                         :: xkl
  real(dp)                                         :: ykl
  real(dp)                                         :: zkl
  real(dp)                                         :: xj
  real(dp)                                         :: yj
  real(dp)                                         :: zj
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
#ifdef TRACE
  call trace_in('reaxffmds')
#endif
!
  t1 = g_cpu_time()
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
  if (status/=0) call outofmemory('reaxFFmds','nboatomptr')
  allocate(nboatomRptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','nboatomRptr')
  allocate(nbodummyptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','nbodummyptr')
!
!  Allocate memory that does not depend on maxneigh
!
  allocate(numattodoptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','numattodoptr')
  allocate(numattodonptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','numattodonptr')
  allocate(numattodoRptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','numattodoRptr')
  allocate(nbos(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','nbos')
  allocate(nbosptr(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','nbosptr')
  allocate(latomdone(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','latomdone')
  allocate(lopanyneigh(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','lopanyneigh')
  allocate(nfreeatom(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','nfreeatom')
  allocate(nneigh(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','nneigh')
  allocate(sumsij(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','sumsij')
  allocate(sum(max(3*numat,14_i4)),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','sum')
  allocate(sum0(max(3*numat,14_i4)),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','sum0')
!
!  Reinitialisation point should maxneigh be increased
!
100 continue
  lmaxneighok = .true.
  if (allocated(rneigh)) then
    deallocate(dsboproddbo,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','dsboproddbo')
    deallocate(nonzeroptr_j2k2,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_j2k2')
    deallocate(nonzeroptr_jk,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_jk')
    deallocate(nonzeroptr_ij2,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_ij2')
    deallocate(nonzeroptr_ik,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_ik')
    deallocate(nonzeroptr_ij,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_ij')
    deallocate(nonzeroptr,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr')
    deallocate(d1BOjk2_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOjk2_pipi')
    deallocate(d1BOjk2_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOjk2_pi')
    deallocate(d1BOjk2_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOjk2_s')
    deallocate(d1BOjk2,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOjk2')
    deallocate(d1BOjk_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOjk_pipi')
    deallocate(d1BOjk_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOjk_pi')
    deallocate(d1BOjk_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOjk_s')
    deallocate(d1BOjk,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOjk')
    deallocate(d1BOik_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOik_pipi')
    deallocate(d1BOik_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOik_pi')
    deallocate(d1BOik_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOik_s')
    deallocate(d1BOik,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOik')
    deallocate(d1BOij2_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOij2_pipi')
    deallocate(d1BOij2_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOij2_pi')
    deallocate(d1BOij2_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOij2_s')
    deallocate(d1BOij2,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOij2')
    deallocate(d1BOi_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOi_pipi')
    deallocate(d1BOi_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOi_pi')
    deallocate(d1BOi_s,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOi_s')
    deallocate(d1BOi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOi')
    deallocate(d1i,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1i')
    deallocate(devaldsbo_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','devaldsbo_neigh')
    deallocate(devalddeltaj_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','devalddeltaj_neigh')
    deallocate(detorsdBOpi_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','detorsdBOpi_neigh')
    deallocate(devaldBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','devaldBO_neigh')
    deallocate(debpendBO_neigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','debpendBO_neigh')
    deallocate(d1BOp_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOp_pipi')
    deallocate(d1BOp_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOp_pi')
    deallocate(d1BOp,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','d1BOp')
    deallocate(deltalplp,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','deltalplp')
    deallocate(deltalp,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','deltalp')
    deallocate(deltap,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','deltap')
    deallocate(delta,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','delta')
    deallocate(BOp_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','BOp_pipi')
    deallocate(BOp_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','BOp_pi')
    deallocate(BOp,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','BOp')
    deallocate(BO_pipi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','BO_pipi')
    deallocate(BO_pi,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','BO_pi')
    deallocate(BO,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','BO')
    deallocate(zneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','zneigh')
    deallocate(yneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','yneigh')
    deallocate(xneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','xneigh')
    deallocate(rneigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','rneigh')
    deallocate(ineigh,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','ineigh')
    deallocate(neighnoRptr,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','neighnoRptr')
    deallocate(neighno,stat=status)
    if (status/=0) call deallocate_error('reaxFFmds','neighno')
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
!  Allocate local memory
!
  allocate(neighno(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','neighno')
  allocate(neighnoRptr(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','neighnoRptr')
  allocate(ineigh(3_i4,maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','ineigh')
  allocate(rneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','rneigh')
  allocate(xneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','xneigh')
  allocate(yneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','yneigh')
  allocate(zneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','zneigh')
  allocate(BO(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','BO')
  allocate(BO_pi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','BO_pi')
  allocate(BO_pipi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','BO_pipi')
  allocate(BOp(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','BOp')
  allocate(BOp_pi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','BOp_pi')
  allocate(BOp_pipi(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','BOp_pipi')
  allocate(delta(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','delta')
  allocate(deltap(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','deltap')
  allocate(deltalp(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','deltalp')
  allocate(deltalplp(numat),stat=status)
  if (status/=0) call outofmemory('reaxFFmds','deltalplp')
  if (lgrad1) then
    allocate(nonzeroptr(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr')
    allocate(nonzeroptr_ij(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_ij')
    allocate(nonzeroptr_ij2(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_ij2')
    allocate(nonzeroptr_ik(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_ik')
    allocate(nonzeroptr_jk(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_jk')
    allocate(nonzeroptr_j2k2(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_j2k2')
    allocate(d1BOp(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOp')
    allocate(d1BOp_pi(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOp_pi')
    allocate(d1BOp_pipi(maxneigh,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOp_pipi')
    allocate(debpendBO_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','debpendBO_neigh')
    allocate(devaldBO_neigh(maxneigh1a),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','devaldBO_neigh')
    allocate(detorsdBOpi_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','detorsdBOpi_neigh')
    allocate(devalddeltaj_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','devalddeltaj_neigh')
    allocate(devaldsbo_neigh(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','devaldsbo_neigh')
    allocate(d1i(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1i')
    allocate(d1BOi(maxneigh2e),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOi')
    allocate(d1BOi_s(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOi_s')
    allocate(d1BOi_pi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOi_pi')
    allocate(d1BOi_pipi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOi_pipi')
    allocate(d1BOij2(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOij2')
    allocate(d1BOij2_s(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOij2_s')
    allocate(d1BOij2_pi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOij2_pi')
    allocate(d1BOij2_pipi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOij2_pipi')
    allocate(d1BOik(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOik')
    allocate(d1BOik_s(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOik_s')
    allocate(d1BOik_pi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOik_pi')
    allocate(d1BOik_pipi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOik_pipi')
    allocate(d1BOjk(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk')
    allocate(d1BOjk_s(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk_s')
    allocate(d1BOjk_pi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk_pi')
    allocate(d1BOjk_pipi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk_pipi')
    allocate(d1BOjk2(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk2')
    allocate(d1BOjk2_s(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk2_s')
    allocate(d1BOjk2_pi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk2_pi')
    allocate(d1BOjk2_pipi(maxneigh2t),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk2_pipi')
    allocate(dsboproddbo(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','dsboproddbo')
  else
    allocate(nonzeroptr(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr')
    allocate(nonzeroptr_ij(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_ij')
    allocate(nonzeroptr_ij2(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_ij2')
    allocate(nonzeroptr_ik(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_ik')
    allocate(nonzeroptr_jk(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_jk')
    allocate(nonzeroptr_j2k2(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','nonzeroptr_j2k2')
    allocate(d1BOp(1,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOp')
    allocate(d1BOp_pi(1,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOp_pi')
    allocate(d1BOp_pipi(1,numat),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOp_pipi')
    allocate(debpendBO_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','debpendBO_neigh')
    allocate(devaldBO_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','devaldBO_neigh')
    allocate(detorsdBOpi_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','detorsdBOpi_neigh')
    allocate(devalddeltaj_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','devalddeltaj_neigh')
    allocate(devaldsbo_neigh(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','devaldsbo_neigh')
    allocate(d1i(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1i')
    allocate(d1BOi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOi')
    allocate(d1BOi_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOi_s')
    allocate(d1BOi_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOi_pi')
    allocate(d1BOi_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOi_pipi')
    allocate(d1BOij2(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOij2')
    allocate(d1BOij2_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOij2_s')
    allocate(d1BOij2_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOij2_pi')
    allocate(d1BOij2_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOij2_pipi')
    allocate(d1BOik(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOik')
    allocate(d1BOik_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOik_s')
    allocate(d1BOik_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOik_pi')
    allocate(d1BOik_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOik_pipi')
    allocate(d1BOjk(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk')
    allocate(d1BOjk_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk_s')
    allocate(d1BOjk_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk_pi')
    allocate(d1BOjk_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk_pipi')
    allocate(d1BOjk2(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk2')
    allocate(d1BOjk2_s(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk2_s')
    allocate(d1BOjk2_pi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk2_pi')
    allocate(d1BOjk2_pipi(1),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','d1BOjk2_pipi')
    allocate(dsboproddbo(maxneigh),stat=status)
    if (status/=0) call outofmemory('reaxFFmds','dsboproddbo')
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
      call stopnow('reaxFFmds')
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
  if (nprocs.gt.1) then
    numattodo = 0
    numattodoRptr(1:numat) = 0
    do i = 1+procid,numat,nprocs
      numattodo = numattodo + 1
      numattodoptr(numattodo) = i
      numattodoRptr(i) = numattodo
    enddo
    ii = 0
    do i = 1,numat
      numattodonptr(i) = ii
      ii = ii + 1
      ii = mod(ii,nprocs)
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
  else
    numattodo = numat
    numattodo_main = numat
    do i = 1,numat
      numattodoptr(i) = i
      numattodonptr(i) = 0_i4
      numattodoRptr(i) = i
    enddo
  endif
!*******************************
!  Sort neighbours into order  *
!******************************* 
  if (lspatialboOK) then
    do i = 1,numat
!               
!  Build pointer
!               
      do nn = 1,nneigh(i)
        nmin = numat + 1 
        do nn2 = nn,nneigh(i) 
          if (neighno(nn2,i).lt.nmin) then
            nmin = neighno(nn2,i)
            nptr = nn2  
          endif
        enddo       
!         
!  Sort quantities
!
        if (nptr.ne.nn) then
          itmp = neighno(nptr,i)
          neighno(nptr,i) = neighno(nn,i)
          neighno(nn,i)  = itmp
          itmp = ineigh(1,nptr,i)
          ineigh(1,nptr,i) = ineigh(1,nn,i)
          ineigh(1,nn,i)  = itmp
          itmp = ineigh(2,nptr,i)
          ineigh(2,nptr,i) = ineigh(2,nn,i)
          ineigh(2,nn,i)  = itmp
          itmp = ineigh(3,nptr,i)
          ineigh(3,nptr,i) = ineigh(3,nn,i)
          ineigh(3,nn,i)  = itmp
          rtmp = rneigh(nptr,i)
          rneigh(nptr,i) = rneigh(nn,i)
          rneigh(nn,i)  = rtmp
          rtmp = xneigh(nptr,i)
          xneigh(nptr,i) = xneigh(nn,i)
          xneigh(nn,i)  = rtmp
          rtmp = yneigh(nptr,i)
          yneigh(nptr,i) = yneigh(nn,i)
          yneigh(nn,i)  = rtmp
          rtmp = zneigh(nptr,i)
          zneigh(nptr,i) = zneigh(nn,i)
          zneigh(nn,i)  = rtmp
        endif  
      enddo         
    enddo
  endif
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
!*************************************************************************
!  Loop over pairs of atoms to compute bond order prime and delta prime  *
!*************************************************************************
!
!  mneigh contains the real value of maxneigh need after removing neighbours
!  whose bond order tolerance is below the allowed threshold
!
  mneigh = 0
  do ii = 1,numattodo
    i = numattodoptr(ii)
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
            call reaxFFbo_sigma(nboij,nspeci,nspecj,rij,BOpij,dBOpijdr,d2BOpijdr2,lgrad1,.false., &
                                .not.lStrict,lBOnonzero)
            if (lBOnonzero) then
              lanyBOnonzero = .true.
              BOp(ni,i) = BOpij
              deltap(i) = deltap(i) + BOpij
              if (lgrad1) then
                rrij = 1.0_dp/rij
                d1BOp(ni,i) = rrij*dBOpijdr
              endif
!
!  Pi
!
              call reaxFFbo_pi(nboij,nspeci,nspecj,rij,BOpij,dBOpijdr,d2BOpijdr2,lgrad1,.false., &
                               .not.lStrict)
              BOp_pi(ni,i) = BOpij
              deltap(i) = deltap(i) + BOpij
              if (lgrad1) then
                d1BOp_pi(ni,i) = rrij*dBOpijdr
              endif
!
!  Pi_pi
!
              call reaxFFbo_pipi(nboij,nspeci,nspecj,rij,BOpij,dBOpijdr,d2BOpijdr2,lgrad1,.false., &
                                 .not.lStrict)
              BOp_pipi(ni,i) = BOpij
              deltap(i) = deltap(i) + BOpij
              if (lgrad1) then
                d1BOp_pipi(ni,i) = rrij*dBOpijdr
              endif
            else
              BOp(ni,i) = 0.0_dp
              BOp_pi(ni,i) = 0.0_dp
              BOp_pipi(ni,i) = 0.0_dp
              if (lgrad1) then
                d1BOp(ni,i) = 0.0_dp
                d1BOp_pi(ni,i) = 0.0_dp
                d1BOp_pipi(ni,i) = 0.0_dp
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
!  Set neighnoRptr
!
  do ii = 1,numattodo
    i = numattodoptr(ii)
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
      if (status/=0) call outofmemory('reaxffmds','StatMPI')
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
      if (status/=0) call deallocate_error('reaxffmds','StatMPI')
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
      if (status/=0) call outofmemory('reaxffmds','nntmp')
      allocate(botmp(2*ntmp),stat=status)
      if (status/=0) call outofmemory('reaxffmds','botmp')
      allocate(StatMPI(MPI_Status_Size,2),stat=status)
      if (status/=0) call outofmemory('reaxffmds','StatMPI')
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
      if (status/=0) call deallocate_error('reaxffmds','StatMPI')
      deallocate(botmp,stat=status)
      if (status/=0) call deallocate_error('reaxffmds','botmp')
      deallocate(nntmp,stat=status)
      if (status/=0) call deallocate_error('reaxffmds','nntmp')
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
  if (nprocs.gt.1) then
!
!  Globalise deltap
!
    sum(1:numat) = 0.0_dp
    do ii = 1,numattodo_main
      i = numattodoptr(ii)
      sum(i) = deltap(i)
    enddo
    deltap(1:numat) = 0.0_dp
    call sumall(sum,deltap,numat,"reaxffmds","deltap")
  endif
!*****************************
!  Loop over pairs of atoms  *
!*****************************
  iloop: do ii = 1,numattodo
    i = numattodoptr(ii)
!
!  If this atom has no ReaxFF terms then there is nothing to do...
!
    if (nbos(i).eq.0) cycle iloop
!
!  Set variables relating to i
!
    nspeci = nbosptr(i)
    nregioni = nregionno(nrelf2a(i))
    nregiontypi = nregiontype(nregioni,ncf)
    lslicei = lsliceatom(nsft + nrelf2a(i))
!
!  Set total number of distances for neighbours of i
!
    nneighi1   = nneigh(i)*(nneigh(i) + 1)/2 
    nneighi2   = nneigh(i) + nneighi1
    nneighi22  = nneighi2*(nneighi2 + 1)/2
!
    nneighi2t  = nneigh(i) + nneighi1 + nneigh(i)*(mneigh+1)*mneigh
    if (lextended) then
      nneighi2e  = nneigh(i) + nneighi1 + nneigh(i)*(2*mneigh+1)*mneigh
    else
      nneighi2e  = nneighi2t
    endif
    nneighi2e2 = nneighi2e*(nneighi2e + 1)/2
    nneighi2t2 = nneighi2t*(nneighi2t + 1)/2
!
!  Initialise derivative storage for neighbours of i
!
    if (lgrad1) then
      d1i(1:nneighi2e) = 0.0_dp
    endif
    if (nprocs.gt.1) then
!---------------------
!  Parallel version  |
!---------------------
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
          call reaxFF_bomds(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                            d1BOp,d1BOp_pi,d1BOp_pipi,BOij_s,BOij_pi,BOij_pipi,d1BOi_s,d1BOi_pi, &
                            d1BOi_pipi,nonzero_ij,nonzero_ij_i,nonzeroptr_ij,lgrad1)
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
          BO_pi(ni,i) = BOij_pi
          BO_pipi(ni,i) = BOij_pipi
          if (lreaxff_ebond.and.lQMMMok.and.ii.le.numattodo_main) then
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
!
!  Derivatives of Bond Order potential energy
!
            if (lgrad1) then
              deijdBO_s    = - scale*reaxFFDe(1,nboij)*expij*(1.0_dp - reaxFFpbe(1,nboij)*reaxFFpbe(2,nboij)*BOijpbe2)
              deijdBO_pi   = - scale*reaxFFDe(2,nboij) 
              deijdBO_pipi = - scale*reaxFFDe(3,nboij) 
              do nk = 1,nonzero_ij
                k = nonzeroptr_ij(nk)
                d1i(k) = d1i(k) + deijdBO_s*d1BOi_s(nk) + deijdBO_pi*d1BOi_pi(nk) + deijdBO_pipi*d1BOi_pipi(nk)
              enddo
            endif
          endif
!
!  End condition section on i or j being associated with moving atom
!
        endif
      enddo
    else
!-------------------
!  Serial version  |
!-------------------
!
!  Loop over neighbours of i (=> j)
!
      do ni = 1,nneigh(i)
        j = neighno(ni,i)
!
!  If j > i then exit loop
!
        if (j.gt.i) exit
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
!  If i = j set scale to be half to correct for double counting
!
          if (i.eq.j) then
            scale = 0.5_dp
          else
            scale = 1.0_dp
          endif
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
          call reaxFF_bomds(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                            d1BOp,d1BOp_pi,d1BOp_pipi,BOij_s,BOij_pi,BOij_pipi,d1BOi_s,d1BOi_pi, &
                            d1BOi_pipi,nonzero_ij,nonzero_ij_i,nonzeroptr_ij,lgrad1)
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
!
!  Derivatives of Bond Order potential energy
!
            if (lgrad1) then
              deijdBO_s    = - scale*reaxFFDe(1,nboij)*expij*(1.0_dp - reaxFFpbe(1,nboij)*reaxFFpbe(2,nboij)*BOijpbe2)
              deijdBO_pi   = - scale*reaxFFDe(2,nboij)
              deijdBO_pipi = - scale*reaxFFDe(3,nboij)
              do nk = 1,nonzero_ij
                k = nonzeroptr_ij(nk)
                d1i(k) = d1i(k) + deijdBO_s*d1BOi_s(nk) + deijdBO_pi*d1BOi_pi(nk) + deijdBO_pipi*d1BOi_pipi(nk)
              enddo
            endif
          endif
!
!  End condition section on i or j being associated with moving atom
!
        endif
      enddo
    endif
!
!  Add derivatives due to neighbours of i
!
    if (lreaxff_ebond.and.lgrad1.and.ii.le.numattodo_main) then
      call d1add(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nbodummyptr,d1i,.true.,.false.,.true.)
    endif
  enddo iloop
!***************************************************
!  Loop over neighbours of atoms to compute delta  *
!***************************************************
  delta(1:numat) = 0.0_dp
  deltalp(1:numat) = 0.0_dp
  deltalplp(1:numat) = 0.0_dp
  ddeltaeddeltai = 0.5_dp
  do ii = 1,numattodo_main
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
  if (nprocs.gt.1) then
!   
!  Global sum of delta terms
!   
    sum(1:numat) = delta(1:numat)
    sum(numat+1:2*numat) = deltalp(1:numat)
    sum(2*numat+1:3*numat) = deltalplp(1:numat)
    call sumall(sum,sum0,3_i4*numat,"reaxffmds","delta")
    delta(1:numat) = sum0(1:numat)
    deltalp(1:numat) = sum0(numat+1:2*numat)
    deltalplp(1:numat) = sum0(2*numat+1:3*numat)
  endif
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
      if (status/=0) call outofmemory('reaxffmds','nntmp')
      allocate(botmp(2*ntmp),stat=status)
      if (status/=0) call outofmemory('reaxffmds','botmp')
      allocate(StatMPI(MPI_Status_Size,2),stat=status)
      if (status/=0) call outofmemory('reaxffmds','StatMPI')
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
      if (status/=0) call deallocate_error('reaxffmds','StatMPI')
      deallocate(botmp,stat=status)
      if (status/=0) call deallocate_error('reaxffmds','botmp')
      deallocate(nntmp,stat=status)
      if (status/=0) call deallocate_error('reaxffmds','nntmp')
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
  do ii = 1,numattodo_main
    i = numattodoptr(ii)
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
          nneighi2e2 = nneighi2e*(nneighi2e + 1)/2
          nneighi2t2 = nneighi2t*(nneighi2t + 1)/2
!
!  Initialise derivative storage for neighbours of i
!
          d1i(1:nneighi2e) = 0.0_dp
        endif
!
!  Lone pair energy
!
        if (lreaxFF_elone) then
          if (abs(reaxFFlp(2,nspeci)).gt.1.0d-12) then
            explp = exp(-reaxFFlam(30)*deltalplp(i))
            trm1lp = 1.0_dp/(1.0_dp + explp)
            esite = reaxFFlp(2,nspeci)*deltalplp(i)*trm1lp
            elone = elone + esite
!
!  Site energy contributions for elone
!
            siteenergy(i) = siteenergy(i) + esite
          endif
        endif
!
!  Bond penalty energy
!
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
          esite = sumover*otrm2*deltalpcorr_i*otrm3
          eover = eover + esite
!
!  Site energy contributions for eover
!
          siteenergy(i) = siteenergy(i) + esite
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
          esite  = - reaxFFuc1(nspeci)*(1.0_dp - expuc1)*utrm1*utrm2
          eunder = eunder + esite
!
!  Site energy contributions for eunder
!
          siteenergy(i) = siteenergy(i) + esite
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
          else
!
!  Old discontinuous form
!
            rintde = dble(int(delta_e))
            drintdeddeltae = 0.0_dp
          endif
!
!  Derivatives of deltalp
!
          trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
          explp2 = exp(-reaxFFlam(29)*trm2lp**2)
          ddeltalpddeltai = 0.5_dp*(drintdeddeltae + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(1.0_dp - drintdeddeltae)) ! 1/2 is to convert from delta_e to delta_i
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
          else
!
!  Old discontinuous form
!
            rintde = dble(int(delta_e))
            drintdeddeltae = 0.0_dp
          endif
!
          trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
          explp2 = exp(-reaxFFlam(29)*trm2lp**2)
          ddeltalplpddeltai = 0.5_dp*(drintdeddeltae + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(1.0_dp - drintdeddeltae)) ! 1/2 is to convert from delta_e to delta_i
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
          else
!
!  Old discontinuous form
!
            rintde = dble(int(delta_e))
            drintdeddeltae = 0.0_dp
          endif
!
!  Derivatives of deltalp
!
          trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
          explp2 = exp(-reaxFFlam(29)*trm2lp**2)
          ddeltalpddeltai = 0.5_dp*(drintdeddeltae + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(1.0_dp - drintdeddeltae)) ! 1/2 is to convert from delta_e to delta_i
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
            else
              call ataper(.false.,BO(ni,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fij,dfijdBOij,d2fijdBOij2, &
                          d3fijdBOij3,lgrad1,.false.,.false.)
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
                else
                  call ataper(.false.,BO(nk,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fik,dfikdBOik,d2fikdBOik2, &
                              d3fikdBOik3,lgrad1,.false.,.false.)
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
                                  d2thetadr2,lgrad1,.false.)
                if (lBOijOK.and.lBOikOK.and.(BO(ni,i)*BO(nk,i).gt.reaxFFatol2)) then
!------------------
!  Valence energy |
!------------------
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
                                     d2f7dBOij2,d2f7dBOijdBOik,d2f7dBOik2,lgrad1,.false.)
                      call reaxFF_f8(nspeci,nspecj,nspeck,nval3,delta(i),f8,df8ddeltai,d2f8ddeltai2,lgrad1,.false.)
                      call reaxFF_theta0(nspeci,nspecj,nspeck,nval3,sbo,theta0,dtheta0dsbo,d2theta0dsbo2,lgrad1,.false.)
!
!  Set angle dependence
!
                      exptheta = exp(-pval2*(theta0 - theta)**2)
                      fangle = (1.0_dp - exptheta)
                      if (lgrad1) then
                        dfangledtheta = - 2.0_dp*pval2*(theta0 - theta)*exptheta
                        dfangledtheta0 = - dfangledtheta
                      endif
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
                        d1i(ni) = d1i(ni) + devaldtheta*dthetadrij
                        d1i(nk) = d1i(nk) + devaldtheta*dthetadrik
                        indjk = nneigh(i) + ni*(ni - 1)/2 + nk
                        d1i(indjk) = d1i(indjk) + devaldtheta*dthetadrjk
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
                    call reaxFF_f9(delta(i),f9,df9ddeltai,d2f9ddeltai2,lgrad1,.false.)
!
!  Evaluate penalty energy
!
                    ppen2 = reaxFFlam(20)
                    exppen1 = exp(-ppen2*(BOij - 2.0_dp)**2)
                    exppen2 = exp(-ppen2*(BOik - 2.0_dp)**2)
                    epentrm = reaxFFpen3(indsjk,nspeci)*f9*exppen1*exppen2
                    esite = fij*fik*epentrm
                    epen = epen + esite
!
!  Site energy for epen
!
                    esite = esite/3.0_dp
                    siteenergy(i) = siteenergy(i) + esite
                    siteenergy(j) = siteenergy(j) + esite
                    siteenergy(k) = siteenergy(k) + esite
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
                    esite = fij*fik*ecoatrm
                    ecoa = ecoa + esite
!
!  Site energy for ecoa
!
                    esite = esite/3.0_dp
                    siteenergy(i) = siteenergy(i) + esite
                    siteenergy(j) = siteenergy(j) + esite
                    siteenergy(k) = siteenergy(k) + esite
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
                    endif
                  endif
                endif
!
              endif
            enddo
!
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
              else                                 
                call ataper(.false.,BO(ni,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fij,dfijdBOij,d2fijdBOij2,d3fijdBOij3, &
                            lgrad1,.false.,.false.)
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
                    else
                      call ataper(.false.,BO(nk,i),reaxFFatol,reaxFFtaperscale*reaxFFatol,fik,dfikdBOik,d2fikdBOik2, &
                                  d3fikdBOik3,lgrad1,.false.,.false.)
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
                        else
                          call ataper(.false.,BO(nl,j),reaxFFatol,reaxFFtaperscale*reaxFFatol,fjl,dfjldBOjl,d2fjldBOjl2, &
                                      d3fjldBOjl3,lgrad1,.false.,.false.)
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
                                                 d2sinkijdr2,lgrad1,.false.)
                            call reaxff_sintheta(rneigh(ni,i),rneigh(nl,j),ril,sinijl,dsinijldr(1),dsinijldr(2),dsinijldr(3), &
                                                 d2sinijldr2,lgrad1,.false.)
!
!  Compute torsion angle
!
                            call reaxff_cosphi(rneigh(nk,i),rjk,rkl,rneigh(ni,i),ril,rneigh(nl,j),cosphi,cosp1d,cosp2d, &
                                               lgrad1,.false.)
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
                                            d2f10dBOijdBOjl,d2f10dBOjl2,lgrad1,.false.)
                            call reaxFF_f11(nspeci,nspecj,delta(i),delta(j),f11,df11ddeltaij,d2f11ddeltaij2,lgrad1,.false.)
!
                            if (lgrad1) then
                              df11ddeltai = df11ddeltaij
                              df11ddeltaj = df11ddeltaij
                            endif
!
                            call reaxFF_f12(BOik,BOij,BOjl,f12,df12dBOik,df12dBOij,df12dBOjl, &
                                            d2f12dBOik2,d2f12dBOikdBOij,d2f12dBOikdBOjl,d2f12dBOij2, &
                                            d2f12dBOijdBOjl,d2f12dBOjl2,lStrict,lgrad1,.false.)
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
!
                            if (lgrad1) then
                              indbojl = ni*maxneigh + nl   ! 2nd shell index for l in shell of j
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
                              do m = 1,6
                                d1i(indcosphi(m)) = d1i(indcosphi(m)) + detorsdcosphi*cosp1d(m)
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
!  Set indices for sines
!
                              indsinkij(1) = ni
                              indsinkij(2) = nk
                              indsinkij(3) = indjk
                              indsinijl(1) = ni
                              indsinijl(2) = indjl
                              indsinijl(3) = indil
!
                              do m = 1,3
                                d1i(indsinkij(m)) = d1i(indsinkij(m)) + detorsdsinkij*dsinkijdr(m)
                                d1i(indsinijl(m)) = d1i(indsinijl(m)) + detorsdsinijl*dsinijldr(m)
                              enddo
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
          else
            ddeltalpcorrddeltai = 1.0_dp - otrm1*ddeltalpddeltai
          endif
!
!  Set up term needed later for derivative of deltalpcorr_i with respect to sumdeltabopi
!
          if (.not.lreaxFFunder(nspeci)) then
            ddeltalpcorrdsumdeltabopi = 0.0_dp
          else
            ddeltalpcorrdsumdeltabopi = (otrm1**2)*expdi*reaxFFlam(6)*reaxFFlam(31)*deltalp(i)
          endif
!---------------------------------------------------------------------------------
!  Set up terms for partial derivatives of eunder                                |
!---------------------------------------------------------------------------------
          dexpuc1ddeltalpcorr = reaxFFlam(7)*expuc1
          dutrm1ddeltalpcorr  = utrm1*utrm1*expuc2*reaxFFoc1(nspeci)
          dutrm2dsumdeltabopi = - utrm2*utrm2*expuc3*reaxFFlam(8)*reaxFFlam(9)
!
!  Compute derivatives of eunder w.r.t. deltalpcorr and sumdeltabopi
!
          deunderddeltalpcorr = - reaxFFuc1(nspeci)*(-dexpuc1ddeltalpcorr*utrm1 + (1.0_dp - expuc1)*dutrm1ddeltalpcorr)*utrm2
          deunderdsumdeltabopi = - reaxFFuc1(nspeci)*(1.0_dp - expuc1)*utrm1*dutrm2dsumdeltabopi
!---------------------------------------------------------------------------------
!  Set up terms for partial derivatives of eover                                 |
!---------------------------------------------------------------------------------
          deoverdsumover = otrm2*deltalpcorr_i*otrm3
          dotrm2ddeltalpcorr = - otrm2**2
          dotrm3ddeltalpcorr = - (otrm3**2)*reaxFFoc1(nspeci)*expoc
          deoverddeltalpcorr = sumover*(otrm2*otrm3 + deltalpcorr_i*(dotrm2ddeltalpcorr*otrm3 + otrm2*dotrm3ddeltalpcorr))
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
!---------------------------------------------------------------------------------
!  Set up terms for partial derivatives of ebpen                                 |
!---------------------------------------------------------------------------------
          debpenddeltai = 0.0_dp
          debpendBO_neigh(1:nneigh(i)) = 0.0_dp
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
                endif
              endif
            enddo
          endif
!---------------------------------------------------------------------------------
!  Sum up terms for partial derivatives of energy                                |
!---------------------------------------------------------------------------------
          detotddeltalpcorr = deunderddeltalpcorr + deoverddeltalpcorr
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
              nneighj2 = nneigh(j)*(nneigh(j)+1)/2
              nneighj2t = nneigh(j) + nneighj2 + nneigh(j)*(mneigh+1)*mneigh
              nneighj2t2 = nneighj2t*(nneighj2t+1)/2
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
              else
                rintde = dble(int(delta_e))
                drintdeddeltae = 0.0_dp
              endif
!
              trm2lp = 2.0_dp*(1.0_dp + delta_e - rintde)
              explp2 = exp(-reaxFFlam(29)*trm2lp**2)
              ddeltalpddeltaj = 0.5_dp*(drintdeddeltae + 4.0_dp*reaxFFlam(29)*trm2lp*explp2*(1.0_dp - drintdeddeltae)) ! 1/2 is to convert from delta_e to delta_i
!
!  Set delta' 
!
              deltapi = deltap(i)
              deltapj = deltap(j)
!***************************************
!  Valid pairwise reaxFF contribution  *
!***************************************
              call reaxFF_bomds(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                                d1BOp,d1BOp_pi,d1BOp_pipi,BOij_s,BOij_pi,BOij_pipi,d1BOi_s,d1BOi_pi, &
                                d1BOi_pipi,nonzero_ij,nonzero_ij_i,nonzeroptr_ij,lgrad1)
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
                dsumdeltabopiddeltaj = (BO_pi(ni,i) + BO_pipi(ni,i))
              else
                dsumdeltabopidbopi = delta(j) - deltalp(j)
                dsumdeltabopiddeltaj = (1.0_dp - ddeltalpddeltaj)*(BO_pi(ni,i) + BO_pipi(ni,i))
              endif
!
!  Build derivatives of deltalpcorr w.r.t. BO pi and pi-pi and deltaj, plus mixed terms
!
              ddeltalpcorrdbopi = ddeltalpcorrdsumdeltabopi*dsumdeltabopidbopi
              ddeltalpcorrddeltaj = ddeltalpcorrdsumdeltabopi*dsumdeltabopiddeltaj
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
              do nk = 1,nonzero_ij
                k = nonzeroptr_ij(nk)
!
!  Derivatives w.r.t. BOij (via deltai)
!
                d1i(k) = d1i(k) + detotdbo_tot*d1BOi(nk)
!
!  Derivatives w.r.t. BOij from angles and torsions via deltai and BOij
!
                d1i(k) = d1i(k) + (devalddeltai + devaldBO_neigh(ni) + devaldsbo_neigh(ni))*d1BOi(nk)
!
!  Derivatives w.r.t. BOij_pi and BOij_pipi
!
                d1i(k) = d1i(k) + detotdbopi_tot*(d1BOi_pi(nk) + d1BOi_pipi(nk))
!
!  Derivatives w.r.t. BOij_pi and BOij_pipi from angles
!
                d1i(k) = d1i(k) + devaldsbo*(d1BOi_pi(nk) + d1BOi_pipi(nk))*sbopi
!
!  Derivatives of etors with respect to BOpi_ij
!
                d1i(k) = d1i(k) + detorsdBOpi_neigh(ni)*d1BOi_pi(nk)
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
!  Calculate total bond order derivatives from sum of components
!
                  do nl = 1,nonzero_jk
                    d1BOjk(nl) = d1BOjk_s(nl) + d1BOjk_pi(nl) + d1BOjk_pipi(nl)
                  enddo
!
                  do nl = 1,nonzero_jk
                    l = nonzeroptr_jk(nl)
                    d1i(l) = d1i(l) + (detotddeltaj_tot + devalddeltaj_neigh(ni))*d1BOjk(nl)
                  enddo
!
!  Derivatives w.r.t. BOjl from torsions 
!
                  indbojk = ni*maxneigh + nj
                  do nl = 1,nonzero_jk
                    l = nonzeroptr_jk(nl)
                    d1i(l) = d1i(l) + devaldBO_neigh(indbojk)*d1BOjk(nl)
                  enddo
                endif
!
!  End of loop over neighbours of j
!
              enddo
!********************************
!  End of delta_j contribution  *
!********************************
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
          call d1add(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nbodummyptr,d1i,.true.,lextended,.true.)
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
  if (lspatialboOK) then
!************************************
!  Spatial decomposition algorithm  *
!************************************
    maxxy = nspcellbo(1)*nspcellbo(2)
    maxx  = nspcellbo(1)
!
!  Compute number of cells to be search for hydrogen bonds
!
    ncellsearchHB(1) = 1 + reaxFFrhtol/rnearestxbo
    ncellsearchHB(2) = 1 + reaxFFrhtol/rnearestybo
    ncellsearchHB(3) = 1 + reaxFFrhtol/rnearestzbo
!
!  Loop over atoms looking for hydrogen
!
    do i = procid+1,numat,nprocs
      if (nbos(i).gt.0.and.nat(i).eq.1) then
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
          if (lgrad1) then
!           
!  Set total number of distances for neighbours of i
!             
            nneighi1   = nneigh(i)*(nneigh(i) + 1)/2
            nneighi2   = nneigh(i) + nneighi1 
            nneighi2t  = nneigh(i) + nneighi1 + nneigh(i)*(mneigh+1)*mneigh
            nneighi22  = nneighi2*(nneighi2 + 1)/2
            nneighi2t2 = nneighi2t*(nneighi2t + 1)/2
!           
!  Initialise derivative storage for neighbours of i
!       
            d1i(1:nneighi2t) = 0.0_dp
          endif
!
!  Loop over atoms that are connected by bond order to i
!
          do ni = 1,nneigh(i)
            j = neighno(ni,i)
            if (nbos(j).gt.0) then
              nspecj = nbosptr(j)
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
                  else
                    call ataper(.false.,BO(ni,i),reaxFFhtol,reaxFFtaperscale*reaxFFhtol,fHB,dfHBdBO,d2fHBdBO2,d3fHBdBO3, &
                                lgrad1,.false.,.false.)
                  endif
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
                    call reaxFF_bomds(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                                      d1BOp,d1BOp_pi,d1BOp_pipi,BOij_s,BOij_pi,BOij_pipi,d1BOi_s,d1BOi_pi, &
                                      d1BOi_pipi,nonzero_ij,nonzero_ij_i,nonzeroptr_ij,lgrad1)
!
!  Calculate total bond order derivatives from sum of components
!
                    do nk = 1,nonzero_ij
                      d1BOi(nk) = d1BOi_s(nk) + d1BOi_pi(nk) + d1BOi_pipi(nk)
                    enddo
                  endif
!
!  Get cell index for j - by definition this should not be a buffer cell
!
                  ind1 = natomcellbo(j)
                  ind2 = ind1 - 1
                  isz = ind2/maxxy
                  ind2 = ind2 - maxxy*isz
                  isy = ind2/maxx
                  isx = ind2 - maxx*isy + 1
                  isy = isy + 1
                  isz = isz + 1
!
!  Set cell search bounds
!
                  nspupper(1) = min(isx+ncellsearchHB(1),nspcellbo(1))
                  nspupper(2) = min(isy+ncellsearchHB(2),nspcellbo(2))
                  nspupper(3) = min(isz+ncellsearchHB(3),nspcellbo(3))
                  nsplower(1) = max(isx-ncellsearchHB(1),1)
                  nsplower(2) = max(isy-ncellsearchHB(2),1)
                  nsplower(3) = max(isz-ncellsearchHB(3),1)
!  
!  Set coordinates of atom j
!           
                  xj = xinboxbo(j)
                  yj = yinboxbo(j)
                  zj = zinboxbo(j)
!  
!  Loop over neighbouring cells
!           
                  do imz = nsplower(3),nspupper(3)
                    do imy = nsplower(2),nspupper(2)
                      do imx = nsplower(1),nspupper(1)
                        indn = (imz-1)*maxxy + (imy-1)*maxx + imx
! 
!  Loop over atoms within neighbouring cells
!
                        nk = nspcellatbo(indn)
                        n1k = nspcellat1ptrbo(indn)
                        kkloop: do kk = 1,nk
                          k = nspcellatptrbo(n1k+kk)
                          kc = nspcellatptrcellbo(n1k+kk)
!
!  Does k have a reaxFF species?
!
                          if (nbos(k).gt.0) then
                            nspeck = nbosptr(k)
!
!  Is this triad valid for a hydrogen bond?
!
                            phb1 = reaxFFhb3(2,nspeck,nspecj,nspeci)
                            if (abs(phb1).lt.1.0d-12) cycle kkloop
!
!  Compute basic interatomic vector
!
                            xkj = xvec2cell(kc) + xinboxbo(k) - xj
                            ykj = yvec2cell(kc) + yinboxbo(k) - yj
                            zkj = zvec2cell(kc) + zinboxbo(k) - zj
!
!  Find valid vectors
!
                            rkj2 = xkj*xkj + ykj*ykj + zkj*zkj
                            if (rkj2.gt.cut2) cycle kkloop
!
!  Exclude j = k case
!
                            if (j.eq.k.and.rkj2.lt.1.0d-12) cycle kkloop
!
!  Calculate distances for i - k
!
                            xki = xkj + xji
                            yki = ykj + yji
                            zki = zkj + zji
                            rki2 = xki*xki + yki*yki + zki*zki
!
!  Exclude i = k case 
!
                            if (i.eq.k.and.rki2.lt.1.0d-12) cycle kkloop
!
                            rki = sqrt(rki2)
                            rkj = sqrt(rkj2)
!
!  Set remaining parameters
!
                            r0hb = reaxFFhb3(1,nspeck,nspecj,nspeci)
                            phb2 = reaxFFhb3(3,nspeck,nspecj,nspeci)
                            phb3 = reaxFFhb3(4,nspeck,nspecj,nspeci)
!  
!  Compute theta for j - i - k
!                     
                            call reaxff_theta(rji,rki,rkj,theta,dthetadrij,dthetadrik,dthetadrjk, &
                                              d2thetadr2,lgrad1,.false.)
!
!  Compute distance taper function
!
                            if (lStrict) then
                              frHB = 1.0_dp
                              dfrHBdr = 0.0_dp
                            else
                              call ataper(.true.,rkj,reaxFFrhtollower,reaxFFrhtol,frHB,dfrHBdr,d2frHBdr2,d3frHBdr3, &
                                          lgrad1,.false.,.false.)
!
!  Convert taper derivatives to correct form
!
                              if (lgrad1) then
                                rrjk = 1.0_dp/rkj
                                dfrHBdr = rrjk*dfrHBdr
                              endif
                            endif
!
!  Compute remaining terms for hydrogen bond energy
!
                            rrik = 1.0_dp/rki
                            rrij = 1.0_dp/rji
                            sinhalftheta  = sin(0.5_dp*theta)
                            sin4halftheta = sinhalftheta**4
                            exphb1 = exp(-phb2*BO(ni,i))
                            exphb2 = exp(-phb3*((r0hb*rrik) + (rki/r0hb) - 2.0_dp))
!
                            ehb = ehb + fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta
!
!  Site energy for ehb
!
                            esite = fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta/3.0_dp
                            siteenergy(i) = siteenergy(i) + esite
                            siteenergy(j) = siteenergy(j) + esite
                            siteenergy(k) = siteenergy(k) + esite
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
!  Set region for k
!
                              nregionk = nregionno(nrelf2a(k))
!
!  Hydrogen bond derivatives with respect to the distance i-j
!
                              dehbdrij = fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrij
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
                                call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,1), &
                                  d2r2dx2(1,1),d2r2dsdx(1,1,1),d2r2ds2(1,1,1),.false.)
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
                              dehbdrik = fHB*frHB*phb1*(1.0_dp - exphb1)*dexphb2drik*sin4halftheta + &
                                         fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrik
!
                              xdrv(i) = xdrv(i) - dehbdrik*xki
                              ydrv(i) = ydrv(i) - dehbdrik*yki
                              zdrv(i) = zdrv(i) - dehbdrik*zki
                              xdrv(k) = xdrv(k) + dehbdrik*xki
                              ydrv(k) = ydrv(k) + dehbdrik*yki
                              zdrv(k) = zdrv(k) + dehbdrik*zki
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
                                call real1strterm(ndim,xki,yki,zki,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,2), &
                                    d2r2dx2(1,2),d2r2dsdx(1,1,2),d2r2ds2(1,1,2),.false.)
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
                              dehbdrjk = fHB*dfrHBdr*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta + &
                                         fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrjk
!
                              xdrv(j) = xdrv(j) - dehbdrjk*xkj
                              ydrv(j) = ydrv(j) - dehbdrjk*ykj
                              zdrv(j) = zdrv(j) - dehbdrjk*zkj
                              xdrv(k) = xdrv(k) + dehbdrjk*xkj
                              ydrv(k) = ydrv(k) + dehbdrjk*ykj
                              zdrv(k) = zdrv(k) + dehbdrjk*zkj
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
                                call real1strterm(ndim,xkj,ykj,zkj,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,3), &
                                  d2r2dx2(1,3),d2r2dsdx(1,1,3),d2r2ds2(1,1,3),.false.)
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
!
!  Derivatives of hydrogen bond energy with respect to bond orders -> add to eval terms
!
                              dehbdBOij = phb1*frHB*(dfHBdBO*(1.0_dp - exphb1) - fHB*dexphb1dBOij)*exphb2*sin4halftheta
!
                              do nl = 1,nonzero_ij
                                l = nonzeroptr_ij(nl)
                                d1i(l) = d1i(l) + dehbdBOij*d1BOi(nl)
                              enddo
                            endif
!
!  End check on whether k is ReaxFF species
!
                          endif
!
!  End of loop over kk
!
                        enddo kkloop
!
!  End loops over neighbouring cells
!
                      enddo
                    enddo
                  enddo
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
!  Add derivatives due to neighbours of i 
!
          if (lgrad1) then
            call d1add(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nbodummyptr,d1i,.true.,.false.,.true.)
          endif
!
!  End of QMMM condition
!
        endif
      endif
    enddo
  else
!***************************
!  Standard N^2 algorithm  *
!***************************
!
!  Loop over atoms looking for hydrogen
!
    do i = procid+1,numat,nprocs
      if (nbos(i).gt.0.and.nat(i).eq.1) then
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
          if (lgrad1) then
!           
!  Set total number of distances for neighbours of i
!             
            nneighi1   = nneigh(i)*(nneigh(i) + 1)/2
            nneighi2   = nneigh(i) + nneighi1 
            nneighi2t  = nneigh(i) + nneighi1 + nneigh(i)*(mneigh+1)*mneigh
            nneighi22  = nneighi2*(nneighi2 + 1)/2
            nneighi2t2 = nneighi2t*(nneighi2t + 1)/2
!           
!  Initialise derivative storage for neighbours of i
!       
            d1i(1:nneighi2t) = 0.0_dp
          endif
!
!  Loop over atoms that are connected by bond order to i
!
          do ni = 1,nneigh(i)
            j = neighno(ni,i)
            if (nbos(j).gt.0) then
              nspecj = nbosptr(j)
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
                  else 
                    call ataper(.false.,BO(ni,i),reaxFFhtol,reaxFFtaperscale*reaxFFhtol,fHB,dfHBdBO,d2fHBdBO2,d3fHBdBO3, &
                                lgrad1,.false.,.false.)
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
                    call reaxFF_bomds(i,j,nspeci,nspecj,ni,nj,deltapi,deltapj,mneigh,nneigh,BOp,BOp_pi,BOp_pipi, &
                                      d1BOp,d1BOp_pi,d1BOp_pipi,BOij_s,BOij_pi,BOij_pipi,d1BOi_s,d1BOi_pi, &
                                      d1BOi_pipi,nonzero_ij,nonzero_ij_i,nonzeroptr_ij,lgrad1)
! 
!  Convert non-zero indices
! 
                    ind0 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 + (ni - 1)*(mneigh+1)*mneigh + ni*mneigh
                    do nl = nonzero_ij_i+1,nonzero_ij
                      nonzeroptr_ij(nl) = ind0 + nonzeroptr_ij(nl)
                    enddo     
!
!  Calculate total bond order derivatives from sum of components
!
                    do nk = 1,nonzero_ij
                      d1BOi(nk) = d1BOi_s(nk) + d1BOi_pi(nk) + d1BOi_pipi(nk)
                    enddo
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
                        else
                          call ataper(.true.,rkj,reaxFFrhtollower,reaxFFrhtol,frHB,dfrHBdr,d2frHBdr2,d3frHBdr3, &
                                      lgrad1,.false.,.false.)
!
!  Convert taper derivatives to correct form
!
                          if (lgrad1) then
                            dfrHBdr = rrjk*dfrHBdr
                          endif
                        endif
!  
!  Compute theta for j - i - k
!                     
                        call reaxff_theta(rji,rki,rkj,theta,dthetadrij,dthetadrik,dthetadrjk, &
                                          d2thetadr2,lgrad1,.false.)
!
!  Compute remaining terms for hydrogen bond energy
!
                        sinhalftheta  = sin(0.5_dp*theta)
                        sin4halftheta = sinhalftheta**4
                        exphb1 = exp(-phb2*BO(ni,i))
                        exphb2 = exp(-phb3*((r0hb*rrik) + (rki/r0hb) - 2.0_dp))
!
                        ehb = ehb + fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta
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
!  Hydrogen bond derivatives with respect to the distance i-j
!
                          dehbdrij = fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrij
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
                            call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,1), &
                              d2r2dx2(1,1),d2r2dsdx(1,1,1),d2r2ds2(1,1,1),.false.)
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
                          dehbdrik = fHB*frHB*phb1*(1.0_dp - exphb1)*dexphb2drik*sin4halftheta + &
                                     fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrik
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
                            call real1strterm(ndim,xki,yki,zki,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,2), &
                              d2r2dx2(1,2),d2r2dsdx(1,1,2),d2r2ds2(1,1,2),.false.)
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
                          dehbdrjk = fHB*dfrHBdr*phb1*(1.0_dp - exphb1)*exphb2*sin4halftheta + &
                                     fHB*frHB*phb1*(1.0_dp - exphb1)*exphb2*dsin4halfthetadtheta*dthetadrjk
!
                          xdrv(j) = xdrv(j) - dehbdrjk*xkj
                          ydrv(j) = ydrv(j) - dehbdrjk*ykj
                          zdrv(j) = zdrv(j) - dehbdrjk*zkj
                          xdrv(k) = xdrv(k) + dehbdrjk*xkj
                          ydrv(k) = ydrv(k) + dehbdrjk*ykj
                          zdrv(k) = zdrv(k) + dehbdrjk*zkj
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
                            call real1strterm(ndim,xkj,ykj,zkj,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,3), &
                              d2r2dx2(1,3),d2r2dsdx(1,1,3),d2r2ds2(1,1,3),.false.)
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
!
!  Derivatives of hydrogen bond energy with respect to bond orders -> add to eval terms
!
                          dehbdBOij = phb1*frHB*(dfHBdBO*(1.0_dp - exphb1) - fHB*dexphb1dBOij)*exphb2*sin4halftheta
!
                          do nl = 1,nonzero_ij
                            l = nonzeroptr_ij(nl)
                            d1i(l) = d1i(l) + dehbdBOij*d1BOi(nl)
                          enddo
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
!  Add derivatives due to neighbours of i 
!
          if (lgrad1) then
            call d1add(i,maxneigh,mneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nbodummyptr,d1i,.true.,.false.,.true.)
          endif
!
!  End of QMMM condition
!
        endif
      endif
    enddo
  endif
  endif   ! End of condition on hydrogen bonding energy being included
!***************************************
!  Twobody VDW & Coulomb contribution  *
!***************************************
  t3 = g_cpu_time()
  if (lreaxFF_ecoul.and.lreaxFFQ) then
    if (literativeQ) then
      call setreaxffQiter(nboatom,nboatomptr,nboatomRptr,nbos,nbosptr,qreaxFF,eself,.false.)
    else
      call setreaxffQ(nboatom,nboatomRptr,nbos,nbosptr,qreaxFF,eself,.false.,.false.)
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
  if (lspatialboOK) then
    derive0self = 0.0_dp
    d2self = 0.0_dp
!************************************
!  Spatial decomposition algorithm  *
!************************************
    maxxy = nspcellbo(1)*nspcellbo(2)
    maxx  = nspcellbo(1)
!
!  Loop over all local spatial cells except buffer regions
!
    do ixyz = 1,ncellpernodebo
      if (.not.lbuffercellbo(ixyz)) then
        ind1 = ncellnodeptrbo(ixyz)
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
        nspupper(1) = min(ix+ncellsearchbo(1),nspcellbo(1))
        nspupper(2) = min(iy+ncellsearchbo(2),nspcellbo(2))
        nspupper(3) = min(iz+ncellsearchbo(3),nspcellbo(3))
        nsplower(1) = max(ix-ncellsearchbo(1),1)
        nsplower(2) = max(iy-ncellsearchbo(2),1)
        nsplower(3) = max(iz-ncellsearchbo(3),1)
!
!  Get number of atoms in this cell
!
        ni = nspcellatbo(ind1)
        n1i = nspcellat1ptrbo(ind1)
!
!  Outer loop over atoms within this cell
!
        do ii = 1,ni
          i = nspcellatptrbo(n1i+ii)
          ib = nboatomptr(i)
!
!  Does i have a reaxFF species or a fixed charge?
!
          nspeci = nbosptr(i)
          lfixQi = lreaxFFqfix(nspeci)
          lreactivei = (nbos(i).gt.0)
!
          if (lreactivei.or.lfixQi) then
            gammai = reaxFFgamma(nspeci)
            ic = nspcellatptrcellbo(n1i+ii)
!
!  Set region for i
!
            nregioni = nregionno(nrelf2a(i))
            nregiontypi = nregiontype(nregioni,ncf)
!
!  Set coordinates of atom i
!
            xi = xinboxbo(i) + xvec2cell(ic)
            yi = yinboxbo(i) + yvec2cell(ic)
            zi = zinboxbo(i) + zvec2cell(ic)
!
!  Loop over neighbouring cells
!         
            do imz = nsplower(3),nspupper(3)
              do imy = nsplower(2),nspupper(2)
                do imx = nsplower(1),nspupper(1)
                  indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!  
!  Loop over atoms within neighbouring cells
!               
                  nj = nspcellatbo(indn)
                  n1j = nspcellat1ptrbo(indn)
                  do jj = 1,nj
                    j = nspcellatptrbo(n1j+jj)
                    jb = nboatomptr(j)
!
!  Does j have a reaxFF species or a fixed charge?
!
                    nspecj = nbosptr(j)
                    lfixQj = lreaxFFqfix(nspecj)
                    lreactivej = (nbos(j).gt.0)
!
                    if (lreactivej.or.lfixQj) then
                      gammaj = reaxFFgamma(nspecj)
                      jc = nspcellatptrcellbo(n1j+jj)
!  
!  Only calculate lower-half triangular interactions
!                 
                      if (j.lt.i.or.(j.eq.i.and.ind1.ne.indn)) then
!  
!  Set coordinate differences and calculate square of distance
!                   
                        xji = xvec2cell(jc) + xinboxbo(j) - xi
                        yji = yvec2cell(jc) + yinboxbo(j) - yi
                        zji = zvec2cell(jc) + zinboxbo(j) - zi
                        rij2 = xji*xji + yji*yji + zji*zji
                        if (rij2.lt.cut2) then
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
                            if (i.eq.j) then
                              scale = 0.5_dp
                            else
                              scale = 1.0_dp
                            endif
!
!  Compute Morse parameters for pair if not explicitly input
!
                            if (nspeci.ge.nspecj) then
                              ind = nspeci*(nspeci - 1)/2 + nspecj
                            else
                              ind = nspecj*(nspecj - 1)/2 + nspeci
                            endif
                            qij = scale*qreaxFF(ib)*qreaxFF(jb)*angstoev
                            rij = sqrt(rij2)
!
!  Compute taper function(s)
!
                            call p7reaxFFvdwtaper(rij,reaxFFcutoffVDW,tp,dtpdr,d2tpdr2,lgrad1,.false.)
                            if (reaxFFcutoffQ.ne.reaxFFcutoffVDW) then
                              call p7reaxFFqtaper(rij,reaxFFcutoffQ,tpQ,dtpQdr,d2tpQdr2,lgrad1,.false.)
                            else
                              tpQ      = tp
                              dtpQdr   = dtpdr
                            endif
                            if (lreaxff_evdw.and.lreactivei.and.lreactivej) then
!
!  Compute f13
!
                              call reaxFF_f13(rij,reaxFFgammaVDW(ind),f13,df13drij,d2f13drij2,lgrad1,.false.)
!
!  Energy terms
!
                              VDWtrm  = reaxFFalphaVDW(ind)*(1.0_dp - f13/reaxFFr0VDW(ind))
                              expVDW1 = exp(0.5_dp*VDWtrm)
                              expVDW2 = expVDW1*expVDW1
!
!  Add VDW energy
!
                              etrm = scale*reaxFFDeVDW(ind)*(expVDW2 - 2.0_dp*expVDW1)
                              evdw = evdw + etrm*tp
!
!  Site energies for evdw
!
                              siteenergy(i) = siteenergy(i) + 0.5_dp*etrm*tp
                              siteenergy(j) = siteenergy(j) + 0.5_dp*etrm*tp
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
                              gam = qij/(rij**3 + gammaij)**third
                              ecoul = ecoul + gam*tpQ
!
!  Site energies for ecoul
!
                              siteenergy(i) = siteenergy(i) + 0.5_dp*gam*tpQ
                              siteenergy(j) = siteenergy(j) + 0.5_dp*gam*tpQ
                            endif
!
                            if (lgrad1) then
                              if (lreaxff_evdw.and.lreactivei.and.lreactivej) then
!
!  Compute derivative of VDW energy
!
                                dVDWtrmdf13 = - reaxFFalphaVDW(ind)/reaxFFr0VDW(ind)
                                dVDWtrmdrij = dVDWtrmdf13*df13drij
                                detrmdrij = scale*reaxFFDeVDW(ind)*(expVDW2 - expVDW1)*dVDWtrmdrij
                                devdwdrij = detrmdrij*tp + etrm*dtpdr
                              else
                                devdwdrij = 0.0_dp
                              endif
!
                              if (lstr) then
                                call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,1), &
                                  d2r2dx2(1,1),d2r2dsdx(1,1,1),d2r2ds2(1,1,1),.false.)
                              endif
!
                              if (lreaxFF_ecoul.and..not.lQMMMelectro) then
!
!  Compute derivative of Coulomb energy
!
                                dgamdrij = - (gam*rij)/(rij*rij2 + gammaij)
                                decouldrij = gam*dtpQdr + dgamdrij*tpQ
                              else
                                decouldrij = 0.0_dp
                              endif
!
!  Get Cartesian components of vector and derivative
!
                              xd = xji*(devdwdrij + decouldrij)
                              yd = yji*(devdwdrij + decouldrij)
                              zd = zji*(devdwdrij + decouldrij)
!             
!  Add derivatives to main arrays
!             
                              xdrv(i) = xdrv(i) - xd
                              ydrv(i) = ydrv(i) - yd
                              zdrv(i) = zdrv(i) - zd
                              xdrv(j) = xdrv(j) + xd          
                              ydrv(j) = ydrv(j) + yd
                              zdrv(j) = zdrv(j) + zd
!
                              if (nregioni.ne.nregionj) then
                                xregdrv(nregioni) = xregdrv(nregioni) - xd
                                yregdrv(nregioni) = yregdrv(nregioni) - yd
                                zregdrv(nregioni) = zregdrv(nregioni) - zd
                                xregdrv(nregionj) = xregdrv(nregionj) + xd          
                                yregdrv(nregionj) = yregdrv(nregionj) + yd
                                zregdrv(nregionj) = zregdrv(nregionj) + zd
                              endif
!
                              if (lstr) then
                                rstrdloc(1:nstrains) = 0.0_dp
                                do kl = 1,nstrains
                                  ns = nstrptr(kl)
                                  rstrdloc(kl) = rstrdloc(kl) + (devdwdrij + decouldrij)*dr2ds(ns,1)
                                enddo
                                do kl = 1,nstrains 
                                  rstrd(kl) = rstrd(kl) + rstrdloc(kl)
                                enddo
                                if (latomicstress) then 
                                  do kl = 1,nstrains
                                    atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*rstrdloc(kl)
                                    atomicstress(kl,j) = atomicstress(kl,j) + 0.5_dp*rstrdloc(kl)
                                  enddo
                                endif
                              endif
                            endif
                          endif
                        endif
                      endif
                    endif
!
!  End loop over j
!
                  enddo 
!               
!  End loops over neighbouring cells
!  
                enddo
              enddo
            enddo
!
!  End check on whether i has a reaxFF species
!
          endif
!
!  End loop over atom i
!
        enddo
!
!  End checks on whether cell is required
!
      endif
!
!  End loop over cells on node
!
    enddo
  else
!***************************
!  Standard N^2 algorithm  *
!***************************
    derive0self = 0.0_dp
    d2self = 0.0_dp
!
!  Set lower bound for i loop
!
    if (ndim.gt.0) then
      imin = 1
    else
      imin = 2
    endif
!
!  Loop over pairs
!
    do ii = imin+procid,nboatom,nprocs
      i = nboatomRptr(ii)
!
!  Does i have a reaxFF species or a fixed charge?
!
      nspeci = nbosptr(i)
      lfixQi = lreaxFFqfix(nspeci)
      lreactivei = (nbos(i).gt.0)
!
      if (lreactivei.or.lfixQi) then
        gammai = reaxFFgamma(nspeci)
!
!  Set region for i
!
        nregioni = nregionno(nrelf2a(i))
        nregiontypi = nregiontype(nregioni,ncf)
!
!  Set upper bound for j loop
!
        if (ndim.gt.0) then
          jmax = ii
        else
          jmax = ii - 1
        endif
!
!  Loop over second atom
!
        jloop: do jj = 1,jmax
          j = nboatomRptr(jj)
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
!  If no distances then cycle
!
            if (nor.eq.0) cycle jloop
!
            if (i.eq.j) then
              scale = 0.5_dp
            else
              scale = 1.0_dp
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
              qij = scale*qreaxFF(ii)*qreaxFF(jj)*angstoev
!
!  Loop over valid distances and calculate contributions
!
              do n = 1,nor
                if (ndim.gt.0) rij2 = dist(n)
                rij = sqrt(rij2)
!
!  Compute taper function(s)
!
                call p7reaxFFvdwtaper(rij,reaxFFcutoffVDW,tp,dtpdr,d2tpdr2,lgrad1,.false.)
                if (reaxFFcutoffQ.ne.reaxFFcutoffVDW) then
                  call p7reaxFFqtaper(rij,reaxFFcutoffQ,tpQ,dtpQdr,d2tpQdr2,lgrad1,.false.)
                else
                  tpQ      = tp
                  dtpQdr   = dtpdr
                endif
                if (lreaxff_evdw.and.lreactivei.and.lreactivej) then
!
!  Compute f13
!
                  call reaxFF_f13(rij,reaxFFgammaVDW(ind),f13,df13drij,d2f13drij2,lgrad1,.false.)
!
!  Energy terms
!
                  VDWtrm  = reaxFFalphaVDW(ind)*(1.0_dp - f13/reaxFFr0VDW(ind))
                  expVDW1 = exp(0.5_dp*VDWtrm)
                  expVDW2 = expVDW1*expVDW1
!
!  Add VDW energy
!
                  etrm = scale*reaxFFDeVDW(ind)*(expVDW2 - 2.0_dp*expVDW1)
                  evdw = evdw + etrm*tp
!
!  Site energies for evdw
!
                  siteenergy(i) = siteenergy(i) + 0.5_dp*etrm*tp
                  siteenergy(j) = siteenergy(j) + 0.5_dp*etrm*tp
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
                  ecoul = ecoul + gam*tpQ
!
!  Site energies for ecoul
!
                  siteenergy(i) = siteenergy(i) + 0.5_dp*gam*tpQ
                  siteenergy(j) = siteenergy(j) + 0.5_dp*gam*tpQ
                endif
!
                if (lgrad1) then
                  if (lreaxff_evdw.and.lreactivei.and.lreactivej) then
!
!  Compute derivatives of VDW energy
!
                    dVDWtrmdf13 = - reaxFFalphaVDW(ind)/reaxFFr0VDW(ind)
                    dVDWtrmdrij = dVDWtrmdf13*df13drij
                    detrmdrij = scale*reaxFFDeVDW(ind)*(expVDW2 - expVDW1)*dVDWtrmdrij
                    devdwdrij = detrmdrij*tp + etrm*dtpdr
                  else
                    devdwdrij = 0.0_dp
                  endif
                  if (lreaxFF_ecoul.and..not.lQMMMelectro) then
!
!  Compute derivative of Coulomb energy
!           
                    dgamdrij = - (gam*rij)/(rij*rij2 + gammaij)
                    decouldrij = gam*dtpQdr + dgamdrij*tpQ
                  else
                    decouldrij = 0.0_dp
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
                  xdrv(j) = xdrv(j) + xd          
                  ydrv(j) = ydrv(j) + yd          
                  zdrv(j) = zdrv(j) + zd          
                  if (nregioni.ne.nregionj) then
                    xregdrv(nregioni) = xregdrv(nregioni) - xd
                    yregdrv(nregioni) = yregdrv(nregioni) - yd
                    zregdrv(nregioni) = zregdrv(nregioni) - zd
                    xregdrv(nregionj) = xregdrv(nregionj) + xd          
                    yregdrv(nregionj) = yregdrv(nregionj) + yd
                    zregdrv(nregionj) = zregdrv(nregionj) + zd
                  endif
                  if (lstr) then
                    call real1strterm(ndim,xji,yji,zji,0.0_dp,0.0_dp,0.0_dp,dr2ds(1,1), &
                      d2r2dx2(1,1),d2r2dsdx(1,1,1),d2r2ds2(1,1,1),.false.)
                    rstrdloc(1:nstrains) = 0.0_dp
                    do kl = 1,nstrains
                      ns = nstrptr(kl)
                      rstrdloc(kl) = rstrdloc(kl) + (devdwdrij + decouldrij)*dr2ds(ns,1)
                    enddo
                    do kl = 1,nstrains 
                      rstrd(kl) = rstrd(kl) + rstrdloc(kl)
                    enddo
                    if (latomicstress) then 
                      do kl = 1,nstrains
                        atomicstress(kl,i) = atomicstress(kl,i) + 0.5_dp*rstrdloc(kl)
                        atomicstress(kl,j) = atomicstress(kl,j) + 0.5_dp*rstrdloc(kl)
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
  if (index(keyword,'verb').ne.0.and.nprocs.gt.1)  then
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
    call sumall(sum0,sum,14_i4,"reaxffmds","energies")
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
  deallocate(nonzeroptr_j2k2,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_j2k2')
  deallocate(nonzeroptr_jk,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_jk')
  deallocate(nonzeroptr_ij2,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_ij2')
  deallocate(nonzeroptr_ik,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_ik')
  deallocate(nonzeroptr_ij,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr_ij')
  deallocate(nonzeroptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nonzeroptr')
  deallocate(d1BOjk2_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOjk2_pipi')
  deallocate(d1BOjk2_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOjk2_pi')
  deallocate(d1BOjk2_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOjk2_s')
  deallocate(d1BOjk2,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOjk2')
  deallocate(d1BOjk_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOjk_pipi')
  deallocate(d1BOjk_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOjk_pi')
  deallocate(d1BOjk_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOjk_s')
  deallocate(d1BOjk,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOjk')
!
  deallocate(dsboproddbo,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','dsboproddbo')
  deallocate(d1BOik_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOik_pipi')
  deallocate(d1BOik_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOik_pi')
  deallocate(d1BOik_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOik_s')
  deallocate(d1BOik,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOik')
  deallocate(d1BOij2_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOij2_pipi')
  deallocate(d1BOij2_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOij2_pi')
  deallocate(d1BOij2_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOij2_s')
  deallocate(d1BOij2,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOij2')
  deallocate(d1BOi_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOi_pipi')
  deallocate(d1BOi_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOi_pi')
  deallocate(d1BOi_s,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOi_s')
  deallocate(d1BOi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOi')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1i')
  deallocate(devaldsbo_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','devaldsbo_neigh')
  deallocate(devalddeltaj_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','devalddeltaj_neigh')
  deallocate(detorsdBOpi_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','detorsdBOpi_neigh')
  deallocate(devaldBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','devaldBO_neigh')
  deallocate(debpendBO_neigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','debpendBO_neigh')
  deallocate(d1BOp_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOp_pipi')
  deallocate(d1BOp_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOp_pi')
  deallocate(d1BOp,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','d1BOp')
  deallocate(deltalplp,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','deltalplp')
  deallocate(deltalp,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','deltalp')
  deallocate(deltap,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','deltap')
  deallocate(delta,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','delta')
  deallocate(BOp_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','BOp_pipi')
  deallocate(BOp_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','BOp_pi')
  deallocate(BOp,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','BOp')
  deallocate(BO_pipi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','BO_pipi')
  deallocate(BO_pi,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','BO_pi')
  deallocate(BO,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','BO')
  deallocate(zneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','zneigh')
  deallocate(yneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','yneigh')
  deallocate(xneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','xneigh')
  deallocate(rneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','rneigh')
  deallocate(ineigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','ineigh')
  deallocate(neighnoRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','neighnoRptr')
  deallocate(neighno,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','neighno')
  deallocate(sum0,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','sum0')
  deallocate(sum,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','sum')
  deallocate(sumsij,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','sumsij')
  deallocate(nneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nneigh')
  deallocate(nfreeatom,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nfreeatom')
  deallocate(lopanyneigh,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','lopanyneigh')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','latomdone')
  deallocate(nbosptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nbosptr')
  deallocate(nbos,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nbos')
  deallocate(numattodoRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','numattodoRptr')
  deallocate(numattodonptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','numattodonptr')
  deallocate(numattodoptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','numattodoptr')
  deallocate(nbodummyptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nbodummyptr')
  deallocate(nboatomRptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nboatomRptr')
  deallocate(nboatomptr,stat=status)
  if (status/=0) call deallocate_error('reaxFFmds','nboatomptr')
!
  t2 = g_cpu_time()
  treaxFF = treaxFF + t2 - t1 - (t4 - t3)
  teem    = teem + t4 - t3
#ifdef TRACE
  call trace_out('reaxffmds')
#endif
!
  return
  end
