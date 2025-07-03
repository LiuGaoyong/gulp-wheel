  subroutine electricfield(efield,lgrad1,lgrad2)
!
!  Subroutine for calculating the energy / forces due to an external
!  electric field. Note that this can only be done in non-periodic
!  directions and that the energy is defined as the integral of the
!  field relative to the initial geometry.
!
!   3/07 Created from force.f
!   6/09 Site energy added
!  11/09 Region derivatives added
!   6/11 Modified to allow for periodic electric field
!   8/11 Possible symmetry removed as electric field is incompatible
!        with this.
!   9/11 Derivatives of variable charge distribution added
!   9/11 Check for ReaxFF incompatibility added
!  10/11 1-D and 2-D direction settings added in mixed fractional-Cartesian
!  11/11 Region-region energy contributions stored
!   2/12 Addressing of nregionno corrected to refer to asymmetric unit
!   4/12 xvir, yvir and zvir removed
!  12/12 Time-dependent field added
!  12/12 Delay and end added for time-dependent field
!  12/12 Modified to allow for multiple time-dependent fields
!   2/18 Trace added
!   5/18 Parallel use of dqdxyz corrected
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
!   7/23 nregionnocfg added
!   8/23 Use of initial coordinates changed for offsets
!   8/23 Charge derivatives removed from first derivatives as these 
!        are not needed due to Hellmann-Feynman
!   8/23 Second derivatives added
!   9/23 Electric field strain first derivatives added
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
!  Julian Gale, CIC, Curtin University, August 2023
!
  use g_constants,    only : pi
  use control,        only : lreaxFF, leem, lgfnff
  use current
  use derivatives,    only : xdrv, ydrv, zdrv, rstrd
  use derivatives,    only : xregdrv, yregdrv, zregdrv
  use derivatives,    only : dqdxyz, derv2
  use energies,       only : siteenergy, eregion2region
  use field
  use general,        only : timesofar
  use m_strain,       only : cartstrterm
  use mdlogic,        only : lmd
  use moldyn,         only : tmdfieldstart, tmdfieldstop
  use parallel,       only : nprocs, procid, node2atom, natomsonnode, atom2node, atom2local
  use symmetry,       only : lstr
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  logical,                       intent(in)    :: lgrad1
  logical,                       intent(in)    :: lgrad2
  real(dp),                      intent(out)   :: efield
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: iloc
  integer(i4)                                  :: indj
  integer(i4)                                  :: ix
  integer(i4)                                  :: iy
  integer(i4)                                  :: iz
  integer(i4)                                  :: ixf
  integer(i4)                                  :: iyf
  integer(i4)                                  :: izf
  integer(i4)                                  :: j
  integer(i4)                                  :: kl
  integer(i4)                                  :: nregioni
  real(dp)                                     :: dxyzds(6,3)
  real(dp)                                     :: d2xyzdsdx(6,3,3)
  real(dp)                                     :: d2xyzds2(6,6,3)
  real(dp)                                     :: esum
  real(dp)                                     :: fieldx
  real(dp)                                     :: fieldy
  real(dp)                                     :: fieldz
  real(dp)                                     :: fnorm
  real(dp)                                     :: tdf
  real(dp)                                     :: tsf
  real(dp)                                     :: twopi
!
!  Initialise integral of force x distance
!
  efield = 0.0_dp
!
!  If number of atoms is zero then return as there is nothing to do
!
  if (numat.eq.0) return
#ifdef TRACE
  call trace_in('electricfield')
#endif
!
!  Static field
!
  if (lfieldcfg(ncf)) then
!
!  For 3-D convert the field direction from cell vector to Cartesian
!
    if (ndim.eq.3) then
      fieldx = fielddirectioncfg(1,ncf)*rv(1,1) + fielddirectioncfg(2,ncf)*rv(1,2) + fielddirectioncfg(3,ncf)*rv(1,3)
      fieldy = fielddirectioncfg(1,ncf)*rv(2,1) + fielddirectioncfg(2,ncf)*rv(2,2) + fielddirectioncfg(3,ncf)*rv(2,3)
      fieldz = fielddirectioncfg(1,ncf)*rv(3,1) + fielddirectioncfg(2,ncf)*rv(3,2) + fielddirectioncfg(3,ncf)*rv(3,3)
      fnorm = fieldx**2 + fieldy**2 + fieldz**2
      fnorm = 1.0_dp/sqrt(fnorm)
      fieldx = fieldcfg(ncf)*fieldx*fnorm
      fieldy = fieldcfg(ncf)*fieldy*fnorm
      fieldz = fieldcfg(ncf)*fieldz*fnorm
    elseif (ndim.eq.2) then
      fieldx = fielddirectioncfg(1,ncf)*rv(1,1) + fielddirectioncfg(2,ncf)*rv(1,2)
      fieldy = fielddirectioncfg(1,ncf)*rv(2,1) + fielddirectioncfg(2,ncf)*rv(2,2)
      fieldz = fielddirectioncfg(3,ncf)
      fnorm = fieldx**2 + fieldy**2 + fieldz**2
      fnorm = 1.0_dp/sqrt(fnorm)
      fieldx = fieldcfg(ncf)*fieldx*fnorm
      fieldy = fieldcfg(ncf)*fieldy*fnorm
      fieldz = fieldcfg(ncf)*fieldz*fnorm
    elseif (ndim.eq.1) then
      fieldx = fielddirectioncfg(1,ncf)*rv(1,1)
      fieldy = fielddirectioncfg(2,ncf)
      fieldz = fielddirectioncfg(3,ncf)
      fnorm = fieldx**2 + fieldy**2 + fieldz**2
      fnorm = 1.0_dp/sqrt(fnorm)
      fieldx = fieldcfg(ncf)*fieldx*fnorm
      fieldy = fieldcfg(ncf)*fieldy*fnorm
      fieldz = fieldcfg(ncf)*fieldz*fnorm
    else
!
!  Find norm of field direction and scale components
!
      fnorm = fielddirectioncfg(1,ncf)**2 + fielddirectioncfg(2,ncf)**2 + fielddirectioncfg(3,ncf)**2
      fnorm = 1.0_dp/sqrt(fnorm)
      fieldx = fieldcfg(ncf)*fielddirectioncfg(1,ncf)*fnorm
      fieldy = fieldcfg(ncf)*fielddirectioncfg(2,ncf)*fnorm
      fieldz = fieldcfg(ncf)*fielddirectioncfg(3,ncf)*fnorm
    endif
  else
    fieldx = 0.0_dp
    fieldy = 0.0_dp
    fieldz = 0.0_dp
  endif
!
!  Time-dependent electric field
!
  if (lmd.and.(ntdfieldcfg(ncf).gt.0)) then
    tsf = timesofar - tmdfieldstart(ncf)
    if (tsf.gt.0.0_dp.and.(tsf.lt.tmdfieldstop(ncf).or.tmdfieldstop(ncf).eq.0.0_dp)) then
      twopi = 2.0_dp*pi
      do i = 1,ntdfieldcfg(ncf)
        tdf = td_fieldcfg(1,i,ncf)*cos(twopi*(tsf*td_fieldcfg(2,i,ncf) + td_fieldcfg(3,i,ncf)))
        fnorm = td_fielddirectioncfg(1,i,ncf)**2 + td_fielddirectioncfg(2,i,ncf)**2 + td_fielddirectioncfg(3,i,ncf)**2
        fnorm = 1.0_dp/sqrt(fnorm)
        fieldx = fieldx + tdf*td_fielddirectioncfg(1,i,ncf)*fnorm
        fieldy = fieldy + tdf*td_fielddirectioncfg(2,i,ncf)*fnorm
        fieldz = fieldz + tdf*td_fielddirectioncfg(3,i,ncf)*fnorm
      enddo
    endif
  endif
!
!  Loop over asymmetric unit performing integral
!
  do i = 1+procid,numat,nprocs
    esum = qa(i)*(fieldx*(xalat(i) - xvoffset) + &
                  fieldy*(yalat(i) - yvoffset) + &
                  fieldz*(zalat(i) - zvoffset))
    efield = efield - esum
!
    nregioni = nregionno(nrelf2a(i))
    eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) - esum
    siteenergy(i) = siteenergy(i) - esum
  enddo
!
  if (lgrad1) then
    do iloc = 1,natomsonnode
      i = node2atom(iloc)
      xdrv(i) = xdrv(i) - fieldx*qa(i)
      ydrv(i) = ydrv(i) - fieldy*qa(i)
      zdrv(i) = zdrv(i) - fieldz*qa(i)
!
      nregioni = nregionno(nrelf2a(i))
      xregdrv(nregioni) = xregdrv(nregioni) - fieldx*qa(i)
      yregdrv(nregioni) = yregdrv(nregioni) - fieldy*qa(i)
      zregdrv(nregioni) = zregdrv(nregioni) - fieldz*qa(i)
    enddo
!
!  Compute strain derivatives, though there may be issues have a field in a periodic system
!
    if (lstr) then
      do iloc = 1,natomsonnode
        i = node2atom(iloc)
        call cartstrterm(ndim,xalat(i),yalat(i),zalat(i),0.0_dp,0.0_dp,0.0_dp,dxyzds, &
                         d2xyzdsdx,d2xyzds2,.false.)
        do kl = 1,nstrains
          rstrd(kl) = rstrd(kl) - qa(i)*(fieldx*dxyzds(kl,1) + &
                                         fieldy*dxyzds(kl,2) + & 
                                         fieldz*dxyzds(kl,3))
        enddo
      enddo
    endif
!
    if (lgrad2.and.(leem.or.lreaxff.or.lgfnff)) then
!
!  Variable charge second derivatives with field
!
      if (nprocs.gt.1) then
        do iloc = 1,natomsonnode
          i = node2atom(iloc)
!  
          ix = 3*(iloc-1) + 1
          iy = ix + 1
          iz = ix + 2
!
          ixf = 3*(i-1) + 1
          iyf = ixf + 1
          izf = ixf + 2
!       
          indj = 0
          do j = 1,numat
            derv2(indj+1,ix) = derv2(indj+1,ix) - fieldx*dqdxyz(indj+1,i) - fieldx*dqdxyz(ixf,j)
            derv2(indj+2,ix) = derv2(indj+2,ix) - fieldx*dqdxyz(indj+2,i) - fieldy*dqdxyz(ixf,j)
            derv2(indj+3,ix) = derv2(indj+3,ix) - fieldx*dqdxyz(indj+3,i) - fieldz*dqdxyz(ixf,j)
            derv2(indj+1,iy) = derv2(indj+1,iy) - fieldy*dqdxyz(indj+1,i) - fieldx*dqdxyz(iyf,j)
            derv2(indj+2,iy) = derv2(indj+2,iy) - fieldy*dqdxyz(indj+2,i) - fieldy*dqdxyz(iyf,j)
            derv2(indj+3,iy) = derv2(indj+3,iy) - fieldy*dqdxyz(indj+3,i) - fieldz*dqdxyz(iyf,j)
            derv2(indj+1,iz) = derv2(indj+1,iz) - fieldz*dqdxyz(indj+1,i) - fieldx*dqdxyz(izf,j)
            derv2(indj+2,iz) = derv2(indj+2,iz) - fieldz*dqdxyz(indj+2,i) - fieldy*dqdxyz(izf,j)
            derv2(indj+3,iz) = derv2(indj+3,iz) - fieldz*dqdxyz(indj+3,i) - fieldz*dqdxyz(izf,j)
            indj = indj + 3
          enddo
        enddo
      else
        do i = 1,numat
          ix = 3*(i-1) + 1
          iy = ix + 1
          iz = ix + 2
!
          indj = 0
          do j = 1,numat
            derv2(indj+1,ix) = derv2(indj+1,ix) - fieldx*dqdxyz(indj+1,i) - fieldx*dqdxyz(ix,j)
            derv2(indj+2,ix) = derv2(indj+2,ix) - fieldx*dqdxyz(indj+2,i) - fieldy*dqdxyz(ix,j)
            derv2(indj+3,ix) = derv2(indj+3,ix) - fieldx*dqdxyz(indj+3,i) - fieldz*dqdxyz(ix,j)
            derv2(indj+1,iy) = derv2(indj+1,iy) - fieldy*dqdxyz(indj+1,i) - fieldx*dqdxyz(iy,j)
            derv2(indj+2,iy) = derv2(indj+2,iy) - fieldy*dqdxyz(indj+2,i) - fieldy*dqdxyz(iy,j)
            derv2(indj+3,iy) = derv2(indj+3,iy) - fieldy*dqdxyz(indj+3,i) - fieldz*dqdxyz(iy,j)
            derv2(indj+1,iz) = derv2(indj+1,iz) - fieldz*dqdxyz(indj+1,i) - fieldx*dqdxyz(iz,j)
            derv2(indj+2,iz) = derv2(indj+2,iz) - fieldz*dqdxyz(indj+2,i) - fieldy*dqdxyz(iz,j)
            derv2(indj+3,iz) = derv2(indj+3,iz) - fieldz*dqdxyz(indj+3,i) - fieldz*dqdxyz(iz,j)
            indj = indj + 3
          enddo
        enddo
      endif
    endif
  endif
#ifdef TRACE
  call trace_out('electricfield')
#endif
!
  return
  end
