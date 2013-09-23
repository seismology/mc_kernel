!====================
module coord_mapping
!====================

  !use data_mesh
  use global_parameters
  !use data_rota
  !use coord_trafo

  implicit none


  public :: rotate_frame_rd
  public :: compute_srdphirdzrd

  private

  contains 
!!!!-------------------------------------------------------!!!
subroutine rotate_frame_rd(npts,srd,phird,zrd,xgd,ygd,zgd,phigr,thetagr)

    implicit none
    integer, intent(in) :: npts
    !< Number of points to rotate
    real(kind=realkind), dimension(npts), intent(in) :: xgd,ygd,zgd 
    !< Coordinates to rotate (in x, y, z)
    real(kind=realkind), intent(in) :: phigr,thetagr
    !< Rotation angles phi and theta
    real(kind=realkind), dimension(npts), intent(out) :: srd,zrd,phird
    !< Rotated coordinates (in s, z, phi)
    real(kind=realkind), dimension(npts) :: xp,yp,zp,xp_cp,yp_cp,zp_cp
    real(kind=realkind) :: phi_cp,rgd,thetagd
    integer :: i
    character(len=55) :: filename
    double precision :: thgr_dble,phgr_dble

    thgr_dble=dble(thetagr)
    phgr_dble=dble(phigr)

    !!first rotation (longitude)
    xp_cp=xgd*dcos(phgr_dble)+ygd*dsin(phgr_dble)
    yp_cp=-xgd*dsin(phgr_dble)+ygd*dcos(phgr_dble)
    zp_cp=zgd

    !second rotation (colat)
    xp=xp_cp*dcos(thgr_dble)-zp_cp*dsin(thgr_dble)
    yp=yp_cp
    zp=xp_cp*dsin(thgr_dble)+zp_cp*dcos(thgr_dble)

    srd=dsqrt(dble(xp)**2+dble(yp)**2)
    zrd=zp
    do i=1,npts
       phi_cp=atan2(yp(i),xp(i))
       if (phi_cp<0.d0) then
          phird(i)=2.d0*pi+phi_cp
       else
          phird(i)=phi_cp
       endif
       if (phgr_dble==0.0 .and. ygd(i)==0.0)  phird(i)=0.
    enddo

    write(6,*)'Done with rotating frame rd.'

end subroutine rotate_frame_rd

!-----------------------------------------------------------------------
subroutine compute_srdphirdzrd(srd,phird,zrd,xgd,ygd,zgd, &
                               xorig,yorig,zorig,phigr,thetagr)
  use coord_trafo, only: xyz2rthetaphi
  real(kind=realkind), intent(out)::  srd,phird,zrd
  real(kind=realkind), intent(in) :: xgd,ygd,zgd
  real(kind=realkind), intent(in) :: xorig,yorig,zorig
  real(kind=realkind), intent(in) :: phigr,thetagr
  real(kind=realkind) :: rgd,phigd,thetagd,alpha
  real(kind=realkind) :: arr

!  thetagd and phigd as computed from
!  xgd(ipt),ygd(ipt),zgd(ipt), using the xyz2rthetaphi routine. 
  call xyz2rthetaphi(rgd,thetagd,phigd,xgd,ygd,zgd)

!  In the case of a source at the NP .AND. the z-component of the receiver, 
!  things are fine.

!  find closest gll point for that specific diffracting point in the D domain
!  for that specific diffracting point

  arr = (xorig**2+yorig**2+zorig**2)**(-.5)

  srd  =arr*( (ygd*zorig-zgd*yorig)**2 + &
             (zgd*xorig-xgd*zorig)**2 + &
             (xgd*yorig-ygd*xorig)**2     )**(.5)
  zrd = arr*(xgd*xorig+ygd*yorig+zgd*zorig)

! phird
  alpha = abs(phigr-phigd)
  phird = atan(sin(alpha)/(sin(thetagr)/tan(thetagd) -cos(thetagr)*cos(alpha)))
  if (abs(phird)<smallval_dble) phird = 0.d0

end subroutine compute_srdphirdzrd

end module coord_mapping
