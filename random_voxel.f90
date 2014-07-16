!=========================================================================================
module voxel

  use global_parameters
  implicit none
  private

  public :: generate_random_points_vox
  public :: get_volume_vox

contains

!-----------------------------------------------------------------------------------------
function generate_random_points_vox(v, n)
!
!  Returns uniform points in a voxel (i.e. non-conforming blocks as defined
!  in Boschi et al. 1998). v1=x, v2=y, v3=z
!

  integer, intent(in)        ::  n ! number of random points
  real(kind=dp), intent(in)  ::  v(3,8) ! the 8 vertices of the input element
  real(kind=dp)              ::  s(3,8) ! the 8 vertices in spherical coordinates
  real(kind=dp)              ::  generate_random_points_vox(3,n) ! out in cart. coord.
  real(kind=dp)              ::  rvec(3,n)

  real(kind=dp)              ::  r_tmp, r_min, r_max
  real(kind=dp)              ::  ph_tmp, ph_min, ph_max
  real(kind=dp)              ::  th_tmp, th_min, th_max
  integer                    ::  i

  ! Generate a matrix with random numbers
  call random_number(rvec)

  ! Convert cartesian to spherical coordinates
  s(1,:)=sqrt ( v(1,:)**2 + v(2,:)**2 + v(3,:)**2 )    ! radius r 
  s(2,:)=acos ( v(3,:) / s(1,:) )                    ! latitude (phi)
  s(3,:)=atan2 ( v(2,:), v(1,:) )                   ! longitude (theta)
     
  ! Determine minimum and maximum r, phi and theta (in radians)
  r_min=minval(s(1,:))
  r_max=maxval(s(1,:))
  ph_min=minval(s(2,:))
  ph_max=maxval(s(2,:))
  th_min=minval(s(3,:))
  th_max=maxval(s(3,:))

  ! Loop over number of random points
  do i=1,n
     r_tmp  = (rvec(1,i)*(r_max**3-r_min**3)+r_min**3)**(1./3.)     
     ph_tmp = acos(cos(ph_min)+(cos(ph_max)-cos(ph_min))*rvec(2,i))
     th_tmp = th_min+abs(th_max-th_min)*rvec(3,i)

     ! Convert back to cartesian
     generate_random_points_vox(1,i)=r_tmp*cos(th_tmp)*sin(ph_tmp)
     generate_random_points_vox(2,i)=r_tmp*sin(th_tmp)*sin(ph_tmp)
     generate_random_points_vox(3,i)=r_tmp*cos(ph_tmp)
  end do

 

end function generate_random_points_vox
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function get_volume_vox(v)
!
! Computes the volume of a voxel
!

  real(kind=dp), intent(in)  ::  v(3,8) ! voxel in sph. coord
  real(kind=dp)              ::  s(3,8) ! the 8 vertices in spherical coordinates
  real(kind=dp)              ::  get_volume_vox

  real(kind=dp)              ::  vol_shell
  real(kind=dp)              ::  frac,area

  real(kind=dp)              ::  r_min, r_max
  real(kind=dp)              ::  ph_min, ph_max
  real(kind=dp)              ::  th_min, th_max


  ! ! Convert cartesian to spherical coordinates
  ! s(1,:)=sqrt ( v(1,:)**2 + v(2,:)**2 + v(3,:)**2 )    ! radius r 
  ! s(2,:)=acos ( v(3,:) / s(1,:) )                    ! latitude (phi)
  ! s(3,:)=atan2 ( v(2,:), v(1,:) )                   ! longitude (theta)

  ! radius r
  s(1,:)=sqrt(v(1,:)**2+v(2,:)**2+v(3,:)**2)     
  ! latitude (phi)
  s(2,:)=asin(v(3,:)/s(1,:))
  ! longitude (theta)
  s(3,:)=atan2(v(2,:),v(1,:))

  ! determine minimum and maximum r, phi and theta (in radians)
  r_min=minval(s(1,:))
  r_max=maxval(s(1,:))
  ph_min=minval(s(2,:))
  ph_max=maxval(s(2,:))
  th_min=minval(s(3,:))
  th_max=maxval(s(3,:))

  vol_shell = (4/3) * pi * ( r_max**3 - r_min**3 )
  frac = abs(th_min-th_max) * abs(sin(ph_min)-sin(ph_max)) / (4*pi);
  get_volume_vox = vol_shell * frac

  

end function get_volume_vox
!-----------------------------------------------------------------------------------------




end module
!=========================================================================================
