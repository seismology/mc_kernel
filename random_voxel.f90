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
!  in Boschi et al. 1998).
!
!  http://mathworld.wolfram.com/SphericalCoordinates.html
!
!

  integer, intent(in)        ::  n           ! number of random points
  real(kind=dp), intent(in)  ::  v(3,8)      ! the 8 vertices of the input element
  real(kind=dp)              ::  generate_random_points_vox(3,n) ! output (cartesian)
  real(kind=dp)              ::  rvec(3,n)   ! random numbers
  real(kind=dp)              ::  sph_rng(6)  ! spherical range: rmin,rmax,phmin,phmax,thmin,thmax
  real(kind=dp)              ::  sph_pnt(3)  ! spherical point: r,ph,th
  real(kind=dp)              ::  u(3)        ! cartesian point: x,y,z
  integer                    ::  i

  ! Generate a matrix with random numbers
  call random_number(rvec)
  call cartesian_to_spherical_range(v,sph_rng)

  ! Loop over number of random points
  do i=1,n

     ! Trigonometric method to find uniformly distributed
     ! random points in a spherical segment
     sph_pnt(1) = (rvec(1,i) * (sph_rng(2)**3-sph_rng(1)**3) &
                  + sph_rng(1)**3)**(1./3.)     
     sph_pnt(2) = acos(cos(sph_rng(3))+(cos(sph_rng(4))& 
                  - cos(sph_rng(3)))*rvec(2,i))
     sph_pnt(3) = sph_rng(5)+abs(sph_rng(6)-sph_rng(5))&
                  * rvec(3,i)

		  

     ! Convert current pont to cartesian
     call spherical_to_cartesian_point(sph_pnt,u)

     ! Return to cartesian coordinates
     generate_random_points_vox(:,i)=u(:)

  end do

 

end function generate_random_points_vox
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function get_volume_vox(v)
!
! Computes the volume of a voxel, by computing the surface area fraction
! of a lat-lon quadrilateral and multiplying with volume of spherical shell
!

  real(kind=dp), intent(in)  ::  v(3,8) ! voxel in sph. coord
  real(kind=dp)              ::  sph_rng(6) ! the 8 vertices in spherical coordinates
  real(kind=dp)              ::  get_volume_vox

  real(kind=dp)              ::  vol_shell
  real(kind=dp)              ::  frac,area

  ! Convert cartesian vertices to spherical range
  call cartesian_to_spherical_range(v,sph_rng)

  ! We need latitude instead of colatitude
  sph_rng(3)=pi/2.-sph_rng(3)
  sph_rng(4)=pi/2.-sph_rng(4)

  ! Compute volume of spherical shell and latitude-longitude quadrilateral
  ! The volume can be computed by multiplying the both
  vol_shell = (4/3) * pi * ( sph_rng(2)**3 &
              - sph_rng(1)**3 )

  frac = abs(sph_rng(5)-sph_rng(6)) &
       * abs(sin(sph_rng(3))-sin(sph_rng(4))) &
       / (4*pi);

  get_volume_vox = vol_shell * frac

end function get_volume_vox
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
pure subroutine cartesian_to_spherical_range ( v, spherical_range )
!
! Converts a voxel given in terms of its 8 vertices in cartesian
! coordinates, to a spherical range, given as minimum and maximum
! colatitude, longitude and radius
!
  implicit none

  real(kind=dp), intent(in)  :: v(3,8) !< The voxel vertices
  real(kind=dp), intent(out) :: spherical_range(6)   !< spherical range output
  real(kind=dp)              :: s(3,8) !< the 8 vertices in spherical coordinates

  ! Convert cartesian to spherical coordinates
  s(1,:) = sqrt   ( v(1,:)**2 + v(2,:)**2 + v(3,:)**2 )  ! radius r 
  s(2,:) = acos   ( v(3,:) / s(1,:) )                    ! colatitude (phi)
                                                         ! latitude = 90 - colatitude
  s(3,:) = atan2  ( v(2,:), v(1,:) )                     ! longitude (theta)
     
  ! Determine minimum and maximum r, phi and theta (in radians)
  spherical_range(1) = minval ( s(1,:) ) ! r_min
  spherical_range(2) = maxval ( s(1,:) ) ! r_max
  spherical_range(3) = minval ( s(2,:) ) ! ph_min
  spherical_range(4) = maxval ( s(2,:) ) ! ph_max
  spherical_range(5) = minval ( s(3,:) ) ! th_min
  spherical_range(6) = maxval ( s(3,:) ) ! th_max

end subroutine cartesian_to_spherical_range
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
pure subroutine spherical_to_cartesian_point ( v, u )
!
! Converts a point in spherical coordinates to a
! point in cartesian coordinates
!
  implicit none

  real(kind=dp), intent(in)  :: v(3) !< Input in spherical coordinates
  real(kind=dp), intent(out) :: u(3) !< Output in cartesian coordinates

  u(1)=v(1)*cos(v(3))*sin(v(2))
  u(2)=v(1)*sin(v(3))*sin(v(2))
  u(3)=v(1)*cos(v(2))

end subroutine spherical_to_cartesian_point
!-----------------------------------------------------------------------------------------



end module
!=========================================================================================
