!******************************************************************************
!
!    This file is part of:
!    MC Kernel: Calculating seismic sensitivity kernels on unstructured meshes
!    Copyright (C) 2016 Simon Staehler, Martin van Driel, Ludwig Auer
!
!    You can find the latest version of the software at:
!    <https://www.github.com/tomography/mckernel>
!
!    MC Kernel is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MC Kernel is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MC Kernel. If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

!=========================================================================================
module polytopes

  use global_parameters
  implicit none
  private

  public :: generate_random_points_tet, generate_random_points_poly
  public :: generate_random_points_ref_tri
  public :: get_volume_tet, get_volume_poly, get_center_tet, get_center_tri
  public :: determinant_4, point_in_triangle, point_in_triangle_3d
  public :: point_in_tetrahedron

  interface determinant_4
     module procedure  :: determinant_4_dp
     module procedure  :: determinant_4_qp
  end interface determinant_4


contains

!-----------------------------------------------------------------------------------------
function generate_random_points_tet(v, npoints, quasirandom)
!  returns uniform points in a tetrahedron.
!  Adapted from J. Burkardt's code collection, referencing:
!
!    Claudio Rocchini, Paolo Cignoni,
!    Generating Random Points in a Tetrahedron,
!    Journal of Graphics Tools,
!    Volume 5, Number 5, 2000, pages 9-12.

  use halton_sequence, only: get_halton
  use simple_routines, only: check_limits

  integer, intent(in)           ::  npoints
  real(kind=dp), intent(in)     ::  v(3,4)
  logical, intent(in), optional ::  quasirandom
  real(kind=dp)                 ::  generate_random_points_tet(3,npoints)

  logical                       ::  invalid_random_number
  real(kind=dp)                 ::  c(4), t
  real(kind=dp)                 ::  coordinates(3,npoints)
  integer                       ::  ipoint

  if (present(quasirandom)) then
    if (quasirandom) then
      call get_halton(coordinates)
    else
      call random_number(coordinates)
    end if
  else
    call random_number(coordinates)
  end if

  invalid_random_number = check_limits(real(coordinates, kind=sp),   &
                                       limits = [0.0, 1.0],          &
                                       array_name = 'random_number')

  if (invalid_random_number) then
    print *, 'ERROR: random number generator returned number outside range [0,1]'
    print *, '       quasirandom: ', quasirandom
  end if

  do ipoint = 1, npoints

     c(1:3) = coordinates(:, ipoint)

     if (sum(c(1:2)) > 1.0d0) then
        c(1) = 1.0d0 - c(1)
        c(2) = 1.0d0 - c(2)
     end if

     if (sum(c(2:3)) > 1.0d0) then
        t = c(3)
        c(3) = 1.0d0 - c(1) - c(2)
        c(2) = 1.0d0 - t
     else if (sum(c(1:3)) > 1.0d0) then
        t = c(3)
        c(3) = sum(c(1:3)) - 1.0d0
        c(1) = 1.0d0 - c(2) - t
     end if

     c(4) = 1.0d0 - sum(c(1:3))
  
     ! c(1:4) are the barycentric coordinates of the point.
  
     generate_random_points_tet(1:3, ipoint) = matmul(v(1:3, 1:4), c(1:4))

  end do

end function generate_random_points_tet
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function generate_random_points_ref_tri( npoints, quasirandom ) result(x)
!< Generate random points in the reference triangle (0,0), (0,1), (1,0)

  use halton_sequence, only: get_halton
  use simple_routines, only: check_limits

  integer, parameter            :: dim_num = 2
  integer, intent(in)           :: npoints  ! Number of points
  logical, intent(in), optional :: quasirandom

  logical                       ::  invalid_random_number
  real(kind=dp)                 :: x(dim_num, npoints)
  real(kind=dp)                 :: r(dim_num, npoints)
  integer                       :: ipoint

  if (present(quasirandom)) then
    if (quasirandom) then
      call get_halton(r)
    else
      call random_number(r)
    end if
  else
    call random_number(r)
  end if

  invalid_random_number = check_limits(real(r, kind=sp),             &
                                       limits = [0.0, 1.0],          &
                                       array_name = 'random_number')

  if (invalid_random_number) then
    print *, 'ERROR: random number generator returned number outside range [0,1]'
    print *, '       quasirandom: ', quasirandom
  end if


  do ipoint = 1, npoints

    if ( 1 < sum (r(:,ipoint)) ) then
      x(:, ipoint) = 1 - (r(:,ipoint)*0.999+0.0005)
    else
      x(:, ipoint) = r(:,ipoint)*0.999+0.0005
    end if
  end do
end function generate_random_points_ref_tri
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function generate_random_points_poly( nv, v, n, quasirandom ) result(x)
!! UNIFORM_IN_POLYGON_MAP maps uniform points into a polygon.
!
!  Discussion:
!
!    If the polygon is regular, or convex, or at least star-shaped,
!    this routine will work.
!
!    This routine assumes that all points between the centroid and
!    any point on the boundary lie within the polygon.
!
!    This routine is valid for spatial dimension DIM_NUM = 2.
!
  use halton_sequence, only: get_halton
  use simple_routines, only: check_limits

  integer, parameter            :: dim_num = 2
  integer, intent(in)           :: n  ! Number of points
  integer, intent(in)           :: nv ! Degree of polynomial
  real(kind=dp), intent(in)     :: v(dim_num, nv)
  logical, intent(in), optional :: quasirandom

  real(kind=dp)                 :: x(dim_num, n)
  logical                       ::  invalid_random_number

  real(kind=dp)                 :: area(nv)
  real(kind=dp)                 :: area_norm(nv)
  real(kind=dp)                 :: area_int(nv)
  real(kind=dp)                 :: area_percent
  real(kind=dp)                 :: centroid(dim_num)
  integer                       :: i, ip1, j, k
  real(kind=dp)                 :: r(2, n)
  real(kind=dp)                 :: t(dim_num,3)

  if (present(quasirandom)) then
    if (quasirandom) then
      call get_halton(r)
    else
      call random_number(r)
    end if
  else
    call random_number(r)
  end if

  invalid_random_number = check_limits(real(r, kind=sp),             &
                                       limits = [0.0, 1.0],          &
                                       array_name = 'random_number')

  if (invalid_random_number) then
    print *, 'ERROR: random number generator returned number outside range [0,1]'
    print *, '       quasirandom: ', quasirandom
  end if


  if (nv == 3) then
    do j = 1, n
      if ( sum(r(:,j)) > 1  ) then
        r(:,j) = 1 - r(:,j)
      end if

      x(:,j) = ( 1 - r(1,j) - r(2,j) ) * v(:,1) &
                   + r(1,j)            * v(:,2) &
                            + r(2,j)   * v(:,3) 
    end do

  else
    !
    !  Find the centroid.
    !
      call polygon_centroid_2d ( nv, dble(v), centroid )
    !
    !  Determine the areas of each triangle.
    !
      do i = 1, nv

        if ( i < nv ) then
          ip1 = i + 1
        else
          ip1 = 1
        end if

        t(1:2,1) = v(1:2,i)
        t(1:2,2) = v(1:2,ip1)
        t(1:2,3) = centroid(1:2)

        call triangle_area_2d ( t, area(i) )

      end do
    !
    !  Normalize the areas.
    !
      area_norm(1:nv) = area(1:nv) / sum ( area(1:nv) )
    !
    !  Replace each area by the sum of itself and all previous ones.
    !
      area_int(1) = area_norm(1)
      do i = 2, nv
        area_int(i) = area_norm(i) + area_int(i-1)
      end do

      do j = 1, n
    !
    !  Choose a triangle at random, based on areas.
    !
        !area_percent = r8_uniform_01 ( seed )
        call random_number(area_percent) 

        do k = 1, nv

          if ( area_percent <= area_int(k) ) then
            i = k
            exit
          end if

        end do
    !
    !  Now choose a point at random in the triangle.
    !
        if ( i < nv ) then
          ip1 = i + 1
        else
          ip1 = 1
        end if

        !call r8vec_uniform_01 ( dim_num, seed, r )

        if ( 1.0D0 < sum (r(:,j)) ) then
          r(:,j) = 1.0D0 - r(:,j)
        end if

        if ((ip1.gt.nv).or.(i.gt.nv)) then
            print *, i, ip1
            print *, area
            print *, area_norm
            print *, area_int
            print *, area_percent
            print *, v(:,1)
            print *, v(:,2)
            print *, v(:,3)
            print *, centroid
        end if
        x(:,j) = ( 1.0D0 - r(1,j) - r(2,j) ) * v(:,i) &
                         + r(1,j)            * v(:,ip1) &
                                  + r(2,j)   * centroid(:)

      end do

  end if ! nv==3

end function generate_random_points_poly
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function determinant_4_dp(a)
! computes the determinant of a 4 by 4 matrix, uses qp internally.
  real(kind=dp), intent(in)  :: a(4,4)
  real(kind=qp)              :: a_qp(4,4)
  real(kind=dp)              :: determinant_4_dp

  a_qp = a

  determinant_4_dp = real(&
      a_qp(1,1) * ( &
        a_qp(2,2) * ( a_qp(3,3) * a_qp(4,4) - a_qp(3,4) * a_qp(4,3) ) &
      - a_qp(2,3) * ( a_qp(3,2) * a_qp(4,4) - a_qp(3,4) * a_qp(4,2) ) &
      + a_qp(2,4) * ( a_qp(3,2) * a_qp(4,3) - a_qp(3,3) * a_qp(4,2) ) ) &
    - a_qp(1,2) * ( &
        a_qp(2,1) * ( a_qp(3,3) * a_qp(4,4) - a_qp(3,4) * a_qp(4,3) ) &
      - a_qp(2,3) * ( a_qp(3,1) * a_qp(4,4) - a_qp(3,4) * a_qp(4,1) ) &
      + a_qp(2,4) * ( a_qp(3,1) * a_qp(4,3) - a_qp(3,3) * a_qp(4,1) ) ) &
    + a_qp(1,3) * ( &
        a_qp(2,1) * ( a_qp(3,2) * a_qp(4,4) - a_qp(3,4) * a_qp(4,2) ) &
      - a_qp(2,2) * ( a_qp(3,1) * a_qp(4,4) - a_qp(3,4) * a_qp(4,1) ) &
      + a_qp(2,4) * ( a_qp(3,1) * a_qp(4,2) - a_qp(3,2) * a_qp(4,1) ) ) &
    - a_qp(1,4) * ( &
        a_qp(2,1) * ( a_qp(3,2) * a_qp(4,3) - a_qp(3,3) * a_qp(4,2) ) &
      - a_qp(2,2) * ( a_qp(3,1) * a_qp(4,3) - a_qp(3,3) * a_qp(4,1) ) &
      + a_qp(2,3) * ( a_qp(3,1) * a_qp(4,2) - a_qp(3,2) * a_qp(4,1) ) ) &
                       , kind = dp)

end function determinant_4_dp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function determinant_4_qp(a)
! computes the determinant of a 4 by 4 matrix with quadruple precision
  real(kind=qp), intent(in)  :: a(4,4)
  real(kind=qp)              :: a_qp(4,4)
  real(kind=qp)              :: determinant_4_qp

  a_qp = a

  determinant_4_qp = &
      a_qp(1,1) * ( &
        a_qp(2,2) * ( a_qp(3,3) * a_qp(4,4) - a_qp(3,4) * a_qp(4,3) ) &
      - a_qp(2,3) * ( a_qp(3,2) * a_qp(4,4) - a_qp(3,4) * a_qp(4,2) ) &
      + a_qp(2,4) * ( a_qp(3,2) * a_qp(4,3) - a_qp(3,3) * a_qp(4,2) ) ) &
    - a_qp(1,2) * ( &
        a_qp(2,1) * ( a_qp(3,3) * a_qp(4,4) - a_qp(3,4) * a_qp(4,3) ) &
      - a_qp(2,3) * ( a_qp(3,1) * a_qp(4,4) - a_qp(3,4) * a_qp(4,1) ) &
      + a_qp(2,4) * ( a_qp(3,1) * a_qp(4,3) - a_qp(3,3) * a_qp(4,1) ) ) &
    + a_qp(1,3) * ( &
        a_qp(2,1) * ( a_qp(3,2) * a_qp(4,4) - a_qp(3,4) * a_qp(4,2) ) &
      - a_qp(2,2) * ( a_qp(3,1) * a_qp(4,4) - a_qp(3,4) * a_qp(4,1) ) &
      + a_qp(2,4) * ( a_qp(3,1) * a_qp(4,2) - a_qp(3,2) * a_qp(4,1) ) ) &
    - a_qp(1,4) * ( &
        a_qp(2,1) * ( a_qp(3,2) * a_qp(4,3) - a_qp(3,3) * a_qp(4,2) ) &
      - a_qp(2,2) * ( a_qp(3,1) * a_qp(4,3) - a_qp(3,3) * a_qp(4,1) ) &
      + a_qp(2,3) * ( a_qp(3,1) * a_qp(4,2) - a_qp(3,2) * a_qp(4,1) ) ) 

end function determinant_4_qp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine polygon_centroid_2d ( n, v, centroid )

!*****************************************************************************80
!
!! POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
!
!  Formula:
!
!    Denoting the centroid coordinates by CENTROID, then
!
!      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
!      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
!
!    Green's theorem states that
!
!      Integral ( Polygon boundary ) ( M dx + N dy ) =
!      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
!
!    Using M = 0 and N = x * x / 2, we get:
!
!      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x * x dy,
!
!    which becomes
!
!      CENTROID(1) = 1/6 Sum ( 1 <= I <= N )
!        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
!
!    where, when I = N, the index "I+1" is replaced by 1.
!
!    A similar calculation gives us a formula for CENTROID(2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerard Bashein, Paul Detmer,
!    Centroid of a Polygon,
!    Graphics Gems IV, edited by Paul Heckbert,
!    AP Professional, 1994.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of sides of the polygonal shape.
!
!    Input, real ( kind = 8 ) V(2,N), the coordinates of the vertices
!    of the shape.
!
!    Output, real ( kind = 8 ) CENTROID(2), the coordinates of the
!    centroid of the shape.
!
  implicit none

  integer, intent(in)        :: n
  real(kind=dp), intent(in)  :: v(2,n)
  real(kind=dp), intent(out) :: centroid(2)

  real(kind=dp)              :: area
  integer                    :: i
  integer                    :: ip1
  real(kind=dp)              :: temp

  area = 0
  centroid(:) = 0

  do i = 1, n

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    temp = ( v(1,i) * v(2,ip1) - v(1,ip1) * v(2,i) )

    area = area + temp

    centroid(:) = centroid(:) + ( v(:,ip1) + v(:,i) ) * temp

  end do

  area = area * 0.5D+00

  centroid(:) = centroid(:) / ( 6.0D+00 * area )

end subroutine polygon_centroid_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine triangle_area_2d ( t, area )
! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
  implicit none

  real(kind=dp), intent(in)  :: t(2,3) !< The triangle vertices
  real(kind=dp), intent(out) :: area   !< area of the triangle

  area = 0.5D+00 * abs ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

end subroutine triangle_area_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function get_volume_tet(v)
! computes the volume of a tetrahedron in 3D.
  real(kind=dp), intent(in)  ::  v(3,4) !< Vertices
  real(kind=dp)              ::  a(4,4)
  real(kind=dp)              ::  get_volume_tet

  a(1:4,1:3) = transpose(v)
  a(1:4,4) = [1, 1, 1, 1]

  get_volume_tet = abs(determinant_4(a)) / 6.0d0

end function get_volume_tet
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Computes the centroid of a tetrahedron
pure function get_center_tet(v)

  real(kind=dp), intent(in) :: v(3,4)
  real(kind=dp)             :: get_center_tet(3)
  integer i

  do i = 1,3
    get_center_tet(i) = sum(v(i,1:4)) / 4.0d0
  end do

end function get_center_tet
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Computes the centroid of a triangle
pure function get_center_tri(v)

  real(kind=dp), intent(in) :: v(3,3)
  real(kind=dp)             :: get_center_tri(3)
  integer i

  do i = 1,3
    get_center_tri(i) = sum(v(i,1:3)) / 3.0d0
  end do

end function get_center_tri
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function get_volume_poly(n, v) result(area)
!! get_volume_poly computes the area of a polygon in 3D.
!
!  Discussion:
!
!    The computation is not valid unless the vertices of the polygon
!    lie in a plane, so that the polygon that is defined is "flat".
!
!    The polygon does not have to be "regular", that is, neither its
!    sides nor its angles need to be equal.
!
  use simple_routines, only  : cross
  implicit none


  integer, parameter        :: dim_num = 3  !< Dimension of the space we're
                                            !! living in
  integer, intent(in)       :: n            !< N, degree of polygon
  real(kind=dp), intent(in) :: v(dim_num,n) !< Vertices of polygon
  
  real(kind=dp)             :: area         !< Output, the area of the polygon
  integer                   :: i
  integer                   :: ip1
  real(kind=dp)             :: normal(dim_num)
  
  normal(:) = 0
  
  do i = 1, n
     if ( i < n ) then
        ip1 = i + 1
     else
        ip1 = 1
     end if
     
     normal(:) = normal(:) + cross(v(:,i), v(:,ip1))
  end do
  
  area = 0.5D+00 * sqrt(sum(normal(:)**2))

end function get_volume_poly
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure function point_in_triangle(r, x)
! Checks whether point x is in the triangle between r1, r2, r3
   real(kind=dp), intent(in)  :: r(2,3)
   real(kind=dp), intent(in)  :: x(:,:)
   real(kind=dp)              :: a, b, c, denominator
   logical                    :: point_in_triangle(size(x,2))
   integer                    :: ipoint, npoints

   npoints = size(x,2)

   do ipoint = 1, npoints
      denominator = ((r(2,2) - r(2,3))*(r(1,1) -      r(1,3)) + (r(1,3) - r(1,2))*(r(2,1)      - r(2,3)))
      a =       dble((r(2,2) - r(2,3))*(x(1,ipoint) - r(1,3)) + (r(1,3) - r(1,2))*(x(2,ipoint) - r(2,3))) / denominator
      b =       dble((r(2,3) - r(2,1))*(x(1,ipoint) - r(1,3)) + (r(1,1) - r(1,3))*(x(2,ipoint) - r(2,3))) / denominator
      c = 1 - a - b

      point_in_triangle(ipoint) = ((0 <= a).and.(a <= 1).and.(0 <= b).and.(b <= 1).and.(0 <= c).and.(c <= 1))
   end do

end function point_in_triangle
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
function point_in_triangle_3d(r, p, isinplane)
  use simple_routines, only  : absreldiff, cross
! Using barycentric coordinates, following 
! http://math.stackexchange.com/questions/4322/check-whether-a-point-is-within-a-3d-triangle
 integer, parameter             :: ndim = 3
 real(kind=dp), intent(in)      :: r(ndim,3)  !< Vertices of the triangle
 real(kind=dp), intent(in)      :: p(:,:)     !< Set of points to check
 logical, intent(out), optional :: isinplane(size(p,2))
 
 real(kind=qp)                  :: a, b, c
 real(kind=dp)                  :: x(3), y(3), k, area
 integer                        :: ipoint, npoints
 logical                        :: point_in_triangle_3d(size(p,2))


 point_in_triangle_3d = .false.
 if (present(isinplane)) isinplane = .false.
 npoints = size(p,2)
 area = norm2(cross( r(:,2)-r(:,1), r(:,3)-r(:,1) )) / 2
 do ipoint = 1, npoints
   ! Step 1:
   ! Check whether point is in the plane spanned by the triangle, by comparing
   ! normal vectors of planes spanned by r1r2 and r1r3 with the one spanned by
   ! Pr2 and Pr3.
   x = cross(r(:,2)-r(:,1),      r(:,3)-r(:,1)      )
   y = cross(r(:,2)-p(:,ipoint), r(:,3)-p(:,ipoint) ) 

   ! Check whether x and y are parallel
   k = norm2(cross(x, y))

   ! If they are not, cycle, and point_in_triangle_3d(ipoint) remains false
   if (abs(k)>1d-10) cycle
   
   if (present(isinplane)) isinplane(ipoint) = .true.

   ! Step 2:
   ! Check whether point is within triangle, based on barycentric coordinates
   ! a,b,c
   a = norm2(cross( r(:,2) - p(:,ipoint), r(:,3) - p(:,ipoint) )) / (2*area)
   b = norm2(cross( r(:,3) - p(:,ipoint), r(:,1) - p(:,ipoint) )) / (2*area)
   c = norm2(cross( r(:,1) - p(:,ipoint), r(:,2) - p(:,ipoint) )) / (2*area) !1 - a - b

   if (any([a,b,c]<=-1d-8).or.any([a,b,c]>1+1d-8).or.(sum([a,b,c])>1+1d-8)) cycle

   point_in_triangle_3d(ipoint) = .true.

 end do

end function point_in_triangle_3d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Tests whether a point is in a tetrahedron, by testing whether it is on the 
! same side of each plane, by testing whether all Di have the same sign.
pure function point_in_tetrahedron(r, p) 
  integer, parameter             :: ndim = 3
  real(kind=dp), intent(in)      :: r(ndim,4)  !< Vertices of the tetrahedron
  real(kind=dp), intent(in)      :: p(:,:)     !< Set of points to check
  logical                        :: point_in_tetrahedron(size(p,2))

  real(kind=qp), dimension(4,4)  :: M0, M1, M2, M3, M4
  real(kind=qp)                  :: D0, D1, D2, D3, D4
  integer                        :: ipoint, npoints


  npoints = size(p,2)
  point_in_tetrahedron = .false.

  M0(1:ndim,:) = r
  M0(4,:)      = 1
  D0           = determinant_4(M0)

  ! Check, whether the tetrahedron is degenerate, i.e. all points are in a plane
  if (abs(D0).lt.epsilon(D0)) then
    point_in_tetrahedron = .false.
    return
  end if

  do ipoint = 1, npoints
    M1(1:ndim,:) = r
    M1(4,:)      = 1
    M1(1:ndim,1) = p(:, ipoint)
    D1           = determinant_4(M1) / D0

    M2(1:ndim,:) = r
    M2(4,:)      = 1
    M2(1:ndim,2) = p(:, ipoint)
    D2           = determinant_4(M2) / D0

    M3(1:ndim,:) = r
    M3(4,:)      = 1
    M3(1:ndim,3) = p(:, ipoint)
    D3           = determinant_4(M3) / D0

    M4(1:ndim,:) = r
    M4(4,:)      = 1
    M4(1:ndim,4) = p(:, ipoint)
    D4           = determinant_4(M4) / D0

    ! Since all Di are divided by D0, they should all be positive now
    ! Since the D? are normalized 1d-10 is enough to take care of fp inaccuracy
    if (any([D1, D2, D3, D4] < -1d-9)) then
      point_in_tetrahedron(ipoint) = .false.

    else
      point_in_tetrahedron(ipoint) = .true.
    end if
  end do

end function point_in_tetrahedron
!------------------------------------------------------------------------------


end module polytopes
!=========================================================================================
