!=========================================================================================
module tetrahedra

  use global_parameters
  implicit none
  private

  public :: generate_random_points_tet, generate_random_points_poly
  public :: get_volume_tet, get_volume_poly
  public :: rmat4_det, point_in_triangle
contains

!-----------------------------------------------------------------------------------------
function generate_random_points_tet(v, n)
!  returns uniform points in a tetrahedron.
!
!    This code is distributed under the GNU LGPL license.
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Claudio Rocchini, Paolo Cignoni,
!    Generating Random Points in a Tetrahedron,
!    Journal of Graphics Tools,
!    Volume 5, Number 5, 2000, pages 9-12.

  integer, intent(in)        ::  n
  real(kind=dp), intent(in)  ::  v(3,4)
  real(kind=dp)              ::  generate_random_points_tet(3,n)

  real(kind=dp)              ::  c(4)
  real(kind=dp)              ::  coordinates(3,n)
  integer                    ::  j
  real(kind=dp )             ::  t

  call random_number(coordinates)
  do j = 1, n

     c(1:3) = dble(coordinates(:, j))

     if ( 1.0D+00 < c(1) + c(2) ) then
        c(1) = 1.0D+00 - c(1)
        c(2) = 1.0D+00 - c(2)
     end if

     if ( 1.0D+00 < c(2) + c(3) ) then
        t = c(3)
        c(3) = 1.0D+00 - c(1) - c(2)
        c(2) = 1.0D+00 - t
     else if ( 1.0D+00 < c(1) + c(2) + c(3) ) then
        t = c(3)
        c(3) = c(1) + c(2) + c(3) - 1.0D+00
        c(1) = 1.0D+00 - c(2) - t
     end if

     c(4) = 1.0D+00 - c(1) - c(2) - c(3)
  
     ! c(1:4) are the barycentric coordinates of the point.
  
     generate_random_points_tet(1:3,j) = matmul( dble(v(1:3,1:4)), dble(c(1:4)) )

  end do

end function generate_random_points_tet
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function generate_random_points_poly( nv, v, n ) result(x)
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
  implicit none

  integer, parameter           :: dim_num = 2
  integer, intent(in)          :: n  ! Number of points
  integer, intent(in)          :: nv ! Degree of polynomial
  real(kind=dp), intent(in)    :: v(dim_num, nv)
  real(kind=dp)                :: x(dim_num, n)

  real(kind=dp)                :: area(nv)
  real(kind=dp)                :: area_norm(nv)
  real(kind=dp)                :: area_int(nv)
  real(kind=dp)                :: area_percent
  real(kind=dp)                :: centroid(dim_num)
  integer                      :: i, ip1, j, k
  real(kind=dp)                :: r(2)
  real(kind=dp)                :: t(dim_num,3)


  if (nv == 3) then
    do j = 1, n
      call random_number(r)

      if ( 1.0D+00 < sum (r(:)) ) then
        r(:) = 1.0D+00 - r(:)
      end if

      x(:,j) = ( 1.0D+00 - r(1) - r(2) ) * v(:,1) &
                         + r(1)          * v(:,2) &
                                + r(2)   * v(:,3) 
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
        call random_number(r)

        if ( 1.0D+00 < sum (r(:)) ) then
          r(:) = 1.0D+00 - r(:)
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
        x(:,j) = ( 1.0D+00 - r(1) - r(2) ) * dble(v(:,i)) &
                           + r(1)          * dble(v(:,ip1)) &
                                  + r(2)   * centroid(:)

      end do

  end if ! nv==3

end function

!-----------------------------------------------------------------------------------------
function rmat4_det ( a )
! RMAT4_DET computes the determinant of a 4 by 4 matrix.
  real(kind=dp), intent(in)  :: a(4,4)
  real(kind=dp)              :: rmat4_det

  rmat4_det = &
      a(1,1) * ( &
        a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
      - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
      + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
    - a(1,2) * ( &
        a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
      - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
      + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
    + a(1,3) * ( &
        a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
      - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
      + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
    - a(1,4) * ( &
        a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
      - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
      + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

end function rmat4_det
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
function get_volume_tet(v)
! TETRA_VOLUME_3D computes the volume of a tetrahedron in 3D.
  real(kind=dp), intent(in)  ::  v(3,4) !< Vertices
  real(kind=dp)              ::  a(4,4)
  real(kind=dp)              ::  get_volume_tet

  a(1:4,1:3) = transpose(v)
  a(1:4,4) = [1, 1, 1, 1]

  get_volume_tet = abs ( rmat4_det ( a ) ) / 6.0d+00

end function get_volume_tet
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_volume_poly(n, v) result(area)
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



function cross(a,b)
! Compute the cross product between vectors a and b
  real(kind=dp)             :: cross(3)
  real(kind=dp), intent(in) :: a(3), b(3)
           
  cross(1) = a(2)*b(3) - a(3)*b(2)
  cross(2) = a(3)*b(1) - a(1)*b(3)
  cross(3) = a(1)*b(2) - b(1)*a(2)

end function cross


function point_in_triangle(r, x)
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

!   function pointInTriangle(x1, y1, x2, y2, x3, y3, x, y:Number):Boolean
!      {
!        var a:Number = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / denominator;
!         var b:Number = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / denominator;
!          var c:Number = 1 - a - b;
!           
!           return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1;
!
end function point_in_triangle
end module
!=========================================================================================
