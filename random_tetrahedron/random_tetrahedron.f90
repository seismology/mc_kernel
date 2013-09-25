!=========================================================================================
module tetrahedra

  integer, parameter :: dp = 8, sp = 4
contains

!-----------------------------------------------------------------------------------------
function generate_random_point(v, n)
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

  implicit none

  integer, intent(in)        ::  n
  real(kind=dp), intent(in)  ::  v(3,4)
  real(kind=dp)              ::  generate_random_point(3,n)

  real(kind=dp)              ::  c(4)
  real(kind=sp)              ::  coordinates(3,n)
  integer                    ::  j
  integer                    ::  seed
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
  
     generate_random_point(1:3,j) = matmul( v(1:3,1:4), dble(c(1:4)) )

  end do

  return

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function rmat4_det ( a )
! RMAT4_DET computes the determinant of a 4 by 4 matrix.
  implicit none
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
function tetra_volume_3d ( v )
! TETRA_VOLUME_3D computes the volume of a tetrahedron in 3D.
  implicit none

  real(kind=dp), intent(in)  ::  v(3,4)
  real(kind=dp)              ::  a(4,4)
  real(kind=dp)              ::  tetra_volume_3d

  a(1:4,1:3) = transpose(v)
  a(1:4,4) = [1, 1, 1, 1]

  tetra_volume_3d = abs ( rmat4_det ( a ) ) / 6.0d+00

end function tetra_volume_3d
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
