module random_points

contains
!*****************************************************************************80
!
!! UNIFORM_IN_TETRAHEDRON returns uniform points in a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2009
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
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V(3,4), the vertices of the tetrahedron.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(3,N), the points.
!
function generate_random_point(v, n)
    implicit none

    integer, intent(in)        ::  n
    real(kind=sp), intent(in)  ::  v(3,4)
    real(kind=sp), intent(out) ::  generate_random_point(3,n)

    real(kind=sp)              ::  c(4)
    real(kind=sp)              ::  coordinates(3,n)
    integer                    ::  j
    integer                    ::  seed
    real(kind=8 )          ::  t

    call random_number(coordinates)
    do j = 1, n

        c(1:3) = dble(coordinates(:, j))
      !call random_number(c(1))
      !call random_number(c(2))
      !call random_number(c(3))
       !call r8vec_uniform_01 ( 3, seed, c )

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
!  
!    C(1:4) are the barycentric coordinates of the point.
!  
        x(1:3,j) = real(matmul ( dble(v(1:3,1:4)), c(1:4) ))

    end do

    return
end subroutine

end module
