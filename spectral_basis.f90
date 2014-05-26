!=========================================================================================
module spectral_basis
    use global_parameters, only: sp, dp, pi

    implicit none
    private

    public :: lagrange_interpol_1D
    public :: lagrange_interpol_2D

    public :: zelegl
    public :: zemngl2

contains

!-----------------------------------------------------------------------------------------
!> computes the Lagrangian interpolation polynomial of a function defined by its values at 
!  a set of collocation points
pure function lagrange_interpol_1D(points, coefficients, x)

  real(dp), intent(in)  :: points(0:)
  real(dp), intent(in)  :: coefficients(0:size(points)-1)
  real(dp), intent(in)  :: x
  real(dp)              :: lagrange_interpol_1D
  real(dp)              :: l_j

  integer               :: j, m, n

  n = size(points) - 1
  lagrange_interpol_1D = 0

  !compare: http://en.wikipedia.org/wiki/Lagrange_polynomial#Definition

  do j=0, n
     l_j = 1
     do m=0, n
        if (m == j) cycle
        l_j = l_j * (x - points(m)) / (points(j) - points(m))
     enddo
     lagrange_interpol_1D = lagrange_interpol_1D + coefficients(j) * l_j
  enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> computes the Lagrangian interpolation polynomial of a function defined by its values at 
!  a set of collocation points in 2D, where the points are a tensorproduct of two sets of
!  points in 1D
pure function lagrange_interpol_2D(points1, points2, coefficients, x1, x2)

  real(dp), intent(in)  :: points1(0:), points2(0:)
  real(dp), intent(in)  :: coefficients(0:size(points1)-1, 0:size(points2)-1)
  real(dp), intent(in)  :: x1, x2
  real(dp)              :: lagrange_interpol_2D
  real(dp)              :: l_i(0:size(points1)-1), l_j(0:size(points2)-1)

  integer               :: i, j, m, n1, n2

  n1 = size(points1) - 1
  n2 = size(points2) - 1

  do i=0, n1
     l_i(i) = 1
     do m=0, n1
        if (m == i) cycle
        l_i(i) = l_i(i) * (x1 - points1(m)) / (points1(i) - points1(m))
     enddo
  enddo

  do j=0, n2
     l_j(j) = 1
     do m=0, n2
        if (m == j) cycle
        l_j(j) = l_j(j) * (x2 - points2(m)) / (points2(j) - points2(m))
     enddo
  enddo

  lagrange_interpol_2D = 0

  do i=0, n1
     do j=0, n2
        lagrange_interpol_2D = lagrange_interpol_2D  + coefficients(i,j) * l_i(i) * l_j(j)
     enddo
  enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> computes the nodes relative to the legendre gauss-lobatto formula
pure function zelegl(n)

  integer, intent(in)   :: n            ! Order of the formula
  real(dp)              :: zelegl(0:n)  ! Vector of the nodes
  real(kind=dp)         :: sn, x, c, etx, dy, d2y, y
  integer               :: i, n2, it

  if (n == 0) then
     ! does this make sense at all?
     zelegl(0) = 0
     return
  endif

  n2 = (n - 1) / 2

  sn = 2 * n - 4 * n2 - 3
  zelegl(0) = -1
  zelegl(n) =  1
  
  if (n  ==  1) return

  zelegl(n2+1) = 0
  x = 0
  call valepo(n, x, y, dy, d2y)
  
  if(n == 2) return

  c  = pi / n

  do i=1, n2
     etx = dcos(c * i)
     do it=1, 8
        call valepo(n, etx, y, dy, d2y)
        etx = etx - dy / d2y
     end do
     zelegl(i) = -etx
     zelegl(n-i) = etx
  end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> computes the value of the legendre polynomial of degree n
!! and its first and second derivatives at a given point
pure subroutine valepo(n, x, y, dy, d2y)

  integer, intent(in)   ::  n   !< degree of the polynomial
  real(dp), intent(in)  ::  x   !< point in which the computation is performed
  real(dp), intent(out) ::  y   !< value of the polynomial in x
  real(dp), intent(out) ::  dy  !< value of the first derivative in x
  real(dp), intent(out) ::  d2y !< value of the second derivative in x
  real(kind=dp)         ::  c1, c2, c4, ym, yp, dym, dyp, d2ym, d2yp
  integer               ::  i

  y   = 1.d0
  dy  = 0.d0
  d2y = 0.d0
  if(n == 0) return

  y   = x
  dy  = 1.d0
  d2y = 0.d0
  if(n == 1) return

  yp   = 1.d0
  dyp  = 0.d0
  d2yp = 0.d0
  do i = 2, n
     c1   = dfloat(i)
     c2   = 2.d0*c1-1.d0
     c4   = c1-1.d0
     ym   = y
     y    = (c2*x*y-c4*yp)/c1
     yp   = ym
     dym  = dy
     dy   = (c2*x*dy-c4*dyp+c2*yp)/c1
     dyp  = dym
     d2ym = d2y
     d2y  = (c2*x*d2y-c4*d2yp+2.d0*c2*dyp)/c1
     d2yp = d2ym
  enddo

end subroutine valepo
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!>   Computes the nodes relative to the modified legendre gauss-lobatto
!!   FORMULA along the s-axis
!!   Relies on computing the eigenvalues of tridiagonal matrix. 
!!   The nodes correspond to the second quadrature formula proposed
!!   by Azaiez et al.  
pure function zemngl2(n)
  
  integer, intent(in)                :: n            !< Order of the formula
  real(dp)                           :: zemngl2(0:n) !< vector of the nodes, et(i), i=0,n.
  real(dp), dimension(n-1)           :: d, e
  integer                            :: i, n2
  real(kind=dp)                      :: x

  if (n  ==  0) return

     n2 = (n-1)/2
     zemngl2(0) = -1.d0
     zemngl2(n) = 1.d0
  if (n  ==  1) return

     zemngl2(n2+1) = 2d-1
     x = 2d-1
  if(n  ==  2) return

  ! Form the matrix diagonals and subdiagonals according to
  ! formulae page 109 of Azaiez, Bernardi, Dauge and Maday.

  do i = 1, n-1
     d(i) = 3d0 / (4d0 * (i + 0.5d0) * (i + 3d0 * 2d0))
  end do

  do i = 1, n-2
     e(i+1) = dsqrt(i**2 + 3d0) / (2 * (i + 3d0 / 2))
  end do

  ! Compute eigenvalues
  call tqli(d, e, n-1)

  ! Sort them in increasing order
  call bubblesort(d, e, n-1)

  zemngl2(1:n-1) = e(1:n-1)

end function zemngl2
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> This routines returns the eigenvalues of the tridiagonal matrix 
!! which diagonal and subdiagonal coefficients are contained in d(1:n) and
!! e(2:n) respectively. e(1) is free. The eigenvalues are returned in array d
pure subroutine tqli(d,e,n)

  integer, intent(in)             :: n
  real(kind=dp), intent(inout)    :: d(n)
  real(kind=dp), intent(inout)    :: e(n)
  integer                         :: i, iter, l, m
  real(kind=dp)                   :: b, c, dd, f, g, p, r, s

  do i = 2, n
    e(i-1) = e(i)
  end do

  e(n) = 0
  do l=1, n
     iter = 0
     iterate: do
     do m = l, n-1
       dd = abs(d(m)) + abs(d(m+1))
       if (abs(e(m)) + dd .eq. dd) exit
     end do

     if( m == l ) exit iterate
     !if( iter == 30 ) stop 'too many iterations in tqli'
     iter = iter + 1
     g = (d(l+1) - d(l)) / (2. * e(l))
     r = pythag(g, 1d0)
     g = d(m) - d(l) + e(l) / (g + sign(r,g))
     s = 1
     c = 1
     p = 0
     do i = m-1,l,-1
        f      = s * e(i)
        b      = c * e(i)
        r      = pythag(f, g)
        e(i+1) = r
        if(r == 0)then
           d(i+1) = d(i+1) - p
           e(m)   = 0
           cycle iterate
        endif
        s      = f / r
        c      = g / r
        g      = d(i+1) - p
        r      = (d(i) - g) * s + 2. * c * b
        p      = s * r
        d(i+1) = g + p
        g      = c * r - b
     end do
     d(l) = d(l) - p
     e(l) = g
     e(m) = 0
     end do iterate
  end do

end subroutine tqli
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> This routine reorders array vin(n) in increasing order and outputs array vout(n).
pure subroutine bubblesort(vin, vout, n)

  integer, intent(in)            :: n
  real(kind=dp), intent(in)      :: vin(n)
  real(kind=dp), intent(out)     :: vout(n)
  integer                        :: rankmax
  integer                        :: rank(n)
  integer                        :: i, j

  rankmax = 1

  do i = 1, n

     rank(i) = 1

     do j = 1, n
        if((vin(i) > vin(j)) .and. (i /= j)) rank(i) = rank(i) + 1
     end do

     rankmax = max(rank(i), rankmax)
     vout(rank(i)) = vin(i)

  end do

end subroutine bubblesort
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function pythag(a,b)

  real(kind=dp), intent(in) :: a, b
  real(kind=dp)             :: absa, absb

  absa = dabs(a)
  absb = dabs(b)

  if(absa > absb) then
     pythag = absa * sqrt(1. + (absb / absa)**2)
  elseif(absb == 0)then
     pythag = 0
  else
     pythag = absb * sqrt(1. + (absa / absb)**2)
  endif

end function pythag
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
