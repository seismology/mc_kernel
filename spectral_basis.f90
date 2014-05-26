!=========================================================================================
module spectral_basis
    use global_parameters, only: sp, dp, pi

    implicit none
    private

    public :: zelegl

contains

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

end module
!=========================================================================================
