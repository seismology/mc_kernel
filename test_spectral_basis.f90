!=========================================================================================
module test_spectral_basis

   use global_parameters
   use spectral_basis
   use ftnunit
   implicit none
   public

contains

!-----------------------------------------------------------------------------------------
subroutine test_lagrange_interpol_1D

  real(dp), allocatable :: points(:), coefficients(:)
  real(dp)              :: x, interpol, interpol_ref

  ! automatic allocation
  points = [-1, 1]
  coefficients = [-1, 1]
  x = 0
  interpol = lagrange_interpol_1D(points, coefficients, x)
  interpol_ref = 0
  call assert_comparable(10 + real(interpol), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 1')

  x = -1
  interpol = lagrange_interpol_1D(points, coefficients, x)
  interpol_ref = -1
  call assert_comparable(10 + real(interpol), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 2')

  x = 1
  interpol = lagrange_interpol_1D(points, coefficients, x)
  interpol_ref = 1
  call assert_comparable(10 + real(interpol), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 3')

  x = 2
  interpol = lagrange_interpol_1D(points, coefficients, x)
  interpol_ref = 2
  call assert_comparable(10 + real(interpol), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 4')

end subroutine test_lagrange_interpol_1D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_lagrange_interpol_2D

  real(dp), allocatable :: points1(:), points2(:), coefficients(:,:)
  real(dp)              :: x1, x2, interpol, interpol_ref

  ! automatic allocation
  points1 = [-1, 1]
  points2 = [-1, 1]
  coefficients = reshape([1, 1, 1, 1], [2,2])
  x1 = 0
  x2 = 0
  interpol = lagrange_interpol_2D(points1, points2, coefficients, x1, x2)
  interpol_ref = 1
  call assert_comparable(10 + real(interpol), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 2d 1')

  coefficients = reshape([-1, 1, -1, 1], [2,2])
  x1 = 0
  x2 = 0
  interpol = lagrange_interpol_2D(points1, points2, coefficients, x1, x2)
  interpol_ref = 0
  call assert_comparable(10 + real(interpol), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 2D 2')

  x1 = 0
  x2 = 0.5
  interpol = lagrange_interpol_2D(points1, points2, coefficients, x1, x2)
  interpol_ref = 0.
  call assert_comparable(10 + real(interpol), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 2D 3')

  x1 = 0.5
  x2 = 0.
  interpol = lagrange_interpol_2D(points1, points2, coefficients, x1, x2)
  interpol_ref = 0.5
  call assert_comparable(10 + real(interpol), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 2D 4')

end subroutine test_lagrange_interpol_2D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_gll_points()

  real(dp), allocatable :: eta(:), eta_ref(:)
  integer               :: n

  n = 1
  ! automatic allocation
  eta = zelegl(n)
  eta_ref = [-1, 1]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=1')

  n = 2
  eta = zelegl(n)
  eta_ref = [-1, 0, 1]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=2')

  ! reference values from:
  ! http://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Lobatto_rules

  ! alternatively in mathematica: 
  ! NSolve[(1 - x*x) D[LegendreP[n, x], x] == 0, x]
  
  n = 3
  eta = zelegl(n)
  eta_ref = [-1d0, -.2d0**.5d0, .2d0**.5d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=3')

  n = 4
  eta = zelegl(n)
  eta_ref = [-1d0, -(3d0/7d0)**.5d0, 0d0, (3d0/7d0)**.5d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=4')

end subroutine test_gll_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_glj_points()

  real(dp), allocatable :: eta(:), eta_ref(:)
  integer               :: n

  ! reference solution from mathematica: 
  ! Sort[NSolve[D[(LegendreP[n, x] + LegendreP[n + 1, x] )/(1 + x), x] == 0, x]]

  n = 1
  ! automatic allocation
  eta = zemngl2(n)
  eta_ref = [-1, 1]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'glj points n=1')

  n = 2
  eta = zemngl2(n)
  eta_ref = [-1d0, 0.2d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'glj points n=2')

  n = 3
  eta = zemngl2(n)
  eta_ref = [-1d0, -0.26120387496374153d0, 0.54691816067802723d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'glj points n=3')

  n = 4
  eta = zemngl2(n)
  eta_ref = [-1d0, -0.50778762955831502d0, 0.13230082077732325d0, &
             0.70882014211432509d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'glj points n=4')

end subroutine test_glj_points
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
