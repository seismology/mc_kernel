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

  real(dp)              :: points(2)
  real(sp)              :: coefficients(2)
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

  real(dp)              :: points1(2), points2(2)
  real(sp)              :: coefficients(2,2)
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
subroutine test_lagrange_interpol_2D_td

  real(dp)              :: points1(2), points2(2)
  real(sp)              :: coefficients(1,2,2)
  real(dp)              :: x1, x2, interpol(1), interpol_ref

  points1 = [-1, 1]
  points2 = [-1, 1]
  coefficients = reshape([1, 1, 1, 1], [1,2,2])
  x1 = 0
  x2 = 0
  interpol = lagrange_interpol_2D_td(points1, points2, coefficients, x1, x2)
  interpol_ref = 1
  call assert_comparable(10 + real(interpol(1)), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 2d 1')

  coefficients = reshape([-1, 1, -1, 1], [1,2,2])
  x1 = 0
  x2 = 0
  interpol = lagrange_interpol_2D_td(points1, points2, coefficients, x1, x2)
  interpol_ref = 0
  call assert_comparable(10 + real(interpol(1)), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 2D 2')

  x1 = 0
  x2 = 0.5
  interpol = lagrange_interpol_2D_td(points1, points2, coefficients, x1, x2)
  interpol_ref = 0.
  call assert_comparable(10 + real(interpol(1)), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 2D 3')

  x1 = 0.5
  x2 = 0.
  interpol = lagrange_interpol_2D_td(points1, points2, coefficients, x1, x2)
  interpol_ref = 0.5
  call assert_comparable(10 + real(interpol(1)), 10 + real(interpol_ref), 1e-7, &
                         'lagrange interpolation 2D 4')

end subroutine test_lagrange_interpol_2D_td
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_gll_points()

  real(dp), allocatable :: eta(:), eta_ref(:)
  integer               :: n

  n = 1
  allocate(eta(n+1), eta_ref(n+1))
  eta = zelegl(n)
  eta_ref = [-1, 1]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=1')

  n = 2
  deallocate(eta, eta_ref)
  allocate(eta(n+1), eta_ref(n+1))
  eta = zelegl(n)
  eta_ref = [-1, 0, 1]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=2')

  ! reference values from:
  ! http://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Lobatto_rules

  ! alternatively in mathematica: 
  ! NSolve[(1 - x*x) D[LegendreP[n, x], x] == 0, x]
  
  n = 3
  deallocate(eta, eta_ref)
  allocate(eta(n+1), eta_ref(n+1))
  eta = zelegl(n)
  eta_ref = [-1d0, -.2d0**.5d0, .2d0**.5d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=3')

  n = 4
  deallocate(eta, eta_ref)
  allocate(eta(n+1), eta_ref(n+1))
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
  allocate(eta(n+1), eta_ref(n+1))
  eta = zemngl2(n)
  eta_ref = [-1, 1]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'glj points n=1')

  n = 2
  deallocate(eta, eta_ref)
  allocate(eta(n+1), eta_ref(n+1))
  eta = zemngl2(n)
  eta_ref = [-1d0, 0.2d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'glj points n=2')

  n = 3
  deallocate(eta, eta_ref)
  allocate(eta(n+1), eta_ref(n+1))
  eta = zemngl2(n)
  eta_ref = [-1d0, -0.26120387496374153d0, 0.54691816067802723d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'glj points n=3')

  n = 4
  deallocate(eta, eta_ref)
  allocate(eta(n+1), eta_ref(n+1))
  eta = zemngl2(n)
  eta_ref = [-1d0, -0.50778762955831502d0, 0.13230082077732325d0, &
             0.70882014211432509d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'glj points n=4')

end subroutine test_glj_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_derivative_tensors()
  
  integer                    :: npol
  real(kind=dp), allocatable :: G0(:), G0_ref(:)
  real(kind=dp), allocatable :: G1(:,:), G1_ref(:)
  real(kind=dp), allocatable :: G2(:,:), G2_ref(:)

  npol = 4
  allocate(G0(0:npol))
  allocate(G1(0:npol,0:npol))
  allocate(G2(0:npol,0:npol))

  allocate(G0_ref(0:npol))
  allocate(G1_ref(1:(npol+1)**2))
  allocate(G2_ref(1:(npol+1)**2))


  G1_ref = &
   [-4.00000000,  6.69583750, -4.63974380,   3.19390655,  -1.25000000, &
    -0.61643892, -1.01582170,  2.49033880,  -1.36116433,   0.50308620, &
     0.16810566, -0.98008018, -0.44157877,   1.80197537,  -0.54842215, &
    -0.10722228,  0.49635028, -1.66964221,  -0.29259955,   1.57311368, &
     0.20000000, -0.87433373,  2.42184663,  -7.49751282,   5.75000000]

  G2_ref = &
   [-5.00000000,  6.75650263, -2.66666675,  1.41016424, -0.50000000, &
    -1.24099028, -0.00000000,  1.74574316, -0.76376259,  0.25900974, &
     0.37500000, -1.33658457,  0.00000000,  1.33658457, -0.37500000, &
    -0.25900974,  0.76376259, -1.74574316, -0.00000000,  1.24099028, &
     0.50000000, -1.41016424,  2.66666675, -6.75650263,  5.00000000]

  G0_ref = [-4.0, 6.69583750, -4.63974380, 3.19390655, -1.25]


  G1 = def_lagrange_derivs_glj(npol, G0)
  G2 = def_lagrange_derivs_gll(npol)

  call assert_comparable_real1d(real(reshape(G1, (/(npol+1)**2/))), &
                                real(G1_ref), &
                                1e-7, 'derivatives tensor, glj')

  call assert_comparable_real1d(real(G0), real(G0_ref), &
                                1e-7, 'derivatives tensor, axial')

  call assert_comparable_real1d(real(reshape(G2, (/(npol+1)**2/))), &
                                real(G2_ref), &
                                1e-7, 'derivatives tensor, gll')

end subroutine test_derivative_tensors
!-----------------------------------------------------------------------------------------


end module
!=========================================================================================
