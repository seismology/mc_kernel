!=========================================================================================
module test_spectral_basis

   use global_parameters
   use spectral_basis
   use ftnunit
   implicit none
   public

contains

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

  ! no independent reference solution, as this seems not to be standard GLJ(0,1)
  ! integration.
 
  n = 1
  ! automatic allocation
  eta = zemngl2(n)
  eta_ref = [-1, 1]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=1')

  n = 2
  eta = zemngl2(n)
  eta_ref = [-1d0, 0.2d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=1')

  n = 3
  eta = zemngl2(n)
  eta_ref = [-1d0, -0.34589528640782047d0, 0.45482385783639190d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=1')

  n = 4
  eta = zemngl2(n)
  eta_ref = [-1d0, -0.50762608041791035d0, 4.6255004736169870d-2, 0.59410917091983573d0, 1d0]
  call assert_comparable_real1d(5 + real(eta), 5 + real(eta_ref), &
                                1e-7, 'gll points n=1')

end subroutine test_glj_points
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
