!=========================================================================================
module test_sem_derivatives

  use global_parameters
  use finite_elem_mapping
  use sem_derivatives
  use spectral_basis
  use ftnunit
  implicit none
  public

contains

!-----------------------------------------------------------------------------------------
subroutine test_gradient()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol
  real(kind=dp), allocatable   :: G1(:,:), G1T(:,:)
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: u(:,:), grad_u(:,:,:), grad_u_ref(:,:,:)

  real(kind=dp), allocatable   :: glj_points(:)
  real(kind=dp), allocatable   :: gll_points(:)

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1

  npol = 1

  G1 = def_lagrange_derivs_glj(npol)
  G2 = def_lagrange_derivs_gll(npol)

  G1T = transpose(G1)
  G2T = transpose(G2)

  glj_points = zemngl2(npol)
  gll_points = zelegl(npol)

  allocate(u(0:npol,0:npol))
  allocate(grad_u_ref(0:npol,0:npol,2))
  
  u = 1
  grad_u_ref = 0
  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'gradient of a constant = 0')

  u(0,:) = 0
  u(1,:) = 1
  grad_u_ref(:,:,1) = 1
  grad_u_ref(:,:,2) = 0
  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'linear in s')

  u(:,0) = 0
  u(:,1) = 1
  grad_u_ref(:,:,1) = 0
  grad_u_ref(:,:,2) = 1
  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'linear in z')

end subroutine 
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
