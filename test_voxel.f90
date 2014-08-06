!=========================================================================================
module test_voxel

  use global_parameters
  use voxel
  use inversion_mesh
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_generate_random_points_vox
  ! Tests whether, indeed all random points lie within a voxel
  
  real(kind=dp), dimension(3,1000)  :: points
  real(kind=dp), dimension(3,8)     :: vertices_A,vertices_B

  ! Strangely shaped test voxel A at the north pole
  vertices_A(:,1) = [-5.291225e+02, -6.479882e-14, 6.047898e+03]
  vertices_A(:,2) = [3.239941e-14,  -5.291225e+02, 6.047898e+03]
  vertices_A(:,3) = [2.276260e-29,  -3.717415e-13, 6.071000e+03]
  vertices_A(:,4) = [-3.717415e-13, -4.552521e-29, 6.071000e+03]
  vertices_A(:,5) = [-5.552692e+02, -6.800087e-14, 6.346756e+03]
  vertices_A(:,6) = [3.400043e-14,  -5.552692e+02, 6.346756e+03]
  vertices_A(:,7) = [2.388742e-29,  -3.901112e-13, 6.371000e+03]
  vertices_A(:,8) = [-3.901112e-13, -4.777485e-29, 6.371000e+03]

  ! "Normal" test voxel B somewhere in the the mantle
  vertices_B(:,1) = [-1.202779e+03, -3.304609e+03, 4.191029e+03]
  vertices_B(:,2) = [-9.101866e+02, -3.396863e+03, 4.191029e+03]
  vertices_B(:,3) = [-1.486218e+03, -3.187204e+03, 4.191029e+03]
  vertices_B(:,4) = [-1.758346e+03, -3.045544e+03, 4.191029e+03]
  vertices_B(:,5) = [-1.268733e+03, -3.485815e+03, 4.420842e+03]
  vertices_B(:,6) = [-9.600963e+02, -3.583128e+03, 4.420842e+03]
  vertices_B(:,7) = [-8.567194e+02, -3.197320e+03, 4.727326e+03]
  vertices_B(:,8) = [-1.132124e+03, -3.110486e+03, 4.727326e+03]

  points = generate_random_points_vox(vertices_A, 1000)
  call assert_true(all(point_in_voxel(vertices_A, points)), 'Random points are in triangle')

  points = generate_random_points_vox(vertices_B, 1000)
  call assert_true(all(point_in_voxel(vertices_B, points)), 'Random points are in triangle')

  ! @TODO
  ! * Check whether density is approximately 
  !   equal in 4 subvoxels of these voxel

  ! ratio_region(1) = real(count(point_in_triangle(new_vertices, points))) / real(npoints)

end subroutine test_generate_random_points_vox
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine test_get_volume_vox
  real(kind=dp), dimension(3,8)  :: vertices_A,vertices_B

  ! Strangely shaped test voxel A at the north pole
  vertices_A(:,1) = [-5.291225e+02, -6.479882e-14, 6.047898e+03]
  vertices_A(:,2) = [3.239941e-14,  -5.291225e+02, 6.047898e+03]
  vertices_A(:,3) = [2.276260e-29,  -3.717415e-13, 6.071000e+03]
  vertices_A(:,4) = [-3.717415e-13, -4.552521e-29, 6.071000e+03]
  vertices_A(:,5) = [-5.552692e+02, -6.800087e-14, 6.346756e+03]
  vertices_A(:,6) = [3.400043e-14,  -5.552692e+02, 6.346756e+03]
  vertices_A(:,7) = [2.388742e-29,  -3.901112e-13, 6.371000e+03]
  vertices_A(:,8) = [-3.901112e-13, -4.777485e-29, 6.371000e+03]

  ! "Normal" test voxel B somewhere in the the mantle
  vertices_B(:,1) = [-1.202779e+03, -3.304609e+03, 4.191029e+03]
  vertices_B(:,2) = [-9.101866e+02, -3.396863e+03, 4.191029e+03]
  vertices_B(:,3) = [-1.486218e+03, -3.187204e+03, 4.191029e+03]
  vertices_B(:,4) = [-1.758346e+03, -3.045544e+03, 4.191029e+03]
  vertices_B(:,5) = [-1.268733e+03, -3.485815e+03, 4.420842e+03]
  vertices_B(:,6) = [-9.600963e+02, -3.583128e+03, 4.420842e+03]
  vertices_B(:,7) = [-8.567194e+02, -3.197320e+03, 4.727326e+03]
  vertices_B(:,8) = [-1.132124e+03, -3.110486e+03, 4.727326e+03]


!  call assert_comparable_real(real(get_volume_vox(vertices_A)), 1./6., 1e-8, &
!       ' Volume of voxel A')

!  call assert_comparable_real(real(get_volume_vox(vertices_B)), 1./6., 1e-8, &
!       ' Volume of voxel A')


end subroutine test_get_volume_vox
!-----------------------------------------------------------------------------------------

! @ TODO
! * Subroutine to test spherical to cartesian 
! * converions and vice versa


end module
!=========================================================================================
