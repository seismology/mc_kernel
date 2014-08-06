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
  real(kind=dp), dimension(3,8)     :: vertices_A,vertices_B,vertices_C
  real(kind=dp)                     :: ratio_subvoxel(8)

  ! "Weird" test voxel A (1 in vox_15l_5deg.dat), at the north pole
  ! Connectivity 1: 1, 2, 3, 4, 5, 6, 7, 8
  ! lamin=85.0 lamax=90.0 lomin=0.0 lomax=90.0 rmin=6071.0 rmax=6371.0
  vertices_A(:,1) = [-5.291225e+02, -6.479882e-14, 6.047898e+03]
  vertices_A(:,2) = [3.239941e-14,  -5.291225e+02, 6.047898e+03]
  vertices_A(:,3) = [2.276260e-29,  -3.717415e-13, 6.071000e+03]
  vertices_A(:,4) = [-3.717415e-13, -4.552521e-29, 6.071000e+03]
  vertices_A(:,5) = [-5.552692e+02, -6.800087e-14, 6.346756e+03]
  vertices_A(:,6) = [3.400043e-14,  -5.552692e+02, 6.346756e+03]
  vertices_A(:,7) = [2.388742e-29,  -3.901112e-13, 6.371000e+03]
  vertices_A(:,8) = [-3.901112e-13, -4.777485e-29, 6.371000e+03]

  ! Normal test voxel B (4193 in vox_15l_5deg.dat), somewhere in the the mantle
  ! Connectivity 4193: 6798, 6799, 6697, 6696, 4641, 4642, 4540, 4539
  ! lamin=50.0 lamax=55.0 lomin=70.0 lomax=75.0 rmin=5471.0 rmax=5771.0
  vertices_B(:,1) = [-1.202779e+03, -3.304609e+03, 4.191029e+03]
  vertices_B(:,2) = [-9.101866e+02, -3.396863e+03, 4.191029e+03]
  vertices_B(:,3) = [-8.121837e+02, -3.031111e+03, 4.481581e+03]
  vertices_B(:,4) = [-1.073272e+03, -2.948790e+03, 4.481581e+03]
  vertices_B(:,5) = [-1.268733e+03, -3.485815e+03, 4.420842e+03]
  vertices_B(:,6) = [-9.600963e+02, -3.583128e+03, 4.420842e+03]
  vertices_B(:,7) = [-8.567194e+02, -3.197320e+03, 4.727326e+03]
  vertices_B(:,8) = [-1.132124e+03, -3.110486e+03, 4.727326e+03]

  ! Test voxel C (957 in vox_15l_5deg.dat), close to the equator
  ! Connectivity 957: 2153, 2155, 2009, 2007, 2154, 2156, 2010, 2008
  ! lamin=0.0 lamax=5.0 lomin=170.0 lomax=175.0 rmin=6071.0 rmax=6371.0
  vertices_C(:,1) = [5.978768e+03, -1.054218e+03, 0.000000e+00]
  vertices_C(:,2) = [6.047898e+03, -5.291225e+02, 0.000000e+00]
  vertices_C(:,3) = [6.024884e+03, -5.271090e+02, 5.291225e+02]
  vertices_C(:,4) = [5.956017e+03, -1.050206e+03, 5.291225e+02]
  vertices_C(:,5) = [6.274210e+03, -1.106313e+03, 0.000000e+00]
  vertices_C(:,6) = [6.346756e+03, -5.552692e+02, 0.000000e+00]
  vertices_C(:,7) = [6.322605e+03, -5.531563e+02, 5.552692e+02]
  vertices_C(:,8) = [6.250335e+03, -1.102103e+03, 5.552692e+02]

  points = generate_random_points_vox(vertices_A, 1000)
  call assert_true(all(point_in_voxel(vertices_A, points)), 'Random points are in voxel A')

  points = generate_random_points_vox(vertices_B, 1000)
  call assert_true(all(point_in_voxel(vertices_B, points)), 'Random points are in voxel B')

  points = generate_random_points_vox(vertices_C, 1000)
  call assert_true(all(point_in_voxel(vertices_C, points)), 'Random points are in voxel C')

  ! @TODO
  ! * Check whether density is approximately 
  !   equal in 4 subvoxels of these voxel

  ! ratio_region(1) = real(count(point_in_triangle(new_vertices, points))) / real(npoints)

end subroutine test_generate_random_points_vox
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine test_get_volume_vox
  real(kind=dp), dimension(3,8)  :: vertices_A,vertices_B,vertices_C
  real(kind=dp)                  :: sph_rng(6) ! the 8 vertices in spherical coordinates
  real(kind=dp)                  :: volume_box

  ! "Weird" test voxel A (1 in vox_15l_5deg.dat), at the north pole
  ! Connectivity 1: 1, 2, 3, 4, 5, 6, 7, 8
  ! lamin=85.0 lamax=90.0 lomin=0.0 lomax=90.0 rmin=6071.0 rmax=6371.0
  vertices_A(:,1) = [-5.291225e+02, -6.479882e-14, 6.047898e+03]
  vertices_A(:,2) = [3.239941e-14,  -5.291225e+02, 6.047898e+03]
  vertices_A(:,3) = [2.276260e-29,  -3.717415e-13, 6.071000e+03]
  vertices_A(:,4) = [-3.717415e-13, -4.552521e-29, 6.071000e+03]
  vertices_A(:,5) = [-5.552692e+02, -6.800087e-14, 6.346756e+03]
  vertices_A(:,6) = [3.400043e-14,  -5.552692e+02, 6.346756e+03]
  vertices_A(:,7) = [2.388742e-29,  -3.901112e-13, 6.371000e+03]
  vertices_A(:,8) = [-3.901112e-13, -4.777485e-29, 6.371000e+03]

  ! Normal test voxel B (4193 in vox_15l_5deg.dat), somewhere in the the mantle
  ! Connectivity 4193: 6798, 6799, 6697, 6696, 4641, 4642, 4540, 4539
  ! lamin=50.0 lamax=55.0 lomin=70.0 lomax=75.0 rmin=5471.0 rmax=5771.0
  vertices_B(:,1) = [-1.202779e+03, -3.304609e+03, 4.191029e+03]
  vertices_B(:,2) = [-9.101866e+02, -3.396863e+03, 4.191029e+03]
  vertices_B(:,3) = [-8.121837e+02, -3.031111e+03, 4.481581e+03]
  vertices_B(:,4) = [-1.073272e+03, -2.948790e+03, 4.481581e+03]
  vertices_B(:,5) = [-1.268733e+03, -3.485815e+03, 4.420842e+03]
  vertices_B(:,6) = [-9.600963e+02, -3.583128e+03, 4.420842e+03]
  vertices_B(:,7) = [-8.567194e+02, -3.197320e+03, 4.727326e+03]
  vertices_B(:,8) = [-1.132124e+03, -3.110486e+03, 4.727326e+03]


  ! Test voxel C (957 in vox_15l_5deg.dat), close to the equator
  ! Connectivity 957: 2153, 2155, 2009, 2007, 2154, 2156, 2010, 2008
  ! lamin=0.0 lamax=5.0 lomin=170.0 lomax=175.0 rmin=6071.0 rmax=6371.0
  vertices_C(:,1) = [5.978768e+03, -1.054218e+03, 0.000000e+00]
  vertices_C(:,2) = [6.047898e+03, -5.291225e+02, 0.000000e+00]
  vertices_C(:,3) = [6.024884e+03, -5.271090e+02, 5.291225e+02]
  vertices_C(:,4) = [5.956017e+03, -1.050206e+03, 5.291225e+02]
  vertices_C(:,5) = [6.274210e+03, -1.106313e+03, 0.000000e+00]
  vertices_C(:,6) = [6.346756e+03, -5.552692e+02, 0.000000e+00]
  vertices_C(:,7) = [6.322605e+03, -5.531563e+02, 5.552692e+02]
  vertices_C(:,8) = [6.250335e+03, -1.102103e+03, 5.552692e+02]

  call assert_comparable_real(real(get_volume_vox(vertices_A)), 6.9412e7, 1e5, &
       ' Volume of voxel A is correct')

  call assert_comparable_real(real(get_volume_vox(vertices_B)), 4.3940e7, 1e5, &
       ' Volume of voxel B is correct')

  ! For test voxel C, which is close to the euqator, the volume of 
  ! a simple box, should be "close" to the true voxel volume
  call cartesian_to_spherical_range(vertices_C,sph_rng)
  volume_box = dabs(sph_rng(1)-sph_rng(2))*(((2*pi*sph_rng(2))/360)*&
               dabs(sph_rng(5)*rad2deg-sph_rng(6)*rad2deg))**2

  call assert_comparable_real(real(get_volume_vox(vertices_C)/volume_box), 0.9524, 1e-3, &
       ' Ratio of simple box to true volume of voxel C is correct')


end subroutine test_get_volume_vox
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine test_coordinate_transformations_vox
  real(kind=dp), dimension(3,8)  :: vertices_A,vertices_B,vertices_C
  real(kind=dp)                  ::  sph_rng(6) 
  real(kind=dp)                  ::  sph_pnt(3,8) 
  real(kind=dp)                  ::  crt_pnt(3) 

  ! "Weird" test voxel A (1 in vox_15l_5deg.dat), at the north pole
  ! Connectivity 1: 1, 2, 3, 4, 5, 6, 7, 8
  ! lamin=85.0 lamax=90.0 lomin=0.0 lomax=90.0 rmin=6071.0 rmax=6371.0
  vertices_A(:,1) = [-5.291225e+02, -6.479882e-14, 6.047898e+03]
  vertices_A(:,2) = [3.239941e-14,  -5.291225e+02, 6.047898e+03]
  vertices_A(:,3) = [2.276260e-29,  -3.717415e-13, 6.071000e+03]
  vertices_A(:,4) = [-3.717415e-13, -4.552521e-29, 6.071000e+03]
  vertices_A(:,5) = [-5.552692e+02, -6.800087e-14, 6.346756e+03]
  vertices_A(:,6) = [3.400043e-14,  -5.552692e+02, 6.346756e+03]
  vertices_A(:,7) = [2.388742e-29,  -3.901112e-13, 6.371000e+03]
  vertices_A(:,8) = [-3.901112e-13, -4.777485e-29, 6.371000e+03]

  ! Normal test voxel B (4193 in vox_15l_5deg.dat), somewhere in the the mantle
  ! Connectivity 4193: 6798, 6799, 6697, 6696, 4641, 4642, 4540, 4539
  ! lamin=50.0 lamax=55.0 lomin=70.0 lomax=75.0 rmin=5471.0 rmax=5771.0
  vertices_B(:,1) = [-1.202779e+03, -3.304609e+03, 4.191029e+03]
  vertices_B(:,2) = [-9.101866e+02, -3.396863e+03, 4.191029e+03]
  vertices_B(:,3) = [-8.121837e+02, -3.031111e+03, 4.481581e+03]
  vertices_B(:,4) = [-1.073272e+03, -2.948790e+03, 4.481581e+03]
  vertices_B(:,5) = [-1.268733e+03, -3.485815e+03, 4.420842e+03]
  vertices_B(:,6) = [-9.600963e+02, -3.583128e+03, 4.420842e+03]
  vertices_B(:,7) = [-8.567194e+02, -3.197320e+03, 4.727326e+03]
  vertices_B(:,8) = [-1.132124e+03, -3.110486e+03, 4.727326e+03]

  ! Test voxel C (957 in vox_15l_5deg.dat), close to the equator
  ! Connectivity 957: 2153, 2155, 2009, 2007, 2154, 2156, 2010, 2008
  ! lamin=0.0 lamax=5.0 lomin=170.0 lomax=175.0 rmin=6071.0 rmax=6371.0
  vertices_C(:,1) = [5.978768e+03, -1.054218e+03, 0.000000e+00]
  vertices_C(:,2) = [6.047898e+03, -5.291225e+02, 0.000000e+00]
  vertices_C(:,3) = [6.024884e+03, -5.271090e+02, 5.291225e+02]
  vertices_C(:,4) = [5.956017e+03, -1.050206e+03, 5.291225e+02]
  vertices_C(:,5) = [6.274210e+03, -1.106313e+03, 0.000000e+00]
  vertices_C(:,6) = [6.346756e+03, -5.552692e+02, 0.000000e+00]
  vertices_C(:,7) = [6.322605e+03, -5.531563e+02, 5.552692e+02]
  vertices_C(:,8) = [6.250335e+03, -1.102103e+03, 5.552692e+02]

  ! Test cartesian_to_spherical_range
  ! Note that cartesian_to_spherical_range returns colatitude, and 
  ! longitude between -180 and 180 when comparing to lat and lon
  ! as given above

  ! Voxel A
  call cartesian_to_spherical_range(vertices_A,sph_rng)
  call assert_comparable_real(real(sph_rng(1)), 6371., 1e-1, ' Voxel A: rmax is correct')
  call assert_comparable_real(real(sph_rng(2)), 6071., 1e-1, ' Voxel A: rmin is correct')
  call assert_comparable_real(90.-real(sph_rng(3)*rad2deg), 90., 1e-3, ' Voxel A: lamax is correct')
  call assert_comparable_real(90.-real(sph_rng(4)*rad2deg), 85., 1e-3, ' Voxel A: lamin is correct')
  call assert_comparable_real(real(sph_rng(6)*rad2deg), -90., 1e-3, ' Voxel A: lomax is correct')
  call assert_comparable_real(real(sph_rng(5)*rad2deg), -180., 1e-3, ' Voxel A: lomin is correct')

  ! Voxel B
  call cartesian_to_spherical_range(vertices_B,sph_rng)
  call assert_comparable_real(real(sph_rng(1)), 5471., 1e-1, ' Voxel B: rmax is correct')
  call assert_comparable_real(real(sph_rng(2)), 5771., 1e-1, ' Voxel B: rmin is correct')
  call assert_comparable_real(90.-real(sph_rng(3)*rad2deg), 55., 1e-3, ' Voxel B: lamax is correct')
  call assert_comparable_real(90.-real(sph_rng(4)*rad2deg), 50., 1e-3, ' Voxel B: lamin is correct')
  call assert_comparable_real(real(sph_rng(6)*rad2deg), -105., 1e-3, ' Voxel B: lomax is correct')
  call assert_comparable_real(real(sph_rng(5)*rad2deg), -110., 1e-3, ' Voxel B: lomin is correct')

  ! Voxel C
  call cartesian_to_spherical_range(vertices_C,sph_rng)
  call assert_comparable_real(real(sph_rng(1)), 6371., 1e-1, ' Voxel C: rmax is correct')
  call assert_comparable_real(real(sph_rng(2)), 6071., 1e-1, ' Voxel C: rmin is correct')
  call assert_comparable_real(90.-real(sph_rng(3)*rad2deg), 5., 1e-3, ' Voxel C: lamax is correct')
  call assert_comparable_real(90.-real(sph_rng(4)*rad2deg), 0., 1e-3, ' Voxel C: lamin is correct')
  call assert_comparable_real(real(sph_rng(6)*rad2deg), -5., 1e-3, ' Voxel C: lomax is correct')
  call assert_comparable_real(real(sph_rng(5)*rad2deg), -10., 1e-3, ' Voxel C: lomin is correct')


  ! Test cartesian_to_spherical_points
  call cartesian_to_spherical_points(vertices_B, sph_pnt)
  call assert_comparable_real(real(sph_pnt(1,1)), 5471., 1e-3, ' r coordinate is correct')
  call assert_comparable_real(90.-real(sph_pnt(2,1)*rad2deg), 50., 1e-3, ' ph coordinate is correct')
  call assert_comparable_real(real(sph_pnt(3,1)*rad2deg), -110., 1e-3, ' th coordinate is correct')

  ! Test spherical_to_cartesian_point
  call spherical_to_cartesian_point(sph_pnt(:,1),crt_pnt)
  call assert_comparable_real(real(crt_pnt(1)), real(vertices_B(1,1)), 1e-3, ' x coordinate is correct')
  call assert_comparable_real(real(crt_pnt(2)), real(vertices_B(2,1)), 1e-3, ' y coordinate is correct')
  call assert_comparable_real(real(crt_pnt(3)), real(vertices_B(3,1)), 1e-3, ' z coordinate is correct')


end subroutine test_coordinate_transformations_vox
!-----------------------------------------------------------------------------------------


end module
!=========================================================================================
