module test_tetrahedra
  use tetrahedra
  use ftnunit
  implicit none
  public
contains

!******************************************************************************
  
subroutine test_generate_random_point

  real(kind=dp), dimension(3,1000)  :: points
  real(kind=dp), dimension(3,4)     :: vertices

  vertices(:,1) = [1, 0, 0]
  vertices(:,2) = [0, 1, 0]
  vertices(:,3) = [0, 0, 1]
  vertices(:,4) = [0, 0, 0]

  points = generate_random_point(vertices, 1000)

  call assert_true(all(sum(points,1)<1), 'Random points are in tetrahedron')

end subroutine test_generate_random_point

!******************************************************************************

subroutine test_rmat4_det
  real(kind=dp), dimension(4,4)    :: matrix
  real(kind=sp)                    :: determinant

  matrix(1,:) = [  -1.564774d-01,   8.314711d-01,   5.844147d-01,    9.189849d-01]
  matrix(2,:) = [   3.114814d-01,  -9.285766d-01,   6.982586d-01,    8.679865d-01]
  matrix(3,:) = [   3.574703d-01,   5.154803d-01,   4.862649d-01,   -2.155460d-01] 
  matrix(4,:) = [   3.109558d-01,  -6.576266d-01,   4.120922d-01,   -9.363343d-01]

  determinant = -6.062465e-01

  call assert_comparable_real(real(rmat4_det(matrix)), determinant, 1e-8, &
                              'Determinant of random matrix')

end subroutine test_rmat4_det

!******************************************************************************

subroutine test_tetra_volume_3d
  real(kind=dp), dimension(3,4)  :: vertices


  vertices(:,1) = [1, 0, 0]
  vertices(:,2) = [0, 1, 0]
  vertices(:,3) = [0, 0, 1]
  vertices(:,4) = [0, 0, 0]

  call assert_comparable_real(real(tetra_volume_3d(vertices)), 1./6., 1e-8, &
                              ' Volume of tetrahedron 1')


  vertices(:,1) = [11, -20, 40]
  vertices(:,2) = [10, -21, 40]
  vertices(:,3) = [10, -20, 41]
  vertices(:,4) = [10, -20, 40]

  call assert_comparable_real(real(tetra_volume_3d(vertices)), 1./6., 1e-8, &
                              ' Volume of tetrahedron 2')

end subroutine








end module
