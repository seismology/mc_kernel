!=========================================================================================
module test_tetrahedra

  use global_parameters
  use tetrahedra
  use inversion_mesh, only: plane_exp_pro2
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_generate_random_point_poly_3

  integer, parameter                  :: npoints = 100000
  real(kind=dp), dimension(2,npoints) :: points
  real(kind=dp), dimension(2,3)       :: vertices, new_vertices
  real(kind=sp)                       :: ratio_region(3)


  call random_number(vertices)
  vertices = (vertices-0.5) * 2

  points = generate_random_points_poly(3, vertices, npoints)

  call assert_true(all(point_in_triangle(vertices, points)), 'Random points are in triangle')

  ! Divide in half, A, B, BC/2+B
  new_vertices(:,1) =  vertices(:,1)
  new_vertices(:,2) =  vertices(:,2)
  new_vertices(:,3) = (vertices(:,3) - vertices(:,2)) * 0.5 + vertices(:,2)
  
  ratio_region(1) = real(count(point_in_triangle(new_vertices, points))) / real(npoints)
  
  ! Divide in half, B, C, CA/2+C
  new_vertices(:,1) =  vertices(:,2)
  new_vertices(:,2) =  vertices(:,3)
  new_vertices(:,3) = (vertices(:,1) - vertices(:,3)) * 0.5 + vertices(:,3)

  ratio_region(2) = real(count(point_in_triangle(new_vertices, points))) / real(npoints)
  
  ! Divide in half, C, A, AB/2+A
  new_vertices(:,1) =  vertices(:,3)
  new_vertices(:,2) =  vertices(:,1)
  new_vertices(:,3) = (vertices(:,2) - vertices(:,1)) * 0.5 + vertices(:,1)

  ratio_region(3) = real(count(point_in_triangle(new_vertices, points))) / real(npoints)
  
  call assert_comparable_real(ratio_region(1), 0.5, 1e-2, 'Correct density in Region 1')
  call assert_comparable_real(ratio_region(2), 0.5, 1e-2, 'Correct density in Region 2')
  call assert_comparable_real(ratio_region(3), 0.5, 1e-2, 'Correct density in Region 3')

end subroutine test_generate_random_point_poly_3
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine test_generate_random_point_poly_4

  integer, parameter                  :: npoints = 100000
  real(kind=dp), dimension(2,npoints) :: points
  real(kind=dp), dimension(2,4)       :: vertices
  real(kind=sp)                       :: ratio_region(5)

  vertices(:,1) = [ 0,  1]
  vertices(:,2) = [ 1,  0]
  vertices(:,3) = [ 0, -1]
  vertices(:,4) = [-1,  0]

  points = generate_random_points_poly(4, vertices, npoints)

  call assert_true(all(sum(abs(points),1)<1), 'Random points are in polygon')

  ! Check for point density in four regions
  ! Region 1: x>0.5
  ratio_region(1) = real(count(points(1,:).ge.0.5))  / real(npoints)
  ! Region 2: x<-0.5
  ratio_region(2) = real(count(points(1,:).le.-0.5)) / real(npoints)
  ! Region 3: y>0.5
  ratio_region(3) = real(count(points(2,:).ge.0.5))  / real(npoints)
  ! Region 4: y<-0.5
  ratio_region(4) = real(count(points(2,:).le.-0.5)) / real(npoints)
  ! Region 5: inner polygon with |x|+|y|<0.5
  ratio_region(5) = real(count(sum(abs(points),1).le.0.5)) / real(npoints)

  call assert_comparable_real(ratio_region(1), 0.125, 5e-2, 'Correct density in Region 1')
  call assert_comparable_real(ratio_region(2), 0.125, 5e-2, 'Correct density in Region 2')
  call assert_comparable_real(ratio_region(3), 0.125, 5e-2, 'Correct density in Region 3')
  call assert_comparable_real(ratio_region(4), 0.125, 5e-2, 'Correct density in Region 4')
  call assert_comparable_real(ratio_region(5), 0.25,  5e-2, 'Correct density in Region 5')

end subroutine test_generate_random_point_poly_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_generate_random_point_tet

  real(kind=dp), dimension(3,1000)  :: points
  real(kind=dp), dimension(3,4)     :: vertices

  vertices(:,1) = [1, 0, 0]
  vertices(:,2) = [0, 1, 0]
  vertices(:,3) = [0, 0, 1]
  vertices(:,4) = [0, 0, 0]

  points = generate_random_points_tet(vertices, 1000)

  call assert_true(all(sum(points,1)<1), 'Random points are in tetrahedron')

end subroutine test_generate_random_point_tet
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_generate_random_point_triangle_space

  integer, parameter                  :: npoints = 100000
  logical                             :: isintriangle(npoints)
  real(kind=dp), dimension(3,npoints) :: points
  real(kind=dp), dimension(2,npoints) :: points_2d
  real(kind=dp), dimension(2,3)       :: vertices, new_vertices
  real(kind=dp)                       :: p_2d(2,3), v_2d(3,0:2), p_ref(3,3)
  integer                             :: nvec, ipoint
  real(kind=sp)                       :: ratio_region(3)

  isintriangle = .false.
  nvec = 2
  p_ref(:,1) = [0, 0, 0]
  p_ref(:,2) = [1, 0, 1]
  p_ref(:,3) = [1, 1, 1]
  call plane_exp_pro2(p_ref   = p_ref,        &
                      npoints = nvec,         &
                      p_3d    = p_ref(:,2:3), &
                      p_2d    = p_2d(:,2:3),    &
                      vec     = v_2d(:,1:2)    )

  points_2d = generate_random_points_poly(3, p_2d, npoints)
   
  v_2d(:,0) = p_ref(:,1)
  do ipoint = 1, npoints
      points(:,ipoint) =   v_2d(:,0)                          &
                         + v_2d(:,1) * points_2d(1,ipoint)    &
                         + v_2d(:,2) * points_2d(2,ipoint)
  end do

  isintriangle = (points(1,:).eq.points(3,:) .and. &
                  points(3,:).le.1           .and. &
                  points(3,:).ge.points(2,:)      )

  call assert_true(isintriangle, 'Points are in triangle')

  ! x==z<0.5, ratio should be 0.25
  ratio_region(1) = real(count(points(1,:)<0.5)) / real(npoints)
  ! y>0.5, ratio should be 0.25
  ratio_region(2) = real(count(points(2,:)>0.5)) / real(npoints)
  ! y<0.5 and x>0.5, ratio should be 0.5
  ratio_region(3) = real(count(points(2,:)<0.5.and.points(1,:)>0.5)) / real(npoints)

  call assert_comparable_real(ratio_region(1), 0.25, 5e-2, 'Correct density in Region 1')
  call assert_comparable_real(ratio_region(2), 0.25, 5e-2, 'Correct density in Region 2')
  call assert_comparable_real(ratio_region(3), 0.5,  5e-2, 'Correct density in Region 3')

end subroutine test_generate_random_point_triangle_space
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_get_volume_poly_3
  real(kind=dp), dimension(3,3)  :: vertices
  real(kind=sp)                  :: area, area_exp
  
  ! Simple triangle in plane
  vertices(:,1) = [0, 0, 0]
  vertices(:,2) = [1, 0, 0]
  vertices(:,3) = [0, 1, 0]
  area = real(get_volume_poly(3,vertices))
  area_exp = 0.5
  call assert_comparable_real(area, area_exp, 1.e-7, 'Area of planar triangle')


  ! Triangle in space
  vertices(:,1) = [0, 0, 0]
  vertices(:,2) = [1, 1, 0]
  vertices(:,3) = [1, 0, 1]

  area = real(get_volume_poly(3,vertices))
  area_exp = sqrt(2.0)*sqrt(6.0)/4.d0
  call assert_comparable_real(area, area_exp, 1.e-7, 'Area of triangle in space')

end subroutine test_get_volume_poly_3
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_get_volume_poly_4
  real(kind=dp), dimension(3,4)  :: vertices, vertices_rot
  real(kind=sp)                  :: area, area_ref

  ! Random values generated with Matlab script unit_tests/rotate_poly.m
  vertices(1,:) = [0.8147236864, 0.9057919371, 0.1269868163, 0.9133758561]
  vertices(2,:) = [0.6323592462, 0.0975404050, 0.2784982189, 0.5468815192]
  vertices(3,:) = [0, 0, 0, 0]
  area_ref = 0.1531723990
  !N = [0.6996987980, 0.7050929803, 0.1151758710]
  !alpha = 349.4134014338
  vertices_rot(1,:) = [0.8263364349, 0.9008053136, 0.1341153648, 0.9216049462]
  vertices_rot(2,:) = [0.6165487444, 0.0851452750, 0.2744938058, 0.5305435831]
  vertices_rot(3,:) = [0.0262420220, 0.1061754819, -0.0187917249, 0.0500267944]
  
  area = get_volume_poly(4,vertices)
  call assert_comparable_real(area, area_ref, 1e-7, 'Area of planar quadrilateral')
  
  area = get_volume_poly(4,vertices_rot)
  call assert_comparable_real(area, area_ref, 1e-7, 'Area of rotated quadrilateral')
  
end subroutine test_get_volume_poly_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_tetra_volume_3d
  real(kind=dp), dimension(3,4)  :: vertices

  vertices(:,1) = [1, 0, 0]
  vertices(:,2) = [0, 1, 0]
  vertices(:,3) = [0, 0, 1]
  vertices(:,4) = [0, 0, 0]

  call assert_comparable_real(real(get_volume_tet(vertices)), 1./6., 1e-8, &
                              ' Volume of tetrahedron 1')


  vertices(:,1) = [11, -20, 40]
  vertices(:,2) = [10, -21, 40]
  vertices(:,3) = [10, -20, 41]
  vertices(:,4) = [10, -20, 40]

  call assert_comparable_real(real(get_volume_tet (vertices)), 1./6., 1e-8, &
                              ' Volume of tetrahedron 2')

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
