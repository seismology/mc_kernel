!=========================================================================================
module test_inversion_mesh

  use global_parameters
  use inversion_mesh
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_mesh_read
  type(inversion_mesh_type)    :: inv_mesh
  integer                      :: npoints, nelems

  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST')

  npoints = inv_mesh%get_nvertices()
  nelems = inv_mesh%get_nelements()

  call assert_equal(npoints, 5, 'number of vertices in vertices.TEST')
  call assert_equal(nelems, 2, 'number of elements in facets.TEST')

  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_dump
  type(inversion_mesh_type)    :: inv_mesh
  integer                      :: myunit, ierr

  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST')

  call inv_mesh%dump_mesh_xdmf('unit_tests/testmesh')

  call assert_file_exists('unit_tests/testmesh.xdmf', 'test xdmf dump')
  call assert_file_exists('unit_tests/testmesh_points.dat', 'test xdmf dump')
  call assert_file_exists('unit_tests/testmesh_grid.dat', 'test xdmf dump')

  ! tidy up
  !open(newunit=myunit, file='unit_tests/testmesh.xdmf', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testmesh_points.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testmesh_grid.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_dump2
  type(inversion_mesh_type)    :: inv_mesh
  integer                      :: myunit, ierr

  call inv_mesh%read_abaqus_mesh('unit_tests/tetrahedron.inp')

  call inv_mesh%dump_mesh_xdmf('unit_tests/testmesh_abaqus')

  call assert_file_exists('unit_tests/testmesh_abaqus.xdmf', 'test xdmf dump')
  call assert_file_exists('unit_tests/testmesh_abaqus_points.dat', 'test xdmf dump')
  call assert_file_exists('unit_tests/testmesh_abaqus_grid.dat', 'test xdmf dump')

  ! tidy up
  open(newunit=myunit, file='unit_tests/testmesh_abaqus.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testmesh_abaqus_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testmesh_abaqus_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable         :: datat(:,:)
  integer                           :: npoints, myunit, ierr

  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST')
  call inv_mesh%init_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat(3,npoints))

  datat(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_data_snap(datat(1,:), 1, 'x')
  call inv_mesh%set_data_snap(datat(2,:), 2, 'x')
  call inv_mesh%set_data_snap(datat(3,:), 3, 'x')

  call inv_mesh%dump_mesh_data_xdmf('unit_tests/testdata')

  call assert_file_exists('unit_tests/testdata.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testdata_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testdata_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testdata_data.dat', 'test xdmf data dump')

  ! tidy up
  open(newunit=myunit, file='unit_tests/testdata.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testdata_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testdata_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  open(newunit=myunit, file='unit_tests/testdata_data.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump2
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable         :: datat(:,:)
  integer                           :: npoints, myunit, ierr

  call inv_mesh%read_abaqus_mesh('unit_tests/circle.inp')
  call inv_mesh%init_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat(3,npoints))

  datat(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_data_snap(datat(1,:), 1, 'x')
  call inv_mesh%set_data_snap(datat(2,:), 2, 'x')
  call inv_mesh%set_data_snap(datat(3,:), 3, 'x')

  call inv_mesh%dump_mesh_data_xdmf('unit_tests/testcircle')

  call assert_file_exists('unit_tests/testcircle.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_data.dat', 'test xdmf data dump')

  ! tidy up
  open(newunit=myunit, file='unit_tests/testcircle.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testcircle_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testcircle_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  open(newunit=myunit, file='unit_tests/testcircle_data.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump3
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable         :: datat(:,:)
  integer                           :: npoints, myunit, ierr

  call inv_mesh%read_abaqus_mesh('unit_tests/circle_quad2.inp')
  call inv_mesh%init_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat(3,npoints))

  datat(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_data_snap(datat(1,:), 1, 'x')
  call inv_mesh%set_data_snap(datat(2,:), 2, 'x')
  call inv_mesh%set_data_snap(datat(3,:), 3, 'x')

  call inv_mesh%dump_mesh_data_xdmf('unit_tests/testcircle_quad')

  call assert_file_exists('unit_tests/testcircle_quad.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_quad_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_quad_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_quad_data.dat', 'test xdmf data dump')

  ! tidy up
  open(newunit=myunit, file='unit_tests/testcircle_quad.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testcircle_quad_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testcircle_quad_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  open(newunit=myunit, file='unit_tests/testcircle_quad_data.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump4
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable         :: datat(:,:)
  integer                           :: npoints, myunit, ierr

  call inv_mesh%read_abaqus_mesh('unit_tests/sphere.inp')
  call inv_mesh%init_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat(3,npoints))

  datat(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_data_snap(datat(1,:), 1, 'x')
  call inv_mesh%set_data_snap(datat(2,:), 2, 'x')
  call inv_mesh%set_data_snap(datat(3,:), 3, 'x')

  call inv_mesh%dump_mesh_data_xdmf('unit_tests/testsphere')

  call assert_file_exists('unit_tests/testsphere.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testsphere_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testsphere_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testsphere_data.dat', 'test xdmf data dump')

  ! tidy up
  open(newunit=myunit, file='unit_tests/testsphere.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testsphere_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testsphere_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  open(newunit=myunit, file='unit_tests/testsphere_data.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump5
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable        :: datat(:,:)
  integer                           :: npoints, myunit, ierr

  call inv_mesh%read_abaqus_mesh('unit_tests/tetrahedron.inp')
  call inv_mesh%init_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat(3,npoints))

  datat(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_data_snap(datat(1,:), 1, 'x')
  call inv_mesh%set_data_snap(datat(2,:), 2, 'x')
  call inv_mesh%set_data_snap(datat(3,:), 3, 'x')

  call inv_mesh%dump_mesh_data_xdmf('unit_tests/testtetrahedrons')

  call assert_file_exists('unit_tests/testtetrahedrons.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testtetrahedrons_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testtetrahedrons_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testtetrahedrons_data.dat', 'test xdmf data dump')

  ! tidy up
  open(newunit=myunit, file='unit_tests/testtetrahedrons.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testtetrahedrons_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testtetrahedrons_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  open(newunit=myunit, file='unit_tests/testtetrahedrons_data.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_valence
  type(inversion_mesh_type)    :: inv_mesh

  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST')

  call assert_equal(inv_mesh%get_valence(1), 2, 'valence of first vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(2), 2, 'valence of first vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(3), 2, 'valence of first vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(4), 1, 'valence of first vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(5), 1, 'valence of first vertex in facets.TEST')

  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_get_connected_elements
  type(inversion_mesh_type)    :: inv_mesh
  integer, allocatable         :: cn_elems(:)

  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST')

  ! uncomment to NOT use automatic allocation
  !allocate(cn_elems(inv_mesh%get_valence(1)))
  cn_elems = inv_mesh%get_connected_elements(1)

  ! somewhat redundant two tests, for now to test ftnunit :)
  call assert_true(.not. any(cn_elems == -1), 'get connected elements')
  call assert_true(cn_elems /= -1, 'get connected elements')

  call assert_equal(inv_mesh%get_connected_elements(1), (/1 ,2/), 'get connected elements')
  call assert_equal(inv_mesh%get_connected_elements(4), (/1/), 'get connected elements')
  call assert_equal(inv_mesh%get_connected_elements(5), (/2/), 'get connected elements')

  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_initialize_mesh
  type(inversion_mesh_type)  :: inv_mesh
  real(kind=dp)              :: vertices(3,3)
  integer                    :: connectivity(1,3)

  vertices(:,1) = [0, 0, 0]
  vertices(:,2) = [1, 0, 0]
  vertices(:,3) = [0, 0, 1]

  connectivity(1,:) = [1, 2, 3]

  call inv_mesh%initialize_mesh('tri', vertices, connectivity)

  !call assert_equal(inv_mesh%get_nelements,             [1],       'get nelements')
  !call assert_equal(inv_mesh%get_connectivity(1),       [1, 2, 3], 'get connectivity')
  !call assert_equal(inv_mesh%get_vertices(2),           [1, 0, 0], 'get vertex 2')
  !call assert_equal(inv_mesh%get_valence(1),            [1],       'get valence')
  !call assert_comparable_real1d(inv_mesh%get_volume(1), [0.5],1e-8,'get volume')

  call inv_mesh%freeme()

end subroutine
!-----------------------------------------------------------------------------------------


end module
!=========================================================================================
