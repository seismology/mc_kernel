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
  !open(newunit=myunit, file='unit_tests/testmesh_abaqus.xdmf', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testmesh_abaqus_points.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testmesh_abaqus_grid.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable        :: datat_node(:,:), datat_cell(:,:)
  integer                           :: npoints, nelements, myunit, ierr, i

  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST')

  call inv_mesh%init_node_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat_node(3,npoints))

  datat_node(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_node_data_snap(datat_node(1,:), 1, 'x')
  call inv_mesh%set_node_data_snap(datat_node(2,:), 2, 'x')
  call inv_mesh%set_node_data_snap(datat_node(3,:), 3, 'x')

  call inv_mesh%dump_node_data_xdmf('unit_tests/testdata')

  call assert_file_exists('unit_tests/testdata.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testdata_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testdata_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testdata_data.dat', 'test xdmf data dump')


  call inv_mesh%init_cell_data(3)

  nelements = inv_mesh%get_nelements()
  allocate(datat_cell(3,nelements))

  do i=1, nelements
    datat_cell(:,i) = inv_mesh%get_volume(i)
  enddo
  datat_cell(2,:) = datat_cell(2,:) + 1
  datat_cell(3,:) = datat_cell(3,:) + 2

  call inv_mesh%set_cell_data_snap(datat_cell(1,:), 1, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(2,:), 2, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(3,:), 3, 'x')

  call inv_mesh%dump_cell_data_xdmf('unit_tests/testcelldata')

  call assert_file_exists('unit_tests/testcelldata.xdmf', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testcelldata_points.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testcelldata_grid.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testcelldata_data.dat', 'test xdmf cell data dump')

  ! tidy up
  !open(newunit=myunit, file='unit_tests/testdata.xdmf', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testdata_points.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testdata_grid.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  !
  !open(newunit=myunit, file='unit_tests/testdata_data.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump2
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable        :: datat_node(:,:), datat_cell(:,:)
  integer                           :: npoints, nelements, myunit, ierr, i

  call inv_mesh%read_abaqus_mesh('unit_tests/circle.inp')
  call inv_mesh%init_node_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat_node(3,npoints))

  datat_node(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_node_data_snap(datat_node(1,:), 1, 'x')
  call inv_mesh%set_node_data_snap(datat_node(2,:), 2, 'x')
  call inv_mesh%set_node_data_snap(datat_node(3,:), 3, 'x')

  call inv_mesh%dump_node_data_xdmf('unit_tests/testcircle')

  call assert_file_exists('unit_tests/testcircle.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_data.dat', 'test xdmf data dump')


  call inv_mesh%init_cell_data(3)

  nelements = inv_mesh%get_nelements()
  allocate(datat_cell(3,nelements))

  do i=1, nelements
    datat_cell(:,i) = inv_mesh%get_volume(i)
  enddo
  datat_cell(2,:) = datat_cell(2,:) + 1
  datat_cell(3,:) = datat_cell(3,:) + 2

  call inv_mesh%set_cell_data_snap(datat_cell(1,:), 1, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(2,:), 2, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(3,:), 3, 'x')

  call inv_mesh%dump_cell_data_xdmf('unit_tests/testcellcircle')

  call assert_file_exists('unit_tests/testcellcircle.xdmf', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testcellcircle_points.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testcellcircle_grid.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testcellcircle_data.dat', 'test xdmf cell data dump')

  ! tidy up
  !open(newunit=myunit, file='unit_tests/testcircle.xdmf', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testcircle_points.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testcircle_grid.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  !
  !open(newunit=myunit, file='unit_tests/testcircle_data.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump3
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable        :: datat_node(:,:), datat_cell(:,:)
  integer                           :: npoints, nelements, myunit, ierr, i

  call inv_mesh%read_abaqus_mesh('unit_tests/circle_quad2.inp')
  call inv_mesh%init_node_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat_node(3,npoints))

  datat_node(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_node_data_snap(datat_node(1,:), 1, 'x')
  call inv_mesh%set_node_data_snap(datat_node(2,:), 2, 'x')
  call inv_mesh%set_node_data_snap(datat_node(3,:), 3, 'x')

  call inv_mesh%dump_node_data_xdmf('unit_tests/testcircle_quad')

  call assert_file_exists('unit_tests/testcircle_quad.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_quad_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_quad_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testcircle_quad_data.dat', 'test xdmf data dump')


  call inv_mesh%init_cell_data(3)

  nelements = inv_mesh%get_nelements()
  allocate(datat_cell(3,nelements))

  do i=1, nelements
    datat_cell(:,i) = inv_mesh%get_volume(i)
  enddo
  datat_cell(2,:) = datat_cell(2,:) + 1
  datat_cell(3,:) = datat_cell(3,:) + 2

  call inv_mesh%set_cell_data_snap(datat_cell(1,:), 1, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(2,:), 2, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(3,:), 3, 'x')

  call inv_mesh%dump_cell_data_xdmf('unit_tests/testcellcircle_quad')

  call assert_file_exists('unit_tests/testcellcircle_quad.xdmf', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testcellcircle_quad_points.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testcellcircle_quad_grid.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testcellcircle_quad_data.dat', 'test xdmf cell data dump')

  ! tidy up
  !open(newunit=myunit, file='unit_tests/testcircle_quad.xdmf', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testcircle_quad_points.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testcircle_quad_grid.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  !
  !open(newunit=myunit, file='unit_tests/testcircle_quad_data.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump4
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable        :: datat_node(:,:), datat_cell(:,:)
  integer                           :: npoints, nelements, myunit, ierr, i

  call inv_mesh%read_abaqus_mesh('unit_tests/sphere.inp')
  call inv_mesh%init_node_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat_node(3,npoints))

  datat_node(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_node_data_snap(datat_node(1,:), 1, 'x')
  call inv_mesh%set_node_data_snap(datat_node(2,:), 2, 'x')
  call inv_mesh%set_node_data_snap(datat_node(3,:), 3, 'x')

  call inv_mesh%dump_node_data_xdmf('unit_tests/testsphere')

  call assert_file_exists('unit_tests/testsphere.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testsphere_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testsphere_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testsphere_data.dat', 'test xdmf data dump')


  call inv_mesh%init_cell_data(3)

  nelements = inv_mesh%get_nelements()
  allocate(datat_cell(3,nelements))

  do i=1, nelements
    datat_cell(:,i) = inv_mesh%get_volume(i)
  enddo
  datat_cell(2,:) = datat_cell(2,:) + 1
  datat_cell(3,:) = datat_cell(3,:) + 2

  call inv_mesh%set_cell_data_snap(datat_cell(1,:), 1, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(2,:), 2, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(3,:), 3, 'x')

  call inv_mesh%dump_cell_data_xdmf('unit_tests/testspherecell')

  call assert_file_exists('unit_tests/testspherecell.xdmf', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testspherecell_points.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testspherecell_grid.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testspherecell_data.dat', 'test xdmf cell data dump')


  ! tidy up
  !open(newunit=myunit, file='unit_tests/testsphere.xdmf', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testsphere_points.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testsphere_grid.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  !
  !open(newunit=myunit, file='unit_tests/testsphere_data.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump5
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable        :: datat_node(:,:), datat_cell(:,:)
  integer                           :: npoints, nelements, myunit, ierr, i

  call inv_mesh%read_abaqus_mesh('unit_tests/tetrahedron.inp')
  call inv_mesh%init_node_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat_node(3,npoints))

  datat_node(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_node_data_snap(datat_node(1,:), 1, 'x')
  call inv_mesh%set_node_data_snap(datat_node(2,:), 2, 'x')
  call inv_mesh%set_node_data_snap(datat_node(3,:), 3, 'x')

  call inv_mesh%dump_node_data_xdmf('unit_tests/testtetrahedrons')

  call assert_file_exists('unit_tests/testtetrahedrons.xdmf', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testtetrahedrons_points.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testtetrahedrons_grid.dat', 'test xdmf data dump')
  call assert_file_exists('unit_tests/testtetrahedrons_data.dat', 'test xdmf data dump')


  call inv_mesh%init_cell_data(3)

  nelements = inv_mesh%get_nelements()
  allocate(datat_cell(3,nelements))

  do i=1, nelements
    datat_cell(:,i) = inv_mesh%get_volume(i)
  enddo
  datat_cell(2,:) = datat_cell(2,:) + 1
  datat_cell(3,:) = datat_cell(3,:) + 2

  call inv_mesh%set_cell_data_snap(datat_cell(1,:), 1, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(2,:), 2, 'x')
  call inv_mesh%set_cell_data_snap(datat_cell(3,:), 3, 'x')

  call inv_mesh%dump_cell_data_xdmf('unit_tests/testtetrahedronscell')

  call assert_file_exists('unit_tests/testtetrahedronscell.xdmf', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testtetrahedronscell_points.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testtetrahedronscell_grid.dat', 'test xdmf cell data dump')
  call assert_file_exists('unit_tests/testtetrahedronscell_data.dat', 'test xdmf cell data dump')


  ! tidy up
  !open(newunit=myunit, file='unit_tests/testtetrahedrons.xdmf', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testtetrahedrons_points.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')

  !open(newunit=myunit, file='unit_tests/testtetrahedrons_grid.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  !
  !open(newunit=myunit, file='unit_tests/testtetrahedrons_data.dat', iostat=ierr)
  !if (ierr == 0) close(myunit, status='delete')
  
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

  integer                    :: nelements_out, connectivity_out(1,3), valence_out
  real(kind=dp)              :: vertex_out(1,3), volume_out

  vertices(:,1) = [0, 0, 0]
  vertices(:,2) = [1, 0, 0]
  vertices(:,3) = [0, 0, 1]

  connectivity(1,:) = [1, 2, 3]

  call inv_mesh%initialize_mesh(1, vertices, connectivity)

  !nelements_out = inv_mesh%get_nelements()
  !call assert_equal(nelements_out,          1,       'get nelements')

  !connectivity_out = inv_mesh%get_element(1)
  !call assert_equal(connectivity_out,       [1, 2, 3], 'get connectivity')

  !vertex_out = inv_mesh%get_vertices(2)
  !call assert_comparable_real1d(vertex_out, [1, 0, 0], 1e-8, 'get vertex 2')

  !valence_out = inv_mesh%get_valence(1)
  !call assert_equal(valence_out,            [1],       'get valence')

  !volume_out = inv_mesh%get_volume(1)
  !call assert_comparable_real1d(volume_out, [0.5],     1e-8, 'get volume')

  call inv_mesh%freeme()

end subroutine
!-----------------------------------------------------------------------------------------


end module
!=========================================================================================
