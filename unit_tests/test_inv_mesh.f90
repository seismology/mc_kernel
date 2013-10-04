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

  call inv_mesh%read_tet_mesh('vertices.TEST', 'facets.TEST')

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

  call inv_mesh%read_tet_mesh('vertices.TEST', 'facets.TEST')

  call inv_mesh%dump_mesh_xdmf('testmesh')

  call assert_file_exists('testmesh.xdmf', 'test xdmf dump')
  call assert_file_exists('testmesh_points.dat', 'test xdmf dump')
  call assert_file_exists('testmesh_grid.dat', 'test xdmf dump')

  ! tidy up
  open(newunit=myunit, file='testmesh.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='testmesh_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='testmesh_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable         :: datat(:,:)
  integer                           :: npoints, myunit, ierr

  call inv_mesh%read_tet_mesh('vertices.TEST', 'facets.TEST')
  call inv_mesh%init_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat(3,npoints))

  datat(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_data_snap(datat(1,:), 1, 'x')
  call inv_mesh%set_data_snap(datat(2,:), 2, 'x')
  call inv_mesh%set_data_snap(datat(3,:), 3, 'x')

  call inv_mesh%dump_mesh_data_xdmf('testdata')

  call assert_file_exists('testdata.xdmf', 'test xdmf data dump')
  call assert_file_exists('testdata_points.dat', 'test xdmf data dump')
  call assert_file_exists('testdata_grid.dat', 'test xdmf data dump')
  call assert_file_exists('testdata_data.dat', 'test xdmf data dump')

  ! tidy up
  open(newunit=myunit, file='testdata.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='testdata_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='testdata_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  open(newunit=myunit, file='testdata_data.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_data_dump2
  type(inversion_mesh_data_type)    :: inv_mesh
  real(kind=sp), allocatable         :: datat(:,:)
  integer                           :: npoints, myunit, ierr

  call inv_mesh%read_abaqus_mesh('circle.inp')
  call inv_mesh%init_data(3)

  npoints = inv_mesh%get_nvertices()
  allocate(datat(3,npoints))

  datat(:,:) = inv_mesh%get_vertices()
  call inv_mesh%set_data_snap(datat(1,:), 1, 'x')
  call inv_mesh%set_data_snap(datat(2,:), 2, 'x')
  call inv_mesh%set_data_snap(datat(3,:), 3, 'x')

  call inv_mesh%dump_mesh_data_xdmf('testcircle')

  call assert_file_exists('testcircle.xdmf', 'test xdmf data dump')
  call assert_file_exists('testcircle_points.dat', 'test xdmf data dump')
  call assert_file_exists('testcircle_grid.dat', 'test xdmf data dump')
  call assert_file_exists('testcircle_data.dat', 'test xdmf data dump')

  ! tidy up
  open(newunit=myunit, file='testcircle.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='testcircle_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='testcircle_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  open(newunit=myunit, file='testcircle_data.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_valence
  type(inversion_mesh_type)    :: inv_mesh

  call inv_mesh%read_tet_mesh('vertices.TEST', 'facets.TEST')

  call assert_equal(inv_mesh%get_valence(1), 2, 'valence of firt vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(2), 2, 'valence of firt vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(3), 2, 'valence of firt vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(4), 1, 'valence of firt vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(5), 1, 'valence of firt vertex in facets.TEST')

  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_get_connected_elements
  type(inversion_mesh_type)    :: inv_mesh
  integer, allocatable         :: cn_elems(:)

  call inv_mesh%read_tet_mesh('vertices.TEST', 'facets.TEST')

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


end module
!=========================================================================================
