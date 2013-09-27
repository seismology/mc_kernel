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

  call assert_equal(npoints, 4, 'number of vertices in vertices.TEST')
  call assert_equal(nelems, 1, 'number of elements in facets.TEST')

  call inv_mesh%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_dump
  type(inversion_mesh_type)    :: inv_mesh
  integer                      :: myunit, ierr

  call inv_mesh%read_tet_mesh('vertices.TEST', 'facets.TEST')

  call inv_mesh%dump_tet_mesh_xdmf('testmesh')

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

  call inv_mesh%dump_tet_mesh_data_xdmf('testdata')

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



!
!  
!  call inv_mesh%freeme()

end module
!=========================================================================================
