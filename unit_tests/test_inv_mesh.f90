module test_inversion_mesh

  use inversion_mesh
  use ftnunit
  implicit none
  public
contains

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

subroutine test_mesh_dump
  type(inversion_mesh_type)    :: inv_mesh
  integer                      :: npoints, nelems

  call inv_mesh%read_tet_mesh('vertices.TEST', 'facets.TEST')

  call inv_mesh%dump_tet_mesh_xdmf('testmesh')

  call assert_file_exists('testmesh.xdmf', 'test xdmf dump')
  call assert_file_exists('testmesh_points.dat', 'test xdmf dump')
  call assert_file_exists('testmesh_grid.dat', 'test xdmf dump')
  
  call inv_mesh%freeme()
end subroutine

end module
