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

  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  npoints = inv_mesh%get_nvertices()
  nelems = inv_mesh%get_nelements()

  call assert_equal(npoints, 5, 'number of vertices in vertices.TEST')
  call assert_equal(nelems, 2, 'number of elements in facets.TEST')

  call inv_mesh%freeme()
end subroutine test_mesh_read
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_dump
  type(inversion_mesh_type)    :: inv_mesh
  integer                      :: myunit, ierr

  ! tidy up
  open(newunit=myunit, file='unit_tests/testmesh.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testmesh_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests/testmesh_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')
  
  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  call inv_mesh%dump_mesh_xdmf('unit_tests_output/testmesh')

  call assert_file_exists('unit_tests_output/testmesh.xdmf', 'test xdmf dump')
  call assert_file_exists('unit_tests_output/testmesh_points.dat', 'test xdmf dump')
  call assert_file_exists('unit_tests_output/testmesh_grid.dat', 'test xdmf dump')

  call inv_mesh%freeme()
end subroutine test_mesh_dump
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_dump2
  type(inversion_mesh_type)    :: inv_mesh
  integer                      :: myunit, ierr

  ! tidy up
  open(newunit=myunit, file='unit_tests_output/testmesh_abaqus.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests_output/testmesh_abaqus_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests_output/testmesh_abaqus_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  call inv_mesh%read_abaqus_mesh('unit_tests/tetrahedron.inp','onvertices')

  call inv_mesh%dump_mesh_xdmf('unit_tests_output/testmesh_abaqus')

  call assert_file_exists('unit_tests_output/testmesh_abaqus.xdmf', 'test xdmf dump')
  call assert_file_exists('unit_tests_output/testmesh_abaqus_points.dat', 'test xdmf dump')
  call assert_file_exists('unit_tests_output/testmesh_abaqus_grid.dat', 'test xdmf dump')
  
  call inv_mesh%freeme()
end subroutine test_mesh_dump2
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mesh_dump3
  type(inversion_mesh_type)    :: inv_mesh
  integer                      :: myunit, ierr

  ! tidy up
  open(newunit=myunit, file='unit_tests_output/testmesh_abaqus_merge.xdmf', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests_output/testmesh_abaqus_merge_points.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  open(newunit=myunit, file='unit_tests_output/testmesh_abaqus_merge_grid.dat', iostat=ierr)
  if (ierr == 0) close(myunit, status='delete')

  call inv_mesh%read_abaqus_mesh('unit_tests/test_merge.inp','onvertices')

  call inv_mesh%dump_mesh_xdmf('unit_tests_output/testmesh_abaqus_merge')

  call assert_file_exists('unit_tests_output/testmesh_abaqus_merge.xdmf', 'test xdmf dump')
  call assert_file_exists('unit_tests_output/testmesh_abaqus_merge_points.dat', 'test xdmf dump')
  call assert_file_exists('unit_tests_output/testmesh_abaqus_merge_grid.dat', 'test xdmf dump')
  
  call inv_mesh%freeme()
end subroutine test_mesh_dump3
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_append_variable
  type (inversion_mesh_variable_type), allocatable  :: variable(:)
  character(len=16), dimension(:), allocatable      :: entry_names_ref_1
  character(len=16), dimension(:), allocatable      :: entry_names_ref_2
  character(len=16), dimension(:), allocatable      :: entry_names_ref_3

  entry_names_ref_1 = ['entry1', 'entry2']
  ! Test initialization of inversion_mesh_data_type variable
  call append_variable(variable, var_name = 'variable_test', &  
                                 nentries = 2,               &
                                 istime   = .false.,         &
                                 entry_names = entry_names_ref_1)

  call assert_true(allocated(variable), 'variable is allocated')
  call assert_equal(size(variable), 1,  'variable has size 1')
  call assert_equal(variable(1)%nentries, 2,  'First datum has has size 2')
  call assert_false(variable(1)%istime,    'variable is not a time')
  call assert_true(variable(1)%entry_names.eq.entry_names_ref_1, &
                   'Entry names are correct: ')

  entry_names_ref_2 = ['entry1', 'entry2', 'entry3', 'entry4']
  ! Test initialization of inversion_mesh_data_type variable
  call append_variable(variable, var_name = 'variable_test_2', &  
                                 nentries = 4,                 &
                                 istime   = .false.,           &
                                 entry_names = entry_names_ref_2)

  call assert_true(allocated(variable), 'variable is allocated')
  call assert_equal(size(variable), 2,  'variable has size 2')
  
  ! First verify that existing entries are not compromised
  call assert_equal(variable(1)%nentries, 2,  'First datum has has size 2')
  call assert_false(variable(1)%istime,    'variable is not a time')
  call assert_true(variable(1)%entry_names.eq.entry_names_ref_1, &
                   'Entry names are correct: ')

  ! Now check the new entries               
  call assert_equal(variable(2)%nentries, 4,  'Second datum has has size 3')
  call assert_false(variable(2)%istime,    'variable is not a time')
  call assert_true(variable(2)%entry_names.eq.entry_names_ref_2, &
                   'Entry names are correct: ')

  entry_names_ref_3 = ['entry1']
  ! Test initialization of inversion_mesh_data_type variable
  call append_variable(variable, var_name = 'variable_test_3', &  
                                 nentries = 1,                 &
                                 istime   = .false.,           &
                                 entry_names = entry_names_ref_3)

  call assert_true(allocated(variable), 'variable is allocated')
  call assert_equal(size(variable), 3,  'variable has size 3')

  ! First verify that existing entries are not compromised
  call assert_equal(variable(1)%nentries, 2,  'First datum has has size 2')
  call assert_false(variable(1)%istime,    'variable is not a time')
  call assert_true(variable(1)%entry_names.eq.entry_names_ref_1, &
                   'Entry names are correct: ')

  call assert_equal(variable(2)%nentries, 4,  'Second datum has has size 3')
  call assert_false(variable(2)%istime,    'variable is not a time')
  call assert_true(variable(2)%entry_names.eq.entry_names_ref_2, &
                   'Entry names are correct: ')

  ! Now check the new entries               
  call assert_equal(variable(3)%nentries, 1,  'Third datum has has size 3')
  call assert_false(variable(3)%istime,    'variable is not a time')
  call assert_true(variable(3)%entry_names.eq.entry_names_ref_3, &
                   'Entry names are correct: ')



end subroutine test_append_variable
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_init_node_data
  use netcdf
  type(inversion_mesh_data_type)    :: inv_mesh
  integer                           :: nf_status, ncid, grp_ncid, dim_len, dimid
  character(len=nf90_max_name)      :: dim_name


  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  call inv_mesh%init_node_data()
  nf_status = nf90_close(ncid = inv_mesh%ncid)

  nf_status = nf90_open(ncid = ncid, path=inv_mesh%filename_tmp_out, mode=NF90_NOWRITE)
  call assert_equal(nf_status, NF90_NOERR, &
                    'Error at nf90_open: '//trim(nf90_strerror(nf_status)))

  nf_status = nf90_inq_ncid(ncid = ncid, name='node_data', grp_ncid = grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, &
                    'Error at Group inquiry: '//trim(nf90_strerror(nf_status)))

  nf_status = nf90_close(ncid = ncid)

end subroutine test_init_node_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_init_mixed_data
  use netcdf
  type(inversion_mesh_data_type)    :: inv_mesh
  integer                           :: nf_status, ncid, grp_ncid, dim_len, dimid
  character(len=nf90_max_name)      :: dim_name


  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  call inv_mesh%init_cell_data()

  call inv_mesh%init_node_data()

  call inv_mesh%init_cell_data()

  call inv_mesh%init_node_data()

  nf_status = nf90_close(ncid = inv_mesh%ncid)

  nf_status = nf90_open(ncid = ncid, path=inv_mesh%filename_tmp_out, mode=NF90_NOWRITE)
  call assert_equal(nf_status, NF90_NOERR, &
                    'Error at nf90_open: '//trim(nf90_strerror(nf_status)))

  nf_status = nf90_inq_ncid(ncid = ncid, name='cell_data', grp_ncid = grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, &
                    'Error at Cell group inquiry: '//trim(nf90_strerror(nf_status)))

  nf_status = nf90_inq_ncid(ncid = ncid, name='node_data', grp_ncid = grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, &
                    'Error at Node group inquiry: '//trim(nf90_strerror(nf_status)))

  nf_status = nf90_close(ncid = ncid)

end subroutine test_init_mixed_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_init_cell_data
  use netcdf
  type(inversion_mesh_data_type)    :: inv_mesh
  integer                           :: nf_status, ncid, grp_ncid, dim_len, dimid
  character(len=nf90_max_name)      :: dim_name


  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  call inv_mesh%init_cell_data()
  nf_status = nf90_close(ncid = inv_mesh%ncid)

  nf_status = nf90_open(ncid = ncid, path=inv_mesh%filename_tmp_out, mode=NF90_NOWRITE)
  call assert_equal(nf_status, NF90_NOERR, &
                    'Error at nf90_open: '//trim(nf90_strerror(nf_status)))

  nf_status = nf90_inq_ncid(ncid = ncid, name='cell_data', grp_ncid = grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, &
                    'Error at Group inquiry: '//trim(nf90_strerror(nf_status)))

  call inv_mesh%freeme()

  nf_status = nf90_close(ncid = ncid)

end subroutine test_init_cell_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_set_mixed_data
  use netcdf
  use nc_routines
  type(inversion_mesh_data_type)    :: inv_mesh
  character(len=nf90_max_name)      :: variable_names(3)
  integer                           :: variable_length(3)

  integer                           :: nf_status, ncid, id, jd, grp_ncid
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), &
                                       testvar3_ref(:,:), testvar4_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_ref_node(:,:), testvar2_ref_node(:,:), &
                                       testvar3_ref_node(:,:), testvar4_ref_node(:,:)
  character(len=nf90_max_name)      :: dim_name
  character(len=16), dimension(:), allocatable  :: entry_names_cell_1, entry_names_cell_2, &
                                                   entry_names_cell_3, entry_names_cell_4
  character(len=16), dimension(:), allocatable  :: entry_names_node_1, entry_names_node_2, &
                                                   entry_names_node_3, entry_names_node_4



  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  variable_names  = ['variable_1', 'variable_2', 'variable_3']
  variable_length = [2, 3, 4]

  allocate(testvar1_ref(2,2)) ! Mesh has two cells
  call random_number(testvar1_ref)
  allocate(testvar2_ref(2,3)) ! Mesh has two cells
  call random_number(testvar2_ref)
  allocate(testvar3_ref(2,4)) ! Mesh has two cells
  call random_number(testvar3_ref)
  allocate(testvar4_ref(2,1)) ! Mesh has two cells
  call random_number(testvar4_ref)

  allocate(testvar1_ref_node(5,2)) ! Mesh has five points
  call random_number(testvar1_ref_node)
  allocate(testvar2_ref_node(5,3)) ! Mesh has five points
  call random_number(testvar2_ref_node)
  allocate(testvar3_ref_node(5,4)) ! Mesh has five points
  call random_number(testvar3_ref_node)
  allocate(testvar4_ref_node(5,1)) ! Mesh has five points
  call random_number(testvar4_ref_node)

  ! CELL DATA
  call inv_mesh%init_cell_data()

  ! Write whole variable at once
  entry_names_cell_1 = ['entry1', 'entry2']
  call inv_mesh%add_cell_variable('cell_variable_1', &
                                  nentries = 2,      &
                                  entry_names = entry_names_cell_1)
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_1',   &
                              values    = testvar1_ref)

  ! Write variable element-wise
  entry_names_cell_2 = ['entry1', 'entry2', 'entry3']
  call inv_mesh%add_cell_variable('cell_variable_2',         &
                                  nentries    = 3,           &
                                  entry_names = entry_names_cell_2)
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_2',          &
                              values    = testvar2_ref(1:1,:),   &
                              ielement  = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_2',          &
                              values    = testvar2_ref(2:2,:),   &
                              ielement  = [2, 2])
  
  ! Just to add confusion, init both once more
  call inv_mesh%init_cell_data()
  call inv_mesh%init_node_data()

  ! Write variable snap-wise
  entry_names_cell_3 = ['entry1', 'entry2', 'entry3', 'entry4']
  call inv_mesh%add_cell_variable('cell_variable_3', &
                                  nentries = 4,      &
                                  entry_names = entry_names_cell_3)
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar3_ref(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar3_ref(:,2:4),   &
                              ientry    = [2, 4])

  ! Add one-entry variable at the end
  entry_names_cell_4 = ['entry1']
  call inv_mesh%add_cell_variable('cell_variable_4', &
                                  nentries = 1,      &
                                  entry_names = entry_names_cell_4)
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_4',          &
                              values    = testvar4_ref(:,:),   &
                              ientry    = [1, 1])

  ! END OF CELL DATA                          

  ! NODE DATA                          
  call inv_mesh%init_node_data()

  ! Write whole variable at once
  entry_names_node_1 = ['entry1', 'entry2']
  call inv_mesh%add_node_variable('node_variable_1', &
                                  nentries = 2,      &
                                  entry_names = entry_names_node_1)
  call inv_mesh%add_node_data(var_name  = 'node_variable_1',   &
                              values    = testvar1_ref_node)

  ! Write variable element-wise
  entry_names_node_2 = ['entry1', 'entry2', 'entry3']
  call inv_mesh%add_node_variable('node_variable_2',         &
                                  nentries    = 3,           &
                                  entry_names = entry_names_node_2)
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar2_ref_node(1:1,:),   &
                              ielement  = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar2_ref_node(2:2,:),   &
                              ielement  = [2, 2])
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar2_ref_node(3:5,:),   &
                              ielement  = [3, 5])

  ! Just to add confusion, init both once more
  call inv_mesh%init_cell_data()
  call inv_mesh%init_node_data()

  ! Write variable snap-wise
  entry_names_node_3 = ['entry1', 'entry2', 'entry3', 'entry4']
  call inv_mesh%add_node_variable('node_variable_3', &
                                  nentries = 4,      &
                                  entry_names = entry_names_node_3)
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar3_ref_node(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar3_ref_node(:,2:4),   &
                              ientry    = [2, 4])

  ! Add one-entry variable at the end
  entry_names_node_4 = ['entry1']
  call inv_mesh%add_node_variable('node_variable_4', &
                                  nentries = 1,      &
                                  entry_names = entry_names_node_4)
  call inv_mesh%add_node_data(var_name  = 'node_variable_4',          &
                              values    = testvar4_ref_node(:,:),   &
                              ientry    = [1, 1])

  ! END OF NODE DATA                           

  ! Test variables in inv_mesh object

  ! First cell variable, 2 entries
  call assert_equal(inv_mesh%variable_cell(1)%nentries, 2, 'First cell variable has has size 2')
  call assert_false(inv_mesh%variable_cell(1)%istime,      'First cell variable is not a time')
  call assert_true(inv_mesh%variable_cell(1)%entry_names.eq.entry_names_cell_1, &
                   'Entry names are correct: ')

  ! Second cell variable, 3 entries
  call assert_equal(inv_mesh%variable_cell(2)%nentries, 3, 'Second cell variable has has size 2')
  call assert_false(inv_mesh%variable_cell(2)%istime,      'Second cell variable is not a time')
  call assert_true(inv_mesh%variable_cell(2)%entry_names.eq.entry_names_cell_2, &
                   'Entry names are correct: ')

  ! Third cell variable, 4 entries
  call assert_equal(inv_mesh%variable_cell(3)%nentries, 4, 'Third cell variable has has size 4')
  call assert_false(inv_mesh%variable_cell(3)%istime,      'Third cell variable is not a time')
  call assert_true(inv_mesh%variable_cell(3)%entry_names.eq.entry_names_cell_3, &
                   'Entry names are correct: ')

  ! Fourth cell variable, 4 entries
  call assert_equal(inv_mesh%variable_cell(4)%nentries, 1, 'Fourth cell variable has has size 4')
  call assert_false(inv_mesh%variable_cell(4)%istime,      'Fourth cell variable is not a time')
  call assert_true(inv_mesh%variable_cell(4)%entry_names.eq.entry_names_cell_4, &
                   'Entry names are correct: ')

  ! First node variable, 2 entries
  call assert_equal(inv_mesh%variable_node(1)%nentries, 2, 'First node variable has has size 2')
  call assert_false(inv_mesh%variable_node(1)%istime,      'First node variable is not a time')
  call assert_true(inv_mesh%variable_node(1)%entry_names.eq.entry_names_node_1, &
                   'Entry names are correct: ')

  ! Second node variable, 3 entries
  call assert_equal(inv_mesh%variable_node(2)%nentries, 3, 'Second node variable has has size 2')
  call assert_false(inv_mesh%variable_node(2)%istime,      'Second node variable is not a time')
  call assert_true(inv_mesh%variable_node(2)%entry_names.eq.entry_names_node_2, &
                   'Entry names are correct: ')

  ! Third node variable, 4 entries
  call assert_equal(inv_mesh%variable_node(3)%nentries, 4, 'Third node variable has has size 4')
  call assert_false(inv_mesh%variable_node(3)%istime,      'Third node variable is not a time')
  call assert_true(inv_mesh%variable_node(3)%entry_names.eq.entry_names_node_3, &
                   'Entry names are correct: ')

  ! Fourth node variable, 4 entries
  call assert_equal(inv_mesh%variable_node(4)%nentries, 1, 'Fourth node variable has has size 4')
  call assert_false(inv_mesh%variable_node(4)%istime,      'Fourth node variable is not a time')
  call assert_true(inv_mesh%variable_node(4)%entry_names.eq.entry_names_node_4, &
                   'Entry names are correct: ')


end subroutine test_set_mixed_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_set_node_data_and_dump()
  use netcdf
  use nc_routines
  type(inversion_mesh_data_type)    :: inv_mesh
  character(len=nf90_max_name)      :: variable_names(3)
  integer                           :: variable_length(3)

  integer                           :: nf_status, ncid, id, jd, grp_ncid
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), testvar3_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_res(:,:), testvar2_res(:,:), testvar3_res(:,:)
  character(len=nf90_max_name)      :: dim_name


  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  variable_names  = ['variable_1', 'variable_2', 'variable_3']
  variable_length = [2, 3, 4]

  allocate(testvar1_ref(5,2)) ! Mesh has five vertices
  call random_number(testvar1_ref)

  allocate(testvar2_ref(5,3)) ! Mesh has five vertices
  call random_number(testvar2_ref)
  
  allocate(testvar3_ref(5,4)) ! Mesh has five vertices
  call random_number(testvar3_ref)

  call inv_mesh%init_node_data()

  ! Write whole variable at once
  call inv_mesh%add_node_variable('variable_1', &
                                  nentries = 2)
  call inv_mesh%add_node_data(var_name  = 'variable_1',   &
                              values    = testvar1_ref)

  ! Write variable element-wise
  call inv_mesh%add_node_variable('variable_2', &
                                  nentries = 3)
  call inv_mesh%add_node_data(var_name  = 'variable_2',          &
                              values    = testvar2_ref(1:1,:),   &
                              ielement  = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'variable_2',          &
                              values    = testvar2_ref(2:2,:),   &
                              ielement  = [2, 2])
  call inv_mesh%add_node_data(var_name  = 'variable_2',          &
                              values    = testvar2_ref(3:5,:),   &
                              ielement  = [3, 5])

  ! Write variable snap-wise
  call inv_mesh%add_node_variable('variable_3', &
                                  nentries = 4)
  call inv_mesh%add_node_data(var_name  = 'variable_3',          &
                              values    = testvar3_ref(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'variable_3',          &
                              values    = testvar3_ref(:,2:4),   &
                              ientry    = [2, 4])

  ! Write file to disk
  call inv_mesh%dump_data_xdmf('unit_tests_output/test_set_node_data')

  call nc_open_for_read(filename = 'unit_tests_output/test_set_node_data.nc', ncid = ncid)
  nf_status = nf90_inq_ncid(ncid=ncid, name='node_data', grp_ncid=grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, 'Group for node data has been created')

  ! Get variable 1
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'variable_1', & 
                         values  = testvar1_res, &
                         limits  = [0.,1.] )

  do id = 1, size(testvar1_ref, 1) 
    call assert_comparable(testvar1_ref(id,:), testvar1_res(id,:),  &
                           1e-7, 'Variable 1 retrieved successfully')
  end do

  ! Get variable 2
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'variable_2', & 
                         values  = testvar2_res, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar2_ref, 1) 
    call assert_comparable(testvar2_ref(id,:), testvar2_res(id,:),  &
                           1e-7, 'Variable 2 retrieved successfully')
  end do


  ! Get variable 3
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'variable_3', & 
                         values  = testvar3_res, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar3_ref, 1) 
    call assert_comparable(testvar3_ref(id,:), testvar3_res(id,:),  &
                           1e-7, 'Variable 3 retrieved successfully')
  end do

  nf_status = nf90_close(ncid = ncid)

end subroutine test_set_node_data_and_dump
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_set_cell_data_and_dump()
  use netcdf
  use nc_routines
  type(inversion_mesh_data_type)    :: inv_mesh
  character(len=nf90_max_name)      :: variable_names(3)
  integer                           :: variable_length(3)

  integer                           :: nf_status, ncid, id, jd, grp_ncid
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), testvar3_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_res(:,:), testvar2_res(:,:), testvar3_res(:,:)
  character(len=nf90_max_name)      :: dim_name


  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  variable_names  = ['variable_1', 'variable_2', 'variable_3']
  variable_length = [2, 3, 4]

  allocate(testvar1_ref(2,2)) ! Mesh has two cells
  call random_number(testvar1_ref)
  allocate(testvar2_ref(2,3)) ! Mesh has two cells
  call random_number(testvar2_ref)
  allocate(testvar3_ref(2,4)) ! Mesh has two cells
  call random_number(testvar3_ref)

  call inv_mesh%init_cell_data()

  ! Write whole variable at once
  call inv_mesh%add_cell_variable('variable_1', &
                                  nentries = 2)
  call inv_mesh%add_cell_data(var_name  = 'variable_1',   &
                              values    = testvar1_ref)

  ! Write variable element-wise
  call inv_mesh%add_cell_variable('variable_2', &
                                  nentries = 3)
  call inv_mesh%add_cell_data(var_name  = 'variable_2',          &
                              values    = testvar2_ref(1:1,:),   &
                              ielement  = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'variable_2',          &
                              values    = testvar2_ref(2:2,:),   &
                              ielement  = [2, 2])

  ! Write variable snap-wise
  call inv_mesh%add_cell_variable('variable_3', &
                                  nentries = 4)
  call inv_mesh%add_cell_data(var_name  = 'variable_3',          &
                              values    = testvar3_ref(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'variable_3',          &
                              values    = testvar3_ref(:,2:4),   &
                              ientry    = [2, 4])

  ! Write file to disk
  call inv_mesh%dump_data_xdmf('unit_tests_output/test_set_cell_data')

  call nc_open_for_read(filename = 'unit_tests_output/test_set_cell_data.nc', ncid = ncid)
  nf_status = nf90_inq_ncid(ncid=ncid, name='cell_data', grp_ncid=grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, 'Group for cell data has been created')

  ! Get variable 1
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'variable_1', & 
                         values  = testvar1_res, &
                         limits  = [0.,1.] )

  do id = 1, size(testvar1_ref, 2) 
    call assert_comparable(testvar1_ref(:,id), testvar1_res(:,id),  &
                           1e-7, 'Variable 1 retrieved successfully')
  end do

  ! Get variable 2
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'variable_2', & 
                         values  = testvar2_res, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar2_ref, 2) 
    call assert_comparable(testvar2_ref(:,id), testvar2_res(:,id),  &
                           1e-7, 'Variable 2 retrieved successfully')
  end do


  ! Get variable 3
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'variable_3', & 
                         values  = testvar3_res, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar3_ref, 2) 
    call assert_comparable(testvar3_ref(:,id), testvar3_res(:,id),  &
                           1e-7, 'Variable 3 retrieved successfully')
  end do

  nf_status = nf90_close(ncid = ncid)

end subroutine test_set_cell_data_and_dump
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_set_mixed_data_and_dump()
  use netcdf
  use nc_routines
  type(inversion_mesh_data_type)    :: inv_mesh
  character(len=nf90_max_name)      :: variable_names(3)
  integer                           :: variable_length(3)

  integer                           :: nf_status, ncid, id, jd, grp_ncid
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), testvar3_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_ref_node(:,:), testvar2_ref_node(:,:), testvar3_ref_node(:,:)
  real(kind=sp), allocatable        :: testvar1_res(:,:), testvar2_res(:,:), testvar3_res(:,:)
  real(kind=sp), allocatable        :: testvar1_res_node(:,:), testvar2_res_node(:,:), testvar3_res_node(:,:)
  character(len=nf90_max_name)      :: dim_name


  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  variable_names  = ['variable_1', 'variable_2', 'variable_3']
  variable_length = [2, 3, 4]

  allocate(testvar1_ref(2,2)) ! Mesh has two cells
  call random_number(testvar1_ref)
  allocate(testvar2_ref(2,3)) ! Mesh has two cells
  call random_number(testvar2_ref)
  allocate(testvar3_ref(2,4)) ! Mesh has two cells
  call random_number(testvar3_ref)

  allocate(testvar1_ref_node(5,2)) ! Mesh has five points
  call random_number(testvar1_ref_node)
  allocate(testvar2_ref_node(5,3)) ! Mesh has five points
  call random_number(testvar2_ref_node)
  allocate(testvar3_ref_node(5,4)) ! Mesh has five points
  call random_number(testvar3_ref_node)

  ! CELL DATA
  call inv_mesh%init_cell_data()

  ! Write whole variable at once
  call inv_mesh%add_cell_variable('cell_variable_1', &
                                  nentries = 2)
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_1',   &
                              values    = testvar1_ref)

  ! Write variable element-wise
  call inv_mesh%add_cell_variable('cell_variable_2', &
                                  nentries = 3)
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_2',          &
                              values    = testvar2_ref(1:1,:),   &
                              ielement  = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_2',          &
                              values    = testvar2_ref(2:2,:),   &
                              ielement  = [2, 2])
  
  ! Just to add confusion, init both once more
  call inv_mesh%init_cell_data()
  call inv_mesh%init_node_data()

  ! Write variable snap-wise
  call inv_mesh%add_cell_variable('cell_variable_3', &
                                  nentries = 4)
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar3_ref(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar3_ref(:,2:4),   &
                              ientry    = [2, 4])

  ! END OF CELL DATA                          

  ! NODE DATA                          
  call inv_mesh%init_node_data()

  ! Write whole variable at once
  call inv_mesh%add_node_variable('node_variable_1', &
                                  nentries = 2)
  call inv_mesh%add_node_data(var_name  = 'node_variable_1',   &
                              values    = testvar1_ref_node)

  ! Write variable element-wise
  call inv_mesh%add_node_variable('node_variable_2', &
                                  nentries = 3)
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar2_ref_node(1:1,:),   &
                              ielement  = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar2_ref_node(2:2,:),   &
                              ielement  = [2, 2])
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar2_ref_node(3:5,:),   &
                              ielement  = [3, 5])

  ! Just to add confusion, init both once more
  call inv_mesh%init_cell_data()
  call inv_mesh%init_node_data()

  ! Write variable snap-wise
  call inv_mesh%add_node_variable('node_variable_3', &
                                  nentries = 4)
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar3_ref_node(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar3_ref_node(:,2:4),   &
                              ientry    = [2, 4])

  ! END OF NODE DATA                           


  ! Write file to disk
  call inv_mesh%dump_data_xdmf('unit_tests_output/test_set_mixed_data')

  call nc_open_for_read(filename = 'unit_tests_output/test_set_mixed_data.nc', ncid = ncid)

  ! CHECK FOR CELL DATA
  nf_status = nf90_inq_ncid(ncid=ncid, name='cell_data', grp_ncid=grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, 'Group for cell data has been created')

  ! Get variable 1
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'cell_variable_1', & 
                         values  = testvar1_res, &
                         limits  = [0.,1.] )

  do id = 1, size(testvar1_ref, 2) 
    call assert_comparable(testvar1_ref(:,id), testvar1_res(:,id),  &
                           1e-7, 'Variable 1 retrieved successfully')
  end do

  ! Get variable 2
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'cell_variable_2', & 
                         values  = testvar2_res, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar2_ref, 2) 
    call assert_comparable(testvar2_ref(:,id), testvar2_res(:,id),  &
                           1e-7, 'Variable 2 retrieved successfully')
  end do


  ! Get variable 3
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'cell_variable_3', & 
                         values  = testvar3_res, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar3_ref, 2) 
    call assert_comparable(testvar3_ref(:,id), testvar3_res(:,id),  &
                           1e-7, 'Variable 3 retrieved successfully')
  end do

  deallocate(testvar1_res)
  deallocate(testvar2_res)
  deallocate(testvar3_res)


  ! CHECK FOR NODE DATA
  nf_status = nf90_inq_ncid(ncid=ncid, name='node_data', grp_ncid=grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, 'Group for node data has been created')

  ! Get variable 1
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'node_variable_1', & 
                         values  = testvar1_res_node, &
                         limits  = [0.,1.] )

  do id = 1, size(testvar1_ref, 2) 
    call assert_comparable(testvar1_ref_node(:,id), testvar1_res_node(:,id),  &
                           1e-7, 'Variable 1 retrieved successfully')
  end do

  ! Get variable 2
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'node_variable_2', & 
                         values  = testvar2_res_node, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar2_ref, 2) 
    call assert_comparable(testvar2_ref_node(:,id), testvar2_res_node(:,id),  &
                           1e-7, 'Variable 2 retrieved successfully')
  end do


  ! Get variable 3
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'node_variable_3', & 
                         values  = testvar3_res_node, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar3_ref, 2) 
    call assert_comparable(testvar3_ref_node(:,id), testvar3_res_node(:,id),  &
                           1e-7, 'Variable 3 retrieved successfully')
  end do

  nf_status = nf90_close(ncid = ncid)

end subroutine test_set_mixed_data_and_dump
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_valence
  type(inversion_mesh_type) :: inv_mesh

  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  call assert_equal(inv_mesh%get_valence(1), 2, 'valence of first vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(2), 2, 'valence of second vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(3), 2, 'valence of third vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(4), 1, 'valence of fourth vertex in facets.TEST')
  call assert_equal(inv_mesh%get_valence(5), 1, 'valence of fifth vertex in facets.TEST')

  call inv_mesh%freeme()

end subroutine test_valence
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_get_connected_elements
  type(inversion_mesh_type) :: inv_mesh
  integer, allocatable      :: cn_elems(:)

  call inv_mesh%read_tet_mesh('unit_tests/vertices.TEST', 'unit_tests/facets.TEST', &
                              'onvertices')

  ! uncomment to NOT use automatic allocation
  allocate(cn_elems(inv_mesh%get_valence(1)))
  cn_elems = inv_mesh%get_connected_elements(1)

  ! somewhat redundant two tests, for now to test ftnunit :)
  call assert_true(.not. any(cn_elems == -1), 'get connected elements')
  call assert_true(cn_elems /= -1, 'get connected elements')
  call assert_equal(inv_mesh%get_connected_elements(1), [1 ,2], 'get connected elements')
  call assert_equal(inv_mesh%get_connected_elements(4), [1], 'get connected elements')
  call assert_equal(inv_mesh%get_connected_elements(5), [2], 'get connected elements')

  call inv_mesh%freeme()

end subroutine test_get_connected_elements
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_random_points_triangle_mesh
  use tetrahedra, only        : point_in_triangle_3d
  type(inversion_mesh_type)  :: inv_mesh
  integer, parameter         :: ndim = 3, nvertices = 4
  real(kind=dp)              :: vertices(ndim, nvertices), offset(ndim)
  integer                    :: connectivity(3, 1)

  integer, parameter         :: npoints = 1000, ntriangle = 10
  integer                    :: nbasisfuncs_per_elem, itriangle, i
  real(kind=dp)              :: points(ndim, npoints)
  logical                    :: isintriangle(npoints)
  logical                    :: isinplane(npoints)

  do itriangle = 1, ntriangle

      call random_number(vertices)
      call random_number(offset)

      do i = 1,3
         vertices(1,:) = (vertices(1,:)-0.5)*2 + (offset(1)-0.5)*20
         vertices(2,:) = (vertices(2,:)-0.5)*2 + (offset(2)-0.5)*20
         vertices(3,:) = (vertices(3,:)-0.5)*2 + (offset(3)-0.5)*20
      end do

      connectivity(:,1) = [1, 2, 3]

      nbasisfuncs_per_elem = 1
      call inv_mesh%initialize_mesh(1, vertices, connectivity, nbasisfuncs_per_elem)
      
      points = inv_mesh%generate_random_points(1, npoints, .true.)
      
      isintriangle = point_in_triangle_3d(vertices, points, isinplane)

      call assert_true(isinplane,    'Points are in plane of triangle')
      call assert_true(isintriangle, 'Points are in triangle')

      call inv_mesh%freeme()
  end do


end subroutine test_random_points_triangle_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_initialize_mesh
  type(inversion_mesh_type)  :: inv_mesh
  real(kind=dp)              :: vertices(3,3)
  integer                    :: connectivity(3,1)

  integer                    :: nelements_out, connectivity_out(1,3)
  integer                    :: nbasisfuncs_per_elem, valence_out
  real(kind=sp)              :: vertex_out(3,3), volume_out

  vertices(:,1) = [0, 0, 0]
  vertices(:,2) = [1, 0, 0]
  vertices(:,3) = [0, 0, 1]

  connectivity(:,1) = [1, 2, 3]

  nbasisfuncs_per_elem = 8

  call inv_mesh%initialize_mesh(1, vertices, connectivity, nbasisfuncs_per_elem)

  nelements_out = inv_mesh%get_nelements()
  call assert_equal(nelements_out,          1,       'get nelements')

  connectivity_out = transpose(inv_mesh%get_connectivity())
  call assert_equal_int1d(connectivity_out(1,:),       [1, 2, 3], 'get connectivity')

  vertex_out = real(inv_mesh%get_vertices(), kind=sp)
  call assert_comparable_real1d(vertex_out(:,2), [1., 0., 0.], 1e-8, 'get vertex 2')

  valence_out = inv_mesh%get_valence(1)
  call assert_equal(valence_out,            1,       'get valence')

  volume_out = real(inv_mesh%get_volume(1), kind=sp)
  call assert_comparable_real(volume_out, 0.5,     1e-8, 'get volume')

  call inv_mesh%freeme()

end subroutine test_initialize_mesh
!-----------------------------------------------------------------------------------------


end module test_inversion_mesh
!=========================================================================================
