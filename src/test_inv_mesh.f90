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

  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
                              'onvertices')

  npoints = inv_mesh%get_nvertices()
  nelems = inv_mesh%get_nelements()

  call assert_equal(npoints, 5, 'number of vertices in vertices.TEST')
  call assert_equal(nelems, 2, 'number of elements in facets.TEST')

  call inv_mesh%freeme()
end subroutine test_mesh_read
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_append_variable
  type (inversion_mesh_variable_type), allocatable  :: variable(:)
  character(len=80)                                 :: entry_names_ref_1(2)
  character(len=80)                                 :: entry_names_ref_2(4)
  character(len=80)                                 :: entry_names_ref_3(1)

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
  !print *, variable(1)%entry_names, entry_names_ref_1 
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
  !print *, variable(1)%entry_names, entry_names_ref_1 
  call assert_true(variable(1)%entry_names.eq.entry_names_ref_1, &
                   'Entry names are correct: ')

  ! Now check the new entries               
  call assert_equal(variable(2)%nentries, 4,  'Second datum has has size 3')
  call assert_false(variable(2)%istime,    'variable is not a time')
  !print *, variable(2)%entry_names, entry_names_ref_2 
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
  !print *, variable(1)%entry_names, entry_names_ref_1 
  call assert_true(variable(1)%entry_names.eq.entry_names_ref_1, &
                   'Entry names are correct: ')

  call assert_equal(variable(2)%nentries, 4,  'Second datum has has size 3')
  call assert_false(variable(2)%istime,    'variable is not a time')
  !print *, variable(2)%entry_names, entry_names_ref_2 
  call assert_true(variable(2)%entry_names.eq.entry_names_ref_2, &
                   'Entry names are correct: ')

  ! Now check the new entries               
  call assert_equal(variable(3)%nentries, 1,  'Third datum has has size 3')
  call assert_false(variable(3)%istime,    'variable is not a time')
  !print *, variable(3)%entry_names, entry_names_ref_3 
  call assert_true(variable(3)%entry_names.eq.entry_names_ref_3, &
                   'Entry names are correct: ')



end subroutine test_append_variable
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_init_node_data
  use netcdf
  type(inversion_mesh_data_type)    :: inv_mesh
  integer                           :: nf_status, ncid, grp_ncid


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
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
  integer                           :: nf_status, ncid, grp_ncid


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
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
  integer                           :: nf_status, ncid, grp_ncid


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
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
  character(len=NF90_MAX_NAME)      :: variable_names(3)
  integer                           :: variable_length(3)

  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), &
                                       testvar3_ref(:,:), testvar4_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_ref_node(:,:), testvar2_ref_node(:,:), &
                                       testvar3_ref_node(:,:), testvar4_ref_node(:,:)
  character(len=16)                 :: entry_names_cell_1(2), entry_names_cell_2(3), &
                                       entry_names_cell_3(4), entry_names_cell_4(1)
  character(len=16)                 :: entry_names_node_1(2), entry_names_node_2(3), &
                                       entry_names_node_3(4), entry_names_node_4(1)



  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
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
  character(len=NF90_MAX_NAME)      :: variable_names(3)
  integer                           :: variable_length(3)

  integer                           :: nf_status, ncid, id, grp_ncid
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), testvar3_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_res(:,:), testvar2_res(:,:), testvar3_res(:,:)


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
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
  call inv_mesh%dump_data_xdmf('./output/test_set_node_data')
  call inv_mesh%free_node_and_cell_data()

  call nc_open_for_read(filename = './output/test_set_node_data.nc', ncid = ncid)
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
  character(len=NF90_MAX_NAME)      :: variable_names(3)
  integer                           :: variable_length(3)

  integer                           :: nf_status, ncid, id, grp_ncid
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), testvar3_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_res(:,:), testvar2_res(:,:), testvar3_res(:,:)


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
                              'onvertices')

  variable_names  = ['variable_1', 'variable_2', 'variable_3']
  variable_length = [2, 3, 4]

  allocate(testvar1_ref(2,variable_length(1))) ! Mesh has two cells
  call random_number(testvar1_ref)
  allocate(testvar2_ref(2,variable_length(2))) ! Mesh has two cells
  call random_number(testvar2_ref)
  allocate(testvar3_ref(2,variable_length(3))) ! Mesh has two cells
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
  call inv_mesh%dump_data_xdmf('./output/test_set_cell_data')

  call nc_open_for_read(filename = './output/test_set_cell_data.nc', ncid = ncid)
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
  character(len=NF90_MAX_NAME)      :: variable_names(3)
  integer                           :: variable_length(3)

  integer                           :: nf_status, ncid, id, grp_ncid
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), testvar3_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_ref_node(:,:), testvar2_ref_node(:,:), testvar3_ref_node(:,:)
  real(kind=sp), allocatable        :: testvar1_res(:,:), testvar2_res(:,:), testvar3_res(:,:)
  real(kind=sp), allocatable        :: testvar1_res_node(:,:), testvar2_res_node(:,:), testvar3_res_node(:,:)


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
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
  call inv_mesh%dump_data_xdmf('./output/test_set_mixed_data')

  call nc_open_for_read(filename = './output/test_set_mixed_data.nc', ncid = ncid)

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
! Test whether two datasets for the same mesh can be written into two separate files
! Needed for the wavefield output, which is written into another file than the kernels.
subroutine test_set_mixed_data_and_dump_into_two_files()
  use netcdf
  use nc_routines
  type(inversion_mesh_data_type)    :: inv_mesh
  character(len=NF90_MAX_NAME)      :: variable_names(3)
  integer                           :: variable_length(6)

  integer                           :: nf_status, ncid, id, grp_ncid

  ! Variables for first output file
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), testvar3_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_ref_node(:,:), testvar2_ref_node(:,:), testvar3_ref_node(:,:)
  real(kind=sp), allocatable        :: testvar1_res(:,:), testvar2_res(:,:), testvar3_res(:,:)
  real(kind=sp), allocatable        :: testvar1_res_node(:,:), testvar2_res_node(:,:), testvar3_res_node(:,:)

  ! Variables for second output file
  real(kind=sp), allocatable        :: testvar4_ref(:,:), testvar5_ref(:,:), testvar6_ref(:,:)
  real(kind=sp), allocatable        :: testvar4_ref_node(:,:), testvar5_ref_node(:,:), testvar6_ref_node(:,:)
  real(kind=sp), allocatable        :: testvar4_res(:,:), testvar5_res(:,:), testvar6_res(:,:)
  real(kind=sp), allocatable        :: testvar4_res_node(:,:), testvar5_res_node(:,:), testvar6_res_node(:,:)


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
                              'onvertices')

  variable_names  = ['variable_1', 'variable_2', 'variable_3']
  variable_length = [2, 3, 4, 5, 6, 7]


  ! Variables for first output file
  allocate(testvar1_ref(2,variable_length(1))) ! Mesh has two cells
  call random_number(testvar1_ref)
  allocate(testvar2_ref(2,variable_length(2))) ! Mesh has two cells
  call random_number(testvar2_ref)
  allocate(testvar3_ref(2,variable_length(3))) ! Mesh has two cells
  call random_number(testvar3_ref)

  allocate(testvar1_ref_node(5,variable_length(1))) ! Mesh has five points
  call random_number(testvar1_ref_node)
  allocate(testvar2_ref_node(5,variable_length(2))) ! Mesh has five points
  call random_number(testvar2_ref_node)
  allocate(testvar3_ref_node(5,variable_length(3))) ! Mesh has five points
  call random_number(testvar3_ref_node)


  ! Variables for second output file
  variable_names  = ['variable_1', 'variable_2', 'variable_3']

  allocate(testvar4_ref(2,variable_length(4))) ! Mesh has two cells
  call random_number(testvar4_ref)
  allocate(testvar5_ref(2,variable_length(5))) ! Mesh has two cells
  call random_number(testvar5_ref)
  allocate(testvar6_ref(2,variable_length(6))) ! Mesh has two cells
  call random_number(testvar6_ref)

  allocate(testvar4_ref_node(5,variable_length(4))) ! Mesh has five points
  call random_number(testvar4_ref_node)
  allocate(testvar5_ref_node(5,variable_length(5))) ! Mesh has five points
  call random_number(testvar5_ref_node)
  allocate(testvar6_ref_node(5,variable_length(6))) ! Mesh has five points
  call random_number(testvar6_ref_node)


  ! Data for first output file
  ! CELL DATA
  call inv_mesh%init_cell_data()

  ! Write whole variable at once
  call inv_mesh%add_cell_variable('cell_variable_1', &
                                  nentries = variable_length(1))
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_1',   &
                              values    = testvar1_ref)

  ! Write variable element-wise
  call inv_mesh%add_cell_variable('cell_variable_2', &
                                  nentries = variable_length(2))
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
                                  nentries = variable_length(3))
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar3_ref(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar3_ref(:,2:variable_length(3)),   &
                              ientry    = [2, variable_length(3)])

  ! END OF CELL DATA                          

  ! NODE DATA                          
  call inv_mesh%init_node_data()

  ! Write whole variable at once
  call inv_mesh%add_node_variable('node_variable_1', &
                                  nentries = variable_length(1))
  call inv_mesh%add_node_data(var_name  = 'node_variable_1',   &
                              values    = testvar1_ref_node)

  ! Write variable element-wise
  call inv_mesh%add_node_variable('node_variable_2', &
                                  nentries = variable_length(2))
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
                                  nentries = variable_length(3))
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar3_ref_node(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar3_ref_node(:,2:variable_length(3)),   &
                              ientry    = [2, variable_length(3)])

  ! END OF NODE DATA                           


  ! Write file to first disk
  call inv_mesh%dump_data_xdmf('./output/test_set_mixed_data_1st_file')

  call inv_mesh%free_node_and_cell_data()

  
  ! Data for second output file
  ! CELL DATA
  call inv_mesh%init_cell_data()

  ! Write whole variable at once
  call inv_mesh%add_cell_variable('cell_variable_1', &
                                  nentries = variable_length(4))
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_1',   &
                              values    = testvar4_ref)

  ! Write variable element-wise
  call inv_mesh%add_cell_variable('cell_variable_2', &
                                  nentries = variable_length(5))
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_2',          &
                              values    = testvar5_ref(1:1,:),   &
                              ielement  = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_2',          &
                              values    = testvar5_ref(2:2,:),   &
                              ielement  = [2, 2])
  
  ! Just to add confusion, init both once more
  call inv_mesh%init_cell_data()
  call inv_mesh%init_node_data()

  ! Write variable snap-wise
  call inv_mesh%add_cell_variable('cell_variable_3', &
                                  nentries = variable_length(6))
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar6_ref(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar6_ref(:,2:variable_length(6)),   &
                              ientry    = [2, variable_length(6)])

  ! END OF CELL DATA                          

  ! NODE DATA                          
  call inv_mesh%init_node_data()

  ! Write whole variable at once
  call inv_mesh%add_node_variable('node_variable_1', &
                                  nentries = variable_length(4))
  call inv_mesh%add_node_data(var_name  = 'node_variable_1',   &
                              values    = testvar4_ref_node)

  ! Write variable element-wise
  call inv_mesh%add_node_variable('node_variable_2', &
                                  nentries = variable_length(5))
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar5_ref_node(1:1,:),   &
                              ielement  = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar5_ref_node(2:2,:),   &
                              ielement  = [2, 2])
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar5_ref_node(3:5,:),   &
                              ielement  = [3, 5])

  ! Just to add confusion, init both once more
  call inv_mesh%init_cell_data()
  call inv_mesh%init_node_data()

  ! Write variable snap-wise
  call inv_mesh%add_node_variable('node_variable_3', &
                                  nentries = variable_length(6))
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar6_ref_node(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar6_ref_node(:,2:variable_length(6)),   &
                              ientry    = [2, variable_length(6)])

  ! END OF NODE DATA                           


  ! Write data to second file
  call inv_mesh%dump_data_xdmf('./output/test_set_mixed_data_2nd_file')

  call inv_mesh%free_node_and_cell_data()


  ! Check first output file
  call nc_open_for_read(filename = './output/test_set_mixed_data_1st_file.nc', ncid = ncid)

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


  ! Check second output file
  call nc_open_for_read(filename = './output/test_set_mixed_data_2nd_file.nc', ncid = ncid)

  ! CHECK FOR CELL DATA
  nf_status = nf90_inq_ncid(ncid=ncid, name='cell_data', grp_ncid=grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, 'Group for cell data from 2nd file has been created')

  ! Get variable 1
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'cell_variable_1', & 
                         values  = testvar4_res, &
                         limits  = [0.,1.] )

  do id = 1, size(testvar4_ref, 2) 
    call assert_comparable(testvar4_ref(:,id), testvar4_res(:,id),  &
                           1e-7, 'Variable 1 retrieved successfully from 2nd file')
  end do

  ! Get variable 2
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'cell_variable_2', & 
                         values  = testvar5_res, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar5_ref, 2) 
    call assert_comparable(testvar5_ref(:,id), testvar5_res(:,id),  &
                           1e-7, 'Variable 2 retrieved successfully from 2nd file')
  end do


  ! Get variable 3
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'cell_variable_3', & 
                         values  = testvar6_res, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar6_ref, 2) 
    call assert_comparable(testvar6_ref(:,id), testvar6_res(:,id),  &
                           1e-7, 'Variable 3 retrieved successfully from 2nd file')
  end do

  deallocate(testvar4_res)
  deallocate(testvar5_res)
  deallocate(testvar6_res)


  ! CHECK FOR NODE DATA
  nf_status = nf90_inq_ncid(ncid=ncid, name='node_data', grp_ncid=grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, 'Group for node data from 2nd file has been created')

  ! Get variable 1
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'node_variable_1', & 
                         values  = testvar4_res_node, &
                         limits  = [0.,1.] )

  do id = 1, size(testvar4_ref, 2) 
    call assert_comparable(testvar4_ref_node(:,id), testvar4_res_node(:,id),  &
                           1e-7, 'Variable 1 retrieved successfully from 2nd file')
  end do

  ! Get variable 2
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'node_variable_2', & 
                         values  = testvar5_res_node, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar5_ref, 2) 
    call assert_comparable(testvar5_ref_node(:,id), testvar5_res_node(:,id),  &
                           1e-7, 'Variable 2 retrieved successfully from 2nd file')
  end do


  ! Get variable 3
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'node_variable_3', & 
                         values  = testvar6_res_node, &
                         limits  = [0.,1.] )
  
  do id = 1, size(testvar6_ref, 2) 
    call assert_comparable(testvar6_ref_node(:,id), testvar6_res_node(:,id),  &
                           1e-7, 'Variable 3 retrieved successfully from 2nd file')
  end do

  nf_status = nf90_close(ncid = ncid)

end subroutine test_set_mixed_data_and_dump_into_two_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_set_node_time_data_and_dump
  use netcdf
  use nc_routines
  type(inversion_mesh_data_type)    :: inv_mesh
  character(len=NF90_MAX_NAME)      :: variable_names(3)
  integer                           :: variable_length(3), chunksizes(2), grp_ncid, id
  integer                           :: ncid, nf_status, varid, dimids(2)
  logical                           :: iscontiguous
  character(len=NF90_MAX_NAME)      :: dim_name

  ! Variables for first output file
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), testvar3_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_res(:,:), testvar2_res(:,:), testvar3_res(:,:)


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
                              'onvertices')

  variable_names  = ['variable_1', 'variable_2', 'variable_3']
  variable_length = [10, 10, 10]

  allocate(testvar1_ref(5,variable_length(1))) ! Mesh has five vertices
  call random_number(testvar1_ref)
  allocate(testvar2_ref(5,variable_length(2))) ! Mesh has five vertices
  call random_number(testvar2_ref)
  allocate(testvar3_ref(5,variable_length(3))) ! Mesh has five vertices
  call random_number(testvar3_ref)


  ! Data for first output file
  ! node DATA
  call inv_mesh%init_node_data(dt = 0.1d0, starttime = 0d0)

  ! Write whole variable at once
  call inv_mesh%add_node_variable('node_'//variable_names(1),    &
                                  nentries = variable_length(1), &
                                  istime   = .true.)
  call inv_mesh%add_node_data(var_name  = 'node_variable_1',   &
                              values    = testvar1_ref)

  ! Write variable element-wise
  call inv_mesh%add_node_variable('node_'//variable_names(2),    &
                                  nentries = variable_length(2), &
                                  istime   = .true.)
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar2_ref(1:1,:),   &
                              ielement  = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'node_variable_2',          &
                              values    = testvar2_ref(2:5,:),   &
                              ielement  = [2, 5])
  
  ! Just to add confusion, init both once more
  call inv_mesh%init_node_data()
  call inv_mesh%init_node_data()

  ! Write variable snap-wise
  call inv_mesh%add_node_variable('node_'//variable_names(3),    &
                                  nentries = variable_length(3), &
                                  istime   = .true.)
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar3_ref(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_node_data(var_name  = 'node_variable_3',          &
                              values    = testvar3_ref(:,2:variable_length(3)),   &
                              ientry    = [2, variable_length(3)])

  ! END OF node DATA                          

  ! Write file to disk
  call inv_mesh%dump_data_xdmf('./output/test_set_node_data_time')
  call inv_mesh%free_node_and_cell_data()

  call nc_open_for_read(filename = './output/test_set_node_data_time.nc', ncid = ncid)
  call check(nf90_inq_ncid(ncid=ncid, name='node_data', grp_ncid=grp_ncid))
  call assert_equal(nf_status, NF90_NOERR, 'Group for node data has been created. Error: '//&
                    NF90_STRERROR(nf_status))

  ! Get variable 1
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'node_'//variable_names(1), &
                         values  = testvar1_res, &
                         limits  = [0.,1.],      &
                         varid   = varid)

  call assert_comparable(testvar1_ref, testvar1_res,  &
                         1e-7, 'Variable 1 retrieved successfully')

  call check(nf90_inquire_variable(grp_ncid, varid,           &
                                   contiguous = iscontiguous, &
                                   dimids     = dimids,       &
                                   chunksizes = chunksizes))

  call check(nf90_inquire_dimension(grp_ncid, dimids(2),      &
                                    name      = dim_name))

  call assert_true((trim(dim_name).eq.'time'), &
    'time dimension is called '''//trim(dim_name)//''' instead of ''time''')
  call assert_true(.not.iscontiguous, 'Time variable should be chunked')
  call assert_equal(chunksizes,  [inv_mesh%get_nvertices(), 1], &
                    'Chunk size in time dimension is one')

  ! Get variable 2
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'node_'//variable_names(2), &
                         values  = testvar2_res, &
                         limits  = [0.,1.],      &
                         varid   = varid)
  
  call assert_comparable(testvar2_ref, testvar2_res,  &
                         1e-7, 'Variable 2 retrieved successfully')

  call check(nf90_inquire_variable(grp_ncid, varid,           &
                                   contiguous = iscontiguous, &
                                   chunksizes = chunksizes,   &
                                   dimids     = dimids))

  call check(nf90_inquire_dimension(grp_ncid, dimids(2),      &
                                    name      = dim_name))

  call assert_true((trim(dim_name).eq.'time'), &
    'time dimension is called '''//trim(dim_name)//''' instead of ''time''')
  call assert_true(.not.iscontiguous, 'Time variable should be chunked')
  call assert_equal(chunksizes,  [inv_mesh%get_nvertices(), 1], &
                    'Chunk size in time dimension is one')

  ! Get variable 3
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'node_'//variable_names(3), &
                         values  = testvar3_res, &
                         limits  = [0.,1.],      &
                         varid   = varid)
  
  call assert_comparable(testvar3_ref, testvar3_res,  &
                         1e-7, 'Variable 3 retrieved successfully')

  call check(nf90_inquire_variable(grp_ncid, varid,           &
                                   contiguous = iscontiguous, &
                                   chunksizes = chunksizes))

  call check(nf90_inquire_dimension(grp_ncid, dimids(2),      &
                                    name      = dim_name))

  call assert_true((trim(dim_name).eq.'time'), &
    'time dimension is called '''//trim(dim_name)//''' instead of ''time''')
  call assert_true(.not.iscontiguous, 'Time variable should be chunked')
  call assert_equal(chunksizes,  [inv_mesh%get_nvertices(), 1], &
                    'Chunk size in time dimension is one')
  call check(nf90_close(ncid = ncid))


end subroutine test_set_node_time_data_and_dump
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_set_cell_time_data_and_dump
  use netcdf
  use nc_routines
  type(inversion_mesh_data_type)    :: inv_mesh
  character(len=NF90_MAX_NAME)      :: variable_names(3)
  integer                           :: variable_length(3), chunksizes(2), grp_ncid, id
  integer                           :: ncid, nf_status, varid, dimids(2)
  logical                           :: iscontiguous
  character(len=NF90_MAX_NAME)      :: dim_name

  ! Variables for first output file
  real(kind=sp), allocatable        :: testvar1_ref(:,:), testvar2_ref(:,:), testvar3_ref(:,:)
  real(kind=sp), allocatable        :: testvar1_res(:,:), testvar2_res(:,:), testvar3_res(:,:)


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
                              'onvertices')

  variable_names  = ['variable_1', 'variable_2', 'variable_3']
  variable_length = [10, 10, 10]

  allocate(testvar1_ref(2,variable_length(1))) ! Mesh has two cells
  call random_number(testvar1_ref)
  allocate(testvar2_ref(2,variable_length(2))) ! Mesh has two cells
  call random_number(testvar2_ref)
  allocate(testvar3_ref(2,variable_length(3))) ! Mesh has two cells
  call random_number(testvar3_ref)


  ! Data for first output file
  ! CELL DATA
  call inv_mesh%init_cell_data(dt = 0.1d0, starttime = 0d0)

  ! Write whole variable at once
  call inv_mesh%add_cell_variable('cell_'//variable_names(1),    &
                                  nentries = variable_length(1), &
                                  istime   = .true.)
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_1',   &
                              values    = testvar1_ref)

  ! Write variable element-wise
  call inv_mesh%add_cell_variable('cell_'//variable_names(2),    &
                                  nentries = variable_length(2), &
                                  istime   = .true.)
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
  call inv_mesh%add_cell_variable('cell_'//variable_names(3),    &
                                  nentries = variable_length(3), &
                                  istime   = .true.)
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar3_ref(:,1:1),   &
                              ientry    = [1, 1])
  call inv_mesh%add_cell_data(var_name  = 'cell_variable_3',          &
                              values    = testvar3_ref(:,2:variable_length(3)),   &
                              ientry    = [2, variable_length(3)])

  ! END OF CELL DATA                          

  ! Write file to disk
  call inv_mesh%dump_data_xdmf('./output/test_set_cell_data_time')
  call inv_mesh%free_node_and_cell_data()

  call nc_open_for_read(filename = './output/test_set_cell_data_time.nc', ncid = ncid)
  nf_status = nf90_inq_ncid(ncid=ncid, name='cell_data', grp_ncid=grp_ncid)
  call assert_equal(nf_status, NF90_NOERR, 'Group for cell data has been created')

  ! Get variable 1
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'cell_'//variable_names(1), &
                         values  = testvar1_res, &
                         limits  = [0.,1.],      &
                         varid   = varid)

  call assert_comparable(testvar1_ref(:,:), testvar1_res(:,:),  &
                         1e-7, 'Variable 1 retrieved successfully')

  nf_status = nf90_inquire_variable(grp_ncid, varid,           &
                                    contiguous = iscontiguous, &
                                    dimids     = dimids,       &
                                    chunksizes = chunksizes)

  nf_status = nf90_inquire_dimension(grp_ncid, dimids(2),      &
                                     name      = dim_name)

  call assert_true((trim(dim_name).eq.'time'), &
    'time dimension is called '''//trim(dim_name)//''' instead of ''time''')
  call assert_true(.not.iscontiguous, 'Time variable should be chunked')
  call assert_equal(chunksizes,  [inv_mesh%get_nelements(), 1], &
                    'Chunk size in time dimension is one')

  ! Get variable 2
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'cell_'//variable_names(2), &
                         values  = testvar2_res, &
                         limits  = [0.,1.],      &
                         varid   = varid)
  
  call assert_comparable(testvar2_ref(:,:), testvar2_res(:,:),  &
                         1e-7, 'Variable 2 retrieved successfully')

  nf_status = nf90_inquire_variable(grp_ncid, varid,           &
                                    contiguous = iscontiguous, &
                                    chunksizes = chunksizes,   &
                                    dimids     = dimids)

  nf_status = nf90_inquire_dimension(grp_ncid, dimids(2),      &
                                     name      = dim_name)

  call assert_true((trim(dim_name).eq.'time'), &
    'time dimension is called '''//trim(dim_name)//''' instead of ''time''')
  call assert_true(.not.iscontiguous, 'Time variable should be chunked')
  call assert_equal(chunksizes,  [inv_mesh%get_nelements(), 1], &
                    'Chunk size in time dimension is one')

  ! Get variable 3
  call nc_getvar_by_name(ncid    = grp_ncid,     &
                         varname = 'cell_'//variable_names(3), &
                         values  = testvar3_res, &
                         limits  = [0.,1.],      &
                         varid   = varid)
  
  call assert_comparable(testvar3_ref(:,:), testvar3_res(:,:),  &
                         1e-7, 'Variable 3 retrieved successfully')

  nf_status = nf90_inquire_variable(grp_ncid, varid,           &
                                    contiguous = iscontiguous, &
                                    chunksizes = chunksizes)

  nf_status = nf90_inquire_dimension(grp_ncid, dimids(2),      &
                                     name      = dim_name)

  call assert_true((trim(dim_name).eq.'time'), &
    'time dimension is called '''//trim(dim_name)//''' instead of ''time''')
  call assert_true(.not.iscontiguous, 'Time variable should be chunked')
  call assert_equal(chunksizes,  [inv_mesh%get_nelements(), 1], &
                    'Chunk size in time dimension is one')
  nf_status = nf90_close(ncid = ncid)


end subroutine test_set_cell_time_data_and_dump
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_valence
  type(inversion_mesh_type) :: inv_mesh

  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
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

  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
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
subroutine test_point_in_element_triangle_mesh
  type(inversion_mesh_type)  :: inv_mesh
  integer, parameter         :: ndim = 3, nvertices = 3, nstep = 10
  real(kind=dp)              :: vertices(ndim, 0:nvertices-1)
  integer                    :: connectivity(3, 1), itest
  integer, parameter         :: ntest = (nstep+3)**2, nbasisfuncs_per_elem = 1
  real(kind=dp)              :: testpoint_in(ndim, ntest)
  logical                    :: is_in_element_one, is_in_element_npoint_ref_temp(ntest)
  logical, allocatable       :: is_in_element_npoint(:), is_in_element_npoint_ref(:)
  integer                    :: ix, iy, ntest_tot

  ! Simple vertex
  vertices(:, 0) = [0, 0, 0]
  vertices(:, 1) = [1, 0, 0]
  vertices(:, 2) = [1, 1, 0]

  itest = 0

  ! Create along each base vector of the tetrahedron points from [-0.1 to 1.1]
  do ix = -1, nstep + 1
    do iy = -1, nstep - ix + 1
      itest = itest + 1
      testpoint_in(:, itest) =   vertices(:,0)                                        & 
                               + (real(ix, kind=dp) / nstep) * (vertices(:,1) - vertices(:,0))  &
                               + (real(iy, kind=dp) / nstep) * (vertices(:,2) - vertices(:,0))

      if (any([ix, iy] == -1)) then
        is_in_element_npoint_ref_temp(itest) = .false.
      elseif (ix+iy>nstep) then
        is_in_element_npoint_ref_temp(itest) = .false.
      else
        is_in_element_npoint_ref_temp(itest) = .true.
      end if

    end do
  end do

  ntest_tot = itest
  allocate(is_in_element_npoint_ref(ntest_tot))
  is_in_element_npoint_ref = is_in_element_npoint_ref_temp(1:ntest_tot)

  connectivity(:,1) = [1, 2, 3]
  call inv_mesh%initialize_mesh(1, vertices, connectivity, nbasisfuncs_per_elem)

  do itest = 1, ntest_tot
    is_in_element_one = inv_mesh%point_in_element(1, testpoint_in(:, itest))
    call assert_true((is_in_element_one .eqv. is_in_element_npoint_ref(itest)), &
                     'Points are correctly found to be in element - 1 point')
    if (is_in_element_one .neqv. is_in_element_npoint_ref(itest)) then
      print *, 'Is inside ', is_in_element_one, ', should be: ', is_in_element_npoint_ref(itest)
      print *, 'P: ', testpoint_in(:, itest)
    end if
  end do

  allocate(is_in_element_npoint(ntest_tot))
  is_in_element_npoint = inv_mesh%point_in_element(1, testpoint_in(:, 1:ntest_tot))
  call assert_true((is_in_element_npoint .eqv. is_in_element_npoint_ref), &
                   'Points are correctly found to be in element - n points')
  call inv_mesh%freeme()


  vertices(:, 0) = [2, 1, 0]
  vertices(:, 1) = [4,-3, 1]
  vertices(:, 2) = [2, 2, 2]

  itest = 0

  ! Create along each base vector of the tetrahedron points from [-0.1 to 1.1]
  do ix = -1, nstep + 1
    do iy = -1, nstep - ix + 1
      itest = itest + 1
      testpoint_in(:, itest) =   vertices(:,0)                                        & 
                               + (real(ix, kind=dp) / nstep) * (vertices(:,1) - vertices(:,0))  &
                               + (real(iy, kind=dp) / nstep) * (vertices(:,2) - vertices(:,0))

      if (any([ix, iy] == -1)) then
        is_in_element_npoint_ref_temp(itest) = .false.
      elseif (ix+iy>nstep) then
        is_in_element_npoint_ref_temp(itest) = .false.
      else
        is_in_element_npoint_ref_temp(itest) = .true.
      end if

    end do
  end do

  ntest_tot = itest
  is_in_element_npoint_ref = is_in_element_npoint_ref_temp(1:ntest_tot)

  connectivity(:,1) = [1, 2, 3]
  call inv_mesh%initialize_mesh(1, vertices, connectivity, nbasisfuncs_per_elem)

  do itest = 1, ntest_tot
    is_in_element_one = inv_mesh%point_in_element(1, testpoint_in(:, itest))
    call assert_true((is_in_element_one .eqv. is_in_element_npoint_ref(itest)), &
                     'Points are correctly found to be in element - 1 point')
    if (is_in_element_one .neqv. is_in_element_npoint_ref(itest)) then
      print *, 'Is inside ', is_in_element_one, ', should be: ', is_in_element_npoint_ref(itest)
      print *, 'P: ', testpoint_in(:, itest)
    end if
  end do

  is_in_element_npoint = inv_mesh%point_in_element(1, testpoint_in(:, 1:ntest_tot))
  call assert_true((is_in_element_npoint .eqv. is_in_element_npoint_ref), &
                   'Points are correctly found to be in element - n points')
  call inv_mesh%freeme()

end subroutine test_point_in_element_triangle_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_point_in_element_tetrahedral_mesh
  type(inversion_mesh_type)  :: inv_mesh
  integer, parameter         :: ndim = 3, nvertices = 4
  real(kind=dp)              :: vertices(ndim, 0:nvertices-1)
  integer                    :: connectivity(4, 1), itest
  integer, parameter         :: nbasisfuncs_per_elem = 1, nstep = 10
  integer, parameter         :: ntest = (nstep+3)**3 / 3
  real(kind=dp)              :: testpoint_in(ndim, ntest)
  logical                    :: is_in_element_one, is_in_element_npoint_ref_temp(ntest)
  logical, allocatable       :: is_in_element_npoint(:), is_in_element_npoint_ref(:)
  integer                    :: ix, iy, iz, ntest_tot

  ! Try simplistic tetrahedron
  vertices(:, 0) = [0, 0, 0]
  vertices(:, 1) = [1, 0, 0]
  vertices(:, 2) = [0, 1, 0]
  vertices(:, 3) = [0, 0, 1]

  itest = 0

  ! Create along each base vector of the tetrahedron points from [-0.1 to 1.1]
  do ix = -1, nstep + 1
    do iy = -1, nstep - ix + 1
      do iz = -1, nstep - ix - iy + 1
        itest = itest + 1
        testpoint_in(:, itest) =   vertices(:,0)                                        & 
                                 + (real(ix, kind=dp) / nstep) * (vertices(:,1) - vertices(:,0))  &
                                 + (real(iy, kind=dp) / nstep) * (vertices(:,2) - vertices(:,0))  &
                                 + (real(iz, kind=dp) / nstep) * (vertices(:,3) - vertices(:,0))

        if (any([ix, iy, iz] == -1)) then
          is_in_element_npoint_ref_temp(itest) = .false.
        elseif (ix+iy+iz>nstep) then
          is_in_element_npoint_ref_temp(itest) = .false.
        else
          is_in_element_npoint_ref_temp(itest) = .true.
        end if

      end do
    end do
  end do

  ntest_tot = itest
  allocate(is_in_element_npoint_ref(ntest_tot))
  is_in_element_npoint_ref = is_in_element_npoint_ref_temp(1:ntest_tot)

  connectivity(:,1) = [1, 2, 3, 4]
  call inv_mesh%initialize_mesh(4, vertices, connectivity, nbasisfuncs_per_elem)

  do itest = 1, ntest_tot
    is_in_element_one = inv_mesh%point_in_element(1, testpoint_in(:, itest))
    call assert_true((is_in_element_one .eqv. is_in_element_npoint_ref(itest)), &
                     'Points are correctly found to be in element - 1 point')
  end do

  allocate(is_in_element_npoint(ntest_tot))
  is_in_element_npoint = inv_mesh%point_in_element(1, testpoint_in(:, 1:ntest_tot))
  call assert_true((is_in_element_npoint .eqv. is_in_element_npoint_ref), &
                   'Points are correctly found to be in element - n points')
  call inv_mesh%freeme()

  ! Try rotated tetrahedron at arbitrary location
  vertices(:, 0) = [2, 1, 0]
  vertices(:, 1) = [4,-3, 1]
  vertices(:, 2) = [2, 2, 2]
  vertices(:, 3) = [2,-8, 1]

  itest = 0
  do ix = 1, nstep
    do iy = 1, nstep - ix
      do iz = 1, nstep - ix - iy
        itest = itest + 1
        testpoint_in(:, itest) =   vertices(:,0)                                        & 
                                 + (real(ix, kind=dp) / nstep) * (vertices(:,1) - vertices(:,0))  &
                                 + (real(iy, kind=dp) / nstep) * (vertices(:,2) - vertices(:,0))  &
                                 + (real(iz, kind=dp) / nstep) * (vertices(:,3) - vertices(:,0))

      end do
    end do
  end do

  ntest_tot = itest

  connectivity(:,1) = [1, 2, 3, 4]
  call inv_mesh%initialize_mesh(4, vertices, connectivity, nbasisfuncs_per_elem)

  do itest = 1, ntest_tot

    is_in_element_one = inv_mesh%point_in_element(1, testpoint_in(:, itest))
    call assert_true(is_in_element_one, 'Points are correctly found to be in element - 1 point')

  end do

  is_in_element_npoint = inv_mesh%point_in_element(1, testpoint_in(:, 1:ntest_tot))
  call assert_true(is_in_element_npoint, 'Points are correctly found to be in element - n points')
  call inv_mesh%freeme()

end subroutine test_point_in_element_tetrahedral_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_random_points_triangle_mesh
  use tetrahedra, only        : point_in_triangle_3d
  use halton_sequence, only   : free_halton
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
      
      ! Using pseudorandom numbers 
      points = inv_mesh%generate_random_points(1, npoints, .false.)
      isintriangle = point_in_triangle_3d(vertices, points, isinplane)

      call assert_true(isinplane,    'Points are in plane of triangle')
      call assert_true(isintriangle, 'Points are in triangle')

      ! Using quasirandom numbers 
      call free_halton()
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

!-----------------------------------------------------------------------------------------
subroutine test_weight
  use halton_sequence, only          : free_halton
  type(inversion_mesh_data_type)    :: inv_mesh
  integer, parameter                :: nrandom = 10
  real(kind=dp)                     :: points(3,4), weight(4), weight_ref(4)
  real(kind=dp)                     :: random_points(3,nrandom), weight_random(nrandom)
  integer                           :: ielement, ivertex
 


  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
                              'onvertices')

  ! Weight should be one at each of the vertices for the respective parameter and 
  ! zero for all the others
  do ielement = 1, inv_mesh%get_nelements()
    points = inv_mesh%get_element(ielement)
    do ivertex = 1, 4 ! Tetrahedral mesh
      weight = inv_mesh%weights(ielement, ivertex, points)
      weight_ref = 0
      weight_ref(ivertex) = 1
      call assert_comparable(weight, weight_ref, 1d-10, 'Weight at vertex is one')
    end do
  end do

  ! Weight should be between one and zero everywhere (inside the element)
  do ielement = 1, inv_mesh%get_nelements()
    random_points = inv_mesh%generate_random_points( ielement, nrandom, .false. )
    do ivertex = 1, 4 ! Tetrahedral mesh
      weight_random = inv_mesh%weights(ielement, ivertex, random_points)
      call assert_true(all(weight_random<1.d0), 'Weight is below one everywhere')
      call assert_true(all(weight_random>0.d0), 'Weight is above zero everywhere')
    end do
  end do

  ! Test with quasi-random numbers
  
  ! Reset Halton number generator (may have been used in earlier tests)
  call free_halton()

  do ielement = 1, inv_mesh%get_nelements()
    random_points = inv_mesh%generate_random_points( ielement, nrandom, .true. )
    do ivertex = 1, 4 ! Tetrahedral mesh
      weight_random = inv_mesh%weights(ielement, ivertex, random_points)
      call assert_true(all(weight_random<1.d0), 'Weight is below one everywhere')
      call assert_true(all(weight_random>0.d0), 'Weight is above zero everywhere')
    end do
  end do

end subroutine test_weight
!-----------------------------------------------------------------------------------------

!!-----------------------------------------------------------------------------------------
!! Test the weight<0 issue occuring on KH's huge mesh
!! Tests the specific element where it was occuring
!subroutine test_weight_issue30
!  use halton_sequence, only          : free_halton
!  type(inversion_mesh_data_type)    :: inv_mesh
!  integer, parameter                :: nrandom = 100000
!  real(kind=dp)                     :: points(3,4), weight(4), weight_ref(4)
!  real(kind=dp)                     :: random_points(3,nrandom), weight_random(nrandom)
!  real(kind=dp)                     :: random_point(3,1), weight_random_single(1)
!  integer                           :: ielement, ivertex, ivertex_test
! 
!
!
!  call inv_mesh%read_tet_mesh('unit_tests/vertices.issue30', &
!                              'unit_tests/facets.issue30', &
!                              'onvertices')
!
!  ! Weight should be one at each of the vertices for the respective parameter and 
!  ! zero for all the others
!  do ielement = 1, inv_mesh%get_nelements()
!    points = inv_mesh%get_element(ielement)
!    do ivertex = 1, 4 ! Tetrahedral mesh
!      weight = inv_mesh%weights(ielement, ivertex, points)
!      weight_ref = 0
!      weight_ref(ivertex) = 1
!      call assert_comparable(weight, weight_ref, 1d-10, 'Weight at vertex is one')
!    end do
!  end do
!
!  ! Weight should be between one and zero everywhere (inside the element)
!  do ielement = 1, inv_mesh%get_nelements()
!    random_points = inv_mesh%generate_random_points( ielement, nrandom, .false. )
!    do ivertex = 1, 4 ! Tetrahedral mesh
!      weight_random = inv_mesh%weights(ielement, ivertex, random_points)
!      call assert_true(all(weight_random<1.d0), 'Weight is below one everywhere')
!      call assert_true(all(weight_random>0.d0), 'Weight is above zero everywhere')
!    end do
!  end do
!
!  ! Weight should be between one and zero everywhere (inside the element)
!  ! Test with quasi-random numbers
!  
!  ! Reset Halton number generator (may have been used in earlier tests)
!  call free_halton()
!
!  do ielement = 1, inv_mesh%get_nelements()
!    random_points = inv_mesh%generate_random_points( ielement, nrandom, .true. )
!    do ivertex = 1, 4 ! Tetrahedral mesh
!      weight_random = inv_mesh%weights(ielement, ivertex, random_points)
!      call assert_true(all(weight_random<1.d0), 'Weight is below one everywhere')
!      call assert_true(all(weight_random>0.d0), 'Weight is above zero everywhere')
!    end do
!  end do
!
!end subroutine test_weight_issue30
!!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Not called by default
! Takes several hours to finish and tests mainly the stability of the Halton sequencies 
! for more than 1e9 random numbers (int4 overflow)
subroutine test_weight_large
  use halton_sequence, only          : free_halton
  type(inversion_mesh_data_type)    :: inv_mesh
  integer, parameter                :: nrandom = 1000
  real(kind=dp)                     :: points(3,4), weight(4), weight_ref(4)
  real(kind=dp)                     :: random_points(3,nrandom), weight_random(nrandom)
  integer                           :: ielement, ivertex
 


  call inv_mesh%read_tet_mesh('Meshes/vertices.fine_20_96', &
                              'Meshes/facets.fine_20_96',   &
                              'onvertices')

  ! Weight should be one at each of the vertices for the respective parameter and 
  ! zero for all the others
  do ielement = 1, inv_mesh%get_nelements()
    points = inv_mesh%get_element(ielement)
    do ivertex = 1, 4 ! Tetrahedral mesh
      weight = inv_mesh%weights(ielement, ivertex, points)
      weight_ref = 0
      weight_ref(ivertex) = 1
      call assert_comparable(1+weight, 1+weight_ref, 1d-10, 'Weight at vertex is one')
    end do
  end do

  ! Weight should be between one and zero everywhere (inside the element)
  do ielement = 1, inv_mesh%get_nelements()
    random_points = inv_mesh%generate_random_points( ielement, nrandom, .false. )
    do ivertex = 1, 4 ! Tetrahedral mesh
      weight_random = inv_mesh%weights(ielement, ivertex, random_points)
      call assert_true(all(weight_random<1.d0), 'Weight is below one everywhere')
      call assert_true(all(weight_random>0.d0), 'Weight is above zero everywhere')
    end do
  end do

  ! Reset Halton number generator (may have been used in earlier tests)
  call free_halton()

  do ielement = 1, inv_mesh%get_nelements()
    random_points = inv_mesh%generate_random_points( ielement, nrandom, .true. )
    do ivertex = 1, 4 ! Tetrahedral mesh
      weight_random = inv_mesh%weights(ielement, ivertex, random_points)
      call assert_true(all(weight_random<1.d0), 'Weight is below one everywhere')
      call assert_true(all(weight_random>0.d0), 'Weight is above zero everywhere')
    end do
  end do

end subroutine test_weight_large
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_integration_in_tetrahedron
  use netcdf
  use montecarlo,                  only: integrated_type, allallconverged, allisconverged
  type(integrated_type)             :: int_test
  type(inversion_mesh_data_type)    :: inv_mesh
  integer, parameter                :: nptperstep = 10000
  real(kind=dp)                     :: values(nptperstep,1), coords(3, nptperstep), integral(1)
  real(kind=dp), parameter          :: pi = 3.141419265d0

  call inv_mesh%read_tet_mesh('./vertices.TEST', './facets.TEST', &
                              'onvertices')


  ! First integrate over constant value in tetrahedron - should result in volume
  call int_test%initialize_montecarlo(nfuncs = 1,                        & 
                                      volume = inv_mesh%get_volume(1),   & 
                                      allowed_error = 1d-2,              &
                                      allowed_relerror = 1d-2)
  
  values = 1
  call int_test%check_montecarlo_integral(values)

  call assert_comparable(int_test%getintegral(), [1.d9/6.d0], 1d-10, &
                         'Integration over constant function gives volume')
  call int_test%freeme()



  ! Integrate over sphere inscribed into tetrahedron
  call int_test%initialize_montecarlo(nfuncs = 1,                        & 
                                      volume = inv_mesh%get_volume(1),   & 
                                      allowed_error = 1d-3,              &
                                      allowed_relerror = 1d-3)
  do while (.not.int_test%areallconverged()) ! Beginning of Monte Carlo loop
      coords = inv_mesh%generate_random_points(1, nptperstep, .false.)
      values = 0

      ! Sphere with radius sqrt(1e3/3), which touches tetrahedral large plane
      where( sum(coords**2, dim=1) < (1d6/3.) ) values(:,1) = 1

      call int_test%check_montecarlo_integral(values)
  end do
 

  integral = int_test%getintegral()
  call assert_comparable(integral, [pi*(1./3.)**(1.5) * 1d9/6.], &
                         1d-3, 'Integral == 1.10076')

end subroutine test_integration_in_tetrahedron
!-----------------------------------------------------------------------------------------

end module test_inversion_mesh
!=========================================================================================
