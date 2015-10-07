module test_nc_routines
  use global_parameters
  use netcdf
  use nc_routines
  use ftnunit

  implicit none

  character(len=128), save         :: nc_filename='./output/test_nc_routines.nc'
  integer                          :: dim_size(3) 
  real(kind=sp), allocatable, save :: testdata_fp(:,:,:)
  integer      , allocatable, save :: testdata_int(:,:,:)

contains

!-----------------------------------------------------------------------------------------
subroutine test_nc_create_testfile
  integer       :: ncid, nmode, idim_1, idim_2, idim_3
  integer       :: nc_dimid(3), nc_varid(6)
  
  
  dim_size = [10, 20, 30]

  allocate(testdata_fp(  dim_size(1), dim_size(2), dim_size(3) ))
  allocate(testdata_int( dim_size(1), dim_size(2), dim_size(3) ))

  do idim_1 = 1, dim_size(1)
    do idim_2 = 1, dim_size(2)
      do idim_3 = 1, dim_size(3)
        testdata_fp(idim_1, idim_2, idim_3) =   idim_1                       &
                                              + dim_size(1) * idim_2         &
                                              + product(dim_size(1:2)) * idim_3
        testdata_int(idim_1, idim_2, idim_3) = - idim_1                      &
                                               - dim_size(1) * idim_2        &
                                               - product(dim_size(1:2)) * idim_3
      end do
    end do
  end do

  nmode = ior(NF90_CLOBBER, NF90_NETCDF4)
  call check(nf90_create(path = nc_filename, cmode = nmode, ncid = ncid))

  call check(nf90_def_dim(ncid, 'dimension_1', dim_size(1), nc_dimid(1)))
  call check(nf90_def_dim(ncid, 'dimension_2', dim_size(2), nc_dimid(2)))
  call check(nf90_def_dim(ncid, 'dimension_3', dim_size(3), nc_dimid(3)))

  ! 1D Float
  call check(nf90_def_var(ncid = ncid, name = 'var_1d_float', xtype = NF90_FLOAT, &
                          dimids = [nc_dimid(1)], varid = nc_varid(1)))

  ! 2D Float
  call check(nf90_def_var(ncid = ncid, name = 'var_2d_float', xtype = NF90_FLOAT, &
                          dimids = [nc_dimid(1:2)], varid = nc_varid(2)))

  ! 3D Float
  call check(nf90_def_var(ncid = ncid, name = 'var_3d_float', xtype = NF90_FLOAT, &
                          dimids = [nc_dimid(1:3)], varid = nc_varid(3)))

  ! 1D Integer
  call check(nf90_def_var(ncid = ncid, name = 'var_1d_int', xtype = NF90_INT, &
                          dimids = [nc_dimid(1)], varid = nc_varid(4)))

  ! 2D Integer
  call check(nf90_def_var(ncid = ncid, name = 'var_2d_int', xtype = NF90_INT, &
                          dimids = [nc_dimid(1:2)], varid = nc_varid(5)))

  ! 3D Integer
  call check(nf90_def_var(ncid = ncid, name = 'var_3d_int', xtype = NF90_INT, &
                          dimids = [nc_dimid(1:3)], varid = nc_varid(6)))


  call check(nf90_enddef(ncid = ncid))

  call check(nf90_put_var(ncid = ncid, varid = nc_varid(1), values = testdata_fp(:,1,1)))
  call check(nf90_put_var(ncid = ncid, varid = nc_varid(2), values = testdata_fp(:,:,1)))
  call check(nf90_put_var(ncid = ncid, varid = nc_varid(3), values = testdata_fp(:,:,:)))

  call check(nf90_put_var(ncid = ncid, varid = nc_varid(4), values = testdata_int(:,1,1)))
  call check(nf90_put_var(ncid = ncid, varid = nc_varid(5), values = testdata_int(:,:,1)))
  call check(nf90_put_var(ncid = ncid, varid = nc_varid(6), values = testdata_int(:,:,:)))

  call check(nf90_close(ncid = ncid))

end subroutine test_nc_create_testfile
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_create_file
  integer           :: ncid, formatnum, status
  character(len=80) :: filename = './output/test_nc_create_file.nc'

  call nc_create_file(filename = filename, ncid = ncid)

  call assert_file_exists(filename, 'Create netcdf file')

  status = nf90_inquire(ncid = ncid, formatNum = formatnum)
  call assert_equal(status, NF90_NOERR, 'Inquire successful: '//nf90_strerror(status))
  call assert_equal(formatnum, NF90_FORMAT_NETCDF4, 'Opened file has NetCDF4 format')

  call nc_close_file(ncid)

  status = nf90_inquire(ncid = ncid, formatNum = formatnum)
  call assert_equal(status, NF90_EBADID, 'File closed again')

end subroutine test_nc_create_file
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_create_group
  integer           :: ncid, formatnum, status, ncid_group, ncid_group_2
  character(len=80) :: filename = './output/test_nc_create_group.nc'

  call nc_create_file(filename = filename, ncid = ncid)
  
  ! Try to create once
  call nc_create_group(ncid = ncid, group_name = 'test', ncid_group = ncid_group)

  ! Try to create same group once more, should work since overwrite is set
  call nc_create_group(ncid = ncid, group_name = 'test', &
                       ncid_group = ncid_group_2, overwrite = .true.)
  call assert_equal(ncid_group, ncid_group_2, 'Second call to nc_create_group results in same ncid')

  call assert_file_exists(filename, 'Create netcdf file')

  status = nf90_inquire(ncid = ncid, formatNum = formatnum)
  call assert_equal(status, NF90_NOERR, 'Inquire successful: '//nf90_strerror(status))
  call assert_equal(formatnum, NF90_FORMAT_NETCDF4, 'Opened file has NetCDF4 format')
  
  status = nf90_inq_ncid(ncid = ncid, name = 'test', grp_ncid = ncid_group)
  call assert_equal(status, NF90_NOERR, 'Inquire group successful: '//nf90_strerror(status))

  call nc_close_file(ncid)

  status = nf90_inquire(ncid = ncid, formatNum = formatnum)
  call assert_equal(status, NF90_EBADID, 'File closed again')

end subroutine test_nc_create_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_open_for_read
  integer           :: ncid, formatnum, status
  character(len=80) :: filename = './output/test_nc_open_for_write'

  ! Create test file
  call nc_create_file(filename = filename, ncid = ncid)
  call nc_close_file(ncid)

  ! Try to open
  call nc_open_for_read(filename = filename, ncid = ncid)

  status = nf90_inquire(ncid = ncid, formatNum = formatnum)

  call assert_equal(status, NF90_NOERR, 'Inquire successful: '//nf90_strerror(status))
  call assert_equal(formatnum, NF90_FORMAT_NETCDF4, 'Opened file has NetCDF4 format')

  status = nf90_redef(ncid)
  call assert_equal(status, -37, 'File opened for read only')

  call check(nf90_close(ncid))

end subroutine test_nc_open_for_read
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_open_for_write
  integer           :: ncid, formatnum, status
  character(len=80) :: filename = './output/test_nc_open_for_write'

  ! Create test file
  call nc_create_file(filename = filename, ncid = ncid)
  call nc_close_file(ncid)

  ! Try to open
  call nc_open_for_write(filename = filename, ncid = ncid)

  status = nf90_inquire(ncid = ncid, formatNum = formatnum)

  call assert_equal(status, NF90_NOERR, 'Inquire successful: '//nf90_strerror(status))
  call assert_equal(formatnum, NF90_FORMAT_NETCDF4, 'Opened file has NetCDF4 format')

  status = nf90_redef(ncid)
  call assert_equal(status, NF90_NOERR, 'File opened for write')

  call check(nf90_close(ncid))

end subroutine test_nc_open_for_write
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_getvar_1d_float
  integer                     :: ncid
  real(kind=sp), allocatable  :: readdata_fp(:)

  call nc_open_for_read(filename = nc_filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_1d_float', values = readdata_fp)

  call assert_comparable(readdata_fp, testdata_fp(:,1,1), 1e-8, 'Retrieve 1D Float')

  call check(nf90_close(ncid))

end subroutine test_nc_getvar_1d_float
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_getvar_2d_float
  integer                     :: ncid, idim_2
  real(kind=sp), allocatable  :: readdata_fp(:,:)

  call nc_open_for_read(filename = nc_filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_2d_float', values = readdata_fp)

  do idim_2 = 1, size(readdata_fp, 2)
    call assert_comparable(readdata_fp(:, idim_2), &
                           testdata_fp(:, idim_2, 1), &
                           1e-8, 'Retrieve 3D Float')
  end do

  call check(nf90_close(ncid))

end subroutine test_nc_getvar_2d_float
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_getvar_3d_float
  integer                     :: ncid, idim_2, idim_3
  real(kind=sp), allocatable  :: readdata_fp(:,:,:)

  call nc_open_for_read(filename = nc_filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_float', values = readdata_fp)

  do idim_2 = 1, size(readdata_fp, 2)
    do idim_3 = 1, size(readdata_fp, 3)
      call assert_comparable(readdata_fp(:, idim_2, idim_3), &
                             testdata_fp(:, idim_2, idim_3), &
                             1e-8, 'Retrieve 3D Float')
    end do
  end do

  call check(nf90_close(ncid))

end subroutine test_nc_getvar_3d_float
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_getvar_1d_int
  integer                     :: ncid
  integer, allocatable        :: readdata_int(:)

  call nc_open_for_read(filename = nc_filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_1d_int', values = readdata_int)

  call assert_equal(readdata_int, testdata_int(:,1,1), 'Retrieve 1D int')

  call check(nf90_close(ncid))

end subroutine test_nc_getvar_1d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_getvar_2d_int
  integer                     :: ncid, idim_2
  integer, allocatable        :: readdata_int(:,:)

  call nc_open_for_read(filename = nc_filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_2d_int', values = readdata_int)

  do idim_2 = 1, size(readdata_int, 2)
    call assert_equal(readdata_int(:, idim_2), & 
                      testdata_int(:, idim_2, 1), &
                      'Retrieve 2D int')
  end do

  call check(nf90_close(ncid))

end subroutine test_nc_getvar_2d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_getvar_3d_int
  integer                     :: ncid, idim_2, idim_3
  integer, allocatable        :: readdata_int(:,:,:)

  call nc_open_for_read(filename = nc_filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_int', values = readdata_int)

  do idim_2 = 1, size(readdata_int, 2)
    do idim_3 = 1, size(readdata_int, 3)
      call assert_equal(readdata_int(:, idim_2, idim_3), & 
                        testdata_int(:, idim_2, idim_3), &
                        'Retrieve 3D int')
    end do
  end do

  call check(nf90_close(ncid))

end subroutine test_nc_getvar_3d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_putvar_1d
  integer                       :: ncid, idim_1, idim_2, idim_3
  character(len=128)            :: filename='./output/test_nc_putvar_1d.nc'
  real(kind=sp), allocatable    :: readdata_fp(:)
  integer,       allocatable    :: readdata_int(:)
  
  dim_size = [10, 20, 30]

  if (allocated(testdata_fp)) deallocate(testdata_fp)
  allocate(testdata_fp(dim_size(1), dim_size(2), dim_size(3) ))

  if (allocated(testdata_int)) deallocate(testdata_int)
  allocate(testdata_int(dim_size(1), dim_size(2), dim_size(3) ))

  do idim_1 = 1, dim_size(1)
    do idim_2 = 1, dim_size(2)
      do idim_3 = 1, dim_size(3)
        testdata_fp(idim_1, idim_2, idim_3) =   idim_1                       &
                                              + dim_size(1) * idim_2         &
                                              + product(dim_size(1:2)) * idim_3
        testdata_int(idim_1, idim_2, idim_3) = idim_1                       &
                                             + dim_size(1) * idim_2         &
                                             + product(dim_size(1:2)) * idim_3
      end do
    end do
  end do

  call nc_create_file(filename = filename, ncid = ncid, overwrite = .true.)

  call nc_putvar_by_name(ncid = ncid, varname = 'var_1d_float', values = testdata_fp(:,1,1))
  call nc_putvar_by_name(ncid = ncid, varname = 'var_1d_int',   values = testdata_int(:,1,1))

  call nc_close_file(ncid = ncid)

  call nc_open_for_read(filename = filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_1d_float', values = readdata_fp)
  call nc_getvar_by_name(ncid = ncid, varname = 'var_1d_int',   values = readdata_int)

  call assert_comparable(readdata_fp(:),       &
                         testdata_fp(:, 1, 1), &
                         1e-8, 'Retrieve 1D Float')
  call assert_equal(readdata_int(:),       &
                    testdata_int(:, 1, 1), &
                    'Retrieve 1D Integer')

  call nc_close_file(ncid = ncid)

end subroutine test_nc_putvar_1d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine test_nc_putvar_2d
  integer                       :: ncid, idim_1, idim_2, idim_3
  character(len=128)            :: filename='./output/test_nc_putvar_2d.nc'
  real(kind=sp), allocatable    :: readdata_fp(:,:)
  integer,       allocatable    :: readdata_int(:,:)
  
  
  dim_size = [10, 20, 30]

  if (allocated(testdata_fp)) deallocate(testdata_fp)
  allocate(testdata_fp(dim_size(1), dim_size(2), dim_size(3) ))

  if (allocated(testdata_int)) deallocate(testdata_int)
  allocate(testdata_int(dim_size(1), dim_size(2), dim_size(3) ))

  do idim_1 = 1, dim_size(1)
    do idim_2 = 1, dim_size(2)
      do idim_3 = 1, dim_size(3)
        testdata_fp(idim_1, idim_2, idim_3) =   idim_1                       &
                                              + dim_size(1) * idim_2         &
                                              + product(dim_size(1:2)) * idim_3
        testdata_int(idim_1, idim_2, idim_3) = idim_1                       &
                                             + dim_size(1) * idim_2         &
                                             + product(dim_size(1:2)) * idim_3
      end do
    end do
  end do

  call nc_create_file(filename = filename, ncid = ncid, overwrite = .true.)

  call nc_putvar_by_name(ncid = ncid, varname = 'var_2d_float', values = testdata_fp(:,:,1))
  call nc_putvar_by_name(ncid = ncid, varname = 'var_2d_int',   values = testdata_int(:,:,1))

  call nc_close_file(ncid = ncid)

  call nc_open_for_read(filename = filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_2d_float', values = readdata_fp)
  call nc_getvar_by_name(ncid = ncid, varname = 'var_2d_int',   values = readdata_int)

  do idim_2 = 1, size(readdata_fp, 2)
    call assert_comparable(readdata_fp(:, idim_2),    &
                           testdata_fp(:, idim_2, 1), &
                           1e-8, 'Retrieve 2D Float')
    call assert_equal(readdata_int(:, idim_2),    &
                      testdata_int(:, idim_2, 1), &
                      'Retrieve 2D Integer')
  end do

  call nc_close_file(ncid = ncid)

end subroutine test_nc_putvar_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine test_nc_putvar_3d
  integer                       :: ncid, idim_1, idim_2, idim_3
  character(len=128)            :: filename='./output/test_nc_putvar_3d.nc'
  real(kind=sp), allocatable    :: readdata_fp(:,:,:)
  integer,       allocatable    :: readdata_int(:,:,:)
  
  
  dim_size = [10, 20, 30]

  if (allocated(testdata_fp)) deallocate(testdata_fp)
  allocate(testdata_fp(dim_size(1), dim_size(2), dim_size(3) ))

  if (allocated(testdata_int)) deallocate(testdata_int)
  allocate(testdata_int(dim_size(1), dim_size(2), dim_size(3) ))

  do idim_1 = 1, dim_size(1)
    do idim_2 = 1, dim_size(2)
      do idim_3 = 1, dim_size(3)
        testdata_fp(idim_1, idim_2, idim_3) = idim_1                       &
                                            + dim_size(1) * idim_2         &
                                            + product(dim_size(1:2)) * idim_3
        testdata_int(idim_1, idim_2, idim_3) = idim_1                       &
                                             + dim_size(1) * idim_2         &
                                             + product(dim_size(1:2)) * idim_3
      end do
    end do
  end do

  call nc_create_file(filename = filename, ncid = ncid, overwrite = .true.)

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_float', values = testdata_fp)
  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_int',   values = testdata_int)

  call nc_close_file(ncid = ncid)

  call nc_open_for_read(filename = filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_float', values = readdata_fp)
  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_int',   values = readdata_int)

  do idim_2 = 1, size(readdata_fp, 2)
    do idim_3 = 1, size(readdata_fp, 3)
      call assert_comparable(readdata_fp(:, idim_2, idim_3), &
                             testdata_fp(:, idim_2, idim_3), &
                             1e-8, 'Retrieve 3D Float')
      call assert_equal(readdata_int(:, idim_2, idim_3), &
                        testdata_int(:, idim_2, idim_3), &
                        'Retrieve 3D Integer')
    end do
  end do

  call nc_close_file(ncid = ncid)

end subroutine test_nc_putvar_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_nc_putvar_1d_into_nd
  integer                       :: ncid, idim_1, idim_2, idim_3
  character(len=128)            :: filename='./output/test_nc_putvar_1d_into_nd.nc'
  real(kind=sp), allocatable    :: readdata_fp(:,:,:)
  real(kind=sp)                 :: array_1d_fp(10)
  integer,       allocatable    :: readdata_int(:,:,:)
  integer                       :: array_1d_int(10)
  
  
  dim_size = [10, 20, 30]

  if (allocated(testdata_fp)) deallocate(testdata_fp)
  allocate(testdata_fp(dim_size(1), dim_size(2), dim_size(3) ))

  if (allocated(testdata_int)) deallocate(testdata_int)
  allocate(testdata_int(dim_size(1), dim_size(2), dim_size(3) ))

  do idim_1 = 1, dim_size(1)
    do idim_2 = 1, dim_size(2)
      do idim_3 = 1, dim_size(3)
        testdata_fp(idim_1, idim_2, idim_3) = idim_1                       &
                                            + dim_size(1) * idim_2         &
                                            + product(dim_size(1:2)) * idim_3
        testdata_int(idim_1, idim_2, idim_3) = idim_1                       &
                                             + dim_size(1) * idim_2         &
                                             + product(dim_size(1:2)) * idim_3
      end do
    end do
  end do

  array_1d_fp  = [-1., -1., -1., -1., -1., -1., -1., -1., -1., -1.]
  array_1d_int = [-1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1]

  ! Test routine nc_putvar_by_name_1d_into_nd and nc_putvar_by_name_1d_into_nd_int along first dimension
  call nc_create_file (filename = filename, ncid = ncid)

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_float', values = testdata_fp)
  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_int',   values = testdata_int)

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_float', values = array_1d_fp, &
                         start = [1, 5, 1], count = [10, 1, 1])

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_int', values = array_1d_int, &
                         start = [1, 5, 1], count = [10, 1, 1])

  call nc_close_file(ncid = ncid)

  call nc_open_for_read(filename = filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_float', values = readdata_fp)
  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_int',   values = readdata_int)

  call assert_comparable(readdata_fp(1:10, 5, 1),   &
                         array_1d_fp,               &
                         1e-8, 'Retrieve 3D Float')
  call assert_equal(readdata_int(1:10, 5, 1), &
                    array_1d_int,             &
                    'Retrieve 3D Integer')

  call nc_close_file(ncid = ncid)

  ! Test routine nc_putvar_by_name_1d_into_nd and nc_putvar_by_name_1d_into_nd_int along second dimension
  call nc_open_for_write(filename = filename, ncid = ncid)

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_float', values = testdata_fp)
  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_int',   values = testdata_int)

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_float', values = array_1d_fp, &
                         start = [3, 1, 1], count = [1, 10, 1])

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_int', values = array_1d_int, &
                         start = [3, 1, 1], count = [1, 10, 1])

  call nc_close_file(ncid = ncid)

  call nc_open_for_read(filename = filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_float', values = readdata_fp)
  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_int',   values = readdata_int)

  call assert_comparable(readdata_fp(3, 1:10, 1),   &
                         array_1d_fp,               &
                         1e-8, 'Retrieve 3D Float')
  call assert_equal(readdata_int(3, 1:10, 1), &
                    array_1d_int,             &
                    'Retrieve 3D Integer')

  call nc_close_file(ncid = ncid)

  ! Test routine nc_putvar_by_name_1d_into_nd and nc_putvar_by_name_1d_into_nd_int along third dimension
  call nc_open_for_write(filename = filename, ncid = ncid)

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_float', values = testdata_fp)
  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_int',   values = testdata_int)

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_float', values = array_1d_fp, &
                         start = [3, 1, 1], count = [1, 1, 10])

  call nc_putvar_by_name(ncid = ncid, varname = 'var_3d_int', values = array_1d_int, &
                         start = [3, 1, 1], count = [1, 1, 10])

  call nc_close_file(ncid = ncid)

  call nc_open_for_read(filename = filename, ncid = ncid)

  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_float', values = readdata_fp)
  call nc_getvar_by_name(ncid = ncid, varname = 'var_3d_int',   values = readdata_int)

  call assert_comparable(readdata_fp(3, 1, 1:10),   &
                         array_1d_fp,               &
                         1e-8, 'Retrieve 3D Float')
  call assert_equal(readdata_int(3, 1, 1:10), &
                    array_1d_int,             &
                    'Retrieve 3D Integer')

  call nc_close_file(ncid = ncid)

end subroutine test_nc_putvar_1d_into_nd
!-----------------------------------------------------------------------------------------

end module test_nc_routines
