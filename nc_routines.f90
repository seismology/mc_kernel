module nc_routines
    use global_parameters, only: sp, dp, lu_out, verbose, myrank
    use commpi,            only: pabort    
    use simple_routines,   only: check_limits
    use netcdf
    implicit none

    interface nc_getvar_by_name
      module procedure  :: nc_getvar_by_name_1d
      module procedure  :: nc_getvar_by_name_2d
      module procedure  :: nc_getvar_by_name_3d
      module procedure  :: nc_getvar_by_name_1d_int
      module procedure  :: nc_getvar_by_name_2d_int
      module procedure  :: nc_getvar_by_name_3d_int
    end interface nc_getvar_by_name

    interface nc_putvar_by_name
      module procedure  :: nc_putvar_by_name_1d
      module procedure  :: nc_putvar_by_name_2d
      module procedure  :: nc_putvar_by_name_3d
      module procedure  :: nc_putvar_by_name_1d_into_nd
      module procedure  :: nc_putvar_by_name_1d_int
      module procedure  :: nc_putvar_by_name_2d_int
      module procedure  :: nc_putvar_by_name_3d_int
      module procedure  :: nc_putvar_by_name_1d_into_nd_int
    end interface nc_putvar_by_name

    interface nc_getvar
      module procedure  :: getvar_real1d
      module procedure  :: getvar_real1d_dble
      module procedure  :: getvar_int1d
      module procedure  :: getvar_real2d
      module procedure  :: getvar_int2d
      module procedure  :: getvar_real3d
      module procedure  :: getvar_int3d
    end interface nc_getvar

    interface nc_putvar
      module procedure  :: putvar_real1d
      module procedure  :: putvar_real2d
      module procedure  :: putvar_real3d
    end interface nc_putvar

contains

!-----------------------------------------------------------------------------------------
subroutine nc_open_for_read(filename, ncid)
   character(len=*), intent(in)  :: filename
   integer, intent(out)          :: ncid
   character(len=512)            :: fmtstring
   integer                       :: status

   status = nf90_open(path     = filename,              &
                      mode     = NF90_NOWRITE,          &
                      ncid     = ncid)

   if (status.ne.NF90_NOERR) then
      fmtstring = "('ERROR: CPU ', I4, ' tried to open file ''', A, ''', " &
                    // "but could not find it')"
      print fmtstring, myrank, trim(filename)
      stop
   end if

end subroutine nc_open_for_read
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_create_file(filename, ncid, cache_size, overwrite)
   character(len=*), intent(in)  :: filename
   integer, intent(out)          :: ncid
   integer, intent(in), optional :: cache_size
   logical, intent(in), optional :: overwrite
   character(len=512)            :: fmtstring
   integer                       :: status, mode
   integer                       :: cache_size_local = 512*2**20 !512 MB default value

   mode = NF90_NETCDF4

   if (present(cache_size)) cache_size_local = cache_size

   if (present(overwrite)) then
     if (overwrite) mode = mode + NF90_CLOBBER
   end if

   status = nf90_create(path       = filename,         &
                        cmode      = mode,             &
                        cache_size = cache_size_local, &
                        ncid       = ncid)

   if (status.eq.NF90_EEXIST) then
     fmtstring = "('ERROR: CPU ', I4, ' tried to create file ''', A, ''', " &
              // "but it already exists and overwrite was not set.')"
     print fmtstring, myrank, trim(filename)
     stop
   elseif (status.ne.NF90_NOERR) then
     fmtstring = "('ERROR: CPU ', I4, ' tried to create file ''', A, ''', " &
              // "but could not. Check permissions.')"
     print fmtstring, myrank, trim(filename)
     stop
   end if
   call check(nf90_enddef(ncid = ncid))

end subroutine nc_create_file
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_create_group(group_name, ncid, ncid_group, overwrite)
   character(len=*), intent(in)  :: group_name
   integer, intent(in)           :: ncid
   logical, intent(in), optional :: overwrite
   integer, intent(inout)        :: ncid_group
   character(len=512)            :: fmtstring
   logical                       :: overwrite_loc = .false.
   integer                       :: status

   if (present(overwrite)) overwrite_loc = overwrite

   status = nf90_inq_ncid(ncid = ncid, name = group_name, grp_ncid = ncid_group)

   select case(status)
   case(NF90_NOERR) 
     if (.not.overwrite_loc) then
       fmtstring = "('ERROR: CPU ', I4, ' tried to create group ''', A, ''', " &
                // "but it already exists.')"
       print fmtstring, myrank, trim(group_name)
       stop
     end if

   case default
     call check(nf90_redef(ncid = ncid))
     status = nf90_def_grp(ncid, trim(group_name), ncid_group)

     if (status.ne.NF90_NOERR) then
       fmtstring = "('ERROR: CPU ', I4, ' tried to create group ''', A, ''', " &
                // "but could not. Check permissions.')"
       print fmtstring, myrank, trim(group_name)
       print *, nf90_strerror(status)
       stop
     end if
     call check(nf90_enddef(ncid = ncid))
   end select

end subroutine nc_create_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_open_for_write(filename, ncid, cache_size)
   character(len=*), intent(in)  :: filename
   integer, intent(out)          :: ncid
   integer, intent(in), optional :: cache_size
   character(len=512)            :: fmtstring
   integer                       :: status, mode
   integer                       :: cache_size_local = 512*2**20 !512 MB default value

   if (present(cache_size)) cache_size_local = cache_size

   status = nf90_open(path       = filename,         &
                      mode       = NF90_WRITE,       &
                      cache_size = cache_size_local, &
                      ncid       = ncid)

   if (status.ne.NF90_NOERR) then
     fmtstring = "('CPU ', I4, ' tried to open file ''', A, ''' for writing, " &
              // "but could not find it. Creates it instead.')"
     print fmtstring, myrank, trim(filename)
     status = nf90_create(path       = filename,         &
                          cmode      = NF90_NETCDF4,   &
                          cache_size = cache_size_local, &
                          ncid       = ncid)

     if (status.ne.NF90_NOERR) then
        fmtstring = "('ERROR: CPU ', I4, ' tried to create file ''', A, ''', " &
                 // "but could not. Check permissions.')"
        print fmtstring, myrank, trim(filename)
        stop
     end if
     call check(nf90_enddef(ncid = ncid))
   end if

end subroutine nc_open_for_write
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_close_file(ncid)
   integer, intent(in)           :: ncid
   integer                       :: status

   status = nf90_close(ncid = ncid)

end subroutine nc_close_file
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvarid(ncid, name, varid, noerr)
    integer, intent(in)           :: ncid
    character(len=*), intent(in)  :: name
    integer, intent(out)          :: varid
    logical, intent(in), optional :: noerr
    integer                       :: status

    status = nf90_inq_varid( ncid  = ncid, &
                             name  = name, &
                             varid = varid )

    if (status.ne.NF90_NOERR) then
      ! Variable not found, return varid -1
      if (present(noerr)) then
        if (noerr) then
          varid = -1
        else
          write(6,100) myrank, trim(name), ncid
          stop
        end if
      else
        write(6,100) myrank, trim(name), ncid
        stop
      end if
      if (verbose>1) then
        write(6,100) myrank, trim(name), ncid
        call flush(6)
      end if

    elseif (verbose>1) then
      ! Variable found
      write(lu_out,101) trim(name), ncid, varid
      call flush(lu_out)

    end if

100 format('ERROR: CPU ', I4, ' could not find variable: ''', A, ''' in NCID', I7)
101 format('    Variable ''', A, ''' found in NCID', I7, ', has ID:', I7)
end subroutine getvarid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getgrpid(ncid, name, grp_ncid)
    integer, intent(in)          :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out)         :: grp_ncid
    integer                      :: status

    status = nf90_inq_ncid( ncid     = ncid, &
                            name     = name, &
                            grp_ncid = grp_ncid )
    if (status.ne.NF90_NOERR) then
        write(6,100) myrank, trim(name), ncid
        call flush(6)
        call pabort
    elseif (verbose>1) then
        write(lu_out,101) trim(name), ncid, grp_ncid
        call flush(lu_out)
    end if
100 format('ERROR: CPU ', I4, ' could not find group: ''', A, ''' in NCID', I7)
101 format('    Group ''', A, ''' found in NCID', I7, ', has ID:', I7)
end subroutine getgrpid
!-----------------------------------------------------------------------------------------
 
!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_1d_int(ncid, varname, values, limits)
!< Looks up a 1D variable name and returns the complete variable
  integer, intent(in)                :: ncid
  character(len=*), intent(in)       :: varname
  integer, allocatable, intent(out)  :: values(:)
  integer, intent(in), optional      :: limits(2)

  integer                            :: variable_id, dimid(1), npoints, variable_type
  integer                            :: status
  integer                            :: limits_loc(2)
  logical                            :: have_limits = .false., out_of_limit = .false.

  if (verbose>1) then
    write(lu_out,"(' Trying to read 1D int variable ', A, '...')") trim(varname)
    call flush(lu_out)
  end if

  call getvarid( ncid  = ncid,      &
                 name  = varname,   &
                 varid = variable_id)

  status = nf90_inquire_variable(ncid   = ncid,           &
                                 varid  = variable_id,    &
                                 xtype  = variable_type,  &
                                 dimids = dimid  )

  status = nf90_inquire_dimension(ncid  = ncid,           &
                                  dimid = dimid(1),       &
                                  len   = npoints)

  allocate(values(npoints))

  select case(variable_type)
  case(NF90_INT)
    call nc_getvar(ncid   = ncid,                 &
                   varid  = variable_id,          &
                   start  = 1,                    &
                   count  = npoints,              &
                   values = values(:)) 
  case default
    write(*,*) 'nc_getvar_by_name_1d_int works only for NF90_INT variables'
    call flush(6)
    call pabort()
  end select

  if (present(limits)) then
     limits_loc = limits
     have_limits = .true.
  else
     status = nf90_get_att(ncid   = ncid,          &
                           varid  = variable_id,   &
                           name   = 'valid_range', &
                           values = limits_loc) 
     if (status.eq.NF90_NOERR) have_limits = .true.
  end if

  if (have_limits) then
     if (minval(values).lt.(minval(limits_loc))) then
        write(*,*) 'ERROR: Value in NetCDF file smaller than limit!'
        write(*,*) 'Variable name: ', trim(varname) 
        write(*,*) 'Element nr   : ', minloc(values)
        write(*,*) 'Element value: ', minval(values)
        write(*,*) 'Lower limit  : ', minval(limits_loc)
        out_of_limit = .true.
     end if
     if (maxval(values).gt.(maxval(limits_loc))) then
        write(*,*) 'ERROR: Value in NetCDF file larger than limit!'
        write(*,*) 'Variable name: ', trim(varname) 
        write(*,*) 'Element nr   : ', maxloc(values)
        write(*,*) 'Element value: ', maxval(values)
        write(*,*) 'Upper limit  : ', maxval(limits_loc)
        out_of_limit = .true.
     end if
     if (out_of_limit) then
       call flush(6)
       call pabort(do_traceback=.false.)
     end if
  end if

end subroutine nc_getvar_by_name_1d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_1d(ncid, varname, values, limits)
!< Looks up a 1D variable name and returns the complete variable
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  real(kind=sp), allocatable, intent(inout)  :: values(:)
  real(kind=sp), intent(in), optional        :: limits(2)

  integer                                    :: variable_id, dimid(1), npoints, variable_type
  integer                                    :: status
  real(kind=sp)                              :: limits_loc(2)
  logical                                    :: have_limits = .false., out_of_limit = .false.


  if (verbose>1) then
    write(lu_out,"(' Trying to read 1D fp variable ', A, '...')") trim(varname)
    call flush(lu_out)
  end if

  call  getvarid( ncid  = ncid,            &
                  name  = varname,   &
                  varid = variable_id)
 
  ! Inquire variable type and dimension ids
  status = nf90_inquire_variable(ncid   = ncid,           &
                                 varid  = variable_id,    &
                                 xtype  = variable_type,  &
                                 dimids = dimid  )

  ! Inquire size of dimension 1                            
  status = nf90_inquire_dimension(ncid  = ncid,           &
                                  dimid = dimid(1),       &
                                  len   = npoints)

  ! Allocate output variable with size of NetCDF variable
  allocate(values(npoints))

  select case(variable_type)
  case(NF90_FLOAT)
    call nc_getvar(ncid   = ncid,                 &
                   varid  = variable_id,          &
                   start  = 1,                    &
                   count  = npoints,              &
                   values = values) 
  case default
    write(*,*) 'nc_getvar_by_name_1d is only implemented for NF90_FLOAT'
    call flush(6)
    call pabort()
  end select

  if (present(limits)) then
     limits_loc = limits
     have_limits = .true.
  else
     status = nf90_get_att(ncid   = ncid,          &
                           varid  = variable_id,   &
                           name   = 'valid_range', &
                           values = limits_loc) 
     if (status.eq.NF90_NOERR) have_limits = .true.
  end if

  if (have_limits) then
     out_of_limit = check_limits(values, limits_loc, varname)
     if (out_of_limit) then
       call flush(6)
       call pabort(do_traceback=.false.)
     end if
  end if

end subroutine nc_getvar_by_name_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_2d_int(ncid, varname, values, limits)
!< Looks up a 2D variable name and returns the complete variable
  integer, intent(in)                :: ncid
  character(len=*), intent(in)       :: varname
  integer, allocatable, intent(out)  :: values(:,:)
  integer, intent(in), optional      :: limits(2)

  integer                            :: variable_id, dimid(2), len_dim1, len_dim2, variable_type
  integer                            :: status
  integer                            :: limits_loc(2)
  logical                            :: have_limits = .false., out_of_limit = .false.

  if (verbose>1) then
    write(lu_out,"(' Trying to read 2D int variable ', A, '...')") trim(varname)
    call flush(lu_out)
  end if

  call  getvarid( ncid  = ncid,            &
                  name  = varname,   &
                  varid = variable_id)

  status = nf90_inquire_variable(ncid   = ncid,           &
                                 varid  = variable_id,    &
                                 xtype  = variable_type,  &
                                 dimids = dimid  )

  status = nf90_inquire_dimension(ncid  = ncid,           &
                                  dimid = dimid(1),       &
                                  len   = len_dim1)

  status = nf90_inquire_dimension(ncid  = ncid,           &
                                  dimid = dimid(2),       &
                                  len   = len_dim2)

  allocate(values(len_dim1, len_dim2))

  select case(variable_type)
  case(NF90_INT)
    call nc_getvar(ncid   = ncid,                 &
                   varid  = variable_id,          &
                   start  = [1, 1],               &
                   count  = [len_dim1, len_dim2], &
                   values = values)
  case default
    write(*,*) 'nc_getvar_by_name_2d_int is only implemented for NF90_INT variables'
    call pabort()
  end select

  if (present(limits)) then
     limits_loc = limits
     have_limits = .true.
  else
     status = nf90_get_att(ncid   = ncid,          &
                           varid  = variable_id,   &
                           name   = 'valid_range', &
                           values = limits_loc) 
     if (status.eq.NF90_NOERR) have_limits = .true.
  end if

  if (have_limits) then
     if (minval(values).lt.(minval(limits_loc))) then
        write(*,*) 'ERROR: Value in NetCDF file smaller than limit!'
        write(*,*) 'Variable name: ', trim(varname) 
        write(*,*) 'Element nr   : ', minloc(values)
        write(*,*) 'Element value: ', minval(values)
        write(*,*) 'Lower limit  : ', minval(limits_loc)
        out_of_limit = .true.
     end if
     if (maxval(values).gt.(maxval(limits_loc))) then
        write(*,*) 'ERROR: Value in NetCDF file larger than limit!'
        write(*,*) 'Variable name: ', trim(varname) 
        write(*,*) 'Element nr   : ', maxloc(values)
        write(*,*) 'Element value: ', maxval(values)
        write(*,*) 'Upper limit  : ', maxval(limits_loc)
        out_of_limit = .true.
     end if
     if (out_of_limit) then
       call flush(6)
       call pabort(do_traceback=.false.)
     end if
  end if

end subroutine nc_getvar_by_name_2d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_2d(ncid, varname, values, limits)
!< Looks up a 2D variable name and returns the complete variable
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  real(kind=sp), allocatable, intent(inout)  :: values(:,:)
  real(kind=sp), intent(in), optional        :: limits(2)

  integer                                    :: variable_id, dimid(2), len_dim1, len_dim2
  integer                                    :: variable_type
  integer                                    :: status
  real(kind=sp)                              :: limits_loc(2)
  logical                                    :: have_limits = .false., out_of_limit = .false.

  if (verbose>1) then
    write(lu_out,"(' Trying to read 2D fp variable ', A, '...')") trim(varname)
    call flush(lu_out)
  end if

  call  getvarid( ncid  = ncid,       &
                  name  = varname,    &
                  varid = variable_id)

  status = nf90_inquire_variable(ncid   = ncid,           &
                                 varid  = variable_id,    &
                                 xtype  = variable_type,  &
                                 dimids = dimid  )

  status = nf90_inquire_dimension(ncid  = ncid,           &
                                  dimid = dimid(1),       &
                                  len   = len_dim1)

  status = nf90_inquire_dimension(ncid  = ncid,           &
                                  dimid = dimid(2),       &
                                  len   = len_dim2)

  allocate(values(len_dim1, len_dim2))

  select case(variable_type)
  case(NF90_FLOAT)
    call nc_getvar(ncid   = ncid,                 &
                   varid  = variable_id,          &
                   start  = [1, 1],               &
                   count  = [len_dim1, len_dim2], &
                   values = values) 
  case default
    write(*,*) 'nc_getvar_by_name_2d is only implemented for NF90_FLOAT variables'
    call flush(6)
    call pabort()
  end select

  if (present(limits)) then
     limits_loc = limits
     have_limits = .true.
  else
     status = nf90_get_att(ncid   = ncid,          &
                           varid  = variable_id,   &
                           name   = 'valid_range', &
                           values = limits_loc) 
     if (status.eq.NF90_NOERR) have_limits = .true.
  end if

  if (have_limits) then
     out_of_limit = check_limits(values, limits_loc, varname)
     if (out_of_limit) then
       call flush(6)
       call pabort(do_traceback=.false.)
     end if
  end if

end subroutine nc_getvar_by_name_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_3d_int(ncid, varname, values, limits)
!< Looks up a 3D variable name and returns the complete variable
  integer, intent(in)                :: ncid
  character(len=*), intent(in)       :: varname
  integer, allocatable, intent(out)  :: values(:,:,:)
  integer, intent(in), optional      :: limits(2)

  integer                            :: variable_id, dimid(3), len_dim(3), variable_type
  integer                            :: status, i_dim
  integer                            :: limits_loc(2)
  logical                            :: have_limits = .false., out_of_limit = .false.

  if (verbose>1) then
    write(lu_out,"(' Trying to read 3D int variable ', A, '...')") trim(varname)
    call flush(lu_out)
  end if

  call  getvarid( ncid  = ncid,            &
                  name  = varname,   &
                  varid = variable_id)


  status = nf90_inquire_variable(ncid   = ncid,           &
                                 varid  = variable_id,    &
                                 xtype  = variable_type,  &
                                 dimids = dimid  )

  do i_dim = 1, 3             
    status = nf90_inquire_dimension(ncid  = ncid,           &
                                    dimid = dimid(i_dim),       &
                                    len   = len_dim(i_dim))
  end do                                


  allocate(values(len_dim(1), len_dim(2), len_dim(3)))

  select case(variable_type)
  case(NF90_INT)
    call nc_getvar(ncid   = ncid,                 &
                   varid  = variable_id,          &
                   start  = [1, 1, 1],            &
                   count  = [len_dim(1), len_dim(2), len_dim(3)], &
                   values = values)
  case default
    write(*,*) 'nc_getvar_by_name_3d_int is only implemented for NF90_INT variables'
    call flush(6)
    call pabort()
  end select

  if (present(limits)) then
     limits_loc = limits
     have_limits = .true.
  else
     status = nf90_get_att(ncid   = ncid,          &
                           varid  = variable_id,   &
                           name   = 'valid_range', &
                           values = limits_loc) 
     if (status.eq.NF90_NOERR) have_limits = .true.
  end if

  if (have_limits) then
     if (minval(values).lt.(minval(limits_loc))) then
        write(*,*) 'ERROR: Value in NetCDF file smaller than limit!'
        write(*,*) 'Variable name: ', trim(varname) 
        write(*,*) 'Element nr   : ', minloc(values)
        write(*,*) 'Element value: ', minval(values)
        write(*,*) 'Lower limit  : ', minval(limits_loc)
        out_of_limit = .true.
     end if
     if (maxval(values).gt.(maxval(limits_loc))) then
        write(*,*) 'ERROR: Value in NetCDF file larger than limit!'
        write(*,*) 'Variable name: ', trim(varname) 
        write(*,*) 'Element nr   : ', maxloc(values)
        write(*,*) 'Element value: ', maxval(values)
        write(*,*) 'Upper limit  : ', maxval(limits_loc)
        out_of_limit = .true.
     end if
     if (out_of_limit) then
       call flush(6)
       call pabort(do_traceback=.false.)
     end if
  end if

end subroutine nc_getvar_by_name_3d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_3d(ncid, varname, values, limits)
!< Looks up a 3D variable name and returns the complete variable
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  real(kind=sp), allocatable, intent(out)    :: values(:,:,:)
  real(kind=sp), intent(in), optional        :: limits(2)

  integer                                    :: variable_id, dimid(3), len_dim(3), variable_type
  integer                                    :: status, i_dim
  real(kind=sp)                              :: limits_loc(2)
  logical                                    :: have_limits = .false., out_of_limit = .false.

  if (verbose>1) then
    write(lu_out,"(' Trying to read 3D fp variable ', A, '...')") trim(varname)
    call flush(lu_out)
  end if

  call  getvarid( ncid  = ncid,   &
                  name  = varname,   &
                  varid = variable_id)


  status = nf90_inquire_variable(ncid   = ncid,           &
                                 varid  = variable_id,    &
                                 xtype  = variable_type,  &
                                 dimids = dimid  )

  do i_dim = 1, 3             
    status = nf90_inquire_dimension(ncid  = ncid,           &
                                    dimid = dimid(i_dim),       &
                                    len   = len_dim(i_dim))
  end do                                


  allocate(values(len_dim(1), len_dim(2), len_dim(3)))

  select case(variable_type)
  case(NF90_FLOAT)
    call nc_getvar(ncid   = ncid,                 &
                   varid  = variable_id,          &
                   start  = [1, 1, 1],            &
                   count  = [len_dim(1), len_dim(2), len_dim(3)], &
                   values = values)
  case default
    write(*,*) 'nc_getvar_by_name_3d get_mesh_variable is only implemented for NF90_FLOAT variables'
    call flush(6)
    call pabort()
  end select

  if (present(limits)) then
     limits_loc = limits
     have_limits = .true.
  else
     status = nf90_get_att(ncid   = ncid,          &
                           varid  = variable_id,   &
                           name   = 'valid_range', &
                           values = limits_loc) 
     if (status.eq.NF90_NOERR) have_limits = .true.
  end if

  if (have_limits) then
     if (minval(values).lt.(minval(limits_loc))) then
        write(*,*) 'ERROR: Value in NetCDF file smaller than limit!'
        write(*,*) 'Variable name: ', trim(varname) 
        write(*,*) 'Element nr   : ', minloc(values)
        write(*,*) 'Element value: ', minval(values)
        write(*,*) 'Lower limit  : ', minval(limits_loc)
        out_of_limit = .true.
     end if
     if (maxval(values).gt.(maxval(limits_loc))) then
        write(*,*) 'ERROR: Value in NetCDF file larger than limit!'
        write(*,*) 'Variable name: ', trim(varname) 
        write(*,*) 'Element nr   : ', maxloc(values)
        write(*,*) 'Element value: ', maxval(values)
        write(*,*) 'Upper limit  : ', maxval(limits_loc)
        out_of_limit = .true.
     end if
     if (out_of_limit) then
       call flush(6)
       call pabort(do_traceback=.false.)
     end if
  end if

end subroutine nc_getvar_by_name_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_putvar_by_name_1d(ncid, varname, values, start, count)
!< Write data into a 1D variable, or create it, if necessary
  integer, parameter                         :: ndim = 1
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  real(kind=sp), intent(in)                  :: values(:)
  integer, intent(in), optional              :: start, count
  character(len=nf90_max_name)               :: dimname
  character(len=80)                          :: fmtstring
  integer                                    :: i, variable_id, dimid(ndim), status
  integer                                    :: ncid_root

  call getvarid(ncid  = ncid,        &
                name  = varname,     &
                varid = variable_id, &
                noerr = .true.)

  ! Get ncid of root group. Used for dimensions, which should always be defined
  ! in the root group (at least according to my taste)
  status = 0
  ncid_root = ncid
  do while (status.eq.NF90_NOERR) !Check whether we are already in root group
    status = nf90_inq_grp_parent(ncid_root, ncid_root) 
  end do

  ! if variable not found              
  if (variable_id==-1) then
    if (present(start).and.present(count)) then
      fmtstring = "('ERROR: Variable with name ', A, ' cannot be created.')"
      print fmtstring, trim(varname)
      print '(A)', "Arguments ''start'' or ''count'' are only allowed, if variable already exists"
      stop
    else
      call check(nf90_redef(ncid = ncid))
      do i = 1, ndim
        write(dimname, '(A,"_",I1)') trim(varname), i
        status = nf90_def_dim(ncid  = ncid,     &
                              name  = dimname,          &
                              len   = size(values, i),  &
                              dimid = dimid(i))
      end do

      status = nf90_def_var( ncid   = ncid,          &
                             name   = trim(varname), &
                             xtype  = NF90_FLOAT,    &
                             dimids = dimid,         &
                             varid  = variable_id) 
      call check(nf90_enddef(ncid = ncid))
    end if
  end if

  if (present(start).and.present(count)) then
    call putvar_real1d(ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = start,       &
                       count  = count)
  else
    call putvar_real1d(ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = 1,           &
                       count  = size(values, 1) )
  end if

end subroutine nc_putvar_by_name_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_putvar_by_name_1d_int(ncid, varname, values, start, count)
!< Write data into a 1D variable, or create it, if necessary
  integer, parameter                         :: ndim = 1
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  integer, intent(in)                        :: values(:)
  integer, intent(in), optional              :: start, count
  character(len=nf90_max_name)               :: dimname
  character(len=80)                          :: fmtstring
  integer                                    :: i, variable_id, dimid(ndim), status
  integer                                    :: ncid_root

  call getvarid(ncid  = ncid,        &
                name  = varname,     &
                varid = variable_id, &
                noerr = .true.)

  ! Get ncid of root group. Used for dimensions, which should always be defined
  ! in the root group (at least according to my taste)
  status = 0
  ncid_root = ncid
  do while (status.eq.NF90_NOERR) !Check whether we are already in root group
    status = nf90_inq_grp_parent(ncid_root, ncid_root) 
  end do

  ! if variable not found              
  if (variable_id==-1) then
    if (present(start).and.present(count)) then
      fmtstring = "('ERROR: Variable with name ', A, ' cannot be created.')"
      print fmtstring, trim(varname)
      print '(A)', "Arguments ''start'' or ''count'' are only allowed, if variable already exists"
      stop
    else
      call check(nf90_redef(ncid = ncid))
      do i = 1, ndim
        write(dimname, '(A,"_",I1)') trim(varname), i
        status = nf90_def_dim(ncid  = ncid,     &
                              name  = dimname,          &
                              len   = size(values, i),  &
                              dimid = dimid(i))
      end do

      status = nf90_def_var( ncid   = ncid,          &
                             name   = trim(varname), &
                             xtype  = NF90_INT,      &
                             dimids = dimid,         &
                             varid  = variable_id) 
      call check(nf90_enddef(ncid = ncid))
    end if
  end if

  if (present(start).and.present(count)) then
    call putvar_int1d( ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = start,       &
                       count  = count)
  else
    call putvar_int1d (ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = 1,           &
                       count  = size(values, 1) )
  end if

end subroutine nc_putvar_by_name_1d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_putvar_by_name_2d(ncid, varname, values, start, count)
!< Write data into a 2D variable, or create it, if necessary
  integer, parameter                         :: ndim = 2
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  real(kind=sp), intent(in)                  :: values(:,:)
  integer, intent(in), optional              :: start(ndim), count(ndim)

  character(len=nf90_max_name)               :: dimname
  character(len=80)                          :: fmtstring
  integer                                    :: i, variable_id, dimid(ndim), status
  integer                                    :: ncid_root

  call getvarid(ncid  = ncid,        &
                name  = varname,     &
                varid = variable_id, &
                noerr = .true.)

  ! Get ncid of root group. Used for dimensions, which should always be defined
  ! in the root group (at least according to my taste)
  status = 0
  ncid_root = ncid
  do while (status.eq.NF90_NOERR) !Check whether we are already in root group
    status = nf90_inq_grp_parent(ncid_root, ncid_root) 
  end do

  ! if variable not found              
  if (variable_id==-1) then
    if (present(start).and.present(count)) then
      fmtstring = "('ERROR: Variable with name ', A, ' cannot be created.')"
      print fmtstring, trim(varname)
      print '(A)', "Arguments ''start'' or ''count'' are only allowed, if variable already exists"
      stop
    else
      call check(nf90_redef(ncid = ncid))
      do i = 1, ndim
        write(dimname, '(A,"_",I1)') trim(varname), i
        status = nf90_def_dim(ncid  = ncid,     &
                              name  = dimname,          &
                              len   = size(values, i),  &
                              dimid = dimid(i))
      end do

      status = nf90_def_var( ncid   = ncid,          &
                             name   = trim(varname), &
                             xtype  = NF90_FLOAT,    &
                             dimids = dimid,         &
                             varid  = variable_id) 
      call check(nf90_enddef(ncid = ncid))
    end if
  end if

  if (present(start).and.present(count)) then
    call putvar_real2d(ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = start,       &
                       count  = count)
  else
    call putvar_real2d(ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = [1, 1],      &
                       count  = [size(values, 1), size(values, 2)])
  end if


end subroutine nc_putvar_by_name_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_putvar_by_name_2d_int(ncid, varname, values, start, count)
!< Write data into a 2D variable, or create it, if necessary
  integer, parameter                         :: ndim = 2
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  integer, intent(in)                        :: values(:,:)
  integer, intent(in), optional              :: start(ndim), count(ndim)

  character(len=nf90_max_name)               :: dimname
  character(len=80)                          :: fmtstring
  integer                                    :: i, variable_id, dimid(ndim), status
  integer                                    :: ncid_root

  call getvarid(ncid  = ncid,        &
                name  = varname,     &
                varid = variable_id, &
                noerr = .true.)

  ! Get ncid of root group. Used for dimensions, which should always be defined
  ! in the root group (at least according to my taste)
  status = 0
  ncid_root = ncid
  do while (status.eq.NF90_NOERR) !Check whether we are already in root group
    status = nf90_inq_grp_parent(ncid_root, ncid_root) 
  end do

  ! if variable not found              
  if (variable_id==-1) then
    if (present(start).and.present(count)) then
      fmtstring = "('ERROR: Variable with name ', A, ' cannot be created.')"
      print fmtstring, trim(varname)
      print '(A)', "Arguments ''start'' or ''count'' are only allowed, if variable already exists"
      stop
    else
      call check(nf90_redef(ncid = ncid))
      do i = 1, ndim
        write(dimname, '(A,"_",I1)') trim(varname), i
        status = nf90_def_dim(ncid  = ncid,     &
                              name  = dimname,          &
                              len   = size(values, i),  &
                              dimid = dimid(i))
      end do

      status = nf90_def_var( ncid   = ncid,          &
                             name   = trim(varname), &
                             xtype  = NF90_INT,      &
                             dimids = dimid,         &
                             varid  = variable_id) 
      call check(nf90_enddef(ncid = ncid))
    end if
  end if

  if (present(start).and.present(count)) then
    call putvar_int2d( ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = start,       &
                       count  = count)
  else
    call putvar_int2d( ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = [1, 1],      &
                       count  = [size(values, 1), size(values, 2)])
  end if


end subroutine nc_putvar_by_name_2d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_putvar_by_name_3d(ncid, varname, values, start, count)
!< Write data into a 3D variable, or create it, if necessary
  integer, parameter                         :: ndim = 3
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  real(kind=sp), intent(in)                  :: values(:,:,:)
  integer, intent(in), optional              :: start(ndim), count(ndim)
  character(len=nf90_max_name)               :: dimname
  character(len=80)                          :: fmtstring
  integer                                    :: i, variable_id, dimid(ndim), status
  integer                                    :: ncid_root

  call getvarid(ncid  = ncid,        &
                name  = varname,     &
                varid = variable_id, &
                noerr = .true.)

  ! Get ncid of root group. Used for dimensions, which should always be defined
  ! in the root group (at least according to my taste)
  status = 0
  ncid_root = ncid
  do while (status.eq.NF90_NOERR) !Check whether we are already in root group
    status = nf90_inq_grp_parent(ncid_root, ncid_root) 
  end do

  ! if variable not found              
  if (variable_id==-1) then
    if (present(start).and.present(count)) then
      fmtstring = "('ERROR: Variable with name ', A, ' cannot be created.')"
      print fmtstring, trim(varname)
      print '(A)', "Arguments ''start'' or ''count'' are only allowed, if variable already exists"
      stop
    else
      call check(nf90_redef(ncid = ncid))
      do i = 1, ndim
        write(dimname, '(A,"_",I1)') trim(varname), i
        status = nf90_def_dim(ncid  = ncid,     &
                              name  = dimname,          &
                              len   = size(values, i),  &
                              dimid = dimid(i))
      end do

      status = nf90_def_var( ncid   = ncid,          &
                             name   = trim(varname), &
                             xtype  = NF90_FLOAT,    &
                             dimids = dimid,         &
                             varid  = variable_id) 
      call check(nf90_enddef(ncid = ncid))
    end if
  end if

  if (present(start).and.present(count)) then
    call putvar_real3d(ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = start,       &
                       count  = count)
  else
    call putvar_real3d(ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = [1, 1, 1],   &
                       count  = [size(values, 1), size(values, 2), size(values, 3)])
  end if


end subroutine nc_putvar_by_name_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_putvar_by_name_3d_int(ncid, varname, values, start, count)
!< Write data into a 3D variable, or create it, if necessary
  integer, parameter                         :: ndim = 3
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  integer, intent(in)                        :: values(:,:,:)
  integer, intent(in), optional              :: start(ndim), count(ndim)
  character(len=nf90_max_name)               :: dimname
  character(len=80)                          :: fmtstring
  integer                                    :: i, variable_id, dimid(ndim), status
  integer                                    :: ncid_root

                  
  call getvarid(ncid  = ncid,        &
                name  = varname,     &
                varid = variable_id, &
                noerr = .true.)

  ! Get ncid of root group. Used for dimensions, which should always be defined
  ! in the root group (at least according to my taste)
  status = 0
  ncid_root = ncid
  do while (status.eq.NF90_NOERR) !Check whether we are already in root group
    status = nf90_inq_grp_parent(ncid_root, ncid_root) 
  end do

  ! if variable not found              
  if (variable_id==-1) then
    if (present(start).and.present(count)) then
      fmtstring = "('ERROR: Variable with name ', A, ' cannot be created.')"
      print fmtstring, trim(varname)
      print '(A)', "Arguments ''start'' or ''count'' are only allowed, if variable already exists"
      stop
    else
      call check(nf90_redef(ncid = ncid))
      do i = 1, ndim
        write(dimname, '(A,"_",I1)') trim(varname), i
        status = nf90_def_dim(ncid  = ncid,     &
                              name  = dimname,          &
                              len   = size(values, i),  &
                              dimid = dimid(i))
      end do

      status = nf90_def_var( ncid   = ncid,          &
                             name   = trim(varname), &
                             xtype  = NF90_INT,      &
                             dimids = dimid,         &
                             varid  = variable_id) 
      call check(nf90_enddef(ncid = ncid))
    end if
  end if

  if (present(start).and.present(count)) then
    call putvar_int3d( ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = start,       &
                       count  = count)
  else
    call putvar_int3d( ncid   = ncid,        &
                       varid  = variable_id, &
                       values = values,      &
                       start  = [1, 1, 1],   &
                       count  = [size(values, 1), size(values, 2), size(values, 3)])
  end if


end subroutine nc_putvar_by_name_3d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_putvar_by_name_1d_into_nd(ncid, varname, values, start, count)
!< Write data into a 3D variable, or create it, if necessary
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  real(kind=sp), intent(in)                  :: values(:)
  integer, intent(in)                        :: start(:), count(:)
  character(len=80)                          :: fmtstring
  integer                                    :: variable_id 

  ! A rather preculiar way of determining that count is larger 1 on all but one dimensio
  ! We can't use the FORTRAN intrinsic 'count' here, well, because we have a variable
  ! of this name. sucks
  if (product(count).ne.maxval(count)) then
    fmtstring = "('ERROR: Variable with name ', A, ' cannot be written.')"
    print fmtstring, trim(varname)
    print '(A)', "Argument ''count'' must be 1 at all but one dimension. "
    print *, 'count: ', count
    print *, 'start: ', start
    stop
  end if

  ! Check whether variable size and values of count are compatible
  if (product(count).ne.size(values, 1)) then
    fmtstring = "('ERROR: Variable with name ', A, ' cannot be written.')"
    print fmtstring, trim(varname)
    print '(A)', "Value of ''count'' must fit size of variable to dump. "
    print *, 'count:          ', count
    print *, 'size(variable): ', start
    stop
  end if
                  
  call getvarid(ncid  = ncid,        &
                name  = varname,     &
                varid = variable_id, &
                noerr = .true.)

  ! if variable not found -> ERROR
  if (variable_id==-1) then
    fmtstring = "('ERROR: Variable with name ', A, ' does not exist.')"
    print fmtstring, trim(varname)
    print '(A)', "Arguments ''start'' or ''count'' are only allowed, if variable already exists"
    stop
  end if

  select case (size(start,1))
  case(2)
    call putvar_real2d(ncid   = ncid,        &
                       varid  = variable_id, &
                       values = reshape(values, count(1:2)),     &
                       start  = start,       &
                       count  = count)
  case(3)
    call putvar_real3d(ncid   = ncid,        &
                       varid  = variable_id, &
                       values = reshape(values, count(1:3)),     &
                       start  = start,       &
                       count  = count)
  end select


end subroutine nc_putvar_by_name_1d_into_nd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_putvar_by_name_1d_into_nd_int(ncid, varname, values, start, count)
!< Write data into a 3D variable, or create it, if necessary
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  integer, intent(in)                        :: values(:)
  integer, intent(in)                        :: start(:), count(:)
  character(len=80)                          :: fmtstring
  integer                                    :: variable_id 

  ! A rather preculiar way of determining that count is larger 1 on all but one dimensio
  ! We can't use the FORTRAN intrinsic 'count' here, well, because we have a variable
  ! of this name. sucks
  if (product(count).ne.maxval(count)) then
    fmtstring = "('ERROR: Variable with name ', A, ' cannot be written.')"
    print fmtstring, trim(varname)
    print '(A)', "Argument ''count'' must be 1 at all but one dimension. "
    print *, 'count: ', count
    print *, 'start: ', start
    stop
  end if

  ! Check whether variable size and values of count are compatible
  if (product(count).ne.size(values, 1)) then
    fmtstring = "('ERROR: Variable with name ', A, ' cannot be written.')"
    print fmtstring, trim(varname)
    print '(A)', "Value of ''count'' must fit size of variable to dump. "
    print *, 'count:          ', count
    print *, 'size(variable): ', start
    stop
  end if

  call getvarid(ncid  = ncid,        &
                name  = varname,     &
                varid = variable_id, &
                noerr = .true.)

  ! if variable not found -> ERROR
  if (variable_id==-1) then
    fmtstring = "('ERROR: Variable with name ', A, ' does not exist.')"
    print fmtstring, trim(varname)
    print '(A)', "Arguments ''start'' or ''count'' are only allowed, if variable already exists"
    stop
  end if

  select case (size(start,1))
  case(2)
    call putvar_int2d( ncid   = ncid,        &
                       varid  = variable_id, &
                       values = reshape(values, count(1:2)),     &
                       start  = start,       &
                       count  = count)
  case(3)
    call putvar_int3d( ncid   = ncid,        &
                       varid  = variable_id, &
                       values = reshape(values, count(1:3)),     &
                       start  = start,       &
                       count  = count)
  end select


end subroutine nc_putvar_by_name_1d_into_nd_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvar_int1d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid, start, count
   integer, intent(out)         :: values(:)
   integer                      :: xtype, ndims, status, dimsize
   integer                      :: dimid(10), status_read
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       call check(status)
   end if

   if (size(values).ne.count) then
       write(*,100) myrank, trim(varname), varid, ncid, size(values), count
       stop
   end if

   status_read = nf90_get_var(ncid   = ncid,           &
                              varid  = varid,          &
                              values = values,         &
                              start  = [start],        &
                              count  = [count] )

                      
   if (status_read.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)
       if (ndims.ne.1) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           call check(status)
       end if
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       status = nf90_inquire_dimension(ncid  = ncid,     &
                                       dimid = dimid(1), &
                                       name  = dimname,  &
                                       len   = dimsize )

       if (start + count - 1 > dimsize) then
           write(*,102) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
           call check(status)
       end if

       write(*,103) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
       print *, trim(nf90_strerror(status_read))
       stop
   
   elseif (verbose>1) then
       write(lu_out,200) myrank, real(count) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 1D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of one')
102 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A)
103 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Read', F10.3, ' MB from 1D variable in NCID', I7, ', with ID:', I7)
end subroutine getvar_int1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvar_real1d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid, start, count
   real, intent(out)            :: values(:)
   integer                      :: xtype, ndims, status, dimsize
   integer                      :: dimid(10), status_read
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       call check(status)
   end if

   if (size(values).ne.count) then
       write(*,100) myrank, trim(varname), varid, ncid, size(values), count
       stop
   end if

   status_read = nf90_get_var(ncid   = ncid,           &
                              varid  = varid,          &
                              values = values,         &
                              start  = [start],        &
                              count  = [count] )

                      
   if (status_read.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)
       if (ndims.ne.1) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           call check(status)
       end if
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       status = nf90_inquire_dimension(ncid  = ncid,     &
                                       dimid = dimid(1), &
                                       name  = dimname,  &
                                       len   = dimsize )

       if (start + count - 1 > dimsize) then
           write(*,102) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
           call check(status)
       end if

       write(*,103) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
       print *, trim(nf90_strerror(status))
   
   elseif (verbose>1) then
       write(lu_out,200) myrank, real(count) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 1D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of one')
102 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A)
103 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Read', F10.3, ' MB from 1D variable in NCID', I7, ', with ID:', I7)
end subroutine getvar_real1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvar_real1d_dble(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid, start, count
   real(kind=dp), intent(out)   :: values(:)
   integer                      :: xtype, ndims, status, dimsize
   integer                      :: dimid(10), status_read
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       call check(status)
   end if

   if (size(values).ne.count) then
       write(*,100) myrank, trim(varname), varid, ncid, size(values), count
       stop
   end if

   status_read = nf90_get_var(ncid   = ncid,           &
                              varid  = varid,          &
                              values = values,         &
                              start  = [start],        &
                              count  = [count] )

                      
   if (status_read.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)
       if (ndims.ne.1) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           call check(status)
       end if
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       status = nf90_inquire_dimension(ncid  = ncid,     &
                                       dimid = dimid(1), &
                                       name  = dimname,  &
                                       len   = dimsize )
       if (start + count - 1 > dimsize) then
           write(*,102) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
           call check(status)
       end if

       write(*,103) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
       print *, trim(nf90_strerror(status_read))
       stop
   
   elseif (verbose>1) then
       write(lu_out,200) myrank, real(count) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 1D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of one')
102 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A)
103 format('ERROR: CPU ', I4, ' could not read 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Read', F10.3, ' MB from 1D variable in NCID', I7, ', with ID:', I7)
end subroutine getvar_real1d_dble
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvar_real2d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(2), count(2)
   real, intent(out)            :: values(:,:)
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10), status_read
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       call check(status)
   end if

   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 2
       if (size(values,idim).ne.count(idim)) then
           write(*,100) myrank, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           call flush(lu_out)
           stop
       end if
   end do

   ! Read data from file
   status_read = nf90_get_var(ncid   = ncid,           &
                              varid  = varid,          &
                              values = values,         &
                              start  = start,          &
                              count  = count )

                      
   ! If an error has occurred, try to find a reason                  
   if (status_read.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims.ne.2) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           call check(status)
       end if

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       do idim = 1, 2
           status = nf90_inquire_dimension(ncid  = ncid,        &
                                           dimid = dimid(idim), &
                                           name  = dimname,     &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim 
               call check(status)
           end if

           ! Otherwise just dump as much information as possible and stop
           write(*,103) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)
           print *, trim(nf90_strerror(status_read))

       end do

       stop
   
   elseif (verbose>1) then
       ! Everything okay
       write(lu_out,200) myrank, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 2D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not read 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not read 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of two')
102 format('ERROR: CPU ', I4, ' could not read 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not read 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Read', F10.3, ' MB into 2D variable in NCID', I7, ', with ID:', I7)
end subroutine getvar_real2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvar_int2d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(2), count(2)
   integer, intent(inout), allocatable :: values(:,:)
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10), status_read
   character(len=nf90_max_name) :: varname, dimname

   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   end if

   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 2
       if (size(values,idim).ne.count(idim)) then
           write(*,100) myrank, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           stop
       end if
   end do

   ! Read data from file
   status_read = nf90_get_var(ncid   = ncid,           &
                              varid  = varid,          &
                              values = values,         &
                              start  = start,          &
                              count  = count )

   ! If an error has occurred, try to find a reason                  
   if (status_read.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims.ne.2) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       end if

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       do idim = 1, 2
           status = nf90_inquire_dimension(ncid  = ncid,        &
                                           dimid = dimid(idim), &
                                           name  = dimname,     &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim 
               print *, trim(nf90_strerror(status_read))
               stop 2
           end if

           ! Otherwise just dump as much information as possible and stop
           write(*,103) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)

       end do

       print *, trim(nf90_strerror(status_read))
       stop
   
   elseif (verbose>1) then
       ! Everything okay
       write(lu_out,200) myrank, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 2D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not read 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not read 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of two')
102 format('ERROR: CPU ', I4, ' could not read 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not read 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Read', F10.3, ' MB into 2D variable in NCID', I7, ', with ID:', I7)
end subroutine getvar_int2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvar_int3d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(3), count(3)
   integer, intent(out)         :: values(count(1), count(2), count(3))
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10), status_read
   character(len=nf90_max_name) :: varname, dimname

   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       print *, trim(nf90_strerror(status))
       call flush(6)
       stop
   end if
   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 3
       if (size(values,idim).ne.count(idim)) then
           write(*,100) myrank, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           call flush(6)
           stop
       end if
   end do

   ! Read data from file
   status_read = nf90_get_var(ncid   = ncid,           &
                              varid  = varid,          &
                              values = values,         &
                              start  = start,          &
                              count  = count )

                      
   ! If an error has occurred, try to find a reason                  
   if (status_read.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims.ne.3) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           call flush(6)
           stop
       end if

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       do idim = 1, 3
           status = nf90_inquire_dimension(ncid  = ncid,        &
                                           dimid = dimid(idim), &
                                           name  = dimname,     &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim 
               call flush(6)
               print *, trim(nf90_strerror(status))
               stop
           end if

           ! Otherwise just dump as much information as possible and stop
           write(*,103) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)
           print *, trim(nf90_strerror(status_read))
           call flush(6)

       end do

       stop
   
   elseif (verbose>1) then
       ! Everything okay
       write(lu_out,200) myrank, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 3D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not read 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not read 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of three')
102 format('ERROR: CPU ', I4, ' could not read 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not read 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Read', F10.3, ' MB from 3D variable in NCID', I7, ', with ID:', I7)
end subroutine getvar_int3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvar_real3d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(3), count(3)
   real, intent(out)            :: values(:,:,:)
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10), status_read
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   end if
   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 3
       if (size(values,idim).ne.count(idim)) then
           write(*,100) myrank, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           stop
       end if
   end do

   ! Read data from file
   status_read = nf90_get_var(ncid   = ncid,      &
                         varid  = varid,          &
                         values = values,         &
                         start  = start,          &
                         count  = count )

                      
   ! If an error has occurred, try to find a reason                  
   if (status_read.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims.ne.3) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       end if

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       do idim = 1, 3
           status = nf90_inquire_dimension(ncid  = ncid,        &
                                           dimid = dimid(idim), &
                                           name  = dimname,     &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim 
               print *, trim(nf90_strerror(status))
               stop
           end if

           ! Otherwise just dump as much information as possible and stop
           write(*,103) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)
           print *, trim(nf90_strerror(status_read))

       end do

       stop
   
   elseif (verbose>1) then
       ! Everything okay
       write(lu_out,200) myrank, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 3D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not read 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not read 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of three')
102 format('ERROR: CPU ', I4, ' could not read 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not read 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Read', F10.3, ' MB from 3D variable in NCID', I7, ', with ID:', I7)
end subroutine getvar_real3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine putvar_real1d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid, start, count
   real, intent(in)             :: values(:)
   integer                      :: xtype, ndims, status, dimsize
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   end if

   if (size(values).ne.count) then
       write(*,100) myrank, trim(varname), varid, ncid, size(values), count
       stop
   end if

   status = nf90_put_var(ncid   = ncid,           &
                         varid  = varid,          &
                         values = values,         &
                         start  = [start],        &
                         count  = [count] )

                      
   if (status.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)
       if (ndims.ne.1) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       end if
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       status = nf90_inquire_dimension(ncid  = ncid,     &
                                       dimid = dimid(1), &
                                       name  = dimname,  &
                                       len   = dimsize )
       if (start + count - 1 > dimsize) then
           write(*,102) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))
           stop
       end if

       write(*,103) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
       print *, trim(nf90_strerror(status))
       stop
   
   elseif (verbose>1) then
       write(lu_out,200) myrank, real(count) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 1D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of one')
102 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A)
103 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Wrote', F10.3, ' MB into 1D variable in NCID', I7, ', with ID:', I7)
end subroutine putvar_real1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine putvar_real2d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(2), count(2)
   real, intent(in)             :: values(:,:)
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   end if
   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 2
       if (size(values,idim).ne.count(idim)) then
           write(*,100) myrank, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           stop
       end if
   end do

   ! Write data to file
   status = nf90_put_var(ncid   = ncid,           &
                         varid  = varid,          &
                         values = values,         &
                         start  = start,          &
                         count  = count )

                      
   ! If an error has occurred, try to find a reason                  
   if (status.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims.ne.2) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       end if

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       do idim = 1, 2
           status = nf90_inquire_dimension(ncid  = ncid,        &
                                           dimid = dimid(idim), &
                                           name  = dimname,     &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim 
               print *, trim(nf90_strerror(status))
               stop
           end if

           ! Otherwise just dump as much information as possible and stop
           write(*,103) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))

       end do

       stop
   
   elseif (verbose>1) then
       ! Everything okay
       write(lu_out,200) myrank, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 2D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of two')
102 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Wrote', F10.3, ' MB into 2D variable in NCID', I7, ', with ID:', I7)
end subroutine putvar_real2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine putvar_real3d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(3), count(3)
   real, intent(in)             :: values(:,:,:)
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   end if
   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 3
       if (size(values,idim).ne.count(idim)) then
           write(*,100) myrank, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           stop
       end if
   end do

   ! Write data to file
   status = nf90_put_var(ncid   = ncid,           &
                         varid  = varid,          &
                         values = values,         &
                         start  = start,          &
                         count  = count )

                      
   ! If an error has occurred, try to find a reason                  
   if (status.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims.ne.3) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       end if

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       do idim = 1, 3
           status = nf90_inquire_dimension(ncid  = ncid,        &
                                           dimid = dimid(idim), &
                                           name  = dimname,     &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim 
               print *, trim(nf90_strerror(status))
               stop
           end if

           ! Otherwise just dump as much information as possible and stop
           write(*,103) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))

       end do

       stop
   
   elseif (verbose>1) then
       ! Everything okay
       write(lu_out,200) myrank, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 3D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of three')
102 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Wrote', F10.3, ' MB into 3D variable in NCID', I7, ', with ID:', I7)
end subroutine putvar_real3d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine putvar_int1d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid, start, count
   integer, intent(in)          :: values(:)
   integer                      :: xtype, ndims, status, dimsize
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   end if

   if (size(values).ne.count) then
       write(*,100) myrank, trim(varname), varid, ncid, size(values), count
       stop
   end if

   status = nf90_put_var(ncid   = ncid,           &
                         varid  = varid,          &
                         values = values,         &
                         start  = [start],        &
                         count  = [count] )

                      
   if (status.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)
       if (ndims.ne.1) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       end if
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       status = nf90_inquire_dimension(ncid  = ncid,     &
                                       dimid = dimid(1), &
                                       name  = dimname,  &
                                       len   = dimsize )
       if (start + count - 1 > dimsize) then
           write(*,102) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))
           stop
       end if

       write(*,103) myrank, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
       print *, trim(nf90_strerror(status))
       stop
   
   elseif (verbose>1) then
       write(lu_out,200) myrank, real(count) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 1D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of one')
102 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A)
103 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Wrote', F10.3, ' MB into 1D variable in NCID', I7, ', with ID:', I7)
end subroutine putvar_int1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine putvar_int2d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(2), count(2)
   integer, intent(in)          :: values(:,:)
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   end if
   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 2
       if (size(values,idim).ne.count(idim)) then
           write(*,100) myrank, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           stop
       end if
   end do

   ! Write data to file
   status = nf90_put_var(ncid   = ncid,           &
                         varid  = varid,          &
                         values = values,         &
                         start  = start,          &
                         count  = count )

                      
   ! If an error has occurred, try to find a reason                  
   if (status.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims.ne.2) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       end if

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       do idim = 1, 2
           status = nf90_inquire_dimension(ncid  = ncid,        &
                                           dimid = dimid(idim), &
                                           name  = dimname,     &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim 
               print *, trim(nf90_strerror(status))
               call pabort()
               stop
           end if

           ! Otherwise just dump as much information as possible and stop
           write(*,103) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))

       end do

       stop
   
   elseif (verbose>1) then
       ! Everything okay
       write(lu_out,200) myrank, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 2D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of two')
102 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Wrote', F10.3, ' MB into 2D variable in NCID', I7, ', with ID:', I7)
end subroutine putvar_int2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine putvar_int3d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(3), count(3)
   integer, intent(in)          :: values(:,:,:)
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid,     &
                                  varid = varid,    &
                                  name  = varname )

   if (status.ne.NF90_NOERR) then
       write(*,99) myrank, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   end if
   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 3
       if (size(values,idim).ne.count(idim)) then
           write(*,100) myrank, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           stop
       end if
   end do

   ! Write data to file
   status = nf90_put_var(ncid   = ncid,           &
                         varid  = varid,          &
                         values = values,         &
                         start  = start,          &
                         count  = count )

                      
   ! If an error has occurred, try to find a reason                  
   if (status.ne.NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid,    &
                                      varid = varid,    &
                                      name  = varname,  &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims.ne.3) then
           write(*,101) myrank, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       end if

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid,     &
                                      varid  = varid,    &
                                      name   = varname,  &
                                      xtype  = xtype,    &
                                      ndims  = ndims,    &
                                      dimids = dimid  )

       do idim = 1, 3
           status = nf90_inquire_dimension(ncid  = ncid,        &
                                           dimid = dimid(idim), &
                                           name  = dimname,     &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim 
               print *, trim(nf90_strerror(status))
               stop
           end if

           ! Otherwise just dump as much information as possible and stop
           write(*,103) myrank, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))

       end do

       stop
   
   elseif (verbose>1) then
       ! Everything okay
       write(lu_out,200) myrank, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(lu_out)
   end if
    
99  format('ERROR: CPU ', I4, ' could not find 3D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of three')
102 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Wrote', F10.3, ' MB into 3D variable in NCID', I7, ', with ID:', I7)
end subroutine putvar_int3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_create_var_by_name(ncid, varname, sizes, dimension_names)
!< Create variable by its size without putting data into it
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  integer, intent(in)                        :: sizes(:)
  character(*), intent(in), optional         :: dimension_names(:)
  character(len=nf90_max_name), allocatable  :: dimensions(:)
  integer                                    :: ndim
  character(len=nf90_max_name)               :: dimname
  character(len=80)                          :: fmtstring
  integer                                    :: i, variable_id, status, ncid_root
  integer, allocatable                       :: dimid(:)

  call getvarid(ncid  = ncid,        &
                name  = varname,     &
                varid = variable_id, &
                noerr = .true.)

  ndim = size(sizes)
  allocate(dimensions(ndim))

  ! Get ncid of root group. Used for dimensions, which should always be defined
  ! in the root group (at least according to my taste)
  status = 0
  ncid_root = ncid
  do while (status.eq.NF90_NOERR) !Check whether we are already in root group
    status = nf90_inq_grp_parent(ncid_root, ncid_root) 
  end do

  if (present(dimension_names)) then
    if (size(dimension_names,1).ne.ndim) then
      write(*,*) 'ERROR: Each dimension must have a name in nc_create_var_by_name'
      write(*,*) '       To fill dimension names automatically, omit dimension_names'
      write(*,*) '       or set the elements to blank'
      stop
    end if

    dimensions = dimension_names
  end if

  allocate(dimid(ndim))

  ! if variable not found              
  if (variable_id==-1) then
    call check(nf90_redef(ncid = ncid_root))
    do i = 1, ndim
      if (trim(dimensions(i)) == '') then
        write(dimname, '(A,"_",I1)') trim(varname), i
        status = nf90_def_dim(ncid  = ncid_root,       &
                              name  = trim(dimname),   &
                              len   = sizes(i),        &
                              dimid = dimid(i))
      else
        status = nf90_inq_dimid(ncid = ncid_root,       &
                                name = dimensions(i),   &
                                dimid = dimid(i))
        if (status.ne.NF90_NOERR) then
          call check(nf90_def_dim(ncid  = ncid_root,           &
                                  name  = trim(dimensions(i)), &
                                  len   = sizes(i),            &
                                  dimid = dimid(i)) )
        end if

      end if
    end do

    status = nf90_def_var( ncid   = ncid,          &
                           name   = trim(varname), &
                           xtype  = NF90_FLOAT,    &
                           dimids = dimid,         &
                           varid  = variable_id) 
    call check(nf90_enddef(ncid = ncid_root))

  else
    fmtstring = "('ERROR: Variable with name ', A, ' cannot be created.')"
    print fmtstring, trim(varname)
    print '(A)', "it already exists"
    stop
  end if

end subroutine nc_create_var_by_name
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Translates netcdf error codes into error messages
subroutine check(status)
  implicit none
  integer, intent (in) :: status

  if(status /= nf90_noerr) then 
     print *, '************************************************************************'
     print *, 'Problem with NetCDF on node', myrank
     print *, trim(nf90_strerror(status))
     print *, '************************************************************************'
     call flush(6)
     call pabort()
  end if
end subroutine check  
!-----------------------------------------------------------------------------------------

end module nc_routines
