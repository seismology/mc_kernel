module nc_routines
    use global_parameters, only: sp, dp, lu_out, verbose, myrank
    use commpi,            only: pabort    
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

    interface nc_getvar
      module procedure  :: getvar_real1d
      module procedure  :: getvar_real1d_dble
      module procedure  :: getvar_int1d
      module procedure  :: getvar_real2d
      module procedure  :: getvar_int2d
      module procedure  :: getvar_real3d
      module procedure  :: getvar_int3d
      !module procedure  :: getvar_real1d_from_3d
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
                      mode     = nf90_nowrite,          &
                      ncid     = ncid)

   if (status.ne.nf90_noerr) then
      fmtstring = "('ERROR: CPU ', I4, ' tried to open file ''', A, ''', " &
                    // "but could not find it')"
      print fmtstring, myrank, trim(filename)
      stop
   end if

end subroutine nc_open_for_read
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvarid(ncid, name, varid)
    integer, intent(in)          :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out)         :: varid
    integer                      :: status

    status = nf90_inq_varid( ncid  = ncid, &
                             name  = name, &
                             varid = varid )
    if (status.ne.NF90_NOERR) then
        write(6,100) myrank, trim(name), ncid
        call pabort 
    elseif (verbose>1) then
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
subroutine nc_getvar_by_name_1d_int(ncid, name, values)
!< Looks up a 1D variable name and returns the complete variable
  integer, intent(in)           :: ncid
  character(len=*), intent(in)  :: name
  integer, allocatable, intent(out)  :: values(:)

  integer                       :: variable_id, dimid(1), npoints, variable_type
  integer                       :: status

  if (verbose>1) then
    write(lu_out,"(' Trying to read 1D variable ', A, '...')") trim(name)
  end if

  call  getvarid( ncid  = ncid,            &
                  name  = name,   &
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
    write(*,*) 'get_mesh_variable is only implemented for NF90_INT variables'
    call pabort()
  end select

end subroutine nc_getvar_by_name_1d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_1d(ncid, name, values, limits)
!< Looks up a 1D variable name and returns the complete variable
  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: name
  real(kind=sp), allocatable, intent(inout)  :: values(:)
  real(kind=sp), intent(in), optional        :: limits(2)

  integer                        :: variable_id, dimid(1), npoints, variable_type
  integer                        :: status
  real(kind=sp)                  :: limits_loc(2)
  logical                        :: have_limits = .false.


  if (verbose>1) then
    write(lu_out,"(' Trying to read 1D variable ', A, '...')") trim(name)
  end if

  call  getvarid( ncid  = ncid,            &
                  name  = name,   &
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
    write(*,*) 'get_mesh_variable is only implemented for NF90_FLOAT and NF90_DBLE variables'
    call pabort()
  end select

  have_limits = .false.
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
        write(*,*) 'Variable name: ', trim(name) 
        write(*,*) 'Element nr   : ', minloc(values)
        write(*,*) 'Element value: ', minval(values)
        write(*,*) 'Lower limit  : ', minval(limits_loc)
     end if
     if (maxval(values).gt.(maxval(limits_loc))) then
        write(*,*) 'ERROR: Value in NetCDF file larger than limit!'
        write(*,*) 'Variable name: ', trim(name) 
        write(*,*) 'Element nr   : ', maxloc(values)
        write(*,*) 'Element value: ', maxval(values)
        write(*,*) 'Upper limit  : ', maxval(limits_loc)
     end if
  end if

end subroutine nc_getvar_by_name_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_2d_int(ncid, name, values)
!< Looks up a 2D variable name and returns the complete variable
  integer, intent(in)           :: ncid
  character(len=*), intent(in)  :: name
  integer, allocatable, intent(out)  :: values(:,:)

  integer                       :: variable_id, dimid(2), len_dim1, len_dim2, variable_type
  integer                       :: status

  if (verbose>1) then
    write(lu_out,"(' Trying to read 2D variable ', A, '...')") trim(name)
  end if

  call  getvarid( ncid  = ncid,            &
                  name  = name,   &
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
    write(*,*) 'get_mesh_variable is only implemented for NF90_INT variables'
    call pabort()
  end select

end subroutine nc_getvar_by_name_2d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_2d(ncid, name, values)
!< Looks up a 2D variable name and returns the complete variable
  integer, intent(in)           :: ncid
  character(len=*), intent(in)  :: name
  real(kind=sp), allocatable, intent(inout)  :: values(:,:)

  integer                       :: variable_id, dimid(2), len_dim1, len_dim2, variable_type
  integer                       :: status

  if (verbose>1) then
    write(lu_out,"(' Trying to read 2D variable ', A, '...')") trim(name)
  end if

  call  getvarid( ncid  = ncid,            &
                  name  = name,   &
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
    write(*,*) 'get_mesh_variable is only implemented for NF90_FLOAT variables'
    call pabort()
  end select

end subroutine nc_getvar_by_name_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_3d_int(ncid, name, values)
!< Looks up a 3D variable name and returns the complete variable
  integer, intent(in)                :: ncid
  character(len=*), intent(in)       :: name
  integer, allocatable, intent(out)  :: values(:,:,:)

  integer                       :: variable_id, dimid(3), len_dim(3), variable_type
  integer                       :: status, i_dim

  if (verbose>1) then
    write(lu_out,"(' Trying to read 3D variable ', A, '...')") trim(name)
  end if

  call  getvarid( ncid  = ncid,            &
                  name  = name,   &
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
    write(*,*) 'get_mesh_variable is only implemented for NF90_INT variables'
    call pabort()
  end select

end subroutine nc_getvar_by_name_3d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_getvar_by_name_3d(ncid, name, values)
!< Looks up a 3D variable name and returns the complete variable
  integer, intent(in)                      :: ncid
  character(len=*), intent(in)             :: name
  real(kind=sp), allocatable, intent(out)  :: values(:,:,:)

  integer                       :: variable_id, dimid(3), len_dim(3), variable_type
  integer                       :: status, i_dim

  if (verbose>1) then
    write(lu_out,"(' Trying to read 3D variable ', A, '...')") trim(name)
  end if

  call  getvarid( ncid  = ncid,            &
                  name  = name,   &
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
    write(*,*) 'get_mesh_variable is only implemented for NF90_FLOAT variables'
    call pabort()
  end select

end subroutine nc_getvar_by_name_3d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine getvar_int1d(ncid, varid, values, start, count)
!< Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid, start, count
   integer, intent(out)         :: values(:)
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

   status = nf90_get_var(ncid   = ncid,           &
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

   status = nf90_get_var(ncid   = ncid,           &
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

   status = nf90_get_var(ncid   = ncid,           &
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
   status = nf90_get_var(ncid   = ncid,           &
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
               stop 2
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
   status = nf90_get_var(ncid   = ncid,           &
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
               stop 2
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
   integer, intent(out)         :: values(:,:,:)
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
   status = nf90_get_var(ncid   = ncid,           &
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
   status = nf90_get_var(ncid   = ncid,           &
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
!> Translates netcdf error codes into error messages
subroutine check(status)
  implicit none
  integer, intent (in) :: status

  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     call pabort
     !call tracebackqq()
  end if
end subroutine check  
!-----------------------------------------------------------------------------------------

end module
