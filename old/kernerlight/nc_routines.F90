module nc_routines
    use global_parameters
    use commpi
    use parameter_types, only: ncparamtype
    use netcdf
    implicit none
    public :: open_netcdf,  close_netcdf, nc_read_summedstraintrace
    public :: nc_read_att_int, nc_read_att_real, nc_read_att_char
    public :: check, getvarid
    private
contains
!-----------------------------------------------------------------------------------------
subroutine close_netcdf(io)

    use parameter_types, only: ioparamtype
    type(ioparamtype), intent(inout)  :: io
    integer                           :: isim
    do isim = 1, io%nsim_fwd
        call check( nf90_close(ncid=io%fwd(isim)%ncid) )
    end do
    do isim = 1, io%nsim_bwd
        call check( nf90_close(ncid=io%bwd(isim)%ncid) )
    end do

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine open_netcdf(io)

    use parameter_types, only: ioparamtype
    type(ioparamtype), intent(inout) :: io
    integer                          :: status, isim
    character(len=200)               :: format20, format21, filename

    format20 = "('Trying to open NetCDF file ', A, ' on CPU ', I5)"
    format21 = "('Succeded,  has NCID ', I6, ', Snapshots group NCID:', I6)"

    do isim = 1, io%nsim_fwd
        ! Forward wavefield
        filename=trim(io%fwd(isim)%meshdir)//'ordered_output.nc4'
        inquire(file=filename, exist=io%fwd(isim)%ordered_output)
        if (.not.io%fwd(isim)%ordered_output) then
            filename = trim(io%fwd(isim)%meshdir)//'Data/axisem_output.nc4'
        end if
        write(6,format20) filename, mynum
        call check( nf90_open(path=filename, &
                              mode=nf90_nowrite, ncid=io%fwd(isim)%ncid) )
        call check(nf90_inq_ncid(io%fwd(isim)%ncid, name="Snapshots", grp_ncid=io%fwd(isim)%snap))
        call check(nf90_inq_ncid(io%fwd(isim)%ncid, name="Surface", grp_ncid=io%fwd(isim)%surf))
        call check(nf90_inq_ncid(io%fwd(isim)%ncid, name="Mesh", grp_ncid=io%fwd(isim)%mesh))
        call nc_read_att_int(io%fwd(isim)%ndumps, 'number of strain dumps', io%fwd(isim))
        write(6,format21) io%fwd(isim)%ncid, io%fwd(isim)%snap 
    end do
        
    ! Backward wavefield

    do isim = 1, io%nsim_bwd
        filename = trim(io%bwd(isim)%meshdir)//'Data/axisem_output.nc4'
        write(6,format20) filename, mynum
        call check( nf90_open(path=filename, &
                              mode=nf90_nowrite, ncid=io%bwd(isim)%ncid) )
        call check(nf90_inq_ncid(io%bwd(isim)%ncid, name="Surface", grp_ncid=io%bwd(isim)%surf))
        call check(nf90_inq_ncid(io%bwd(isim)%ncid, name="Snapshots", grp_ncid=io%bwd(isim)%snap))
        call check(nf90_inq_ncid(io%bwd(isim)%ncid, name="Mesh", grp_ncid=io%bwd(isim)%mesh))
        call nc_read_att_int(io%bwd(isim)%ndumps, 'number of strain dumps', io%bwd(isim))
        write(6,format21) io%bwd(isim)%ncid, io%bwd(isim)%snap 
    end do

    call flush(6) 

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine check(status)
! Translates netcdf error codes into error messages

  implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop 2
  end if
end subroutine check  
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
        write(6,100) mynum, trim(name), ncid
        stop
    elseif (mynum.eq.0) then
        write(6,101) trim(name), ncid, varid
    end if
100 format('ERROR: CPU ', I4, ' could not find variable: ''', A, ''' in NCID', I7)
101 format('Variable ''', A, ''' found in NCID', I7, ', has ID:', I7)
    end subroutine

!-----------------------------------------------------------------------------------------
subroutine nc_read_summedstraintrace(u, azim_factor, ipt_mesh, nc)
    type(ncparamtype), intent(in)                :: nc
    integer, dimension(:), allocatable           :: ipt_mesh
    real(kind=8), dimension(:,:), intent(inout)  :: u
    real, dimension(:), intent(in)               :: azim_factor
    !include 'mesh_params_kernel.h'
    real, dimension(nc%ndumps)                   :: utemp
    integer                                      :: ncvarid_straintrace, ncid_snap
    integer                                      :: ipt, npts


    npts = size(u,2)
    call check( nf90_inq_varid( ncid  = nc%snap, &
                                name  = "straintrace", &
                                varid = ncvarid_straintrace) )
    do ipt = 1, npts

        call check( nf90_get_var( ncid   = nc%snap,  &
                                  varid  = ncvarid_straintrace, &
                                  start  = (/1,ipt_mesh(ipt)/),  &
                                  count  = (/nc%ndumps,1/), &
                                  values = utemp )) 
        u(:,ipt) = u(:,ipt) + dble(azim_factor(ipt) * utemp)
    end do
end subroutine


!-----------------------------------------------------------------------------------------

subroutine nc_read_surfcoords(nsurf, surftheta, nc)
    integer, intent(out)    :: nsurf
    type(ncparamtype), intent(in) :: nc
    real, allocatable      :: surftheta(:)
    integer                :: ncvarid_theta, nc_surf_dimid

    call check( nf90_inq_varid(ncid=nc%surf, name="elem_theta", &
                              varid=ncvarid_theta) )
    call check( nf90_get_att(ncid = nc%surf, &
                             name = 'nsurfelem', &
                             varid = NF90_GLOBAL, &
                             values = nsurf) ) 

    call check( nf90_get_var(ncid  = nc%surf, &
                             varid = ncvarid_theta, &
                             values = surftheta )) 
end subroutine

!-----------------------------------------------------------------------------------------
subroutine nc_read_surf(recel, u, v, nc, srctype)
    type(ncparamtype), intent(in) :: nc
    integer, intent(in)                    :: recel
    real, intent(out), dimension(nc%ndumps,3) :: u, v
    character(len=10)                      :: srctype
    integer                                :: ncvarid_surfdisp, ncvarid_surfvelo
    call check( nf90_inq_varid(ncid=nc%surf, name="displacement", &
                               varid=ncvarid_surfdisp) )
    call check( nf90_inq_varid(ncid=nc%surf, name="velocity", &
                               varid=ncvarid_surfvelo) )
    if (srctype=='monopole') then
        call check( nf90_get_var(ncid=nc%surf, varid=ncvarid_surfdisp, &
                                 start = (/1,1,recel/), &
                                 count=(/nc%ndumps,2,1/), &
                                 values = u(:,1:2) )) 
        call check( nf90_get_var(ncid=nc%surf, varid=ncvarid_surfvelo, &
                                 start = (/1,1,recel/), &
                                 count=(/nc%ndumps,2,1/), &
                                 values = v(:,1:2) )) 
    else
        call check( nf90_get_var(ncid=nc%surf, varid=ncvarid_surfdisp, &
                                 start = (/1,1,recel/), &
                                 count=(/nc%ndumps,3,1/), &
                                 values = u )) 
        call check( nf90_get_var(ncid=nc%surf, varid=ncvarid_surfvelo, &
                                 start = (/1,1,recel/), &
                                 count=(/nc%ndumps,3,1/), &
                                 values = v )) 
    end if

end subroutine nc_read_surf

!-----------------------------------------------------------------------------------------
!> Read NetCDF attribute of type Integer
subroutine nc_read_att_int(attribute_value, attribute_name, nc)
  character(len=*),  intent(in)     :: attribute_name
  integer, intent(out)              :: attribute_value
  type(ncparamtype), intent(in)     :: nc
  integer                           :: status

  status = nf90_get_att(nc%ncid, NF90_GLOBAL, attribute_name, attribute_value)
  if (status.ne.NF90_NOERR) then
      write(6,*) 'Could not find attribute ', attribute_name, ' in NetCDF file ', &
                 nc%meshdir, '/Data/axisem_output.nc4'
      stop
  end if
end subroutine nc_read_att_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Read NetCDF attribute of type Character
subroutine nc_read_att_char(attribute_value, attribute_name, nc)
  character(len=*),  intent(in)     :: attribute_name
  character(len=*), intent(out)     :: attribute_value
  type(ncparamtype), intent(in)     :: nc
  integer                           :: status

  status = nf90_get_att(nc%ncid, NF90_GLOBAL, attribute_name, attribute_value)
  if (status.ne.NF90_NOERR) then
      write(6,*) 'Could not find attribute ', attribute_name, ' in NetCDF file ', &
                 nc%meshdir, '/Data/axisem_output.nc4'
      stop
  end if
end subroutine nc_read_att_char
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Read NetCDF attribute of type Real
subroutine nc_read_att_real(attribute_value, attribute_name, nc)
  character(len=*),  intent(in)     :: attribute_name
  real, intent(out)                 :: attribute_value
  type(ncparamtype), intent(in)     :: nc
  integer                           :: status

  status = nf90_get_att(nc%ncid, NF90_GLOBAL, attribute_name, attribute_value)
  if (status.ne.NF90_NOERR) then
      write(6,*) 'Could not find attribute ', attribute_name, ' in NetCDF file ', &
                 nc%meshdir, '/Data/axisem_output.nc4'
      stop
  end if
end subroutine nc_read_att_real
!-----------------------------------------------------------------------------------------


end module
