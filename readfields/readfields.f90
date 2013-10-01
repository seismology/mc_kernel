module readfields
    use kdtree2_module
    !use netcdf
    include 'netcdf.h'

    integer, parameter                    :: sp = 4, dp = 8

    ! Should be moved to somewhere else later

    type meshtype
        real(kind=sp), allocatable        :: s(:), z(:)
        integer                           :: npoints
    end type

    type ncparamtype
        integer                           :: ncid
        integer                           :: snap, surf, mesh  ! Group IDs
        integer                           :: straintrace, dev_strain ! Variable IDs
        character(len=200)                :: meshdir
        integer                           :: ndumps
        logical                           :: ordered_output
    end type

    type netcdf_type
        private

        integer                           :: nsim_fwd, nsim_bwd
        type(ncparamtype), allocatable    :: fwd(:)
        type(ncparamtype), allocatable    :: bwd(:)

        type(kdtree2), pointer            :: fwdtree, bwdtree
        type(kdtree2_result), allocatable :: nextpoint(:)

        type(meshtype)                    :: fwdmesh, bwdmesh

        logical                           :: params_set
        logical                           :: files_open
        logical                           :: meshes_read
        logical                           :: kdtree_built

        contains 
            procedure, pass               :: set_params
            procedure, pass               :: open_files
            procedure, pass               :: read_meshes
            procedure, pass               :: build_kdtree
            procedure, pass               :: load_fw_points
            procedure, pass               :: load_bw_points

    end type

contains

subroutine set_params(this, parameters)
    class(netcdf_type)            :: this
    class(parameter_type)         :: parameters


    this%nsim_fwd = parameters%nsim_fwd
    this%nsim_bwd =  1 !parameters%nsim_bwd
    
    allocate( this%fwd(this%nsim_fwd) )
    allocate( this%bwd(this%nsim_bwd) )

    select case(this%nsim_fwd)
    case(1)
        this%fwd(1)%meshdir = parameters%dir_fwdmesh

    case(2) 
        this%fwd(1)%meshdir = trim(parameters%fwd_dir)//'PZ/'
        this%fwd(2)%meshdir = trim(parameters%fwd_dir)//'PX/'

    case(4)
        this%fwd(1)%meshdir = trim(parameters%fwd_dir)//'MZZ/'
        this%fwd(2)%meshdir = trim(parameters%fwd_dir)//'MXX_P_MYY/'
        this%fwd(3)%meshdir = trim(parameters%fwd_dir)//'MXZ_MYZ/'
        this%fwd(4)%meshdir = trim(parameters%fwd_dir)//'MXY_MXX_M_MYY/'
    end select
    
    do isim = 1, this%nsim_bwd
        this%bwd(isim)%meshdir = parameters%dir_bwdmesh
    end do

    this%params_set = .true.

end subroutine


!-------------------------------------------------------------------------------
subroutine open_files(this)

    class(netcdf_type)               :: this
    integer                          :: status, isim
    character(len=200)               :: format20, format21, filename

    if (.not.this%params_set) then
        print *, 'ERROR in open_files(): Parameters have to be set first'
        print *, 'Call set_params before open_files()'
        stop
    end if

    format20 = "('Trying to open NetCDF file ', A, ' on CPU ', I5)"
    format21 = "('Succeded,  has NCID ', I6, ', Snapshots group NCID:', I6)"

    do isim = 1, this%nsim_fwd
        ! Forward wavefield
        filename=trim(this%fwd(isim)%meshdir)//'ordered_output.nc4'
        
        inquire(file=filename, exist=this%fwd(isim)%ordered_output)
        if (.not.this%fwd(isim)%ordered_output) then
            filename = trim(this%fwd(isim)%meshdir)//'Data/axisem_output.nc4'
        end if
        
        write(6,format20) filename, mynum
        status = nf90_open(       path     = filename,              &
                                  mode     = nf90_nowrite,          &
                                  ncid     = this%fwd(isim)%ncid)

        call check(nf90_inq_ncid( ncid     = this%fwd(isim)%ncid,   &
                                  name     = "Snapshots",           &
                                  grp_ncid = this%fwd(isim)%snap))
        call check(nf90_inq_varid(ncid     = this%fwd(isim)%snap,   &
                                  name     = "straintrace",         &
                                  varid    = this%fwd(isim)%straintrace) )
        
        call check(nf90_inq_ncid(ncid      = this%fwd(isim)%ncid,   &
                                 name      = "Surface",             &
                                 grp_ncid  = this%fwd(isim)%surf))

        call check(nf90_inq_ncid(ncid      = this%fwd(isim)%ncid,   &
                                 name      = "Mesh",                &
                                 grp_ncid  = this%fwd(isim)%mesh))

        call nc_read_att_int(    this%fwd(isim)%ndumps,             &
                                 'number of strain dumps',          &
                                 this%fwd(isim))
        write(6,format21) this%fwd(isim)%ncid, this%fwd(isim)%snap 
    end do
        

    do isim = 1, this%nsim_bwd
        ! Backward wavefield
        filename=trim(this%bwd(isim)%meshdir)//'ordered_output.nc4'
        
        inquire(file=filename, exist=this%bwd(isim)%ordered_output)
        if (.not.this%bwd(isim)%ordered_output) then
            filename = trim(this%bwd(isim)%meshdir)//'Data/axisem_output.nc4'
        end if

        write(6,format20) filename, mynum
        call check( nf90_open(    path     = filename,              &
                                  mode     = nf90_nowrite,          &
                                  ncid     = this%bwd(isim)%ncid) )

        call check(nf90_inq_ncid( ncid     = this%bwd(isim)%ncid,   &
                                  name     = "Snapshots",           &
                                  grp_ncid = this%bwd(isim)%snap))

        call check(nf90_inq_varid(ncid     = this%bwd(isim)%snap,   &
                                  name     = "straintrace",         &
                                  varid    = this%bwd(isim)%straintrace) )
        
        call check(nf90_inq_ncid( ncid     = this%bwd(isim)%ncid,   &
                                  name     = "Surface",             &
                                  grp_ncid = this%bwd(isim)%surf))
        call check(nf90_inq_ncid( ncid     = this%bwd(isim)%ncid,   &
                                  name     = "Mesh",                &
                                  grp_ncid = this%bwd(isim)%mesh))
        call nc_read_att_int(     this%bwd(isim)%ndumps,            &
                                  'number of strain dumps',         &
                                  this%bwd(isim))
        write(6,format21) this%bwd(isim)%ncid, this%bwd(isim)%snap 
    end do

    call flush(6) 
    this%files_open = .true.

end subroutine


!-------------------------------------------------------------------------------
function load_fw_points(this, coordinates, sourceparams)
    class(netcdf_type)            :: this
    real(kind=dp), intent(in)     :: coordinates(:,:)
    type(src_param_type)          :: source_params
    real(kind=sp), allocatable    :: load_fw_points(:,:)
    real(kind=sp)                 :: utemp(this%fwd(1)%ndumps)
    type(kdtree2_result), pointer :: nextpoint(:)

    integer                       :: npoints, pointid(size(coordinates,2))
    
    npoints = size(coordinates,2)

    allocate(load_fw_points(this%fwd(1)%ndumps, npoints))
    load_fw_points(:,:) = 0.0
    
    ! Rotate points to FWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,           &
                          coordinates(1,:), coordinates(2,:), coordinates(3,:), &
                          phi, theta)


    do ipoint = 1, npoints
        
        call kdtree2_n_nearest( this%fwdtree,                           &
                                [rotmesh_s(ipoint), rotmesh_z(ipoint)], &
                                nn = 1,                                 &
                                results = nextpoint )
        
        pointid(ipoint) = nextpoint(1)%idx

        do isim = 1, this%nsim_fwd
            call check( nf90_get_var( ncid   = this%fwd(isim)%snap,        & 
                                      varid  = this%fwd(isim)%straintrace, &
                                      start  = [1, pointid(ipoint)],       &
                                      count  = [this%fwd%ndumps, 1],       &
                                      values = utemp) )
            load_fw_points(:, ipoint) = load_fw_points(:,ipoint) &
                                        + azim_factor(rotmesh_phi(ipoint), &
                                                      source_params%mij, isim) &
                                        * utemp
        end do !isim
        
    end do !ipoint

end function load_fw_points

!-------------------------------------------------------------------------------
function load_bw_points(this, coordinates, receiver)
    class(netcdf_type)            :: this
    real(kind=dp), intent(in)     :: coordinates(:,:)
    type(rec_param_type)          :: receiver
    real(kind=sp), allocatable    :: load_bw_points(:,:)
    real(kind=sp)                 :: utemp(this%bwd(1)%ndumps)
    type(kdtree2_result), pointer :: nextpoint(:)

    integer                       :: npoints , pointid(size(coordinates,2))

    npoints = size(coordinates,2)

    allocate(load_bw_points(this%bwd(1)%ndumps,npoints))
    load_bw_points(:,:) = 0.0
    
    ! Rotate points to BWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,           &
                          coordinates(1,:), coordinates(2,:), coordinates(3,:), &
                          phi, theta)


    do ipoint = 1, npoints
        
        call kdtree2_n_nearest( this%bwdtree,                           &
                                [rotmesh_s(ipoint), rotmesh_z(ipoint)], &
                                nn = 1,                                 &
                                results = nextpoint )
        
        pointid(ipoint) = nextpoint(1)%idx

        call check( nf90_get_var( ncid   = this%bwd(1)%snap,           & 
                                  varid  = this%bwd(1)%straintrace,    &
                                  start  = [1, pointid(ipoint)],       &
                                  count  = [this%bwd%ndumps, 1],       &
                                  values = utemp) )

        select case(receiver%component)
        case('Z')
            load_bw_points(:,ipoint) = utemp
        case('R')
            load_bw_points(:,ipoint) = cos(rotmesh_phi) * utemp
        case('T')
            load_bw_points(:,ipoint) = - sin(rotmesh_phi) * utemp 
        end select

    end do !ipoint

end function load_bw_points

!-------------------------------------------------------------------------------
function azim_factor(phi, mij, isim)
    real(kind=sp), intent(in)    :: phi
    real(kind=sp), intent(in)    :: mij(6)
    integer, intent(in)          :: isim
    real(kind=sp)                :: azim_factor


    select case(isim)
    case(1) ! Mzz
       azim_factor = Mij(1)
       
    case(2) ! Mxx+Myy
       azim_factor = 0.5*(Mij(2)+Mij(3))
       
    case(3) ! dipole
       azim_factor = Mij(4)*cos(phi) + Mij(5)*sin(phi)
       
    case(4) ! quadrupole
       azim_factor = 0.5*(Mij(2)-Mij(3))*cos(2.d0*phi) + Mij(6)*sin(2.d0*phi)

    case default
       write(6,*)mynum,'unknown number of simulations',isim
    end select

end function

!-------------------------------------------------------------------------------
subroutine build_kdtree(this)
    class(netcdf_type)    :: this

    if (.not.this%meshes_read) then
        print *, 'ERROR in build_kdtree(): Meshes have not been read yet'
        print *, 'Call read_meshes() before build_kdtree!'
        stop
    end if

    ! KDtree in forward field
    this%fwdtree => kdtree2_create(reshape([this%fwdmesh%s, this%fwdmesh%z], &
                                           [2, this%fwdmesh%npoints]    ),&
                                   dim = 2, &
                                   sort = .false.,&
                                   rearrange = .true.)

    ! KDtree in backward field
    this%bwdtree => kdtree2_create(reshape([this%bwdmesh%s, this%bwdmesh%z], &
                                           [2, this%bwdmesh%npoints]    ),&
                                   dim = 2, &
                                   sort = .false.,&
                                   rearrange = .true.)

end subroutine build_kdtree

!-------------------------------------------------------------------------------

subroutine read_meshes(this)
    use netcdf

    class(netcdf_type)             :: this

    integer                        :: ncvarid_mesh_s, ncvarid_mesh_z
   
    if (.not.this%files_open) then
        print *, 'ERROR in read_meshes(): Files have not been opened!'
        print *, 'Call open_files() before read_meshes()'
        stop
    end if

    ! Forward SEM mesh
    if (mynum.eq.0) write(6,*) 'Read SEM mesh from first forward simulation'
    
    call check( nf90_get_att(ncid = this%fwd(1)%ncid, &
                             name = 'npoints', &
                             varid = NF90_GLOBAL, &
                             values = this%fwdmesh%npoints) ) 

    allocate(this%fwdmesh%s(this%fwdmesh%npoints))
    allocate(this%fwdmesh%z(this%fwdmesh%npoints))
    
    call  getvarid( ncid  = this%fwd(1)%mesh, &
                    name  = "mesh_S", &
                    varid = ncvarid_mesh_s) 
    call  getvarid( ncid  = this%fwd(1)%mesh, &
                    name  = "mesh_Z", &
                    varid = ncvarid_mesh_z) 

    call check( nf90_get_var(ncid   = this%fwd(1)%mesh, &
                             varid  = ncvarid_mesh_s, &
                             values = this%fwdmesh%s )) 
    call check( nf90_get_var(ncid   = this%fwd(1)%mesh, &
                             varid  = ncvarid_mesh_z, &
                             values = this%fwdmesh%z )) 
   
    ! Backward SEM mesh                     
    if (mynum.eq.0) write(6,*) 'Read SEM mesh from first backward simulation'
    
    call check( nf90_get_att(ncid = this%bwd(1)%ncid, &
                             name = 'npoints', &
                             varid = NF90_GLOBAL, &
                             values = this%bwdmesh%npoints) ) 

    allocate(this%bwdmesh%s(this%bwdmesh%npoints))
    allocate(this%bwdmesh%z(this%bwdmesh%npoints))
    
    call  getvarid( ncid  = this%bwd(1)%mesh, &
                    name  = "mesh_S", &
                    varid = ncvarid_mesh_s) 
    call  getvarid( ncid  = this%bwd(1)%mesh, &
                    name  = "mesh_Z", &
                    varid = ncvarid_mesh_z) 

    call check( nf90_get_var(ncid   = this%bwd(1)%mesh, &
                             varid  = ncvarid_mesh_s, &
                             values = this%bwdmesh%s )) 
    call check( nf90_get_var(ncid   = this%bwd(1)%mesh, &
                             varid  = ncvarid_mesh_z, &
                             values = this%bwdmesh%z )) 

    this%meshes_read = .true.

end subroutine
 

!-------------------------------------------------------------------------------
subroutine check(status)
! Translates netcdf error codes into error messages

  implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop 2
  end if
end subroutine check  
!-------------------------------------------------------------------------------

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



!-------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------


end module
