module readfields
    use global_parameters, only            : sp, dp, pi
    use type_parameter, only               : src_param_type, rec_param_type, parameter_type
    use buffer, only                       : buffer_type
    use netcdf
    use kdtree2_module                     

    implicit none

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
        type(buffer_type)                 :: buffer
    end type

    type netcdf_type
        private

        integer                           :: nsim_fwd, nsim_bwd
        type(ncparamtype), allocatable    :: fwd(:)
        type(ncparamtype), allocatable    :: bwd(:)

        type(kdtree2), pointer            :: fwdtree, bwdtree
        !type(kdtree2_result), allocatable :: nextpoint(:)

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
            procedure, pass               :: close_files

    end type

    integer, parameter                    :: mynum = 0

contains

subroutine set_params(this, parameters)
    class(netcdf_type)            :: this
    class(parameter_type)         :: parameters
    integer                       :: isim


    this%nsim_fwd = parameters%nsim_fwd
    this%nsim_bwd =  1 !parameters%nsim_bwd
    
    allocate( this%fwd(this%nsim_fwd) )
    allocate( this%bwd(this%nsim_bwd) )

    select case(this%nsim_fwd)
    case(1)
        this%fwd(1)%meshdir = parameters%dir_fwdmesh

    case(2) 
        this%fwd(1)%meshdir = trim(parameters%dir_fwdmesh)//'/PZ/'
        this%fwd(2)%meshdir = trim(parameters%dir_fwdmesh)//'/PX/'

    case(4)
        this%fwd(1)%meshdir = trim(parameters%dir_fwdmesh)//'/MZZ/'
        this%fwd(2)%meshdir = trim(parameters%dir_fwdmesh)//'/MXX_P_MYY/'
        this%fwd(3)%meshdir = trim(parameters%dir_fwdmesh)//'/MXZ_MYZ/'
        this%fwd(4)%meshdir = trim(parameters%dir_fwdmesh)//'/MXY_MXX_M_MYY/'
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
        
        write(6,format20) trim(filename), mynum
        status = nf90_open(       path     = filename,              &
                                  mode     = nf90_nowrite,          &
                                  ncid     = this%fwd(isim)%ncid)

        call check(nf90_inq_ncid( ncid     = this%fwd(isim)%ncid,   &
                                  name     = "Snapshots",           &
                                  grp_ncid = this%fwd(isim)%snap))

        
        call check(nf90_inq_varid(ncid     = this%fwd(isim)%snap,   &
                                  name     = "straintrace",         &
                                  varid    = this%fwd(isim)%straintrace) )
        write(6,format21) this%fwd(isim)%ncid, this%fwd(isim)%snap 
        
        if (this%fwd(isim)%ordered_output) then
            filename = trim(this%fwd(isim)%meshdir)//'Data/axisem_output.nc4'
       
            write(6,format20) trim(filename), mynum
            status = nf90_open(   path     = filename,              &
                                  mode     = nf90_nowrite,          &
                                  ncid     = this%fwd(isim)%ncid)

        end if
        
        call check(nf90_inq_ncid(ncid      = this%fwd(isim)%ncid,   &
                                 name      = "Surface",             &
                                 grp_ncid  = this%fwd(isim)%surf))

        call check(nf90_inq_ncid(ncid      = this%fwd(isim)%ncid,   &
                                 name      = "Mesh",                &
                                 grp_ncid  = this%fwd(isim)%mesh))

        call nc_read_att_int(    this%fwd(isim)%ndumps,             &
                                 'number of strain dumps',          &
                                 this%fwd(isim))
    end do
        

    do isim = 1, this%nsim_bwd
        ! Backward wavefield
        filename=trim(this%bwd(isim)%meshdir)//'ordered_output.nc4'
        
        inquire(file=filename, exist=this%bwd(isim)%ordered_output)
        if (.not.this%bwd(isim)%ordered_output) then
            filename = trim(this%bwd(isim)%meshdir)//'Data/axisem_output.nc4'
        end if

        write(6,format20) trim(filename), mynum
        call check( nf90_open(    path     = filename,              &
                                  mode     = nf90_nowrite,          &
                                  ncid     = this%bwd(isim)%ncid) )

        call check(nf90_inq_ncid( ncid     = this%bwd(isim)%ncid,   &
                                  name     = "Snapshots",           &
                                  grp_ncid = this%bwd(isim)%snap))

        call check(nf90_inq_varid(ncid     = this%bwd(isim)%snap,   &
                                  name     = "straintrace",         &
                                  varid    = this%bwd(isim)%straintrace) )
        write(6,format21) this%bwd(isim)%ncid, this%bwd(isim)%snap 
        
        if (this%bwd(isim)%ordered_output) then
            filename = trim(this%bwd(isim)%meshdir)//'Data/axisem_output.nc4'
       
            write(6,format20) trim(filename), mynum
            status = nf90_open(   path     = filename,              &
                                  mode     = nf90_nowrite,          &
                                  ncid     = this%bwd(isim)%ncid)

        end if

        call check(nf90_inq_ncid( ncid     = this%bwd(isim)%ncid,   &
                                  name     = "Surface",             &
                                  grp_ncid = this%bwd(isim)%surf))
        call check(nf90_inq_ncid( ncid     = this%bwd(isim)%ncid,   &
                                  name     = "Mesh",                &
                                  grp_ncid = this%bwd(isim)%mesh))
        call nc_read_att_int(     this%bwd(isim)%ndumps,            &
                                  'number of strain dumps',         &
                                  this%bwd(isim))
    end do

    call flush(6) 
    this%files_open = .true.

    do isim = 1, this%nsim_fwd
       status = this%fwd(isim)%buffer%init(100, this%fwd(isim)%ndumps)
    end do
    do isim = 1, this%nsim_bwd
       status = this%bwd(isim)%buffer%init(100, this%bwd(isim)%ndumps)
    end do

end subroutine open_files

!-------------------------------------------------------------------------------
subroutine close_files(this)
    class(netcdf_type)   :: this
    integer              :: status, isim

    do isim = 1, this%nsim_fwd
       status = nf90_close(this%fwd(isim)%ncid)
       status = this%fwd(isim)%buffer%freeme()
    end do
    deallocate(this%fwd)
    do isim = 1, this%nsim_bwd
       status = nf90_close(this%bwd(isim)%ncid)
       status = this%bwd(isim)%buffer%freeme()
    end do
    deallocate(this%bwd)

end subroutine close_files

!-------------------------------------------------------------------------------
function load_fw_points(this, coordinates, source_params)
    use sorting, only                               : mergesort_3
    class(netcdf_type)                             :: this
    real(kind=dp), intent(in)                      :: coordinates(:,:)
    type(src_param_type)                           :: source_params
    real(kind=sp)                                  :: load_fw_points(this%fwd(1)%ndumps, size(coordinates,2))
    real(kind=sp)                                  :: utemp(this%fwd(1)%ndumps)
    type(kdtree2_result), allocatable              :: nextpoint(:)

    integer                                        :: npoints
    integer, dimension(size(coordinates,2))        :: pointid, idx
    integer                                        :: ipoint, isim, status
    real(kind=sp), dimension(size(coordinates,2))  :: rotmesh_s, rotmesh_phi, rotmesh_z
    
    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_fw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       stop 2
    end if
    npoints = size(coordinates,2)

    load_fw_points(:,:) = 0.0
    
    ! Rotate points to FWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,           &
                          coordinates(1,:)*1d3, coordinates(2,:)*1d3, coordinates(3,:)*1d3, &
                          source_params%lon, source_params%colat)

    allocate(nextpoint(1))
    do ipoint = 1, npoints
        call kdtree2_n_nearest( this%fwdtree,                           &
                                [rotmesh_s(ipoint), rotmesh_z(ipoint)], &
                                nn = 1,                                 &
                                results = nextpoint )
        
        pointid(ipoint) = nextpoint(1)%idx
        !print *, 'Original coordinates: ', coordinates(:,ipoint)
        !print *, 'Coordinates:    ', rotmesh_s(ipoint), rotmesh_z(ipoint), ', next pointid: ', pointid(ipoint)
        !print *, 'CO of SEM point:', this%fwdmesh%s(pointid(ipoint)), this%fwdmesh%z(pointid(ipoint))

        idx(ipoint) = ipoint

    end do

    !call mergesort_3(pointid, idx)

    do ipoint = 1, npoints
        

        isim = 1
        !print *, 'Azim factors: '
        do isim = 1, this%nsim_fwd
           !write(*,*) 'Reading point ', pointid(ipoint), ' (', ipoint, ') with ', &
           !           this%fwd(isim)%ndumps, ' values'
            status = this%fwd(isim)%buffer%get(pointid(ipoint), utemp)
            if (status.ne.0) then
               write(*,*) 'Did not find point', ipoint, ' in buffer, rereading'
               if (this%fwd(isim)%ordered_output) then
                  call check( nf90_get_var( ncid   = this%fwd(isim)%snap,        & 
                                            varid  = this%fwd(isim)%straintrace, &
                                            start  = [1, pointid(ipoint)],       &
                                            count  = [this%fwd(isim)%ndumps, 1], &
                                            values = utemp) )
               else
                  call check( nf90_get_var( ncid   = this%fwd(isim)%snap,        & 
                                            varid  = this%fwd(isim)%straintrace, &
                                            start  = [pointid(ipoint), 1],       &
                                            count  = [1, this%fwd(isim)%ndumps], &
                                            values = utemp) )
               end if
               status = this%fwd(isim)%buffer%put(pointid(ipoint), utemp)
            else
               write(*,*) 'Found point', ipoint, ' (',pointid(ipoint),') in buffer!'
            end if
            load_fw_points(:, (ipoint)) = load_fw_points(:,ipoint) &
                                        + azim_factor(rotmesh_phi(ipoint), &
                                                      source_params%mij, isim) &
                                        * utemp
            !load_fw_points(:, ipoint) = utemp
            !print *, 'isim', isim, '; azim. factor:', azim_factor(rotmesh_phi(ipoint),&
            !                                                      source_params%mij, isim)
        end do !isim
        
    end do !ipoint

end function load_fw_points

!-------------------------------------------------------------------------------
function load_bw_points(this, coordinates, receiver)
    use sorting, only                               : mergesort_3
    class(netcdf_type)                             :: this
    real(kind=dp), intent(in)                      :: coordinates(:,:)
    type(rec_param_type)                           :: receiver
    real(kind=sp)                                  :: load_bw_points(this%bwd(1)%ndumps, size(coordinates,2))
    real(kind=sp)                                  :: utemp(this%bwd(1)%ndumps)
    type(kdtree2_result), allocatable              :: nextpoint(:)

    integer, dimension(size(coordinates,2))        :: pointid, idx
    integer                                        :: npoints
    integer                                        :: ipoint, status
    real(kind=sp), dimension(size(coordinates,2))  :: rotmesh_s, rotmesh_phi, rotmesh_z

    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_bw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       stop 2
    end if
    npoints = size(coordinates,2)

    !allocate(load_bw_points(this%bwd(1)%ndumps,npoints))
    load_bw_points(:,:) = 0.0
    
    ! Rotate points to BWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,           &
                          coordinates(1,:), coordinates(2,:), coordinates(3,:), &
                          receiver%lon, receiver%colat)

    allocate(nextpoint(1))
    do ipoint = 1, npoints
        call kdtree2_n_nearest( this%bwdtree,                           &
                                [rotmesh_s(ipoint), rotmesh_z(ipoint)], &
                                nn = 1,                                 &
                                results = nextpoint )
        pointid(ipoint) = nextpoint(1)%idx
    end do
    
    !call mergesort_3(pointid, idx)


    do ipoint = 1, npoints
        
        status = this%bwd(1)%buffer%get(pointid(ipoint), utemp)
        if (status.ne.0) then
           write(*,*) 'Did not find point', ipoint, ' in buffer, rereading'
           if (this%bwd(1)%ordered_output) then
              call check( nf90_get_var( ncid   = this%bwd(1)%snap,        & 
                                        varid  = this%bwd(1)%straintrace, &
                                        start  = [1, pointid(ipoint)],       &
                                        count  = [this%bwd(1)%ndumps, 1], &
                                        values = utemp) )
           else
              call check( nf90_get_var( ncid   = this%bwd(1)%snap,        & 
                                        varid  = this%bwd(1)%straintrace, &
                                        start  = [pointid(ipoint), 1],       &
                                        count  = [1, this%bwd(1)%ndumps], &
                                        values = utemp) )
           end if
           status = this%bwd(1)%buffer%put(pointid(ipoint), utemp)
        else
           write(*,*) 'Found point', ipoint, ' (',pointid(ipoint),') in buffer!'
        end if


        select case(receiver%component)
        case('Z')
            load_bw_points(:,(ipoint)) =                      utemp
        case('R')
            load_bw_points(:,(ipoint)) =   cos(rotmesh_phi) * utemp
        case('T')
            load_bw_points(:,(ipoint)) = - sin(rotmesh_phi) * utemp 
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
    class(netcdf_type)         :: this
    real(kind=sp), allocatable :: mesh(:,:)
    integer                    :: ipoint

    if (.not.this%meshes_read) then
        print *, 'ERROR in build_kdtree(): Meshes have not been read yet'
        print *, 'Call read_meshes() before build_kdtree!'
        stop
    end if

    allocate(mesh(2, this%fwdmesh%npoints))
    mesh = transpose(reshape([this%fwdmesh%s, this%fwdmesh%z], [this%fwdmesh%npoints, 2]))

    !do ipoint = 1, this%fwdmesh%npoints
    !   write(41,*) mesh(:,ipoint)
    !   write(42,*) this%fwdmesh%s(ipoint), this%fwdmesh%z(ipoint)
    !end do
    write(*,*) 'Building forward KD-Tree'
    ! KDtree in forward field
    !this%fwdtree => kdtree2_create(reshape([this%fwdmesh%s, this%fwdmesh%z], &
    !                                       [2, this%fwdmesh%npoints]    ),&
    this%fwdtree => kdtree2_create(mesh,              &
                                   dim = 2,           &
                                   sort = .true.,     &
                                   rearrange = .true.)

    ! KDtree in backward field
    write(*,*) 'Building backward KD-Tree'
    mesh = transpose(reshape([this%bwdmesh%s, this%bwdmesh%z], [this%bwdmesh%npoints, 2]))
    this%bwdtree => kdtree2_create(mesh,              &
                                   dim = 2,           &
                                   sort = .true.,     &
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
    
    call nc_read_att_int(this%fwdmesh%npoints, &
                         'npoints', &
                         this%fwd(1))

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
    
    call nc_read_att_int(this%bwdmesh%npoints, &
                         'npoints', &
                         this%bwd(1))

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
      write(6,*) 'Could not find attribute ', trim(attribute_name)
      write(6,*) ' in NetCDF file ', trim(nc%meshdir), '/Data/axisem_output.nc4'
      write(6,*) ' with NCID: ', nc%ncid
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
      write(6,*) 'Could not find attribute ', trim(attribute_name)
      write(6,*) ' in NetCDF file ', trim(nc%meshdir), '/Data/axisem_output.nc4'
      write(6,*) ' with NCID: ', nc%ncid
      stop 2
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
      write(6,*) 'Could not find attribute ', trim(attribute_name)
      write(6,*) ' in NetCDF file ', trim(nc%meshdir), '/Data/axisem_output.nc4'
      write(6,*) ' with NCID: ', nc%ncid
      stop 2
  end if
end subroutine nc_read_att_real
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine rotate_frame_rd(npts, srd,phird,zrd, xgd,ygd,zgd, phigr,thetagr)

    implicit none
    integer, intent(in)                         :: npts
    !< Number of points to rotate
    
    real(kind=dp), dimension(npts), intent(in)  :: xgd, ygd, zgd 
    !< Coordinates to rotate (in x, y, z)
    
    real(kind=dp), intent(in)                   :: phigr, thetagr
    !< Rotation angles phi and theta
    
    real(kind=sp), dimension(npts), intent(out) :: srd, zrd, phird
    !< Rotated coordinates (in s, z, phi)

    real(kind=dp), dimension(npts)              :: xp, yp, zp, xp_cp, yp_cp, zp_cp
    real(kind=dp)                               :: phi_cp, rgd, thetagd
    integer                                     :: ipt


    !!first rotation (longitude)
    xp_cp =  xgd * dcos(phigr) + ygd * dsin(phigr)
    yp_cp = -xgd * dsin(phigr) + ygd * dcos(phigr)
    zp_cp =  zgd

    !second rotation (colat)
    xp = xp_cp * dcos(thetagr) - zp_cp * dsin(thetagr)
    yp = yp_cp 
    zp = xp_cp * dsin(thetagr) + zp_cp * dcos(thetagr)

    srd = dsqrt(xp*xp + yp*yp)
    zrd = zp
    do ipt = 1, npts
       phi_cp = datan2(yp(ipt), xp(ipt))
       if (phi_cp<0.d0) then
          phird(ipt) = 2.d0 * pi + phi_cp
       else
          phird(ipt) = phi_cp
       endif
       if (phigr==0.0 .and. ygd(ipt)==0.0)  phird(ipt)=0.
    enddo

    write(6,*)'Done with rotating frame rd.'

end subroutine rotate_frame_rd
!-------------------------------------------------------------------------------



end module
