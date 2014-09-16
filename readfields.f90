!=========================================================================================
module readfields
    use global_parameters, only            : sp, dp, pi, deg2rad, rad2deg, verbose, lu_out, &
                                             myrank, id_buffer, id_netcdf, id_rotate,       &
                                             id_load_strain, id_kdtree, id_calc_strain,     &
                                             id_find_point_fwd, id_find_point_bwd, id_lagrange

    use source_class,      only            : src_param_type
    use receiver_class,    only            : rec_param_type
    use buffer,            only            : buffer_type
    use clocks_mod,        only            : tick
    use commpi,            only            : pabort
    use nc_routines,       only            : getgrpid, getvarid, nc_open_for_read, nc_getvar, &
                                             check

    use rotations
    use netcdf
    use kdtree2_module                     

    implicit none
    private
    public                                 :: semdata_type

    type meshtype
        real(kind=sp), allocatable         :: s(:), z(:)            ! Coordinates of all GLL points
        real(kind=sp), allocatable         :: s_mp(:), z_mp(:)      ! Coordinates of element midpoints
        integer, allocatable               :: corner_point_ids(:,:) ! (4,nelem)
        integer, allocatable               :: eltype(:)             ! (nelem)
        integer                            :: npoints, nelem
        real(kind=dp), allocatable         :: theta(:)
        integer                            :: nsurfelem
    end type

    type ncparamtype
        integer                            :: ncid
        integer                            :: file_version
        real(kind=dp)                      :: planet_radius
        real(kind=dp)                      :: rmin, rmax
        real(kind=dp)                      :: colatmin, colatmax
        integer                            :: snap, surf, mesh, seis  ! Group IDs
        integer                            :: strainvarid(6)          ! Variable IDs
        integer                            :: displvarid(3)           ! Variable IDs
        integer                            :: seis_disp, seis_velo    ! Variable IDs
        integer                            :: stf_varid               ! Variable IDs
        integer                            :: stf_d_varid             ! Variable IDs
        integer                            :: fem_mesh_varid          ! Variable IDs
        integer                            :: sem_mesh_varid          ! Variable IDs
        integer                            :: eltype_varid            ! Variable IDs
        integer                            :: axis_varid              ! Variable IDs
        integer                            :: mesh_s_varid            ! Variable IDs
        integer                            :: mesh_z_varid            ! Variable IDs
        integer                            :: mesh_rho_varid          ! Variable IDs
        integer                            :: mesh_lambda_varid       ! Variable IDs
        integer                            :: mesh_mu_varid           ! Variable IDs
        integer                            :: mesh_vs_varid           ! Variable IDs
        integer                            :: mesh_vp_varid           ! Variable IDs
        integer                            :: chunk_gll
        integer                            :: count_error_pointoutside
        character(len=200)                 :: meshdir
        character(len=12)                  :: dump_type
        integer                            :: ndumps, nseis, ngll, npol
        integer                            :: source_shift_samples    
        real(kind=dp)                      :: source_shift_t
        character(len=10)                  :: source_type
        character(len=10)                  :: excitation_type
        real(kind=sp), allocatable         :: stf(:), stf_d(:)
        type(buffer_type)                  :: buffer_strain
        type(buffer_type)                  :: buffer_disp
        type(buffer_type)                  :: buffer
        real(kind=dp)                      :: dt
        real(kind=dp)                      :: amplitude
        real(kind=dp), public, allocatable :: G1(:,:), G1T(:,:)
        real(kind=dp), public, allocatable :: G2(:,:), G2T(:,:)
        real(kind=dp), public, allocatable :: G0(:)
        real(kind=dp), public, allocatable :: gll_points(:), glj_points(:)
    end type

    type semdata_type
        private

        integer, public                    :: nsim_fwd, nsim_bwd
        type(ncparamtype), allocatable     :: fwd(:)
        type(ncparamtype), allocatable     :: bwd(:)

        type(kdtree2), pointer             :: fwdtree, bwdtree
        type(meshtype)                     :: fwdmesh, bwdmesh

        logical                            :: params_set
        logical                            :: files_open
        logical                            :: meshes_read
        logical                            :: kdtree_built
        
        character(len=32)                  :: strain_type  !< full tensor or straintrace
        integer                            :: ndim          !< Number of dimensions which has to be read to calculate 
                                                            !! Kernel on parameter model_param

        real(kind=dp), public              :: dt
        integer,       public              :: ndumps, decimate_factor
        integer,       public              :: nseis 
        integer,       public              :: npol
        real(kind=dp), public, allocatable :: G1(:,:), G1T(:,:)
        real(kind=dp), public, allocatable :: G2(:,:), G2T(:,:)
        real(kind=dp), public, allocatable :: G0(:)
        real(kind=dp), public, allocatable :: gll_points(:), glj_points(:)
        real(kind=dp), public              :: windowlength
        real(kind=dp), public              :: timeshift_fwd, timeshift_bwd
        real(kind=dp), public              :: amplitude_fwd, amplitude_bwd
        real(kind=dp), public, allocatable :: veloseis(:,:), dispseis(:,:)
        real(kind=dp), public, allocatable :: stf_fwd(:), stf_bwd(:)
        real(kind=dp), public, allocatable :: stf_d_fwd(:), stf_d_bwd(:)
        integer                            :: strain_buffer_size
        integer                            :: displ_buffer_size
        character(len=12)                  :: dump_type
         
        real(kind=dp), dimension(3,3)      :: rot_mat, trans_rot_mat

        contains 
            procedure, pass                :: get_ndim 
            procedure, pass                :: set_params
            procedure, pass                :: open_files
            procedure, pass                :: check_consistency
            procedure, pass                :: read_meshes
            procedure, pass                :: build_kdtree
            procedure, pass                :: load_fw_points
            procedure, pass                :: load_bw_points
            procedure, pass                :: close_files
            procedure, pass                :: load_seismogram

    end type
contains

!-----------------------------------------------------------------------------------------
function get_ndim(this)
    class(semdata_type)            :: this
    integer                        :: get_ndim
    if (.not.this%params_set) then
        print *, 'ERROR in get_ndim(): Parameters have to be set first'
        print *, 'Call set_params before get_ndim()'
        call pabort
    end if
    get_ndim = this%ndim
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine set_params(this, fwd_dir, bwd_dir, strain_buffer_size, displ_buffer_size, strain_type)
    class(semdata_type)            :: this
    character(len=512), intent(in) :: fwd_dir, bwd_dir
    integer,            intent(in) :: strain_buffer_size
    integer,            intent(in) :: displ_buffer_size
    character(len=32),   intent(in), optional :: strain_type
    character(len=512)             :: dirnam
    logical                        :: moment=.false., force=.false., single=.false.

    this%strain_buffer_size = strain_buffer_size
    this%displ_buffer_size = displ_buffer_size

    dirnam = trim(fwd_dir)//'/MZZ/simulation.info'
    write(lu_out,*) 'Inquiring: ', trim(dirnam)
    inquire( file = trim(dirnam), exist = moment)

    dirnam = trim(fwd_dir)//'/PZ/simulation.info'
    write(lu_out,*) 'Inquiring: ', trim(dirnam)
    inquire( file = trim(dirnam), exist = force)

    dirnam = trim(fwd_dir)//'/simulation.info'
    write(lu_out,*) 'Inquiring: ', trim(dirnam)
    inquire( file = trim(dirnam), exist = single)

    if (moment) then
       this%nsim_fwd = 4
       write(lu_out,*) 'Forward simulation was ''moment'' source'
    elseif (force) then
       this%nsim_fwd = 2
       write(lu_out,*) 'Forward simulation was ''forces'' source'
    elseif (single) then
       this%nsim_fwd = 1
       write(lu_out,*) 'Forward simulation was ''single'' source'
    else 
       this%nsim_fwd = 0
       write(*,*) 'ERROR: Forward rundir does not seem to be an axisem rundirectory'
       call pabort
    end if

    moment = .false.
    force  = .false.
    single = .false.

    dirnam = trim(bwd_dir)//'/MZZ/simulation.info'
    write(lu_out,*) 'Inquiring: ', trim(dirnam)
    inquire( file = trim(dirnam), exist = moment)

    dirnam = trim(bwd_dir)//'/PZ/simulation.info'
    write(lu_out,*) 'Inquiring: ', trim(dirnam)
    inquire( file = trim(dirnam), exist = force)

    dirnam = trim(bwd_dir)//'/simulation.info'
    write(lu_out,*) 'Inquiring: ', trim(dirnam)
    inquire( file = trim(dirnam), exist = single)

    if (moment) then
       this%nsim_bwd = 4
       write(lu_out,*) 'Backward simulation was ''moment'' source'
       write(lu_out,*) 'This is not implemented yet!'
       call pabort()
    elseif (force) then
       this%nsim_bwd = 2
       write(lu_out,*) 'Backward simulation was ''forces'' source'
    elseif (single) then
       this%nsim_bwd = 1
       write(lu_out,*) 'Backward simulation was ''single'' source'
    else 
       this%nsim_bwd = 0
       write(lu_out,*) 'WARNING: Backward rundir does not seem to be an axisem rundirectory'
       write(lu_out,*) 'continuing anyway, as this is default in db mode'
    end if

    allocate( this%fwd(this%nsim_fwd) )
    allocate( this%bwd(this%nsim_bwd) )

    select case(this%nsim_fwd)
    case(1)    ! Single
        this%fwd(1)%meshdir = fwd_dir//'/'

    case(2)    ! Forces
        this%fwd(1)%meshdir = trim(fwd_dir)//'/PZ/'
        this%fwd(2)%meshdir = trim(fwd_dir)//'/PX/'

    case(4)    ! Moment
        this%fwd(1)%meshdir = trim(fwd_dir)//'/MZZ/'
        this%fwd(2)%meshdir = trim(fwd_dir)//'/MXX_P_MYY/'
        this%fwd(3)%meshdir = trim(fwd_dir)//'/MXZ_MYZ/'
        this%fwd(4)%meshdir = trim(fwd_dir)//'/MXY_MXX_M_MYY/'
    end select
    
    select case(this%nsim_bwd)
    case(1)    ! Single
        this%bwd(1)%meshdir = bwd_dir//'/'

    case(2)    ! Forces
        this%bwd(1)%meshdir = trim(bwd_dir)//'/PZ/'
        this%bwd(2)%meshdir = trim(bwd_dir)//'/PX/'

    case(4)    ! Moment
        this%bwd(1)%meshdir = trim(bwd_dir)//'/MZZ/'
        this%bwd(2)%meshdir = trim(bwd_dir)//'/MXX_P_MYY/'
        this%bwd(3)%meshdir = trim(bwd_dir)//'/MXZ_MYZ/'
        this%bwd(4)%meshdir = trim(bwd_dir)//'/MXY_MXX_M_MYY/'

    end select
    
    if (present(strain_type)) then
       this%strain_type = strain_type
    else
       this%strain_type = 'straintensor_full'
    end if

    select case(trim(this%strain_type))
    case('straintensor_trace')
       this%ndim = 1
    case('straintensor_full')
       this%ndim = 6
    case default
        print *, 'ERROR in set_params(): unknown straintensor output format '//this%strain_type
        call pabort
    end select
    write(lu_out, *) 'Straintensor output variant: ', trim(this%strain_type), &
                     ', Dimension of wavefields: ', this%ndim

    call flush(lu_out)
    this%params_set = .true.

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine open_files(this)

    class(semdata_type)              :: this
    integer                          :: status, isim, chunks(2), deflev
    character(len=200)               :: format20, format21, filename
    character(len=11)                :: nc_strain_varnamelist(6)
    character(len=11)                :: nc_displ_varnamelist(3)
    real(kind=sp)                    :: temp
    integer                          :: istrainvar, idisplvar

    if (.not.this%params_set) then
        print *, 'ERROR in open_files(): Parameters have to be set first'
        print *, 'Call set_params before open_files()'
        call pabort
    end if

    nc_strain_varnamelist = ['strain_dsus', 'strain_dsuz', 'strain_dpup', &
                             'strain_dsup', 'strain_dzup', 'straintrace']
           
    nc_displ_varnamelist  = ['disp_s     ', 'disp_p     ', 'disp_z     ']

    format20 = "('  Trying to open NetCDF file ', A, ' on CPU ', I5)"
    format21 = "('  Succeded,  has NCID ', I6, ', Snapshots group NCID: ', I6)"

    do isim = 1, this%nsim_fwd

        ! Forward wavefield
        filename=trim(this%fwd(isim)%meshdir)//'/Data/ordered_output.nc4'
        
        if (verbose>0) write(lu_out,format20) trim(filename), myrank
        call nc_open_for_read(    filename = filename,              &
                                  ncid     = this%fwd(isim)%ncid) 

        call nc_read_att_int(this%fwd(isim)%file_version, 'file version', this%fwd(isim))

        call nc_read_att_dble(this%fwd(isim)%planet_radius, 'planet radius', this%fwd(isim))
        this%fwd(isim)%planet_radius = this%fwd(isim)%planet_radius * 1d3

        call nc_read_att_dble(this%fwd(isim)%rmin, 'kernel wavefield rmin', this%fwd(isim))
        call nc_read_att_dble(this%fwd(isim)%rmax, 'kernel wavefield rmax', this%fwd(isim))
        this%fwd(isim)%rmin = this%fwd(isim)%rmin * 1d3
        this%fwd(isim)%rmax = this%fwd(isim)%rmax * 1d3

        call nc_read_att_dble(this%fwd(isim)%colatmin, 'kernel wavefield colatmin', this%fwd(isim))
        call nc_read_att_dble(this%fwd(isim)%colatmax, 'kernel wavefield colatmax', this%fwd(isim))

        call nc_read_att_char(this%fwd(isim)%dump_type, &
                              'dump type (displ_only, displ_velo, fullfields)', &
                               this%fwd(isim))

        call nc_read_att_char(this%fwd(isim)%source_type, &
                              'source type', &
                               this%fwd(isim))

        call nc_read_att_char(this%fwd(isim)%excitation_type, &
                              'excitation type', &
                               this%fwd(isim))

        call getgrpid(  ncid     = this%fwd(isim)%ncid,   &
                        name     = "Snapshots",           &
                        grp_ncid = this%fwd(isim)%snap)


        
        if (trim(this%fwd(isim)%dump_type) == 'displ_only') then
            do idisplvar = 1, 3
                status = nf90_inq_varid(ncid  = this%fwd(isim)%snap,                  &
                                        name  = nc_displ_varnamelist(idisplvar),      &
                                        varid = this%fwd(isim)%displvarid(idisplvar)) 
                
                if (status.ne.NF90_NOERR) then
                    this%fwd(isim)%displvarid(idisplvar) = -1
                    if (idisplvar == 1) then
                        print *, 'Did not find variable ''disp_s'' in NetCDF file'
                        call pabort
                    end if
                end if
            end do
            call check(nf90_inquire_variable(ncid       = this%fwd(isim)%snap,   &
                                             varid      = this%fwd(isim)%displvarid(1), &
                                             chunksizes = chunks, &
                                             deflate_level = deflev) )

            call nc_read_att_int(this%fwd(isim)%npol, 'npol', this%fwd(isim))



        elseif (trim(this%fwd(isim)%dump_type) == 'fullfields') then
            do istrainvar = 1, 6
                status = nf90_inq_varid(ncid  = this%fwd(isim)%snap,                  &
                                        name  = nc_strain_varnamelist(istrainvar),    &
                                        varid = this%fwd(isim)%strainvarid(istrainvar)) 
                
                if (status.ne.NF90_NOERR) then
                    this%fwd(isim)%strainvarid(istrainvar) = -1
                    if (istrainvar == 6) then
                        print *, 'Did not find variable ''straintrace'' in NetCDF file'
                        call pabort
                    end if
                end if
            end do
            call check(nf90_inquire_variable(ncid       = this%fwd(isim)%snap,   &
                                             varid      = this%fwd(isim)%strainvarid(6), &
                                             chunksizes = chunks, &
                                             deflate_level = deflev) )


        else
           print *, 'ERROR: dump_type ', this%fwd(isim)%dump_type, ' not implemented!'
           call pabort
        endif



        write(lu_out, "('  FWD SIM:', I2, ', Chunksizes:', 2(I7), ', deflate level: ', I2)") &
              isim, chunks, deflev

        this%fwd(isim)%chunk_gll = chunks(1)

        if (verbose>0) write(lu_out,format21) this%fwd(isim)%ncid, this%fwd(isim)%snap 
        
        call getgrpid(           ncid      = this%fwd(isim)%ncid,   &
                                 name      = "Surface",             &
                                 grp_ncid  = this%fwd(isim)%surf)


        call getvarid(           ncid      = this%fwd(isim)%surf,   &
                                 name      = "displacement",         &
                                 varid     = this%fwd(isim)%seis_disp) 
        
        call getvarid(           ncid      = this%fwd(isim)%surf,   &
                                 name      = "velocity",            &
                                 varid     = this%fwd(isim)%seis_velo) 


        call getgrpid(           ncid      = this%fwd(isim)%ncid,   &
                                 name      = "Mesh",                &
                                 grp_ncid  = this%fwd(isim)%mesh)

        call nc_read_att_int(    this%fwd(isim)%ndumps,             &
                                 'number of strain dumps',          &
                                 this%fwd(isim))

        call nc_read_att_dble(   this%fwd(isim)%dt,               &
                                 'strain dump sampling rate in sec', &
                                 this%fwd(isim))

        call nc_read_att_int(    this%fwd(isim)%nseis,             &
                                 'length of seismogram  in time samples', &
                                 this%fwd(isim))

        call nc_read_att_real(   temp, &
                                 'source shift factor in sec',     &
                                 this%fwd(isim))
        this%fwd(isim)%source_shift_t = real(temp, kind=dp)
        
        call nc_read_att_int(    this%fwd(isim)%source_shift_samples,    &
                                 'source shift factor for deltat_coarse',     &
                                 this%fwd(isim))
        
        call nc_read_att_real(   temp,      &
                                 'scalar source magnitude',     &
                                 this%fwd(isim))
        this%fwd(isim)%amplitude = real(temp, kind=dp)

        !call check(nf90_inq_ncid( ncid     = this%fwd(isim)%ncid,   &
        !                          name     = "Seismograms",         &
        !                          grp_ncid = this%fwd(isim)%seis))

        call getvarid(            ncid     = this%fwd(isim)%surf,   &
                                  name     = "stf_dump",            &
                                  varid    = this%fwd(isim)%stf_varid)

        call getvarid(            ncid     = this%fwd(isim)%surf,   &
                                  name     = "stf_d_dump",            &
                                  varid    = this%fwd(isim)%stf_d_varid)        

        allocate( this%fwd(isim)%stf( this%fwd(isim)%ndumps ) )
        call check(nf90_get_var(  ncid   = this%fwd(isim)%surf,   &
                                  varid  = this%fwd(isim)%stf_varid, &
                                  values = this%fwd(isim)%stf  ))

        allocate( this%fwd(isim)%stf_d( this%fwd(isim)%ndumps ) )
        call check(nf90_get_var(  ncid   = this%fwd(isim)%surf,   &
                                  varid  = this%fwd(isim)%stf_d_varid, &
                                  values = this%fwd(isim)%stf_d  ))
        
    end do
        
    call flush(lu_out)

    do isim = 1, this%nsim_bwd
        ! Backward wavefield
        filename=trim(this%bwd(isim)%meshdir)//'/Data/ordered_output.nc4'
        
        if (verbose>0) write(lu_out,format20) trim(filename), myrank
        call nc_open_for_read(filename = filename,              &
                              ncid     = this%bwd(isim)%ncid) 

        call nc_read_att_int(this%bwd(isim)%file_version, 'file version', this%bwd(isim))

        call nc_read_att_dble(this%bwd(isim)%planet_radius, 'planet radius', this%bwd(isim))
        this%bwd(isim)%planet_radius = this%bwd(isim)%planet_radius * 1d3

        call nc_read_att_dble(this%bwd(isim)%rmin, 'kernel wavefield rmin', this%bwd(isim))
        call nc_read_att_dble(this%bwd(isim)%rmax, 'kernel wavefield rmax', this%bwd(isim))
        this%fwd(isim)%rmin = this%fwd(isim)%rmin * 1d3
        this%fwd(isim)%rmax = this%fwd(isim)%rmax * 1d3

        call nc_read_att_dble(this%bwd(isim)%colatmin, 'kernel wavefield colatmin', this%bwd(isim))
        call nc_read_att_dble(this%bwd(isim)%colatmax, 'kernel wavefield colatmax', this%bwd(isim))

        call nc_read_att_char(this%bwd(isim)%dump_type, &
                              'dump type (displ_only, displ_velo, fullfields)', &
                               this%bwd(isim))

        call nc_read_att_char(this%bwd(isim)%source_type, &
                              'source type', &
                               this%bwd(isim))

        call nc_read_att_char(this%bwd(isim)%excitation_type, &
                              'excitation type', &
                               this%bwd(isim))

        call getgrpid( ncid     = this%bwd(isim)%ncid,   &
                       name     = "Snapshots",           &
                       grp_ncid = this%bwd(isim)%snap)

        if (trim(this%bwd(isim)%dump_type) == 'displ_only') then
            do idisplvar = 1, 3
                status = nf90_inq_varid(ncid  = this%bwd(isim)%snap,                  &
                                        name  = nc_displ_varnamelist(idisplvar),      &
                                        varid = this%bwd(isim)%displvarid(idisplvar)) 
                
                if (status.ne.NF90_NOERR) then
                    this%bwd(isim)%displvarid(idisplvar) = -1
                    if (idisplvar == 1) then
                        print *, 'Did not find variable ''disp_s'' in NetCDF file'
                        call pabort
                    end if
                end if
            end do
            call check(nf90_inquire_variable(ncid       = this%bwd(isim)%snap,   &
                                             varid      = this%bwd(isim)%displvarid(1), &
                                             chunksizes = chunks, &
                                             deflate_level = deflev) )

            call nc_read_att_int(this%bwd(isim)%npol, 'npol', this%bwd(isim))

        elseif (trim(this%bwd(isim)%dump_type) == 'fullfields') then
            do istrainvar = 1, 6            
                status = nf90_inq_varid(ncid  = this%bwd(isim)%snap,                  &
                                        name  = nc_strain_varnamelist(istrainvar),    &
                                        varid = this%bwd(isim)%strainvarid(istrainvar)) 
    
                if (status.ne.NF90_NOERR) then
                    this%bwd(isim)%strainvarid(istrainvar) = -1
                    if (istrainvar == 6.) then
                        print *, 'Did not find variable ''straintrace'' in NetCDF file'
                        call pabort
                    end if
                end if
            end do
    
            call check(nf90_inquire_variable(ncid       = this%bwd(isim)%snap,   &
                                             varid      = this%bwd(isim)%strainvarid(6), &
                                             chunksizes = chunks, &
                                             deflate_level = deflev) )
        else
           print *, 'ERROR: dump_type ', this%bwd(isim)%dump_type, ' not implemented!'
           call pabort
        endif

        write(lu_out, "('BWD SIM:', I2, ', Chunksizes:', 2(I7), ', deflate level: ', I2)") &
              isim, chunks, deflev

        this%bwd(isim)%chunk_gll = chunks(1)
        
        if (verbose>0) write(lu_out,format21) this%bwd(isim)%ncid, this%bwd(isim)%snap 
        
        call getgrpid( ncid     = this%bwd(isim)%ncid,   &
                       name     = "Surface",             &
                       grp_ncid = this%bwd(isim)%surf)
        call getgrpid( ncid     = this%bwd(isim)%ncid,   &
                       name     = "Mesh",                &
                       grp_ncid = this%bwd(isim)%mesh)
        !call getgrpid( ncid     = this%bwd(isim)%ncid,   &
        !               name     = "Seismograms",         &
        !               grp_ncid = this%bwd(isim)%seis)

        call nc_read_att_int(     this%bwd(isim)%ndumps,            &
                                  'number of strain dumps',         &
                                  this%bwd(isim))

        call nc_read_att_int(     this%bwd(isim)%nseis,             &
                                  'length of seismogram  in time samples', &
                                  this%bwd(isim))

        call nc_read_att_dble(    this%bwd(isim)%dt,                &
                                  'strain dump sampling rate in sec', &
                                  this%bwd(isim))

        call nc_read_att_real(   temp, &
                                 'source shift factor in sec',     &
                                 this%bwd(isim))
        this%bwd(isim)%source_shift_t = real(temp, kind=dp)

        call nc_read_att_int(    this%bwd(isim)%source_shift_samples,    &
                                 'source shift factor for deltat_coarse',     &
                                 this%bwd(isim))
        
        call nc_read_att_real(   temp,      &
                                 'scalar source magnitude',     &
                                 this%bwd(isim))
        this%bwd(isim)%amplitude = real(temp, kind=dp)
        
        ! call getvarid(            ncid     = this%bwd(isim)%surf,   &
        !                           name     = "stf_dump",            &
        !                           varid    = this%bwd(isim)%stf_varid)



        call getvarid(            ncid     = this%bwd(isim)%surf,   &
                                  name     = "stf_dump",            &
                                  varid    = this%bwd(isim)%stf_varid)

        call getvarid(            ncid     = this%bwd(isim)%surf,   &
                                  name     = "stf_d_dump",            &
                                  varid    = this%bwd(isim)%stf_d_varid)


        
        allocate( this%bwd(isim)%stf( this%bwd(isim)%ndumps ) )
        call check(nf90_get_var(  ncid   = this%bwd(isim)%surf,   &
                                  varid  = this%bwd(isim)%stf_varid, &
                                  values = this%bwd(isim)%stf  ))

        allocate( this%bwd(isim)%stf_d( this%bwd(isim)%ndumps ) )
        call check(nf90_get_var(  ncid   = this%bwd(isim)%surf,   &
                                  varid  = this%bwd(isim)%stf_d_varid, &
                                  values = this%bwd(isim)%stf_d  ))

        
    end do


    call flush(lu_out)
    call this%check_consistency()


    call flush(6) 
    this%files_open = .true.

    !@TODO memory could be used more efficient for monopole sources in the buffers
    !Initialize Buffers. 
    select case(trim(this%dump_type))
    case('displ_only')
       do isim = 1, this%nsim_fwd
          status = this%fwd(isim)%buffer_disp%init(this%displ_buffer_size, this%fwd(isim)%ndumps, 3)
          select case(this%strain_type)
          case('straintensor_trace')
            status = this%fwd(isim)%buffer_strain%init(this%strain_buffer_size,      &
                                                       this%fwd(isim)%ndumps, &
                                                       this%fwd(isim)%npol+1,   &
                                                       this%fwd(isim)%npol+1)
          case('straintensor_full')
            status = this%fwd(isim)%buffer_strain%init(this%strain_buffer_size,      &
                                                       this%fwd(isim)%ndumps, &
                                                       this%fwd(isim)%npol+1,   &
                                                       this%fwd(isim)%npol+1,   &
                                                       6)
          end select
          this%fwd(isim)%count_error_pointoutside = 0
       end do

       do isim = 1, this%nsim_bwd
          status = this%bwd(isim)%buffer_disp%init(this%displ_buffer_size, this%bwd(isim)%ndumps, 3)
          select case(this%strain_type)
          case('straintensor_trace')
            status = this%bwd(isim)%buffer_strain%init(this%strain_buffer_size,      &
                                                       this%bwd(isim)%ndumps, &
                                                       this%bwd(isim)%npol+1,   &
                                                       this%bwd(isim)%npol+1)
          case('straintensor_full')
            status = this%bwd(isim)%buffer_strain%init(this%strain_buffer_size,      &
                                                       this%bwd(isim)%ndumps, &
                                                       this%bwd(isim)%npol+1,   &
                                                       this%bwd(isim)%npol+1,   &
                                                       6)
          end select
          this%bwd(isim)%count_error_pointoutside = 0
       end do
    case('fullfields')
       do isim = 1, this%nsim_fwd
          status = this%fwd(isim)%buffer%init(this%strain_buffer_size, this%fwd(isim)%ndumps, this%ndim)
       end do

       do isim = 1, this%nsim_bwd
          status = this%bwd(isim)%buffer%init(this%strain_buffer_size, this%bwd(isim)%ndumps, this%ndim)
       end do
    case default
       print *, 'Unknown dump type in solver'
       call pabort()
    end select


    call flush(lu_out)

end subroutine open_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine close_files(this)
    class(semdata_type)   :: this
    integer              :: status, isim

    ! Destroy kdtree
    if (this%kdtree_built) then
      call kdtree2_destroy(this%fwdtree)
      call kdtree2_destroy(this%bwdtree)
    end if

    ! Free buffers
    select case(trim(this%dump_type))
    case('fullfields')
      do isim = 1, this%nsim_fwd
         status = nf90_close(this%fwd(isim)%ncid)
         if (verbose>0) then
            write(lu_out,'(A,I1,A,F9.6)') ' Buffer efficiency fwd(', isim, '): ',  &
                                     this%fwd(isim)%buffer%efficiency()
         end if
         status = this%fwd(isim)%buffer%freeme()
      end do

      do isim = 1, this%nsim_bwd
         status = nf90_close(this%bwd(isim)%ncid)
         if (verbose>0) then
            write(lu_out,'(A,F9.6)') ' Buffer efficiency bwd   : ', & 
                                this%bwd(isim)%buffer%efficiency()
         end if
         status = this%bwd(isim)%buffer%freeme()
      end do

    case('displ_only')
      do isim = 1, this%nsim_fwd
         status = nf90_close(this%fwd(isim)%ncid)
         if (verbose>0) then
            write(lu_out,'(A,I1,A,F9.6)') ' Strain buffer efficiency fwd(', isim, '): ',  &
                                     this%fwd(isim)%buffer_strain%efficiency()
            write(lu_out,'(A,I1,A,F9.6)') ' Displ. buffer efficiency fwd(', isim, '): ',  &
                                     this%fwd(isim)%buffer_disp%efficiency()
         end if
         status = this%fwd(isim)%buffer_strain%freeme()
         status = this%fwd(isim)%buffer_disp%freeme()
      end do
      if (this%nsim_fwd > 0) &
         write(lu_out,'(A,I8)') ' Points outside of element (fwd): ', &
                                this%fwd(1)%count_error_pointoutside

      do isim = 1, this%nsim_bwd
         status = nf90_close(this%bwd(isim)%ncid)
         if (verbose>0) then
            write(lu_out,'(A,I1,A,F9.6)') ' Strain buffer efficiency bwd(', isim, '): ',  &
                                     this%bwd(isim)%buffer_strain%efficiency()
            write(lu_out,'(A,I1,A,F9.6)') ' Displ. buffer efficiency bwd(', isim, '): ',  &
                                     this%bwd(isim)%buffer_disp%efficiency()
         end if
         status = this%bwd(isim)%buffer_strain%freeme()
         status = this%bwd(isim)%buffer_disp%freeme()
      end do
      if (this%nsim_bwd > 0) &
         write(lu_out,'(A,I8)') ' Points outside of element (bwd): ', &
                                this%bwd(1)%count_error_pointoutside
    end select

    deallocate(this%fwd)
    deallocate(this%bwd)

    call flush(lu_out)

end subroutine close_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine check_consistency(this)
    !< Checks consistency of the wavefield dumps
    class(semdata_type)    :: this
    integer                :: isim
    real(kind=dp)          :: dt_agreed
    character(len=512)     :: fmtstring, fmtstring_stf
    character(len=12)      :: dump_type_agreed
    integer                :: ndumps_agreed, nseis_agreed, npol_agreed
    real(kind=dp)          :: source_shift_agreed_fwd, source_shift_agreed_bwd
    real(kind=dp)          :: amplitude_agreed_fwd, amplitude_agreed_bwd
    real(kind=dp), allocatable  :: stf_agreed_fwd(:), stf_d_agreed_fwd(:)
    real(kind=dp), allocatable  :: stf_agreed_bwd(:), stf_d_agreed_bwd(:)

    if (this%nsim_fwd > 0) then
       allocate(stf_agreed_fwd(this%fwd(1)%ndumps))
       allocate(stf_d_agreed_fwd(this%fwd(1)%ndumps))
    endif
    if (this%nsim_bwd > 0) then
       allocate(stf_agreed_bwd(this%fwd(1)%ndumps))
       allocate(stf_d_agreed_bwd(this%fwd(1)%ndumps))
    endif

    ! Check whether the dump_type is the same in all files
    dump_type_agreed = this%fwd(1)%dump_type
    

    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1)' 
    do isim = 1, this%nsim_fwd
       if (dump_type_agreed /= this%fwd(isim)%dump_type) then
          write(*,fmtstring) 'dump_type', isim
          call pabort
       end if
    end do   

    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1)' 

    do isim = 1, this%nsim_bwd
       if (dump_type_agreed /= this%bwd(isim)%dump_type) then
          write(*,fmtstring) 'dump_type', isim
          call pabort
       end if
    end do

    this%dump_type = dump_type_agreed

    ! Check whether the sampling period is the same in all files
    dt_agreed = this%fwd(1)%dt
    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the others")' 
    do isim = 1, this%nsim_fwd
       if (dt_agreed.ne.this%fwd(isim)%dt) then
          write(*,fmtstring) 'dt', isim, dt_agreed, this%fwd(isim)%dt
          call pabort
       end if
    end do

    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the forward case")' 


    do isim = 1, this%nsim_bwd
       if (dt_agreed.ne.this%bwd(isim)%dt) then
          write(*,fmtstring) 'dt', isim, dt_agreed, this%bwd(isim)%dt
          call pabort
       end if
    end do

    this%dt = dt_agreed
 
    ! Check whether npol is the same in all files
    if (trim(this%dump_type) == 'displ_only') then
        npol_agreed = this%fwd(1)%npol

        fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                    '  in simulation ", I1, " (", I7,") vs ", I7, " in the others")' 
        do isim = 1, this%nsim_fwd
           if (npol_agreed.ne.this%fwd(isim)%npol) then
              write(*,fmtstring) 'npol', isim, npol_agreed, this%fwd(isim)%npol
              call pabort
           end if
        end do

        fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                    '  in simulation ", I1, "(",I7,"s) vs ", I7, " in the forward case")' 

        do isim = 1, this%nsim_bwd
           if (npol_agreed.ne.this%bwd(isim)%npol) then
              write(*,fmtstring) 'npol', isim, npol_agreed, this%bwd(isim)%npol
              call pabort
           end if
        end do

        this%npol = npol_agreed
    endif


    ! Check whether the number of dumps (time samples) is the same in all files
    ndumps_agreed = this%fwd(1)%ndumps
    nseis_agreed  = this%fwd(1)%nseis

    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",I7,") vs ", I7, " in the others")' 
    do isim = 1, this%nsim_fwd
       if (ndumps_agreed.ne.this%fwd(isim)%ndumps) then
          write(*,fmtstring) 'ndumps', isim, ndumps_agreed, this%fwd(isim)%ndumps
          call pabort
       end if
       if (nseis_agreed.ne.this%fwd(isim)%nseis) then
          write(*,fmtstring) 'nseis', isim, nseis_agreed, this%fwd(isim)%nseis
          call pabort
       end if
    end do

    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",I7,"s) vs ", I7, " in the forward case")' 

    do isim = 1, this%nsim_bwd
       if (ndumps_agreed.ne.this%bwd(isim)%ndumps) then
          write(*,fmtstring) 'ndumps', isim, ndumps_agreed, this%bwd(isim)%ndumps
          call pabort
       end if
       if (nseis_agreed.ne.this%bwd(isim)%nseis) then
          write(*,fmtstring) 'nseis', isim, nseis_agreed, this%bwd(isim)%nseis
          call pabort
       end if
    end do

    this%ndumps = ndumps_agreed
    this%windowlength = ndumps_agreed * dt_agreed

    ! Check whether the source time shift and stf are the same in all files
    source_shift_agreed_fwd = this%fwd(1)%source_shift_t
    stf_agreed_fwd = this%fwd(1)%stf
    stf_d_agreed_fwd = this%fwd(1)%stf_d
    amplitude_agreed_fwd = this%fwd(1)%amplitude

    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the others")' 
    fmtstring_stf = '("Inconsistency in forward simulations: ", A, " is different \'// &
                    '  in simulation ", I1, " vs the others")' 
    do isim = 1, this%nsim_fwd
       if (source_shift_agreed_fwd.ne.this%fwd(isim)%source_shift_t) then
          write(*,fmtstring) 'source time shift', isim, source_shift_agreed_fwd, &
                             this%fwd(isim)%source_shift_t
          call pabort
       end if
       if (any(abs(stf_agreed_fwd - this%fwd(isim)%stf).gt.1e-10)) then
           write(*,fmtstring) 'stf', isim
           call pabort
       end if
       if (any(abs(stf_d_agreed_fwd - this%fwd(isim)%stf_d).gt.1e-10)) then
           write(*,fmtstring) 'stf_d', isim
           call pabort
       end if
       if (amplitude_agreed_fwd.ne.this%fwd(isim)%amplitude) then
          write(*,fmtstring) 'source amplitude', isim, amplitude_agreed_fwd, &
                             this%fwd(isim)%amplitude
          call pabort
       end if
    end do

    this%timeshift_fwd = real(source_shift_agreed_fwd, kind=dp)
    allocate(this%stf_fwd(ndumps_agreed))
    allocate(this%stf_d_fwd(ndumps_agreed))
    this%stf_fwd = real(stf_agreed_fwd, kind=dp)
    this%stf_d_fwd = real(stf_d_agreed_fwd, kind=dp)
    this%amplitude_fwd = real(amplitude_agreed_fwd, kind=dp)

    if (this%nsim_bwd > 0) then

       source_shift_agreed_bwd = this%bwd(1)%source_shift_t
       stf_agreed_bwd = this%bwd(1)%stf
       stf_d_agreed_bwd = this%bwd(1)%stf_d 
       amplitude_agreed_bwd = this%bwd(1)%amplitude
       fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                   '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the others")' 
       fmtstring_stf = '("Inconsistency in backward simulations: ", A, " is different \'// &
                       '  in simulation ", I1, " vs the others")' 

       do isim = 1, this%nsim_bwd
          if (source_shift_agreed_bwd.ne.this%bwd(isim)%source_shift_t) then
             write(*,fmtstring) 'source time shift', isim, source_shift_agreed_bwd, &
                                this%bwd(isim)%source_shift_t
             call pabort
          end if
          if (any(abs(stf_agreed_bwd - this%bwd(isim)%stf).gt.1e-10)) then
              write(*,fmtstring) 'stf', isim
              call pabort
          end if
          if (any(abs(stf_d_agreed_bwd - this%bwd(isim)%stf_d).gt.1e-10)) then
              write(*,fmtstring) 'stf_d', isim
              call pabort
          end if
          if (amplitude_agreed_bwd.ne.this%bwd(isim)%amplitude) then
             write(*,fmtstring) 'source amplitude', isim, amplitude_agreed_bwd, &
                                this%bwd(isim)%amplitude
             call pabort
          end if
       end do

       this%timeshift_bwd = real(source_shift_agreed_bwd, kind=dp)
       allocate(this%stf_bwd(ndumps_agreed))
       allocate(this%stf_d_bwd(ndumps_agreed))
       this%stf_bwd = real(stf_agreed_bwd, kind=dp)
       this%stf_d_bwd = real(stf_d_agreed_bwd, kind=dp)
       this%amplitude_bwd = real(amplitude_agreed_bwd, kind=dp)

    endif

    this%dt = dt_agreed
    this%decimate_factor = nseis_agreed / ndumps_agreed
    this%nseis  = ndumps_agreed * this%decimate_factor        

    call flush(lu_out)

end subroutine check_consistency
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_fw_points(this, coordinates, source_params, coeffs)
    use finite_elem_mapping, only      : inside_element

    class(semdata_type)               :: this
    real(kind=dp), intent(in)         :: coordinates(:,:)
    type(src_param_type), intent(in)  :: source_params
    real(kind=dp)                     :: load_fw_points(this%ndumps, this%ndim, &
                                                        size(coordinates,2))

    real(kind=dp), intent(out), optional :: coeffs(3,size(coordinates,2)) ! for now: vp, vs, rho

    type(kdtree2_result), allocatable :: nextpoint(:)
    integer                           :: npoints, nnext_points
    integer                           :: pointid(size(coordinates,2))
    integer                           :: ipoint, inext_point, isim, iclockold, icp
    integer                           :: corner_point_ids(4), eltype(1), axis_int(1)
    logical                           :: axis
    integer, allocatable              :: gll_point_ids(:,:)
    integer                           :: id_elem
    real(kind=dp)                     :: corner_points(4,2)
    real(kind=dp)                     :: rotmesh_s(size(coordinates,2)), rotmesh_s_buff
    real(kind=dp)                     :: rotmesh_phi(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_z(size(coordinates,2))
    real(kind=dp)                     :: utemp(this%ndumps, this%ndim)
    real(kind=dp)                     :: coeff_buff(1)
    real(kind=dp)                     :: xi, eta

    
    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_fw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       call pabort 
    end if
    npoints = size(coordinates,2)
    
    if (trim(this%dump_type) == 'displ_only') then
        nnext_points = 6 ! 6, because this is the maximum valence in the mesh
        allocate(gll_point_ids(0:this%npol, 0:this%npol))
    else
        nnext_points = 1
    endif

    
    ! Rotate points to FWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,   &
                          coordinates*1d3,                              &
                          source_params%lon, source_params%colat)

    allocate(nextpoint(nnext_points))
    load_fw_points(:,:,:) = 0.0
    do ipoint = 1, npoints
        ! map points from outside earth to the surface:
        if (rotmesh_s(ipoint)**2 + rotmesh_z(ipoint)**2 > this%fwd(1)%planet_radius**2) then
           rotmesh_s_buff = rotmesh_s(ipoint) &
                               / (rotmesh_s(ipoint)**2 + rotmesh_z(ipoint)**2)**0.5d0 &
                               * this%fwd(1)%planet_radius
           rotmesh_z(ipoint) = rotmesh_z(ipoint) &
                               / (rotmesh_s(ipoint)**2 + rotmesh_z(ipoint)**2)**0.5d0 &
                               * this%fwd(1)%planet_radius
           rotmesh_s(ipoint) = rotmesh_s_buff
        endif

        iclockold = tick()
        call kdtree2_n_nearest( this%fwdtree,                           &
                                real([rotmesh_s(ipoint), rotmesh_z(ipoint)]), &
                                nn = nnext_points,                            &
                                results = nextpoint )
        iclockold = tick(id=id_kdtree, since=iclockold)
        
        pointid(ipoint) = nextpoint(1)%idx

        if (trim(this%dump_type) == 'displ_only') then
            iclockold = tick()
            do inext_point = 1, nnext_points
                corner_point_ids = this%fwdmesh%corner_point_ids(:, nextpoint(inext_point)%idx)
                eltype = this%fwdmesh%eltype(nextpoint(inext_point)%idx)
                
                do icp = 1, 4
                    corner_points(icp, 1) = this%fwdmesh%s(corner_point_ids(icp)+1)
                    corner_points(icp, 2) = this%fwdmesh%z(corner_point_ids(icp)+1)
                enddo                        
                ! test point to be inside, if so, exit
                if (inside_element(rotmesh_s(ipoint), rotmesh_z(ipoint), &
                                   corner_points, eltype(1), xi=xi, eta=eta, &
                                   tolerance=1d-3)) then
                    if (verbose > 1) then
                       write(6,*) 'coordinates= ', coordinates(:,ipoint)
                       write(6,*) 's, z       = ', rotmesh_s(ipoint), rotmesh_z(ipoint)
                       write(6,*) 'eltype     = ', eltype
                       write(6,*) 'xi, eta    = ', xi, eta
                       write(6,*) 'element id = ', nextpoint(inext_point)%idx
                    endif
                    exit
                endif
            enddo

            if (inext_point > nnext_points) then
               write(6,*) 'ERROR: element not found. (fwd)'
               write(6,*) '       Probably outside depth/distance range in the netcdf file?'
               write(6,*) '       Try increasing nnext_points in case this problem persists'
               write(6,*) 'coordinates= ', coordinates(:,ipoint)
               write(6,*) 's, z       = ', rotmesh_s(ipoint), rotmesh_z(ipoint)
               write(6,*) 'radius     = ', norm2([rotmesh_s(ipoint), rotmesh_z(ipoint)])
               do icp = 1, 4
                 write(6,*) 'cp: ', icp, ', s: ',   corner_points(icp, 1) 
                 write(6,*) 'cp: ', icp, ', z: ',   corner_points(icp, 2)
               enddo                        
               write(6,*) 'eltype     = ', eltype
               write(6,*) 'xi, eta    = ', xi, eta
               write(6,*) 'element id = ', nextpoint(inext_point)%idx
               call pabort(do_traceback = .false.)
               this%fwd(1)%count_error_pointoutside = this%fwd(1)%count_error_pointoutside + 1
               cycle
            endif

            id_elem = nextpoint(inext_point)%idx
         
            ! get gll points of spectral element
            gll_point_ids = -1
            if (verbose > 1) &
                write(6,*) 'element id = ', id_elem !nextpoint(inext_point)%idx
            call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
                                    varid  = this%fwd(1)%sem_mesh_varid, &
                                    start  = [1, 1, id_elem], &
                                    count  = [this%npol+1, this%npol+1, 1], &
                                    values = gll_point_ids))
            if (verbose > 1) &
                write(6,*) 'gll_point_ids = ', gll_point_ids(:,0)

            call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
                                    varid  = this%fwd(1)%axis_varid, &
                                    start  = [id_elem], &
                                    count  = [1], &
                                    values = axis_int))

            if (axis_int(1) == 1) then
               axis = .true.
            elseif (axis_int(1) == 0) then
               axis = .false.
            else
               call pabort
            endif

            if (verbose > 1) &
               write(6,*) 'axis = ', axis

            iclockold = tick(id=id_find_point_fwd, since=iclockold)
        endif
    
        ! Load model coefficients vp, vs and rho at point ipoint
        ! @ TODO : We need to store anisotropic model parameters
        !          in the wavefield netcdf files
        ! @ TODO : Ludwig is not sure whether [nextpoint(1)%idx]
        !          is the correct way to load coeffs from sem mesh
        ! Load coefficient vp
        call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
             varid  = this%fwd(1)%mesh_vp_varid,             &
             start  = [pointid(ipoint)],                     &
             count  = [1],                                   &
             values = coeff_buff))
        coeffs(1,ipoint) = coeff_buff(1)/1e3 ! convert to km/s
        ! Load coefficient vs
        call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
             varid  = this%fwd(1)%mesh_vs_varid,             &
             start  = [pointid(ipoint)],                     & 
             count  = [1],                                   &
             values = coeff_buff))
        coeffs(2,ipoint) = coeff_buff(1)/1e3 ! convert to km/s
        ! Load coefficient rho
        call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
             varid  = this%fwd(1)%mesh_rho_varid,            &
             start  = [pointid(ipoint)],                     & 
             count  = [1],                                   &
             values = coeff_buff))
        coeffs(3,ipoint) = coeff_buff(1)*1e9 ! convert to kg/(km^3)

        select case(trim(this%strain_type))
        case('straintensor_trace')    
           
           do isim = 1, this%nsim_fwd
              if (trim(this%dump_type) == 'displ_only') then
                 utemp = load_strain_point_interp(this%fwd(isim), gll_point_ids, &
                      xi, eta, this%strain_type, &
                      corner_points, eltype(1), axis, &
                      id_elem = id_elem) !/ this%fwd(isim)%amplitude
              else
                 utemp = load_strain_point(this%fwd(isim),      &
                      pointid(ipoint),     &
                      this%strain_type)  !/ this%fwd(isim)%amplitude
              endif
              
              iclockold = tick()
              
              load_fw_points(:, :, ipoint) = load_fw_points(:,:,ipoint)                   &
                   + utemp * azim_factor(rotmesh_phi(ipoint),     &
                   source_params%mij, isim, 1) 
              iclockold = tick(id=id_rotate, since=iclockold)
           end do !isim

        case('straintensor_full')

           do isim = 1, this%nsim_fwd

              if (trim(this%dump_type) == 'displ_only') then
                 utemp = load_strain_point_interp(this%fwd(isim), gll_point_ids, &
                      xi, eta, this%strain_type, &
                      corner_points, eltype(1), axis, &
                      id_elem = id_elem) !/ this%fwd(isim)%amplitude
              else
                 utemp = load_strain_point(this%fwd(isim),      &
                      pointid(ipoint),     &
                      this%strain_type)  !/ this%fwd(isim)%amplitude
              endif
              
              iclockold = tick()

              load_fw_points(:,1,ipoint) = load_fw_points(:,1,ipoint) &
                    + utemp(:,1) * azim_factor(rotmesh_phi(ipoint), source_params%mij, isim, 1) 
              load_fw_points(:,2,ipoint) = load_fw_points(:,2,ipoint) &
                    + utemp(:,2) * azim_factor(rotmesh_phi(ipoint), source_params%mij, isim, 1) 
              load_fw_points(:,3,ipoint) = load_fw_points(:,3,ipoint) &
                    + utemp(:,3) * azim_factor(rotmesh_phi(ipoint), source_params%mij, isim, 1) 
              load_fw_points(:,4,ipoint) = load_fw_points(:,4,ipoint) &
                    + utemp(:,4) * azim_factor(rotmesh_phi(ipoint), source_params%mij, isim, 2) 
              load_fw_points(:,5,ipoint) = load_fw_points(:,5,ipoint) &
                    + utemp(:,5) * azim_factor(rotmesh_phi(ipoint), source_params%mij, isim, 1) 
              load_fw_points(:,6,ipoint) = load_fw_points(:,6,ipoint) &
                    + utemp(:,6) * azim_factor(rotmesh_phi(ipoint), source_params%mij, isim, 2) 
              
              iclockold = tick(id=id_rotate, since=iclockold)
              
           end do !isim

           load_fw_points(:,:,ipoint) = rotate_symm_tensor_voigt_src_to_xyz(load_fw_points(:,:,ipoint), &
                                          source_params%lon, this%ndumps)

           load_fw_points(:,:,ipoint) = rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(load_fw_points(:,:,ipoint), &
                                          source_params%lon, source_params%colat, this%ndumps)



        end select

    end do !ipoint

end function load_fw_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine load_seismogram(this, receivers, src)
!< This function loads a seismogram 
   class(semdata_type)      :: this
   type(rec_param_type)     :: receivers(:)
   type(src_param_type)     :: src
   real(kind=dp)            :: seismogram_disp(this%ndumps)
   real(kind=dp)            :: seismogram_velo(this%ndumps)
   real(kind=sp)            :: utemp(this%ndumps,1,1)
   real(kind=dp)            :: Mij_scale(6), mij_prefact(4)
   integer                  :: reccomp, isurfelem, irec, isim, nrec

   if (.not.this%meshes_read) then
       print *, 'ERROR in load_seismogram(): Meshes have not been read yet'
       print *, 'Call read_meshes() before load_seismogram!'
       call pabort
   end if
      
   nrec = size(receivers)
   allocate(this%veloseis(this%ndumps, nrec))
   allocate(this%dispseis(this%ndumps, nrec))
   
   Mij_scale = src%mij / this%fwd(1)%amplitude

   write(lu_out, '(A, ES11.3)') '  Forward simulation amplitude: ', this%fwd(1)%amplitude
 
   do irec = 1, nrec
      write(lu_out, '(A,F8.4,A,F8.4)') '  Receiver theta: ', receivers(irec)%theta/deg2rad, &
                                       ', phi: ', receivers(irec)%phi/deg2rad
      write(lu_out, '(A)')             '                  Mij     Mij_scaled'
      write(lu_out, '(A,2(ES11.3))')   '  Mrr:       ', src%mij(1), mij_scale(1) 
      write(lu_out, '(A,2(ES11.3))')   '  Mtt:       ', src%mij(2), mij_scale(2)
      write(lu_out, '(A,2(ES11.3))')   '  Mpp:       ', src%mij(3), mij_scale(3)
      write(lu_out, '(A,2(ES11.3))')   '  Mrt:       ', src%mij(4), mij_scale(4)
      write(lu_out, '(A,2(ES11.3))')   '  Mrp:       ', src%mij(5), mij_scale(5)
      write(lu_out, '(A,2(ES11.3))')   '  Mtp:       ', src%mij(6), mij_scale(6)


      !print '(A,6(F8.5,/))',     'Mij_scale: ', Mij_scale
      select case(receivers(irec)%component)
      case('Z')
         mij_prefact(1) = Mij_scale(1)
         mij_prefact(2) = Mij_scale(2) + Mij_scale(3)
         mij_prefact(3) =   Mij_scale(4) * cos(receivers(irec)%phi) &
                          + Mij_scale(5) * sin(receivers(irec)%phi)
         mij_prefact(4) =  (Mij_scale(2) - Mij_scale(3)) * cos(2. * receivers(irec)%phi)  &
                          +           2. * Mij_scale(6) * sin(2. * receivers(irec)%phi) 
         reccomp = 1
      case('T')
         mij_prefact(1) = 0.0
         mij_prefact(2) = 0.0
         mij_prefact(3) = - Mij_scale(4) * sin(receivers(irec)%phi) &
                          + Mij_scale(5) * cos(receivers(irec)%phi)
         mij_prefact(4) =  (Mij_scale(3) - Mij_scale(2)) * sin(2. * receivers(irec)%phi) &
                          +           2. * Mij_scale(6)  * cos(2. * receivers(irec)%phi)
         reccomp = 2
      case('R')
         mij_prefact(1) = Mij_scale(1)
         mij_prefact(2) = Mij_scale(2) + Mij_scale(3)
         mij_prefact(3) =   Mij_scale(4) * cos(receivers(irec)%phi) &
                          + Mij_scale(5) * sin(receivers(irec)%phi)
         mij_prefact(4) =  (Mij_scale(2) - Mij_scale(3)) * cos(2. * receivers(irec)%phi)  &
                          +           2. * Mij_scale(6)  * sin(2. * receivers(irec)%phi) 
         reccomp = 3
      case default
         print *, 'ERROR: Unknown receiver component: ', receivers(irec)%component
         call pabort
      end select
      
      isurfelem = minloc( abs(this%fwdmesh%theta*deg2rad - receivers(irec)%theta), 1 )
      write(lu_out,'(A,F8.4,A,I5,A,F8.4)') &
                'Receiver with theta ', receivers(irec)%theta/deg2rad, &
                                    ' has element ', isurfelem, &
                                    ' with theta: ', this%fwdmesh%theta(isurfelem)
      
      seismogram_disp = 0.0
      seismogram_velo = 0.0

      write(lu_out,'(A,4(E12.4))') 'Mij prefactors: ', mij_prefact
      
      do isim = 1, this%nsim_fwd
         write(lu_out,'(A,I1,A,I5,A,I2,A,I6)') &
                'Sim: ', isim, ' Read element', isurfelem, &
                                               ', component: ', reccomp, ', no of samples:', this%ndumps
         ! Displacement seismogram
         !call check( nf90_get_var( ncid   = this%fwd(isim)%surf,        & 
         !                          varid  = this%fwd(isim)%seis_disp,   &
         !                          start  = [1, reccomp, isurfelem],    &
         !                          count  = [this%ndumps, 1, 1],        &
         !                          values = utemp) )
      
         call nc_getvar( ncid   = this%fwd(isim)%surf,        & 
                         varid  = this%fwd(isim)%seis_disp,   &
                         start  = [1, reccomp, isurfelem],    &
                         count  = [this%ndumps, 1, 1],        &
                         values = utemp) 
      
         seismogram_disp = real(utemp(:,1,1), kind=dp) * mij_prefact(isim) + seismogram_disp

         ! Velocity seismogram
         !call check( nf90_get_var( ncid   = this%fwd(isim)%surf,        & 
         !                          varid  = this%fwd(isim)%seis_velo,   &
         !                          start  = [1, reccomp, isurfelem],    &
         !                          count  = [this%ndumps, 1, 1],        &
         !                          values = utemp) )
      
         call nc_getvar( ncid   = this%fwd(isim)%surf,        & 
                         varid  = this%fwd(isim)%seis_velo,   &
                         start  = [1, reccomp, isurfelem],    &
                         count  = [this%ndumps, 1, 1],        &
                         values = utemp) 
      
         seismogram_velo = real(utemp(:,1,1), kind=dp) * mij_prefact(isim) + seismogram_velo

      end do

      this%dispseis(:, irec) = seismogram_disp(1:this%ndumps)
      this%veloseis(:, irec) = seismogram_velo(1:this%ndumps)

   end do
  
   call flush(lu_out)


end subroutine load_seismogram
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_bw_points(this, coordinates, receiver)
    use finite_elem_mapping, only      : inside_element

    class(semdata_type)                :: this
    real(kind=dp), intent(in)          :: coordinates(:,:)
    type(rec_param_type)               :: receiver
    real(kind=dp)                      :: load_bw_points(this%ndumps, this%ndim, &
                                                         size(coordinates,2))

    type(kdtree2_result), allocatable :: nextpoint(:)
    integer                           :: npoints, nnext_points, id_elem
    integer                           :: pointid(size(coordinates,2))
    integer                           :: ipoint, inext_point, iclockold, icp
    integer                           :: corner_point_ids(4), eltype(1), axis_int(1)
    logical                           :: axis
    integer, allocatable              :: gll_point_ids(:,:)
    real(kind=dp)                     :: corner_points(4,2)
    real(kind=dp)                     :: rotmesh_s(size(coordinates,2)), rotmesh_s_buff
    real(kind=dp)                     :: rotmesh_phi(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_z(size(coordinates,2))
    real(kind=dp)                     :: utemp(this%ndumps, this%ndim)
    real(kind=dp)                     :: xi, eta

    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_bw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       call pabort
    end if
    npoints = size(coordinates,2)

    if (trim(this%dump_type) == 'displ_only') then
        nnext_points = 6 ! 6, because this is the maximum valence in the mesh
        allocate(gll_point_ids(0:this%npol, 0:this%npol))
    else
        nnext_points = 1
    endif

    ! Rotate points to BWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,   &
                          coordinates*1e3,                              &
                          receiver%lon, receiver%colat)

    allocate(nextpoint(nnext_points))
    load_bw_points(:,:,:) = 0.0
    do ipoint = 1, npoints
        ! map points from outside earth to the surface:
        if (rotmesh_s(ipoint)**2 + rotmesh_z(ipoint)**2 > this%bwd(1)%planet_radius**2) then

           rotmesh_s_buff = rotmesh_s(ipoint) &
                               / (rotmesh_s(ipoint)**2 + rotmesh_z(ipoint)**2)**0.5d0 &
                               * this%bwd(1)%planet_radius
           rotmesh_z(ipoint) = rotmesh_z(ipoint) &
                               / (rotmesh_s(ipoint)**2 + rotmesh_z(ipoint)**2)**0.5d0 &
                               * this%bwd(1)%planet_radius
           rotmesh_s(ipoint) = rotmesh_s_buff
        endif
        iclockold = tick()
        call kdtree2_n_nearest( this%bwdtree,                           &
                                real([rotmesh_s(ipoint), rotmesh_z(ipoint)], kind=sp), &
                                nn = nnext_points,                            &
                                results = nextpoint )
        iclockold = tick(id=id_kdtree, since=iclockold)
    
        pointid(ipoint) = nextpoint(1)%idx
        if (trim(this%dump_type) == 'displ_only') then
            iclockold = tick()
          
            do inext_point = 1, nnext_points
                ! get cornerpoints of finite element
                corner_point_ids = this%bwdmesh%corner_point_ids(:, nextpoint(inext_point)%idx)
                eltype = this%bwdmesh%eltype(nextpoint(inext_point)%idx)
                
                do icp = 1, 4
                    corner_points(icp, 1) = this%bwdmesh%s(corner_point_ids(icp)+1)
                    corner_points(icp, 2) = this%bwdmesh%z(corner_point_ids(icp)+1)
                enddo                        

                ! test point to be inside, if so, exit
                if (inside_element(rotmesh_s(ipoint), rotmesh_z(ipoint), &
                                   corner_points, eltype(1), xi=xi, eta=eta, &
                                   tolerance=1d-3)) then
                    if (verbose > 1) then
                       write(6,*) 'coordinates= ', coordinates(:,ipoint)
                       write(6,*) 's, z       = ', rotmesh_s(ipoint), rotmesh_z(ipoint)
                       write(6,*) 'eltype     = ', eltype
                       write(6,*) 'xi, eta    = ', xi, eta
                       write(6,*) 'element id = ', nextpoint(inext_point)%idx
                    endif
                    exit
                endif
            enddo

            if (inext_point > nnext_points) then
               write(6,*) 'ERROR: element not found. (bwd)'
               write(6,*) '       Probably outside depth/distance range in the netcdf file?'
               write(6,*) '       Try increasing nnext_points in case this problem persists'
               write(6,*) 'radius     = ', norm2([rotmesh_s(ipoint), rotmesh_z(ipoint)])
               this%bwd(1)%count_error_pointoutside = this%bwd(1)%count_error_pointoutside + 1
               cycle
               call pabort(do_traceback = .false.)
            endif

            id_elem = nextpoint(inext_point)%idx

            ! get gll points of spectral element
            gll_point_ids = -1
            if (verbose > 1) &
                write(6,*) 'element id = ', nextpoint(inext_point)%idx
            call check(nf90_get_var(ncid   = this%bwd(1)%mesh,   &
                                    varid  = this%bwd(1)%sem_mesh_varid, &
                                    start  = [1, 1, nextpoint(inext_point)%idx], &
                                    count  = [this%npol+1, this%npol+1, 1], &
                                    values = gll_point_ids))
            if (verbose > 1) &
                write(6,*) 'gll_point_ids = ', gll_point_ids(:,0)

            call check(nf90_get_var(ncid   = this%bwd(1)%mesh,   &
                                    varid  = this%bwd(1)%axis_varid, &
                                    start  = [nextpoint(inext_point)%idx], &
                                    count  = [1], &
                                    values = axis_int))

            if (axis_int(1) == 1) then
               axis = .true.
            elseif (axis_int(1) == 0) then
               axis = .false.
            else
               call pabort
            endif

            if (verbose > 1) &
               write(6,*) 'axis = ', axis

            iclockold = tick(id=id_find_point_bwd, since=iclockold)
        endif
    
        select case(receiver%component)
        case('Z')
            if (trim(this%dump_type) == 'displ_only') then
               utemp = load_strain_point_interp(this%bwd(1), gll_point_ids,     &
                                                xi, eta, this%strain_type,      &
                                                corner_points, eltype(1), axis, &
                                                id_elem = id_elem)
            else
               utemp = load_strain_point(this%bwd(1), pointid(ipoint), this%strain_type)
            endif
            load_bw_points(:,:,ipoint) &
                =                               utemp / this%bwd(1)%amplitude
        case('R')
            if (trim(this%dump_type) == 'displ_only') then
               utemp = load_strain_point_interp(this%bwd(2), gll_point_ids,     &
                                                xi, eta, this%strain_type,      &
                                                corner_points, eltype(1), axis, &
                                                id_elem = id_elem)
            else
               utemp = load_strain_point(this%bwd(2), pointid(ipoint), this%strain_type)
            endif
            load_bw_points(:,:,ipoint) &
                =   dcos(rotmesh_phi(ipoint)) * utemp / this%bwd(2)%amplitude
        case('T')
            if (trim(this%dump_type) == 'displ_only') then
               utemp = load_strain_point_interp(this%bwd(2), gll_point_ids,     &
                                                xi, eta, this%strain_type,      &
                                                corner_points, eltype(1), axis, &
                                                id_elem = id_elem)
            else
               utemp = load_strain_point(this%bwd(2), pointid(ipoint), this%strain_type)
            endif
            load_bw_points(:,:,ipoint) &
                = - dsin(rotmesh_phi(ipoint)) * utemp / this%bwd(2)%amplitude 
        end select

        ! only need to rotate in case of vs
        if (this%strain_type.eq.'straintensor_full') then

           load_bw_points(:,:,ipoint) = rotate_symm_tensor_voigt_src_to_xyz(load_bw_points(:,:,ipoint), &
                                          receiver%lon, this%ndumps)

           load_bw_points(:,:,ipoint) = rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(load_bw_points(:,:,ipoint), &
                                          receiver%lon, receiver%colat, this%ndumps)
        end if

    end do !ipoint


end function load_bw_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_strain_point(sem_obj, pointid, strain_type)

    type(ncparamtype), intent(in)   :: sem_obj
    integer, intent(in)             :: pointid
    character(len=*), intent(in)    :: strain_type
    real(kind=dp), allocatable      :: load_strain_point(:,:)
    real(kind=dp), allocatable      :: strain_buff(:,:)

    integer                         :: start_chunk, iread, gll_to_read
    integer                         :: iclockold, iclockold_total
    integer                         :: status, istrainvar
    real(kind=sp), allocatable      :: utemp(:,:)
    real(kind=sp), allocatable      :: utemp_chunk(:,:,:)

    if (trim(sem_obj%dump_type) /= 'fullfields') then
        write(6,*) 'ERROR: trying to read strain from a file that was not'
        write(6,*) '       written with dump_type "fullfields"'
        write(6,*) sem_obj%dump_type
        call pabort(do_traceback=.false.)
    endif

    iclockold_total = tick()

    select case(strain_type)
    case('straintensor_trace')
        allocate(load_strain_point(sem_obj%ndumps, 1))
        allocate(utemp(sem_obj%ndumps, 1))
        allocate(utemp_chunk(sem_obj%chunk_gll, sem_obj%ndumps, 1))

        iclockold = tick()
        status = sem_obj%buffer%get(pointid, utemp)
        iclockold = tick(id=id_buffer, since=iclockold)

        if (status.ne.0) then
           start_chunk = ((pointid-1) / sem_obj%chunk_gll) * sem_obj%chunk_gll + 1

           ! Only read to last point, not further
           gll_to_read = min(sem_obj%chunk_gll, sem_obj%ngll + 1 - start_chunk)
           iclockold = tick()
           call nc_getvar( ncid   = sem_obj%snap,           & 
                           varid  = sem_obj%strainvarid(6), &
                           start  = [start_chunk, 1],       &
                           count  = [gll_to_read, sem_obj%ndumps], &
                           values = utemp_chunk(1:gll_to_read, :, 1)) 

           iclockold = tick(id=id_netcdf, since=iclockold)

           do iread = 0, sem_obj%chunk_gll - 1
               status = sem_obj%buffer%put(start_chunk + iread, utemp_chunk(iread+1,:,:))
           end do

           iclockold = tick(id=id_buffer, since=iclockold)

           load_strain_point(:,1) = real(utemp_chunk(pointid-start_chunk+1,:,1), kind=dp)
        else
           load_strain_point(:,1) = real(utemp(:,1), kind=dp)
        end if

    case('straintensor_full')
        allocate(utemp(sem_obj%ndumps, 6))
        allocate(utemp_chunk(sem_obj%chunk_gll, sem_obj%ndumps, 6))
        allocate(strain_buff(sem_obj%ndumps, 6))

        status = sem_obj%buffer%get(pointid, utemp)
        if (status.ne.0) then
            start_chunk = ((pointid-1) / sem_obj%chunk_gll) * sem_obj%chunk_gll + 1
            
            ! Only read to last point, not further
            gll_to_read = min(sem_obj%chunk_gll, sem_obj%ngll + 1 - start_chunk)
            do istrainvar = 1, 6

                if (sem_obj%strainvarid(istrainvar).eq.-1) then
                    strain_buff(:, istrainvar) = 0
                    cycle ! For monopole source which does not have this component.
                endif

                iclockold = tick()

                call nc_getvar( ncid   = sem_obj%snap,           & 
                                varid  = sem_obj%strainvarid(istrainvar), &
                                start  = [start_chunk, 1],       &
                                count  = [gll_to_read, sem_obj%ndumps], &
                                values = utemp_chunk(1:gll_to_read, :, istrainvar)) 

                iclockold = tick(id=id_netcdf, since=iclockold)
                strain_buff(:,istrainvar) &
                     = real(utemp_chunk(pointid-start_chunk+1, :, istrainvar), kind=dp)

            end do
            do iread = 0, gll_to_read - 1
                status = sem_obj%buffer%put(start_chunk + iread, utemp_chunk(iread+1,:,:))
            end do
        else
           strain_buff(:,:) = real(utemp, kind=dp)
        endif

        allocate(load_strain_point(sem_obj%ndumps, 6))
        ! transform strain to voigt mapping
        ! from:
        ! ['strain_dsus', 'strain_dsuz', 'strain_dpup', &
        !  'strain_dsup', 'strain_dzup', 'straintrace']
        ! to:
        ! dsus, dpup, dzuz, dzup, dsuz, dsup
        load_strain_point(:,1) = strain_buff(:,1)
        load_strain_point(:,2) = strain_buff(:,3)
        load_strain_point(:,3) = strain_buff(:,6) - strain_buff(:,1) - strain_buff(:,3)
        load_strain_point(:,4) = -strain_buff(:,5)
        load_strain_point(:,5) = strain_buff(:,2)
        load_strain_point(:,6) = -strain_buff(:,4)

    end select

    iclockold_total = tick(id=id_load_strain, since=iclockold)
end function load_strain_point
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_strain_point_interp(sem_obj, pointids, xi, eta, strain_type, nodes, &
                                  element_type, axis, id_elem)
    !< Calculates strain in element given by pointids, nodes. 
    !! Strain is then interpolated to the point defined by xi, eta.
    !! If parameter id_elem is present, it checks strain buffer, whether strain for this
    !! element has been calculated before. In this case, only the interpolation is done
    !! (order of magnitude faster)
    use sem_derivatives
    use spectral_basis, only : lagrange_interpol_2D_td

    type(ncparamtype), intent(in)   :: sem_obj
    integer,           intent(in)   :: pointids(0:sem_obj%npol, 0:sem_obj%npol) 
                                                    !< ID of GLL/GLI points in element of
                                                    !! interest
    real(kind=dp),     intent(in)   :: xi, eta      !< Coordinates at which to interpolate
                                                    !! strain
    character(len=*),  intent(in)   :: strain_type  !< Model parameter (decides on 
                                                    !! straintrace (vp) or full tensor (vs)
    real(kind=dp),     intent(in)   :: nodes(4,2)   !< Coordinates of element corner
                                                    !! points
    integer,           intent(in)   :: element_type !< Element type in the solver
    logical,           intent(in)   :: axis         !< Axis element or not 
    integer, optional, intent(in)   :: id_elem   !< ID of element to interpolate strain in
                                                 !! Giving this argument activates the
                                                 !! strain buffer. Omitting it restores
                                                 !! classic behaviour 
    real(kind=dp),     allocatable  :: load_strain_point_interp(:,:)

    logical                         :: use_strainbuffer
    integer                         :: start_chunk, iread, gll_to_read
    integer                         :: iclockold, status, idisplvar
    real(kind=sp), allocatable      :: utemp_chunk(:,:,:)
    real(kind=sp), allocatable      :: ubuff(:,:)
    real(kind=sp)                   :: utemp(1:sem_obj%ndumps, &
                                             0:sem_obj%npol,   &
                                             0:sem_obj%npol,   &
                                             3)
    real(kind=sp)                   :: strain(1:sem_obj%ndumps, &
                                              0:sem_obj%npol,   &
                                              0:sem_obj%npol,   &
                                              6)
    real(kind=sp)                   :: straintrace(1:sem_obj%ndumps, &
                                                   0:sem_obj%npol,   &
                                                   0:sem_obj%npol)
    real(kind=dp), allocatable      :: G(:,:), GT(:,:)
    real(kind=dp), allocatable      :: col_points_xi(:), col_points_eta(:)
    integer                         :: ipol, jpol, i, iclockold_total

    iclockold_total = tick()

    use_strainbuffer = present(id_elem)

    if (trim(sem_obj%dump_type) /= 'displ_only') then
        write(6,*) 'ERROR: trying to read interpolated strain from a file that was not'
        write(6,*) '       written with dump_type "displ_only"'
        call pabort
    endif

    if (axis) then
        G  = sem_obj%G2
        GT = sem_obj%G1T
        col_points_xi  = sem_obj%glj_points
        col_points_eta = sem_obj%gll_points
    else
        G  = sem_obj%G2
        GT = sem_obj%G2T
        col_points_xi  = sem_obj%gll_points
        col_points_eta = sem_obj%gll_points
    endif 


    if (use_strainbuffer) then

        if(id_elem<=0) then
            print *, 'id_elem is zero or smaller: ', id_elem
            call pabort()
        end if

        iclockold = tick()

        select case(strain_type)
        case('straintensor_trace')
            status = sem_obj%buffer_strain%get(id_elem, straintrace)
        case('straintensor_full')
            status = sem_obj%buffer_strain%get(id_elem, strain)
        case default
            status = - 1
        end select
        iclockold = tick(id=id_buffer, since=iclockold)
    else
        status = - 1
    end if

    if (status.ne.0) then

      allocate(utemp_chunk(sem_obj%chunk_gll, sem_obj%ndumps, 3))
      allocate(ubuff(sem_obj%ndumps, 3))

      ! load displacements from all GLL points
      do ipol = 0, sem_obj%npol
         do jpol = 0, sem_obj%npol

            iclockold = tick()
            status = sem_obj%buffer_disp%get(pointids(ipol,jpol), ubuff(:,:))
            iclockold = tick(id=id_buffer, since=iclockold)
            if (status.ne.0) then
               start_chunk &
                    = (pointids(ipol,jpol) / sem_obj%chunk_gll) * sem_obj%chunk_gll + 1
               ! Only read to last point, not further
               gll_to_read = min(sem_obj%chunk_gll, sem_obj%ngll + 1 - start_chunk)

               do idisplvar = 1, 3

                   if (sem_obj%displvarid(idisplvar).eq.-1) then
                       utemp(:, ipol, jpol, idisplvar) = 0
                       cycle ! For monopole source which does not have this component.
                   endif

                   iclockold = tick()
                   call check(nf90_get_var( ncid   = sem_obj%snap,           & 
                                            varid  = sem_obj%displvarid(idisplvar), &
                                            start  = [start_chunk, 1],       &
                                            count  = [gll_to_read, sem_obj%ndumps], &
                                            values = utemp_chunk(1:gll_to_read, :, idisplvar)))

                   iclockold = tick(id=id_netcdf, since=iclockold)
                   utemp(:,ipol,jpol, idisplvar) &
                        = utemp_chunk(pointids(ipol,jpol) - start_chunk + 2,:,idisplvar)
               enddo
               do iread = 0, sem_obj%chunk_gll - 1
                   status = sem_obj%buffer_disp%put(start_chunk + iread - 1, &
                                                    utemp_chunk(iread+1,:,:) )
               end do
            else
               utemp(:,ipol,jpol,:) = ubuff(:,:)
            endif
         enddo
      enddo

      iclockold = tick()
      select case(strain_type)
      case('straintensor_trace')
          ! compute straintrace
          if (sem_obj%excitation_type == 'monopole') then
              straintrace = real(straintrace_monopole(real(utemp, kind=dp), G, GT, col_points_xi, &
                                                 col_points_eta, sem_obj%npol, &
                                                 sem_obj%ndumps, nodes, element_type, axis), kind=sp)

          elseif (sem_obj%excitation_type == 'dipole') then
              straintrace = real(straintrace_dipole(real(utemp, kind=dp), G, GT, col_points_xi, &
                                               col_points_eta, sem_obj%npol, sem_obj%ndumps, &
                                               nodes, element_type, axis), kind=sp)

          elseif (sem_obj%excitation_type == 'quadpole') then
              straintrace = real(straintrace_quadpole(real(utemp, kind=dp), G, GT, col_points_xi, &
                                                 col_points_eta, sem_obj%npol, &
                                                 sem_obj%ndumps, nodes, element_type, axis), kind=sp)
          else
              print *, 'ERROR: unknown excitation_type: ', sem_obj%excitation_type
              call pabort
          endif

          iclockold = tick(id=id_calc_strain, since=iclockold)
          if (use_strainbuffer) & 
              status = sem_obj%buffer_strain%put(id_elem, straintrace)

      case('straintensor_full')
          ! compute full strain tensor
          if (sem_obj%excitation_type == 'monopole') then
              strain = real(strain_monopole(real(utemp, kind=dp), G, GT, col_points_xi, &
                                       col_points_eta, sem_obj%npol, sem_obj%ndumps, nodes, &
                                       element_type, axis), kind=sp)

          elseif (sem_obj%excitation_type == 'dipole') then
              strain = real(strain_dipole(real(utemp, kind=dp), G, GT, col_points_xi, &
                                     col_points_eta, sem_obj%npol, sem_obj%ndumps, nodes, &
                                     element_type, axis), kind=sp)

          elseif (sem_obj%excitation_type == 'quadpole') then
              strain = real(strain_quadpole(real(utemp, kind=dp), G, GT, col_points_xi, &
                                       col_points_eta, sem_obj%npol, sem_obj%ndumps, nodes, &
                                       element_type, axis), kind=sp)
          else
              print *, 'ERROR: unknown excitation_type: ', sem_obj%excitation_type
              call pabort
          endif
          
          iclockold = tick(id=id_calc_strain, since=iclockold)

          if (use_strainbuffer) & 
              status = sem_obj%buffer_strain%put(id_elem, strain)

      end select
    
    endif ! Element not found in buffer
    
    select case(strain_type)
    case('straintensor_trace')
        allocate(load_strain_point_interp(sem_obj%ndumps, 1))
        load_strain_point_interp(:, 1) &
            = lagrange_interpol_2D_td(col_points_xi, col_points_eta, &
                                      real(straintrace(:,:,:), kind=dp), xi, eta)

        iclockold = tick(id=id_lagrange, since=iclockold)

    case('straintensor_full')
        allocate(load_strain_point_interp(sem_obj%ndumps, 6))
        do i = 1, 6
            load_strain_point_interp(:, i) &
                = lagrange_interpol_2D_td(col_points_xi, col_points_eta, &
                                          real(strain(:,:,:,i), kind=dp), xi, eta)
        enddo

        iclockold = tick(id=id_lagrange, since=iclockold)

        !@TODO for consistency with SOLVER output
        load_strain_point_interp(:, 4) = - load_strain_point_interp(:, 4) 
        load_strain_point_interp(:, 6) = - load_strain_point_interp(:, 6) 

    end select

    iclockold_total = tick(id=id_load_strain, since=iclockold_total)

end function load_strain_point_interp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine build_kdtree(this)
    class(semdata_type)        :: this
    real(kind=sp), allocatable :: mesh(:,:)

    if (.not.this%meshes_read) then
        print *, 'ERROR in build_kdtree(): Meshes have not been read yet'
        print *, 'Call read_meshes() before build_kdtree!'
        call pabort
    end if

    write(lu_out,*) ' Reshaping mesh variables'
    call flush(lu_out)

    if (trim(this%dump_type) == 'displ_only') then
        allocate(mesh(2, this%fwdmesh%nelem)) ! midpoints only
        mesh = transpose(reshape([this%fwdmesh%s_mp, this%fwdmesh%z_mp], &
                                 [this%fwdmesh%nelem, 2]))
    else
        allocate(mesh(2, this%fwdmesh%npoints))
        mesh = transpose(reshape([this%fwdmesh%s, this%fwdmesh%z],       &
                                 [this%fwdmesh%npoints, 2]))
    endif


    write(lu_out,*) ' Building forward KD-Tree'
    call flush(lu_out)
    ! KDtree in forward field
    this%fwdtree => kdtree2_create(mesh,              &
                                   dim = 2,           &
                                   sort = .true.,     &
                                   rearrange = .true.)
    
    deallocate(mesh)                           

    ! KDtree in backward field
    if (this%nsim_bwd > 0) then

        if (trim(this%dump_type) == 'displ_only') then
            allocate(mesh(2, this%bwdmesh%nelem)) ! midpoints only
            mesh = transpose(reshape([this%bwdmesh%s_mp, this%bwdmesh%z_mp], &
                                     [this%bwdmesh%nelem, 2]))
        else
            allocate(mesh(2, this%bwdmesh%npoints))
            mesh = transpose(reshape([this%bwdmesh%s, this%bwdmesh%z],       &
                                     [this%bwdmesh%npoints, 2]))
        endif

        write(lu_out,*) ' Building backward KD-Tree'
        call flush(lu_out)
        this%bwdtree => kdtree2_create(mesh,              &
                                       dim = 2,           &
                                       sort = .true.,     &
                                       rearrange = .true.)
        deallocate(mesh)                           
    endif

    call flush(lu_out)

end subroutine build_kdtree
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_meshes(this)
    use netcdf
    use spectral_basis, only : zelegl, zemngl2, &
                               def_lagrange_derivs_gll, def_lagrange_derivs_glj

    class(semdata_type)        :: this
    integer                    :: ncvarid_mesh_s, ncvarid_mesh_z
    integer                    :: surfdimid, ncvarid_theta
    integer                    :: isim
    logical                    :: mesherror
    real(kind=sp), allocatable :: theta(:)
   
    if (.not.this%files_open) then
        print *, 'ERROR in read_meshes(): Files have not been opened!'
        print *, 'Call open_files() before read_meshes()'
        call pabort
    end if

    ! Forward SEM mesh
    write(lu_out,*) '  Read SEM mesh from first forward simulation'
    
    call nc_read_att_int(this%fwdmesh%npoints, 'npoints', this%fwd(1))
    
    do isim = 1, this%nsim_fwd
       this%fwd(isim)%ngll = this%fwdmesh%npoints
    end do

    if (trim(this%dump_type) == 'displ_only') then
        call nc_read_att_int(this%fwdmesh%nelem, 'nelem_kwf_global', this%fwd(1))
        
        allocate(this%fwdmesh%s_mp(this%fwdmesh%nelem))
        allocate(this%fwdmesh%z_mp(this%fwdmesh%nelem))

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "fem_mesh",            &
                      varid = this%fwd(1)%fem_mesh_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "sem_mesh",            &
                      varid = this%fwd(1)%sem_mesh_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "eltype",              &
                      varid = this%fwd(1)%eltype_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "axis",              &
                      varid = this%fwd(1)%axis_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_S",              &
                      varid = this%fwd(1)%mesh_s_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_Z",              &
                      varid = this%fwd(1)%mesh_z_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_rho",              &
                      varid = this%fwd(1)%mesh_rho_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_vp",              &
                      varid = this%fwd(1)%mesh_vp_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_vs",              &
                      varid = this%fwd(1)%mesh_vs_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_lambda",              &
                      varid = this%fwd(1)%mesh_lambda_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_mu",              &
                      varid = this%fwd(1)%mesh_mu_varid)
           
        call  getvarid( ncid  = this%fwd(1)%mesh, &
                        name  = "mp_mesh_S", &
                        varid = ncvarid_mesh_s) 

        call  getvarid( ncid  = this%fwd(1)%mesh, &
                        name  = "mp_mesh_Z", &
                        varid = ncvarid_mesh_z) 

        call nc_getvar(ncid   = this%fwd(1)%mesh,       &
                       varid  = ncvarid_mesh_s,         &
                       start  = 1,                    &
                       count  = this%fwdmesh%nelem, &
                       values = this%fwdmesh%s_mp ) 

        call nc_getvar(ncid   = this%fwd(1)%mesh,       &
                       varid  = ncvarid_mesh_z,         &
                       start  = 1,                    &
                       count  = this%fwdmesh%nelem, &
                       values = this%fwdmesh%z_mp ) 

!       Only executed, if compiled for the KERNER. This block is skipped for the RDBM!
        allocate(this%fwdmesh%corner_point_ids(4, this%fwdmesh%nelem))
        allocate(this%fwdmesh%s(this%fwdmesh%npoints))
        allocate(this%fwdmesh%z(this%fwdmesh%npoints))
        allocate(this%fwdmesh%eltype(this%fwdmesh%nelem))

        call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
                                varid  = this%fwd(1)%fem_mesh_varid, &
                                start  = [1, 1], &
                                count  = [4, this%fwdmesh%nelem], &
                                values = this%fwdmesh%corner_point_ids))
        
        call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
                                varid  = this%fwd(1)%eltype_varid, &
                                start  = [1], &
                                count  = [this%fwdmesh%nelem], &
                                values = this%fwdmesh%eltype))
        
        call nc_getvar(ncid   = this%fwd(1)%mesh,           &
                       varid  = this%fwd(1)%mesh_s_varid,   &
                       start  = 1,                          &
                       count  = this%fwdmesh%npoints,       &
                       values = this%fwdmesh%s ) 

        call nc_getvar(ncid   = this%fwd(1)%mesh,           &
                       varid  = this%fwd(1)%mesh_z_varid,   &
                       start  = 1,                          &
                       count  = this%fwdmesh%npoints,       &
                       values = this%fwdmesh%z ) 
    
    else
        allocate(this%fwdmesh%s(this%fwdmesh%npoints))
        allocate(this%fwdmesh%z(this%fwdmesh%npoints))

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_rho",              &
                      varid = this%fwd(1)%mesh_rho_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_vp",              &
                      varid = this%fwd(1)%mesh_vp_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_vs",              &
                      varid = this%fwd(1)%mesh_vs_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_lambda",              &
                      varid = this%fwd(1)%mesh_lambda_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_mu",              &
                      varid = this%fwd(1)%mesh_mu_varid)
           
        call  getvarid( ncid  = this%fwd(1)%mesh, &
                        name  = "mesh_S", &
                        varid = ncvarid_mesh_s) 

        call  getvarid( ncid  = this%fwd(1)%mesh, &
                        name  = "mesh_Z", &
                        varid = ncvarid_mesh_z) 

        call nc_getvar(ncid   = this%fwd(1)%mesh,       &
                       varid  = ncvarid_mesh_s,         &
                       start  = 1,                    &
                       count  = this%fwdmesh%npoints, &
                       values = this%fwdmesh%s ) 

        call nc_getvar(ncid   = this%fwd(1)%mesh,       &
                       varid  = ncvarid_mesh_z,         &
                       start  = 1,                    &
                       count  = this%fwdmesh%npoints, &
                       values = this%fwdmesh%z ) 
    endif

   
   
    ! Backward SEM mesh                     
    if (this%nsim_bwd > 0) then
        write(lu_out,*) 'Read SEM mesh from first backward simulation'
        
        call nc_read_att_int(this%bwdmesh%npoints, 'npoints', this%bwd(1))
        
        do isim = 1, this%nsim_bwd
           this%bwd(isim)%ngll = this%bwdmesh%npoints
        end do

        if (trim(this%dump_type) == 'displ_only') then
            call nc_read_att_int(this%bwdmesh%nelem, 'nelem_kwf_global', this%fwd(1))
            
            allocate(this%bwdmesh%s_mp(this%fwdmesh%nelem))
            allocate(this%bwdmesh%z_mp(this%fwdmesh%nelem))

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "fem_mesh",            &
                          varid = this%bwd(1)%fem_mesh_varid)

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "sem_mesh",            &
                          varid = this%bwd(1)%sem_mesh_varid)

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "eltype",              &
                          varid = this%bwd(1)%eltype_varid)

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "axis",              &
                          varid = this%bwd(1)%axis_varid)

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "mesh_S",              &
                          varid = this%bwd(1)%mesh_s_varid)

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "mesh_Z",              &
                          varid = this%bwd(1)%mesh_z_varid)

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "mesh_rho",              &
                          varid = this%bwd(1)%mesh_rho_varid)

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "mesh_lambda",              &
                          varid = this%bwd(1)%mesh_lambda_varid)

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "mesh_mu",              &
                          varid = this%bwd(1)%mesh_mu_varid)
            
            call  getvarid( ncid  = this%bwd(1)%mesh, &
                            name  = "mp_mesh_S", &
                            varid = ncvarid_mesh_s) 
            call  getvarid( ncid  = this%bwd(1)%mesh, &
                            name  = "mp_mesh_Z", &
                            varid = ncvarid_mesh_z) 

            call nc_getvar(ncid   = this%bwd(1)%mesh, &
                           varid  = ncvarid_mesh_s,   &
                           start  = 1,                &
                           count  = this%bwdmesh%nelem, &
                           values = this%bwdmesh%s_mp ) 
            call nc_getvar(ncid   = this%bwd(1)%mesh, &
                           varid  = ncvarid_mesh_z,   &
                           start  = 1,                &
                           count  = this%bwdmesh%nelem, &
                           values = this%bwdmesh%z_mp ) 

!           Only executed, if compiled for the KERNER. This block is skipped for the RDBM!
            allocate(this%bwdmesh%corner_point_ids(4, this%bwdmesh%nelem))
            allocate(this%bwdmesh%s(this%bwdmesh%npoints))
            allocate(this%bwdmesh%z(this%bwdmesh%npoints))
            allocate(this%bwdmesh%eltype(this%bwdmesh%nelem))

            call check(nf90_get_var(ncid   = this%bwd(1)%mesh,   &
                                    varid  = this%bwd(1)%fem_mesh_varid, &
                                    start  = [1, 1], &
                                    count  = [4, this%bwdmesh%nelem], &
                                    values = this%bwdmesh%corner_point_ids))
            
            call check(nf90_get_var(ncid   = this%bwd(1)%mesh,   &
                                    varid  = this%bwd(1)%eltype_varid, &
                                    start  = [1], &
                                    count  = [this%bwdmesh%nelem], &
                                    values = this%bwdmesh%eltype))
            
            call nc_getvar(ncid   = this%bwd(1)%mesh,            &
                           varid  = this%bwd(1)%mesh_s_varid,    &
                           start  = 1,                           &
                           count  = this%bwdmesh%npoints,        &
                           values = this%bwdmesh%s ) 
            call nc_getvar(ncid   = this%bwd(1)%mesh,            &
                           varid  = this%bwd(1)%mesh_z_varid,    &
                           start  = 1,                           &
                           count  = this%bwdmesh%npoints,        &
                           values = this%bwdmesh%z ) 

        else

            allocate(this%bwdmesh%s(this%bwdmesh%npoints))
            allocate(this%bwdmesh%z(this%bwdmesh%npoints))
            
            call  getvarid( ncid  = this%bwd(1)%mesh, &
                            name  = "mesh_S", &
                            varid = ncvarid_mesh_s) 
            call  getvarid( ncid  = this%bwd(1)%mesh, &
                            name  = "mesh_Z", &
                            varid = ncvarid_mesh_z) 

            call nc_getvar(ncid   = this%bwd(1)%mesh, &
                           varid  = ncvarid_mesh_s,   &
                           start  = 1,                &
                           count  = this%bwdmesh%npoints, &
                           values = this%bwdmesh%s ) 
            call nc_getvar(ncid   = this%bwd(1)%mesh, &
                           varid  = ncvarid_mesh_z,   &
                           start  = 1,                &
                           count  = this%bwdmesh%npoints, &
                           values = this%bwdmesh%z ) 
        endif
    endif

    ! Read surface element theta
    ! Forward mesh
    call check( nf90_inq_dimid( ncid  = this%fwd(1)%surf, &
                                name  = 'surf_elems',     &
                                dimid = surfdimid) ) 

    call check( nf90_inquire_dimension(ncid  = this%fwd(1)%surf,        & 
                                       dimid = surfdimid,               &
                                       len   = this%fwdmesh%nsurfelem) )

    call  getvarid( ncid  = this%fwd(1)%surf, &
                    name  = "elem_theta",     &
                    varid = ncvarid_theta) 

    allocate( this%fwdmesh%theta(this%fwdmesh%nsurfelem) )
    allocate( theta(this%fwdmesh%nsurfelem) )
    call nc_getvar( ncid   = this%fwd(1)%surf,   &
                    varid  = ncvarid_theta,          &
                    start  = 1,                      & 
                    count  = this%fwdmesh%nsurfelem, &
                    values = theta                    )
    this%fwdmesh%theta = real(theta, kind=dp)
    
    ! Backward mesh
    if (this%nsim_bwd > 0) then

       call check( nf90_inq_dimid( ncid  = this%bwd(1)%surf, &
                                   name  = 'surf_elems',     &
                                   dimid = surfdimid) ) 

       call check( nf90_inquire_dimension(ncid  = this%bwd(1)%surf,        & 
                                          dimid = surfdimid,               &
                                          len   = this%bwdmesh%nsurfelem) )

       call  getvarid(ncid  = this%bwd(1)%surf, &
                      name  = "elem_theta",     &
                      varid = ncvarid_theta) 

      allocate( this%bwdmesh%theta(this%bwdmesh%nsurfelem) )
      ! sure that fwd is correct here??
      call nc_getvar( ncid   = this%bwd(1)%surf,   &
                      varid  = ncvarid_theta,          &
                      start  = 1,                      & 
                      count  = this%bwdmesh%nsurfelem, &
                      values = theta                    )
      this%bwdmesh%theta = real(theta, kind=dp)
   endif
                             
    ! define terms needed to compute gradient
    if (trim(this%dump_type) == 'displ_only') then
        allocate(this%G1(0:this%npol,0:this%npol))
        allocate(this%G1T(0:this%npol,0:this%npol))
        allocate(this%G2(0:this%npol,0:this%npol))
        allocate(this%G2T(0:this%npol,0:this%npol))
        allocate(this%G0(0:this%npol))

        allocate(this%gll_points(0:this%npol))
        allocate(this%glj_points(0:this%npol))

        this%gll_points = zelegl(this%npol)
        this%glj_points = zemngl2(this%npol)

        this%G1 = def_lagrange_derivs_glj(this%npol, this%G0)
        this%G2 = def_lagrange_derivs_gll(this%npol)

        this%G1T = transpose(this%G1)
        this%G2T = transpose(this%G2)

        do isim = 1, this%nsim_fwd
           allocate(this%fwd(isim)%gll_points(0:this%npol))
           allocate(this%fwd(isim)%glj_points(0:this%npol))

           this%fwd(isim)%gll_points = this%gll_points
           this%fwd(isim)%glj_points = this%glj_points

           allocate(this%fwd(isim)%G1(0:this%npol,0:this%npol))
           allocate(this%fwd(isim)%G1T(0:this%npol,0:this%npol))
           allocate(this%fwd(isim)%G2(0:this%npol,0:this%npol))
           allocate(this%fwd(isim)%G2T(0:this%npol,0:this%npol))
           allocate(this%fwd(isim)%G0(0:this%npol))

           this%fwd(isim)%G1 = this%G1
           this%fwd(isim)%G2 = this%G2
           this%fwd(isim)%G1T = this%G1T
           this%fwd(isim)%G2T = this%G2T
           this%fwd(isim)%G0 = this%G0
        end do

        do isim = 1, this%nsim_bwd
           allocate(this%bwd(isim)%gll_points(0:this%npol))
           allocate(this%bwd(isim)%glj_points(0:this%npol))
           
           this%bwd(isim)%gll_points = this%gll_points
           this%bwd(isim)%glj_points = this%glj_points
           
           allocate(this%bwd(isim)%G1(0:this%npol,0:this%npol))
           allocate(this%bwd(isim)%G1T(0:this%npol,0:this%npol))
           allocate(this%bwd(isim)%G2(0:this%npol,0:this%npol))
           allocate(this%bwd(isim)%G2T(0:this%npol,0:this%npol))
           allocate(this%bwd(isim)%G0(0:this%npol))
           
           this%bwd(isim)%G1 = this%G1
           this%bwd(isim)%G2 = this%G2
           this%bwd(isim)%G1T = this%G1T
           this%bwd(isim)%G2T = this%G2T
           this%bwd(isim)%G0 = this%G0
        end do

    endif

    ! Mesh sanity checks
    mesherror = .false.
    if (allocated(this%fwdmesh%s)) then
       if (maxval(abs(this%fwdmesh%s)) > 1e32) then
          write(*,*) 'Maximum value of S in the forward mesh is unreasonably large'
          write(*,*) 'maxval(S): ', this%fwdmesh%s(maxloc(abs(this%fwdmesh%s))), ' m'
          mesherror = .true.
       end if
    end if
    if (allocated(this%fwdmesh%z)) then
       if (maxval(abs(this%fwdmesh%z)) > 1e32) then
          write(*,*) 'Maximum value of Z in the forward mesh is unreasonably large'
          write(*,*) 'maxval(Z): ', this%fwdmesh%z(maxloc(abs(this%fwdmesh%z))), ' m'
          mesherror = .true.
       end if
    end if
    if (allocated(this%fwdmesh%s_mp)) then
       if (maxval(abs(this%fwdmesh%s_mp)) > 1e32) then
          write(*,*) 'Maximum value of S_mp in the forward mesh is unreasonably large'
          write(*,*) 'maxval(S_mp): ', this%fwdmesh%s_mp(maxloc(abs(this%fwdmesh%s_mp))), ' m'
          mesherror = .true.
       end if
    end if
    if (allocated(this%fwdmesh%z_mp)) then
       if (maxval(abs(this%fwdmesh%z_mp)) > 1e32) then
          write(*,*) 'Maximum value of Z_mp in the forward mesh is unreasonably large'
          write(*,*) 'maxval(Z_mp): ', this%fwdmesh%z_mp(maxloc(abs(this%fwdmesh%z_mp))), ' m'
          mesherror = .true.
       end if
    end if
    if (maxval(this%fwdmesh%theta).gt.180.0) then
       write(*,*) 'Maximum value of theta in the backward mesh is larger than 180°'
       write(*,*) 'maxval(theta): ', this%fwdmesh%theta(maxloc(abs(this%fwdmesh%theta)))
       write(*,*) 'maxloc(theta): ', maxloc(abs(this%fwdmesh%theta))
       mesherror = .true.
    end if

    if (this%nsim_bwd > 0) then
       if (allocated(this%bwdmesh%s)) then
          if (maxval(abs(this%bwdmesh%s)) > 1e32) then
             write(*,*) 'Maximum value of S in the backward mesh is unreasonably large'
             write(*,*) 'maxval(S): ', this%bwdmesh%s(maxloc(abs(this%bwdmesh%s))), ' m'
             mesherror = .true.
          end if
       end if
       if (allocated(this%bwdmesh%z)) then
          if (maxval(abs(this%bwdmesh%z)) > 1e32) then
             write(*,*) 'Maximum value of Z in the backward mesh is unreasonably large'
             write(*,*) 'maxval(Z): ', this%bwdmesh%z(maxloc(abs(this%bwdmesh%z))), ' m'
             mesherror = .true.
          end if
       end if
       if (allocated(this%bwdmesh%s_mp)) then
          if (maxval(abs(this%bwdmesh%s_mp)) > 1e32) then
             write(*,*) 'Maximum value of S_mp in the backward mesh is unreasonably large'
             write(*,*) 'maxval(S_mp): ', this%bwdmesh%s_mp(maxloc(abs(this%bwdmesh%s_mp))), ' m'
             mesherror = .true.
          end if
       end if
       if (allocated(this%bwdmesh%z_mp)) then
          if (maxval(abs(this%bwdmesh%z_mp)) > 1e32) then
             write(*,*) 'Maximum value of Z_mp in the backward mesh is unreasonably large'
             write(*,*) 'maxval(Z_mp): ', this%bwdmesh%z_mp(maxloc(abs(this%bwdmesh%z_mp))), ' m'
             mesherror = .true.
          end if
       end if
       if (maxval(this%bwdmesh%theta).gt.180.0) then
          write(*,*) 'Maximum value of theta in the backward mesh is larger than 180°'
          write(*,*) 'maxval(theta): ', this%bwdmesh%theta(maxloc(abs(this%bwdmesh%theta)))
          write(*,*) 'maxloc(theta): ', maxloc(abs(this%bwdmesh%theta))
          mesherror = .true.
       end if
    endif


    if (mesherror) then
       write(*,*) 'ERROR: One or more mesh errors found!'
       call pabort
    end if

                                   
   this%meshes_read = .true.

   write(lu_out, *) 'Forward and backward SEM mesh reading succeeded'
   call flush(lu_out)

end subroutine read_meshes
!-----------------------------------------------------------------------------------------
 
!-----------------------------------------------------------------------------------------
!> Read NetCDF attribute of type Integer
subroutine nc_read_att_int(attribute_value, attribute_name, nc)
  character(len=*),  intent(in)     :: attribute_name
  integer, intent(out)              :: attribute_value
  type(ncparamtype), intent(in)     :: nc
  integer                           :: status

  status = nf90_get_att(nc%ncid, NF90_GLOBAL, attribute_name, attribute_value)
  if (status.ne.NF90_NOERR) then
      write(6,*) 'Could not find attribute ', trim(attribute_name)
      write(6,*) ' in NetCDF file ', trim(nc%meshdir), '/Data/ordered_output.nc4'
      write(6,*) ' with NCID: ', nc%ncid
      call pabort
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
      write(6,*) 'Could not find attribute ', trim(attribute_name)
      write(6,*) ' in NetCDF file ', trim(nc%meshdir), '/Data/ordered_output.nc4'
      write(6,*) ' with NCID: ', nc%ncid
      call pabort 
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
      write(6,*) 'Could not find attribute ', trim(attribute_name)
      write(6,*) ' in NetCDF file ', trim(nc%meshdir), '/Data/ordered_output.nc4'
      write(6,*) ' with NCID: ', nc%ncid
      call pabort
  end if
end subroutine nc_read_att_real
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Read NetCDF attribute of type Double
subroutine nc_read_att_dble(attribute_value, attribute_name, nc)
  character(len=*),  intent(in)     :: attribute_name
  real(kind=dp), intent(out)        :: attribute_value
  type(ncparamtype), intent(in)     :: nc
  integer                           :: status

  status = nf90_get_att(nc%ncid, NF90_GLOBAL, attribute_name, attribute_value)
  if (status.ne.NF90_NOERR) then
      write(6,*) 'Could not find attribute ', trim(attribute_name)
      write(6,*) ' in NetCDF file ', trim(nc%meshdir), '/Data/ordered_output.nc4'
      write(6,*) ' with NCID: ', nc%ncid
      call pabort
  end if
end subroutine nc_read_att_dble
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
