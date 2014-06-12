!=========================================================================================
module readfields
    use global_parameters, only            : sp, dp, pi, deg2rad, rad2deg, verbose, lu_out, &
                                             myrank, id_buffer, id_netcdf, id_rotate
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
    public                                :: semdata_type

    type meshtype
        real(kind=sp), allocatable        :: s(:), z(:)
        integer                           :: npoints, nelem
        real(kind=dp), allocatable        :: theta(:)
        integer                           :: nsurfelem
    end type

    type ncparamtype
        integer                           :: ncid
        integer                           :: snap, surf, mesh, seis  ! Group IDs
        integer                           :: strainvarid(6)          ! Variable IDs
        integer                           :: displvarid(3)           ! Variable IDs
        integer                           :: seis_disp, seis_velo    ! Variable IDs
        integer                           :: stf_varid               ! Variable IDs
        integer                           :: fem_mesh_varid          ! Variable IDs
        integer                           :: sem_mesh_varid          ! Variable IDs
        integer                           :: eltype_varid            ! Variable IDs
        integer                           :: mesh_s_varid            ! Variable IDs
        integer                           :: mesh_z_varid            ! Variable IDs
        integer                           :: chunk_gll
        character(len=200)                :: meshdir
        character(len=12)                 :: dump_type
        integer                           :: ndumps, nseis, ngll, npol
        integer                           :: source_shift_samples    
        real(kind=dp)                     :: source_shift_t
        real(kind=sp), allocatable        :: stf(:)
        type(buffer_type), allocatable    :: buffer(:)
        real(kind=dp)                     :: dt
        real(kind=dp)                     :: amplitude
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
        
        character(len=4)                   :: model_param   !< Parameter for which to calculate kernel
        integer                            :: ndim          !< Number of dimensions which has to be read to calculate 
                                                            !! Kernel on parameter model_param

        real(kind=dp), public              :: dt
        integer,       public              :: ndumps, decimate_factor
        integer,       public              :: nseis 
        integer,       public              :: npol
        real(kind=dp), public              :: windowlength
        real(kind=dp), public              :: timeshift_fwd, timeshift_bwd
        real(kind=dp), public, allocatable :: veloseis(:,:), dispseis(:,:)
        real(kind=dp), public, allocatable :: stf_fwd(:), stf_bwd(:)
        integer                            :: buffer_size
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
            procedure, pass                :: load_fw_points_rdbm
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
subroutine set_params(this, fwd_dir, bwd_dir, buffer_size, model_param)
    class(semdata_type)            :: this
    character(len=512), intent(in) :: fwd_dir, bwd_dir
    integer,            intent(in) :: buffer_size
    character(len=4),   intent(in), optional :: model_param
    character(len=512)             :: dirnam
    integer                        :: isim
    logical                        :: moment=.false., force=.false., single=.false.

    this%buffer_size = buffer_size

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
       write(lu_out,*) 'ERROR: Forward rundir does not seem to be an axisem rundirectory'
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
       write(lu_out,*) 'Backword simulation was ''moment'' source'
       write(lu_out,*) 'This is not implemented yet!'
       call pabort()
    elseif (force) then
       this%nsim_bwd = 2
       write(lu_out,*) 'Backword simulation was ''forces'' source'
       write(lu_out,*) 'This is not implemented yet!'
    elseif (single) then
       this%nsim_bwd = 1
       write(lu_out,*) 'Backword simulation was ''single'' source'
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
    
    if (present(model_param)) then
       this%model_param = model_param
    else
       this%model_param = 'vp'
    end if

    select case(trim(this%model_param))
    case('vp')
       this%ndim = 1
    case('vs')
       this%ndim = 6
    case default
        print *, 'ERROR in set_params(): unknown model param '//this%model_param
        call pabort
    end select
    write(lu_out, *) 'Model parameter: ', trim(this%model_param), &
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

        call nc_read_att_char(this%fwd(isim)%dump_type, &
                              'dump type (displ_only, displ_velo, fullfields)', &
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
        
        allocate( this%fwd(isim)%stf( this%fwd(isim)%ndumps ) )
        call check(nf90_get_var(  ncid   = this%fwd(isim)%surf,   &
                                  varid  = this%fwd(isim)%stf_varid, &
                                  values = this%fwd(isim)%stf  ))
        
    end do
        
    call flush(lu_out)

    do isim = 1, this%nsim_bwd
        ! Backward wavefield
        filename=trim(this%bwd(isim)%meshdir)//'/Data/ordered_output.nc4'
        
        if (verbose>0) write(lu_out,format20) trim(filename), myrank
        call nc_open_for_read(filename = filename,              &
                              ncid     = this%bwd(isim)%ncid) 

        call nc_read_att_char(this%bwd(isim)%dump_type, &
                              'dump type (displ_only, displ_velo, fullfields)', &
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
        
        call getvarid(            ncid     = this%bwd(isim)%surf,   &
                                  name     = "stf_dump",            &
                                  varid    = this%bwd(isim)%stf_varid)
        
        allocate( this%bwd(isim)%stf( this%bwd(isim)%ndumps ) )
        call check(nf90_get_var(  ncid   = this%bwd(isim)%surf,   &
                                  varid  = this%bwd(isim)%stf_varid, &
                                  values = this%bwd(isim)%stf  ))
        
    end do


    call flush(lu_out)
    call this%check_consistency()

    call flush(6) 
    this%files_open = .true.

    do isim = 1, this%nsim_fwd
        allocate(this%fwd(isim)%buffer(this%ndim))
        do istrainvar = 1, this%ndim
           status = this%fwd(isim)%buffer(istrainvar)%init(this%buffer_size,     &
                                                           this%fwd(isim)%ndumps)
        end do
    end do

    do isim = 1, this%nsim_bwd
        allocate(this%bwd(isim)%buffer(this%ndim))
        do istrainvar = 1, this%ndim
           status = this%bwd(isim)%buffer(istrainvar)%init(this%buffer_size,     &
                                                           this%bwd(isim)%ndumps)
        end do
    end do
    call flush(lu_out)

end subroutine open_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine close_files(this)
    class(semdata_type)   :: this
    integer              :: status, isim, istrainvar

    do isim = 1, this%nsim_fwd
       status = nf90_close(this%fwd(isim)%ncid)
       if (verbose>0) then
          write(lu_out,'(A,I1,A,F9.6)') ' Buffer efficiency fwd(', isim, '): ',  &
                                   this%fwd(isim)%buffer(1)%efficiency()
       end if
       do istrainvar = 1, this%ndim
           status = this%fwd(isim)%buffer(istrainvar)%freeme()
       end do
    end do
    deallocate(this%fwd)

    do isim = 1, this%nsim_bwd
       status = nf90_close(this%bwd(isim)%ncid)
       if (verbose>0) then
       write(lu_out,'(A,F9.6)') ' Buffer efficiency bwd   : ', & 
                           this%bwd(isim)%buffer(1)%efficiency()
       end if
       do istrainvar = 1, this%ndim
           status = this%bwd(isim)%buffer(istrainvar)%freeme()
       end do
    end do
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
    real(kind=dp)          :: stf_agreed_fwd(this%fwd(1)%ndumps)
    real(kind=dp)          :: stf_agreed_bwd(this%bwd(1)%ndumps)

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
    npol_agreed = this%fwd(1)%npol

    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",I7,") vs ", I7, " in the others")' 
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
    end do

    this%timeshift_fwd = real(source_shift_agreed_fwd, kind=dp)
    allocate(this%stf_fwd(ndumps_agreed))
    this%stf_fwd = real(stf_agreed_fwd, kind=dp)
    
    if (this%nsim_bwd > 0) then
       source_shift_agreed_bwd = this%bwd(1)%source_shift_t
       stf_agreed_bwd = this%bwd(1)%stf
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
       end do
       this%timeshift_bwd = real(source_shift_agreed_bwd, kind=dp)
       allocate(this%stf_bwd(ndumps_agreed))
       this%stf_bwd = real(stf_agreed_bwd, kind=dp)
    endif


    this%dt = dt_agreed
    this%decimate_factor = nseis_agreed / ndumps_agreed
    this%nseis  = ndumps_agreed * this%decimate_factor        

    call flush(lu_out)

end subroutine check_consistency
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_fw_points(this, coordinates, source_params)
    use simple_routines, only          : mult2d_1d

    class(semdata_type)               :: this
    real(kind=dp), intent(in)         :: coordinates(:,:)
    type(src_param_type), intent(in)  :: source_params
    type(kdtree2_result), allocatable :: nextpoint(:)
    real(kind=dp)                     :: load_fw_points(this%ndumps, this%ndim, &
                                                        size(coordinates,2))

    integer                           :: npoints
    integer                           :: pointid(size(coordinates,2))
    integer                           :: ipoint, isim, iclockold
    real(kind=dp)                     :: rotmesh_s(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_phi(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_z(size(coordinates,2))
    real(kind=dp)                     :: utemp(this%ndumps, this%ndim)
    
    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_fw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       call pabort 
    end if
    npoints = size(coordinates,2)

    
    ! Rotate points to FWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,   &
                          coordinates*1d3,                              &
                          source_params%lon, source_params%colat)

    allocate(nextpoint(1))
    do ipoint = 1, npoints
        call kdtree2_n_nearest( this%fwdtree,                           &
                                real([rotmesh_s(ipoint), rotmesh_z(ipoint)]), &
                                nn = 1,                                 &
                                results = nextpoint )
        
        pointid(ipoint) = nextpoint(1)%idx
    end do

    load_fw_points(:,:,:) = 0.0
    
    do ipoint = 1, npoints
    
       do isim = 1, this%nsim_fwd
            utemp = load_strain_point(this%fwd(isim),      &
                                      pointid(ipoint),     &
                                      this%model_param)

            iclockold = tick()
            select case(trim(this%model_param))
            !if (this%model_param.eq.'vs') then
            case('vp')
                load_fw_points(:, :, ipoint) = load_fw_points(:,:,ipoint)                   &
                                             + utemp * azim_factor(rotmesh_phi(ipoint),     &
                                                                   source_params%mij, isim, 1) 
            case('vs')
                !@TODO: utemp needs to be summed with azimfactors first before
                !       beeing rotated. I'd suggest summation first, because
                !       rotation is not for free.
                load_fw_points(:, :, ipoint) = load_fw_points(:, :, ipoint)                 &
                                             + rotate_straintensor(utemp,                   &
                                                                   rotmesh_phi(ipoint),     &
                                                                   source_params%mij, isim) 
            end select
            iclockold = tick(id=id_rotate, since=iclockold)
        end do !isim


        !read(*,*)
        !write(1002,'(7(ES15.6))') rotmesh_s(ipoint), rotmesh_z(ipoint), rotmesh_phi(ipoint), &
        !              azim_factor(rotmesh_phi(ipoint), source_params%mij, 1), &
        !              azim_factor(rotmesh_phi(ipoint), source_params%mij, 2), &
        !              azim_factor(rotmesh_phi(ipoint), source_params%mij, 3), &
        !              azim_factor(rotmesh_phi(ipoint), source_params%mij, 4)
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
    class(semdata_type)                :: this
    real(kind=dp), intent(in)          :: coordinates(:,:)
    type(rec_param_type)               :: receiver
    real(kind=dp)                      :: load_bw_points(this%ndumps, this%ndim, &
                                                         size(coordinates,2))
    real(kind=sp)                      :: utemp(this%ndumps, this%ndim)
    type(kdtree2_result), allocatable  :: nextpoint(:)

    integer                            :: pointid(size(coordinates,2))
    integer                            :: npoints
    integer                            :: ipoint
    real(kind=dp)                      :: rotmesh_s(size(coordinates,2))
    real(kind=dp)                      :: rotmesh_phi(size(coordinates,2))
    real(kind=dp)                      :: rotmesh_z(size(coordinates,2))

    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_bw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       call pabort 
    end if
    npoints = size(coordinates,2)

    ! Rotate points to BWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,   &
                          coordinates*1e3,                              &
                          receiver%lon, receiver%colat)

    allocate(nextpoint(1))
    do ipoint = 1, npoints
        call kdtree2_n_nearest( this%bwdtree,                           &
                                real([rotmesh_s(ipoint), rotmesh_z(ipoint)]), &
                                nn = 1,                                 &
                                results = nextpoint )
    
        pointid(ipoint) = nextpoint(1)%idx
    end do
    
    load_bw_points(:,:,:) = 0.0
    
    do ipoint = 1, npoints
        
        select case(receiver%component)
        case('Z')
            utemp = load_strain_point(this%bwd(1), pointid(ipoint), this%model_param)
            load_bw_points(:,:,ipoint) &
                =                               utemp / this%bwd(1)%amplitude
        case('R')
            utemp = load_strain_point(this%bwd(2), pointid(ipoint), this%model_param)
            load_bw_points(:,:,ipoint) &
                =   dcos(rotmesh_phi(ipoint)) * utemp / this%bwd(2)%amplitude
        case('T')
            utemp = load_strain_point(this%bwd(2), pointid(ipoint), this%model_param)
            load_bw_points(:,:,ipoint) &
                = - dsin(rotmesh_phi(ipoint)) * utemp / this%bwd(2)%amplitude 
        end select

        if (this%model_param.eq.'vs') then
            load_bw_points(:,:,ipoint) &
                = rotate_straintensor(load_bw_points(:,:,ipoint), &
                                      rotmesh_phi(ipoint),        &
                                      real([1, 1, 1, 0, 0, 0], kind=dp), 1)
        end if
    end do !ipoint

end function load_bw_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_fw_points_rdbm(this, source_params, reci_source_params, component)
    use simple_routines, only          : mult2d_1d
    use finite_elem_mapping, only      : inside_element

    class(semdata_type)               :: this
    type(src_param_type), intent(in)  :: source_params(:)
    type(src_param_type), intent(in)  :: reci_source_params
    character(len=1), intent(in)      :: component
    !real(kind=dp)                     :: load_fw_points_rdbm(this%ndumps, this%ndim, &
    !                                                         size(source_params))
    real(kind=dp)                     :: load_fw_points_rdbm(this%ndumps, 1, &
                                                             size(source_params))

    type(kdtree2_result), allocatable :: nextpoint(:)
    integer                           :: npoints, nnext_points
    integer                           :: pointid(size(source_params))
    integer                           :: ipoint, inext_point, isim, iclockold, i, icp
    integer                           :: corner_point_ids(4), eltype(1)
    integer, allocatable              :: gll_point_ids(:,:)
    real(kind=dp)                     :: corner_points(4,2)
    real(kind=dp)                     :: cps(1), cpz(1)
    real(kind=dp)                     :: rotmesh_s(size(source_params))
    real(kind=dp)                     :: rotmesh_phi(size(source_params))
    real(kind=dp)                     :: rotmesh_z(size(source_params))
    real(kind=dp)                     :: utemp(this%ndumps, this%ndim)
    real(kind=dp)                     :: coordinates(size(source_params),3)
    real(kind=dp)                     :: mij_buff(6)
    real(kind=dp)                     :: xi, eta

    if (trim(this%dump_type) == 'displ_only') then
        nnext_points = 6 ! 6, because this is the maximum valence in the mesh
        allocate(gll_point_ids(0:this%npol, 0:this%npol))
    else
        nnext_points = 1
    endif
    
    npoints = size(source_params)

    do ipoint = 1, npoints
        coordinates(ipoint,1) = source_params(ipoint)%x
        coordinates(ipoint,2) = source_params(ipoint)%y
        coordinates(ipoint,3) = source_params(ipoint)%z
    enddo
    
    ! Rotate points to FWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z, coordinates * 1d3, &
                          reci_source_params%lon, reci_source_params%colat)


    allocate(nextpoint(nnext_points))
    do ipoint = 1, npoints
        call kdtree2_n_nearest( this%fwdtree, &
                                real([rotmesh_s(ipoint), rotmesh_z(ipoint)]), &
                                nn = nnext_points, &
                                results = nextpoint )
        
        pointid(ipoint) = nextpoint(1)%idx

        if (verbose > 0) &
            write(6,*) 'nearest point', nextpoint(1)%dis**.5, nextpoint(1)%idx, &
                       (rotmesh_s(ipoint)**2  + rotmesh_z(ipoint)**2)**.5, &
                       (this%fwdmesh%s(pointid(ipoint))**2  &
                            + this%fwdmesh%z(pointid(ipoint))**2)**.5 

        
        if (trim(this%dump_type) == 'displ_only') then
            do inext_point = 1, nnext_points
                ! get cornerpoints of finite element
                call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
                                        varid  = this%fwd(1)%fem_mesh_varid, &
                                        start  = [1, nextpoint(inext_point)%idx], &
                                        count  = [4, 1], &
                                        values = corner_point_ids))
                
                call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
                                        varid  = this%fwd(1)%eltype_varid, &
                                        start  = [nextpoint(inext_point)%idx], &
                                        count  = [1], &
                                        values = eltype))
                
                do icp = 1, 4
                    call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
                                            varid  = this%fwd(1)%mesh_s_varid, &
                                            start  = [corner_point_ids(icp) + 1], &
                                            count  = [1], &
                                            values = cps))
                    
                    call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
                                            varid  = this%fwd(1)%mesh_z_varid, &
                                            start  = [corner_point_ids(icp) + 1], &
                                            count  = [1], &
                                            values = cpz))

                    corner_points(icp, 1) = cps(1)
                    corner_points(icp, 2) = cpz(1)
                enddo                        
                ! test point to be inside, if so, exit
                if (inside_element(rotmesh_s(ipoint), rotmesh_z(ipoint), &
                                   corner_points, eltype(1), xi=xi, eta=eta)) then
                    write(6,*) 'eltype     = ', eltype
                    write(6,*) 'xi, eta    = ', xi, eta
                    write(6,*) 'element id = ', nextpoint(inext_point)%idx
                    exit
                endif
            enddo

            if (inext_point >= nnext_points) then
               write(6,*) 'ERROR: element not found. '
               write(6,*) '       Try increasing nnext_points in case this problem persists'
               call pabort
            endif

            ! get gll points of spectral element
            gll_point_ids = -1
            write(6,*) 'npol = ', this%npol
            write(6,*) 'element id = ', nextpoint(inext_point)%idx
            call check(nf90_get_var(ncid   = this%fwd(1)%mesh,   &
                                    varid  = this%fwd(1)%sem_mesh_varid, &
                                    start  = [1, 1, nextpoint(inext_point)%idx], &
                                    count  = [this%npol+1, this%npol+1, 1], &
                                    values = gll_point_ids))
            write(6,*) 'gll_point_ids = ', gll_point_ids(:,0)
        endif
    end do

    load_fw_points_rdbm(:,:,:) = 0.0
    
    do ipoint = 1, npoints
    
       if (this%model_param == 'vp') then
          select case(component)
          case('Z')
               isim = 1
               utemp = load_strain_point(this%fwd(isim),      &
                                         pointid(ipoint),     &
                                         this%model_param)

               load_fw_points_rdbm(:, :, ipoint) = utemp

          case('R')
               isim = 2
               utemp = load_strain_point(this%fwd(isim),      &
                                         pointid(ipoint),     &
                                         this%model_param)

               load_fw_points_rdbm(:, :, ipoint) = utemp 

          case('T')
               load_fw_points_rdbm(:, :, ipoint) = 0

          case('N')
               isim = 2
               utemp = load_strain_point(this%fwd(isim),      &
                                         pointid(ipoint),     &
                                         this%model_param)

               load_fw_points_rdbm(:, :, ipoint) = &
                       - utemp * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 1d0, 0d0/), isim, 1) 

          case('E')
               isim = 2
               utemp = load_strain_point(this%fwd(isim),      &
                                         pointid(ipoint),     &
                                         this%model_param)

               load_fw_points_rdbm(:, :, ipoint) = &
                       utemp * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 0d0, 1d0/), isim, 1) 

          case default
               write(6,*) 'component "', component, '" unknown or not yet implemented'
               call pabort
          end select
       elseif (this%model_param == 'vs') then
          select case(component)
          case('Z')
               isim = 1
               utemp = load_strain_point(this%fwd(isim),      &
                                         pointid(ipoint),     &
                                         this%model_param)

               ! rotate source mt to global cartesian system
               mij_buff = rotate_symm_tensor_voigt_xyz_src_to_xyz_earth( &
                                source_params(ipoint)%mij_voigt, &
                                source_params(ipoint)%lon, &
                                source_params(ipoint)%colat)

               ! rotate source mt to receiver cartesian system
               mij_buff = rotate_symm_tensor_voigt_xyz_earth_to_xyz_src( &
                                mij_buff, &
                                reci_source_params%lon, &
                                reci_source_params%colat)

               ! rotate source mt to receiver s,phi,z system
               mij_buff = rotate_symm_tensor_voigt_xyz_to_src(mij_buff, rotmesh_phi(ipoint))

               mij_buff = mij_buff / this%fwd(isim)%amplitude

               load_fw_points_rdbm(:, :, ipoint) = 0
               
               do i = 1, 3
                  load_fw_points_rdbm(:, 1, ipoint) = &
                        load_fw_points_rdbm(:, 1, ipoint) &
                            + mij_buff(i) * utemp(:,i)
               enddo 

               ! components 4-6 need a factor of two because of voigt mapping
               ! without factor of two in the strain
               i = 5
               load_fw_points_rdbm(:, 1, ipoint) = &
                     load_fw_points_rdbm(:, 1, ipoint) &
                         + 2 * mij_buff(i) * utemp(:,i)

          case('N')
               isim = 2
               utemp = load_strain_point(this%fwd(isim),      &
                                         pointid(ipoint),     &
                                         this%model_param)

               ! rotate source mt to global cartesian system
               mij_buff = rotate_symm_tensor_voigt_xyz_src_to_xyz_earth( &
                                source_params(ipoint)%mij_voigt, &
                                source_params(ipoint)%lon, &
                                source_params(ipoint)%colat)

               ! rotate source mt to receiver cartesian system
               mij_buff = rotate_symm_tensor_voigt_xyz_earth_to_xyz_src( &
                                mij_buff, &
                                reci_source_params%lon, &
                                reci_source_params%colat)

               ! rotate source mt to receiver s,phi,z system
               mij_buff = rotate_symm_tensor_voigt_xyz_to_src(mij_buff, rotmesh_phi(ipoint))

               mij_buff = mij_buff / this%fwd(isim)%amplitude

               load_fw_points_rdbm(:, :, ipoint) = 0

               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(1) * utemp(:,1) &
                        * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 1d0, 0d0/), isim, 1) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(2) * utemp(:,2) &
                        * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 1d0, 0d0/), isim, 1) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(3) * utemp(:,3) &
                        * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 1d0, 0d0/), isim, 1) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(4) * utemp(:,4) &
                        * 2 * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 1d0, 0d0/), isim, 2) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(5) * utemp(:,5) &
                        * 2 * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 1d0, 0d0/), isim, 1) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(6) * utemp(:,6) &
                        * 2 * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 1d0, 0d0/), isim, 2) 
               
               !@TODO not sure why we need the - sign here. Might be because N
               !      is in negative theta direction
               load_fw_points_rdbm(:, 1, ipoint) = - load_fw_points_rdbm(:, 1, ipoint)

          case('E')
               isim = 2
               utemp = load_strain_point(this%fwd(isim),      &
                                         pointid(ipoint),     &
                                         this%model_param)

               ! rotate source mt to global cartesian system
               mij_buff = rotate_symm_tensor_voigt_xyz_src_to_xyz_earth( &
                                source_params(ipoint)%mij_voigt, &
                                source_params(ipoint)%lon, &
                                source_params(ipoint)%colat)

               ! rotate source mt to receiver cartesian system
               mij_buff = rotate_symm_tensor_voigt_xyz_earth_to_xyz_src( &
                                mij_buff, &
                                reci_source_params%lon, &
                                reci_source_params%colat)

               ! rotate source mt to receiver s,phi,z system
               mij_buff = rotate_symm_tensor_voigt_xyz_to_src(mij_buff, rotmesh_phi(ipoint))

               mij_buff = mij_buff / this%fwd(isim)%amplitude

               load_fw_points_rdbm(:, :, ipoint) = 0
               
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(1) * utemp(:,1) &
                        * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 0d0, 1d0/), isim, 1) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(2) * utemp(:,2) &
                        * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 0d0, 1d0/), isim, 1) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(3) * utemp(:,3) &
                        * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 0d0, 1d0/), isim, 1) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(4) * utemp(:,4) &
                        * 2 * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 0d0, 1d0/), isim, 2) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(5) * utemp(:,5) &
                        * 2 * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 0d0, 1d0/), isim, 1) 
               load_fw_points_rdbm(:, 1, ipoint) &
                    = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(6) * utemp(:,6) &
                        * 2 * azim_factor_bw(rotmesh_phi(ipoint), (/0d0, 0d0, 1d0/), isim, 2) 

          case default
               write(6,*) 'component "', component, '" unknown or not yet implemented'
               call pabort
          end select
         
       else
          call pabort
       endif

    end do !ipoint

end function load_fw_points_rdbm
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_strain_point(sem_obj, pointid, model_param)

    type(ncparamtype), intent(in)   :: sem_obj
    integer, intent(in)             :: pointid
    character(len=*), intent(in)    :: model_param
    real(kind=dp), allocatable      :: load_strain_point(:,:)
    real(kind=dp), allocatable      :: strain_buff(:,:)

    integer                         :: start_chunk, iread, gll_to_read
    integer                         :: iclockold, status, istrainvar
    real(kind=sp)                   :: utemp(sem_obj%ndumps)
    real(kind=sp)                   :: utemp_chunk(sem_obj%chunk_gll, sem_obj%ndumps)

    if (trim(sem_obj%dump_type) /= 'fullfields') then
        write(6,*) 'ERROR: trying to read strain from a file that was not'
        write(6,*) '       written with dump_type "fullfields"'
        call pabort
    endif

    select case(model_param)
    case('vp')
        allocate(load_strain_point(sem_obj%ndumps, 1))

        !iclockold = tick()
        status = sem_obj%buffer(1)%get(pointid, utemp)
        !iclockold = tick(id=id_buffer, since=iclockold)

        if (status.ne.0) then
           start_chunk = ((pointid-1) / sem_obj%chunk_gll) * sem_obj%chunk_gll + 1

           ! Only read to last point, not further
           gll_to_read = min(sem_obj%chunk_gll, sem_obj%ngll + 1 - start_chunk)
           !iclockold = tick()
           call nc_getvar( ncid   = sem_obj%snap,           & 
                           varid  = sem_obj%strainvarid(6), &
                           start  = [start_chunk, 1],       &
                           count  = [gll_to_read, sem_obj%ndumps], &
                           values = utemp_chunk(1:gll_to_read, :)) 

           !iclockold = tick(id=id_netcdf, since=iclockold)

           do iread = 0, sem_obj%chunk_gll - 1
               status = sem_obj%buffer(1)%put(start_chunk + iread, utemp_chunk(iread+1,:))
           end do
           !iclockold = tick(id=id_buffer, since=iclockold)
           load_strain_point(:,1) = real(utemp_chunk(pointid-start_chunk+1, :), kind=dp)
        else
           load_strain_point(:,1) = real(utemp, kind=dp)
        end if

    case('vs')
        allocate(strain_buff(sem_obj%ndumps, 6))

        do istrainvar = 1, 6

            if (sem_obj%strainvarid(istrainvar).eq.-1) then
                strain_buff(:, istrainvar) = 0
                cycle ! For monopole source which does not have this component.
            endif

            !iclockold = tick()
            status = sem_obj%buffer(istrainvar)%get(pointid, utemp)
            !iclockold = tick(id=id_buffer, since=iclockold)
            
            if (status.ne.0) then
               start_chunk = ((pointid-1) / sem_obj%chunk_gll) * sem_obj%chunk_gll + 1
               
               ! Only read to last point, not further
               gll_to_read = min(sem_obj%chunk_gll, sem_obj%ngll + 1 - start_chunk)

               !iclockold = tick()
               call nc_getvar( ncid   = sem_obj%snap,           & 
                               varid  = sem_obj%strainvarid(istrainvar), &
                               start  = [start_chunk, 1],       &
                               count  = [gll_to_read, sem_obj%ndumps], &
                               values = utemp_chunk(1:gll_to_read, :)) 

               !iclockold = tick(id=id_netcdf, since=iclockold)
               do iread = 0, gll_to_read - 1
                   status = sem_obj%buffer(istrainvar)%put(start_chunk + iread, &
                                                           utemp_chunk(iread+1,:))
               end do
               !iclockold = tick(id=id_buffer, since=iclockold)
               strain_buff(:,istrainvar) &
                    = real(utemp_chunk(pointid-start_chunk+1, :), kind=dp)
            else
               strain_buff(:,istrainvar) = real(utemp, kind=dp)
            end if

        end do

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

end function load_strain_point
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine build_kdtree(this)
    class(semdata_type)        :: this
    real(kind=sp), allocatable :: mesh(:,:)
    integer                    :: npoints

    if (.not.this%meshes_read) then
        print *, 'ERROR in build_kdtree(): Meshes have not been read yet'
        print *, 'Call read_meshes() before build_kdtree!'
        call pabort
    end if

    write(lu_out,*) ' Reshaping mesh variables'
    call flush(lu_out)

    if (trim(this%dump_type) == 'displ_only') then
        npoints = this%fwdmesh%nelem ! midpoints only
    else
        npoints = this%fwdmesh%npoints
    endif

    allocate(mesh(2, npoints))
    mesh = transpose(reshape([this%fwdmesh%s, this%fwdmesh%z], [npoints, 2]))

    print *, npoints

    write(lu_out,*) ' Building forward KD-Tree'
    call flush(lu_out)
    ! KDtree in forward field
    this%fwdtree => kdtree2_create(mesh,              &
                                   dim = 2,           &
                                   sort = .true.,     &
                                   rearrange = .true.)
    
    deallocate(mesh)                           

    if (trim(this%dump_type) == 'displ_only') then
        npoints = this%bwdmesh%nelem ! midpoints only
    else
        npoints = this%bwdmesh%npoints
    endif

    ! KDtree in backward field
    if (this%nsim_bwd > 0) then
        write(lu_out,*) ' Building backward KD-Tree'
        call flush(lu_out)
        allocate(mesh(2, npoints))
        mesh = transpose(reshape([this%bwdmesh%s, this%bwdmesh%z], [npoints, 2]))
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
        
        allocate(this%fwdmesh%s(this%fwdmesh%nelem))
        allocate(this%fwdmesh%z(this%fwdmesh%nelem))

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
                      name  = "mesh_S",              &
                      varid = this%fwd(1)%mesh_s_varid)

        call getvarid(ncid  = this%fwd(1)%mesh,   &
                      name  = "mesh_Z",              &
                      varid = this%fwd(1)%mesh_z_varid)
           
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
                       values = this%fwdmesh%s ) 
        call nc_getvar(ncid   = this%fwd(1)%mesh,       &
                       varid  = ncvarid_mesh_z,         &
                       start  = 1,                    &
                       count  = this%fwdmesh%nelem, &
                       values = this%fwdmesh%z ) 
    
    else
        allocate(this%fwdmesh%s(this%fwdmesh%npoints))
        allocate(this%fwdmesh%z(this%fwdmesh%npoints))
           
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
            
            allocate(this%bwdmesh%s(this%fwdmesh%nelem))
            allocate(this%bwdmesh%z(this%fwdmesh%nelem))

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
                          name  = "mesh_S",              &
                          varid = this%bwd(1)%mesh_s_varid)

            call getvarid(ncid  = this%bwd(1)%mesh,   &
                          name  = "mesh_Z",              &
                          varid = this%bwd(1)%mesh_z_varid)
           
            
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
                           values = this%bwdmesh%s ) 
            call nc_getvar(ncid   = this%bwd(1)%mesh, &
                           varid  = ncvarid_mesh_z,   &
                           start  = 1,                &
                           count  = this%bwdmesh%nelem, &
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
                             

    ! Mesh sanity checks
    mesherror = .false.
    if (maxval(abs(this%fwdmesh%s)) > 1e32) then
       write(*,*) 'Maximum value of S in the forward mesh is unreasonably large'
       write(*,*) 'maxval(S): ', this%fwdmesh%s(maxloc(abs(this%fwdmesh%s))), ' m'
       mesherror = .true.
    end if
    if (maxval(abs(this%fwdmesh%z)) > 1e32) then
       write(*,*) 'Maximum value of Z in the forward mesh is unreasonably large'
       write(*,*) 'maxval(Z): ', this%fwdmesh%z(maxloc(abs(this%fwdmesh%z))), ' m'
       mesherror = .true.
    end if
    if (maxval(this%fwdmesh%theta).gt.180.0) then
       write(*,*) 'Maximum value of theta in the backward mesh is larger than 180°'
       write(*,*) 'maxval(theta): ', this%fwdmesh%theta(maxloc(abs(this%fwdmesh%theta)))
       write(*,*) 'maxloc(theta): ', maxloc(abs(this%fwdmesh%theta))
       mesherror = .true.
    end if

    if (this%nsim_bwd > 0) then
        if (maxval(abs(this%bwdmesh%s)) > 1e32) then
           write(*,*) 'Maximum value of S in the backward mesh is unreasonably large'
           write(*,*) 'maxval(S): ', this%bwdmesh%s(maxloc(abs(this%bwdmesh%s))), ' m'
           mesherror = .true.
        end if
        if (maxval(abs(this%bwdmesh%z)) > 1e32) then
           write(*,*) 'Maximum value of Z in the backward mesh is unreasonably large'
           write(*,*) 'maxval(Z): ', this%bwdmesh%z(maxloc(abs(this%bwdmesh%z))), ' m'
           mesherror = .true.
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
