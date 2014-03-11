module readfields
    use global_parameters, only            : sp, dp, pi, deg2rad, verbose, lu_out, myrank
    !use type_parameter, only               : src_param_type, rec_param_type, parameter_type
    use source_class,          only        : src_param_type
    use receiver_class,        only        : rec_param_type
    use buffer, only                       : buffer_type
    use netcdf
    use kdtree2_module                     

    implicit none

    type meshtype
        real(kind=sp), allocatable        :: s(:), z(:)
        integer                           :: npoints
        real(kind=dp), allocatable        :: theta(:)
        integer                           :: nsurfelem
    end type

    type ncparamtype
        integer                           :: ncid
        integer                           :: snap, surf, mesh, seis  ! Group IDs
        integer                           :: straintrace, dev_strain ! Variable IDs
        integer                           :: seis_disp, seis_velo    ! Variable IDs
        character(len=200)                :: meshdir
        integer                           :: ndumps, nseis
        integer                           :: source_shift_samples    
        real(kind=dp)                     :: source_shift_t
        !logical                           :: ordered_output
        type(buffer_type)                 :: buffer
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
        real(kind=dp), public              :: dt
        integer,       public              :: ndumps, decimate_factor
        integer,       public              :: nseis ! ndumps * decimate_factor
        real(kind=dp), public              :: windowlength
        real(kind=dp), public              :: timeshift_fwd, timeshift_bwd
        real(kind=dp), public, allocatable :: veloseis(:,:), dispseis(:,:)
         
        real(kind=dp), dimension(3,3)      :: rot_mat, trans_rot_mat

        contains 
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

!-------------------------------------------------------------------------------
subroutine set_params(this, fwd_dir, bwd_dir)
    class(semdata_type)            :: this
    character(len=512), intent(in) :: fwd_dir, bwd_dir
    character(len=512)             :: dirnam
    integer                        :: isim
    logical                        :: moment=.false., force=.false.

    dirnam = trim(fwd_dir)//'/MZZ/simulation.info'
    write(lu_out,*) 'Inquiring: ', trim(dirnam)
    inquire( file = trim(dirnam), exist = moment)

    dirnam = trim(fwd_dir)//'/PZ/simulation.info'
    write(lu_out,*) 'Inquiring: ', trim(dirnam)
    inquire( file = trim(dirnam), exist = force)
    if (moment) then
       this%nsim_fwd = 4
       write(lu_out,*) 'Forward simulation was ''moment'' source'
    elseif (force) then
       this%nsim_fwd = 2
       write(lu_out,*) 'Forward simulation was ''forces'' source'
    else 
       this%nsim_fwd = 1
       write(lu_out,*) 'Forward simulation was ''single'' source'
    end if

    !this%nsim_fwd = parameters%nsim_fwd
    this%nsim_bwd =  1 !parameters%nsim_bwd
    
    allocate( this%fwd(this%nsim_fwd) )
    allocate( this%bwd(this%nsim_bwd) )

    select case(this%nsim_fwd)
    case(1)
        this%fwd(1)%meshdir = fwd_dir//'/'

    case(2) 
        this%fwd(1)%meshdir = trim(fwd_dir)//'/PZ/'
        this%fwd(2)%meshdir = trim(fwd_dir)//'/PX/'

    case(4)
        this%fwd(1)%meshdir = trim(fwd_dir)//'/MZZ/'
        this%fwd(2)%meshdir = trim(fwd_dir)//'/MXX_P_MYY/'
        this%fwd(3)%meshdir = trim(fwd_dir)//'/MXZ_MYZ/'
        this%fwd(4)%meshdir = trim(fwd_dir)//'/MXY_MXX_M_MYY/'
    end select
    
    do isim = 1, this%nsim_bwd
        this%bwd(isim)%meshdir = trim(bwd_dir)//'/'
    end do

    this%params_set = .true.

end subroutine
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine open_files(this)

    class(semdata_type)              :: this
    integer                          :: status, isim, chunks(2), deflev
    character(len=200)               :: format20, format21, filename
    real(kind=sp)                    :: temp

    if (.not.this%params_set) then
        print *, 'ERROR in open_files(): Parameters have to be set first'
        print *, 'Call set_params before open_files()'
        stop
    end if

    format20 = "('  Trying to open NetCDF file ', A, ' on CPU ', I5)"
    format21 = "('  Succeded,  has NCID ', I6, ', Snapshots group NCID: ', I6)"

    do isim = 1, this%nsim_fwd
        ! Forward wavefield
        filename=trim(this%fwd(isim)%meshdir)//'/ordered_output.nc4'
        
        !inquire(file=filename, exist=this%fwd(isim)%ordered_output)
        !if (.not.this%fwd(isim)%ordered_output) then
        !    filename = trim(this%fwd(isim)%meshdir)//'/Data/axisem_output.nc4'
        !end if
        
        if (verbose>0) write(6,format20) trim(filename), myrank
        call nc_open_for_read(    filename = filename,              &
                                  ncid     = this%fwd(isim)%ncid) 

        call getgrpid(  ncid     = this%fwd(isim)%ncid,   &
                        name     = "Snapshots",           &
                        grp_ncid = this%fwd(isim)%snap)

        
        call getvarid(  ncid     = this%fwd(isim)%snap,   &
                        name     = "straintrace",         &
                        varid    = this%fwd(isim)%straintrace) 
        call check(nf90_inquire_variable(ncid       = this%fwd(isim)%snap,   &
                                         varid      = this%fwd(isim)%straintrace, &
                                         chunksizes = chunks, &
                                         deflate_level = deflev) )

        write(lu_out, "('FWD SIM:', I2, ', Chunksizes:', 2(I7), ', deflate level: ', I2)") &
              isim, chunks, deflev
        if (verbose>0) write(lu_out,format21) this%fwd(isim)%ncid, this%fwd(isim)%snap 
        
        !if (this%fwd(isim)%ordered_output) then
        !    filename = trim(this%fwd(isim)%meshdir)//'/Data/axisem_output.nc4'
       
        !    if (verbose>0) write(6,format20) trim(filename), myrank
        !    call nc_open_for_read(filename = filename,              &
        !                          ncid     = this%fwd(isim)%ncid) 

        !end if
        
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
        ! Was a hack because dt was written out wrong in earlier AxiSEM versions
        !this%fwd(isim)%dt = 1.7064171433448792
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

        !call getvarid(            ncid     = this%fwd(isim)%surf,   &
        !                          name     = "stf_seis",         &
        !                          varid    = varid_stf)
        !
        !call check(nf90_get_var(  ncid   = this%fwd(isim)%seis, &
        !                          varid  = this%fwd(isim)%stf
    end do
        

    do isim = 1, this%nsim_bwd
        ! Backward wavefield
        filename=trim(this%bwd(isim)%meshdir)//'/ordered_output.nc4'
        
        !inquire(file=filename, exist=this%bwd(isim)%ordered_output)
        !if (.not.this%bwd(isim)%ordered_output) then
        !    filename = trim(this%bwd(isim)%meshdir)//'/Data/axisem_output.nc4'
        !end if

        if (verbose>0) write(6,format20) trim(filename), myrank
        call nc_open_for_read(filename = filename,              &
                              ncid     = this%bwd(isim)%ncid) 

        call getgrpid( ncid     = this%bwd(isim)%ncid,   &
                       name     = "Snapshots",           &
                       grp_ncid = this%bwd(isim)%snap)

        call getvarid(ncid     = this%bwd(isim)%snap,   &
                      name     = "straintrace",         &
                      varid    = this%bwd(isim)%straintrace) 
        if (verbose>0) write(6,format21) this%bwd(isim)%ncid, this%bwd(isim)%snap 
        
        !if (this%bwd(isim)%ordered_output) then
        !    filename = trim(this%bwd(isim)%meshdir)//'/Data/axisem_output.nc4'
       
        !    if (verbose>0) write(6,format20) trim(filename), myrank
        !    call nc_open_for_read(filename = filename,              &
        !                          ncid     = this%bwd(isim)%ncid) 

        !end if

        call getgrpid( ncid     = this%bwd(isim)%ncid,   &
                       name     = "Surface",             &
                       grp_ncid = this%bwd(isim)%surf)
        call getgrpid( ncid     = this%bwd(isim)%ncid,   &
                       name     = "Mesh",                &
                       grp_ncid = this%bwd(isim)%mesh)
        !call getgrpid( ncid     = this%bwd(isim)%ncid,   &
        !               name     = "Seismograms",         &
        !               grp_ncid = this%bwd(isim)%seis)
        !call check(nf90_inq_ncid( ncid     = this%bwd(isim)%ncid,   &
        !                          name     = "Surface",             &
        !                          grp_ncid = this%bwd(isim)%surf))
        !call check(nf90_inq_ncid( ncid     = this%bwd(isim)%ncid,   &
        !                          name     = "Mesh",                &
        !                          grp_ncid = this%bwd(isim)%mesh))
        !call check(nf90_inq_ncid( ncid     = this%bwd(isim)%ncid,   &
        !                          name     = "Seismograms",         &
        !                          grp_ncid = this%bwd(isim)%seis))

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
        
    end do


    call this%check_consistency()

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
    class(semdata_type)   :: this
    integer              :: status, isim

    do isim = 1, this%nsim_fwd
       status = nf90_close(this%fwd(isim)%ncid)
       if (verbose>0) then
          write(lu_out,'(A,I1,A,F9.6)') ' Buffer efficiency fwd(', isim, '): ',  &
                                   this%fwd(isim)%buffer%efficiency()
       end if
       status = this%fwd(isim)%buffer%freeme()
    end do
    deallocate(this%fwd)
    do isim = 1, this%nsim_bwd
       status = nf90_close(this%bwd(isim)%ncid)
       if (verbose>0) then
       write(lu_out,'(A,F9.6)') ' Buffer efficiency bwd   : ', & 
                           this%bwd(1)%buffer%efficiency()
       end if
       status = this%bwd(isim)%buffer%freeme()
    end do
    deallocate(this%bwd)

end subroutine close_files

!------------------------------------------------------------------------------
subroutine check_consistency(this)
    !< Checks consistency of the wavefield dumps
    class(semdata_type)    :: this
    integer                :: isim
    real(kind=dp)          :: dt_agreed
    character(len=512)     :: fmtstring
    integer                :: ndumps_agreed, nseis_agreed
    real(kind=dp)          :: source_shift_agreed_fwd, source_shift_agreed_bwd

    ! Check whether the sampling period is the same in all files
    dt_agreed = this%fwd(1)%dt
    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the others")' 
    do isim = 1, this%nsim_fwd
       if (dt_agreed.ne.this%fwd(isim)%dt) then
          write(*,fmtstring) 'dt', isim, dt_agreed, this%fwd(isim)%dt
          stop
       end if
    end do

    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the forward case")' 

    do isim = 1, this%nsim_bwd
       if (dt_agreed.ne.this%bwd(isim)%dt) then
          write(*,fmtstring) 'dt', isim, dt_agreed, this%bwd(isim)%dt
          stop
       end if
    end do

    this%dt = dt_agreed

    ndumps_agreed = this%fwd(1)%ndumps
    nseis_agreed  = this%fwd(1)%nseis


    ! Check whether the number of dumps (time samples) is the same in all files
    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",I7,") vs ", I7, " in the others")' 
    do isim = 1, this%nsim_fwd
       if (ndumps_agreed.ne.this%fwd(1)%ndumps) then
          write(*,fmtstring) 'ndumps', isim, ndumps_agreed, this%fwd(isim)%ndumps
          stop
       end if
       if (nseis_agreed.ne.this%fwd(1)%nseis) then
          write(*,fmtstring) 'nseis', isim, nseis_agreed, this%fwd(isim)%nseis
          stop
       end if
    end do


    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",I7,"s) vs ", I7, " in the forward case")' 

    do isim = 1, this%nsim_bwd
       if (ndumps_agreed.ne.this%bwd(1)%ndumps) then
          write(*,fmtstring) 'ndumps', isim, ndumps_agreed, this%bwd(isim)%ndumps
          stop
       end if
       if (nseis_agreed.ne.this%bwd(1)%nseis) then
          write(*,fmtstring) 'nseis', isim, nseis_agreed, this%bwd(isim)%nseis
          stop
       end if
    end do

    this%ndumps = ndumps_agreed
    this%windowlength = ndumps_agreed * dt_agreed

    ! Check whether the source time shift is the same in all files
    source_shift_agreed_fwd = this%fwd(1)%source_shift_t
    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the others")' 
    do isim = 1, this%nsim_fwd
       if (source_shift_agreed_fwd.ne.this%fwd(isim)%source_shift_t) then
          write(*,fmtstring) 'source time shift', isim, source_shift_agreed_fwd, &
                             this%fwd(isim)%source_shift_t
          stop
       end if
    end do

    source_shift_agreed_bwd = this%bwd(1)%source_shift_t
    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the others")' 

    do isim = 1, this%nsim_bwd
       if (source_shift_agreed_bwd.ne.this%bwd(isim)%source_shift_t) then
          write(*,fmtstring) 'source time shift', isim, source_shift_agreed_bwd, &
                             this%bwd(isim)%source_shift_t
          stop
       end if
    end do

    this%timeshift_fwd = real(source_shift_agreed_fwd, kind=dp)
    this%timeshift_bwd = real(source_shift_agreed_bwd, kind=dp)

    this%dt = dt_agreed
    this%decimate_factor = nseis_agreed / ndumps_agreed
    this%nseis  = ndumps_agreed * this%decimate_factor        

end subroutine

!-------------------------------------------------------------------------------
function load_fw_points(this, coordinates, source_params)
    use source_class, only             : src_param_type
    class(semdata_type)               :: this
    real(kind=dp), intent(in)         :: coordinates(:,:)
    type(src_param_type), intent(in)  :: source_params
    real(kind=dp)                     :: load_fw_points(this%fwd(1)%ndumps, &
                                                        size(coordinates,2))
    real(kind=sp)                     :: utemp(this%fwd(1)%ndumps)
    type(kdtree2_result), allocatable :: nextpoint(:)

    integer                           :: npoints
    integer                           :: pointid(size(coordinates,2))
    integer                           :: ipoint, isim, status
    real(kind=dp)                     :: rotmesh_s(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_phi(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_z(size(coordinates,2))
    
    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_fw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       stop 2
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
        !print *, 'Original coordinates: ', coordinates(:,ipoint)
        !print *, 'Coordinates:    ', rotmesh_s(ipoint), rotmesh_z(ipoint), ', next pointid: ', pointid(ipoint)
        !print *, 'CO of SEM point:', this%fwdmesh%s(pointid(ipoint)), this%fwdmesh%z(pointid(ipoint))

        write(1000, '(7(ES15.6), I8)') coordinates(:, ipoint), rotmesh_s(ipoint), rotmesh_z(ipoint), &
                                       this%fwdmesh%s(pointid(ipoint)), this%fwdmesh%z(pointid(ipoint)), pointid(ipoint)

    end do

    load_fw_points(:,:) = 0.0
    
    do ipoint = 1, npoints
        
        !print *, 'Azim factors: '
        do isim = 1, this%nsim_fwd
           !write(*,*) 'Reading point ', pointid(ipoint), ' (', ipoint, ') with ', &
           !           this%fwd(isim)%ndumps, ' values'
            status = this%fwd(isim)%buffer%get(pointid(ipoint), utemp)
            if (status.ne.0) then
               !write(*,*) 'Did not find point', ipoint, ' in buffer, rereading'
               !if (this%fwd(isim)%ordered_output) then
               !   call check( nf90_get_var( ncid   = this%fwd(isim)%snap,        & 
               !                             varid  = this%fwd(isim)%straintrace, &
               !                             start  = [1, pointid(ipoint)],       &
               !                             count  = [this%fwd(isim)%ndumps, 1], &
               !                             values = utemp) )
               !else
                  call check( nf90_get_var( ncid   = this%fwd(isim)%snap,        & 
                                            varid  = this%fwd(isim)%straintrace, &
                                            start  = [pointid(ipoint), 1],       &
                                            count  = [1, this%fwd(isim)%ndumps], &
                                            values = utemp) )
               !end if
               status = this%fwd(isim)%buffer%put(pointid(ipoint), utemp)
            !else
            !   write(*,*) 'Found point', ipoint, ' (',pointid(ipoint),') in buffer!'
            end if
            load_fw_points(:, ipoint) = load_fw_points(:,ipoint) &
                                      + azim_factor(rotmesh_phi(ipoint), &
                                                    source_params%mij, isim) &
                                      * real(utemp, kind=dp)
            !load_fw_points(:, ipoint) = utemp
            !print *, 'isim', isim, '; azim. factor:', azim_factor(rotmesh_phi(ipoint),&
            !                                                      source_params%mij, isim)
        end do !isim
        !read(*,*)
        write(1002,'(7(ES15.6))') rotmesh_s(ipoint), rotmesh_z(ipoint), rotmesh_phi(ipoint), &
                      azim_factor(rotmesh_phi(ipoint), source_params%mij, 1), &
                      azim_factor(rotmesh_phi(ipoint), source_params%mij, 2), &
                      azim_factor(rotmesh_phi(ipoint), source_params%mij, 3), &
                      azim_factor(rotmesh_phi(ipoint), source_params%mij, 4)
    end do !ipoint

end function load_fw_points

!-------------------------------------------------------------------------------
subroutine load_seismogram(this, receivers, src)
!< This function loads a seismogram 
   class(semdata_type)      :: this
   type(rec_param_type)     :: receivers(:)
   type(src_param_type)     :: src
   real(kind=dp)            :: seismogram_disp(this%ndumps)
   real(kind=dp)            :: seismogram_velo(this%ndumps)
   real(kind=sp)            :: utemp(this%ndumps)
   real(kind=dp)            :: Mij_scale(6), mij_prefact(4)
   integer                  :: reccomp, isurfelem, irec, isim, nrec

   if (.not.this%meshes_read) then
       print *, 'ERROR in load_seismogram(): Meshes have not been read yet'
       print *, 'Call read_meshes() before build_kdtree!'
       stop
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
         stop
      end select
      
      isurfelem = minloc( abs(this%fwdmesh%theta*deg2rad - receivers(irec)%theta), 1 )
      write(lu_out,'(A,F8.4,A,I5,A,F8.4)') 'Receiver with theta ', receivers(irec)%theta/deg2rad, &
                                    ' has element ', isurfelem, &
                                    ' with theta: ', this%fwdmesh%theta(isurfelem)
      
      seismogram_disp = 0.0
      seismogram_velo = 0.0

      write(lu_out,'(A,4(E12.4))') 'Mij prefactors: ', mij_prefact
      
      do isim = 1, this%nsim_fwd
         write(lu_out,'(A,I1,A,I5,A,I2,A,I6)') 'Sim: ', isim, ' Read element', isurfelem, &
                                               ', component: ', reccomp, ', no of samples:', this%ndumps
         ! Displacement seismogram
         call check( nf90_get_var( ncid   = this%fwd(isim)%surf,        & 
                                   varid  = this%fwd(isim)%seis_disp,   &
                                   start  = [1, reccomp, isurfelem],    &
                                   count  = [this%ndumps, 1, 1],        &
                                   values = utemp) )
      
         seismogram_disp = real(utemp, kind=dp) * mij_prefact(isim) + seismogram_disp

         ! Velocity seismogram
         call check( nf90_get_var( ncid   = this%fwd(isim)%surf,        & 
                                   varid  = this%fwd(isim)%seis_velo,   &
                                   start  = [1, reccomp, isurfelem],    &
                                   count  = [this%ndumps, 1, 1],        &
                                   values = utemp) )
      
         seismogram_velo = real(utemp, kind=dp) * mij_prefact(isim) + seismogram_velo
      end do

      this%dispseis(:, irec) = seismogram_disp(1:this%ndumps)
      this%veloseis(:, irec) = seismogram_velo(1:this%ndumps)
   end do
  

end subroutine load_seismogram

!-------------------------------------------------------------------------------
function load_bw_points(this, coordinates, receiver)
    class(semdata_type)                            :: this
    real(kind=dp), intent(in)                      :: coordinates(:,:)
    type(rec_param_type)                           :: receiver
    real(kind=dp)                                  :: load_bw_points(this%ndumps, size(coordinates,2))
    real(kind=sp)                                  :: utemp(this%ndumps)
    type(kdtree2_result), allocatable              :: nextpoint(:)

    integer, dimension(size(coordinates,2))        :: pointid, idx
    integer                                        :: npoints
    integer                                        :: ipoint, status
    real(kind=dp)                     :: rotmesh_s(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_phi(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_z(size(coordinates,2))

    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_bw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       stop 2
    end if
    npoints = size(coordinates,2)

    !allocate(load_bw_points(this%bwd(1)%ndumps,npoints))
    
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
        !print *, 'Original coordinates: ', coordinates(:,ipoint)
        !print *, 'Coordinates:    ', rotmesh_s(ipoint), rotmesh_z(ipoint), ', next pointid: ', pointid(ipoint)
        !print *, 'CO of SEM point:', this%bwdmesh%s(pointid(ipoint)), this%bwdmesh%z(pointid(ipoint))

        write(1001, '(7(ES15.6), I8)') coordinates(:, ipoint), rotmesh_s(ipoint), rotmesh_z(ipoint), &
                                       this%bwdmesh%s(pointid(ipoint)), this%bwdmesh%z(pointid(ipoint)), pointid(ipoint)

    end do
    
    load_bw_points(:,:) = 0.0
    
    do ipoint = 1, npoints
        
        status = this%bwd(1)%buffer%get(pointid(ipoint), utemp)
        if (status.ne.0) then
           !write(*,*) 'Did not find point', ipoint, ' in buffer, rereading'
           !if (this%bwd(1)%ordered_output) then
           !   call check( nf90_get_var( ncid   = this%bwd(1)%snap,        & 
           !                             varid  = this%bwd(1)%straintrace, &
           !                             start  = [1, pointid(ipoint)],    &
           !                             count  = [this%bwd(1)%ndumps, 1], &
           !                             values = utemp) )
           !else
              call check( nf90_get_var( ncid   = this%bwd(1)%snap,        & 
                                        varid  = this%bwd(1)%straintrace, &
                                        start  = [pointid(ipoint), 1],    &
                                        count  = [1, this%bwd(1)%ndumps], &
                                        values = utemp) )
           !end if
           status = this%bwd(1)%buffer%put(pointid(ipoint), utemp)
        !else
        !   write(*,*) 'Found point', ipoint, ' (',pointid(ipoint),') in buffer!'
        end if


        select case(receiver%component)
        case('Z')
            load_bw_points(:,(ipoint)) =                       real(utemp, kind=dp) / this%bwd(1)%amplitude
        case('R')
            load_bw_points(:,(ipoint)) =   dcos(rotmesh_phi) * real(utemp, kind=dp) / this%bwd(1)%amplitude
        case('T')
            load_bw_points(:,(ipoint)) = - dsin(rotmesh_phi) * real(utemp, kind=dp) / this%bwd(1)%amplitude 
        end select

        

    end do !ipoint

end function load_bw_points

!-------------------------------------------------------------------------------
function azim_factor(phi, mij, isim)
    real(kind=dp), intent(in)    :: phi
    real(kind=dp), intent(in)    :: mij(6)
    integer, intent(in)          :: isim
    real(kind=dp)                :: azim_factor


    select case(isim)
    case(1) ! Mzz
       azim_factor = Mij(1)
       
    case(2) ! Mxx+Myy
       azim_factor = 0.5*(Mij(2)+Mij(3))
       
    case(3) ! dipole
       azim_factor = Mij(4)*dcos(phi) + Mij(5)*dsin(phi)
       
    case(4) ! quadrupole
       azim_factor = 0.5*(Mij(2)-Mij(3))*dcos(2.d0*phi) + Mij(6)*dsin(2.d0*phi)

    case default
       write(6,*) myrank,': unknown number of simulations',isim
    end select

end function

!-------------------------------------------------------------------------------
subroutine build_kdtree(this)
    class(semdata_type)        :: this
    real(kind=sp), allocatable :: mesh(:,:)
    integer                    :: ipoint

    if (.not.this%meshes_read) then
        print *, 'ERROR in build_kdtree(): Meshes have not been read yet'
        print *, 'Call read_meshes() before build_kdtree!'
        stop
    end if

    allocate(mesh(2, this%fwdmesh%npoints))
    mesh = transpose(reshape([this%fwdmesh%s, this%fwdmesh%z], [this%fwdmesh%npoints, 2]))

    !print *, 'Maximum s:', maxval(this%fwdmesh%s)
    !print *, 'Minimum s:', minval(this%fwdmesh%s)
    !print *, 'Maximum z:', maxval(this%fwdmesh%z)
    !print *, 'Minimum z:', minval(this%fwdmesh%z)

    write(lu_out,*) ' Building forward KD-Tree'
    ! KDtree in forward field
    this%fwdtree => kdtree2_create(mesh,              &
                                   dim = 2,           &
                                   sort = .true.,     &
                                   rearrange = .true.)
    
    deallocate(mesh)                           

    ! KDtree in backward field
    write(lu_out,*) ' Building backward KD-Tree'
    allocate(mesh(2, this%bwdmesh%npoints))
    mesh = transpose(reshape([this%bwdmesh%s, this%bwdmesh%z], [this%bwdmesh%npoints, 2]))
    this%bwdtree => kdtree2_create(mesh,              &
                                   dim = 2,           &
                                   sort = .true.,     &
                                   rearrange = .true.)
    deallocate(mesh)                           

end subroutine build_kdtree

!-------------------------------------------------------------------------------

subroutine read_meshes(this)
    use netcdf
    class(semdata_type)    :: this
    integer                :: ncvarid_mesh_s, ncvarid_mesh_z
    integer                :: surfdimid, ncvarid_theta
    integer                :: ielem
   
    if (.not.this%files_open) then
        print *, 'ERROR in read_meshes(): Files have not been opened!'
        print *, 'Call open_files() before read_meshes()'
        stop
    end if

    ! Forward SEM mesh
    write(lu_out,*) '  Read SEM mesh from first forward simulation'
    
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

    call check( nf90_get_var(ncid   = this%fwd(1)%mesh,       &
                             varid  = ncvarid_mesh_s,         &
                             start  = [1],                    &
                             count  = [this%fwdmesh%npoints], &
                             values = this%fwdmesh%s )) 
    call check( nf90_get_var(ncid   = this%fwd(1)%mesh,       &
                             varid  = ncvarid_mesh_z,         &
                             start  = [1],                    &
                             count  = [this%fwdmesh%npoints], &
                             values = this%fwdmesh%z )) 
   
    ! Backward SEM mesh                     
    write(lu_out,*) 'Read SEM mesh from first backward simulation'
    
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
                             varid  = ncvarid_mesh_s,   &
                             start  = [1],                    &
                             count  = [this%bwdmesh%npoints], &
                             values = this%bwdmesh%s )) 
    call check( nf90_get_var(ncid   = this%bwd(1)%mesh, &
                             varid  = ncvarid_mesh_z,   &
                             start  = [1],                    &
                             count  = [this%bwdmesh%npoints], &
                             values = this%bwdmesh%z )) 

   
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
   call check( nf90_get_var( ncid   = this%fwd(1)%surf,    &
                             varid  = ncvarid_theta,       &
                             values = this%fwdmesh%theta) )

   
   ! Backward mesh
   call check( nf90_inq_dimid( ncid  = this%bwd(1)%surf, &
                               name  = 'surf_elems',     &
                               dimid = surfdimid) ) 

   call check( nf90_inquire_dimension(ncid  = this%bwd(1)%surf,        & 
                                      dimid = surfdimid,               &
                                      len   = this%fwdmesh%nsurfelem) )

   call  getvarid(ncid  = this%bwd(1)%surf, &
                  name  = "elem_theta",     &
                  varid = ncvarid_theta) 

   allocate( this%bwdmesh%theta(this%fwdmesh%nsurfelem) )
   call check( nf90_get_var( ncid   = this%bwd(1)%surf,    &
                             varid  = ncvarid_theta,       &
                             values = this%bwdmesh%theta) )
                             

                                   
    this%meshes_read = .true.

end subroutine read_meshes
 

!-------------------------------------------------------------------------------
subroutine check(status)
! Translates netcdf error codes into error messages

  implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop 2
     !call tracebackqq()
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
        write(6,100) myrank, trim(name), ncid
        stop
    elseif (myrank.eq.0) then
        if (verbose>0) write(lu_out,101) trim(name), ncid, varid
    end if
100 format('ERROR: CPU ', I4, ' could not find variable: ''', A, ''' in NCID', I7)
101 format('  Variable ''', A, ''' found in NCID', I7, ', has ID:', I7)
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine nc_open_for_read(filename, ncid)
   character(len=*), intent(in)  :: filename
   integer, intent(out)          :: ncid
   character(len=512)            :: fmtstring
   integer                       :: status

   status = nf90_open(path     = filename,              &
                      mode     = nf90_nowrite,          &
                      ncid     = ncid)

   if (status.ne.nf90_noerr) then
      fmtstring = "('ERROR: CPU ', I4, ' tried to open file ''', A, ''', but could not find it')"
      print fmtstring, myrank, trim(filename)
      stop
   end if

end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
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
        stop
    elseif (verbose>1) then
        write(lu_out,101) trim(name), ncid, grp_ncid
        call flush(lu_out)
    end if
100 format('ERROR: CPU ', I4, ' could not find group: ''', A, ''' in NCID', I7)
101 format('    Group ''', A, ''' found in NCID', I7, ', has ID:', I7)
end subroutine getgrpid
!-------------------------------------------------------------------------------

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
!> Read NetCDF attribute of type Double
subroutine nc_read_att_dble(attribute_value, attribute_name, nc)
  character(len=*),  intent(in)     :: attribute_name
  real(kind=dp), intent(out)        :: attribute_value
  type(ncparamtype), intent(in)     :: nc
  integer                           :: status

  status = nf90_get_att(nc%ncid, NF90_GLOBAL, attribute_name, attribute_value)
  if (status.ne.NF90_NOERR) then
      write(6,*) 'Could not find attribute ', trim(attribute_name)
      write(6,*) ' in NetCDF file ', trim(nc%meshdir), '/Data/axisem_output.nc4'
      write(6,*) ' with NCID: ', nc%ncid
      stop 2
  end if
end subroutine nc_read_att_dble
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine rotate_frame_rd(npts, srd, phird, zrd, rgd, phigr, thetagr)

    implicit none
    integer, intent(in)                            :: npts
    !< Number of points to rotate
    
    real(kind=dp), dimension(3, npts), intent(in)  :: rgd
    !< Coordinates to rotate (in x, y, z)
    
    real(kind=dp), intent(in)                      :: phigr, thetagr
    !< Rotation angles phi and theta
    
    real(kind=dp), dimension(npts), intent(out)    :: srd, zrd, phird
    !< Rotated coordinates (in s, z, phi)

    real(kind=dp), dimension(npts)                 :: xp, yp, zp
    real(kind=dp), dimension(npts)                 :: xp_cp, yp_cp, zp_cp
    real(kind=dp)                                  :: phi_cp
    integer                                        :: ipt

    !!first rotation (longitude)
    xp_cp =  rgd(1,:) * dcos(phigr) + rgd(2,:) * dsin(phigr)
    yp_cp = -rgd(1,:) * dsin(phigr) + rgd(2,:) * dcos(phigr)
    zp_cp =  rgd(3,:)

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
    enddo
end subroutine rotate_frame_rd
!-------------------------------------------------------------------------------

!=============================================================================
end module
