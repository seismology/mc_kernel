!******************************************************************************
!
!    This file is part of:
!    MC Kernel: Calculating seismic sensitivity kernels on unstructured meshes
!    Copyright (C) 2016 Simon Staehler, Martin van Driel, Ludwig Auer
!
!    You can find the latest version of the software at:
!    <https://www.github.com/tomography/mckernel>
!
!    MC Kernel is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MC Kernel is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MC Kernel. If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

!=========================================================================================
module readfields
    use global_parameters, only            : sp, dp, pi, deg2rad, rad2deg, verbose, lu_out, &
                                             myrank, long,                                  &
                                             id_buffer, id_netcdf, id_rotate,               &
                                             id_load_strain, id_kdtree, id_calc_strain,     &
                                             id_find_point_fwd, id_find_point_bwd, id_lagrange, &
                                             id_allocate, id_allocate_2

    use source_class,      only            : src_param_type
    use receiver_class,    only            : rec_param_type
    use buffer,            only            : buffer_type
    use clocks_mod,        only            : tick
    use commpi,            only            : pabort
    use nc_routines,       only            : getgrpid, getvarid, nc_open_for_read,  &
                                             nc_getvar, nc_getvar_by_name, check

    use receivers_rdbm,    only            : receivers_rdbm_type

    use rotations,         only            : azim_factor, azim_factor_bw,                   &
                                             rotate_symm_tensor_voigt_src_to_xyz,           &
                                             rotate_symm_tensor_voigt_xyz_src_to_xyz_earth, &
                                             rotate_symm_tensor_voigt_xyz_earth_to_xyz_src, &
                                             rotate_symm_tensor_voigt_xyz_to_src,           &    
                                             rotate_frame_rd
    use kdtree2_module,    only            : kdtree2

    use interpolate_mesh,  only            : parameter_interpolator

    implicit none
    private
    public                                 :: semdata_type, meshtype, get_chunk_bounds, dampen_field, nelem_to_read_max

    integer, parameter                     :: min_file_version = 3
    integer                                :: npoints_max      = 10
    integer, protected, save               :: nelem_to_read_max      !< How many elements to read 
                                                                     !! for the merged database case.
    integer, parameter                     :: nmodel_parameters_sem_file = 6 !< For the anisotropic
                                                                             !! case. Will increase
                                                                             !! for attenuation
    real(kind=sp), allocatable, save       :: u_batch(:,:,:,:,:)
    real(kind=sp), allocatable, save       :: utemp_fwd(:,:,:,:,:)
    real(kind=sp), allocatable, save       :: utemp_bwd(:,:,:,:,:)
    real(kind=dp), allocatable, save       :: strain_fwd(:,:,:,:,:,:)
    real(kind=dp), allocatable, save       :: strain_bwd(:,:,:,:,:,:)
    real(kind=dp), allocatable, save       :: straintrace_fwd(:,:,:,:,:)
    real(kind=dp), allocatable, save       :: straintrace_bwd(:,:,:,:,:)

    type meshtype
        real(kind=sp), allocatable         :: s(:), z(:)            !< Coordinates of all GLL points
        real(kind=sp), allocatable         :: s_mp(:), z_mp(:)      !< Coordinates of element midpoints
        integer, allocatable               :: corner_point_ids(:,:) !< (4,nelem)
        integer, allocatable               :: eltype(:)             !< (nelem)
        type(parameter_interpolator)       :: vp, vs, rho           !< Model parameters
        type(parameter_interpolator)       :: lambda, mu            !< Elastic parameters
        type(parameter_interpolator)       :: phi, xi, eta          !< Anisotropic parameters
        real(kind=sp)                      :: vp_max, vs_max        !< Maximum values of vp, vs in the model
        integer, allocatable               :: isaxis(:)             !< Is this point at the axis?
        integer, allocatable               :: gll_point_ids(:,:,:)  !< IDs of GLL points for this element
        integer                            :: npoints, nelem
    end type

    type ncparamtype
        integer                            :: ncid
        integer                            :: file_version
        real(kind=dp)                      :: planet_radius
        real(kind=dp)                      :: rmin, rmax
        real(kind=dp)                      :: colatmin, colatmax
        integer                            :: snap, mesh, seis        ! Group IDs
        integer                            :: strainvarid(6)          ! Variable IDs
        integer                            :: displvarid(3)           ! Variable IDs
        integer                            :: mergedvarid             ! Variable IDs
        integer                            :: chunk_gll
        integer                            :: count_error_pointoutside
        character(len=200)                 :: meshdir
        integer                            :: ndumps, nseis, ngll, npol
        integer                            :: source_shift_samples    
        real(kind=dp)                      :: source_shift_t
        character(len=10)                  :: source_type
        character(len=10)                  :: stf_type
        character(len=10)                  :: excitation_type
        real(kind=sp), allocatable         :: stf(:), stf_d(:)
        type(buffer_type)                  :: buffer_strain
        type(buffer_type)                  :: buffer_disp
        type(buffer_type)                  :: buffer
        real(kind=dp)                      :: dt
        real(kind=dp)                      :: amplitude
        real(kind=dp)                      :: source_depth
        real(kind=dp), public, allocatable :: G1(:,:), G1T(:,:)
        real(kind=dp), public, allocatable :: G2(:,:), G2T(:,:)
        real(kind=dp), public, allocatable :: G0(:)
        real(kind=dp), public, allocatable :: gll_points(:), glj_points(:)
        logical                            :: merged = .false.
        integer                            :: nsim_merged = 0
        integer                            :: npoints, nelem
        logical                            :: parallel_read
    end type

    type semdata_type
        private

        integer, public                    :: nsim_fwd, nsim_bwd
        integer, public                    :: nfiles_fwd, nfiles_bwd
        type(ncparamtype), allocatable, public     :: fwd(:)
        type(ncparamtype), allocatable     :: bwd(:)

        type(kdtree2), pointer, private    :: fwdtree, bwdtree        !< Contain all points
        type(kdtree2), pointer, private    :: fwdtree_mp, bwdtree_mp  !< Contain only midpoints
        type(meshtype)                     :: fwdmesh, bwdmesh

        logical, private                   :: params_set   = .false.
        logical, private                   :: files_open   = .false.
        logical, private                   :: meshes_read  = .false.
        logical, private                   :: kdtree_built = .false.
        
        character(len=32)                  :: strain_type  !< full tensor or straintrace
        integer                            :: ndim     !< Number of dimensions which has to be read to calculate 
                                                       !! Kernel on parameter model_param

        logical, private                   :: merged_fwd, merged_bwd !< Whether this object has 
                                                                     !! its data stored in a 
                                                                     !! merged database file
        logical, private                   :: parallel_read

        real(kind=dp), public              :: dt
        integer,       public              :: ndumps, decimate_factor
        integer,       public              :: nseis 
        integer,       public              :: npol
        real(kind=dp), public, allocatable :: G1(:,:), G1T(:,:)
        real(kind=dp), public, allocatable :: G2(:,:), G2T(:,:)
        real(kind=dp), public, allocatable :: G0(:)
        real(kind=dp), public, allocatable :: gll_points(:), glj_points(:)
        real(kind=dp), public              :: windowlength
        real(kind=dp), public              :: desired_source_depth
        real(kind=dp), public              :: timeshift_fwd, timeshift_bwd
        real(kind=dp), public              :: amplitude_fwd, amplitude_bwd
        real(kind=dp), public, allocatable :: seis(:,:) 
        real(kind=dp), public, allocatable :: stf_fwd(:), stf_bwd(:)
        real(kind=dp), public, allocatable :: stf_d_fwd(:), stf_d_bwd(:)
        integer                            :: strain_buffer_size
        integer                            :: displ_buffer_size
         
        real(kind=dp), dimension(3,3)      :: rot_mat, trans_rot_mat

        contains 
            procedure, pass                :: get_ndim 
            procedure, pass                :: get_mesh 
            procedure, pass                :: get_v_max
            procedure, pass                :: set_params
            procedure, pass                :: open_files
            procedure, pass                :: close_files
            procedure, pass                :: check_consistency
            procedure, pass                :: read_meshes
            procedure, pass, private       :: build_kdtree
            procedure, pass                :: load_fw_points
            procedure, pass                :: load_fw_points_rdbm
            procedure, pass                :: load_bw_points
            procedure, pass                :: load_model_coeffs
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
function get_mesh(this, fwd_or_bwd)
   class(semdata_type)           :: this
   character(len=3), intent(in)  :: fwd_or_bwd
   type(meshtype)                :: get_mesh

   select case(fwd_or_bwd)
   case('fwd')
      get_mesh = this%fwdmesh
   case('bwd')
      get_mesh = this%bwdmesh
   case default
      write(*,*) 'ERROR: get_mesh can only get "fwd" or "bwd" mesh!'
      call pabort(do_traceback=.false.)
   end select
end function get_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_v_max(this)
  ! Return the maximum p and s velocities of the forward mesh
  class(semdata_type)            :: this
  real(kind=sp)                  :: get_v_max(2)

  get_v_max(1) = this%fwdmesh%vp_max
  get_v_max(2) = this%fwdmesh%vs_max

end function get_v_max
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine set_params(this, fwd_dir, bwd_dir, strain_buffer_size, displ_buffer_size, &
                      strain_type, desired_source_depth, npoints, parallel_read)
    class(semdata_type)            :: this
    character(len=512), intent(in) :: fwd_dir, bwd_dir
    integer,            intent(in) :: strain_buffer_size
    integer,            intent(in) :: displ_buffer_size
    integer,            intent(in) :: npoints
    character(len=*),   intent(in) :: strain_type
    logical,            intent(in) :: parallel_read
    real(kind=dp),      intent(in) :: desired_source_depth

    this%strain_buffer_size = strain_buffer_size
    this%displ_buffer_size = displ_buffer_size

    this%desired_source_depth = desired_source_depth

    npoints_max = npoints

    ! Get simulation type of forward simulation and whether it was stored in 
    ! separate directories or in one merged file
    call get_simulation_type(fwd_dir, this%merged_fwd, &
                             this%nsim_fwd, this%nfiles_fwd)
    allocate(this%fwd(this%nfiles_fwd))
    if (this%merged_fwd) this%fwd(1)%nsim_merged = this%nsim_fwd
    call set_simulation_paths(this%fwd, fwd_dir, this%nfiles_fwd)

    ! Get simulation type of backward simulation and whether it was stored in 
    ! separate directories or in one merged file
    call get_simulation_type(bwd_dir, this%merged_bwd, &
                             this%nsim_bwd, this%nfiles_bwd)
    allocate(this%bwd(this%nfiles_bwd))
    if (this%merged_bwd) this%bwd(1)%nsim_merged = this%nsim_bwd
    call set_simulation_paths(this%bwd, bwd_dir, this%nfiles_bwd)

    this%strain_type = strain_type

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

    this%parallel_read = parallel_read
    write(lu_out, *) 'Using Parallel NetCDF4: ', this%parallel_read

    call flush(lu_out)
    this%params_set = .true.

end subroutine set_params
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine open_files(this)
    class(semdata_type)              :: this
    integer                          :: ifile
    character(len=200)               :: filename

    if (.not.this%params_set) then
        print *, 'ERROR in open_files(): Parameters have to be set first'
        print *, 'Call set_params before open_files()'
        call pabort
    end if

    do ifile = 1, this%nfiles_fwd
        ! Forward wavefield
        if (this%merged_fwd) then
          filename=trim(this%fwd(ifile)%meshdir)//'/merged_instaseis_db.nc4'
        else
          filename=trim(this%fwd(ifile)%meshdir)//'/Data/ordered_output.nc4'
        end if

        call open_file_read_varids_and_attributes(this%fwd(ifile),   &
                                                  filename,          &
                                                  this%merged_fwd,   &
                                                  this%parallel_read)
    end do
        
    call flush(lu_out)

    do ifile = 1, this%nfiles_bwd
        ! Backward wavefield
        if (this%merged_bwd) then
          filename=trim(this%bwd(ifile)%meshdir)//'/merged_instaseis_db.nc4'
        else
          filename=trim(this%bwd(ifile)%meshdir)//'/Data/ordered_output.nc4'
        end if

        call open_file_read_varids_and_attributes(this%bwd(ifile),   &
                                                  filename,          &
                                                  this%merged_bwd,   &
                                                  this%parallel_read)
        
    end do


    call flush(lu_out)
    call this%check_consistency()

    this%files_open = .true.

    !@TODO memory could be used more efficient for monopole sources in the buffers
    !Initialize Buffers. 
    do ifile = 1, this%nfiles_fwd
      if (this%merged_fwd) then
        call init_merged_buffer(this%fwd(ifile), this%displ_buffer_size, &
                                this%strain_buffer_size, this%strain_type)
      else
        call init_disp_only_buffer(this%fwd(ifile), this%displ_buffer_size, &
                                   this%strain_buffer_size, this%strain_type)
      end if
      this%fwd(ifile)%count_error_pointoutside = 0
    end do

    do ifile = 1, this%nfiles_bwd
      if (this%merged_bwd) then
        call init_merged_buffer(this%bwd(ifile), this%displ_buffer_size, &
                                this%strain_buffer_size, this%strain_type)
      else
        call init_disp_only_buffer(this%bwd(ifile), this%displ_buffer_size, &
                                   this%strain_buffer_size, this%strain_type)
      end if
      this%bwd(ifile)%count_error_pointoutside = 0
    end do

    call flush(lu_out)

end subroutine open_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_merged_buffer(nc_obj, displ_buffer_size, strain_buffer_size, &
                                 strain_type)
  type(ncparamtype)             :: nc_obj 
  integer, intent(in)           :: strain_buffer_size, displ_buffer_size
  character(len=*), intent(in)  :: strain_type
  integer                       :: status

  status = nc_obj%buffer_disp%init(displ_buffer_size, &
                                   'real',            &
                                   nc_obj%ndumps,     &
                                   nc_obj%npol+1,     &
                                   nc_obj%npol+1,     &
                                   nint(nc_obj%nsim_merged*2.5))

  ! Should be at least one, even for small buffers
  ! Select such that not every read overwrites the complete buffer
  ! 20 is a heuristic value
  nelem_to_read_max = max(int(displ_buffer_size / 20), 1) 

  if (.not.allocated(u_batch)) then
    allocate(u_batch(nc_obj%ndumps,                  &
                     0:nc_obj%npol,                  &
                     0:nc_obj%npol,                  &
                     1:10,                           &
                     nelem_to_read_max))
    u_batch = 0                 
  end if

  select case(strain_type)
  case('straintensor_trace')
    status = nc_obj%buffer_strain%init(strain_buffer_size, &
                                       'dble',             &
                                       nc_obj%ndumps,      &
                                       nc_obj%npol+1,      &
                                       nc_obj%npol+1,      &
                                       nc_obj%nsim_merged)
    select case(nc_obj%nsim_merged)
    case(4)
      allocate(straintrace_fwd(1:nc_obj%ndumps, &
                               0:nc_obj%npol,   &
                               0:nc_obj%npol,   &
                               1:nc_obj%nsim_merged, &
                               1:npoints_max))
      straintrace_fwd = 0
    case(2)
      allocate(straintrace_bwd(1:nc_obj%ndumps, &
                               0:nc_obj%npol,   &
                               0:nc_obj%npol,   &
                               1:nc_obj%nsim_merged, &
                               1:npoints_max))
      straintrace_bwd = 0
    case default 
      print *, 'Unknown value for nsim_merged: ', nc_obj%nsim_merged
      call pabort()
    end select
  case('straintensor_full')
    status = nc_obj%buffer_strain%init(strain_buffer_size, &
                                       'dble',             &
                                       nc_obj%ndumps,      &
                                       nc_obj%npol+1,      &
                                       nc_obj%npol+1,      &
                                       6,                  &
                                       nc_obj%nsim_merged)
  case default 
    print *, 'Unknown value for strain_type:', strain_type
    call pabort()
  end select

  select case(nc_obj%nsim_merged)
  case(4)
    allocate(strain_fwd(1:nc_obj%ndumps, &
                        0:nc_obj%npol,   &
                        0:nc_obj%npol,   &
                        6,               &
                        1:nc_obj%nsim_merged, &
                        1:npoints_max))
    allocate(utemp_fwd(1:nc_obj%ndumps, &
                       0:nc_obj%npol,   &
                       0:nc_obj%npol,   &
                       1:10,            &
                       1:npoints_max))
    strain_fwd = 0
    utemp_fwd = 0
  case(2)
    allocate(strain_bwd(1:nc_obj%ndumps, &
                        0:nc_obj%npol,   &
                        0:nc_obj%npol,   &
                        6,               &
                        1:nc_obj%nsim_merged, &
                        1:npoints_max))
    allocate(utemp_bwd(1:nc_obj%ndumps, &
                       0:nc_obj%npol,   &
                       0:nc_obj%npol,   &
                       1:5,             &
                       1:npoints_max))
    strain_bwd = 0
    utemp_bwd = 0
  end select
end subroutine init_merged_buffer
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_disp_only_buffer(nc_obj, displ_buffer_size, strain_buffer_size, &
                                 strain_type)
  type(ncparamtype)             :: nc_obj 
  integer, intent(in)           :: strain_buffer_size, displ_buffer_size
  character(len=*), intent(in)  :: strain_type
  integer                       :: status

  status = nc_obj%buffer_disp%init(displ_buffer_size, 'real', nc_obj%ndumps, 3)
  select case(strain_type)
  case('straintensor_trace')
    status = nc_obj%buffer_strain%init(strain_buffer_size, &
                                       'dble',             &
                                       nc_obj%ndumps,      &
                                       nc_obj%npol+1,      &
                                       nc_obj%npol+1)
  case('straintensor_full')
    status = nc_obj%buffer_strain%init(strain_buffer_size, &
                                       'dble',             &
                                       nc_obj%ndumps,      &
                                       nc_obj%npol+1,      &
                                       nc_obj%npol+1,      &
                                       6)
  end select

end subroutine init_disp_only_buffer
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine open_file_read_varids_and_attributes(nc_obj, filename, merged, parallel_read)
  use netcdf, only              : nf90_inq_varid, nf90_inquire_variable, &
                                  nf90_get_var, NF90_NOERR
  use commpi, only              : MPI_COMM_NODE, mpi_info
  type(ncparamtype)             :: nc_obj 
  character(len=*), intent(in)  :: filename
  logical, intent(in)           :: merged, parallel_read
  integer                       :: idisplvar
  integer                       :: status, chunks(2), deflev
  integer                       :: stf_grp_id
  character(len=200)            :: format20, format21
  character(len=11)             :: nc_strain_varnamelist(6)
  character(len=11)             :: nc_displ_varnamelist(3)
  real(kind=sp)                 :: temp

  nc_strain_varnamelist = ['strain_dsus', 'strain_dsuz', 'strain_dpup', &
                           'strain_dsup', 'strain_dzup', 'straintrace']
          
  nc_displ_varnamelist  = ['disp_s     ', 'disp_p     ', 'disp_z     ']

  format20 = "('  Trying to open NetCDF file ', A, ' on CPU ', I5)"
  format21 = "('  Succeded,  has NCID ', I6, ', Snapshots group NCID: ', I6)"

  if (verbose>0) write(lu_out,format20) trim(filename), myrank
  if (parallel_read) then
    call nc_open_for_read(filename = filename,              &
                          ncid     = nc_obj%ncid,           &
                          comm     = MPI_COMM_NODE,         &
                          info     = mpi_info) 
  else
    call nc_open_for_read(filename = filename,              &
                          ncid     = nc_obj%ncid) 
  end if

  nc_obj%parallel_read = parallel_read

  call nc_read_att_int(nc_obj%file_version, 'file version', nc_obj)

  if (nc_obj%file_version < min_file_version) then
     print *, 'ERROR: AxiSEM NetCDF file too old. '
     print *, 'Filename: ', trim(nc_obj%meshdir)//'/Data/ordered_output.nc4'
     print *, 'Minimum file version: ', min_file_version, &
              ', found: ', nc_obj%file_version
            
     call pabort(do_traceback=.false.)
  endif

  call nc_read_att_char(nc_obj%stf_type, 'source time function', nc_obj)

  call nc_read_att_dble(nc_obj%source_depth, 'source depth in km', nc_obj)

  call nc_read_att_dble(nc_obj%planet_radius, 'planet radius', nc_obj)
  nc_obj%planet_radius = nc_obj%planet_radius * 1d3

  call nc_read_att_dble(nc_obj%rmin, 'kernel wavefield rmin', nc_obj)
  call nc_read_att_dble(nc_obj%rmax, 'kernel wavefield rmax', nc_obj)
  nc_obj%rmin = nc_obj%rmin * 1d3
  nc_obj%rmax = nc_obj%rmax * 1d3

  call nc_read_att_dble(nc_obj%colatmin, 'kernel wavefield colatmin', nc_obj)
  call nc_read_att_dble(nc_obj%colatmax, 'kernel wavefield colatmax', nc_obj)

  call nc_read_att_char(nc_obj%source_type, &
                        'source type', &
                         nc_obj)

  call nc_read_att_char(nc_obj%excitation_type, &
                        'excitation type', &
                         nc_obj)

  call nc_read_att_int(nc_obj%npol, 'npol', &
                       nc_obj)

  if (.not.merged) then                     
    call getgrpid(  ncid     = nc_obj%ncid,   &
                    name     = "Snapshots",           &
                    grp_ncid = nc_obj%snap)

    do idisplvar = 1, 3
        status = nf90_inq_varid(ncid  = nc_obj%snap,                  &
                                name  = nc_displ_varnamelist(idisplvar),      &
                                varid = nc_obj%displvarid(idisplvar)) 
        
        if (status.ne.NF90_NOERR) then
            nc_obj%displvarid(idisplvar) = -1
            if (idisplvar == 1) then
                print *, 'Did not find variable ''disp_s'' in NetCDF file'
                call pabort
            end if
        end if
    end do
    call check(nf90_inquire_variable(ncid       = nc_obj%snap,   &
                                     varid      = nc_obj%displvarid(1), &
                                     chunksizes = chunks, &
                                     deflate_level = deflev) )


    write(lu_out, "('  File', A, ', Chunksizes:', 2(I7), ', deflate level: ', I2)") &
          trim(filename), chunks, deflev

    nc_obj%chunk_gll = chunks(1)

    if (verbose>0) write(lu_out,format21) nc_obj%ncid, nc_obj%snap 
    
    stf_grp_id = nc_obj%snap

  else ! Merged snapshot file

    status = nf90_inq_varid(ncid  = nc_obj%ncid,         &
                            name  = 'merged_snapshots',       &
                            varid = nc_obj%mergedvarid)

    stf_grp_id = nc_obj%ncid

  end if

  call nc_read_att_int(    nc_obj%ndumps,             &
                           'number of strain dumps',  &
                           nc_obj)

  call nc_getvar_by_name(  ncid    = stf_grp_id,    &
                           varname = 'stf_dump',    &
                           limits  = [0., 1e36],    &
                           values  = nc_obj%stf,    &
                           collective = parallel_read)

  call nc_getvar_by_name(  ncid    = stf_grp_id,    &
                           varname = 'stf_d_dump',  &
                           limits  = [-1e36, 1e36],    &
                           values  = nc_obj%stf_d ,    &
                           collective = parallel_read)

  call getgrpid(           ncid      = nc_obj%ncid,   &
                           name      = "Mesh",                &
                           grp_ncid  = nc_obj%mesh)

  call nc_read_att_dble(   nc_obj%dt,               &
                           'strain dump sampling rate in sec', &
                           nc_obj)

  call nc_read_att_int(    nc_obj%nseis,             &
                           'length of seismogram  in time samples', &
                           nc_obj)

  call nc_read_att_real(   temp, &
                           'source shift factor in sec',     &
                           nc_obj)
  nc_obj%source_shift_t = real(temp, kind=dp)
  
  call nc_read_att_int(    nc_obj%source_shift_samples,    &
                           'source shift factor for deltat_coarse',     &
                           nc_obj)
  
  call nc_read_att_real(   temp,      &
                           'scalar source magnitude',     &
                           nc_obj)
  nc_obj%amplitude = real(temp, kind=dp)


end subroutine open_file_read_varids_and_attributes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine close_files(this)
    use kdtree2_module, only : kdtree2_destroy
    use netcdf,         only : nf90_close
    class(semdata_type)     :: this
    integer                 :: status, ifile

    ! Destroy kdtree
    if (this%kdtree_built) then
      call kdtree2_destroy(this%fwdtree)
      call kdtree2_destroy(this%bwdtree)
      call kdtree2_destroy(this%fwdtree_mp)
      call kdtree2_destroy(this%bwdtree_mp)
    end if

    ! Free buffers

      do ifile = 1, this%nfiles_fwd
         status = nf90_close(this%fwd(ifile)%ncid)
         if (verbose>0) then
            write(lu_out,'(A,I1,A,F9.6)') ' Strain buffer efficiency fwd(', ifile, '): ',  &
                                     this%fwd(ifile)%buffer_strain%efficiency()
            write(lu_out,'(A,I1,A,F9.6)') ' Displ. buffer efficiency fwd(', ifile, '): ',  &
                                     this%fwd(ifile)%buffer_disp%efficiency()
         end if
         status = this%fwd(ifile)%buffer_strain%freeme()
         status = this%fwd(ifile)%buffer_disp%freeme()
      end do
      write(lu_out,'(A,I8)') ' Points outside of element (fwd): ', &
                             this%fwd(1)%count_error_pointoutside

      do ifile = 1, this%nfiles_bwd
         status = nf90_close(this%bwd(ifile)%ncid)
         if (verbose>0) then
            write(lu_out,'(A,I1,A,F9.6)') ' Strain buffer efficiency bwd(', ifile, '): ',  &
                                     this%bwd(ifile)%buffer_strain%efficiency()
            write(lu_out,'(A,I1,A,F9.6)') ' Displ. buffer efficiency bwd(', ifile, '): ',  &
                                     this%bwd(ifile)%buffer_disp%efficiency()
         end if
         status = this%bwd(ifile)%buffer_strain%freeme()
         status = this%bwd(ifile)%buffer_disp%freeme()
      end do
      write(lu_out,'(A,I8)') ' Points outside of element (bwd): ', &
                             this%bwd(1)%count_error_pointoutside

    deallocate(this%fwd)
    deallocate(this%bwd)

    if (allocated(u_batch)) deallocate(u_batch)
    if (allocated(utemp_fwd)) deallocate(utemp_fwd)
    if (allocated(utemp_bwd)) deallocate(utemp_bwd)
    if (allocated(strain_fwd)) deallocate(strain_fwd)
    if (allocated(strain_bwd)) deallocate(strain_bwd)
    if (allocated(straintrace_fwd)) deallocate(straintrace_fwd)
    if (allocated(straintrace_bwd)) deallocate(straintrace_bwd)

    call flush(lu_out)

end subroutine close_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine check_consistency(this)
    !< Checks consistency of the wavefield dumps
    !! and write agreed values into the semdata_type object "this"
    class(semdata_type)    :: this
    integer                :: ifile
    real(kind=dp)          :: dt_agreed
    character(len=512)     :: fmtstring, fmtstring_stf
    integer                :: ndumps_agreed, nseis_agreed, npol_agreed
    real(kind=dp)          :: source_shift_agreed_fwd, source_shift_agreed_bwd
    real(kind=dp)          :: amplitude_agreed_fwd, amplitude_agreed_bwd
    real(kind=dp), allocatable  :: stf_agreed_fwd(:), stf_d_agreed_fwd(:)
    real(kind=dp), allocatable  :: stf_agreed_bwd(:), stf_d_agreed_bwd(:)

    write(lu_out, *) 'Checking consistency of wavefield files...'
100 format('  Parameter ', A, ' has value ', I8)
101 format('  Parameter ', A, ' has value ', E15.8)

    allocate(stf_agreed_fwd(this%fwd(1)%ndumps))
    allocate(stf_d_agreed_fwd(this%fwd(1)%ndumps))
    allocate(stf_agreed_bwd(this%fwd(1)%ndumps))
    allocate(stf_d_agreed_bwd(this%fwd(1)%ndumps))

    ! Check whether the STF in AxiSEM was correct
    do ifile = 1, this%nfiles_fwd
      if (trim(this%fwd(ifile)%stf_type).ne.'gauss_0') then
        print *, 'ERROR: Invalid AxiSEM source time function: ', this%fwd(ifile)%stf_type
        print *, '       Please run AxiSEM with ''gauss_0'' to ensure correct units'
        print *, '       and avoid aliasing'
      end if
    end do
    do ifile = 1, this%nfiles_bwd
      if (trim(this%bwd(ifile)%stf_type).ne.'gauss_0') then
        print *, 'ERROR: Invalid AxiSEM source time function: ', this%bwd(ifile)%stf_type
        print *, '       Please run AxiSEM with ''gauss_0'' to ensure correct units'
        print *, '       and avoid aliasing'
      end if
    end do

    ! Check whether the depth in CMTSOLUTION is consistent with the depth of the AxiSEM fwd run
    do ifile = 1, this%nfiles_fwd
      if (this%fwd(ifile)%source_depth.ne.this%desired_source_depth) then
        print *, 'ERROR: Source depth in CMTSOLUTION and AxiSEM fwd run are inconsistent!'
        print *, '       Depth in CMTSOLUTION: ', this%desired_source_depth
        print *, '       Depth in AxiSEM run:  ', this%fwd(ifile)%source_depth
        stop
      end if
    end do


    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1)' 
    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1)' 

    ! Check whether the sampling period is the same in all files
    dt_agreed = this%fwd(1)%dt
    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the others")' 
    do ifile = 1, this%nfiles_fwd
       if (dt_agreed.ne.this%fwd(ifile)%dt) then
          write(*,fmtstring) 'dt', ifile, dt_agreed, this%fwd(ifile)%dt
          call pabort
       end if
    end do

    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the forward case")' 


    do ifile = 1, this%nfiles_bwd
       if (dt_agreed.ne.this%bwd(ifile)%dt) then
          write(*,fmtstring) 'dt', ifile, dt_agreed, this%bwd(ifile)%dt
          call pabort
       end if
    end do

    this%dt = dt_agreed
    write(lu_out, 101) 'dt', this%dt
 
    ! Check whether npol is the same in all files
    npol_agreed = this%fwd(1)%npol

    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, " (", I7,") vs ", I7, " in the others")' 
    do ifile = 1, this%nfiles_fwd
       if (npol_agreed.ne.this%fwd(ifile)%npol) then
          write(*,fmtstring) 'npol', ifile, npol_agreed, this%fwd(ifile)%npol
          call pabort
       end if
    end do

    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",I7,"s) vs ", I7, " in the forward case")' 

    do ifile = 1, this%nfiles_bwd
       if (npol_agreed.ne.this%bwd(ifile)%npol) then
          write(*,fmtstring) 'npol', ifile, npol_agreed, this%bwd(ifile)%npol
          call pabort
       end if
    end do

    this%npol = npol_agreed
    write(lu_out, 100) 'npol', this%npol


    ! Check whether the number of dumps (time samples) is the same in all files
    ndumps_agreed = this%fwd(1)%ndumps
    nseis_agreed  = this%fwd(1)%nseis

    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",I7,") vs ", I7, " in the others")' 
    do ifile = 1, this%nfiles_fwd
       if (ndumps_agreed.ne.this%fwd(ifile)%ndumps) then
          write(*,fmtstring) 'ndumps', ifile, ndumps_agreed, this%fwd(ifile)%ndumps
          call pabort
       end if
       if (nseis_agreed.ne.this%fwd(ifile)%nseis) then
          write(*,fmtstring) 'nseis', ifile, nseis_agreed, this%fwd(ifile)%nseis
          call pabort
       end if
    end do

    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",I7,"s) vs ", I7, " in the forward case")' 

    do ifile = 1, this%nfiles_bwd
       if (ndumps_agreed.ne.this%bwd(ifile)%ndumps) then
          write(*,fmtstring) 'ndumps', ifile, ndumps_agreed, this%bwd(ifile)%ndumps
          call pabort
       end if
       if (nseis_agreed.ne.this%bwd(ifile)%nseis) then
          write(*,fmtstring) 'nseis', ifile, nseis_agreed, this%bwd(ifile)%nseis
          call pabort
       end if
    end do

    this%ndumps = ndumps_agreed
    this%windowlength = ndumps_agreed * dt_agreed
    write(lu_out, 100) 'ndumps', this%ndumps
    write(lu_out, 101) 'windowlength', this%windowlength

    ! Check whether the source time shift and stf are the same in all files
    source_shift_agreed_fwd = this%fwd(1)%source_shift_t
    stf_agreed_fwd = this%fwd(1)%stf
    stf_d_agreed_fwd = this%fwd(1)%stf_d
    amplitude_agreed_fwd = this%fwd(1)%amplitude

    fmtstring = '("Inconsistency in forward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the others")' 
    fmtstring_stf = '("Inconsistency in forward simulations: ", A, " is different \'// &
                    '  in simulation ", I1, " vs the others")' 
    !do ifile = 1, this%nfiles_fwd
     !  if (source_shift_agreed_fwd.ne.this%fwd(ifile)%source_shift_t) then
      !    write(*,fmtstring) 'source time shift', ifile, source_shift_agreed_fwd, &
       !                      this%fwd(ifile)%source_shift_t
        !  call pabort
      ! end if
      ! if (any(abs(stf_agreed_fwd - this%fwd(ifile)%stf).gt.1e-10)) then
       !    write(*,fmtstring) 'stf', ifile
        !   call pabort
      ! end if
      ! if (any(abs(stf_d_agreed_fwd - this%fwd(ifile)%stf_d).gt.1e-10)) then
       !    write(*,fmtstring) 'stf_d', ifile
        !   call pabort
       !end if
       !if (amplitude_agreed_fwd.ne.this%fwd(ifile)%amplitude) then
        !  write(*,fmtstring) 'source amplitude', ifile, amplitude_agreed_fwd, &
                             !this%fwd(ifile)%amplitude
         ! call pabort
      ! end if
   ! end do

    this%timeshift_fwd = real(source_shift_agreed_fwd, kind=dp)
    allocate(this%stf_fwd(ndumps_agreed))
    allocate(this%stf_d_fwd(ndumps_agreed))
    this%stf_fwd = real(stf_agreed_fwd, kind=dp)
    this%stf_d_fwd = real(stf_d_agreed_fwd, kind=dp)
    this%amplitude_fwd = real(amplitude_agreed_fwd, kind=dp)


    source_shift_agreed_bwd = this%bwd(1)%source_shift_t
    stf_agreed_bwd = this%bwd(1)%stf
    stf_d_agreed_bwd = this%bwd(1)%stf_d 
    amplitude_agreed_bwd = this%bwd(1)%amplitude
    fmtstring = '("Inconsistency in backward simulations: ", A, " is different \'// &
                '  in simulation ", I1, "(",F9.4,"s) vs ", F9.4, " in the others")' 
    fmtstring_stf = '("Inconsistency in backward simulations: ", A, " is different \'// &
                    '  in simulation ", I1, " vs the others")' 

    do ifile = 1, this%nfiles_bwd
       if (source_shift_agreed_bwd.ne.this%bwd(ifile)%source_shift_t) then
          write(*,fmtstring) 'source time shift', ifile, source_shift_agreed_bwd, &
                             this%bwd(ifile)%source_shift_t
          call pabort
       end if
       if (any(abs(stf_agreed_bwd - this%bwd(ifile)%stf).gt.1e-10)) then
           write(*,fmtstring) 'stf', ifile
           call pabort
       end if
       if (any(abs(stf_d_agreed_bwd - this%bwd(ifile)%stf_d).gt.1e-10)) then
           write(*,fmtstring) 'stf_d', ifile
           call pabort
       end if
       if (amplitude_agreed_bwd.ne.this%bwd(ifile)%amplitude) then
          write(*,fmtstring) 'source amplitude', ifile, amplitude_agreed_bwd, &
                             this%bwd(ifile)%amplitude
          call pabort
       end if
    end do

    this%timeshift_bwd = real(source_shift_agreed_bwd, kind=dp)
    allocate(this%stf_bwd(ndumps_agreed))
    allocate(this%stf_d_bwd(ndumps_agreed))
    this%stf_bwd = real(stf_agreed_bwd, kind=dp)
    this%stf_d_bwd = real(stf_d_agreed_bwd, kind=dp)
    this%amplitude_bwd = real(amplitude_agreed_bwd, kind=dp)

    this%dt = dt_agreed
    this%decimate_factor = nseis_agreed / ndumps_agreed
    this%nseis  = ndumps_agreed * this%decimate_factor        

    write(lu_out, '(A)') ' all simulation parameters agree and have been set'

    call flush(lu_out)

end subroutine check_consistency
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Loads the model coefficients for a selected coordinate 
function load_model_coeffs(this, coordinates_xyz) result(model)
   use background_model, only         : backgroundmodel_type
   use kdtree2_module, only           : kdtree2_result, kdtree2_n_nearest

   class(semdata_type)               :: this
   real(kind=dp), intent(in)         :: coordinates_xyz(:,:)
   type(backgroundmodel_type)        :: model

   real(kind=dp)                     :: coordinates_r(size(coordinates_xyz,2))
   real(kind=sp)                     :: coeffs(nmodel_parameters_sem_file, size(coordinates_xyz,2)) 

   if (.not.this%kdtree_built) then
      print *, 'ERROR: KDTree is not built yet. Call build_kdtree before loading points!'
      call pabort()
   end if

   coordinates_r = norm2(coordinates_xyz, dim=1)

   coeffs = get_model_coeffs(this, coordinates_r)

   ! Combine 6 mesh values to get the 12 parameters of backgroundmodel.f90
   call model%combine(coeffs)

end function load_model_coeffs
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Gets the model coefficients for a selected point
function get_model_coeffs(this, r) result(coeffs)
   class(semdata_type), intent(in) :: this
   real(kind=dp), intent(in)       :: r(:)
   real(kind=sp)                   :: coeffs(nmodel_parameters_sem_file, size(r))
   
   ! Load model coefficients vp, vs and rho at point ipoint
   ! Load coefficient vp
   coeffs(1,:) = this%fwdmesh%vp%get(r)
   ! Load coefficient vs
   coeffs(2,:) = this%fwdmesh%vs%get(r)
   ! Load coefficient rho
   coeffs(3,:) = this%fwdmesh%rho%get(r)
   ! Load coefficient phi
   coeffs(4,:) = this%fwdmesh%phi%get(r)
   ! Load coefficient xi
   coeffs(5,:) = this%fwdmesh%xi%get(r)
   ! Load coefficient eta
   coeffs(6,:) = this%fwdmesh%eta%get(r)

end function get_model_coeffs
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine load_seismogram(this, rec_in, src_in)
!< This function loads a seismogram via the reciprocity database mode
   class(semdata_type)               :: this

   type(src_param_type)              :: src_in
   type(rec_param_type)              :: rec_in(:)

   type(src_param_type), allocatable :: src_rdbm(:)
   type(receivers_rdbm_type)         :: rec_rdbm

   integer                           :: nrec
   integer                           :: irec
 
   real(kind=dp), allocatable        :: seismogram(:,:,:)

   if (.not.this%kdtree_built) then
      print *, 'ERROR: KDTree is not built yet. Call build_kdtree before loading points!'
      call pabort()
   end if

   nrec=size(rec_in)

   ! just need 1 source, this is just a hook
   allocate(src_rdbm(1))
   src_rdbm(1) = src_in

   call rec_rdbm%create_reci_sources(rec_in)

   allocate(this%seis(this%ndumps, nrec))

   allocate(seismogram(this%ndumps, 1, 1)) ! last 1 means 1 source

   do irec=1,nrec

      seismogram = this%load_fw_points_rdbm(src_rdbm, rec_rdbm%reci_sources(irec), &
                                            rec_in(irec)%component)

      this%seis(:, irec) = seismogram(:,1,1)  
   
   end do

   deallocate(seismogram)


end subroutine load_seismogram
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine find_element(mesh, tree, s, z, xi, eta, corner_points, id_elem, eltype, &
                        gll_point_ids, axis)
  use kdtree2_module, only           : kdtree2_result, kdtree2_n_nearest, kdtree2
  use finite_elem_mapping, only      : inside_element
  type(meshtype)                    :: mesh
  type(kdtree2), pointer            :: tree
  real(kind=dp), intent(inout)      :: s, z
  real(kind=dp), intent(out)        :: xi, eta, corner_points(4,2)
  integer,       intent(out)        :: id_elem, eltype, gll_point_ids(0:, 0:)
  logical,       intent(out)        :: axis

  integer                           :: inext_point, icp
  integer                           :: corner_point_ids(4)
  integer, parameter                :: nnext_points = 12
  type(kdtree2_result), allocatable :: nextpoint(:)

  ! Find the twelve closest midpoints first
  allocate(nextpoint(nnext_points))
  call kdtree2_n_nearest( tree,                   &
                          real([s, z], kind=sp),  &
                          nn = nnext_points,      &
                          results = nextpoint )

  ! Check, whether point is in any of the twelve closest elements
  do inext_point = 1, nnext_points
      ! get cornerpoints of finite element
      corner_point_ids = mesh%corner_point_ids(:, nextpoint(inext_point)%idx)
      eltype = mesh%eltype(nextpoint(inext_point)%idx)
      
      do icp = 1, 4
          corner_points(icp, 1) = mesh%s(corner_point_ids(icp)+1)
          corner_points(icp, 2) = mesh%z(corner_point_ids(icp)+1)
      enddo                        

      ! test point to be inside, if so, exit
      if (inside_element(s, z, &
                         corner_points, eltype, xi=xi, eta=eta, &
                         tolerance=1d-3)) then
          if (verbose > 1) then
             write(6,*) 's, z       = ', s, z
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
     write(6,*) 'radius     = ', norm2([s, z])
     !this%bwd(1)%count_error_pointoutside = this%bwd(1)%count_error_pointoutside + 1
     !cycle
     call pabort(do_traceback = .false.)
  endif

  id_elem = nextpoint(inext_point)%idx

  ! get gll points of spectral element
  gll_point_ids = -1
  if (verbose > 1) &
      write(6,*) 'element id = ', nextpoint(inext_point)%idx
  
  ! gll_point_ids starts at 0 in NetCDF file
  gll_point_ids = mesh%gll_point_ids(:,:,id_elem) + 1
  if (verbose > 1) &
      write(6,*) 'gll_point_ids = ', gll_point_ids(:,0)


  if (mesh%isaxis(id_elem) == 1) then
     axis = .true.
  elseif (mesh%isaxis(id_elem) == 0) then
     axis = .false.
  else
     call pabort
  endif

  if (verbose > 1) &
     write(6,*) 'axis = ', axis

end subroutine find_element
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_fw_points(this, coordinates, source_params, model)
    use finite_elem_mapping, only      : inside_element
    use background_model, only         : backgroundmodel_type
    use simple_routines, only          : check_NaN
    use kdtree2_module, only           : kdtree2_result, kdtree2_n_nearest
    use rotations, only                : azim_factor_nsim

    class(semdata_type)               :: this
    real(kind=dp), intent(in)         :: coordinates(:,:)
    type(src_param_type), intent(in)  :: source_params
    real(kind=dp)                     :: load_fw_points(this%ndumps, this%ndim, &
                                                        size(coordinates,2))

    type(backgroundmodel_type), intent(out), optional :: model

    type(kdtree2_result), allocatable :: nextpoint(:)
    integer                           :: npoints, nnext_points
    integer                           :: ipoint, isim
    integer(kind=long)                :: iclockold
    integer                           :: eltype(size(coordinates,2))
    logical                           :: axis(size(coordinates,2))
    integer, allocatable              :: gll_point_ids(:,:,:)
    integer                           :: id_elem(size(coordinates,2))
    real(kind=dp)                     :: corner_points(4,2,size(coordinates,2))
    real(kind=dp)                     :: rotmesh_s(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_phi(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_z(size(coordinates,2))
    real(kind=dp)                     :: utemp(this%ndumps, this%ndim, &
                                               this%nsim_fwd, size(coordinates,2))
    real(kind=sp), allocatable        :: coeffs(:,:)
    real(kind=dp)                     :: xi(size(coordinates,2)), eta(size(coordinates,2))
    real(kind=dp)                     :: az_1(4, size(coordinates,2))
    real(kind=dp)                     :: az_2(4, size(coordinates,2))
    real(kind=dp)                     :: dist_from_core_square


    if (.not.this%kdtree_built) then
       print *, 'ERROR: KDTree is not built yet. Call build_kdtree before loading points!'
       call pabort()
    end if
    
    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_fw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       call pabort 
    end if
    npoints = size(coordinates,2)
    
    nnext_points = 6 ! 6, because this is the maximum valence in the mesh
    allocate(gll_point_ids(0:this%npol, 0:this%npol, npoints))

    allocate(coeffs(6,npoints))

    ! Rotate points to FWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,   &
                          coordinates,                                  &
                          source_params%lon, source_params%colat)

    if (present(model)) then
      coeffs = get_model_coeffs(this, norm2(coordinates, dim=1))
    end if

    allocate(nextpoint(nnext_points))
    load_fw_points(:,:,:) = 0.0
    iclockold = tick()
    do ipoint = 1, npoints
        ! map points from outside earth to the surface:
        dist_from_core_square = rotmesh_s(ipoint)**2 + rotmesh_z(ipoint)**2
        if (dist_from_core_square > this%fwd(1)%planet_radius**2) then
           rotmesh_s(ipoint) = rotmesh_s(ipoint) &
                               / sqrt(dist_from_core_square) * this%fwd(1)%planet_radius
           rotmesh_z(ipoint) = rotmesh_z(ipoint) &
                               / sqrt(dist_from_core_square) * this%fwd(1)%planet_radius
        endif

        ! Find the twelve closest midpoints first
        call find_element(mesh    = this%fwdmesh,                    &
                          tree    = this%fwdtree_mp,                 &
                          s       = rotmesh_s(ipoint),               &
                          z       = rotmesh_z(ipoint),               &
                          xi      = xi(ipoint),                      &
                          eta     = eta(ipoint),                     &
                          corner_points = corner_points(:,:,ipoint), &
                          gll_point_ids = gll_point_ids(:,:,ipoint), &
                          id_elem = id_elem(ipoint),                 &
                          eltype  = eltype(ipoint),                  &
                          axis    = axis(ipoint)  )

    end do
    iclockold = tick(id=id_find_point_fwd, since=iclockold)

    az_1 = azim_factor_nsim(rotmesh_phi, source_params%mij, 1) 
    az_2 = azim_factor_nsim(rotmesh_phi, source_params%mij, 2) 

    if (this%merged_fwd) then
      call load_strain_point_merged_fwd( &
                              sem_obj       = this%fwd(1),               &
                              xi            = xi,                &
                              eta           = eta,               &
                              strain_type   = this%strain_type,          &
                              nodes         = corner_points, &
                              element_type  = eltype,            &
                              axis          = axis,              &
                              id_elem       = id_elem,           &
                              strain_result = utemp)
       utemp = utemp / this%fwd(1)%amplitude
    else
      do ipoint = 1, npoints
        do isim = 1, this%nsim_fwd
          utemp(:,:,isim,ipoint) = load_strain_point_interp( &
                                     this%fwd(isim), &
                                     gll_point_ids(:,:,ipoint),  &
                                     xi(ipoint), eta(ipoint), this%strain_type,      &
                                     corner_points(:,:,ipoint), eltype(ipoint), axis(ipoint), &
                                     id_elem = id_elem(ipoint))              &
                             / this%fwd(isim)%amplitude
        end do 
      end do
    end if

    do ipoint = 1, npoints

        select case(trim(this%strain_type))
        case('straintensor_trace')    

          ! Set NaNs to zero
          where(utemp.ne.utemp) utemp = 0.0
          
          iclockold = tick()

          do isim = 1, this%nsim_fwd
            load_fw_points(:, :, ipoint) = load_fw_points(:,:,ipoint)          &
                 + utemp(:, :, isim, ipoint) * az_1(isim, ipoint)
          end do
          iclockold = tick(id=id_rotate, since=iclockold)

        case('straintensor_full')

           iclockold = tick()

           do isim = 1, this%nsim_fwd  
             load_fw_points(:,1,ipoint) = load_fw_points(:,1,ipoint) &
                   + utemp(:,1,isim,ipoint) * az_1(isim, ipoint)               
             load_fw_points(:,2,ipoint) = load_fw_points(:,2,ipoint) &
                   + utemp(:,2,isim,ipoint) * az_1(isim, ipoint)
             load_fw_points(:,3,ipoint) = load_fw_points(:,3,ipoint) &
                   + utemp(:,3,isim,ipoint) * az_1(isim, ipoint)
             load_fw_points(:,4,ipoint) = load_fw_points(:,4,ipoint) &
                   + utemp(:,4,isim,ipoint) * az_2(isim, ipoint)
             load_fw_points(:,5,ipoint) = load_fw_points(:,5,ipoint) &
                   + utemp(:,5,isim,ipoint) * az_1(isim, ipoint)
             load_fw_points(:,6,ipoint) = load_fw_points(:,6,ipoint) &
                   + utemp(:,6,isim,ipoint) * az_2(isim, ipoint)
           end do 
           iclockold = tick(id=id_rotate, since=iclockold)

           load_fw_points(:,:,ipoint) = rotate_symm_tensor_voigt_src_to_xyz( &
                                          load_fw_points(:,:,ipoint),        &
                                          source_params%lon, this%ndumps    )

           load_fw_points(:,:,ipoint) = rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(        &
                                          load_fw_points(:,:,ipoint),                         &
                                          source_params%lon, source_params%colat, this%ndumps)

        end select

    end do !ipoint

    if (present(model)) call model%combine(coeffs)

end function load_fw_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_bw_points(this, coordinates, receiver)
    use finite_elem_mapping, only      : inside_element
    use simple_routines, only          : check_NaN
    use kdtree2_module, only           : kdtree2_result, kdtree2_n_nearest

    class(semdata_type)               :: this
    real(kind=dp), intent(in)         :: coordinates(:,:)
    type(rec_param_type)              :: receiver
    real(kind=dp)                     :: load_bw_points(this%ndumps, this%ndim, &
                                                         size(coordinates,2))

    type(kdtree2_result), allocatable :: nextpoint(:)
    integer                           :: npoints, nnext_points, isim
    integer                           :: ipoint
    integer(kind=long)                :: iclockold
    integer                           :: eltype(size(coordinates,2))
    logical                           :: axis(size(coordinates,2))
    integer, allocatable              :: gll_point_ids(:,:,:)
    integer                           :: id_elem(size(coordinates,2))
    real(kind=dp)                     :: corner_points(4,2,size(coordinates,2))
    real(kind=dp)                     :: rotmesh_s(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_phi(size(coordinates,2))
    real(kind=dp)                     :: rotmesh_z(size(coordinates,2))
    real(kind=dp)                     :: utemp_nsim(this%ndumps, this%ndim, &
                                                    this%nsim_bwd, size(coordinates,2))
    real(kind=dp)                     :: utemp(this%ndumps, this%ndim, size(coordinates,2))
    real(kind=dp)                     :: xi(size(coordinates,2)), eta(size(coordinates,2))
    real(kind=dp)                     :: az_1, az_2, az(6, size(coordinates,2))
    real(kind=dp)                     :: dist_from_core_square

    
    if (.not.this%kdtree_built) then
       print *, 'ERROR: KDTree is not built yet. Call build_kdtree before loading points!'
       call pabort()
    end if

    if (size(coordinates,1).ne.3) then
       write(*,*) ' Error in load_bw_points: input variable coordinates has to be a '
       write(*,*) ' 3 x npoints array'
       call pabort
    end if
    npoints = size(coordinates,2)

    nnext_points = 6 ! 6, because this is the maximum valence in the mesh
    allocate(gll_point_ids(0:this%npol, 0:this%npol, size(coordinates,2)))

    ! Rotate points to BWD coordinate system
    call rotate_frame_rd( npoints, rotmesh_s, rotmesh_phi, rotmesh_z,   &
                          coordinates,                                  &
                          receiver%lon, receiver%colat)

    allocate(nextpoint(nnext_points))
    load_bw_points(:,:,:) = 0.0
    do ipoint = 1, npoints
        iclockold = tick()
        ! map points from outside earth to the surface:
        dist_from_core_square = rotmesh_s(ipoint)**2 + rotmesh_z(ipoint)**2
        if (dist_from_core_square > this%fwd(1)%planet_radius**2) then
           rotmesh_s(ipoint) = rotmesh_s(ipoint) &
                               / sqrt(dist_from_core_square) * this%fwd(1)%planet_radius
           rotmesh_z(ipoint) = rotmesh_z(ipoint) &
                               / sqrt(dist_from_core_square) * this%fwd(1)%planet_radius
        endif

        call find_element(mesh    = this%bwdmesh, &
                          tree    = this%bwdtree_mp, &
                          s       = rotmesh_s(ipoint),               &
                          z       = rotmesh_z(ipoint),               &
                          xi      = xi(ipoint),                      &
                          eta     = eta(ipoint),                     &
                          corner_points = corner_points(:,:,ipoint), &
                          gll_point_ids = gll_point_ids(:,:,ipoint), &
                          id_elem = id_elem(ipoint),                 &
                          eltype  = eltype(ipoint),                  &
                          axis    = axis(ipoint)  )
        iclockold = tick(id=id_find_point_bwd, since=iclockold)

        select case(receiver%component)
        case('Z')
          isim = 1
          az   = 1

        case('R') 
          isim = 2
          az_1 = azim_factor_bw(rotmesh_phi(ipoint), [0d0, 1d0, 0d0], isim, 1)
          az_2 = azim_factor_bw(rotmesh_phi(ipoint), [0d0, 1d0, 0d0], isim, 2)
          az(1, ipoint) = az_1
          az(2, ipoint) = az_1
          az(3, ipoint) = az_1
          az(4, ipoint) = az_2
          az(5, ipoint) = az_1 * 2
          az(6, ipoint) = az_2 * 2

        case('T') 
          isim = 2
          az_1 = azim_factor_bw(rotmesh_phi(ipoint), [0d0, 0d0, 1d0], isim, 1)
          az_2 = azim_factor_bw(rotmesh_phi(ipoint), [0d0, 0d0, 1d0], isim, 2)
          az(1, ipoint) = az_1
          az(2, ipoint) = az_1
          az(3, ipoint) = az_1
          az(4, ipoint) = az_2
          az(5, ipoint) = az_1 * 2
          az(6, ipoint) = az_2 * 2

        end select

    end do


    if (this%merged_bwd) then
      call load_strain_point_merged_bwd(this%bwd(1),                &
                                        xi, eta, this%strain_type,  &
                                        corner_points, eltype,      &
                                        axis, id_elem,              &
                                        strain_result = utemp_nsim)              
      utemp = utemp_nsim(:,:,isim,:) 
    else
      do ipoint = 1, npoints
        utemp(:,:,ipoint) = load_strain_point_interp(this%bwd(isim), &
                                         gll_point_ids(:,:,ipoint),  &
                                         xi(ipoint), eta(ipoint), this%strain_type,      &
                                         corner_points(:,:,ipoint), eltype(ipoint), &
                                         axis(ipoint), id_elem(ipoint))              
      end do
    end if


    do ipoint = 1, npoints

        load_bw_points(:, 1, ipoint) = utemp(:,1, ipoint) * az(1, ipoint)

        if (this%strain_type.eq.'straintensor_full') then
          load_bw_points(:, 2, ipoint) = utemp(:,2,ipoint) * az(2, ipoint)
          load_bw_points(:, 3, ipoint) = utemp(:,3,ipoint) * az(3, ipoint)
          load_bw_points(:, 4, ipoint) = utemp(:,4,ipoint) * az(4, ipoint)
          load_bw_points(:, 5, ipoint) = utemp(:,5,ipoint) * az(5, ipoint)
          load_bw_points(:, 6, ipoint) = utemp(:,6,ipoint) * az(6, ipoint)
        end if
        
        load_bw_points(:,:,ipoint) = load_bw_points(:,:,ipoint) &
                                     / this%bwd(1)%amplitude

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
function load_fw_points_rdbm(this, source_params, reci_source_params, component)
    use finite_elem_mapping, only       : inside_element
    use kdtree2_module, only            : kdtree2_result, kdtree2_n_nearest

    class(semdata_type)                     :: this
    type(src_param_type), intent(in)        :: source_params(:)
    type(src_param_type), intent(in)        :: reci_source_params
    character(len=1), intent(in)            :: component
    real(kind=dp), allocatable              :: load_fw_points_rdbm(:,:,:)

    type(kdtree2_result), allocatable :: nextpoint(:)
    integer                           :: npoints, nnext_points, id_elem
    integer                           :: ipoint, isim, i
    integer                           :: eltype(1)
    logical                           :: axis
    integer, allocatable              :: gll_point_ids(:,:)
    real(kind=dp)                     :: corner_points(4,2,1)
    real(kind=dp)                     :: rotmesh_s(size(source_params))
    real(kind=dp)                     :: rotmesh_phi(size(source_params))
    real(kind=dp)                     :: rotmesh_z(size(source_params))
    real(kind=dp)                     :: utemp_nsim(this%ndumps, 6, this%nsim_bwd, 1)
    real(kind=dp)                     :: utemp(this%ndumps, 6)
    real(kind=dp)                     :: coordinates(3,size(source_params))
    real(kind=dp)                     :: mij_buff(6)
    real(kind=dp)                     :: xi, eta

    if (.not.this%kdtree_built) then
       print *, 'ERROR: KDTree is not built yet. Call build_kdtree before loading points!'
       call pabort()
    end if

    nnext_points = 6 ! 6, because this is the maximum valence in the mesh
    allocate(gll_point_ids(0:this%npol, 0:this%npol))

    allocate(load_fw_points_rdbm(this%ndumps, 1, size(source_params)))
    load_fw_points_rdbm(:,:,:) = 0.0
   
    coordinates(:, 1) = source_params(1)%get_r(this%fwd(1)%planet_radius)
    
    ! Rotate points to FWD coordinate system
    call rotate_frame_rd( 1, rotmesh_s, rotmesh_phi, rotmesh_z, coordinates, &
                          reci_source_params%lon, reci_source_params%colat)

    allocate(nextpoint(nnext_points))
    ipoint = 1

    call find_element(mesh    = this%fwdmesh, &
                      tree    = this%fwdtree_mp, &
                      s       = rotmesh_s(ipoint), &
                      z       = rotmesh_z(ipoint),  &
                      xi      = xi,                 &
                      eta     = eta,                &
                      corner_points = corner_points(:,:,1), &
                      gll_point_ids = gll_point_ids, &
                      id_elem = id_elem,             &
                      eltype  = eltype(1),              &
                      axis    = axis  )

    ! rotate source mt to global cartesian system
    mij_buff = rotate_symm_tensor_voigt_xyz_src_to_xyz_earth( &
                     source_params(ipoint)%mij_voigt, &
                     source_params(ipoint)%lon, &
                     source_params(ipoint)%colat)

    ! rotate source mt to receiver cartesian system
    mij_buff = rotate_symm_tensor_voigt_xyz_earth_to_xyz_src( &
                     mij_buff, reci_source_params%lon, reci_source_params%colat)

    ! rotate source mt to receiver s,phi,z system
    mij_buff = rotate_symm_tensor_voigt_xyz_to_src(mij_buff, rotmesh_phi(ipoint))
    mij_buff = mij_buff / this%bwd(1)%amplitude

    select case(component)
    case('Z')
      isim = 1
    case('R')
      isim = 2
    case('T')
      isim = 2
    case default
      write(6,*) 'component "', component, '" unknown or not yet implemented'
      call pabort
    end select

    if (this%merged_bwd) then
      call load_strain_point_merged_bwd(this%bwd(1),                    &
                                        [xi], [eta], 'straintensor_full',   &
                                        corner_points, eltype, [axis], &
                                        [id_elem], use_buffer=.false.,      &
                                        strain_result = utemp_nsim)
      utemp = utemp_nsim(:,:,isim,1)
    else
      utemp = load_strain_point_interp_seismogram(this%bwd(isim), gll_point_ids, &
                                       xi, eta, &
                                       corner_points, eltype(1), axis)
    end if

    load_fw_points_rdbm(:, :, ipoint) = 0

    select case(component)
    case('Z')
         do i = 1, 3
            load_fw_points_rdbm(:, 1, ipoint) = &
                  load_fw_points_rdbm(:, 1, ipoint) + mij_buff(i) * utemp(:,i)
         enddo 

         ! components 4-6 need a factor of two because of voigt mapping
         ! without factor of two in the strain
         i = 5
         load_fw_points_rdbm(:, 1, ipoint) = &
               load_fw_points_rdbm(:, 1, ipoint) + 2 * mij_buff(i) * utemp(:,i)

    case('R')
         load_fw_points_rdbm(:, 1, ipoint) &
              = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(1) * utemp(:,1)
         load_fw_points_rdbm(:, 1, ipoint) &
              = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(2) * utemp(:,2)
         load_fw_points_rdbm(:, 1, ipoint) &
              = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(3) * utemp(:,3) 
         load_fw_points_rdbm(:, 1, ipoint) &
              = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(5) * utemp(:,5) * 2 
         
         !@TODO not sure why we need the - sign here. Might be because N
         !      is in negative theta direction
         load_fw_points_rdbm(:, 1, ipoint) = - load_fw_points_rdbm(:, 1, ipoint)


    case('T')
         load_fw_points_rdbm(:, 1, ipoint) &
              = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(4) * utemp(:,4) * 2
         load_fw_points_rdbm(:, 1, ipoint) &
              = load_fw_points_rdbm(:, 1, ipoint) + mij_buff(6) * utemp(:,6) * 2 
        
    end select

end function load_fw_points_rdbm
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function load_strain_point_interp_seismogram(sem_obj, pointids, xi, eta, nodes, &
                                             element_type, axis) !, id_elem)
    !< Calculates strain in element given by pointids, nodes. 
    !! This routine is specifically for seismogram retrieval and does not use the buffer
    !! Since seismograms are loaded only once at startup, this does not hurt performance
    !! much and since seismometers are usually not collocated, the buffer would not be
    !! of much use anyway.
    !! Strain is then interpolated to the point defined by xi, eta.
    use sem_derivatives
    use spectral_basis, only : lagrange_interpol_2D_td
    use simple_routines, only        : check_limits

    type(ncparamtype), intent(in)   :: sem_obj
    integer,           intent(in)   :: pointids(0:sem_obj%npol, 0:sem_obj%npol) 
                                                    !< ID of GLL/GLI points in element of
                                                    !! interest
    real(kind=dp),     intent(in)   :: xi, eta      !< Coordinates at which to interpolate
                                                    !! strain
    real(kind=dp),     intent(in)   :: nodes(4,2)   !< Coordinates of element corner
                                                    !! points
    integer,           intent(in)   :: element_type !< Element type in the solver
    logical,           intent(in)   :: axis         !< Axis element or not 
    real(kind=dp)                   :: load_strain_point_interp_seismogram(sem_obj%ndumps,6)

    integer                         :: start_chunk, gll_to_read
    integer                         :: idisplvar
    real(kind=sp), allocatable      :: utemp_chunk(:,:,:)
    real(kind=dp)                   :: utemp(1:sem_obj%ndumps, &
                                             0:sem_obj%npol,   &
                                             0:sem_obj%npol,   &
                                             3)
    real(kind=dp)                   :: strain(1:sem_obj%ndumps, &
                                              0:sem_obj%npol,   &
                                              0:sem_obj%npol,   &
                                              6)
    real(kind=dp)                   :: G( 0:sem_obj%npol, 0:sem_obj%npol)
    real(kind=dp)                   :: GT(0:sem_obj%npol, 0:sem_obj%npol)
    real(kind=dp)                   :: col_points_xi(0:sem_obj%npol), col_points_eta(0:sem_obj%npol)
    integer                         :: ipol, jpol, i
    integer(kind=long)              :: iclockold, iclockold_total
    logical                         :: strain_nan

    iclockold_total = tick()

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


      allocate(utemp_chunk(sem_obj%chunk_gll, sem_obj%ndumps, 3))

      ! load displacements from all GLL points
      do ipol = 0, sem_obj%npol
         do jpol = 0, sem_obj%npol

               call get_chunk_bounds(pointid     = pointids(ipol, jpol), &
                                     chunksize   = sem_obj%chunk_gll,    &
                                     npoints     = sem_obj%ngll,         &
                                     start_chunk = start_chunk,          &   
                                     count_chunk = gll_to_read )

               do idisplvar = 1, 3

                   if (sem_obj%displvarid(idisplvar).eq.-1) then
                       utemp(:, ipol, jpol, idisplvar) = 0
                       cycle ! For monopole source which does not have this component.
                   endif

                   iclockold = tick()

                   call nc_getvar( ncid   = sem_obj%snap,                  & 
                                   varid  = sem_obj%displvarid(idisplvar), &
                                   start  = [start_chunk, 1],              &
                                   count  = [gll_to_read, sem_obj%ndumps], &
                                   values = utemp_chunk(1:gll_to_read, :, idisplvar))

                   strain_nan = check_limits(utemp_chunk(1:gll_to_read,:,idisplvar), &
                                             array_name='strain')

                   ! Set NaNs in utemp_chunk to zero
                   where (utemp_chunk.ne.utemp_chunk) utemp_chunk=0.0

                   !print *, 'suceeded'
                   !call flush(6)
                   iclockold = tick(id=id_netcdf, since=iclockold)
                   utemp(:,ipol,jpol, idisplvar) &
                        = utemp_chunk(pointids(ipol,jpol) - start_chunk + 1,:,idisplvar)
               enddo

!               do iread = 0, sem_obj%chunk_gll - 1
!                   status = sem_obj%buffer_disp%put(start_chunk + iread, &
!                                                    utemp_chunk(iread+1,:,:) )
!               end do
!            else
!               utemp(:,ipol,jpol,:) = real(ubuff(:,:), kind=dp)
!            endif
         enddo
      enddo

      iclockold = tick()
!      select case(strain_type)
!
!      case('straintensor_full')
          ! compute full strain tensor
          if (sem_obj%excitation_type == 'monopole') then
              strain = strain_monopole(utemp(:,:,:,1:3:2), G, GT, col_points_xi, &
                                       col_points_eta, sem_obj%npol, sem_obj%ndumps, nodes, &
                                       element_type, axis)

          elseif (sem_obj%excitation_type == 'dipole') then
              strain = strain_dipole(utemp, G, GT, col_points_xi, &
                                     col_points_eta, sem_obj%npol, sem_obj%ndumps, nodes, &
                                     element_type, axis)

          elseif (sem_obj%excitation_type == 'quadpole') then
              strain = strain_quadpole(utemp, G, GT, col_points_xi, &
                                       col_points_eta, sem_obj%npol, sem_obj%ndumps, nodes, &
                                       element_type, axis)
          else
              print *, 'ERROR: unknown excitation_type: ', sem_obj%excitation_type
              call pabort
          endif
          
          iclockold = tick(id=id_calc_strain, since=iclockold)

        do i = 1, 6
            load_strain_point_interp_seismogram(:, i) &
                = lagrange_interpol_2D_td(col_points_xi, col_points_eta, &
                                          real(strain(:,:,:,i), kind=dp), xi, eta)
        enddo

        iclockold = tick(id=id_lagrange, since=iclockold)

        !@TODO for consistency with SOLVER output
        load_strain_point_interp_seismogram(:, 4) = - load_strain_point_interp_seismogram(:, 4) 
        load_strain_point_interp_seismogram(:, 6) = - load_strain_point_interp_seismogram(:, 6) 

    !end select

    iclockold_total = tick(id=id_load_strain, since=iclockold_total)

end function load_strain_point_interp_seismogram
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
    use simple_routines, only        : check_limits

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
    integer                         :: idisplvar
    real(kind=sp), allocatable      :: utemp_chunk(:,:,:)
    real(kind=sp), allocatable      :: ubuff(:,:)
    real(kind=dp)                   :: utemp(1:sem_obj%ndumps, &
                                             0:sem_obj%npol,   &
                                             0:sem_obj%npol,   &
                                             3)
    real(kind=dp)                   :: strain(1:sem_obj%ndumps, &
                                              0:sem_obj%npol,   &
                                              0:sem_obj%npol,   &
                                              6)
    real(kind=dp)                   :: straintrace(1:sem_obj%ndumps, &
                                                   0:sem_obj%npol,   &
                                                   0:sem_obj%npol)
    real(kind=dp)                   :: G( 0:sem_obj%npol, 0:sem_obj%npol)
    real(kind=dp)                   :: GT(0:sem_obj%npol, 0:sem_obj%npol)
    real(kind=dp)                   :: col_points_xi(0:sem_obj%npol)
    real(kind=dp)                   :: col_points_eta(0:sem_obj%npol)
    integer                         :: ipol, jpol, i, status
    integer(kind=long)              :: iclockold, iclockold_total
    logical                         :: strain_nan

    iclockold_total = tick()

    use_strainbuffer = present(id_elem)

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
      utemp_chunk = 0
      allocate(ubuff(sem_obj%ndumps, 3))
      ubuff = 0

      ! load displacements from all GLL points
      do ipol = 0, sem_obj%npol
         do jpol = 0, sem_obj%npol

            iclockold = tick()
            status = sem_obj%buffer_disp%get(pointids(ipol,jpol), ubuff(:,:))
            iclockold = tick(id=id_buffer, since=iclockold)
            if (status.ne.0) then

               call get_chunk_bounds(pointid     = pointids(ipol, jpol), &
                                     chunksize   = sem_obj%chunk_gll,    &
                                     npoints     = sem_obj%ngll,         &
                                     start_chunk = start_chunk,          &
                                     count_chunk = gll_to_read )

               do idisplvar = 1, 3

                   if (sem_obj%displvarid(idisplvar).eq.-1) then
                       utemp(:, ipol, jpol, idisplvar) = 0
                       cycle ! For monopole source which does not have this component.
                   endif

                   iclockold = tick()
                   !print *, 'Trying to read data!'
                   !print *, 'pointid:      ', pointids(ipol, jpol)
                   !print *, 'start_chunk:  ', start_chunk
                   !print *, 'gll_to_read:  ', gll_to_read
                   !print *, 'last_element: ', start_chunk + gll_to_read - 1
                   !call flush(6)

                   call nc_getvar( ncid   = sem_obj%snap,                  & 
                                   varid  = sem_obj%displvarid(idisplvar), &
                                   start  = [start_chunk, 1],              &
                                   count  = [gll_to_read, sem_obj%ndumps], &
                                   values = utemp_chunk(1:gll_to_read, :, idisplvar))

                   strain_nan = check_limits(utemp_chunk(1:gll_to_read,:,idisplvar), &
                                             array_name='strain')

                   ! Set NaNs in utemp_chunk to zero
                   where (utemp_chunk.ne.utemp_chunk) utemp_chunk=0.0

                   !print *, 'suceeded'
                   !call flush(6)
                   iclockold = tick(id=id_netcdf, since=iclockold)
                   utemp(:,ipol,jpol, idisplvar) &
                        = real(utemp_chunk(pointids(ipol,jpol) - start_chunk + 1, &
                                           :,                                     &
                                           idisplvar),                            &
                               kind=dp)
               enddo

               do iread = 0, sem_obj%chunk_gll - 1
                   status = sem_obj%buffer_disp%put(start_chunk + iread, &
                                                    utemp_chunk(iread+1,:,:) )
               end do
            else
               utemp(:,ipol,jpol,:) = real(ubuff(:,:), kind=dp)
            endif
         enddo
      enddo

      iclockold = tick()
      select case(strain_type)
      case('straintensor_trace')
          ! compute straintrace
          if (sem_obj%excitation_type == 'monopole') then
              straintrace = straintrace_monopole(utemp(:,:,:,1:3:2), G, GT, col_points_xi, &
                                                 col_points_eta, sem_obj%npol, &
                                                 sem_obj%ndumps, nodes, element_type, axis)

          elseif (sem_obj%excitation_type == 'dipole') then
              straintrace = straintrace_dipole(utemp, G, GT, col_points_xi, &
                                               col_points_eta, sem_obj%npol, sem_obj%ndumps, &
                                               nodes, element_type, axis)

          elseif (sem_obj%excitation_type == 'quadpole') then
              straintrace = straintrace_quadpole(utemp, G, GT, col_points_xi, &
                                                 col_points_eta, sem_obj%npol, &
                                                 sem_obj%ndumps, nodes, element_type, axis)
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
              strain = strain_monopole(utemp(:,:,:,1:3:2), G, GT, col_points_xi, &
                                       col_points_eta, sem_obj%npol, sem_obj%ndumps, nodes, &
                                       element_type, axis)

          elseif (sem_obj%excitation_type == 'dipole') then
              strain = strain_dipole(utemp, G, GT, col_points_xi, &
                                     col_points_eta, sem_obj%npol, sem_obj%ndumps, nodes, &
                                     element_type, axis)

          elseif (sem_obj%excitation_type == 'quadpole') then
              strain = strain_quadpole(utemp, G, GT, col_points_xi, &
                                       col_points_eta, sem_obj%npol, sem_obj%ndumps, nodes, &
                                       element_type, axis)
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
                                      straintrace, xi, eta)

        iclockold = tick(id=id_lagrange, since=iclockold)

    case('straintensor_full')
        allocate(load_strain_point_interp(sem_obj%ndumps, 6))
        do i = 1, 6
            load_strain_point_interp(:, i) &
                = lagrange_interpol_2D_td(col_points_xi, col_points_eta, &
                                          strain(:,:,:,i), xi, eta)
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
subroutine load_strain_point_merged_fwd(sem_obj, xi, eta, strain_type, nodes, &
                                        element_type, axis, id_elem, use_buffer, strain_result)
    use spectral_basis, only         : lagrange_interpol_2D_td
    use sem_derivatives, only        : straintrace_merged_fwd, strain_merged_fwd
    type(ncparamtype)               :: sem_obj
    real(kind=dp),     intent(in)   :: xi(:), eta(:)      !< Coordinates at which to interpolate
                                                    !! strain
    character(len=*),  intent(in)   :: strain_type  !< Model parameter (decides on 
                                                    !! straintrace (vp) or full tensor (vs)
    real(kind=dp),     intent(in)   :: nodes(:,:,:) !< Coordinates of element corner
                                                    !! points
    integer,           intent(in)   :: element_type(:) !< Element type in the solver
    logical,           intent(in)   :: axis(:)         !< Axis element or not 
    integer,           intent(in)   :: id_elem(:)      !< ID of element to interpolate strain in
    logical, optional, intent(in)   :: use_buffer
    real(kind=dp),     intent(out)  :: strain_result(:,:,:,:)

    integer                         :: ipoint, npoints

    real(kind=dp)                   :: G( 0:sem_obj%npol, 0:sem_obj%npol, size(xi,1))
    real(kind=dp)                   :: GT(0:sem_obj%npol, 0:sem_obj%npol, size(xi,1))
    real(kind=dp)                   :: col_points_xi(0:sem_obj%npol, size(xi,1))
    real(kind=dp)                   :: col_points_eta(0:sem_obj%npol, size(xi,1))
    integer(kind=long)              :: iclockold_total, iclockold
    integer                         :: status_strain(size(xi,1))
    integer                         :: status_disp = -1
    integer                         :: temp
    logical                         :: use_buffer_loc
    logical                         :: to_be_read_from_disk(size(xi,1))

    iclockold_total = tick()
    iclockold = tick()

    npoints = size(xi,1)
    status_strain = -1

    if (present(use_buffer)) then
      use_buffer_loc = use_buffer
    else
      use_buffer_loc = .true.
    end if

    !select case(strain_type)
    !case('straintensor_trace')
    !  straintrace_fwd = 0
    !case('straintensor_full')
    !  strain_fwd = 0            
    !end select
    !utemp_fwd = 0

    iclockold = tick(id=id_allocate, since=iclockold)
    if (use_buffer_loc) then
      select case(strain_type)
      case('straintensor_trace')
        do ipoint = 1, npoints
          status_strain(ipoint) = sem_obj%buffer_strain%get(id_elem(ipoint), &
                                                            straintrace_fwd(:,:,:,:,ipoint))
        end do
      case('straintensor_full')
        do ipoint = 1, npoints
          status_strain(ipoint) = sem_obj%buffer_strain%get(id_elem(ipoint), &
                                                            strain_fwd(:,:,:,:,:,ipoint))
        end do
      end select
    else
      status_strain = - 1
    end if

    to_be_read_from_disk = .true.
    do ipoint = 1, npoints
      ! Try displacement buffer
      iclockold = tick()
      status_disp = sem_obj%buffer_disp%get(id_elem(ipoint), &
                                            utemp_fwd(:,:,:,:,ipoint))
      ! If status was zero (successful), then this does not have to be read
      to_be_read_from_disk(ipoint) = ((status_disp.ne.0) .and.      & 
                                      (status_strain(ipoint).ne.0))
    end do

    iclockold = tick(id=id_buffer, since=iclockold)

    ! If there's anything to be read from disk, do it
    if (any(to_be_read_from_disk)) then
      call read_points_from_disk(sem_obj, id_elem, &
                                 to_be_read_from_disk, utemp_fwd)
    end if
      
    iclockold = tick()

    !if (any(utemp_fwd.ne.utemp_fwd)) print *, 'Found NaN in utemp_fwd!'

    set_element_matrices: do ipoint = 1, npoints
      G(:,:,ipoint)  = sem_obj%G2
      col_points_eta(:, ipoint) = sem_obj%gll_points
      if (axis(ipoint)) then
        GT(:,:,ipoint) = sem_obj%G1T
        col_points_xi(:, ipoint)  = sem_obj%glj_points
      else
        GT(:,:,ipoint) = sem_obj%G2T
        col_points_xi(:, ipoint)  = sem_obj%gll_points
      endif 
    end do set_element_matrices

    iclockold = tick(id=id_allocate_2, since=iclockold)

    select case(strain_type)
    case('straintensor_trace')
      ! Compute just trace of strain (for vp or lambda kernels)
      do concurrent (ipoint = 1: npoints, status_strain(ipoint).eq.-1)
          straintrace_fwd(:,:,:,:,ipoint) = straintrace_merged_fwd(                 &
                                           u = real(utemp_fwd(:,:,:,:,ipoint), kind=dp), &
                                           G = G(:,:,ipoint),                   &
                                           GT = GT(:,:,ipoint),                 &
                                           xi = col_points_xi(:,ipoint),        &          
                                           eta = col_points_eta(:,ipoint),      &
                                           npol = sem_obj%npol,                 &
                                           nsamp = sem_obj%ndumps,              &
                                           nodes = nodes(:,:,ipoint),           &
                                           element_type = element_type(ipoint), &
                                           axial = axis(ipoint))
      end do

    case('straintensor_full')
      ! compute full strain tensor
      do concurrent (ipoint = 1: npoints, status_strain(ipoint).eq.-1)
        strain_fwd(:,:,:,:,:,ipoint) = strain_merged_fwd(                   &
                               u = real(utemp_fwd(:,:,:,:,ipoint), kind=dp),  &
                               G = G(:,:,ipoint),                   &
                               GT = GT(:,:,ipoint),                 &
                               xi = col_points_xi(:,ipoint),        &          
                               eta = col_points_eta(:,ipoint),      &
                               npol = sem_obj%npol,                 &
                               nsamp = sem_obj%ndumps,              &
                               nodes = nodes(:,:,ipoint),           &
                               element_type = element_type(ipoint), &
                               axial = axis(ipoint))
                             end do
    end select

      !select case(strain_type)
      !case('straintensor_trace')
      !  if (any(straintrace.ne.straintrace)) print *, 'Found NaN in straintrace, ipoint:', ipoint, status_strain(ipoint)
      !case('straintensor_full')
      !  if (any(strain.ne.strain)) print *, 'Found NaN in strain, ipoint:', ipoint, status_strain(ipoint)
      !end select
    iclockold = tick(id=id_calc_strain, since=iclockold)

    select case(strain_type)
    case('straintensor_trace')
      do concurrent (ipoint = 1:npoints)
          strain_result(:, 1, :, ipoint)                      &
              = lagrange_interpol_2D_td(col_points_xi(:,ipoint),            &
                                        col_points_eta(:,ipoint),           &
                                        straintrace_fwd(:, :, :, :, ipoint), &
                                        xi(ipoint),                         &
                                        eta(ipoint))
      end do 
    case('straintensor_full')
      do concurrent (ipoint = 1:npoints)
        strain_result(:, :, :, ipoint)                    &
            = lagrange_interpol_2D_td(col_points_xi(:,ipoint),          &
                                      col_points_eta(:,ipoint),         &
                                      strain_fwd(:, :, :, :, :, ipoint), &
                                      xi(ipoint),                       &
                                      eta(ipoint))
      enddo
    end select
    iclockold = tick(id=id_lagrange, since=iclockold)



    if (trim(strain_type).eq.'straintensor_full') then
      !@TODO for consistency with SOLVER output
      strain_result(:, 4, :, :) = - strain_result(:, 4, :, :) 
      strain_result(:, 6, :, :) = - strain_result(:, 6, :, :) 
    end if

    iclockold = tick()
    if (use_buffer_loc) then
      select case(strain_type)
      case('straintensor_trace')
        do ipoint = 1, npoints
          temp = sem_obj%buffer_strain%put(id_elem(ipoint), &
                                           straintrace_fwd(:,:,:,:,ipoint))
        end do
      case('straintensor_full')
        do ipoint = 1, npoints
          temp = sem_obj%buffer_strain%put(id_elem(ipoint), &
                                           strain_fwd(:,:,:,:,:,ipoint))
        end do
      end select
    end if
    iclockold = tick(id=id_buffer, since=iclockold)
    iclockold_total = tick(id=id_load_strain, since=iclockold_total)


end subroutine load_strain_point_merged_fwd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine load_strain_point_merged_bwd(sem_obj, xi, eta, strain_type, nodes, &
                                        element_type, axis, id_elem, use_buffer, strain_result)
    use spectral_basis, only         : lagrange_interpol_2D_td
    use sem_derivatives, only        : straintrace_merged_bwd, strain_merged_bwd
    type(ncparamtype)               :: sem_obj
    real(kind=dp),     intent(in)   :: xi(:), eta(:)      !< Coordinates at which to interpolate
                                                    !! strain
    character(len=*),  intent(in)   :: strain_type  !< Model parameter (decides on 
                                                    !! straintrace (vp) or full tensor (vs)
    real(kind=dp),     intent(in)   :: nodes(:,:,:) !< Coordinates of element corner
                                                    !! points
    integer,           intent(in)   :: element_type(:) !< Element type in the solver
    logical,           intent(in)   :: axis(:)         !< Axis element or not 
    integer,           intent(in)   :: id_elem(:)      !< ID of element to interpolate strain in
    logical, optional, intent(in)   :: use_buffer
    real(kind=dp),     intent(out)  :: strain_result(:,:,:,:)

    integer                         :: ipoint, npoints

    real(kind=dp)                   :: G( 0:sem_obj%npol, 0:sem_obj%npol, size(xi,1))
    real(kind=dp)                   :: GT(0:sem_obj%npol, 0:sem_obj%npol, size(xi,1))
    real(kind=dp)                   :: col_points_xi(0:sem_obj%npol, size(xi,1))
    real(kind=dp)                   :: col_points_eta(0:sem_obj%npol, size(xi,1))
    integer(kind=long)              :: iclockold_total, iclockold
    integer                         :: status_strain(size(xi,1))
    integer                         :: status_disp = -1
    integer                         :: temp
    logical                         :: use_buffer_loc
    logical                         :: to_be_read_from_disk(size(xi,1))

    iclockold_total = tick()
    iclockold = tick()

    npoints = size(xi,1)
    status_strain = -1

    if (present(use_buffer)) then
      use_buffer_loc = use_buffer
    else
      use_buffer_loc = .true.
    end if

    !select case(strain_type)
    !case('straintensor_trace')
    !  straintrace_bwd = 0
    !case('straintensor_full')
    !  strain_bwd = 0            
    !end select
    !utemp_fwd = 0

    iclockold = tick(id=id_allocate, since=iclockold)
    if (use_buffer_loc) then
      select case(strain_type)
      case('straintensor_trace')
        do ipoint = 1, npoints
          status_strain(ipoint) = sem_obj%buffer_strain%get(id_elem(ipoint), &
                                                            straintrace_bwd(:,:,:,:,ipoint))
        end do
      case('straintensor_full')
        do ipoint = 1, npoints
          status_strain(ipoint) = sem_obj%buffer_strain%get(id_elem(ipoint), &
                                                            strain_bwd(:,:,:,:,:,ipoint))
        end do
      end select
    else
      status_strain = - 1
    end if

    to_be_read_from_disk = .true.
    do ipoint = 1, npoints
      ! Try displacement buffer
      iclockold = tick()
      status_disp = sem_obj%buffer_disp%get(id_elem(ipoint), &
                                            utemp_bwd(:,:,:,:,ipoint))
      ! If status was zero (successful), then this does not have to be read
      to_be_read_from_disk(ipoint) = ((status_disp.ne.0) .and.      & 
                                      (status_strain(ipoint).ne.0))
    end do

    iclockold = tick(id=id_buffer, since=iclockold)

    ! If there's anything to be read from disk, do it
    if (any(to_be_read_from_disk)) then
      call read_points_from_disk(sem_obj, id_elem, &
                                 to_be_read_from_disk, utemp_bwd)
    end if
      
    iclockold = tick()

    set_element_matrices: do ipoint = 1, npoints
      G(:,:,ipoint)  = sem_obj%G2
      col_points_eta(:, ipoint) = sem_obj%gll_points
      if (axis(ipoint)) then
        GT(:,:,ipoint) = sem_obj%G1T
        col_points_xi(:, ipoint)  = sem_obj%glj_points
      else
        GT(:,:,ipoint) = sem_obj%G2T
        col_points_xi(:, ipoint)  = sem_obj%gll_points
      endif 
    end do set_element_matrices

    iclockold = tick(id=id_allocate_2, since=iclockold)

    select case(strain_type)
    case('straintensor_trace')
      ! Compute just trace of strain (for vp or lambda kernels)
      do concurrent (ipoint = 1: npoints, status_strain(ipoint).eq.-1)
          straintrace_bwd(:,:,:,:,ipoint) = straintrace_merged_bwd(                 &
                                           u = real(utemp_bwd(:,:,:,:,ipoint), kind=dp), &
                                           G = G(:,:,ipoint),                   &
                                           GT = GT(:,:,ipoint),                 &
                                           xi = col_points_xi(:,ipoint),        &          
                                           eta = col_points_eta(:,ipoint),      &
                                           npol = sem_obj%npol,                 &
                                           nsamp = sem_obj%ndumps,              &
                                           nodes = nodes(:,:,ipoint),           &
                                           element_type = element_type(ipoint), &
                                           axial = axis(ipoint))
      end do

    case('straintensor_full')
      ! compute full strain tensor
      do concurrent (ipoint = 1: npoints, status_strain(ipoint).eq.-1)
        strain_bwd(:,:,:,:,:,ipoint) = strain_merged_bwd(                   &
                               u = real(utemp_bwd(:,:,:,:,ipoint), kind=dp),  &
                               G = G(:,:,ipoint),                   &
                               GT = GT(:,:,ipoint),                 &
                               xi = col_points_xi(:,ipoint),        &          
                               eta = col_points_eta(:,ipoint),      &
                               npol = sem_obj%npol,                 &
                               nsamp = sem_obj%ndumps,              &
                               nodes = nodes(:,:,ipoint),           &
                               element_type = element_type(ipoint), &
                               axial = axis(ipoint))
      end do                         
    end select

      !select case(strain_type)
      !case('straintensor_trace')
      !  if (any(straintrace.ne.straintrace)) print *, 'Found NaN in straintrace, ipoint:', ipoint, status_strain(ipoint)
      !case('straintensor_full')
      !  if (any(strain.ne.strain)) print *, 'Found NaN in strain, ipoint:', ipoint, status_strain(ipoint)
      !end select
    iclockold = tick(id=id_calc_strain, since=iclockold)

    select case(strain_type)
    case('straintensor_trace')
      do concurrent (ipoint = 1:npoints)
          strain_result(:, 1, :, ipoint)                      &
              = lagrange_interpol_2D_td(col_points_xi(:,ipoint),            &
                                        col_points_eta(:,ipoint),           &
                                        straintrace_bwd(:, :, :, :, ipoint), &
                                        xi(ipoint),                         &
                                        eta(ipoint))
      end do 
    case('straintensor_full')
      do concurrent (ipoint = 1:npoints)
        strain_result(:, :, :, ipoint)                    &
            = lagrange_interpol_2D_td(col_points_xi(:,ipoint),          &
                                      col_points_eta(:,ipoint),         &
                                      strain_bwd(:, :, :, :, :, ipoint), &
                                      xi(ipoint),                       &
                                      eta(ipoint))
      enddo
    end select
    iclockold = tick(id=id_lagrange, since=iclockold)



    if (trim(strain_type).eq.'straintensor_full') then
      !@TODO for consistency with SOLVER output
      strain_result(:, 4, :, :) = - strain_result(:, 4, :, :) 
      strain_result(:, 6, :, :) = - strain_result(:, 6, :, :) 
    end if

    iclockold = tick()
    if (use_buffer_loc) then
      select case(strain_type)
      case('straintensor_trace')
        do ipoint = 1, npoints
          temp = sem_obj%buffer_strain%put(id_elem(ipoint), &
                                           straintrace_bwd(:,:,:,:,ipoint))
        end do
      case('straintensor_full')
        do ipoint = 1, npoints
          temp = sem_obj%buffer_strain%put(id_elem(ipoint), &
                                           strain_bwd(:,:,:,:,:,ipoint))
        end do
      end select
    end if
    iclockold = tick(id=id_buffer, since=iclockold)
    iclockold_total = tick(id=id_load_strain, since=iclockold_total)


end subroutine load_strain_point_merged_bwd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!!subroutine load_strain_point_merged(sem_obj, xi, eta, strain_type, nodes, &
!!                                  element_type, axis, id_elem, use_buffer, strain_result)
!!    use spectral_basis, only         : lagrange_interpol_2D_td
!!    use sem_derivatives, only        : straintrace_merged_fwd, strain_merged_fwd
!!    use sem_derivatives, only        : straintrace_merged_bwd, strain_merged_bwd
!!    type(ncparamtype)               :: sem_obj
!!    real(kind=dp),     intent(in)   :: xi(:), eta(:)      !< Coordinates at which to interpolate
!!                                                    !! strain
!!    character(len=*),  intent(in)   :: strain_type  !< Model parameter (decides on 
!!                                                    !! straintrace (vp) or full tensor (vs)
!!    real(kind=dp),     intent(in)   :: nodes(:,:,:) !< Coordinates of element corner
!!                                                    !! points
!!    integer,           intent(in)   :: element_type(:) !< Element type in the solver
!!    logical,           intent(in)   :: axis(:)         !< Axis element or not 
!!    integer,           intent(in)   :: id_elem(:)      !< ID of element to interpolate strain in
!!    logical, optional, intent(in)   :: use_buffer
!!    real(kind=dp),     intent(out)  :: strain_result(:,:,:,:)
!!
!!    integer                         :: ipoint, npoints
!!
!!    real(kind=dp), allocatable      :: strain(:,:,:,:,:,:)
!!
!!    real(kind=dp), allocatable      :: straintrace(:,:,:,:,:)
!!
!!    real(kind=dp)                   :: G( 0:sem_obj%npol, 0:sem_obj%npol, size(xi,1))
!!    real(kind=dp)                   :: GT(0:sem_obj%npol, 0:sem_obj%npol, size(xi,1))
!!    real(kind=dp)                   :: col_points_xi(0:sem_obj%npol, size(xi,1))
!!    real(kind=dp)                   :: col_points_eta(0:sem_obj%npol, size(xi,1))
!!    real(kind=sp)                   :: utemp(1:sem_obj%ndumps, &
!!                                             0:sem_obj%npol,   &
!!                                             0:sem_obj%npol,   &
!!                                             1:nint(sem_obj%nsim_merged*2.5),     &
!!                                             size(xi,1)) 
!!    integer(kind=long)              :: iclockold_total, iclockold
!!    integer                         :: status_strain(size(xi,1))
!!    integer                         :: status_disp = -1
!!    integer                         :: ndirection, i, isim
!!    integer                         :: ielem_read, nelem_to_read, temp
!!    logical                         :: use_buffer_loc
!!    logical                         :: to_be_read_from_disk(size(xi,1))
!!    integer                         :: ielem_in_batch(size(xi,1))
!!    integer                         :: first_elem, last_elem
!!
!!    iclockold_total = tick()
!!    iclockold = tick()
!!
!!    npoints = size(xi,1)
!!    status_strain = -1
!!
!!    if (present(use_buffer)) then
!!      use_buffer_loc = use_buffer
!!    else
!!      use_buffer_loc = .true.
!!    end if
!!
!!
!!    if (sem_obj%nsim_merged==4) then
!!      ndirection = 10
!!    elseif (sem_obj%nsim_merged==2) then
!!      ndirection = 5
!!    end if
!!
!!
!!    select case(strain_type)
!!    case('straintensor_trace')
!!      allocate(straintrace(1:sem_obj%ndumps, &
!!                           0:sem_obj%npol,   &
!!                           0:sem_obj%npol,   &
!!                           1:sem_obj%nsim_merged, &
!!                           1:npoints))
!!      straintrace = 0
!!    case('straintensor_full')
!!      allocate(strain(1:sem_obj%ndumps, &
!!                      0:sem_obj%npol,   &
!!                      0:sem_obj%npol,   &
!!                      6,                &
!!                      1:sem_obj%nsim_merged, &
!!                      1:npoints))
!!      strain = 0            
!!    end select
!!
!!    iclockold = tick(id=id_allocate, since=iclockold)
!!    if (use_buffer_loc) then
!!      select case(strain_type)
!!      case('straintensor_trace')
!!        do ipoint = 1, npoints
!!          status_strain(ipoint) = sem_obj%buffer_strain%get(id_elem(ipoint), &
!!                                                            straintrace(:,:,:,:,ipoint))
!!        end do
!!      case('straintensor_full')
!!        do ipoint = 1, npoints
!!          status_strain(ipoint) = sem_obj%buffer_strain%get(id_elem(ipoint), &
!!                                                            strain(:,:,:,:,:,ipoint))
!!        end do
!!      end select
!!    else
!!      status_strain = - 1
!!    end if
!!
!!    to_be_read_from_disk = .true.
!!    do ipoint = 1, npoints
!!      ! Try displacement buffer
!!      iclockold = tick()
!!      status_disp = sem_obj%buffer_disp%get(id_elem(ipoint), &
!!                                            utemp(:,:,:,:,ipoint))
!!      ! If status was zero (successful), then this does not have to be read
!!      to_be_read_from_disk(ipoint) = ((status_disp.ne.0) .and.      & 
!!                                      (status_strain(ipoint).ne.0))
!!    end do
!!
!!    iclockold = tick(id=id_buffer, since=iclockold)
!!
!!    ! If there's anything to be read from disk, do it
!!    if (any(to_be_read_from_disk)) then
!!      call read_points_from_disk(sem_obj, id_elem, to_be_read_from_disk, utemp)
!!    end if
!!      
!!    iclockold = tick()
!!
!!    !if (any(utemp.ne.utemp)) print *, 'Found NaN in utemp!'
!!
!!    set_element_matrices: do ipoint = 1, npoints
!!      G(:,:,ipoint)  = sem_obj%G2
!!      col_points_eta(:, ipoint) = sem_obj%gll_points
!!      if (axis(ipoint)) then
!!        GT(:,:,ipoint) = sem_obj%G1T
!!        col_points_xi(:, ipoint)  = sem_obj%glj_points
!!      else
!!        GT(:,:,ipoint) = sem_obj%G2T
!!        col_points_xi(:, ipoint)  = sem_obj%gll_points
!!      endif 
!!    end do set_element_matrices
!!
!!    iclockold = tick(id=id_allocate_2, since=iclockold)
!!
!!    select case(strain_type)
!!    case('straintensor_trace')
!!      ! Compute just trace of strain (for vp or lambda kernels)
!!      select case(sem_obj%nsim_merged)
!!      case(2)
!!        do concurrent (ipoint = 1: npoints, status_strain(ipoint).eq.-1)
!!            straintrace(:,:,:,:,ipoint) = straintrace_merged_bwd(                 &
!!                                             u = real(utemp(:,:,:,:,ipoint), kind=dp), &
!!                                             G = G(:,:,ipoint),                   &
!!                                             GT = GT(:,:,ipoint),                 &
!!                                             xi = col_points_xi(:,ipoint),        &          
!!                                             eta = col_points_eta(:,ipoint),      &
!!                                             npol = sem_obj%npol,                 &
!!                                             nsamp = sem_obj%ndumps,              &
!!                                             nodes = nodes(:,:,ipoint),           &
!!                                             element_type = element_type(ipoint), &
!!                                             axial = axis(ipoint))
!!          !end if
!!        end do
!!      case(4)
!!        do concurrent (ipoint = 1: npoints, status_strain(ipoint).eq.-1)
!!            straintrace(:,:,:,:,ipoint) = straintrace_merged_fwd(                 &
!!                                             u = real(utemp(:,:,:,:,ipoint), kind=dp), &
!!                                             G = G(:,:,ipoint),                   &
!!                                             GT = GT(:,:,ipoint),                 &
!!                                             xi = col_points_xi(:,ipoint),        &          
!!                                             eta = col_points_eta(:,ipoint),      &
!!                                             npol = sem_obj%npol,                 &
!!                                             nsamp = sem_obj%ndumps,              &
!!                                             nodes = nodes(:,:,ipoint),           &
!!                                             element_type = element_type(ipoint), &
!!                                             axial = axis(ipoint))
!!          !end if
!!        end do
!!      end select
!!
!!    case('straintensor_full')
!!      ! compute full strain tensor
!!      select case(sem_obj%nsim_merged)
!!      case(2)
!!        do concurrent (ipoint = 1: npoints, status_strain(ipoint).eq.-1)
!!          strain(:,:,:,:,:,ipoint) = strain_merged_bwd(                   &
!!                                 u = real(utemp(:,:,:,:,ipoint), kind=dp),  &
!!                                 G = G(:,:,ipoint),                   &
!!                                 GT = GT(:,:,ipoint),                 &
!!                                 xi = col_points_xi(:,ipoint),        &          
!!                                 eta = col_points_eta(:,ipoint),      &
!!                                 npol = sem_obj%npol,                 &
!!                                 nsamp = sem_obj%ndumps,              &
!!                                 nodes = nodes(:,:,ipoint),           &
!!                                 element_type = element_type(ipoint), &
!!                                 axial = axis(ipoint))
!!        end do
!!      case(4)
!!        do concurrent (ipoint = 1: npoints, status_strain(ipoint).eq.-1)
!!          strain(:,:,:,:,:,ipoint) = strain_merged_fwd(                   &
!!                                 u = real(utemp(:,:,:,:,ipoint), kind=dp),  &
!!                                 G = G(:,:,ipoint),                   &
!!                                 GT = GT(:,:,ipoint),                 &
!!                                 xi = col_points_xi(:,ipoint),        &          
!!                                 eta = col_points_eta(:,ipoint),      &
!!                                 npol = sem_obj%npol,                 &
!!                                 nsamp = sem_obj%ndumps,              &
!!                                 nodes = nodes(:,:,ipoint),           &
!!                                 element_type = element_type(ipoint), &
!!                                 axial = axis(ipoint))
!!        end do
!!      end select
!!    end select
!!
!!      !select case(strain_type)
!!      !case('straintensor_trace')
!!      !  if (any(straintrace.ne.straintrace)) print *, 'Found NaN in straintrace, ipoint:', ipoint, status_strain(ipoint)
!!      !case('straintensor_full')
!!      !  if (any(strain.ne.strain)) print *, 'Found NaN in strain, ipoint:', ipoint, status_strain(ipoint)
!!      !end select
!!    iclockold = tick(id=id_calc_strain, since=iclockold)
!!
!!    select case(strain_type)
!!    case('straintensor_trace')
!!      do concurrent (ipoint = 1:npoints)
!!          strain_result(:, 1, :, ipoint)                      &
!!              = lagrange_interpol_2D_td(col_points_xi(:,ipoint),            &
!!                                        col_points_eta(:,ipoint),           &
!!                                        straintrace(:, :, :, :, ipoint), &
!!                                        xi(ipoint),                         &
!!                                        eta(ipoint))
!!      end do 
!!    case('straintensor_full')
!!      do concurrent (ipoint = 1:npoints)
!!        strain_result(:, :, :, ipoint)                    &
!!            = lagrange_interpol_2D_td(col_points_xi(:,ipoint),          &
!!                                      col_points_eta(:,ipoint),         &
!!                                      strain(:, :, :, :, :, ipoint), &
!!                                      xi(ipoint),                       &
!!                                      eta(ipoint))
!!      enddo
!!    end select
!!    iclockold = tick(id=id_lagrange, since=iclockold)
!!
!!
!!
!!    if (trim(strain_type).eq.'straintensor_full') then
!!      !@TODO for consistency with SOLVER output
!!      strain_result(:, 4, :, :) = - strain_result(:, 4, :, :) 
!!      strain_result(:, 6, :, :) = - strain_result(:, 6, :, :) 
!!    end if
!!
!!    iclockold = tick()
!!    if (use_buffer_loc) then
!!      select case(strain_type)
!!      case('straintensor_trace')
!!        do ipoint = 1, npoints
!!          temp = sem_obj%buffer_strain%put(id_elem(ipoint), &
!!                                           straintrace(:,:,:,:,ipoint))
!!        end do
!!      case('straintensor_full')
!!        do ipoint = 1, npoints
!!          temp = sem_obj%buffer_strain%put(id_elem(ipoint), &
!!                                           strain(:,:,:,:,:,ipoint))
!!        end do
!!      end select
!!    end if
!!    iclockold = tick(id=id_buffer, since=iclockold)
!!    iclockold_total = tick(id=id_load_strain, since=iclockold_total)
!!
!!
!!end subroutine load_strain_point_merged
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_points_from_disk(sem_obj, id_elem, to_be_read_from_disk, u)
  use netcdf, only              : nf90_get_var ! HACK, create 5D wrapper in nc_routines 
  type(ncparamtype)            :: sem_obj
  integer, intent(in)          :: id_elem(:)
  logical, intent(inout)       :: to_be_read_from_disk(:)
  real(kind=sp), intent(inout) :: u(:,:,:,:,:)

  integer                      :: first_elem, last_elem, nelem_to_read
  integer                      :: ielem_read, ipoint, status_disp
  integer                      :: npoints, ndirection
  integer                      :: ielem_in_batch(size(id_elem,1))
  integer(kind=long)           :: iclockold

  ndirection = size(u, 4)
  npoints = size(id_elem)

  first_elem = minval(id_elem, mask=to_be_read_from_disk)
  last_elem  = maxval(id_elem, mask=to_be_read_from_disk) 
  nelem_to_read = last_elem - first_elem + 1

  !write(lu_out, *) 'nelem: ', nelem_to_read, 'first: ', first_elem, 'last: ', last_elem
  !write(*, *) 'nelem: ', nelem_to_read, 'first: ', first_elem, 'last: ', last_elem, ' ids:', id_elem
  u_batch = 0

  nread_larger_max: if (nelem_to_read < nelem_to_read_max) then
    !write(lu_out, *) 'reading all at once'
    !write(*, *) 'reading all at once'
    do ipoint = 1, npoints
      if (to_be_read_from_disk(ipoint)) then
        ielem_in_batch(ipoint) = id_elem(ipoint) + 1 - first_elem
      end if
    end do

    iclockold = tick()
    call check(nf90_get_var(ncid   = sem_obj%ncid,                      & 
                            varid  = sem_obj%mergedvarid,               &
                            start  = [1, 1, 1, 1, first_elem],          &
                            count  = [sem_obj%ndumps, sem_obj%npol+1,   &
                                      sem_obj%npol+1, ndirection,       &
                                      nelem_to_read],                   &
                            values = u_batch(:,:,:,1:ndirection,1:nelem_to_read)))
    iclockold = tick(id=id_netcdf, since=iclockold)

    do ielem_read = 1, nelem_to_read
      status_disp = sem_obj%buffer_disp%put(first_elem + ielem_read - 1, &
                                            u_batch(:,                   &
                                                    :,                   &
                                                    :,                   &
                                                    1:ndirection,        &
                                                    ielem_read))
    end do
    iclockold = tick(id=id_buffer, since=iclockold)

    do ipoint = 1, npoints
      if (to_be_read_from_disk(ipoint)) then
        u(:,:,:,:,ipoint) = u_batch(:,:,:,1:ndirection,ielem_in_batch(ipoint))
      end if
    end do

  else 
    !write(lu_out, *) 'reading one at a time'
    !write(*, *) 'reading one at a time'
    nelem_to_read = npoints
    iclockold = tick()

    do ipoint = 1, npoints
      iclockold = tick()
      call check(nf90_get_var(ncid   = sem_obj%ncid,                      & 
                              varid  = sem_obj%mergedvarid,               &
                              start  = [1, 1, 1, 1, id_elem(ipoint)],     &
                              count  = [sem_obj%ndumps, sem_obj%npol+1,   &
                                        sem_obj%npol+1, ndirection,       &
                                        1],                               &
                              values = u_batch(1:sem_obj%ndumps,          &
                                               0:sem_obj%npol,            &
                                               0:sem_obj%npol,            &
                                               1:ndirection,              &
                                               1)))
      iclockold = tick(id=id_netcdf, since=iclockold)
      status_disp = sem_obj%buffer_disp%put(id_elem(ipoint),             &
                                            u_batch(:,                   &
                                                    :,                   &
                                                    :,                   &
                                                    1:ndirection,        &
                                                    1))
      iclockold = tick(id=id_buffer, since=iclockold)

      if (to_be_read_from_disk(ipoint)) then
        u(:,:,:,:,ipoint) = u_batch(:, :, :, 1:ndirection, 1)
      end if

    end do

  end if nread_larger_max

end subroutine read_points_from_disk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine build_kdtree(this)
    use kdtree2_module, only    : kdtree2_create, kdtree2_destroy
    class(semdata_type)        :: this
    real(kind=sp), allocatable :: mesh(:,:)

    !if (.not.this%meshes_read) then
    !    print *, 'ERROR in build_kdtree(): Meshes have not been read yet'
    !    print *, 'Call read_meshes() before build_kdtree!'
    !    call pabort
    !end if

    ! Destroy kdtree
    if (this%kdtree_built) then
      if (verbose>0) then 
         print *, 'WARNING in build_kdtree(): Meshes have already been built'
         print *, 'Destroying the old trees...'
      end if
      call kdtree2_destroy(this%fwdtree)
      call kdtree2_destroy(this%bwdtree)
    end if

    allocate(mesh(2, this%fwdmesh%npoints))
    mesh = transpose(reshape([this%fwdmesh%s, this%fwdmesh%z],       &
                             [this%fwdmesh%npoints, 2]))


    write(lu_out,*) ' Building forward KD-Tree'
    call flush(lu_out)
    ! KDtree in forward field
    this%fwdtree => kdtree2_create(mesh,              &
                                   dim = 2,           &
                                   sort = .true.,     &
                                   rearrange = .true.)
    deallocate(mesh)                           


    write(lu_out,*) ' Building forward midpoint-only KD-Tree'
    allocate(mesh(2, this%fwdmesh%nelem)) ! midpoints only
    mesh = transpose(reshape([this%fwdmesh%s_mp, this%fwdmesh%z_mp], &
                             [this%fwdmesh%nelem, 2]))
    ! KDtree in forward midpoint field
    this%fwdtree_mp => kdtree2_create(mesh,              &
                                      dim = 2,           &
                                      sort = .true.,     &
                                      rearrange = .true.)
    deallocate(mesh)                           

    

    ! KDtree in backward field

    allocate(mesh(2, this%bwdmesh%npoints))
    mesh = transpose(reshape([this%bwdmesh%s, this%bwdmesh%z],       &
                             [this%bwdmesh%npoints, 2]))

    write(lu_out,*) ' Building backward KD-Tree'
    call flush(lu_out)
    this%bwdtree => kdtree2_create(mesh,              &
                                   dim = 2,           &
                                   sort = .true.,     &
                                   rearrange = .true.)
    deallocate(mesh)                           

    write(lu_out,*) ' Building backward midpoint-only KD-Tree'
    allocate(mesh(2, this%bwdmesh%nelem)) ! midpoints only
    mesh = transpose(reshape([this%bwdmesh%s_mp, this%bwdmesh%z_mp], &
                             [this%bwdmesh%nelem, 2]))
    ! KDtree in backward midpoint field
    this%bwdtree_mp => kdtree2_create(mesh,              &
                                      dim = 2,           &
                                      sort = .true.,     &
                                      rearrange = .true.)
    deallocate(mesh)

    call flush(lu_out)

    this%kdtree_built = .true.

end subroutine build_kdtree
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_meshes(this)
   use netcdf

   class(semdata_type)        :: this
   integer                    :: isim
   
   if (.not.this%files_open) then
       print *, 'ERROR in read_meshes(): Files have not been opened!'
       print *, 'Call open_files() before read_meshes()'
       call pabort
   end if

   ! Forward SEM mesh
   write(lu_out,*) '  Read SEM mesh from first forward simulation'
   
   call nc_read_att_int(this%fwdmesh%npoints, 'npoints', this%fwd(1))
   call nc_read_att_int(this%fwdmesh%nelem, 'nelem_kwf_global', this%fwd(1))
   write(lu_out, *) 'Mesh has ', this%fwdmesh%npoints, ' points, ', &
                                 this%fwdmesh%nelem, ' elements'
   
   do isim = 1, this%nfiles_fwd
      this%fwd(isim)%ngll = this%fwdmesh%npoints
      this%fwd(isim)%nelem = this%fwdmesh%nelem
   end do

   call cache_mesh(this%fwd(1)%mesh, this%fwdmesh) 
   
   ! Backward SEM mesh                     
   write(lu_out,*) 'Read SEM mesh from first backward simulation'
   
   call nc_read_att_int(this%bwdmesh%npoints, 'npoints', this%bwd(1))
   call nc_read_att_int(this%bwdmesh%nelem, 'nelem_kwf_global', this%bwd(1))
   write(lu_out, *) 'Mesh has ', this%fwdmesh%npoints, ' points, ', &
                                 this%fwdmesh%nelem, ' elements'
 
   do isim = 1, this%nfiles_bwd
      this%bwd(isim)%ngll = this%bwdmesh%npoints
      this%bwd(isim)%nelem = this%bwdmesh%nelem
   end do
   
   call cache_mesh(this%bwd(1)%mesh, this%bwdmesh) 
                             
   ! define terms needed to compute gradient
   call calc_gradient_terms(this)

   ! Build KDTree
   write(lu_out, *) 'Build KD-Trees'
   call flush(lu_out)
   call this%build_kdtree()

   ! Load mesh model 
   write(lu_out, *) 'Load forward mesh model parameters and create interpolation objects'
   call flush(lu_out)
   call load_model_parameter(this%fwd(1)%mesh, this%fwdmesh, this%fwdtree, & 
                             this%fwd(1)%planet_radius)

   write(lu_out, *) 'Load backward mesh model parameters and create interpolation objects'
   call flush(lu_out)
   call load_model_parameter(this%bwd(1)%mesh, this%bwdmesh, this%bwdtree, & 
                             this%bwd(1)%planet_radius)

   this%meshes_read = .true.

   write(lu_out, *) 'Forward and backward SEM mesh reading succeeded'
   call flush(lu_out)

end subroutine read_meshes
!-----------------------------------------------------------------------------------------
 
!-----------------------------------------------------------------------------------------
!> Read and cache mesh variables
subroutine cache_mesh(ncid, mesh)
  use mpi_nc_routines,      only : mpi_getvar_by_name
  use commpi,               only : MPI_COMM_NODE
  integer, intent(in)           :: ncid
  type(meshtype)                :: mesh

  write(lu_out, '(A)', advance='no') 'Reading mesh parameters...'
  call flush(lu_out)

  call mpi_getvar_by_name(ncid   = ncid,          &
                         varname = 'mesh_S',      &
                         limits  = [0., 1e9],     & 
                         comm    = MPI_COMM_NODE, &
                         values  = mesh%s   )

  !write(lu_out,*) 'sizeof mesh%s: ', sizeof(mesh%s)

              
  call mpi_getvar_by_name(ncid   = ncid,          &
                         varname   = 'mesh_Z',      &
                         limits = [-1e9, 1e9],   & 
                         comm      = MPI_COMM_NODE,    &
                         values = mesh%z   )

  !write(lu_out,*) 'sizeof mesh%s: ', sizeof(mesh%s)
              

  
  call mpi_getvar_by_name(ncid   = ncid,         &
                         varname   = 'eltype',     &
                         limits = [0, 3],       &
                         comm      = MPI_COMM_NODE,    &
                         values = mesh%eltype)
  !write(lu_out,*) 'sizeof mesh%eltype: ', sizeof(mesh%eltype)

  call mpi_getvar_by_name(ncid   = ncid,         &
                         varname   = 'axis',       &
                         limits = [0, 1],       &
                         comm      = MPI_COMM_NODE,    &
                         values = mesh%isaxis)
  !write(lu_out,*) 'sizeof mesh%isaxis: ', sizeof(mesh%isaxis)

  call mpi_getvar_by_name(ncid   = ncid,         &
                         varname   = 'mp_mesh_S',  &
                         limits = [0., 1e9],    & 
                         comm      = MPI_COMM_NODE,    &
                         values = mesh%s_mp )
  !write(lu_out,*) 'sizeof mesh%s_mp: ', sizeof(mesh%s_mp)
              
  call mpi_getvar_by_name(ncid   = ncid,         &
                         varname   = 'mp_mesh_Z',  &
                         limits = [-1e9, 1e9],  & 
                         comm      = MPI_COMM_NODE,    &
                         values = mesh%z_mp )
  !write(lu_out,*) 'sizeof mesh%z_mp: ', sizeof(mesh%z_mp)

  call mpi_getvar_by_name(ncid   = ncid,         &
                         varname   = 'fem_mesh',   &
                         limits = [0, size(mesh%s)-1], &
                         comm      = MPI_COMM_NODE,    &
                         values = mesh%corner_point_ids )
  !write(lu_out,*) 'sizeof mesh%corner_point_ids: ', sizeof(mesh%corner_point_ids)

  call mpi_getvar_by_name(ncid   = ncid,         &
                         varname   = 'sem_mesh',   &
                         limits = [0, size(mesh%s)-1], &
                         comm      = MPI_COMM_NODE,    &
                         values = mesh%gll_point_ids)
  !write(lu_out,*) 'sizeof mesh%gll_point_ids: ', sizeof(mesh%gll_point_ids)


  write(lu_out, *) ' done'

end subroutine cache_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine load_model_parameter(ncid, mesh, tree, radius)
  use mpi_nc_routines,   only : mpi_getvar_by_name
  use commpi, only            : MPI_COMM_NODE
  use interpolate_mesh, only  : create_interpolator
  type(meshtype)             :: mesh
  integer, intent(in)        :: ncid
  real(kind=dp), intent(in)  :: radius
  type(kdtree2), pointer     :: tree
  real(kind=sp), allocatable :: param_tmp(:)

  call mpi_getvar_by_name(ncid   = ncid,          &
                          varname   = 'mesh_vp',     &
                          limits = [0.0, 2e4],    & 
                          comm   = MPI_COMM_NODE, &
                          values = param_tmp  )
  mesh%vp = create_interpolator(param_tmp, tree, radius)
  mesh%vp_max = maxval(param_tmp)
  deallocate(param_tmp)
              
  call mpi_getvar_by_name(ncid   = ncid,          &
                          varname   = 'mesh_vs',     &
                          limits = [0.0, 2e4],    & 
                          comm   = MPI_COMM_NODE, &
                          values = param_tmp)
  mesh%vs = create_interpolator(param_tmp, tree, radius)
  mesh%vs_max = maxval(param_tmp)
  deallocate(param_tmp)
              
  call mpi_getvar_by_name(ncid   = ncid,          &
                          varname   = 'mesh_rho',    &
                          limits = [0.0, 2e4],    & 
                          comm   = MPI_COMM_NODE, &
                          values = param_tmp)
  mesh%rho = create_interpolator(param_tmp, tree, radius)
  deallocate(param_tmp)
              
  call mpi_getvar_by_name(ncid   = ncid,          &
                          varname   = 'mesh_lambda', &
                          limits = [1e9, 1e15],   & 
                          comm   = MPI_COMM_NODE, &
                          values = param_tmp)
  mesh%lambda = create_interpolator(param_tmp, tree, radius)
  deallocate(param_tmp)
              
  call mpi_getvar_by_name(ncid   = ncid,          &
                          varname   = 'mesh_mu',     &
                          limits = [0.0, 1e12],   & 
                          comm   = MPI_COMM_NODE, &
                          values = param_tmp)
  mesh%mu = create_interpolator(param_tmp, tree, radius)
  deallocate(param_tmp)
              
  call mpi_getvar_by_name(ncid   = ncid,          &
                          varname   = 'mesh_phi',    &
                          limits = [0.0, 3.0],    & 
                          comm   = MPI_COMM_NODE, &
                          values = param_tmp)
  mesh%phi = create_interpolator(param_tmp, tree, radius)
  deallocate(param_tmp)
              
  call mpi_getvar_by_name(ncid   = ncid,          &
                          varname   = 'mesh_xi',  &
                          limits = [0.0, 3.0],    & 
                          comm   = MPI_COMM_NODE, &
                          values = param_tmp)
  mesh%xi = create_interpolator(param_tmp, tree, radius)
  deallocate(param_tmp)
              
  call mpi_getvar_by_name(ncid   = ncid,          &
                          varname   = 'mesh_eta',    &
                          limits = [0.0, 1e12],   & 
                          comm   = MPI_COMM_NODE, &
                          values = param_tmp)
  mesh%eta = create_interpolator(param_tmp, tree, radius)
  deallocate(param_tmp)

end subroutine load_model_parameter 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Calculates the terms G1, G2, G1T, G2T needed to compute gradients and assigns them
!! to the variables of type ncparamtype of forward and backward field
subroutine calc_gradient_terms(sem_var)
  use spectral_basis, only : zelegl, zemngl2, &
                             def_lagrange_derivs_gll, def_lagrange_derivs_glj
  class(semdata_type)     :: sem_var

  integer                 :: ifile

  write(lu_out, '(" Calculating gradient terms for npol = ", I1, "...")') sem_var%npol

  allocate(sem_var%G1(0:sem_var%npol,0:sem_var%npol))
  allocate(sem_var%G1T(0:sem_var%npol,0:sem_var%npol))
  allocate(sem_var%G2(0:sem_var%npol,0:sem_var%npol))
  allocate(sem_var%G2T(0:sem_var%npol,0:sem_var%npol))
  allocate(sem_var%G0(0:sem_var%npol))

  allocate(sem_var%gll_points(0:sem_var%npol))
  allocate(sem_var%glj_points(0:sem_var%npol))

  sem_var%gll_points = zelegl(sem_var%npol)
  sem_var%glj_points = zemngl2(sem_var%npol)

  sem_var%G1 = def_lagrange_derivs_glj(sem_var%npol, sem_var%G0)
  sem_var%G2 = def_lagrange_derivs_gll(sem_var%npol)

  sem_var%G1T = transpose(sem_var%G1)
  sem_var%G2T = transpose(sem_var%G2)

  do ifile = 1, sem_var%nfiles_fwd
     allocate(sem_var%fwd(ifile)%gll_points(0:sem_var%npol))
     allocate(sem_var%fwd(ifile)%glj_points(0:sem_var%npol))

     sem_var%fwd(ifile)%gll_points = sem_var%gll_points
     sem_var%fwd(ifile)%glj_points = sem_var%glj_points

     allocate(sem_var%fwd(ifile)%G1(0:sem_var%npol,0:sem_var%npol))
     allocate(sem_var%fwd(ifile)%G1T(0:sem_var%npol,0:sem_var%npol))
     allocate(sem_var%fwd(ifile)%G2(0:sem_var%npol,0:sem_var%npol))
     allocate(sem_var%fwd(ifile)%G2T(0:sem_var%npol,0:sem_var%npol))
     allocate(sem_var%fwd(ifile)%G0(0:sem_var%npol))

     sem_var%fwd(ifile)%G1 = sem_var%G1
     sem_var%fwd(ifile)%G2 = sem_var%G2
     sem_var%fwd(ifile)%G1T = sem_var%G1T
     sem_var%fwd(ifile)%G2T = sem_var%G2T
     sem_var%fwd(ifile)%G0 = sem_var%G0
  end do

  do ifile = 1, sem_var%nfiles_bwd
     allocate(sem_var%bwd(ifile)%gll_points(0:sem_var%npol))
     allocate(sem_var%bwd(ifile)%glj_points(0:sem_var%npol))
     
     sem_var%bwd(ifile)%gll_points = sem_var%gll_points
     sem_var%bwd(ifile)%glj_points = sem_var%glj_points
     
     allocate(sem_var%bwd(ifile)%G1(0:sem_var%npol,0:sem_var%npol))
     allocate(sem_var%bwd(ifile)%G1T(0:sem_var%npol,0:sem_var%npol))
     allocate(sem_var%bwd(ifile)%G2(0:sem_var%npol,0:sem_var%npol))
     allocate(sem_var%bwd(ifile)%G2T(0:sem_var%npol,0:sem_var%npol))
     allocate(sem_var%bwd(ifile)%G0(0:sem_var%npol))
     
     sem_var%bwd(ifile)%G1 = sem_var%G1
     sem_var%bwd(ifile)%G2 = sem_var%G2
     sem_var%bwd(ifile)%G1T = sem_var%G1T
     sem_var%bwd(ifile)%G2T = sem_var%G2T
     sem_var%bwd(ifile)%G0 = sem_var%G0
  end do
end subroutine calc_gradient_terms
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Calculates start and count to read one chunk (and not exceed size of variable)
subroutine get_chunk_bounds(pointid, chunksize, npoints, start_chunk, count_chunk)
  integer, intent(in)       :: pointid      !< ID of point to read
  integer, intent(in)       :: chunksize    !< Chunk size of variable
  integer, intent(in)       :: npoints      !< Size of variable
  integer, intent(out)      :: start_chunk  !< Start of chunk in which pointid is
  integer, intent(out)      :: count_chunk  !< Size of chunk in which pointid is 
                                            !! Normally == chunksize, but should not be larger
                                            !! than npoints
  integer                   :: ichunk       !< Number of chunk to read (starting from 0)

  if ((pointid > npoints).or.(pointid<1)) then
    print *, 'ERROR: Requesting chunk for point ', pointid
    print *, '        variable bounds are 0 and ', npoints
    call pabort()
  end if

  ! Chunk 1 ranges from 1 to chunksize, 
  ! Chunk i from (i-1)*chunksize+1 to i*chunksize
  ichunk = (pointid - 1) / chunksize  ! Integer division

  start_chunk = ichunk * chunksize + 1
  
  count_chunk = min(chunksize, npoints - start_chunk + 1)
                                            
end subroutine get_chunk_bounds
!-----------------------------------------------------------------------------------------
 
!-----------------------------------------------------------------------------------------
!> Dampen field variable around a central point
subroutine dampen_field(field, r_points, r_src_rec, r_max)
  real(kind=dp), intent(inout)      :: field(:,:,:)    !< Variable to dampen
  real(kind=dp), intent(in)         :: r_points(:,:)   !< Locations of points
  real(kind=dp), intent(in)         :: r_src_rec(3)    !< Location of damping center
  real(kind=dp), intent(in)         :: r_max           !< Distance at which damping starts

  real(kind=dp)                     :: dist
  integer                           :: ipoint, npoints 
  
  ! Only damp, if r_max is larger zero
  if (r_max > 0.0d0) then
    npoints = size(field,1)
    do ipoint = 1, npoints
      dist = norm2(r_points(:,ipoint) - r_src_rec)
      if (dist<r_max) then
        field(ipoint,:,:) = field(ipoint,:,:) * dist / r_max
      end if
    end do
  end if

end subroutine dampen_field
!-----------------------------------------------------------------------------------------
 
!-----------------------------------------------------------------------------------------
subroutine get_simulation_type(sim_dir, merged, nsim, nfiles)
  ! Get simulation type of backward simulation and whether it was stored in 
  ! separate directories or in one merged file
  ! Convention for nsim: If merged=false (separate databases for each run of moment source
  !                      then nsim is equal to the rank of nc_obj
  !                      If merged=true, rank of nc_obj is always one

  character(len=*), intent(in)   :: sim_dir
  logical, intent(out)           :: merged
  integer, intent(out)           :: nsim
  integer, intent(out)           :: nfiles
  logical                        :: moment=.false., force=.false.
  logical                        :: single=.false.
  character(len=512)             :: dirnam

  dirnam = trim(sim_dir)//'/merged_instaseis_db.nc4'
  write(lu_out,*) 'Inquiring: ', trim(dirnam)
  inquire( file = trim(dirnam), exist = merged)
  
  dirnam = trim(sim_dir)//'/MZZ/Data/ordered_output.nc4'
  write(lu_out,*) 'Inquiring: ', trim(dirnam)
  inquire( file = trim(dirnam), exist = moment)

  dirnam = trim(sim_dir)//'/PZ/Data/ordered_output.nc4'
  write(lu_out,*) 'Inquiring: ', trim(dirnam)
  inquire( file = trim(dirnam), exist = force)

  dirnam = trim(sim_dir)//'/Data/ordered_output.nc4'
  write(lu_out,*) 'Inquiring: ', trim(dirnam)
  inquire( file = trim(dirnam), exist = single)
 
  if (moment) then
     nsim = 4
     nfiles = 4
     write(lu_out,*) 'Simulation was ''moment'' source'
  elseif (force) then
     nsim = 2
     nfiles = 2
     write(lu_out,*) 'Simulation was ''forces'' source'
  elseif (single) then
     nsim = 1
     nfiles = 1
     write(lu_out,*) 'Simulation was ''single'' source'
  elseif (merged) then
     nsim = get_nsim_from_merged_file(trim(sim_dir) &
                                      //'/merged_instaseis_db.nc4')
     nfiles = 1
     write(lu_out,*) 'Simulation in merged file, nsim=', nsim
  else 
     write(*,*) 'ERROR: Run directory (as set in inparam_basic)'
     write(*,*) trim(sim_dir)
     write(*,*) 'does not seem to be an axisem rundirectory'
     call pabort(do_traceback=.false.)
  end if

end subroutine get_simulation_type
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine set_simulation_paths(nc_obj, sim_dir, nfiles)
  type(ncparamtype)            :: nc_obj(:)
  character(len=*), intent(in) :: sim_dir
  integer, intent(in)          :: nfiles 

  select case(nfiles)
  case(1)    ! Single or merged
      nc_obj(1)%meshdir = sim_dir//'/'

  case(2)    ! Forces
      nc_obj(1)%meshdir = trim(sim_dir)//'/PZ/'
      nc_obj(2)%meshdir = trim(sim_dir)//'/PX/'

  case(4)    ! Moment
      nc_obj(1)%meshdir = trim(sim_dir)//'/MZZ/'
      nc_obj(2)%meshdir = trim(sim_dir)//'/MXX_P_MYY/'
      nc_obj(3)%meshdir = trim(sim_dir)//'/MXZ_MYZ/'
      nc_obj(4)%meshdir = trim(sim_dir)//'/MXY_MXX_M_MYY/'
  end select

end subroutine set_simulation_paths
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_nsim_from_merged_file(filename)
  use netcdf,     only               : nf90_get_att, NF90_GLOBAL, NF90_NOERR  
  use nc_routines, only             : nc_close_file, nc_open_for_read
  
  character(len=*), intent(in)      :: filename
  
  integer                           :: ncid, status

  call nc_open_for_read(filename = filename, ncid = ncid) 
  status = nf90_get_att(ncid, NF90_GLOBAL, 'nsim', get_nsim_from_merged_file)
  call nc_close_file(ncid)

end function get_nsim_from_merged_file
!-----------------------------------------------------------------------------------------
 
!-----------------------------------------------------------------------------------------
!> Read NetCDF attribute of type Integer
subroutine nc_read_att_int(attribute_value, attribute_name, nc)
  use netcdf,     only               : nf90_get_att, NF90_GLOBAL, NF90_NOERR  
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
  use netcdf,     only               : nf90_get_att, NF90_GLOBAL, NF90_NOERR  
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
  use netcdf,     only               : nf90_get_att, NF90_GLOBAL, NF90_NOERR  
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
  use netcdf,     only               : nf90_get_att, NF90_GLOBAL, NF90_NOERR  
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
