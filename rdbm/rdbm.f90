!=========================================================================================
program rdbm

  use readfields,           only : semdata_type
  use commpi
  use global_parameters
  use source_class
  use resampling
  use filtering,            only : filter_type
  use type_parameter
  use receivers_rdbm
  use inversion_mesh
  use fft, only : taperandzeropad
  use clocks_mod
  use progressbar_mod

  implicit none

  type(semdata_type)                  :: sem_data
  type(parameter_type)                :: parameters
  type(receivers_rdbm_type)           :: receivers
  type(src_param_type), allocatable   :: sources(:)
  type(progressbar)                   :: bar
  type(inversion_mesh_data_type)      :: mesh

  character(len=512)                  :: bwd_dir
  character(len=4)                    :: model_param
  real(kind=dp), allocatable          :: fw_field(:,:,:)
  integer                             :: i, j
  integer                             :: lu_seis
  character(len=10)                   :: fname

  type(resampling_type)               :: resamp
  type(filter_type)                   :: filter
  character(len=32)                   :: filtername
  real(kind=dp)                       :: filter_df
  integer                             :: filter_nomega
  real(kind=dp), allocatable          :: fw_field_res(:,:)
  real(kind=dp), allocatable          :: mu(:)

  integer                             :: iclockold

  real(kind=dp), dimension(:), allocatable      :: T
  real(kind=dp)                                 :: dt_out

  call start_clock()
  iclockold = tick()

  ! read input files
  call parameters%read_parameters('inparam_basic')

  if (parameters%receiver_file_type == 'stations') then
     call receivers%read_stations_file(fname=parameters%receiver_file)
  else if (parameters%receiver_file_type == 'colatlon') then
     call receivers%read_receiver_dat(fname=parameters%receiver_file)
  else if (parameters%receiver_file_type == 'abaqus') then
     call mesh%read_abaqus_mesh(parameters%receiver_file, 'onvertices')
     call receivers%init_xyz(mesh%get_vertices_dp())
  else
     write(6,*) 'ERROR: unknown receiver file type'
     call pabort
  endif

  ! open heavy data files and initialize meshes
  bwd_dir = ''

  if (trim(parameters%source_type) == 'explosion') then
     model_param = 'vp'
  else
     model_param = 'vs'
  endif

  call sem_data%set_params(parameters%sim_dir, bwd_dir, &
                           parameters%buffer_size, model_param)
  call sem_data%open_files()
  call sem_data%read_meshes()
  call sem_data%build_kdtree()


  ! initialize sources
  if (trim(parameters%source_file_type) == 'cmtsolution') then

     allocate(sources(1))
     call sources(1)%read_cmtsolution(fname=trim(parameters%source_file_name))
     sources(1)%shift_time_sample = (sources(1)%shift_time - sem_data%timeshift_fwd) / sem_data%dt

  else if (trim(parameters%source_file_type) == 'cmtsolutions') then

     allocate(sources(parameters%nsources))
     do i = 1, parameters%nsources
        call sources(i)%read_cmtsolution(fname=trim(parameters%source_files(i)))
        sources(i)%shift_time_sample = (sources(i)%shift_time - sem_data%timeshift_fwd) / sem_data%dt
     enddo

  else if (trim(parameters%source_file_type) == 'srf') then

     write(6,*) 'reading srf file'

     call read_srf(parameters%source_file_name, sources, nsources=parameters%nsources)

     write(6,*) parameters%nsources
     do i = 1, parameters%nsources
        sources(i)%shift_time_sample = (sources(i)%shift_time - sem_data%timeshift_fwd) / sem_data%dt
     enddo
     allocate(mu(parameters%nsources))
  else
     write(6,*) 'ERROR: unkown source file type ', parameters%source_file_type
     call pabort
  endif

  ! initialize filter
  if (parameters%apply_filter) then
     filtername = 'filter1'
     filter_nomega = sem_data%ndumps + 1
     filter_df = 0.5d0 / sem_data%dt / filter_nomega
     call filter%create(filtername, filter_df, filter_nomega, parameters%filterclass, parameters%filterfreqs)
  endif

  ! initialize resampling
  if (parameters%resample .or. parameters%time_shift) then
     allocate(fw_field_res(parameters%nsamp * 2, parameters%nsources))
     if (parameters%apply_filter) then
        call resamp%init(sem_data%ndumps * 2, parameters%nsamp * 2, parameters%nsources, filter=filter)
     else
        call resamp%init(sem_data%ndumps * 2, parameters%nsamp * 2, parameters%nsources)
     endif
  else
     parameters%nsamp = sem_data%ndumps
     allocate(fw_field_res(parameters%nsamp, parameters%nsources))
  endif

  if (parameters%receiver_file_type == 'abaqus') then
     call mesh%init_node_data(parameters%nsamp)
  endif

  allocate(fw_field(sem_data%ndumps, 1, parameters%nsources))

  ! initialize time traces
  allocate(T(1:parameters%nsamp))
  dt_out = sem_data%dt * sem_data%ndumps / parameters%nsamp

  do i = 1, parameters%nsamp
     T(i) = dt_out * (i - 1)
  end do

  ! initialization done
  iclockold = tick(id=id_init, since=iclockold)

  write(6,*) 'Initiatlization Done.'
  write(6,*) 'Working on ', receivers%num_rec, ' receivers'

  do i=1, receivers%num_rec

     !if (receivers%reci_sources(i)%lond < -30 .or. receivers%reci_sources(i)%lond > 30 &
     !      .or. receivers%reci_sources(i)%latd < -30 .or. receivers%reci_sources(i)%latd > 30) then
     !   cycle
     !endif

     ! load data from file
     iclockold = tick()
     if (trim(parameters%source_file_type) == 'srf') then
        fw_field = sem_data%load_fw_points_rdbm(sources, receivers%reci_sources(i), &
                                                parameters%component, mu=mu)
        do j = 1, parameters%nsources
           ! correcting shear modulus used during reading the source
           fw_field(:,:,j) = fw_field(:,:,j) * mu(j) / 32d9
        enddo
     else
        fw_field = sem_data%load_fw_points_rdbm(sources, receivers%reci_sources(i), &
                                                parameters%component)
     endif
     iclockold = tick(id=id_load, since=iclockold)

     ! resample
     if (parameters%resample .and. .not. parameters%time_shift) then
        iclockold = tick()
        call resamp%resample(taperandzeropad(fw_field(:,1,:), ntaper=10, &
                                             ntimes=sem_data%ndumps * 2, &
                                             end_only=.true.), &
                             fw_field_res)
        iclockold = tick(id=id_resamp, since=iclockold)
     elseif (parameters%time_shift) then
        iclockold = tick()
        call resamp%resample_timeshift(taperandzeropad(fw_field(:,1,:), ntaper=10, &
                                                       ntimes=sem_data%ndumps * 2, &
                                                       end_only=.true.), &
                                       fw_field_res, sources=sources)
        iclockold = tick(id=id_resamp, since=iclockold)
     else
        fw_field_res = fw_field(:,1,:)
     endif

     !! ouput to file
     iclockold = tick()
     if (parameters%receiver_file_type == 'abaqus') then
        !TODO: for now summing over all sources, would it be useful not to?
        call mesh%set_node_data_trace(real(sum(fw_field_res(1:parameters%nsamp,:), dim=2), kind=sp), i)
     else
        write(fname,'("seis_",I0.5)') i
        open(newunit=lu_seis, file=fname)

        if (trim(parameters%mode) == 'tomography') then
           do j = 1, parameters%nsamp
              write(lu_seis,*) T(j), fw_field_res(j,:)
           enddo
        elseif (trim(parameters%mode) == 'finitesource') then
           do j = 1, parameters%nsamp
              write(lu_seis,*) T(j), sum(fw_field_res(j,:))
           enddo
        else
           write(6,*) 'ERROR: unknown mode "', trim(parameters%mode), '"'
           call pabort
        endif
        close(lu_seis)
     endif
     iclockold = tick(id=id_out, since=iclockold)

     call bar%printbar(100 * i / receivers%num_rec)
  enddo

  if (parameters%receiver_file_type == 'abaqus') then
     !TODO: read filename from inparam
     call mesh%dump_node_data_xdmf('testdata')
  endif

  ! finish
  call sem_data%close_files
  call end_clock()

contains

!-----------------------------------------------------------------------------------------
subroutine start_clock
  !
  ! Driver routine to start the timing, using the clocks_mod module.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  use global_parameters
  use clocks_mod, only : clock_id, clocks_init
  
  implicit none
  
  character(len=8)  :: mydate
  character(len=10) :: mytime

  call date_and_time(mydate,mytime) 
  write(lu_out,11) mydate(5:6), mydate(7:8), mydate(1:4), mytime(1:2), mytime(3:4)

11 format('     Kerner started on ', A2,'/',A2,'/',A4,' at ', A2,'h ',A2,'min',/)

  call clocks_init(0)

  id_init   = clock_id('initialization')
  id_load   = clock_id('loading seismograms from file')
  id_resamp = clock_id('resampling')
  id_out    = clock_id('output')

end subroutine start_clock
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine end_clock 
  !
  ! Wapper routine to end timing and display clock informations.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  use global_parameters, only : lu_out
  use clocks_mod,        only : clocks_exit

  implicit none

  write(lu_out,*)
  write(lu_out,"(10x,'Summary of timing measurements:')")
  write(lu_out,*)

  call clocks_exit(0)

  write(lu_out,*)

end subroutine end_clock
!-----------------------------------------------------------------------------------------

end program
!=========================================================================================
