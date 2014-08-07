!=========================================================================================
program rdbm

  use readfields, only : semdata_type
  use commpi
  use global_parameters
  use source_class
  use resampling
  use type_parameter
  use receivers_rdbm
  use fft, only : taperandzeropad
  use clocks_mod

  implicit none

  type(semdata_type)                  :: sem_data
  type(parameter_type)                :: parameters
  type(receivers_rdbm_type)           :: receivers
  type(src_param_type), allocatable   :: sources(:)

  character(len=512)                  :: bwd_dir
  character(len=4)                    :: model_param
  real(kind=dp), allocatable          :: fw_field(:,:,:)
  integer                             :: i, j
  integer                             :: lu_seis
  character(len=8)                    :: fname

  type(resampling_type)               :: resamp
  real(kind=dp), allocatable          :: fw_field_res(:,:)

  integer                             :: iclockold

  real(kind=dp), dimension(:), allocatable      :: T
  real(kind=dp)                                 :: dt_out

  call start_clock()
  iclockold = tick()

  ! read input files
  call parameters%read_parameters('inparam_basic')

  if (parameters%receiver_file_type == 'stations') then
     call receivers%read_stations_file()
  else if (parameters%receiver_file_type == 'colatlon') then
     call receivers%read_receiver_dat()
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

  ! initialize resampling
  if (.not. parameters%resample) then
     parameters%nsamp = sem_data%ndumps
     allocate(fw_field_res(parameters%nsamp, parameters%nsources))
  else
     allocate(fw_field_res(parameters%nsamp * 2, parameters%nsources))
     call resamp%init(sem_data%ndumps * 2, parameters%nsamp * 2, parameters%nsources)
  endif

  allocate(fw_field(sem_data%ndumps, 1, parameters%nsources))

  ! initialize sources
  allocate(sources(parameters%nsources))
  do i = 1, parameters%nsources
     call sources(i)%read_cmtsolution(fname=trim(parameters%source_files(i)))
  enddo

  ! initialize time traces
  allocate(T(1:parameters%nsamp))
  dt_out = sem_data%dt * sem_data%ndumps / parameters%nsamp

  do i = 1, parameters%nsamp
     T(i) = dt_out * (i - 1)
  end do

  ! initialization done
  iclockold = tick(id=id_init, since=iclockold)

  do i=1, receivers%num_rec
     ! load data from file
     iclockold = tick()
     fw_field = sem_data%load_fw_points_rdbm(sources, receivers%reci_sources(i), &
                                             parameters%component)
     iclockold = tick(id=id_load, since=iclockold)

     ! resample
     if (parameters%resample) then
        iclockold = tick()
        call resamp%resample(taperandzeropad(fw_field(:,1,:), ntaper=0, &
                                             ntimes=sem_data%ndumps * 2), &
                             fw_field_res)
        iclockold = tick(id=id_resamp, since=iclockold)
     else
        fw_field_res = fw_field(:,1,:)
     endif

     ! ouput to file
     iclockold = tick()
     write(fname,'("seis_",I0.3)') i
     open(newunit=lu_seis, file=fname)
     do j = 1, parameters%nsamp
        write(lu_seis,*) T(j), fw_field_res(j,:)
     enddo
     close(lu_seis)
     iclockold = tick(id=id_out, since=iclockold)
  enddo

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
