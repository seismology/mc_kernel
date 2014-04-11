!=========================================================================================
program rdbm

  use readfields, only : semdata_type
  use commpi
  use global_parameters
  use source_class
  use resampling
  use type_parameter
  use receivers_rdbm

  implicit none

  type(semdata_type)                  :: sem_data
  !type(src_param_type)                :: reci_source
  type(parameter_type)                :: parameters
  type(receivers_rdbm_type)           :: receivers

  character(len=512)                  :: bwd_dir
  character(len=4)                    :: model_param
  real(kind=dp), allocatable          :: source_coordinates(:,:)
  real(kind=dp)                       :: x, y, r, th
  real(kind=dp), allocatable          :: fw_field(:,:,:)
  integer                             :: i, j

  type(resampling_type)               :: resamp
  real(kind=dp), allocatable          :: fw_field_res(:,:)
  integer                             :: ntraces

  real(kind=dp), dimension(:), allocatable      :: T
  real(kind=dp)                                 :: dt_out

  call parameters%read_parameters('inparam_basic')

  call receivers%read_stations_file()

  ntraces = 1

  write(*,*) '***************************************************************'
  write(*,*) ' Initialize and open AxiSEM wavefield files'
  write(*,*) '***************************************************************'

  bwd_dir = ''

  model_param = 'vp'
  call sem_data%set_params(parameters%sim_dir, bwd_dir, parameters%buffer_size, model_param)
  call sem_data%open_files()
  call sem_data%read_meshes()
  call sem_data%build_kdtree()

  !call reci_source%init(90d0, 0d0, (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0 /))


  allocate(fw_field(sem_data%ndumps, 1, ntraces))
  allocate(fw_field_res(parameters%nsamp, ntraces))
  allocate(source_coordinates(3, ntraces))

  r = 5000.
  th = 0.

  do i = 1, ntraces
     x = sin(th / 180 * pi) * r
     y = cos(th / 180 * pi) * r
     source_coordinates(:,i)  = (/0d0, x, y/)
     r = r - 5
  enddo

  allocate(T(1:parameters%nsamp))

  dt_out = sem_data%dt * sem_data%ndumps / parameters%nsamp

  do i = 1, parameters%nsamp
     T(i) = dt_out * (i - 1)
  end do

  call resamp%init(sem_data%ndumps, parameters%nsamp, ntraces)

  do i=1, receivers%num_rec
     fw_field = sem_data%load_fw_points_rdbm(source_coordinates, receivers%reci_sources(i))

     call resamp%resample(fw_field(:,1,:), fw_field_res)

     !do i = 1, sem_data%ndumps
     !   write(111,*) real(i-1) / (sem_data%ndumps -1), fw_field(i,1,:)
     !enddo

     do j = 1, parameters%nsamp
        write(110+i,*) T(j), fw_field_res(j,:)
     enddo
  enddo

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

  id_fft         = clock_id('FFT routines')
  id_fwd         = clock_id('Reading fwd field')
  id_netcdf      = clock_id(' - NetCDF routines')
  id_rotate      = clock_id(' - Rotate fields')
  id_buffer      = clock_id(' - Buffer routines')
  id_bwd         = clock_id('Reading bwd field')
  id_mc          = clock_id('Monte Carlo routines')
  id_filter_conv = clock_id('Filtering and convolution')
  id_inv_mesh    = clock_id('Inversion mesh routines')
  id_kernel      = clock_id('Kernel routines')
  id_init        = clock_id('Initialization per task')
  id_mpi         = clock_id('MPI communication with Master')

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
