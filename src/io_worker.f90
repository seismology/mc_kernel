!=========================================================================================
module ioworker_mod

  use global_parameters,           only: sp, dp, long, pi, deg2rad, &
                                         init_random_seed, myrank, lu_out
  implicit none
  private

  integer, parameter                   :: DIETAG = 999

  public                               :: do_ioworker, stop_ioworker

contains

!-----------------------------------------------------------------------------------------
subroutine do_ioworker()
  use global_parameters,           only: DIETAG, id_mpi, id_read_params, id_init_fft, id_int_hetero
  use inversion_mesh,              only: inversion_mesh_data_type
  use readfields,                  only: semdata_type
  use heterogeneities,             only: hetero_type
  use type_parameter,              only: parameter_type
  use fft,                         only: rfft_type, taperandzeropad
  use clocks_mod,                  only: tick

  implicit none
  type(parameter_type)                :: parameters
  type(inversion_mesh_data_type)      :: inv_mesh
  type(semdata_type)                  :: sem_data
  type(rfft_type)                     :: fft_data
  type(hetero_type)                   :: het_model

  integer                             :: ndumps, ntimes, nomega
  integer                             :: ikernel, ierror
  integer(kind=long)                  :: iclockold
  real(kind=dp)                       :: df
  integer                             :: itask
  character(len=255)                  :: fmtstring
  integer                             :: nptperstep

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Read input files for parameters, source and receivers'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  iclockold = tick()

  call parameters%read_parameters()
  call parameters%read_source()
  call parameters%read_receiver()

  nptperstep = parameters%npoints_per_step

  ! Master and slave part ways here for some time. 
  ! Master reads in the inversion mesh, slaves initialize the FFT

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Initialize and open AxiSEM wavefield files'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  call sem_data%set_params(parameters%fwd_dir,     &
                           parameters%bwd_dir,     &
                           parameters%strain_buffer_size, & 
                           parameters%displ_buffer_size, & 
                           parameters%strain_type_fwd,    &
                           parameters%source%depth)

  call sem_data%open_files()
  call sem_data%read_meshes()

  call sem_data%load_seismogram_rdbm(parameters%receiver, parameters%source)

  ndumps = sem_data%ndumps

  iclockold = tick(id=id_read_params, since=iclockold)

  if (parameters%int_over_hetero) then
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Initialize heterogeneity structure'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)

    call het_model%load_het_rtpv(parameters%hetero_file)
    call het_model%build_hetero_kdtree()

    iclockold = tick(id=id_int_hetero, since=iclockold)
  end if

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Initialize FFT'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  call fft_data%init(ndumps, sem_data%get_ndim(), nptperstep, sem_data%dt, &
                     fftw_plan=parameters%fftw_plan)

  ntimes = fft_data%get_ntimes()
  nomega = fft_data%get_nomega()
  df     = fft_data%get_df()
  fmtstring = '(A, I8, A, I8)'
  write(lu_out,fmtstring) '  ntimes: ',  ntimes,     '  , nfreq: ', nomega
  fmtstring = '(A, F8.3, A, F8.3, A)'
  write(lu_out,fmtstring) '  dt:     ', sem_data%dt, ' s, df:    ', df*1000, ' mHz'

  iclockold = tick(id=id_init_fft, since=iclockold)

  ! Master and slave synchronize again

  iclockold = tick()
  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Define filters'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  call parameters%read_filter(nomega, df)

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Define kernels'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  call parameters%read_kernel(sem_data, parameters%filter)

  iclockold = tick(id=id_read_params, since=iclockold)

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') 'Starting my job as an IO worker'
  write(lu_out,'(A)') '***************************************************************'
  flush(lu_out)
  call loop_ioworker(sem_data)

  write(lu_out,'(A)')
  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Close AxiSEM wavefield files'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)
  call sem_data%close_files()

end subroutine do_ioworker
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine loop_ioworker(fields)
  use commpi,            only        : MPI_COMM_NODE
  use global_parameters, only        : nproc_node
  use readfields,        only        : semdata_type, load_single_point_from_file
# ifndef include_mpi
  use mpi
# endif
# ifdef include_mpi
  include 'mpif.h'
# endif

  
  type(semdata_type), intent(in)    :: fields
  integer                           :: mpistatus(MPI_STATUS_SIZE)
  integer                           :: pointids(0:fields%npol, 0:fields%npol)
  integer                           :: rank_sender, field_tag, ierror, myrank_loc
  real(kind=sp)                     :: u(fields%ndumps, fields%npol+1, fields%npol+1, 3)
  integer                           :: t2, t1, ticks_per_sec
  logical                           :: alldone(nproc_node) 
  integer                           :: npts, iproc
  logical                           :: youvegotmail

  alldone(:) = .false.

  call MPI_COMM_RANK( MPI_COMM_NODE, myrank_loc, ierror )
  write(lu_out,*) 'Expecting requests on rank ', myrank_loc, ' on node: ', MPI_COMM_NODE
  flush(lu_out)


  ! Number of points to read
  npts = (fields%npol+1)**2

  receive_io_requests: do 

    probe_slaves: do 
      iproc = iproc + 1
      if (iproc > nproc_node-1) iproc = 1

      ! Non-blocking MPI_probe. Checks whether rank 'iproc' requested data
      call MPI_IProbe(iproc,             &
                      MPI_ANY_TAG,       &
                      MPI_COMM_NODE,     &
                      youvegotmail,      &
                      MPI_STATUS_IGNORE, &
                      ierror)
      if (youvegotmail) exit
    enddo probe_slaves

    ! Receive request from any of the other workers in MPI communicator MPI_COMM_NODE
    call MPI_Recv(pointids,         & ! message buffer
                  npts,             & ! one data item
                  MPI_INTEGER,      & ! data item is an integer
                  iproc,            & ! receive from rank 'iproc'
                  MPI_ANY_TAG,      & ! any type of message
                  MPI_COMM_NODE,    & ! communicator for this node
                  mpistatus,        & ! info about the received message
                  ierror)

    write(lu_out, '(A,I2,A,I3)', advance='no') 'Received tag', field_tag, ' from rank:', rank_sender

    !call system_clock( t1, ticks_per_sec)
    if (ierror.ne.MPI_SUCCESS) then
      print *, 'MPI_Recv error on IO worker: ', ierror
    end if
    rank_sender = mpistatus(MPI_SOURCE)
    field_tag = mpistatus(MPI_TAG)

    select case(field_tag)
    case(1,2,3,4)
      call load_single_point_from_file(fields%fwd(field_tag), pointids, u)
    case(5,6,7,8)
      call load_single_point_from_file(fields%bwd(field_tag-4), pointids, u)
    case(DIETAG) ! the DIETAG
      alldone(rank_sender) = .true.
      write(lu_out, *) 'Received die tag from rank:', rank_sender, alldone
    case default
      print *, 'Invalid field tag: ', field_tag
      stop
    end select
    !call system_clock( t2, ticks_per_sec)
    !write(lu_out,'(A, F8.6, A)'), ', took ', 1.e3/real(ticks_per_sec)*(t2-t1), ' ms to answer'

    if (field_tag.ne.DIETAG) then
      ! Send the same worker the loaded time series (blocking)
      call MPI_RSend(u,                     & ! message buffer
                     3*fields%ndumps*npts,  & ! three dimensions per time step
                     MPI_REAL,              & ! data item is a single-precision float
                     rank_sender,           & ! to who we just received from
                     field_tag,             & ! user chosen message tag
                     MPI_COMM_NODE,         & ! communicator for this node
                     ierror)

       if (ierror.ne.MPI_SUCCESS) then
         print *, 'MPI_RSend error on IO worker: ', ierror
       end if

    else ! This worker is finished, check whether all are
      if (all(alldone)) exit
    end if

  end do receive_io_requests

end subroutine loop_ioworker
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> This routine is called by all slaves after they are done to send the IO worker the 
!! DIE tag (-1). If he has received the DIE tag once by all slaves, he shuts himself down.
subroutine stop_ioworker()
  use commpi, only   : MPI_COMM_NODE
# ifndef include_mpi
  use mpi
# endif
# ifdef include_mpi
  include 'mpif.h'
# endif
  integer           :: ierror
  integer           :: tmp(1)

  tmp = 1

  ! Tell IO worker to stop waiting for tasks (unless we're the IO worker)
  call MPI_SSend(tmp(1),                & ! message buffer
                 1,                     & ! One point ID
                 MPI_INTEGER,           & ! data item is a integer
                 0,                     & ! to rank zero, the io-worker
                 DIETAG,                & ! The DIE tag
                 MPI_COMM_NODE,         & ! Node communicator
                 ierror)

end subroutine stop_ioworker
!-----------------------------------------------------------------------------------------

end module ioworker_mod
!=========================================================================================
