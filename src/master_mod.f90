!=========================================================================================
module master_module
  implicit none
  private
  integer, parameter         :: sp = selected_real_kind(6, 37)
  integer, parameter         :: dp = selected_real_kind(15, 307)
  public :: do_master
contains

!-----------------------------------------------------------------------------------------
subroutine do_master(mpi_communicator)
# ifndef include_mpi
  use mpi
# endif
  use master_slave_parameters
  use work_type_mod
  use master_queue

#ifdef include_mpi
  include 'mpif.h'
#endif

  integer               :: mpi_communicator  ! The MPI communicator that engulfs
                                             ! the master and all the slaves
  integer               :: nslaves, rank, ierror
  integer, allocatable  :: output(:,:), sendrequest(:), work_done(:)
  integer               :: mpistatus(MPI_STATUS_SIZE)
  integer               :: itask, ntasks, ioutput, iclock, iclock_ref, ticks_per_sec
  integer               :: itask_result
  real(kind=dp)         :: time
  character(len=64)     :: fmtstring

  call init_queue(ntasks)

  allocate(output(ntasks,2))
  output = -1
  
  ! Find out how many processes there are in the default communicator
  call MPI_COMM_SIZE(mpi_communicator, nslaves, ierror)
  nslaves = nslaves - 1 ! the master does not work
  allocate(sendrequest(nslaves))

  ! Format string for plotting status of slaves
  allocate(work_done(nslaves))
  work_done = 0
  write(fmtstring,"('(',I5,'('' '', I5),'' | '' F6.2,''%'', F10.1, ''sec'')')") nslaves + 1

  if (nslaves > ntasks) then
    write(6,*) 'ERROR: more slaves than tasks'
    stop
  elseif (nslaves < 1) then
    write(6,*) 'ERROR: need at least 1 slave'
    stop
  endif

  write(6, '(A,I8,A,I8)') 'Master initialized, ntasks :', ntasks, ', nslaves:', nslaves

  call system_clock(count=iclock_ref, count_rate=ticks_per_sec)

  itask = 0
  ! Seed the slaves; send one unit of work to each slave.
  do rank=1, nslaves

    ! Find the next item of work to do
    itask = itask + 1

    ! fill sendbuffer
    call get_next_task(itask)

    ! Send it to each rank (blocking)
    call MPI_Send(wt,                & ! message buffer
                  1,                 & ! one data item
                  wt%mpitype,        & ! data item is an integer
                  rank,              & ! destination process rank
                  WORKTAG,           & ! user chosen message tag
                  mpi_communicator,  & ! default communicator
                  sendrequest(rank), &
                  ierror)
  enddo

  ! Loop over tasks until there is no more work to be done
  ioutput = 0
  do itask=nslaves+1, ntasks

    ! Receive results from any !sic! slave (blocking!)
    ioutput = ioutput + 1
    call MPI_Recv(wt,               & ! message buffer
                  1,                & ! one data item
                  wt%mpitype,       & ! data item is an integer
                  MPI_ANY_SOURCE,   & ! receive from any sender
                  MPI_ANY_TAG,      & ! any type of message
                  mpi_communicator, & ! default communicator
                  mpistatus,        & ! info about the received message
                  ierror)

    itask_result = wt%itask              
    call extract_receive_buffer(itask_result, mpistatus(MPI_SOURCE))
    
    ! fill sendbuffer
    call get_next_task(itask)

    ! Send the same slave some more work to do (blocking)
    call MPI_Send(wt,               & ! message buffer
                  1,                & ! one data item
                  wt%mpitype,       & ! data item is an integer
                  mpistatus(MPI_SOURCE), & ! to who we just received from
                  WORKTAG,          & ! user chosen message tag
                  mpi_communicator, & ! default communicator
                  sendrequest(mpistatus(MPI_SOURCE)), &
                  ierror)

    ! Some more stuff before waiting for next result
    ! 1. Plot status of slaves to stdout
    ! 2. Save intermediate results to NetCDF file

    ! Plot status of slaves
    work_done(mpistatus(MPI_SOURCE)) = work_done(mpistatus(MPI_SOURCE)) + 1          
    
    !Get time elapsed since start
    call system_clock(count=iclock)
    time = real(iclock-iclock_ref, kind=dp) / real(ticks_per_sec, kind=dp)

    write(6,fmtstring) work_done, sum(work_done), real(sum(work_done)) / real(ntasks) * 100., time

    ! Save intermediate results to NetCDF file
    call dump_intermediate(itask_result)

  enddo

  ! There's no more work to be distributed, so receive all the outstanding results 
  ! from the slaves (blocking, so when this loop is finished, work is done and
  ! results received in the buffer!).
  do rank=1, nslaves
    ioutput = ioutput + 1
    call MPI_Recv(wt,              & ! message buffer
                  1,               & ! one data item
                  wt%mpitype,      & ! data item is an integer
                  MPI_ANY_SOURCE,  & ! receive from any sender
                  MPI_ANY_TAG,     & ! any type of message
                  mpi_communicator,& ! default communicator
                  mpistatus,       & ! info about the received message
                  ierror)
    
    call extract_receive_buffer(wt%itask, mpistatus(MPI_SOURCE))

    ! Send this processor the die tag
    call MPI_Send(0,               & !
                  0,               & ! empty message
                  MPI_INTEGER,     & !
                  mpistatus(MPI_SOURCE), & ! to who we just received from
                  DIETAG,          & ! the tag conatains the actual information
                  mpi_communicator,& ! default communicator
                  sendrequest(mpistatus(MPI_SOURCE)), &
                  ierror)
    
    ! Plot status of slaves
    work_done(mpistatus(MPI_SOURCE)) = work_done(mpistatus(MPI_SOURCE)) + 1          
    
    !Get time elapsed since start
    call system_clock(count=iclock)
    time = real(iclock-iclock_ref, kind=dp) / real(ticks_per_sec, kind=dp)

    write(6, fmtstring) work_done, sum(work_done), real(sum(work_done)) / real(ntasks) * 100., time

  enddo

  call finalize()

  call delete_intermediate()

end subroutine do_master
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
