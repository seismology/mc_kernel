!=========================================================================================
module parameters
  implicit none
  public
  integer, parameter :: WORKTAG = 1
  integer, parameter :: DIETAG  = 2
end module
!=========================================================================================

!=========================================================================================
module work_type_mod
  use global_parameters
  implicit none
  private

  public :: init_work_type
  public :: work_type
  public :: ntotal_kernel
  public :: work_mpitype

  integer, protected :: ntotal_kernel, ndimensions, nvertices
  integer, protected :: work_mpitype

  type work_type
     sequence           ! force the derived type to be stored contiguously
     ! can't make the init function a member of work_type, because then it is
     ! impossible to use the seqeuence keyword
     real(kind=dp), allocatable :: vertices(:,:)
     real(kind=dp), allocatable :: kernel_values(:,:)
  end type

contains

!-----------------------------------------------------------------------------------------
function init_work_type(nkern, ndim, nverts)
  use mpi
  
  type(work_type)       :: init_work_type
  integer, intent(in)   :: ndim, nverts, nkern
  integer               :: ierr
  integer, allocatable  :: oldtypes(:), blocklengths(:), offsets(:)
  integer, parameter    :: nblocks = 1

  ntotal_kernel = nkern
  nvertices = nverts
  ndimensions = ndim

  allocate(init_work_type%vertices(ndimensions, nvertices))
  allocate(init_work_type%kernel_values(ntotal_kernel, nvertices))
  init_work_type%vertices = 0
  init_work_type%kernel_values = 0
  
  allocate(oldtypes(nblocks))
  allocate(blocklengths(nblocks))
  allocate(offsets(nblocks))

  blocklengths(1) = ndimensions * nvertices + ntotal_kernel * nvertices
  oldtypes(1) = MPI_DOUBLE_PRECISION
  offsets(1) = 0

  call MPI_TYPE_STRUCT(nblocks, blocklengths, offsets, oldtypes, work_mpitype, ierr)
  call MPI_TYPE_COMMIT(work_mpitype, ierr)

end function init_work_type
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================

!=========================================================================================
module master_module
  implicit none
  private
  public :: master
contains

!-----------------------------------------------------------------------------------------
subroutine master()
  use mpi
  use parameters

  integer               :: nslaves, rank, ierror
  integer, allocatable  :: tasks(:), output(:), sendrequest(:)
  integer               :: mpistatus(MPI_STATUS_SIZE)
  integer               :: itask, ntasks, ioutput

  ! initialize work
  ntasks = 10
  allocate(tasks(ntasks))
  do itask=1, ntasks
     tasks(itask) = itask
  enddo

  allocate(output(ntasks))
  output = -1
  
  ! Find out how many processes there are in the default communicator
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nslaves, ierror)
  nslaves = nslaves - 1 ! the master does not work
  allocate(sendrequest(nslaves))

  if (nslaves > ntasks) then
     write(6,*) 'ERROR: more slaves then tasks'
     stop
  endif

  itask = 0
  ! Seed the slaves; send one unit of work to each slave.
  do rank=1, nslaves

    ! Find the next item of work to do
    itask = itask + 1

    ! Send it to each rank (nonblocking)
    call MPI_ISend(tasks(itask),      & ! message buffer
                   1,                 & ! one data item
                   MPI_INT,           & ! data item is an integer
                   rank,              & ! destination process rank
                   WORKTAG,           & ! user chosen message tag
                   MPI_COMM_WORLD,    & ! default communicator
                   sendrequest(rank), &
                   ierror)
  enddo

  ! Loop over tasks until there is no more work to be done
  ioutput = 0
  do itask=nslaves+1, ntasks

    ! Receive results from any !sic! slave (blocking!)
    ioutput = ioutput + 1
    call MPI_Recv(output(ioutput),  & ! message buffer
                  1,                & ! one data item
                  MPI_INT,          & ! data item is an integer
                  MPI_ANY_SOURCE,   & ! receive from any sender
                  MPI_ANY_TAG,      & ! any type of message
                  MPI_COMM_WORLD,   & ! default communicator
                  mpistatus,        & ! info about the received message
                  ierror)

    ! Send the same slave some more work to do (nonblocking)
    call MPI_ISend(tasks(itask),     & ! message buffer
                   1,                & ! one data item
                   MPI_INT,          & ! data item is an integer
                   mpistatus(MPI_SOURCE), & ! to who we just received from
                   WORKTAG,          & ! user chosen message tag
                   MPI_COMM_WORLD,   & ! default communicator
                   sendrequest(mpistatus(MPI_SOURCE)), &
                   ierror)
  enddo

  ! There's no more work to be distributed, so receive all the outstanding results 
  ! from the slaves (blocking, so when this loop is finished, work is done and
  ! results received in the buffer!).
  do rank=1, nslaves
    ioutput = ioutput + 1
    call MPI_Recv(output(ioutput), & ! message buffer
                  1,               & ! one data item
                  MPI_INT,         & ! data item is an integer
                  MPI_ANY_SOURCE,  & ! receive from any sender
                  MPI_ANY_TAG,     & ! any type of message
                  MPI_COMM_WORLD,  & ! default communicator
                  mpistatus,       & ! info about the received message
                  ierror)
  enddo

  write(6,*) output

  ! Tell all the slaves to exit by sending an empty message with the DIETAG.
  do rank=1, nslaves
    call MPI_Send(0,               & !
                  0,               & ! empty message
                  MPI_INT,         & !
                  rank,            & ! destination
                  DIETAG,          & ! the tag conatains the actual information
                  MPI_COMM_WORLD,  & ! default communicator
                  sendrequest(rank), &
                  ierror)
  enddo

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================

!=========================================================================================
module slave_module
  implicit none
  private
  public :: slave
contains

!-----------------------------------------------------------------------------------------
subroutine slave()
  use mpi
  use parameters
  integer   :: task
  integer   :: output
  integer   :: mpistatus(MPI_STATUS_SIZE), ierror

  do while (.true.)
    ! Receive a message from the master
    call MPI_Recv(task, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, ierror)

    ! Check the tag of the received message. If no more work to do, exit loop
    ! and return to main programm
    if (mpistatus(MPI_TAG) == DIETAG) exit

    ! Do the work
    output = work(task)

    ! Send the result back
    call MPI_Send(output, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, ierror)
  enddo
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function work(task)
  integer, intent(in) :: task
  
  work = task ** 2
end function
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================

!=========================================================================================
program master_slave
  use mpi
  use master_module
  use slave_module
  use work_type_mod

  implicit none
  integer               :: myrank, ierror
  integer               :: nkern = 1
  type(work_type)       :: wt

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)

  wt = init_work_type(nkern, 3, 4)
  write(6,*) 'mpi_worktype = ', work_mpitype

  if (myrank == 0) then
     print *, 'MASTER'
     call master()
  else
     print *, 'SLAVE', myrank
     call slave()
  endif

  call MPI_FINALIZE(ierror)
end
!=========================================================================================
