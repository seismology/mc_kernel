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
  public :: wt          ! essentially implementing the singleton pattern by
                        ! providing a single variable of the private type work_type

  type work_type
     sequence           ! force the derived type to be stored contiguously
     ! can't make the init function a member of work_type, because then it is
     ! impossible to use the seqeuence keyword
     ! maybe can get around this by defining one block for each variable in the
     ! derived type and dropping the sequence statement
     integer :: ntotal_kernel
     integer :: ndimensions
     integer :: nvertices
     integer :: mpitype
     real(kind=dp), allocatable :: vertices(:,:)
     real(kind=dp), allocatable :: kernel_values(:,:)
  end type

  type(work_type) :: wt

contains

!-----------------------------------------------------------------------------------------
subroutine init_work_type(nkern, ndim, nverts)
  use mpi
  
  integer, intent(in)   :: ndim, nverts, nkern
  integer               :: ierr, i
  integer, allocatable  :: oldtypes(:), blocklengths(:)
  integer(kind=MPI_ADDRESS_KIND), allocatable  :: offsets(:)
  integer, parameter    :: nblocks = 3

  wt%ntotal_kernel = nkern
  wt%nvertices = nverts
  wt%ndimensions = ndim

  allocate(wt%vertices(1:wt%ndimensions, 1:wt%nvertices))
  allocate(wt%kernel_values(wt%ntotal_kernel, wt%nvertices))
  wt%vertices = 0
  wt%kernel_values = 0

  ! define blocks for the mpi type. NB: it seems to be necessary to define one
  ! block per array, otherwise having segfaults.
  allocate(oldtypes(nblocks))
  allocate(blocklengths(nblocks))
  allocate(offsets(nblocks))

  blocklengths(1) = 4
  blocklengths(2) = wt%ndimensions * wt%nvertices
  blocklengths(3) = wt%ntotal_kernel * wt%nvertices

  oldtypes(1) = MPI_INTEGER
  oldtypes(2) = MPI_DOUBLE_PRECISION
  oldtypes(3) = MPI_DOUBLE_PRECISION

  ! find memory offsets, more stable then computing with MPI_TYPE_EXTEND
  call MPI_GET_ADDRESS(wt%ntotal_kernel, offsets(1), ierr)
  call MPI_GET_ADDRESS(wt%vertices,      offsets(2), ierr)
  call MPI_GET_ADDRESS(wt%kernel_values, offsets(3), ierr)

  ! make offsets relative
  do i=2, size(offsets)
      offsets(i) = offsets(i) - offsets(1)
  end do
  offsets(1) = 0

  call MPI_TYPE_CREATE_STRUCT(nblocks, blocklengths, offsets, oldtypes, wt%mpitype, ierr)
  call MPI_TYPE_COMMIT(wt%mpitype, ierr)

end subroutine init_work_type
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
  use work_type_mod

  integer               :: nslaves, rank, ierror
  integer, allocatable  :: tasks(:), output(:,:), sendrequest(:)
  integer               :: mpistatus(MPI_STATUS_SIZE)
  integer               :: itask, ntasks, ioutput

  ! initialize work
  ntasks = 20
  allocate(tasks(ntasks))
  do itask=1, ntasks
     tasks(itask) = itask
  enddo

  allocate(output(ntasks,2))
  output = -1
  
  ! Find out how many processes there are in the default communicator
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nslaves, ierror)
  nslaves = nslaves - 1 ! the master does not work
  allocate(sendrequest(nslaves))

  if (nslaves > ntasks) then
     write(6,*) 'ERROR: more slaves then tasks'
     stop
  elseif (nslaves < 1) then
     write(6,*) 'ERROR: need at least 1 slave'
     stop
  endif

  itask = 0
  ! Seed the slaves; send one unit of work to each slave.
  do rank=1, nslaves

     ! Find the next item of work to do
     itask = itask + 1

     ! fill sendbuffer
     wt%ntotal_kernel = tasks(itask)
     wt%ndimensions = 10

     ! Send it to each rank (blocking)
     call MPI_Send(wt,                & ! message buffer
                   1,                 & ! one data item
                   wt%mpitype,        & ! data item is an integer
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
     call MPI_Recv(wt,               & ! message buffer
                   1,                & ! one data item
                   wt%mpitype,       & ! data item is an integer
                   MPI_ANY_SOURCE,   & ! receive from any sender
                   MPI_ANY_TAG,      & ! any type of message
                   MPI_COMM_WORLD,   & ! default communicator
                   mpistatus,        & ! info about the received message
                   ierror)

     ! extract from receive buffer
     output(ioutput,1) = wt%ntotal_kernel
     output(ioutput,2) = wt%ndimensions

     
     ! fill sendbuffer
     wt%ntotal_kernel = tasks(itask)
     wt%ndimensions = 10

     ! Send the same slave some more work to do (blocking)
     call MPI_Send(wt,               & ! message buffer
                   1,                & ! one data item
                   wt%mpitype,       & ! data item is an integer
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
     call MPI_Recv(wt,              & ! message buffer
                   1,               & ! one data item
                   wt%mpitype,      & ! data item is an integer
                   MPI_ANY_SOURCE,  & ! receive from any sender
                   MPI_ANY_TAG,     & ! any type of message
                   MPI_COMM_WORLD,  & ! default communicator
                   mpistatus,       & ! info about the received message
                   ierror)
     
     ! extract from receive buffer
     output(ioutput,1) = wt%ntotal_kernel
     output(ioutput,2) = wt%ndimensions
  enddo

  do itask=1, ntasks
     write(6,*) output(itask,2), output(itask,1)
  enddo

  ! Tell all the slaves to exit by sending an empty message with the DIETAG.
  do rank=1, nslaves
     call MPI_Send(0,               & !
                   0,               & ! empty message
                   MPI_INTEGER,     & !
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
  use work_type_mod
  
  integer   :: task
  integer   :: output
  integer   :: mpistatus(MPI_STATUS_SIZE), ierror

  do while (.true.)
    ! Receive a message from the master
    call MPI_Recv(wt, 1, wt%mpitype, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, ierror)

    ! Check the tag of the received message. If no more work to do, exit loop
    ! and return to main programm
    if (mpistatus(MPI_TAG) == DIETAG) exit

    task = wt%ntotal_kernel

    ! Do the work
    output = work(task)
    
    wt%ntotal_kernel = output
    wt%ndimensions = task
    wt%vertices = task

    ! Send the result back
    call MPI_Send(wt, 1, wt%mpitype, 0, 0, MPI_COMM_WORLD, ierror)
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

!!=========================================================================================
!program master_slave
!  use mpi
!  use master_module
!  use slave_module
!  use work_type_mod
!
!  implicit none
!  integer               :: myrank, ierror
!  integer               :: nkern
!  
!  nkern = 10
!
!  call MPI_INIT(ierror)
!  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
!
!  call init_work_type(nkern, 3, 4)
!
!  if (myrank == 0) then
!     print *, 'MASTER'
!     call master()
!  else
!     print *, 'SLAVE', myrank
!     call slave()
!  endif
!
!  call MPI_FINALIZE(ierror)
!end
!!=========================================================================================
