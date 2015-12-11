!=========================================================================================
module commpi
  
  ! Wrapper routines to invoke the MPI library. 
  ! in case you have problems with the mpi module, you might try to use the
  ! include below, in which case you will have to specify the location in the 
  ! Makefile or copy to the build directory!

# ifndef include_mpi
  use mpi
# endif
  use global_parameters, only: sp, dp, master, myrank, nproc, firstslave, testing
  implicit none

# ifdef include_mpi
  include 'mpif.h'
# endif

  private 

  integer, protected    :: MPI_COMM_NODE, MPI_COMM_MASTER_SLAVES, MPI_INFO

  interface pbroadcast
    module procedure pbroadcast_float
    module procedure pbroadcast_float_arr
    module procedure pbroadcast_dble
    module procedure pbroadcast_dble_arr
    module procedure pbroadcast_char
    module procedure pbroadcast_char_arr
    module procedure pbroadcast_int
    module procedure pbroadcast_int_arr
    module procedure pbroadcast_log
  end interface pbroadcast

  public  :: pbroadcast
  public  :: pbroadcast_float, pbroadcast_float_arr
  public  :: pbroadcast_dble, pbroadcast_dble_arr
  public  :: pbroadcast_char, pbroadcast_char_arr
  public  :: pbroadcast_int, pbroadcast_int_arr
  public  :: pbroadcast_log
  public  :: ppinit, pbarrier, pbarrier_node, ppend, pabort, ppsplit
  public  :: MPI_COMM_NODE, MPI_COMM_MASTER_SLAVES, mpi_info

contains

!-----------------------------------------------------------------------------------------
subroutine ppinit
!< Start message-passing interface, assigning the total number of processors 
!! nproc and each processor with its local number mynum=0,...,nproc-1.

  use global_parameters, only  : set_lu_out, set_myrank, set_nproc, & 
                                 set_master, set_firstslave, testing
  integer                     :: ierror, lu_out_loc, myrank_loc, nproc_loc
  character(len=11)           :: fnam
  
  call MPI_INIT( ierror)
  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank_loc, ierror )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc_loc, ierror )

  mpi_info = MPI_INFO_NULL

  call set_myrank(myrank_loc)
  call set_nproc(nproc_loc)

  if ((nproc_loc < 2).and. (.not. testing)) then
    print "('ERROR: Number of CPUs has to be larger one, is:', I5)", nproc_loc
    stop
  end if

  if (myrank == 0) then
      call set_master(.true.)
      call set_firstslave(.false.)
      if (.not.testing) then
        call set_lu_out(6)
      end if
  else
      call set_master(.false.)
      !if (myrank == 1) then
      !  call set_firstslave(.true.)
      !else
      !  call set_firstslave(.false.)
      !end if
      if (.not.testing) then
        write(fnam,"('OUTPUT_', I4.4)") myrank
        open(newunit=lu_out_loc, file=fnam, status='replace')
        call set_lu_out(lu_out_loc)
      end if
  end if

  print *, 'I have rank ', myrank

  
end subroutine ppinit
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Uses MPI_Comm_split to create new communicators spanning only one node. It assumes that 
!! the master uses one node on its own and then assigns NSLAVES_PER_NODE slaves to each
!! node. 
subroutine ppsplit(nslaves_per_node_in)
  use global_parameters, only    : myrank, myrank_node, nproc_node, &
                                   set_myrank_node, set_nproc_node, &
                                   set_myrank_master_slaves, set_nproc_master_slaves, &
                                   set_firstslave, set_ioworker,    &
                                   lu_out, ioworker, dist_io

  integer, intent(in), optional :: nslaves_per_node_in !< How many slaves per node?
                                                       !! If this argument is missing, no slave
                                                       !! gets the IOworker flag and we only create
                                                       !! a MPI communicator for all slaves.
  integer                       :: nslaves_per_node
  integer                       :: mynode, myrank_node_loc, nproc_node_loc, ierror
  integer                       :: master_or_slave !< 0 if this rank is an IO worker (no real slave
                                                   !! 1 if it is the master or a non-IO slave
  integer                       :: myrank_master_slaves_loc, nproc_master_slaves_loc

  if (present(nslaves_per_node_in)) then
    nslaves_per_node = nslaves_per_node_in
  else
    nslaves_per_node = nproc - 1
  end if

  if (myrank==0) then
    ! I am the master, I have a node on my own.
    mynode = 0
  else
    ! I am not the master, I have to share a node
    mynode = int((myrank-1)/nslaves_per_node) + 1
  end if

  ! Create new communicator MPI_COMM_NODE, which contains all slaves on this node
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, mynode, myrank, MPI_COMM_NODE, ierror)
  call MPI_COMM_RANK(MPI_COMM_NODE, myrank_node_loc, ierror )
  call MPI_COMM_SIZE(MPI_COMM_NODE, nproc_node_loc, ierror )

  call set_myrank_node(myrank_node_loc)
  call set_nproc_node(nproc_node_loc)

  print *, 'I am rank ', myrank, ' on node ', mynode, ', myrank_node: ', myrank_node, ', nproc_node', nproc_node
  call pbarrier()

  if ((myrank>0.and.myrank_node==0).and.(dist_io)) then
    write(lu_out,*) ' RANK: ', myrank, ' is a IO-worker for this run'
    call set_ioworker(.true.)
    call set_firstslave(.true.)

  elseif (myrank>0.and.myrank_node==0) then
    call set_firstslave(.true.)
    call set_ioworker(.false.)

  else
    call set_firstslave(.false.)
    call set_ioworker(.false.)

  end if

  ! Create a new communicator MPI_COMM_MASTER_SLAVES, which contains the master and
  ! all slaves that are not IO workers. This is used to send tasks around (the IO
  ! workers should not receive any tasks)
  if (ioworker) then
    master_or_slave = 0
  else
    master_or_slave = 1
  end if

  call MPI_COMM_SPLIT(MPI_COMM_WORLD, master_or_slave, myrank, &
                      MPI_COMM_MASTER_SLAVES, ierror)
  call MPI_COMM_RANK(MPI_COMM_MASTER_SLAVES, myrank_master_slaves_loc, ierror )
  call MPI_COMM_SIZE(MPI_COMM_MASTER_SLAVES, nproc_master_slaves_loc, ierror )
  call set_myrank_master_slaves(myrank_master_slaves_loc)
  call set_nproc_master_slaves(nproc_master_slaves_loc)

end subroutine ppsplit
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_char(input_char, input_proc, comm)

  integer, intent(in)           :: input_proc
  character(*), intent(inout)   :: input_char
  integer, intent(in), optional :: comm
  integer                       :: mpi_comm
  integer                       :: ierror
  logical                       :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (present(comm)) then
    mpi_comm = comm
  else
    mpi_comm = MPI_COMM_WORLD
  end if

  if (isinitialized) then
    call mpi_bcast(input_char, len(input_char), MPI_CHARACTER, input_proc, &
                   mpi_comm, ierror)
    call mpi_barrier(mpi_comm, ierror)
  else
    input_char = input_char
  end if

end subroutine pbroadcast_char
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_char_arr(input_char, input_proc, comm)

  integer, intent(in)          :: input_proc
  character(*), intent(inout)  :: input_char(:)
  integer, intent(in), optional :: comm
  integer                       :: mpi_comm
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (present(comm)) then
    mpi_comm = comm
  else
    mpi_comm = MPI_COMM_WORLD
  end if

  if (isinitialized) then
    call mpi_bcast(input_char, size(input_char), MPI_CHARACTER, input_proc, &
                   mpi_comm, ierror)
    call mpi_barrier(mpi_comm, ierror)
  else
    input_char = input_char
  end if

end subroutine pbroadcast_char_arr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_log(input_log, input_proc, comm)

  integer, intent(in)          :: input_proc
  logical, intent(inout)       :: input_log
  integer, intent(in), optional :: comm
  integer                       :: mpi_comm
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (present(comm)) then
    mpi_comm = comm
  else
    mpi_comm = MPI_COMM_WORLD
  end if

  if (isinitialized) then
    call mpi_bcast(input_log, 1, MPI_LOGICAL, input_proc, mpi_comm, ierror)
    call mpi_barrier(mpi_comm, ierror)
  else
    input_log = input_log
  end if

end subroutine pbroadcast_log
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_int(input_int, input_proc, comm)

  integer, intent(in)          :: input_proc
  integer, intent(inout)       :: input_int
  integer, intent(in), optional :: comm
  integer                       :: mpi_comm
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (present(comm)) then
    mpi_comm = comm
  else
    mpi_comm = MPI_COMM_WORLD
  end if

  if (isinitialized) then
    call mpi_bcast(input_int, 1, MPI_INTEGER, input_proc, mpi_comm, ierror)
    call mpi_barrier(mpi_comm, ierror)
  else 
    input_int = input_int
  end if

end subroutine pbroadcast_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_int_arr(input_int, input_proc, comm)

  integer, intent(in)          :: input_proc
  integer, intent(inout)       :: input_int(:)
  integer, intent(in), optional :: comm
  integer                       :: mpi_comm
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (present(comm)) then
    mpi_comm = comm
  else
    mpi_comm = MPI_COMM_WORLD
  end if

  if (isinitialized) then
    call mpi_bcast(input_int, size(input_int), MPI_INTEGER, input_proc, &
                   mpi_comm, ierror)
    call mpi_barrier(mpi_comm, ierror)
  else 
    input_int = input_int
  end if

end subroutine pbroadcast_int_arr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_dble(input_dble, input_proc, comm)

  integer, intent(in)          :: input_proc
  real(kind=dp), intent(inout) :: input_dble
  integer, intent(in), optional :: comm
  integer                       :: mpi_comm
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (present(comm)) then
    mpi_comm = comm
  else
    mpi_comm = MPI_COMM_WORLD
  end if

  if (isinitialized) then
    call mpi_bcast(input_dble, 1, MPI_DOUBLE_PRECISION, input_proc, &
                   mpi_comm, ierror)
    call mpi_barrier(mpi_comm, ierror)
  else 
    input_dble = input_dble
  end if

end subroutine pbroadcast_dble
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_dble_arr(input_dble, input_proc, comm)

  integer, intent(in)          :: input_proc
  real(kind=dp), intent(inout) :: input_dble(:)
  integer, intent(in), optional :: comm
  integer                       :: mpi_comm
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (present(comm)) then
    mpi_comm = comm
  else
    mpi_comm = MPI_COMM_WORLD
  end if

  if (isinitialized) then
    call mpi_bcast(input_dble, size(input_dble), MPI_DOUBLE_PRECISION, &
                   input_proc, mpi_comm, ierror)
    call mpi_barrier(mpi_comm, ierror)
  else 
    input_dble = input_dble
  end if

end subroutine pbroadcast_dble_arr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_float(input_float, input_proc, comm)

  integer, intent(in)          :: input_proc
  real(kind=sp), intent(inout) :: input_float
  integer, intent(in), optional :: comm
  integer                       :: mpi_comm
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (present(comm)) then
    mpi_comm = comm
  else
    mpi_comm = MPI_COMM_WORLD
  end if

  if (isinitialized) then
    call mpi_bcast(input_float, 1, MPI_REAL, input_proc, &
                   mpi_comm, ierror)
    call mpi_barrier(mpi_comm, ierror)
  else 
    input_float = input_float
  end if

end subroutine pbroadcast_float
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_float_arr(input_float, input_proc, comm)

  integer, intent(in)          :: input_proc
  real(kind=sp), intent(inout) :: input_float(:)
  integer, intent(in), optional :: comm
  integer                       :: mpi_comm
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (present(comm)) then
    mpi_comm = comm
  else
    mpi_comm = MPI_COMM_WORLD
  end if

  if (isinitialized) then
    call mpi_bcast(input_float, size(input_float), MPI_REAL, &
                   input_proc, mpi_comm, ierror)
    call mpi_barrier(mpi_comm, ierror)
  else 
    input_float = input_float
  end if

end subroutine pbroadcast_float_arr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbarrier
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  end if

end subroutine pbarrier
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbarrier_node
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     CALL MPI_BARRIER(MPI_COMM_NODE, ierror)
  end if

end subroutine pbarrier_node
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine ppend
!< Calls MPI_FINALIZE
  integer :: ierror

  call MPI_FINALIZE(ierror)

end subroutine ppend
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pabort(do_traceback)
!< Calls MPI_ABORT
!! By default, does also a traceback. Can be suppressed by the argument do_traceback=.false
# if defined(__INTEL_COMPILER)
  use ifcore, only: tracebackqq ! Ifort helper routines
# endif  
  logical, intent(in), optional :: do_traceback
  integer                       :: ierror
  logical                       :: isinitialized, do_traceback_loc
  character(len=512)            :: msg

  call flush(6)
  print *, 'Processor ', myrank, ' has found an error and aborts computation'

  do_traceback_loc = .true.
  if (present(do_traceback)) do_traceback_loc = do_traceback

  !if (do_traceback_loc) then
#   if defined(__GFORTRAN__)
#   if (__GNUC_MINOR__>=8)
     call backtrace
#   endif
#   endif
#   if defined(__INTEL_COMPILER)
     write(msg,"('Processor ', I5, ' found an error in line:')") myrank
     call flush(6)
     call tracebackqq(string=msg, status=ierror)
#   endif
  !end if

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     call flush(6)
     call MPI_ABORT(MPI_COMM_WORLD, 0, ierror)
  else
     stop
  end if

end subroutine pabort
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
