!===============
module commpi
!===============
  
  ! Wrapper routines to invoke the MPI library. 
  ! This routine is the sole place for parallel interactions. 
  ! in case you have problems with the mpi module, you might try to use the
  ! include below, in which case you will have to specify the location in the 
  ! Makefile or copy to the build directory!
  use mpi
  use global_parameters, only: dp
  implicit none

  integer :: myrank, nproc
  logical :: master

  public  :: pbroadcast_dble, pbroadcast_char, pbroadcast_log
  public  :: pbroadcast_int_arr, pbroadcast_int, ppinit

contains

!----------------------------------------------------------------------------
subroutine ppinit
!< Start message-passing interface, assigning the total number of processors 
!! nproc and each processor with its local number mynum=0,...,nproc-1.

  integer :: ierror
  
  call MPI_INIT( ierror)
  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierror )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierror )

  if (myrank == 0) then
      master = .true.
  else
      master = .false.
  end if
  
end subroutine ppinit
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_char(input_char,input_proc)

  integer, intent(in)           :: input_proc
  character(*), intent(inout)   :: input_char
  integer                       :: ierror


  call mpi_bcast(input_char, len(input_char), MPI_CHARACTER, input_proc, &
                 MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)

end subroutine pbroadcast_char
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_log(input_log,input_proc)

  integer, intent(in)    :: input_proc
  logical, intent(inout) :: input_log
  integer                :: ierror

  call mpi_bcast(input_log, 1, MPI_LOGICAL, input_proc, MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)

end subroutine pbroadcast_log
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_int(input_int,input_proc)

  integer, intent(in)    :: input_proc
  integer, intent(inout) :: input_int
  integer                :: ierror

  call mpi_bcast(input_int, 1, MPI_INTEGER, input_proc, MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)

end subroutine pbroadcast_int
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_int_arr(input_int, input_proc)

  integer, intent(in)    :: input_proc
  integer, intent(inout) :: input_int(:)
  integer                :: ierror

  call mpi_bcast(input_int, size(input_int), MPI_INTEGER, input_proc, &
                 MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)

end subroutine pbroadcast_int_arr
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_dble(input_dble,input_proc)

  integer, intent(in)          :: input_proc
  real(kind=dp), intent(inout) :: input_dble
  integer                      :: ierror

  call mpi_bcast(input_dble, 1, MPI_DOUBLE_PRECISION, input_proc, &
                 MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)

end subroutine pbroadcast_dble
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pbroadcast_dble_arr(input_dble,input_proc)

  integer, intent(in)          :: input_proc
  real(kind=dp), intent(inout) :: input_dble(:)
  integer                      :: ierror

  call mpi_bcast(input_dble, size(input_dble), MPI_DOUBLE_PRECISION, &
                 input_proc, MPI_COMM_WORLD, ierror)
  call mpi_barrier(MPI_COMM_WORLD, ierror)

end subroutine pbroadcast_dble_arr
!=============================================================================

end module
