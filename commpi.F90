!===============
module commpi
!===============
  
  ! Wrapper routines to invoke the MPI library. 
  ! This routine is the sole place for parallel interactions. 
  ! in case you have problems with the mpi module, you might try to use the
  ! include below, in which case you will have to specify the location in the 
  ! Makefile or copy to the build directory!
  use mpi
  use global_parameters, only: dp, lu_out, master, myrank, nproc, firstslave
  implicit none

  private 

  public  :: pbroadcast_dble, pbroadcast_dble_arr, pbroadcast_char, pbroadcast_log
  public  :: pbroadcast_int_arr, pbroadcast_int, ppinit, pbarrier, ppend, pabort

contains

!----------------------------------------------------------------------------
subroutine ppinit
!< Start message-passing interface, assigning the total number of processors 
!! nproc and each processor with its local number mynum=0,...,nproc-1.

  integer            :: ierror
  character(len=10)  :: fnam
  
  call MPI_INIT( ierror)
  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierror )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierror )

  firstslave = .false. 

  if (myrank == 0) then
      master = .true.
      lu_out = 6
  else
      master = .false.
      if (myrank == 1) firstslave = .true.
      write(fnam,"('OUTPUT_', I3.3)") myrank
      open(newunit=lu_out, file=fnam, status='replace')
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

!-----------------------------------------------------------------------------
subroutine pbarrier
  integer :: ierror
 
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

end subroutine pbarrier
!=============================================================================

!-----------------------------------------------------------------------------
subroutine ppend
!< Calls MPI_FINALIZE
  integer :: ierror

  call MPI_FINALIZE(ierror)

end subroutine ppend
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pabort
!< Calls MPI_ABORT
  integer :: ierror

  print *, 'Processor ', myrank, ' has found an error and aborts computation'

#if defined(__GFORTRAN__)
#if (__GNUC_MINOR__>=8)
  call backtrace
#endif
#endif
#if defined(__INTEL_COMPILER)
  call tracebackqq()
#endif

  call MPI_ABORT(MPI_COMM_WORLD, 0, ierror)

  stop

end subroutine pabort
!=============================================================================
end module
