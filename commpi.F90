!=========================================================================================
module commpi
  
  ! Wrapper routines to invoke the MPI library. 
  ! This routine is the sole place for parallel interactions. 
  ! in case you have problems with the mpi module, you might try to use the
  ! include below, in which case you will have to specify the location in the 
  ! Makefile or copy to the build directory!

# ifndef include_mpi
  use mpi
# endif
  use global_parameters, only: dp, master, myrank, nproc, firstslave
  implicit none

# ifdef include_mpi
  include 'mpif.h'
# endif

  private 

  public  :: pbroadcast_dble, pbroadcast_dble_arr
  public  :: pbroadcast_char, pbroadcast_char_arr
  public  :: pbroadcast_int, pbroadcast_int_arr
  public  :: ppinit, pbarrier, ppend, pabort
  public  :: pbroadcast_log

contains

!-----------------------------------------------------------------------------------------
subroutine ppinit
!< Start message-passing interface, assigning the total number of processors 
!! nproc and each processor with its local number mynum=0,...,nproc-1.

  use global_parameters, only  : set_lu_out, set_myrank, set_nproc
  integer                     :: ierror, lu_out_loc, myrank_loc, nproc_loc
  character(len=11)           :: fnam
  
  call MPI_INIT( ierror)
  call MPI_COMM_RANK( MPI_COMM_WORLD, myrank_loc, ierror )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc_loc, ierror )

  call set_myrank(myrank_loc)
  call set_nproc(nproc_loc)

  firstslave = .false. 

  if (myrank == 0) then
      master = .true.
      call set_lu_out(6)
  else
      master = .false.
      if (myrank == 1) firstslave = .true.
      write(fnam,"('OUTPUT_', I4.4)") myrank
      open(newunit=lu_out_loc, file=fnam, status='replace')
      call set_lu_out(lu_out_loc)
  end if

  
end subroutine ppinit
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_char(input_char,input_proc)

  integer, intent(in)          :: input_proc
  character(*), intent(inout)  :: input_char
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     call mpi_bcast(input_char, len(input_char), MPI_CHARACTER, input_proc, &
                    MPI_COMM_WORLD, ierror)
     call mpi_barrier(MPI_COMM_WORLD, ierror)
  else
     input_char = input_char
  end if

end subroutine pbroadcast_char
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_char_arr(input_char,input_proc)

  integer, intent(in)          :: input_proc
  character(*), intent(inout)  :: input_char(:)
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     call mpi_bcast(input_char, size(input_char), MPI_CHARACTER, input_proc, &
                    MPI_COMM_WORLD, ierror)
     call mpi_barrier(MPI_COMM_WORLD, ierror)
  else
     input_char = input_char
  end if

end subroutine pbroadcast_char_arr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_log(input_log,input_proc)

  integer, intent(in)          :: input_proc
  logical, intent(inout)       :: input_log
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     call mpi_bcast(input_log, 1, MPI_LOGICAL, input_proc, MPI_COMM_WORLD, ierror)
     call mpi_barrier(MPI_COMM_WORLD, ierror)
  else
     input_log = input_log
  end if

end subroutine pbroadcast_log
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_int(input_int,input_proc)

  integer, intent(in)          :: input_proc
  integer, intent(inout)       :: input_int
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     call mpi_bcast(input_int, 1, MPI_INTEGER, input_proc, MPI_COMM_WORLD, ierror)
     call mpi_barrier(MPI_COMM_WORLD, ierror)
  else 
     input_int = input_int
  end if

end subroutine pbroadcast_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_int_arr(input_int, input_proc)

  integer, intent(in)          :: input_proc
  integer, intent(inout)       :: input_int(:)
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     call mpi_bcast(input_int, size(input_int), MPI_INTEGER, input_proc, &
                    MPI_COMM_WORLD, ierror)
     call mpi_barrier(MPI_COMM_WORLD, ierror)
  else 
     input_int = input_int
  end if

end subroutine pbroadcast_int_arr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_dble(input_dble,input_proc)

  integer, intent(in)          :: input_proc
  real(kind=dp), intent(inout) :: input_dble
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     call mpi_bcast(input_dble, 1, MPI_DOUBLE_PRECISION, input_proc, &
                    MPI_COMM_WORLD, ierror)
     call mpi_barrier(MPI_COMM_WORLD, ierror)
  else 
     input_dble = input_dble
  end if

end subroutine pbroadcast_dble
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbroadcast_dble_arr(input_dble,input_proc)

  integer, intent(in)          :: input_proc
  real(kind=dp), intent(inout) :: input_dble(:)
  integer                      :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     call mpi_bcast(input_dble, size(input_dble), MPI_DOUBLE_PRECISION, &
                    input_proc, MPI_COMM_WORLD, ierror)
     call mpi_barrier(MPI_COMM_WORLD, ierror)
  else 
     input_dble = input_dble
  end if

end subroutine pbroadcast_dble_arr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pbarrier
  integer :: ierror
  logical                      :: isinitialized

  call MPI_Initialized(isinitialized, ierror)

  if (isinitialized) then
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  end if

end subroutine pbarrier
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
