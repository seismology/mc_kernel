
!===============
module commpi
!===============
!
! Wrapper routines to invoke the MPI library. 
! This routine is the sole place for parallel interactions. 
! In other words, the routines are called by wrappers that contain an 
! if (nproc>1) statement such that a serial version 
! (i.e. without MPI libraries) of this code *shall* run after merely taking 
! out this module commpi. Consequently, this module is 'public' 
! (every routine here needs to be called by external wrapper)
!
! WARNING: Need to make sure this is consistent for other issues like 
! source and receiver locations in global/local reference.
!
use global_parameters

! in case you have problems with the mpi module, you might try to use the
! include below, in which case you will have to specify the location in the 
! Makefile or copy to the build directory!
use mpi
implicit none
!include 'mpif.h'

integer          :: nproc, mynum, ierror
character(len=4) :: appmynum, appnproc

public 

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------------
subroutine pinit
!
! Start message-passing interface, assigning the total number of processors 
! nproc and each processor with its local number mynum=0,...,nproc-1.
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call MPI_INIT( ierror)
  call MPI_COMM_RANK( MPI_COMM_WORLD, mynum, ierror )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierror )


end subroutine pinit
!=============================================================================

!-----------------------------------------------------------------------------
subroutine pend

  call MPI_FINALIZE(ierror)

end subroutine pend
!=============================================================================

!!-----------------------------------------------------------------------------
!subroutine abort
! 
!  CALL MPI_ABORT(MPI_COMM_WORLD, IERROR)
!
!end subroutine abort
!!=============================================================================

!-----------------------------------------------------------------------------
subroutine barrier
 
  CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

end subroutine barrier
!=============================================================================
end module
