!=========================================================================================
module commpi
  ! a cleaned up version of the commpi module to compile without mpi, only provides
  ! functionality of pabort function to produce a backtrace
  
  implicit none

  private 

  public  :: pabort

contains

!-----------------------------------------------------------------------------------------
subroutine pabort

#if defined(__GFORTRAN__)
#if (__GNUC_MINOR__>=8)
  call backtrace
#endif
#endif
#if defined(__INTEL_COMPILER)
  call tracebackqq()
#endif

  stop

end subroutine pabort
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
