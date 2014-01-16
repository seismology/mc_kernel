!=========================================================================================
module master_queue
  implicit none
  private
  
  public :: init_work, get_next_task, extract_receive_buffer

  integer, allocatable  :: tasks(:)

contains

!-----------------------------------------------------------------------------------------
subroutine init_work(ntasks)
  
  integer, intent(out)  :: ntasks

  integer               :: itask
  
  ! initialize work
  ntasks = 20
  allocate(tasks(ntasks))
  do itask=1, ntasks
     tasks(itask) = itask
  enddo

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_next_task(itask)
  
  use work_mod
  integer, intent(in)   :: itask

  wt%ntotal_kernel = tasks(itask)
  wt%ndimensions = 10

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine extract_receive_buffer(ioutput, output)
  
  use work_mod
  integer, intent(in)       :: ioutput
  integer, intent(inout)    :: output(:,:)

  output(ioutput,1) = wt%ntotal_kernel
  output(ioutput,2) = wt%ndimensions

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
