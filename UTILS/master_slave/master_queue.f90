!=========================================================================================
module master_queue
  implicit none
  private
  
  public :: init_queue, get_next_task, extract_receive_buffer

  integer, allocatable  :: tasks(:)

contains

!-----------------------------------------------------------------------------------------
subroutine init_queue(ntasks)
! anything that should be done before starting the loop over the work. for now,
! the number of tasks is fixed here

  integer, intent(out)  :: ntasks
  integer               :: itask
  
  ! initialize tasks
  ntasks = 20
  allocate(tasks(ntasks))
  do itask=1, ntasks
     tasks(itask) = itask
  enddo

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_next_task(itask)
! put a new piece of work in the send buffer
  
  use work_type_mod
  integer, intent(in)   :: itask

  wt%ntotal_kernel = tasks(itask)
  wt%ndimensions = 10

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine extract_receive_buffer(ioutput, output)
! extract information received back from a slave

  use work_type_mod
  integer, intent(in)       :: ioutput
  integer, intent(inout)    :: output(:,:)

  output(ioutput,1) = wt%ntotal_kernel
  output(ioutput,2) = wt%ndimensions

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
