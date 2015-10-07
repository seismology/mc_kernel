!=========================================================================================
module slave_work
  use global_parameters
  use work_type_mod

  implicit none
  private
  
  public :: init_work
  public :: work

contains

!-----------------------------------------------------------------------------------------
subroutine init_work()
! put anything that should happen before receiving any task from the master here  

end subroutine init_work
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine work()
! actual work to be done based on the messages received from the master

  wt%ndimensions = wt%ntotal_kernel
  wt%vertices = wt%ntotal_kernel
  wt%ntotal_kernel = wt%ntotal_kernel ** 2

end subroutine work
!-----------------------------------------------------------------------------------------


end module
!=========================================================================================
