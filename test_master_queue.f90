!=========================================================================================
module test_master_queue

  use global_parameters
  use master_queue
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
! These routines are difficult to test one by one, so one big test for all
subroutine test_master_all
  integer       ::  ntasks

  testing = .false.
  call init_queue(ntasks, 'unit_tests/inparam_test')

  call finalize

end subroutine test_master_all
!-----------------------------------------------------------------------------------------

end module test_master_queue

