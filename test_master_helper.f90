!=========================================================================================
module test_master_helper

  use global_parameters
  use master_helper
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_create_tasks
  logical               :: completed(10)
  integer               :: nelems_per_task, ntasks, ntasks_ref
  integer, allocatable  :: elems_in_task(:,:), elems_in_task_ref(:,:)

  ! Test 1: No element completed yet, 5 elems per task. 
  !         Result should be 2 tasks, with elems 1:5 and 6:10

  completed = .false.
  nelems_per_task = 5

  ntasks_ref = 2
  allocate(elems_in_task_ref(ntasks_ref, nelems_per_task))
  elems_in_task_ref(1,:) = [1, 2, 3, 4,  5]
  elems_in_task_ref(2,:) = [6, 7, 8, 9, 10]

  call create_tasks(completed, nelems_per_task, ntasks, elems_in_task)
  
  call assert_equal(ntasks, ntasks_ref, 'Create 2 tasks')
  call assert_equal(elems_in_task(1,:), elems_in_task_ref(1,:), 'First task: 1:5')
  call assert_equal(elems_in_task(2,:), elems_in_task_ref(2,:), 'Second task: 6:10')
  deallocate(elems_in_task)
  deallocate(elems_in_task_ref)

  ! Test 2: First three elements completed, 5 elems per task. 
  !         Result should be 2 tasks, with elems 4:8 and [9,10,-1,-1,-1]

  completed = .false.
  completed(1:3) = .true.
  nelems_per_task = 5

  ntasks_ref = 2
  allocate(elems_in_task_ref(ntasks_ref, nelems_per_task))
  elems_in_task_ref(1,:) = [4, 5, 6, 7,  8]
  elems_in_task_ref(2,:) = [9, 10, -1, -1, -1]

  call create_tasks(completed, nelems_per_task, ntasks, elems_in_task)
  
  call assert_equal(ntasks, ntasks_ref, 'Create 2 tasks')
  call assert_equal(elems_in_task(1,:), elems_in_task_ref(1,:), 'First task: 4:8')
  call assert_equal(elems_in_task(2,:), elems_in_task_ref(2,:), 'Second task: 9:10, then -1')
  deallocate(elems_in_task)
  deallocate(elems_in_task_ref)

  ! Test 3: Elements 2, 4, 9 completed, 5 elems per task. 
  !         Result should be 2 tasks, with elems [1,3,5,6,7] and [8,10,-1,-1,-1]

  completed = .false.
  completed([2,4,9]) = .true.
  nelems_per_task = 5

  ntasks_ref = 2
  allocate(elems_in_task_ref(ntasks_ref, nelems_per_task))
  elems_in_task_ref(1,:) = [1, 3, 5, 6, 7]
  elems_in_task_ref(2,:) = [8, 10, -1, -1, -1]

  call create_tasks(completed, nelems_per_task, ntasks, elems_in_task)
  
  call assert_equal(ntasks, ntasks_ref, 'Create 2 tasks')
  call assert_equal(elems_in_task(1,:), elems_in_task_ref(1,:), 'First task: 4:8')
  call assert_equal(elems_in_task(2,:), elems_in_task_ref(2,:), 'Second task: 9:10, then -1')
  deallocate(elems_in_task)
  deallocate(elems_in_task_ref)


  ! Test 4: Elements 2, 4, 9 completed, 2 elems per task. 
  !         Result should be 4 tasks, with elems [1,3], [5,6], [7,8] and [10,-1]

  completed = .false.
  completed([2,4,9]) = .true.
  nelems_per_task = 2

  ntasks_ref = 4
  allocate(elems_in_task_ref(ntasks_ref, nelems_per_task))
  elems_in_task_ref(1,:) = [1, 3]
  elems_in_task_ref(2,:) = [5, 6]
  elems_in_task_ref(3,:) = [7, 8]
  elems_in_task_ref(4,:) = [10, -1]

  call create_tasks(completed, nelems_per_task, ntasks, elems_in_task)
  
  call assert_equal(ntasks, ntasks_ref, 'Create 4 tasks')
  call assert_equal(elems_in_task(1,:), elems_in_task_ref(1,:), 'First task:  [1,3]')
  call assert_equal(elems_in_task(2,:), elems_in_task_ref(2,:), 'Second task: [5,6]')
  call assert_equal(elems_in_task(3,:), elems_in_task_ref(3,:), 'Third task:  [7,8]')
  call assert_equal(elems_in_task(4,:), elems_in_task_ref(4,:), 'Fourth task: [10,-1]')
  deallocate(elems_in_task)
  deallocate(elems_in_task_ref)


end subroutine test_create_tasks
!-----------------------------------------------------------------------------------------

end module test_master_helper
