!=========================================================================================
module test_kernel

  use global_parameters
  use kernel
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------

subroutine test_kernel_cut_timewindow
   integer, parameter         :: ntimes = 100
   real(kind=sp)              :: timeseries(ntimes), dataseries(ntimes,1)
   real(kind=sp)              :: dt = 0.1
   real(kind=sp), allocatable :: cut_dataseries(:,:), cut_dataseries_ref(:)
   integer                    :: i

   do i=1, ntimes
      timeseries(i) = i * dt
      dataseries(i,1) = i 
   end do

   cut_dataseries = cut_timewindow(timeseries, dataseries, [1.45, 2.01])
   
   allocate(cut_dataseries_ref(6))
   cut_dataseries_ref = [15, 16, 17, 18, 19, 20]

   call assert_comparable_real1d(cut_dataseries(:,1), cut_dataseries_ref, 1e-9, &
                                 'cut_dataseries == cut_dataseries_ref')

   cut_dataseries = cut_timewindow(timeseries, dataseries, [7.95, 8.05])
  
   deallocate(cut_dataseries_ref)
   allocate(cut_dataseries_ref(1))
   cut_dataseries_ref = [80]

   call assert_comparable_real1d(cut_dataseries(:,1), cut_dataseries_ref, 1e-9, &
                                 'cut_dataseries == cut_dataseries_ref')

end subroutine test_kernel_cut_timewindow


end module test_kernel
