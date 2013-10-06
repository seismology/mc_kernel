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

   ! Time window with normal length
   do i=1, ntimes
      timeseries(i) = i * dt
      dataseries(i,1) = i 
   end do

   
   allocate(cut_dataseries_ref(6))
   cut_dataseries_ref = [15, 16, 17, 18, 19, 20]

   cut_dataseries = cut_timewindow(timeseries, dataseries, [1.45, 2.01])

   call assert_comparable_real1d(cut_dataseries(:,1), cut_dataseries_ref, 1e-9, &
                                 'cut_dataseries == cut_dataseries_ref')

   deallocate(cut_dataseries_ref)
  
   allocate(cut_dataseries_ref(6))
   cut_dataseries_ref = [15, 16, 17, 18, 19, 20]

   cut_dataseries = cut_timewindow(timeseries, dataseries, [1.45, 2.00])

   call assert_comparable_real1d(cut_dataseries(:,1), cut_dataseries_ref, 1e-9, &
                                 'cut_dataseries == cut_dataseries_ref')

   deallocate(cut_dataseries_ref)

   allocate(cut_dataseries_ref(6))
   cut_dataseries_ref = [15, 16, 17, 18, 19]

   cut_dataseries = cut_timewindow(timeseries, dataseries, [1.45, 1.99])

   call assert_comparable_real1d(cut_dataseries(:,1), cut_dataseries_ref, 1e-9, &
                                 'cut_dataseries == cut_dataseries_ref')

   deallocate(cut_dataseries_ref)

   ! Time window with length 1
   allocate(cut_dataseries_ref(1))
   cut_dataseries_ref = [80]

   cut_dataseries = cut_timewindow(timeseries, dataseries, [7.95, 8.05])

   call assert_comparable_real1d(cut_dataseries(:,1), cut_dataseries_ref, 1e-9, &
                                 'cut_dataseries == cut_dataseries_ref')


   ! Time window with negative starting time

   do i=1, ntimes
      timeseries(i) = i * dt - 5
      dataseries(i,1) = i - 50
   end do

   deallocate(cut_dataseries_ref)
   allocate(cut_dataseries_ref(6))
   
   cut_dataseries_ref = [-3, -2, -1, 0, 1, 2]

   cut_dataseries = cut_timewindow(timeseries, dataseries, [-0.39, 0.29])

   call assert_comparable_real1d(cut_dataseries(:,1), cut_dataseries_ref, 1e-9, &
                                 'cut_dataseries == cut_dataseries_ref')

end subroutine test_kernel_cut_timewindow

!-------------------------------------------------------------------------------
subroutine test_kernel_init

   type(kernelspec_type)       :: kernel
   type(filter_type), target   :: gabor

   character(len=32)           :: filtername

   ! Just some arbitrary filter to create
   filtername = 'Gabor'
   call gabor%create(1.0, 256, filtername, [5.0, 0.5, 0., 0.])


   call kernel%init(time_window     = [1.0, 2.0], &
                    filter          = gabor,      &
                    misfit_type     = 'CC  ',     &  
                    model_parameter = 'vp  ')

   call assert_true(kernel%isinitialized(), 'Kernel initialized')
   call kernel%freeme()
   
   call assert_false(kernel%isinitialized(), 'Kernel not initialized')

end subroutine test_kernel_init

!-------------------------------------------------------------------------------

end module test_kernel
