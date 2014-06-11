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
   real(kind=dp)              :: timeseries(ntimes), dataseries(ntimes)
   real(kind=dp)              :: dt = 0.1d0
   real(kind=dp), allocatable :: cut_dataseries(:), cut_dataseries_ref(:)
   integer                    :: i

   ! Time window with normal length
   do i=1, ntimes
      timeseries(i) = i * dt
      dataseries(i) = i 
      ! print *, timeseries(i), dataseries(i)
   end do

   
   ! Test 1: Time window bounds are between samples
   allocate(cut_dataseries_ref(6))
   cut_dataseries_ref = [15, 16, 17, 18, 19, 20]

   call cut_timewindow(timeseries, dataseries, [1.45d0, 2.01d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 1)')

   deallocate(cut_dataseries_ref)
  
   ! Test 2: Upper time window bound is exactly on sample
   !         Correct behaviour would be that this sample is included in the resulting time window.
   !         This test may fail due to floating point comparison, however, 2 should have an exact 
   !         floating point representation 
   allocate(cut_dataseries_ref(6))
   cut_dataseries_ref = [15, 16, 17, 18, 19, 20]

   call cut_timewindow(timeseries, dataseries, [1.45d0, 2.0d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 2)')

   deallocate(cut_dataseries_ref)

   ! Test 3: Upper time window bound is narrowly (but within real(8) resolution) below sample
   allocate(cut_dataseries_ref(5))
   cut_dataseries_ref = [15, 16, 17, 18, 19]

   call cut_timewindow(timeseries, dataseries, [1.45d0, 1.9999999d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 3)')


   deallocate(cut_dataseries_ref)

   ! Test 4: Time window with length 1
   allocate(cut_dataseries_ref(1))
   cut_dataseries_ref = [80]

   call cut_timewindow(timeseries, dataseries, [7.95d0, 8.05d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 4)' )



   ! Test 5: Time window with negative starting time
   do i=1, ntimes
      timeseries(i) = i * dt - 5
      dataseries(i) = i - 50
   end do

   deallocate(cut_dataseries_ref)
   allocate(cut_dataseries_ref(6))
   
   cut_dataseries_ref = [-3, -2, -1, 0, 1, 2]

   call cut_timewindow(timeseries, dataseries, [-0.39d0, 0.29d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 5)')


end subroutine test_kernel_cut_timewindow

!-------------------------------------------------------------------------------
subroutine test_kernel_init

   type(kernelspec_type)       :: kernel
   type(filter_type), target   :: gabor

   character(len=32)           :: filtername, filterclass
   real(kind=dp)               :: veloseis(256)

   call assert_false(kernel%isinitialized(), 'Kernel not initialized')
   ! Just some arbitrary filter to create
   filtername  = 'Gabor'
   filterclass = 'Gabor'
   call gabor%create(filtername, 1.0d0, 257, filterclass, [5.0d0, 0.5d0, 0.d0, 0.d0])

   veloseis = 0.0

   call kernel%init(name            = 'Testkernel      ',&
                    time_window     = [1.0d0, 2.0d0],    &
                    filter          = gabor,             &
                    misfit_type     = 'CC  ',            &  
                    model_parameter = 'vp  ',            &
                    seis            = veloseis,          &
                    dt              = 0.1d0,             &
                    timeshift_fwd   = 1.0d0   )

   call assert_true(kernel%isinitialized(), 'Kernel initialized')
   call kernel%freeme()
   
   call assert_false(kernel%isinitialized(), 'Kernel not initialized')

end subroutine test_kernel_init

!-------------------------------------------------------------------------------

end module test_kernel
