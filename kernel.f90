module kernel
use global_parameters, only: sp, dp
use filtering,         only: filter_type
implicit none
    type kernelspec_type
        real, dimension(2)                   :: time_window
        real, dimension(4)                   :: corner_frequencies
        integer                              :: filter_type
        character(len=4)                     :: misfit_type
        character(len=4)                     :: model_parameter
        !integer                              :: receiver_index
        !integer                              :: src_index
        !type(rec_param_type), pointer        :: receiver
        type(filter_type), pointer           :: filter
        logical                              :: initialized = .false.

        contains 
           procedure, pass                   :: calc_misfit_kernel
           procedure, pass                   :: init 
           procedure, pass                   :: isinitialized
           procedure, pass                   :: freeme
    end type


contains 

subroutine init(this, time_window, filter, misfit_type, model_parameter)
   class(kernelspec_type)                  :: this
   real(kind=sp), dimension(2), intent(in) :: time_window
   type(filter_type), target, intent(in)   :: filter
   character(len=4), intent(in)            :: misfit_type
   character(len=4), intent(in)            :: model_parameter

   if(this%initialized) then
      stop 'This kernel is already initialized'
   end if
   this%time_window       = time_window
   this%filter            => filter
   this%misfit_type       = misfit_type
   this%model_parameter   = model_parameter
   this%initialized       = .true.

end subroutine

logical function isinitialized(this)
   class(kernelspec_type)                   :: this
   isinitialized = this%initialized
end function

subroutine freeme(this)
   class(kernelspec_type)                   :: this
   this%filter => NULL()
   this%initialized = .false.
end subroutine

function calc_misfit_kernel(this, t, timeseries, veloseis)
   !< This routine can take multiple time series and calculate the kernel at each
   class(kernelspec_type)                  :: this
   real(kind=dp), intent(in)               :: timeseries(:,:), t(:), veloseis(:)
   real(kind=sp)                           :: calc_misfit_kernel(size(timeseries,2))
   real(kind=sp)                           :: dt
   integer                                 :: lenseis, itrace, ntrace

   lenseis = size(veloseis,1)
   ntrace  = size(timeseries,2)

   select case(trim(this%misfit_type))
   case('CC')
      dt = t(2) - t(1)
      do itrace = 1, ntrace
         
         calc_misfit_kernel(itrace) = sum(  cut_timewindow(t,                               &
                                                           timeseries(:, itrace),           &
                                                           this%time_window)                & 
                                          * cut_timewindow(t(1:lenseis),                    &
                                                           veloseis,  &
                                                           this%time_window), 1)            &
                                      * dt
         !print *, 'Kernel value: ', calc_misfit_kernel(itrace), ', dt:', dt
      end do
   end select

end function

function cut_timewindow(t, f, timewindow)
   real(kind=dp), intent(in)  :: t(:), f(:)
   real(kind=sp), intent(in)  :: timewindow(2)
   real(kind=sp), allocatable :: cut_timewindow(:)
   real(kind=sp), allocatable :: cut_timewindow_temp(:)
   integer                    :: ntimes, i, iintimewindow
   integer                    :: ntraces
   real(kind=sp)              :: dt
   
   ntimes = size(t)
   dt = t(2) - t(1)
   if(timewindow(1).lt.t(1)) then
      write(*,*) 'Time window starts before beginning of time series'
      write(*,*) 'time window:', timewindow, '; t(1):', t(1)
      stop
   end if
   if(timewindow(2).gt.t(ntimes)) then
      write(*,*) 'Time window ends after beginning of time series'
      write(*,*) 'time window:', timewindow, '; t(ntimes):', t(ntimes)
      stop
   end if

   if ((timewindow(2) - timewindow(1)).le.0) then
      write(*,*) 'length of time window is negative'
      write(*,*) 'Beginning: ', timewindow(1), '; end: ', timewindow(2)
      stop
   end if

   allocate(cut_timewindow_temp(ntimes))
   iintimewindow = 0
   do i = 1, ntimes
      if (t(i).le.timewindow(2).and.t(i).ge.timewindow(1)) then
         iintimewindow = iintimewindow + 1
         cut_timewindow_temp(iintimewindow) = f(i)
      end if
      if (t(i).gt.timewindow(2)) exit
   end do

   if (allocated(cut_timewindow)) deallocate(cut_timewindow)
   allocate(cut_timewindow(iintimewindow))

   cut_timewindow(:) = cut_timewindow_temp(1:iintimewindow)
   

end function

end module
