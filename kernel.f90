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
   integer                                 :: lenseis

   lenseis = size(veloseis,1)

   select case(trim(this%misfit_type))
   case('CC')
      dt = t(2) - t(1)
      calc_misfit_kernel = sum(  cut_timewindow(t, timeseries,          &
                                                this%time_window)       & 
                               * cut_timewindow(t(1:lenseis),          &
                                                reshape(veloseis, [lenseis,1]),              &
                                                this%time_window), 1)  &
                           * dt
   end select

end function

function cut_timewindow(t, f, timewindow)
   real(kind=dp), intent(in)  :: t(:), f(:,:)
   real(kind=sp), intent(in)  :: timewindow(2)
   real(kind=sp), allocatable :: cut_timewindow(:,:)
   real(kind=sp), allocatable :: cut_timewindow_temp(:,:)
   integer                    :: ntimes, i, iintimewindow
   integer                    :: ntraces
   real(kind=sp)              :: dt
   
   ntimes = size(t)
   ntraces = size(f,2)
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

   allocate(cut_timewindow_temp(ntimes, ntraces))
   iintimewindow = 0
   do i = 1, ntimes
      if (t(i).le.timewindow(2).and.t(i).ge.timewindow(1)) then
         iintimewindow = iintimewindow + 1
         cut_timewindow_temp(iintimewindow,:) = f(i,:)
      end if
      if (t(i).gt.timewindow(2)) exit
   end do

   if (allocated(cut_timewindow)) deallocate(cut_timewindow)
   allocate(cut_timewindow(iintimewindow, ntraces))

   cut_timewindow(:,:) = cut_timewindow_temp(1:iintimewindow, :)
   

end function

end module
