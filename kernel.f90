module kernel
use global_parameters, only: sp, dp
use filtering,         only: filter_type
implicit none
    type kernelspec_type
        real, dimension(2)                   :: time_window
        real, allocatable, dimension(:)      :: veloseis          ! Velocity seismogram 
                                                                  ! in the time window of 
                                                                  ! this kernel
        real(kind=sp)                        :: normalization
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

!-------------------------------------------------------------------------------
subroutine init(this, time_window, filter, misfit_type, model_parameter, &
                veloseis, dt)
   use fft, only: rfft_type, taperandzeropad
   class(kernelspec_type)                  :: this
   real(kind=sp), intent(in)               :: time_window(2)
   type(filter_type), target, intent(in)   :: filter
   real(kind=sp), intent(in)               :: veloseis(:)
   real(kind=dp), intent(in)               :: dt
   
   type(rfft_type)                         :: fft_data
   complex(kind=dp), allocatable           :: veloseis_fd(:,:)
   real(kind=sp), allocatable              :: veloseis_td(:,:)
   real(kind=dp), allocatable              :: veloseis_filtered(:,:), t(:)
   character(len=4), intent(in)            :: misfit_type
   character(len=4), intent(in)            :: model_parameter

   integer                                 :: ntimes, ntimes_ft, nomega

   if(this%initialized) then
      stop 'This kernel is already initialized'
   end if
   this%time_window       = time_window + 15.0
   this%filter            => filter
   this%misfit_type       = misfit_type
   this%model_parameter   = model_parameter
   this%initialized       = .true.


   ! Save velocity seismogram in the timewindow of this kernel
   ntimes = size(veloseis,1)
   !print *, 'Length of velocity seismogram', ntimes
   allocate(veloseis_td(ntimes, 1))
   veloseis_td(:,1) = veloseis

   call fft_data%init(ntimes, 1, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   allocate(veloseis_fd(nomega,1))
   allocate(veloseis_filtered(ntimes_ft,1))

   call fft_data%rfft(taperandzeropad(veloseis_td, ntimes_ft), veloseis_fd)
   veloseis_fd = this%filter%apply_2d(veloseis_fd)
   call fft_data%irfft(veloseis_fd, veloseis_filtered)

   allocate(t(ntimes_ft))
   t = fft_data%get_t()

   call fft_data%freeme()

   this%veloseis = cut_timewindow(t,                      &
                                  veloseis_filtered(:,1), &
                                  this%time_window)

   !print *, 'Length of time window', size(this%veloseis,1)
                                 
   this%normalization = sum(this%veloseis**2)
   print *, 'Normalization coefficient: ', this%normalization

end subroutine init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
logical function isinitialized(this)
   class(kernelspec_type)                   :: this
   isinitialized = this%initialized
end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine freeme(this)
   class(kernelspec_type)                   :: this
   this%filter => NULL()
   this%initialized = .false.
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function calc_misfit_kernel(this, t, timeseries)
   !< This routine can take multiple time series and calculate the kernel at each
   class(kernelspec_type)                  :: this
   real(kind=dp), intent(in)               :: timeseries(:,:), t(:)
   real(kind=sp)                           :: calc_misfit_kernel(size(timeseries,2))
   real(kind=sp)                           :: dt
   integer                                 :: lenseis, itrace, ntrace

   ntrace  = size(timeseries,2)

   select case(trim(this%misfit_type))
   case('CC')
      dt = t(2) - t(1)
      do itrace = 1, ntrace
         
         calc_misfit_kernel(itrace) = sum(  cut_timewindow(t,                               &
                                                           timeseries(:, itrace),           &
                                                           this%time_window)                & 
                                          * this%veloseis, 1  )                             &
                                      * dt / this%normalization
         !print *, 'Kernel value: ', calc_misfit_kernel(itrace), ', dt:', dt
      end do
   end select

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
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

   if (iintimewindow.eq.0) then
       print *, 'Time window: ', timewindow, ', t(1):', t(1), ', t(ntimes):', t(ntimes)
       stop 'Time window length was zero'
   end if

   if (allocated(cut_timewindow)) then
       print *, 'Warning: Time window was already allocated'
       deallocate(cut_timewindow)
   end if
   allocate(cut_timewindow(iintimewindow))

   cut_timewindow(:) = cut_timewindow_temp(1:iintimewindow)
   

end function
!-------------------------------------------------------------------------------

end module
