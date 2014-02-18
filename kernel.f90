module kernel
use global_parameters,                    only: sp, dp, verbose
use filtering,                            only: filter_type
implicit none
    type kernelspec_type
        character(len=16)                    :: name
        real(kind=sp), dimension(2)          :: time_window
        real(kind=sp), allocatable           :: veloseis(:)       ! Velocity seismogram 
                                                                  ! in the time window of 
                                                                  ! this kernel
        real(kind=dp), allocatable           :: t(:)
        real(kind=dp)                        :: dt
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
           procedure, pass                   :: init 
           procedure, pass                   :: calc_misfit_kernel
           procedure, pass                   :: isinitialized
           procedure, pass                   :: freeme
           procedure, pass                   :: apply_filter
    end type


contains 

!-------------------------------------------------------------------------------
subroutine init(this, name, time_window, filter, misfit_type, model_parameter, &
                veloseis, dt, timeshift_fwd)
   use fft,       only                      : rfft_type, taperandzeropad
   use filtering, only                      : timeshift
   class(kernelspec_type)                  :: this
   character(len=16), intent(in)           :: name
   real(kind=sp), intent(in)               :: time_window(2)
   type(filter_type), target, intent(in)   :: filter
   character(len=4), intent(in)            :: misfit_type
   character(len=4), intent(in)            :: model_parameter
   real(kind=sp), intent(in)               :: veloseis(:)
   real(kind=dp), intent(in)               :: dt
   real(kind=sp), intent(in)               :: timeshift_fwd
   
   type(rfft_type)                         :: fft_data
   complex(kind=dp), allocatable           :: veloseis_fd(:,:)
   real(kind=sp),    allocatable           :: veloseis_td(:,:), t_cut(:)
   real(kind=dp),    allocatable           :: veloseis_filtered(:,:)
   character(len=32)                       :: fmtstring

   integer                                 :: ntimes, ntimes_ft, nomega, isample

   if(this%initialized) then
      stop 'This kernel is already initialized'
   end if
   this%name              = name
   this%time_window       = time_window 
   this%filter            => filter
   this%misfit_type       = misfit_type
   this%model_parameter   = model_parameter
   this%initialized       = .true.
   this%dt                = dt

   if (verbose>0) then
      print *, '  ---------------------------------------------------------'
      print '(2(A))',         '   Initialize kernel ', this%name
      print *, '  ---------------------------------------------------------'
      print '(A,2(F8.1))',    '   Time window:  ', this%time_window
      print '(2(A))',         '   Misfit type:  ', this%misfit_type
      print '(2(A))',         '   Model param:  ', this%model_parameter
      print '(A,L)',          '   Initialized:  ', this%initialized
      print '(2(A))',         '   Filter type:  ', this%filter%filterclass
      print '(A,4(F8.3))',    '   Filter freq:  ', this%filter%frequencies
      print '(A,F8.3)',       '   Time shift :  ', timeshift_fwd
   end if
   !print '(A,I6)',         'len(Veloseis):', size(veloseis,1) 


   ! Save velocity seismogram in the timewindow of this kernel

   ! Allocate temporary variable
   ntimes = size(veloseis,1)
   allocate(veloseis_td(ntimes, 1))
   veloseis_td(:,1) = veloseis

   ! Initialize FFT for the velocity seismogram
   call fft_data%init(ntimes, 1, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   fmtstring = '(A, I8, A, I8)'
   print fmtstring, '  ntimes: ',  ntimes,     '  , nfreq: ', nomega

   allocate(veloseis_fd(nomega,1))
   allocate(veloseis_filtered(ntimes_ft,1))
   allocate(this%t(ntimes_ft))
   this%t = fft_data%get_t()


   ! FFT, timeshift and filter the velocity seismogram
   call fft_data%rfft(taperandzeropad(veloseis_td, ntimes_ft), veloseis_fd)
   call timeshift(veloseis_fd, fft_data%get_f(), timeshift_fwd)
   veloseis_fd = this%filter%apply_2d(veloseis_fd)
   call fft_data%irfft(veloseis_fd, veloseis_filtered)
   call fft_data%freeme()

   ! Cut timewindow of the kernel from the velocity seismogram
   call cut_timewindow(this%t,                 &
                       veloseis_filtered(:,1), &
                       this%time_window,       &
                       this%veloseis )
   call cut_timewindow(this%t,                 &
                       this%t,                 &
                       this%time_window,       &
                       t_cut )

   if (abs(sum(this%veloseis**2)).lt.1.e-20) then
       this%normalization = 0
   else
       this%normalization = sum(this%veloseis**2)
   end if
   
   if (verbose>0) then
      print *, '  Normalization coefficient: ', this%normalization
      print *, '  Length of seismogram: ', size(this%veloseis,1), ' samples'
      print *, '  ---------------------------------------------------------'
      print *, ''
   end if

   open(unit=100,file='seismogram_raw_'//trim(this%name), action='write')
   do isample = 1, size(veloseis,1)
      write(100,*) this%t(isample), veloseis(isample)
   end do
   close(100)

   open(unit=100,file='seismogram_'//trim(this%name), action='write')
   do isample = 1, size(veloseis_filtered,1)
      write(100,*) this%t(isample), veloseis_filtered(isample,1)
   end do
   close(100)

   open(unit=100,file='seismogram_cut_'//trim(this%name), action='write')
   do isample = 1, size(this%veloseis,1)
      write(100,*) t_cut(isample), this%veloseis(isample)
   end do
   close(100)

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
function calc_misfit_kernel(this, timeseries)
   !< This routine can take multiple time series and calculate the kernel at each
   class(kernelspec_type)                  :: this
   real(kind=dp), intent(in)               :: timeseries(:,:)
   real(kind=sp)                           :: calc_misfit_kernel(size(timeseries,2))
   real(kind=sp), allocatable              :: cut_timeseries(:)
   real(kind=sp)                           :: dt
   integer                                 :: lenseis, itrace, ntrace

   ntrace  = size(timeseries,2)

   select case(trim(this%misfit_type))
   case('CC')
      do itrace = 1, ntrace
         call cut_timewindow(this%t,                 &
                             timeseries(:, itrace),  &
                             this%time_window,       &
                             cut_timeseries)
         
         calc_misfit_kernel(itrace) = sum( cut_timeseries * this%veloseis, 1 ) &
                                      * this%dt / this%normalization
      end do
   end select

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine cut_timewindow(t, x, timewindow, cut_tw)
   real(kind=dp), intent(in)  :: t(:), x(:)
   real(kind=sp), intent(in)  :: timewindow(2)
   real(kind=sp), allocatable :: cut_tw(:)
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
         cut_timewindow_temp(iintimewindow) = x(i)
      end if
      if (t(i).gt.timewindow(2)) exit
   end do

   if (iintimewindow.eq.0) then
       print *, 'Time window: ', timewindow, ', t(1):', t(1), ', t(ntimes):', t(ntimes)
       stop 'Time window length was zero'
   end if

   !if (allocated(cut_timewindow)) then
   if (allocated(cut_tw)) then
       !print *, 'Warning: Time window was already allocated'
       deallocate(cut_tw)
   end if
   allocate(cut_tw(iintimewindow))

   cut_tw(:) = cut_timewindow_temp(1:iintimewindow)
   

end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function apply_filter(this, freq_series)
   class(kernelspec_type)           :: this
   complex(kind=dp), intent(in)     :: freq_series(:,:)
   complex(kind=dp)                 :: apply_filter(size(freq_series,1), &
                                                    size(freq_series,2))


   apply_filter = this%filter%apply_2d(freq_series)
   apply_filter = this%filter%apply_2d(apply_filter)

end function apply_filter
!-------------------------------------------------------------------------------




end module
