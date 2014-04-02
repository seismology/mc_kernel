module kernel
use global_parameters,                    only: sp, dp, verbose, lu_out
use filtering,                            only: filter_type
implicit none

    type kernelspec_type
        character(len=32), public            :: name
        real(kind=dp), dimension(2)          :: time_window
        real(kind=dp), allocatable           :: seis(:)           ! Seismogram (Velocity or displacement)
                                                                  ! in the time window of 
                                                                  ! this kernel
        real(kind=dp), allocatable           :: t(:)
        real(kind=dp)                        :: dt
        real(kind=dp)                        :: normalization
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
           procedure, pass                   :: apply_filter_2d
           procedure, pass                   :: apply_filter_3d
           generic                           :: apply_filter => apply_filter_2d, apply_filter_3d
    end type

    !interface apply_filter
    !   module procedure                      :: apply_filter_3d
    !   module procedure                      :: apply_filter_2d
    !end interface apply_filter

contains 

!-------------------------------------------------------------------------------
subroutine init(this, name, time_window, filter, misfit_type, model_parameter, &
                seis, dt, timeshift_fwd)
   use fft,       only                      : rfft_type, taperandzeropad
   use filtering, only                      : timeshift
   class(kernelspec_type)                  :: this
   character(len=16), intent(in)           :: name
   real(kind=dp), intent(in)               :: time_window(2)
   type(filter_type), target, intent(in)   :: filter
   character(len=4), intent(in)            :: misfit_type
   character(len=4), intent(in)            :: model_parameter
   real(kind=dp), intent(in)               :: seis(:)
   real(kind=dp), intent(in)               :: dt
   real(kind=dp), intent(in)               :: timeshift_fwd
   
   type(rfft_type)                         :: fft_data
   complex(kind=dp), allocatable           :: seis_fd(:,:)
   real(kind=dp),    allocatable           :: seis_td(:,:), t_cut(:)
   real(kind=dp),    allocatable           :: seis_filtered(:,:)
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
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,'(2(A))')         '   Initialize kernel ', this%name
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,'(A,2(F8.1))')    '   Time window:  ', this%time_window
      write(lu_out,'(2(A))')         '   Misfit type:  ', this%misfit_type
      write(lu_out,'(2(A))')         '   Model param:  ', this%model_parameter
      write(lu_out,'(A,L)')          '   Initialized:  ', this%initialized
      write(lu_out,'(2(A))')         '   Filter type:  ', this%filter%filterclass
      write(lu_out,'(A,4(F8.3))')    '   Filter freq:  ', this%filter%frequencies
      write(lu_out,'(A,F8.3)')       '   Time shift :  ', timeshift_fwd
   end if


   ! Save seismogram in the timewindow of this kernel

   ! Allocate temporary variable
   ntimes = size(seis, 1)
   allocate(seis_td(ntimes, 1))
   seis_td(:,1) = seis

   ! Initialize FFT for the seismogram
   call fft_data%init(ntimes, 1, 1, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   fmtstring = '(A, I8, A, I8)'
   write(lu_out,fmtstring) '  ntimes: ',  ntimes,     '  , nfreq: ', nomega

   allocate(seis_fd(nomega, 1))
   allocate(seis_filtered(ntimes_ft, 1))
   allocate(this%t(ntimes_ft))
   this%t = fft_data%get_t()


   ! FFT, timeshift and filter the seismogram
   call fft_data%rfft(taperandzeropad(seis_td, ntimes_ft), seis_fd)
   call timeshift(seis_fd, fft_data%get_f(), timeshift_fwd)
   seis_fd = this%filter%apply_2d(seis_fd)
   call fft_data%irfft(seis_fd, seis_filtered)
   call fft_data%freeme()

   ! Cut timewindow of the kernel from the seismogram
   call cut_timewindow(this%t,             &
                       seis_filtered(:,1), &
                       this%time_window,   &
                       this%seis )
   call cut_timewindow(this%t,             &
                       this%t,             &
                       this%time_window,   &
                       t_cut )

   if (abs(sum(this%seis**2)).lt.1.e-30) then
       this%normalization = 0
   else
       this%normalization = 1./sum(this%seis**2)
   end if
   
   if (verbose>0) then
      write(lu_out,*) '  Normalization coefficient: ', this%normalization
      write(lu_out,*) '  Length of seismogram: ', size(this%seis,1), ' samples'
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,*) ''
   end if

   open(unit=100,file='seismogram_raw_'//trim(this%name), action='write')
   do isample = 1, size(seis,1)
      write(100,*) this%t(isample), seis(isample)
   end do
   close(100)

   open(unit=100,file='seismogram_'//trim(this%name), action='write')
   do isample = 1, size(seis_filtered,1)
      write(100,*) this%t(isample), seis_filtered(isample,1)
   end do
   close(100)

   open(unit=100,file='seismogram_cut_'//trim(this%name), action='write')
   do isample = 1, size(this%seis,1)
      write(100,*) t_cut(isample), this%seis(isample)
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
   real(kind=dp)                           :: calc_misfit_kernel(size(timeseries,2))
   real(kind=dp), allocatable              :: cut_timeseries(:)
   real(kind=dp)                           :: dt
   integer                                 :: lenseis, itrace, ntrace, lu_errorlog
   character(len=64)                       :: fmtstring

   ntrace = size(timeseries,2)

   select case(trim(this%misfit_type))
   case('CC')
      do itrace = 1, ntrace
         call cut_timewindow(this%t,                 &
                             timeseries(:, itrace),  &
                             this%time_window,       &
                             cut_timeseries)
         
         calc_misfit_kernel(itrace) = integrate( cut_timeseries * this%seis, this%dt ) !&
                                      !* this%normalization
      end do

   case('AM')
      do itrace = 1, ntrace
         call cut_timewindow(this%t,                 &
                             timeseries(:, itrace),  &
                             this%time_window,       &
                             cut_timeseries)
         
         calc_misfit_kernel(itrace) = integrate( cut_timeseries * this%seis, this%dt )! &
                                      !* this%normalization
      end do
   end select

   if (any(calc_misfit_kernel.ne.calc_misfit_kernel)) then ! Kernel is NaN
       print *, 'ERROR Kernel has value NaN!'
       print *, 'Values of Kernel:'
       write(fmtstring, "('(',I6,'(ES11.3))')") ntrace
       print fmtstring, calc_misfit_kernel
       print *, 'Value of normalization:'
       print *, this%normalization
       print *, 'Convolved time series is written into "errorlog_timeseries",'
       print *, 'seismogram is written into "errorlog_seismogram"'
       open(newunit=lu_errorlog, file='errorlog_timeseries', status='replace')
       write(lu_errorlog, *) cut_timeseries
       close(lu_errorlog)
       open(newunit=lu_errorlog, file='errorlog_seismogram', status='replace')
       write(lu_errorlog, *) this%seis
       close(lu_errorlog)
       stop
   end if

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine cut_timewindow(t, x, timewindow, cut_tw)
   real(kind=dp), intent(in)  :: t(:), x(:)
   real(kind=dp), intent(in)  :: timewindow(2)
   real(kind=dp), allocatable :: cut_tw(:)
   real(kind=dp), allocatable :: cut_timewindow_temp(:)
   integer                    :: ntimes, i, iintimewindow
   integer                    :: ntraces
   real(kind=dp)              :: dt
   
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
function apply_filter_3d(this, freq_series)
   class(kernelspec_type)           :: this
   complex(kind=dp), intent(in)     :: freq_series(:,:,:)
   complex(kind=dp)                 :: apply_filter_3d(size(freq_series,1), &
                                                       size(freq_series,2), &
                                                       size(freq_series,3))


   apply_filter_3d = this%filter%apply_3d(freq_series)
   apply_filter_3d = this%filter%apply_3d(apply_filter_3d)

end function apply_filter_3d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function apply_filter_2d(this, freq_series)
   class(kernelspec_type)           :: this
   complex(kind=dp), intent(in)     :: freq_series(:,:)
   complex(kind=dp)                 :: apply_filter_2d(size(freq_series,1), &
                                                       size(freq_series,2))


   apply_filter_2d = this%filter%apply_2d(freq_series)
   apply_filter_2d = this%filter%apply_2d(apply_filter_2d)

end function apply_filter_2d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
pure function integrate(timeseries, dt)
    real(kind=dp), intent(in)   :: timeseries(:), dt
    real(kind=dp)               :: integrate
    integer                     :: npoints

    npoints = size(timeseries, 1)

    ! Trapezoidal rule: I = dt/2 * (f(x1) + 2f(x2) + ... + 2f(xN-1) + f(xN))
    integrate = timeseries(1) + timeseries(npoints) + 2*sum(timeseries(2:npoints-1), 1)
    integrate = integrate * dt * 0.5

end function integrate
!-------------------------------------------------------------------------------




end module
