!=========================================================================================
module resampling
  use global_parameters
  use commpi,               only : pabort
  use fft,                  only : rfft_type, taperandzeropad
  use filtering,            only : filter_type

  implicit none

  private
  public :: resampling_type

  type :: resampling_type
     private
     type(rfft_type)            :: fft_in, fft_out
     logical                    :: initialized = .false.
     integer                    :: ntimes_in, ntimes_out, ntraces
     real(kind=dp), allocatable :: data_out_buf(:,:)
     logical                    :: apply_filter = .false.
     type(filter_type)          :: filter
     contains
     procedure, pass            :: init
     procedure, pass            :: resample
     procedure, pass            :: resample_timeshift
     procedure, pass            :: freeme
     procedure, pass            :: get_ntraces
     procedure, pass            :: get_ntimes_in
     procedure, pass            :: get_ntimes_out
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine init(this, ntimes_in, ntimes_out, ntraces, filter)!, measure)

  class(resampling_type)        :: this
  integer, intent(in)           :: ntimes_in, ntimes_out, ntraces
  type(filter_type), optional   :: filter
  !logical, intent(in)       :: measure

  ! hardcoding ndim=1 for now, may be dummy variable in the future
  integer                   :: ndim 
  ndim = 1

  call this%fft_in%init(ntimes_in, ndim, ntraces, nfft=ntimes_in)!, measure=measure)
  call this%fft_out%init(ntimes_out, ndim, ntraces, nfft=ntimes_out)!, measure=measure)

  this%ntimes_in = ntimes_in
  this%ntimes_out = ntimes_out
  this%ntraces = ntraces

  allocate(this%data_out_buf(this%fft_out%get_ntimes(), ntraces))

  if(present(filter)) then
     this%filter = filter
     this%apply_filter = .true.
  else
     this%apply_filter = .false.
  endif

  this%initialized = .true.
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> resampling routine for rational number resampling, inspriation taken from
!! scipy.signal.resample
subroutine resample(this, data_in, data_out)

  class(resampling_type)            :: this
  real(kind=dp), intent(in)         :: data_in(:,:)
  real(kind=dp), intent(out)        :: data_out(:,:)

  complex(kind=dp), allocatable     :: dataf_in(:,:), dataf_out(:,:)
  integer                           :: nomega_min, i

  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing resampling type that is not initialized'
     call pabort 
  end if

  if (any(shape(data_in) /= [this%ntimes_in, this%ntraces])) then
     write(*,*) 'ERROR: shape mismatch in first argument - ' &
                    // 'shape does not match the plan for resampling'
     write(*,*) 'is: ', shape(data_in), '; should be: [', this%ntimes_in, &
                    ', ', this%ntraces, ']'
     call pabort
  end if

  if (any(shape(data_out) &
        /= [this%ntimes_out, this%ntraces])) then
     write(*,*) 'ERROR: shape mismatch in second argument - ' &
                    // 'shape does not match the plan for resampling'
     write(*,*) 'is: ', shape(data_out), '; should be: [', this%ntimes_out, &
                    ', ', this%ntraces, ']'
     call pabort
  end if

  allocate(dataf_in(1:this%fft_in%get_nomega(), 1:this%fft_in%get_ntraces()))
  allocate(dataf_out(1:this%fft_out%get_nomega(), 1:this%fft_out%get_ntraces()))

  call this%fft_in%rfft(taperandzeropad(data_in, this%fft_in%get_ntimes(), ntaper=0), dataf_in)
  
  nomega_min = min(this%fft_in%get_nomega(), this%fft_out%get_nomega())

  if (this%apply_filter) then
     do i = 1, this%fft_in%get_ntraces()
        dataf_in(:,i) = dataf_in(:,i) * this%filter%transferfunction
     enddo
  endif

  dataf_out(1:nomega_min,:) = dataf_in(1:nomega_min,:)

  if (nomega_min < this%fft_out%get_nomega()) &
     dataf_out(nomega_min+1:,:) = 0

  call this%fft_out%irfft(dataf_out, this%data_out_buf)

  data_out(:,:) = this%data_out_buf(1:this%ntimes_out,:) &
                    * real(this%ntimes_out, kind=dp) / real(this%ntimes_in, kind=dp)

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> resampling routine for rational number resampling, inspriation taken from
!! scipy.signal.resample
subroutine resample_timeshift(this, data_in, data_out, sources, ts_fwd_sample)

  use source_class

  class(resampling_type)                :: this
  real(kind=dp), intent(in)             :: data_in(:,:)
  real(kind=dp), intent(out)            :: data_out(:,:)
  class(src_param_type)                 :: sources(:)
  real(kind=dp), intent(in), optional   :: ts_fwd_sample

  complex(kind=dp), allocatable     :: dataf_in(:,:), dataf_out(:,:)
  complex(kind=dp), allocatable     :: shift_fd(:)
  integer                           :: nomega_min, i

  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing resampling type that is not initialized'
     call pabort 
  end if

  if (any(shape(data_in) /= [this%ntimes_in, this%ntraces])) then
     write(*,*) 'ERROR: shape mismatch in first argument - ' &
                    // 'shape does not match the plan for resampling'
     write(*,*) 'is: ', shape(data_in), '; should be: [', this%ntimes_in, &
                    ', ', this%ntraces, ']'
     call pabort
  end if

  if (any(shape(data_out) &
        /= [this%ntimes_out, this%ntraces])) then
     write(*,*) 'ERROR: shape mismatch in second argument - ' &
                    // 'shape does not match the plan for resampling'
     write(*,*) 'is: ', shape(data_out), '; should be: [', this%ntimes_out, &
                    ', ', this%ntraces, ']'
     call pabort
  end if

  if (size(sources) /= this%ntraces) then
     write(*,*) 'ERROR: shape mismatch for argument sources'
     call pabort
  end if

  allocate(dataf_in(1:this%fft_in%get_nomega(), 1:this%fft_in%get_ntraces()))
  allocate(dataf_out(1:this%fft_out%get_nomega(), 1:this%fft_out%get_ntraces()))
  allocate(shift_fd(1:this%fft_in%get_nomega()))

  call this%fft_in%rfft(taperandzeropad(data_in, this%fft_in%get_ntimes(), ntaper=0), dataf_in)
  
  ! time shift (dt in fft is 1!!)
  do i = 1, this%fft_in%get_ntraces()
     if (present(ts_fwd_sample)) then
        shift_fd = exp(- this%fft_in%get_f() * 2 * pi &
                       * (sources(i)%shift_time_sample - ts_fwd_sample) * cmplx(0, 1))
     else
        shift_fd = exp(- this%fft_in%get_f() * 2 * pi &
                       * sources(i)%shift_time_sample * cmplx(0, 1))
     endif
     dataf_in(:,i) = dataf_in(:,i) * shift_fd
     if (this%apply_filter) &
        dataf_in(:,i) = dataf_in(:,i) * this%filter%transferfunction
  enddo

  ! resample
  nomega_min = min(this%fft_in%get_nomega(), this%fft_out%get_nomega())
  dataf_out(1:nomega_min,:) = dataf_in(1:nomega_min,:)

  if (nomega_min < this%fft_out%get_nomega()) &
     dataf_out(nomega_min+1:,:) = 0

  call this%fft_out%irfft(dataf_out, this%data_out_buf)

  data_out(:,:) = this%data_out_buf(1:this%ntimes_out,:) &
                    * real(this%ntimes_out, kind=dp) / real(this%ntimes_in, kind=dp)

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine freeme(this)
  class(resampling_type) :: this

  call this%fft_in%freeme
  call this%fft_out%freeme
  this%initialized = .false.
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntimes_in(this)
  class(resampling_type)        :: this
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing resampling type that is not initialized'
     call pabort 
  end if
  get_ntimes_in = this%ntimes_in
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntimes_out(this)
  class(resampling_type)        :: this
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing resampling type that is not initialized'
     call pabort 
  end if
  get_ntimes_out = this%ntimes_out
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntraces(this)
  class(resampling_type)        :: this
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing resampling type that is not initialized'
     call pabort 
  end if
  get_ntraces = this%ntraces
end function
!-----------------------------------------------------------------------------------------


end module resampling
!=========================================================================================
