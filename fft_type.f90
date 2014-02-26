!=========================================================================================
module fft
  use global_parameters
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f'

  private
  public :: rfft_type
  public :: taperandzeropad

  type :: rfft_type
     private
     integer(8)                 :: plan_fft, plan_ifft
     integer                    :: nomega, ntimes, ntraces
     real(kind=dp)              :: dt, df
     logical                    :: initialized = .false.
     integer                    :: fftw_mode = FFTW_ESTIMATE
     real(kind=dp), allocatable :: t(:), f(:)
     contains
     procedure, pass            :: get_nomega
     procedure, pass            :: get_ntimes
     procedure, pass            :: get_ntraces
     procedure, pass            :: get_initialized
     procedure, pass            :: get_f
     procedure, pass            :: get_df
     procedure, pass            :: get_t
     procedure, pass            :: init
     procedure, pass            :: rfft
     procedure, pass            :: irfft
     procedure, pass            :: freeme
  end type

contains

!-----------------------------------------------------------------------------------------
integer function get_nomega(this)
  class(rfft_type)        :: this
  if (.not. this%initialized) &
     stop 'ERROR: accessing fft type that is not initialized'
  get_nomega = this%nomega
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntimes(this)
  class(rfft_type)        :: this
  if (.not. this%initialized) &
     stop 'ERROR: accessing fft type that is not initialized'
  get_ntimes = this%ntimes
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntraces(this)
  class(rfft_type)        :: this
  if (.not. this%initialized) &
     stop 'ERROR: accessing fft type that is not initialized'
  get_ntraces = this%ntraces
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
logical function get_initialized(this)
  class(rfft_type)        :: this
  get_initialized = this%initialized
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_df(this)
   class(rfft_type) :: this
   real(kind=dp)    :: get_df
   if (.not. this%initialized) &
      stop 'ERROR: accessing fft type that is not initialized'
   get_df = this%df
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_f(this)
   class(rfft_type) :: this
   real(kind=dp)    :: get_f(this%nomega)
   if (.not. this%initialized) &
      stop 'ERROR: accessing fft type that is not initialized'
   get_f(:) = this%f(:)
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_t(this)
   class(rfft_type) :: this
   real(kind=dp)    :: get_t(this%ntimes)
   if (.not. this%initialized) &
      stop 'ERROR: accessing fft type that is not initialized'
   get_t(:) = this%t(:)
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init(this, ntimes_in, ntraces, dt)
  class(rfft_type)                     :: this
  integer, intent(in)                  :: ntimes_in, ntraces
  real(kind=dp), intent(in), optional  :: dt
  integer                              :: ntimes_np2, nomega_np2, i
  integer                              :: rank, istride, ostride, np2

  real(kind=dp), dimension(:,:), allocatable       :: datat
  complex(kind=dp), dimension(:,:), allocatable    :: dataf

  np2 = nextpow2(ntimes_in)
  nomega_np2 = np2 + 1
  ntimes_np2 = np2 * 2

  this%ntimes = ntimes_np2
  this%nomega = nomega_np2
  this%ntraces = ntraces
 
  ! determine spectral resolution
  if (present(dt)) then
     this%dt = dt
  else
     this%dt = 1.0
  end if

  this%df    = 1. / (this%dt*ntimes_np2)
  allocate( this%f( this%nomega ) )
  allocate( this%t( this%ntimes ) )
  
  do i = 1, this%nomega
     this%f(i) = (i-1) * this%df
  end do
  do i = 1, this%ntimes
     this%t(i) = i * this%dt
  end do

  ! temporaryly allocating work arrays, are needed in generation of the plan
  allocate(datat(1:ntimes_np2, 1:ntraces))
  allocate(dataf(1:nomega_np2, 1:ntraces))

  rank    = 1
  istride = 1
  ostride = 1

  call dfftw_plan_many_dft_r2c(this%plan_fft, rank, ntimes_np2, ntraces, datat, &
                               ntraces, istride, ntimes_np2, dataf, &
                               ntraces, ostride, nomega_np2, this%fftw_mode)

  call dfftw_plan_many_dft_c2r(this%plan_ifft, rank, ntimes_np2, ntraces, dataf, &
                               ntraces, istride, nomega_np2, datat, &
                               ntraces, ostride, ntimes_np2, this%fftw_mode)


  this%initialized = .true.

  deallocate(datat)
  deallocate(dataf)

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine rfft(this, datat, dataf)
  class(rfft_type) :: this
  real(kind=dp), intent(in)        :: datat(:,:)
  complex(kind=dp), intent(out)    :: dataf(:,:)

  if (.not. this%initialized) &
     stop 'ERROR: trying fft on a plan that was not initialized'

  if (any(shape(datat) /= (/this%ntimes, this%ntraces/))) then
     write(*,*) 'ERROR: shape mismatch in first argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(datat), '; should be: (', this%ntimes, ', ', this%ntraces, ')'
     stop
  end if

  if (any(shape(dataf) /= (/this%nomega, this%ntraces/))) then
     write(*,*) 'ERROR: shape mismatch in second argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(dataf), '; should be: (', this%nomega, ', ', this%ntraces, ')'
     stop
  end if

  ! use specific interfaces including the buffer arrays to make sure the
  ! compiler does not change the order of execution (stupid, but known issue)
  call dfftw_execute_dft_r2c(this%plan_fft, datat, dataf)
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine irfft(this, dataf, datat)
  class(rfft_type)              :: this
  complex(kind=dp), intent(in)  :: dataf(:,:)
  real(kind=dp), intent(out)    :: datat(:,:)

  if (.not. this%initialized) then
     stop 'ERROR: trying inverse fft on a plan that was not initialized'
  endif

  if (any(shape(dataf) /= (/this%nomega, this%ntraces/))) then
     write(*,*) 'ERROR: shape mismatch in first argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(dataf), '; should be: (', this%nomega, ', ', this%ntraces, ')'
     stop 
  end if

  if (any(shape(datat) /= (/this%ntimes, this%ntraces/))) then
     write(*,*) 'ERROR: shape mismatch in second argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(datat), '; should be: (', this%ntimes, ', ', this%ntraces, ')'
     stop 
  end if

  ! use specific interfaces including the buffer arrays to make sure the
  ! compiler does not change the order of execution (stupid, but known issue)
  call dfftw_execute_dft_c2r(this%plan_ifft, dataf, datat)

  ! normalization, see
  ! http://www.fftw.org/doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html
  datat = datat / this%ntimes
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine freeme(this)
  class(rfft_type) :: this

  call dfftw_destroy_plan(this%plan_fft)
  call dfftw_destroy_plan(this%plan_ifft)
  this%initialized = .false.
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure integer function nextpow2(n)
  integer, intent(in) :: n

  nextpow2 = 2
  do while (nextpow2 < n) 
     nextpow2 = nextpow2 * 2
  end do
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function taperandzeropad(array, ntimes)
  real(kind=dp), intent(in)     :: array(:,:)
  integer,       intent(in)     :: ntimes
  real(kind=dp)                 :: taperandzeropad(ntimes, size(array,2))
  real(kind=dp), allocatable    :: win(:,:)

  integer                       :: ntaper ! Taper length in samples
  real, parameter               :: D=3.3  ! Decay constant
  integer                       :: ntimes_in ! Length of incoming time series
  integer                       :: i
  
  ntimes_in = size(array,1)
  if (ntimes<ntimes_in) then
     write(*,*) 'For zeropadding to work, ntimes has to be equal or larger ', & 
                'than the length of the input array'
     stop
  end if

  ntaper   = ntimes_in / 4

  allocate(win(ntimes_in, size(array,2)))
  win = 1
  do i = 1, ntaper
     win(i,:) = exp( -(D * (ntaper-i+1) / ntaper)**2 )
  end do
  do i = 0, ntaper-1
     win(ntimes_in - i,:) = exp( -(D * (ntaper-i) / ntaper)**2 )
  end do
  taperandzeropad(:,:) = 0
  taperandzeropad(1:size(array,1), :) = array * win
  deallocate(win)

  
end function taperandzeropad
!-----------------------------------------------------------------------------------------

end module fft
!=========================================================================================
