!=========================================================================================
module fft
  use global_parameters
  use simple_routines, only      : mult2d_1d, mult3d_1d
  use commpi,          only      : pabort
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f'

  private
  public :: rfft_type
  public :: taperandzeropad

  type :: rfft_type
     private
     integer(8)                 :: plan_fft, plan_ifft
     integer(8)                 :: plan_fft_1d, plan_ifft_1d
     integer                    :: nomega, ntimes, ntraces, ndim, ntraces_fft
     real(kind=dp)              :: dt, df
     logical                    :: initialized = .false.
     integer                    :: fftw_mode = FFTW_ESTIMATE
     real(kind=dp), allocatable :: t(:), f(:)
     contains
     procedure, pass            :: get_nomega
     procedure, pass            :: get_ntimes
     procedure, pass            :: get_ntraces
     procedure, pass            :: get_ndim
     procedure, pass            :: get_initialized
     procedure, pass            :: get_f
     procedure, pass            :: get_df
     procedure, pass            :: get_t
     procedure, pass            :: init
     procedure, pass            :: rfft_md
     procedure, pass            :: irfft_md
     procedure, pass            :: rfft_1d
     procedure, pass            :: irfft_1d
     generic                    :: rfft => rfft_md, rfft_1d
     generic                    :: irfft => irfft_md, irfft_1d
     procedure, pass            :: freeme
  end type

  interface taperandzeropad
     module procedure           :: taperandzeropad_1d
     module procedure           :: taperandzeropad_md
  end interface
contains

!-----------------------------------------------------------------------------------------
integer function get_nomega(this)
  class(rfft_type)        :: this
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing fft type that is not initialized'
     call pabort 
  end if
  get_nomega = this%nomega
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntimes(this)
  class(rfft_type)        :: this
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing fft type that is not initialized'
     call pabort 
  end if
  get_ntimes = this%ntimes
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntraces(this)
  class(rfft_type)        :: this
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing fft type that is not initialized'
     call pabort 
  end if
  get_ntraces = this%ntraces
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ndim(this)
  class(rfft_type)        :: this
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing fft type that is not initialized'
     call pabort 
  end if
  get_ndim = this%ndim
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
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing fft type that is not initialized'
     call pabort 
  end if
  get_df = this%df
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_f(this)
  class(rfft_type) :: this
  real(kind=dp)    :: get_f(this%nomega)
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing fft type that is not initialized'
     call pabort 
  end if
  get_f(:) = this%f(:)
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_t(this)
  class(rfft_type) :: this
  real(kind=dp)    :: get_t(this%ntimes)
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing fft type that is not initialized'
     call pabort 
  end if
  get_t(:) = this%t(:)
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init(this, ntimes_in, ndim, ntraces, dt, measure)
!< This routines initializes a FFTw plan for 1D DFTs. 
!! The time series is supposed to have length ntimes_in, and is stored along the first 
!! dimension of a three-dimensional array. The other two dimensions are ndim and ntraces.
!! Internally, this array will always be converted to two-dimensional array of dimensions
!! [ntimes_in, ndim*ntraces]
  class(rfft_type)                     :: this
  integer, intent(in)                  :: ntimes_in, ntraces, ndim
  real(kind=dp), intent(in), optional  :: dt
  logical, intent(in), optional        :: measure 
  integer                              :: ntimes_np2, nomega_np2, i, ntraces_fft
  integer                              :: rank, istride, ostride, np2

  real(kind=dp), dimension(:,:), allocatable       :: datat
  complex(kind=dp), dimension(:,:), allocatable    :: dataf

  np2 = nextpow2(ntimes_in)
  nomega_np2 = np2 + 1
  ntimes_np2 = np2 * 2

  ntraces_fft = ntraces * ndim

  this%ntimes      = ntimes_np2
  this%nomega      = nomega_np2
  this%ndim        = ndim
  this%ntraces     = ntraces 
  this%ntraces_fft = ntraces_fft 

  this%fftw_mode = FFTW_ESTIMATE
  if (present(measure)) then
      if (measure) then
          this%fftw_mode = FFTW_MEASURE
      end if
  end if


 
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
  ! Creating plan for multi-dimensional FFT (is needed for forward FT)
  allocate(datat(1:ntimes_np2, 1:ntraces_fft))
  allocate(dataf(1:nomega_np2, 1:ntraces_fft))

  rank    = 1
  istride = 1
  ostride = 1

  call dfftw_plan_many_dft_r2c(this%plan_fft, rank, ntimes_np2, ntraces_fft, datat, &
                               ntraces_fft, istride, ntimes_np2, dataf, &
                               ntraces_fft, ostride, nomega_np2, this%fftw_mode)

  call dfftw_plan_many_dft_c2r(this%plan_ifft, rank, ntimes_np2, ntraces_fft, dataf, &
                               ntraces_fft, istride, nomega_np2, datat, &
                               ntraces_fft, ostride, ntimes_np2, this%fftw_mode)


  deallocate(datat)
  deallocate(dataf)

  ! Creating plan for one-dimensional FFT (is needed in backward FT)
  allocate(datat(1:ntimes_np2, 1:ntraces))
  allocate(dataf(1:nomega_np2, 1:ntraces))

  rank    = 1
  istride = 1
  ostride = 1

  call dfftw_plan_many_dft_r2c(this%plan_fft_1d, rank, ntimes_np2, ntraces, datat, &
                               ntraces, istride, ntimes_np2, dataf, &
                               ntraces, ostride, nomega_np2, this%fftw_mode)

  call dfftw_plan_many_dft_c2r(this%plan_ifft_1d, rank, ntimes_np2, ntraces, dataf, &
                               ntraces, istride, nomega_np2, datat, &
                               ntraces, ostride, ntimes_np2, this%fftw_mode)


  deallocate(datat)
  deallocate(dataf)

  this%initialized = .true.

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine rfft_1d(this, datat, dataf)
  class(rfft_type) :: this
  real(kind=dp), intent(in)        :: datat(:,:)
  complex(kind=dp), intent(out)    :: dataf(:,:)
  !real(kind=dp)                    :: datat(this%ntimes, this%ntraces)
  !complex(kind=dp)                 :: dataf(this%nomega, this%ntraces)

  if (.not. this%initialized) then
     write(*,*) 'ERROR: trying fft on a plan that was not initialized'
     call pabort 
  end if

  if (any(shape(datat) /= [this%ntimes, this%ntraces])) then
     write(*,*) 'ERROR: shape mismatch in first argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(datat), '; should be: [', this%ntimes, ', ', this%ntraces, ']'
     call pabort
  end if

  if (any(shape(dataf) /= [this%nomega, this%ntraces])) then
     write(*,*) 'ERROR: shape mismatch in second argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(dataf), '; should be: [', this%nomega, ', ', this%ntraces, ']'
     call pabort
  end if

  ! use specific interfaces including the buffer arrays to make sure the
  ! compiler does not change the order of execution (stupid, but known issue)
  call dfftw_execute_dft_r2c(this%plan_fft_1d, datat, dataf)

end subroutine rfft_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine irfft_1d(this, dataf, datat)
  class(rfft_type)              :: this
  complex(kind=dp), intent(in)  :: dataf(:,:)
  real(kind=dp), intent(out)    :: datat(:,:)

  if (.not. this%initialized) then
     write(*,*) 'ERROR: trying fft on a plan that was not initialized'
     call pabort 
  end if

  if (any(shape(dataf) /= (/this%nomega, this%ntraces/))) then
     write(*,*) 'ERROR: shape mismatch in first argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(dataf), '; should be: [', this%nomega, ', ', this%ntraces, ']'
     call pabort 
  end if

  if (any(shape(datat) /= (/this%ntimes, this%ntraces/))) then
     write(*,*) 'ERROR: shape mismatch in second argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(datat), '; should be: [', this%ntimes, ', ', this%ntraces, ']'
     call pabort 
  end if

  ! use specific interfaces including the buffer arrays to make sure the
  ! compiler does not change the order of execution (stupid, but known issue)
  call dfftw_execute_dft_c2r(this%plan_ifft, dataf, datat)
  ! normalization, see
  ! http://www.fftw.org/doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html
  datat = datat / this%ntimes

end subroutine irfft_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine rfft_md(this, datat_in, dataf_out)
  class(rfft_type) :: this
  real(kind=dp), intent(in)        :: datat_in(:,:,:)
  complex(kind=dp), intent(out)    :: dataf_out(:,:,:)
  real(kind=dp)                    :: datat(this%ntimes, this%ntraces_fft)
  complex(kind=dp)                 :: dataf(this%nomega, this%ntraces_fft)

  if (.not. this%initialized) then
     write(*,*) 'ERROR: trying fft on a plan that was not initialized'
     call pabort 
  end if

  if (any(shape(datat_in) /= [this%ntimes, this%ndim, this%ntraces])) then
     write(*,*) 'ERROR: shape mismatch in first argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(datat_in), '; should be: [', this%ntimes, ', ', this%ndim, ', ', this%ntraces, ']'
     call pabort
  end if

  if (any(shape(dataf_out) /= [this%nomega, this%ndim, this%ntraces])) then
     write(*,*) 'ERROR: shape mismatch in second argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(dataf_out), '; should be: [', this%nomega, ', ', this%ndim, ', ', this%ntraces, ']'
     call pabort
  end if

  ! use specific interfaces including the buffer arrays to make sure the
  ! compiler does not change the order of execution (stupid, but known issue)
  datat = reshape(datat_in, [this%ntimes, this%ntraces_fft])
  call dfftw_execute_dft_r2c(this%plan_fft, datat, dataf)
  dataf_out = reshape(dataf, [this%nomega, this%ndim, this%ntraces])

end subroutine rfft_md
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine irfft_md(this, dataf_in, datat_out)
  class(rfft_type)              :: this
  complex(kind=dp), intent(in)  :: dataf_in(:,:,:)
  real(kind=dp), intent(out)    :: datat_out(:,:,:)
  real(kind=dp)                 :: datat(this%ntimes, this%ntraces_fft)
  complex(kind=dp)              :: dataf(this%nomega, this%ntraces_fft)


  if (.not. this%initialized) then
     write(*,*) 'ERROR: trying fft on a plan that was not initialized'
     call pabort 
  end if

  if (any(shape(dataf_in) /= (/this%nomega, this%ndim, this%ntraces/))) then
     write(*,*) 'ERROR: shape mismatch in first argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(dataf_in), '; should be: [', this%nomega, ', ', this%ndim, ', ', this%ntraces, ']'
     call pabort 
  end if

  if (any(shape(datat_out) /= (/this%ntimes, this%ndim, this%ntraces/))) then
     write(*,*) 'ERROR: shape mismatch in second argument - shape does not match the plan for fftw'
     write(*,*) 'is: ', shape(datat_out), '; should be: [', this%ntimes, ', ', this%ndim, ', ', this%ntraces, ']'
     call pabort 
  end if

  ! use specific interfaces including the buffer arrays to make sure the
  ! compiler does not change the order of execution (stupid, but known issue)
  dataf = reshape(dataf_in, [this%nomega, this%ntraces_fft])
  call dfftw_execute_dft_c2r(this%plan_ifft, dataf, datat)
  ! normalization, see
  ! http://www.fftw.org/doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html
  datat = datat / this%ntimes

  datat_out = reshape(datat, [this%ntimes, this%ndim, this%ntraces]) 

end subroutine irfft_md
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
function taperandzeropad_1d(array, ntimes, ntaper)
  real(kind=dp), intent(in)     :: array(:,:)
  integer,       intent(in)     :: ntimes
  integer, intent(in), optional :: ntaper !< Taper length in samples, 
                                          !! default is 1/4 of signal length
  real(kind=dp)                 :: taperandzeropad_1d(ntimes, size(array,2))
  real(kind=dp), allocatable    :: win(:)

  real, parameter               :: D=3.3  ! Decay constant
  integer                       :: ntaper_loc
  integer                       :: ntimes_in ! Length of incoming time series
  integer                       :: i
  
  ntimes_in = size(array,1)
  if (ntimes<ntimes_in) then
     write(*,*) 'For zeropadding to work, ntimes has to be equal or larger ', & 
                'than the length of the input array'
     call pabort
  end if

  if (present(ntaper)) then
     ntaper_loc = ntaper
  else
     ntaper_loc = ntimes_in / 4
  end if

  allocate(win(ntimes_in))
  win = 1
  do i = 1, ntaper_loc
     win(i) = exp( -(D * (ntaper_loc-i+1) / ntaper_loc)**2 )
  end do
  do i = 0, ntaper_loc-1
     win(ntimes_in - i) = exp( -(D * (ntaper_loc-i) / ntaper_loc)**2 )
  end do
  taperandzeropad_1d(:,:) = 0
  taperandzeropad_1d(1:size(array,1),:) = mult2d_1d(array, win)
  deallocate(win)

  
end function taperandzeropad_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function taperandzeropad_md(array, ntimes, ntaper)
  real(kind=dp), intent(in)     :: array(:,:,:)
  integer,       intent(in)     :: ntimes
  integer, intent(in), optional :: ntaper !< Taper length in samples, 
                                          !! default is 1/4 of signal length
  real(kind=dp)                 :: taperandzeropad_md(ntimes, size(array,2), size(array,3))
  real(kind=dp), allocatable    :: win(:)

  real, parameter               :: D=3.3  ! Decay constant
  integer                       :: ntaper_loc
  integer                       :: ntimes_in ! Length of incoming time series
  integer                       :: i
  
  ntimes_in = size(array,1)
  if (ntimes<ntimes_in) then
     write(*,*) 'For zeropadding to work, ntimes has to be equal or larger ', & 
                'than the length of the input array'
     call pabort
  end if

  if (present(ntaper)) then
     ntaper_loc = ntaper
  else
     ntaper_loc = ntimes_in / 4
  end if

  allocate(win(ntimes_in))
  win = 1
  do i = 1, ntaper_loc
     win(i) = exp( -(D * (ntaper_loc-i+1) / ntaper_loc)**2 )
  end do
  do i = 0, ntaper_loc-1
     win(ntimes_in - i) = exp( -(D * (ntaper_loc-i) / ntaper_loc)**2 )
  end do
  taperandzeropad_md(:,:,:) = 0
  taperandzeropad_md(1:size(array,1), :, :) = mult3d_1d(array, win)
  deallocate(win)

  
end function taperandzeropad_md
!-----------------------------------------------------------------------------------------

end module fft
!=========================================================================================
