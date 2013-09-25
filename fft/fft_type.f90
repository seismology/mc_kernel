!=========================================================================================
module fft
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  private
  public :: rfft_type

  type :: rfft_type
     private
     integer(8)  :: plan_fft, plan_ifft
     integer     :: nomega, ntimes, ntraces
     logical     :: initialized = .false.
     integer     :: fftw_mode = FFTW_ESTIMATE
     contains
     procedure, pass :: get_nomega
     procedure, pass :: get_ntimes
     procedure, pass :: get_ntraces
     procedure, pass :: get_initialized
     procedure, pass :: init
     procedure, pass :: rfft
     procedure, pass :: irfft
     procedure, pass :: freeme
  end type

contains

!-----------------------------------------------------------------------------------------
integer function get_nomega(this)
  class(rfft_type)        :: this
  get_nomega = this%nomega
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntimes(this)
  class(rfft_type)        :: this
  get_ntimes = this%ntimes
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntraces(this)
  class(rfft_type)        :: this
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
subroutine init(this, ntimes_in, ntraces)
  class(rfft_type)        :: this
  integer, intent(in)     :: ntimes_in, ntraces
  integer                 :: ntimes_np2, nomega_np2
  integer                 :: rank, istride, ostride, np2

  real(kind=8), dimension(:,:), allocatable       :: datat
  complex(kind=8), dimension(:,:), allocatable    :: dataf

  np2 = nextpow2(ntimes_in)
  nomega_np2 = np2 + 1
  ntimes_np2 = np2 * 2

  this%ntimes = ntimes_np2
  this%nomega = nomega_np2
  this%ntraces = ntraces

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
  real(kind=8), intent(in)        :: datat(:,:)
  complex(kind=8), intent(out)    :: dataf(:,:)

  if (.not. this%initialized) &
     stop 'ERROR: trying fft on a plan that was not initialized'

  if (any(shape(datat) /= (/this%ntimes, this%ntraces/))) &
     stop 'ERROR: shape mismatch in first argument - shape does not match the plan for fftw'

  if (any(shape(dataf) /= (/this%nomega, this%ntraces/))) &
     stop 'ERROR: shape mismatch in second argument - shape does not match the plan for fftw'

  ! use specific interfaces including the buffer arrays to make sure the
  ! compiler does not change the order of execution (stupid, but known issue)
  call dfftw_execute_dft_r2c(this%plan_fft, datat, dataf)
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine irfft(this, dataf, datat)
  class(rfft_type) :: this
  complex(kind=8), intent(in)     :: dataf(:,:)
  real(kind=8), intent(out)       :: datat(:,:)

  if (.not. this%initialized) then
     stop 'ERROR: trying inverse fft on a plan that was not initialized'
  endif

  if (any(shape(dataf) /= (/this%nomega, this%ntraces/))) &
     stop 'ERROR: shape mismatch in first argument - shape does not match the plan for fftw'

  if (any(shape(datat) /= (/this%ntimes, this%ntraces/))) &
     stop 'ERROR: shape mismatch in second argument - shape does not match the plan for fftw'

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
  this%initialized = .true.
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

end module fft
!=========================================================================================
