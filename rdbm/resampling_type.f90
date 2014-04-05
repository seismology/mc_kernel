!=========================================================================================
module resampling
  use global_parameters
  use commpi,               only : pabort
  use fft,                  only : rfft_type

  implicit none

  private
  public :: resampling_type

  type :: resampling_type
     private
     type(rfft_type)            :: fft_in, fft_out
     logical                    :: initialized = .false.
     contains
     procedure, pass            :: init
     procedure, pass            :: resample
     procedure, pass            :: freeme
     procedure, pass            :: get_ntraces
     procedure, pass            :: get_ntimes_in
     procedure, pass            :: get_ntimes_out
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine init(this, ntimes_in, ntimes_out, ntraces, measure)

  class(resampling_type)    :: this
  integer, intent(in)       :: ntimes_in, ntimes_out, ntraces
  logical, intent(in)       :: measure

  ! hardcoding ndim=1 for now, may be dummy variable in the future
  integer                   :: ndim 
  ndim = 1

  call this%fft_in%init(ntimes_in, ndim, ntraces, measure=measure)
  call this%fft_out%init(ntimes_out, ndim, ntraces, measure=measure)

  this%initialized = .true.
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> resampling routine for rational number resampling, inspriation taken from
!! scipy.signal.resample
subroutine resample(this, data_in, data_out)

  class(resampling_type)            :: this
  real(kind=dp), intent(in)         :: data_in(:,:), data_out(:,:)

  complex(kind=dp), allocatable     :: dataf_in(:,:), dataf_out(:,:)
  integer                           :: nomega_min

  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing resampling type that is not initialized'
     call pabort 
  end if

  if (any(shape(data_in) /= [this%fft_in%get_ntimes(), this%fft_in%get_ntraces()])) then
     write(*,*) 'ERROR: shape mismatch in first argument - ' &
                    // 'shape does not match the plan for resampling'
     write(*,*) 'is: ', shape(data_in), '; should be: [', this%fft_in%get_ntimes(), &
                    ', ', this%fft_in%get_ntraces(), ']'
     call pabort
  end if

  if (any(shape(data_out) &
        /= [this%fft_out%get_ntimes(), this%fft_out%get_ntraces()])) then
     write(*,*) 'ERROR: shape mismatch in second argument - ' &
                    // 'shape does not match the plan for resampling'
     write(*,*) 'is: ', shape(data_out), '; should be: [', this%fft_out%get_ntimes(), &
                    ', ', this%fft_out%get_ntraces(), ']'
     call pabort
  end if

  allocate(dataf_in(1:this%fft_in%get_nomega(), 1:this%fft_in%get_ntraces()))
  allocate(dataf_out(1:this%fft_out%get_nomega(), 1:this%fft_out%get_ntraces()))

  call this%fft_in%rfft(data_in, dataf_in)
  
  nomega_min = min(this%fft_in%get_nomega(), this%fft_out%get_nomega())

  dataf_out(1:nomega_min,:) = dataf_in(1:nomega_min,:)

  if (nomega_min < this%fft_out%get_nomega()) &
     dataf_out(nomega_min+1:,:) = 0

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
  get_ntimes_in = this%fft_in%get_ntimes()
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntimes_out(this)
  class(resampling_type)        :: this
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing resampling type that is not initialized'
     call pabort 
  end if
  get_ntimes_out = this%fft_out%get_ntimes()
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntraces(this)
  class(resampling_type)        :: this
  if (.not. this%initialized) then
     write(*,*) 'ERROR: accessing resampling type that is not initialized'
     call pabort 
  end if
  get_ntraces = this%fft_in%get_ntraces()
end function
!-----------------------------------------------------------------------------------------


end module resampling
!=========================================================================================
