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

  ! hardcoding ndim for now, may be dummy variable in the future
  integer                   :: ndim 
  ndim = 1

  call this%fft_in%init(ntimes_in, ndim, ntraces, measure=measure)
  call this%fft_out%init(ntimes_out, ndim, ntraces, measure=measure)

  this%initialized = .true.
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
