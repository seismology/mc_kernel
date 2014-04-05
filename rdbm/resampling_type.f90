!=========================================================================================
module resampling
  use global_parameters
  use fft

  implicit none

  private
  public :: resampling_type

  type :: resampling_type
     private
     type(rfft_type)            :: fft_in, fft_out
     contains
     procedure, pass            :: init
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
end subroutine
!-----------------------------------------------------------------------------------------

end module resampling
!=========================================================================================
