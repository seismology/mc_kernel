  module data_fft
  use global_parameters, only : realkind
  implicit none
  public
!-------------------------------------------------------------------------------------
  integer :: ntimes ! number of discrete time steps
  integer :: nomega ! number of discrete frequencies
  integer :: npad   ! number of extra discrete timesteps for padding 
!
  real(kind=8), dimension(:), allocatable :: times
  real(kind=8) :: timewidth
  real(kind=8), dimension(:), allocatable :: omega, period, frequency

! plan for 1D example
  real(kind=8), dimension(:), allocatable   :: uphys1d
  complex(kind=8),dimension(:), allocatable ::  uspec1d,v0spec
  integer*8 plan_fftf1d,plan_fftb1d ! plans for the one-dimensional forward and backward ffts

! plan for full kernel
  integer*8  plan_fftf_ker, plan_fftb_ker

! filtering
  complex(kind=8),allocatable :: specfilt(:)
  real(kind=8) :: filt_power
  integer :: filter_what
  character(len=3) :: filter_type
  real(kind=8) :: filter_period,filter_period_low,filter_period_hi
  real(kind=8) :: fmin,fmax
  real(kind=8) :: omegamin,omegamax
  
!-------------------------------------------------------------------------------------
  end module data_fft
