!=========================================================================================
module test_lanczos

   use global_parameters, only: sp, dp, pi
   use lanczos
   use ftnunit
   implicit none
   public
contains

!-----------------------------------------------------------------------------------------
subroutine test_lanczos_resample

  integer, parameter        :: ntimes = 1000, ncoarsen = 4
  real(kind=dp), parameter  :: thalf  = 100.d0
  real(kind=dp)             :: data_in(ntimes)
  real(kind=dp), allocatable:: data_resampled(:)
  real(kind=dp), allocatable:: data_in_coarse(:)
  integer                   :: x

  data_in(:) = 0
  data_in(ntimes/4+1:ntimes) = [(2.d0 * exp(-(x/thalf)**2)*x/thalf**2, x=1, ntimes-ntimes/4 )]

  data_in_coarse = data_in(1:ntimes:ncoarsen)

  data_resampled = lanczos_resample(si      = data_in_coarse,          &
                                    dt_old  = real(ncoarsen, kind=dp), &
                                    dt_new  = 1.0d0,                   &
                                    a       = 12)


  call assert_equal(size(data_in), size(data_resampled), 'Size of resampled array is as expected')

  call assert_comparable(data_in+1d0, data_resampled+1d0, 1d-4, &
                         'Resampled data is equal input data')
  

end subroutine test_lanczos_resample
!-----------------------------------------------------------------------------------------

end module test_lanczos      
!=========================================================================================
