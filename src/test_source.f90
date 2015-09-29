!=========================================================================================
module test_source

   use global_parameters, only: dp, sp, pi
   use source_class
   use ftnunit
   implicit none
   public
contains

!!-----------------------------------------------------------------------------------------
!subroutine test_resample_stf
!   character(len=256)   :: srf_file
!   integer              :: npoints, nsources, i, lu
!   type(src_param_type), allocatable :: sources(:)
!   real(kind=dp), allocatable        :: stf_ref(:)
!
!   srf_file = 'unit_tests/standard_rupture_fomat_testfile.srf'
!   call read_srf(srf_file, sources, npoints=npoints, nsources=nsources)
!
!   open(newunit=lu, file='unit_tests/srf_stf_in')
!   write(lu,'(2E20.10)') (sources(1)%stf_dt * (i-1), sources(1)%stf(i), i=1, size(sources(1)%stf))
!   close(lu)
!
!   call sources(1)%resample_stf(sources(1)%stf_dt, size(sources(1)%stf))
!
!   call assert_comparable_real1d(real(sources(1)%stf, sp) + 1, &
!                                 real(sources(1)%stf_resampled, sp) + 1, 1e-5, &
!                                 'in = out')
!
!   call sources(1)%resample_stf(sources(1)%stf_dt / 5.1, size(sources(1)%stf) * 6)
!
!   allocate(stf_ref(size(sources(1)%stf) * 6))
!   open(newunit=lu, file='unit_tests/test_stf_resample')
!   read(lu,*) stf_ref
!   close(lu)
!
!   call assert_comparable_real1d(real(sources(1)%stf_resampled, sp) + 1, &
!                                 real(stf_ref, sp) + 1, 1e-5, &
!                                 'odd sampling')
!
!end subroutine
!!-----------------------------------------------------------------------------------------
!
!!-----------------------------------------------------------------------------------------
!subroutine test_fft_stf
!   character(len=256)   :: srf_file
!   integer              :: npoints, nsources, i, lu
!   type(src_param_type), allocatable :: sources(:)
!   real(kind=dp), allocatable        :: stf_ref(:)
!   real(kind=dp), allocatable        :: stf_fwd(:)
!
!   srf_file = 'unit_tests/standard_rupture_fomat_testfile.srf'
!   call read_srf(srf_file, sources, npoints=npoints, nsources=nsources)
!
!   do i=1, nsources
!      call sources(i)%resample_stf(sources(i)%stf_dt, 16)
!   enddo
!
!   ! only fft and compare to ref data
!   call fft_stf(sources)
!
!   !open(newunit=lu, file='unit_tests/test_stf_fft')
!   !write(lu,*) abs(sources(1)%stf_fd)
!   !close(lu)
!
!   allocate(stf_ref(17))
!   open(newunit=lu, file='unit_tests/test_stf_fft')
!   read(lu,*) stf_ref
!   close(lu)
!
!   call assert_comparable_real1d(real(abs(sources(1)%stf_fd), sp), &
!                                 real(stf_ref, sp), 1e-5, &
!                                 'stf fft')
!
!   ! stf_fwd = dirac should have no effect in deconvolution
!   allocate(stf_fwd(16))
!   stf_fwd = 0
!   stf_fwd(1) = 1d0 / sources(1)%stf_dt_resampled
!
!   call fft_stf(sources, stf_fwd)
!
!   call assert_comparable_real1d(real(abs(sources(1)%stf_reconv_fd), sp), &
!                                 real(stf_ref, sp), 1e-4, &
!                                 'stf fft dirac')
!
!
!   ! stf_fwd = desired stf -> reconv stf = 1
!   stf_fwd = sources(1)%stf_resampled
!
!   call fft_stf(sources, stf_fwd)
!   stf_ref = 1
!
!   call assert_comparable_real1d(real(abs(sources(1)%stf_reconv_fd), sp), &
!                                 real(stf_ref, sp), 1e-5, &
!                                 'stf fft, fwd = desired')
!
!   ! stf_fwd = dirac + have_stf = false -> no effect
!   stf_fwd = 0
!   stf_fwd(1) = 1d0 / sources(1)%stf_dt_resampled
!
!   do i=1, nsources
!      sources(1)%have_stf = .false.
!      call sources(i)%resample_stf(sources(i)%stf_dt, 16)
!   enddo
!
!   call fft_stf(sources, stf_fwd)
!   stf_ref = 1
!
!   call assert_comparable_real1d(real(abs(sources(1)%stf_reconv_fd), sp), &
!                                 real(stf_ref, sp), 1e-5, &
!                                 'stf_fwd = dirac, have_stf = .false.')
!
!end subroutine
!!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
