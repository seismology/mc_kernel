!=========================================================================================
module test_source

   use global_parameters, only: dp, sp, pi
   use source_class
   use ftnunit
   implicit none
   public
contains

!-----------------------------------------------------------------------------------------
subroutine test_read_srf
   character(len=256)   :: srf_file
   integer              :: npoints, nsources, i
   type(src_param_type), allocatable :: sources(:)

   srf_file = 'unit_tests/standard_rupture_fomat_testfile.srf'
   call read_srf(srf_file, sources, npoints=npoints, nsources=nsources)

   call assert_equal(npoints, 10, 'Number of points')
   call assert_equal(nsources, 8, 'Number of sources')

   !do i=1, nsources
   !   write(6,'(6e10.2)') sources(i)%mij_voigt(:)
   !enddo

   call assert_comparable_real1d(real(sources(1)%mij_voigt, sp) + 1e16, &
                                 [0., 0., 0., 0., 0., -0.64e16] + 1e16, 1e-5, &
                                 'strike = 0, rake = 0, dip = 90')

   call assert_comparable_real1d(real(sources(2)%mij_voigt, sp) + 1e16, &
                                 [0.64e16, -0.64e16, 0., 0., 0., 0.] + 1e16, 1e-5, &
                                 'strike = 0, rake = -45, dip = 90')

   call assert_comparable_real1d(real(sources(3)%mij_voigt, sp) + 1e16, &
                                 [0., 0., 0., 0., -0.64e16, 0.] + 1e16, 1e-5, &
                                 'strike = 0, rake = 0, dip = 0')

   do i=1, nsources
      call assert_comparable(real(sum(sources(i)%mij_voigt(1:3)), sp) + 1e22, 1e22, 1e-5, &
                             'trace of DC is zero')
   enddo

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_resample_stf
   character(len=256)   :: srf_file
   integer              :: npoints, nsources, i, lu
   type(src_param_type), allocatable :: sources(:)
   real(kind=dp), allocatable        :: stf_ref(:)

   srf_file = 'unit_tests/standard_rupture_fomat_testfile.srf'
   call read_srf(srf_file, sources, npoints=npoints, nsources=nsources)

   open(newunit=lu, file='unit_tests/srf_stf_in')
   write(lu,'(2E20.10)') (sources(1)%stf_dt * (i-1), sources(1)%stf(i), i=1, size(sources(1)%stf))
   close(lu)

   call sources(1)%resample_stf(sources(1)%stf_dt, size(sources(1)%stf))

   call assert_comparable_real1d(real(sources(1)%stf, sp) + 1, &
                                 real(sources(1)%stf_resampled, sp) + 1, 1e-5, &
                                 'in = out')

   call sources(1)%resample_stf(sources(1)%stf_dt / 5.1, size(sources(1)%stf) * 6)

   allocate(stf_ref(size(sources(1)%stf) * 6))
   open(newunit=lu, file='unit_tests/test_stf_resample')
   read(lu,*) stf_ref
   close(lu)

   call assert_comparable_real1d(real(sources(1)%stf_resampled, sp) + 1, &
                                 real(stf_ref, sp) + 1, 1e-5, &
                                 'odd sampling')

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
