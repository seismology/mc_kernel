!=========================================================================================
module test_source

   use global_parameters, only: sp, pi
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

   call assert_comparable_real1d(real(sources(1)%mij_voigt, sp) + 1e22, &
                                 [0., 0., 0., 0., 0., -0.64e22] + 1e22, 1e-5, &
                                 'strike = 0, rake = 0, dip = 90')

   call assert_comparable_real1d(real(sources(2)%mij_voigt, sp) + 1e22, &
                                 [0.64e22, -0.64e22, 0., 0., 0., 0.] + 1e22, 1e-5, &
                                 'strike = 0, rake = -45, dip = 90')

   call assert_comparable_real1d(real(sources(3)%mij_voigt, sp) + 1e22, &
                                 [0., 0., 0., 0., -0.64e22, 0.] + 1e22, 1e-5, &
                                 'strike = 0, rake = 0, dip = 0')

   do i=1, nsources
      call assert_comparable(real(sum(sources(i)%mij_voigt(1:3)), sp) + 1e22, 1e22, 1e-5, &
                             'trace of DC is zero')
   enddo

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
