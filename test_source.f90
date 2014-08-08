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
   integer              :: npoints, nsources
   type(src_param_type), allocatable :: sources(:)

   srf_file = 'standard_rupture_fomat_testfile.srf'
   call read_srf(srf_file, sources, npoints=npoints, nsources=nsources)

   call assert_equal(npoints, 10, 'Number of points')
   call assert_equal(nsources, 8, 'Number of sources')

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
