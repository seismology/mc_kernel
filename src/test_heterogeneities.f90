!=========================================================================================
module test_heterogeneities

  use global_parameters
  use heterogeneities
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_load_het_rtpv 

  type(hetero_type)     :: het_model
  character(len=64)     :: filename
  integer               :: ix, iy, iz
  real(kind=dp)         :: coeffs_ref, coeffs(3,1), coords(3,1)

  filename = './test_het.rtpv'

  call het_model%load_het_rtpv(filename)

  !call assert_equal(het_model%nhet, 8, 'Correct number of points in heterogeneity file')

  call het_model%build_hetero_kdtree()

  do ix = 0, 1
    do iy = 0, 1
      do iz = 0, 1
        coords(:,1) = real([ix, iy, iz], kind=dp) * 1d3
        coeffs = het_model%load_model_coeffs(coords, [1.0d0]) 
        coeffs_ref = ix + iy*2 + iz*4
        print *, coeffs(1,1), coeffs_ref, ix, iy, iz
        call assert_comparable(coeffs(1,1),   coeffs_ref, 1d-10, 'vp is #index')
        call assert_comparable(coeffs(2,1),  -coeffs_ref, 1d-10, 'vs is -#index')
        call assert_comparable(coeffs(3,1), 2*coeffs_ref, 1d-10, 'rho is 2*#index')
      end do
    end do
  end do

end subroutine test_load_het_rtpv 
!-----------------------------------------------------------------------------------------

end module test_heterogeneities
!=========================================================================================
