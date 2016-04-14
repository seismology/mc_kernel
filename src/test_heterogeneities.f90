!******************************************************************************
!
!    This file is part of:
!    MC Kernel: Calculating seismic sensitivity kernels on unstructured meshes
!    Copyright (C) 2016 Simon Staehler, Martin van Driel, Ludwig Auer
!
!    You can find the latest version of the software at:
!    <https://www.github.com/tomography/mckernel>
!
!    MC Kernel is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MC Kernel is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MC Kernel. If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

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
  integer               :: ix, iy, iz, ipoint
  integer, parameter    :: npoints = 8
  real(kind=dp)         :: coeffs_ref, coeffs(3,npoints), coords(3,npoints)
  real(kind=dp)         :: weights(npoints)

  filename = './test_het.rtpv'

  call het_model%load_het_rtpv(filename)

  call het_model%build_hetero_kdtree()

  ipoint = 1 
  do ix = 0, 1
    do iy = 0, 1
      do iz = 0, 1
        coords(:,ipoint) = real([ix, iy, iz], kind=dp) * 1d3
        ipoint = ipoint + 1
      end do
    end do
  end do

  weights = 1
  
  coeffs = het_model%load_model_coeffs(coords, weights) 

  ipoint = 1 
  do ix = 0, 1
    do iy = 0, 1
      do iz = 0, 1
        coeffs_ref = ix + iy*2 + iz*4
        call assert_comparable(coeffs(1,ipoint),   coeffs_ref, 1d-10, 'vp is #index')
        call assert_comparable(coeffs(2,ipoint),  -coeffs_ref, 1d-10, 'vs is -#index')
        call assert_comparable(coeffs(3,ipoint), 2*coeffs_ref, 1d-10, 'rho is 2*#index')
        ipoint = ipoint + 1
      end do
    end do
  end do

  ! Repeat test, but with individual weights

  call random_number(weights)
  
  coeffs = het_model%load_model_coeffs(coords, weights) 

  ipoint = 1 
  do ix = 0, 1
    do iy = 0, 1
      do iz = 0, 1
        coeffs_ref = ix + iy*2 + iz*4
        call assert_comparable(coeffs(1,ipoint),   coeffs_ref*weights(ipoint), &
                               1d-10, 'vp is #index (random weight)')
        call assert_comparable(coeffs(2,ipoint),  -coeffs_ref*weights(ipoint), &
                               1d-10, 'vs is -#index (random weight)')
        call assert_comparable(coeffs(3,ipoint), 2*coeffs_ref*weights(ipoint), &
                               1d-10, 'rho is 2*#index (random weight)')
        ipoint = ipoint + 1
      end do
    end do
  end do

end subroutine test_load_het_rtpv 
!-----------------------------------------------------------------------------------------

end module test_heterogeneities
!=========================================================================================
