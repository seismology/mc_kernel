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
  real(kind=dp)             :: data_resampled(ntimes)
  real(kind=dp)             :: data_in_coarse(250)
  integer                   :: x

  data_in(:) = 0
  data_in(ntimes/4+1:ntimes) = [(2.d0 * exp(-(x/thalf)**2)*x/thalf**2, x=1, ntimes-ntimes/4 )]

  !allocate(data_in_coarse(ceiling(real(ntimes)/ncoarsen)))
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
