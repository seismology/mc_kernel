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

module interpolate_mesh

  use global_parameters, only             : sp, dp
  implicit none
  private
  type parameter_interpolator
      real(kind=sp), allocatable         :: values(:)
      integer                            :: ndepth
      real(kind=dp)                      :: dr
      contains 
        private
        procedure, pass                  :: init
        procedure, pass, public          :: get
  end type

  public parameter_interpolator, create_interpolator

contains

!-----------------------------------------------------------------------------------------
function create_interpolator(param_tmp, tree, radius) result(interpolator)
  use kdtree2_module, only          :  kdtree2_result, kdtree2_n_nearest, kdtree2
  real(kind=sp), intent(in)         :: param_tmp(:)
  type(kdtree2), pointer            :: tree
  real(kind=dp), intent(in)         :: radius
  type(parameter_interpolator)      :: interpolator

  real(kind=sp), parameter          :: dz = 0.1d3 ! 100m
  real(kind=sp)                     :: coordinates_sz(2)
  real(kind=sp), allocatable        :: values(:)
  integer                           :: idepth, ndepth, pointid
  type(kdtree2_result), allocatable :: nextpoint(:)
  integer, parameter                :: nnext_points = 1

  ndepth = int(radius / dz) + 1
  allocate(values(ndepth))

  ! nextpoint has to be allocatable in kdtree module
  allocate(nextpoint(nnext_points))

  do idepth = 1, ndepth
    coordinates_sz(1) = 0
    coordinates_sz(2) = (idepth - 1) * dz
    call kdtree2_n_nearest( tree,                &
                            coordinates_sz,      &  
                            nn = nnext_points,   &
                            results = nextpoint )

    pointid = nextpoint(1)%idx
    values(idepth) = param_tmp(pointid)
  end do

  call interpolator%init(values, dz)  

end function create_interpolator
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init(this, values, dr)
  class(parameter_interpolator)       :: this 
  real(kind=sp), intent(in)           :: values(:)
  real(kind=sp)                       :: dr

  this%ndepth = size(values)
  allocate(this%values(this%ndepth))
  this%values = values

  this%dr = dr

end subroutine init
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function get(this, r) result(values)
  class(parameter_interpolator), intent(in)   :: this
  real(kind=dp), intent(in)                   :: r(:)
  real(kind=sp)                               :: values(size(r))

  integer                                     :: idx(size(r))

  idx = min(max(int(r / this%dr + 1), 1), this%ndepth)

  values = this%values(idx)

end function get
!-----------------------------------------------------------------------------------------

end module interpolate_mesh
