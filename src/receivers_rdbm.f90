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
module receivers_rdbm

  use global_parameters
  use source_class,         only : src_param_type
  use receiver_class,      only : rec_param_type
  use commpi,               only : pabort

  implicit none

  private
  public   :: receivers_rdbm_type

  type receivers_rdbm_type
     type(src_param_type), allocatable      :: reci_sources(:)
     integer                                :: num_rec
     logical                                :: initialized = .false.
     contains
     procedure, pass                        :: create_reci_sources
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine create_reci_sources(this, rec_in)
  class(receivers_rdbm_type)       :: this
  type(rec_param_type), intent(in) :: rec_in(:)

  integer                          :: i,nrec

  nrec = size(rec_in)

  if (this%initialized) then
     write(6,*) 'ERROR: trying to initialize aleady initalized rdbm receiver object'
     call pabort
  endif

  allocate(this%reci_sources(nrec))

  do i=1,nrec
     call this%reci_sources(i)%init(rec_in(i)%latd, rec_in(i)%lond, (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0 /), 0d0)
  enddo


  this%initialized = .true.

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
