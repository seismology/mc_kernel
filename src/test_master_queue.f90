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
module test_master_queue

  use global_parameters
  use master_queue
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
! These routines are difficult to test one by one, so one big test for all
subroutine test_master_all
  integer       ::  ntasks

  testing = .false.
  call init_queue(ntasks, './inparam_test')

  call finalize

end subroutine test_master_all
!-----------------------------------------------------------------------------------------

end module test_master_queue

