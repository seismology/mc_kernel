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

module master_helper
  use global_parameters, only   : verbose, sp, lu_out
  implicit none

  private
  public create_tasks

  contains

!-----------------------------------------------------------------------------------------
!> Create tasks. Takes into account, which elements have already been computed in
!!               an earlier run
subroutine create_tasks(completed, nelems_per_task, ntasks, elems_in_task)
  logical, intent(in)               :: completed(:)       !< The size of completed sets the 
                                                          !! total number of elements (completed
                                                          !! or not)
  integer, intent(in)               :: nelems_per_task    !< How many elements per task?
  integer, intent(out)              :: ntasks             !< Returns number of incomplete tasks
  integer, allocatable, intent(out) :: elems_in_task(:,:) !< Maps elements to tasks

  integer                           :: itask, iel, iel_in_task, nelems, ncompleted
  character(len=30)                 :: fmtstring
  
  nelems = size(completed)

  fmtstring = '(A, I8, A, I8, A, I8)'
  ! Calculate number of tasks
  ntasks = ceiling(real(count(.not.completed), kind=sp) / nelems_per_task)
  ncompleted = count(completed)

  write(lu_out,fmtstring) ' nelements: ',  nelems, ', ntasks: ', ntasks, ', ncompleted: ', ncompleted
   
  allocate(elems_in_task(ntasks, nelems_per_task))

  iel = 0
  do itask = 1, ntasks
    !print *, 'itask: ', itask
    !call flush(6)
    iel_in_task = 1
    !do iel = 1, parameters%nelems_per_task
    do while (iel_in_task<=nelems_per_task)
      !iel = iel_in_task + (itask-1) * nelems_per_task
      iel = iel + 1
      if (iel <= nelems) then
          if (completed(iel)) cycle
          elems_in_task(itask, iel_in_task) = iel
      else
          elems_in_task(itask, iel_in_task) = -1
      end if
      !print *, 'iel_in_task: ', iel_in_task, ', iel: ', iel
      !call flush(6)
      iel_in_task = iel_in_task + 1
    end do
  enddo

end subroutine create_tasks  
!-----------------------------------------------------------------------------------------


end module master_helper
