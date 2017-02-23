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
module plot_wavefields

  ! Collect routines for wavefield plotting
  use inversion_mesh,              only: inversion_mesh_data_type
  use type_parameter,              only: parameter_type
  use global_parameters,           only: sp, dp, lu_out, long, verbose
  implicit none
  private

  public init_wavefield_file, write_wavefield_data

contains

!-----------------------------------------------------------------------------------------
subroutine init_wavefield_file(inv_mesh, parameters, dt, ndumps)
  type(inversion_mesh_data_type) :: inv_mesh
  type(parameter_type)           :: parameters
  real(kind=dp), intent(in)      :: dt         
  integer,       intent(in)      :: ndumps
  integer, parameter             :: ndim = 6
  integer                        :: ikernel, idim
  character(len=2), parameter   :: dim_name(ndim) = ['tt', 'pp', 'rr', 'pr', 'tr', 'tp']
  character(len=32)             :: fmtstring
  character(len=512)            :: var_name

  call inv_mesh%init_cell_data(starttime = 0.d0, dt = dt)

  do ikernel = 1, parameters%nkernel

    fmtstring = '("fw_", A, "_straintrace")'
    write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name)
    call inv_mesh%add_cell_variable(var_name, nentries=ndumps, istime=.true.)

    fmtstring = '("bw_", A, "_straintrace")'
    write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name)
    call inv_mesh%add_cell_variable(var_name, nentries=ndumps, istime=.true.)

    fmtstring = '("conv_", A, "_")'
    write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name)
    call inv_mesh%add_cell_variable(var_name, nentries=ndumps, istime=.true.)

  end do

end subroutine init_wavefield_file
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_wavefield_data(inv_mesh, parameters, idx_elems, &
                                fw_field, bw_field, conv_field)
  type(inversion_mesh_data_type) :: inv_mesh
  type(parameter_type)           :: parameters
  integer, intent(in)            :: idx_elems(:)                             
  real(kind=sp), intent(in)      :: fw_field(:,:,:,:), bw_field(:,:,:,:), conv_field(:,:,:)

  real(kind=sp)                  :: val_out(size(fw_field, 4), size(fw_field, 1))
  integer                        :: nelements, ielement
  integer, parameter             :: ndim = 6
  integer                        :: ikernel, idim
  character(len=2), parameter    :: dim_name(ndim) = ['tt', 'pp', 'rr', 'pr', 'tr', 'tp']
  character(len=32)              :: fmtstring
  character(len=512)             :: var_name

  do ikernel = 1, parameters%nkernel

    fmtstring = '("fw_", A, "_straintrace")'
    write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name)

    if (ndim==6) then ! full strain tensor, reduce to straintrace
      val_out(:, :) = transpose(sum(fw_field(:, [1, 2, 3], ikernel, :), dim=2))
    !  val_out(:, :) = transpose(fw_field(:, idim, ikernel, :))
    else
      val_out(:, :) = transpose(fw_field(:, 1, ikernel, :))
    end if
    call write_variable(inv_mesh, var_name, val_out, idx_elems)

    fmtstring = '("bw_", A, "_straintrace")'
    write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name)

    if (ndim==6) then ! full strain tensor, reduce to straintrace
      val_out(:, :) = transpose(sum(bw_field(:, [1, 2, 3], ikernel, :), dim=2))
    else
      val_out(:, :) = transpose(bw_field(:, 1, ikernel, :))
    end if
    call write_variable(inv_mesh, var_name, val_out, idx_elems)


    fmtstring = '("conv_", A, "_")'
    write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name)
    val_out(:, :) = transpose(conv_field(:, ikernel, :))
    call write_variable(inv_mesh, var_name, val_out, idx_elems)

  end do

end subroutine write_wavefield_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_variable(inv_mesh, var_name, val_out, idx_elems)
  ! Call add_cell_data for variable var_name. 
  ! Mainly a wrapper to include the timing routines and check whether only a 
  ! part of val_out should be written to disk.
  use simple_routines, only       : first_occurence
  type(inversion_mesh_data_type) :: inv_mesh
  character(len=512), intent(in) :: var_name
  real(kind=sp),      intent(in) :: val_out(:,:)
  integer,            intent(in) :: idx_elems(:)

  integer                        :: nelements
  integer(kind=long)             :: iclock, iclock_ref, ticks_per_sec, ilast
  real(kind=dp)                  :: waittime

  nelements = size(idx_elems, 1)

  call system_clock(count=iclock_ref, count_rate=ticks_per_sec)

  ! The last task is not full, find last element
  ! If first_occurence==-1, all elements are legit
  ilast = first_occurence(idx_elems, -1) - 1

  if (ilast<0) then ! All elements are legit
    call inv_mesh%add_cell_data(var_name = var_name,   &
                                values   = val_out,    &
                                ielement = [idx_elems(1), idx_elems(nelements)])
  else
    call inv_mesh%add_cell_data(var_name = var_name,             &
                                values   = val_out(1:ilast, :),  &
                                ielement = [idx_elems(1), idx_elems(ilast)])
  end if

  call system_clock(count=iclock)
  waittime = real(iclock-iclock_ref, kind=dp) / real(ticks_per_sec, kind=dp)

  if (verbose>1) then
    write(lu_out,*) 'Variable ', trim(var_name), ': took ', waittime, ' seconds'
    call flush(lu_out)
  end if
end subroutine write_variable
!-----------------------------------------------------------------------------------------

end module plot_wavefields
!=========================================================================================
