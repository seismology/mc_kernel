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
module test_worktype

  use global_parameters
  use work_type_mod
  use ftnunit
  use commpi
# ifndef include_mpi
  use mpi
# endif
  implicit none

# ifdef include_mpi
  include 'mpif.h'
# endif
  public
contains


!------------------------------------------------------------------------------
subroutine test_init_work_type

  integer           :: nkernel = 2
  integer           :: nelems_per_task = 3
  integer           :: nvertices_per_task = 10
  integer           :: nvertices_per_elem = 4
  integer           :: nbasisfuncs_per_elem = 1
  integer           :: nmodel_parameters = 6
  integer           :: nmodel_parameters_hetero = 6

  logical           :: plot_wavefields = .true.
  integer           :: ndumps = 100
  integer           :: ndim = 6
  real              :: dt = 1e-3

  call init_work_type(nkernel                  = nkernel                  , &
                      nelems_per_task          = nelems_per_task          , &
                      nvertices                = nvertices_per_task       , &
                      nvertices_per_elem       = nvertices_per_elem       , &
                      nbasisfuncs_per_elem     = nbasisfuncs_per_elem     , &
                      nmodel_parameters        = nmodel_parameters        , &
                      nmodel_parameters_hetero = nmodel_parameters_hetero , &
                      plot_wavefields          = plot_wavefields          , &
                      ndumps                   = ndumps                   , &
                      ndim                     = ndim                     , &
                      dt                       = dt                       )
  


  ! Check values of parameters in resulting worktype
  call assert_equal(wt%ntotal_kernel           , nkernel, 'nkernel')
  call assert_equal(wt%nelems_per_task         , nelems_per_task, 'nelems_per_task')
  call assert_equal(wt%nvertices               , nvertices_per_task, 'nvertices')
  call assert_equal(wt%nvertices_per_elem      , nvertices_per_elem, 'nvertices_per_elem')
  call assert_equal(wt%nbasisfuncs_per_elem    , nbasisfuncs_per_elem, 'nbasisfuncs_per_elem')
  call assert_equal(wt%nmodel_parameters       , nmodel_parameters, 'nmodel_parameters')
  call assert_equal(wt%nmodel_parameters_hetero, nmodel_parameters_hetero, 'nmodel_parameters_hetero')

  call assert_equal(wt%ndumps                  , ndumps, 'ndumps')
  call assert_equal(wt%ndim                    , ndim, 'ndim')
  call assert_equal(wt%dt                      , nint(dt*1.e6) ,'dt in microseconds')


  ! Check size of arrays in resulting worktype
  call assert_equal(size(wt%connectivity),     nvertices_per_elem*nelems_per_task, &
                    'size wt%connectivity')
  call assert_equal(size(wt%vertices),         3 * nvertices_per_task, &
                    'size wt%vertices')
  call assert_equal(size(wt%kernel_values),    nelems_per_task * nkernel * nbasisfuncs_per_elem, &
                    'size wt%kernel_values')
  call assert_equal(size(wt%kernel_variance),  nelems_per_task * nkernel * nbasisfuncs_per_elem, &
                    'size wt%kernel_variance')
  call assert_equal(size(wt%niterations),      nelems_per_task * nkernel, &
                    'size wt%niterations')
  call assert_equal(size(wt%computation_time), nelems_per_task, &
                    'size wt%computation_time')
  call assert_equal(size(wt%model),            nelems_per_task * nmodel_parameters * nbasisfuncs_per_elem, &
                    'size wt%model')
  call assert_equal(size(wt%hetero_model),     nelems_per_task * nmodel_parameters_hetero * nbasisfuncs_per_elem, &
                    'size wt%hetero_model')

  call assert_equal(size(wt%fw_field),         ndumps * ndim * nkernel * nbasisfuncs_per_elem * nelems_per_task, &
                    'size wt%fw_field')
  call assert_equal(size(wt%bw_field),         ndumps * ndim * nkernel * nbasisfuncs_per_elem * nelems_per_task, &
                    'size wt%bw_field')
  call assert_equal(size(wt%conv_field),       ndumps * 1 * nkernel * nbasisfuncs_per_elem * nelems_per_task, &
                    'size wt%conv_field')
                    

  call free_work_type()                 

  call ppend()
                    
end subroutine test_init_work_type
!------------------------------------------------------------------------------
end module test_worktype

