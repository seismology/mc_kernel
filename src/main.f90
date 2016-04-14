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

program kerner_code

#ifndef include_mpi
    use mpi
#endif
    use commpi,                      only: ppinit, pbroadcast_int, ppend, pabort, pbarrier,&
                                           pbroadcast_log, ppsplit
    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, init_random_seed, &
                                           master, lu_out, myrank

    use simple_routines,             only: lowtrim

    use inversion_mesh,              only: inversion_mesh_data_type
    use readfields,                  only: semdata_type
    use type_parameter,              only: parameter_type
    use ftnunit,                     only: runtests_init, runtests, runtests_final
    use unit_tests,                  only: test_all
    use slave_mod,                   only: do_slave
    use master_module,               only: do_master
    use work_type_mod,               only: init_work_type
    use background_model,            only: nmodel_parameters
    use heterogeneities,             only: nmodel_parameters_hetero

    implicit none

#ifdef include_mpi
    include 'mpif.h'
#endif

    type(inversion_mesh_data_type)      :: inv_mesh
    type(parameter_type)                :: parameters
    type(semdata_type)                  :: sem_data

    integer                             :: nvertices_per_elem
    integer                             :: nvertices_per_task
    integer                             :: nbasisfuncs_per_elem
    integer                             :: nbasisfuncs_per_task
    logical                             :: plot_wavefields
    integer                             :: ndim
    integer                             :: ndumps
    real(kind=sp)                       :: dt

    verbose = 0

    call init_random_seed()

    call runtests_init
    call runtests( test_all )
    call runtests_final

    verbose = 1

    call ppinit()

    call ppsplit()
    !if (master) print *, 'Rank: ', myrank, ' is master!', MPI_COMM_NODE
    !if (firstslave) print *, 'Rank: ', myrank, ' is first slave!', MPI_COMM_NODE
    !if (.not.(master.or.ioworker.or.firstslave)) print *, 'Rank: ', myrank, ' is a slave!', MPI_COMM_NODE

    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' MPI communication initialized, I have rank', myrank
    write(lu_out,*) '***************************************************************'
   
    if (master) then
      call start_clock_master()
    else
      call start_clock_slave()
    end if

    call parameters%read_parameters()
    call parameters%read_source()
    call parameters%read_receiver()
    
    if (master) then
        ! Get type of mesh and number of vertices per element
        
        select case(lowtrim(parameters%mesh_file_type))
        case('tetrahedral') 
            nvertices_per_elem = 4
            nbasisfuncs_per_elem = 4
        case('abaqus') 
            call inv_mesh%read_abaqus_meshtype(parameters%mesh_file,parameters%int_type)
            nbasisfuncs_per_elem = inv_mesh%nbasisfuncs_per_elem
            nvertices_per_elem = inv_mesh%nvertices_per_elem
            call inv_mesh%freeme()
        case default
            write(*,*) 'Unknown mesh type: ', trim(parameters%mesh_file_type)
            call pabort(do_traceback=.false.)
        end select

        ! Get number of wavefield dumps (samples)
        call sem_data%set_params(parameters%fwd_dir,            &
                                 parameters%bwd_dir,            &
                                 parameters%strain_buffer_size, & 
                                 parameters%displ_buffer_size,  & 
                                 parameters%strain_type_fwd,    &
                                 parameters%source%depth,       &
                                 1,                             &
                                 .false.)

        call sem_data%open_files()

        ndumps = sem_data%ndumps
        ndim   = sem_data%get_ndim()
        dt     = real(sem_data%dt, kind=sp) !SP because work type has only 4Byte entries
        call sem_data%close_files()

    end if

    call pbroadcast_int(nbasisfuncs_per_elem, 0)
    call pbroadcast_int(nvertices_per_elem, 0)

    nvertices_per_task = parameters%nelems_per_task * nvertices_per_elem
    nbasisfuncs_per_task = parameters%nelems_per_task * nbasisfuncs_per_elem

    plot_wavefields = parameters%plot_wavefields

    call pbroadcast_int(nbasisfuncs_per_task, 0)
    call pbroadcast_int(nvertices_per_task, 0)
    call pbroadcast_int(ndumps, 0)
    call pbroadcast_int(ndim, 0)
    call pbroadcast_log(plot_wavefields, 0)

    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Initialize MPI work type'
    write(lu_out,*) '***************************************************************'


    call init_work_type(nkernel                  = parameters%nkernel,          &
                        nelems_per_task          = parameters%nelems_per_task,  &
                        nvertices                = nvertices_per_task,          &
                        nvertices_per_elem       = nvertices_per_elem,          &
                        nbasisfuncs_per_elem     = nbasisfuncs_per_elem,        &
                        nmodel_parameters        = nmodel_parameters,           &
                        nmodel_parameters_hetero = nmodel_parameters_hetero,    &
                        plot_wavefields          = plot_wavefields,             &
                        ndumps                   = ndumps,                      &
                        ndim                     = ndim,                        &
                        dt                       = dt                          )

    

    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Master and slave part ways'
    write(lu_out,*) '***************************************************************'
    if (master) then
       call do_master()
    else
       call do_slave()
    endif

    call end_clock()   

    ! Wait for all other threads to arrive here before finalizing
    call pbarrier()
    call ppend()
  
    write(lu_out,*)
    write(lu_out,*) ' Finished!'
contains


!-----------------------------------------------------------------------------
!> Driver routine to start the timing for the master, using the clocks_mod module.
subroutine start_clock_master
  use global_parameters
  use clocks_mod, only : clock_id, clocks_init
  
  implicit none
  
  character(len=8)  :: mydate
  character(len=10) :: mytime

  call date_and_time(mydate,mytime) 
  write(lu_out,11) mydate(5:6), mydate(7:8), mydate(1:4), mytime(1:2), mytime(3:4)

11 format('     Kerner started on ', A2,'/',A2,'/',A4,' at ', A2,'h ',A2,'min',/)


  call clocks_init(0)

  id_read_params     = clock_id('Read parameters')
  id_create_tasks    = clock_id('Create tasks for slaves')
  id_get_next_task   = clock_id('Get next task for Slave')
  id_extract         = clock_id('Extract receive buffer')
  id_mpi             = clock_id('MPI communication with Slaves')
  id_dump            = clock_id('Dump intermediate results')
  id_write_kernel    = clock_id('Write Kernels to disk')
  id_mult_kernel     = clock_id('Multiply Kernels with Model')


end subroutine start_clock_master
!=============================================================================

!-----------------------------------------------------------------------------
!> Driver routine to start the timing for the slave, using the clocks_mod module.
subroutine start_clock_slave
  use global_parameters
  use clocks_mod, only : clock_id, clocks_init
  
  implicit none
  
  character(len=8)  :: mydate
  character(len=10) :: mytime

  call date_and_time(mydate,mytime) 
  write(lu_out,11) mydate(5:6), mydate(7:8), mydate(1:4), mytime(1:2), mytime(3:4)

11 format('     Kerner started on ', A2,'/',A2,'/',A4,' at ', A2,'h ',A2,'min',/)


  call clocks_init(0)

  id_read_params     = clock_id('Read parameters')
  id_init_fft        = clock_id('Initialize FFT plan')
  id_fft             = clock_id('FFT routines')
  id_bwd             = clock_id('Reading bwd field')
  id_fwd             = clock_id('Reading fwd field')
  id_kdtree          = clock_id(' - KD-tree lookup (only fwd)')
  id_find_point_fwd  = clock_id(' - Find next GLL point (only fwd)')
  id_find_point_bwd  = clock_id(' - Find next GLL point (only bwd)')
  id_load_strain     = clock_id(' - Load_strain (fwd and bwd)')
  id_netcdf          = clock_id(' - - NetCDF routines')
  id_rotate          = clock_id(' - - Rotate fields')
  id_buffer          = clock_id(' - - Buffer routines')
  id_calc_strain     = clock_id(' - - Calculate strain')
  id_lagrange        = clock_id(' - - Lagrange interpolation')
  id_mc              = clock_id('Monte Carlo routines')
  id_filter_conv     = clock_id('Filtering and convolution')
  id_inv_mesh        = clock_id('Inversion mesh routines')
  id_int_model       = clock_id('Integrate model paramaters')
  id_int_hetero      = clock_id('Integrate heterogeneity model')
  id_kernel          = clock_id('Kernel routines')
  id_init            = clock_id('Initialization per task')
  id_mpi             = clock_id('MPI communication with Master')
  id_out             = clock_id('Write wavefields to disk')
  id_finalize        = clock_id('Finalization of output files')
  id_element         = clock_id('Time spent for one element')
  id_allocate        = clock_id('Allocation in load_strain')
  id_allocate_2      = clock_id('Allocation in load_strain 2')


end subroutine start_clock_slave
!=============================================================================

!-----------------------------------------------------------------------------
subroutine end_clock 
  !
  ! Wapper routine to end timing and display clock informations.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  use clocks_mod, only : clocks_exit

  implicit none

  write(lu_out,*)
  write(lu_out,"(10x,'Summary of timing measurements:')")
  write(lu_out,*)

  call clocks_exit(0)

  write(lu_out,*)

end subroutine end_clock
!=============================================================================
end program
