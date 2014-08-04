program kerner_code

    use mpi
    use commpi,                      only: ppinit, pbroadcast_int, ppend
    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, init_random_seed, &
                                           master, lu_out, myrank, id_fft, id_fwd, id_bwd

    use inversion_mesh,              only: inversion_mesh_data_type
    use type_parameter,              only: parameter_type
    use ftnunit,                     only: runtests_init, runtests, runtests_final
    use unit_tests,                  only: test_all
    use slave_mod,                   only: do_slave
    use master_module,               only: do_master
    use work_type_mod,               only: init_work_type
    use plot_wavefields_mod,         only: plot_wavefields

    implicit none
    type(inversion_mesh_data_type)      :: inv_mesh
    type(parameter_type)                :: parameters

    integer                             :: nvertices_per_elem
    integer                             :: nvertices_per_task
    integer                             :: nbasisfuncs_per_elem
    integer                             :: nbasisfuncs_per_task
    integer                             :: ierror
    

    verbose = 0

    call init_random_seed()

    call runtests_init
    call runtests( test_all )
    call runtests_final

    verbose = 1


    call ppinit
    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' MPI communication initialized, I have rank', myrank
    write(lu_out,*) '***************************************************************'
   
    if (.not.master) call start_clock()

    call parameters%read_parameters()
    call parameters%read_source()
    call parameters%read_receiver()
    


    select case(trim(parameters%whattodo))
    case('integratekernel')

       if (master) then
           ! Get type of mesh and number of vertices per element
           if (trim(parameters%mesh_file).eq.'Karin') then
               nvertices_per_elem = 4
               nbasisfuncs_per_elem = 4
           else
               call inv_mesh%read_abaqus_meshtype(parameters%mesh_file,parameters%inttype)
               nbasisfuncs_per_elem = inv_mesh%nbasisfuncs_per_elem
               nvertices_per_elem = inv_mesh%nvertices_per_elem
               call inv_mesh%freeme()
           end if
       end if

       call pbroadcast_int(nbasisfuncs_per_elem, 0)
       call pbroadcast_int(nvertices_per_elem, 0)

       nvertices_per_task = parameters%nelems_per_task * nvertices_per_elem
       nbasisfuncs_per_task = parameters%nelems_per_task * nbasisfuncs_per_elem

       call pbroadcast_int(nbasisfuncs_per_task, 0)
       call pbroadcast_int(nvertices_per_task, 0)

       write(lu_out,*) '***************************************************************'
       write(lu_out,*) ' Initialize MPI work type'
       write(lu_out,*) '***************************************************************'


        call init_work_type(nkernel            = parameters%nkernel,         &
                            nelems_per_task    = parameters%nelems_per_task, &
                            nvertices          = nvertices_per_task,   &
                            nvertices_per_elem = nvertices_per_elem,          &
                            nbasisfuncs_per_elem = nbasisfuncs_per_elem)
      

        write(lu_out,*) '***************************************************************'
        write(lu_out,*) ' Master and slave part ways'
        write(lu_out,*) '***************************************************************'
        if (master) then
           call do_master()
        else
           call do_slave()
        endif
 
        call ppend()
        if (.not.master) call end_clock()   
  
    case('plot_wavefield')
        if (master) then
            call start_clock()
            call plot_wavefields()
            call end_clock()
        else
            print *, 'Nothing to do on rank ', myrank
        end if
    end select


    write(lu_out,*)
    write(lu_out,*) ' Finished!'
contains


!-----------------------------------------------------------------------------
subroutine start_clock
  !
  ! Driver routine to start the timing, using the clocks_mod module.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  use global_parameters
  use clocks_mod, only : clock_id, clocks_init
  
  implicit none
  
  character(len=8)  :: mydate
  character(len=10) :: mytime

  call date_and_time(mydate,mytime) 
  write(lu_out,11) mydate(5:6), mydate(7:8), mydate(1:4), mytime(1:2), mytime(3:4)

11 format('     Kerner started on ', A2,'/',A2,'/',A4,' at ', A2,'h ',A2,'min',/)


  call clocks_init(0)

  !id_read = clock_id('read_parameters')
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
  id_kernel          = clock_id('Kernel routines')
  id_init            = clock_id('Initialization per task')
  id_mpi             = clock_id('MPI communication with Master')


end subroutine start_clock
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
