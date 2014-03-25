program kerner

    use mpi
    use commpi,                      only: ppinit, pbroadcast_int
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

    integer                             :: nelems, nvertices, nvertices_per_elem, ndimensions
    integer                             :: nvertices_per_task
    integer                             :: ierror, ntasks, iclockold
    character(len=64)                   :: fmtstring
    

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
            else
                call inv_mesh%read_abaqus_meshtype(parameters%mesh_file)
                nvertices_per_elem = inv_mesh%nvertices_per_elem
            end if
        end if
        call pbroadcast_int(nvertices_per_elem, 0)
        nvertices_per_task = parameters%nelems_per_task * nvertices_per_elem
        call pbroadcast_int(nvertices_per_task, 0)

        write(lu_out,*) '***************************************************************'
        write(lu_out,*) ' Initialize MPI work type'
        write(lu_out,*) '***************************************************************'

        call init_work_type(nkernel            = parameters%nkernel,         &
                            nelems_per_task    = parameters%nelems_per_task, &
                            nvertices          = nvertices_per_task,   &
                            nvertices_per_elem = nvertices_per_elem)
      

        write(lu_out,*) '***************************************************************'
        write(lu_out,*) ' Master and slave part ways'
        write(lu_out,*) '***************************************************************'
        if (master) then
           !print *, 'MASTER'
           call do_master()
        else
           !print *, 'SLAVE', myrank
           call do_slave()
        endif
      
        call MPI_FINALIZE(ierror)

    case('plot_wavefield')
        if (master) then
            call plot_wavefields()
        else
            print *, 'Nothing to do on rank ', myrank
        end if
    end select

    if (.not.master) call end_clock()

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
  id_fft         = clock_id('FFT routines')
  id_fwd         = clock_id('Reading fwd field')
  id_netcdf      = clock_id(' - NetCDF routines')
  id_rotate      = clock_id(' - Rotate fields')
  id_buffer      = clock_id(' - Buffer routines')
  id_bwd         = clock_id('Reading bwd field')
  id_mc          = clock_id('Monte Carlo routines')
  id_filter_conv = clock_id('Filtering and convolution')
  id_inv_mesh    = clock_id('Inversion mesh routines')
  id_kernel      = clock_id('Kernel routines')
  id_init        = clock_id('Initialization per task')
  id_mpi         = clock_id('MPI communication with Master')
  !idold05 = clock_id('monte_carlo_routines')
  !idold06 = clock_id('bkgrdmodel_testing')
  !idold07 = clock_id('get_global no loop')
  !idold08 = clock_id('glob-slob/flob numbering')
  !idold09 = clock_id('get_global in loop')
  !idold11 = clock_id('create_pdb')
  !
  !idold12 = clock_id('define_glocal_numbering')
  !idold13 = clock_id('define_sflocal_numbering')
  !idold14 = clock_id('generate_serendipity_per_proc')


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
