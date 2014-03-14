program kerner

    use mpi
    use commpi,                      only: ppinit, pbroadcast_int
    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, init_random_seed, &
                                           master, lu_out, myrank

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
    integer                             :: ierror, ntasks
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
    
    call parameters%read_parameters('inparam_basic')
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

    write(lu_out,*)
    write(lu_out,*) ' Finished!'

end program
