program kerner

    use mpi
    use commpi,                      only: ppinit, master
    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, init_random_seed, &
                                           myrank, lu_out

    use inversion_mesh,              only: inversion_mesh_data_type
    use work_type_mod
    use type_parameter,              only: parameter_type
    use ftnunit,                     only: runtests_init, runtests, runtests_final
    use unit_tests,                  only: test_all

    implicit none
    type(inversion_mesh_data_type)      :: inv_mesh
    type(parameter_type)                :: parameters

    integer                             :: nelems, nvertices, nvertices_per_elem, ndimensions
    integer                             :: ierror, ntasks
    character(len=64)                   :: fmtstring
    

    verbose = 0

    call init_random_seed()

    call runtests_init
    call runtests( test_all )
    call runtests_final

    verbose = 1


    write(*,*) '***************************************************************'
    write(*,*) ' Initialize MPI communication'
    write(*,*) '***************************************************************'
    call ppinit

    
    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Read input files for parameters, source and receivers'
    write(lu_out,*) '***************************************************************'
    call parameters%read_parameters('inparam_basic')
    call parameters%read_source()
    call parameters%read_receiver()
    

    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Read inversion mesh'
    write(lu_out,*) '***************************************************************'
    !call inv_mesh%read_tet_mesh('vertices.USA10', 'facets.USA10')
    call inv_mesh%read_abaqus_mesh(parameters%mesh_file)

    nvertices = inv_mesh%get_nvertices()
    nvertices_per_elem = inv_mesh%nvertices_per_elem
    nelems    = inv_mesh%get_nelements()
    fmtstring = '(A, I8, A, I8)'
    write(lu_out,fmtstring) '  nvertices: ',  nvertices, ', nelems: ', nelems
    

    select case(trim(parameters%whattodo))
    case('integratekernel')

        write(lu_out,*) '***************************************************************'
        write(lu_out,*) ' Initialize MPI work type'
        write(lu_out,*) '***************************************************************'
        ntasks = parameters%nelems_per_task
        call init_work_type(nkern = parameters%nkernel, &
                            nelems_per_task = parameters%nelems_per_task, &
                            nvertices_per_elem = inv_mesh%nvertices_per_elem)
      

        write(lu_out,*) '***************************************************************'
        write(lu_out,*) ' Master and slave part ways'
        write(lu_out,*) '***************************************************************'
        if (master) then
           print *, 'MASTER'
           call do_master(parameters, inv_mesh)
        else
           print *, 'SLAVE', myrank
           call do_slave(parameters, inv_mesh)
        endif
      
        call MPI_FINALIZE(ierror)

    case('plot_wavefield')
        stop 'Closed due to roadworks'
    end select

    write(lu_out,*)
    write(lu_out,*) ' Finished!'

    
contains

!=========================================================================================
subroutine do_slave(parameters, inv_mesh)
    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, &
                                           init_random_seed, DIETAG
    use inversion_mesh,              only: inversion_mesh_data_type
    use readfields,                  only: semdata_type
    use type_parameter,              only: parameter_type
    use fft,                         only: rfft_type, taperandzeropad
    use mpi

    implicit none
    type(inversion_mesh_data_type), intent(in)      :: inv_mesh
    type(parameter_type), intent(in)                :: parameters
    type(semdata_type)                  :: sem_data
    type(rfft_type)                     :: fft_data

    integer                             :: npoints, ndumps, nelems, ntimes, nomega
    integer                             :: ikernel, nvertices
    integer                             :: mpistatus(MPI_STATUS_SIZE)
    real(kind=dp)                       :: df
    character(len=32)                   :: filtername
    character(len=64)                   :: fmtstring
    integer                             :: nptperstep
    real(kind=dp), allocatable          :: workresult(:,:,:)

    nptperstep = parameters%npoints_per_step

    nvertices = inv_mesh%nvertices_per_elem
    ! nelems    = inv_mesh%get_nelements()

    !write(*,*) '***************************************************************'
    !write(*,*) ' Initialize and open AxiSEM wavefield files'
    !write(*,*) '***************************************************************'
    call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir)
    call sem_data%open_files()
    call sem_data%read_meshes()
    call sem_data%build_kdtree()

    call sem_data%load_seismogram(parameters%receiver, parameters%source)

    ndumps = sem_data%ndumps

    !print *, 'What to do: ', trim(parameters%whattodo)

    !!write(*,*) '***************************************************************'
    !write(*,*) ' Initialize FFT'
    !write(*,*) '***************************************************************'
    call fft_data%init(ndumps, nptperstep, sem_data%dt)
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    df     = fft_data%get_df()
    fmtstring = '(A, I8, A, I8)'
    write(lu_out,fmtstring) '  ntimes: ',  ntimes,     '  , nfreq: ', nomega
    fmtstring = '(A, F8.3, A, F8.3, A)'
    write(lu_out,fmtstring) '  dt:     ', sem_data%dt, ' s, df:    ', df*1000, ' mHz'

    !write(*,*) '***************************************************************'
    !write(*,*) ' Define filters'
    !write(*,*) '***************************************************************'
    call parameters%read_filter(nomega, df)

    !write(*,*) '***************************************************************'
    !write(*,*) ' Define kernels'
    !write(*,*) '***************************************************************'
    call parameters%read_kernel(sem_data, parameters%filter)

    do 
       ! Receive a message from the master
       call MPI_Recv(wt, 1, wt%mpitype, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, ierror)

       ! Check the tag of the received message. If no more work to do, exit loop
       ! and return to main programm
       if (mpistatus(MPI_TAG) == DIETAG) exit

       !print *, wt%ielements
       write(lu_out,*) '***************************************************************'

       allocate(workresult(parameters%nkernel, nvertices, wt%nelems_per_task))
       workresult = slave_work(parameters, sem_data, inv_mesh, fft_data, wt%ielements)
       wt%kernel_values = workresult
       deallocate(workresult)

       ! Send the result back
       call MPI_Send(wt, 1, wt%mpitype, 0, 0, MPI_COMM_WORLD, ierror)
          
    end do

    write(lu_out,*)
    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Free memory of Kernelspecs'
    write(lu_out,*) '***************************************************************'
    do ikernel = 1, parameters%nkernel
       call parameters%kernel(ikernel)%freeme() 
    end do
   
    !case('plot_wavefield')
    !    stop 'Closed due to construction work'
    !   write(*,*) '***************************************************************'
    !   write(*,*) 'Initialize output file'
    !   write(*,*) '***************************************************************'
    !   call inv_mesh%init_data(parameters%nkernel)

    !   write(*,*) '***************************************************************'
    !   write(*,*) ' Initialize FFT'
    !   write(*,*) '***************************************************************'
    !   call fft_data%init(ndumps, nvertices, sem_data%dt)
    !   ntimes = fft_data%get_ntimes()
    !   nomega = fft_data%get_nomega()
    !   df     = fft_data%get_df()
    !   fmtstring = '(A, I8, A, I8)'
    !   print fmtstring, '  ntimes: ',  ntimes,     '  , nfreq: ', nomega
    !   fmtstring = '(A, F8.3, A, F8.3, A)'
    !   print fmtstring, '  dt:     ', sem_data%dt, ' s, df:    ', df*1000, ' mHz'

    !   print *, 'Initialize XDMF file'
    !   allocate(co_points(3, nvertices))
    !   co_points = inv_mesh%get_vertices()
    !   call inv_mesh%init_data(ndumps*3)

    !   write(*,*) ' Read in forward field'
    !   allocate(fw_field(ndumps, nvertices))
    !   fw_field = sem_data%load_fw_points(dble(co_points), parameters%source)

    !   
    !   ! Dump forward field to XDMF file
    !   do idump = 1, ndumps
    !       if (mod(idump, 100)==0) write(*,*) ' Passing dump ', idump, ' to inversion mesh datatype'
    !       !Test of planar wave , works
    !       !fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
    !       !bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)
    !       call inv_mesh%set_data_snap(real(fw_field(idump,:), kind=sp), idump, 'fwd_wavefield')
    !   end do
    !   write(*,*) ' FFT forward field'
    !   allocate(fw_field_fd(nomega, nvertices))
    !   call fft_data%rfft(taperandzeropad(fw_field, ntimes), fw_field_fd)
    !   deallocate(fw_field)
    !   call timeshift( fw_field_fd, fft_data%get_f(), sem_data%timeshift_fwd )

    !   write(*,*) ' Read in backward field'
    !   allocate(bw_field(ndumps, nvertices))
    !   bw_field = sem_data%load_bw_points(dble(co_points), parameters%receiver(1))
    !   ! Dump backward field to XDMF file
    !   do idump = 1, ndumps
    !       if (mod(idump, 100)==0) write(*,*) ' Passing dump ', idump, ' to inversion mesh datatype'
    !       !Test of planar wave , works
    !       !fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
    !       !bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)
    !       call inv_mesh%set_data_snap(real(bw_field(idump,:), kind=sp), idump+ndumps, 'bwd_wavefield')
    !   end do
    !   write(*,*) ' FFT backward field'
    !   allocate(bw_field_fd  (nomega, nvertices))
    !   call fft_data%rfft(taperandzeropad(bw_field, ntimes), bw_field_fd)
    !   deallocate(bw_field)
    !   call timeshift( bw_field_fd, fft_data%get_f(), sem_data%timeshift_bwd )

    !   write(*,*) ' Convolve wavefields'
    !   allocate(conv_field(ntimes, nvertices))
    !   allocate(conv_field_fd(nomega, nvertices))
    !   conv_field_fd = fw_field_fd * bw_field_fd
    !   call fft_data%irfft(conv_field_fd, conv_field)
    !   deallocate(conv_field_fd)
    !!   write(*,*) ' Apply filter'
    !!   fw_field_fd = gabor20%apply_2d(fw_field_fd)

    !!   write(*,*) ' iFFT product of fields'
    !!   allocate(conv_field   (ntimes, nptperstep))
    !!   call fft_data%irfft(fw_field_fd*bw_field_fd, conv_field)
    !!
    !!   !write(*,*) ' Write pickpoint to disc'
    !!   !write(41,*) fw_field(:,10000)
    !!   !write(42,*) bw_field(:,10000)
    !!   !write(43,*) conv_field(:,10000)
    !!   !read(*,*)


    !   do idump = 1, ndumps
    !      if (mod(idump, 100)==0) write(*,*) ' Passing dump ', idump, ' of convolved wavefield'
    !      call inv_mesh%set_data_snap(real(conv_field(idump,:)), idump+ndumps*2, 'field_convolved')
    !   end do 

    !!   write(*,*)
    !!   write(*,*) ' Writing data to disk'
    !   call inv_mesh%dump_mesh_data_xdmf('wavefield')
    !   call inv_mesh%freeme()

    !end select

    write(lu_out,*)
    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Free memory of inversion mesh datatype'
    write(lu_out,*) '***************************************************************'
    call inv_mesh%freeme

    write(lu_out,*)
    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Close AxiSEM wavefield files'
    write(lu_out,*) '***************************************************************'
    call sem_data%close_files()

    write(lu_out,*)
    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Free memory of FFT datatype'
    write(lu_out,*) '***************************************************************'
    call fft_data%freeme()

end subroutine do_slave
!=========================================================================================


!=========================================================================================
function slave_work(parameters, sem_data, inv_mesh, fft_data, elementlist) result(results)

    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, myrank

    use inversion_mesh,              only: inversion_mesh_data_type
    use readfields,                  only: semdata_type
    use type_parameter,              only: parameter_type
    use fft,                         only: rfft_type, taperandzeropad
    use filtering,                   only: timeshift
    use montecarlo,                  only: integrated_type, allallconverged, allisconverged

    type(inversion_mesh_data_type), intent(in)    :: inv_mesh
    type(parameter_type),           intent(in)    :: parameters
    type(semdata_type),             intent(in)    :: sem_data
    type(rfft_type),                intent(in)    :: fft_data
    integer,                        intent(in)    :: elementlist(:)

    type(integrated_type), allocatable  :: int_kernel(:)

    real(kind=dp),    allocatable       :: results(:,:,:)
    real(kind=dp),    allocatable       :: element_points(:,:,:)
    real(kind=dp),    allocatable       :: co_points(:,:), K_x(:,:), veloseis(:)
    real(kind=dp),    allocatable       :: fw_field(:,:)
    real(kind=dp),    allocatable       :: bw_field(:,:)
    complex(kind=dp), allocatable       :: fw_field_fd(:,:)
    complex(kind=dp), allocatable       :: bw_field_fd(:,:)
    complex(kind=dp), allocatable       :: conv_field_fd(:,:), conv_field_fd_filt(:,:)
    real(kind=dp),    allocatable       :: conv_field(:,:)
    real(kind=dp),    allocatable       :: random_points(:,:), kernelvalue(:,:), kernelvalue_vertex(:,:)
    integer,          allocatable       :: niterations(:,:)
    !real(kind=dp),    allocatable       :: volume(:)
    real(kind=dp)                       :: volume

    !real(kind=sp), allocatable          :: co_element(:,:) !(3,4)

    integer                             :: idump, ipoint, ielement, irec, iel, ivertex, ikernel
    integer                             :: nptperstep, ndumps, ntimes, nomega, nelements

    nptperstep = parameters%npoints_per_step
    ndumps = sem_data%ndumps
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    nvertices = inv_mesh%nvertices_per_elem
    nelements = size(elementlist)


    allocate(results(parameters%nkernel, nvertices, nelements))

    !write(*,*) '***************************************************************'
    !write (*,*) 'Initialize Kernel variables'
    !write(*,*) '***************************************************************'
    !allocate(niterations(parameters%nkernel, nelems))
    !niterations = 0
    !open(unit=lu_iterations, file='niterations.txt', action='write')

    allocate(kernelvalue(nptperstep, parameters%nkernel))
    allocate(kernelvalue_vertex(nptperstep, parameters%nkernel))

    !allocate(co_element(3, nvertices))
    
    allocate(fw_field(ndumps, nptperstep))
    allocate(bw_field(ndumps, nptperstep))

    allocate(conv_field   (ntimes, nptperstep))
    allocate(fw_field_fd  (nomega, nptperstep))
    allocate(bw_field_fd  (nomega, nptperstep))
    allocate(conv_field_fd(nomega, nptperstep))
    allocate(conv_field_fd_filt(nomega, nptperstep))

    !allocate(volume(nelements))
    volume = 0.0

    allocate(niterations(parameters%nkernel, nelements))
    niterations = 0

    !write(*,*) '***************************************************************'
    !write(*,*) ' Start loop over elements'
    !write(*,*) '***************************************************************'

    allocate(random_points(3, nptperstep))
    allocate(int_kernel(inv_mesh%nvertices_per_elem))


    do iel = 1, nelements 
        ielement = elementlist(iel)
        if (ielement==-1) cycle

        !co_element = inv_mesh%get_element(ielement)
        volume = inv_mesh%get_volume(ielement)

        
        ! Initialize Monte Carlo integral for this element
        do ivertex = 1, inv_mesh%nvertices_per_elem
            call int_kernel(ivertex)%initialize_montecarlo(parameters%nkernel, &
                                                           volume,   &
                                                           parameters%allowed_error) 
        end do
      
        do while (.not.allallconverged(int_kernel)) ! Beginning of Monte Carlo loop
           random_points = inv_mesh%generate_random_points( ielement, nptperstep )
           
           ! Stop MC integration in this element after max_iter iterations
           if (any(niterations(:, iel)>parameters%max_iter)) exit 
           
           ! Load, FT and timeshift forward field
           fw_field = sem_data%load_fw_points( random_points, parameters%source )
           call fft_data%rfft( taperandzeropad(fw_field, ntimes), fw_field_fd )
           call timeshift( fw_field_fd, fft_data%get_f(), sem_data%timeshift_fwd )
         
           do irec = 1, parameters%nrec
                
              if (allallconverged(int_kernel,[(ikernel, ikernel =                                 & 
                                                        parameters%receiver(irec)%firstkernel,    &
                                                        parameters%receiver(irec)%lastkernel) ])) &
                                              cycle
          
              ! Load, FT and timeshift backward field
              bw_field = sem_data%load_bw_points( random_points, parameters%receiver(irec) )
              call fft_data%rfft( taperandzeropad(bw_field, ntimes), bw_field_fd )
              call timeshift( bw_field_fd, fft_data%get_f(), sem_data%timeshift_bwd )

              ! Convolve forward and backward fields
              conv_field_fd = fw_field_fd * bw_field_fd

              do ikernel = parameters%receiver(irec)%firstkernel, &
                           parameters%receiver(irec)%lastkernel 

                 ! If this kernel is already converged, go to the next one
                 if (allisconverged(int_kernel, ikernel)) cycle
                 niterations(ikernel, iel) = niterations(ikernel, iel) + 1

                 ! Apply Filter 
                 conv_field_fd_filt = parameters%kernel(ikernel)%apply_filter(conv_field_fd)

                 ! Backward FFT
                 call fft_data%irfft(conv_field_fd_filt, conv_field)

                 ! Calculate Scalar kernel from convolved time traces
                 kernelvalue(:,ikernel) = &
                     parameters%kernel(ikernel)%calc_misfit_kernel(conv_field)
              end do ! ikernel

           end do ! irec
           ! Check for convergence
           
           do ivertex = 1, inv_mesh%nvertices_per_elem
               do ikernel = 1, parameters%nkernel
                   kernelvalue_vertex(:, ikernel) = kernelvalue(:, ikernel) &
                                                  * inv_mesh%weights(ielement, ivertex, random_points)
               end do
               call int_kernel(ivertex)%check_montecarlo_integral(kernelvalue_vertex)
           end do

           ! Print convergence info
           write(fmtstring,"('(I6, I6, ES10.3, A, L1, ', I2, '(ES11.3), A, ', I2, '(ES10.3))')") &
                           parameters%nkernel, parameters%nkernel 
           write(lu_out,fmtstring) myrank, ielement, volume, ' Converged? ',                     &
                                   int_kernel(1)%areallconverged(), int_kernel(1)%getintegral(), &
                                   ' +- ', sqrt(int_kernel(1)%getvariance())

        end do ! Kernel coverged

        do ivertex = 1, inv_mesh%nvertices_per_elem
            !wt%kernel_values(:, ivertex, iel) = int_kernel(ivertex)%getintegral()
            results(:, ivertex, iel) = int_kernel(ivertex)%getintegral()
        end do

    end do ! iel

end function slave_work
!=========================================================================================

!=========================================================================================
subroutine do_master(parameters, inv_mesh)
    use commpi, only                             : nproc, master
    use global_parameters 
    type(inversion_mesh_data_type), intent(in)  :: inv_mesh
    type(parameter_type), intent(in)            :: parameters
    integer                                     :: ikernel, ielement, iel, ivertex
    real(kind=dp),    allocatable               :: K_x(:,:)
    integer,          allocatable               :: connectivity(:,:)
    integer                                     :: nslaves, rank, ierror
    integer, allocatable                        :: tasks(:), output(:,:), sendrequest(:), work_done(:)
    integer                                     :: mpistatus(MPI_STATUS_SIZE)
    integer                                     :: itask, ntasks, ioutput
    integer, allocatable                        :: elems_in_task(:,:)
    character(len=64)                           :: fmtstring

    allocate(connectivity(inv_mesh%nvertices_per_elem, nelems))
    connectivity = inv_mesh%get_connectivity()
    allocate(K_x(nvertices, parameters%nkernel))
    K_x = 0.0


    ! Parallelization is implemented only for kernel calculation, not for wavefield plotting
    select case(trim(parameters%whattodo) )
    case('integratekernel')


       write(lu_out,*) '***************************************************************'
       write(lu_out,*) ' Define filters'
       write(lu_out,*) '***************************************************************'
       call parameters%read_filter()

       write(lu_out,*) '***************************************************************'
       write(lu_out,*) ' Define kernels'
       write(lu_out,*) '***************************************************************'
       call parameters%read_kernel()


       nelems    = inv_mesh%get_nelements()
       fmtstring = '(A, I8, A, I8)'
       ! Calculate number of tasks
       ntasks = inv_mesh%get_nelements() / parameters%nelems_per_task + 1
       print fmtstring, '  nelements: ',  inv_mesh%get_nelements(), ', ntasks: ', ntasks
        
       allocate(tasks(ntasks))
       allocate(elems_in_task(ntasks, parameters%nelems_per_task))
       do itask = 1, ntasks
          tasks(itask) = itask
          do iel = 1, parameters%nelems_per_task
              if (iel + (itask-1) * parameters%nelems_per_task <= nelems) then
                  elems_in_task(itask, iel) = iel + (itask-1) * parameters%nelems_per_task
              else
                  elems_in_task(itask, iel) = -1
              end if
          end do
       enddo
        
       ! Find out how many processes there are in the default communicator
       nslaves = nproc - 1 ! the master does not work
       allocate(sendrequest(nslaves))
       allocate(work_done(nslaves))
       work_done = 0

       if (nslaves > ntasks) then
          write(6,*) 'ERROR: more slaves then tasks'
          stop
       elseif (nslaves < 1) then
          write(6,*) 'ERROR: need at least 1 slave'
          stop
       endif



       itask = 0
       ! Seed the slaves; send one unit of work to each slave.
       do rank=1, nslaves

          ! Find the next item of work to do
          itask = itask + 1

          ! fill sendbuffer
          wt%ielements     = elems_in_task(itask,:)

          ! Send it to each rank (blocking)
          call MPI_Send(wt,                & ! message buffer
                        1,                 & ! one data item
                        wt%mpitype,        & ! data item is an integer
                        rank,              & ! destination process rank
                        WORKTAG,           & ! user chosen message tag
                        MPI_COMM_WORLD,    & ! default communicator
                        sendrequest(rank), &
                        ierror)
       enddo

       ! Loop over tasks until there is no more work to be done
       ioutput = 0
       do itask=nslaves+1, ntasks

          ! Receive results from any !sic! slave (blocking!)
          ioutput = ioutput + 1
          call MPI_Recv(wt,               & ! message buffer
                        1,                & ! one data item
                        wt%mpitype,       & ! data item is an integer
                        MPI_ANY_SOURCE,   & ! receive from any sender
                        MPI_ANY_TAG,      & ! any type of message
                        MPI_COMM_WORLD,   & ! default communicator
                        mpistatus,        & ! info about the received message
                        ierror)

          ! extract from receive buffer
          do iel = 1, parameters%nelems_per_task
              ielement = wt%ielements(iel)
              if (ielement.eq.-1) cycle
              do ivertex = 1, inv_mesh%nvertices_per_elem
                 K_x(connectivity(ivertex, ielement),:) = K_x(connectivity(ivertex, ielement),:) &
                                                          + wt%kernel_values(:, ivertex, iel)
              end do
          end do
          
          ! fill sendbuffer
          wt%ielements     = elems_in_task(itask,:)

          ! Send the same slave some more work to do (blocking)
          call MPI_Send(wt,               & ! message buffer
                        1,                & ! one data item
                        wt%mpitype,       & ! data item is an integer
                        mpistatus(MPI_SOURCE), & ! to who we just received from
                        WORKTAG,          & ! user chosen message tag
                        MPI_COMM_WORLD,   & ! default communicator
                        sendrequest(mpistatus(MPI_SOURCE)), &
                        ierror)

          work_done(mpistatus(MPI_SOURCE)) = work_done(mpistatus(MPI_SOURCE)) + 1          
          write(fmtstring,"('(',I3,'(I5, '' ''), F8.3''%'')')") nslaves + 1
          write(lu_out,fmtstring) work_done, sum(work_done), real(sum(work_done)) / real(ntasks) * 100.

       enddo
       write(lu_out,*) '***************************************************************'
       write(lu_out,*) ' All work distributed'
       write(lu_out,*) '***************************************************************'

       ! There's no more work to be distributed, so receive all the outstanding results 
       ! from the slaves (blocking, so when this loop is finished, work is done and
       ! results received in the buffer!).
       do rank=1, nslaves
          ioutput = ioutput + 1
          call MPI_Recv(wt,              & ! message buffer
                        1,               & ! one data item
                        wt%mpitype,      & ! data item is an integer
                        MPI_ANY_SOURCE,  & ! receive from any sender
                        MPI_ANY_TAG,     & ! any type of message
                        MPI_COMM_WORLD,  & ! default communicator
                        mpistatus,       & ! info about the received message
                        ierror)
          
          ! extract from receive buffer
          do iel = 1, parameters%nelems_per_task
              ielement = wt%ielements(iel)
              if (ielement.eq.-1) cycle
              do ivertex = 1, inv_mesh%nvertices_per_elem
                 K_x(connectivity(ivertex, ielement),:) = K_x(connectivity(ivertex, ielement),:) &
                                                          + wt%kernel_values(:, ivertex, iel)
              end do
          end do

          work_done(mpistatus(MPI_SOURCE)) = work_done(mpistatus(MPI_SOURCE)) + 1          
          write(fmtstring,"('(',I3,'(I5, '' ''), F8.3''%'')')") nslaves + 1
          write(lu_out,fmtstring) work_done, sum(work_done), real(sum(work_done)) / real(ntasks) * 100.
          !write(fmtstring,"('(',I3,'(I5, '' ''))')") nslaves + 1
          !print fmtstring, work_done, sum(work_done)

       enddo

!       do itask=1, ntasks
!          write(6,*) output(itask,2), output(itask,1)
!       enddo
!
       write(lu_out,*) '***************************************************************'
       write(lu_out,*) ' All work collected, shutting down the slaves'
       write(lu_out,*) '***************************************************************'

       ! Tell all the slaves to exit by sending an empty message with the DIETAG.
       do rank=1, nslaves
          call MPI_Send(0,               & !
                        0,               & ! empty message
                        MPI_INTEGER,     & !
                        rank,            & ! destination
                        DIETAG,          & ! the tag conatains the actual information
                        MPI_COMM_WORLD,  & ! default communicator
                        sendrequest(rank), &
                        ierror)
       enddo

       write(lu_out,*) '***************************************************************'
       write(lu_out,*) 'Initialize output file'
       write(lu_out,*) '***************************************************************'
       call inv_mesh%init_data(parameters%nkernel)



       ! Save big kernel variable to disk
       !if ((mod(ielement, 1000)==0).or.(ielement==nelems)) then
       write(lu_out,*) 'Write Kernel to disk'
       do ikernel = 1, parameters%nkernel
          call inv_mesh%set_data_snap(real(K_x(:,ikernel), kind=sp), &
                                      ikernel, parameters%kernel(ikernel)%name )
       end do
       !call inv_mesh%set_data_snap(real(volume(:), kind=sp), &
       !                               parameters%nkernel + 1, 'volume' )

       call inv_mesh%dump_mesh_data_xdmf('gaborkernel')
       !end if

       !write(lu_iterations,*) ielement, niterations(:, ielement)
    end select


end subroutine do_master
!=========================================================================================



end program
