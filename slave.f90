module slave_mod

    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, init_random_seed, &
                                           myrank, lu_out
    use work_type_mod
    use mpi
    implicit none

    type slave_result_type  
      real(kind=dp), allocatable :: kernel_values(:,:,:)
      real(kind=dp), allocatable :: kernel_errors(:,:,:)
      integer,       allocatable :: niterations(:,:)
    end type

contains

!=========================================================================================
subroutine do_slave() !, inv_mesh)
    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, &
                                           init_random_seed, DIETAG
    use inversion_mesh,              only: inversion_mesh_data_type
    use readfields,                  only: semdata_type
    use type_parameter,              only: parameter_type
    use fft,                         only: rfft_type, taperandzeropad
    use mpi

    implicit none
    !type(parameter_type), intent(in)    :: parameters
    type(parameter_type)                :: parameters
    type(inversion_mesh_data_type)      :: inv_mesh
    type(semdata_type)                  :: sem_data
    type(rfft_type)                     :: fft_data
    type(slave_result_type)             :: slave_result

    integer                             :: npoints, ndumps, nelems, ntimes, nomega
    integer                             :: ikernel, ierror
    integer                             :: mpistatus(MPI_STATUS_SIZE)
    real(kind=dp)                       :: df
    real(kind=dp), allocatable          :: K_x(:,:)
    integer                             :: iel, ielement, ivertex
    character(len=32)                   :: filtername
    character(len=64)                   :: fmtstring, fnam
    integer                             :: nptperstep
    real(kind=dp), allocatable          :: workresult(:,:,:)

    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Read input files for parameters, source and receivers'
    write(lu_out,*) '***************************************************************'
    call parameters%read_parameters('inparam_basic')
    call parameters%read_source()
    call parameters%read_receiver()


    !write(lu_out,*) '***************************************************************'
    !write(lu_out,*) ' Read inversion mesh'
    !write(lu_out,*) '***************************************************************'
    !call inv_mesh%read_tet_mesh('vertices.USA10', 'facets.USA10')
    !call inv_mesh%read_abaqus_mesh(parameters%mesh_file)


    nptperstep = parameters%npoints_per_step
    !nvertices = inv_mesh%nvertices_per_elem
    ! nelems    = inv_mesh%get_nelements()

    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Initialize and open AxiSEM wavefield files'
    write(lu_out,*) '***************************************************************'
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

    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Define filters'
    write(lu_out,*) '***************************************************************'
    call parameters%read_filter(nomega, df)

    write(lu_out,*) '***************************************************************'
    write(lu_out,*) ' Define kernels'
    write(lu_out,*) '***************************************************************'
    call parameters%read_kernel(sem_data, parameters%filter)

    do 
       ! Receive a message from the master
       call MPI_Recv(wt, 1, wt%mpitype, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, ierror)

       call inv_mesh%initialize_mesh(wt%ielement_type, wt%vertices, wt%connectivity)
       ! Check the tag of the received message. If no more work to do, exit loop
       ! and return to main programm
       if (mpistatus(MPI_TAG) == DIETAG) exit

       write(lu_out,*) '***************************************************************'

       slave_result = slave_work(parameters, sem_data, inv_mesh, fft_data)

       wt%kernel_values = slave_result%kernel_values
       wt%kernel_errors = slave_result%kernel_errors
       wt%niterations   = slave_result%niterations


       !! Write out my part of the mesh
       !call inv_mesh%init_node_data(parameters%nkernel)
       !
       !if (.not.allocated(K_x)) allocate(K_x(inv_mesh%get_nvertices(), parameters%nkernel))
       !K_x = 0
       !do iel = 1, parameters%nelems_per_task
       !    ielement = iel  
       !    if (ielement.eq.-1) cycle
       !    do ivertex = 1, inv_mesh%nvertices_per_elem
       !       K_x(wt%connectivity(ivertex, ielement),:) = K_x(wt%connectivity(ivertex, ielement),:) &
       !                                                + wt%kernel_values(:, ivertex, iel)
       !    end do
       !end do

       !write(fnam, "('partial_kernels/kernel_part_data_', I5.5)") wt%itask
       !do ikernel = 1, parameters%nkernel
       !   call inv_mesh%set_node_data_snap(real(K_x(:,ikernel), kind=sp), &
       !                               ikernel, parameters%kernel(ikernel)%name )
       !end do

       !call inv_mesh%dump_node_data_xdmf(fnam)

       ! Finished writing my part of the mesh
       call inv_mesh%freeme

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
function slave_work(parameters, sem_data, inv_mesh, fft_data) result(slave_result)

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
    !integer,                        intent(in)    :: elementlist(:)
    type(slave_result_type)                       :: slave_result

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
    real(kind=dp)                       :: volume

    character(len=64)                   :: fmtstring
    integer                             :: idump, ipoint, ielement, irec, iel, ivertex, ikernel
    integer                             :: nptperstep, ndumps, ntimes, nomega, nelements
    integer                             :: nvertices_per_elem

    nptperstep = parameters%npoints_per_step
    ndumps = sem_data%ndumps
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    nvertices_per_elem = inv_mesh%nvertices_per_elem
    nelements = inv_mesh%get_nelements()


    allocate(slave_result%kernel_values(parameters%nkernel, nvertices_per_elem, nelements))
    allocate(slave_result%kernel_errors(parameters%nkernel, nvertices_per_elem, nelements))
    allocate(slave_result%niterations(parameters%nkernel, nelements))

    allocate(kernelvalue(nptperstep, parameters%nkernel))
    allocate(kernelvalue_vertex(nptperstep, parameters%nkernel))

    
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


    do ielement = 1, nelements 

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
           if (any(niterations(:, ielement)>parameters%max_iter)) exit 
           
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
                 niterations(ikernel, ielement) = niterations(ikernel, ielement) + 1

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
            slave_result%kernel_values(:, ivertex, ielement) = int_kernel(ivertex)%getintegral()
            slave_result%kernel_errors(:, ivertex, ielement) = sqrt(int_kernel(ivertex)%getvariance())
            slave_result%niterations(:,ielement)             = niterations(:,ielement)
            call int_kernel(ivertex)%freeme()
        end do

    end do ! iel

end function slave_work
!=========================================================================================

















end module
