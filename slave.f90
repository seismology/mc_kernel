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
subroutine do_slave() 
    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, &
                                           init_random_seed, DIETAG, id_mpi
    use inversion_mesh,              only: inversion_mesh_data_type
    use readfields,                  only: semdata_type
    use type_parameter,              only: parameter_type
    use fft,                         only: rfft_type, taperandzeropad
    use mpi
    use clocks_mod,                  only: tick

    implicit none
    type(parameter_type)                :: parameters
    type(inversion_mesh_data_type)      :: inv_mesh
    type(semdata_type)                  :: sem_data
    type(rfft_type)                     :: fft_data
    type(slave_result_type)             :: slave_result

    integer                             :: ndumps, ntimes, nomega
    integer                             :: ikernel, ierror, iclockold
    integer                             :: mpistatus(MPI_STATUS_SIZE)
    real(kind=dp)                       :: df
    real(kind=dp), allocatable          :: K_x(:,:)
    integer                             :: ielement, itask
    character(len=255)                  :: fmtstring, fnam
    integer                             :: nptperstep

    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Read input files for parameters, source and receivers'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)

    call parameters%read_parameters()
    call parameters%read_source()
    call parameters%read_receiver()

    nptperstep = parameters%npoints_per_step

    ! Master and slave part ways here for some time. 
    ! Master reads in the inversion mesh, slaves initialize the FFT

    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Initialize and open AxiSEM wavefield files'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)

    call sem_data%set_params(parameters%fwd_dir,     &
                             parameters%bwd_dir,     &
                             parameters%buffer_size, & 
                             parameters%model_param)

    call sem_data%open_files()
    call sem_data%read_meshes()
    call sem_data%build_kdtree()

    call sem_data%load_seismogram(parameters%receiver, parameters%source)

    ndumps = sem_data%ndumps

    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Initialize FFT'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)

    call fft_data%init(ndumps, sem_data%get_ndim(), nptperstep, sem_data%dt, &
                       fftw_plan=parameters%fftw_plan)
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    df     = fft_data%get_df()
    fmtstring = '(A, I8, A, I8)'
    write(lu_out,fmtstring) '  ntimes: ',  ntimes,     '  , nfreq: ', nomega
    fmtstring = '(A, F8.3, A, F8.3, A)'
    write(lu_out,fmtstring) '  dt:     ', sem_data%dt, ' s, df:    ', df*1000, ' mHz'

    ! Master and slave synchronize again

    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Define filters'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)

    call parameters%read_filter(nomega, df)

    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Define kernels'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)

    call parameters%read_kernel(sem_data, parameters%filter)
    
    itask = 0
    do 
       write(lu_out,'(A)') '***************************************************************'
       write(lu_out,'(A)') ' Waiting for tasks from master rank'
       write(lu_out,'(A)') '***************************************************************'
       call flush(lu_out)


       ! Receive a message from the master
       iclockold = tick()
       call MPI_Recv(wt, 1, wt%mpitype, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, ierror)
       iclockold = tick(id=id_mpi, since=iclockold)

       ! Check the tag of the received message. If no more work to do, exit loop
       ! and return to main program
       if (mpistatus(MPI_TAG) == DIETAG) exit

       itask = itask + 1
       write(lu_out,'(A)') '***************************************************************'
       write(lu_out,'(A, I3, A)') ' Received task # ', itask, ', going to work'
       write(lu_out,'(A)') '***************************************************************'
       call flush(lu_out)


       call inv_mesh%initialize_mesh(wt%ielement_type, wt%vertices, &
                                     wt%connectivity, wt%nbasisfuncs_per_elem)


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
       iclockold = tick()
       call MPI_Send(wt, 1, wt%mpitype, 0, 0, MPI_COMM_WORLD, ierror)
       iclockold = tick(id=id_mpi, since=iclockold)
          
    end do

    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A, I3, A)') ' All work done. Did ', itask, ' tasks in total'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)


    write(lu_out,'(A)')
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Free memory of Kernelspecs'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)
    do ikernel = 1, parameters%nkernel
       call parameters%kernel(ikernel)%freeme() 
    end do
   
    write(lu_out,'(A)')
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Free memory of inversion mesh datatype'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)
    call inv_mesh%freeme

    write(lu_out,'(A)')
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Close AxiSEM wavefield files'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)
    call sem_data%close_files()

    write(lu_out,'(A)')
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Free memory of FFT datatype'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)
    call fft_data%freeme()

end subroutine do_slave
!=========================================================================================


!=========================================================================================
function slave_work(parameters, sem_data, inv_mesh, fft_data) result(slave_result)

    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, myrank, &
                                           id_fwd, id_bwd, id_fft, id_mc, id_filter_conv, &
                                           id_inv_mesh, id_kernel, id_init

    use inversion_mesh,              only: inversion_mesh_data_type
    use readfields,                  only: semdata_type
    use type_parameter,              only: parameter_type
    use fft,                         only: rfft_type, taperandzeropad
    use filtering,                   only: timeshift_type
    use montecarlo,                  only: integrated_type, allallconverged, allisconverged
    use clocks_mod,                  only: tick

    type(inversion_mesh_data_type), intent(in)    :: inv_mesh
    type(parameter_type),           intent(in)    :: parameters
    type(semdata_type),             intent(in)    :: sem_data
    type(rfft_type),                intent(in)    :: fft_data
    type(slave_result_type)                       :: slave_result
    type(timeshift_type)                          :: timeshift_fwd, timeshift_bwd

    type(integrated_type), allocatable  :: int_kernel(:)

    real(kind=dp),    allocatable       :: fw_field(:,:,:)
    real(kind=dp),    allocatable       :: bw_field(:,:,:)
    complex(kind=dp), allocatable       :: fw_field_fd(:,:,:)
    complex(kind=dp), allocatable       :: bw_field_fd(:,:,:)
    complex(kind=dp), allocatable       :: conv_field_fd(:,:), conv_field_fd_filt(:,:)
    real(kind=dp),    allocatable       :: conv_field(:,:)
    real(kind=dp),    allocatable       :: random_points(:,:), kernelvalue(:,:), kernelvalue_weighted(:,:)
    integer,          allocatable       :: niterations(:,:)
    real(kind=dp)                       :: volume

    character(len=256)                  :: fmtstring
    integer                             :: ielement, irec, ivertex, ikernel, ibasisfunc
    integer                             :: nptperstep, ndumps, ntimes, nomega, nelements
    integer                             :: nvertices_per_elem
    integer                             :: nbasisfuncs_per_elem
    integer                             :: iclockold
    integer                             :: ndim
    integer, parameter                  :: taper_length = 2 ! This is the bare minimum. It does 
                                                            ! not produce artifacts yet, at least 
                                                            ! no obvious ones
    
    iclockold = tick()

    ndim = sem_data%get_ndim()
    nptperstep = parameters%npoints_per_step
    ndumps = sem_data%ndumps
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    nbasisfuncs_per_elem = inv_mesh%nbasisfuncs_per_elem
    nelements = inv_mesh%get_nelements()


    allocate(slave_result%kernel_values(parameters%nkernel, nbasisfuncs_per_elem, nelements))
    allocate(slave_result%kernel_errors(parameters%nkernel, nbasisfuncs_per_elem, nelements))
    allocate(slave_result%niterations(parameters%nkernel, nelements))

    allocate(kernelvalue_weighted(nptperstep, parameters%nkernel))
    allocate(kernelvalue(nptperstep, parameters%nkernel))
    
    allocate(fw_field(ndumps, ndim, nptperstep))
    allocate(bw_field(ndumps, ndim, nptperstep))

    allocate(conv_field        (ntimes, nptperstep))
    allocate(fw_field_fd       (nomega, ndim, nptperstep))
    allocate(bw_field_fd       (nomega, ndim, nptperstep))
    allocate(conv_field_fd     (nomega, nptperstep))
    allocate(conv_field_fd_filt(nomega, nptperstep))

    volume = 0.0

    allocate(niterations(parameters%nkernel, nelements))
    niterations = 0

    if (.not.parameters%deconv_stf) then
      call timeshift_fwd%init(fft_data%get_f(), sem_data%timeshift_fwd)
      call timeshift_bwd%init(fft_data%get_f(), sem_data%timeshift_bwd)
    end if
    

    !write(*,*) '***************************************************************'
    !write(*,*) ' Start loop over elements'
    !write(*,*) '***************************************************************'

    allocate(random_points(3, nptperstep))
    allocate(int_kernel(nbasisfuncs_per_elem))

    iclockold = tick(id=id_init, since=iclockold)

    do ielement = 1, nelements 

        iclockold = tick()
        volume = inv_mesh%get_volume(ielement)

        ! Initialize Monte Carlo integral for this element
        iclockold = tick(id=id_inv_mesh)
        do ibasisfunc = 1, nbasisfuncs_per_elem
            call int_kernel(ibasisfunc)%initialize_montecarlo(parameters%nkernel, &
                                                           volume,   &
                                                           parameters%allowed_error
                                                           parameters%allowed_relative_error) 
        end do
      
        do while (.not.allallconverged(int_kernel)) ! Beginning of Monte Carlo loop
           iclockold = tick(id=id_mc)

           random_points = inv_mesh%generate_random_points( ielement, nptperstep )
           iclockold = tick(id=id_inv_mesh)
           
           ! Stop MC integration in this element after max_iter iterations
           if (any(niterations(:, ielement)>parameters%max_iter)) exit 
           
           ! Load forward field
           iclockold = tick()
           fw_field = sem_data%load_fw_points( random_points, parameters%source )
           iclockold = tick(id=id_fwd, since=iclockold)
           
           ! FFT of forward field
           call fft_data%rfft( taperandzeropad(fw_field, ntimes, taper_length), fw_field_fd )
           iclockold = tick(id=id_fft, since=iclockold)

           ! Timeshift forward field
           ! Not necessary, if STF is deconvolved
           if (.not.parameters%deconv_stf) then
              call timeshift_fwd%apply(fw_field_fd)
              iclockold = tick(id=id_filter_conv, since=iclockold)
           end if

         
           do irec = 1, parameters%nrec
                
              if (allallconverged(int_kernel,[(ikernel, ikernel =                                 & 
                                                        parameters%receiver(irec)%firstkernel,    &
                                                        parameters%receiver(irec)%lastkernel) ])) &
                                              cycle
          
              iclockold = tick(id=id_mc, since=iclockold)

              ! Load backward field
              bw_field = sem_data%load_bw_points( random_points, parameters%receiver(irec) )
              iclockold = tick(id=id_bwd, since=iclockold)

              ! FFT of backward field
              call fft_data%rfft( taperandzeropad(bw_field, ntimes, taper_length), bw_field_fd )
              iclockold = tick(id=id_fft, since=iclockold)

              ! Timeshift backward field
              ! Not necessary, if STF is deconvolved
              if (.not.parameters%deconv_stf) then
                 call timeshift_bwd%apply(bw_field_fd)
              end if

              ! Convolve forward and backward fields
              ! The summation over dimension 2 is necessary for vs kernels
              conv_field_fd = sum(fw_field_fd * bw_field_fd, 2)
              iclockold = tick(id=id_filter_conv, since=iclockold)

              do ikernel = parameters%receiver(irec)%firstkernel, &
                           parameters%receiver(irec)%lastkernel 

                 ! If this kernel is already converged, go to the next one
                 if (allisconverged(int_kernel, ikernel)) cycle
                 niterations(ikernel, ielement) = niterations(ikernel, ielement) + 1
                 iclockold = tick(id=id_mc, since=iclockold)

                 ! Apply Filter 
                 conv_field_fd_filt = parameters%kernel(ikernel)%apply_filter(conv_field_fd)
                 iclockold = tick(id=id_filter_conv, since=iclockold)

                 ! iFFT of multiplied spectra
                 call fft_data%irfft(conv_field_fd_filt, conv_field)
                 iclockold = tick(id=id_fft, since=iclockold)

                 ! Calculate Scalar kernel from convolved time traces
                 kernelvalue(:,ikernel) = &
                     parameters%kernel(ikernel)%calc_misfit_kernel(conv_field)
                 iclockold = tick(id=id_kernel, since=iclockold)
              end do ! ikernel
              iclockold = tick(id=id_mc, since=iclockold)

           end do ! irec
           iclockold = tick(id=id_mc, since=iclockold)


           ! Check for convergence           
           do ibasisfunc = 1, nbasisfuncs_per_elem !inv_mesh%nvertices_per_elem
              iclockold = tick()
              do ikernel = 1, parameters%nkernel
                 kernelvalue_weighted(:, ikernel) = kernelvalue(:, ikernel) &
                      * inv_mesh%weights(ielement, ibasisfunc, random_points)
              end do
              iclockold = tick(id=id_kernel, since=iclockold)
              call int_kernel(ibasisfunc)%check_montecarlo_integral(kernelvalue_weighted)
              iclockold = tick(id=id_mc, since=iclockold)
           end do        
              
           ! Print convergence info
           write(fmtstring,"('(I6, I6, ES10.3, A, I0.5, A, I0.5, ', I4, '(ES11.3), A, ', I4, '(ES10.3))')") &
                           parameters%nkernel, parameters%nkernel 
           write(lu_out,fmtstring) myrank, ielement, volume, ' Converged? ',                            &
                                   int_kernel(1)%countconverged(), '/', parameters%nkernel,             &
                                   int_kernel(1)%getintegral(), ' +- ', sqrt(int_kernel(1)%getvariance())

        end do ! End of Monte Carlo loop
    
        if (any(niterations(:, ielement)>parameters%max_iter)) then
           fmtstring = "('Max number of iterations reached. ', I5 , ' kernels were not converged.')"
           write(lu_out, fmtstring) parameters%nkernel - int_kernel(1)%countconverged()
        else
           fmtstring = "('All kernels converged after', I5, ' iterations')"
           write(lu_out, fmtstring) maxval(niterations(:, ielement))
        end if


        do ibasisfunc = 1, nbasisfuncs_per_elem
           
           slave_result%kernel_values(:, ibasisfunc, ielement) = int_kernel(ibasisfunc)%getintegral()
           slave_result%kernel_errors(:, ibasisfunc, ielement) = sqrt(int_kernel(ibasisfunc)%getvariance())
           slave_result%niterations(:,ielement)             = niterations(:,ielement)

           call int_kernel(ibasisfunc)%freeme()
        end do

    end do ! ielement

    call timeshift_fwd%freeme()
    call timeshift_bwd%freeme()
end function slave_work
!=========================================================================================




end module
