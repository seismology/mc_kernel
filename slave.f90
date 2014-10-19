!=========================================================================================
module slave_mod

  use global_parameters,           only: sp, dp, pi, deg2rad, &
                                         init_random_seed, myrank, lu_out
  use work_type_mod
# ifndef include_mpi
  use mpi
# endif
  implicit none

# ifdef include_mpi
  include 'mpif.h'
# endif

  type slave_result_type  
    real(kind=dp), allocatable :: kernel_values(:,:,:)
    real(kind=dp), allocatable :: kernel_variance(:,:,:)
    integer,       allocatable :: niterations(:,:)
    real(kind=dp), allocatable :: model(:,:,:)
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine do_slave() 
    use global_parameters,           only: DIETAG, id_mpi
    use inversion_mesh,              only: inversion_mesh_data_type
    use readfields,                  only: semdata_type
    use type_parameter,              only: parameter_type
    use fft,                         only: rfft_type, taperandzeropad
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
    integer                             :: itask
    character(len=255)                  :: fmtstring
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
                             parameters%strain_buffer_size, & 
                             parameters%displ_buffer_size, & 
                             parameters%strain_type_fwd)

    call sem_data%open_files()
    call sem_data%read_meshes()
    call sem_data%build_kdtree()

!    call sem_data%load_seismogram(parameters%receiver, parameters%source)
    call sem_data%load_seismogram_rdbm(parameters%receiver, parameters%source)

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

       ! Check the tag of the received message. If no more work to do, exit loop
       ! and return to main program
       if (mpistatus(MPI_TAG) == DIETAG) exit
       iclockold = tick(id=id_mpi, since=iclockold)

       itask = itask + 1
       write(lu_out,'(A)') '***************************************************************'
       write(lu_out,'(A, I3, A)') ' Received task # ', itask, ', going to work'
       write(lu_out,'(A)') '***************************************************************'
       call flush(lu_out)


       call inv_mesh%initialize_mesh(wt%ielement_type, wt%vertices, &
                                     wt%connectivity, wt%nbasisfuncs_per_elem)


       slave_result = slave_work(parameters, sem_data, inv_mesh, fft_data)

       wt%kernel_values   = slave_result%kernel_values
       wt%kernel_variance = slave_result%kernel_variance
       wt%niterations     = slave_result%niterations
       wt%model           = slave_result%model

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
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function slave_work(parameters, sem_data, inv_mesh, fft_data) result(slave_result)

    use global_parameters,           only: sp, dp, pi, deg2rad, myrank, &
                                           id_fwd, id_bwd, id_fft, id_mc, id_filter_conv, &
                                           id_inv_mesh, id_kernel, id_init

    use inversion_mesh,              only: inversion_mesh_data_type
    use readfields,                  only: semdata_type
    use type_parameter,              only: parameter_type
    use fft,                         only: rfft_type, taperandzeropad
    use filtering,                   only: timeshift_type
    use montecarlo,                  only: integrated_type, allallconverged, allisconverged
    use kernel,                      only: calc_basekernel, calc_physical_kernels
    use backgroundmodel,             only: backgroundmodel_type, nmodel_parameters, &
                                           weight_parameters
    use clocks_mod,                  only: tick

    type(inversion_mesh_data_type), intent(in)    :: inv_mesh
    type(parameter_type),           intent(in)    :: parameters
    type(semdata_type),             intent(in)    :: sem_data
    type(rfft_type),                intent(in)    :: fft_data
    type(slave_result_type)                       :: slave_result

    type(timeshift_type)                :: timeshift_fwd, timeshift_bwd
    type(integrated_type), allocatable  :: int_kernel(:), int_model(:)
    type(backgroundmodel_type)          :: bg_model

    real(kind=dp),    allocatable       :: fw_field(:,:,:)
    real(kind=dp),    allocatable       :: bw_field(:,:,:)
    complex(kind=dp), allocatable       :: fw_field_fd(:,:,:)
    complex(kind=dp), allocatable       :: bw_field_fd(:,:,:)
    complex(kind=dp), allocatable       :: conv_field_fd(:,:), conv_field_fd_filt(:,:)
    real(kind=dp),    allocatable       :: conv_field(:,:)
    real(kind=dp),    allocatable       :: random_points(:,:) 
    real(kind=dp),    allocatable       :: kernelvalue_basekers(:,:,:)
    real(kind=dp),    allocatable       :: kernelvalue_weighted(:,:,:) 
    real(kind=dp),    allocatable       :: kernelvalue_physical(:,:) ! this is linear combs of ...
    integer,          allocatable       :: niterations(:,:)
    real(kind=dp),    allocatable       :: model_random_points(:,:) 

    real(kind=dp)                       :: volume

    character(len=256)                  :: fmtstring
    integer                             :: ielement, irec, ikernel, ibasisfunc, ibasekernel
    integer                             :: nptperstep, ndumps, ntimes, nomega, nelements
    integer                             :: nbasisfuncs_per_elem, nbasekernels
    integer                             :: iclockold
    integer                             :: ndim
    integer, parameter                  :: taper_length = 10! This is the bare minimum. It does 
                                                            ! not produce artifacts yet, at least 
                                                            ! no obvious ones
    integer, parameter                  :: nptperstep_model = 100 !< Points per MC iteration for 
                                                                  !! Integration of model parameters
                                                                  !! Larger than for kernels, since
                                                                  !! evaluation is very cheap

    iclockold = tick()

    ndim = sem_data%get_ndim()
    nptperstep = parameters%npoints_per_step
    ndumps = sem_data%ndumps
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    nbasisfuncs_per_elem = inv_mesh%nbasisfuncs_per_elem
    nelements = inv_mesh%get_nelements()
    nbasekernels = 6

    allocate(slave_result%kernel_values(parameters%nkernel, nbasisfuncs_per_elem, nelements))
    allocate(slave_result%kernel_variance(parameters%nkernel, nbasisfuncs_per_elem, nelements))
    allocate(slave_result%niterations(parameters%nkernel, nelements))
    allocate(slave_result%model(nmodel_parameters, nbasisfuncs_per_elem, nelements))
    
    allocate(fw_field(ndumps, ndim, nptperstep))
    allocate(bw_field(ndumps, ndim, nptperstep))

    allocate(fw_field_fd       (nomega, ndim, nptperstep))
    allocate(bw_field_fd       (nomega, ndim, nptperstep))

    allocate(conv_field        (ntimes, nptperstep))
    allocate(conv_field_fd     (nomega, nptperstep))
    allocate(conv_field_fd_filt(nomega, nptperstep))

    allocate(kernelvalue_basekers(nptperstep, parameters%nkernel, nbasekernels))
    allocate(kernelvalue_weighted(nptperstep, parameters%nkernel, nbasekernels))
    allocate(kernelvalue_physical(nptperstep, parameters%nkernel))

    allocate(model_random_points(nptperstep_model, nmodel_parameters))

    volume = 0.0

    allocate(niterations(parameters%nkernel, nelements))
    niterations = 0

    if (.not.parameters%deconv_stf) then
      call timeshift_fwd%init_ts(fft_data%get_f(), sem_data%timeshift_fwd)
      call timeshift_bwd%init_ts(fft_data%get_f(), sem_data%timeshift_bwd)
    end if
    
!    write(*,*) '***************************************************************'
!    write(*,*) ' Start loop over elements'
!    write(*,*) '***************************************************************'

    allocate(random_points(3, nptperstep))
    allocate(int_kernel(nbasisfuncs_per_elem))
    allocate(int_model(nbasisfuncs_per_elem))

    iclockold = tick(id=id_init, since=iclockold)

    do ielement = 1, nelements 


        ! Get volume of current element
        iclockold = tick()
        volume = inv_mesh%get_volume(ielement)
        iclockold = tick(id=id_inv_mesh)

        ! Calculate integrated model parameters for this element
        do ibasisfunc = 1, nbasisfuncs_per_elem
           call int_model(ibasisfunc)%initialize_montecarlo(nfuncs = nmodel_parameters,   & 
                                                            volume = volume,              &
                                                            allowed_error = 1d-8,         &
                                                            allowed_relerror = 1d-2)
        end do
        do while (.not.allallconverged(int_model))
           random_points = inv_mesh%generate_random_points( ielement, nptperstep_model, &
                                                               parameters%quasirandom)
           do ibasisfunc = 1, nbasisfuncs_per_elem
              model_random_points = weight_parameters(sem_data%load_model_coeffs(random_points), &
                                                      inv_mesh%weights(ielement,                 &  
                                                                       ibasisfunc,               &
                                                                       random_points) )

              call int_model(ibasisfunc)%check_montecarlo_integral(model_random_points)
           end do
        end do

        ! Initialize basis kernel Monte Carlo integrals for current element
        iclockold = tick()
        do ibasisfunc = 1, nbasisfuncs_per_elem
           call int_kernel(ibasisfunc)%initialize_montecarlo(parameters%nkernel, &
                                                             volume,                         &
                                                             parameters%allowed_error,       &
                                                             parameters%allowed_relative_error) 
        end do

        do while (.not.allallconverged(int_kernel)) ! Beginning of Monte Carlo loop


           iclockold = tick(id=id_mc)

           ! Generate random points
           random_points = inv_mesh%generate_random_points( ielement, nptperstep, &
                                                            parameters%quasirandom)
           iclockold = tick(id=id_inv_mesh)
           
           ! Stop MC integration in this element after max_iter iterations
           if (any(niterations(:, ielement)>parameters%max_iter)) exit 

           ! Load forward field
           iclockold = tick()
           fw_field = sem_data%load_fw_points( random_points, parameters%source, & 
                                               model = bg_model)
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
         
           ! Loop over all receivers in receiver input file
           do irec = 1, parameters%nrec
              
              if (allallconverged(int_kernel,[(ikernel, ikernel =                         &
                                               parameters%receiver(irec)%firstkernel,     &
                                               parameters%receiver(irec)%lastkernel) ]))  &
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

              iclockold = tick(id=id_filter_conv, since=iclockold)
                 
              ! Loop over all base kernels
              ! Calculate the waveform kernel in the basic parameters 
              ! (lambda, mu, rho, a, b, c)  
              do ibasekernel = 1, nbasekernels
                 ! Check whether any actual kernel on this receiver needs this base kernel
                 if (.not.parameters%receiver(irec)%needs_basekernel(ibasekernel)) cycle

                 ! Calculate base kernel numero ibasekernel
                 ! TODO: It might be possible that the BWD field is only a trace, to
                 !       save loading the full tensor, if there are only VP kernels at
                 !       this receiver. 
                 !       However, in the moment this is not implemented in readfields.f90
                 conv_field_fd = calc_basekernel(ibasekernel,                           &
                                                 parameters%strain_type_fwd,            &
                                                 parameters%strain_type_fwd,            &
                                                 !parameters%receiver(irec)%strain_type, &
                                                 fw_field_fd, bw_field_fd)
                  

                 iclockold = tick(id=id_filter_conv, since=iclockold)
                 
                 do ikernel = parameters%receiver(irec)%firstkernel, &
                              parameters%receiver(irec)%lastkernel 
                             
                       ! Check whether kernel (ikernel) actually needs the base kernel 
                       ! (ibasekernel), otherwise cycle
                       if (.not.parameters%kernel(ikernel)%needs_basekernel(ibasekernel)) cycle
                        

                       ! If this kernel is already converged, go to the next one
                       if (allisconverged(int_kernel, ikernel)) cycle
                       niterations(ikernel, ielement) = niterations(ikernel, ielement) + 1
                       iclockold = tick(id=id_mc, since=iclockold)

                       ! Apply Filter 
                       conv_field_fd_filt = parameters%kernel(ikernel)%apply_filter(conv_field_fd)
                       iclockold = tick(id=id_filter_conv, since=iclockold)

                       ! Apply Inverse FFT
                       call fft_data%irfft(conv_field_fd_filt, conv_field)
                       iclockold = tick(id=id_fft, since=iclockold)

                       ! Calculate scalar misfit base kernels from convolved time traces
                       kernelvalue_basekers(:,ikernel,ibasekernel) =         &
                                            parameters%kernel(ikernel)       &
                                            %calc_misfit_kernel(conv_field)
                                                                                           
                       iclockold = tick(id=id_kernel, since=iclockold)

                 end do ! ikernel

              end do ! ibasekernel 
              iclockold = tick(id=id_mc, since=iclockold)

           end do ! irec
           iclockold = tick(id=id_mc, since=iclockold)

           ! Check for convergence           
           do ibasisfunc = 1, nbasisfuncs_per_elem
              iclockold = tick()

              ! Compute weighted base kernels
              do ikernel = 1, parameters%nkernel

                 do ibasekernel = 1, nbasekernels
                    ! Check whether kernel (ikernel) actually needs the base kernel 
                    ! (ibasekernel), otherwise cycle
                    if (.not.parameters%kernel(ikernel)%needs_basekernel(ibasekernel)) cycle
                    kernelvalue_weighted(:, ikernel, ibasekernel) =    &
                         kernelvalue_basekers(:, ikernel, ibasekernel) &
                         * inv_mesh%weights(ielement, ibasisfunc, random_points)                 
                 end do
                 
                 ! Compute the kernels for the actual physical parameters of interest
                 kernelvalue_physical(:, ikernel) = calc_physical_kernels(        &
                                      parameters%kernel(ikernel)%model_parameter, &
                                      kernelvalue_weighted(:, ikernel, :),        &
                                      bg_model = bg_model)

              end do

              ! Check if the summed kernels for each basis function have converged
              iclockold = tick(id=id_kernel, since=iclockold)
              call int_kernel(ibasisfunc)%check_montecarlo_integral(kernelvalue_physical)
              iclockold = tick(id=id_mc, since=iclockold)
           end do        
              
           ! Print convergence info
           if (parameters%detailed_convergence) then
              write(fmtstring,"('(I6, I6, ES10.3, A, I0.5, A, I0.5, ', I4, '(ES11.3), A, ', I4, '(ES10.3))')") &
                              parameters%nkernel, parameters%nkernel 
              write(lu_out,fmtstring) myrank, ielement, volume, ' Converged? ',                            &
                                      int_kernel(1)%countconverged(), '/', parameters%nkernel,             &
                                      int_kernel(1)%getintegral(), ' +- ', sqrt(int_kernel(1)%getvariance())
           end if

        end do ! End of Monte Carlo loop
    
        if (any(niterations(:, ielement)>parameters%max_iter)) then
           fmtstring = "('Element', I6, ': Max number of iterations reached. ', I5 , ' kernels were not converged.')"
           write(lu_out, fmtstring) ielement, parameters%nkernel - int_kernel(1)%countconverged()
        else
           fmtstring = "('Element', I6, ': All kernels converged after', I5, ' iterations')"
           write(lu_out, fmtstring) ielement, maxval(niterations(:, ielement))
        end if

        ! Pass over results to master
        do ibasisfunc = 1, nbasisfuncs_per_elem          
           slave_result%kernel_values(:, ibasisfunc, ielement)   = int_kernel(ibasisfunc)%getintegral()
           slave_result%kernel_variance(:, ibasisfunc, ielement) = int_kernel(ibasisfunc)%getvariance()
           slave_result%niterations(:,ielement)                  = niterations(:,ielement)
           slave_result%model(:, ibasisfunc, ielement)           = int_model(ibasisfunc)%getintegral()
           call int_kernel(ibasisfunc)%freeme()
        end do

    end do ! ielement

    call timeshift_fwd%freeme()
    call timeshift_bwd%freeme()
end function slave_work
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
