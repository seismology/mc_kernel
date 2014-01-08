program kerner

    use global_parameters,           only: sp, dp, pi, deg2rad, verbose, init_random_seed
    use inversion_mesh,              only: inversion_mesh_data_type
    use readfields,                  only: semdata_type
    use type_parameter,              only: parameter_type
    use tetrahedra,                  only: tetra_volume_3d, generate_random_point
    use fft,                         only: rfft_type, taperandzeropad
    use filtering,                   only: timeshift
    use montecarlo,                  only: integrated_type
    !use kernel,                      only: kernelspec_type

    use ftnunit,                     only: runtests_init, runtests, runtests_final
    use unit_tests,                  only: test_all

    implicit none
    type(inversion_mesh_data_type)      :: inv_mesh
    type(parameter_type)                :: parameters
    type(semdata_type)                  :: sem_data
    type(rfft_type)                     :: fft_data
    !type(filter_type)                   :: gabor40, gabor20, gabor10
    type(integrated_type)               :: int_kernel
    !type(kernelspec_type), allocatable  :: kernelspec(:)

    integer                             :: npoints, nelems, ntimes, nomega
    integer                             :: idump, ipoint, ielement, ndumps, irec
    integer                             :: nkernel, ikernel, ivertex, nvertices
    real(kind=sp),    allocatable       :: element_points(:,:,:)
    real(kind=sp),    allocatable       :: co_points(:,:), K_x(:,:), veloseis(:)
    real(kind=sp),    allocatable       :: fw_field(:,:)
    real(kind=sp),    allocatable       :: bw_field(:,:)
    complex(kind=dp), allocatable       :: fw_field_fd(:,:)
    complex(kind=dp), allocatable       :: bw_field_fd(:,:)
    complex(kind=dp), allocatable       :: conv_field_fd(:,:), conv_field_fd_filt(:,:)
    real(kind=dp),    allocatable       :: conv_field(:,:)
    real(kind=dp),    allocatable       :: random_points(:,:), kernelvalue(:,:)
    integer,          allocatable       :: connectivity(:,:), niterations(:,:)
    real(kind=sp)                       :: co_element(3,4)
    real(kind=dp)                       :: volume
    real(kind=dp)                       :: df
    character(len=32)                   :: filtername, whattodo
    character(len=64)                   :: fmtstring

    integer, parameter                  :: nptperstep = 25, lu_iterations = 400

    verbose = 0

    call init_random_seed()

    call runtests_init
    call runtests( test_all )
    call runtests_final

    verbose = 1

    write(*,*) '***************************************************************'
    write(*,*) ' Read input files for parameters, source and receivers'
    write(*,*) '***************************************************************'
    call parameters%read_parameters('inparam_basic')
    call parameters%read_source()
    call parameters%read_receiver()


    write(*,*) '***************************************************************'
    write(*,*) ' Initialize and open AxiSEM wavefield files'
    write(*,*) '***************************************************************'
    call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir)
    call sem_data%open_files()
    call sem_data%read_meshes()
    call sem_data%build_kdtree()

    call sem_data%load_seismogram(parameters%receiver, parameters%source)

    ndumps = sem_data%ndumps


    write(*,*) '***************************************************************'
    write(*,*) ' Read inversion mesh'
    write(*,*) '***************************************************************'
    call inv_mesh%read_tet_mesh('vertices.USA10', 'facets.USA10')

    nvertices = inv_mesh%get_nvertices()
    nelems    = inv_mesh%get_nelements()
    fmtstring = '(A, I8, A, I8)'
    print fmtstring, '  nvertices: ',  nvertices, ', nelems: ', nelems

    write(*,*) '***************************************************************'
    write(*,*) ' Initialize FFT'
    write(*,*) '***************************************************************'
    call fft_data%init(ndumps, nptperstep, sem_data%dt)
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    df     = fft_data%get_df()
    fmtstring = '(A, I8, A, I8)'
    print fmtstring, '  ntimes: ',  ntimes,     '  , nfreq: ', nomega
    fmtstring = '(A, F8.3, A, F8.3, A)'
    print fmtstring, '  dt:     ', sem_data%dt, ' s, df:    ', df*1000, ' mHz'


    write(*,*) '***************************************************************'
    write(*,*) ' Define filters'
    write(*,*) '***************************************************************'
    
    call parameters%read_filter(nomega, df)


    write(*,*) '***************************************************************'
    write(*,*) ' Define kernels'
    write(*,*) '***************************************************************'

    call parameters%read_kernel(sem_data, parameters%filter)



    !allocate(parameters%receiver(1))

    !!ANMO       IU       34.9459   -106.4572    1720.0   100.0
    !! Distance 112.577
    !call parameters%receiver(1)%init(lat       = 34.9459 , &
    !                                 lon       = -106.4572, &
    !                                 component = 'Z')
    !! Station Pitcairn PTCN
    !call parameters%receiver(1)%init(lat       = -25.08 , &
    !                                 lon       = 229.890, &
    !                                 component = 'Z')




    ! Test parameters set

    !write(*,*) '***************************************************************'
    !write(*,*) ' Define kernelspecs'
    !write(*,*) '***************************************************************'
    !nkernel = 08
    !allocate(kernelspec(nkernel))
    !allocate(veloseis(ndumps))
    !veloseis = sem_data%load_seismogram(parameters%receiver(1), parameters%source)

    !call kernelspec(1)%init(name            = 'P_40s           ',     &
    !                        time_window     = [675.0, 735.0],         &
    !                        filter          = gabor40,                &
    !                        misfit_type     = 'CC  ',                 &  
    !                        model_parameter = 'vp  ',                 &
    !                        veloseis        = veloseis,               &
    !                        dt              = sem_data%dt,            &
    !                        timeshift_fwd   = sem_data%timeshift_fwd)

    !call kernelspec(2)%init(name            = 'P_20s           ',     &
    !                        time_window     = [675.0, 705.0],         &
    !                        filter          = gabor20,                &
    !                        misfit_type     = 'CC  ',                 &  
    !                        model_parameter = 'vp  ',                 &
    !                        veloseis        = veloseis,               &
    !                        dt              = sem_data%dt,            &
    !                        timeshift_fwd   = sem_data%timeshift_fwd)

    !call kernelspec(3)%init(name            = 'P_10s           ',     &
    !                        time_window     = [675.0, 690.0],         &
    !                        filter          = gabor10,                &
    !                        misfit_type     = 'CC  ',                 &  
    !                        model_parameter = 'vp  ',                 &
    !                        veloseis        = veloseis,               &
    !                        dt              = sem_data%dt,            &
    !                        timeshift_fwd   = sem_data%timeshift_fwd)

    !call kernelspec(4)%init(name            = 'PP_10s          ',     &
    !                        time_window     = [837.00,0852.00],       &
    !                        filter          = gabor10,                &
    !                        misfit_type     = 'CC  ',                 &  
    !                        model_parameter = 'vp  ',                 &
    !                        veloseis        = veloseis,               &
    !                        dt              = sem_data%dt,            &
    !                        timeshift_fwd   = sem_data%timeshift_fwd)

    !call kernelspec(5)%init(name            = 'PP_20s          ',     &
    !                        time_window     = [837.00, 0867.0],       &
    !                        filter          = gabor20,                &
    !                        misfit_type     = 'CC  ',                 &  
    !                        model_parameter = 'vp  ',                 &
    !                        veloseis        = veloseis,               &
    !                        dt              = sem_data%dt,            &
    !                        timeshift_fwd   = sem_data%timeshift_fwd)

    !call kernelspec(6)%init(name            = 'PP_40s          ',     &
    !                        time_window     = [837.00, 877.00],       &
    !                        filter          = gabor40,                &
    !                        misfit_type     = 'CC  ',                 &  
    !                        model_parameter = 'vp  ',                 &
    !                        veloseis        = veloseis,               &
    !                        dt              = sem_data%dt,            &
    !                        timeshift_fwd   = sem_data%timeshift_fwd)

    !call kernelspec(7)%init(name            = 'PPP_10s         ',     &
    !                        time_window     = [944.00, 0959.0],       &
    !                        filter          = gabor10,                &
    !                        misfit_type     = 'CC  ',                 &  
    !                        model_parameter = 'vp  ',                 &
    !                        veloseis        = veloseis,               &
    !                        dt              = sem_data%dt,            &
    !                        timeshift_fwd   = sem_data%timeshift_fwd)

    !call kernelspec(8)%init(name            = 'PPP_20s         ',     &
    !                        time_window     = [944.00, 974.00],       &
    !                        filter          = gabor20,                &
    !                        misfit_type     = 'CC  ',                 &  
    !                        model_parameter = 'vp  ',                 &
    !                        veloseis        = veloseis,               &
    !                        dt              = sem_data%dt,            &
    !                        timeshift_fwd   = sem_data%timeshift_fwd)

!    call kernelspec(9)%init(name            = 'PPP_20s         ',     &
!                            time_window     = [2215.0, 2245.0],       &
!                            filter          = gabor20,                &
!                            misfit_type     = 'CC  ',                 &  
!                            model_parameter = 'vp  ',                 &
!                            veloseis        = veloseis,               &
!                            dt              = sem_data%dt,            &
!                            timeshift_fwd   = sem_data%timeshift_fwd)
!
!    call kernelspec(10)%init(name           = 'PPP_40s         ',     &
!                            time_window     = [2215.0, 2275.0],       &
!                            filter          = gabor40,                &
!                            misfit_type     = 'CC  ',                 &  
!                            model_parameter = 'vp  ',                 &
!                            veloseis        = veloseis,               &
!                            dt              = sem_data%dt,            &
!                            timeshift_fwd   = sem_data%timeshift_fwd)

    whattodo = 'integratekernel'

    select case(trim(whattodo))
    case('integratekernel')
       write(*,*) '***************************************************************'
       write(*,*) 'Initialize output file'
       write(*,*) '***************************************************************'
       call inv_mesh%init_data(nkernel)
       
       write(*,*) '***************************************************************'
       write (*,*) 'Initialize Kernel variables'
       write(*,*) '***************************************************************'
       allocate(niterations(nkernel, nelems))
       niterations = 0
       open(unit=lu_iterations, file='niterations.txt', action='write')
       allocate(K_x(nvertices, nkernel))
       allocate(kernelvalue(nptperstep, nkernel))
       allocate(connectivity(4, nelems))
       
       allocate(fw_field(ndumps, nptperstep))
       allocate(bw_field(ndumps, nptperstep))

       allocate(conv_field   (ntimes, nptperstep))
       allocate(fw_field_fd  (nomega, nptperstep))
       allocate(bw_field_fd  (nomega, nptperstep))
       allocate(conv_field_fd(nomega, nptperstep))
       allocate(conv_field_fd_filt(nomega, nptperstep))

       write(*,*) '***************************************************************'
       write(*,*) ' Loading Connectivity'
       write(*,*) '***************************************************************'
       connectivity = inv_mesh%get_connectivity()


       write(*,*) '***************************************************************'
       write(*,*) ' Start loop over elements'
       write(*,*) '***************************************************************'

       do ielement = 1, nelems

          co_element = inv_mesh%get_element(ielement)
          volume = tetra_volume_3d(dble(co_element))

          ! Omit elements in the core
          if (all( sum(co_element**2, 1).lt.2890.0**2 )) then
             print '(I6,E11.3,A)', ielement, volume, ' in core, skipping'
             print '(4(" ", F6.1))', sqrt(sum(co_element**2, 1))
             cycle
          end if
          
          ! Initialize Monte Carlo integral for this element
          call int_kernel%initialize_montecarlo(nkernel, volume, parameters%allowed_error) 
         
          do while (.not.int_kernel%areallconverged()) ! Beginning of Monte Carlo loop
             random_points = generate_random_point( dble(co_element), nptperstep )
             
             ! Load, FT and timeshift forward field
             fw_field = sem_data%load_fw_points( random_points, parameters%source )
             call fft_data%rfft( taperandzeropad(fw_field, ntimes), fw_field_fd )
             call timeshift( fw_field_fd, fft_data%get_f(), sem_data%timeshift_fwd )
           
             do irec = 1, parameters%nrec

                if (int_kernel%areallconverged([ (ikernel, ikernel =                        &
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
                             parameters%receiver(irec)%lastkernel !, nkernel
                   if (int_kernel%isconverged(ikernel)) cycle
                   !niterations(ikernel, ielement) = niterations(ikernel, ielement) + 1

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
             call int_kernel%check_montecarlo_integral(kernelvalue)
             
             ! Print convergence info
             fmtstring = '(I6, E10.3, A, L1, 8(E11.3), A, 8(E10.3))'
             print fmtstring, ielement, volume, ' Converged? ', int_kernel%areallconverged(), &
                              int_kernel%getintegral(), ' +- ', sqrt(int_kernel%getvariance())

          end do ! Kernel coverged
         
          ! Save integral values to the big kernel variable
          do ivertex = 1, 4
             K_x(connectivity(ivertex, ielement),:) = K_x(connectivity(ivertex, ielement),:) + int_kernel%getintegral() / volume
          end do

          call int_kernel%freeme()

          ! Save big kernel variable to disk
          if (mod(ielement, 100)==0) then
             write(*,*) 'Write Kernel to disk'
             do ikernel = 1, nkernel
                call inv_mesh%set_data_snap(K_x(:,ikernel), ikernel, parameters%kernel(ikernel)%name )
             end do

             call inv_mesh%dump_mesh_data_xdmf('gaborkernel')
          end if

          write(lu_iterations,*) ielement, niterations(:, ielement)
       end do ! ielement
       close(lu_iterations)
   
    !case('plot_wavefield')

    !   print *, 'Initialize XDMF file'
    !   allocate(co_points(3, nvertices))
    !   co_points = inversion_mesh%get_vertices()
    !   call inversion_mesh%init_data(ndumps*2 + ntimes)

    !   write(*,*) ' Read in forward field'
    !   allocate(fw_field(ndumps, nvertices))
    !   fw_field = sem_data%load_fw_points(dble(co_points), parameters%source)
    !   
    !   ! Dump forward field to XDMF file
    !   do idump = 1, ndumps
    !       write(*,*) ' Passing dump ', idump, ' to inversion mesh datatype'
    !       !Test of planar wave , works
    !       !fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
    !       !bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)
    !       call inversion_mesh%set_data_snap(fw_field(idump,:), idump, 'fwd_wavefield')
    !   end do
    !   write(*,*) ' FFT forward field'
    !   allocate(fw_field_fd  (nomega, nptperstep))
    !   call fft_data%rfft(taperandzeropad(fw_field, ntimes), fw_field_fd)
    !   deallocate(fw_field)

    !   write(*,*) ' Read in backward field'
    !   allocate(bw_field(ndumps, nvertices))
    !   bw_field = sem_data%load_bw_points(dble(co_points), parameters%receiver(1))
    !   ! Dump backward field to XDMF file
    !   do idump = 1, ndumps
    !       write(*,*) ' Passing dump ', idump, ' to inversion mesh datatype'
    !       !Test of planar wave , works
    !       !fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
    !       !bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)
    !       call inversion_mesh%set_data_snap(bw_field(idump,:), idump+ndumps, 'bwd_wavefield')
    !   end do
    !   write(*,*) ' FFT backward field'
    !   allocate(bw_field_fd  (nomega, nptperstep))
    !   call fft_data%rfft(taperandzeropad(bw_field, ntimes), bw_field_fd)
    !   deallocate(bw_field)

    !   write(*,*) ' Apply filter'
    !   fw_field_fd = gabor20%apply_2d(fw_field_fd)

    !   write(*,*) ' iFFT product of fields'
    !   allocate(conv_field   (ntimes, nptperstep))
    !   call fft_data%irfft(fw_field_fd*bw_field_fd, conv_field)
    !
    !   !write(*,*) ' Write pickpoint to disc'
    !   !write(41,*) fw_field(:,10000)
    !   !write(42,*) bw_field(:,10000)
    !   !write(43,*) conv_field(:,10000)
    !   !read(*,*)


    !   do idump = 1, ntimes
    !      write(*,*) ' Passing dump ', idump, ' of convolved wavefield'
    !      call inversion_mesh%set_data_snap(real(conv_field(idump,:)), idump+ndumps*2, 'field_convolved')
    !   end do 

    !   write(*,*)
    !   write(*,*) ' Writing data to disk'
    !   call inversion_mesh%dump_mesh_data_xdmf('wavefield')
    !   call inversion_mesh%freeme()

    end select

    write(*,*)
    write(*,*) '***************************************************************'
    write(*,*) ' Free memory of Kernelspecs'
    write(*,*) '***************************************************************'
    do ikernel = 1, nkernel
       call parameters%kernel(ikernel)%freeme() 
    end do

    write(*,*)
    write(*,*) '***************************************************************'
    write(*,*) ' Free memory of inversion mesh datatype'
    write(*,*) '***************************************************************'
    call inv_mesh%freeme

    write(*,*)
    write(*,*) '***************************************************************'
    write(*,*) ' Close AxiSEM wavefield files'
    write(*,*) '***************************************************************'
    call sem_data%close_files()

    write(*,*)
    write(*,*) '***************************************************************'
    write(*,*) ' Free memory of FFT datatype'
    write(*,*) '***************************************************************'
    call fft_data%freeme()

    write(*,*)
    write(*,*) ' Finished!'

    
contains


!subroutine make_receiver_kernel_mapping(receiver_info, kernelspec)
!    type(kernelspec_type), target, intent(in), dimension(:)   :: kernelspec
!    type(receiver_type), intent(inout), dimension(:)  :: receiver_info
!
!
!    nreceiver = size(receiver_info)
!    
!    do ikernel = 1, nkernel
!        do ireceiver = 1, nreceiver
!            if (ireceiver.eq.kernelspec%receiver_index) then
!                receiver_info(ireceiver)%nkernel = receiver_info(ireceiver)%nkernel + 1
!            end if
!        end do
!    end do
!
!    do ireceiver = 1, nreceiver
!        allocate(receiver_info(ireceiver)%kernel(receiver_info(ireceiver)%nkernel))
!    end do
!
!    do ireceiver = 1, nreceiver
!        ikernel_rec = 0
!        do ikernel = 1, nkernel
!            if (ireceiver.eq.kernelspec%receiver_index) then
!                ikernel_rec = ikernel_rec + 1
!                receiver_info(ireceiver)%kernel(ikernel_rec) => kernelspec(ikernel)
!                kernelspec(ikernel)%receiver => receiver_info(ireceiver)
!            end if 
!        end do
!    end do
!
!end subroutine
!
!
!subroutine calc_kernel(inv_points, receiver_info, parameters)
!
!          
!    !type(kernelspec_type), intent(in), dimension(:) :: kernelspec
!    type(netcdf_type), intent(in)                   :: sem_data
!    type(receiver_type), intent(in), dimension(:)   :: receiver_info
!    type(parameter_type), intent(in)                :: parameters
!    type(integrated_type)                           :: int_kernel
!
!    real(kind=dp), dimension(4), intent(in)         :: inv_points
!
!    type(kernelspec_type)                           :: kernel
!    type(sem_type)                                  :: sem_data
!    real(kind=dp), dimension(100)                   :: random_points
!    integer                                         :: ireceiver, nreceiver, ikernel, nkernel
!    logical, dimension(:), allocatable              :: converged
!
!    nreceiver = size(receiver_info)
!    ! Has to be called before: 
!    ! call sem_data%set_params
!    ! call sem_data%read_meshes
!    ! call sem_data%build_tree
!
!
!    do ireceiver = 1, nreceiver
!
!        nkernel = receiver_info(ireceiver)%nkernel
!
!        call int_kernel%initialize_montecarlo(nkernel, volume, parameters%allowed_error) 
!
!        seismograms = sem_data%load_seismograms(parameters%receiver(ireceiver))
!
!        do while (.not.int_kernel%areallconverged()) ! Beginning of Monte Carlo loop
!
!            random_points = generate_random_point(inv_points, n)
!
!            !fw_fields = load_fw_fields(random_point, parameters%netcdf_info)
!            fw_fields = sem_data%load_fw_points(random_points, parameters%source)
!
!            summed_fw_field_fd = FFT(summed_fw_field, parameters%FFT_plan)
!
!            bw_fields = sem_data%load_bw_points(random_points, parameters%receiver(ireceiver))
!            !bw_field = load_bw_field(random_point, receiver_info(ireceiver), parameters%netcdf_info)
!
!            bw_field_fd = FFT(bw_field, parameters%FFT_plan)
!
!            wavefield_kernel_fd = convolve(summed_fw_field_fd, bw_field_fd)
!
!            do ikernel = 1, nkernel
!                
!                kernelspec = parameters%receiver(ireceiver)%kernel(ikernel)
!
!                if (int_kernel%isconverged(ikernel)) cycle
!                
!                filtered_wavefield_kernel_fd = filter(wavefield_kernel_fd, kernelspec%filter)
!                filtered_wavefield_kernel_td = iFFT(filtered_wavefield_kernel_fd, parameters%FFT_plan)
!
!                cut_filtered_wavefield_kernel_td = cut_time_window(filtered_wavefield_kernel_td, kernel(ikernel)%time_window)
!                kernelvalues(ikernel) = calc_misfit_kernel(filtered_wavefield_kernel_td, seismograms, kernelspec%misfit)
!
!
!            end do ! iKernel
!            call int_kernel%check_montecarlo_integral(kernelvalues)
!
!        end do ! Convergence
!
!        ! Fill up G-matrix
!
!    end do !Receivers
!
!end subroutine
!

end program
