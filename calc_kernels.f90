program kerner

use global_parameters,           only: sp, dp, pi, deg2rad
use inversion_mesh,              only: inversion_mesh_data_type
use readfields,                  only: semdata_type
use type_parameter,              only: parameter_type
use tetrahedra,                  only: tetra_volume_3d, generate_random_point
use fft,                         only: rfft_type, taperandzeropad
use filtering,                   only: filter_type
use montecarlo,                  only: integrated_type
use kernel,                      only: kernelspec_type

    implicit none
    type(inversion_mesh_data_type)  :: inversion_mesh
    type(parameter_type)            :: parameters
    type(semdata_type)              :: sem_data
    type(rfft_type)                 :: fft_data
    type(filter_type)               :: gabor40, gabor20 
    type(integrated_type)           :: int_kernel
    type(kernelspec_type)           :: kernelspec(2)

    integer                         :: npoints, nelems, ntimes, nomega
    integer                         :: idump, ipoint, ielement, ndumps
    integer                         :: nkernel, ikernel, ivertex, nvertices
    real(kind=sp), allocatable      :: element_points(:,:,:)
    real(kind=sp), allocatable      :: co_points(:,:), K_x(:,:), veloseis(:)
    real(kind=sp), allocatable      :: fw_field(:,:)
    real(kind=sp), allocatable      :: bw_field(:,:)
    complex(kind=dp), allocatable   :: fw_field_fd(:,:)
    complex(kind=dp), allocatable   :: bw_field_fd(:,:)
    complex(kind=dp), allocatable   :: conv_field_fd(:,:), conv_field_fd_filt(:,:)
    real(kind=dp), allocatable      :: conv_field(:,:)
    real(kind=dp), allocatable      :: random_points(:,:), kernelvalue(:,:)
    real(kind=sp)                   :: co_element(3,4)
    real(kind=dp)                   :: volume
    real(kind=sp)                   :: df
    integer, allocatable            :: connectivity(:,:)
    character(len=32)               :: filtername, whattodo

    integer, parameter              :: nptperstep = 10


    call inversion_mesh%read_tet_mesh('vertices.USA10', 'facets.USA10')

    nvertices = inversion_mesh%get_nvertices()
    nelems  = inversion_mesh%get_nelements()

    !allocate(element_points(3,4,nelems))
    !
    !element_points = inversion_mesh%get_elements()
    allocate(co_points(3,nptperstep))
    co_points = inversion_mesh%get_vertices()

    ! Set testing parameters
    parameters%allowed_error = 1e-15

    call parameters%source%init(lat = 30.d0,   &
                                lon = -90.d0,  &
                                mij = dble([0., 0., 0., 0., 1., 0.] ))

    allocate(parameters%receiver(1))

    parameters%receiver(1)%component = 'Z'

    parameters%receiver(1)%latd   = 90
    parameters%receiver(1)%lond   = 0
    parameters%receiver(1)%colatd = 90 - parameters%receiver(1)%latd

    parameters%receiver(1)%colat = parameters%receiver(1)%colatd * deg2rad
    parameters%receiver(1)%lon   = parameters%receiver(1)%lond   * deg2rad

    call parameters%receiver(1)%rotate_receiver(parameters%source)

    parameters%nsim_fwd = 4
    parameters%nsim_bwd = 1

    !parameters%dir_fwdmesh = 'wavefield_example/fwd/'
    !parameters%dir_bwdmesh = 'wavefield_example/bwd/'
    parameters%dir_fwdmesh = 'wavefield_example_10s/fwd'
    parameters%dir_bwdmesh = 'wavefield_example_10s/bwd/'

    ! Test parameters set

    call sem_data%set_params(parameters)
    call sem_data%open_files()
    call sem_data%read_meshes()
    call sem_data%build_kdtree()

    ndumps = sem_data%ndumps

    allocate(fw_field(ndumps, nptperstep))
    allocate(bw_field(ndumps, nptperstep))

    call fft_data%init(ndumps, nptperstep, sem_data%dt)
    ntimes = fft_data%get_ntimes()
    nomega = fft_data%get_nomega()
    df     = fft_data%get_df()
    print *, 'ntimes: ', ntimes, ', nomega: ', nomega
    print *, 'dt: ', sem_data%dt, ', df: ', df

    allocate(conv_field   (ntimes, nptperstep))
    allocate(fw_field_fd  (nomega, nptperstep))
    allocate(bw_field_fd  (nomega, nptperstep))
    allocate(conv_field_fd(nomega, nptperstep))
    allocate(conv_field_fd_filt(nomega, nptperstep))

    write(*,*) ' Define filters and kernelspecs'
    filtername = 'Gabor'
    call gabor40%create(df, nomega, filtername, [40.0, 0.5, 0., 0.])

    call kernelspec(1)%init(time_window     = [590.0, 650.0], &
                            filter          = gabor40,        &
                            misfit_type     = 'CC  ',         &  
                            model_parameter = 'vp  ')

    call gabor20%create(df, nomega, filtername, [20.0, 0.5, 0., 0.])

    call kernelspec(2)%init(time_window     = [590.0, 620.0], &
                            filter          = gabor20,        &
                            misfit_type     = 'CC  ',         &  
                            model_parameter = 'vp  ')

    whattodo = 'integratekernel'
    nkernel = 2
    select case(trim(whattodo))
    case('integratekernel')
       write(*,*) 'Initialize output file'
       call inversion_mesh%init_data(nkernel)
       
       write (*,*) 'Initialize Kernel variables'
       allocate(K_x(nvertices, nkernel))
       allocate(kernelvalue(nptperstep, nkernel))
       allocate(connectivity(4, nelems))
       allocate(veloseis(ndumps))
       
       write(*,*) 'Loading Connectivity'
       connectivity = inversion_mesh%get_connectivity()
   
       !do ivertex = 1, nvertices
       !   nelems  = inversion_mesh%get_valence(ivertex)
       !   elems   = inversion_mesh%get_connected_elements(ivertex)
       !   volume  = 0.0
       !   allocate(random_points(3, npoints * nvertices))
       !   do ielement = 1, nelems
       !      co_element = inversion_mesh%get_element(elems(ielement))
       !      volume = volume + tetra_volume_3d(dble(co_element))

       !      random_points((ielement-1)*npoints+1:ielement*npoints) = &
       !                   generate_random_point(dble(co_element), npoints)
       !   end do

       veloseis = sem_data%load_seismogram(parameters%receiver(1), parameters%source)

       write (*,*) 'Start loop over elements'
       do ielement = 1, nelems
          co_element = inversion_mesh%get_element(ielement)
          volume = volume + tetra_volume_3d(dble(co_element))
          call int_kernel%initialize_montecarlo(nkernel, volume, parameters%allowed_error) 
         
          do while (.not.int_kernel%areallconverged()) ! Beginning of Monte Carlo loop
             random_points =  generate_random_point(dble(co_element), nptperstep)
             fw_field = sem_data%load_fw_points(random_points, parameters%source)
             bw_field = sem_data%load_bw_points(random_points, parameters%receiver(1))
             call fft_data%rfft(taperandzeropad(fw_field, ntimes), fw_field_fd)
             call fft_data%rfft(taperandzeropad(bw_field, ntimes), bw_field_fd)
             
             conv_field_fd = fw_field_fd * bw_field_fd

             do ikernel = 1, nkernel
                conv_field_fd_filt = kernelspec(ikernel)%filter%apply_2d(conv_field_fd)

                call fft_data%irfft(conv_field_fd_filt, conv_field)

                kernelvalue(:,ikernel) = kernelspec(ikernel)%calc_misfit_kernel(fft_data%get_t(), &
                                                                                conv_field,       &
                                                                                dble(veloseis))
             end do
             call int_kernel%check_montecarlo_integral(kernelvalue)
             print *, 'Converged ? ', int_kernel%areallconverged(), int_kernel%getintegral(), &
                      '+-', sqrt(int_kernel%getvariance())

          end do
         
          do ivertex = 1, 4
             K_x(connectivity(ivertex, ielement),:) = K_x(connectivity(ivertex, ielement),:) + int_kernel%getintegral()
          end do
          if (mod(ielement, 100)==0) then
             write(*,*) 'Write Kernel to disk'
             call inversion_mesh%set_data_snap(K_x(:,1), 1, 'gabor_20')
             call inversion_mesh%set_data_snap(K_x(:,2), 2, 'gabor_40')
             call inversion_mesh%dump_mesh_data_xdmf('gaborkernel')
          end if
       end do
   
    case('plot_wavefield')

       write(*,*) ' Read in forward field'
       fw_field = sem_data%load_fw_points(dble(co_points), parameters%source)
       write(*,*) ' Read in backward field'
       bw_field = sem_data%load_bw_points(dble(co_points), parameters%receiver(1))
       write(*,*) ' FFT forward field'
       call fft_data%rfft(taperandzeropad(fw_field, ntimes), fw_field_fd)
       write(*,*) ' FFT forward field'
       call fft_data%rfft(taperandzeropad(bw_field, ntimes), bw_field_fd)

       write(*,*) ' Apply filter'
       fw_field_fd = gabor20%apply_2d(fw_field_fd)

       write(*,*) ' iFFT product of fields'
       call fft_data%irfft(fw_field_fd*bw_field_fd, conv_field)
    
       write(*,*) ' Write pickpoint to disc'
       write(41,*) fw_field(:,10000)
       write(42,*) bw_field(:,10000)
       write(43,*) conv_field(:,10000)
       read(*,*)

       ! Dump to XDMF file
       call inversion_mesh%init_data(ndumps*2 + ntimes)
       do idump = 1, ndumps
           write(*,*) ' Passing dump ', idump, ' to inversion mesh datatype'
           !Test of planar wave , works
           !fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
           !bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)
           call inversion_mesh%set_data_snap(fw_field(idump,:), idump, 'fwd_wavefield')
           call inversion_mesh%set_data_snap(bw_field(idump,:), idump+ndumps, 'bwd_wavefield')
       end do

       do idump = 1, ntimes
          write(*,*) ' Passing dump ', idump, ' of convolved wavefield'
          call inversion_mesh%set_data_snap(real(conv_field(idump,:)), idump+ndumps*2, 'field_convolved')
       end do 

       write(*,*)
       write(*,*) ' Writing data to disk'
       call inversion_mesh%dump_mesh_data_xdmf('wavefield')

    end select
    write(*,*)
    write(*,*) ' Free memory of inversion mesh datatype'
    call inversion_mesh%freeme

    call sem_data%close_files()

    write(*,*)
    write(*,*) ' Finished!'

    
    return




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
