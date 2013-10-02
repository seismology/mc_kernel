program kerner

use global_parameters
use inversion_mesh
use readfields

    implicit none
    type(inversion_mesh_data_type)  :: inversion_mesh
    !type(inversion_mesh_type)       :: inversion_mesh
    type(parameter_type)            :: parameters
    type(netcdf_type)               :: sem_data

    integer                         :: npoints, nelems, ndumps
    integer                         :: idump, ipoint
    real(kind=sp), allocatable      :: element_points(:,:,:)
    real(kind=sp), allocatable      :: co_points(:,:)
    real(kind=sp), allocatable      :: fw_field(:,:)
    real(kind=sp), allocatable      :: bw_field(:,:)


    call inversion_mesh%read_tet_mesh('vertices.USA10', 'facets.USA10')

    npoints = inversion_mesh%get_nvertices()
    nelems  = inversion_mesh%get_nelements()
    ndumps  = 200
    
    !allocate(element_points(3,4,nelems))
    !
    !element_points = inversion_mesh%get_elements()
    allocate(co_points(3,npoints))
    co_points = inversion_mesh%get_vertices()

    parameters%allowed_error = 1e-3

    parameters%source%lat = 50
    parameters%source%lon = -90
    parameters%source%colat = 90 - parameters%source%lat
    parameters%source%mij = [0, 0, 0, 0, 1, 0]

    allocate(parameters%receiver(1))
    parameters%receiver(1)%component = 'Z'

    parameters%receiver(1)%lat = -30
    parameters%receiver(1)%colat = 90 - parameters%receiver(1)%lat
    parameters%receiver(1)%lon = 120


    parameters%nsim_fwd = 4
    parameters%nsim_bwd = 1

    parameters%dir_fwdmesh = 'wavefield_example/fwd/'
    parameters%dir_bwdmesh = 'wavefield_example/bwd/'

    call sem_data%set_params(parameters)
    call sem_data%open_files()
    call sem_data%read_meshes()
    call sem_data%build_kdtree()

    allocate(fw_field(ndumps, npoints))
    allocate(bw_field(ndumps, npoints))

    !fw_field = sem_data%load_fw_points(dble(co_points), parameters%source)
    !bw_field = sem_data%load_bw_points(dble(co_points), parameters%receiver(1))
    
    call inversion_mesh%init_data(ndumps)
    do idump = 1, ndumps
        write(*,*) ' Passing dump ', idump, ' to inversion mesh datatype'
        !Test of planar wave , works
        fw_field(idump,:) = sin(co_points(1,:)/1000 + idump*0.1)
        bw_field(idump,:) = sin(co_points(2,:)/1000 + idump*0.1)
        call inversion_mesh%set_data_snap(fw_field(idump,:), idump, 'fwd_wavefield')
        call inversion_mesh%set_data_snap(bw_field(idump,:), idump, 'bwd_wavefield')
    end do
    write(*,*) ' Writing data to disk'
    call inversion_mesh%dump_tet_mesh_data_xdmf('wavefield')

    write(*,*) ' Free memory of inversion mesh datatype'
    call inversion_mesh%freeme



    
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
