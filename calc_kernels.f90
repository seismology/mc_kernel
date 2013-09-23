type receiver_type
    real(kind=dp)                        :: lat, lon
    real(kind=dp)                        :: theta, phi
    integer                              :: nkernel
    pointer, dimension(:), allocatable   :: kernel
end type

type kernelspec_type
    real, dimension(2) :: time_window
    real, dimension(4) :: corner_frequencies
    integer            :: filter_type
    character          :: misfit_type
    character          :: model_parameter
    integer            :: receiver_index
    integer            :: src_index
    pointer            :: receiver
    pointer            :: filter
end type

type parameter_type
    type(kdtree_type)  :: fw_tree, bw_tree
    type(netcdf_type)  :: netcdf_info
    type(fft_type)     :: FFT_plan

end type

type(kdtree_type)
    type(kdtree2_result), allocatable    :: nextpoint(:)
    type(kdtree2),pointer                :: tree
end type
    


subroutine make_receiver_kernel_mapping(receiver_info, kernelspec)
    type(kernelspec_type), target, intent(in), dimension(:)   :: kernelspec
    type(receiver_type), intent(inout), dimension(:)  :: receiver_info


    nreceiver = size(receiver_info)
    
    do ikernel = 1, nkernel
        do ireceiver = 1, nreceiver
            if (ireceiver.eq.kernelspec%receiver_index) then
                receiver_info(ireceiver)%nkernel = receiver_info(ireceiver)%nkernel + 1
            end if
        end do
    end do

    do ireceiver = 1, nreceiver
        allocate(receiver_info(ireceiver)%kernel(receiver_info(ireceiver)%nkernel))
    end do

    do ireceiver = 1, nreceiver
        ikernel_rec = 0
        do ikernel = 1, nkernel
            if (ireceiver.eq.kernelspec%receiver_index) then
                ikernel_rec = ikernel_rec + 1
                receiver_info(ireceiver)%kernel(ikernel_rec) => kernelspec(ikernel)
                kernelspec(ikernel)%receiver => receiver_info(ireceiver)
            end if 
        end do
    end do

end subroutine


subroutine calc_kernel(inv_points, receiver_info, parameters)

          
    !type(kernelspec_type), intent(in), dimension(:) :: kernelspec
    type(receiver_type), intent(in), dimension(:)   :: receiver_info
    type(parameter_type), intent(in)                :: parameters

    real(kind=dp), dimension(4), intent(in)         :: inv_points

    type(kernelspec_type)                           :: kernel
    real(kind=dp), dimension(100)                   :: random_points
    integer                                         :: ireceiver, nreceiver, ikernel, nkernel
    logical, dimension(:), allocatable              :: converged

    nreceiver = size(receiver_info)

    do ireceiver = 1, nreceiver

        nkernel = receiver_info(ireceiver)%nkernel

        allocate(converged(nkernel))
        converged = .false.

        seismograms = load_synthetic_seismograms(receiver_info, parameters%netcdf_info)

        do while any(.not.converged) ! Beginning of Monte Carlo loop

            random_point = generate_random_point(inv_points, n)

            fw_fields = load_fw_fields(random_point, parameters%netcdf_info)
            summed_fw_field = sum_over_6_src_comp(fw_fields, src_info%Mij, receiver_info(ireceiver))
            summed_fw_field_fd = FFT(summed_fw_field, parameters%FFT_plan)

            bw_field = load_bw_field(random_point, receiver_info(ireceiver), parameters%netcdf_info)

            bw_field_fd = FFT(bw_field, parameters%FFT_plan)

            wavefield_kernel_fd = convolve(summed_fw_field_fd, bw_field_fd)

            do ikernel = 1, nkernel
                
                kernelspec = receiver_info(ireceiver)%kernel(ikernel)

                if (converged(ikernel)) cycle
                
                filtered_wavefield_kernel_fd = filter(wavefield_kernel_fd, kernelspec%filter)
                filtered_wavefield_kernel_td = iFFT(filtered_wavefield_kernel_fd, parameters%FFT_plan)

                call cut_time_window(filtered_wavefield_kernel_td, kernel(ikernel)%time_window)
                kernel(ikern) = calc_misfit_kernel(filtered_wavefield_kernel_td, kernelspec%misfit)

                integrated_kernel(ikernel) = add_random_point_to_integral(kernel(ikernel), converged(ikernel))

            end do ! iKernel

        end do ! Convergence

        ! Fill up G-matrix

    end do !Receivers

end subroutine
