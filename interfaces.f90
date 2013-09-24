

function load_synthetic_seismograms(receiver_info, netcdf_info)
    type(receiver_type), intent(in)                       :: receiver_info
    type(netcdf_type), intent(in)                         :: netcdf_info
    real(kind=dp), dimension(:,3)                         :: load_synthetic_seismograms
end function

function load_fw_fields(coordinate, netcdf_info)
    real(kind=sp), dimension(:,:), intent(in)             :: coordinate
    type(netcdf_type) , intent(in)                        :: netcdf_info
    real(kind=sp), dimension(:,:,:), intent(out)          :: load_fw_fields
end function

function load_bw_fields(coordinate, receiver_info, netcdf_info)
    real(kind=sp), dimension(:,:), intent(in)             :: coordinate
    type(receiver_type), intent(in)                       :: receiver_info
    type(netcdf_type) , intent(in)                        :: netcdf_info
    real(kind=sp), dimension(:,:,:), intent(out)          :: load_fw_fields
end function

function sum_over_6_src_comp(fields, mij, receiver_info)
    real(kind=sp), dimension(:,:,:), intent(in)           :: fields
    real(kind=sp), dimension(6), intent(in)               :: mij
    type(receiver_type), intent(in)                       :: receiver_info
    real(kind=sp), dimension(:,:), intent(out)            :: sum_over_6_src_comp
end function

function FFT(field, fft_plan)
    real(kind=sp), dimension(:,:), intent(in)             :: field
    type(fft_type), intent(in)                            :: fft_plan
    real(kind=sp), dimension(size(field)), intent(out)    :: FFT
end function

function iFFT(field, fft_plan)
    real(kind=sp), dimension(:,:), intent(in)             :: field
    type(fft_type), intent(in)                            :: fft_plan
    real(kind=sp), dimension(size(field)), intent(out)    :: iFFT
end function

function convolve(onefield, anotherfield)
    real(kind=sp), dimension(:,:), intent(in)             :: onefield, anotherfield
    real(kind=sp), dimension(size(onefield)), intent(out) :: convolve
end function

function filter(timeseries, filterseries)
    real(kind=sp), dimension(:,:), intent(in)             :: timeseries
    real(kind=sp), dimension(:), intent(in)               :: filterseries
    real(kind=sp), dimension(size(timeseries)), intent(out) :: filter
end function

function cut_time_window(timeseries, time_window) 
    real(kind=sp), dimension(:,:), intent(in)             :: timeseries
    real(kind=sp), dimension(2)                           :: time_window
    real(kind=sp), allocatable, dimension(:,:), intent(out) :: cut_time_window
end function

function calc_misfit_kernel(wavefield_kernel, misfit_criterion)
    real(kind=sp), dimension(:,:). intent(in)             :: wavefield_kernel
    character, (len=16), intent(in)                       :: misfit_criterion
    real(kind=sp), dimension(:), intent(out)              :: calc_misfit_kernel
end function

