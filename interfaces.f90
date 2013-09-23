function load_synthetic_seismograms(receiver_info, netcdf_info)
    type(receiver_type), intent(in) :: receiver_info
    type(netcdf_type), intent(in)   :: netcdf_info
    real(kind=dp), dimension(:,3)   :: load_synthetic_seismograms
end function


function load_fw_fields(coordinate, netcdf_info)
    real(kind=sp), dimension(3,100), intent(in)    :: coordinate
    type(netcdf_type) , intent(in)                 :: netcdf_info
    real(kind=sp), dimension(3,100,:), intent(out) :: load_fw_fields
end function

function load_bw_fields(coordinate, receiver_info, netcdf_info)
    real(kind=sp), dimension(3,100), intent(in)    :: coordinate
    type(receiver_type), intent(in)                :: receiver_info
    type(netcdf_type) , intent(in)                 :: netcdf_info
    real(kind=sp), dimension(3,100,:), intent(out) :: load_fw_fields
end function

function FFT(field, fft_plan)
    real(kind=sp), dimension(:,:), intent(in)      :: field
    type(fft_type), intent(in)                     :: fft_plan
    real(kind=sp), dimension(:,:), intent(out)     :: FFT
end function

function iFFT(field, fft_plan)
    real(kind=sp), dimension(:,:), intent(in)      :: field
    type(fft_type), intent(in)                     :: fft_plan
    real(kind=sp), dimension(:,:), intent(out)     :: iFFT
end function
