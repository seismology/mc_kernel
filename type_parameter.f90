module type_parameter
    use global_parameters,               only : sp, dp
    
    type src_param_type
        real(kind=sp)                        :: mij(6)
        real(kind=dp)                        :: colat, lat, lon
        real(kind=dp)                        :: colatd, latd, lond
    end type

    type rec_param_type
        character(len=1)                     :: component
        real(kind=dp)                        :: colat, lat, lon
        real(kind=dp)                        :: colatd, latd, lond
        integer                              :: nkernel
        type(kernelspec_type), pointer       :: kernel(:)
    end type

    type parameter_type
        type(src_param_type)                 :: source
        type(rec_param_type), allocatable    :: receiver(:)

        real(kind=sp)                        :: allowed_error
        character(len=512)                   :: dir_fwdmesh
        character(len=512)                   :: dir_bwdmesh
        integer                              :: nsim_fwd, nsim_bwd
    end type

    type kernelspec_type
        real, dimension(2)                   :: time_window
        real, dimension(4)                   :: corner_frequencies
        integer                              :: filter_type
        character                            :: misfit_type
        character                            :: model_parameter
        integer                              :: receiver_index
        integer                              :: src_index
        type(rec_param_type), pointer        :: receiver
        !pointer                              :: filter
    end type

end module
