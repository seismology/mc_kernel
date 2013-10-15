module type_parameter
    use global_parameters,               only : sp, dp, pi, deg2rad
    use source, only: src_param_type
    use receiver, only: rec_param_type
    implicit none    

    type parameter_type
        type(src_param_type)                 :: source
        type(rec_param_type), allocatable    :: receiver(:)

        real(kind=sp)                        :: allowed_error
        character(len=512)                   :: dir_fwdmesh
        character(len=512)                   :: dir_bwdmesh
        integer                              :: nsim_fwd, nsim_bwd
    end type

!    type kernelspec_type
!        real, dimension(2)                   :: time_window
!        real, dimension(4)                   :: corner_frequencies
!        integer                              :: filter_type
!        character                            :: misfit_type
!        character                            :: model_parameter
!        integer                              :: receiver_index
!        integer                              :: src_index
!        type(rec_param_type), pointer        :: receiver
!        !pointer                              :: filter
!    end type

end module
!------------------------------------------------------------------------------
