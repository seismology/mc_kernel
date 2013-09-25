module montecarlo
    
    implicit none
    integer, parameter                           :: dp = 8, sp = 4

    type                                         :: integrated_type
        private
        integer                                  :: nfuncs
        real(kind=dp), dimension(:), allocatable :: fsum, f2sum
        real(kind=dp), dimension(:), allocatable :: integral
        real(kind=dp), dimension(:), allocatable :: allowed_variance
        real(kind=dp), dimension(:), allocatable :: variance
        real(kind=dp)                            :: volume
        logical, dimension(:), allocatable       :: converged
        integer                                  :: nmodels
        contains
        procedure, pass                      :: check_montecarlo_integral
        procedure, pass                      :: initialize_montecarlo
        procedure, pass                      :: areallconverged
        procedure, pass                      :: getintegral, getvariance

    end type

contains

subroutine check_montecarlo_integral(this, func_values)
    class(integrated_type)                    :: this
    real(kind=dp), dimension(:,:), intent(in) :: func_values !(npoints, nfuncs)

    integer                                   :: npoints, ifuncs

    
    npoints = size(func_values, 1)
    this%nmodels = this%nmodels + npoints
    do ifuncs = 1, this%nfuncs
        if(this%converged(ifuncs)) cycle
        ! https://en.wikipedia.org/wiki/Monte_Carlo_Integration
        this%fsum(ifuncs)     = this%fsum(ifuncs)  + sum(func_values(:,ifuncs))
        this%f2sum(ifuncs)    = this%f2sum(ifuncs) + sum(func_values(:,ifuncs))**2

        this%integral(ifuncs) = this%fsum(ifuncs) * &
                                this%volume / this%nmodels
        this%variance(ifuncs) = this%volume**2 / this%nmodels *    &
                                (  (this%f2sum(ifuncs) / this%nmodels)     &
                                 - (this%fsum(ifuncs)  / this%nmodels)**2  )

        if (this%variance(ifuncs) < this%allowed_variance(ifuncs)) then
            this%converged(ifuncs) = .true.
        end if
    end do

end subroutine check_montecarlo_integral

subroutine initialize_montecarlo(this, nfuncs, volume, allowed_error)
    class(integrated_type), intent(inout) :: this 
    integer, intent(in)                   :: nfuncs
    real(kind=dp), intent(in)             :: volume, allowed_error

    if(allocated(this%fsum)) deallocate(this%fsum)
    if(allocated(this%f2sum)) deallocate(this%f2sum)
    if(allocated(this%converged)) deallocate(this%converged)
    if(allocated(this%integral)) deallocate(this%integral)
    if(allocated(this%variance)) deallocate(this%variance)
    if(allocated(this%allowed_variance)) deallocate(this%allowed_variance)
    
    allocate(this%fsum(nfuncs))
    allocate(this%f2sum(nfuncs))
    allocate(this%converged(nfuncs))
    allocate(this%integral(nfuncs))
    allocate(this%variance(nfuncs))
    allocate(this%allowed_variance(nfuncs))

    this%volume           = volume
    this%nfuncs           = nfuncs
    this%fsum             = 0
    this%f2sum            = 0
    this%allowed_variance = allowed_error ** 2
    this%converged        = .false.
    this%nmodels          = 0

end subroutine initialize_montecarlo

function areallconverged(this)
    class(integrated_type)  :: this
    logical                 :: areallconverged

    areallconverged = all(this%converged)
end function

function getintegral(this)
    class(integrated_type)                :: this
    real(kind=dp), dimension(this%nfuncs) :: getintegral

    getintegral = this%integral

end function

function isconverged(this, ifunc)
    class(integrated_type)                :: this
    integer                               :: ifunc
    logical                               :: isconverged

    isconverged = this%converged(ifunc)
end function

function getvariance(this)
    class(integrated_type)                :: this
    real(kind=dp), dimension(this%nfuncs) :: getvariance

    getvariance = this%variance

end function
end module
