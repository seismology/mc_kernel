module montecarlo
    
    implicit none
    integer, parameter                           :: dp = 8, sp = 4

    type                                         :: integrated_type
        private
        integer                                  :: nfuncs
        real(kind=dp), dimension(:), allocatable :: fsum, f2sum
        real(kind=dp), dimension(:), allocatable :: integral
        real(kind=dp), dimension(:), allocatable :: allowed_error
        real(kind=dp), dimension(:), allocatable :: variance
        real(kind=dp)                            :: volume
        logical, dimension(:), allocatable       :: converged
        integer                                  :: nmodels
        contains
        procedure, pass                      :: check_montecarlo_integral
        procedure, pass                      :: initialize_montecarlo

    end type

contains

subroutine check_montecarlo_integral(this, func_values)
    class(integrated_type)                    :: this
    real(kind=sp), dimension(:,:), intent(in) :: func_values

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

        if (this%variance(ifuncs) < this%allowed_error(ifuncs)) then
            this%converged(ifuncs) = .true.
        end if
    end do

end subroutine check_montecarlo_integral

subroutine initialize_montecarlo(this, nfuncs, volume)
    class(integrated_type), intent(inout) :: this 
    integer, intent(in)                   :: nfuncs
    real(kind=dp), intent(in)             :: volume

    deallocate(this%fsum)
    deallocate(this%f2sum)
    deallocate(this%converged)
    deallocate(this%variance)
    deallocate(this%allowed_error)
    
    allocate(this%fsum(nfuncs))
    allocate(this%f2sum(nfuncs))
    allocate(this%converged(nfuncs))
    allocate(this%variance(nfuncs))
    allocate(this%allowed_error(nfuncs))

    this%volume        = volume
    this%nfuncs        = nfuncs
    this%fsum          = 0
    this%f2sum         = 0
    this%allowed_error = 1d-10
    this%converged     = .false.
    this%nmodels       = 0

end subroutine initialize_montecarlo

function areallconverged(this)
    class(integrated_type)  :: this
    logical                 :: areallconverged

    areallconverged = all(this%converged)
end function

end module
