module montecarlo
    
    use global_parameters
    implicit none

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
        logical                                  :: isinitialized = .false.
        
     contains
        procedure, pass                          :: check_montecarlo_integral
        procedure, pass                          :: initialize_montecarlo
        procedure, pass                          :: areallconverged
        procedure, pass                          :: getintegral
        procedure, pass                          :: getvariance
        procedure, pass                          :: isconverged
        procedure, pass                          :: freeme

    end type

    public                                       :: allallconverged
    public                                       :: allisconverged
contains

!-------------------------------------------------------------------------------
subroutine check_montecarlo_integral(this, func_values)
    class(integrated_type)                    :: this
    real(kind=dp), dimension(:,:), intent(in) :: func_values !(npoints, nfuncs)

    integer                                   :: npoints, ifuncs
    
    if (.not.this%isinitialized) then
       stop 'Initialize this MC type first'
    end if
    npoints = size(func_values, 1)
    this%nmodels = this%nmodels + npoints
    do ifuncs = 1, this%nfuncs
        if(this%converged(ifuncs)) cycle
        ! https://en.wikipedia.org/wiki/Monte_Carlo_Integration
        this%fsum(ifuncs)     = this%fsum(ifuncs)  + sum(func_values(:,ifuncs))
        this%f2sum(ifuncs)    = this%f2sum(ifuncs) + sum(func_values(:,ifuncs)**2)

        this%integral(ifuncs) = this%fsum(ifuncs) * &
                                this%volume / this%nmodels
        this%variance(ifuncs) = this%volume**2 / this%nmodels *    &
                                (  (this%f2sum(ifuncs) / this%nmodels)     &
                                 - (this%fsum(ifuncs)  / this%nmodels)**2  )

        if (this%variance(ifuncs) < this%allowed_variance(ifuncs)) then
           this%converged(ifuncs) = .true.
           ! Extra condition: Consider not converged, if error is more than twice 
           ! the integral value and the integral value is more than 1E-3 of allowed 
           ! error.
           if (this%variance(ifuncs)         > 2. * abs(this%integral(ifuncs))  .and. &
               this%variance(ifuncs) * 100.0 > this%allowed_variance(ifuncs) ) then
               this%converged(ifuncs) = .false.
           end if
        end if
    end do

end subroutine check_montecarlo_integral
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine initialize_montecarlo(this, nfuncs, volume, allowed_error)
    class(integrated_type), intent(inout) :: this 
    integer, intent(in)                   :: nfuncs
    real(kind=dp), intent(in)             :: volume
    real(kind=dp), intent(in)             :: allowed_error

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
    this%isinitialized    = .true.

end subroutine initialize_montecarlo
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine freeme(this)
    class(integrated_type), intent(inout) :: this 

    deallocate(this%fsum)
    deallocate(this%f2sum)
    deallocate(this%converged)
    deallocate(this%integral)
    deallocate(this%variance)
    deallocate(this%allowed_variance)
    this%isinitialized = .false.
    
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function getintegral(this)
    class(integrated_type)                :: this
    real(kind=dp), dimension(this%nfuncs) :: getintegral

    if (.not.this%isinitialized) then
       stop 'Initialize this MC type first'
    end if
    getintegral = this%integral

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function getvariance(this)
    class(integrated_type)                :: this
    real(kind=dp), dimension(this%nfuncs) :: getvariance

    if (.not.this%isinitialized) then
       stop 'Initialize this MC type first'
    end if
    getvariance = this%variance

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function areallconverged(this, ikernels)
    class(integrated_type)        :: this
    logical                       :: areallconverged
    integer, optional, intent(in) :: ikernels(:)

    if (.not.this%isinitialized) then
       stop 'Initialize this MC type first'
    end if

    if (present(ikernels)) then
       areallconverged = all(this%converged(ikernels))
    else 
       areallconverged = all(this%converged)
    end if

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function isconverged(this, ifunc)
    class(integrated_type)                :: this
    integer                               :: ifunc
    logical                               :: isconverged

    if (.not.this%isinitialized) then
       stop 'Initialize this MC type first'
    end if
    isconverged = this%converged(ifunc)

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function allallconverged(int_object, ikernels)
    type(integrated_type), intent(in)   :: int_object(:)
    logical                              :: allallconverged
    logical                              :: anynotconverged
    integer, optional, intent(in)        :: ikernels(:)
    integer                              :: ivertex, nvertices

    nvertices = size(int_object,1)
    anynotconverged = .false.

    do ivertex = 1, nvertices
        if (present(ikernels)) then
           anynotconverged = .not.(all(int_object(ivertex)%converged(ikernels))).or.anynotconverged
        else 
           anynotconverged = .not.(all(int_object(ivertex)%converged(:))).or.anynotconverged
        end if
    end do
    allallconverged = .not.anynotconverged

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function allisconverged(int_object, ikernel)
    type(integrated_type), intent(in)   :: int_object(:)
    logical                              :: allisconverged
    logical                              :: anynotconverged
    integer, optional, intent(in)        :: ikernel
    integer                              :: ivertex, nvertices

    nvertices = size(int_object)
    anynotconverged = .false.
    do ivertex = 1, nvertices
        anynotconverged = .not.(int_object(ivertex)%converged(ikernel)).or.anynotconverged
    end do
    allisconverged = .not.anynotconverged

end function
!-------------------------------------------------------------------------------
end module
