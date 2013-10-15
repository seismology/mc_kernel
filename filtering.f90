module filtering
   use global_parameters
   implicit none
   type filter_type

       complex(kind=dp), allocatable   :: transferfunction(:)
       integer                         :: nomega               !< Number of frequencies
       real(kind=sp), allocatable      :: f(:)
       logical                         :: initialized = .false.
       character(len=32)               :: filterclass
       real(kind=sp)                   :: frequencies(4)


   contains

       procedure, pass   :: create
       procedure, pass   :: deleteme
       procedure, pass   :: apply_1d
       procedure, pass   :: apply_2d
       procedure, pass   :: isinitialized
       procedure, pass   :: get_transferfunction

   end type


contains

! -----------------------------------------------------------------------------
subroutine create(this, domega, nomega, filterclass, frequencies)
    class(filter_type)              :: this
    integer, intent(in)             :: nomega
    real(kind=sp), intent(in)       :: domega, frequencies(4)
    integer                         :: iomega
    character(len=32), intent(in)   :: filterclass
    character(len=32)               :: fmtstring

    fmtstring = '("Class: ", A, ", freq: ", 4(F8.4))'

    if (this%initialized) then
       write(*,*) 'This filter is already initialized as'
       write(*, fmtstring) trim(this%filterclass), this%frequencies
       write(*,*) 'delete it first before using'
       stop
    end if
    allocate(this%transferfunction(nomega))
    allocate(this%f(nomega))

    this%transferfunction  = 0.0
    this%nomega            = nomega
    this%filterclass       = filterclass
    this%frequencies       = frequencies

    do iomega = 1, nomega
       this%f(iomega) = domega * (iomega-1)
    end do

    select case(trim(filterclass))
    case('Gabor')
       this%transferfunction = loggabor(this%f, frequencies(1), frequencies(2))
    case default
       stop 'Unknown filter type'
    end select

    open(10, file='filterresponse.txt', action='write')
    do iomega = 1, nomega
       write(10,*), this%f(iomega), real(this%transferfunction(iomega)), imag(this%transferfunction(iomega))
    end do
    close(10)
    this%initialized = .true.
end subroutine

! -----------------------------------------------------------------------------
subroutine deleteme(this)
    class(filter_type)              :: this

    this%filterclass = ''
    deallocate(this%transferfunction)
    deallocate(this%f)
    this%initialized = .false.
end subroutine

! -----------------------------------------------------------------------------
function apply_1d(this, freq_series)
   class(filter_type)               :: this
   complex(kind=dp), intent(in)     :: freq_series(:)
   complex(kind=dp)                 :: apply_1d(size(freq_series))
   if (.not.this%initialized) then
      write(*,*) 'Filter is not initialized yet'
      stop
   end if
   if (size(freq_series, 1).ne.this%nomega) then
      write(*,*) 'Filter length: ', this%nomega, ', data length: ', size(freq_series, 1)
      stop
   end if
   apply_1d = freq_series * this%transferfunction
end function apply_1d 

! -----------------------------------------------------------------------------
function apply_2d(this, freq_series)
   class(filter_type)               :: this
   complex(kind=dp), intent(in)     :: freq_series(:,:)
   complex(kind=dp)                 :: apply_2d(size(freq_series,1), size(freq_series,2))
   integer                          :: itrace

   if (.not.this%initialized) then
      write(*,*) 'Filter is not initialized yet'
      stop
   end if
   if (size(freq_series, 1).ne.this%nomega) then
      write(*,*) 'Filter length: ', this%nomega, ', data length: ', size(freq_series, 1)
      stop
   end if
   do itrace = 1, (size(freq_series, 2))
      apply_2d(:,itrace) = freq_series(:,itrace) * this%transferfunction(:)
   end do

end function apply_2d

! -----------------------------------------------------------------------------
function get_transferfunction(this)
    class(filter_type)           :: this
    real(kind=sp)                :: get_transferfunction(this%nomega, 3)

    if (.not.this%initialized) then
       write(*,*) 'Filter is not initialized yet'
       stop
    end if
    get_transferfunction(:,1) = this%f(:)
    get_transferfunction(:,2) = real(this%transferfunction(:))
    get_transferfunction(:,3) = imag(this%transferfunction(:))

end function

! -----------------------------------------------------------------------------
function isinitialized(this)
   class(filter_type)            :: this
   logical                       :: isinitialized

   isinitialized = this%initialized
end function



!-----------------------------------------------------------------------------!
!                 BLOCK WITH SEPARATE FILTER ROUTINES                         !
!-----------------------------------------------------------------------------!
function loggabor(f, p_c, sigma)             !< Log-Gabor filter as employed by K. Sigloch
    real(kind=sp), intent(in)       :: f(:)
    real(kind=sp), intent(in)       :: p_c   !< Center period of the filter
    real(kind=sp), intent(in)       :: sigma !< (fixed) ratio of sigma divided by p_c, 
                                             !! where sigma is the standard
                                             !! deviation of the Gaussian that describes the
                                             !! log-Gabor filter in the *time* domain, 
                                             !! and fc is the filter's center frequency.
                                             !! Hence larger sigmaIfc means narrower bandwidth
                                             !! sigmaIfc = .5 means the one-sigma interval i
                                             !! equals one octave

    real(kind=sp)                   :: loggabor(size(f))
    real(kind=sp)                   :: f_c

    f_c = 1 / p_c
    loggabor = exp( -(log(f/f_c))**2 / ( 2 * log(sigma)**2) )
    ! G(f,k) = exp( -(ln(f/f_k))^2 / (2*ln(sigmaIfc)^2)  )

end function loggabor
end module
