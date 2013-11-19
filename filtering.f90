module filtering
   use global_parameters
   implicit none
   type filter_type

       character(len=32)               :: name
       complex(kind=dp), allocatable   :: transferfunction(:)
       integer                         :: nfreq                !< Number of frequencies
       real(kind=dp), allocatable      :: f(:)
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
!> Create a filter object with the specified parameters
subroutine create(this, name, dfreq, nfreq, filterclass, frequencies)
    class(filter_type)              :: this
    integer, intent(in)             :: nfreq
    real(kind=dp), intent(in)       :: dfreq, frequencies(4)
    integer                         :: ifreq
    character(len=32), intent(in)   :: name
    character(len=32), intent(in)   :: filterclass
    character(len=32)               :: fmtstring
    character(len=64)               :: fnam

    fmtstring = '("Class: ", A, ", freq: ", 4(F8.4))'

    if (this%initialized) then
       write(*,*) 'ERROR: This filter is already initialized as'
       write(*, fmtstring) trim(this%filterclass), this%frequencies
       write(*,*) 'delete it first before using'
       stop
    end if
    allocate(this%transferfunction(nfreq))
    allocate(this%f(nfreq))

    this%name              = name
    this%transferfunction  = 0.0
    this%nfreq             = nfreq
    this%filterclass       = filterclass
    this%frequencies       = frequencies

    do ifreq = 1, nfreq
       this%f(ifreq) = dfreq * (ifreq-1)
    end do

    select case(trim(filterclass))
    case('Gabor')
       this%transferfunction = loggabor(this%f, frequencies(1), frequencies(2))
    case default
       print *, 'ERROR: Unknown filter type: ', trim(filterclass)
       stop
    end select

20  format('filterresponse_', A, '_', 2(F08.3))
    write(fnam,20) trim(filterclass), frequencies(1:2)
    open(10, file=trim(fnam), action='write')
    do ifreq = 1, nfreq
       write(10,*), this%f(ifreq), real(this%transferfunction(ifreq)), imag(this%transferfunction(ifreq))
    end do
    close(10)
    this%initialized = .true.
end subroutine

! -----------------------------------------------------------------------------
!> Delete this filter and free the memory
subroutine deleteme(this)
    class(filter_type)              :: this

    this%filterclass = ''
    deallocate(this%transferfunction)
    deallocate(this%f)
    this%initialized = .false.
end subroutine

! -----------------------------------------------------------------------------
!> Apply this filter to one trace (in the frequency domain)
function apply_1d(this, freq_series)
   class(filter_type)               :: this
   complex(kind=dp), intent(in)     :: freq_series(:)
   complex(kind=dp)                 :: apply_1d(size(freq_series))
   if (.not.this%initialized) then
      write(*,*) 'ERROR: Filter is not initialized yet'
      stop
   end if
   if (size(freq_series, 1).ne.this%nfreq) then
      write(*,*) 'ERROR: Filter length: ', this%nfreq, ', data length: ', size(freq_series, 1)
      stop
   end if
   apply_1d = freq_series * this%transferfunction
end function apply_1d 

! -----------------------------------------------------------------------------
!> Apply this filter to multiple traces (in the frequency domain)
function apply_2d(this, freq_series)
   class(filter_type)               :: this
   complex(kind=dp), intent(in)     :: freq_series(:,:)
   complex(kind=dp)                 :: apply_2d(size(freq_series,1), size(freq_series,2))
   integer                          :: itrace

   if (.not.this%initialized) then
      write(*,*) 'ERROR: Filter is not initialized yet'
      stop
   end if
   if (size(freq_series, 1).ne.this%nfreq) then
      write(*,*) 'ERROR: Filter length: ', this%nfreq, ', data length: ', size(freq_series, 1)
      stop
   end if
   do itrace = 1, (size(freq_series, 2))
      apply_2d(:,itrace) = freq_series(:,itrace) * this%transferfunction(:)
   end do

end function apply_2d

! -----------------------------------------------------------------------------
!> Returns the transfer function of this filter
function get_transferfunction(this)
    class(filter_type)           :: this
    real(kind=sp)                :: get_transferfunction(this%nfreq, 3)

    if (.not.this%initialized) then
       write(*,*) 'ERROR: Filter is not initialized yet'
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

! ------------------------------------------------------------------------------
!> Apply a timeshift of dtshift on the frequency domain traces in field
subroutine timeshift(field, freq, dtshift)

   complex(kind=dp), intent(inout)  :: field(:,:) !< Frequency domain traces to apply 
                                                  !! time shift on. 
                                                  !! Dimension: nfreq x ntraces
   real(kind=dp),    intent(in)     :: freq(:)    !< 1D array with the frequencies of 
                                                  !! the entries in field
                                                  !! Dimension: nfreq
   real(kind=sp),    intent(in)     :: dtshift    !< Time shift to apply (in seconds)
   
   integer                          :: ipoint
   complex(kind=dp), allocatable    :: shift_fd(:)

   allocate( shift_fd(size(field,1)) )
   shift_fd = exp( 2*pi * cmplx(0., 1.) * freq(:) * dtshift)

   do ipoint = 1, size(field,2)
      field(:, ipoint) = field(:, ipoint) * shift_fd
   end do

end subroutine timeshift

!-----------------------------------------------------------------------------!
!                 BLOCK WITH SEPARATE FILTER ROUTINES                         !
!-----------------------------------------------------------------------------!
function loggabor(f, p_c, sigma)             !< Log-Gabor filter as employed by K. Sigloch
    real(kind=dp), intent(in)       :: f(:)  !< Frequency array
    real(kind=dp), intent(in)       :: p_c   !< Center period of the filter
    real(kind=dp), intent(in)       :: sigma !< (fixed) ratio of sigma divided by p_c, 
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
