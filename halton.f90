!=========================================================================================
module halton_sequence
  use global_parameters
  implicit none
  private
  
  public :: init_halton
  public :: get_halton      
  public :: free_halton
                            

  integer, save, allocatable :: step(:)   !< stores current step of the sequence
  integer, save, allocatable :: seed(:)   !< initial seed (step that first sequence member 
                                          !! corresponds to)
  integer, save, allocatable :: base(:)   !< base of the halton sequence
  integer, save, allocatable :: leap(:)   !< how many sequence members to leap over
  integer, save              :: ndim

  integer, parameter         :: dim_max = 10
  integer, parameter         :: primes(dim_max) = [2, 3, 5, 7, 11, 13, 17, 19, 21, 23]

  logical, save              :: initialized = .false.
  
  
contains 

!-----------------------------------------------------------------------------------------
subroutine init_halton(ndim_in, seed_in, base_in, leap_in)
  !< Initializes a ndimensional Halton sequence

  integer, intent(in)            :: ndim_in           !< Dimensionality of Halton sequence
  integer, intent(in), optional  :: seed_in(ndim_in)  !< Seed for each dimension. Basically
                                                      !! the step with which to start the 
                                                      !! sequence
  integer, intent(in), optional  :: base_in(ndim_in)  !< Bases of the sequence. Has to be 
                                                      !! ndim (different) prime numbers
                                                      !! Preferrably small ones
  integer, intent(in), optional  :: leap_in(ndim_in)  !< Leap, i.e. number of sequence members
                                                      !! to omit after each call
  integer                        :: i_dim

  ndim = ndim_in

  allocate(seed(ndim))
  allocate(step(ndim))
  allocate(base(ndim))
  allocate(leap(ndim))


  if (present(seed_in)) then
     seed = seed_in
  else 
     seed = 20
  end if

  if (present(base_in)) then
     base = base_in
     do i_dim = 1, ndim
        if (.not.any(base(i_dim).eq.primes)) then
           write(*,'(A)') 'ERROR in Halton sequence initialization: '
           write(*,'(A,I3,A,I3,A)') 'Base number ', base(i_dim), '(dimension: ', i_dim, ') '
           write(*,'(A)') 'is not a prime number'
           stop
        end if
     end do
  else
     base(1:ndim) = primes(1:ndim)
  end if

  if (present(leap_in)) then
     leap = leap_in
  else 
     leap = 2
  end if
  
  step = 0

  initialized = .true.

end subroutine init_halton
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_halton(r)
!< Gets n values from a ndimensional Halton sequence. The dimensionality should be the 
!! same as in the initialization. If the sequence had not been initialized before, it is 
!! initialized here. If it had been initialized with a different dimensionality, an error 
!! occurs.
  implicit none

  real(kind=dp), intent(out) :: r(:,:) ! (ndim, npoints)
  real(kind=8)               :: digit(size(r,2))
  real(kind=8)               :: base_inv(size(r,2))
  integer                    :: seed2(size(r,2)), i_dim, j
  integer                    :: npoints

  npoints = size(r,2)

  if (.not.initialized) call init_halton(ndim_in = size(r,1))

  if (size(r,1).ne.ndim) then
     write(*,'(A,I3,A)') 'Halton sequence had been initialized for dimension ', ndim, ',' 
     write(*,'(A,I3,A)') 'but a series with number ', size(r,1), ' is requested.'
     stop
  end if 

  r(:,:) = 0

  do i_dim = 1, ndim 

    do j = 1, npoints
      seed2(j) = seed(i_dim) + ( step(i_dim) + j - 1 ) * leap(i_dim)
    end do

    base_inv = 1.0D+00 / real ( base(i_dim), kind = 8 )

    do while ( any ( seed2(:) /= 0 ) )
      digit(:) = mod ( seed2(:), base(i_dim) )
      r(i_dim,:) = r(i_dim,:) + real ( digit(:), kind = 8 ) * base_inv
      base_inv = base_inv / real ( base(i_dim), kind = 8 )
      seed2(:) = seed2(:) / base(i_dim)
    end do

  end do
  
  step = step + npoints

end subroutine get_halton
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine free_halton

  if (initialized) then
    ndim = -1
    deallocate(seed)
    deallocate(step)
    deallocate(base)
    deallocate(leap)
  end if

  initialized = .false.

end subroutine free_halton
!-----------------------------------------------------------------------------------------

end module halton_sequence
!=========================================================================================
