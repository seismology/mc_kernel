module simple_routines
   use global_parameters, only: sp, dp

   implicit none
   private
   public mult2d_1d, mult3d_1d, absreldiff, check_NaN, cross
   public to_lower, lowtrim, firstderiv

   interface mult2d_1d
      module procedure  :: mult2d_1d_dble
      module procedure  :: mult2d_1d_cmplx
   end interface mult2d_1d

   interface mult3d_1d
      module procedure  :: mult3d_1d_dble
      module procedure  :: mult3d_1d_cmplx
   end interface mult3d_1d

   interface absreldiff
      module procedure  :: absreldiff_real
      module procedure  :: absreldiff_dble
   end interface absreldiff

   interface check_NaN
     module procedure   :: check_NaN_1d_sp
     module procedure   :: check_NaN_1d_dp
     module procedure   :: check_NaN_2d_sp
     module procedure   :: check_NaN_2d_dp
   end interface check_NaN

   interface cross
     module procedure   :: cross_sp
     module procedure   :: cross_dp
   end interface cross
  

   contains

!------------------------------------------------------------------------------
function mult2d_1d_dble(A, B) result(C)
! Multiply 1D array B along first dimension of 3D array A
   real(kind=dp), intent(in) :: A(:,:), B(:)
   real(kind=dp)             :: C(size(A,1), size(A,2))

   C = A * spread(B, 2, size(A,2))
end function mult2d_1d_dble
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
function mult2d_1d_cmplx(A, B) result(C)
! Multiply 1D array B along first dimension of 3D array A
   complex(kind=dp), intent(in) :: A(:,:), B(:)
   complex(kind=dp)             :: C(size(A,1), size(A,2))

   C = A * spread(B, 2, size(A,2))
end function mult2d_1d_cmplx
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
function mult3d_1d_dble(A, B) result(C)
! Multiply 1D array B along first dimension of 3D array A
   real(kind=dp), intent(in) :: A(:,:,:), B(:)
   real(kind=dp)             :: C(size(A,1), size(A,2), size(A,3))

   C = A * spread(spread(B, 2, size(A,2)), 3, size(A,3))
end function mult3d_1d_dble
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
function mult3d_1d_cmplx(A, B) result(C)
! Multiply 1D array B along first dimension of 3D array A
   complex(kind=dp), intent(in) :: A(:,:,:), B(:)
   complex(kind=dp)             :: C(size(A,1), size(A,2), size(A,3))

   C = A * spread(spread(B, 2, size(A,2)), 3, size(A,3))
end function mult3d_1d_cmplx
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
function to_lower(strIn) result(strOut)
!< Converts string to lowercase, adapted from http://www.star.le.ac.uk/~cgp/fortran.html
    implicit none

    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut
    integer                      :: i,j

    do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("A") .and. j<=iachar("Z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))+32)
        else
            strOut(i:i) = strIn(i:i)
        end if
    end do

end function to_lower
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
function lowtrim(strIn) result(strOut)
!< Converts string to lowercase and trims it
    implicit none

    character(len=*), intent(in)   :: strIn
    character(len=len_trim(strIn)) :: strOut

    strOut = trim(to_lower(strIn))

end function lowtrim 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
function firstderiv(timeseries) 
!< Calculates the first derivative of timeseries, using a compact stencil
real(kind=dp), intent(in) :: timeseries(:)
real(kind=dp)             :: firstderiv(size(timeseries,1))
integer                   :: ntimes

ntimes = size(timeseries, 1)
firstderiv       = 0
firstderiv(2:ntimes-1) = (timeseries(3:ntimes) - timeseries(1:ntimes-2)) / 2

end function firstderiv
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure real(kind=sp) function absreldiff_real(x1,x2)
!< Calculate the absolute relative difference between two numbers. Relative to
!! X1 or the non-zero one, that is. Single precision version
  real(kind=sp), intent(in) :: x1, x2
  real(kind=sp)             :: denominator

  denominator = max(abs(x1), abs(x2))

  if (denominator == 0) then
     absreldiff_real = 0
  else
     absreldiff_real = abs(x1-x2) / denominator
  end if

end function absreldiff_real
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure real(kind=dp) function absreldiff_dble(x1,x2)
!< Calculate the absolute relative difference between two numbers. Relative to
!! X1 or the non-zero one, that is. Double precision version
  real(kind=dp), intent(in) :: x1, x2
  real(kind=dp)             :: denominator

  denominator = max(abs(x1), abs(x2))

  if (denominator == 0) then
     absreldiff_dble = 0
  else
     absreldiff_dble = dabs(x1-x2) / denominator
  end if

end function absreldiff_dble
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure function cross_sp(a,b)
!< Compute the cross product between vectors a and b
  real(kind=sp)             :: cross_sp(3)
  real(kind=sp), intent(in) :: a(3), b(3)
           
  cross_sp(1) = a(2)*b(3) - a(3)*b(2)
  cross_sp(2) = a(3)*b(1) - a(1)*b(3)
  cross_sp(3) = a(1)*b(2) - a(2)*b(1)

end function cross_sp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure function cross_dp(a,b)
!< Compute the cross product between vectors a and b
  real(kind=dp)             :: cross_dp(3)
  real(kind=dp), intent(in) :: a(3), b(3)
           
  cross_dp(1) = a(2)*b(3) - a(3)*b(2)
  cross_dp(2) = a(3)*b(1) - a(1)*b(3)
  cross_dp(3) = a(1)*b(2) - a(2)*b(1)

end function cross_dp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Checks, whether any value of 'array' is NaN and returns its location
!> 1D single precision version
pure subroutine check_NaN_1d_sp(array, isnan, nan_loc)
  real(kind=sp), intent(in)     :: array(:)    !< 1D array to check
  logical,      intent(out)     :: isnan       !< Is any value of array NaN  
  integer,      intent(out)     :: nan_loc     !< Location of first NaN in array
  integer                       :: i1

  nan_loc = 0
  isnan = .false.

  if (any(array.ne.array)) then
    do i1 = 1, size(array,1)
      if (array(i1).ne.array(i1)) then
        nan_loc = i1
        isnan = .true.
      end if
    end do
  end if

end subroutine check_NaN_1d_sp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Checks, whether any value of 'array' is NaN and returns its location
!> 1D double precision version
pure subroutine check_NaN_1d_dp(array, isnan, nan_loc)
  real(kind=dp), intent(in)     :: array(:)    !< 1D array to check
  logical,      intent(out)     :: isnan       !< Is any value of array NaN  
  integer,      intent(out)     :: nan_loc     !< Location of first NaN in array
  integer                       :: i1

  nan_loc = 0
  isnan = .false.

  if (any(array.ne.array)) then
    do i1 = 1, size(array,1)
      if (array(i1).ne.array(i1)) then
        nan_loc = i1
        isnan = .true.
      end if
    end do
  end if

end subroutine check_NaN_1d_dp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Checks, whether any value of 'array' is NaN and returns its location
!> 2D single precision version
pure subroutine check_NaN_2d_sp(array, isnan, nan_loc)
  real(kind=sp), intent(in)     :: array(:,:)  !< 2D array to check
  logical,      intent(out)     :: isnan       !< Is any value of array NaN  
  integer,      intent(out)     :: nan_loc(2)  !< Location of first NaN in array
  integer                       :: i1, i2

  nan_loc = [0, 0]
  isnan = .false.

  if (any(array.ne.array)) then
    do i1 = 1, size(array,1)
      do i2 = 1, size(array,2)
        if (array(i1,i2).ne.array(i1,i2)) then
          nan_loc = [i1, i2]
          isnan = .true.
          exit
        end if
      end do
      if (isnan) exit
    end do
  end if

end subroutine check_NaN_2d_sp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Checks, whether any value of 'array' is NaN and returns its location
!> 2D double precision version
pure subroutine check_NaN_2d_dp(array, isnan, nan_loc)
  real(kind=dp), intent(in)     :: array(:,:)  !< 2D array to check
  logical,      intent(out)     :: isnan       !< Is any value of array NaN  
  integer,      intent(out)     :: nan_loc(2)  !< Location of first NaN in array
  integer                       :: i1, i2

  nan_loc = [0, 0]
  isnan = .false.

  if (any(array.ne.array)) then
    do i1 = 1, size(array,1)
      do i2 = 1, size(array,2)
        if (array(i1,i2).ne.array(i1,i2)) then
          nan_loc = [i1, i2]
          isnan = .true.
          exit
        end if
      end do
      if (isnan) exit
    end do
  end if

end subroutine check_NaN_2d_dp
!------------------------------------------------------------------------------

end module simple_routines
