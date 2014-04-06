module simple_routines
   use global_parameters, only: dp

   implicit none

   interface mult2d_1d
      module procedure  :: mult2d_1d_dble
      module procedure  :: mult2d_1d_cmplx
   end interface mult2d_1d

   interface mult3d_1d
      module procedure  :: mult3d_1d_dble
      module procedure  :: mult3d_1d_cmplx
   end interface mult3d_1d

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
function firstderiv(timeseries) 
!< Calculates the first derivative of timeseries, using a compact stencil
real(kind=dp), intent(in) :: timeseries(:)
real(kind=dp)             :: firstderiv(size(timeseries,1))
integer                   :: ntimes

ntimes = size(timeseries, 1)
firstderiv       = 0
firstderiv(2:99) = (timeseries(3:100) - timeseries(1:98)) / 2

end function firstderiv
!------------------------------------------------------------------------------

end module simple_routines
