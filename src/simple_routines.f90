!******************************************************************************
!
!    This file is part of:
!    MC Kernel: Calculating seismic sensitivity kernels on unstructured meshes
!    Copyright (C) 2016 Simon Staehler, Martin van Driel, Ludwig Auer
!
!    You can find the latest version of the software at:
!    <https://www.github.com/tomography/mckernel>
!
!    MC Kernel is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MC Kernel is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MC Kernel. If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

!=========================================================================================
module simple_routines
   use global_parameters, only: sp, dp, qp, verbose

   implicit none
   private
   public mult2d_1d, mult3d_1d, absreldiff, check_NaN, cross
   public to_lower, lowtrim, check_limits
   public cumsum_trapezoidal

   interface mult2d_1d
      module procedure  :: mult2d_1d_dble
      module procedure  :: mult2d_1d_cmplx
   end interface mult2d_1d

   interface mult3d_1d
      module procedure  :: mult3d_1d_dble
      module procedure  :: mult3d_1d_cmplx
   end interface mult3d_1d

   interface absreldiff
      module procedure  :: absreldiff_sp
      module procedure  :: absreldiff_dp
      module procedure  :: absreldiff_qp
   end interface absreldiff

   interface cross
     module procedure   :: cross_sp
     module procedure   :: cross_dp
     module procedure   :: cross_qp
   end interface cross

   interface cumsum_trapezoidal
     module procedure   :: cumsum_trapezoidal_1d
     module procedure   :: cumsum_trapezoidal_2d
     module procedure   :: cumsum_trapezoidal_3d
   end interface cumsum_trapezoidal

   interface check_NaN
     module procedure   :: check_NaN_1d_sp
     module procedure   :: check_NaN_1d_dp
     module procedure   :: check_NaN_2d_sp
     module procedure   :: check_NaN_2d_dp
   end interface check_NaN

   interface check_limits
     module procedure   :: checklim_1d
     module procedure   :: checklim_2d
     module procedure   :: checklim_3d
     module procedure   :: checklim_1d_int
     module procedure   :: checklim_2d_int
     module procedure   :: checklim_3d_int
   end interface check_limits

   contains

!------------------------------------------------------------------------------
pure function mult2d_1d_dble(A, B) result(C)
! Multiply 1D array B along first dimension of 3D array A
   real(kind=dp), intent(in) :: A(:,:), B(:)
   real(kind=dp)             :: C(size(A,1), size(A,2))

   C = A * spread(B, 2, size(A,2))
end function mult2d_1d_dble
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure function mult2d_1d_cmplx(A, B) result(C)
! Multiply 1D array B along first dimension of 3D array A
   complex(kind=dp), intent(in) :: A(:,:), B(:)
   complex(kind=dp)             :: C(size(A,1), size(A,2))

   C = A * spread(B, 2, size(A,2))
end function mult2d_1d_cmplx
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
pure function mult3d_1d_dble(A, B) result(C)
! Multiply 1D array B along first dimension of 3D array A
   real(kind=dp), intent(in) :: A(:,:,:), B(:)
   real(kind=dp)             :: C(size(A,1), size(A,2), size(A,3))

   C = A * spread(spread(B, 2, size(A,2)), 3, size(A,3))
end function mult3d_1d_dble
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure function mult3d_1d_cmplx(A, B) result(C)
! Multiply 1D array B along first dimension of 3D array A
   complex(kind=dp), intent(in) :: A(:,:,:), B(:)
   complex(kind=dp)             :: C(size(A,1), size(A,2), size(A,3))

   C = A * spread(spread(B, 2, size(A,2)), 3, size(A,3))
end function mult3d_1d_cmplx
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure function to_lower(strIn) result(strOut)
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
pure function lowtrim(strIn) result(strOut)
!< Converts string to lowercase and trims it
    implicit none

    character(len=*), intent(in)   :: strIn
    character(len=len_trim(strIn)) :: strOut

    strOut = trim(to_lower(strIn))

end function lowtrim 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure real(kind=sp) function absreldiff_sp(x1,x2)
!< Calculate the absolute relative difference between two numbers. Relative to
!! X1 or the non-zero one, that is. Single precision version
  real(kind=sp), intent(in) :: x1, x2
  real(kind=sp)             :: denominator

  denominator = max(abs(x1), abs(x2))

  if (denominator == 0) then
     absreldiff_sp = 0
  else
     absreldiff_sp = abs(x1-x2) / denominator
  end if

end function absreldiff_sp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure real(kind=dp) function absreldiff_dp(x1,x2)
!< Calculate the absolute relative difference between two numbers. Relative to
!! X1 or the non-zero one, that is. Double precision version
  real(kind=dp), intent(in) :: x1, x2
  real(kind=dp)             :: denominator

  denominator = max(abs(x1), abs(x2))

  if (denominator == 0) then
     absreldiff_dp = 0
  else
     absreldiff_dp = dabs(x1-x2) / denominator
  end if

end function absreldiff_dp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure real(kind=qp) function absreldiff_qp(x1,x2)
!< Calculate the absolute relative difference between two numbers. Relative to
!! X1 or the non-zero one, that is. Quadruple precision version
  real(kind=qp), intent(in) :: x1, x2
  real(kind=qp)             :: denominator

  denominator = max(abs(x1), abs(x2))

  if (denominator == 0) then
     absreldiff_qp = 0
  else
     absreldiff_qp = abs(x1-x2) / denominator
  end if

end function absreldiff_qp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure function cross_sp(a,b)
!< Compute the cross product between vectors a and b
  real(kind=sp)             :: cross_sp(3)
  real(kind=sp), intent(in) :: a(3), b(3)
  real(kind=qp)             :: a_qp(3), b_qp(3)

  a_qp = a
  b_qp = b
           
  cross_sp(1) = real(a_qp(2)*b_qp(3) - a_qp(3)*b_qp(2), kind=sp)
  cross_sp(2) = real(a_qp(3)*b_qp(1) - a_qp(1)*b_qp(3), kind=sp)
  cross_sp(3) = real(a_qp(1)*b_qp(2) - a_qp(2)*b_qp(1), kind=sp)

end function cross_sp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure function cross_dp(a,b)
!< Compute the cross product between vectors a and b
  real(kind=dp)             :: cross_dp(3)
  real(kind=dp), intent(in) :: a(3), b(3)
  real(kind=qp)             :: a_qp(3), b_qp(3)

  a_qp = a
  b_qp = b
           
  cross_dp(1) = real(a_qp(2)*b_qp(3) - a_qp(3)*b_qp(2), kind=dp)
  cross_dp(2) = real(a_qp(3)*b_qp(1) - a_qp(1)*b_qp(3), kind=dp)
  cross_dp(3) = real(a_qp(1)*b_qp(2) - a_qp(2)*b_qp(1), kind=dp)
                                                
end function cross_dp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
pure function cross_qp(a,b)
!< Compute the cross product between vectors a and b
  real(kind=qp)             :: cross_qp(3)
  real(kind=qp), intent(in) :: a(3), b(3)

  cross_qp(1) = a(2)*b(3) - a(3)*b(2)
  cross_qp(2) = a(3)*b(1) - a(1)*b(3)
  cross_qp(3) = a(1)*b(2) - a(2)*b(1)
                                                
end function cross_qp
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Integrate timeseries once using the trapezoidal rule. 1D version
pure function cumsum_trapezoidal_1d(timeseries, dt) result(cumsum)
    real(kind=dp), intent(in)   :: timeseries(:), dt
    real(kind=dp)               :: cumsum(size(timeseries))
    integer                     :: npoints, ipoint

    npoints = size(timeseries, 1)

    ! Trapezoidal rule: I = dt/2 * (f(x1) + 2f(x2) + ... + 2f(xN-1) + f(xN))
    cumsum(1) = timeseries(1)
    if (npoints>1) then
      cumsum(2) = cumsum(1) + timeseries(2)
      if (npoints>2) then
        do ipoint = 3, npoints
          cumsum(ipoint) = cumsum(ipoint-1) + timeseries(ipoint-1) + timeseries(ipoint) 
          !timeseries(1) + timeseries(ipoint) + 2*sum(timeseries(2:ipoint-1), 1)
        end do
      end if
    end if
    cumsum = cumsum * dt * 0.5

end function cumsum_trapezoidal_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Integrate timeseries once using the trapezoidal rule. 2D version
!! Integration is done along first dimension.
!! @TODO: This is inefficient memory access and should be optimized
pure function cumsum_trapezoidal_2d(timeseries, dt) result(cumsum)
    real(kind=dp), intent(in)   :: timeseries(:,:), dt
    real(kind=dp)               :: cumsum(size(timeseries, 1), &
                                          size(timeseries, 2))
    integer                     :: npoints, ipoint

    npoints = size(timeseries, 1)

    ! Trapezoidal rule: I = dt/2 * (f(x1) + 2f(x2) + ... + 2f(xN-1) + f(xN))
    cumsum(1,:) = timeseries(1,:)
    if (npoints>1) then
      cumsum(2,:) = cumsum(1,:) + timeseries(2,:)
      if (npoints>2) then
        do ipoint = 3, npoints
          cumsum(ipoint,:) = cumsum(ipoint-1,:) &
                           + timeseries(ipoint-1,:) + timeseries(ipoint,:) 
          !timeseries(1) + timeseries(ipoint) + 2*sum(timeseries(2:ipoint-1), 1)
        end do
      end if
    end if
    cumsum = cumsum * dt * 0.5

end function cumsum_trapezoidal_2d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Integrate timeseries once using the trapezoidal rule. 3D version
!! Integration is done along first dimension.
!! @TODO: This is inefficient memory access and should be optimized
pure function cumsum_trapezoidal_3d(timeseries, dt) result(cumsum)
    real(kind=dp), intent(in)   :: timeseries(:,:,:), dt
    real(kind=dp)               :: cumsum(size(timeseries, 1), &
                                          size(timeseries, 2), &
                                          size(timeseries, 3))
    integer                     :: npoints, ipoint

    npoints = size(timeseries, 1)

    ! Trapezoidal rule: I = dt/2 * (f(x1) + 2f(x2) + ... + 2f(xN-1) + f(xN))
    cumsum(1,:,:) = timeseries(1,:,:)
    if (npoints>1) then
      cumsum(2,:,:) = cumsum(1,:,:) + timeseries(2,:,:)
      if (npoints>2) then
        do ipoint = 3, npoints
          cumsum(ipoint,:,:) = cumsum(ipoint-1,:,:) & 
                             + timeseries(ipoint-1,:,:) + timeseries(ipoint,:,:) 
          !timeseries(1) + timeseries(ipoint) + 2*sum(timeseries(2:ipoint-1), 1)
        end do
      end if
    end if
    cumsum = cumsum * dt * 0.5

end function cumsum_trapezoidal_3d
!-------------------------------------------------------------------------------


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

!-----------------------------------------------------------------------------------------
!> Checks whether any value of 'array' is outside of 'limits' or NaN
!! This function does not use polymorphism, since it is buggy and generally stinks
function checklim_1d_int(array, limits, array_name, ntoosmall, ntoolarge) result(out_of_limit)
  integer, intent(in)                    ::  array(:)
  integer, intent(in), optional          ::  limits(2)
  character(len=*), intent(in), optional ::  array_name
  integer, intent(out), optional         ::  ntoosmall, ntoolarge

  logical                                ::  out_of_limit
  integer                                ::  limits_loc(2)
  integer                                ::  ntoosmall_loc, ntoolarge_loc
  logical                                ::  toosmall(size(array)), toolarge(size(array))
  integer                                ::  i

  out_of_limit = .false.

  if (present(limits)) then
    limits_loc(1) = limits(1)
    limits_loc(2) = limits(2)
  else
    limits_loc(1) = -huge(array)
    limits_loc(2) = huge(array)
  end if

  toosmall = (.not.array.ge.limits_loc(1)) ! This catches NaNs as well, which give .false.
                                       ! for every binary comparison.

  toolarge = (.not.array.le.limits_loc(2)) ! This catches NaNs as well, which give .false.
                                       ! for every binary comparison.

  ntoosmall_loc = count(toosmall)
  if (ntoosmall_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array smaller than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Lower limit  :              ', minval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoosmall_loc    
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do i = 1, size(array, 1)
          if (.not.array(i).ge.(limits_loc(1))) then
              write(*,*) i, array(i)
          end if
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  ntoolarge_loc = count(toolarge)
  if (ntoolarge_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array larger than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Upper limit  :              ', maxval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoolarge_loc
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do i = 1, size(array, 1)
          if (.not.array(i).le.(limits_loc(2))) then
            write(*,*) i, array(i)
          end if
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  if (present(ntoosmall)) ntoosmall = ntoosmall_loc
  if (present(ntoolarge)) ntoolarge = ntoolarge_loc

end function checklim_1d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Checks whether any value of 'array' is outside of 'limits_loc' or NaN (2D version)
function checklim_2d_int(array, limits, array_name, ntoosmall, ntoolarge) result(out_of_limit)
  integer, intent(in)                    ::  array(:,:)
  integer, intent(in), optional          ::  limits(2)
  character(len=*), intent(in), optional ::  array_name
  integer, intent(out), optional         ::  ntoosmall, ntoolarge

  logical                                ::  out_of_limit
  integer                                ::  limits_loc(2)
  integer                                ::  ntoosmall_loc, ntoolarge_loc
  logical                                ::  toosmall(size(array,1), size(array,2))
  logical                                ::  toolarge(size(array,1), size(array,2))
  integer                                ::  i, j

  out_of_limit = .false.

  if (present(limits)) then
    limits_loc(1) = limits(1)
    limits_loc(2) = limits(2)
  else
    limits_loc(1) = -huge(array)
    limits_loc(2) = huge(array)
  end if

  toosmall = (.not.array.ge.limits_loc(1)) ! This catches NaNs as well, which give .false.
                                       ! for every binary comparison.

  toolarge = (.not.array.le.limits_loc(2)) ! This catches NaNs as well, which give .false.
                                       ! for every binary comparison.

  ntoosmall_loc = count(toosmall)
  if (ntoosmall_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array smaller than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Lower limit  :              ', minval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoosmall_loc
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do j = 1, size(array, 2)
          do i = 1, size(array, 1)
            if (.not.array(i,j).ge.(limits_loc(1))) then
              write(*,*) i, j, array(i,j)
            end if
          end do
        end do
      end if 
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  ntoolarge_loc = count(toolarge)
  if (ntoolarge_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array larger than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Upper limit  :              ', maxval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoolarge_loc
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do j = 1, size(array, 2)
          do i = 1, size(array, 1)
            if (.not.array(i,j).le.(limits_loc(2))) then
              write(*,*) i, j, array(i,j)
            end if
          end do
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  if (present(ntoosmall)) ntoosmall = ntoosmall_loc
  if (present(ntoolarge)) ntoolarge = ntoolarge_loc

end function checklim_2d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Checks whether any value of 'array' is outside of 'limits_loc' or NaN
function checklim_3d_int(array, limits, array_name, ntoosmall, ntoolarge) result(out_of_limit)
  integer, intent(in)                    ::  array(:,:,:)
  integer, intent(in), optional          ::  limits(2)
  character(len=*), intent(in), optional ::  array_name
  integer, intent(out), optional         ::  ntoosmall, ntoolarge

  logical                                ::  out_of_limit
  integer                                ::  limits_loc(2)
  integer                                ::  ntoosmall_loc, ntoolarge_loc
  logical, allocatable                   ::  toosmall(:,:,:), toolarge(:,:,:)
  integer                                ::  i, j, k

  out_of_limit = .false.

  if (present(limits)) then
    limits_loc(1) = limits(1)
    limits_loc(2) = limits(2)
  else
    limits_loc(1) = -huge(array)
    limits_loc(2) = huge(array)
  end if

  allocate(toosmall(size(array,1), size(array,2), size(array,3)))
  allocate(toolarge(size(array,1), size(array,2), size(array,3)))
  toosmall = (.not.array.ge.limits_loc(1)) ! This catches NaNs as well, which give .false.
                                           ! for every binary comparison.

  toolarge = (.not.array.le.limits_loc(2)) ! This catches NaNs as well, which give .false.
                                           ! for every binary comparison.

  ntoosmall_loc = count(toosmall)
  if (ntoosmall_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array smaller than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Lower limit  :              ', minval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoosmall_loc
      write(*,*) 'Details for first violating element: follow'
      if (verbose>1) then
        do k = 1, size(array, 3)
          do j = 1, size(array, 2)
            do i = 1, size(array, 1)
              if (.not.array(i,j,k).ge.(limits_loc(1))) then
                write(*,*) i, j, k, array(i,j,k)
              end if
            end do
          end do
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  ntoolarge_loc = count(toolarge)
  if (ntoolarge_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array larger than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Upper limit  :              ', maxval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoolarge_loc
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do k = 1, size(array, 3)
          do j = 1, size(array, 2)
            do i = 1, size(array, 1)
              if (.not.array(i,j,k).le.(limits_loc(2))) then
                write(*,*) i, j, k, array(i,j,k)
              end if
            end do
          end do
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  if (present(ntoosmall)) ntoosmall = ntoosmall_loc
  if (present(ntoolarge)) ntoolarge = ntoolarge_loc

end function checklim_3d_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Checks whether any value of 'array' is outside of 'limits' or NaN
function checklim_1d(array, limits, array_name, ntoosmall, ntoolarge) result(out_of_limit)
  real(kind=sp), intent(in)              ::  array(:)
  real(kind=sp), intent(in), optional    ::  limits(2)
  character(len=*), intent(in), optional ::  array_name
  integer, intent(out), optional         ::  ntoosmall, ntoolarge

  logical                                ::  out_of_limit
  real(kind=sp)                          ::  limits_loc(2)
  integer                                ::  ntoosmall_loc, ntoolarge_loc
  logical                                ::  toosmall(size(array)), toolarge(size(array))
  integer                                ::  i

  out_of_limit = .false.

  if (present(limits)) then
    limits_loc(1) = limits(1) !- tiny(array)
    limits_loc(2) = limits(2) !+ tiny(array)
  else
    limits_loc(1) = -huge(array)
    limits_loc(2) = huge(array)
  end if

  toosmall = (.not.array.ge.limits_loc(1)) ! This catches NaNs as well, which give .false.
                                       ! for every binary comparison.

  toolarge = (.not.array.le.limits_loc(2)) ! This catches NaNs as well, which give .false.
                                       ! for every binary comparison.

  ntoosmall_loc = count(toosmall)
  if (ntoosmall_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array smaller than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Lower limit  :              ', minval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoosmall_loc
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do i = 1, size(array, 1)
          if (.not.array(i).ge.(limits_loc(1))) then
              write(*,*) i, array(i)
          end if
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  ntoolarge_loc = count(toolarge)
  if (ntoolarge_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array larger than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Upper limit  :              ', maxval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoolarge_loc
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do i = 1, size(array, 1)
          if (.not.array(i).le.(limits_loc(2))) then
            write(*,*) i, array(i)
          end if
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  if (present(ntoosmall)) ntoosmall = ntoosmall_loc
  if (present(ntoolarge)) ntoolarge = ntoolarge_loc

end function checklim_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Checks whether any value of 'array' is outside of 'limits_loc' or NaN (2D version)
function checklim_2d(array, limits, array_name, ntoosmall, ntoolarge) result(out_of_limit)
  real(kind=sp), intent(in)              ::  array(:,:)
  real(kind=sp), intent(in), optional    ::  limits(2)
  character(len=*), intent(in), optional ::  array_name
  integer, intent(out), optional         ::  ntoosmall, ntoolarge

  logical                                ::  out_of_limit
  real(kind=sp)                          ::  limits_loc(2)
  integer                                ::  ntoosmall_loc, ntoolarge_loc
  logical                                ::  toosmall(size(array,1), size(array,2))
  logical                                ::  toolarge(size(array,1), size(array,2))
  integer                                ::  i, j

  out_of_limit = .false.

  if (present(limits)) then
    limits_loc(1) = limits(1) !- tiny(array)
    limits_loc(2) = limits(2) !+ tiny(array)
  else
    limits_loc(1) = -huge(array)
    limits_loc(2) = huge(array)
  end if

  toosmall = (.not.array.ge.limits_loc(1)) ! This catches NaNs as well, which give .false.
                                       ! for every binary comparison.

  toolarge = (.not.array.le.limits_loc(2)) ! This catches NaNs as well, which give .false.
                                       ! for every binary comparison.

  ntoosmall_loc = count(toosmall)
  if (ntoosmall_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array smaller than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Lower limit  :              ', minval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoosmall_loc
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do j = 1, size(array, 2)
          do i = 1, size(array, 1)
            if (.not.array(i,j).ge.(limits_loc(1))) then
              write(*,*) i, j, array(i,j)
            end if
          end do
        end do
      end if 
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  ntoolarge_loc = count(toolarge)
  if (ntoolarge_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array larger than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Upper limit  :              ', maxval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoolarge_loc
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do j = 1, size(array, 2)
          do i = 1, size(array, 1)
            if (.not.array(i,j).le.(limits_loc(2))) then
              write(*,*) i, j, array(i,j)
            end if
          end do
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  if (present(ntoosmall)) ntoosmall = ntoosmall_loc
  if (present(ntoolarge)) ntoolarge = ntoolarge_loc

end function checklim_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Checks whether any value of 'array' is outside of 'limits_loc' or NaN
function checklim_3d(array, limits, array_name, ntoosmall, ntoolarge) result(out_of_limit)
  real(kind=sp), intent(in)              ::  array(:,:,:)
  real(kind=sp), intent(in), optional    ::  limits(2)
  character(len=*), intent(in), optional ::  array_name
  integer, optional, intent(out)         ::  ntoosmall, ntoolarge

  logical                                ::  out_of_limit
  real(kind=sp)                          ::  limits_loc(2)
  integer                                ::  ntoosmall_loc, ntoolarge_loc
  logical, allocatable                   ::  toosmall(:,:,:), toolarge(:,:,:)
  integer                                ::  i, j, k

  out_of_limit = .false.

  if (present(limits)) then
    limits_loc(1) = limits(1) !- tiny(array)
    limits_loc(2) = limits(2) !+ tiny(array)
  else
    limits_loc(1) = -huge(array)
    limits_loc(2) = huge(array)
  end if

  allocate(toosmall(size(array,1), size(array,2), size(array,3)))
  allocate(toolarge(size(array,1), size(array,2), size(array,3)))

  toosmall = (.not.array.ge.limits_loc(1)) ! This catches NaNs as well, which give .false.
                                           ! for every binary comparison.

  toolarge = (.not.array.le.limits_loc(2)) ! This catches NaNs as well, which give .false.
                                           ! for every binary comparison.

  ntoosmall_loc = count(toosmall)
  if (ntoosmall_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array smaller than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Lower limit  :              ', minval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoosmall_loc
      write(*,*) 'Details for first violating element: follow'
      if (verbose>1) then
        do k = 1, size(array, 3)
          do j = 1, size(array, 2)
            do i = 1, size(array, 1)
              if (.not.array(i,j,k).ge.(limits_loc(1))) then
                write(*,*) i, j, k, array(i,j,k)
              end if
            end do
          end do
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  ntoolarge_loc = count(toolarge)
  if (ntoolarge_loc>0) then
    if (verbose>0) then
      write(*,*) '**********************************************************************'
      write(*,*) 'ERROR: Value in array larger than limit!'
      write(*,*) 'Variable name:              ', trim(array_name) 
      write(*,*) 'Upper limit  :              ', maxval(limits_loc)
      write(*,*) 'Number of violating values: ', ntoolarge_loc
      if (verbose>1) then
        write(*,*) 'Details for first violating element: follow'
        do k = 1, size(array, 3)
          do j = 1, size(array, 2)
            do i = 1, size(array, 1)
              if (.not.array(i,j,k).le.(limits_loc(2))) then
                write(*,*) i, j, k, array(i,j,k)
              end if
            end do
          end do
        end do
      end if
      write(*,*) '**********************************************************************'
      call flush(6)
    end if
    out_of_limit = .true.
  end if

  if (present(ntoosmall)) ntoosmall = ntoosmall_loc
  if (present(ntoolarge)) ntoolarge = ntoolarge_loc

end function checklim_3d
!-----------------------------------------------------------------------------------------

end module simple_routines
!=========================================================================================
