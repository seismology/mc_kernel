!=========================================================================================
module sem_derivatives
  use global_parameters,      only : dp
  use commpi,                 only : pabort
  use finite_elem_mapping,    only : inv_jacobian

  implicit none
  private

  public :: axisym_gradient

  interface axisym_gradient
    module procedure  :: axisym_gradient
    module procedure  :: axisym_gradient_td
  end interface

  interface mxm
    module procedure  :: mxm
    module procedure  :: mxm_atd
    module procedure  :: mxm_btd
  end interface

contains

!-----------------------------------------------------------------------------------------
function axisym_gradient_td(f, G, GT, xi, eta, npol, nsamp, nodes, element_type)
  ! Computes the axisymmetric gradient of scalar field f
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: f(1:nsamp,0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: axisym_gradient_td(1:nsamp,0:npol,0:npol,1:2)

  real(kind=dp)                 :: inv_j_npol(0:npol,0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm1(1:nsamp,0:npol,0:npol)
  real(kind=dp)                 :: mxm2(1:nsamp,0:npol,0:npol)

  do ipol = 0, npol
     do jpol = 0, npol
        inv_j_npol(ipol,jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
     enddo
  enddo

!        | dxi  / ds  dxi  / dz |
! J^-1 = |                      |
!        | deta / ds  deta / dz |

  mxm1 = mxm(GT,f)
  mxm2 = mxm(f,G)

  do jpol = 0, npol
     do ipol = 0, npol
        axisym_gradient_td(:,ipol,jpol,1) =   &
                inv_j_npol(ipol,jpol,1,1) * mxm1(:,ipol,jpol) &
              + inv_j_npol(ipol,jpol,2,1) * mxm2(:,ipol,jpol)
        axisym_gradient_td(:,ipol,jpol,2) =   &
                inv_j_npol(ipol,jpol,1,2) * mxm1(:,ipol,jpol) &
              + inv_j_npol(ipol,jpol,2,2) * mxm2(:,ipol,jpol)
     enddo
  enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function axisym_gradient(f, G, GT, xi, eta, npol, nodes, element_type)
  ! Computes the axisymmetric gradient of scalar field f
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: f(0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: axisym_gradient(0:npol,0:npol,1:2)

  real(kind=dp)                 :: inv_j_npol(0:npol,0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm1(0:npol,0:npol)
  real(kind=dp)                 :: mxm2(0:npol,0:npol)

  do ipol = 0, npol
     do jpol = 0, npol
        inv_j_npol(ipol,jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
     enddo
  enddo

!        | dxi  / ds  dxi  / dz |
! J^-1 = |                      |
!        | deta / ds  deta / dz |

  mxm1 = mxm(GT,f)
  mxm2 = mxm(f,G)
  axisym_gradient(:,:,1) =   inv_j_npol(:,:,1,1) * mxm1 &
                           + inv_j_npol(:,:,2,1) * mxm2
  axisym_gradient(:,:,2) =   inv_j_npol(:,:,1,2) * mxm1 &
                           + inv_j_npol(:,:,2,2) * mxm2

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! a is time dependent
pure function mxm_atd(a, b)

  real(kind=dp), intent(in)  :: a(1:,0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_atd(1:size(a,1), 0:size(a,2)-1,0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_atd = 0

  do j = 0, size(b,2) -1
     do i = 0, size(a,2) -1
        do k = 0, size(a,3) -1
           mxm_atd(:,i,j) = mxm_atd(:,i,j) + a(:,i,k) * b(k,j)
        enddo
     end do
  end do 

end function mxm_atd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! b is time dependent
pure function mxm_btd(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(1:,0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_btd(1:size(b,1),0:size(a,1)-1,0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_btd = 0

  do j = 0, size(b,2) -1
     do i = 0, size(a,1) -1
        do k = 0, size(a,2) -1
           mxm_btd(:,i,j) = mxm_btd(:,i,j) + a(i,k) * b(:,k,j)
        enddo
     end do
  end do 

end function mxm_btd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
pure function mxm(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm(0:size(a,1)-1,0:size(b,2)-1)    !< Result
  integer                    :: i, j

  do j = 0, size(b,2) -1
     do i = 0, size(a,1) -1
        mxm(i,j) = sum(a(i,:) * b(:,j))
     end do
  end do 

end function mxm
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
