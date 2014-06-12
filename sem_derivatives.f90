!=========================================================================================
module sem_derivatives
    use global_parameters,      only : dp
    use commpi,                 only : pabort
    use finite_elem_mapping,    only : inv_jacobian

    implicit none
    private

    public :: axisym_gradient

contains

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
!! Size is fixed to npol x npol
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
