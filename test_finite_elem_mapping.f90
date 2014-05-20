!=========================================================================================
module test_finite_elem_mapping

   use global_parameters
   use finite_elem_mapping
   use ftnunit
   implicit none
   public

contains

!-----------------------------------------------------------------------------------------
subroutine test_mapping_subpar_xieta_to_sz()

   real(kind=dp)    :: xi, eta, sz(2), sz_ref(2)
   real(kind=dp)    :: nodes(4,2)

! 4 - - - - - - - 3
! |       ^       |
! |   eta |       |
! |       |       |
! |        --->   |
! |        xi     |
! |               |
! |               |
! 1 - - - - - - - 2
    
   nodes(1,:) = [-1,-1]
   nodes(2,:) = [ 1,-1]
   nodes(3,:) = [ 1, 1]
   nodes(4,:) = [-1, 1]

   xi = 0
   eta = 0
   sz_ref = [xi, eta]
   sz = mapping_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(sz), 1 + real(sz_ref), &
                                 1e-7, 'ref to ref is identity, [0 ,0]')

   call random_number(xi)
   call random_number(eta)
   xi  =  xi * 2 - 1
   eta = eta * 2 - 1
   sz_ref = [xi, eta]
   sz = mapping_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(sz), 1 + real(sz_ref), &
                                 1e-7, 'ref to ref is identity, random xi, eta')

   nodes(1,:) = [0,0]
   nodes(2,:) = [2,0]
   nodes(3,:) = [2,2]
   nodes(4,:) = [0,2]

   call random_number(xi)
   call random_number(eta)
   xi  =  xi * 2 - 1
   eta = eta * 2 - 1
   sz_ref = [xi+1, eta+1]
   sz = mapping_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(sz), 1 + real(sz_ref), &
                                 1e-7, 'pure translation, random xi, eta')
   
   nodes(1,:) = [ 0, 0]
   nodes(2,:) = [10, 0]
   nodes(3,:) = [10,20]
   nodes(4,:) = [ 0,20]

   call random_number(xi)
   call random_number(eta)
   xi  =  xi * 2 - 1
   eta = eta * 2 - 1
   sz_ref = [xi * 5 + 5, eta * 10 + 10]
   sz = mapping_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(sz), 1 + real(sz_ref), &
                                 1e-7, 'stretching + translation, random xi, eta')
   
   nodes(2,:) = [ 0, 0]
   nodes(3,:) = [10, 0]
   nodes(4,:) = [10,20]
   nodes(1,:) = [ 0,20]

   call random_number(xi)
   call random_number(eta)
   xi  =  xi * 2 - 1
   eta = eta * 2 - 1
   sz_ref = [eta * 5 + 5, -xi * 10 + 10]
   sz = mapping_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(sz), 1 + real(sz_ref), &
                                 1e-7, 'rotation+ stretching + translation, random xi, eta')

end subroutine test_mapping_subpar_xieta_to_sz
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_jacobian_subpar()

   real(kind=dp)    :: xi, eta, jacobian(2,2), jacobian_ref(2,2)
   real(kind=dp)    :: nodes(4,2)

!            | ds / dxi  ds / deta |
! jacobian = |                     |
!            | dz / dxi  dz / deta |

   nodes(1,:) = [-1,-1]
   nodes(2,:) = [ 1,-1]
   nodes(3,:) = [ 1, 1]
   nodes(4,:) = [-1, 1]

   call random_number(xi)
   call random_number(eta)
   jacobian_ref(1,:) = [1,0]
   jacobian_ref(2,:) = [0,1]
   jacobian = jacobian_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(reshape(jacobian, [4])), &
                                 1 + real(reshape(jacobian_ref, [4])), &
                                 1e-7, 'ref to ref, jacobian is identity')

   nodes(2,:) = [-1,-1]
   nodes(3,:) = [ 1,-1]
   nodes(4,:) = [ 1, 1]
   nodes(1,:) = [-1, 1]

   call random_number(xi)
   call random_number(eta)
   jacobian_ref(1,:) = [ 0,1]
   jacobian_ref(2,:) = [-1,0]
   jacobian = jacobian_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(reshape(jacobian, [4])), &
                                 1 + real(reshape(jacobian_ref, [4])), &
                                 1e-7, 'pure rotation')
   
   nodes(1,:) = [0,0]
   nodes(2,:) = [2,0]
   nodes(3,:) = [2,2]
   nodes(4,:) = [0,2]

   call random_number(xi)
   call random_number(eta)
   jacobian_ref(1,:) = [1,0]
   jacobian_ref(2,:) = [0,1]
   jacobian = jacobian_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(reshape(jacobian, [4])), &
                                 1 + real(reshape(jacobian_ref, [4])), &
                                 1e-7, 'pure translation, jacobian is identity')

   nodes(1,:) = [ 0, 0]
   nodes(2,:) = [10, 0]
   nodes(3,:) = [10,20]
   nodes(4,:) = [ 0,20]

   call random_number(xi)
   call random_number(eta)
   jacobian_ref(1,:) = [5, 0]
   jacobian_ref(2,:) = [0,10]
   jacobian = jacobian_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(reshape(jacobian, [4])), &
                                 1 + real(reshape(jacobian_ref, [4])), &
                                 1e-7, 'continuous stretching, jacobian is constant')

end subroutine test_jacobian_subpar
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_inv_jacobian_subpar()

!                | dxi  / ds  dxi  / dz |
! inv_jacobian = |                      |
!                | deta / ds  deta / dz |

   real(kind=dp)    :: xi, eta, inv_jacobian(2,2), inv_jacobian_ref(2,2)
   real(kind=dp)    :: nodes(4,2)

   nodes(1,:) = [-1,-1]
   nodes(2,:) = [ 1,-1]
   nodes(3,:) = [ 1, 1]
   nodes(4,:) = [-1, 1]

   call random_number(xi)
   call random_number(eta)
   inv_jacobian_ref(1,:) = [1,0]
   inv_jacobian_ref(2,:) = [0,1]
   inv_jacobian = inv_jacobian_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(reshape(inv_jacobian, [4])), &
                                 1 + real(reshape(inv_jacobian_ref, [4])), &
                                 1e-7, 'ref to ref, inv_jacobian is identity')

   nodes(2,:) = [-1,-1]
   nodes(3,:) = [ 1,-1]
   nodes(4,:) = [ 1, 1]
   nodes(1,:) = [-1, 1]

   call random_number(xi)
   call random_number(eta)
   inv_jacobian_ref(1,:) = [0,-1]
   inv_jacobian_ref(2,:) = [1, 0]
   inv_jacobian = inv_jacobian_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(reshape(inv_jacobian, [4])), &
                                 1 + real(reshape(inv_jacobian_ref, [4])), &
                                 1e-7, 'pure rotation')
   
   nodes(1,:) = [0,0]
   nodes(2,:) = [2,0]
   nodes(3,:) = [2,2]
   nodes(4,:) = [0,2]

   call random_number(xi)
   call random_number(eta)
   inv_jacobian_ref(1,:) = [1,0]
   inv_jacobian_ref(2,:) = [0,1]
   inv_jacobian = inv_jacobian_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(reshape(inv_jacobian, [4])), &
                                 1 + real(reshape(inv_jacobian_ref, [4])), &
                                 1e-7, 'pure translation, inv_jacobian is identity')

   nodes(1,:) = [ 0, 0]
   nodes(2,:) = [10, 0]
   nodes(3,:) = [10,20]
   nodes(4,:) = [ 0,20]

   call random_number(xi)
   call random_number(eta)
   inv_jacobian_ref(1,:) = [0.2d0, 0.0d0]
   inv_jacobian_ref(2,:) = [0.0d0, 0.1d0]
   inv_jacobian = inv_jacobian_subpar(xi, eta, nodes)

   call assert_comparable_real1d(1 + real(reshape(inv_jacobian, [4])), &
                                 1 + real(reshape(inv_jacobian_ref, [4])), &
                                 1e-7, 'continuous stretching, inv_jacobian is constant')

end subroutine test_inv_jacobian_subpar
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
