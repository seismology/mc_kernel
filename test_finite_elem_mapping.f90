!=========================================================================================
module test_finite_elem_mapping

   use global_parameters
   use finite_elem_mapping
   use ftnunit
   implicit none
   public

contains

!-----------------------------------------------------------------------------------------
subroutine test_mapping_xieta_to_sz()

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

end subroutine test_mapping_xieta_to_sz
!-----------------------------------------------------------------------------------------
end module
!=========================================================================================
