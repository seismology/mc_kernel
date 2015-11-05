!=========================================================================================
module test_simple_routines

  use global_parameters
  use simple_routines
  use ftnunit
  implicit none
  public
contains


!------------------------------------------------------------------------------
subroutine test_mult2d_1d
  integer, parameter     :: npoint = 1000, ndim2 = 3
  integer                :: idim2
  real(kind=dp)          :: vect_1d(npoint), vect_2d(npoint, ndim2)
  real(kind=dp)          :: vect_res(npoint, ndim2), vect_ref(npoint, ndim2)

  call random_number(vect_1d) 
  call random_number(vect_2d)

  vect_res = vect_2d
  call mult2d_1d(vect_res, vect_1d)

  do idim2 = 1, ndim2
     vect_ref(:, idim2) = vect_1d(:) * vect_2d(:, idim2)
     call assert_comparable(vect_ref(:, idim2), vect_res(:, idim2), 1d-15, &
                            'mult2d_1d gives same result as loop over dimensions')
  end do
     
end subroutine test_mult2d_1d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_mult3d_1d
  integer, parameter     :: npoint = 1000, ndim2 = 3, ndim3 = 2
  integer                :: idim2, idim3
  real(kind=dp)          :: vect_1d(npoint), vect_3d(npoint, ndim2, ndim3)
  real(kind=dp)          :: vect_res(npoint, ndim2, ndim3), vect_ref(npoint, ndim2, ndim3)

  call random_number(vect_1d) 
  call random_number(vect_3d)

  vect_res = vect_3d
  call mult3d_1d(vect_res, vect_1d)

  do idim2 = 1, ndim2
     do idim3 = 1, ndim3
        vect_ref(:, idim2, idim3) = vect_1d(:) * vect_3d(:, idim2, idim3)
        call assert_comparable(vect_ref(:, idim2, idim3), vect_res(:, idim2, idim3), &
                               1d-15,                                                &
                               'mult2d_1d gives same result as loop over dimensions')
     end do
  end do
     
end subroutine test_mult3d_1d
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine test_mult2d_1d_complex
  integer, parameter     :: npoint = 1000, ndim2 = 3
  integer                :: idim2
  real(kind=dp)          :: vect_1d_real(npoint), vect_2d_real(npoint, ndim2)
  real(kind=dp)          :: vect_1d_imag(npoint), vect_2d_imag(npoint, ndim2)
  complex(kind=dp)       :: vect_1d(npoint), vect_2d(npoint, ndim2)
  complex(kind=dp)       :: vect_res(npoint, ndim2), vect_ref(npoint, ndim2)

  call random_number(vect_1d_real) 
  call random_number(vect_2d_real)
  call random_number(vect_1d_imag) 
  call random_number(vect_2d_imag)

  vect_2d = cmplx(vect_2d_real, vect_2d_imag, kind=dp)
  vect_1d = cmplx(vect_1d_real, vect_1d_imag, kind=dp)

  vect_res = vect_2d
  call mult2d_1d(vect_res, vect_1d)

  do idim2 = 1, ndim2
     vect_ref(:, idim2) = vect_1d(:) * vect_2d(:, idim2)
     call assert_comparable(abs(vect_ref(:, idim2)), abs(vect_res(:, idim2)), 1d-15, &
                            'mult2d_1d gives same result as loop over dimensions')
     call assert_comparable(aimag(vect_ref(:, idim2)), aimag(vect_res(:, idim2)), 1d-15, &
                            'mult2d_1d gives same result as loop over dimensions')
  end do
     
end subroutine test_mult2d_1d_complex
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_mult3d_1d_complex
  integer, parameter     :: npoint = 1000, ndim2 = 3, ndim3 = 2
  integer                :: idim2, idim3
  real(kind=dp)          :: vect_1d_real(npoint), vect_3d_real(npoint, ndim2, ndim3)
  real(kind=dp)          :: vect_1d_imag(npoint), vect_3d_imag(npoint, ndim2, ndim3)
  complex(kind=dp)       :: vect_1d(npoint), vect_3d(npoint, ndim2, ndim3)
  complex(kind=dp)       :: vect_res(npoint, ndim2, ndim3), vect_ref(npoint, ndim2, ndim3)

  call random_number(vect_1d_real) 
  call random_number(vect_3d_real)
  call random_number(vect_1d_imag) 
  call random_number(vect_3d_imag)

  vect_3d = cmplx(vect_3d_real, vect_3d_imag, kind=dp)
  vect_1d = cmplx(vect_1d_real, vect_1d_imag, kind=dp)

  vect_res = vect_3d
  call mult3d_1d(vect_res, vect_1d)

  do idim2 = 1, ndim2
     do idim3 = 1, ndim3
        vect_ref(:, idim2, idim3) = vect_1d(:) * vect_3d(:, idim2, idim3)
        call assert_comparable(abs(vect_ref(:, idim2, idim3)), abs(vect_res(:, idim2, idim3)), &
                               1d-15,                                                &
                               'mult3d_1d gives same result as loop over dimensions')
        call assert_comparable(aimag(vect_ref(:, idim2, idim3)), aimag(vect_res(:, idim2, idim3)), &
                               1d-15,                                                &
                               'mult3d_1d gives same result as loop over dimensions')
     end do
  end do
     
end subroutine test_mult3d_1d_complex
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_to_lower
  character(len=40)    :: string, string_ref

  string = 'LiEsT iRgEnDJEMANd dieSEN CODE? 146 &%$ '
  string_ref = 'liest irgendjemand diesen code? 146 &%$ '

  call assert_true(to_lower(string)==string_ref,  & 
                   'to_lower transforms all uppercase letters, but leaves special ' // &
                   'characters and numbers intact')

end subroutine test_to_lower
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_lowtrim
  character(len=48)    :: string, string_ref

  string = 'LiEsT iRgEnDJEMANd dieSEN CODE? 146 &%$%        '
  string_ref = 'liest irgendjemand diesen code? 146 &%$%'

  call assert_true(lowtrim(string)==string_ref,  & 
                   'trimlow transforms all uppercase letters, but leaves special ' // &
                   'characters and numbers intact and trims')

end subroutine test_lowtrim
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine test_absreldiff_sp
  real(kind=sp)    :: number1, number2, diff
  
  ! Test 1: Compare 1 with 0.999
  number1 = 1
  number2 = 1 + 1d-03
  
  diff = absreldiff(number1, number2)
  call assert_comparable(diff, 1e-03, 1e-03, '1: Relative difference is computed correctly')
  
  ! Test 2: Compare random number with itself times 0.999
  call random_number(number1)
  number2 = number1 - 1e-03*number1

  diff = absreldiff(number1, number2)
  call assert_comparable(diff, 1e-03, 1e-03, '2: Relative difference is computed correctly')
  
  ! Test 3: Compare negative random number with itself times 0.999
  call random_number(number1)
  number1 = - number1              
  number2 =  number1 - 1e-03*number1

  diff = absreldiff(number1, number2)
  call assert_comparable(diff, 1e-03, 1e-03, '3: Relative difference is computed correctly')
  
  ! Test 3: Compare number with itself
  call random_number(number1)
  number2 = number1

  diff = absreldiff(number1, number2)
  call assert_comparable(diff, 0e0, 1e-20, '4: Relative difference is computed correctly')

end subroutine test_absreldiff_sp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_absreldiff_dp
  real(kind=dp)    :: number1, number2, diff
  
  ! Test 1: Compare 1 with 0.999999999999999999
  number1 = 1
  number2 = 1 + 1d-10
  
  diff = absreldiff(number1, number2)
  call assert_comparable(diff, 1d-10, 1d-06, '1: Relative difference is computed correctly')
  
  ! Test 2: Compare random number with itself times 0.9999999999 
  call random_number(number1)
  number2 = number1 - 1d-10*number1

  diff = absreldiff(number1, number2)
  call assert_comparable(diff, 1d-10, 1d-06, '2: Relative difference is computed correctly')
  
  ! Test 3: Compare negative random number with itself times 0.9999999999 
  call random_number(number1)
  number1 = - number1              
  number2 =  number1 - 1d-07*number1

  diff = absreldiff(number1, number2)
  call assert_comparable(diff, 1d-07, 1d-06, '3: Relative difference is computed correctly')
  
  ! Test 3: Compare number with itself
  call random_number(number1)
  number2 = number1

  diff = absreldiff(number1, number2)
  call assert_comparable(diff, 0d0, 1d-20, '4: Relative difference is computed correctly')

end subroutine test_absreldiff_dp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_absreldiff_qp
  real(kind=qp)    :: number1, number2, diff
  
  ! Test 1: Compare 1 with 0.999999999999999999
  number1 = 1
  number2 = 1 + 1d-10
  
  diff = absreldiff(number1, number2)
  call assert_comparable(real(diff, kind=dp), 1d-10, 1d-06, &
                         '1: Relative difference is computed correctly')
  
  ! Test 2: Compare random number with itself times 0.9999999999 
  call random_number(number1)
  number2 = number1 - 1d-10*number1

  diff = absreldiff(number1, number2)
  call assert_comparable(real(diff, kind=dp), 1d-10, 1d-06, &
                         '2: Relative difference is computed correctly')
  
  ! Test 3: Compare negative random number with itself times 0.9999999999 
  call random_number(number1)
  number1 = - number1              
  number2 =  number1 - 1d-07*number1

  diff = absreldiff(number1, number2)
  call assert_comparable(real(diff, kind=dp), 1d-07, 1d-06, &
                         '3: Relative difference is computed correctly')
  
  ! Test 3: Compare number with itself
  call random_number(number1)
  number2 = number1

  diff = absreldiff(number1, number2)
  call assert_comparable(real(diff, kind=dp), 0d0, 1d-20, &
                         '4: Relative difference is computed correctly')

end subroutine test_absreldiff_qp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_cross_sp
  real(kind=sp), dimension(3) :: v1, v2, prod, prod_ref

  ! Test 1: Cross product of pre-defined vectors
  v1 = [3, -2, 9]
  v2 = [-1, 8, 2]

  prod = cross(v1, v2)
  prod_ref =  [-76, -15, 22]

  call assert_comparable(prod, prod_ref, 1e-6, 'Cross product correct')

  ! Test 2: Cross product of random vector with parallel one
  call random_number(v1)
  v2 = 2 * v1

  prod = cross(v1, v2)
  prod_ref = [0, 0, 0]
  call assert_comparable(prod, prod_ref, 1e-6, 'Cross product of parallel vectors is zero')
                
  
  ! Test 3: Cross product of random 2D vector with orthogonal one
  call random_number(v1)
  v1(3) = 0

  v2(1) = v1(2)
  v2(2) = -v1(1)
  v2(3) = 0

  prod = cross(v1, v2)

  call assert_comparable(norm2(prod), norm2(v1)*norm2(v2), 1e-6, &
                         'Length of cross product of orthogonal vectors is product of their length')


  ! Test 4: Cross product of random vectors to check whether v1xv2 is orthogonal to v1,v2
  
  call random_number(v1)
  call random_number(v2)

  prod = cross(v1, v2)

  call assert_comparable(1+dot_product(v1, prod), 1e0, 1e-6, &
                         'Cross product is orthogonal to vector 1')
  call assert_comparable(1+dot_product(v2, prod), 1e0, 1e-6, &
                         'Cross product is orthogonal to vector 2')
end subroutine test_cross_sp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_cross_dp
  real(kind=dp), dimension(3) :: v1, v2, prod, prod_ref

  ! Test 1: Cross product of pre-defined vectors
  v1 = [3, -2, 9]
  v2 = [-1, 8, 2]

  prod = cross(v1, v2)
  prod_ref =  [-76, -15, 22]

  call assert_comparable(prod, prod_ref, 1d-7, 'Cross product correct')

  ! Test 2: Cross product of random vector with parallel one
  call random_number(v1)
  v2 = 2 * v1

  prod = cross(v1, v2)
  prod_ref = [0, 0, 0]
  call assert_comparable(prod, prod_ref, 1d-7, 'Cross product of parallel vectors is zero')
                
  
  ! Test 3: Cross product of random 2D vector with orthogonal one
  call random_number(v1)
  v1(3) = 0

  v2(1) = v1(2)
  v2(2) = -v1(1)
  v2(3) = 0

  prod = cross(v1, v2)

  call assert_comparable(norm2(prod), norm2(v1)*norm2(v2), 1d-7, &
                         'Length of cross product of orthogonal vectors is product of their length')


  ! Test 4: Cross product of random vectors to check whether v1xv2 is orthogonal to v1,v2
  
  call random_number(v1)
  call random_number(v2)

  prod = cross(v1, v2)

  call assert_comparable(1+dot_product(v1, prod), 1d0, 1d-10, &
                         'Cross product is orthogonal to vector 1')
  call assert_comparable(1+dot_product(v2, prod), 1d0, 1d-10, &
                         'Cross product is orthogonal to vector 2')
end subroutine test_cross_dp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_cross_qp
  real(kind=qp), dimension(3) :: v1, v2
  real(kind=dp), dimension(3) :: prod, prod_ref

  ! Test 1: Cross product of pre-defined vectors
  v1 = [3, -2, 9]
  v2 = [-1, 8, 2]

  prod = cross(v1, v2)
  prod_ref =  [-76, -15, 22]

  call assert_comparable(prod, prod_ref, 1d-7, 'Cross product correct')

  ! Test 2: Cross product of random vector with parallel one
  call random_number(v1)
  v2 = 2 * v1

  prod = cross(v1, v2)
  prod_ref = [0, 0, 0]
  call assert_comparable(prod, prod_ref, 1d-7, 'Cross product of parallel vectors is zero')
                
  
  ! Test 3: Cross product of random 2D vector with orthogonal one
  call random_number(v1)
  v1(3) = 0

  v2(1) = v1(2)
  v2(2) = -v1(1)
  v2(3) = 0

  prod = cross(v1, v2)

  call assert_comparable(norm2(prod), norm2(real(v1, kind=dp))*norm2(real(v2, kind=dp)), 1d-15, &
                         'Length of cross product of orthogonal vectors is product of their length')

  ! Test 4: Cross product of random vectors to check whether v1xv2 is orthogonal to v1,v2
  
  call random_number(v1)
  call random_number(v2)

  prod = cross(v1, v2)

  call assert_comparable(real(1+dot_product(v1, prod), kind=dp), 1d0, 1d-15, &
                         'Cross product is orthogonal to vector 1')
  call assert_comparable(real(1+dot_product(v2, prod), kind=dp), 1d0, 1d-15, &
                         'Cross product is orthogonal to vector 2')
end subroutine test_cross_qp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_checklim_1d_int
  integer           ::  array(10)
  integer           ::  ntoosmall, ntoolarge
  logical           ::  out_of_limit

  array = 1

  ! Test without limits. Should always be in limits then
  out_of_limit = check_limits(array     = array)
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test without ntoosmall
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [0, 2] )
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test with all arguments
  out_of_limit = check_limits(array     = array,     &
                              limits    = [0,2],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

  ! One is too large, one too small
  array(:) = 1
  array(3) = -1
  array(7) = 3
  out_of_limit = check_limits(array     = array,     &
                              limits    = [0,2],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')


  ! All are too large
  array(:) = 3
  out_of_limit = check_limits(array     = array,     &
                              limits    = [0,2],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 10, 'one is too large')

  ! Limits are the same
  array(:) = 1
  array(3) = 0
  array(7) = 2
  out_of_limit = check_limits(array     = array,     &
                              limits    = [1,1],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')

  ! Limits are the same - all within
  array = 1
  out_of_limit = check_limits(array     = array,     &
                              limits    = [1,1],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

end subroutine test_checklim_1d_int
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_checklim_2d_int
  integer           ::  array(10, 10)
  integer           ::  ntoosmall, ntoolarge
  logical           ::  out_of_limit

  array(:,:) = 1

  ! Test without limits. Should always be in limits then
  out_of_limit = check_limits(array     = array)
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test without ntoosmall
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [0, 2] )
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test with all arguments
  out_of_limit = check_limits(array     = array,     &
                              limits    = [0,2],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

  ! One is too large, one too small
  array(:,:) = 1
  array(3,1) = -1
  array(7,6) = 3
  out_of_limit = check_limits(array     = array,     &
                              limits    = [0,2],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')


  ! All are too large
  array(:,:) = 3
  out_of_limit = check_limits(array     = array,     &
                              limits    = [0,2],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'all are out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 100, 'one is too large')

  ! Limits are the same
  array(:,:) = 1
  array(3,8) = 0
  array(1,6) = 2
  out_of_limit = check_limits(array     = array,     &
                              limits    = [1,1],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')

  ! Limits are the same - all within
  array(:,:) = 1
  out_of_limit = check_limits(array     = array,     &
                              limits    = [1,1],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

end subroutine test_checklim_2d_int
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_checklim_3d_int
  integer           ::  array(10,10,10)
  integer           ::  ntoosmall, ntoolarge
  logical           ::  out_of_limit

  array(:,:,:) = 1

  ! Test without limits. Should always be in limits then
  out_of_limit = check_limits(array     = array)
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test without ntoosmall
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [0, 2] )
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test with all arguments
  out_of_limit = check_limits(array     = array,     &
                              limits    = [0,2],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

  ! One is too large, one too small
  array(:,:,:) = 1
  array(3,8,2) = -1
  array(1,7,6) = 3
  out_of_limit = check_limits(array     = array,     &
                              limits    = [0,2],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')


  ! All are too large
  array(:,:,:) = 3
  out_of_limit = check_limits(array     = array,     &
                              limits    = [0,2],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 1000, 'one is too large')

  ! Limits are the same
  array(:,:,:) = 1
  array(3,8,1) = 0
  array(7,2,9) = 2
  out_of_limit = check_limits(array     = array,     &
                              limits    = [1,1],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')

  ! Limits are the same - all within
  array(:,:,:) = 1
  out_of_limit = check_limits(array     = array,     &
                              limits    = [1,1],     &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

end subroutine test_checklim_3d_int
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_checklim_1d
  real(kind=sp)     ::  array(10)
  integer           ::  ntoosmall, ntoolarge
  logical           ::  out_of_limit

  array(:) = 1.0

  ! Test without limits. Should always be in limits then
  out_of_limit = check_limits(array     = array)
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test without ntoosmall
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)])
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test with all arguments
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)], &
                              ntoosmall = ntoosmall,                            &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

  ! One is too large, one too small
  array(:) = 1.0
  array(3) = 1.0 - 2*epsilon(1.0)
  array(7) = 1.0 + 2*epsilon(1.0)
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)], &
                              ntoosmall = ntoosmall,                            &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')

  ! All are too large
  array(:) = 2.0
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)], &
                              ntoosmall = ntoosmall,                            &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 10, '10 (all) are too large')

  ! Limits are the same
  array(:) = 1.0
  array(3) = 1.0 - 2*epsilon(1.0)
  array(7) = 1.0 + 2*epsilon(1.0)
  out_of_limit = check_limits(array     = array,      &
                              limits    = [1.0, 1.0], &
                              ntoosmall = ntoosmall,  &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')

  ! Limits are the same - all within
  array(:) = 1.0
  out_of_limit = check_limits(array     = array,     &
                              limits    = [1.0,1.0], &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

end subroutine test_checklim_1d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_checklim_2d
  real(kind=sp)     ::  array(10,10)
  integer           ::  ntoosmall, ntoolarge
  logical           ::  out_of_limit

  array(:,:) = 1.0

  ! Test without limits. Should always be in limits then
  out_of_limit = check_limits(array     = array)
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test without ntoosmall
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)])
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test with all arguments
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)], &
                              ntoosmall = ntoosmall,                            &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

  ! One is too large, one too small
  array(:,:) = 1.0
  array(3,2) = 1.0 - 2*epsilon(1.0)
  array(7,9) = 1.0 + 2*epsilon(1.0)
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)], &
                              ntoosmall = ntoosmall,                            &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')

  ! All are too large
  array(:,:) = 2.0
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)], &
                              ntoosmall = ntoosmall,                            &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 100, '100 (all) are too large')

  ! Limits are the same
  array(:,:) = 1.0
  array(3,3) = 1.0 - 2*epsilon(1.0)
  array(7,1) = 1.0 + 2*epsilon(1.0)
  out_of_limit = check_limits(array     = array,      &
                              limits    = [1.0, 1.0], &
                              ntoosmall = ntoosmall,  &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')

  ! Limits are the same - all within
  array(:,:) = 1.0
  out_of_limit = check_limits(array     = array,     &
                              limits    = [1.0,1.0], &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

end subroutine test_checklim_2d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_checklim_3d
  real(kind=sp)     ::  array(10,10,10)
  integer           ::  ntoosmall, ntoolarge
  logical           ::  out_of_limit

  array(:,:,:) = 1.0

  ! Test without limits. Should always be in limits then
  out_of_limit = check_limits(array     = array)
  call assert_false(out_of_limit, 'none is out of limit')

  ! Test without ntoosmall
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)])
  call assert_false(out_of_limit, 'none is out of limit')


  ! Test with all arguments
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)], &
                              ntoosmall = ntoosmall,                            &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

  ! One is too large, one too small
  array(:,:,:) = 1.0
  array(3,9,8) = 1.0 - 2*epsilon(1.0)
  array(7,2,3) = 1.0 + 2*epsilon(1.0)
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)], &
                              ntoosmall = ntoosmall,                            &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')

  ! All are too large
  array(:,:,:) = 2.0
  out_of_limit = check_limits(array     = array,                                &
                              limits    = [1.0-epsilon(1.0), 1.0+epsilon(1.0)], &
                              ntoosmall = ntoosmall,                            &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 1000, '1000 (all) are too large')

  ! Limits are the same
  array(:,:,:) = 1.0
  array(3,8,1) = 1.0 - 2*epsilon(1.0)
  array(7,2,1) = 1.0 + 2*epsilon(1.0)
  out_of_limit = check_limits(array     = array,      &
                              limits    = [1.0, 1.0], &
                              ntoosmall = ntoosmall,  &
                              ntoolarge = ntoolarge )

  call assert_true(out_of_limit, 'some are out of limit')
  call assert_equal(ntoosmall, 1, 'one is too small')
  call assert_equal(ntoolarge, 1, 'one is too large')

  ! Limits are the same - all within
  array(:,:,:) = 1.0
  out_of_limit = check_limits(array     = array,     &
                              limits    = [1.0,1.0], &
                              ntoosmall = ntoosmall, &
                              ntoolarge = ntoolarge )

  call assert_false(out_of_limit, 'none is out of limit')
  call assert_equal(ntoosmall, 0, 'none is too small')
  call assert_equal(ntoolarge, 0, 'none is too large')

end subroutine test_checklim_3d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_cumsum_trapezoidal_1d
  real(kind=dp), allocatable :: data_in(:), data_out(:), data_ref(:)
  real(kind=dp)              :: dt, misfit
  integer                    :: lu_ref, lu_res, len_ref, isample

  open(newunit=lu_ref, file='./cumsum_1d_ref', action='read')
  read(lu_ref, *) len_ref, dt
  allocate(data_in(len_ref))
  allocate(data_ref(len_ref))
  allocate(data_out(len_ref))

  do isample = 1, len_ref
    read(lu_ref, *) data_in(isample), data_ref(isample)
  end do
  close(lu_ref)

  data_out = cumsum_trapezoidal(data_in, dt)

  misfit = norm2(data_out -  &
                 data_ref) / &
           norm2(data_ref)
  call assert_true(misfit<1d-6, 'misfit is small')

  open(newunit=lu_res, file='./output/cumsum_1d_res', action='write')
  do isample = 1, len_ref
    write(lu_res, *) data_ref(isample), data_out(isample)
  end do
  close(lu_res)

end subroutine test_cumsum_trapezoidal_1d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_cumsum_trapezoidal_2d
  real(kind=dp), allocatable :: data_in(:,:), data_out(:,:), data_ref(:,:)
  real(kind=dp)              :: dt, misfit
  integer                    :: lu_ref, lu_res, len_ref, isample

  open(newunit=lu_ref, file='./cumsum_2d_ref', action='read')
  read(lu_ref, *) len_ref, dt
  allocate(data_in(len_ref, 2))
  allocate(data_ref(len_ref, 2))
  allocate(data_out(len_ref, 2))

  do isample = 1, len_ref
    read(lu_ref, *) data_in(isample, :), data_ref(isample, :)
  end do
  close(lu_ref)

  data_out = cumsum_trapezoidal(data_in, dt)

  misfit = norm2(data_out -  &
                 data_ref) / &
           norm2(data_ref)
  call assert_true(misfit<1d-6, 'misfit is small')

  open(newunit=lu_res, file='./output/cumsum_2d_res', action='write')
  do isample = 1, len_ref
    write(lu_res, *) data_ref(isample, :), data_out(isample, :)
  end do
  close(lu_res)

end subroutine test_cumsum_trapezoidal_2d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_cumsum_trapezoidal_3d
  real(kind=dp), allocatable :: data_in(:,:,:), data_out(:,:,:), data_ref(:,:,:)
  real(kind=dp)              :: dt, misfit
  integer                    :: lu_ref, lu_res, len_ref, isample

  open(newunit=lu_ref, file='./cumsum_3d_ref', action='read')
  read(lu_ref, *) len_ref, dt
  allocate(data_in(len_ref, 2, 2))
  allocate(data_ref(len_ref, 2, 2))
  allocate(data_out(len_ref, 2, 2))

  do isample = 1, len_ref
    read(lu_ref, *) data_in(isample, :, :), data_ref(isample, :, :)
  end do
  close(lu_ref)

  data_out = cumsum_trapezoidal(data_in, dt)

  misfit = norm2(data_out -  &
                 data_ref) / &
           norm2(data_ref)
  call assert_true(misfit<1d-6, 'misfit is small')

  open(newunit=lu_res, file='./output/cumsum_3d_res', action='write')
  do isample = 1, len_ref
    write(lu_res, *) data_ref(isample, :, :), data_out(isample, :, :)
  end do
  close(lu_res)

end subroutine test_cumsum_trapezoidal_3d
!------------------------------------------------------------------------------

end module test_simple_routines
!=========================================================================================
