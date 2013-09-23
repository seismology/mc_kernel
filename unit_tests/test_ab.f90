module test_ab

  use ab
  use ftnunit
  implicit none
  public
contains

subroutine test_trivial
  call assert_equal(4, a_plus_b(1,3), 'a_plus_b(1,3) = 4')
end subroutine

subroutine test_all
  call test(test_trivial, 'TEST a_plus_b')
end subroutine

end module
