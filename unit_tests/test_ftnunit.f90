program test_ftnunit

  use ftnunit
  use ab
  use test_ab
  implicit none

  call runtests_init
  call runtests( test_all )
  call runtests_final

end program
