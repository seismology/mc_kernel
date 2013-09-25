program test_ftnunit

  use ftnunit
  use test_mc
  implicit none

  call runtests_init
  call runtests( test_all )
  call runtests_final

end program
