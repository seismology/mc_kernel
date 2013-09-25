program test_ftnunit

  use ftnunit
  use test_montecarlo
  use test_fft
  implicit none

  call runtests_init
  call runtests( test_all )
  call runtests_final


  contains

  subroutine test_all
    call test(test_unit_hexagon, 'TEST MC_unit_hexagon')
    call test(test_fft_dirac, 'TEST FFT_dirac')
  end subroutine

end program
