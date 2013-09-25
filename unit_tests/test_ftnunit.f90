program test_ftnunit

  use ftnunit
  use test_montecarlo
  use test_fft
  use test_tetrahedra
  implicit none

  call runtests_init
  call runtests( test_all )
  call runtests_final


  contains

  subroutine test_all
    call test(test_unit_hexagon, 'TEST MC_unit_hexagon')
    call test(test_fft_dirac, 'TEST FFT_dirac')
    call test(test_fft_inverse, 'TEST FFT_inverse')
    call test(test_generate_random_point, 'TEST Random points in Tetrahedra')
    call test(test_rmat4_det, 'TEST Matrix determinant')
    call test(test_tetra_volume_3d, 'TEST Tetrahedron volume')
  end subroutine

end program
