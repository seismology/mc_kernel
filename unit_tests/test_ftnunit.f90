!=========================================================================================
program test_ftnunit

  use ftnunit
  use test_montecarlo
  use test_fft
  use test_tetrahedra
  use test_inversion_mesh
  implicit none

  call runtests_init
  call runtests( test_all )
  call runtests_final

contains

!-----------------------------------------------------------------------------------------
subroutine test_all
  ! test_montecarlo
  call test(test_unit_hexagon, 'TEST MC_unit_hexagon')
  call test(test_sphere_in_tetrahedron, 'TEST MC sphere in tetrahedron')

  ! test_fft
  call test(test_fft_dirac, 'TEST FFT_dirac')
  call test(test_fft_inverse, 'TEST FFT_inverse')

  ! test_tetrahedra
  call test(test_generate_random_point, 'TEST Random points in Tetrahedra')
  call test(test_rmat4_det, 'TEST Matrix determinant')
  call test(test_tetra_volume_3d, 'TEST Tetrahedron volume')

  ! test_inversion_mesh
  call test(test_mesh_read, 'TEST reading tetrahedral mesh')
  call test(test_mesh_dump, 'TEST reading/dumping tetrahedral mesh')
end subroutine
!-----------------------------------------------------------------------------------------

end program
!=========================================================================================
