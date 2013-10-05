!=========================================================================================
program test_ftnunit

  use ftnunit
  use test_montecarlo
  use test_fft
  use test_tetrahedra
  use test_inversion_mesh
  use test_buffer
  use test_filter
  implicit none

  call runtests_init
  call runtests( test_all )
  call runtests_final

contains

!-----------------------------------------------------------------------------------------
subroutine test_all
  ! test_montecarlo
  write(6,'(/,a)') 'TEST MONTECARLO MODULE'
  call test(test_unit_hexagon, 'MC_unit_hexagon')
  call test(test_sphere_in_tetrahedron, 'MC sphere in tetrahedron')

  ! test_fft
  write(6,'(/,a)') 'TEST FFT MODULE'
  call test(test_fft_dirac, 'FFT_dirac')
  call test(test_fft_sine, 'FFT_sine')
  call test(test_fft_inverse, 'FFT_inverse')
  call test(test_fft_convolve, 'FFT_convolve')
  call test(test_fft_taperandzeropad, 'FFT_taperandzeropad')

  ! test filter
  write(6,'(/,a)') 'TEST FILTER MODULE'
  call test(test_filter_gabor_response, 'Test Gabor filter')

  ! test_tetrahedra
  write(6,'(/,a)') 'TEST TETRAHEDRON MODULE'
  call test(test_generate_random_point, 'Random points in Tetrahedra')
  call test(test_rmat4_det, 'Matrix determinant')
  call test(test_tetra_volume_3d, 'Tetrahedron volume')

  ! test_inversion_mesh
  write(6,'(/,a)') 'TEST INVERSION MESH MODULE'
  call test(test_mesh_read, 'reading tetrahedral mesh')
  call test(test_mesh_dump, 'reading/dumping tetrahedral mesh')
  call test(test_mesh_data_dump, 'reading/dumping tetrahedral mesh with data')
  call test(test_mesh_data_dump2, &
            'reading/dumping triangular mesh from abaqus file with data')
  call test(test_mesh_data_dump3, &
            'reading/dumping quadrilateral mesh from abaqus file with data')
  call test(test_mesh_data_dump4, &
            'reading/dumping hexahedral mesh from abaqus file with data')
  call test(test_valence, 'computation of valence')
  call test(test_get_connected_elements, 'get connected elements')

  ! test_buffer
  write(6,'(/,a)') 'TEST BUFFER MODULE'
  call test(test_buffer_storage, 'put data into the buffer')
  call test(test_buffer_retrieval, 'get data back from the buffer')
end subroutine
!-----------------------------------------------------------------------------------------

end program
!=========================================================================================
