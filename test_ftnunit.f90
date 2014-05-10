!=========================================================================================
module unit_tests

  use ftnunit, only: test
  use global_parameters, only: lu_out
  use test_montecarlo
  use test_fft
  use test_tetrahedra
  use test_inversion_mesh
  use test_buffer
  use test_filter
  use test_kernel
  use test_readfields
  use test_resampling

  implicit none

contains
!-----------------------------------------------------------------------------------------
subroutine test_all

  call init_output()

  ! test_readfields
  write(6,'(/,a)') 'TEST READFIELDS MODULE'
  call test(test_readfields_rotate, 'Field rotation')
  call test(test_readfields_rotate_straintensor, 'Strainfield rotation')
  call test(test_readfields_rotate_straintensor_voigt, 'Strainfield rotation voigt')
  call test(test_readfields_set_params, 'Set SEM file params')
  call test(test_readfields_open_files, 'Open SEM file')
  call test(test_readfields_load_seismogram, 'Load seismogram')
   
  ! test_montecarlo
  write(6,'(/,a)') 'TEST MONTECARLO MODULE'
  call test(test_mc_meanandvariance, 'MC mean and variance')
  call test(test_mc_unit_hexagon, 'MC unit hexagon')
  call test(test_mc_sphere_in_tetrahedron, 'MC sphere in tetrahedron')

  ! test_fft
  write(6,'(/,a)') 'TEST FFT MODULE'
  call test(test_fft_dirac, 'FFT_dirac')
  call test(test_fft_sine, 'FFT_sine')
  call test(test_fft_inverse, 'FFT_inverse')
  call test(test_fft_convolve, 'FFT_convolve')
  call test(test_fft_taperandzeropad, 'FFT_taperandzeropad')

  ! test_resampling
  write(6,'(/,a)') 'TEST RESAMPLING MODULE'
  call test(test_resampling_const, 'RESAMPLING_const')
  call test(test_resampling_const_ntraces, 'RESAMPLING_const_ntraces')
  call test(test_resampling_triangle, 'RESAMPLING_triangle')

  ! test filter
  write(6,'(/,a)') 'TEST FILTER MODULE'
  call test(test_filter_gabor_response, 'Test Gabor filter')
  call test(test_filter_timeshift, 'Test time shifting routine')
  
  ! test kernel
  write(6,'(/,a)') 'TEST KERNEL MODULE'
  call test(test_kernel_init, 'Test Kernel initialization')
  call test(test_kernel_cut_timewindow, 'Test Time window cutting')

  ! test_tetrahedra
  write(6,'(/,a)') 'TEST TETRAHEDRON MODULE'
  call test(test_generate_random_point_poly_3, 'Random points in Triangle')
  call test(test_generate_random_point_poly_4, 'Random points in Quadrilateral')
  call test(test_generate_random_point_tet, 'Random points in Tetrahedra')
  call test(test_generate_random_point_triangle_space, 'Random points in Triangle in Space')
  call test(test_rmat4_det, 'Matrix determinant')
  call test(test_tetra_volume_3d, 'Tetrahedron volume')
  call test(test_get_volume_poly_3, 'Triangle area')
  call test(test_get_volume_poly_4, 'Quadrilateral area')

  ! test_inversion_mesh
  write(6,'(/,a)') 'TEST INVERSION MESH MODULE'
  call test(test_mesh_read, 'reading tetrahedral mesh')
  call test(test_mesh_dump, 'reading/dumping tetrahedral mesh')
  call test(test_mesh_dump2, 'reading/dumping tetrahedral mesh from abaqus')
  call test(test_mesh_data_dump, 'reading/dumping tetrahedral mesh with data')
  call test(test_mesh_data_dump2, &
            'reading/dumping triangular mesh from abaqus file with data')
  call test(test_mesh_data_dump3, &
            'reading/dumping quadrilateral mesh from abaqus file with data')
  call test(test_mesh_data_dump4, &
            'reading/dumping hexahedral mesh from abaqus file with data')
  call test(test_mesh_data_dump5, &
            'reading/dumping tetrahedral mesh from abaqus file with data')
  call test(test_valence, 'computation of valence')
  call test(test_get_connected_elements, 'get connected elements')
  call test(test_initialize_mesh, 'initialize mesh')

  ! test_buffer
  write(6,'(/,a)') 'TEST BUFFER MODULE'
  call test(test_buffer_storage, 'put data into the buffer')
  call test(test_buffer_retrieval, 'get data back from the buffer')
  call test(test_buffer_overwrite, 'buffer gets overwritten after time')

  call finish_output()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_output()
  character(len=11) :: fnam

  fnam = 'OUTPUT_test'
  open(newunit=lu_out, file=fnam, status='unknown', position='append')
  write(lu_out,*) '*********************************************************************'
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine finish_output()
  close(lu_out)
end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
