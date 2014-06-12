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
  use test_rotations
  use test_resampling
  use test_finite_elem_mapping
  use test_spectral_basis
  use test_sem_derivatives

  implicit none

contains
!-----------------------------------------------------------------------------------------
subroutine test_all

  call init_output()

  ! test sem derivatives
  write(6,'(/,a)') 'TEST SEM DERIVATIVE MODULE'
  call test(test_gradient, 'SEM gradient')

  ! test spectral basis functions
  write(6,'(/,a)') 'TEST SPECTRAL BASIS MODULE'
  call test(test_lagrange_interpol_1D, 'lagrangian interpolation in 1D')
  call test(test_lagrange_interpol_2D, 'lagrangian interpolation in 2D')
  call test(test_lagrange_interpol_2D_td, 'lagrangian interpolation in 2D time dep')

  call test(test_gll_points, 'gll points')
  call test(test_glj_points, 'glj points')

  call test(test_derivative_tensors, 'derivative tensors')

  ! test_finite_elem_mapping
  write(6,'(/,a)') 'TEST FINITE ELEMENT MODULE'
  call test(test_inside_element, 'inside element?')

  call test(test_mapping_semino_xieta_to_sz, 'semino mapping (xi,eta) to (s,z)')
  call test(test_inv_mapping_semino_sz_to_xieta, 'inv semino mapping (s,z) to (xi,eta)')
  call test(test_jacobian_semino, 'jacobian semino')
  call test(test_inv_jacobian_semino, 'inverse jacobian semino')
  
  call test(test_mapping_semiso_xieta_to_sz, 'semiso mapping (xi,eta) to (s,z)')
  call test(test_inv_mapping_semiso_sz_to_xieta, 'inv semiso mapping (s,z) to (xi,eta)')
  call test(test_jacobian_semiso, 'jacobian semiso')
  call test(test_inv_jacobian_semiso, 'inverse jacobian semiso')

  call test(test_mapping_spheroidal_xieta_to_sz, 'spheroidal mapping (xi,eta) to (s,z)')
  call test(test_inv_mapping_spheroidal_sz_to_xieta, 'inv spheroidal mapping (s,z) to (xi,eta)')
  call test(test_jacobian_spheroidal, 'jacobian spheroidal')
  call test(test_inv_jacobian_spheroid, 'inverse jacobian spheroid')

  call test(test_mapping_subpar_xieta_to_sz, 'subpar mapping (xi,eta) to (s,z)')
  call test(test_inv_mapping_subpar_sz_to_xieta, 'inv subpar mapping (s,z) to (xi,eta)')
  call test(test_jacobian_subpar, 'jacobian subpar')
  call test(test_inv_jacobian_subpar, 'inverse jacobian subpar')


  ! test_readfields
  write(6,'(/,a)') 'TEST READFIELDS MODULE'
  call test(test_readfields_set_params, 'Set SEM file params')
  call test(test_readfields_open_files, 'Open SEM file')
  call test(test_readfields_load_seismogram, 'Load seismogram')

  ! test_rotations
  write(6,'(/,a)') 'TEST ROTATIONS MODULE'
  call test(test_readfields_rotate, 'Field rotation')
  call test(test_readfields_rotate_straintensor, 'Strainfield rotation')
  call test(test_readfields_rotate_straintensor_voigt, 'Strainfield rotation voigt')
  call test(test_rotate_symm_tensor_voigt_src_to_xyz_1d, 'symm tensor rotation src to xyz - 1d')
  call test(test_rotate_symm_tensor_voigt_src_to_xyz_2d, 'symm tensor rotation src to xyz - 2d')
  call test(test_rotate_symm_tensor_voigt_xyz_to_src_1d, 'symm tensor rotation xyz to src - 1d')
  call test(test_rotate_symm_tensor_voigt_xyz_to_src_2d, 'symm tensor rotation xyz to src - 2d')
  call test(test_rotate_symm_tensor_voigt_xyz_src_to_xyz_earth_1d, 'symm tensor rotation xyz src to earth - 1d')
  call test(test_rotate_symm_tensor_voigt_xyz_src_to_xyz_earth_2d, 'symm tensor rotation xyz src to earth - 2d')
  call test(test_rotate_symm_tensor_voigt_xyz_earth_to_xyz_src_1d, 'symm tensor rotation xyz earth to src - 1d')
  call test(test_rotate_symm_tensor_voigt_xyz_earth_to_xyz_src_2d, 'symm tensor rotation xyz earth to src - 2d')
   
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
  call test(test_buffer_storage_1d, 'put 1d data into the buffer')
  call test(test_buffer_storage_2d, 'put 2d data into the buffer')
  call test(test_buffer_storage_3d, 'put 3d data into the buffer')
  call test(test_buffer_retrieval_1d, 'get 1d data back from the buffer')
  call test(test_buffer_retrieval_2d, 'get 2d data back from the buffer')
  call test(test_buffer_retrieval_3d, 'get 3d data back from the buffer')
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
