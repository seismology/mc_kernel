!******************************************************************************
!
!    This file is part of:
!    MC Kernel: Calculating seismic sensitivity kernels on unstructured meshes
!    Copyright (C) 2016 Simon Staehler, Martin van Driel, Ludwig Auer
!
!    You can find the latest version of the software at:
!    <https://www.github.com/tomography/mckernel>
!
!    MC Kernel is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MC Kernel is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MC Kernel. If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

!=========================================================================================
module unit_tests

  use ftnunit, only: test
  use global_parameters, only: lu_out, verbose, set_master, set_master
  use commpi, only: ppinit, ppsplit
  use test_source
  use test_mc_integration
  use test_halton_sequence
  use test_lanczos
  use test_fft_type
  use test_polytopes 
  use test_voxel
  use test_inversion_mesh
  use test_buffer
  use test_filter
  use test_kernel
  use test_nc_routines
  use test_readfields
  use test_background_model
  use test_heterogeneities
  use test_rotations
  use test_finite_elem_mapping
  use test_spectral_basis
  use test_sem_derivatives
  use test_simple_routines
  use test_type_parameter
  use test_master_queue
  use test_master_helper
  use test_worktype

  implicit none

contains
!-----------------------------------------------------------------------------------------
subroutine test_all
  integer   :: oldverbose

  verbose = 1
  call set_master(.true.)
  testing = .true.

  call init_output()

  call ppinit()
  call ppsplit(1)

  ! test source routines
  !write(6,'(/,a)') 'TEST SOURCE MODULE'
  !call test(test_resample_stf, 'resample stf')
  !call test(test_fft_stf, 'fft stf')

  ! test simple routines
  write(6,'(/,a)') 'TEST SIMPLE ROUTINES MODULE'
  call test(test_mult2d_1d, 'Multiply 2D array with 1D array')
  call test(test_mult3d_1d, 'Multiply 3D array with 1D array')
  call test(test_mult2d_1d_complex, 'Multiply complex 2D array with 1D array')
  call test(test_mult3d_1d_complex, 'Multiply complex 3D array with 1D array')
  call test(test_to_lower, 'Transform string to lowercase')
  call test(test_lowtrim,  'Transform string to lowercase and trim')
  call test(test_absreldiff_sp, 'Calculate absolute relative difference (single precision)   ')
  call test(test_absreldiff_dp, 'Calculate absolute relative difference (double precision)   ')
  call test(test_absreldiff_qp, 'Calculate absolute relative difference (quadruple precision)')
  call test(test_cross_sp, 'Cross product (single precision)')
  call test(test_cross_dp, 'Cross product (double precision)')
  call test(test_cross_qp, 'Cross product (quadruple precision)')
  oldverbose = verbose; verbose = 0 !These tests produce too much output otherwise
  call test(test_checklim_1d_int, 'Check limits 1D (integer)')
  call test(test_checklim_2d_int, 'Check limits 2D (integer)')
  call test(test_checklim_3d_int, 'Check limits 3D (integer)')
  call test(test_checklim_1d,     'Check limits 1D (float)')
  call test(test_checklim_2d,     'Check limits 2D (float)')
  call test(test_checklim_3d,     'Check limits 3D (float)')
  call test(test_cumsum_trapezoidal_1d, 'Cumulative sum (integration), 1D')
  call test(test_cumsum_trapezoidal_2d, 'Cumulative sum (integration), 2D')
  call test(test_cumsum_trapezoidal_3d, 'Cumulative sum (integration), 3D')
  verbose = oldverbose

  ! test sem derivatives
  write(6,'(/,a)') 'TEST SEM DERIVATIVE MODULE'
  call test(test_strain_monopole, 'SEM strain monopole')
  call test(test_strain_monopole_td, 'SEM strain monopole time dep')
  call test(test_strain_dipole, 'SEM strain dipole')
  call test(test_strain_dipole_td, 'SEM strain dipole time dep')
  call test(test_strain_quadpole, 'SEM strain quadpole')
  call test(test_strain_quadpole_td, 'SEM strain quadpole time dep')
  call test(test_f_over_s, 'SEM f_over_s')
  call test(test_f_over_s_td, 'SEM f_over_s time dep')
  call test(test_dsdf, 'SEM dsdf')
  call test(test_td_dsdf, 'SEM dsdfi time dependent')
  call test(test_gradient, 'SEM gradient')
  call test(test_gradient2, 'SEM gradient, npol=4')
  call test(test_td_gradient, 'time dependent SEM gradient')

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

  ! test_background_models
  write(6,'(/,a)') 'TEST BACKGROUND_MODELS MODULE'
  call test(test_background_model_combine, 'Combine model parameters')
  call test(test_background_model_weight,  'Weight model parameters')
  call test(test_background_model_get_parameter_names, 'Get model parameter names')

  ! test_heterogeneities
  write(6,'(/,a)') 'TEST HETEROGENEITIES MODULE'
  call test(test_load_het_rtpv, 'Load RTPV file')

  ! test_nc_routines
  write(6,'(/,a)') 'TEST NC_ROUTINES MODULE'
  call test_nc_create_testfile()
  call test(test_nc_create_file, 'Create NetCDF file')
  call test(test_nc_create_group, 'Create Group')
  call test(test_nc_open_for_read, 'Open NetCDF file for reading')
  call test(test_nc_open_for_write, 'Open NetCDF file for writing')
  call test(test_nc_getvar_1d_dble, 'Read 1D Double by name')
  call test(test_nc_getvar_1d_float, 'Read 1D Float by name')
  call test(test_nc_getvar_2d_float, 'Read 2D Float by name')
  call test(test_nc_getvar_3d_float, 'Read 3D Float by name')
  call test(test_nc_getvar_1d_int, 'Read 1D Integer by name')
  call test(test_nc_getvar_2d_int, 'Read 2D Integer by name')
  call test(test_nc_getvar_3d_int, 'Read 3D Integer by name')
  call test(test_nc_putvar_1d, 'Write 1D Float by name')
  call test(test_nc_putvar_2d, 'Write 2D Float by name')
  call test(test_nc_putvar_3d, 'Write 3D Float by name')
  call test(test_nc_putvar_1d_into_nd, 'Write 1D slices into 3D variable by name')

  ! test_readfields
  write(6,'(/,a)') 'TEST READFIELDS MODULE'
  call test(test_readfields_set_params, 'Set SEM file params')
  call test(test_readfields_open_files, 'Open SEM file')
  call test(test_readfields_load_straintrace_merged, 'Read straintrace from old and merged db')
  call test(test_readfields_load_strain_merged, 'Read strain from old and merged db')
  call test(test_readfields_load_model_coeffs, 'Read Model coefficients from SEM mesh')
  call test(test_load_seismograms_rdbm_Z, 'Read seismograms (Z-component)')
  call test(test_load_seismograms_rdbm_R, 'Read seismograms (R-component)')
  call test(test_load_seismograms_rdbm_T, 'Read seismograms (T-component)')
  call test(test_load_seismograms_rdbm_merged_Z, 'Read seismograms from merged db (Z-component)')
  call test(test_load_seismograms_rdbm_merged_R, 'Read seismograms from merged db (R-component)')
  call test(test_load_seismograms_rdbm_merged_T, 'Read seismograms from merged db (T-component)')
  call test(test_dampen_field, 'Dampen fields around source')
  call test(test_get_chunk_bounds, 'Get_chunk_bounds')

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
   
  ! test_mc_integration
  write(6,'(/,a)') 'TEST MC_INTEGRATION MODULE'
  call test(test_mc_meanandvariance, 'MC mean and variance')
  call test(test_mc_unit_hexagon, 'MC unit hexagon')
  call test(test_mc_sphere_in_tetrahedron, 'MC sphere in tetrahedron')

  ! test_halton_sequence
  write(6,'(/,a)') 'TEST HALTON SEQUENCE MODULE'
  call test(test_init_halton, 'Init_Halton sequence')
  call test(test_get_halton,  'Get_Halton sequence')

  ! test_lanczos
  write(6,'(/,a)') 'TEST LANCZOS RESAMPLING MODULE'
  call test(test_lanczos_resample, 'Lanczos resampling')

  ! test_fft_type
  write(6,'(/,a)') 'TEST FFT MODULE'
  call test(test_fft_dirac, 'FFT_dirac')
  call test(test_fft_sine, 'FFT_sine')
  call test(test_fft_inverse, 'FFT_inverse')
  call test(test_fft_inverse_md, 'FFT_inverse_multiple_traces')
  call test(test_fft_inverse_taz, 'FFT_inverse_taperandzeropad')
  call test(test_fft_convolve, 'FFT_convolve')
  call test(test_fft_taperandzeropad, 'FFT_taperandzeropad')
  call test(test_fft_parseval, 'FFT_Parseval_theorem')

  ! test filter
  write(6,'(/,a)') 'TEST FILTER MODULE'
  call test(test_filter_ident, 'Test Identical filter')
  call test(test_filter_gabor_response, 'Test Gabor filter')
  call test(test_filter_butterworth_lp_response, 'Test Butterworth LP filter')
  call test(test_filter_butterworth_hp_response, 'Test Butterworth HP filter')
  call test(test_filter_butterworth_bp_response, 'Test Butterworth BP filter')
  call test(test_filter_timeshift, 'Test time shifting routine')
  
  ! test kernel
  write(6,'(/,a)') 'TEST KERNEL MODULE'
  call test(test_kernel_init, 'Test Kernel initialization')
  call test(test_integrate_parseval_real, 'Test Parseval integration (real numbers)')
  call test(test_integrate_parseval_complex, 'Test Parseval integration (complex numbers)')
  call test(test_kernel_cut_timewindow, 'Test Time window cutting')
  call test(test_tabulate_kernel, 'Test base kernel tabulation')

  ! test type_parameter
  write(6,'(/,a)') 'TEST TYPE_PARAMETER MODULE'
  call test(test_parameter_reading, 'Reading in all parameters')

  ! test_polytopes
  write(6,'(/,a)') 'TEST POLYTOPES MODULE'
  call test(test_generate_random_point_poly_3, 'Random points in Triangle')
  call test(test_generate_random_point_poly_4, 'Random points in Quadrilateral')
  call test(test_generate_random_point_tet, 'Random points in Tetrahedra')
  call test(test_generate_random_point_triangle_space, 'Random points in Triangle in Space')
  call test(test_generate_random_point_triangle, 'Random points in reference triangle')
  call test(test_determinant, 'Matrix determinant')
  call test(test_tetra_volume_3d, 'Tetrahedron volume')
  call test(test_get_volume_poly_3, 'Triangle area')
  call test(test_get_volume_poly_4, 'Quadrilateral area')
  call test(test_point_in_triangle_3d, 'Test for test whether Point in triangle')
  call test(test_generate_random_point_poly_3_quasi, 'Quasi-Random points in Triangle')
  call test(test_generate_random_point_poly_4_quasi, 'Quasi-Random points in Quadrilateral')
  call test(test_generate_random_point_tet_quasi, 'Quasi-Random points in Tetrahedra')
  call test(test_generate_random_point_triangle_space_quasi, 'Quasi-Random points in Triangle in Space')
  call test(test_generate_random_point_triangle_quasi, 'Quasi-Random points in reference triangle')

  ! test_voxel
  write(6,'(/,a)') 'TEST VOXEL MODULE'
  call test(test_generate_random_points_vox, 'Random points in Voxel')
  call test(test_get_volume_vox, 'Volume of voxel')
  call test(test_coordinate_transformations_vox, 'Voxel related coordinate transformations')

  ! test_inversion_mesh
  write(6,'(/,a)') 'TEST INVERSION MESH MODULE'
  call test(test_mesh_read, 'reading tetrahedral mesh')
  ! Explicit mesh dumping does not exist in NetCDF output anymore
  !call test(test_mesh_dump, 'reading/dumping tetrahedral mesh')
  !call test(test_mesh_dump2, 'reading/dumping tetrahedral mesh from abaqus')
  !call test(test_mesh_dump3, 'reading/dumping tetrahedral mesh from abaqus with multiple element blocks')
  !call test(test_mesh_sort, 'sorting a mesh and dumping to xdmf')
  call test(test_append_variable, 'Append data to variable type')
  call test(test_init_node_data, 'Initialize node data')
  call test(test_init_cell_data, 'Initialize cell data')
  call test(test_init_mixed_data, 'Initialize mixed data')
  call test(test_set_node_data_and_dump, 'Set node data')
  call test(test_set_cell_data_and_dump, 'Set cell data')
  call test(test_set_cell_time_data_and_dump, 'Set cell data on time dimension')
  call test(test_set_node_time_data_and_dump, 'Set node data on time dimension')
  call test(test_set_mixed_data_and_dump, 'Set mixed data')
  call test(test_set_mixed_data_and_dump_into_two_files, 'Dump mixed data into two files')
  call test(test_set_mixed_data, 'Set mixed data with entry names')
  call test(test_valence, 'computation of valence')
  call test(test_get_connected_elements, 'get connected elements')
  call test(test_initialize_mesh, 'initialize mesh')
  call test(test_point_in_element_triangle_mesh, 'point_in_element for triangle mesh')
  call test(test_point_in_element_tetrahedral_mesh, 'point_in_element for tetrahedral mesh')
  call test(test_random_points_triangle_mesh, &
            'generate random numbers on triangular mesh')
  call test(test_weight, 'weight function (hat functions)')
  !call test(test_weight_issue30, 'weight function (specific case of issue 30)')
  !!call test(test_weight_large, 'weight function on a very large mesh (hat functions)')
  call test(test_integration_in_tetrahedron, 'MC-integrate in tetrahedral mesh element')

  ! test_buffer
  write(6,'(/,a)') 'TEST BUFFER MODULE'
  call test(test_buffer_storage_1d_sp, 'put 1d data into the buffer (single)')
  call test(test_buffer_storage_2d_sp, 'put 2d data into the buffer (single)')
  call test(test_buffer_storage_3d_sp, 'put 3d data into the buffer (single)')
  call test(test_buffer_storage_4d_sp, 'put 4d data into the buffer (single)')
  call test(test_buffer_storage_5d_sp, 'put 5d data into the buffer (single)')
  call test(test_buffer_retrieval_1d_sp, 'get 1d data back from the buffer (single)')
  call test(test_buffer_retrieval_2d_sp, 'get 2d data back from the buffer (single)')
  call test(test_buffer_retrieval_3d_sp, 'get 3d data back from the buffer (single)')
  call test(test_buffer_retrieval_4d_sp, 'get 4d data back from the buffer (single)')
  call test(test_buffer_retrieval_5d_sp, 'get 5d data back from the buffer (single)')
  call test(test_buffer_overwrite_sp, 'buffer gets overwritten after time')
  ! double precision buffer tests
  call test(test_buffer_storage_1d_dp, 'put 1d data into the buffer (double)')
  call test(test_buffer_storage_2d_dp, 'put 2d data into the buffer (double)')
  call test(test_buffer_storage_3d_dp, 'put 3d data into the buffer (double)')
  call test(test_buffer_storage_4d_dp, 'put 4d data into the buffer (double)')
  call test(test_buffer_storage_5d_dp, 'put 5d data into the buffer (double)')
  call test(test_buffer_retrieval_1d_dp, 'get 1d data back from the buffer (double)')
  call test(test_buffer_retrieval_2d_dp, 'get 2d data back from the buffer (double)')
  call test(test_buffer_retrieval_3d_dp, 'get 3d data back from the buffer (double)')
  call test(test_buffer_retrieval_4d_dp, 'get 4d data back from the buffer (double)')
  call test(test_buffer_retrieval_5d_dp, 'get 5d data back from the buffer (double)')
  call test(test_buffer_overwrite_dp, 'buffer gets overwritten after time')

  ! Test work type
  call test(test_init_work_type, 'work type initialization')

  ! Test master queue
  ! @TODO: Test does not work, should be checked
  !write(6,'(/,a)') 'TEST MASTER QUEUE'
  !call test(test_master_all, 'Master queue init_queue and finalize')
  
  ! Test master helper
  write(6,'(/,a)') 'TEST MASTER HELPER'
  call test(test_create_tasks, 'Create tasks from ''completed'' array')

  call finish_output()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_output()
  character(len=11) :: fnam
  integer           :: lu_out_local

  fnam = 'OUTPUT_test'
  open(newunit=lu_out_local, file=fnam, status='unknown', position='append')
  call set_lu_out(lu_out_local)
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
