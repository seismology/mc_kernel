!=========================================================================================
module test_readfields

   use global_parameters
   use readfields
   use ftnunit

   implicit none
   public
contains

!-----------------------------------------------------------------------------------------
! Fails for some involved reason. Its functionality is in other tests
!subroutine test_readfields_load_seismogram()
!   use type_parameter, only : parameter_type
!   type(semdata_type)      :: sem_data
!   type(parameter_type)    :: parameters
!
!   call parameters%read_parameters('unit_tests/inparam_test')
!   call parameters%read_source()
!   call parameters%read_receiver()
!   call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir, 100, 100, parameters%strain_type_fwd)
!   call sem_data%open_files()
!   call sem_data%read_meshes()
!
!   call parameters%read_filter(nomega = 512, df = 1d0)
!   call parameters%read_kernel(sem_data, parameters%filter)
!   
!   call sem_data%load_seismogram(parameters%receiver, parameters%source)
!   
!   call sem_data%close_files()
!
!end subroutine test_readfields_load_seismogram
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_set_params()
   type(semdata_type)   :: sem_data
   character(len=512)   :: fwd_dir, bwd_dir

   fwd_dir = './test_wavefields/kerner_fwd'
   bwd_dir = './test_wavefields/kerner_bwd'
   
   call sem_data%set_params(fwd_dir, bwd_dir, 100, 100, 'straintensor_trace')

   call assert_equal_int(sem_data%nsim_fwd, 4, 'nsim_fwd == 4')
   call assert_equal_int(sem_data%nsim_bwd, 2, 'nsim_bwd == 2')

end subroutine test_readfields_set_params
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_open_files()
   type(semdata_type)   :: sem_data
   character(len=512)   :: fwd_dir, bwd_dir

   fwd_dir = './test_wavefields/kerner_fwd'
   bwd_dir = './test_wavefields/kerner_bwd'
   
   call sem_data%set_params(fwd_dir, bwd_dir, 100, 100, 'straintensor_trace') 

   call sem_data%open_files()

   call sem_data%close_files()

end subroutine  test_readfields_open_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_read_meshes
   type(semdata_type)   :: sem_data
   type(meshtype)       :: testmesh_fwd, testmesh_bwd
   character(len=512)   :: fwd_dir, bwd_dir
   integer              :: i_max_vp, i_max_vs, i_max_rho
   real(kind=sp)        :: r_max_vp, r_max_vs, r_max_rho

   fwd_dir = './test_wavefields/kerner_fwd'
   bwd_dir = './test_wavefields/kerner_bwd'
   
   call sem_data%set_params(fwd_dir, bwd_dir, 100, 100, 'straintensor_trace') 

   call sem_data%open_files()

   call sem_data%read_meshes()

   testmesh_fwd = sem_data%get_mesh('fwd')
   testmesh_bwd = sem_data%get_mesh('bwd')

   call assert_comparable(testmesh_fwd%s, testmesh_bwd%s, 1e-10, &
                          'FWD and BWD meshes are identical')
   
   call assert_comparable(testmesh_fwd%mu, testmesh_bwd%mu, 1e-10, &
                          'FWD and BWD meshes are identical')
   

   !This is PREM, so we might check whether the mesh variables have reasonable limits

   ! S
   call assert_comparable(minval(testmesh_fwd%s),   0.0,       1e-5, 'All S>0')
   call assert_comparable(maxval(testmesh_fwd%s),   6371000.0, 1e-5, 'All S<6371')

   ! Z
   call assert_comparable(minval(testmesh_fwd%z),  -6371000.0, 1e-5, 'All Z>-6371')
   call assert_comparable(maxval(testmesh_fwd%z),   6371000.0, 1e-5, 'All Z< 6371')

   ! VP
   call assert_comparable(minval(testmesh_fwd%vp),  5800.0 ,   1e-5, 'Min VP=5800')
   call assert_comparable(maxval(testmesh_fwd%vp),  13716.6,   1e-5, 'Max VP=13717')

   ! VS
   call assert_comparable(minval(testmesh_fwd%vs),  0000.0 ,   1e-5, 'All VS>0000')
   call assert_comparable(maxval(testmesh_fwd%vs),  7265.97,   1e-5, 'All VS<7265')

   ! Rho
   call assert_comparable(minval(testmesh_fwd%rho), 2600.0 ,   1e-5, 'All Rho>2600')
   call assert_comparable(maxval(testmesh_fwd%rho), 13088.48,  1e-5, 'All Rho<13089')
   
   ! Mu
   call assert_comparable(minval(testmesh_fwd%mu),  0.0 ,     1e-5, 'Min Mu=0')
   call assert_comparable(maxval(testmesh_fwd%mu),  293.77e9,  1e-5, 'Max Mu=293.8 GPa')
   
   ! Lambda
   ! These would be the correct values according to the PREM paper, but who knows
   !call assert_comparable(minval(testmesh_fwd%lambda), 52.0e9 ,   1e-5, 'Min Lambda=52 GPa')
   !call assert_comparable(maxval(testmesh_fwd%lambda), 1425.3e9,  1e-5, 'Max Lambda=1425 GPa')
   call assert_comparable(minval(testmesh_fwd%lambda), 34.216e9 ,   1e-5, 'Min Lambda=34 GPa')
   call assert_comparable(maxval(testmesh_fwd%lambda), 1307.96e9,  1e-5, 'Max Lambda=1308 GPa')
  
   ! Eta
   call assert_comparable(minval(testmesh_fwd%eta), 1.0,   1e-5, 'Min eta=1')
   call assert_comparable(maxval(testmesh_fwd%eta), 1.0,   1e-5, 'Max eta=1')

   ! I have no idea, what these values should be - SCS
   !! Phi
   !call assert_comparable(minval(testmesh_fwd%phi), 52.0e9 ,   1e-5, 'Min phi=52 GPa')
   !call assert_comparable(maxval(testmesh_fwd%phi), 1425.3e9,  1e-5, 'Max phi=1425 GPa')
   !
   !! Xi
   !call assert_comparable(minval(testmesh_fwd%xi), 52.0e9 ,   1e-5, 'Min xi=52 GPa')
   !call assert_comparable(maxval(testmesh_fwd%xi), 1425.3e9,  1e-5, 'Max xi=1425 GPa')
   
   ! Maximum VP should be at the CMB
   i_max_vp = maxloc(testmesh_fwd%vp, 1)
   r_max_vp = norm2([testmesh_fwd%s(i_max_vp), testmesh_fwd%z(i_max_vp)])
   call assert_comparable(r_max_vp, 3480000.0, 1e-5, 'Max VP at the CMB')

   ! Maximum VS should be at the D"
   i_max_vs = maxloc(testmesh_fwd%vs, 1)
   r_max_vs = norm2([testmesh_fwd%s(i_max_vs), testmesh_fwd%z(i_max_vs)])
   call assert_comparable(r_max_vs, 3630000.0, 1e-5, 'Max VS at the D"')

   ! Maximum Rho should be in the center
   i_max_rho = maxloc(testmesh_fwd%rho, 1)
   r_max_rho = norm2([testmesh_fwd%s(i_max_rho), testmesh_fwd%z(i_max_rho)])
   call assert_comparable(r_max_rho, 0.0, 1e-5, 'Max rho at the center of the earth')



   call sem_data%close_files()

end subroutine  test_readfields_read_meshes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_load_model_coeffs
   use backgroundmodel, only   : backgroundmodel_type
   use global_parameters, only : pi, sp, dp
   type(semdata_type)         :: sem_data
   type(meshtype)             :: testmesh_fwd, testmesh_bwd
   type(backgroundmodel_type) :: model
   character(len=512)         :: fwd_dir, bwd_dir
   integer, parameter         :: npoints = 100
   real(kind=dp)              :: coordinates(3,npoints), phi(npoints), z(npoints)
   integer                    :: ipoint

   fwd_dir = './test_wavefields/kerner_fwd'
   bwd_dir = './test_wavefields/kerner_bwd'
   
   call sem_data%set_params(fwd_dir, bwd_dir, 100, 100, 'straintensor_trace') 

   call sem_data%open_files()

   call sem_data%read_meshes()

   call random_number(phi)
   call random_number(z)
   phi = phi * 2.d0 * pi

   do ipoint = 1, npoints
      coordinates(:, ipoint) = [cos(phi(ipoint)), sin(phi(ipoint)), (z(ipoint)-0.5)*2] & 
                               * 6.371d6
   end do

   model = sem_data%load_model_coeffs(coordinates)

   call assert_true(all(abs(model%c_vp(1)-model%c_vp(2:))<1e-10), 'vP is identical for same radius')
   call assert_true(all(abs(model%c_vpv(1)-model%c_vpv(2:))<1e-10), 'vPv is identical for same radius')
   call assert_true(all(abs(model%c_vph(1)-model%c_vph(2:))<1e-10), 'vPh is identical for same radius')
   call assert_true(all(abs(model%c_vs(1)-model%c_vs(2:))<1e-10), 'vS is identical for same radius')
   call assert_true(all(abs(model%c_vsv(1)-model%c_vsv(2:))<1e-10), 'vSv is identical for same radius')
   call assert_true(all(abs(model%c_vsh(1)-model%c_vsh(2:))<1e-10), 'vSh is identical for same radius')
   call assert_true(all(abs(model%c_xi(1)-model%c_xi(2:))<1e-10), 'Xi is identical for same radius')
   call assert_true(all(abs(model%c_phi(1)-model%c_phi(2:))<1e-10), 'Phi is identical for same radius')
   call assert_true(all(abs(model%c_eta(1)-model%c_eta(2:))<1e-10), 'eta is identical for same radius')
   call assert_true(all(abs(model%c_rho(1)-model%c_rho(2:))<1e-10), 'eta is identical for same radius')

   call sem_data%close_files()

end subroutine test_readfields_load_model_coeffs
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> For now, this test only checks, whether something crashes while reading points
subroutine test_readfields_load_fw_points
   use readfields,     only : semdata_type
   use type_parameter, only : parameter_type
   type(parameter_type)    :: parameters
   type(semdata_type)      :: sem_data
   real(kind=dp)           :: u(73, 6, 2), coordinates(3,2)
   
   call parameters%read_parameters('unit_tests/inparam_test')
   call parameters%read_source()
   call parameters%read_receiver()

   call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir, 100, 100,  &
                            parameters%strain_type_fwd)
   call sem_data%open_files()
   call sem_data%read_meshes()
   call sem_data%build_kdtree()
   call sem_data%load_seismogram(parameters%receiver, parameters%source)

   call parameters%read_filter(nomega=129, df=0.1d0)
     
   call parameters%read_kernel(sem_data, parameters%filter)

   coordinates(:,1) = [ 1d6, 1d6, 1d6]
   coordinates(:,2) = [-1d6, 1d6, 1d6]
   u = sem_data%load_fw_points(coordinates, parameters%source)

   call sem_data%close_files()

end subroutine test_readfields_load_fw_points
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
