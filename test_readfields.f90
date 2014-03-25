!=========================================================================================
module test_readfields

   use global_parameters
   use readfields
   use ftnunit
   implicit none
   public
contains

!------------------------------------------------------------------------------
subroutine test_readfields_load_seismogram()
   use type_parameter, only : parameter_type
   type(semdata_type)      :: sem_data
   type(parameter_type)    :: parameters

   call parameters%read_parameters('unit_tests/inparam_test')
   call parameters%read_source()
   call parameters%read_receiver()
   call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir, 100)
   call sem_data%open_files()
   call sem_data%read_meshes()
   call sem_data%load_seismogram(parameters%receiver, parameters%source)
   
   !print *, sem_data%veloseis
   
   
   call sem_data%close_files()


end subroutine test_readfields_load_seismogram
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_readfields_set_params()
   type(semdata_type)   :: sem_data
   character(len=512)   :: fwd_dir, bwd_dir

   fwd_dir = './wavefield/fwd'
   bwd_dir = './wavefield/bwd'
   
   call sem_data%set_params(fwd_dir, bwd_dir, 100) 

   call assert_equal_int(sem_data%nsim_fwd, 4, 'nsim_fwd == 4')
   call assert_equal_int(sem_data%nsim_bwd, 1, 'nsim_bwd == 1')

end subroutine test_readfields_set_params
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_readfields_open_files()
   type(semdata_type)   :: sem_data
   character(len=512)   :: fwd_dir, bwd_dir

   fwd_dir = './wavefield/fwd'
   bwd_dir = './wavefield/bwd'
   
   call sem_data%set_params(fwd_dir, bwd_dir, 100) 

   call sem_data%open_files()

   call sem_data%close_files()

end subroutine  test_readfields_open_files
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
subroutine test_readfields_rotate

   integer, parameter :: npoints = 3
   real(kind=dp)      :: r_in(3,npoints)
   real(kind=dp)      :: phi_rot, theta_rot
   real(kind=dp)      :: s_out(npoints), z_out(npoints), phi_out(npoints)
   integer            :: itheta, iphi
   character(len=16)  :: fmtstring
 
   r_in(:,1) = [1, 0, 0]
   r_in(:,2) = [0, 1, 0]
   r_in(:,3) = [0, 0, 1]
   phi_rot = 0

   ! Tests start here   
   r_in(:,1) = [1, 0, 0]
   r_in(:,2) = [0, 1, 0]
   r_in(:,3) = [0, 0, 1]
   phi_rot = 0
   theta_rot = 0
 
   call rotate_frame_rd(npoints, s_out, phi_out, z_out, r_in, phi_rot, theta_rot)

   call assert_comparable_real1d(1+real(z_out), 1+[0.0, 0.0, 1.0],              &
                                 1e-7, 'Z no rotation')
   call assert_comparable_real1d(1+real(s_out), 1+[1.0, 1.0, 0.0], 1e-7,        &
                                 'S no rotation')
   call assert_comparable_real1d(1+real(phi_out), 1+real([0.0d0, pi/2, 0.0d0]), &
                                 1e-7, 'Phi no rotation')
 
   phi_rot = 0
   theta_rot = 90*deg2rad
 
   call rotate_frame_rd(npoints, s_out, phi_out, z_out, r_in, phi_rot, theta_rot)

   call assert_comparable_real1d(real(1+z_out), 1+[1.0, 0.0, 0.0],   & 
                                 1e-6, 'Z pi/2 rotation')
   call assert_comparable_real1d(real(1+s_out), 1+[0.0, 1.0, 1.0],   & 
                                 1e-6, 'S pi/2 rotation')
   call assert_comparable_real1d(real(phi_out), real([pi, pi/2, pi]),& 
                                 1e-6, 'Phi pi/2 rotation')
 
   phi_rot = 0
   theta_rot = -90*deg2rad
 
   call rotate_frame_rd(npoints, s_out, phi_out, z_out, r_in, phi_rot, theta_rot)

   call assert_comparable_real1d(real(1+z_out), 1+[-1.0, 0.0, 0.0],          &
                                 1e-6, 'Z -pi/2 rotation')
   call assert_comparable_real1d(real(1+s_out), 1+[0.0, 1.0, 1.0],           &
                                 1e-6, 'S -pi/2 rotation')
   call assert_comparable_real1d(real(1+phi_out), 1+real([pi, pi/2, 0.0d0]), &
                                 1e-6, 'Phi -pi/2 rotation')

   phi_rot = 0
   theta_rot = 180*deg2rad
 
   call rotate_frame_rd(npoints, s_out, phi_out, z_out, r_in, phi_rot, theta_rot)

   call assert_comparable_real1d(real(1+z_out), 1+[0.0, 0.0, -1.0],             &
                                 1e-6, 'Z 2pi rotation')
   call assert_comparable_real1d(real(1+s_out), 1+[1.0, 1.0, .0],               &
                                 1e-6, 'S 2pi rotation')
   call assert_comparable_real1d(real(1+phi_out), 1+real([pi, pi/2, 0.0d0]),    &
                                 1e-6, 'Phi 2pi rotation')

   phi_rot = 0
   theta_rot = 45*deg2rad
 
   call rotate_frame_rd(npoints, s_out, phi_out, z_out, r_in, phi_rot, theta_rot)

   call assert_comparable_real1d(real(1+z_out), 1+[sqrt(0.5), 0.0, sqrt(0.5)],  &
                                 1e-6, 'Z pi/4 rotation')
   call assert_comparable_real1d(real(1+s_out), 1+[sqrt(0.5), 1.0, sqrt(0.5)],  &
                                 1e-6, 'S pi/4 rotation')
   call assert_comparable_real1d(real(1+phi_out), 1+real([0.0d0, pi/2, pi]),    &
                                 1e-6, 'Phi pi/4 rotation')

   phi_rot = -45 * deg2rad
   theta_rot = 0
 
   call rotate_frame_rd(npoints, s_out, phi_out, z_out, r_in, phi_rot, theta_rot)

   call assert_comparable_real1d(real(1+z_out),   1+[0.0, 0.0, 1.0],             &
                                 1e-7, 'Z phi=pi/4 rotation')
   call assert_comparable_real1d(real(1+s_out),   1+[1.0, 1.0, 0.0],             &
                                 1e-7, 'S phi=pi/4 rotation')
   call assert_comparable_real1d(real(1+phi_out), 1+real([pi/4, 3*pi/4, 0.0d0]), &
                                 1e-7, 'Phi phi=pi/4 rotation')

end subroutine
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_readfields_rotate_straintensor()
   real(kind=dp)        :: straintensor(1,6)
   real(kind=sp)        :: straintensor_rot(1,6)
   real(kind=sp)        :: straintensor_rot_ref(1,6)
   real(kind=dp)        :: mij(6)
   real(kind=dp)        :: phi
   integer              :: isim

   !! 'strain_dsus', 'strain_dsuz', 'strain_dpup', &
   !! 'strain_dsup', 'strain_dzup', 'straintrace']

   ! Explosion source, pure diagonal strain
   straintensor(1,:) = [1, 0, 1, 0, 0, 3]
   mij = [1, 1, 1, 0, 0, 0] / 2.d0

   ! Azimuth zero
   phi = 0
   straintensor_rot = 0
   do isim = 1, 4
       straintensor_rot = rotate_straintensor(straintensor, phi, mij, isim) &
                        + straintensor_rot
   enddo
   straintensor_rot_ref(1,:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(straintensor_rot(1,:), straintensor_rot_ref(1,:), &
                                 1e-7, 'Rotation of explosion source, phi = 0')

   ! Arbitrary azimuth
   call random_number(phi)
   straintensor_rot = 0
   do isim = 1, 4
       straintensor_rot = rotate_straintensor(straintensor, phi, mij, isim) &
                        + straintensor_rot
   enddo
   straintensor_rot_ref(1,:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(straintensor_rot(1,:), straintensor_rot_ref(1,:), &
                                 1e-7, 'Rotation of explosion source, arbitrary phi')


end subroutine test_readfields_rotate_straintensor
!------------------------------------------------------------------------------
end module
!-----------------------------------------------------------------------------------------
