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
module test_rotations

   use global_parameters
   use ftnunit
   use rotations

   implicit none
   public
contains

!-----------------------------------------------------------------------------------------
subroutine test_readfields_rotate

   integer, parameter :: npoints = 3
   real(kind=dp)      :: r_in(3,npoints)
   real(kind=dp)      :: phi_rot, theta_rot
   real(kind=dp)      :: s_out(npoints), z_out(npoints), phi_out(npoints)
 
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
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_rotate_straintensor()
   real(kind=dp)        :: straintensor(1,6)
   real(kind=sp)        :: straintensor_rot(1,6)
   real(kind=sp)        :: straintensor_rot_ref(1,6)
   real(kind=dp)        :: mij(6)  !Mrr Mtt Mpp Mrt Mrp Mtp
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
       straintensor_rot = real(rotate_straintensor(straintensor, phi, mij, isim), kind=sp) &
                        + straintensor_rot
   enddo
   straintensor_rot_ref(1,:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(straintensor_rot(1,:), straintensor_rot_ref(1,:), &
                                 1e-7, 'Rotation of explosion source, phi = 0')

   ! Arbitrary azimuth
   call random_number(phi)
   straintensor_rot = 0
   do isim = 1, 4
       straintensor_rot = real(rotate_straintensor(straintensor, phi, mij, isim), kind=sp) &
                        + straintensor_rot
   enddo
   straintensor_rot_ref(1,:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(straintensor_rot(1,:), straintensor_rot_ref(1,:), &
                                 1e-7, 'Rotation of explosion source, arbitrary phi')


end subroutine test_readfields_rotate_straintensor
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_rotate_straintensor_voigt()
   real(kind=dp)        :: straintensor(1,6)
   real(kind=sp)        :: straintensor_rot(1,6)
   real(kind=sp)        :: straintensor_rot_ref(1,6)
   real(kind=dp)        :: mij(6)  ! Mtt Mpp Mrr Mrp Mrt Mtp
   real(kind=dp)        :: phi
   integer              :: isim

   ! 'strain_dsus' 
   ! 'strain_dpup' 
   ! 'straintrace' - 'strain_dsus' - 'strain_dpup'
   ! 'strain_dzup'
   ! 'strain_dsuz'
   ! 'strain_dsup'

   ! Explosion source, pure diagonal strain
   straintensor(1,:) = [1, 1, 1, 0, 0, 0]
   mij = [1, 1, 1, 0, 0, 0] / 2.d0

   ! Azimuth zero
   phi = 0
   straintensor_rot = 0
   do isim = 1, 4
       straintensor_rot = real(rotate_straintensor_voigt(straintensor, phi, mij, isim), kind=sp) &
                        + straintensor_rot
   enddo
   straintensor_rot_ref(1,:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(straintensor_rot(1,:), straintensor_rot_ref(1,:), &
                                 1e-7, 'Rotation of explosion source, phi = 0')

   ! Arbitrary azimuth
   call random_number(phi)
   straintensor_rot = 0
   do isim = 1, 4
       straintensor_rot = real(rotate_straintensor_voigt(straintensor, phi, mij, isim), kind=sp) &
                        + straintensor_rot
   enddo
   straintensor_rot_ref(1,:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(straintensor_rot(1,:), straintensor_rot_ref(1,:), &
                                 1e-7, 'Rotation of explosion source, arbitrary phi')


end subroutine test_readfields_rotate_straintensor_voigt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_rotate_symm_tensor_voigt_src_to_xyz_1d
   real(kind=dp)        :: symm_tensor(6)
   real(kind=sp)        :: symm_tensor_rot(6)
   real(kind=sp)        :: symm_tensor_rot_ref(6)
   real(kind=dp)        :: phi

   ! explosion - Azimuth zero
   phi = 0
   symm_tensor(:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! explosion - Arbitrary azimuth
   call random_number(phi)
   symm_tensor(:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')

   ! swap component 1 and 2
   phi = 90 * deg2rad
   symm_tensor(:) = [1, 0, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [0, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')

   ! swap component 4 and 5
   phi = 90 * deg2rad
   symm_tensor(:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, -1, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')

   ! Arbitrary azimuth - inverse test
   call random_number(phi)
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor = rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi)
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, 4, 5, 6]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')

   ! DC 45
   phi = 45 * deg2rad
   symm_tensor(:) = [1, -1, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 0, 1]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 45')

   ! Mxx 45
   phi = 45 * deg2rad
   symm_tensor(:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [.5, .5, 0., 0., 0., .5]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Mxx, phi = 45')

   ! 180
   phi = 180 * deg2rad
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, -4, -5, 6]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '180')

end subroutine test_rotate_symm_tensor_voigt_src_to_xyz_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_rotate_symm_tensor_voigt_src_to_xyz_2d
   real(kind=dp)        :: symm_tensor(2,6)
   real(kind=sp)        :: symm_tensor_rot(2,6)
   real(kind=sp)        :: symm_tensor_rot_ref(6)
   real(kind=dp)        :: phi

   ! explosion - Azimuth zero
   phi = 0
   symm_tensor(1,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! explosion - Arbitrary azimuth
   call random_number(phi)
   symm_tensor(1,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')

   ! swap component 1 and 2
   phi = 90 * deg2rad
   symm_tensor(1,:) = [1, 0, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')

   ! swap component 4 and 5
   phi = 90 * deg2rad
   symm_tensor(1,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor(2,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, -1, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')

   ! Arbitrary azimuth - inverse test
   call random_number(phi)
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor = rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi, 2)
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, 4, 5, 6]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')

   ! DC 45
   phi = 45 * deg2rad
   symm_tensor(1,:) = [1, -1, 0, 0, 0, 0]
   symm_tensor(2,:) = [1, -1, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 0, 1]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 45')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 45')

   ! Mxx 45
   phi = 45 * deg2rad
   symm_tensor(1,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [.5, .5, 0., 0., 0., .5]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Mxx, phi = 45')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Mxx, phi = 45')

   ! 180
   phi = 180 * deg2rad
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, -4, -5, 6]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '180')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '180')

end subroutine test_rotate_symm_tensor_voigt_src_to_xyz_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_rotate_symm_tensor_voigt_xyz_to_src_1d
   real(kind=dp)        :: symm_tensor(6)
   real(kind=sp)        :: symm_tensor_rot(6)
   real(kind=sp)        :: symm_tensor_rot_ref(6)
   real(kind=dp)        :: phi

   ! explosion Azimuth zero
   phi = 0
   symm_tensor(:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! explosion Arbitrary azimuth
   call random_number(phi)
   symm_tensor(:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')

   ! swap component 1 and 2
   phi = 90 * deg2rad
   symm_tensor(:) = [1, 0, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [0, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')

   ! swap component 4 and 5
   phi = 90 * deg2rad
   symm_tensor(:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 1, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')

   ! Arbitrary azimuth - inverse test
   call random_number(phi)
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor = rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi)
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, 4, 5, 6]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')

   ! DC 45
   phi = 45 * deg2rad
   symm_tensor(:) = [1, -1, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 0, -1]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 45')

   ! Mxx 45
   phi = 45 * deg2rad
   symm_tensor(:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [.5, .5, 0., 0., 0., -.5]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Mxx, phi = 45')

   ! 180
   phi = 180 * deg2rad
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, -4, -5, 6]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '180')

end subroutine test_rotate_symm_tensor_voigt_xyz_to_src_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_rotate_symm_tensor_voigt_xyz_to_src_2d
   real(kind=dp)        :: symm_tensor(2,6)
   real(kind=sp)        :: symm_tensor_rot(2,6)
   real(kind=sp)        :: symm_tensor_rot_ref(6)
   real(kind=dp)        :: phi

   ! explosion Azimuth zero
   phi = 0
   symm_tensor(1,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! explosion Arbitrary azimuth
   call random_number(phi)
   symm_tensor(1,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')

   ! swap component 1 and 2
   phi = 90 * deg2rad
   symm_tensor(1,:) = [1, 0, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')

   ! swap component 4 and 5
   phi = 90 * deg2rad
   symm_tensor(1,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor(2,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 1, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')

   ! Arbitrary azimuth - inverse test
   call random_number(phi)
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor = rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi, 2)
   symm_tensor_rot = real(rotate_symm_tensor_voigt_src_to_xyz(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, 4, 5, 6]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')

   ! DC 45
   phi = 45 * deg2rad
   symm_tensor(1,:) = [1, -1, 0, 0, 0, 0]
   symm_tensor(2,:) = [1, -1, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 0, -1]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 45')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 45')

   ! Mxx 45
   phi = 45 * deg2rad
   symm_tensor(1,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [.5, .5, 0., 0., 0., -.5]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Mxx, phi = 45')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Mxx, phi = 45')

   ! 180
   phi = 180 * deg2rad
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_to_src(symm_tensor, phi, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, -4, -5, 6]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '180')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '180')

end subroutine test_rotate_symm_tensor_voigt_xyz_to_src_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_rotate_symm_tensor_voigt_xyz_src_to_xyz_earth_1d
   real(kind=dp)        :: symm_tensor(6)
   real(kind=sp)        :: symm_tensor_rot(6)
   real(kind=sp)        :: symm_tensor_rot_ref(6)
   real(kind=dp)        :: phi, theta

   ! explosion Azimuth zero
   phi = 0
   theta = 0
   symm_tensor(:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! explosion Arbitrary azimuth
   call random_number(phi)
   call random_number(theta)
   symm_tensor(:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')

   ! Arbitrary angles - inverse test
   call random_number(phi)
   call random_number(theta)
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor = rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta)
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, 4, 5, 6]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')

   ! swap component 1 and 2
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(:) = [1, 0, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [0, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')

   ! swap component 4 and 5
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, -1, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')

   ! DC
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, -1, 0] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 90, theta = 0')
   
   ! dipole
   phi = 45 * deg2rad
   theta = 0
   symm_tensor(:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 0, 0, 0, 1] / 2.

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Dipole, phi = 45, theta = 0')

   ! DC
   phi = 0
   theta = 90 * deg2rad
   symm_tensor(:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 0, 1] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')
   
   ! Dipole
   phi = 0
   theta = 45 * deg2rad
   symm_tensor(:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 0, 1, 0, -1, 0] / 2.

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')



   ! theta = 90: x > z, y > y, z > -x
   phi = 0
   theta = 90 * deg2rad
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [3, 2, 1, -6, -5, 4] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! theta = 90, phi = 90: x > -y, y > z, z > -x
   phi = 90 * deg2rad
   theta = 90 * deg2rad
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [2, 3, 1, -5, 6, -4] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')


end subroutine test_rotate_symm_tensor_voigt_xyz_src_to_xyz_earth_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_rotate_symm_tensor_voigt_xyz_src_to_xyz_earth_2d
   real(kind=dp)        :: symm_tensor(2,6)
   real(kind=sp)        :: symm_tensor_rot(2,6)
   real(kind=sp)        :: symm_tensor_rot_ref(6)
   real(kind=dp)        :: phi, theta

   ! explosion Azimuth zero
   phi = 0
   theta = 0
   symm_tensor(1,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! explosion Arbitrary azimuth
   call random_number(phi)
   call random_number(theta)
   symm_tensor(1,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')

   ! Arbitrary angles - inverse test
   call random_number(phi)
   call random_number(theta)
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor = rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2)
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, 4, 5, 6]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')

   ! swap component 1 and 2
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(1,:) = [1, 0, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')

   ! swap component 4 and 5
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(1,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor(2,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, -1, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')

   ! DC
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(1,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor(2,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, -1, 0] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 90, theta = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 90, theta = 0')
   
   ! dipole
   phi = 45 * deg2rad
   theta = 0
   symm_tensor(1,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 0, 0, 0, 1] / 2.

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Dipole, phi = 45, theta = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Dipole, phi = 45, theta = 0')

   ! DC
   phi = 0
   theta = 90 * deg2rad
   symm_tensor(1,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor(2,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 0, 1] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')
   
   ! Dipole
   phi = 0
   theta = 45 * deg2rad
   symm_tensor(1,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 0, 1, 0, -1, 0] / 2.

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')



   ! theta = 90: x > z, y > y, z > -x
   phi = 0
   theta = 90 * deg2rad
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [3, 2, 1, -6, -5, 4] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! theta = 90, phi = 90: x > -y, y > z, z > -x
   phi = 90 * deg2rad
   theta = 90 * deg2rad
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [2, 3, 1, -5, 6, -4] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')


end subroutine test_rotate_symm_tensor_voigt_xyz_src_to_xyz_earth_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_rotate_symm_tensor_voigt_xyz_earth_to_xyz_src_1d
   real(kind=dp)        :: symm_tensor(6)
   real(kind=sp)        :: symm_tensor_rot(6)
   real(kind=sp)        :: symm_tensor_rot_ref(6)
   real(kind=dp)        :: phi, theta

   ! explosion Azimuth zero
   phi = 0
   theta = 0
   symm_tensor(:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! explosion Arbitrary azimuth
   call random_number(phi)
   call random_number(theta)
   symm_tensor(:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')

   ! Arbitrary angles - inverse test
   call random_number(phi)
   call random_number(theta)
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor = rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta)
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, 4, 5, 6]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')

   ! swap component 1 and 2
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(:) = [1, 0, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [0, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')

   ! swap component 4 and 5
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 1, 0]

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')

   ! DC
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 1, 0] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 90, theta = 0')
   
   ! dipole
   phi = 45 * deg2rad
   theta = 0
   symm_tensor(:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 0, 0, 0, -1] / 2.

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Dipole, phi = 45, theta = 0')

   ! DC
   phi = 0
   theta = 90 * deg2rad
   symm_tensor(:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 0, -1] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')

   ! Dipole
   phi = 0
   theta = 45 * deg2rad
   symm_tensor(:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [1, 0, 1, 0, 1, 0] / 2.

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')

   
   
   ! theta = 90: x > -z, y > y, z > x
   phi = 0
   theta = 90 * deg2rad
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [3, 2, 1, 6, -5, -4] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '123456, phi = 0, theta = 90')

   ! theta = 90, phi = 90: x > -z, y > -x, z > y
   phi = 90 * deg2rad
   theta = 90 * deg2rad
   symm_tensor(:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta), kind=sp)
   symm_tensor_rot_ref(:) = [3, 1, 2, -6, -4, 5] 

   call assert_comparable_real1d(symm_tensor_rot(:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '123456, phi = 90, theta = 90')

end subroutine test_rotate_symm_tensor_voigt_xyz_earth_to_xyz_src_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_rotate_symm_tensor_voigt_xyz_earth_to_xyz_src_2d
   real(kind=dp)        :: symm_tensor(2,6)
   real(kind=sp)        :: symm_tensor_rot(2,6)
   real(kind=sp)        :: symm_tensor_rot_ref(6)
   real(kind=dp)        :: phi, theta

   ! explosion Azimuth zero
   phi = 0
   theta = 0
   symm_tensor(1,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = 0')

   ! explosion Arbitrary azimuth
   call random_number(phi)
   call random_number(theta)
   symm_tensor(1,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 1, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Rotation of isotropic tensor, phi = random')

   ! Arbitrary angles - inverse test
   call random_number(phi)
   call random_number(theta)
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor = rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2)
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_src_to_xyz_earth(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 2, 3, 4, 5, 6]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'inverse test, phi = random')

   ! swap component 1 and 2
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(1,:) = [1, 0, 1, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 1, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 1, 1, 0, 0, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 1 and 2, phi = 90')

   ! swap component 4 and 5
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(1,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor(2,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 1, 0]

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'swap 4 and 5, phi = 90')

   ! DC
   phi = 90 * deg2rad
   theta = 0
   symm_tensor(1,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor(2,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 1, 0] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 90, theta = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 90, theta = 0')
   
   ! dipole
   phi = 45 * deg2rad
   theta = 0
   symm_tensor(1,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 1, 0, 0, 0, -1] / 2.

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Dipole, phi = 45, theta = 0')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'Dipole, phi = 45, theta = 0')

   ! DC
   phi = 0
   theta = 90 * deg2rad
   symm_tensor(1,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor(2,:) = [0, 0, 0, 1, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [0, 0, 0, 0, 0, -1] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')

   ! Dipole
   phi = 0
   theta = 45 * deg2rad
   symm_tensor(1,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor(2,:) = [1, 0, 0, 0, 0, 0]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [1, 0, 1, 0, 1, 0] / 2.

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, 'DC, phi = 0, theta = 90')

   
   
   ! theta = 90: x > -z, y > y, z > x
   phi = 0
   theta = 90 * deg2rad
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [3, 2, 1, 6, -5, -4] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '123456, phi = 0, theta = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '123456, phi = 0, theta = 90')

   ! theta = 90, phi = 90: x > -z, y > -x, z > y
   phi = 90 * deg2rad
   theta = 90 * deg2rad
   symm_tensor(1,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor(2,:) = [1, 2, 3, 4, 5, 6]
   symm_tensor_rot = real(rotate_symm_tensor_voigt_xyz_earth_to_xyz_src(symm_tensor, phi, theta, 2), kind=sp)
   symm_tensor_rot_ref(:) = [3, 1, 2, -6, -4, 5] 

   call assert_comparable_real1d(symm_tensor_rot(1,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '123456, phi = 90, theta = 90')
   call assert_comparable_real1d(symm_tensor_rot(2,:) + 1, symm_tensor_rot_ref(:) + 1, &
                                 1e-7, '123456, phi = 90, theta = 90')

end subroutine test_rotate_symm_tensor_voigt_xyz_earth_to_xyz_src_2d
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
