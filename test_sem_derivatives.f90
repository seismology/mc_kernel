!=========================================================================================
module test_sem_derivatives

  use global_parameters
  use finite_elem_mapping
  use sem_derivatives
  use spectral_basis
  use ftnunit
  implicit none
  public

contains

!-----------------------------------------------------------------------------------------
subroutine test_strain_monopole()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol, nsamp, ipol, jpol
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: G1(:,:), G1T(:,:)
  real(kind=dp), allocatable   :: u(:,:,:), strain(:,:,:), strain_ref(:,:,:)

  real(kind=dp), allocatable   :: gll_points(:), glj_points(:)
  real(kind=dp), allocatable   :: s(:,:), z(:,:), sz(:,:,:)
  logical                      :: axial
  real(kind=dp)                :: phi, r1, r2

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1
  axial = .true.

  npol = 4

  G1 = def_lagrange_derivs_glj(npol)
  G2 = def_lagrange_derivs_gll(npol)
  allocate(G1T(0:npol,0:npol))
  allocate(G2T(0:npol,0:npol))
  G1T = transpose(G1)
  G2T = transpose(G2)

  allocate(gll_points(0:npol))
  allocate(glj_points(0:npol))
  gll_points(0:npol) = zelegl(npol)
  glj_points(0:npol) = zemngl2(npol)

  allocate(s(0:npol,0:npol))
  allocate(z(0:npol,0:npol))
  allocate(sz(0:npol,0:npol,2))

  allocate(u(0:npol,0:npol,3))
  allocate(strain(0:npol,0:npol,6))
  allocate(strain_ref(0:npol,0:npol,6))

  nodes(1,:) = [1,2]
  nodes(2,:) = [2,1]
  nodes(3,:) = [3,2]
  nodes(4,:) = [2,3]

  ! linear element
  element_type = 1
  axial = .false.

  do ipol=0, npol 
     do jpol=0, npol 
        sz(ipol, jpol,:) = mapping(glj_points(ipol), gll_points(jpol), &
                                   nodes, element_type)
        s = sz(:,:,1)
        z = sz(:,:,2)
     enddo
  enddo

  u(:,:,1) = s**2 * z
  u(:,:,2) = 0
  u(:,:,3) = s * z**2
  
  strain_ref(:,:,1) = 2 * s * z
  strain_ref(:,:,2) = s * z
  strain_ref(:,:,3) = 2 * s * z
  strain_ref(:,:,4) = 0
  strain_ref(:,:,5) = (s**2 + z**2) / 2
  strain_ref(:,:,6) = 0

  strain = strain_monopole(u, G2, G1T, glj_points, gll_points, npol, nodes, &
                           element_type, axial)
  
  call assert_comparable_real1d(real(reshape(strain, (/(npol+1)**2 * 6/))), &
                                real(reshape(strain_ref, (/(npol+1)**2 * 6/))), &
                                1e-7, 'monopole strain non-axial, linear element')

  phi = 0.1d0
  r1 = 10d0
  r2 = 11d0

  nodes(1,:) = [0d0, r1]
  nodes(2,:) = [r1 * dsin(phi), r1 * dcos(phi)]
  nodes(3,:) = [r2 * dsin(phi), r2 * dcos(phi)]
  nodes(4,:) = [0d0, r2]

  ! speroidal element
  element_type = 0
  axial = .true.

  do ipol=0, npol 
     do jpol=0, npol 
        sz(ipol, jpol,:) = mapping(glj_points(ipol), gll_points(jpol), &
                                   nodes, element_type)
        s = sz(:,:,1)
        z = sz(:,:,2)
     enddo
  enddo

  u(:,:,1) = s**2 * z
  u(:,:,2) = 0
  u(:,:,3) = s * z**2
  
  strain_ref(:,:,1) = 2 * s * z
  strain_ref(:,:,2) = s * z
  strain_ref(:,:,3) = 2 * s * z
  strain_ref(:,:,4) = 0
  strain_ref(:,:,5) = (s**2 + z**2) / 2
  strain_ref(:,:,6) = 0

  strain = strain_monopole(u, G2, G1T, glj_points, gll_points, npol, nodes, &
                           element_type, axial)
  
  call assert_comparable_real1d(10 + real(reshape(strain, (/(npol+1)**2 * 6/))), &
                                10 + real(reshape(strain_ref, (/(npol+1)**2 * 6/))), &
                                1e-5, 'monopole strain axial, spheroidal element')

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_strain_dipole()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol, nsamp, ipol, jpol
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: G1(:,:), G1T(:,:)
  real(kind=dp), allocatable   :: u(:,:,:), strain(:,:,:), strain_ref(:,:,:)

  real(kind=dp), allocatable   :: gll_points(:), glj_points(:)
  real(kind=dp), allocatable   :: s(:,:), z(:,:), sz(:,:,:)
  logical                      :: axial
  real(kind=dp)                :: phi1, phi2, r1, r2

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1
  axial = .true.

  npol = 4

  G1 = def_lagrange_derivs_glj(npol)
  G2 = def_lagrange_derivs_gll(npol)
  allocate(G1T(0:npol,0:npol))
  allocate(G2T(0:npol,0:npol))
  G1T = transpose(G1)
  G2T = transpose(G2)

  allocate(gll_points(0:npol))
  allocate(glj_points(0:npol))
  gll_points(0:npol) = zelegl(npol)
  glj_points(0:npol) = zemngl2(npol)

  allocate(s(0:npol,0:npol))
  allocate(z(0:npol,0:npol))
  allocate(sz(0:npol,0:npol,2))

  allocate(u(0:npol,0:npol,3))
  allocate(strain(0:npol,0:npol,6))
  allocate(strain_ref(0:npol,0:npol,6))

  nodes(1,:) = [1,2]
  nodes(2,:) = [2,1]
  nodes(3,:) = [3,2]
  nodes(4,:) = [2,3]

  ! linear element
  element_type = 1
  axial = .false.

  do ipol=0, npol 
     do jpol=0, npol 
        sz(ipol, jpol,:) = mapping(glj_points(ipol), gll_points(jpol), &
                                   nodes, element_type)
        s = sz(:,:,1)
        z = sz(:,:,2)
     enddo
  enddo

  u(:,:,1) = (s**2 * z + s**2 + z**2) / 2d0
  u(:,:,2) = (s**2 * z - s**2 - z**2) / 2d0
  u(:,:,3) = s * z**2
  
  strain_ref(:,:,1) = 2 * s * z
  strain_ref(:,:,2) = s * z - s - z**2 / s
  strain_ref(:,:,3) = 2 * s * z
  strain_ref(:,:,4) = - (z**2 + 2 * z) / 2
  strain_ref(:,:,5) = (s**2 + z**2) / 2
  strain_ref(:,:,6) = (z**2 / s - s * z - s) / 2

  strain = strain_dipole(u, G2, G1T, glj_points, gll_points, npol, nodes, &
                           element_type, axial)

  call assert_comparable_real1d(real(reshape(strain, (/(npol+1)**2 * 6/))), &
                                real(reshape(strain_ref, (/(npol+1)**2 * 6/))), &
                                1e-7, 'dipole strain non-axial, linear element')

  phi1 = 0.1d0
  phi2 = 0.2d0
  r1 = 10d0
  r2 = 11d0

  nodes(1,:) = [r1 * dsin(phi1), r1 * dcos(phi1)]
  nodes(2,:) = [r1 * dsin(phi2), r1 * dcos(phi2)]
  nodes(3,:) = [r2 * dsin(phi2), r2 * dcos(phi2)]
  nodes(4,:) = [r2 * dsin(phi1), r2 * dcos(phi1)]

  ! speroidal element
  element_type = 0
  axial = .false.

  do ipol=0, npol 
     do jpol=0, npol 
        sz(ipol, jpol,:) = mapping(glj_points(ipol), gll_points(jpol), &
                                   nodes, element_type)
        s = sz(:,:,1)
        z = sz(:,:,2)
     enddo
  enddo

  ! +- basis!
  u(:,:,1) = (s**2 * z + s**2 + z**2) / 2d0
  u(:,:,2) = (s**2 * z - s**2 - z**2) / 2d0
  u(:,:,3) = s * z**2
  
  strain_ref(:,:,1) = 2 * s * z
  strain_ref(:,:,2) = s * z - s - z**2 / s
  strain_ref(:,:,3) = 2 * s * z
  strain_ref(:,:,4) = - (z**2 + 2 * z) / 2
  strain_ref(:,:,5) = (s**2 + z**2) / 2
  strain_ref(:,:,6) = (z**2 / s - s * z - s) / 2

  strain = strain_dipole(u, G2, G1T, glj_points, gll_points, npol, nodes, &
                           element_type, axial)

  call assert_comparable_real1d(10 + real(reshape(strain, (/(npol+1)**2 * 6/))), &
                                10 + real(reshape(strain_ref, (/(npol+1)**2 * 6/))), &
                                1e-5, 'dipole strain axial, spheroidal element')

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_strain_quadpole()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol, nsamp, ipol, jpol
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: G1(:,:), G1T(:,:)
  real(kind=dp), allocatable   :: u(:,:,:), strain(:,:,:), strain_ref(:,:,:)

  real(kind=dp), allocatable   :: gll_points(:), glj_points(:)
  real(kind=dp), allocatable   :: s(:,:), z(:,:), sz(:,:,:)
  logical                      :: axial
  real(kind=dp)                :: phi1, phi2, r1, r2

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1
  axial = .true.

  npol = 4

  G1 = def_lagrange_derivs_glj(npol)
  G2 = def_lagrange_derivs_gll(npol)
  allocate(G1T(0:npol,0:npol))
  allocate(G2T(0:npol,0:npol))
  G1T = transpose(G1)
  G2T = transpose(G2)

  allocate(gll_points(0:npol))
  allocate(glj_points(0:npol))
  gll_points(0:npol) = zelegl(npol)
  glj_points(0:npol) = zemngl2(npol)

  allocate(s(0:npol,0:npol))
  allocate(z(0:npol,0:npol))
  allocate(sz(0:npol,0:npol,2))

  allocate(u(0:npol,0:npol,3))
  allocate(strain(0:npol,0:npol,6))
  allocate(strain_ref(0:npol,0:npol,6))

  nodes(1,:) = [1,2]
  nodes(2,:) = [2,1]
  nodes(3,:) = [3,2]
  nodes(4,:) = [2,3]

  ! linear element
  element_type = 1
  axial = .false.

  do ipol=0, npol 
     do jpol=0, npol 
        sz(ipol, jpol,:) = mapping(gll_points(ipol), gll_points(jpol), &
                                   nodes, element_type)
        s = sz(:,:,1)
        z = sz(:,:,2)
     enddo
  enddo

  u(:,:,1) = s**2 * z
  u(:,:,2) = s**2 + z**2
  u(:,:,3) = s * z**2
  
  strain_ref(:,:,1) = 2 * s * z
  strain_ref(:,:,2) = s * z - 2 * s - 2 * z**2 / s
  strain_ref(:,:,3) = 2 * s * z
  strain_ref(:,:,4) = - z**2 - z
  strain_ref(:,:,5) = (s**2 + z**2) / 2
  strain_ref(:,:,6) = (z**2 / s - 2 * s * z - s) / 2

  strain = strain_quadpole(u, G2, G2T, gll_points, gll_points, npol, nodes, &
                           element_type, axial)

  call assert_comparable_real1d(real(reshape(strain, (/(npol+1)**2 * 6/))), &
                                real(reshape(strain_ref, (/(npol+1)**2 * 6/))), &
                                1e-7, 'quadpole strain non-axial, linear element')

  phi1 = 0.1d0
  phi2 = 0.2d0
  r1 = 10d0
  r2 = 11d0

  nodes(1,:) = [r1 * dsin(phi1), r1 * dcos(phi1)]
  nodes(2,:) = [r1 * dsin(phi2), r1 * dcos(phi2)]
  nodes(3,:) = [r2 * dsin(phi2), r2 * dcos(phi2)]
  nodes(4,:) = [r2 * dsin(phi1), r2 * dcos(phi1)]

  ! speroidal element
  element_type = 0
  axial = .false.

  do ipol=0, npol 
     do jpol=0, npol 
        sz(ipol, jpol,:) = mapping(glj_points(ipol), gll_points(jpol), &
                                   nodes, element_type)
        s = sz(:,:,1)
        z = sz(:,:,2)
     enddo
  enddo

  u(:,:,1) = s**2 * z
  u(:,:,2) = s**2 + z**2
  u(:,:,3) = s * z**2
  
  strain_ref(:,:,1) = 2 * s * z
  strain_ref(:,:,2) = s * z - 2 * s - 2 * z**2 / s
  strain_ref(:,:,3) = 2 * s * z
  strain_ref(:,:,4) = - z**2 - z
  strain_ref(:,:,5) = (s**2 + z**2) / 2
  strain_ref(:,:,6) = (z**2 / s - 2 * s * z - s) / 2

  strain = strain_quadpole(u, G2, G1T, glj_points, gll_points, npol, nodes, &
                           element_type, axial)

  call assert_comparable_real1d(10 + real(reshape(strain, (/(npol+1)**2 * 6/))), &
                                10 + real(reshape(strain_ref, (/(npol+1)**2 * 6/))), &
                                1e-5, 'quadpole strain axial, spheroidal element')

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_f_over_s_td()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol, nsamp, jpol
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: u(:,:,:), u_over_s(:,:,:), u_over_s_ref(:,:,:)

  real(kind=dp), allocatable   :: gll_points(:)
  logical                      :: axial

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1
  axial = .true.

  npol = 4
  nsamp = 2

  G2 = def_lagrange_derivs_gll(npol)
  allocate(G2T(0:npol,0:npol))
  G2T = transpose(G2)
  gll_points = zelegl(npol)

  allocate(u(1:nsamp,0:npol,0:npol))
  allocate(u_over_s(1:nsamp,0:npol,0:npol))
  allocate(u_over_s_ref(1:nsamp,0:npol,0:npol))

  do jpol=0, npol 
     u(1,:,jpol) = (gll_points + 1) / 2d0
  enddo
  u_over_s_ref(1,0:npol,0:npol) = 1

  do jpol=0, npol 
     u(2,:,jpol) = (gll_points + 1)
  enddo
  u_over_s_ref(2,0:npol,0:npol) = 2

  u_over_s = f_over_s(u, G2, G2T, gll_points, gll_points, npol, nsamp, nodes, &
                      element_type, axial)
  
  call assert_comparable_real1d(1 + real(reshape(u_over_s, (/(npol+1)**2/))), &
                                1 + real(reshape(u_over_s_ref, (/(npol+1)**2/))), &
                                1e-7, 'f over s axial')

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_f_over_s()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol, jpol
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: u(:,:), u_over_s(:,:), u_over_s_ref(:,:)

  real(kind=dp), allocatable   :: gll_points(:)
  logical                      :: axial

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1
  axial = .true.

  npol = 4

  G2 = def_lagrange_derivs_gll(npol)
  allocate(G2T(0:npol,0:npol))
  G2T = transpose(G2)
  gll_points = zelegl(npol)

  allocate(u(0:npol,0:npol))
  allocate(u_over_s(0:npol,0:npol))
  allocate(u_over_s_ref(0:npol,0:npol))

  do jpol=0, npol 
     u(:,jpol) = (gll_points + 1) / 2d0
  enddo
  u_over_s_ref(0:npol,0:npol) = 1

  u_over_s(0:npol,0:npol) = f_over_s(u, G2, G2T, gll_points, gll_points, npol, nodes, &
                                     element_type, axial)
  
  call assert_comparable_real1d(1 + real(reshape(u_over_s, (/(npol+1)**2/))), &
                                1 + real(reshape(u_over_s_ref, (/(npol+1)**2/))), &
                                1e-7, 'f over s axial')
  
  nodes(1,:) = [1,0]
  nodes(2,:) = [2,0]
  nodes(3,:) = [2,1]
  nodes(4,:) = [1,1]

  axial = .false.
  
  do jpol=0, npol 
     u(:,jpol) = (gll_points + 1) / 2d0 + 1
  enddo
  u_over_s_ref(0:npol,0:npol) = 1

  u_over_s(0:npol,0:npol) = f_over_s(u, G2, G2T, gll_points, gll_points, npol, nodes, &
                                     element_type, axial)

  call assert_comparable_real1d(1 + real(reshape(u_over_s, (/(npol+1)**2/))), &
                                1 + real(reshape(u_over_s_ref, (/(npol+1)**2/))), &
                                1e-7, 'f over s non axial')

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_dsdf()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: G1(:,:), G1T(:,:)
  real(kind=dp), allocatable   :: u(:,:), dsdu(:), dsdu_ref(:)

  real(kind=dp), allocatable   :: glj_points(:)
  real(kind=dp), allocatable   :: gll_points(:)

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1

  npol = 1

  G1 = def_lagrange_derivs_glj(npol)
  G2 = def_lagrange_derivs_gll(npol)
  allocate(G1T(0:npol,0:npol))
  allocate(G2T(0:npol,0:npol))
  G1T = transpose(G1)
  G2T = transpose(G2)

  glj_points = zemngl2(npol)
  gll_points = zelegl(npol)

  allocate(u(0:npol,0:npol))
  allocate(dsdu_ref(0:npol))
  
  u = 1
  dsdu_ref = 0
  dsdu = dsdf_axis(u, G2, G1T, glj_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(dsdu), 1 + real(dsdu_ref), &
                                1e-7, 'gradient of a constant = 0')

  u(0,:) = 0
  u(1,:) = 1
  dsdu_ref = 1
  dsdu = dsdf_axis(u, G2, G1T, glj_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(dsdu), 1 + real(dsdu_ref), &
                                1e-7, 'gradient of a constant = 0')

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [2,1]
  nodes(4,:) = [0,1]

  u(0,:) = 0
  u(1,:) = 1
  dsdu_ref = [1d0, 0.5d0]
  dsdu = dsdf_axis(u, G2, G1T, glj_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(dsdu), 1 + real(dsdu_ref), &
                                1e-7, 'gradient of a constant = 0')


end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_td_dsdf()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol, nsamp
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: u(:,:,:), dsdu(:,:), dsdu_ref(:,:)

  real(kind=dp), allocatable   :: gll_points(:)

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1

  npol = 1
  nsamp = 2

  G2 = def_lagrange_derivs_gll(npol)
  allocate(G2T(0:npol,0:npol))
  G2T = transpose(G2)

  gll_points = zelegl(npol)

  allocate(u(1:nsamp,0:npol,0:npol))
  allocate(dsdu_ref(1:nsamp,0:npol))
  
  u(1,0,:) = 0
  u(1,1,:) = 1
  dsdu_ref(1,:) = 1
  u(2,0,:) = 0
  u(2,1,:) = -1
  dsdu_ref(2,:) = -1
  dsdu = dsdf_axis(u, G2, G2T, gll_points, gll_points, npol, nsamp, nodes, element_type)
  call assert_comparable_real1d(1 + real(dsdu(1,:)), 1 + real(dsdu_ref(1,:)), &
                                1e-7, 'td dsdf')

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_gradient()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: u(:,:), grad_u(:,:,:), grad_u_ref(:,:,:)

  real(kind=dp), allocatable   :: gll_points(:)

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1

  npol = 1

  G2 = def_lagrange_derivs_gll(npol)
  allocate(G2T(0:npol,0:npol))
  G2T = transpose(G2)

  gll_points = zelegl(npol)

  allocate(u(0:npol,0:npol))
  allocate(grad_u_ref(0:npol,0:npol,2))
  
  u = 1
  grad_u_ref = 0
  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'gradient of a constant = 0')

  u(0,:) = 0
  u(1,:) = 1
  grad_u_ref(:,:,1) = 1
  grad_u_ref(:,:,2) = 0
  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'linear in s')

  u(:,0) = 0
  u(:,1) = 1
  grad_u_ref(:,:,1) = 0
  grad_u_ref(:,:,2) = 1
  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'linear in z')

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_gradient2()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol
  real(kind=dp), allocatable   :: G1(:,:), G1T(:,:)
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: u(:,:), grad_u(:,:,:), grad_u_ref(:,:,:)

  real(kind=dp), allocatable   :: glj_points(:)
  real(kind=dp), allocatable   :: gll_points(:)
  integer                      :: ipol

  nodes(1,:) = [-1,-1]
  nodes(2,:) = [1,-1]
  nodes(3,:) = [1,1]
  nodes(4,:) = [-1,1]

  ! linear element
  element_type = 1

  npol = 4

  G1 = def_lagrange_derivs_glj(npol)
  G2 = def_lagrange_derivs_gll(npol)

  allocate(G1T(0:npol,0:npol))
  allocate(G2T(0:npol,0:npol))
  G1T = transpose(G1)
  G2T = transpose(G2)

  glj_points = zemngl2(npol)
  gll_points = zelegl(npol)

  allocate(u(0:npol,0:npol))
  allocate(grad_u_ref(0:npol,0:npol,2))

  do ipol = 0, npol
     u(:,ipol) = gll_points
  enddo
  grad_u_ref(:,:,1) = 1
  grad_u_ref(:,:,2) = 0
  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'nonaxial, npol = 4, linera in s')
  
  do ipol = 0, npol
     u(ipol,:) = gll_points
  enddo
  grad_u_ref(:,:,1) = 0
  grad_u_ref(:,:,2) = 1
  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'nonaxial, npol = 4, linera in z')

  do ipol = 0, npol
     u(:,ipol) = glj_points
  enddo
  grad_u_ref(:,:,1) = 1
  grad_u_ref(:,:,2) = 0
  grad_u = axisym_gradient(u, G2, G1T, glj_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'axial, npol = 4, linera in s')

  do ipol = 0, npol
     u(ipol,:) = gll_points
  enddo
  grad_u_ref(:,:,1) = 0
  grad_u_ref(:,:,2) = 1
  grad_u = axisym_gradient(u, G2, G1T, glj_points, gll_points, npol, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/(npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/(npol+1)**2 * 2/))), &
                                1e-7, 'axial, npol = 4, linera in z')

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_td_gradient()

  real(kind=dp)                :: nodes(4,2)
  integer                      :: element_type

  integer                      :: npol, nsamp
  real(kind=dp), allocatable   :: G1(:,:), G1T(:,:)
  real(kind=dp), allocatable   :: G2(:,:), G2T(:,:)
  real(kind=dp), allocatable   :: u(:,:,:), grad_u(:,:,:,:), grad_u_ref(:,:,:,:)

  real(kind=dp), allocatable   :: glj_points(:)
  real(kind=dp), allocatable   :: gll_points(:)

  nodes(1,:) = [0,0]
  nodes(2,:) = [1,0]
  nodes(3,:) = [1,1]
  nodes(4,:) = [0,1]

  ! linear element
  element_type = 1

  npol = 1
  nsamp = 2

  G1 = def_lagrange_derivs_glj(npol)
  G2 = def_lagrange_derivs_gll(npol)

  allocate(G1T(0:npol,0:npol))
  allocate(G2T(0:npol,0:npol))
  G1T = transpose(G1)
  G2T = transpose(G2)

  glj_points = zemngl2(npol)
  gll_points = zelegl(npol)

  allocate(u(1:nsamp,0:npol,0:npol))
  allocate(grad_u_ref(1:nsamp,0:npol,0:npol,2))
  
  u = 1
  grad_u_ref = 0
  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nsamp, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/nsamp * (npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/nsamp * (npol+1)**2 * 2/))), &
                                1e-7, 'gradient of a constant = 0')

  u(1,0,:) = 0
  u(1,1,:) = 1
  grad_u_ref(1,:,:,1) = 1
  grad_u_ref(1,:,:,2) = 0
  u(2,:,0) = 0
  u(2,:,1) = 1
  grad_u_ref(2,:,:,1) = 0
  grad_u_ref(2,:,:,2) = 1

  grad_u = axisym_gradient(u, G2, G2T, gll_points, gll_points, npol, nsamp, nodes, element_type)
  call assert_comparable_real1d(1 + real(reshape(grad_u, (/nsamp * (npol+1)**2 * 2/))), &
                                1 + real(reshape(grad_u_ref, (/nsamp * (npol+1)**2 * 2/))), &
                                1e-7, '1: linear in s, 2: linear in z')

end subroutine 
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
