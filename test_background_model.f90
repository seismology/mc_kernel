!=========================================================================================
module test_background_models

  use global_parameters
  use backgroundmodel
  use ftnunit

  implicit none
  public
contains
!-----------------------------------------------------------------------------------------
subroutine test_background_models_combine
  type(backgroundmodel_type) :: bm
  real(kind=sp), allocatable :: coeffs(:,:)
  real(kind=sp)              :: ref_vs (1)
  real(kind=sp)              :: ref_vp (1)
  real(kind=sp)              :: ref_rho(1)
  real(kind=sp)              :: ref_vsh(1)
  real(kind=sp)              :: ref_vsv(1)
  real(kind=sp)              :: ref_vph(1)
  real(kind=sp)              :: ref_vpv(1)
  real(kind=sp)              :: ref_eta(1)
  real(kind=sp)              :: ref_phi(1)
  real(kind=sp)              :: ref_xi(1)

  allocate(coeffs(6,1))
  coeffs(:,1) = [4500.0, 3000.0, 2000.0, 1.1, 0.9, 1.0]

  ref_vs  = 3000.0
  ref_vp  = 4500.0
  ref_rho = 2000.0
  ref_eta = 1
  ref_phi = 1.1
  ref_xi  = 0.9
  ref_vsh = ref_vs
  ref_vsv = sqrt(ref_vsh**2/ref_xi)
  ref_vph = ref_vp
  ref_vpv = sqrt(ref_vph**2*ref_phi)

  call bm%combine(coeffs)

  call assert_comparable(bm%c_vs,  ref_vs,  1e-7, 'VS correct')
  call assert_comparable(bm%c_vp,  ref_vp,  1e-7, 'VP correct')
  call assert_comparable(bm%c_rho, ref_rho, 1e-7, 'Rho correct')
  call assert_comparable(bm%c_vsh, ref_vsh, 1e-7, 'VSh correct')
  call assert_comparable(bm%c_vsv, ref_vsv, 1e-7, 'VSv correct')
  call assert_comparable(bm%c_vph, ref_vph, 1e-7, 'VPh correct')
  call assert_comparable(bm%c_vpv, ref_vpv, 1e-7, 'VPv correct')
  call assert_comparable(bm%c_eta, ref_eta, 1e-7, 'Eta correct')
  call assert_comparable(bm%c_phi, ref_phi, 1e-7, 'Phi correct')
  call assert_comparable(bm%c_xi,  ref_xi,  1e-7, 'Xi correct')

  deallocate(coeffs)
  allocate(coeffs(6,2))
  coeffs(:,1) = [11000.0, 7000.0, 12000.0, 1.2, 0.7, 2.0]
  coeffs(:,2) = [11000.0, 7000.0, 12000.0, 1.2, 0.7, 2.0]

  ref_vs  = 7000.0
  ref_vp  = 11000.0
  ref_rho = 12000.0
  ref_eta = 2
  ref_phi = 1.2
  ref_xi  = 0.7
  ref_vsh = ref_vs
  ref_vsv = sqrt(ref_vsh**2/ref_xi)
  ref_vph = ref_vp
  ref_vpv = sqrt(ref_vph**2*ref_phi)

  call bm%combine(coeffs)

  call assert_comparable(bm%c_vs(1),  ref_vs(1),  1e-7, 'VS correct')
  call assert_comparable(bm%c_vp(1),  ref_vp(1),  1e-7, 'VP correct')
  call assert_comparable(bm%c_rho(1), ref_rho(1), 1e-7, 'Rho correct')
  call assert_comparable(bm%c_vsh(1), ref_vsh(1), 1e-7, 'VSh correct')
  call assert_comparable(bm%c_vsv(1), ref_vsv(1), 1e-7, 'VSv correct')
  call assert_comparable(bm%c_vph(1), ref_vph(1), 1e-7, 'VPh correct')
  call assert_comparable(bm%c_vpv(1), ref_vpv(1), 1e-7, 'VPv correct')
  call assert_comparable(bm%c_eta(1), ref_eta(1), 1e-7, 'Eta correct')
  call assert_comparable(bm%c_phi(1), ref_phi(1), 1e-7, 'Phi correct')
  call assert_comparable(bm%c_xi(1),  ref_xi(1),  1e-7, 'Xi correct')
  

end subroutine test_background_models_combine
!-----------------------------------------------------------------------------------------

end module test_background_models
!=========================================================================================
