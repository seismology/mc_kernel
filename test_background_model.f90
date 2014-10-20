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
  real(kind=dp), allocatable :: params(:,:)
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

!-----------------------------------------------------------------------------------------
subroutine test_background_models_weight
  type(backgroundmodel_type) :: bm
  real(kind=sp), allocatable :: coeffs(:,:)
  real(kind=dp), allocatable :: params(:,:)
  real(kind=dp)              :: ref_vs (1)
  real(kind=dp)              :: ref_vp (1)
  real(kind=dp)              :: ref_rho(1)
  real(kind=dp)              :: ref_vsh(1)
  real(kind=dp)              :: ref_vsv(1)
  real(kind=dp)              :: ref_vph(1)
  real(kind=dp)              :: ref_vpv(1)
  real(kind=dp)              :: ref_eta(1)
  real(kind=dp)              :: ref_phi(1)
  real(kind=dp)              :: ref_xi(1)

  allocate(coeffs(6,2))
  coeffs(:,1) = [4500.0, 3000.0, 2000.0, 1.1, 0.9, 1.0]
  coeffs(:,2) = [4500.0, 3000.0, 2000.0, 1.1, 0.9, 1.0]

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

  allocate(params(10, 2))
  params = bm%weight([1.d0, 3.d0])

  ! First weight was 1
  call assert_comparable(params( 1,1), ref_vp(1),  1d-7, 'VP correct')
  call assert_comparable(params( 2,1), ref_vs(1),  1d-7, 'VS correct')
  call assert_comparable(params( 3,1), ref_rho(1), 1d-7, 'Rho correct')
  call assert_comparable(params( 4,1), ref_vph(1), 1d-7, 'VPh correct')
  call assert_comparable(params( 5,1), ref_vpv(1), 1d-7, 'VPv correct')
  call assert_comparable(params( 6,1), ref_vsh(1), 1d-7, 'VSh correct')
  call assert_comparable(params( 7,1), ref_vsv(1), 1d-7, 'VSv correct')
  call assert_comparable(params( 8,1), ref_eta(1), 1d-7, 'Eta correct')
  call assert_comparable(params( 9,1), ref_phi(1), 1d-7, 'Phi correct')
  call assert_comparable(params(10,1), ref_xi(1),  1d-7, 'Xi correct')

  ! Second weight is 3
  call assert_comparable(params( 1,2), 3.d0*ref_vp(1),  1d-7, 'VP correct')
  call assert_comparable(params( 2,2), 3.d0*ref_vs(1),  1d-7, 'VS correct')
  call assert_comparable(params( 3,2), 3.d0*ref_rho(1), 1d-7, 'Rho correct')
  call assert_comparable(params( 4,2), 3.d0*ref_vph(1), 1d-7, 'VPh correct')
  call assert_comparable(params( 5,2), 3.d0*ref_vpv(1), 1d-7, 'VPv correct')
  call assert_comparable(params( 6,2), 3.d0*ref_vsh(1), 1d-7, 'VSh correct')
  call assert_comparable(params( 7,2), 3.d0*ref_vsv(1), 1d-7, 'VSv correct')
  call assert_comparable(params( 8,2), 3.d0*ref_eta(1), 1d-7, 'Eta correct')
  call assert_comparable(params( 9,2), 3.d0*ref_phi(1), 1d-7, 'Phi correct')
  call assert_comparable(params(10,2), 3.d0*ref_xi(1),  1d-7, 'Xi correct')

end subroutine test_background_models_weight
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_background_models_get_parameter_names()
  type(backgroundmodel_type)    :: bg_model
  character(len=3), allocatable :: parameter_names(:)
  character(len=3), parameter   :: parameter_names_ref(10) =            &
                                  ['vp ', 'vs ', 'rho', 'vph', 'vpv',  &
                                   'vsh', 'vsv', 'eta', 'phi', 'xi ']
  
  parameter_names = bg_model%get_parameter_names()

  call assert_true(parameter_names==parameter_names_ref, 'Model parameter names are correct')

end subroutine test_background_models_get_parameter_names
!-----------------------------------------------------------------------------------------

end module test_background_models
!=========================================================================================
