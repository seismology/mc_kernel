!=========================================================================================
module backgroundmodel
  use global_parameters, only: sp, dp

  type backgroundmodel_type
    real(kind=sp), allocatable :: c_vp(:)  ! velocities at each point in m/s
    real(kind=sp), allocatable :: c_vs(:)  ! velocities at each point in m/s
    real(kind=sp), allocatable :: c_rho(:) ! densities at each point in kg/(m^3)

    real(kind=sp), allocatable :: c_vsh(:)
    real(kind=sp), allocatable :: c_vsv(:)
    real(kind=sp), allocatable :: c_vph(:)
    real(kind=sp), allocatable :: c_vpv(:)
    real(kind=sp), allocatable :: c_eta(:)
    real(kind=sp), allocatable :: c_phi(:)
    real(kind=sp), allocatable :: c_xi(:)
    contains 
      procedure, pass            :: recombine
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine recombine(this, coeffs)
  class(backgroundmodel_type) :: this
  real(kind=sp), intent(in)   :: coeffs(:,:)
  integer                     :: npoints


  npoints = size(coeffs, 2)

  allocate(this%c_vp (npoints))
  allocate(this%c_vs (npoints))
  allocate(this%c_rho(npoints))
  allocate(this%c_vsh(npoints))
  allocate(this%c_vsv(npoints))
  allocate(this%c_vph(npoints))
  allocate(this%c_vpv(npoints))
  allocate(this%c_eta(npoints))

  ! Recombine this coefficients for chosen parameterization
  this%c_vp  = coeffs(1,:)
  this%c_vs  = coeffs(2,:)
  this%c_rho = coeffs(3,:)

  this%c_phi = coeffs(4,:) 
  this%c_xi  = coeffs(5,:) 
  this%c_eta = coeffs(6,:) 
  
  ! from axisem:
  !lambda = rho * (vphtmp**2 - two * vshtmp**2)
  !mu = rho * vshtmp**2
  !xi  = vsh**2 / vsv**2
  !phi = vpv**2 / vph**2

  ! TODO: double check correctness. atm should be consistent with axisem
  ! references:
  !     Nolet(2008), Eq. (16.2)
  !     MvD [Anisotropy Notes, p. 13.4]
  this%c_vsh = this%c_vs
  this%c_vph = this%c_vp

  this%c_vsv = this%c_vsh / sqrt(xi)
  this%c_vpv = this%c_vph * sqrt(phi)
  !this%c_eta = 1.d0

end subroutine 
!-----------------------------------------------------------------------------------------

end module backgroundmodel
!=========================================================================================
