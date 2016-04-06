!=========================================================================================
module background_model
  use global_parameters, only: sp, dp

  implicit none

  integer, parameter           :: nmodel_parameters = 12 !< Number of basic model parameters
                                                         !! which are availabe in backgroundmodel_type

  character(len=3), parameter  :: parameter_name(nmodel_parameters) =  &
                                  ['vp ', 'vs ', 'rho', 'vph', 'vpv', 'vsh', &
                                   'vsv', 'eta', 'phi', 'xi ', 'lam', 'mu ']

  type backgroundmodel_type

    real(kind=sp), allocatable :: c_vp(:)
    real(kind=sp), allocatable :: c_vs(:)
    real(kind=sp), allocatable :: c_rho(:)
    real(kind=sp), allocatable :: c_vph(:)
    real(kind=sp), allocatable :: c_vpv(:)
    real(kind=sp), allocatable :: c_vsh(:)
    real(kind=sp), allocatable :: c_vsv(:)
    real(kind=sp), allocatable :: c_eta(:)
    real(kind=sp), allocatable :: c_phi(:)
    real(kind=sp), allocatable :: c_xi(:)
    real(kind=sp), allocatable :: c_lam(:)
    real(kind=sp), allocatable :: c_mu(:)

    real(kind=sp)              :: max_vp, max_vs

    contains 
      procedure, pass          :: init
      procedure, pass          :: combine
      procedure, pass          :: weight
  end type

contains

!-----------------------------------------------------------------------------------------
!> Just allocate the members of backgroundmodel_type to size npoints
subroutine init(this, npoints)
  class(backgroundmodel_type) :: this
  integer, intent(in)         :: npoints

  if (allocated(this%c_vp)) then
    if (size(this%c_vp).ne.npoints) then
      deallocate(this%c_vp )
      deallocate(this%c_vs )
      deallocate(this%c_rho)
      deallocate(this%c_vsh)
      deallocate(this%c_vsv)
      deallocate(this%c_vph)
      deallocate(this%c_vpv)
      deallocate(this%c_eta)
      deallocate(this%c_phi)
      deallocate(this%c_xi)
      deallocate(this%c_lam)
      deallocate(this%c_mu)
      allocate(this%c_vp (npoints))
      allocate(this%c_vs (npoints))
      allocate(this%c_rho(npoints))
      allocate(this%c_vsh(npoints))
      allocate(this%c_vsv(npoints))
      allocate(this%c_vph(npoints))
      allocate(this%c_vpv(npoints))
      allocate(this%c_eta(npoints))
      allocate(this%c_phi(npoints))
      allocate(this%c_xi(npoints))
      allocate(this%c_lam(npoints))
      allocate(this%c_mu(npoints))
    end if
  else
    allocate(this%c_vp (npoints))
    allocate(this%c_vs (npoints))
    allocate(this%c_rho(npoints))
    allocate(this%c_vsh(npoints))
    allocate(this%c_vsv(npoints))
    allocate(this%c_vph(npoints))
    allocate(this%c_vpv(npoints))
    allocate(this%c_eta(npoints))
    allocate(this%c_phi(npoints))
    allocate(this%c_xi(npoints))
    allocate(this%c_lam(npoints))
    allocate(this%c_mu(npoints))
  end if

  this%max_vp = maxval(this%c_vp)
  this%max_vs = maxval(this%c_vs)

end subroutine init
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Allocate the members of backgroundmodel_type to size npoints and fill with the
!! values given in coeffs (these are the values stored in the SEM mesh for an anisotropic
!! run). Where necessary, these values are combined to calculate the ones not stored 
!! directly in the file (like vpv...)
subroutine combine(this, coeffs)
  class(backgroundmodel_type) :: this
  real(kind=sp), intent(in)   :: coeffs(:,:)
  integer                     :: npoints

  npoints = size(coeffs, 2)

  if (allocated(this%c_vp)) then
    if (size(this%c_vp).ne.npoints) then
      deallocate(this%c_vp )
      deallocate(this%c_vs )
      deallocate(this%c_rho)
      deallocate(this%c_vsh)
      deallocate(this%c_vsv)
      deallocate(this%c_vph)
      deallocate(this%c_vpv)
      deallocate(this%c_eta)
      deallocate(this%c_phi)
      deallocate(this%c_xi)
      deallocate(this%c_lam)
      deallocate(this%c_mu)
      allocate(this%c_vp (npoints))
      allocate(this%c_vs (npoints))
      allocate(this%c_rho(npoints))
      allocate(this%c_vsh(npoints))
      allocate(this%c_vsv(npoints))
      allocate(this%c_vph(npoints))
      allocate(this%c_vpv(npoints))
      allocate(this%c_eta(npoints))
      allocate(this%c_phi(npoints))
      allocate(this%c_xi(npoints))
      allocate(this%c_lam(npoints))
      allocate(this%c_mu(npoints))
    end if
  else
    allocate(this%c_vp (npoints))
    allocate(this%c_vs (npoints))
    allocate(this%c_rho(npoints))
    allocate(this%c_vsh(npoints))
    allocate(this%c_vsv(npoints))
    allocate(this%c_vph(npoints))
    allocate(this%c_vpv(npoints))
    allocate(this%c_eta(npoints))
    allocate(this%c_phi(npoints))
    allocate(this%c_xi(npoints))
    allocate(this%c_lam(npoints))
    allocate(this%c_mu(npoints))
  end if

  ! Recombine this coefficients for chosen parameterization
  this%c_vp  = coeffs(1,:)
  this%c_vs  = coeffs(2,:)
  this%c_rho = coeffs(3,:)

  this%c_phi = coeffs(4,:) 
  this%c_xi  = coeffs(5,:) 
  this%c_eta = coeffs(6,:) 
  
  this%c_lam = coeffs(3,:) * (coeffs(1,:)**2 - 2 * coeffs(2,:)**2)
  this%c_mu  = coeffs(3,:) *  coeffs(2,:)**2
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

  this%c_vsv = this%c_vsh / sqrt(coeffs(5,:))
  this%c_vpv = this%c_vph * sqrt(coeffs(4,:))

end subroutine combine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_parameter_names()
  character(len=3)            :: get_parameter_names(nmodel_parameters)

  get_parameter_names = parameter_name
end function get_parameter_names
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_parameter_index(model_parameter)
  use commpi, only   : pabort
  character(len=3)  :: model_parameter
  integer           :: get_parameter_index

  integer           :: iparam

   ! Determine the index of the model parameter in the list defined in backgroundmodel.f90
   do iparam = 1, nmodel_parameters
     if (model_parameter == parameter_name(iparam)) then
       get_parameter_index = iparam
       exit
     end if
   end do
   if (iparam == nmodel_parameters + 1) then
     print '("ERROR: Unknown model parameter for kernel", A)', &
       trim(model_parameter)
     print '("Available options: ", 10(A3))', parameter_name
     call pabort(do_traceback=.false.)
   end if

end function get_parameter_index
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function weight(this, weights) result(all_coeffs)
  class(backgroundmodel_type) :: this
  real(kind=dp), intent(in)   :: weights(:)
  real(kind=dp)               :: all_coeffs(nmodel_parameters, size(weights,1))

  integer                     :: npoint, ipoint

  npoint = size(weights)

  all_coeffs(1,  :) = this%c_vp (:) 
  all_coeffs(2,  :) = this%c_vs (:) 
  all_coeffs(3,  :) = this%c_rho(:) 
  all_coeffs(4,  :) = this%c_vph(:) 
  all_coeffs(5,  :) = this%c_vpv(:) 
  all_coeffs(6,  :) = this%c_vsh(:) 
  all_coeffs(7,  :) = this%c_vsv(:) 
  all_coeffs(8,  :) = this%c_eta(:) 
  all_coeffs(9,  :) = this%c_phi(:) 
  all_coeffs(10, :) = this%c_xi (:) 
  all_coeffs(11, :) = this%c_lam(:) 
  all_coeffs(12, :) = this%c_mu (:) 

  do ipoint = 1, npoint
     all_coeffs(:, ipoint)  = all_coeffs(:, ipoint) * weights(ipoint)
  end do

end function weight
!-----------------------------------------------------------------------------------------

end module background_model
!=========================================================================================
