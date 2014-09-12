!=========================================================================================
module backgroundmodel
use global_parameters, only: sp, dp

 type backgroundmodel_type
   real(kind=sp), allocatable :: c_vp(:)  ! contains velocities at each point in km/s
   real(kind=sp), allocatable :: c_vs(:)  ! contains velocities at each point in km/s
   real(kind=sp), allocatable :: c_rho(:) ! contains velocities at each point in kg/(km^3)

  !@ TODO : For now assume an isotropic background model
  !         need to store xi in the netcdf files
   real(kind=sp), allocatable :: c_vsh(:)
   real(kind=sp), allocatable :: c_vsv(:)
   real(kind=sp), allocatable :: c_vph(:)
   real(kind=sp), allocatable :: c_vpv(:)
   real(kind=sp), allocatable :: c_eta(:)
   contains 
     procedure, pass            :: recombine
 end type

 contains

!-------------------------------------------------------------------------------
 subroutine recombine(this, coeffs)
   class(backgroundmodel_type) :: this
   real(kind=dp), intent(in)  :: coeffs(:,:)
   integer                    :: npoints


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
   this%c_vp  = coeffs(1,:) ! contains velocities at each point in km/s
   this%c_vs  = coeffs(2,:) ! contains velocities at each point in km/s
   this%c_rho = coeffs(3,:) ! contains velocities at each point in kg/(km^3)

   ! @ TODO : For now assume an isotropic background this
   !          need to store xi in the netcdf files
   if (size(coeffs, 1) == 3) then !isotropic case
     this%c_vsh = this%c_vs
     this%c_vsv = this%c_vs
     this%c_vph = this%c_vp
     this%c_vpv = this%c_vp
     this%c_eta = 1.d0
   else
     !! SOMETHING
     stop
   end if
 end subroutine recombine
!-------------------------------------------------------------------------------

end module backgroundmodel
