

module source
    use global_parameters,               only : sp, dp, pi, deg2rad

    type src_param_type
        real(kind=sp)                        :: mij(6)
        real(kind=dp)                        :: colat, lat, lon
        real(kind=dp)                        :: colatd, latd, lond
        real(kind=dp), dimension(3,3)        :: rot_mat, trans_rot_mat
        contains
           procedure, pass                   :: init
           procedure, pass                   :: def_rot_matrix
    end type
contains

subroutine init(this, lat, lon, mij)
   class(src_param_type)      :: this

   real(kind=dp), intent(in)  :: lat, lon, mij(6)

   this%latd   = lat
   this%lond   = lon
   this%colatd = 90 - this%latd
   this%colat  = this%colatd * deg2rad
   this%lon    = this%lond   * deg2rad

   this%mij    = mij
   call this%def_rot_matrix()

end subroutine

subroutine def_rot_matrix(this)
! This function defines the rotation matrix to rotate coordinates to the 
! forward source system. Taken from AxiSEM solver.
  class(src_param_type)          :: this
  real(kind=dp)                  :: srccolat, srclon

  srccolat = this%colat
  srclon   = this%lon

  ! This is the rotation matrix of Nissen-Meyer, Dahlen, Fournier, GJI 2007.
  this%rot_mat(1,1) =  dcos(srccolat) * dcos(srclon)
  this%rot_mat(2,1) =  dcos(srccolat) * dsin(srclon)
  this%rot_mat(3,1) = -dsin(srccolat)
  this%rot_mat(1,2) = -dsin(srclon)
  this%rot_mat(2,2) =  dcos(srclon)
  this%rot_mat(3,2) =  0.d0
  this%rot_mat(1,3) =  dsin(srccolat) * dcos(srclon)
  this%rot_mat(2,3) =  dsin(srccolat) * dsin(srclon)
  this%rot_mat(3,3) =  dcos(srccolat)

  this%trans_rot_mat = transpose(this%rot_mat)

end subroutine def_rot_matrix
end module

!------------------------------------------------------------------------------

module receiver
    
    use global_parameters,               only : sp, dp, pi, deg2rad
    implicit none

    type rec_param_type
        character(len=1)                     :: component
        real(kind=dp)                        :: colat, lat, lon
        real(kind=dp)                        :: colatd, latd, lond
        real(kind=dp)                        :: theta, phi  !< receiver coordinates in the source system
        integer                              :: nkernel
        contains
           procedure, pass                   :: rotate_receiver
           procedure, pass                   :: init
        !type(kernelspec_type), pointer       :: kernel(:)
    end type

contains

subroutine init(this, lat, lon, component)
   class(rec_param_type)        :: this
   real(kind=sp), intent(in)    :: lat, lon
   character(len=1), intent(in) :: component

   this%component = component

   this%latd   = lat
   this%lond   = lon
   this%colatd = 90 - this%latd

   this%colat = this%colatd * deg2rad
   this%lon   = this%lond   * deg2rad

end subroutine

subroutine rotate_receiver(this, trans_rot_mat)
   class(rec_param_type)      :: this
   real(kind=dp), intent(in)  :: trans_rot_mat(3,3)

   real(kind=dp)              :: x_vec(3), x_vec_rot(3)
   real(kind=dp)              :: r_r

   x_vec(1) = dsin(this%colat) * dcos(this%lon)
   x_vec(2) = dsin(this%colat) * dsin(this%lon)
   x_vec(3) = dcos(this%colat)

   x_vec_rot=matmul(trans_rot_mat,x_vec)
   
   r_r = dsqrt(x_vec_rot(1)**2 + x_vec_rot(2)**2 + x_vec_rot(3)**2)
   this%theta = dacos( x_vec_rot(3) / r_r )
   
   if (x_vec_rot(2) >= 0.d0) then
      this%phi =          dacos( x_vec_rot(1) / (r_r * dsin(this%theta)) )
   else
      this%phi = 2 * pi - dacos( x_vec_rot(1) / (r_r * dsin(this%theta)) )
   end if

end subroutine
end module
