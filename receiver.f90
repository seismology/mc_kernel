!=========================================================================================
module receiver_class
    
    use global_parameters,  only : sp, dp, pi, deg2rad, rad2deg
    use kernel,             only : kernelspec_type
    implicit none

    type rec_param_type
        character(len=1)               :: component
        character(len=16)              :: name
        real(kind=dp)                  :: colat,  lat,  lon  !< receiver coordinates
                                                             !! in the earth system
                                                             !! (radians)
        real(kind=dp)                  :: colatd, latd, lond !< receiver coordinates
                                                             !! in the earth system
                                                             !! (degrees)
        real(kind=dp)                  :: theta, phi         !< receiver coordinates 
                                                             !! in the source system
        real(kind=dp), dimension(3,3)  :: rot_mat, trans_rot_mat
        integer                        :: nkernel
        integer                        :: firstkernel, lastkernel
        type(kernelspec_type), pointer :: kernel(:)
        contains
           procedure, pass      :: rotate_receiver
           procedure, pass      :: init
    end type

contains

!-----------------------------------------------------------------------------------------
subroutine init(this, name, lat, lon, component, nkernel, firstkernel, lastkernel)
   class(rec_param_type)         :: this
   character(len=16), intent(in) :: name
   real(kind=dp), intent(in)     :: lat, lon
   character(len=1), intent(in)  :: component
   integer, intent(in)           :: nkernel, firstkernel, lastkernel

   this%name        = name
   this%component   = component

   this%latd        = lat
   this%lond        = lon
   this%colatd      = 90 - this%latd

   this%colat       = this%colatd * deg2rad
   this%lon         = this%lond   * deg2rad
   this%lat         = this%latd   * deg2rad

   this%nkernel     = nkernel
   this%firstkernel = firstkernel
   this%lastkernel  = lastkernel

   allocate(this%kernel(this%nkernel))

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine rotate_receiver(this, trans_rot_mat)
   implicit none
   class(rec_param_type)      :: this
   real(kind=dp), intent(in)  :: trans_rot_mat(3,3)

   real(kind=dp)              :: x_vec(3), x_vec_rot(3), cosphi

   x_vec(1) = dsin(this%colat) * dcos(this%lon)
   x_vec(2) = dsin(this%colat) * dsin(this%lon)
   x_vec(3) = dcos(this%colat)

   x_vec_rot=matmul(trans_rot_mat,x_vec)
   !print *, 'Rotation matrix:'
   !print '(3(E9.2))', trans_rot_mat
   !print *, 'X-vector original'
   !print '(3(E9.2))', x_vec
   !print *, 'X-vector rotated'
   !print '(3(E9.2))', x_vec_rot

   this%theta = dacos( x_vec_rot(3) )
   
   cosphi = x_vec_rot(1) / dsin(this%theta)
   if (cosphi.gt. 1.0d0)  cosphi =  1.0d0
   if (cosphi.lt.-1.0d0)  cosphi = -1.0d0

   if (x_vec_rot(2) >= 0.d0) then
      this%phi =          dacos( cosphi ) 
   else
      this%phi = 2 * pi - dacos( cosphi )
   end if

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
