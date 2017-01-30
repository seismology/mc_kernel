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
        real(kind=dp), dimension(3)    :: r                  !< receiver location 
                                                             !! in cartesian coordinates
                                                             !! in meters!
        real(kind=dp)                  :: theta, phi         !< receiver coordinates 
                                                             !! in the source system
        real(kind=dp), dimension(3,3)  :: rot_mat, trans_rot_mat
        integer                        :: nkernel
        integer                        :: istf   !< Which stf group does this receiver 
                                                 !! belong to?
        integer                        :: firstkernel, lastkernel
        logical                        :: needs_basekernel(6)
        character(len=32)              :: strain_type
        type(kernelspec_type), pointer :: kernel(:)
        real(kind=dp)                  :: t_last_p=0, t_last_s=0 
        contains
           procedure, pass      :: rotate_receiver
           procedure, pass      :: init
           procedure, pass      :: get_r
    end type

contains

!-----------------------------------------------------------------------------------------
subroutine init(this, name, lat, lon, component, nkernel, firstkernel, lastkernel, istf)
   class(rec_param_type)         :: this
   character(len=16), intent(in) :: name
   real(kind=dp), intent(in)     :: lat, lon
   character(len=1), intent(in)  :: component
   integer, intent(in)           :: nkernel, firstkernel, lastkernel
   integer, intent(in)           :: istf

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

   this%istf = istf

   allocate(this%kernel(this%nkernel))

end subroutine init
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_r(this, planet_radius) result(r)
   class(rec_param_type)                :: this
   real(kind=dp), intent(in), optional  :: planet_radius
   real(kind=dp)                        :: radius
   real(kind=dp)                        :: r(3)

   if (present(planet_radius)) then
     radius = planet_radius
   else 
     radius = 6371d3
   end if

   r(1) = cos(this%lat) * cos(this%lon) * radius
   r(2) = cos(this%lat) * sin(this%lon) * radius
   r(3) = sin(this%lat)                 * radius

end function get_r
!----------------------------------------------------------------------------------------
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

!-----------------------------------------------------------------------------------------
subroutine freeme(this)
   class(rec_param_type)         :: this

   deallocate(this%kernel)

end subroutine freeme
!-----------------------------------------------------------------------------------------
end module
!=========================================================================================
