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
module sem_derivatives
  use global_parameters,      only : dp
  use finite_elem_mapping,    only : inv_jacobian

  implicit none
  private

  public :: strain_merged_fwd
  public :: strain_merged_bwd
  public :: strain_monopole
  public :: strain_dipole
  public :: strain_quadpole
  public :: straintrace_merged_fwd
  public :: straintrace_merged_bwd
  public :: straintrace_monopole
  public :: straintrace_dipole
  public :: straintrace_quadpole
  public :: axisym_gradient
  public :: dsdf_axis
  public :: f_over_s

  interface strain_monopole
    module procedure  :: strain_monopole
    module procedure  :: strain_monopole_td
  end interface

  interface straintrace_monopole
    module procedure  :: straintrace_monopole
    module procedure  :: straintrace_monopole_td
  end interface

  interface strain_dipole
    module procedure  :: strain_dipole
    module procedure  :: strain_dipole_td
  end interface

  interface straintrace_dipole
    module procedure  :: straintrace_dipole
    module procedure  :: straintrace_dipole_td
  end interface

  interface strain_quadpole
    module procedure  :: strain_quadpole
    module procedure  :: strain_quadpole_td
  end interface

  interface straintrace_quadpole
    module procedure  :: straintrace_quadpole
    module procedure  :: straintrace_quadpole_td
  end interface

  interface f_over_s
    module procedure  :: f_over_s
    module procedure  :: f_over_s_td
  end interface

  interface dsdf_axis
    module procedure  :: dsdf_axis
    module procedure  :: dsdf_axis_td
  end interface

  interface axisym_gradient
    module procedure  :: axisym_gradient
    module procedure  :: axisym_gradient_td
  end interface

  interface mxm
    module procedure  :: mxm
    module procedure  :: mxm_atd
    module procedure  :: mxm_btd
  end interface

  interface mxm_ipol0
    module procedure  :: mxm_ipol0
    module procedure  :: mxm_ipol0_atd
    module procedure  :: mxm_ipol0_btd
  end interface

contains

!--MERGED-------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
pure function strain_merged_fwd(u, G, GT, xi, eta, npol, nsamp, &
                       nodes, element_type, axial) result(strain)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol,10)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: strain(1:nsamp,0:npol,0:npol,6,4)

  ! MZZ
  strain(:,:,:,:,1) = strain_monopole(u(:,:,:,1:2), G, GT, xi, eta, npol, &
                                      nsamp, nodes, element_type, axial)
  ! MXX
  strain(:,:,:,:,2) = strain_monopole(u(:,:,:,3:4), G, GT, xi, eta, npol, &
                                      nsamp, nodes, element_type, axial)
  ! MXZ
  strain(:,:,:,:,3) = strain_dipole(u(:,:,:,5:7), G, GT, xi, eta, npol, &
                                    nsamp, nodes, element_type, axial)
  ! MXY
  strain(:,:,:,:,4) = strain_quadpole(u(:,:,:,8:10), G, GT, xi, eta, npol, &
                                      nsamp, nodes, element_type, axial)

end function strain_merged_fwd
!-----------------------------------------------------------------------------------------
pure function strain_merged_bwd(u, G, GT, xi, eta, npol, nsamp, &
                       nodes, element_type, axial) result(strain)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol,5)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: strain(1:nsamp,0:npol,0:npol,6,2)

  ! PZ
  strain(:,:,:,:,1) = strain_monopole(u(:,:,:,1:2), G, GT, xi, eta, npol, &
                                      nsamp, nodes, element_type, axial)
  ! PX
  strain(:,:,:,:,2) = strain_dipole(u(:,:,:,3:5), G, GT, xi, eta, npol, &
                                    nsamp, nodes, element_type, axial)

end function strain_merged_bwd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function strain_merged(u, G, GT, xi, eta, npol, nsamp, nsim, &
                       nodes, element_type, axial) result(strain)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp, nsim
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol, nint(nsim*2.5))
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: strain(1:nsamp,0:npol,0:npol,6,nsim)

  select case(nsim)
  case(4) 
    ! MZZ
    strain(:,:,:,:,1) = strain_monopole(u(:,:,:,1:2), G, GT, xi, eta, npol, &
                                        nsamp, nodes, element_type, axial)
    ! MXX
    strain(:,:,:,:,2) = strain_monopole(u(:,:,:,3:4), G, GT, xi, eta, npol, &
                                        nsamp, nodes, element_type, axial)
    ! MXZ
    strain(:,:,:,:,3) = strain_dipole(u(:,:,:,5:7), G, GT, xi, eta, npol, &
                                      nsamp, nodes, element_type, axial)
    ! MXY
    strain(:,:,:,:,4) = strain_quadpole(u(:,:,:,8:10), G, GT, xi, eta, npol, &
                                        nsamp, nodes, element_type, axial)

  case(2)
    ! PZ
    strain(:,:,:,:,1) = strain_monopole(u(:,:,:,1:2), G, GT, xi, eta, npol, &
                                        nsamp, nodes, element_type, axial)
    ! PX
    strain(:,:,:,:,2) = strain_dipole(u(:,:,:,3:5), G, GT, xi, eta, npol, &
                                      nsamp, nodes, element_type, axial)
  end select

end function strain_merged
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function straintrace_merged_fwd(u, G, GT, xi, eta, npol, nsamp, nodes, &
                            element_type, axial) result(straintrace)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol,10)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: straintrace(1:nsamp,0:npol,0:npol,4)

  ! MZZ
  straintrace(:,:,:,1) = straintrace_monopole(u(:,:,:,1:2), G, GT, xi, eta, npol, &
                                              nsamp, nodes, element_type, axial)
  ! MXX
  straintrace(:,:,:,2) = straintrace_monopole(u(:,:,:,3:4), G, GT, xi, eta, npol, &
                                              nsamp, nodes, element_type, axial)
  ! MXZ
  straintrace(:,:,:,3) = straintrace_dipole(u(:,:,:,5:7), G, GT, xi, eta, npol, &
                                            nsamp, nodes, element_type, axial)
  ! MXY
  straintrace(:,:,:,4) = straintrace_quadpole(u(:,:,:,8:10), G, GT, xi, eta, npol, &
                                              nsamp, nodes, element_type, axial)

end function straintrace_merged_fwd
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
pure function straintrace_merged_bwd(u, G, GT, xi, eta, npol, nsamp, nodes, &
                            element_type, axial) result(straintrace)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol,5)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: straintrace(1:nsamp,0:npol,0:npol,2)

  ! PZ
  straintrace(:,:,:,1) = straintrace_monopole(u(:,:,:,1:2), G, GT, xi, eta, npol, &
                                              nsamp, nodes, element_type, axial)
  ! PX 
  straintrace(:,:,:,2) = straintrace_dipole(u(:,:,:,3:5), G, GT, xi, eta, npol, &
                                            nsamp, nodes, element_type, axial)

end function straintrace_merged_bwd
!-----------------------------------------------------------------------------------------

!--MONOPOLE-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function strain_monopole_td(u, G, GT, xi, eta, npol, nsamp, nodes, element_type, axial)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol, 2)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: strain_monopole_td(1:nsamp,0:npol,0:npol,6)
  
  real(kind=dp)                 :: grad_buff1(1:nsamp,0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff2(1:nsamp,0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,:,1), G, GT, xi, eta, npol, nsamp, &
                               nodes, element_type)
  
  ! 1: dsuz, 2: dzuz
  grad_buff2 = axisym_gradient(u(:,:,:,2), G, GT, xi, eta, npol, nsamp, &
                               nodes, element_type)

  strain_monopole_td(:,:,:,1) = grad_buff1(:,:,:,1)
  strain_monopole_td(:,:,:,2) = f_over_s(u(:,:,:,1), G, GT, xi, eta, npol, nsamp, &
                                         nodes, element_type, axial)
  strain_monopole_td(:,:,:,3) = grad_buff2(:,:,:,2)
  strain_monopole_td(:,:,:,4) = 0
  strain_monopole_td(:,:,:,5) = (grad_buff1(:,:,:,2) + grad_buff2(:,:,:,1)) / 2d0
  strain_monopole_td(:,:,:,6) = 0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function strain_monopole(u, G, GT, xi, eta, npol, nodes, element_type, axial)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: u(0:npol,0:npol, 2)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: strain_monopole(0:npol,0:npol,6)
  
  real(kind=dp)                 :: grad_buff1(0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff2(0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,1), G, GT, xi, eta, npol, nodes, element_type)
  
  ! 1: dsuz, 2: dzuz
  grad_buff2 = axisym_gradient(u(:,:,2), G, GT, xi, eta, npol, nodes, element_type)

  strain_monopole(:,:,1) = grad_buff1(:,:,1)
  strain_monopole(:,:,2) = f_over_s(u(:,:,1), G, GT, xi, eta, npol, nodes, &
                                    element_type, axial)
  strain_monopole(:,:,3) = grad_buff2(:,:,2)
  strain_monopole(:,:,4) = 0
  strain_monopole(:,:,5) = (grad_buff1(:,:,2) + grad_buff2(:,:,1)) / 2d0
  strain_monopole(:,:,6) = 0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function straintrace_monopole_td(u, G, GT, xi, eta, npol, nsamp, nodes, element_type, &
                                 axial)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol, 2)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: straintrace_monopole_td(1:nsamp,0:npol,0:npol)
  
  real(kind=dp)                 :: grad_buff1(1:nsamp,0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff2(1:nsamp,0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,:,1), G, GT, xi, eta, npol, nsamp, &
                               nodes, element_type)
  
  ! 1: dsuz, 2: dzuz
  grad_buff2 = axisym_gradient(u(:,:,:,2), G, GT, xi, eta, npol, nsamp, &
                               nodes, element_type)

  straintrace_monopole_td(:,:,:) &
        = grad_buff1(:,:,:,1)  &
        + f_over_s(u(:,:,:,1), G, GT, xi, eta, npol, nsamp, nodes, element_type, axial) &
        + grad_buff2(:,:,:,2)

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function straintrace_monopole(u, G, GT, xi, eta, npol, nodes, element_type, axial)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: u(0:npol,0:npol, 2)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: straintrace_monopole(0:npol,0:npol)
  
  real(kind=dp)                 :: grad_buff1(0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff2(0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,1), G, GT, xi, eta, npol, nodes, element_type)
  
  ! 1: dsuz, 2: dzuz
  grad_buff2 = axisym_gradient(u(:,:,2), G, GT, xi, eta, npol, nodes, element_type)

  straintrace_monopole(:,:) &
        = grad_buff1(:,:,1) &
        + f_over_s(u(:,:,1), G, GT, xi, eta, npol, nodes, element_type, axial) &
        + grad_buff2(:,:,2)

end function
!-----------------------------------------------------------------------------------------

!--DIPOLE---------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function strain_dipole_td(u, G, GT, xi, eta, npol, nsamp, nodes, element_type, axial)
  ! Computes the strain tensor for displacement u excited bz a dipole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol, 3)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: strain_dipole_td(1:nsamp,0:npol,0:npol,6)
  
  real(kind=dp)                 :: grad_buff1(1:nsamp,0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff2(1:nsamp,0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff3(1:nsamp,0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,:,1), G, GT, xi, eta, npol, nsamp, &
                               nodes, element_type)
  
  ! 1: dsup, 2: dzup
  grad_buff2 = axisym_gradient(u(:,:,:,2), G, GT, xi, eta, npol, nsamp, &
                               nodes, element_type)

  ! 1: dsuz, 2: dzuz
  grad_buff3 = axisym_gradient(u(:,:,:,3), G, GT, xi, eta, npol, nsamp, nodes, &
                               element_type)

  strain_dipole_td(:,:,:,1) = grad_buff1(:,:,:,1)
  strain_dipole_td(:,:,:,2) = f_over_s(u(:,:,:,1) - u(:,:,:,2), G, GT, xi, eta, npol, &
                                       nsamp, nodes, element_type, axial)
  strain_dipole_td(:,:,:,3) = grad_buff3(:,:,:,2)
  strain_dipole_td(:,:,:,4) = - 0.5d0 * (f_over_s(u(:,:,:,3), G, GT, xi, eta, npol, &
                                                  nsamp, nodes, element_type, axial) &
                                         + grad_buff2(:,:,:,2))
  strain_dipole_td(:,:,:,5) = (grad_buff1(:,:,:,2) + grad_buff3(:,:,:,1)) / 2d0
  strain_dipole_td(:,:,:,6) = - f_over_s((u(:,:,:,1) - u(:,:,:,2)) / 2, G, GT, xi, eta, &
                                         npol, nsamp, nodes, element_type, axial) &
                              - grad_buff2(:,:,:,1) / 2d0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function strain_dipole(u, G, GT, xi, eta, npol, nodes, element_type, axial)
  ! Computes the strain tensor for displacement u excited bz a dipole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: u(0:npol,0:npol, 3)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: strain_dipole(0:npol,0:npol,6)
  
  real(kind=dp)                 :: grad_buff1(0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff2(0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff3(0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,1), G, GT, xi, eta, npol, nodes, &
                               element_type)
  
  ! 1: dsup, 2: dzup
  grad_buff2 = axisym_gradient(u(:,:,2), G, GT, xi, eta, npol, nodes, &
                               element_type)

  ! 1: dsuz, 2: dzuz
  grad_buff3 = axisym_gradient(u(:,:,3), G, GT, xi, eta, npol, nodes, element_type)

  strain_dipole(:,:,1) = grad_buff1(:,:,1)
  strain_dipole(:,:,2) = f_over_s(u(:,:,1) - u(:,:,2), G, GT, xi, eta, npol, nodes, &
                                  element_type, axial)
  strain_dipole(:,:,3) = grad_buff3(:,:,2)
  strain_dipole(:,:,4) = - 0.5d0 * (f_over_s(u(:,:,3), G, GT, xi, eta, npol, nodes, &
                                             element_type, axial) &
                                    + grad_buff2(:,:,2))
  strain_dipole(:,:,5) = (grad_buff1(:,:,2) + grad_buff3(:,:,1)) / 2d0
  strain_dipole(:,:,6) = - f_over_s((u(:,:,1) - u(:,:,2)) / 2, G, GT, xi, eta, npol, nodes, &
                                  element_type, axial) &
                         - grad_buff2(:,:,1) / 2d0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function straintrace_dipole_td(u, G, GT, xi, eta, npol, nsamp, nodes, element_type, axial)
  ! Computes the straintrace tensor for displacement u excited by a dipole source
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol, 3)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: straintrace_dipole_td(1:nsamp,0:npol,0:npol)
  
  real(kind=dp)                 :: grad_buff1(1:nsamp,0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff3(1:nsamp,0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,:,1), G, GT, xi, eta, npol, nsamp, &
                               nodes, element_type)
  
  ! 1: dsuz, 2: dzuz
  grad_buff3 = axisym_gradient(u(:,:,:,3), G, GT, xi, eta, npol, nsamp, nodes, &
                               element_type)

  straintrace_dipole_td(:,:,:) = grad_buff1(:,:,:,1) &
                                 + f_over_s(u(:,:,:,1) - u(:,:,:,2), G, GT, xi, eta, &
                                            npol, nsamp, nodes, element_type, axial) &
                                 + grad_buff3(:,:,:,2)

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function straintrace_dipole(u, G, GT, xi, eta, npol, nodes, element_type, axial)
  ! Computes the straintrace tensor for displacement u excited by a dipole source
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: u(0:npol,0:npol, 3)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: straintrace_dipole(0:npol,0:npol)
  
  real(kind=dp)                 :: grad_buff1(0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff3(0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,1), G, GT, xi, eta, npol, nodes, &
                               element_type)
  
  ! 1: dsuz, 2: dzuz
  grad_buff3 = axisym_gradient(u(:,:,3), G, GT, xi, eta, npol, nodes, element_type)

  straintrace_dipole(:,:) = grad_buff1(:,:,1) &
                            + f_over_s(u(:,:,1) - u(:,:,2), G, GT, xi, eta, npol, nodes, &
                                       element_type, axial) &
                            + grad_buff3(:,:,2)

end function
!-----------------------------------------------------------------------------------------

!--QUADRUPOLE-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function strain_quadpole_td(u, G, GT, xi, eta, npol, nsamp, nodes, element_type, axial)
  ! Computes the strain tensor for displacement u excited bz a quadpole_td source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol, 3)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: strain_quadpole_td(1:nsamp,0:npol,0:npol,6)
  
  real(kind=dp)                 :: grad_buff1(1:nsamp,0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff2(1:nsamp,0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff3(1:nsamp,0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,:,1), G, GT, xi, eta, npol, nsamp, nodes, &
                               element_type)
  
  ! 1: dsup, 2: dzup
  grad_buff2 = axisym_gradient(u(:,:,:,2), G, GT, xi, eta, npol, nsamp, nodes, &
                               element_type)

  ! 1: dsuz, 2: dzuz
  grad_buff3 = axisym_gradient(u(:,:,:,3), G, GT, xi, eta, npol, nsamp, nodes, &
                               element_type)

  strain_quadpole_td(:,:,:,1) = grad_buff1(:,:,:,1)
  strain_quadpole_td(:,:,:,2) = f_over_s(u(:,:,:,1) - 2 * u(:,:,:,2), G, GT, xi, eta, &
                                         npol, nsamp, nodes, element_type, axial)
  strain_quadpole_td(:,:,:,3) = grad_buff3(:,:,:,2)
  strain_quadpole_td(:,:,:,4) = - f_over_s(u(:,:,:,3), G, GT, xi, eta, npol, nsamp, &
                                           nodes, element_type, axial) &
                                - grad_buff2(:,:,:,2) / 2d0
  strain_quadpole_td(:,:,:,5) = (grad_buff1(:,:,:,2) + grad_buff3(:,:,:,1)) / 2d0
  strain_quadpole_td(:,:,:,6) = f_over_s(0.5d0 * u(:,:,:,2) - u(:,:,:,1), G, GT, xi, &
                                         eta, npol, nsamp, nodes, element_type, axial) &
                                - grad_buff2(:,:,:,1) / 2d0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function strain_quadpole(u, G, GT, xi, eta, npol, nodes, element_type, axial)
  ! Computes the strain tensor for displacement u excited bz a quadpole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: u(0:npol,0:npol, 3)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: strain_quadpole(0:npol,0:npol,6)
  
  real(kind=dp)                 :: grad_buff1(0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff2(0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff3(0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,1), G, GT, xi, eta, npol, nodes, element_type)
  
  ! 1: dsup, 2: dzup
  grad_buff2 = axisym_gradient(u(:,:,2), G, GT, xi, eta, npol, nodes, element_type)

  ! 1: dsuz, 2: dzuz
  grad_buff3 = axisym_gradient(u(:,:,3), G, GT, xi, eta, npol, nodes, element_type)

  strain_quadpole(:,:,1) = grad_buff1(:,:,1)
  strain_quadpole(:,:,2) = f_over_s(u(:,:,1) - 2 * u(:,:,2), G, GT, xi, eta, npol, nodes, &
                                    element_type, axial)
  strain_quadpole(:,:,3) = grad_buff3(:,:,2)
  strain_quadpole(:,:,4) = - f_over_s(u(:,:,3), G, GT, xi, eta, npol, nodes, &
                                      element_type, axial) &
                           - grad_buff2(:,:,2) / 2d0
  strain_quadpole(:,:,5) = (grad_buff1(:,:,2) + grad_buff3(:,:,1)) / 2d0
  strain_quadpole(:,:,6) = f_over_s(0.5d0 * u(:,:,2) - u(:,:,1), G, GT, xi, eta, npol, &
                                    nodes, element_type, axial)  &
                           - grad_buff2(:,:,1) / 2d0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function straintrace_quadpole_td(u, G, GT, xi, eta, npol, nsamp, nodes, element_type, axial)
  ! Computes the straintrace tensor for displacement u excited bz a quadpole_td source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: u(1:nsamp,0:npol,0:npol, 3)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: straintrace_quadpole_td(1:nsamp,0:npol,0:npol)
  
  real(kind=dp)                 :: grad_buff1(1:nsamp,0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff3(1:nsamp,0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,:,1), G, GT, xi, eta, npol, nsamp, nodes, &
                               element_type)
  
  ! 1: dsuz, 2: dzuz
  grad_buff3 = axisym_gradient(u(:,:,:,3), G, GT, xi, eta, npol, nsamp, nodes, &
                               element_type)

  straintrace_quadpole_td(:,:,:) &
      =   grad_buff1(:,:,:,1) &
        + f_over_s(u(:,:,:,1) - 2 * u(:,:,:,2), G, GT, xi, eta, npol, nsamp, nodes, &
                   element_type, axial) &
        + grad_buff3(:,:,:,2)

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function straintrace_quadpole(u, G, GT, xi, eta, npol, nodes, element_type, axial)
  ! Computes the straintrace tensor for displacement u excited bz a quadpole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: u(0:npol,0:npol, 3)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: straintrace_quadpole(0:npol,0:npol)
  
  real(kind=dp)                 :: grad_buff1(0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff3(0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,1), G, GT, xi, eta, npol, nodes, element_type)

  ! 1: dsuz, 2: dzuz
  grad_buff3 = axisym_gradient(u(:,:,3), G, GT, xi, eta, npol, nodes, element_type)

  straintrace_quadpole(:,:) &
      = grad_buff1(:,:,1) &
        + f_over_s(u(:,:,1) - 2 *u(:,:,2), G, GT, xi, eta, npol, nodes, &
                   element_type, axial) &
        + grad_buff3(:,:,2)

end function
!-----------------------------------------------------------------------------------------

!--GENERAL--------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function f_over_s_td(f, G, GT, xi, eta, npol, nsamp, nodes, element_type, axial)
  ! Computes the f / s
  ! needs G and GT for l'hospitals rule to compute f/s = df/ds at the axis s = 0
  
  use finite_elem_mapping, only : mapping
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: f(nsamp, 0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: f_over_s_td(nsamp, 0:npol,0:npol)

  integer                       :: ipol, jpol
  real(kind=dp)                 :: sz(0:npol,0:npol,1:2)

  do jpol=0, npol
     do ipol=0, npol
        sz(ipol, jpol,:) =  mapping(xi(ipol), eta(jpol), nodes, element_type)
     enddo
  enddo

  if (.not. axial) then
     do jpol=0, npol
        do ipol=0, npol
           f_over_s_td(:,ipol,jpol) = f(:,ipol,jpol) / sz(ipol,jpol,1)
        enddo
     enddo
  else
     do jpol=0, npol
        do ipol=1, npol
           f_over_s_td(:,ipol,jpol) = f(:,ipol,jpol) / sz(ipol,jpol,1)
        enddo
     enddo
     f_over_s_td(:,0,:) = dsdf_axis(f, G, GT, xi, eta, npol, nsamp, nodes, element_type)
  endif

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function f_over_s(f, G, GT, xi, eta, npol, nodes, element_type, axial)
  ! Computes the f / s
  ! needs G and GT for l'hospitals rule to compute f/s = df/ds at the axis s = 0
  
  use finite_elem_mapping, only : mapping
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: f(0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  logical, intent(in)           :: axial
  real(kind=dp)                 :: f_over_s(0:npol,0:npol)

  integer                       :: ipol, jpol
  real(kind=dp)                 :: sz(0:npol,0:npol,1:2)

  do ipol=0, npol
     do jpol=0, npol
        sz(ipol, jpol,:) =  mapping(xi(ipol), eta(jpol), nodes, element_type)
     enddo
  enddo

  if (.not. axial) then
     f_over_s(:,:) = f(:,:) / sz(:,:,1)
  else
     f_over_s(1:npol,0:npol) = f(1:npol,0:npol) / sz(1:npol,0:npol,1)
     f_over_s(0,:) = dsdf_axis(f, G, GT, xi, eta, npol, nodes, element_type)
  endif

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function dsdf_axis_td(f, G, GT, xi, eta, npol, nsamp, nodes, element_type)
  ! Computes the axisymmetric gradient of scalar field f
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: f(1:nsamp,0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: dsdf_axis_td(1:nsamp,0:npol)

  real(kind=dp)                 :: inv_j_npol(0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm_ipol0_1(1:nsamp,0:npol)
  real(kind=dp)                 :: mxm_ipol0_2(1:nsamp,0:npol)

  ipol = 0
  do jpol = 0, npol
     inv_j_npol(jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
  enddo

  mxm_ipol0_1 = mxm_ipol0(GT,f)
  mxm_ipol0_2 = mxm_ipol0(f,G)

  do jpol = 0, npol
     dsdf_axis_td(:,jpol) =   inv_j_npol(jpol,1,1) * mxm_ipol0_1(:,jpol) &
                            + inv_j_npol(jpol,2,1) * mxm_ipol0_2(:,jpol)
  enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function dsdf_axis(f, G, GT, xi, eta, npol, nodes, element_type)
  ! Computes the partial derivative of scalar field f for ipol = 0
  ! needed for l'hospitals rule to compute f/s = df/ds at the axis s = 0
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: f(0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: dsdf_axis(0:npol)

  real(kind=dp)                 :: inv_j_npol(0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm_ipol0_1(0:npol)
  real(kind=dp)                 :: mxm_ipol0_2(0:npol)

  ipol = 0
  do jpol = 0, npol
     inv_j_npol(jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
  enddo

  mxm_ipol0_1 = mxm_ipol0(GT,f)
  mxm_ipol0_2 = mxm_ipol0(f,G)
  dsdf_axis(:) =   inv_j_npol(:,1,1) * mxm_ipol0_1(:) &
                 + inv_j_npol(:,2,1) * mxm_ipol0_2(:)

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function axisym_gradient_td(f, G, GT, xi, eta, npol, nsamp, nodes, element_type)
  ! Computes the axisymmetric gradient of scalar field f
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: f(1:nsamp,0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: axisym_gradient_td(1:nsamp,0:npol,0:npol,1:2)

  real(kind=dp)                 :: inv_j_npol(0:npol,0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm1(1:nsamp,0:npol,0:npol)
  real(kind=dp)                 :: mxm2(1:nsamp,0:npol,0:npol)

  do ipol = 0, npol
     do jpol = 0, npol
        inv_j_npol(ipol,jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
     enddo
  enddo

!        | dxi  / ds  dxi  / dz |
! J^-1 = |                      |
!        | deta / ds  deta / dz |

  mxm1 = mxm(GT,f)
  mxm2 = mxm(f,G)

  do jpol = 0, npol
     do ipol = 0, npol
        axisym_gradient_td(:,ipol,jpol,1) =   &
                inv_j_npol(ipol,jpol,1,1) * mxm1(:,ipol,jpol) &
              + inv_j_npol(ipol,jpol,2,1) * mxm2(:,ipol,jpol)
        axisym_gradient_td(:,ipol,jpol,2) =   &
                inv_j_npol(ipol,jpol,1,2) * mxm1(:,ipol,jpol) &
              + inv_j_npol(ipol,jpol,2,2) * mxm2(:,ipol,jpol)
     enddo
  enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function axisym_gradient(f, G, GT, xi, eta, npol, nodes, element_type)
  ! Computes the axisymmetric gradient of scalar field f
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: f(0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: axisym_gradient(0:npol,0:npol,1:2)

  real(kind=dp)                 :: inv_j_npol(0:npol,0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm1(0:npol,0:npol)
  real(kind=dp)                 :: mxm2(0:npol,0:npol)

  do ipol = 0, npol
     do jpol = 0, npol
        inv_j_npol(ipol,jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
     enddo
  enddo

!        | dxi  / ds  dxi  / dz |
! J^-1 = |                      |
!        | deta / ds  deta / dz |

  mxm1 = mxm(GT,f)
  mxm2 = mxm(f,G)
  axisym_gradient(:,:,1) =   inv_j_npol(:,:,1,1) * mxm1 &
                           + inv_j_npol(:,:,2,1) * mxm2
  axisym_gradient(:,:,2) =   inv_j_npol(:,:,1,2) * mxm1 &
                           + inv_j_npol(:,:,2,2) * mxm2

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! a is time dependent
pure function mxm_atd(a, b)

  real(kind=dp), intent(in)  :: a(1:,0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_atd(1:size(a,1), 0:size(a,2)-1,0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_atd = 0

  do j = 0, size(b,2) -1
     do i = 0, size(a,2) -1
        do k = 0, size(a,3) -1
           mxm_atd(:,i,j) = mxm_atd(:,i,j) + a(:,i,k) * b(k,j)
        enddo
     end do
  end do 

end function mxm_atd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! b is time dependent
pure function mxm_btd(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(1:,0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_btd(1:size(b,1),0:size(a,1)-1,0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_btd = 0

  do j = 0, size(b,2) -1
     do i = 0, size(a,1) -1
        do k = 0, size(a,2) -1
           mxm_btd(:,i,j) = mxm_btd(:,i,j) + a(i,k) * b(:,k,j)
        enddo
     end do
  end do 

end function mxm_btd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
pure function mxm(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm(0:size(a,1)-1,0:size(b,2)-1)    !< Result
  integer                    :: i, j

  do j = 0, size(b,2) -1
     do i = 0, size(a,1) -1
        mxm(i,j) = sum(a(i,:) * b(:,j))
     end do
  end do 

end function mxm
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c, but only computes the component (0,:) of c
!! a is time dependent
pure function mxm_ipol0_atd(a, b)

  real(kind=dp), intent(in)  :: a(1:,0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_ipol0_atd(1:size(a,1), 0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_ipol0_atd = 0
  i = 0

  do j = 0, size(b,2) -1
     do k = 0, size(a,3) -1
        mxm_ipol0_atd(:,j) = mxm_ipol0_atd(:,j) + a(:,i,k) * b(k,j)
     enddo
  end do 

end function mxm_ipol0_atd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c, but only computes the component (0,:) of c
!! b is time dependent
pure function mxm_ipol0_btd(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(1:,0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_ipol0_btd(1:size(b,1),0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_ipol0_btd = 0

  i = 0
  do j = 0, size(b,2) -1
     do k = 0, size(a,2) -1
        mxm_ipol0_btd(:,j) = mxm_ipol0_btd(:,j) + a(i,k) * b(:,k,j)
     enddo
  end do 

end function mxm_ipol0_btd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c, but only computes the component (0,:) of c
pure function mxm_ipol0(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_ipol0(0:size(b,2)-1)    !< Result
  integer                    :: i, j

  i = 0
  do j = 0, size(b,2) -1
     mxm_ipol0(j) = sum(a(i,:) * b(:,j))
  end do 

end function mxm_ipol0
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
