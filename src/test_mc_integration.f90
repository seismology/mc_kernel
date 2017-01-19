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
module test_mc_integration

  use global_parameters
  use mc_integration, only: integrated_type
  use ftnunit,        only: assert_comparable_real1d, assert_true, assert_false, &
                            assert_equal
  implicit none
  public

contains

!-----------------------------------------------------------------------------------------
subroutine test_mc_meanandvariance
    type(integrated_type)            :: mc_integral 
    integer, parameter               :: N = 10000
    real(kind=dp)                    :: values(N,1), volume, randa(N), randb(N)
    real(kind=sp)                    :: variance(1), integral(1), mean_an, variance_an 


    ! Test with discrete set of values, where we know the analytical solution
    volume = 10
    values(1:5, 1) = [2, 2, -3, 1, 8]
    
    call mc_integral%initialize_montecarlo(1, volume, 1d-3)
    call mc_integral%check_montecarlo_integral(values(1:5, 1:1))

    integral = real(mc_integral%getintegral(), kind=sp)
    call assert_comparable_real1d(integral, [20.0], 1.e-7, 'Integral == 20.0')
    variance = real(mc_integral%getvariance(), kind=sp)
    call assert_comparable_real1d(variance, [310.0], 1.e-7, 'Integral error == 310')

    call mc_integral%freeme()


    ! Generate set of values with mean=5 and std=0
    
    volume = 1.0
    values(:,1) = 5.0

    call mc_integral%initialize_montecarlo(1, volume, 1d-3)
    call mc_integral%check_montecarlo_integral(values)

    integral = real(mc_integral%getintegral(), kind=sp)
    call assert_comparable_real1d(integral, [5.0], 1.e-7, 'Mean == 5.0')
    variance = real(mc_integral%getvariance(), kind=sp)
    call assert_comparable_real1d(variance, [0.0], 1.e-7, 'Variance == 0.0')

    call mc_integral%freeme()



    ! Generate normal distributed values with mean=1 and std=2
    volume = 1.0
    call random_number(randa)
    call random_number(randb)
    values(:,1) = sqrt(-2*log(randa)) * cos(2*pi*randb) * 2.0 + 1.0
    
    call mc_integral%initialize_montecarlo(nfuncs=1, &
                                           volume=volume, &
                                           allowed_error=1d-3, &
                                           allowed_relerror=1d-2)
    call mc_integral%check_montecarlo_integral(values)

    mean_an = real(sum(values) / N, kind=sp)
    variance_an = real(sum((values-mean_an)**2) / (N-1), kind=sp)

    integral = real(mc_integral%getintegral(), kind=sp)
    call assert_comparable_real1d(integral, [mean_an], 1.e-4, 'Mean == 1.0')
    variance = real(mc_integral%getvariance(), kind=sp)
    call assert_comparable_real1d(sqrt(variance), [sqrt(variance_an/N)], &
                                  1.e-4, 'Error == 2/sqrt(N)')

    call mc_integral%freeme()

end subroutine test_mc_meanandvariance
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mc_isconverged()
    type(integrated_type)            :: mc_integral 
    integer, parameter               :: N = 10000
    real(kind=dp)                    :: values(N,1), volume, randa(N), randb(N)
    real(kind=dp)                    :: mean_an, variance_an 
    integer                          :: i
    real(kind=dp)                    :: allowed_error 

    volume = 1.0
    allowed_error = 1d-1

    ! Produces normal distributed values with mean 1 and std 2
    call random_number(randa)
    call random_number(randb)
    values(:,1) = sqrt(-2*log(randa)) * cos(2*pi*randb) * 2.0 + 1.0
    
    call mc_integral%initialize_montecarlo(nfuncs=1, &
                                           volume=volume, &
                                           allowed_error=allowed_error)
    variance_an = 1e10

    i = 1
    do while (variance_an > allowed_error**2)                                    
      call assert_false(mc_integral%isconverged(1), 'recognized as not converged')
      call mc_integral%check_montecarlo_integral(reshape(values(i, :), [1, 1]))
      if (i>1) then
        mean_an = sum(values(1:i, 1)) / i
        variance_an = sum((values(1:i, 1)-mean_an)**2) / (i * (i-1.d0))
        ! print *, i, mean_an, mc_integral%getintegral(), variance_an, &
        !          mc_integral%getvariance(), mc_integral%isconverged(1)
      end if
      i = i + 1
    end do
    call assert_true(mc_integral%isconverged(1), 'recognized as converged')

    call mc_integral%freeme()

end subroutine test_mc_isconverged
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mc_areallconverged()
    type(integrated_type)            :: mc_integral 
    integer, parameter               :: N = 10000, nfuncs = 10
    real(kind=dp)                    :: values(N,nfuncs), randa(N, nfuncs), randb(N, nfuncs)
    real(kind=dp)                    :: mean_an(nfuncs), variance_an(nfuncs) 
    integer                          :: i, j, nconverged
    real(kind=dp)                    :: allowed_error, volume
    character(len=80)                :: teststr
    logical                          :: conv(nfuncs)

    volume = 1.0
    allowed_error = 1d-1
    conv = .false.

    ! Produces normal distributed values with mean 1 and std 2
    call random_number(randa)
    call random_number(randb)
    values = sqrt(-2*log(randa)) * cos(2*pi*randb) * 2.0 + 1.0
    
    call mc_integral%initialize_montecarlo(nfuncs=nfuncs, &
                                           volume=volume, &
                                           allowed_error=allowed_error)
    variance_an = 1e10

    i = 1
    do while (any(variance_an > allowed_error**2))                                    
      call mc_integral%check_montecarlo_integral(values(i:i, :))

      if (i>2) then
        mean_an = sum(values(1:i, :), dim=1) / i
        
        do j = 1, nfuncs
          if (.not.conv(j)) then
            variance_an(j) = sum((values(1:i, j) - mean_an(j))**2) / (i * (i-1))
            if (variance_an(j) <= allowed_error**2) conv(j) = .true.
          end if
        end do
        nconverged = count(conv) 
        call assert_equal(mc_integral%countconverged(), nconverged, &
                          'count converged kernels')
        
        ! do j = 1, nfuncs
        !   print *, i, j, mean_an(j), mc_integral%getintegral(j), variance_an(j), &
        !            mc_integral%getvariance(j), mc_integral%isconverged(j)
        !   if (mc_integral%isconverged(j).neqv.conv(j)) then
        !     write(teststr, '("not equal", I2)') j
        !     write(*,*) variance_an(j), mc_integral%getvariance(j)
        !     call assert_true(.false., teststr)
        !   end if
        ! end do
      end if
      i = i + 1
    end do
    call assert_true(mc_integral%areallconverged(), 'recognized as all converged')

    call mc_integral%freeme()

end subroutine test_mc_areallconverged
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mc_unit_hexagon
    type(integrated_type)                      :: mc_integral 
    integer                                    :: iiter
    real(kind=dp), dimension(100,3)            :: coords
    real(kind=dp), dimension(100,1)            :: values
    real(kind=sp), dimension(1)                :: integral
    real(kind=dp), parameter                   :: pitominusthreehalf &
                                                     = 0.063493635934240969389d0
    real(kind=dp), parameter, dimension(2)     :: bounds = [-0.2d0, 0.2d0]
    real(kind=dp)                              :: volume

    volume = (bounds(2) - bounds(1))**3

    call mc_integral%initialize_montecarlo(1, volume, 5d-4)

    iiter = 0

    do while (.not.mc_integral%areallconverged()) ! Beginning of Monte Carlo loop
        iiter = iiter + 1
        !print *, 'Iteration ', iiter 
        call random_number(coords)
        coords = coords * (bounds(2)-bounds(1)) + bounds(1)
        values(:,1) = exp(-sum(coords**2,2)/2) * pitominusthreehalf
        call mc_integral%check_montecarlo_integral(values)
        !print *, 'Value of the integral:', mc_integral%getintegral(), &
        !         '; standard dev.:', sqrt(mc_integral%getvariance())
    end do
 

    integral = real(mc_integral%getintegral(), kind=sp)
    call assert_comparable_real1d(1+integral, 1+[3.98343E-003], 1e-4, 'Integral == reference value')
end subroutine test_mc_unit_hexagon
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mc_sphere_in_tetrahedron
    use polytopes,  only                        : generate_random_points_tet
    type(integrated_type)                      :: mc_integral 
    integer                                    :: iiter
    real(kind=dp), dimension(100,3)            :: coords
    real(kind=dp), dimension(100,1)            :: values
    real(kind=sp), dimension(1)                :: integral
    real(kind=dp)                              :: volume
    real(kind=dp), dimension(3,4)              :: vertices
    real(kind=sp), parameter                   :: pi = 3.141419265

    vertices(:,1) = [1, 0, 0]
    vertices(:,2) = [0, 1, 0]
    vertices(:,3) = [0, 0, 1]
    vertices(:,4) = [0, 0, 0]
    volume = 1./6.

    call mc_integral%initialize_montecarlo(1, volume, 2d-4)

    iiter = 0

    do while (.not.mc_integral%areallconverged()) ! Beginning of Monte Carlo loop
        iiter = iiter + 1
        !print *, 'Iteration ', iiter 
        coords = transpose(generate_random_points_tet(vertices, 100))
        values = 0

        ! Sphere with radius sqrt(1/3)
        where( sum(coords**2, dim=2) < (1./3.) ) values(:,1) = 1

        call mc_integral%check_montecarlo_integral(values)
        !print *, 'Value of the integral:', mc_integral%getintegral(), &
        !         '; std. dev:', sqrt(mc_integral%getvariance())
    end do
 

    integral = real(mc_integral%getintegral(), kind=sp)
    call assert_comparable_real1d(1+integral, 1+[pi*(1./3.)**(1.5)*real(volume)], &
                                  1e-3, 'Integral == 1.10076')
end subroutine test_mc_sphere_in_tetrahedron
!-----------------------------------------------------------------------------------------


end module
!=========================================================================================
