!=========================================================================================
module test_montecarlo

  use global_parameters
  use montecarlo, only: integrated_type
  use ftnunit, only   : test, assert_comparable_real1d
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_unit_hexagon
    type(integrated_type)                    :: mc_integral 
    integer                                   :: ipoint, iiter
    real(kind=dp), dimension(100,3)            :: coords
    real(kind=dp), dimension(100,1)            :: values
    real(kind=sp), dimension(1)                :: integral
    real(kind=dp), parameter                   :: pitominusthreehalf &
                                                     = 0.063493635934240969389d0
    real(kind=dp), parameter, dimension(2)     :: bounds = [-0.2d0, 0.2d0]
    real(kind=dp)                              :: volume

    volume = (bounds(2) - bounds(1))**3

    !print *, 'In test_unit_hexagon'

    ! Integrate a gaussian function with sigma = 1 over [-5,5]^3
    
    call mc_integral%initialize_montecarlo(1, volume, 1d-3)

    iiter = 0

    do while (.not.mc_integral%areallconverged()) ! Beginning of Monte Carlo loop
        iiter = iiter + 1
        !print *, 'Iteration ', iiter 
        call random_number(coords)
        coords = coords * (bounds(2)-bounds(1)) + bounds(1)
        values(:,1) = exp(-sum(coords**2,2)/2) * pitominusthreehalf
        call mc_integral%check_montecarlo_integral(values)
        !print *, 'Value of the integral:', mc_integral%getintegral(), &
        !         '; Variance:', mc_integral%getvariance()
    end do
 

    integral = mc_integral%getintegral()
    call assert_comparable_real1d(integral, [3.98343E-003], 1e-3, 'Integral == 1')
end subroutine test_unit_hexagon
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_sphere_in_tetrahedron
    use tetrahedra, only                       : generate_random_point
    type(integrated_type)                     :: mc_integral 
    integer                                   :: ipoint, iiter
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

    !print *, 'In test_unit_hexagon'

    ! Integrate a gaussian function with sigma = 1 over [-5,5]^3
    
    call mc_integral%initialize_montecarlo(1, volume, 1d-3)

    iiter = 0

    do while (.not.mc_integral%areallconverged()) ! Beginning of Monte Carlo loop
        iiter = iiter + 1
        !print *, 'Iteration ', iiter 
        coords = transpose(generate_random_point(vertices, 100))
        values = 0

        ! Sphere with radius sqrt(1/3)
        where( sum(coords**2, dim=2) < (1./3.) ) values(:,1) = 1

        call mc_integral%check_montecarlo_integral(values)
        !print *, 'Value of the integral:', mc_integral%getintegral(), &
        !         '; Variance:', mc_integral%getvariance()
    end do
 

    integral = mc_integral%getintegral()
    call assert_comparable_real1d(integral, [pi*(1./3.)**(1.5)*real(volume)], &
                                  1e-2, 'Integral == 1')
end subroutine test_sphere_in_tetrahedron
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
