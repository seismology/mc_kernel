!=========================================================================================
module test_montecarlo

  use global_parameters
  use montecarlo, only: integrated_type
  use ftnunit,    only: test, assert_comparable_real1d
  implicit none
  public

contains

!-----------------------------------------------------------------------------------------
subroutine test_mc_meanandvariance
    type(integrated_type)            :: mc_integral 
    integer, parameter               :: N = 10000
    real(kind=dp)                    :: values(N,1), volume, randa(N), randb(N)
    real(kind=sp)                    :: variance(1), integral(1), mean_an, variance_an 

    volume = 1.0
    

    ! Generate values with mean=5 and std=0
    
    values(:,1) = 5.0

    call mc_integral%initialize_montecarlo(1, volume, 1e-3)
    call mc_integral%check_montecarlo_integral(values)

    integral = mc_integral%getintegral()
    call assert_comparable_real1d(integral, [5.0], 1.e-7, 'Mean == 5.0')
    variance = mc_integral%getvariance()
    call assert_comparable_real1d(variance, [0.0], 1.e-7, 'Variance == 0.0')

    call mc_integral%freeme()



    ! Generate normal distributed values with mean=1 and std=2
    
    call random_number(randa)
    call random_number(randb)
    ! Produces normal distributed values with mean 1 and std 2
    values(:,1) = sqrt(-2*log(randa)) * cos(2*pi*randb) * 2.0 + 1.0
    
    call mc_integral%initialize_montecarlo(1, volume, 1e-3)
    call mc_integral%check_montecarlo_integral(values)

    mean_an = sum(values) / N
    variance_an = sum((values-mean_an)**2) / (N-1)
    !print *, 'Analytical mean     : ', mean_an
    !print *, 'Analytical variance : ', variance_an

    integral = mc_integral%getintegral()
    call assert_comparable_real1d(integral, [mean_an], 1.e-4, 'Mean == 1.0')
    variance = mc_integral%getvariance()
    call assert_comparable_real1d(sqrt(variance), [sqrt(variance_an/N)], 1.e-4, 'Error == 2/sqrt(N)')

    call mc_integral%freeme()

end subroutine test_mc_meanandvariance
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mc_unit_hexagon
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

    call mc_integral%initialize_montecarlo(1, volume, 5e-4)

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
 

    integral = mc_integral%getintegral()
    call assert_comparable_real1d(1+integral, 1+[3.98343E-003], 1e-4, 'Integral == reference value')
end subroutine test_mc_unit_hexagon
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_mc_sphere_in_tetrahedron
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

    call mc_integral%initialize_montecarlo(1, volume, 2e-4)

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
        !         '; std. dev:', sqrt(mc_integral%getvariance())
    end do
 

    integral = mc_integral%getintegral()
    call assert_comparable_real1d(1+integral, 1+[pi*(1./3.)**(1.5)*real(volume)], &
                                  1e-3, 'Integral == 1.10076')
end subroutine test_mc_sphere_in_tetrahedron
!-----------------------------------------------------------------------------------------


end module
!=========================================================================================
