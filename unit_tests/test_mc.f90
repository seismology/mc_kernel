module test_montecarlo

  use montecarlo, only: integrated_type
  use ftnunit, only   : test, assert_comparable_real1d
  implicit none
  public
contains

subroutine test_unit_hexagon
    type(integrated_type)                    :: mc_integral 
    integer                                   :: ipoint, iiter
    real(kind=8), dimension(100,3)            :: coords
    real(kind=8), dimension(100,1)            :: values
    real(kind=4), dimension(1)                :: integral
    real(kind=8), parameter                   :: pitominusthreehalf = 0.063493635934240969389d0
    real(kind=8), parameter, dimension(2)     :: bounds = [-0.2d0, 0.2d0]
    real(kind=8)                              :: volume

    volume = (bounds(2) - bounds(1))**3

    print *, 'In test_unit_hexagon'

    ! Integrate a gaussian function with sigma = 1 over [-5,5]^3
    
    call mc_integral%initialize_montecarlo(1, volume, 1d-3)

    iiter = 0

    do while (.not.mc_integral%areallconverged()) ! Beginning of Monte Carlo loop
        iiter = iiter + 1
        print *, 'Iteration ', iiter 
        call random_number(coords)
        coords = coords * (bounds(2)-bounds(1)) + bounds(1)
        values(:,1) = exp(-sum(coords**2,2)/2) * pitominusthreehalf
        call mc_integral%check_montecarlo_integral(values)
        print *, 'Value of the integral:', mc_integral%getintegral(), &
                 '; Variance:', mc_integral%getvariance()
    end do
 

    integral = mc_integral%getintegral()
    call assert_comparable_real1d(integral, [1.0], 1e-5, 'Integral == 1')
end subroutine

end module
