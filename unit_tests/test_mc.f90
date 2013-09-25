module test_mc

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

    print *, 'In test_unit_hexagon'

    ! Integrate a gaussian function with sigma = 1 over [-5,5]^3
    
    call mc_integral%initialize_montecarlo(1, 1000d0, 1d-5)

    iiter = 0

    do while (.not.mc_integral%areallconverged()) ! Beginning of Monte Carlo loop
        iiter = iiter + 1
        print *, 'Iteration ', iiter 
        call random_number(coords)
        coords = coords * 10 - 5
        values(:,1) = exp(-sum(coords**2,2)/2) * pitominusthreehalf
        call mc_integral%check_montecarlo_integral(values)
        print *, 'Value of the integral:', mc_integral%getintegral(), &
                 '; Variance:', mc_integral%getvariance()
    end do
 

    integral = mc_integral%getintegral()
    call assert_comparable_real1d(integral, [1.0], 1e-5, 'Integral == 1')
end subroutine

subroutine test_all
  call test(test_unit_hexagon, 'TEST MC_unit_hexagon')
end subroutine

end module
