!=========================================================================================
module test_halton_sequence

   use global_parameters
   use halton_sequence
   use ftnunit

   implicit none
   public
contains

!-----------------------------------------------------------------------------------------
subroutine test_init_halton

  ! Just to see whether it crashes
  call init_halton(ndim_in = 2,       &
                   seed_in = [1, 1],  &
                   base_in = [2, 3],  &
                   leap_in = [1, 1] )

  call free_halton()

end subroutine test_init_halton
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_get_halton
  real(kind=dp)  :: halton_data(2, 5)
  real(kind=dp)  :: halton_ref(2, 5)

  halton_ref(:,1) = [1.d0/2.d0, 1.d0/3.d0] 
  halton_ref(:,2) = [1.d0/4.d0, 2.d0/3.d0] 
  halton_ref(:,3) = [3.d0/4.d0, 1.d0/9.d0] 
  halton_ref(:,4) = [1.d0/8.d0, 4.d0/9.d0] 
  halton_ref(:,5) = [5.d0/8.d0, 7.d0/9.d0]
   
  ! Just to see whether the automatic initialization works
  call get_halton(halton_data)
  call free_halton()

  ! Get the first 10 members of the Halton sequence (2,3) and check against reference 
  ! solution.

  call init_halton(ndim_in = 2,       &
                   seed_in = [1, 1],  &
                   base_in = [2, 3],  &
                   leap_in = [1, 1] )

  call get_halton(halton_data)

  call assert_comparable(halton_data(1,:), halton_ref(1,:), 1d-7, &
                         'First 5 members of Halton sequence 2 are correct')

  call assert_comparable(halton_data(2,:), halton_ref(2,:), 1d-7, &
                         'First 5 members of Halton sequence 3 are correct')

end subroutine test_get_halton
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!subroutine test_halton_density
!  real(kind=dp)  :: halton_data(2, 5)
!
!  call init_halton(ndim_in = 2,       &
!                   seed_in = [1, 1],  &
!                   base_in = [2, 3],  &
!                   leap_in = [1, 1] )
!
!  call get_halton(halton_data)
!
!  density_middle = count(abs(halton_data(1,:)-0.5)<0.1)
!
!
!end subroutine test_halton_density
!-----------------------------------------------------------------------------------------

end module test_halton_sequence
