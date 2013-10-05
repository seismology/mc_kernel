!=========================================================================================
module test_buffer

  use global_parameters
  use buffer
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------

subroutine test_buffer_storage
   type(buffer_type)                   :: buffer
   integer, parameter                  :: nbuffer = 100
   integer, parameter                  :: lenbuffer = 1000
   real(kind=sp), dimension(lenbuffer) :: stufftostore
   integer                             :: iwrite, status

   status = buffer%init(nbuffer, lenbuffer)
   call assert_equal(status, 0, 'Buffer initialisation succeeds')

   do iwrite = 1, nbuffer + 1
      call random_number(stufftostore)
      status = buffer%put(iwrite, stufftostore)
      call assert_equal(status, 0, 'Buffer put succeeds')
   end do

   status = buffer%freeme()
   call assert_equal(status, 0, 'Buffer freeme succeeds')

end subroutine

subroutine test_buffer_retrieval
   type(buffer_type)                   :: buffer
   integer, parameter                  :: nbuffer = 100
   integer, parameter                  :: lenbuffer = 1000
   real(kind=sp), dimension(lenbuffer) :: stufftostore, stuffretrieved
   integer                             :: iwrite, status

   status = buffer%init(nbuffer, lenbuffer)

   call random_number(stufftostore)

   ! After ten other writes to the buffer, the data should still be there
   status = buffer%put(1, stufftostore)
   call assert_equal(status, 0, 'Buffer put succeeds')
   do iwrite = 1, nbuffer/10
      status = buffer%put(iwrite, stufftostore*iwrite)
   end do
   status = buffer%get(1, stuffretrieved)
   call assert_equal(status, 0, 'Buffer get succeeds')
   call assert_comparable_real1d(stufftostore, stuffretrieved, &
                                 1e-7, 'Buffer gives data back (5% chance of random failure)')
  
   ! After 1000 other writes to the buffer, the data should be overwritten
   status = buffer%put(1, stufftostore)
   call assert_equal(status, 0, 'Buffer put succeeds')
   do iwrite = 1, nbuffer*10
      status = buffer%put(iwrite, stufftostore*iwrite)
   end do
   
   status = buffer%get(1, stuffretrieved)
   call assert_equal(status, -1, 'Buffer overwritten after reasonable time')

end subroutine test_buffer_retrieval


end module
