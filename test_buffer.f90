!=========================================================================================
module test_buffer

  use global_parameters
  use buffer
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_buffer_storage_1d
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
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_buffer_storage_2d
   type(buffer_type)                   :: buffer
   integer, parameter                  :: nbuffer = 100
   integer, parameter                  :: lenbuffer1 = 1000, lenbuffer2 = 2
   real(kind=sp)                       :: stufftostore(lenbuffer1, lenbuffer2)
   integer                             :: iwrite, status

   status = buffer%init(nbuffer, lenbuffer1, lenbuffer2)
   call assert_equal(status, 0, 'Buffer initialisation succeeds')

   do iwrite = 1, nbuffer + 1
      call random_number(stufftostore)
      status = buffer%put(iwrite, stufftostore)
      call assert_equal(status, 0, 'Buffer put succeeds')
   end do

   status = buffer%freeme()
   call assert_equal(status, 0, 'Buffer freeme succeeds')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_buffer_storage_3d
   type(buffer_type)                   :: buffer
   integer, parameter                  :: nbuffer = 100
   integer, parameter                  :: lenbuffer1 = 1000, lenbuffer2 = 2, &
                                          lenbuffer3 = 2
   real(kind=sp)                       :: stufftostore(lenbuffer1, lenbuffer2, lenbuffer3)
   integer                             :: iwrite, status

   status = buffer%init(nbuffer, lenbuffer1, lenbuffer2, lenbuffer3)
   call assert_equal(status, 0, 'Buffer initialisation succeeds')

   do iwrite = 1, nbuffer + 1
      call random_number(stufftostore)
      status = buffer%put(iwrite, stufftostore)
      call assert_equal(status, 0, 'Buffer put succeeds')
   end do

   status = buffer%freeme()
   call assert_equal(status, 0, 'Buffer freeme succeeds')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_buffer_storage_4d
   type(buffer_type)    :: buffer
   integer, parameter   :: nbuffer = 100
   integer, parameter   :: lenbuffer1 = 1000, lenbuffer2 = 2, &
                           lenbuffer3 = 2, lenbuffer4 = 2
   real(kind=sp)        :: stufftostore(lenbuffer1, lenbuffer2, lenbuffer3, lenbuffer4)
   integer              :: iwrite, status

   status = buffer%init(nbuffer, lenbuffer1, lenbuffer2, lenbuffer3, lenbuffer4)
   call assert_equal(status, 0, 'Buffer initialisation succeeds')

   do iwrite = 1, nbuffer + 1
      call random_number(stufftostore)
      status = buffer%put(iwrite, stufftostore)
      call assert_equal(status, 0, 'Buffer put succeeds')
   end do

   status = buffer%freeme()
   call assert_equal(status, 0, 'Buffer freeme succeeds')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_buffer_retrieval_1d
   type(buffer_type)                   :: buffer
   integer, parameter                  :: nbuffer = 100
   integer, parameter                  :: lenbuffer = 1000
   real(kind=sp), dimension(lenbuffer) :: stufftostore, stuffretrieved
   integer                             :: iwrite, status

   status = buffer%init(nbuffer, lenbuffer)

   stuffretrieved = 0.0
   call random_number(stufftostore)

   ! After ten other writes to the buffer, the data should still be there
   status = buffer%put(1, stufftostore)
   call assert_equal(status, 0, 'Buffer put succeeds')
   do iwrite = 2, nbuffer/10
      status = buffer%put(iwrite, stufftostore*iwrite)
      call assert_equal(status, 0, 'Buffer put succeeds')
   end do
   status = buffer%get(1, stuffretrieved)
   call assert_equal(status, 0, 'Buffer get succeeds')
   call assert_comparable_real1d(stufftostore, stuffretrieved, &
                                 1e-7, 'Buffer gives correct data back')
   status = buffer%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_buffer_retrieval_2d
   type(buffer_type)                   :: buffer
   integer, parameter                  :: nbuffer = 100
   integer, parameter                  :: lenbuffer1 = 1000, lenbuffer2 = 2
   real(kind=sp)                       :: stufftostore(lenbuffer1, lenbuffer2), &
                                          stuffretrieved(lenbuffer1, lenbuffer2)
   integer                             :: iwrite, status, ibuff

   status = buffer%init(nbuffer, lenbuffer1, lenbuffer2)

   stuffretrieved = 0.0
   call random_number(stufftostore)

   ! After ten other writes to the buffer, the data should still be there
   status = buffer%put(1, stufftostore)
   call assert_equal(status, 0, 'Buffer put succeeds')

   do iwrite = 2, nbuffer/10
      status = buffer%put(iwrite, stufftostore*iwrite)
      call assert_equal(status, 0, 'Buffer put succeeds')
   end do

   status = buffer%get(1, stuffretrieved)
   call assert_equal(status, 0, 'Buffer get succeeds')
   do ibuff = 1, lenbuffer2
      call assert_comparable_real1d(stufftostore(:,ibuff), stuffretrieved(:,ibuff), &
                                    1e-7, 'Buffer gives correct data back')
   enddo

   status = buffer%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_buffer_retrieval_3d
   type(buffer_type)                   :: buffer
   integer, parameter                  :: nbuffer = 100
   integer, parameter                  :: lenbuffer1 = 1000, lenbuffer2 = 2, &
                                          lenbuffer3 = 2
   real(kind=sp)                       :: stufftostore(lenbuffer1, lenbuffer2, lenbuffer3), &
                                          stuffretrieved(lenbuffer1, lenbuffer2, lenbuffer3)
   integer                             :: iwrite, status, ibuff, jbuff

   status = buffer%init(nbuffer, lenbuffer1, lenbuffer2, lenbuffer3)

   stuffretrieved = 0.0
   call random_number(stufftostore)

   ! After ten other writes to the buffer, the data should still be there
   status = buffer%put(1, stufftostore)
   call assert_equal(status, 0, 'Buffer put succeeds')

   do iwrite = 2, nbuffer/10
      status = buffer%put(iwrite, stufftostore*iwrite)
      call assert_equal(status, 0, 'Buffer put succeeds')
   end do

   status = buffer%get(1, stuffretrieved)
   call assert_equal(status, 0, 'Buffer get succeeds')
   do ibuff = 1, lenbuffer2
      do jbuff = 1, lenbuffer3
         call assert_comparable_real1d(stufftostore(:,ibuff, jbuff), stuffretrieved(:,ibuff, jbuff), &
                                       1e-7, 'Buffer gives correct data back')
      enddo
   enddo

   status = buffer%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_buffer_retrieval_4d
   type(buffer_type)    :: buffer
   integer, parameter   :: nbuffer = 100
   integer, parameter   :: lenbuffer1 = 1000, lenbuffer2 = 2, &
                           lenbuffer3 = 2, lenbuffer4 = 2
   real(kind=sp)        :: stufftostore(lenbuffer1, lenbuffer2, lenbuffer3, lenbuffer4), &
                           stuffretrieved(lenbuffer1, lenbuffer2, lenbuffer3, lenbuffer4)
   integer              :: iwrite, status, ibuff, jbuff, kbuff

   status = buffer%init(nbuffer, lenbuffer1, lenbuffer2, lenbuffer3, lenbuffer4)

   stuffretrieved = 0.0
   call random_number(stufftostore)

   ! After ten other writes to the buffer, the data should still be there
   status = buffer%put(1, stufftostore)
   call assert_equal(status, 0, 'Buffer put succeeds')

   do iwrite = 2, nbuffer/10
      status = buffer%put(iwrite, stufftostore*iwrite)
      call assert_equal(status, 0, 'Buffer put succeeds')
   end do

   status = buffer%get(1, stuffretrieved)
   call assert_equal(status, 0, 'Buffer get succeeds')
   do ibuff = 1, lenbuffer2
      do jbuff = 1, lenbuffer3
         do kbuff = 1, lenbuffer4
            call assert_comparable_real1d(stufftostore(:,ibuff,jbuff,kbuff), &
                                          stuffretrieved(:,ibuff,jbuff,kbuff), &
                                          1e-7, 'Buffer gives correct data back')
         enddo
      enddo
   enddo

   status = buffer%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_buffer_overwrite
   type(buffer_type)                   :: buffer
   integer, parameter                  :: nbuffer = 100
   integer, parameter                  :: lenbuffer = 1000
   real(kind=sp), dimension(lenbuffer) :: stufftostore, stuffretrieved
   integer                             :: iwrite, status

   status = buffer%init(nbuffer, lenbuffer)

   stuffretrieved = 0.0
   call random_number(stufftostore)
  
   ! After 1000 other writes to the buffer, the data should be overwritten
   stuffretrieved = 0.0
   status = buffer%put(1, stufftostore)
   call assert_equal(status, 0, 'Buffer put succeeds')
   do iwrite = 2, nbuffer*10
      status = buffer%put(iwrite, stufftostore*iwrite)
      call assert_equal(status, 0, 'Buffer put succeeds')
   end do
   
   status = buffer%get(1, stuffretrieved)
   call assert_equal(status, -1, 'Buffer overwritten after reasonable time')

   status = buffer%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
