!=========================================================================================
module slave_module
  implicit none
  private
  public :: slave
contains

!-----------------------------------------------------------------------------------------
subroutine slave()
  use mpi
  use master_slave_parameters
  use work_type_mod
  use slave_work
  
  integer   :: mpistatus(MPI_STATUS_SIZE), ierror

  call init_work()

  do while (.true.)
    ! Receive a message from the master
    call MPI_Recv(wt, 1, wt%mpitype, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, ierror)

    ! Check the tag of the received message. If no more work to do, exit loop
    ! and return to main programm
    if (mpistatus(MPI_TAG) == DIETAG) exit

    ! Do the work
    call work()

    ! Send the result back
    call MPI_Send(wt, 1, wt%mpitype, 0, 0, MPI_COMM_WORLD, ierror)
  enddo
end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
