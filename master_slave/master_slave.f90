!=========================================================================================
program master_slave
  use mpi
  use master_module
  use slave_module
  use work_mod

  implicit none
  integer               :: myrank, ierror
  integer               :: nkern
  
  nkern = 10

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)

  call init_work_type(nkern, 3, 4)

  if (myrank == 0) then
     print *, 'MASTER'
     call master()
  else
     print *, 'SLAVE', myrank
     call slave()
  endif

  call MPI_FINALIZE(ierror)
end
!=========================================================================================
