program struct

  use mpi
  
  implicit none

  type :: mytype
     sequence           ! force the derived type to be stored contiguously
     integer            :: myint
     integer            :: myint2
     double precision   :: mydp
     double precision   :: mydp2
  end type
  
  integer               :: ierr, mynum, nproc, sendrequest
  integer               :: mpistatus(MPI_STATUS_SIZE)
  integer, allocatable  :: oldtypes(:), blocklengths(:), offsets(:)
  integer               :: newtype
  integer, parameter    :: nblocks = 2

  type(mytype)          :: sbuf, rbuf
  
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc,  ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mynum, ierr)
 
  allocate(oldtypes(nblocks))
  allocate(blocklengths(nblocks))
  allocate(offsets(nblocks))

  blocklengths(1) = 2
  blocklengths(2) = 2

  oldtypes(1) = MPI_INTEGER
  oldtypes(2) = MPI_DOUBLE_PRECISION

  ! find memory offsets, more stable then computing with MPI_TYPE_EXTEND
  call MPI_ADDRESS(sbuf%myint, offsets(1), ierr)
  call MPI_ADDRESS(sbuf%mydp,  offsets(2), ierr)

  ! make relative
  offsets(2) = offsets(2) - offsets(1)
  offsets(1) = 0

  call MPI_TYPE_STRUCT(nblocks, blocklengths, offsets, oldtypes, newtype, ierr)
  call MPI_TYPE_COMMIT(newtype, ierr)

  if (mynum == 0) then
     sbuf%myint  = 10
     sbuf%myint2 = 20
     sbuf%mydp  = 1.1d0
     sbuf%mydp2 = 2.2d0
     write(6,*) sbuf

     call MPI_Send(sbuf, 1, newtype, 1, 1, MPI_COMM_WORLD, sendrequest, ierr)
  elseif (mynum == 1) then
     call MPI_Recv(rbuf, 1, newtype, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &
                   mpistatus, ierr)
      write(6,*) rbuf
  endif

  call MPI_FINALIZE(ierr)

end program
