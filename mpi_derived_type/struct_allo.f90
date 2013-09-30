program struct

  use mpi
  
  implicit none

  type :: mytype
     sequence           ! force the derived type to be stored contiguously
     integer            :: myint1
     integer            :: myint2
     integer            :: myint3
     integer            :: myint4
     double precision, allocatable   :: mydp1(:)
     double precision, allocatable   :: mydp2(:)
  end type
  
  integer               :: ierr, mynum, nproc, sendrequest, i
  integer               :: mpistatus(MPI_STATUS_SIZE)
  integer, allocatable  :: oldtypes(:), blocklengths(:)
  integer(kind=MPI_ADDRESS_KIND), allocatable :: offsets(:)
  integer               :: newtype
  integer, parameter    :: nblocks = 3, n = 5

  type(mytype)          :: buf
  
  allocate(buf%mydp1(n))
  allocate(buf%mydp2(n))

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc,  ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mynum, ierr)
 
  allocate(oldtypes(nblocks))
  allocate(blocklengths(nblocks))
  allocate(offsets(nblocks))

  blocklengths(1) = 4
  blocklengths(2) = n
  blocklengths(3) = n

  oldtypes(1) = MPI_INTEGER
  oldtypes(2) = MPI_DOUBLE_PRECISION
  oldtypes(3) = MPI_DOUBLE_PRECISION

  ! find memory offsets, more stable then computing with MPI_TYPE_EXTEND
  call MPI_GET_ADDRESS(buf%myint1,  offsets(1), ierr)
  call MPI_GET_ADDRESS(buf%mydp1,   offsets(2), ierr)
  call MPI_GET_ADDRESS(buf%mydp2,   offsets(3), ierr)

  ! make relative
  do i = 2, size(offsets)
     offsets(i) = offsets(i) - offsets(1)
  end do
  offsets(1) = 0

  write(6,*) offsets

  call MPI_TYPE_CREATE_STRUCT(nblocks, blocklengths, offsets, oldtypes, newtype, ierr)
  call MPI_TYPE_COMMIT(newtype, ierr)

  if (mynum == 0) then
     buf%myint1  = 1
     buf%myint2  = 2
     buf%myint3  = 3
     buf%myint4  = 4
     buf%mydp1  = 1.1d0
     buf%mydp2  = 2.1d0
     write(6,*) buf%mydp1
     write(6,*) buf%mydp2

     call MPI_Send(buf, 1, newtype, 1, 1, MPI_COMM_WORLD, sendrequest, ierr)
  elseif (mynum == 1) then
     call MPI_Recv(buf, 1, newtype, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &
                   mpistatus, ierr)
      write(6,*) buf%myint1
      write(6,*) buf%myint2
      write(6,*) buf%myint3
      write(6,*) buf%myint4
      write(6,*) buf%mydp1
      write(6,*) buf%mydp2
  endif

  call MPI_FINALIZE(ierr)

end program
