!=========================================================================================
module work_type_mod
  use global_parameters
  implicit none
  private
  
  public :: init_work_type
  public :: wt          ! essentially implementing the singleton pattern by
                        ! providing a single variable of the private type work_type

  type work_type
     sequence           ! force the derived type to be stored contiguously
     ! can't make the init function a member of work_type, because then it is
     ! impossible to use the seqeuence keyword
     ! maybe can get around this by defining one block for each variable in the
     ! derived type and dropping the sequence statement
     integer :: ntotal_kernel
     integer :: ndimensions
     integer :: nvertices
     integer :: mpitype
     real(kind=dp), allocatable :: vertices(:,:)
     real(kind=dp), allocatable :: kernel_values(:,:)
  end type

  type(work_type) :: wt

contains

!-----------------------------------------------------------------------------------------
subroutine init_work_type(nkern, ndim, nverts)
! set up the type and the mpitype that is used for communication

  use mpi
  
  integer, intent(in)   :: ndim, nverts, nkern
  integer               :: ierr, i
  integer, allocatable  :: oldtypes(:), blocklengths(:)
  integer(kind=MPI_ADDRESS_KIND), allocatable  :: offsets(:)
  integer, parameter    :: nblocks = 3

  wt%ntotal_kernel = nkern
  wt%nvertices = nverts
  wt%ndimensions = ndim

  allocate(wt%vertices(1:wt%ndimensions, 1:wt%nvertices))
  allocate(wt%kernel_values(wt%ntotal_kernel, wt%nvertices))
  wt%vertices = 0
  wt%kernel_values = 0

  ! define blocks for the mpi type. NB: it seems to be necessary to define one
  ! block per array, otherwise having segfaults.
  allocate(oldtypes(nblocks))
  allocate(blocklengths(nblocks))
  allocate(offsets(nblocks))

  blocklengths(1) = 4
  blocklengths(2) = wt%ndimensions * wt%nvertices
  blocklengths(3) = wt%ntotal_kernel * wt%nvertices

  oldtypes(1) = MPI_INTEGER
  oldtypes(2) = MPI_DOUBLE_PRECISION
  oldtypes(3) = MPI_DOUBLE_PRECISION

  ! find memory offsets, more stable then computing with MPI_TYPE_EXTEND
  call MPI_GET_ADDRESS(wt%ntotal_kernel, offsets(1), ierr)
  call MPI_GET_ADDRESS(wt%vertices,      offsets(2), ierr)
  call MPI_GET_ADDRESS(wt%kernel_values, offsets(3), ierr)

  ! make offsets relative
  do i=2, size(offsets)
      offsets(i) = offsets(i) - offsets(1)
  end do
  offsets(1) = 0

  call MPI_TYPE_CREATE_STRUCT(nblocks, blocklengths, offsets, oldtypes, wt%mpitype, ierr)
  call MPI_TYPE_COMMIT(wt%mpitype, ierr)

end subroutine init_work_type
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
