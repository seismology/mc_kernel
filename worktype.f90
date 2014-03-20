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
     integer                    :: ntotal_kernel
     integer                    :: nelems_per_task
     integer                    :: nvertices
     integer                    :: nvertices_per_elem
     integer                    :: mpitype
     integer                    :: itask
     integer                    :: ielement_type !1-tet, 2-quad, 3-tri, 4-hex
     integer, allocatable       :: connectivity(:,:)
     real(kind=dp), allocatable :: vertices(:,:)
     real(kind=dp), allocatable :: kernel_values(:,:,:)
     real(kind=dp), allocatable :: kernel_errors(:,:,:)
     integer, allocatable       :: niterations(:,:)
  end type

  type(work_type) :: wt

contains

!-----------------------------------------------------------------------------------------
subroutine init_work_type(nkernel, nelems_per_task, nvertices, nvertices_per_elem)

  use mpi
  
  integer, intent(in)   :: nkernel, nelems_per_task, nvertices, nvertices_per_elem
  integer               :: ierr, i
  integer, allocatable  :: oldtypes(:), blocklengths(:)
  integer(kind=MPI_ADDRESS_KIND), allocatable  :: offsets(:)
  integer, parameter    :: nblocks = 6
  character(len=64)     :: fmtstring

  wt%ntotal_kernel      = nkernel
  wt%nelems_per_task    = nelems_per_task
  wt%nvertices_per_elem = nvertices_per_elem
  wt%nvertices          = nvertices

  fmtstring = '(A32, I5)'
  write(lu_out, fmtstring), 'nkernel:', wt%ntotal_kernel
  write(lu_out, fmtstring), 'nelems_per_task:', wt%nelems_per_task   
  write(lu_out, fmtstring), 'nvertices_per_elem:', wt%nvertices_per_elem 
  write(lu_out, fmtstring), 'nvertices:', wt%nvertices          

  allocate(wt%connectivity(wt%nvertices_per_elem, wt%nelems_per_task))
  allocate(wt%vertices(3, wt%nvertices))
  allocate(wt%kernel_values(wt%ntotal_kernel, wt%nvertices_per_elem, wt%nelems_per_task))
  allocate(wt%kernel_errors(wt%ntotal_kernel, wt%nvertices_per_elem, wt%nelems_per_task))
  allocate(wt%niterations(wt%ntotal_kernel, wt%nelems_per_task))

  wt%connectivity  = 0
  wt%vertices      = 0
  wt%kernel_values = 0
  wt%kernel_errors = 0
  wt%niterations   = 0

  ! define blocks for the mpi type. NB: it seems to be necessary to define one
  ! block per array, otherwise having segfaults.
  allocate(oldtypes(nblocks))
  allocate(blocklengths(nblocks))
  allocate(offsets(nblocks))

  blocklengths(1) = 7 ! Variable sizes and itask
  blocklengths(2) = wt%nelems_per_task * wt%nvertices_per_elem ! connectivity
  blocklengths(3) = wt%nvertices * 3                  ! vertices
  blocklengths(4) = wt%ntotal_kernel * wt%nvertices_per_elem * wt%nelems_per_task !kernel_values
  blocklengths(5) = wt%ntotal_kernel * wt%nvertices_per_elem * wt%nelems_per_task !kernel_errors
  blocklengths(6) = wt%ntotal_kernel * wt%nelems_per_task                         !niterations

  oldtypes(1) = MPI_INTEGER            ! 
  oldtypes(2) = MPI_INTEGER            ! connectivity
  oldtypes(3) = MPI_DOUBLE_PRECISION   ! vertices 
  oldtypes(4) = MPI_DOUBLE_PRECISION   ! kernel_values 
  oldtypes(5) = MPI_DOUBLE_PRECISION   ! kernel_errors
  oldtypes(6) = MPI_INTEGER            ! kernel_errors

  ! find memory offsets, more stable then computing with MPI_TYPE_EXTEND
  call MPI_GET_ADDRESS(wt%ntotal_kernel, offsets(1), ierr)
  call MPI_GET_ADDRESS(wt%connectivity,  offsets(2), ierr)
  call MPI_GET_ADDRESS(wt%vertices,      offsets(3), ierr)
  call MPI_GET_ADDRESS(wt%kernel_values, offsets(4), ierr)
  call MPI_GET_ADDRESS(wt%kernel_errors, offsets(5), ierr)
  call MPI_GET_ADDRESS(wt%niterations,   offsets(6), ierr)

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
