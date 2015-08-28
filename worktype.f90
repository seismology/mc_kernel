!=========================================================================================
module work_type_mod
  use global_parameters
  implicit none
  private
  
  public :: init_work_type
  public :: free_work_type
  public :: wt          ! essentially implementing the singleton pattern by
                        ! providing a single variable of the private type work_type

  type work_type
     sequence           ! force the derived type to be stored contiguously
     ! can't make the init function a member of work_type, because then it is
     ! impossible to use the seqeuence keyword
     ! maybe can get around this by defining one block for each variable in the
     ! derived type and dropping the sequence statement
     integer                    :: ntotal_kernel = -1
     integer                    :: nelems_per_task
     integer                    :: nvertices
     integer                    :: nvertices_per_elem
     integer                    :: nbasisfuncs_per_elem
     integer                    :: nmodel_parameters
     integer                    :: nmodel_parameters_hetero
     integer                    :: mpitype
     integer                    :: itask
     integer                    :: ielement_type !1-tet, 2-quad, 3-tri, 4-hex, 5-vox
     integer                    :: ndumps
     integer                    :: ndim
     integer                    :: dt   !< dt in microseconds
     integer, allocatable       :: connectivity(:,:)
     real(kind=dp), allocatable :: vertices(:,:)
     real(kind=dp), allocatable :: kernel_values(:,:,:)
     real(kind=dp), allocatable :: kernel_variance(:,:,:)
     integer, allocatable       :: niterations(:,:)
     real(kind=dp), allocatable :: computation_time(:)
     real(kind=dp), allocatable :: model(:,:,:)
     real(kind=dp), allocatable :: hetero_model(:,:,:)
     real(kind=dp), allocatable :: fw_field(:,:,:,:,:)
     real(kind=dp), allocatable :: bw_field(:,:,:,:,:)
     real(kind=dp), allocatable :: conv_field(:,:,:,:,:)

  end type

  type(work_type) :: wt

contains

!-----------------------------------------------------------------------------------------
subroutine init_work_type(nkernel, nelems_per_task, nvertices, nvertices_per_elem, &
                          nbasisfuncs_per_elem, nmodel_parameters, nmodel_parameters_hetero, &
                          plot_wavefields, ndumps, ndim, dt)

# ifndef include_mpi
  use mpi
# endif
  
# ifdef include_mpi
  include 'mpif.h'
# endif

  integer, intent(in)       :: nkernel, nelems_per_task, nvertices, nvertices_per_elem
  integer, intent(in)       :: nbasisfuncs_per_elem, nmodel_parameters, nmodel_parameters_hetero
  logical, intent(in)       :: plot_wavefields
  integer, intent(in)       :: ndumps, ndim
  real(kind=sp), intent(in) :: dt

  integer               :: ierr, i
  integer, allocatable  :: oldtypes(:), blocklengths(:)
  integer(kind=MPI_ADDRESS_KIND), allocatable  :: offsets(:)
  integer, parameter    :: nblocks = 12
  character(len=64)     :: fmtstring

  if (wt%ntotal_kernel .ne. -1) then
    print *, 'ERROR: Trying to init work type that is already initialized!'
    stop
  end if

  wt%ntotal_kernel             = nkernel
  wt%nelems_per_task           = nelems_per_task
  wt%nvertices                 = nvertices
  wt%nvertices_per_elem        = nvertices_per_elem
  wt%nbasisfuncs_per_elem      = nbasisfuncs_per_elem
  wt%nmodel_parameters         = nmodel_parameters
  wt%nmodel_parameters_hetero  = nmodel_parameters_hetero

  wt%ndumps                    = ndumps
  wt%ndim                      = ndim
  wt%dt                        = nint(dt*1.e6) !dt in microseconds

  fmtstring = '(A32, I5)'
  write(lu_out, fmtstring) 'nkernel:', wt%ntotal_kernel
  write(lu_out, fmtstring) 'nelems_per_task:', wt%nelems_per_task   
  write(lu_out, fmtstring) 'nvertices:', wt%nvertices          
  write(lu_out, fmtstring) 'nvertices_per_elem:', wt%nvertices_per_elem 
  write(lu_out, fmtstring) 'nbasisfuncs_per_elem:', wt%nbasisfuncs_per_elem 
  write(lu_out, fmtstring) 'nmodel_parameters:', wt%nmodel_parameters
  write(lu_out, fmtstring) 'nmodel_parameters_hetero:', wt%nmodel_parameters_hetero


  allocate(wt%connectivity(wt%nvertices_per_elem, wt%nelems_per_task))
  allocate(wt%vertices(3, wt%nvertices))
  allocate(wt%kernel_values(wt%ntotal_kernel, wt%nbasisfuncs_per_elem, wt%nelems_per_task))
  allocate(wt%kernel_variance(wt%ntotal_kernel, wt%nbasisfuncs_per_elem, wt%nelems_per_task))
  allocate(wt%niterations(wt%ntotal_kernel, wt%nelems_per_task))
  allocate(wt%computation_time(wt%nelems_per_task))
  allocate(wt%model(wt%nmodel_parameters, wt%nbasisfuncs_per_elem, wt%nelems_per_task))
  allocate(wt%hetero_model(wt%nmodel_parameters_hetero, wt%nbasisfuncs_per_elem, wt%nelems_per_task))

  if (plot_wavefields) then
    allocate(wt%fw_field(wt%ndumps, wt%ndim, wt%ntotal_kernel, wt%nbasisfuncs_per_elem, wt%nelems_per_task))
    allocate(wt%bw_field(wt%ndumps, wt%ndim, wt%ntotal_kernel, wt%nbasisfuncs_per_elem, wt%nelems_per_task))
    allocate(wt%conv_field(wt%ndumps, 1, wt%ntotal_kernel, wt%nbasisfuncs_per_elem, wt%nelems_per_task))
  else
    allocate(wt%fw_field(1, 1, 1, 1, 1))
    allocate(wt%bw_field(1, 1, 1, 1, 1))
    allocate(wt%conv_field(1, 1, 1, 1, 1))
  end if

  wt%connectivity    = 0
  wt%vertices        = 0
  wt%kernel_values   = 0
  wt%kernel_variance = 0
  wt%niterations     = 0
  wt%model           = 0
  wt%hetero_model    = 0
  wt%fw_field        = 0
  wt%bw_field        = 0
  wt%conv_field      = 0

  ! define blocks for the mpi type. NB: it seems to be necessary to define one
  ! block per array, otherwise having segfaults.
  allocate(oldtypes(nblocks))
  allocate(blocklengths(nblocks))
  allocate(offsets(nblocks))

  blocklengths(1) = 13 ! variable sizes and itask
  blocklengths(2) = wt%nelems_per_task * wt%nvertices_per_elem ! connectivity
  blocklengths(3) = wt%nvertices * 3                           ! vertices
  blocklengths(4) = wt%ntotal_kernel * wt%nbasisfuncs_per_elem * wt%nelems_per_task !kernel_values
  blocklengths(5) = wt%ntotal_kernel * wt%nbasisfuncs_per_elem * wt%nelems_per_task !kernel_variance
  blocklengths(6) = wt%ntotal_kernel * wt%nelems_per_task                           !niterations
  blocklengths(7) = wt%nelems_per_task                                              !computation_time
  blocklengths(8) = wt%nmodel_parameters * wt%nbasisfuncs_per_elem &
                    * wt%nelems_per_task                                            !model
  blocklengths(9) = wt%nmodel_parameters_hetero * wt%nbasisfuncs_per_elem &
                    * wt%nelems_per_task                                            !model

  if (plot_wavefields) then
    blocklengths(10) = wt%ndumps * wt%ndim * wt%ntotal_kernel * &       
                       wt%nbasisfuncs_per_elem * wt%nelems_per_task                 !Forward field
    blocklengths(11) = wt%ndumps * wt%ndim * wt%ntotal_kernel * &       
                       wt%nbasisfuncs_per_elem * wt%nelems_per_task                 !Backward field
    blocklengths(12) = wt%ndumps * 1 * wt%ntotal_kernel * &       
                       wt%nbasisfuncs_per_elem * wt%nelems_per_task                 !Convolved field
  else
    blocklengths(10:12) = 1
  end if

  oldtypes(1)  = MPI_INTEGER            ! all variable sizes and itask
  oldtypes(2)  = MPI_INTEGER            ! connectivity
  oldtypes(3)  = MPI_DOUBLE_PRECISION   ! vertices 
  oldtypes(4)  = MPI_DOUBLE_PRECISION   ! kernel_values 
  oldtypes(5)  = MPI_DOUBLE_PRECISION   ! kernel_variance
  oldtypes(6)  = MPI_INTEGER            ! niterations
  oldtypes(7)  = MPI_DOUBLE_PRECISION   ! computation_time
  oldtypes(8)  = MPI_DOUBLE_PRECISION   ! model 
  oldtypes(9)  = MPI_DOUBLE_PRECISION   ! heterogeneity model 
  oldtypes(10) = MPI_DOUBLE_PRECISION   ! forward field
  oldtypes(11) = MPI_DOUBLE_PRECISION   ! backward field
  oldtypes(12) = MPI_DOUBLE_PRECISION   ! convolved field

  ! find memory offsets, more stable then computing with MPI_TYPE_EXTEND
  call MPI_GET_ADDRESS(wt%ntotal_kernel,    offsets(1),  ierr)
  call MPI_GET_ADDRESS(wt%connectivity,     offsets(2),  ierr)
  call MPI_GET_ADDRESS(wt%vertices,         offsets(3),  ierr)
  call MPI_GET_ADDRESS(wt%kernel_values,    offsets(4),  ierr)
  call MPI_GET_ADDRESS(wt%kernel_variance,  offsets(5),  ierr)
  call MPI_GET_ADDRESS(wt%niterations,      offsets(6),  ierr)
  call MPI_GET_ADDRESS(wt%computation_time, offsets(7),  ierr)
  call MPI_GET_ADDRESS(wt%model,            offsets(8),  ierr)
  call MPI_GET_ADDRESS(wt%hetero_model,     offsets(9),  ierr)
  call MPI_GET_ADDRESS(wt%fw_field,         offsets(10), ierr)
  call MPI_GET_ADDRESS(wt%bw_field,         offsets(11), ierr)
  call MPI_GET_ADDRESS(wt%conv_field,       offsets(12), ierr)

  ! make offsets relative
  do i=2, size(offsets)
      offsets(i) = offsets(i) - offsets(1)
  end do
  offsets(1) = 0

  call MPI_TYPE_CREATE_STRUCT(nblocks, blocklengths, offsets, oldtypes, wt%mpitype, ierr)
  call MPI_TYPE_COMMIT(wt%mpitype, ierr)


end subroutine init_work_type
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Free the variables of a worktype
subroutine free_work_type()

  wt%ntotal_kernel             = -1
  wt%nelems_per_task           = -1
  wt%nvertices                 = -1
  wt%nvertices_per_elem        = -1
  wt%nbasisfuncs_per_elem      = -1
  wt%nmodel_parameters         = -1
  wt%nmodel_parameters_hetero  = -1

  wt%ndumps                    = -1
  wt%ndim                      = -1
  wt%dt                        = -1

  deallocate(wt%connectivity)
  deallocate(wt%vertices)
  deallocate(wt%kernel_values)
  deallocate(wt%kernel_variance)
  deallocate(wt%niterations)
  deallocate(wt%computation_time)
  deallocate(wt%model)
  deallocate(wt%hetero_model)
  deallocate(wt%fw_field)
  deallocate(wt%bw_field)
  deallocate(wt%conv_field)

end subroutine free_work_type
!-----------------------------------------------------------------------------------------
end module
!=========================================================================================
