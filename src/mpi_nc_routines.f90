!=========================================================================================
module mpi_nc_routines

# ifndef include_mpi
  use mpi,               only   : MPI_COMM_RANK
# endif
  use commpi,            only   : pbroadcast
  use global_parameters, only   : sp, dp
  use nc_routines,       only   : nc_getvar_by_name

  implicit none
# ifdef include_mpi
  include 'mpif.h'
# endif

  public                       :: mpi_getvar_by_name

  interface mpi_getvar_by_name
    module procedure           :: mpi_getvar_by_name_sp_1d
    module procedure           :: mpi_getvar_by_name_sp_2d
    module procedure           :: mpi_getvar_by_name_int_1d
    module procedure           :: mpi_getvar_by_name_int_2d
    module procedure           :: mpi_getvar_by_name_int_3d
  end interface mpi_getvar_by_name
  

contains
!-----------------------------------------------------------------------------------------
subroutine mpi_getvar_by_name_sp_1d(ncid, varname, limits, values, comm)

  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  real(kind=sp), allocatable, intent(inout)  :: values(:)
  real(kind=sp), intent(in), optional        :: limits(2)
  integer, intent(in)                        :: comm

  integer                                    :: nelem, myrank, ierror

  call MPI_COMM_RANK(comm, myrank, ierror )

  if (myrank==0) then
    call nc_getvar_by_name(ncid    = ncid,     &
                           varname = varname,  &
                           limits  = limits,   & 
                           values  = values )
    nelem = size(values)
    call pbroadcast(nelem, 0, comm)
  else
    call pbroadcast(nelem, 0, comm)
    allocate(values(nelem))
  end if

  call pbroadcast(values, 0, comm)

end subroutine mpi_getvar_by_name_sp_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine mpi_getvar_by_name_sp_2d(ncid, varname, limits, values, comm)

  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  real(kind=sp), allocatable, intent(inout)  :: values(:,:)
  real(kind=sp), allocatable                 :: values_tmp(:)  !< 1D variable for broadcasting
  real(kind=sp), intent(in), optional        :: limits(2)
  integer, intent(in)                        :: comm

  integer                                    :: nelem_1, nelem_2, myrank, ierror

  call MPI_COMM_RANK(comm, myrank, ierror )

  if (myrank==0) then
    call nc_getvar_by_name(ncid    = ncid,     &
                           varname = varname,  &
                           limits  = limits,   & 
                           values  = values )
    nelem_1 = size(values, 1)
    nelem_2 = size(values, 2)
    call pbroadcast(nelem_1, 0, comm)
    call pbroadcast(nelem_2, 0, comm)
  else
    call pbroadcast(nelem_1, 0, comm)
    call pbroadcast(nelem_2, 0, comm)
    allocate(values(nelem_1, nelem_2))
  end if

  ! Create 1D array for broadcasting
  allocate(values_tmp(nelem_1*nelem_2))
  if (myrank==0) values_tmp = reshape(values, shape=[nelem_1*nelem_2])

  call pbroadcast(values_tmp, 0, comm)

  values = reshape(values_tmp, shape=[nelem_1, nelem_2])

end subroutine mpi_getvar_by_name_sp_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine mpi_getvar_by_name_int_1d(ncid, varname, limits, values, comm)

  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  integer, allocatable, intent(inout)        :: values(:)
  integer, intent(in), optional              :: limits(2)
  integer, intent(in)                        :: comm

  integer                                    :: nelem, myrank, ierror

  call MPI_COMM_RANK(comm, myrank, ierror )

  if (myrank==0) then
    call nc_getvar_by_name(ncid    = ncid,     &
                           varname = varname,  &
                           limits  = limits,   & 
                           values  = values )
    nelem = size(values)
    call pbroadcast(nelem, 0, comm)
  else
    call pbroadcast(nelem, 0, comm)
    allocate(values(nelem))
  end if

  call pbroadcast(values, 0, comm)

end subroutine mpi_getvar_by_name_int_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine mpi_getvar_by_name_int_2d(ncid, varname, limits, values, comm)

  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  integer, allocatable, intent(inout)        :: values(:,:)
  integer, allocatable                       :: values_tmp(:)  !< 1D variable for broadcasting
  integer, intent(in), optional              :: limits(2)
  integer, intent(in)                        :: comm

  integer                                    :: nelem_1, nelem_2, myrank, ierror

  call MPI_COMM_RANK(comm, myrank, ierror )

  if (myrank==0) then
    call nc_getvar_by_name(ncid    = ncid,     &
                           varname = varname,  &
                           limits  = limits,   & 
                           values  = values )
    nelem_1 = size(values, 1)
    nelem_2 = size(values, 2)
    call pbroadcast(nelem_1, 0, comm)
    call pbroadcast(nelem_2, 0, comm)
  else
    call pbroadcast(nelem_1, 0, comm)
    call pbroadcast(nelem_2, 0, comm)
    allocate(values(nelem_1, nelem_2))
  end if

  ! Create 1D array for broadcasting
  allocate(values_tmp(nelem_1*nelem_2))
  if (myrank==0) values_tmp = reshape(values, shape=[nelem_1*nelem_2])

  call pbroadcast(values_tmp, 0, comm)

  values = reshape(values_tmp, shape=[nelem_1, nelem_2])

end subroutine mpi_getvar_by_name_int_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine mpi_getvar_by_name_int_3d(ncid, varname, limits, values, comm)

  integer, intent(in)                        :: ncid
  character(len=*), intent(in)               :: varname
  integer, allocatable, intent(inout)        :: values(:,:,:)
  integer, allocatable                       :: values_tmp(:)  !< 1D variable for broadcasting
  integer, intent(in), optional              :: limits(2)
  integer, intent(in)                        :: comm

  integer                                    :: nelem_1, nelem_2, nelem_3
  integer                                    :: myrank, ierror

  call MPI_COMM_RANK(comm, myrank, ierror )

  if (myrank==0) then
    call nc_getvar_by_name(ncid    = ncid,     &
                           varname = varname,  &
                           limits  = limits,   & 
                           values  = values )
    nelem_1 = size(values, 1)
    nelem_2 = size(values, 2)
    nelem_3 = size(values, 3)
    call pbroadcast(nelem_1, 0, comm)
    call pbroadcast(nelem_2, 0, comm)
    call pbroadcast(nelem_3, 0, comm)
  else
    call pbroadcast(nelem_1, 0, comm)
    call pbroadcast(nelem_2, 0, comm)
    call pbroadcast(nelem_3, 0, comm)
    allocate(values(nelem_1, nelem_2, nelem_3))
  end if

  ! Create 1D array for broadcasting
  allocate(values_tmp(nelem_1*nelem_2*nelem_3))
  if (myrank==0) values_tmp = reshape(values, shape=[nelem_1*nelem_2*nelem_3])

  call pbroadcast(values_tmp, 0, comm)

  values = reshape(values_tmp, shape=[nelem_1, nelem_2, nelem_3])

end subroutine mpi_getvar_by_name_int_3d
!-----------------------------------------------------------------------------------------

end module mpi_nc_routines
!=========================================================================================
