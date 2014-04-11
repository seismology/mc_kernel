!=========================================================================================
module receivers_rdbm

  use global_parameters
  use source_class

  implicit none

  type receivers_rdbm_type
     type(src_param_type), allocatable    :: reci_sources(:)
     contains
     procedure, pass                      :: read_stations_file
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine read_stations_file(this)
  class(receivers_rdbm_type)    :: this
  integer                       :: ierror, lu_stations, num_rec
  integer                       :: i
  character(len=128)            :: junk
  character(len=20)             :: station_name
  character(len=20)             :: network_name
  real(kind=dp)                 :: reclat, reclon, recelevation, recbury

  if (verbose > 1) &
     write(6,*)'  reading receiver latitudes and longitudes from STATIONS...'

  ! count number of receivers first
  open(newunit=lu_stations, file='STATIONS', iostat=ierror, status='old', &
       action='read', position='rewind')

  num_rec = 0
  do while(ierror == 0)
     read(lu_stations, *, iostat=ierror) junk
     if(ierror == 0) num_rec = num_rec + 1
  enddo

  close(lu_stations)

  if (verbose > 0) &
     write(6,*)'  ...counted number of stations:', num_rec

  allocate(this%reci_sources(num_rec))

  open(newunit=lu_stations, file='STATIONS', iostat=ierror, status='old', &
       action='read', position='rewind')

  do i=1, num_rec
     read(lu_stations,*) station_name, network_name, reclat, reclon, recelevation, recbury
     if (verbose > 0) &
        write(6,*) station_name, network_name, reclat, reclon, recelevation, recbury
     call this%reci_sources(i)%init(reclat, reclon, (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0 /))
  enddo

  close(lu_stations)

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
