!=========================================================================================
module receivers_rdbm

  use global_parameters
  use source_class
  use commpi,               only : pabort

  implicit none

  private
  public   :: receivers_rdbm_type

  type receivers_rdbm_type
     type(src_param_type), allocatable      :: reci_sources(:)
     integer                                :: num_rec
     logical                                :: initialized = .false.
     contains
     procedure, pass                        :: init_xyz
     procedure, pass                        :: read_stations_file
     procedure, pass                        :: read_receiver_dat
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine init_xyz(this, xyz)
  class(receivers_rdbm_type)    :: this
  real(kind=dp), intent(in)     :: xyz(:,:)

  integer                       :: i

  if (this%initialized) then
     write(6,*) 'ERROR: trying to initialize aleady initalized rdbm receiver object'
     call pabort
  endif

  this%num_rec = size(xyz, dim=2)

  do i=1, this%num_rec
     call this%reci_sources(i)%init_xyz(xyz(1,i), xyz(2,i), xyz(3,i), &
                                        (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0 /))
  enddo

  this%initialized = .true.

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_stations_file(this)
  class(receivers_rdbm_type)    :: this
  integer                       :: ierror, lu_stations
  integer                       :: i
  character(len=128)            :: junk
  character(len=4)              :: station_name
  character(len=2)              :: network_name
  real(kind=dp)                 :: reclat, reclon, recelevation, recbury

  if (this%initialized) then
     write(6,*) 'ERROR: trying to initialize aleady initalized rdbm receiver object'
     call pabort
  endif

  if (verbose > 1) &
     write(6,*) '  reading receiver latitudes and longitudes from STATIONS...'

  ! count number of receivers first
  open(newunit=lu_stations, file='STATIONS', iostat=ierror, status='old', &
       action='read', position='rewind')

  this%num_rec = 0
  do while(ierror == 0)
     read(lu_stations, *, iostat=ierror) junk
     if(ierror == 0) this%num_rec = this%num_rec + 1
  enddo

  close(lu_stations)

  if (verbose > 0) &
     write(6,*)'  ...counted number of stations:', this%num_rec

  allocate(this%reci_sources(this%num_rec))

  open(newunit=lu_stations, file='STATIONS', iostat=ierror, status='old', &
       action='read', position='rewind')

  do i=1, this%num_rec
     read(lu_stations,*) station_name, network_name, reclat, reclon, recelevation, recbury
     if (verbose > 0) &
        write(6,*) station_name, network_name, reclat, reclon, recelevation, recbury
     if (recelevation .ne. 0 .or. recbury .ne. 0) &
        write(6,*) 'WARNING: receivers are assumed to be at the earth surface, '//&
                   'but elevation or bury is nonzero in STATIONS file for station '//&
                   station_name
     ! TODO: hardcoded receiver depth to 0 for now
     call this%reci_sources(i)%init(reclat, reclon, &
                                    (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0 /), 0d0)
  enddo

  close(lu_stations)

  this%initialized = .true.

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_receiver_dat(this)
  class(receivers_rdbm_type)    :: this
  integer                       :: ierror, lu_stations
  integer                       :: i
  real(kind=dp)                 :: reccolat, reclon

  if (this%initialized) then
     write(6,*) 'ERROR: trying to initialize aleady initalized rdbm receiver object'
     call pabort
  endif

  if (verbose > 1) &
     write(6,*) '  reading receiver latitudes and longitudes from receivers.dat...'

  ! count number of receivers first
  open(newunit=lu_stations, file='receivers.dat', iostat=ierror, status='old', &
       action='read', position='rewind')

  read(lu_stations, *, iostat=ierror) this%num_rec

  allocate(this%reci_sources(this%num_rec))

  do i=1, this%num_rec
     read(lu_stations,*) reccolat, reclon
     ! TODO: hardcoded receiver depth to 0 for now
     call this%reci_sources(i)%init(90 - reccolat, reclon, (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0 /), 0d0)
  enddo

  close(lu_stations)

  this%initialized = .true.

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
