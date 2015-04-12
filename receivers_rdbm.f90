!=========================================================================================
module receivers_rdbm

  use global_parameters
  use source_class,         only : src_param_type
  use receiver_class,      only : rec_param_type
  use commpi,               only : pabort

  implicit none

  private
  public   :: receivers_rdbm_type

  type receivers_rdbm_type
     type(src_param_type), allocatable      :: reci_sources(:)
     integer                                :: num_rec
     logical                                :: initialized = .false.
     contains
     procedure, pass                        :: create_reci_sources
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine create_reci_sources(this, rec_in)
  class(receivers_rdbm_type)       :: this
  type(rec_param_type), intent(in) :: rec_in(:)

  integer                          :: ierror
  character(len=512)               :: receivers_file
  integer                          :: i,nrec
  real(kind=dp)                    :: reccolat, reclon

  nrec = size(rec_in)

  if (this%initialized) then
     write(6,*) 'ERROR: trying to initialize aleady initalized rdbm receiver object'
     call pabort
  endif

  allocate(this%reci_sources(nrec))

  do i=1,nrec
     call this%reci_sources(i)%init(rec_in(i)%latd, rec_in(i)%lond, (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0 /), 0d0)
  enddo


  this%initialized = .true.

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
