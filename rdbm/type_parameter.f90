!=========================================================================================
module type_parameter

  use global_parameters
  use source_class,       only : src_param_type
  use receiver_class,     only : rec_param_type
  use commpi,             only : pabort

  implicit none    

  type parameter_type
      type(src_param_type), allocatable     :: source
      type(rec_param_type), allocatable     :: receiver(:)
      integer                               :: nsrc, nrec

      character(len=512)                    :: sim_dir
      character(len=16)                     :: mode
      character(len=16)                     :: source_file_type
      character(len=512)                    :: source_file_wildcard
      character(len=512)                    :: source_file_name
      character(len=512), allocatable       :: source_files(:)
      integer                               :: nsources
      character(len=10)                     :: source_type
      character(len=8)                      :: receiver_file_type
      character(len=512)                    :: receiver_file
      character(len=1)                      :: component
      integer                               :: nsim_fwd
      logical                               :: resample
      integer                               :: nsamp
      logical                               :: apply_filter
      character(len=32)                     :: filterclass
      real(kind=dp)                         :: filterfreqs(4)
      logical                               :: time_shift
      integer                               :: buffer_size
      integer                               :: verbose
      logical                               :: parameters_read = .false.
      contains
         procedure, pass                    :: read_parameters
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine read_parameters(this, input_file_in)

  class(parameter_type)           :: this
  character(len=*), intent(in), optional :: input_file_in
  character(len=256)              :: input_file
  integer                         :: lu_inparam_basic, ioerr, narg, lu, i
  character(len=256)              :: line
  character(len=256)              :: keyword, keyvalue
  character(len=512)              :: fname_buf


  if (present(input_file_in)) then
      input_file = input_file_in
  else
      narg = command_argument_count()
      if (narg<1) then
          print *, 'ERROR: Argument "parameter file" is missing'
          call pabort
      end if
      call get_command_argument(1, input_file)
  end if

  write(lu_out,'(" Reading ", A, "...")') trim(input_file)
  open(newunit=lu_inparam_basic, file=trim(input_file), status='old', &
       action='read', iostat=ioerr)
  if (ioerr /= 0) then
     print *, 'ERROR: Check input file ''', trim(input_file), '''! Is it still there?' 
     call pabort
  end if

  do
     read(lu_inparam_basic, fmt='(a256)', iostat=ioerr) line
     if (ioerr < 0) exit
     if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
  
     read(line,*) keyword, keyvalue 
   
     parameter_to_read : select case(trim(keyword))

     case('SIM_DIR')
        this%sim_dir = keyvalue

     case('MODE')
        this%mode = keyvalue
     
     case('RECFILE_TYPE')
         this%receiver_file_type = keyvalue

     case('RECFILE_NAME')
        this%receiver_file = keyvalue

     case('SRCFILE_TYPE')
        this%source_file_type = keyvalue

     case('SRCFILE_NAME')
        this%source_file_name = keyvalue

     case('SRCFILE_WILDCARD')
        this%source_file_wildcard = keyvalue

     case('SRC_TYPE')
        this%source_type = keyvalue
     
     case('COMPONENT')
        read(keyvalue, *) this%component

     case('RESAMPLE')
        read(keyvalue, *) this%resample

     case('FILTER')
        read(keyvalue, *) this%apply_filter

     case('FILTER_DEF')
        read(keyvalue, *) this%filterclass, this%filterfreqs

     case('NSAMP')
        read(keyvalue, *) this%nsamp

     case('TIMESHIFT')
        read(keyvalue, *) this%time_shift

     case('NETCDF_BUFFER_SIZE')
        read(keyvalue, *) this%buffer_size
     
     case('VERBOSITY')
        read(keyvalue, *) this%verbose

     case default
        write(6,*) 'ERROR: unknown keyword in inparam_basic "', trim(keyword), '"' 
        call pabort

     end select parameter_to_read
     print "('  ', A32,' ',A)", keyword, trim(keyvalue)
  end do
  close(lu_inparam_basic)

  verbose = this%verbose

  if (trim(this%source_file_type) == 'cmtsolutions') then

     ! write filenames of CMTSOLUTION files to cmtsolutions.tmp
     call system('ls '//this%source_file_wildcard//' > cmtsolutions.tmp')

     ! count number of sources
     this%nsources = 0
     open(newunit=lu, file='cmtsolutions.tmp', action="read")
     do
        read(lu, fmt='(a)', iostat=ioerr) fname_buf
        if (ioerr/=0) exit
        this%nsources = this%nsources + 1
     end do

     if (this%verbose > 0) then
        write(6,*) '  '
        write(6,*) "nsources: " , this%nsources
     endif

     ! read filenames
     rewind(lu)
     allocate(this%source_files(this%nsources))
     if (this%verbose > 0) write(6,*) 'source files:'
     do i=1, this%nsources
        read(lu, fmt='(a)', iostat=ioerr) this%source_files(i)
        if (this%verbose > 0) write(6,*) '  ', trim(this%source_files(i))
     end do
     if (this%verbose > 0) write(6,*) '  '

     close(lu)
  else if (trim(this%source_file_type) == 'cmtsolution') then
     this%nsources = 1
  endif

  this%parameters_read = .true.

end subroutine read_parameters
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
