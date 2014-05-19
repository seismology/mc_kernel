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
      character(len=512)                    :: source_file
      character(len=10)                     :: source_type
      character(len=512)                    :: receiver_file
      character(len=8)                      :: receiver_file_type
      character(len=1)                      :: component
      integer                               :: nsim_fwd
      integer                               :: nsamp
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
  integer                         :: lu_inparam_basic, ioerr, narg
  character(len=256)              :: line
  character(len=256)              :: keyword, keyvalue


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

     case('NSAMP')
        read(keyvalue, *) this%nsamp

     case('SRCFILE_TYPE')
        this%source_file = keyvalue

     case('SRC_TYPE')
        this%source_type = keyvalue
     
     case('RECFILE_TYPE')
         this%receiver_file_type = keyvalue

     case('NETCDF_BUFFER_SIZE')
        read(keyvalue, *) this%buffer_size
     
     case('COMPONENT')
        read(keyvalue, *) this%component

     case('VERBOSITY')
        read(keyvalue, *) this%verbose

     end select parameter_to_read
     print "('  ', A32,' ',A)", keyword, trim(keyvalue)
  end do
  close(lu_inparam_basic)

  verbose = this%verbose

  this%parameters_read = .true.

end subroutine read_parameters
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
