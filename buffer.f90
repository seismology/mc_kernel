module buffer
   use global_parameters, only:                 sp
   implicit none
   type buffer_type
      private
      integer, allocatable               :: idx(:)
      real(kind=sp), allocatable         :: val(:,:)
      integer                            :: nbuffer, nvalues
      logical                            :: initialized = .false.
      integer                            :: naccess, nhit

      contains
         procedure, pass                 :: init
         procedure, pass                 :: freeme
         procedure, pass                 :: get
         procedure, pass                 :: put
         procedure, pass                 :: efficiency
   end type

contains

function init(this, nbuffer, nvalues)
    class(buffer_type)      :: this
    integer, intent(in)     :: nbuffer, nvalues
    integer                 :: init

    init = -1
    allocate(this%val(nvalues, nbuffer))
    allocate(this%idx(nbuffer))

    this%nvalues = nvalues
    this%nbuffer = nbuffer

    this%naccess = 0
    this%nhit    = 0

    this%initialized = .true.
    init = 0
end function

function freeme(this)
    class(buffer_type)       :: this
    integer                  :: freeme

    freeme = -1
    deallocate(this%val)
    deallocate(this%idx)

    this%nvalues = 0
    this%nbuffer = 0 

    this%initialized = .false.
    freeme = 0
end function

function get(this, iindex, values)
    class(buffer_type)       :: this
    integer, intent(in)     :: iindex
    real, intent(out)       :: values(this%nvalues)
    integer                 :: get !< status value, 0=success
    integer                 :: ibuffer

    if (.not.this%initialized) then
       stop "Buffer has not been initialized"
    end if
    this%naccess = this%naccess + 1
    get = -1
    do ibuffer = 1, this%nbuffer
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val(:,ibuffer) 
       this%nhit = this%nhit + 1
       get = 0 
       exit
    end do

 end function

 function put(this, iindex, values)
    class(buffer_type)   :: this
    integer, intent(in)  :: iindex
    real, intent(in)     :: values(this%nvalues)
    integer              :: put
    real                 :: randtemp
    integer              :: ibuffer

    put = -1
    call random_number(randtemp)
    ibuffer = int(randtemp*this%nbuffer)+1

    this%idx(ibuffer) = iindex
    this%val(:,ibuffer) = values

    put = 0

 end function

 function efficiency(this)
    class(buffer_type)  :: this
    real(kind=sp)       :: efficiency
  
    efficiency = real(this%nhit)/real(this%naccess)
 end function


end module buffer
