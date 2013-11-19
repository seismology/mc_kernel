module buffer
   use global_parameters, only     : sp, dp
   implicit none
   type buffer_type
      private
      integer, allocatable        :: idx(:)
      real(kind=sp), allocatable  :: val(:,:)
      integer                     :: nbuffer, nvalues
      logical                     :: initialized = .false.
      integer                     :: naccess, nhit

      contains
         procedure, pass          :: init
         procedure, pass          :: freeme
         procedure, pass          :: get
         procedure, pass          :: put
         procedure, pass          :: efficiency
   end type

contains

function init(this, nbuffer, nvalues)
    class(buffer_type)      :: this
    integer, intent(in)     :: nbuffer  !< How many elements should the buffer be 
                                        !! able to store? 
    integer, intent(in)     :: nvalues  !< How many values should one buffer store
                                        !! i.e. how many samples of a time trace.
    integer                 :: init     !< Return value, 0=Success

    init = -1
    allocate(this%val(nvalues, nbuffer))
    allocate(this%idx(nbuffer))

    this%nvalues = nvalues
    this%nbuffer = nbuffer

    this%naccess = 0
    this%nhit    = 0

    this%idx     = -1
    this%val     = 0.0

    this%initialized = .true.
    init = 0
end function

function freeme(this)
    class(buffer_type)       :: this
    integer                  :: freeme    !< Return value, 0=Success


    freeme = -1
    deallocate(this%val)
    deallocate(this%idx)

    this%nvalues = 0
    this%nbuffer = 0 

    this%initialized = .false.
    freeme = 0
end function

function get(this, iindex, values)
    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues)
    integer                     :: get      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       stop "Buffer has not been initialized"
    end if
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       stop
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
    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues) !< Values to store
    integer                   :: put    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       stop
    end if
    if (any(this%idx==iindex)) then
       print *, 'Buffer with this index', iindex, ' already exists'
       put = -1
    else
       call random_number(randtemp)
       ibuffer = int(randtemp*this%nbuffer)+1
       this%idx(ibuffer) = iindex
       this%val(:,ibuffer) = values
       put = 0
    end if

 end function

 function efficiency(this)
    class(buffer_type)  :: this
    real(kind=sp)       :: efficiency
  
    efficiency = real(this%nhit)/real(this%naccess)
 end function


end module buffer
