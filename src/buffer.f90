!=========================================================================================
module buffer

   use global_parameters, only     : sp, dp, lu_out, verbose
   use commpi,            only     : pabort

   implicit none
   private
   public :: buffer_type

   type buffer_type
      private
      integer, allocatable        :: idx(:)
      real(kind=sp), allocatable  :: val_1d_sp(:,:)
      real(kind=sp), allocatable  :: val_2d_sp(:,:,:)
      real(kind=sp), allocatable  :: val_3d_sp(:,:,:,:)
      real(kind=sp), allocatable  :: val_4d_sp(:,:,:,:,:)
      real(kind=sp), allocatable  :: val_5d_sp(:,:,:,:,:,:)
      real(kind=dp), allocatable  :: val_1d_dp(:,:)
      real(kind=dp), allocatable  :: val_2d_dp(:,:,:)
      real(kind=dp), allocatable  :: val_3d_dp(:,:,:,:)
      real(kind=dp), allocatable  :: val_4d_dp(:,:,:,:,:)
      real(kind=dp), allocatable  :: val_5d_dp(:,:,:,:,:,:)
      integer                     :: nbuffer, nput
      integer                     :: ndim,  nvalues_dim1, nvalues_dim2, &
                                     nvalues_dim3, nvalues_dim4, nvalues_dim5
      logical                     :: initialized = .false.
      integer                     :: naccess, nhit

      contains
         procedure, pass          :: init
         procedure, pass          :: freeme
         procedure, pass          :: efficiency

         procedure, pass, private :: get_1d_sp
         procedure, pass, private :: get_2d_sp
         procedure, pass, private :: get_3d_sp
         procedure, pass, private :: get_4d_sp
         procedure, pass, private :: get_5d_sp
         procedure, pass, private :: get_1d_dp
         procedure, pass, private :: get_2d_dp
         procedure, pass, private :: get_3d_dp
         procedure, pass, private :: get_4d_dp
         procedure, pass, private :: get_5d_dp
         generic, public          :: get => get_1d_sp, get_2d_sp, get_3d_sp,  &
                                            get_4d_sp, get_5d_sp,             &
                                            get_1d_dp, get_2d_dp, get_3d_dp,  &
                                            get_4d_dp, get_5d_dp

         procedure, pass, private :: put_1d_sp
         procedure, pass, private :: put_2d_sp
         procedure, pass, private :: put_3d_sp
         procedure, pass, private :: put_4d_sp
         procedure, pass, private :: put_5d_sp
         procedure, pass, private :: put_1d_dp
         procedure, pass, private :: put_2d_dp
         procedure, pass, private :: put_3d_dp
         procedure, pass, private :: put_4d_dp
         procedure, pass, private :: put_5d_dp
         generic, public          :: put => put_1d_sp, put_2d_sp, put_3d_sp,  &
                                            put_4d_sp, put_5d_sp,             &
                                            put_1d_dp, put_2d_dp, put_3d_dp,  &
                                            put_4d_dp, put_5d_dp
   end type
    

contains

!-----------------------------------------------------------------------------------------
function init(this, nbuffer, value_type, nvalues_dim1, nvalues_dim2, nvalues_dim3, &
              nvalues_dim4, nvalues_dim5)
    use simple_routines, only: to_lower
    class(buffer_type)      :: this
    integer, intent(in)     :: nbuffer  !< How many elements should the buffer be 
                                        !! able to store? 
    character(len=4)        :: value_type
    integer, intent(in)     :: nvalues_dim1  !< How many values should one buffer store
                                        !! i.e. how many samples of a time trace.
    integer, intent(in), optional :: nvalues_dim2, nvalues_dim3, nvalues_dim4, nvalues_dim5
    integer                 :: init     !< Return value, 0=Success
    integer                 :: buffer_size

    init = -1

    if (.not. present(nvalues_dim2) .and. .not. present(nvalues_dim3) &
        .and. .not. present(nvalues_dim4) .and. .not. present(nvalues_dim5)) then
      
      select case(to_lower(value_type))
      case('real')
        allocate(this%val_1d_sp(nvalues_dim1, nbuffer))
        this%val_1d_sp = 0
      case('dble')
        allocate(this%val_1d_dp(nvalues_dim1, nbuffer))
        this%val_1d_dp = 0
      end select

      this%ndim = 1
      this%nvalues_dim1 = nvalues_dim1
      this%nvalues_dim2 = 0
      this%nvalues_dim3 = 0
      this%nvalues_dim4 = 0
      buffer_size = nbuffer * nvalues_dim1

    elseif (present(nvalues_dim2) .and. .not. present(nvalues_dim3) &
            .and. .not. present(nvalues_dim4) .and. .not. present(nvalues_dim5)) then
      select case(to_lower(value_type))
      case('real')
        allocate(this%val_2d_sp(nvalues_dim1, nvalues_dim2, nbuffer))
        this%val_2d_sp = 0
      case('dble')
        allocate(this%val_2d_dp(nvalues_dim1, nvalues_dim2, nbuffer))
        this%val_2d_dp = 0
      end select

      this%ndim = 2
      this%nvalues_dim1 = nvalues_dim1
      this%nvalues_dim2 = nvalues_dim2
      this%nvalues_dim3 = 0
      this%nvalues_dim4 = 0
      this%nvalues_dim5 = 0
      buffer_size = nbuffer * nvalues_dim1 * nvalues_dim2

    elseif (present(nvalues_dim2) .and. present(nvalues_dim3) &
            .and. .not. present(nvalues_dim4) .and. .not. present(nvalues_dim5)) then
      select case(to_lower(value_type))
      case('real')
        allocate(this%val_3d_sp(nvalues_dim1, nvalues_dim2, nvalues_dim3, nbuffer))
        this%val_3d_sp = 0
      case('dble')
        allocate(this%val_3d_dp(nvalues_dim1, nvalues_dim2, nvalues_dim3, nbuffer))
        this%val_3d_dp = 0
      end select

      this%ndim = 3
      this%nvalues_dim1 = nvalues_dim1
      this%nvalues_dim2 = nvalues_dim2
      this%nvalues_dim3 = nvalues_dim3
      this%nvalues_dim4 = 0
      this%nvalues_dim5 = 0
      buffer_size = nbuffer * nvalues_dim1 * nvalues_dim2 * nvalues_dim3

    elseif (present(nvalues_dim2) .and. present(nvalues_dim3) &
            .and. present(nvalues_dim4) .and. .not. present(nvalues_dim5)) then
      select case(to_lower(value_type))
      case('real')
        allocate(this%val_4d_sp(nvalues_dim1, nvalues_dim2, nvalues_dim3, nvalues_dim4, nbuffer))
        this%val_4d_sp = 0
      case('dble')
        allocate(this%val_4d_dp(nvalues_dim1, nvalues_dim2, nvalues_dim3, nvalues_dim4, nbuffer))
        this%val_4d_dp = 0
      end select

      this%ndim = 4
      this%nvalues_dim1 = nvalues_dim1
      this%nvalues_dim2 = nvalues_dim2
      this%nvalues_dim3 = nvalues_dim3
      this%nvalues_dim4 = nvalues_dim4
      this%nvalues_dim5 = 0
      buffer_size = nbuffer * nvalues_dim1 * nvalues_dim2 * nvalues_dim3 * nvalues_dim4

    elseif (present(nvalues_dim2) .and. present(nvalues_dim3) &
            .and. present(nvalues_dim4) .and. present(nvalues_dim5)) then
      select case(to_lower(value_type))
      case('real')
        allocate(this%val_5d_sp(nvalues_dim1, nvalues_dim2, nvalues_dim3, nvalues_dim4, &
                                nvalues_dim5, nbuffer))
        this%val_5d_sp = 0
      case('dble')
        allocate(this%val_5d_dp(nvalues_dim1, nvalues_dim2, nvalues_dim3, nvalues_dim4, &
                                nvalues_dim5, nbuffer))
        this%val_5d_dp = 0
      end select
        this%ndim = 5
        this%nvalues_dim1 = nvalues_dim1
        this%nvalues_dim2 = nvalues_dim2
        this%nvalues_dim3 = nvalues_dim3
        this%nvalues_dim4 = nvalues_dim4
        this%nvalues_dim5 = nvalues_dim5
        buffer_size = nbuffer * nvalues_dim1 * nvalues_dim2 * nvalues_dim3 * nvalues_dim4 &
                      * nvalues_dim5

    else
        write(6,*) 'ERROR: illegal combination of optional arguments'
        call pabort
    endif

    if (verbose > 0) then
       write(lu_out, '(A,I5,A,I5,I5,I5,I5,A)') ' Initialize buffer with ', nbuffer, &
                    ' memories for (', this%nvalues_dim1, this%nvalues_dim2, &
                    this%nvalues_dim3, this%nvalues_dim4, ') values'
       write(lu_out,'(A,I12,A,F8.1,A)') 'this is a total of ', buffer_size, ' numbers and uses a size of', &
                                        buffer_size * 4d0 / 1024d0**2, ' MB in memory'
    endif
    
    allocate(this%idx(nbuffer))

    this%nbuffer = nbuffer

    this%naccess = 0
    this%nhit    = 0
    this%nput    = 0

    this%idx     = -1

    this%initialized = .true.
    init = 0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function freeme(this)

    class(buffer_type)       :: this
    integer                  :: freeme    !< Return value, 0=Success

    freeme = -1
    if (allocated(this%val_1d_sp)) deallocate(this%val_1d_sp)
    if (allocated(this%val_2d_sp)) deallocate(this%val_2d_sp)
    if (allocated(this%val_3d_sp)) deallocate(this%val_3d_sp)
    if (allocated(this%val_4d_sp)) deallocate(this%val_4d_sp)
    if (allocated(this%val_5d_sp)) deallocate(this%val_5d_sp)
    if (allocated(this%val_1d_dp)) deallocate(this%val_1d_dp)
    if (allocated(this%val_2d_dp)) deallocate(this%val_2d_dp)
    if (allocated(this%val_3d_dp)) deallocate(this%val_3d_dp)
    if (allocated(this%val_4d_dp)) deallocate(this%val_4d_dp)
    if (allocated(this%val_5d_dp)) deallocate(this%val_5d_dp)
    deallocate(this%idx)

    this%nvalues_dim1 = 0
    this%nvalues_dim2 = 0
    this%nvalues_dim3 = 0
    this%nvalues_dim4 = 0
    this%nvalues_dim5 = 0
    this%nbuffer = 0 

    this%initialized = .false.
    freeme = 0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_1d_sp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues_dim1)
    integer                     :: get_1d_sp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if

    if (this%ndim /= 1) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_1d_sp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle
       values = this%val_1d_sp(:,ibuffer) 
       this%nhit = this%nhit + 1
       get_1d_sp = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_2d_sp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2)
    integer                     :: get_2d_sp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if
    
    if (this%ndim /= 2) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_2d_sp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_2d_sp(:,:,ibuffer) 
       this%nhit = this%nhit + 1
       get_2d_sp = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_3d_sp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                          this%nvalues_dim3)
    integer                     :: get_3d_sp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if
    
    if (this%ndim /= 3) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_3d_sp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_3d_sp(:,:,:,ibuffer) 
       this%nhit = this%nhit + 1
       get_3d_sp = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_4d_sp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                          this%nvalues_dim3, this%nvalues_dim4)
    integer                     :: get_4d_sp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if
    
    if (this%ndim /= 4) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_4d_sp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_4d_sp(:,:,:,:,ibuffer)
       this%nhit = this%nhit + 1
       get_4d_sp = 0 
       exit
    end do

end function get_4d_sp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_5d_sp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                          this%nvalues_dim3, this%nvalues_dim4, &
                                          this%nvalues_dim5)
    integer                     :: get_5d_sp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if
    
    if (this%ndim /= 5) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_5d_sp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_5d_sp(:,:,:,:,:,ibuffer)
       this%nhit = this%nhit + 1
       get_5d_sp = 0 
       exit
    end do

end function get_5d_sp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_1d_sp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues_dim1) !< Values to store
    integer                   :: put_1d_sp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 1) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_1d_sp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_1d_sp(:,ibuffer) = values
       this%nput = this%nput + 1
       put_1d_sp = 0
    end if

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_2d_sp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2) !< Values to store
    integer                   :: put_2d_sp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 2) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_2d_sp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_2d_sp(:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_2d_sp = 0
    end if

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_3d_sp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                        this%nvalues_dim3) !< Values to store
    integer                   :: put_3d_sp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 3) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_3d_sp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_3d_sp(:,:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_3d_sp = 0
    end if

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_4d_sp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                        this%nvalues_dim3, this%nvalues_dim4) !< Values to store
    integer                   :: put_4d_sp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 4) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_4d_sp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_4d_sp(:,:,:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_4d_sp = 0
    end if

end function put_4d_sp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_5d_sp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                        this%nvalues_dim3, this%nvalues_dim4, &
                                        this%nvalues_dim5) !< Values to store
    integer                   :: put_5d_sp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 5) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_5d_sp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_5d_sp(:,:,:,:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_5d_sp = 0
    end if

end function put_5d_sp
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Double precision routines
!-----------------------------------------------------------------------------------------
function get_1d_dp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=dp), intent(out)  :: values(this%nvalues_dim1)
    integer                     :: get_1d_dp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if

    if (this%ndim /= 1) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_1d_dp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle
       values = this%val_1d_dp(:,ibuffer) 
       this%nhit = this%nhit + 1
       get_1d_dp = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_2d_dp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=dp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2)
    integer                     :: get_2d_dp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if
    
    if (this%ndim /= 2) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_2d_dp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_2d_dp(:,:,ibuffer) 
       this%nhit = this%nhit + 1
       get_2d_dp = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_3d_dp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=dp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                          this%nvalues_dim3)
    integer                     :: get_3d_dp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if
    
    if (this%ndim /= 3) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_3d_dp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_3d_dp(:,:,:,ibuffer) 
       this%nhit = this%nhit + 1
       get_3d_dp = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_4d_dp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=dp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                          this%nvalues_dim3, this%nvalues_dim4)
    integer                     :: get_4d_dp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if
    
    if (this%ndim /= 4) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_4d_dp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_4d_dp(:,:,:,:,ibuffer)
       this%nhit = this%nhit + 1
       get_4d_dp = 0 
       exit
    end do

end function get_4d_dp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_5d_dp(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=dp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                          this%nvalues_dim3, this%nvalues_dim4, &
                                          this%nvalues_dim5)
    integer                     :: get_5d_dp      !< status value, 0=success
    integer                     :: ibuffer

    if (.not.this%initialized) then
       write(*, '(A)') "ERROR: Buffer has not been initialized"
       call pabort 
    end if
    
    if (this%ndim /= 5) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if
    
    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if
    
    this%naccess = this%naccess + 1
    get_5d_dp = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_5d_dp(:,:,:,:,:,ibuffer)
       this%nhit = this%nhit + 1
       get_5d_dp = 0 
       exit
    end do

end function get_5d_dp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_1d_dp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=dp), intent(in) :: values(this%nvalues_dim1) !< Values to store
    integer                   :: put_1d_dp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 1) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_1d_dp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_1d_dp(:,ibuffer) = values
       this%nput = this%nput + 1
       put_1d_dp = 0
    end if

end function put_1d_dp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_2d_dp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=dp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2) !< Values to store
    integer                   :: put_2d_dp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 2) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_2d_dp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_2d_dp(:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_2d_dp = 0
    end if

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_3d_dp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=dp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                        this%nvalues_dim3) !< Values to store
    integer                   :: put_3d_dp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 3) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_3d_dp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_3d_dp(:,:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_3d_dp = 0
    end if

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_4d_dp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=dp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                        this%nvalues_dim3, this%nvalues_dim4) !< Values to store
    integer                   :: put_4d_dp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 4) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_4d_dp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_4d_dp(:,:,:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_4d_dp = 0
    end if

end function put_4d_dp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_5d_dp(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=dp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                        this%nvalues_dim3, this%nvalues_dim4, &
                                        this%nvalues_dim5) !< Values to store
    integer                   :: put_5d_dp    !< Return value, 0=Success
    real(kind=dp)             :: randtemp
    integer                   :: ibuffer

    if (iindex<0) then
       write(*,*) 'ERROR: Buffer index must be larger zero, is: ', iindex
       call pabort
    end if

    if (this%ndim /= 5) then
       write(*,*) 'ERROR: wrong rank of argument "values", buffer was initialized with ndim = ', this%ndim
       call pabort
    end if

    if (any(this%idx==iindex)) then
       put_5d_dp = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_5d_dp(:,:,:,:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_5d_dp = 0
    end if

end function put_5d_dp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function efficiency(this)
    class(buffer_type)  :: this
    real(kind=sp)       :: efficiency
  
    if (this%naccess > 0) then
       efficiency = real(this%nhit) / real(this%naccess)
    else
       efficiency = real(this%nhit)
    endif

end function
!-----------------------------------------------------------------------------------------

end module buffer
!=========================================================================================
