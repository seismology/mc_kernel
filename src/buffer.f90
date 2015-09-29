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
      real(kind=sp), allocatable  :: val_1d(:,:)
      real(kind=sp), allocatable  :: val_2d(:,:,:)
      real(kind=sp), allocatable  :: val_3d(:,:,:,:)
      real(kind=sp), allocatable  :: val_4d(:,:,:,:,:)
      integer                     :: nbuffer, nput
      integer                     :: ndim,  nvalues_dim1, nvalues_dim2, &
                                     nvalues_dim3, nvalues_dim4
      logical                     :: initialized = .false.
      integer                     :: naccess, nhit

      contains
         procedure, pass          :: init
         procedure, pass          :: freeme
         procedure, pass          :: efficiency

         procedure, pass, private :: get_1d
         procedure, pass, private :: get_2d
         procedure, pass, private :: get_3d
         procedure, pass, private :: get_4d
         generic, public          :: get => get_1d, get_2d, get_3d, get_4d

         procedure, pass, private :: put_1d
         procedure, pass, private :: put_2d
         procedure, pass, private :: put_3d
         procedure, pass, private :: put_4d
         generic, public          :: put => put_1d, put_2d, put_3d, put_4d
   end type
    

contains

!-----------------------------------------------------------------------------------------
function init(this, nbuffer, nvalues_dim1, nvalues_dim2, nvalues_dim3, nvalues_dim4)

    class(buffer_type)      :: this
    integer, intent(in)     :: nbuffer  !< How many elements should the buffer be 
                                        !! able to store? 
    integer, intent(in)     :: nvalues_dim1  !< How many values should one buffer store
                                        !! i.e. how many samples of a time trace.
    integer, intent(in), optional :: nvalues_dim2, nvalues_dim3, nvalues_dim4
    integer                 :: init     !< Return value, 0=Success
    integer                 :: buffer_size

    init = -1

    if (.not. present(nvalues_dim2) .and. .not. present(nvalues_dim3) &
            .and. .not. present(nvalues_dim4)) then
        allocate(this%val_1d(nvalues_dim1, nbuffer))
        this%ndim = 1
        this%nvalues_dim1 = nvalues_dim1
        this%nvalues_dim2 = 0
        this%nvalues_dim3 = 0
        this%nvalues_dim4 = 0
        this%val_1d = 0
        buffer_size = nbuffer * nvalues_dim1

    elseif (present(nvalues_dim2) .and. .not. present(nvalues_dim3) &
            .and. .not. present(nvalues_dim4)) then
        allocate(this%val_2d(nvalues_dim1, nvalues_dim2, nbuffer))
        this%ndim = 2
        this%nvalues_dim1 = nvalues_dim1
        this%nvalues_dim2 = nvalues_dim2
        this%nvalues_dim3 = 0
        this%nvalues_dim4 = 0
        this%val_2d = 0
        buffer_size = nbuffer * nvalues_dim1 * nvalues_dim2

    elseif (present(nvalues_dim2) .and. present(nvalues_dim3) &
            .and. .not. present(nvalues_dim4)) then
        allocate(this%val_3d(nvalues_dim1, nvalues_dim2, nvalues_dim3, nbuffer))
        this%ndim = 3
        this%nvalues_dim1 = nvalues_dim1
        this%nvalues_dim2 = nvalues_dim2
        this%nvalues_dim3 = nvalues_dim3
        this%nvalues_dim4 = 0
        this%val_3d = 0
        buffer_size = nbuffer * nvalues_dim1 * nvalues_dim2 * nvalues_dim3

    elseif (present(nvalues_dim2) .and. present(nvalues_dim3) &
            .and. present(nvalues_dim4)) then
        allocate(this%val_4d(nvalues_dim1, nvalues_dim2, nvalues_dim3, nvalues_dim4, nbuffer))
        this%ndim = 4
        this%nvalues_dim1 = nvalues_dim1
        this%nvalues_dim2 = nvalues_dim2
        this%nvalues_dim3 = nvalues_dim3
        this%nvalues_dim4 = nvalues_dim4
        this%val_4d = 0
        buffer_size = nbuffer * nvalues_dim1 * nvalues_dim2 * nvalues_dim3 * nvalues_dim4

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
    if (allocated(this%val_1d)) deallocate(this%val_1d)
    if (allocated(this%val_2d)) deallocate(this%val_2d)
    if (allocated(this%val_3d)) deallocate(this%val_3d)
    if (allocated(this%val_4d)) deallocate(this%val_4d)
    deallocate(this%idx)

    this%nvalues_dim1 = 0
    this%nvalues_dim2 = 0
    this%nvalues_dim3 = 0
    this%nvalues_dim4 = 0
    this%nbuffer = 0 

    this%initialized = .false.
    freeme = 0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_1d(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues_dim1)
    integer                     :: get_1d      !< status value, 0=success
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
    get_1d = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_1d(:,ibuffer) 
       this%nhit = this%nhit + 1
       get_1d = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_2d(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2)
    integer                     :: get_2d      !< status value, 0=success
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
    get_2d = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_2d(:,:,ibuffer) 
       this%nhit = this%nhit + 1
       get_2d = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_3d(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                          this%nvalues_dim3)
    integer                     :: get_3d      !< status value, 0=success
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
    get_3d = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_3d(:,:,:,ibuffer) 
       this%nhit = this%nhit + 1
       get_3d = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_4d(this, iindex, values)

    class(buffer_type)          :: this
    integer, intent(in)         :: iindex   !< Index under which the value was stored.
                                            !! E.g. the index of the 
                                            !! point whose values are stored.
    real(kind=sp), intent(out)  :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                          this%nvalues_dim3, this%nvalues_dim4)
    integer                     :: get_4d      !< status value, 0=success
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
    get_4d = -1

    do ibuffer = 1, min(this%nput, this%nbuffer)
       if (this%idx(ibuffer).ne.iindex) cycle

       values = this%val_4d(:,:,:,:,ibuffer)
       this%nhit = this%nhit + 1
       get_4d = 0 
       exit
    end do

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_1d(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues_dim1) !< Values to store
    integer                   :: put_1d    !< Return value, 0=Success
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
       put_1d = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_1d(:,ibuffer) = values
       this%nput = this%nput + 1
       put_1d = 0
    end if

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_2d(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2) !< Values to store
    integer                   :: put_2d    !< Return value, 0=Success
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
       put_2d = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_2d(:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_2d = 0
    end if

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_3d(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                        this%nvalues_dim3) !< Values to store
    integer                   :: put_3d    !< Return value, 0=Success
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
       put_3d = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_3d(:,:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_3d = 0
    end if

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function put_4d(this, iindex, values)

    class(buffer_type)        :: this
    integer, intent(in)       :: iindex      !< Index under which the values can later
                                             !! be accessed
    real(kind=sp), intent(in) :: values(this%nvalues_dim1, this%nvalues_dim2, &
                                        this%nvalues_dim3, this%nvalues_dim4) !< Values to store
    integer                   :: put_4d    !< Return value, 0=Success
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
       put_4d = -1
    else
       if (this%nput < this%nbuffer) then
          ibuffer = this%nput + 1
       else
          call random_number(randtemp)
          ibuffer = int(randtemp*this%nbuffer) + 1
       endif
       this%idx(ibuffer) = iindex
       this%val_4d(:,:,:,:,ibuffer) = values
       this%nput = this%nput + 1
       put_4d = 0
    end if

end function
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
