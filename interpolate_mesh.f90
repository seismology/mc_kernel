module interpolate_mesh

  use global_parameters, only             : sp, dp
  implicit none
  private
  type parameter_interpolator
      real(kind=sp), allocatable         :: values(:)
      integer                            :: ndepth
      real(kind=dp)                      :: dr
      contains 
        private
        procedure, pass                  :: init
        procedure, pass, public          :: get
  end type

  public parameter_interpolator, create_interpolator

contains

!-----------------------------------------------------------------------------------------
function create_interpolator(param_tmp, tree, radius) result(interpolator)
  use kdtree2_module, only          :  kdtree2_result, kdtree2_n_nearest, kdtree2
  real(kind=sp), intent(in)         :: param_tmp(:)
  type(kdtree2), pointer            :: tree
  real(kind=dp), intent(in)         :: radius
  type(parameter_interpolator)      :: interpolator

  real(kind=dp), parameter          :: dz = 0.1d3 ! 100m
  real(kind=sp)                     :: coordinates_sz(2)
  real(kind=sp), allocatable        :: values(:)
  integer                           :: idepth, ndepth, pointid
  type(kdtree2_result), allocatable :: nextpoint(:)
  integer, parameter                :: nnext_points = 1

  ndepth = radius / dz + 1
  allocate(values(ndepth))

  ! nextpoint has to be allocatable in kdtree module
  allocate(nextpoint(nnext_points))

  do idepth = 1, ndepth
    coordinates_sz(1) = 0
    coordinates_sz(2) = (idepth - 1) * dz
    call kdtree2_n_nearest( tree,                &
                            coordinates_sz,      &  
                            nn = nnext_points,   &
                            results = nextpoint )

    pointid = nextpoint(1)%idx
    values(idepth) = param_tmp(pointid)
  end do

  call interpolator%init(values, dz)  

end function create_interpolator
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init(this, values, dr)
  class(parameter_interpolator)       :: this 
  real(kind=sp), intent(in)           :: values(:)
  real(kind=dp)                       :: dr

  this%ndepth = size(values)
  allocate(this%values(this%ndepth))
  this%values = values

  this%dr = dr

end subroutine init
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get(this, r) result(values)
  class(parameter_interpolator)   :: this
  real(kind=dp), intent(in)       :: r(:)
  real(kind=dp), allocatable      :: values(:)

  integer                         :: idx(size(r))

  idx = r/this%dr + 1

  ! Some meshes have points outside of the earth
  ! Set values to outermost layer in this case
  where (idx.gt.this%ndepth) idx=this%ndepth

  if (any(idx.le.0)) then
    print *, 'Depth out of range!'
    stop
  end if

  values = this%values(idx)

end function get
!-----------------------------------------------------------------------------------------

end module interpolate_mesh
