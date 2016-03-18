!=========================================================================================
module heterogeneities

  use global_parameters
  use kdtree2_module
  use commpi, only   : pabort

  implicit none

  integer, parameter           :: nmodel_parameters_hetero = 3   !< Number of basic model parameters
                                                                 !! which may be used in these tests
                                                                 !! So far: vp, vs, rho

  character(len=3), parameter  :: parameter_name_het(nmodel_parameters_hetero) =  &
                                    ['vp ', 'vs ', 'rho']          ! to identify them for integration
  character(len=6), parameter  :: parameter_name_het_store(nmodel_parameters_hetero) =  &
                                    ['dlnvp ', 'dlnvs ', 'dlnrho'] ! how to store in netcdf

  type hetero_type
    private
    type(kdtree2), pointer, private  :: heterotree
    logical, private                 :: kdtree_built = .false.
    integer                          :: nhet
    real(kind=sp), allocatable       :: dlnhet(:,:)
    real(kind=sp), allocatable       :: coords(:,:)

    contains
      procedure, pass   :: load_het_rtpv
      procedure, pass   :: build_hetero_kdtree
      procedure, pass   :: load_model_coeffs
      procedure, nopass   :: get_hetparam_index
  end type hetero_type

contains


!-----------------------------------------------------------------------------------------
!> subroutine to load heterogeneities in rtpv format    
subroutine load_het_rtpv(this, het_file)

  use voxel, only : spherical_to_cartesian_point 
  character(len=*), intent(in)   :: het_file
  class(hetero_type)             :: this

  integer :: iinput_hetero
  integer :: i,ierr

  real(kind=dp) :: r,t,p    ! theta and phi (in degree)
  real(kind=dp) :: sph(3) ! sph point in r,t,p format, r is in [km]
  real(kind=dp) :: car(3) ! car point in x,y,z format

  open(newunit=iinput_hetero, file=trim(het_file), status='old', &
       action='read', iostat=ierr)

  if (ierr.ne.0) then
    print "('Heterogeneity file ', A, ' is missing')", trim(het_file)
    stop
  end if

  read(iinput_hetero,*) this%nhet
  allocate(this%coords(3,this%nhet))
  allocate(this%dlnhet(nmodel_parameters_hetero,this%nhet))

  ! read line by line
  do i=1,this%nhet
    ! r t p dlnvp dlnvs dlnrho
    read(iinput_hetero,*) r,t,p,this%dlnhet(:,i)

    if (abs(t)>90.0d0) then
      print "('ERROR in heterogeneity file ', A, ', l: ', I8)", trim(het_file), i
      print "('Value for theta (col. 2) is outside limits (-90,90): F9.2')", t
      stop
    end if

    sph(1) = r * 1000.d0 ! convert to [m]
    sph(2) = (90.d0 - t) * deg2rad ! convert to colatitude and rad
    sph(3) = p * deg2rad ! convert to rad
    ! routine takes r[m],phi[rad],theta[rad]         
    call spherical_to_cartesian_point(sph,car)
    this%coords(1,i) = real(car(1), kind=sp) ! x
    this%coords(2,i) = real(car(2), kind=sp) ! y
    this%coords(3,i) = real(car(3), kind=sp) ! z
  end do
  close(iinput_hetero)

end subroutine load_het_rtpv
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> sets up a 3d kdtree for the heterogeneity model
subroutine build_hetero_kdtree(this)
  use kdtree2_module, only    : kdtree2_create, kdtree2_destroy     
  class(hetero_type)         :: this
  ! Destroy kdtree
  if (this%kdtree_built) then
    if (verbose>0) then 
      print *, 'WARNING in build_kdtree(): Meshes have already been built'
      print *, 'Destroying the old trees...'
    end if
    call kdtree2_destroy(this%heterotree)
  end if      
  ! KDtree in forward field
  this%heterotree => kdtree2_create(this%coords,       &
                                    dim = 3,           &
                                    sort = .true.,     &
                                    rearrange = .true.) 

  this%kdtree_built = .true.

end subroutine build_hetero_kdtree
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> Loads the model coefficients for a selected coordinate from the kdtree
function load_model_coeffs(this, coordinates_xyz, weights) result(weighted_coeffs)

  use kdtree2_module, only           : kdtree2_result, kdtree2_n_nearest     

  class(hetero_type)                :: this
  real(kind=dp), intent(in)         :: coordinates_xyz(:,:)      
  real(kind=dp), intent(in)         :: weights(:)

  type(kdtree2_result), allocatable :: nextpoint(:)
  integer, parameter                :: nnext_points = 1
  integer                           :: npoints, ipoint
  integer                           :: pointid
  real(kind=dp)                     :: coeffs(nmodel_parameters_hetero, size(weights,1))
  real(kind=dp)                     :: weighted_coeffs(nmodel_parameters_hetero, size(weights,1))

  if (.not.this%kdtree_built) then
    print *, 'ERROR: KDTree is not built yet. Call build_kdtree before loading points!'
    call pabort()
  end if

  if (size(weights,1).ne.size(coordinates_xyz,2)) then
    print *, 'ERROR in load_model_coeffs: Number of weights is not equal number of points'
    print *, 'size(weights,1)         : ', size(weights,1)
    print *, 'size(coordinates_xyz,2) : ', size(coordinates_xyz,2)
  end if

  ! nextpoint has to be allocatable in kdtree module
  allocate(nextpoint(nnext_points))      
  npoints = size(coordinates_xyz, 2)

  ! load points from hetero tree
  do ipoint = 1, npoints
    ! Find nearest point
    call kdtree2_n_nearest( this%heterotree,                  &
                            real(coordinates_xyz(:, ipoint)), &
                            nn = nnext_points,                &
                            results = nextpoint )
    pointid = nextpoint(1)%idx         
    ! Get 6 model coefficients of nearest point from Mesh
    coeffs(:, ipoint) = this%dlnhet(:, pointid)                  
  end do ! ipoint

  ! apply weights
  do ipoint = 1, npoints
    weighted_coeffs(:, ipoint) = coeffs(:, ipoint) * weights(ipoint)
  end do

end function load_model_coeffs
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_hetparam_index(model_parameter)
  use commpi, only   : pabort
  character(len=3)  :: model_parameter
  integer           :: get_hetparam_index      
  integer           :: iparam

  ! Determine the index of the model parameter in the list defined in background_model.f90
  do iparam = 1, nmodel_parameters_hetero
    if (model_parameter == parameter_name_het(iparam)) then
      get_hetparam_index = iparam
      exit
    end if
  end do
  if (iparam == nmodel_parameters_hetero + 1) then
    print '("WARNING: Unknown hetero parameter: ", A)', &
      trim(model_parameter)
    !print '("Available options: ", 10(A3))', parameter_name_het
    print *, "This kernel will not be integrated over at the end"
    get_hetparam_index = -1 !call pabort(do_traceback=.false.)
  end if

end function get_hetparam_index
!-----------------------------------------------------------------------------------------

end module heterogeneities
!=========================================================================================
