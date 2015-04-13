!=========================================================================================
module inversion_mesh

  use global_parameters, only: sp, dp, lu_out

  use tetrahedra,        only: get_volume_tet,                  &
                               get_volume_poly,                 &
                               get_center_tet,                  &
                               generate_random_points_tet,      &
                               generate_random_points_poly,     &
                               generate_random_points_ref_tri

  use voxel,             only: get_volume_vox,                  &
                               get_center_vox,                  &
                               generate_random_points_vox
                               
  use commpi,            only: pabort

  implicit none

  private
  public :: inversion_mesh_type
  public :: inversion_mesh_data_type
  public :: inversion_mesh_variable_type
  public :: plane_exp_pro2
  public :: append_variable

  type                              :: inversion_mesh_variable_type
    character(len=80)               :: var_name
    integer                         :: nentries
    logical                         :: istime
    character(len=80), allocatable  :: entry_names(:)
  end type


  type :: inversion_mesh_type
     private
     integer, public                    :: nvertices_per_elem
     integer, public                    :: nbasisfuncs_per_elem
     integer                            :: nelements, nvertices
     integer, allocatable               :: connectivity(:,:)
     integer, allocatable               :: block_id(:)
     real(kind=dp), allocatable         :: vertices(:,:)
     real(kind=dp), allocatable         :: v_2d(:,:,:), p_2d(:,:,:)
     real(kind=dp), allocatable         :: abinv(:,:,:)
     character(len=4)                   :: element_type
     logical                            :: initialized = .false.
     contains
     procedure, pass :: get_element_type
     procedure, pass :: get_nelements
     procedure, pass :: get_element
     procedure, pass :: get_nvertices
     procedure, pass :: get_nbasisfuncs
     procedure, pass :: get_vertices
     procedure, pass :: get_vertices_dp
     procedure, pass :: get_valence
     procedure, pass :: get_connected_elements
     procedure, pass :: get_connectivity
     procedure, pass :: read_tet_mesh
     procedure, pass :: read_abaqus_mesh
     procedure, pass :: read_abaqus_meshtype
     procedure, pass :: tree_sort
     procedure, pass :: dump_mesh_xdmf
     procedure, pass :: get_volume
     procedure, pass :: get_center
     procedure, pass :: get_model_coeff
     procedure, pass :: generate_random_points
     procedure, pass :: make_2d_vectors
     procedure, pass :: weights
     procedure, pass :: initialize_mesh
     procedure, pass :: freeme
     procedure       :: init_weight_tet_mesh
  end type

  type, extends(inversion_mesh_type)    :: inversion_mesh_data_type
     private
     integer                            :: ntimes_node
     integer                            :: ntimes_cell
     real(kind=sp), allocatable         :: datat_node(:,:)
     real(kind=sp), allocatable         :: datat_cell(:,:)
     integer, allocatable               :: group_id_node(:)
     integer, allocatable               :: group_id_cell(:)
     character(len=32), allocatable     :: data_group_names_node(:)
     character(len=32), allocatable     :: data_group_names_cell(:)
     integer                            :: ngroups_node
     integer                            :: ngroups_cell
     integer, public                    :: ncid      = -1
     integer                            :: ncid_cell = -1
     integer                            :: ncid_node = -1
     integer                            :: nvar_node = 0
     integer                            :: nvar_cell = 0
     integer                            :: ntimes    = 1
     real(kind=dp)                      :: dt, starttime
     character(len=21), public          :: filename_tmp_out
     type(inversion_mesh_variable_type), allocatable, public    :: variable_cell(:)
     type(inversion_mesh_variable_type), allocatable, public    :: variable_node(:)

     contains
     procedure, pass :: get_ntimes_node
     procedure, pass :: get_ntimes_cell

     ! for xdmf format dumps
     procedure, pass :: init_node_data
     procedure, pass :: init_cell_data
     procedure, pass :: add_node_variable
     procedure, pass :: add_node_data
     procedure, pass :: add_cell_variable
     procedure, pass :: add_cell_data
     procedure, pass :: dump_data_xdmf

     procedure, pass :: free_node_and_cell_data

     ! for csr format dumps
     procedure, pass :: dump_node_data_csr
     procedure, pass :: dump_cell_data_csr
  
     ! for ascii format dumps
     procedure, pass :: dump_node_data_ascii
     procedure, pass :: dump_cell_data_ascii

  end type

contains

!-----------------------------------------------------------------------------------------
integer function get_element_type(this)
  class(inversion_mesh_type)        :: this

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  select case(this%element_type)
  case('tri')
     get_element_type = 1
  case('quad')
     get_element_type = 2
  case('hex')
     get_element_type = 3
  case('tet')
     get_element_type = 4
  case('vox')
     get_element_type = 5
  case default
     get_element_type = -1
     call pabort
  end select
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_nelements(this)
  class(inversion_mesh_type)        :: this

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  get_nelements = this%nelements
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_element(this, ielement)
  class(inversion_mesh_type)        :: this
  real(kind=dp)                     :: get_element(3,this%nvertices_per_elem)
  integer, intent(in)               :: ielement
  integer                           :: ivert

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  do ivert= 1, this%nvertices_per_elem
     get_element(:,ivert) = this%vertices(:, this%connectivity(ivert,ielement)+1) 
  enddo
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_nvertices(this)
  class(inversion_mesh_type)        :: this

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  get_nvertices = this%nvertices
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_nbasisfuncs(this, int_type)
  class(inversion_mesh_type)        :: this
  character(len=32)                 :: int_type
  
  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  select case(trim(int_type))
  case('onvertices')
     get_nbasisfuncs = this%nvertices
  case('volumetric')
     get_nbasisfuncs = this%nelements
  case default
     write(*,'(A)') 'ERROR: unknown int_type caused crash in get_nbasisfuncs!'
     get_nbasisfuncs = -1
     call pabort
  end select

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_vertices(this)
  class(inversion_mesh_type)        :: this
  real(kind=sp)                     :: get_vertices(3,this%nvertices)

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  get_vertices = this%vertices
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_vertices_dp(this)
  class(inversion_mesh_type)        :: this
  real(kind=dp)                     :: get_vertices_dp(3,this%nvertices)

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  get_vertices_dp = this%vertices
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_valence(this, ivert)
  class(inversion_mesh_type)        :: this
  integer                           :: get_valence
  integer, intent(in)               :: ivert
  
  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  get_valence = count(this%connectivity == ivert - 1) ! -1 because counting in
                                                      ! the mesh starts from 0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_connected_elements(this, ivert)
  class(inversion_mesh_type)        :: this
  integer, allocatable              :: get_connected_elements(:)
  integer, intent(in)               :: ivert
  integer                           :: ielement, ct
  
  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  allocate(get_connected_elements(this%get_valence(ivert)))
  get_connected_elements = -1

  ct = 1
  do ielement=1, this%nelements
     if (any(this%connectivity(:,ielement) == ivert - 1)) then
        get_connected_elements(ct) = ielement
        ct = ct + 1
     endif
  enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_connectivity(this)
  class(inversion_mesh_type)        :: this
  integer                           :: get_connectivity(this%nvertices_per_elem, &
                                                        this%nelements)
  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if
  get_connectivity = this%connectivity + 1
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntimes_node(this)
  class(inversion_mesh_data_type)        :: this

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if
  get_ntimes_node = this%ntimes_node
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntimes_cell(this)
  class(inversion_mesh_data_type)        :: this

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if
  get_ntimes_cell = this%ntimes_cell
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_volume(this, ielement)
  class(inversion_mesh_type)        :: this
  integer, intent(in)               :: ielement
  real(kind=dp)                     :: get_volume
  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  get_volume = 0

  select case(this%element_type)
  case('tet')
     get_volume = get_volume_tet(this%get_element(ielement))
  case('quad')
     get_volume = get_volume_poly(4, this%get_element(ielement))
  case('tri')
     get_volume = get_volume_poly(3, this%get_element(ielement))
  case('vox')
     get_volume = get_volume_vox(this%get_element(ielement))
  case('hex')
     !get_volume = get_volume_hex(this%get_element(ielement))
  end select

end function get_volume
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_center(this, ielement)
  class(inversion_mesh_type)        :: this
  integer, intent(in)               :: ielement
  real(kind=dp)                     :: get_center(3)
  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort
  end if

  get_center = 0

  select case(this%element_type)
  case('tet')
     get_center = get_center_tet(this%get_element(ielement))
!  case('quad')
!    get_center = get_center_poly(4, this%get_element(ielement))
!  case('tri')
!    get_center = get_center_poly(3, this%get_element(ielement))
  case('vox')
     get_center = get_center_vox(this%get_element(ielement))
!  case('hex')
!    get_center = get_center_hex(this%get_element(ielement))
  case default
     write(*,'(A,A,A)') 'ERROR: get_center for element type ', this%element_type, ' not implemented'
     call pabort 
  end select

end function get_center
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_model_coeff(this, ielement, modelcoeffid, int_type)
!
!< Depending on the int_type (volumetric or vertex based mode) returns
! model coefficient with id modelcoeffid at the center or the vertices of an element
!
  class(inversion_mesh_type)        :: this
  integer, intent(in)               :: ielement
  character(len=*), intent(in)      :: modelcoeffid
  character(len=*), intent(in)      :: int_type
  real(kind=dp), allocatable        :: get_model_coeff(:)
  real(kind=dp), allocatable        :: vertices(:,:)
  real(kind=dp)                     :: center(3) 
  integer                           :: ivertex

  ! @ TODO: This routine is not quite finished but may become useful
  ! if one wants to export model coefficients for plotting

  select case(trim(int_type))
  case('onvertices')  
     allocate(get_model_coeff(this%nvertices))
     allocate(vertices(3,this%nvertices))
     vertices = this%get_element(ielement)
     do ivertex = 1,this%nvertices           
        get_model_coeff(ivertex) = 0.d0 ! @ TODO: load_coefficient(vertices(:,1),modelcoeffs_id(icoeff))
     end do
  case('volumetric')
     allocate(get_model_coeff(1))
     center = this%get_center(ielement)
     get_model_coeff(1) = 0.d0 ! @ TODO: load_coefficient(center,modelcoeffs_id(icoeff))
  end select
 
end function get_model_coeff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function weights(this, ielem, ivertex, points)
!< Calculates the weight that a value at location "points" has on a kernel on vertex 
!! "ivertex" of element "ielement". Quadrature rules come in here   
!! Based on Nolet, Breviary, 12.2, p.225
  use tetrahedra, only            : point_in_tetrahedron
  class(inversion_mesh_type)     :: this
  integer, intent(in)            :: ielem, ivertex
  real(kind=dp), intent(in)      :: points(:,:)
  real(kind=dp)                  :: weights(size(points,2))
  logical                        :: isonplane(size(points,2)), isintetrahedron(size(points,2))
  logical                        :: isdegenerate
  real(kind=dp)                  :: dx, dy, dz
  integer                        :: ipoint, i

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  if (this%nbasisfuncs_per_elem==1) then !Integration on cell
    weights = 1

  else !Integration on vertices
    select case(this%element_type) 
    case('tet')
       do ipoint = 1, size(points, 2)
         dx = points(1, ipoint) - this%vertices(1, this%connectivity(1, ielem) + 1)
         dy = points(2, ipoint) - this%vertices(2, this%connectivity(1, ielem) + 1)
         dz = points(3, ipoint) - this%vertices(3, this%connectivity(1, ielem) + 1)

         weights(ipoint) =   this%abinv(ivertex, 1, ielem) * dx &
                           + this%abinv(ivertex, 2, ielem) * dy &
                           + this%abinv(ivertex, 3, ielem) * dz &
                           + this%abinv(ivertex, 4, ielem) * 1

       end do

       if (any(weights<-1d-9)) then
         isintetrahedron = point_in_tetrahedron(this%vertices(:, this%connectivity(1, ielem) + 1), &
                                                points(:, :),                                      &
                                                isonplane, isdegenerate)

         print *, 'ERROR: weight is smaller zero! Check whether point is outside of element'
         print *, '       Element (local numbering): ', ielem, ', ivertex: ', ivertex
         print *, '       Vertices:  '
         do i = 1, 4
           print *, '       ', i, this%vertices(1:3, this%connectivity(i, 1) + 1)
         end do
         print *, '       Volume:   ', this%get_volume(ielem)
         print *, '       Degener:  ', isdegenerate
         print *, ''

         do ipoint = 1, size(points, 2)
           dx = points(1, ipoint) - this%vertices(1, this%connectivity(1, ielem) + 1)
           dy = points(2, ipoint) - this%vertices(2, this%connectivity(1, ielem) + 1)
           dz = points(3, ipoint) - this%vertices(3, this%connectivity(1, ielem) + 1)

           print *, '       Point:    ', ipoint
           print *, '                 ', points(:, ipoint)
           print *, '       isintet?  ', isintetrahedron(ipoint)
           print *, '       isonsurf? ', isonplane(ipoint)
           print *, '       weight:   ', weights(ipoint)
           print *, '       dx,dy,dz: ', dx, dy, dz
           print *, ''
         end do
         call pabort()
       end if
       if (any(weights>1.d0+1d-9)) then
         print *, 'ERROR: weight is larger one! Check whether point is outside of element'
         print *, '       Element (local numbering): ', ielem, ', ivertex: ', ivertex
         print *, '       Vertices:'
         do i = 1, 4
           print *, '       ', i, this%vertices(1:3, this%connectivity(i, 1) + 1)
         end do
         print *, '       Volume: ', this%get_volume(ielem)
         print *, ''
         do ipoint = 1, size(points, 2)
           dx = points(1, ipoint) - this%vertices(1, this%connectivity(1, ielem) + 1)
           dy = points(2, ipoint) - this%vertices(2, this%connectivity(1, ielem) + 1)
           dz = points(3, ipoint) - this%vertices(3, this%connectivity(1, ielem) + 1)

           print *, '       Point:    ', ipoint
           print *, '                 ', points(:, ipoint)
           print *, '       weight:   ', weights(ipoint)
           print *, '       dx,dy,dz: ', dx, dy, dz
           print *, ''
         end do
         call pabort()
       end if

    case default
        ! Least sophisticated version possible
        weights = 1
    end select
  end if

end function weights
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function generate_random_points(this, ielement, npoints, quasirandom) result(points)
  class(inversion_mesh_type)     :: this
  integer, intent(in)            :: ielement, npoints
  logical, intent(in)            :: quasirandom
  real(kind=dp)                  :: points(3, npoints)
  real(kind=dp)                  :: points2d(2, npoints)
  integer                        :: ipoint

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if

  select case(this%element_type)
  case('tet')
     points = generate_random_points_tet(this%get_element(ielement), npoints,  &
                                         quasirandom )
  case('quad')
     points2d = generate_random_points_poly(4, this%p_2d(:,:,ielement), &
                                            npoints, quasirandom )
     do ipoint = 1, npoints
        points(:,ipoint) =   this%v_2d(:,0,ielement)                         &
                           + this%v_2d(:,1,ielement) * points2d(1,ipoint)    &
                           + this%v_2d(:,2,ielement) * points2d(2,ipoint)
     end do
  case('tri')
     points2d = generate_random_points_ref_tri(npoints, quasirandom)
     do ipoint = 1, npoints
        points(:,ipoint) =   this%v_2d(:,0,ielement)                         &
                           + this%v_2d(:,1,ielement) * points2d(1,ipoint)    &
                           + this%v_2d(:,2,ielement) * points2d(2,ipoint)
     end do
  case('hex')
     ! points  = generate_random_points_hex( this%get_element(ielement), npoints)
     call pabort()
  case('vox')
     points = generate_random_points_vox( this%get_element(ielement), npoints)
  case default
     call pabort()
  end select

end function generate_random_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_tet_mesh(this, filename_vertices, filename_connectivity, int_type)
  class(inversion_mesh_type)        :: this
  character(len=*), intent(in)      :: filename_vertices
  character(len=*), intent(in)      :: filename_connectivity
  character(len=*), intent(in)      :: int_type                  !< Integrate on nodes or cells? 
  integer                           :: iinput_vertices, iinput_connectivity
  integer                           :: i, ierr, ndim

  if (this%initialized) then
     write(*,'(A)') 'ERROR: Trying to read mesh into inversion mesh type that is already initialized'
     call pabort(do_traceback = .false.)
  end if

  this%nvertices_per_elem = 4
  this%element_type = 'tet'
  
  ! read vertices
  open(newunit=iinput_vertices, file=trim(filename_vertices), status='old', &
       action='read', iostat=ierr)
  if ( ierr /= 0 ) then
     write(*,*) 'ERROR: Could not open vertex file: ', trim(filename_vertices)
     write(*,*) 'ierr :', ierr
     call pabort(do_traceback = .false.)
  endif

  read(iinput_vertices,*) ndim ! The dimensionality of the mesh can be changed here
                               ! 3: means tetrahedral meshes
                               ! 2: means probably triangles. This should not be used 
                               !    anymore and is not implemented here, at least for 
                               !    now. SCS September 14

  if (ndim.ne.3) then
     write(*,*) 'ERROR: the dimension of Sigloch/Nolet style meshes should be three (3)'
     write(*,*) '       Defined in the first line of the vertex files'
     call pabort(do_traceback = .false.)
  end if    

  read(iinput_vertices,*) this%nvertices

  allocate(this%vertices(3,this%nvertices))
  do i=1, this%nvertices
     read(iinput_vertices,*) this%vertices(:,i)
  enddo

  close(iinput_vertices)

  ! read connectivity (tetrahedral elements)
  open(newunit=iinput_connectivity, file=trim(filename_connectivity), status='old', &
       action='read', iostat=ierr)
  if ( ierr /= 0 ) then
     write(*,*) 'ERROR: Could not open connectivity file: ', trim(filename_connectivity)
     write(*,*) 'ierr :', ierr
     call pabort(do_traceback = .false.)
  endif

  read(iinput_connectivity,*) this%nelements

  allocate(this%connectivity(this%nvertices_per_elem,this%nelements))
  do i=1, this%nelements
     read(iinput_connectivity,*) this%connectivity(:,i)
  enddo
  close(iinput_connectivity)

  allocate(this%block_id(this%nelements))
  this%block_id = 1
  
  write(lu_out,'(A)')       ' Tetrahedral mesh read in...'
  write(lu_out,'(A,A)')     '  vertex file       : ', trim(filename_vertices)
  write(lu_out,'(A,A)')     '  connectivity file : ', trim(filename_connectivity)
  write(lu_out,'(A,I10,A)') '  contains ', this%nvertices, ' vertices '
  write(lu_out,'(A,I10,A)') '  and      ', this%nelements, ' elements '
  write(lu_out,'(A)')       ''


  ! Convert locations to meters, if necessary
  if (maxval(abs(this%vertices))<1D4) then
    write(lu_out, *) 'Assuming the vertex locations are in km!'
    this%vertices = this%vertices * 1D3 ! Vertices are set in km
  end if

  select case(trim(int_type))
  case('onvertices')
     this%nbasisfuncs_per_elem = this%nvertices_per_elem    
  case('volumetric')
     this%nbasisfuncs_per_elem = 1     
  end select

  ! Initialize weigths and calculate base vectors for each element
  call this%init_weight_tet_mesh()

  this%initialized = .true.
end subroutine read_tet_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine initialize_mesh(this, ielem_type, vertices, connectivity, nbasisfuncs_per_elem) 
  class(inversion_mesh_type)        :: this
  integer                          :: ielem_type !1-tet, 2-quad, 3-tri, 4-hex
  integer                          :: nbasisfuncs_per_elem
  real(kind=dp),    intent(in)     :: vertices(:,:)
  integer,          intent(in)     :: connectivity(:,:)
  character(len=255)               :: fmtstring

  if (this%initialized) then
     write(*,'(A)') 'ERROR: Trying to initialize inversion mesh type that is already initialized'
     call pabort 
  end if
  
  this%nvertices = size(vertices,2)
  this%nelements = size(connectivity,2)

  fmtstring = "('  Initialize mesh with ', I5, ' vertices and ', I5, ' elements')"
  write(lu_out,fmtstring) this%nvertices, this%nelements

  select case(ielem_type)
  case(1) ! tri
     this%nvertices_per_elem = 3
     this%element_type = 'tri'
  case(2) ! quad
     this%nvertices_per_elem = 4
     this%element_type = 'quad'
  case(3) ! hex
     this%nvertices_per_elem = 8
     this%element_type = 'hex'
  case(4) ! tet
     this%nvertices_per_elem = 4
     this%element_type = 'tet'
  case(5) ! vox
     this%nvertices_per_elem = 8
     this%element_type = 'vox'
  case default
     write(6,*) 'ERROR: Initializing with elementtype ', ielem_type, &
                'not yet implemented'
     call pabort
  end select

  if (this%nvertices_per_elem.ne.size(connectivity,1)) then
      write(*,*) 'ERROR at initialize_mesh:'
      write(*,*) 'Wrong number of vertices per element for type ', trim(this%element_type)
      write(*, '(A,I5,A,I2)') 'is: ', size(connectivity,1), ', should be: ', this%nvertices_per_elem
      call pabort
  end if

  ! how many basis functions per elem
  this%nbasisfuncs_per_elem = nbasisfuncs_per_elem    

  ! prepare arrays
  allocate(this%vertices(3,this%nvertices))
  allocate(this%connectivity(this%nvertices_per_elem,this%nelements))
  allocate(this%block_id(this%nelements))

  this%vertices = vertices
  this%connectivity = connectivity - 1
  this%block_id = 1

  if (this%element_type.eq.'tri'.or.this%element_type.eq.'quad') then
     call this%make_2d_vectors
  elseif (this%element_type.eq.'tet') then
     call this%init_weight_tet_mesh()
  end if

  this%initialized = .true.

end subroutine initialize_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_abaqus_meshtype(this, filename, int_type)
  class(inversion_mesh_type)        :: this
  character(len=*), intent(in)      :: filename
  character(len=*), intent(in)      :: int_type
  integer                           :: iinput
  integer                           :: ierr, ct
  character(len=128)                :: line
  character(len=16)                 :: elem_type

  if (this%initialized) then
     write(*,'(A)') 'ERROR: Trying to initialize inversion mesh type that is already initialized'
     call pabort 
  end if

  ! open file 
  open(newunit=iinput, file=trim(filename), status='old', &
       action='read', iostat=ierr)
  if ( ierr /= 0 ) then
     write(*,*) 'ERROR: Could not open file: ', trim(filename)
     write(*,*) 'ierr :', ierr
     call pabort
  endif

  ! go to vertex block
  do
    read(iinput,*) line
    if (index(line, '*NODE') == 1) exit
  enddo

  ! scan number of vertices
  ct = 0
  do
    read(iinput,*) line
    if (index(line, '**') == 1) exit
    ct = ct + 1
  enddo
  this%nvertices = ct
  
  ! go to element block
  do
    read(iinput,*) line, elem_type
    if (index(line, '*ELEMENT') == 1) exit
  enddo
 
  ! get element type
  elem_type = trim(elem_type(6:))

  select case (trim(elem_type))
  case('STRI3')
     this%nvertices_per_elem = 3
     this%element_type = 'tri'
  case('S4R')
     this%nvertices_per_elem = 4
     this%element_type = 'quad'
  case('C3D8R')
     this%nvertices_per_elem = 8
     this%element_type = 'hex'
  case('C3D4')
     this%nvertices_per_elem = 4
     this%element_type = 'tet'
  case('VOX')
     this%nvertices_per_elem = 8
     this%element_type = 'vox'
  case default
     write(6,*) 'ERROR: reading abaqus file with elementtype ', trim(elem_type), &
                'not yet implemented'
     call pabort
  end select
  close(iinput)

  select case(trim(int_type))
  case('onvertices')
     this%nbasisfuncs_per_elem = this%nvertices_per_elem    
  case('volumetric')
     this%nbasisfuncs_per_elem = 1     
  end select

end subroutine read_abaqus_meshtype
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_abaqus_mesh(this, filename, int_type)
  class(inversion_mesh_type)        :: this
  character(len=*), intent(in)      :: filename
  character(len=*), intent(in)      :: int_type
  integer                           :: iinput
  integer                           :: i, ierr, ct, ct_block
  character(len=128)                :: line, line_buff
  character(len=16)                 :: elem_type, elem_type_buff

  if (this%initialized) then
     write(*,'(A)') 'ERROR: Trying to initialize inversion mesh type that is already initialized'
     call pabort 
  end if

  write(lu_out, *) 'Reading Mesh file from ', trim(filename)

  ! open file 
  open(newunit=iinput, file=trim(filename), status='old', &
       action='read', iostat=ierr)
  if ( ierr /= 0 ) then
     write(*,*) 'ERROR: Could not open file: ', trim(filename)
     write(*,*) 'ierr :', ierr
     call pabort
  endif

  ! go to vertex block
  do
    read(iinput,*) line
    if (index(line, '*NODE') == 1) exit
  enddo

  ! scan number of vertices
  ct = 0
  do
    read(iinput,*) line
    if (index(line, '**') == 1) exit
    ct = ct + 1
  enddo
  this%nvertices = ct
  
  ! go to element block
  do
    read(iinput,*) line, elem_type
    if (index(line, '*ELEMENT') == 1) exit
  enddo
 
  ! get element type
  elem_type = trim(elem_type(6:))

  select case (trim(elem_type))
  case('STRI3')
     this%nvertices_per_elem = 3
     this%element_type = 'tri'
  case('S4R')
     this%nvertices_per_elem = 4
     this%element_type = 'quad'
  case('C3D8R')
     this%nvertices_per_elem = 8
     this%element_type = 'hex'
  case('C3D4')
     this%nvertices_per_elem = 4
     this%element_type = 'tet'
  case('VOX')
     this%nvertices_per_elem = 8
     this%element_type = 'vox'
  case default
     write(6,*) 'ERROR: reading abaqus file with elementtype ', trim(elem_type), &
                'not yet implemented'
     call pabort
  end select

  ! scan number of elements
  ct = 0
  do
    read(iinput,'(A128)') line
    if (index(line, '*ELEMENT') == 1) then
        read(line,*) line_buff, elem_type_buff
        elem_type_buff = trim(elem_type_buff(6:))
        if (trim(elem_type_buff) /= trim(elem_type)) then
           write(6,*) 'ERROR: reading abaqus file with multiple elementtypes (', &
                      trim(elem_type), ', ', trim(elem_type_buff), &
                      ') not yet implemented'
           call pabort
        endif
        cycle
    endif
    if (index(line, '**') == 1) exit
    ct = ct + 1
  enddo
  this%nelements = ct

  ! prepare arrays
  allocate(this%vertices(3,this%nvertices))
  allocate(this%connectivity(this%nvertices_per_elem,this%nelements))
  allocate(this%block_id(this%nelements))
  this%block_id = 1

  ! reset to beginning of file
  rewind(iinput)
  
  ! go to node block
  do
    read(iinput,*) line
    if (index(line, '*NODE') == 1) exit
  enddo
 
  ! read vertices
  do i=1, this%nvertices
    read(iinput,*) line, this%vertices(:,i)
  enddo

  ! go to element block
  do
    read(iinput,*) line, elem_type
    if (index(line, '*ELEMENT') == 1) exit
  enddo

  ! read elements
  ct = 0
  ct_block = 1
  do 
     read(iinput,'(A128)') line_buff
     if (index(line_buff, '*ELEMENT') == 1) then
        ct_block = ct_block + 1
        cycle
     endif
     ct = ct + 1
     read(line_buff,*) line, this%connectivity(:,ct)
     this%block_id(ct) = ct_block
     if (ct == this%nelements) exit
  enddo

  close(iinput)

  ! abaqus starts counting at 1
  this%connectivity(:,:) = this%connectivity(:,:) - 1
  
  ! Convert locations to meters, if necessary
  if (maxval(abs(this%vertices))<1D4) then
    write(lu_out, *) 'Assuming the vertex locations are in km!'
    this%vertices = this%vertices * 1D3 ! Vertices are set in km
  end if

  ! Initialize weigths and calculate base vectors for each element
  if (this%element_type.eq.'tri'.or.this%element_type.eq.'quad') then
     call this%make_2d_vectors
  elseif (this%element_type.eq.'tet') then
     call this%init_weight_tet_mesh()
  end if

  write(lu_out, '("  Mesh type: ", A)')  trim(this%element_type)
  write(lu_out, '("  nvertices: ", I8)') this%nvertices
  write(lu_out, '("  nelements: ", I8)') this%nelements

  ! how many basis functions per elem
  select case(trim(int_type))
  case('onvertices')
     this%nbasisfuncs_per_elem = this%nvertices_per_elem    
  case('volumetric')
     this%nbasisfuncs_per_elem = 1     
  end select

  this%initialized = .true.
end subroutine read_abaqus_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_weight_tet_mesh(this)
  class(inversion_mesh_type)       :: this
  real(kind=dp)                    :: x1, y1, z1, ab(4,4)
  integer                          :: ielem, ivertex

  if (this%initialized) then ! BUG REMOVED LA JUNE 15, removed .not., 
     write(*,'(A)') 'ERROR: Trying to initialize inversion mesh type that is already initialized'
     call pabort 
  end if


  allocate(this%abinv(4,4,this%nelements))
  do ielem = 1, this%nelements
     x1 = this%vertices(1, this%connectivity(1,ielem)+1)
     y1 = this%vertices(2, this%connectivity(1,ielem)+1)
     z1 = this%vertices(3, this%connectivity(1,ielem)+1)


     do ivertex = 1, 4
        ab(1, ivertex) = this%vertices(1, this%connectivity(ivertex,ielem)+1) - x1
        ab(2, ivertex) = this%vertices(2, this%connectivity(ivertex,ielem)+1) - y1
        ab(3, ivertex) = this%vertices(3, this%connectivity(ivertex,ielem)+1) - z1
        ab(4, ivertex) = 1

     end do

     this%abinv(:,:,ielem) = invert(ab, 4)

  end do

end subroutine init_weight_tet_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function invert(A, nrows) result(A_inv)
!< Inverts a square matrix of dimension nrows using LAPACK routines
!use lapack95

  integer, intent(in)           :: nrows 
  real(kind=dp), intent(in)     :: A(nrows, nrows)
  real(kind=dp)                 :: A_inv(nrows, nrows)
  integer                       :: info
  real(kind=dp), allocatable    :: work(:)
  integer, allocatable          :: ipiv(:)

  external                      :: dgetrf, dgetri
  ! This is done in LAPACK routines
  A_inv = A
  allocate(ipiv(nrows))
  allocate(work(nrows))
  
  call dgetrf( nrows, nrows, A_inv, nrows, ipiv, info )               

  if (info.ne.0) then
    write(*,*) 'Matrix is numerically singular'
    call pabort()
  end if
  
  ! Inverse der LU-faktorisierten Matrix A        
  call dgetri( nrows, A_inv, nrows, ipiv, work, nrows, info )
  if (info.ne.0) then
     write(*,*) 'Error in matrix inversion'
     call pabort() 
  end if

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine make_2d_vectors(this)
  class(inversion_mesh_type)        :: this
  integer                           :: nvec, ielem

  select case(trim(this%element_type))
  case('tri')
     nvec = 2
     
     allocate(this%v_2d(3,  0:2, this%nelements))
     allocate(this%p_2d(2, nvec, this%nelements))

     do ielem = 1, this%nelements
        this%v_2d(:,0,ielem) = this%vertices(:, this%connectivity(1,ielem)+1)
        this%v_2d(:,1,ielem) = this%vertices(:, this%connectivity(2,ielem)+1) - &
                               this%vertices(:, this%connectivity(1,ielem)+1)
        this%v_2d(:,2,ielem) = this%vertices(:, this%connectivity(3,ielem)+1) - &
                               this%vertices(:, this%connectivity(1,ielem)+1)
     end do
    
     this%p_2d = 0

  case('quad') 
     nvec = 3
  
     allocate(this%v_2d(3,  0:2, this%nelements))
     allocate(this%p_2d(2, nvec, this%nelements))

     do ielem = 1, this%nelements
        this%v_2d(:,0,ielem) = this%vertices(:, this%connectivity(1,ielem)+1)
        call plane_exp_pro2(p_ref   = this%vertices(:, this%connectivity(1:3, ielem)+1),      &
                            npoints = nvec,                                                   &
                            p_3d    = this%vertices(:, this%connectivity(2:nvec+1, ielem)+1), &
                            p_2d    = this%p_2d(:,:,ielem),                                   &
                            vec     = this%v_2d(:,1:2,ielem))
     end do
  case default
     call pabort
  end select
end subroutine make_2d_vectors
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine plane_exp_pro2 ( p_ref, npoints, p_3d, p_2d, vec )
! PLANE_EXP_PRO2 produces 2D coordinates of points that lie in a plane, in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is:
!
!      the plane through P1, P2 and P3.
!
!    The first thing to do is to compute two orthonormal vectors V1 and
!    V2, so that any point P that lies in the plane may be written as
!
!      P = P1 + alpha * V1 + beta * V2
!
!    The vector V1 lies in the direction P2-P1, and V2 lies in
!    the plane, is orthonormal to V1, and has a positive component
!    in the direction of P3-P1.
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Input, integer ( kind = 4 ) N, the number of points to project.
!
!    Input, real ( kind = 8 ) P(3,N), are the Cartesian
!    coordinates of points which lie on the plane spanned by the
!    three points.  These points are not checked to ensure that
!    they lie on the plane.
!
!    Output, real ( kind = 8 ) PP(2,N), the "in-plane"
!    coordinates of the points.  

  integer, parameter         :: dim_num = 3
  integer, intent(in)        :: npoints
  real(kind=dp), intent(in)  :: p_3d(dim_num,npoints)
  real(kind=dp), intent(in)  :: p_ref(dim_num,3)
  real(kind=dp), intent(out) :: p_2d(2, npoints)
  real(kind=dp), intent(out) :: vec(dim_num, 2)
  
  real(kind=dp)              :: dot
  integer                    :: i

  !Compute the two basis vectors for the affine plane.
  vec(:,1) = p_ref(:,2) - p_ref(:,1)
  vec(:,1) = vec(:,1) / norm2(vec(:,1))
  
  vec(:,2) = p_ref(:,3) - p_ref(:,1)
  dot = dot_product ( vec(:,1), vec(:,2) )
  vec(:,2) = vec(:,2) - dot * vec(:,1)
  vec(:,2) = vec(:,2) / norm2(vec(:,2))
  
  !Now decompose each point.
  do i = 1, npoints
     p_2d(1,i) = dot_product ( p_3d(:,i) - p_ref(:,1), vec(:,1) )
     p_2d(2,i) = dot_product ( p_3d(:,i) - p_ref(:,2), vec(:,2) )
  end do
  
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine tree_sort(this)
!< resorts the mesh such that consecutive elements are more likely to be close to each other
!  based on a KD-tree traversal, hence does the sorting in a split second for meshes with 
!  10^5 elements

  use kdtree2_module
  class(inversion_mesh_type)        :: this
  type(kdtree2), pointer            :: tree
  real(kind=sp)                     :: midpoints(3, this%nelements)
  integer                           :: connectivity_sorted(this%nvertices_per_elem, this%nelements)
  integer                           :: i

  ! compute element midpoints
  do i = 1, this%nelements
     midpoints(:,i) = this%get_center(i)
  enddo

  ! build a kdtree on element midpoints
  tree => kdtree2_create(midpoints, dim = 3, sort = .true., rearrange = .true.)

  ! traverse the tree
  do i = 1, this%nelements
     connectivity_sorted(:,i) = this%connectivity(:, tree%ind(i))
  enddo

  ! update connectivity
  this%connectivity(:,:) = connectivity_sorted(:,:)

  call kdtree2_destroy(tree)
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_mesh_xdmf(this, filename)
  class(inversion_mesh_type)        :: this
  character(len=*), intent(in)      :: filename
  integer                           :: iinput_xdmf, iinput_heavy_data
  character(len=16)                 :: xdmf_elem_type
  character(len=512)                :: filename_np

  if (.not. this%initialized) then
     write(*,*) 'ERROR: trying to dump a non initialized mesh'
     call pabort 
  end if

  ! relative filename for xdmf content
  filename_np = trim(filename(index(filename, '/', back=.true.)+1:))

  select case(this%element_type)
  case('tri')
     xdmf_elem_type = 'Triangle'
  case('tet')
     xdmf_elem_type = 'Tetrahedron'
  case('quad')
     xdmf_elem_type = 'Quadrilateral'
  case('hex')
     xdmf_elem_type = 'Hexahedron'
  case('vox')
     xdmf_elem_type = 'Hexahedron'
  case default
     write(6,*) 'ERROR: xdmf dumping for element type ', this%element_type, &
                ' not implemented'
     call pabort
  end select

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  write(iinput_xdmf, 732) trim(xdmf_elem_type), this%nelements, this%nelements, &
                      this%nvertices_per_elem, 'binary', &
                      trim(filename_np)//'_grid.dat', this%nvertices, 'binary', &
                      trim(filename_np)//'_points.dat', &
                      this%nelements, trim(filename_np)//'_block.dat'
  close(iinput_xdmf)

732 format(&    
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/&
    '  <Grid GridType="Uniform">',/&
    '    <Topology TopologyType="', A, '" NumberOfElements="',i10,'">',/&
    '      <DataItem Dimensions="',i10, i3, '" NumberType="Int" Format="',A,'">',/&
    '        ',A,/&
    '      </DataItem>',/&
    '    </Topology>',/&
    '    <Geometry GeometryType="XYZ">',/&
    '      <DataItem Dimensions="',i10,' 3" NumberType="Float" Format="',A,'">',/&
    '        ',A/&
    '      </DataItem>',/&
    '    </Geometry>',/&
    '    <Attribute Name="blockid" AttributeType="Scalar" Center="Cell">',/&
    '      <DataItem Dimensions="',i10,'" NumberType="Int" Format="binary">',/&
    '        ',A,/&
    '      </DataItem>',/&
    '    </Attribute>',/&
    '  </Grid>',/&
    '</Domain>',/&
    '</Xdmf>')

  ! VERTEX data
  open(newunit=iinput_heavy_data, file=trim(filename)//'_points.dat', access='stream', &
      status='replace', form='unformatted', convert='little_endian')
  write(iinput_heavy_data) real(this%vertices, kind=sp)
  close(iinput_heavy_data)

  ! CONNECTIVITY data
  open(newunit=iinput_heavy_data, file=trim(filename)//'_grid.dat', access='stream', &
      status='replace', form='unformatted', convert='little_endian')
  write(iinput_heavy_data) this%connectivity
  close(iinput_heavy_data)

  ! BLOCK data
  open(newunit=iinput_heavy_data, file=trim(filename)//'_block.dat', access='stream', &
      status='replace', form='unformatted', convert='little_endian')
  write(iinput_heavy_data) this%block_id
  close(iinput_heavy_data)
  
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine free_node_and_cell_data(this)
  class(inversion_mesh_data_type)   :: this

  if (allocated(this%datat_node))            deallocate(this%datat_node)
  if (allocated(this%datat_cell))            deallocate(this%datat_cell)
  if (allocated(this%data_group_names_node)) deallocate(this%data_group_names_node)
  if (allocated(this%group_id_node))         deallocate(this%group_id_node)
  if (allocated(this%data_group_names_cell)) deallocate(this%data_group_names_cell)
  if (allocated(this%group_id_cell))         deallocate(this%group_id_cell)

end subroutine free_node_and_cell_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine freeme(this)
  class(inversion_mesh_type)        :: this

  if (allocated(this%vertices)) deallocate(this%vertices)
  if (allocated(this%connectivity)) deallocate(this%connectivity)

  this%nvertices = -1
  this%nelements = -1
  this%nvertices_per_elem = -1
  this%nbasisfuncs_per_elem = -1
  this%element_type = ''

  if (allocated(this%v_2d)) deallocate(this%v_2d)
  if (allocated(this%p_2d)) deallocate(this%p_2d)
  if (allocated(this%abinv)) deallocate(this%abinv)
  if (allocated(this%block_id)) deallocate(this%block_id)

  select type (this)
  type is (inversion_mesh_data_type)
     if (allocated(this%datat_node))            deallocate(this%datat_node)
     if (allocated(this%datat_cell))            deallocate(this%datat_cell)
     if (allocated(this%data_group_names_node)) deallocate(this%data_group_names_node)
     if (allocated(this%group_id_node))         deallocate(this%group_id_node)
     if (allocated(this%data_group_names_cell)) deallocate(this%data_group_names_cell)
     if (allocated(this%group_id_cell))         deallocate(this%group_id_cell)
  end select

  this%initialized = .false.

end subroutine freeme
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_node_data(this, dt, starttime)
# if defined(__INTEL_COMPILER)
use ifport, only: getpid ! For ifort, this module needs to be loaded, 
                         ! for gfortran getpid is a GNU extension
# endif  
  use nc_routines, only                   : nc_create_file, nc_create_group
  class(inversion_mesh_data_type)        :: this
  real(kind=dp), intent(in), optional    :: dt, starttime

  character(len=80)                      :: filename

  ! Check whether the NetCDF output file has already been created and opened
  if (this%ncid.eq.-1) then ! If it has not, do so
    this%nvar_node = 0 
    write(this%filename_tmp_out, "('netcdf_out_', I6.6, '.tmp')") getpid()
    call nc_create_file(filename = this%filename_tmp_out,  &
                        ncid = this%ncid, overwrite=.true.) 
  end if

  call nc_create_group(ncid = this%ncid, group_name = 'node_data', &
                       ncid_group = this%ncid_node, overwrite = .true.)
  
  if (present(dt)) then
    this%dt = dt
  else
    this%dt = 1
  end if
  if (present(starttime)) then
    this%starttime = starttime
  else
    this%starttime = 0
  end if

end subroutine init_node_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_cell_data(this, dt, starttime)
# if defined(__INTEL_COMPILER)
use ifport, only: getpid ! For ifort, this module needs to be loaded, 
                         ! for gfortran getpid is a GNU extension
# endif  
  use nc_routines, only                   : nc_create_file, nc_create_group
  class(inversion_mesh_data_type)        :: this
  real(kind=dp), intent(in), optional    :: dt, starttime

  character(len=80)                      :: filename

  ! Check whether the NetCDF output file has already been created and opened
  if (this%ncid.eq.-1) then ! If it has not, do so
    this%nvar_cell = 0 
    write(this%filename_tmp_out, "('netcdf_out_', I6.6, '.tmp')") getpid()
    call nc_create_file(filename = this%filename_tmp_out,  &
                        ncid = this%ncid, overwrite=.true.) 
  end if

  call nc_create_group(ncid = this%ncid, group_name = 'cell_data', &
                       ncid_group = this%ncid_cell, overwrite = .true.)
  
  if (present(dt)) then
    this%dt = dt
  else
    this%dt = 1
  end if
  if (present(starttime)) then
    this%starttime = starttime
  else
    this%starttime = 0
  end if

end subroutine init_cell_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine add_node_variable(this, var_name, nentries, entry_names, istime)
  use nc_routines, only                      : nc_create_var_by_name
  class(inversion_mesh_data_type)           :: this
  character(len=*), intent(in)              :: var_name
  integer, intent(in)                       :: nentries
  character(len=*), intent(in), optional    :: entry_names(:)
  logical, intent(in), optional             :: istime

  integer                                   :: ivar, ientry
  logical                                   :: istime_loc = .false.
  character(len=80), allocatable            :: entry_names_loc(:)

  if (this%ncid_node==-1) then
    print *, 'ERROR: This mesh has not been initialized for node data yet'
  end if

  if (present(entry_names).and.present(istime)) then
    if (istime) then
      print *, 'ERROR: Cell variable cannot have entry names, if it is a time variable'
      stop
    end if
  end if

  if (present(entry_names)) then
    entry_names_loc = entry_names
  else
    allocate(entry_names_loc(nentries))
    entry_names_loc(:) = var_name
  end if

  if (present(istime)) istime_loc = istime

  call nc_create_var_by_name(ncid    = this%ncid_node,                  &
                             varname = trim(var_name),                  &
                             sizes   = [this%nvertices, nentries],      &
                             dimension_names = ['Vertices', '        '])

  if (istime_loc) this%ntimes = max(this%ntimes, nentries)

  call append_variable(this%variable_node, var_name, nentries, istime_loc, entry_names_loc)

  this%nvar_node = this%nvar_node + 1

end subroutine add_node_variable
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine add_cell_variable(this, var_name, nentries, entry_names, istime)
  use nc_routines, only                      : nc_create_var_by_name
  class(inversion_mesh_data_type)           :: this
  character(len=*), intent(in)              :: var_name
  integer, intent(in)                       :: nentries
  character(len=*), intent(in), optional    :: entry_names(:)
  logical, intent(in), optional             :: istime

  integer                                   :: ivar
  logical                                   :: istime_loc = .false.
  character(len=80), allocatable            :: entry_names_loc(:)

  if (this%ncid_cell==-1) then
    print *, 'ERROR: This mesh has not been initialized for cell data yet'
  end if

  if (present(entry_names).and.present(istime)) then
    if (istime) then
      print *, 'ERROR: Cell variable cannot have entry names, if it is a time variable'
      stop
    end if
  end if

  if (present(entry_names)) then
    entry_names_loc = entry_names
  else
    allocate(entry_names_loc(nentries))
    entry_names_loc(:) = trim(var_name)
  end if

  if (present(istime)) istime_loc = istime

  call nc_create_var_by_name(ncid    = this%ncid_cell,                  &
                             varname = trim(var_name),                  &
                             sizes   = [this%nelements, nentries],      &
                             dimension_names = ['Elements', '        '])

  if (istime_loc) this%ntimes = max(this%ntimes, nentries)

  call append_variable(this%variable_cell, var_name, nentries, istime_loc, entry_names_loc)

  this%nvar_cell = this%nvar_cell + 1

end subroutine add_cell_variable
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Increase size of inversion_mesh_variable_type by one. Fortran style
subroutine append_variable(variable, var_name, nentries, istime, entry_names)
  type(inversion_mesh_variable_type), intent(inout), allocatable  :: variable(:)
  character(len=*),                   intent(in)                  :: var_name
  integer,                            intent(in)                  :: nentries
  character(len=*),                   intent(in)                  :: entry_names(:)
  logical,                            intent(in)                  :: istime

  type (inversion_mesh_variable_type), allocatable  :: variable_temp(:)
  integer                                           :: ivar, nvar, ientry

  if (nentries.ne.size(entry_names)) then
    print *, 'ERROR in append_variable. Wrong number of entry names'
    print *, '      nentries:          ', nentries
    print *, '      size(entry_names): ', size(entry_names)
    print *, '      entry_names:       '
    do ientry = 1, size(entry_names)
      print *, '      ', ientry, trim(entry_names(ientry))
    end do
    call pabort()
  end if

  if (.not.allocated(variable)) then
    ! If it did not exist yet, create only
    allocate(variable(1))
    nvar = 1
  else
    nvar = size(variable)
    ! Copy 'variable' to 'variable_temp'
    allocate(variable_temp(nvar))
    do ivar = 1, nvar
      allocate(variable_temp(ivar)%entry_names(variable(ivar)%nentries)) 
      variable_temp(ivar)%var_name    = variable(ivar)%var_name
      variable_temp(ivar)%nentries    = variable(ivar)%nentries
      variable_temp(ivar)%istime      = variable(ivar)%istime
      variable_temp(ivar)%entry_names = variable(ivar)%entry_names
    end do

    ! Increase size of 'variable' by one
    deallocate(variable)
    allocate(variable(nvar+1))

    ! Copy 'variable_temp' to 'variable'
    do ivar = 1, nvar
      allocate(variable(ivar)%entry_names(nentries)) 
      variable(ivar)%var_name    = variable_temp(ivar)%var_name
      variable(ivar)%nentries    = variable_temp(ivar)%nentries
      variable(ivar)%istime      = variable_temp(ivar)%istime
      variable(ivar)%entry_names = variable_temp(ivar)%entry_names
    end do
    deallocate(variable_temp)

    nvar = nvar + 1
  end if

  ! Assign values to new entry of variable
  variable(nvar)%var_name       = var_name
  variable(nvar)%nentries       = nentries
  variable(nvar)%istime         = istime
  allocate(variable(nvar)%entry_names(nentries)) 
  variable(nvar)%entry_names(:) = entry_names

end subroutine append_variable
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine add_node_data(this, var_name, values, ielement, ientry)
  use nc_routines, only                      : nc_putvar_by_name
  class(inversion_mesh_data_type)           :: this
  character(len=*), intent(in)              :: var_name
  real(kind=sp)                             :: values(:,:)
  integer, intent(in), optional             :: ielement(2), ientry(2)

  integer                                   :: ielement_loc(2), ientry_loc(2), start(2), count(2)
  integer                                   :: ivar

  if (this%ncid_node==-1) then
     write(*,*) 'ERROR: trying to write node data without initialization!'
     call pabort 
  end if

  ! Find variable
  do ivar = 1, this%nvar_node
    if (trim(var_name).eq.trim(this%variable_node(ivar)%var_name)) exit
  end do
  
  ! If variable not found
  if (ivar == this%nvar_node + 1) then
    write(*,*) 'ERROR: Could not find variable in node data list'
    write(*,*) '       var_name: ', trim(var_name)
    write(*,*) '       Available variables for node data:'
    do ivar = 1, this%nvar_node
      write(*,*) '       ', trim(this%variable_node(ivar)%var_name)
    end do
    call pabort ()
  end if

  ! If ielement is not present, assume that plot ranges over all elements
  ! (Full snapshot)
  if (present(ielement)) then
    ielement_loc = ielement
  else
    ielement_loc(1) = 1
    ielement_loc(2) = this%nvertices
  end if

  ! If ientry is not present, assume that plot ranges over all entries for 
  ! this variable (Full time series or all kernels)
  if (present(ientry)) then
    ientry_loc = ientry
  else
    ientry_loc(1) = 1
    ientry_loc(2) = this%variable_node(ivar)%nentries
  end if

  start = [ielement_loc(1),                   ientry_loc(1)]
  count = [ielement_loc(2)-ielement_loc(1)+1, ientry_loc(2)-ientry_loc(1)+1]

  if (size(values,1).ne.count(1) .or. size(values,2).ne.count(2)) then
    write(*,*) 'ERROR: Wrong dimension of input to add_node_data!'
    write(*,*) '       var_name:     ', trim(var_name)
    write(*,*) '       nvertices:    ', this%nvertices
    write(*,*) '       nentries:     ', this%variable_node(ivar)%nentries
    write(*,*) '       size(values): ', size(values,1), size(values,2) 
    write(*,*) '       count:        ', count
    call pabort()
  end if
    
  call nc_putvar_by_name(ncid    = this%ncid_node,       &
                         varname = var_name,             &
                         values  = values,               &
                         start   = start,                &
                         count   = count)

end subroutine add_node_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine add_cell_data(this, var_name, values, ielement, ientry)
  use nc_routines, only                      : nc_putvar_by_name
  class(inversion_mesh_data_type)           :: this
  character(len=*), intent(in)              :: var_name
  real(kind=sp), optional, intent(in)       :: values(:,:)
  integer, intent(in), optional             :: ielement(2), ientry(2)

  integer                                   :: ielement_loc(2), ientry_loc(2), start(2), count(2)
  integer                                   :: ivar

  if (this%ncid_cell==-1) then
     write(*,*) 'ERROR: trying to write cell data without initialization!'
     call pabort 
  end if

  ! Find variable
  do ivar = 1, this%nvar_cell
    if (trim(var_name).eq.trim(this%variable_cell(ivar)%var_name)) exit
  end do
  
  ! If variable not found
  if (ivar == this%nvar_cell + 1) then
    write(*,*) 'ERROR: Could not find variable in cell data list'
    write(*,*) '       var_name: ', trim(var_name)
    write(*,*) '       Available variables for cell data:'
    do ivar = 1, this%nvar_cell
      write(*,*) '       ', trim(this%variable_cell(ivar)%var_name)
    end do
    call pabort ()
  end if

  ! If ielement is not present, assume that plot ranges over all elements
  ! (Full snapshot)
  if (present(ielement)) then
    ielement_loc = ielement
  else
    ielement_loc(1) = 1
    ielement_loc(2) = this%nelements
  end if

  ! If ientry is not present, assume that plot ranges over all entries for 
  ! this variable (Full time series or all kernels)
  if (present(ientry)) then
    ientry_loc = ientry
  else
    ientry_loc(1) = 1
    ientry_loc(2) = this%variable_cell(ivar)%nentries
  end if

  start = [ielement_loc(1),                   ientry_loc(1)]
  count = [ielement_loc(2)-ielement_loc(1)+1, ientry_loc(2)-ientry_loc(1)+1]

  if (size(values,1).ne.count(1) .or. size(values,2).ne.count(2)) then
    write(*,*) 'ERROR: Wrong dimension of input to add_cell_data!'
    write(*,*) '       var_name:     ', trim(var_name)
    write(*,*) '       nelements:    ', this%nelements
    write(*,*) '       nentries:     ', this%variable_cell(ivar)%nentries
    write(*,*) '       size(values): ', size(values,1), size(values,2) 
    write(*,*) '       count:        ', count
    call pabort()
  end if
    
  call nc_putvar_by_name(ncid    = this%ncid_cell,       &
                         varname = var_name,             &
                         values  = values,               &
                         start   = start,                &
                         count   = count)

end subroutine add_cell_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_cell_data_csr( this, data_kernel, nkernels, threshold, filename)
  class(inversion_mesh_data_type)           :: this
  real(kind=sp), intent(in)                 :: data_kernel(:,:)
  character(len=*), intent(in)              :: filename
  integer, intent(in)                       :: nkernels
  real(kind=dp)                             :: threshold
  integer                                   :: iinput_ind
  integer                                   :: iinput_val
  integer                                   :: iinput_pnt
  integer                                   :: ielement
  integer                                   :: ikernel
  integer                                   :: ientry

  open(newunit=iinput_val,file=trim(filename)//'.val',&
       access='direct',recl=4,form='unformatted')
  open(newunit=iinput_ind,file=trim(filename)//'.ind',&
       access='direct',recl=4,form='unformatted')
  open(newunit=iinput_pnt,file=trim(filename)//'.pnt')

  ientry = 0
  do ikernel=1,nkernels
     do ielement=1,this%nelements
        if ( abs(data_kernel(ielement,ikernel)).gt.threshold ) then
           ientry = ientry + 1
           write(iinput_val,rec=ientry) data_kernel(ielement,ikernel)
           write(iinput_ind,rec=ientry) ielement
        end if
     end do
     write(iinput_pnt,*) ientry
  end do

end subroutine dump_cell_data_csr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_node_data_csr( this, data_kernel, nkernels, threshold, filename)
  class(inversion_mesh_data_type)           :: this
  real(kind=sp), intent(in)                 :: data_kernel(:,:)
  character(len=*), intent(in)              :: filename
  integer, intent(in)                       :: nkernels
  real(kind=dp)                             :: threshold
  integer                                   :: iinput_ind
  integer                                   :: iinput_val
  integer                                   :: iinput_pnt
  integer                                   :: ikernel
  integer                                   :: ivertex
  integer                                   :: ientry

  open(newunit=iinput_val,file=trim(filename)//'.val',&
       access='direct',recl=4,form='unformatted')
  open(newunit=iinput_ind,file=trim(filename)//'.ind',&
       access='direct',recl=4,form='unformatted')
  open(newunit=iinput_pnt,file=trim(filename)//'.pnt')

  ientry = 0
  do ikernel=1,nkernels
     do ivertex=1,this%nvertices
        if ( abs(data_kernel(ivertex,ikernel)).gt.threshold ) then
           ientry = ientry + 1
           write(iinput_val,rec=ientry) data_kernel(ivertex,ikernel)
           write(iinput_ind,rec=ientry) ivertex
        end if
     end do
     write(iinput_pnt) ientry
  end do


end subroutine dump_node_data_csr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_cell_data_ascii( this, data_kernel, nkernels, threshold, filename)
  class(inversion_mesh_data_type)           :: this
  real(kind=sp), intent(in)                 :: data_kernel(:,:)
  character(len=*), intent(in)              :: filename
  integer, intent(in)                       :: nkernels
  real(kind=dp)                             :: threshold
  integer                                   :: iinput_ind
  integer                                   :: iinput_val
  integer                                   :: iinput_pnt
  integer                                   :: ielement
  integer                                   :: ikernel
  integer                                   :: ientry

  open(newunit=iinput_val,file=trim(filename)//'.val')
  open(newunit=iinput_ind,file=trim(filename)//'.ind')
  open(newunit=iinput_pnt,file=trim(filename)//'.pnt')

  ientry = 0
  do ikernel=1,nkernels
     do ielement=1,this%nelements
        if ( abs(data_kernel(ielement,ikernel)).gt.threshold ) then
           ientry = ientry + 1
           write(iinput_val, *) data_kernel(ielement,ikernel)
           write(iinput_ind, *) ielement
        end if
     end do
     write(iinput_pnt,*) ientry
  end do

end subroutine dump_cell_data_ascii
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_node_data_ascii( this, data_kernel, nkernels, threshold, filename)
  class(inversion_mesh_data_type)           :: this
  real(kind=sp), intent(in)                 :: data_kernel(:,:)
  character(len=*), intent(in)              :: filename
  integer, intent(in)                       :: nkernels
  real(kind=dp)                             :: threshold
  integer                                   :: iinput_ind
  integer                                   :: iinput_val
  integer                                   :: iinput_pnt
  integer                                   :: ikernel
  integer                                   :: ivertex
  integer                                   :: ientry

  open(newunit=iinput_val,file=trim(filename)//'.val')
  open(newunit=iinput_ind,file=trim(filename)//'.ind')
  open(newunit=iinput_pnt,file=trim(filename)//'.pnt')

  ientry = 0
  do ikernel=1,nkernels
     do ivertex=1,this%nvertices
        if ( abs(data_kernel(ivertex,ikernel)).gt.threshold ) then
           ientry = ientry + 1
           write(iinput_val, '(" ", E15.8)', advance='no') data_kernel(ivertex,ikernel)
           write(iinput_ind, '(" ", I15)', advance='no') ivertex
        end if
     end do
     write(iinput_pnt, *) ientry
     write(iinput_val, '(A)', advance='yes') 
     write(iinput_ind, '(A)', advance='yes') 
  end do


end subroutine dump_node_data_ascii
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_data_xdmf(this, filename)
  use nc_routines, only              : nc_close_file, nc_putvar_by_name
  class(inversion_mesh_data_type)   :: this
  character(len=*), intent(in)      :: filename

  integer                           :: iinput_xdmf, iinput_heavy_data
  integer                           :: itime, ivar, ientry, exitstat = 1
  character(len=16)                 :: xdmf_elem_type
  character(len=512)                :: filename_np, filename_nc
  character(len=512)                :: cmdmsg, sys_cmd

  if (.not. this%initialized) then
     write(*,*) 'ERROR: trying to dump a non initialized mesh'
     call pabort 
  end if
  
  if (.not. allocated(this%datat_cell) .and.this%ncid==-1) then
     write(*,*) 'ERROR: no data to dump available'
     call pabort 
  end if

  ! relative filename for xdmf content
  filename_np = trim(filename(index(filename, '/', back=.true.)+1:))
  ! relative filename of NetCDF file
  filename_nc = trim(filename_np)//'.nc'


  select case(this%element_type)
  case('tri')
     xdmf_elem_type = 'Triangle'
  case('tet')
     xdmf_elem_type = 'Tetrahedron'
  case('quad')
     xdmf_elem_type = 'Quadrilateral'
  case('hex')
     xdmf_elem_type = 'Hexahedron'
  case('vox')
     xdmf_elem_type = 'Hexahedron'
  case default
     write(6,*) 'ERROR: xmdf dumping for element type ', this%element_type, &
                ' not implemented'
     call pabort
  end select

  ! XML header
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  ! start xdmf file, write header
  write(iinput_xdmf, 733) this%nelements, this%nvertices_per_elem, 'hdf', &
                          trim(filename_nc)//':/grid', &
                          this%nvertices, 'hdf', trim(filename_nc)//':/points', &
                          this%nelements, 'hdf', trim(filename_nc)//':/block'

  ! First part: Variables, which are not time variables and therefore 
  !             have only one time step

  ! create new snapshot in the temporal collection
  write(iinput_xdmf, 7341) 'grid', this%starttime, &
                           trim(xdmf_elem_type), this%nelements, &
                           "'", "'", "'", "'", "'", "'"

  do ivar = 1, this%nvar_node
    do ientry = 1, this%variable_node(ivar)%nentries
      ! Time variables will be written into the next time steps in the the 
      ! second part. For all other variables, all entries are written into
      ! the first time step
      if (this%variable_node(ivar)%istime.and.ientry>1) cycle
      write(iinput_xdmf, 73421) trim(this%variable_node(ivar)%entry_names(ientry)), &
                               this%nvertices, ientry-1, this%nvertices,            &
                               this%variable_node(ivar)%nentries, this%nvertices,   &
                               trim(filename_nc),                                   &
                               trim(this%variable_node(ivar)%var_name)
    end do
  end do 

  do ivar = 1, this%nvar_cell
    do ientry = 1, this%variable_cell(ivar)%nentries
      ! Time variables will be written into the next time steps in the the 
      ! second part. For all other variables, all entries are written into
      ! the first time step
      if (this%variable_cell(ivar)%istime.and.ientry>1) cycle
      write(iinput_xdmf, 73422) trim(this%variable_cell(ivar)%entry_names(ientry)), &
                               this%nelements, ientry-1, this%nelements,            &
                               this%variable_cell(ivar)%nentries, this%nelements,   &
                               trim(filename_nc),                                   &
                               trim(this%variable_cell(ivar)%var_name)
    end do
  end do
  ! Finish first snapshot
  write(iinput_xdmf, 7343)


  ! Second part: All other time steps, but only for variables that have
  !              data for this time step.
  do itime = 2, this%ntimes 

    ! create new snapshot in the temporal collection
    write(iinput_xdmf, 7341) 'grid', this%dt * itime + this%starttime, &
                             trim(xdmf_elem_type), this%nelements, &
                             "'", "'", "'", "'", "'", "'"

    do ivar = 1, this%nvar_node
      ! Only write if this is a time variable and we still have entries 
      ! for this time step. Unless, of course, this is the first time 
      ! step, where every entry has to be written.
      if (itime.le.this%variable_node(ivar)%nentries .and. &
          this%variable_node(ivar)%istime) then
        write(iinput_xdmf, 73421) trim(this%variable_node(ivar)%entry_names(ientry)), &
                                 this%nvertices, itime-1, this%nvertices,             &
                                 this%variable_node(ivar)%nentries, this%nvertices,   &
                                 trim(filename_nc),                                   &
                                 trim(this%variable_node(ivar)%var_name)
      end if
    end do 

    do ivar = 1, this%nvar_cell
      if (itime.le.this%variable_cell(ivar)%nentries .and. &
          this%variable_cell(ivar)%istime) then
        write(iinput_xdmf, 73422) trim(this%variable_cell(ivar)%entry_names(ientry)), &
                                 this%nelements, itime-1, this%nelements,             &
                                 this%variable_cell(ivar)%nentries, this%nvertices,   &
                                 trim(filename_nc),                                   &
                                 trim(this%variable_cell(ivar)%var_name)
      end if
    end do
    write(iinput_xdmf, 7343)

  enddo

  ! finish xdmf file
  write(iinput_xdmf, 736)
  close(iinput_xdmf)

733 format(&    
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/,/&
    '<DataItem Name="grid" Dimensions="',i10, i3, '" NumberType="Int" Format="',A,'">',/&
    '  ', A,/&
    '</DataItem>',/,/&
    '<DataItem Name="points" Dimensions="',i10,' 3" NumberType="Float" Format="',A,'">',/&
    '  ', A,/&
    '</DataItem>',/,/&
    '<DataItem Name="block" Dimensions="',i10'" NumberType="Int" Format="',A,'">',/&
    '  ', A,/&
    '</DataItem>',/,/&
    '<Grid Name="CellsTime" GridType="Collection" CollectionType="Temporal">',/)

7341 format(&    
    '    <Grid Name="', A,'" GridType="Uniform">',/&
    '        <Time Value="',F8.2,'" />',/&
    '        <Topology TopologyType="', A, '" NumberOfElements="',i10,'">',/&
    '            <DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'grid', A,']" />',/&
    '        </Topology>',/&
    '        <Geometry GeometryType="XYZ">',/&
    '            <DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '        </Geometry>',/&
    '        <Attribute Name="blockid" AttributeType="Scalar" Center="Cell">',/&
    '            <DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'block', A,']" />',/&
    '        </Attribute>')

73421 format(&    
    '        <Attribute Name="', A,'" AttributeType="Scalar" Center="Node">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/node_data/', A ,/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>')

73422 format(&    
    '        <Attribute Name="', A,'" AttributeType="Scalar" Center="Cell">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/cell_data/', A ,/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>')

7343 format(&    
    '    </Grid>',/)

736 format(&    
    '</Grid>',/&
    '</Domain>',/&
    '</Xdmf>')

  ! VERTEX data
  call nc_putvar_by_name(ncid    = this%ncid,                   &
                         varname = 'points',                    &
                         values  = real(this%vertices, kind=sp))

  ! CONNECTIVITY data
  call nc_putvar_by_name(ncid    = this%ncid,        &
                         varname = 'grid',           &
                         values  = this%connectivity)

  ! BLOCK data
  call nc_putvar_by_name(ncid    = this%ncid,    &
                         varname = 'block',      &
                         values  = this%block_id)

  ! Finalize and close NetCDF file
  call nc_close_file(this%ncid)

  ! Move temporary NetCDF file to the requested location
  sys_cmd = 'mv '//this%filename_tmp_out//' '//trim(filename)//'.nc'
  call execute_command_line(command = sys_cmd, wait = .true., exitstat = exitstat, &
                            cmdmsg = cmdmsg)
  if (exitstat.ne.0) then
    write(*,*) 'ERROR: Renaming of the temporary output file failed with status', exitstat
    write(*,*) '       and message:'
    write(*,*) cmdmsg
    call pabort()
  end if

end subroutine dump_data_xdmf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
# if defined(__INTEL_COMPILER)
! Apparently, ifort does not support execute_command_line yet
subroutine execute_command_line(command, wait, exitstat, cmdmsg)
  use ifport
  character(len=*), intent(in)   :: command
  logical, intent(in)            :: wait
  integer, intent(out)           :: exitstat
  character(len=80), intent(out) :: cmdmsg

  exitstat = system(command)
  cmdmsg = 'Ifort-specific replacement for execute_command_line, with meaningless output msg'
end subroutine execute_command_line
# endif
!-----------------------------------------------------------------------------------------
end module
!=========================================================================================
