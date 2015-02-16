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
  public :: plane_exp_pro2

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
     integer                            :: ncid

     contains
     procedure, pass :: get_ntimes_node
     procedure, pass :: get_ntimes_cell

     ! for xdmf format dumps
     procedure, pass :: init_node_data
     procedure, pass :: init_cell_data
     procedure, pass :: set_node_data_snap
     procedure, pass :: set_cell_data_snap
     procedure, pass :: set_node_data_trace
     procedure, pass :: set_cell_data_trace
     procedure, pass :: dump_node_data_xdmf
     procedure, pass :: dump_cell_data_xdmf

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
  case('quad')
!    get_center = get_center_poly(4, this%get_element(ielement))
  case('tri')
!    get_center = get_center_poly(3, this%get_element(ielement))
  case('vox')
     get_center = get_center_vox(this%get_element(ielement))
  case('hex')
!    get_center = get_center_hex(this%get_element(ielement))
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
  class(inversion_mesh_type)     :: this
  integer, intent(in)            :: ielem, ivertex
  real(kind=dp), intent(in)      :: points(:,:)
  real(kind=dp)                  :: weights(size(points,2))
  real(kind=dp)                  :: dx, dy, dz
  integer                        :: ipoint

  if (.not. this%initialized) then
     write(*,'(A)') 'ERROR: accessing inversion mesh type that is not initialized'
     call pabort 
  end if
  select case(this%element_type) 

  case('tet')
     do ipoint = 1, size(points, 2)
        dx = points(1, ipoint) - this%vertices(1, this%connectivity(1,ielem)+1)
        dy = points(2, ipoint) - this%vertices(2, this%connectivity(1,ielem)+1)
        dz = points(3, ipoint) - this%vertices(3, this%connectivity(1,ielem)+1)

        weights(ipoint) =   this%abinv(ivertex, 1, ielem) * dx &
                          + this%abinv(ivertex, 2, ielem) * dy &
                          + this%abinv(ivertex, 3, ielem) * dz
     end do
  case default
      ! Least sophisticated version possible
      weights = 1
  end select

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
  character(len=*), intent(in)      :: int_type
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
subroutine init_node_data(this, ntimes_node, default_name)
  class(inversion_mesh_data_type)        :: this
  integer, intent(in)                    :: ntimes_node
  character(len=*), intent(in), optional :: default_name

  this%ntimes_node = ntimes_node

  allocate(this%datat_node(this%nvertices, ntimes_node))
  this%datat_node = 0

  allocate(this%data_group_names_node(ntimes_node))

  if (present(default_name)) then
     this%data_group_names_node = trim(default_name)
  else
     this%data_group_names_node = 'node_data'
  endif

  allocate(this%group_id_node(ntimes_node))
  this%group_id_node = 1
  this%ngroups_node = 0

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_cell_data(this, ntimes_cell, default_name, netcdf_in)
  use nc_routines, only                   : nc_open_for_write, nc_createvar_by_name
  class(inversion_mesh_data_type)        :: this
  integer, intent(in)                    :: ntimes_cell
  character(len=*), intent(in), optional :: default_name
  logical, optional                      :: netcdf_in
  character(len=80)                      :: filename

  this%ntimes_cell = ntimes_cell

  allocate(this%data_group_names_cell(ntimes_cell))

  if (present(default_name)) then
     this%data_group_names_cell = trim(default_name)
  else
     this%data_group_names_cell = 'cell_data'
  endif

  if (present(netcdf_in)) then
    if (netcdf_in) then
      filename = 'blubberlutsch.nc'
      call nc_open_for_write(filename = filename, ncid = this%ncid) 
      call nc_createvar_by_name(ncid = this%ncid, &
                                varname = 'data', &
                                sizes   = [this%nelements, ntimes_cell])
      
    end if
  else

    allocate(this%datat_cell(this%nelements, ntimes_cell))
    this%datat_cell = 0

  end if


  allocate(this%group_id_cell(ntimes_cell))
  this%group_id_cell = 1
  this%ngroups_cell = 0

end subroutine
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
subroutine set_node_data_snap(this, data_snap, isnap, data_name)
  class(inversion_mesh_data_type)           :: this
  real(kind=sp), intent(in)                 :: data_snap(:)
  integer, intent(in)                       :: isnap
  character(len=*), intent(in), optional    :: data_name
  integer                                   :: i
  logical                                   :: name_exists

  if (.not. allocated(this%datat_node)) then
     write(*,*) 'ERROR: trying to write node data without initialization!'
     call pabort 
  end if

  if (size(data_snap) /= this%nvertices) then
     write(*,*) 'size(data_snap):', size(data_snap), '; this%nvertices', this%nvertices
     write(*,*) 'data_name:', trim(data_name), '; isnap:', isnap
     write(*,*) 'ERROR: wrong dimensions of input data_snap for writing vertex data'
     call pabort 
  end if

  if (isnap > size(this%datat_node,2)) then
    write(*,*) 'ERROR: XDMF file was initialized for ', size(this%datat_node,2), ' snaps.'
    write(*,*) '       Trying to dump snap with isnap ', isnap
    call pabort()
  end if

  this%datat_node(:,isnap) = data_snap(:)

  if (present(data_name)) then
     name_exists = .false.
     do i=1, this%ngroups_node
        if (this%data_group_names_node(i) == data_name) then
           name_exists = .true.
           exit
        endif
     enddo
     if (name_exists) then
        this%group_id_node(isnap) = i
     else
        this%ngroups_node = this%ngroups_node + 1
        this%group_id_node(isnap) = this%ngroups_node
        this%data_group_names_node(this%ngroups_node) = data_name
     endif
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine set_cell_data_snap(this, data_snap, isnap, data_name)
  class(inversion_mesh_data_type)           :: this
  real(kind=sp), intent(in)                 :: data_snap(:)
  integer, intent(in)                       :: isnap
  character(len=*), intent(in), optional    :: data_name
  integer                                   :: i
  logical                                   :: name_exists

  if (.not. allocated(this%datat_cell)) then
     write(*,*) 'ERROR: trying to write cell data without initialization!'
     call pabort 
  end if

  if (size(data_snap) /= this%nelements) then
     write(*,*) 'size(data_snap):', size(data_snap), '; this%nvertices', this%nvertices
     write(*,*) 'data_name:', trim(data_name), '; isnap:', isnap
     write(*,*) 'ERROR: wrong dimensions of input data_snap for writing cell data'
     call pabort 
  end if

  if (isnap > size(this%datat_cell, 2)) then
    write(*,*) 'ERROR: XDMF file was initialized for ', size(this%datat_cell,2), ' snaps.'
    write(*,*) '       Trying to dump snap with isnap ', isnap
    call pabort()
  end if

  this%datat_cell(:,isnap) = data_snap(:)

  if (present(data_name)) then
     name_exists = .false.
     do i=1, this%ngroups_cell
        if (this%data_group_names_cell(i) == data_name) then
           name_exists = .true.
           exit
        endif
     enddo
     if (name_exists) then
        this%group_id_cell(isnap) = i
     else
        this%ngroups_cell = this%ngroups_cell + 1
        this%group_id_cell(isnap) = this%ngroups_cell
        this%data_group_names_cell(this%ngroups_cell) = data_name
     endif
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine set_node_data_trace(this, data_trace, itrace, data_name)
  class(inversion_mesh_data_type)           :: this
  real(kind=sp), intent(in)                 :: data_trace(:)
  integer, intent(in)                       :: itrace
  character(len=*), intent(in), optional    :: data_name

  if (.not. allocated(this%datat_node)) then
     write(*,*) 'ERROR: trying to write node data without initialization!'
     call pabort 
  end if

  if (size(data_trace) /= this%ntimes_node) then
     write(*,*) 'ERROR: wrong dimensions of input data_trace for writing vertex data'
     write(*,*) 'size(data_trace):', size(data_trace), '; this%ntimes_node', this%ntimes_node
     write(*,*) 'data_name:', trim(data_name), '; itrace:', itrace
     call pabort 
  end if

  if (itrace < 1 .or. itrace > this%nvertices) then
     write(*,*) 'ERROR: index "itrace" out of bounds'
     call pabort 
  end if

  this%datat_node(itrace,:) = data_trace(:)

  if (present(data_name)) then
     this%data_group_names_node = trim(data_name)
  else
     this%data_group_names_node = 'node_data'
  endif

  this%group_id_node = 1
  this%ngroups_node = 1

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine set_cell_data_trace(this, data_trace, itrace, data_name)
  use nc_routines, only                      : nc_putvar_by_name
  class(inversion_mesh_data_type)           :: this
  real(kind=sp), intent(in)                 :: data_trace(:)
  integer, intent(in)                       :: itrace
  character(len=*), intent(in), optional    :: data_name

  if (.not. allocated(this%datat_cell)) then
     write(*,*) 'ERROR: trying to write cell data without initialization!'
     call pabort 
  end if

  if (size(data_trace) /= this%ntimes_cell) then
     write(*,*) 'ERROR: wrong dimensions of input data_trace for writing cell data'
     write(*,*) 'size(data_trace):', size(data_trace), '; this%ntimes_cell', this%ntimes_cell
     write(*,*) 'data_name:', trim(data_name), '; itrace:', itrace
     call pabort 
  end if

  if (itrace < 1 .or. itrace > this%nelements) then
     write(*,*) 'ERROR: index "itrace" out of bounds'
     call pabort 
  end if

  if (this%ncid==-1) then !Opened for binary output
    this%datat_cell(itrace,:) = data_trace(:)
    
  else
    call nc_putvar_by_name(this%ncid,                     &
                           varname = 'data',              &
                           values = data_trace,           &
                           start  = [itrace, 1],          &
                           count  = [1, this%ntimes_cell])

  end if

  if (present(data_name)) then
     this%data_group_names_cell = trim(data_name)
  else
     this%data_group_names_cell = 'cell_data'
  endif

  this%group_id_cell = 1
  this%ngroups_cell = 1

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_cell_data_xdmf(this, filename)
  use nc_routines, only              : nc_close_file
  class(inversion_mesh_data_type)   :: this
  character(len=*), intent(in)      :: filename
  integer                           :: iinput_xdmf, iinput_heavy_data
  integer                           :: i, igroup, itime, isnap, n
  character(len=16)                 :: xdmf_elem_type
  character(len=512)                :: filename_np

  if (.not. this%initialized) then
     write(*,*) 'ERROR: trying to dump a non initialized mesh'
     call pabort 
  end if
  
  if (.not. allocated(this%datat_cell)) then
     write(*,*) 'ERROR: no data to dump available'
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
     write(6,*) 'ERROR: xmdf dumping for element type ', this%element_type, &
                ' not implemented'
     call pabort
  end select

  ! XML header
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  ! start xdmf file, write header
  write(iinput_xdmf, 733) this%nelements, this%nvertices_per_elem, 'binary', &
                          trim(filename_np)//'_grid.dat', &
                          this%nvertices, 'binary', trim(filename_np)//'_points.dat', &
                          this%nelements, 'binary', trim(filename_np)//'_block.dat'

  i = 1
  itime = 1

  ! loop over all data
  do while(i <= this%ntimes_cell)

     ! create new snapshot in the temporal collection
     write(iinput_xdmf, 7341) 'grid', dble(itime), trim(xdmf_elem_type), this%nelements, &
                              "'", "'", "'", "'", "'", "'"

     igroup = 1
     ! loop over groups, that have more items than itime
     do while(igroup <= this%ngroups_cell)

        if (count(this%group_id_cell == igroup) >= itime) then

           ! find the itime'th item in the group
           n = 0
           do isnap=1, this%ntimes_cell
              if (this%group_id_cell(isnap) == igroup) n = n + 1
              if (n == itime) exit
           enddo

           ! write attribute
           ! @TODO: Change numbers to nelements
           if (this%ncid==-1) then
             write(iinput_xdmf, 7342) trim(this%data_group_names_cell(igroup)), &
                                      this%nelements, isnap-1, this%nelements,  &
                                      this%ntimes_cell,                         &
                                      this%nelements, trim(filename_np)//'_data.dat'
           else
             write(iinput_xdmf, 8342) trim(this%data_group_names_cell(igroup)), &
                                      this%nelements, isnap-1, this%nelements,  &
                                      this%ntimes_cell,                         &
                                      this%nelements, 'blubberlutsch.nc'
           end if
           i = i + 1
        endif
        igroup = igroup + 1
     enddo
     ! finish snapshot
     write(iinput_xdmf, 7343)
     itime = itime + 1
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

7342 format(&    
    '        <Attribute Name="', A,'" AttributeType="Scalar" Center="Cell">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary">',/&
    '                   ', A,/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>')

8342 format(&    
    '        <Attribute Name="', A,'" AttributeType="Scalar" Center="Cell">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="hdf">',/&
    '                   ', A, ':/data' ,/&
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

  ! CELL data
  if (this%ncid==-1) then
    open(newunit=iinput_heavy_data, file=trim(filename)//'_data.dat', access='stream', &
        status='replace', form='unformatted', convert='little_endian')
    do i=1, this%ntimes_cell
       write(iinput_heavy_data) real(this%datat_cell(:,i), kind=sp)
    enddo
    close(iinput_heavy_data)
  else
    call nc_close_file(this%ncid)
  end if

end subroutine dump_cell_data_xdmf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_node_data_xdmf(this, filename)
  class(inversion_mesh_data_type)   :: this
  character(len=*), intent(in)      :: filename
  integer                           :: iinput_xdmf, iinput_heavy_data
  integer                           :: i, igroup, itime, isnap, n
  character(len=16)                 :: xdmf_elem_type
  character(len=512)                :: filename_np

  if (.not. this%initialized) then
     write(*,*) 'ERROR: trying to dump a non initialized mesh'
     call pabort 
  end if
  
  if (.not. allocated(this%datat_node)) then
     write(*,*) 'ERROR: no data to dump available'
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
     write(6,*) 'ERROR: xmdf dumping for element type ', this%element_type, &
                ' not implemented'
     call pabort
  end select

  ! XML header
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  ! start xdmf file, write header
  write(iinput_xdmf, 733) this%nelements, this%nvertices_per_elem, 'binary', &
                          trim(filename_np)//'_grid.dat', &
                          this%nvertices, 'binary', trim(filename_np)//'_points.dat', &
                          this%nelements, 'binary', trim(filename_np)//'_block.dat'

  i = 1
  itime = 1

  
  ! loop over all data
  do while(i <= this%ntimes_node)
     ! create new snapshot in the temporal collection
     write(iinput_xdmf, 7341) 'grid', dble(itime), trim(xdmf_elem_type), this%nelements, &
                              "'", "'", "'", "'", "'", "'"

     igroup = 1
     ! loop over groups, that have more items then itime
     do while(igroup <= this%ngroups_node)

        if (count(this%group_id_node == igroup) >= itime) then

           ! find the itime'th item in the group
           n = 0
           do isnap=1, this%ntimes_node
              if (this%group_id_node(isnap) == igroup) n = n + 1
              if (n == itime) exit
           enddo

           ! write attribute
           write(iinput_xdmf, 7342) trim(this%data_group_names_node(igroup)), &
                                    this%nvertices, isnap-1, this%nvertices,  &
                                    this%ntimes_node,                         &
                                    this%nvertices, trim(filename_np)//'_data.dat'
           i = i + 1
        endif
        igroup = igroup + 1
     enddo
     ! finish snapshot
     write(iinput_xdmf, 7343)
     itime = itime + 1
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

7342 format(&    
    '        <Attribute Name="', A,'" AttributeType="Scalar" Center="Node">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary">',/&
    '                   ', A,/&
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

  ! VERTEX data
  open(newunit=iinput_heavy_data, file=trim(filename)//'_data.dat', access='stream', &
      status='replace', form='unformatted', convert='little_endian')
  do i=1, this%ntimes_node
     write(iinput_heavy_data) real(this%datat_node(:,i), kind=sp)
  enddo
  close(iinput_heavy_data)

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
