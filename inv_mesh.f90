!=========================================================================================
module inversion_mesh

  use global_parameters, only: sp, dp
  use tetrahedra,        only: get_volume_tet, get_volume_poly, &
                               generate_random_points_tet, generate_random_points_poly
  implicit none
  private
  public :: inversion_mesh_type
  public :: inversion_mesh_data_type
  public :: plane_exp_pro2

  type :: inversion_mesh_type
     private
     integer, public                    :: nvertices_per_elem
     integer                            :: nelements, nvertices
     integer, allocatable               :: connectivity(:,:)
     real(kind=dp), allocatable         :: vertices(:,:)
     real(kind=dp), allocatable         :: v_2d(:,:,:), p_2d(:,:,:)
     character(len=4)                   :: element_type
     logical                            :: initialized = .false.
     contains
     procedure, pass :: get_nelements
     procedure, pass :: get_element
     procedure, pass :: get_nvertices
     procedure, pass :: get_vertices
     procedure, pass :: get_valence
     procedure, pass :: get_connected_elements
     procedure, pass :: get_connectivity
     procedure, pass :: read_tet_mesh
     procedure, pass :: read_abaqus_mesh
     procedure, pass :: dump_mesh_xdmf
     procedure, pass :: get_volume
     procedure, pass :: generate_random_points
     procedure, pass :: make_2d_vectors
     procedure, pass :: freeme
  end type

  type, extends(inversion_mesh_type)    :: inversion_mesh_data_type
     private
     integer                            :: ntimes
     real(kind=sp), allocatable         :: datat(:,:)
     integer, allocatable               :: group_id(:)
     character(len=16), allocatable     :: data_group_names(:)
     integer                            :: ngroups
     contains
     procedure, pass :: get_ntimes
     procedure, pass :: init_data
     procedure, pass :: set_data_snap
     procedure, pass :: dump_mesh_data_xdmf
  end type

contains

!-----------------------------------------------------------------------------------------
integer function get_nelements(this)
  class(inversion_mesh_type)        :: this
  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'
  get_nelements = this%nelements
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_element(this, ielement)
  class(inversion_mesh_type)        :: this
  real(kind=dp)                     :: get_element(3,this%nvertices_per_elem)
  integer, intent(in)               :: ielement
  integer                           :: ivert

  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'

  do ivert=1, this%nvertices_per_elem
     get_element(:,ivert) = this%vertices(:, this%connectivity(ivert,ielement)+1) 
  enddo
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_nvertices(this)
  class(inversion_mesh_type)        :: this
  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'
  get_nvertices = this%nvertices
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_vertices(this)
  class(inversion_mesh_type)        :: this
  real(kind=sp)                     :: get_vertices(3,this%nvertices)

  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'
  get_vertices = this%vertices
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_valence(this, ivert)
  class(inversion_mesh_type)        :: this
  integer                           :: get_valence
  integer, intent(in)               :: ivert
  
  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'

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
  
  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'

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
  integer                           :: get_connectivity(this%nvertices_per_elem, this%nelements)
  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'
  get_connectivity = this%connectivity + 1
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_ntimes(this)
  class(inversion_mesh_data_type)        :: this
  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh data type that is not initialized'
  get_ntimes = this%ntimes
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_volume(this, ielement)
  class(inversion_mesh_type)        :: this
  integer, intent(in)               :: ielement
  real(kind=dp)                     :: get_volume

  select case(this%element_type)
  case('tet')
     get_volume = get_volume_tet(this%get_element(ielement))
  case('quad')
     get_volume = get_volume_poly(4, this%get_element(ielement))
  case('tri')
     get_volume = get_volume_poly(3, this%get_element(ielement))
  case('hex')
     !get_volume = get_volume_hex(this%get_element(ielement))
  end select

end function get_volume
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function generate_random_points(this, ielement, npoints) result(points)
  class(inversion_mesh_type)     :: this
  integer, intent(in)            :: ielement, npoints
  real(kind=dp)                  :: points(3, npoints)
  real(kind=dp)                  :: points2d(2, npoints)
  integer                        :: ipoint

  select case(this%element_type)
  case('tet')
     points = generate_random_points_tet(this%get_element(ielement), npoints)
  case('quad')
     !generate_random_points = generate_random_points_poly(4, this%get_element(ielement), &
     !                                                     npoints)
     points2d = generate_random_points_poly(4, this%p_2d(:,:,ielement), npoints)
     do ipoint = 1, npoints
        points(:,ipoint) =   this%v_2d(:,0,ielement)                         &
                           + this%v_2d(:,1,ielement) * points2d(1,ipoint)    &
                           + this%v_2d(:,2,ielement) * points2d(2,ipoint)
     end do
  case('tri')
     points2d = generate_random_points_poly(3, this%p_2d(:,:,ielement), npoints)
     do ipoint = 1, npoints
        points(:,ipoint) =   this%v_2d(:,0,ielement)                         &
                           + this%v_2d(:,1,ielement) * points2d(1,ipoint)    &
                           + this%v_2d(:,2,ielement) * points2d(2,ipoint)
     end do
  case('hex')
     !generate_random_points = generate_random_points_hex( this%get_element(ielement), &
     !                                                     ielement, npoints)
  end select

end function generate_random_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_tet_mesh(this, filename_vertices, filename_connectivity)
  class(inversion_mesh_type)        :: this
  character(len=*), intent(in)      :: filename_vertices
  character(len=*), intent(in)      :: filename_connectivity
  integer                           :: iinput_vertices, iinput_connectivity
  integer                           :: i, ierr

  this%nvertices_per_elem = 4
  this%element_type = 'tet'
  
  ! read vertices
  open(newunit=iinput_vertices, file=trim(filename_vertices), status='old', &
       action='read', iostat=ierr)
  if ( ierr /= 0 ) then
     write(*,*) 'ERROR: Could not open file: ', trim(filename_vertices)
     stop
  endif

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
     write(*,*) 'ERROR: Could not open file: ', trim(filename_connectivity)
     stop
  endif

  read(iinput_connectivity,*) this%nelements

  allocate(this%connectivity(this%nvertices_per_elem,this%nelements))
  do i=1, this%nelements
     read(iinput_connectivity,*) this%connectivity(:,i)
  enddo
  close(iinput_connectivity)

  this%initialized = .true.
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_abaqus_mesh(this, filename)
  class(inversion_mesh_type)        :: this
  character(len=*), intent(in)      :: filename
  integer                           :: iinput
  integer                           :: i, ierr, ct
  character(len=128)                :: line
  character(len=16)                 :: elem_type

  ! open file 
  open(newunit=iinput, file=trim(filename), status='old', &
       action='read', iostat=ierr)
  if ( ierr /= 0 ) then
     write(*,*) 'ERROR: Could not open file: ', trim(filename)
     stop
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
  case default
     write(6,*) 'ERROR: reading abaqus file with elementtype ', trim(elem_type), &
                'not yet implemented'
     stop
  end select

  ! scan number of elements
  ct = 0
  do
    read(iinput,*) line
    if (index(line, '**') == 1) exit
    ct = ct + 1
  enddo
  this%nelements = ct

  ! prepare arrays
  allocate(this%vertices(3,this%nvertices))
  allocate(this%connectivity(this%nvertices_per_elem,this%nelements))

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
  do i=1, this%nelements
     read(iinput,*) line, this%connectivity(:,i)
  enddo

  close(iinput)

  ! abaqus starts counting at 1
  this%connectivity(:,:) = this%connectivity(:,:) - 1

  
  if (this%element_type.eq.'tri'.or.this%element_type.eq.'quad') then
     call this%make_2d_vectors
  end if


  this%initialized = .true.
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine make_2d_vectors(this)
  class(inversion_mesh_type)        :: this
  integer                           :: nvec, ielem

  select case(trim(this%element_type))
  case('tri')
     nvec = 2
  case('quad') 
     nvec = 3
  end select
  
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
end subroutine make_2d_vectors
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine plane_exp_pro2 ( p_ref, npoints, p_3d, p_2d, vec )

!*****************************************************************************80
!
!! PLANE_EXP_PRO2 produces 2D coordinates of points that lie in a plane, in 3D.
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
!
  implicit none
  
  integer, parameter         :: dim_num = 3
  integer, intent(in)        :: npoints
  real(kind=dp), intent(in)  :: p_3d(dim_num,npoints)
  real(kind=dp), intent(in)  :: p_ref(dim_num,3)
  real(kind=dp), intent(out) :: p_2d(2, npoints)
  real(kind=dp), intent(out) :: vec(dim_num, 2)
  
  real(kind=dp)              :: dot
  integer                    :: i
  !
  !  Compute the two basis vectors for the affine plane.
  !
  vec(:,1) = p_ref(:,2) - p_ref(:,1)
  
  !call vector_unit_nd ( dim_num, v1 )
  vec(:,1) = vec(:,1) / norm2(vec(:,1))
  
  vec(:,2) = p_ref(:,3) - p_ref(:,1)
  
  dot = dot_product ( vec(:,1), vec(:,2) )
  
  vec(:,2) = vec(:,2) - dot * vec(:,1)
  
  !call vector_unit_nd ( dim_num, v2 )
  vec(:,2) = vec(:,2) / norm2(vec(:,2))
  !
  !  Now decompose each point.
  !
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
  integer                           :: i
  character(len=16)                 :: xdmf_elem_type
  character(len=512)                :: filename_np

  if (.not. this%initialized) &
     stop 'ERROR: trying to dump a non initialized mesh'

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
  case default
     write(6,*) 'ERROR: xmdf dumping for element type ', this%element_type, &
                ' not implemented'
     stop
  end select

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  write(iinput_xdmf, 732) trim(xdmf_elem_type), this%nelements, this%nelements, &
                      this%nvertices_per_elem, 'binary', &
                      trim(filename_np)//'_grid.dat', this%nvertices, 'binary', &
                      trim(filename_np)//'_points.dat'
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
  
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine freeme(this)
  class(inversion_mesh_type)        :: this

  if (allocated(this%vertices)) deallocate(this%vertices)
  if (allocated(this%connectivity)) deallocate(this%connectivity)

  select type (this)
  type is (inversion_mesh_data_type)
     if (allocated(this%datat)) deallocate(this%datat)
  end select

  this%initialized = .false.

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_data(this, ntimes)
  class(inversion_mesh_data_type)   :: this
  integer, intent(in)               :: ntimes

  this%ntimes = ntimes

  allocate(this%datat(this%nvertices, ntimes))
  this%datat = 0

  allocate(this%data_group_names(ntimes))
  this%data_group_names = 'data'

  allocate(this%group_id(ntimes))
  this%group_id = 1

  this%ngroups = 0
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine set_data_snap(this, data_snap, isnap, data_name)
  class(inversion_mesh_data_type)           :: this
  real(kind=sp), intent(in)                 :: data_snap(:)
  integer, intent(in)                       :: isnap
  character(len=*), intent(in), optional    :: data_name
  integer                                   :: i
  logical                                   :: name_exists

  if (.not. allocated(this%datat)) &
     stop 'ERROR: trying to write data without initialization!'

  if (size(data_snap) /= this%nvertices) then
     write(*,*) 'size(data_snap):', size(data_snap), '; this%nvertices', this%nvertices
     stop 'ERROR: wrong dimensions of input data_snap for writing vertex data'
  end if

  this%datat(:,isnap) = data_snap(:)

  if (present(data_name)) then
     name_exists = .false.
     do i=1, this%ngroups
        if (this%data_group_names(i) == data_name) then
           name_exists = .true.
           exit
        endif
     enddo
     if (name_exists) then
        this%group_id(isnap) = i
     else
        this%ngroups = this%ngroups + 1
        this%group_id(isnap) = this%ngroups
        this%data_group_names(this%ngroups) = data_name
     endif
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_mesh_data_xdmf(this, filename)
  class(inversion_mesh_data_type)   :: this
  character(len=*), intent(in)      :: filename
  integer                           :: iinput_xdmf, iinput_heavy_data
  integer                           :: i, igroup, itime, isnap, n
  character(len=16)                 :: xdmf_elem_type
  character(len=512)                :: filename_np

  if (.not. this%initialized) &
     stop 'ERROR: trying to dump a non initialized mesh'
  
  if (.not. allocated(this%datat)) &
     stop 'ERROR: no data to dump available'

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
  case default
     write(6,*) 'ERROR: xmdf dumping for element type ', this%element_type, &
                ' not implemented'
     stop
  end select

  ! XML header
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  ! start xdmf file, write header
  write(iinput_xdmf, 733) this%nelements, this%nvertices_per_elem, 'binary', &
                          trim(filename_np)//'_grid.dat', &
                          this%nvertices, 'binary', trim(filename_np)//'_points.dat'

  i = 1
  itime = 1

  ! loop over all data
  do while(i <= this%ntimes)
     ! create new snapshot in the temporal collection
     write(iinput_xdmf, 7341) 'grid', dble(itime), trim(xdmf_elem_type), this%nelements, &
                              "'", "'", "'", "'"

     igroup = 1
     ! loop over groups, that have more items then itime
     do while(igroup <= this%ngroups)
        if (count(this%group_id == igroup) >= itime) then
           ! find the itime'th item in the group
           n = 0
           do isnap=1, this%ntimes
              if (this%group_id(isnap) == igroup) n = n + 1
              if (n == itime) exit
           enddo

           ! write attribute
           write(iinput_xdmf, 7342) this%data_group_names(igroup), &
                                    this%nvertices, isnap-1, this%nvertices, this%ntimes, &
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
    '</DataItem>',/&
    '<DataItem Name="points" Dimensions="',i10,' 3" NumberType="Float" Format="',A,'">',/&
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
    '        </Geometry>')

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

  ! VERTEX data
  open(newunit=iinput_heavy_data, file=trim(filename)//'_data.dat', access='stream', &
      status='replace', form='unformatted', convert='little_endian')
  do i=1, this%ntimes
     write(iinput_heavy_data) real(this%datat(:,i), kind=sp)
  enddo
  close(iinput_heavy_data)

end subroutine
!-----------------------------------------------------------------------------------------
end module
!=========================================================================================
