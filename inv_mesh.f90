!=========================================================================================
module inversion_mesh

  use global_parameters, only: sp, dp
  use type_parameter
  implicit none
  private
  public :: inversion_mesh_type
  public :: inversion_mesh_data_type

  type :: inversion_mesh_type
     private
     integer                            :: nelements, nvertices, nvertices_per_elem
     integer, allocatable               :: connectivity(:,:)
     real(kind=sp), allocatable         :: vertices(:,:)
     character(len=4)                   :: element_type
     logical                            :: initialized = .false.
     contains
     procedure, pass :: get_nelements
     procedure, pass :: get_element
     procedure, pass :: get_nvertices
     procedure, pass :: get_vertices
     procedure, pass :: get_valence
     procedure, pass :: get_connected_elements
     procedure, pass :: read_tet_mesh
     procedure, pass :: read_abaqus_mesh
     procedure, pass :: dump_mesh_xdmf
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
     procedure, pass :: dump_tet_mesh_data_xdmf
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
  real(kind=sp)                     :: get_element(3,this%nvertices_per_elem)
  integer, intent(in)               :: ielement
  integer                           :: ivert

  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'

  do ivert=1, this%nvertices_per_elem
     get_element(:,ivert) = this%vertices(:, this%connectivity(ivert,ielement)) 
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
  real(kind=sp)                     :: get_connectivity(this%nvertices,this%nelements)
  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'
  get_connectivity = this%connectivity
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

  this%initialized = .true.
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_mesh_xdmf(this, filename)
  class(inversion_mesh_type)        :: this
  character(len=*), intent(in)      :: filename
  integer                           :: iinput_xdmf, iinput_heavy_data
  integer                           :: i
  character(len=16)                 :: xdmf_elem_type

  if (.not. this%initialized) &
     stop 'ERROR: trying to dump a non initialized mesh'

  select case(this%element_type)
  case('tri')
     xdmf_elem_type = 'Triangle'
  case('tet')
     xdmf_elem_type = 'Tetrahedron'
  case default
     write(6,*) 'ERROR: xmdf dumping for element type ', this%element_type, &
                ' not implemented'
     stop
  end select

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  write(iinput_xdmf, 732) trim(xdmf_elem_type), this%nelements, this%nelements, &
                      this%nvertices_per_elem, 'binary', &
                      trim(filename)//'_grid.dat', this%nvertices, 'binary', &
                      trim(filename)//'_points.dat'
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
  write(iinput_heavy_data) this%vertices
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
        write(*,*) 'Trying ', trim(this%data_group_names(i))
        if (this%data_group_names(i) == data_name) then
           name_exists = .true.
           write(*,*) 'Name ', trim(data_name), ' exists as no. ', i
           exit
        endif
     enddo
     if (name_exists) then
        this%group_id(isnap) = i
     else
        write(*,*) 'Did not find name ', trim(data_name)
        this%ngroups = this%ngroups + 1
        this%group_id(isnap) = this%ngroups
        this%data_group_names(this%ngroups) = data_name
     endif
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_tet_mesh_data_xdmf(this, filename)
  class(inversion_mesh_data_type)   :: this
  character(len=*), intent(in)      :: filename
  integer                           :: iinput_xdmf, iinput_heavy_data
  integer                           :: i, igroup, itime, isnap, n

  if (.not. this%initialized) &
     stop 'ERROR: trying to dump a non initialized mesh'

  if (.not. allocated(this%datat)) &
     stop 'ERROR: no data to dump available'
  
  if (this%nvertices_per_elem /= 4) &
     stop 'ERROR: nvertices_per_elem /= 4, so this seems not to be a tet mesh'

  ! XML header
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  ! start xdmf file, write header
  write(iinput_xdmf, 733) this%nelements, 'binary', trim(filename)//'_grid.dat', &
                          this%nvertices, 'binary', trim(filename)//'_points.dat'

  i = 1
  itime = 1

  ! loop over all data
  do while(i <= this%ntimes)
     ! create new snapshot in the temporal collection
     write(iinput_xdmf, 7341) 'grid', dble(itime), this%nelements, "'", "'", "'", "'"

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
                                    this%nvertices, trim(filename)//'_data.dat'
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
    '<DataItem Name="grid" Dimensions="',i10,' 4" NumberType="Int" Format="',A,'">',/&
    '  ', A,/&
    '</DataItem>',/&
    '<DataItem Name="points" Dimensions="',i10,' 3" NumberType="Float" Format="',A,'">',/&
    '  ', A,/&
    '</DataItem>',/,/&
    '<Grid Name="CellsTime" GridType="Collection" CollectionType="Temporal">',/)

7341 format(&    
    '    <Grid Name="', A,'" GridType="Uniform">',/&
    '        <Time Value="',F8.2,'" />',/&
    '        <Topology TopologyType="Tetrahedron" NumberOfElements="',i10,'">',/&
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
  write(iinput_heavy_data) this%vertices
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
     write(iinput_heavy_data) this%datat(:,i)
  enddo
  close(iinput_heavy_data)

  
end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
