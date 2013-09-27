!=========================================================================================
module inversion_mesh

  implicit none
  private
  public :: inversion_mesh_type
  public :: inversion_mesh_data_type

  integer, parameter :: sp = 4, dp = 8

  type :: inversion_mesh_type
     private
     integer                            :: nelements, nvertices
     integer, allocatable               :: connectivity(:,:)
     real(kind=sp), allocatable         :: vertices(:,:)
     logical                            :: initialized = .false.
     contains
     procedure, pass :: get_nelements
     procedure, pass :: get_element
     procedure, pass :: get_elements
     procedure, pass :: get_nvertices
     procedure, pass :: get_vertices
     procedure, pass :: read_tet_mesh
     procedure, pass :: dump_tet_mesh_xdmf
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
  real(kind=sp)                     :: get_element(3,4)
  integer, intent(in)               :: ielement

  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'

  get_element(:,1) = this%vertices(:, this%connectivity(1,ielement)) 
  get_element(:,2) = this%vertices(:, this%connectivity(2,ielement)) 
  get_element(:,3) = this%vertices(:, this%connectivity(3,ielement)) 
  get_element(:,4) = this%vertices(:, this%connectivity(4,ielement)) 
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_elements(this)
  class(inversion_mesh_type)        :: this
  real(kind=sp)                     :: get_elements(3,4,this%nelements)
  integer                           :: ielement
  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'

  do ielement = 1, this%nelements
      get_elements(:,1,ielement) = this%vertices(:, this%connectivity(1,ielement)) 
      get_elements(:,2,ielement) = this%vertices(:, this%connectivity(2,ielement)) 
      get_elements(:,3,ielement) = this%vertices(:, this%connectivity(3,ielement)) 
      get_elements(:,4,ielement) = this%vertices(:, this%connectivity(4,ielement)) 
  end do
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
function get_connectivity(this)
  class(inversion_mesh_type)        :: this
  real(kind=sp)                     :: get_connectivity(4,this%nelements)
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

  allocate(this%connectivity(4,this%nelements))
  do i=1, this%nelements
     read(iinput_connectivity,*) this%connectivity(:,i)
  enddo

  close(iinput_connectivity)

  this%initialized = .true.
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_tet_mesh_xdmf(this, filename)
  class(inversion_mesh_type)        :: this
  character(len=*), intent(in)      :: filename
  integer                           :: iinput_xdmf, iinput_heavy_data
  integer                           :: i

  if (.not. this%initialized) &
     stop 'ERROR: trying to dump a non initialized mesh'

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  write(iinput_xdmf, 732) this%nelements, this%nelements, 'binary', &
                      trim(filename)//'_grid.dat', this%nvertices, 'binary', &
                      trim(filename)//'_points.dat'
  close(iinput_xdmf)

732 format(&    
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/&
    '  <Grid GridType="Uniform">',/&
    '    <Topology TopologyType="Tetrahedron" NumberOfElements="',i10,'">',/&
    '      <DataItem Dimensions="',i10,' 4" NumberType="Int" Format="',A,'">',/&
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

  this%ngroups = 1
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

  if (size(data_snap) /= this%nvertices) &
     stop 'ERROR: wrong dimensions of input data_snap for writing vertex data'

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
subroutine dump_tet_mesh_data_xdmf(this, filename)
  class(inversion_mesh_data_type)   :: this
  character(len=*), intent(in)      :: filename
  integer                           :: iinput_xdmf, iinput_heavy_data
  integer                           :: i, igroup, itime, isnap, n

  if (.not. this%initialized) &
     stop 'ERROR: trying to dump a non initialized mesh'

  if (.not. allocated(this%datat)) &
     stop 'ERROR: no data to dump available'

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
