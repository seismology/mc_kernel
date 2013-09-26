!=========================================================================================
module inversion_mesh

  private
  public :: inversion_mesh_type

  integer, parameter    :: sp = 4, dp = 8

  type :: inversion_mesh_type
     private
     integer                        :: nelements, nvertices
     integer, allocatable           :: connectivity(:,:)
     real(kind=sp), allocatable     :: vertices(:,:)
     logical                        :: initialized = .false.
     contains
     procedure, pass :: get_nelements
     procedure, pass :: get_nvertices
     procedure, pass :: read_tet_mesh
     procedure, pass :: dump_tet_mesh_xdmf
     procedure, pass :: freeme
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
integer function get_nvertices(this)
  class(inversion_mesh_type)        :: this
  if (.not. this%initialized) &
     stop 'ERROR: accessing inversion mesh type that is not initialized'
  get_nvertices = this%nvertices
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_tet_mesh(this, filename_vertices, filename_connectivity)
  class(inversion_mesh_type)        :: this
  character(len=*), intent(in)      :: filename_vertices
  character(len=*), intent(in)      :: filename_connectivity
  integer                           :: iinput_vertices, iinput_connectivity
  integer                           :: i

  ! read vertices
  open(newunit=iinput_vertices, file=trim(filename_vertices), status='old', action='read')

  read(iinput_vertices,*) this%nvertices

  allocate(this%vertices(3,this%nvertices))
  do i=1, this%nvertices
     read(iinput_vertices,*) this%vertices(:,i)
  enddo

  close(iinput_vertices)

  ! read connectivity (tetrahedral elements)
  open(newunit=iinput_connectivity, file=trim(filename_connectivity), status='old', &
       action='read')

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
  deallocate(this%vertices)
  deallocate(this%connectivity)
  this%initialized = .false.
end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================

program test_inversion_mesh
  use inversion_mesh
  type(inversion_mesh_type) :: inv_mesh

  !call inv_mesh%read_tet_mesh('vertices.TEST', 'facets.TEST')
  call inv_mesh%read_tet_mesh('vertices.USA10', 'facets.USA10')

  call inv_mesh%dump_tet_mesh_xdmf('testmesh')
  
  call inv_mesh%freeme()

end program
