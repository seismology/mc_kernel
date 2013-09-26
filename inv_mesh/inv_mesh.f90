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
     !procedure, pass :: dump_mesh_xdmf
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
  character(len=*)                  :: filename_vertices
  character(len=*)                  :: filename_connectivity
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

  call inv_mesh%read_tet_mesh('vertices.TEST', 'facets.TEST')
  call inv_mesh%freeme()

end program
