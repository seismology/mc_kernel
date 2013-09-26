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
     contains
     procedure, pass :: get_nelements
     procedure, pass :: get_nvertices
  end type

contains

!-----------------------------------------------------------------------------------------
integer function get_nelements(this)
  class(inversion_mesh_type)        :: this
  get_nelements = this%nelements
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_nvertices(this)
  class(inversion_mesh_type)        :: this
  get_nvertices = this%nvertices
end function
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
