!=========================================================================================
module master_queue
  use global_parameters,           only: sp, dp, lu_out
  use inversion_mesh,              only: inversion_mesh_data_type
  use type_parameter,              only: parameter_type
  implicit none
  private
  
  public :: init_queue, get_next_task, extract_receive_buffer, finalize

  type(inversion_mesh_data_type)      :: inv_mesh
  type(parameter_type)                :: parameters
  integer, allocatable                :: tasks(:), elems_in_task(:,:)
  real(kind=dp), allocatable          :: K_x(:,:), Err(:,:)
  integer,       allocatable          :: connectivity(:,:)

contains

!-----------------------------------------------------------------------------------------
subroutine init_queue(ntasks)
! anything that should be done before starting the loop over the work. for now,
! the number of tasks is fixed here

  use work_type_mod, only    : wt
  integer, intent(out)  :: ntasks
  integer               :: itask, nelems, iel
  character(len=64)     :: fmtstring
  
  write(lu_out,*) '***************************************************************'
  write(lu_out,*) ' Read input files for parameters, source and receivers'
  write(lu_out,*) '***************************************************************'
  call parameters%read_parameters('inparam_basic')
  call parameters%read_source()
  call parameters%read_receiver()


  write(lu_out,*) '***************************************************************'
  write(lu_out,*) ' Read inversion mesh'
  write(lu_out,*) '***************************************************************'
  !call inv_mesh%read_tet_mesh('vertices.USA10', 'facets.USA10')
  call inv_mesh%read_abaqus_mesh(parameters%mesh_file)
  
  wt%ielement_type = inv_mesh%get_element_type()
  !select case(inv_mesh%element_type)
  !case('tri')
  !   wt%ielement_type = 1
  !case('quad')
  !   wt%ielement_type = 2
  !case('hex')
  !   wt%ielement_type = 3
  !case('tet')
  !   wt%ielement_type = 4
  !end select

  
  nelems    = inv_mesh%get_nelements()
  allocate(connectivity(inv_mesh%nvertices_per_elem, nelems))
  connectivity = inv_mesh%get_connectivity()
  
  write(lu_out,*) '***************************************************************'
  write(lu_out,*) ' Define filters'
  write(lu_out,*) '***************************************************************'
  call parameters%read_filter()

  write(lu_out,*) '***************************************************************'
  write(lu_out,*) ' Define kernels'
  write(lu_out,*) '***************************************************************'
  call parameters%read_kernel()

  
  fmtstring = '(A, I8, A, I8)'
  ! Calculate number of tasks
  ntasks = ceiling(real(inv_mesh%get_nelements()) / parameters%nelems_per_task)
  print fmtstring, '  nelements: ',  nelems, ', ntasks: ', ntasks
   
  allocate(tasks(ntasks))
  allocate(elems_in_task(ntasks, parameters%nelems_per_task))
  do itask = 1, ntasks
     tasks(itask) = itask
     do iel = 1, parameters%nelems_per_task
         if (iel + (itask-1) * parameters%nelems_per_task <= nelems) then
             elems_in_task(itask, iel) = iel + (itask-1) * parameters%nelems_per_task
         else
             elems_in_task(itask, iel) = -1
         end if
     end do
  enddo
  
  allocate(K_x(inv_mesh%get_nvertices(), parameters%nkernel))
  K_x = 0.0
  allocate(Err(inv_mesh%get_nvertices(), parameters%nkernel))
  Err = 0.0

  write(lu_out,*) '***************************************************************'
  write(lu_out,*) ' Starting to distribute the work'
  write(lu_out,*) '***************************************************************'
  
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_next_task(itask)
! put a new piece of work in the send buffer
  
  use work_type_mod, only    : wt
  integer, intent(in)   :: itask
  integer               :: iel, ielement, ivert
  integer, allocatable  :: ivertex(:)

  allocate(ivertex(inv_mesh%nvertices_per_elem))

  !wt%ielements = elems_in_task(itask, :)
  wt%itask = itask

  do iel = 1, parameters%nelems_per_task
      ielement = elems_in_task(itask, iel)
      if (ielement.eq.-1) cycle
      ivertex = [( (iel-1) * inv_mesh%nvertices_per_elem + ivert, &
                   ivert = 1, inv_mesh%nvertices_per_elem )]
      wt%vertices(:, ivertex) = inv_mesh%get_element(ielement)
      wt%connectivity(:, iel) = ivertex
  end do

  !wt%ntotal_kernel = tasks(itask)
  !wt%ndimensions = 10

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine extract_receive_buffer(itask)
! extract information received back from a slave

  use work_type_mod, only    : wt
  integer, intent(in)       :: itask
  integer                   :: iel, ielement, ivertex

  ! extract from receive buffer
  do iel = 1, parameters%nelems_per_task
      !ielement = wt%ielements(iel)
      ielement = elems_in_task(itask, iel)
      if (ielement.eq.-1) cycle
      do ivertex = 1, inv_mesh%nvertices_per_elem
         K_x(connectivity(ivertex, ielement),:) = K_x(connectivity(ivertex, ielement),:) &
                                                  + wt%kernel_values(:, ivertex, iel)
         Err(connectivity(ivertex, ielement),:) = Err(connectivity(ivertex, ielement),:) &
                                                  + wt%kernel_errors(:, ivertex, iel)    &
                                                  / abs(wt%kernel_values(:, ivertex, iel))

      end do
  end do
  !output(ioutput,1) = wt%ntotal_kernel
  !output(ioutput,2) = wt%ndimensions

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine finalize()
  integer       :: ikernel, ivertex

  write(lu_out,*) '***************************************************************'
  write(lu_out,*) 'Initialize output file'
  write(lu_out,*) '***************************************************************'
  call inv_mesh%init_node_data(parameters%nkernel*2)

  ! Save big kernel variable to disk
  write(lu_out,*) 'Write Kernel to disk'
  
  do ikernel = 1, parameters%nkernel
     call inv_mesh%set_node_data_snap(real(K_x(:,ikernel), kind=sp), &
                                 ikernel, 'K_x_'//parameters%kernel(ikernel)%name )
     call inv_mesh%set_node_data_snap(real(Err(:,ikernel), kind=sp), &
                                 ikernel+parameters%nkernel,         &
                                 'Err_'//parameters%kernel(ikernel)%name )
  end do

  call inv_mesh%dump_node_data_xdmf('gaborkernel')
  call inv_mesh%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
