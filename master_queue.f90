!=========================================================================================
module master_queue
  use global_parameters,           only: sp, dp, lu_out
  use inversion_mesh,              only: inversion_mesh_data_type
  use type_parameter,              only: parameter_type
  use simple_routines,             only: lowtrim
  implicit none
  private
  
  public :: init_queue, get_next_task, extract_receive_buffer, finalize

  type(inversion_mesh_data_type)      :: inv_mesh
  type(parameter_type)                :: parameters
  integer, allocatable                :: tasks(:), elems_in_task(:,:)
  real(kind=dp), allocatable          :: K_x(:,:), Var(:,:), Bg_model(:,:)
  integer,       allocatable          :: connectivity(:,:)
  integer,       allocatable          :: niterations(:,:), element_proc(:)

contains

!-----------------------------------------------------------------------------------------
subroutine init_queue(ntasks)
! anything that should be done before starting the loop over the work. for now,
! the number of tasks is fixed here

  use work_type_mod, only    : wt
  integer, intent(out)      :: ntasks
  integer                   :: itask, nelems, iel
  character(len=64)         :: fmtstring
  
  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Read input files for parameters, source and receivers'
  write(lu_out,'(A)') '***************************************************************'
  call parameters%read_parameters()
  call parameters%read_source()
  call parameters%read_receiver()

  ! Master and slave part ways here for some time. 
  ! Master reads in the inversion mesh, slaves initialize the FFT

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Read inversion mesh'
  write(lu_out,'(A)') '***************************************************************'
  select case(lowtrim(parameters%mesh_file_type))
  case('tetrahedral') 
    call inv_mesh%read_tet_mesh(parameters%mesh_file_vert, parameters%mesh_file_face)
  case('abaqus')
    call inv_mesh%read_abaqus_mesh(parameters%mesh_file, parameters%int_type)
  end select
  
  wt%ielement_type = inv_mesh%get_element_type()
  print *, 'Inversion mesh type: ', wt%ielement_type
  
  nelems    = inv_mesh%get_nelements()
  allocate(connectivity(inv_mesh%nvertices_per_elem, nelems))
  connectivity = inv_mesh%get_connectivity()
  
  ! Master and slave synchronize again

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Define filters'
  write(lu_out,'(A)') '***************************************************************'
  call parameters%read_filter()

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Define kernels'
  write(lu_out,'(A)') '***************************************************************'
  call parameters%read_kernel()

  allocate(niterations(parameters%nkernel, nelems))
  allocate(element_proc(nelems))
  
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
  
  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Allocate variables to store result'
  write(lu_out,'(A)') '***************************************************************'
  allocate(K_x(inv_mesh%get_nbasisfuncs(parameters%int_type), parameters%nkernel))
  K_x = 0.0
  allocate(Var(inv_mesh%get_nbasisfuncs(parameters%int_type), parameters%nkernel))
  Var = 0.0
  allocate(Bg_model(inv_mesh%get_nbasisfuncs(parameters%int_type),  &
                                             parameters%nmodel_parameter))
  Bg_Model = 0.0

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Starting to distribute the work'
  write(lu_out,'(A)') '***************************************************************'
  
end subroutine init_queue
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_next_task(itask)
! put a new piece of work in the send buffer
  
  use work_type_mod, only    : wt
  integer, intent(in)   :: itask
  integer               :: iel, ielement, ivert
  integer, allocatable  :: ivertex(:)

  allocate(ivertex(inv_mesh%nvertices_per_elem))

  wt%itask = itask

  do iel = 1, parameters%nelems_per_task
      ielement = elems_in_task(itask, iel)
      if (ielement.eq.-1) cycle
      ivertex = [( (iel-1) * inv_mesh%nvertices_per_elem + ivert, &
                   ivert = 1, inv_mesh%nvertices_per_elem )]
      wt%vertices(:, ivertex) = inv_mesh%get_element(ielement)
      wt%connectivity(:, iel) = ivertex
  end do


end subroutine get_next_task
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine extract_receive_buffer(itask, irank)
! extract information received back from a slave

  use work_type_mod, only    : wt
  integer, intent(in)       :: itask, irank
  integer                   :: iel, ielement, ibasisfunc

  ! extract from receive buffer
  do iel = 1, parameters%nelems_per_task
      ielement = elems_in_task(itask, iel)
      if (ielement.eq.-1) cycle
      do ibasisfunc = 1, inv_mesh%nbasisfuncs_per_elem

         ! are we in volumetric or vertex mode?
         select case(trim(parameters%int_type))
         case('onvertices')         
            K_x(connectivity(ibasisfunc, ielement),:) = K_x(connectivity(ibasisfunc, ielement),:) &
                                                        + wt%kernel_values(:, ibasisfunc, iel)
            Var(connectivity(ibasisfunc, ielement),:) = Var(connectivity(ibasisfunc, ielement),:) &
                                                        + wt%kernel_variance(:, ibasisfunc, iel) 
            Bg_Model(connectivity(ibasisfunc, ielement),:) = Bg_Model(connectivity(ibasisfunc, ielement),:) &
                                                        + wt%model(:, ibasisfunc, iel) 
         case('volumetric')
            K_x(ielement,:)      = K_x(ielement,:) &
                                   + wt%kernel_values(:, ibasisfunc, iel) 
            Var(ielement,:)      = Var(ielement,:) &
                                   + wt%kernel_variance(:, ibasisfunc, iel)   
            Bg_Model(ielement,:) = Bg_Model(ielement,:) &
                                   + wt%model(:, ibasisfunc, iel)   
         end select


      end do
      niterations(:,ielement)  = wt%niterations(:,iel)
      element_proc(ielement) = irank 
  end do


end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine finalize()
  integer                       :: ikernel, imodel_parameter
  real(kind=sp), allocatable    :: rel_error(:)
  character(len=64)             :: xdmf_varname
  real(kind=dp)                 :: total_time

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') 'Initialize output file'
  write(lu_out,'(A)') '***************************************************************'

  ! Save big kernel variable to disk
  write(lu_out,*) 'Write Kernel to disk'

  ! Select dump format
  select case(trim(parameters%dump_type))

  ! Save kernels in xdmf format
  case ('xdmf')
  
     ! Distiguish between volumetric and node-based mode
     select case(trim(parameters%int_type))
     case ('onvertices')
        call inv_mesh%init_node_data(parameters%nkernel*2+parameters%nmodel_parameter)
        allocate(rel_error(size(Var,1)))

        do ikernel = 1, parameters%nkernel
           xdmf_varname = 'K_x_'//parameters%kernel(ikernel)%name
           call inv_mesh%set_node_data_snap(real(K_x(:,ikernel), kind=sp), &
                                            ikernel,                       &     
                                            trim(xdmf_varname))

           rel_error = real(sqrt(Var(:,ikernel))/abs(K_x(:,ikernel)), kind=sp)
           xdmf_varname = 'err_'//parameters%kernel(ikernel)%name 
           call inv_mesh%set_node_data_snap(rel_error,                          &
                                            ikernel+parameters%nkernel,         &
                                            trim(xdmf_varname))
                                           
        end do
        
        do imodel_parameter = 1, parameters%nmodel_parameter
           xdmf_varname = 'Bg_model_'//parameters%bgmodel_parameter_names(imodel_parameter)
           call inv_mesh%set_node_data_snap(real(Bg_Model(:, imodel_parameter), kind=sp), &
                                            imodel_parameter+parameters%nkernel*2,  &
                                            trim(xdmf_varname))
        end do

        call inv_mesh%dump_node_data_xdmf(trim(parameters%output_file)//'_kernel')

     case ('volumetric')
        call inv_mesh%init_cell_data(parameters%nkernel*2+parameters%nmodel_parameter)
        allocate(rel_error(size(Var,1)))

        do ikernel = 1, parameters%nkernel
           xdmf_varname = 'K_x_'//parameters%kernel(ikernel)%name
           call inv_mesh%set_cell_data_snap(real(K_x(:,ikernel), kind=sp), &
                                            ikernel,                       &
                                            trim(xdmf_varname))

           rel_error = real(sqrt(Var(:,ikernel))/abs(K_x(:,ikernel)), kind=sp)
           xdmf_varname = 'err_'//parameters%kernel(ikernel)%name
           call inv_mesh%set_cell_data_snap(rel_error,                   &
                                            ikernel+parameters%nkernel,  &
                                            trim(xdmf_varname))
        end do
        
        do imodel_parameter = 1, parameters%nmodel_parameter
           xdmf_varname = 'Bg_model_'//parameters%bgmodel_parameter_names(imodel_parameter)
           call inv_mesh%set_cell_data_snap(real(Bg_Model(:, imodel_parameter), kind=sp), &
                                            imodel_parameter+parameters%nkernel*2,  &
                                            trim(xdmf_varname))
        end do

        call inv_mesh%dump_cell_data_xdmf(trim(parameters%output_file)//'_kernel')
     end select
  
  ! Save kernels in Yale-style csr format
  case ('csr')

     select case(trim(parameters%int_type))
     case ('onvertices')

        call inv_mesh%dump_node_data_csr ( real(K_x(:,:), kind=sp),  &
                                           parameters%nkernel,       &
                                           parameters%allowed_error, & 
                                           trim(parameters%output_file)//'_kernel')


     case ('volumetric')

        call inv_mesh%dump_cell_data_csr ( real(K_x(:,:), kind=sp),  &
                                           parameters%nkernel,       &
                                           parameters%allowed_error, & 
                                           trim(parameters%output_file)//'_kernel')


     end select

  ! Save kernels in ASCII format
  case ('ascii')

     select case(trim(parameters%int_type))
     case ('onvertices')

        call inv_mesh%dump_node_data_ascii ( real(K_x(:,:), kind=sp),  &
                                           parameters%nkernel,         &
                                           parameters%allowed_error,   & 
                                           trim(parameters%output_file)//'_kernel')


     case ('volumetric')

        call inv_mesh%dump_cell_data_ascii ( real(K_x(:,:), kind=sp),  &
                                           parameters%nkernel,         &
                                           parameters%allowed_error,   & 
                                           trim(parameters%output_file)//'_kernel')

     end select


  end select

  call inv_mesh%free_node_and_cell_data()

  ! Save mesh partition and convergence information
  write(lu_out,*) 'Write mesh partition and convergence to disk'
  call inv_mesh%init_cell_data(parameters%nkernel + 1)
  call inv_mesh%set_cell_data_snap(real(element_proc, kind=sp), 1,  &
                                   'element_proc')
  do ikernel = 1, parameters%nkernel
      call inv_mesh%set_cell_data_snap(real(niterations(ikernel, :), kind=sp), 1+ikernel,&
                                       'nit_'//parameters%kernel(ikernel)%name)
  end do 
  call inv_mesh%dump_cell_data_xdmf(trim(parameters%output_file)//'_kernel_stat')

  call inv_mesh%freeme()

  ! Multiply kernels with model
  do ikernel = 1, parameters%nkernel
    total_time = sum(K_x(:, ikernel) * Bg_model(:, ikernel))
    print '(A,": ",E15.5," s")', parameters%kernel(ikernel)%name, total_time
  end do

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
