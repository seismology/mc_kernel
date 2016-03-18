!=========================================================================================
module master_queue
  use global_parameters,           only: sp, dp, lu_out, long, verbose
  use inversion_mesh,              only: inversion_mesh_data_type
  use type_parameter,              only: parameter_type
  use simple_routines,             only: lowtrim
  use nc_routines,                 only: nc_create_file, nc_close_file, &
                                         nc_putvar_by_name, nc_getvar_by_name, &
                                         nc_open_for_read
  use background_model,            only: nmodel_parameters, parameter_name
  use heterogeneities,             only: nmodel_parameters_hetero,parameter_name_het, &
                                         parameter_name_het_store
  use master_helper,               only: create_tasks
  implicit none
  private
  
  public :: init_queue, get_next_task, extract_receive_buffer, finalize, & 
            dump_intermediate, delete_intermediate

  type(inversion_mesh_data_type), save :: inv_mesh
  type(parameter_type),           save :: parameters
  integer, allocatable,           save :: elems_in_task(:,:)
  real(kind=dp), allocatable,     save :: K_x(:,:), Var(:,:), Bg_Model(:,:), Het_Model(:,:)
  real(kind=sp), allocatable,     save :: fw_field(:,:,:,:), bw_field(:,:,:,:), conv_field(:,:,:)
  integer,       allocatable,     save :: connectivity(:,:)
  integer,       allocatable,     save :: niterations(:,:), element_proc(:)
  real(kind=dp), allocatable,     save :: computation_time(:)
  integer(kind=long),             save :: iclockold_mpi
  integer,                        save :: ncid_intermediate

contains

!-----------------------------------------------------------------------------------------
subroutine init_queue(ntasks, inparam_file)
! anything that should be done before starting the loop over the work. for now,
! the number of tasks is fixed here
! Reads the intermediate results file and returns the boolean array 'completed', which
! tells, which values have already been completed.

  use clocks_mod,    only                   : tick
  use global_parameters, only               : id_read_params, id_create_tasks
  use work_type_mod, only                   : wt
  integer, intent(out)                     :: ntasks
  character(len=*), intent(in), optional   :: inparam_file
  integer                                  :: nelems, nbasisfuncs
  integer(kind=long)                       :: iclockold
  character(len=64)                        :: filename
  logical                                  :: intermediate_results_exist
  logical, allocatable                     :: completed(:)
  
  iclockold = tick()

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Read input files for parameters, source and receivers'
  write(lu_out,'(A)') '***************************************************************'
  if (present(inparam_file)) then !only relevant for unit tests
    call parameters%read_parameters(trim(inparam_file))
  else
    call parameters%read_parameters()
  end if
  call parameters%read_source()
  call parameters%read_receiver()

  ! Master and slave part ways here for some time. 
  ! Master reads in the inversion mesh, slaves initialize the FFT

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Read inversion mesh'
  write(lu_out,'(A)') '***************************************************************'
  select case(lowtrim(parameters%mesh_file_type))
  case('tetrahedral') 
    call inv_mesh%read_tet_mesh(parameters%mesh_file_vert, &
                                parameters%mesh_file_face, &
                                parameters%int_type)
  case('abaqus')
    call inv_mesh%read_abaqus_mesh(parameters%mesh_file, &
                                   parameters%int_type)
  end select
  
  wt%ielement_type = inv_mesh%get_element_type()
  print *, 'Inversion mesh type: ', wt%ielement_type

  if (parameters%sort_mesh_elements) then
    print *, 'Sorting inversion mesh'
    call inv_mesh%tree_sort()
  end if
  
  nelems = inv_mesh%get_nelements()
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

  iclockold = tick(id=id_read_params, since=iclockold)

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Allocate variables to store result'
  write(lu_out,'(A)') '***************************************************************'
  ! Number of vertices or number of cells, depending on mesh type
  nbasisfuncs = inv_mesh%get_nbasisfuncs(parameters%int_type)
  allocate(K_x(nbasisfuncs, parameters%nkernel))
  allocate(Var(nbasisfuncs, parameters%nkernel))
  allocate(Bg_Model(nbasisfuncs, nmodel_parameters))
  allocate(Het_Model(nbasisfuncs, nmodel_parameters_hetero))
  K_x(:,:) = 0.0
  Var(:,:) = 0.0
  Bg_Model(:,:) = 0.0
  Het_Model(:,:) = 0.0

  if (parameters%plot_wavefields) then
    allocate(fw_field(nbasisfuncs, wt%ndumps, wt%ndim, &
                      parameters%nkernel))
    allocate(bw_field(nbasisfuncs, wt%ndumps, wt%ndim, &
                      parameters%nkernel))
    allocate(conv_field(nbasisfuncs, wt%ndumps, parameters%nkernel))
    fw_field = 0.0
    bw_field = 0.0
    conv_field = 0.0
  end if

  allocate(niterations(nelems, parameters%nkernel))
  allocate(computation_time(nelems))
  allocate(element_proc(nelems))
  niterations(:,:)     = -1
  computation_time(:)  = -1
  element_proc(:)      = -1

  ! Allocate variable to store whether task has already been completed
  allocate(completed(nelems))
  completed(:) = .false.
  
  if (parameters%create_intermediate) then
    ! Check for intermediate results and load them if available
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Check for intermediate results from previous run'
    write(lu_out,'(A)') '***************************************************************'
    filename = 'intermediate_results.nc'
    inquire( file = trim(filename), exist = intermediate_results_exist)

    ! If available, read (sets global variables)
    if (intermediate_results_exist) call read_intermediate(filename, completed)
  end if

  ! Create mapping array elems_in_task that contains which element belongs to which
  ! task. Uses the boolean array completed as input
  call create_tasks(completed, parameters%nelems_per_task, ntasks, elems_in_task)
  iclockold = tick(id=id_create_tasks, since=iclockold)

  if (parameters%create_intermediate) then
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Create file for intermediate results '
    write(lu_out,'(A)') '***************************************************************'
    call create_intermediate(filename)
  end if

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Starting to distribute the work'
  write(lu_out,'(A)') '***************************************************************'
  
end subroutine init_queue
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_next_task(itask)
! put a new piece of work in the send buffer
  use clocks_mod,    only    : tick
  use work_type_mod, only    : wt
  use global_parameters, only: id_get_next_task
  integer, intent(in)       :: itask
  integer                   :: iel, ielement, ivert
  integer(kind=long)        :: iclockold
  integer, allocatable      :: ivertex(:)

  iclockold = tick()
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
  
  iclockold_mpi = tick(id=id_get_next_task, since=iclockold)

end subroutine get_next_task
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine extract_receive_buffer(itask, irank)
! extract information received back from a slave

  use clocks_mod,    only    : tick
  use global_parameters, only: id_extract, id_mpi
  use work_type_mod, only    : wt
  integer, intent(in)       :: itask, irank
  integer                   :: iel, ielement, ibasisfunc, ipoint
  real(kind=dp)             :: valence
  integer(kind=long)        :: iclockold

  if (.not.iclockold_mpi==-1) then
    iclockold = tick(id=id_mpi, since=iclockold_mpi)
    iclockold_mpi = -1
  else
    iclockold = tick()
  end if

  ! extract from receive buffer
  ! iel = element number within task (in the local mesh of the slave)
  ! ielement = element number within the global mesh of the master
  do iel = 1, parameters%nelems_per_task
    ielement = elems_in_task(itask, iel)
    if (ielement.eq.-1) cycle
    do ibasisfunc = 1, inv_mesh%nbasisfuncs_per_elem

      ! are we in volumetric or vertex mode?
      select case(trim(parameters%int_type))
      case('onvertices')         
        ipoint = connectivity(ibasisfunc, ielement)
        valence = 1.d0/real(inv_mesh%nbasisfuncs_per_elem, kind=dp)
      case('volumetric')
        ipoint = ielement
        valence = 1.d0
      end select

      K_x(ipoint, :)       = K_x(ipoint,:)       + wt%kernel_values(:, ibasisfunc, iel) * valence
      Var(ipoint, :)       = Var(ipoint,:)       + wt%kernel_variance(:, ibasisfunc, iel)  * valence
      if (parameters%int_over_background) then
        Bg_Model(ipoint, :)  = Bg_Model(ipoint,:)  + wt%model(:, ibasisfunc, iel) * valence
      end if
      if (parameters%int_over_hetero) then
        Het_Model(ipoint, :) = Het_Model(ipoint,:) + wt%hetero_model(:, ibasisfunc, iel)  * valence
      end if
      if (parameters%plot_wavefields) then
        fw_field(ipoint, :, :, :) = fw_field(ipoint, :, :, :) + real(wt%fw_field(:, :, :, ibasisfunc, iel) * valence, kind=sp)
        bw_field(ipoint, :, :, :) = bw_field(ipoint, :, :, :) + real(wt%bw_field(:, :, :, ibasisfunc, iel) * valence, kind=sp)
        conv_field(ipoint, :, :)  = conv_field(ipoint, :, :)  + real(wt%conv_field(:, 1, :, ibasisfunc, iel) * valence, kind=sp)
      end if


    end do
    niterations(ielement,:)     = wt%niterations(:,iel)
    computation_time(ielement)  = wt%computation_time(iel)
    element_proc(ielement)      = irank 
  end do

  iclockold = tick(id=id_extract, since=iclockold)

end subroutine extract_receive_buffer
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine finalize()
  use work_type_mod, only        : wt
  integer                       :: ikernel
  real(kind=dp)                 :: total_time
  real(kind=dp)                 :: delay_time
  integer                       :: imodel, idim
  character(len=32)             :: fmtstring
  character(len=512)            :: var_name
  character(len=2), parameter   :: dim_name(6) = ['tt', 'pp', 'rr', 'pr', 'tr', 'tp']

  if (parameters%create_intermediate) then
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') 'Finalize intermediate file'
    write(lu_out,'(A)') '***************************************************************'
    call nc_close_file(ncid = ncid_intermediate)
  end if

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
       call inv_mesh%init_node_data()

       call inv_mesh%add_node_variable(var_name    = 'kernel',           &
                                       nentries    = parameters%nkernel, &
                                       entry_names = [('K_x_'//parameters%kernel(ikernel)%name, &
                                                       ikernel = 1, parameters%nkernel)] )

       call inv_mesh%add_node_variable(var_name    = 'error',            &
                                       nentries    = parameters%nkernel, &
                                       entry_names = [('err_'//parameters%kernel(ikernel)%name, &
                                                       ikernel = 1, parameters%nkernel)] )

                                                    
       call inv_mesh%add_node_data(var_name = 'kernel',                      &
                                   values   = real(K_x, kind=sp) )
       call inv_mesh%add_node_data(var_name = 'error',                       &
                                   values   = real(sqrt(Var/K_x), kind=sp) )

       if (parameters%int_over_background) then
         call inv_mesh%add_node_variable(var_name    = 'model',                                  &
                                         nentries    = nmodel_parameters,              &
                                         entry_names = parameter_name )

         call inv_mesh%add_node_data(var_name = 'model',                       &
                                     values   = real(Bg_Model, kind=sp) )
       end if
       if (parameters%int_over_hetero) then
         call inv_mesh%add_node_variable(var_name    = 'hetero',                                  &
                                         nentries    = nmodel_parameters_hetero,              &
                                         entry_names = parameter_name_het_store )

         call inv_mesh%add_node_data(var_name = 'hetero',                       &
                                     values   = real(Het_Model, kind=sp) )
       end if
                                     
    case ('volumetric')
       call inv_mesh%init_cell_data()

       call inv_mesh%add_cell_variable(var_name    = 'kernel',           &
                                       nentries    = parameters%nkernel, &
                                       entry_names = [('K_x_'//parameters%kernel(ikernel)%name, &
                                                       ikernel = 1, parameters%nkernel)] )

       call inv_mesh%add_cell_variable(var_name    = 'error',            &
                                       nentries    = parameters%nkernel, &
                                       entry_names = [('err_'//parameters%kernel(ikernel)%name, &
                                                       ikernel = 1, parameters%nkernel)] )

                                                    
       call inv_mesh%add_cell_data(var_name = 'kernel',                      &
                                   values   = real(K_x, kind=sp) )
       call inv_mesh%add_cell_data(var_name = 'error',                       &
                                   values   = real(sqrt(Var/K_x), kind=sp) )

       if (parameters%int_over_background) then
         call inv_mesh%add_cell_variable(var_name    = 'model',                                  &
                                         nentries    = nmodel_parameters,              &
                                         entry_names = parameter_name )

         call inv_mesh%add_cell_data(var_name = 'model',                       &
                                     values   = real(Bg_Model, kind=sp) )
       end if
       if (parameters%int_over_hetero) then
         call inv_mesh%add_cell_variable(var_name    = 'hetero',                                  &
                                         nentries    = nmodel_parameters_hetero,              &
                                         entry_names = parameter_name_het_store )
         call inv_mesh%add_cell_data(var_name = 'hetero',                       &
                                     values   = real(Het_Model, kind=sp) )
       end if

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


  ! Save mesh partition and convergence information
  write(lu_out,*) 'Write mesh partition and convergence to disk'

  call inv_mesh%init_cell_data()
  call inv_mesh%add_cell_variable(var_name    = 'niterations',      &
                                  nentries    = parameters%nkernel, &
                                  entry_names = [('nit_'//parameters%kernel(ikernel)%name, &
                                                  ikernel = 1, parameters%nkernel)] )
  call inv_mesh%add_cell_variable(var_name    = 'computation_time', &
                                  nentries    = 1 )                 
  call inv_mesh%add_cell_variable(var_name    = 'element_proc', &
                                  nentries    = 1 )

  call inv_mesh%add_cell_data(var_name = 'niterations',       &
                              values   = real(niterations, kind=sp))
  call inv_mesh%add_cell_data(var_name = 'computation_time',  &
                              values   = reshape(real(computation_time, kind=sp), [size(element_proc,1), 1]))
  call inv_mesh%add_cell_data(var_name = 'element_proc',      &
                              values   = reshape(real(element_proc, kind=sp), [size(element_proc,1), 1]))


  call inv_mesh%dump_data_xdmf(trim(parameters%output_file)//'_kernel')
      
  call inv_mesh%free_node_and_cell_data()

  if (parameters%plot_wavefields) then
    write(lu_out,*) 'Write wavefields to disk'
    select case(trim(parameters%int_type))
    case ('volumetric')
      call inv_mesh%init_cell_data(starttime = 0.d0, dt = real(wt%dt*1d-6, kind=dp))

      do ikernel = 1, parameters%nkernel
        do idim = 1, wt%ndim
          fmtstring = '("fw_", A, "_", A2)'
          write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name), dim_name(idim)
          call inv_mesh%add_cell_variable(var_name, nentries=wt%ndumps, istime=.true.)
          call inv_mesh%add_cell_data(var_name = var_name,                            &
                                      values   = fw_field(:, :, idim, ikernel))
                                      
          fmtstring = '("bw_", A, "_", A2)'
          write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name), dim_name(idim)
          call inv_mesh%add_cell_variable(var_name, nentries=wt%ndumps, istime=.true.)
          call inv_mesh%add_cell_data(var_name = var_name,                            &
                                      values   = bw_field(:, :, idim, ikernel))
                                      
        end do

        fmtstring = '("conv_", A, "_")'
        write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name)
        call inv_mesh%add_cell_variable(var_name, nentries=wt%ndumps, istime=.true.)
        call inv_mesh%add_cell_data(var_name = var_name,                            &
                                    values   = conv_field(:, :, ikernel))

      end do

    case ('onvertices')
      call inv_mesh%init_node_data(starttime = 0.d0, dt = real(wt%dt*1d-6, kind=dp))

      do ikernel = 1, parameters%nkernel
        do idim = 1, wt%ndim
          fmtstring = '("fw_", A, "_", A2)'
          write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name), dim_name(idim)
          call inv_mesh%add_node_variable(var_name, nentries=wt%ndumps, istime=.true.)
          call inv_mesh%add_node_data(var_name = var_name,                            &
                                      values   = fw_field(:, :, idim, ikernel))
                                      
          fmtstring = '("bw_", A, "_", A2)'
          write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name), dim_name(idim)
          call inv_mesh%add_node_variable(var_name, nentries=wt%ndumps, istime=.true.)
          call inv_mesh%add_node_data(var_name = var_name,                            &
                                      values   = bw_field(:, :, idim, ikernel))
                                      
        end do

        fmtstring = '("conv_", A, "_", I1)'
        write(var_name, fmtstring) trim(parameters%kernel(ikernel)%name), 1
        call inv_mesh%add_node_variable(var_name, nentries=wt%ndumps, istime=.true.)
        call inv_mesh%add_node_data(var_name = var_name,                            &
                                    values   = conv_field(:, :, ikernel))

      end do

    end select
    call inv_mesh%dump_data_xdmf(trim(parameters%output_file)//'_wavefield')
    call inv_mesh%free_node_and_cell_data()
  end if
    


  call inv_mesh%freeme()


  if (parameters%int_over_background) then
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') 'Minimum and maximum values of background model parameters'
    write(lu_out,'(A)') '***************************************************************'
    do imodel = 1, nmodel_parameters
      print *, imodel, parameter_name(imodel), minval(Bg_Model(:,imodel)), maxval(Bg_Model(:,imodel))
    end do
  end if

  if (parameters%int_over_hetero) then
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') 'Minimum and maximum values of 3D model parameters'
    write(lu_out,'(A)') '***************************************************************'
    do imodel = 1, nmodel_parameters_hetero
      print *, imodel, parameter_name_het(imodel), minval(Het_Model(:,imodel)), maxval(Het_Model(:,imodel))
    end do
  end if

  ! Multiply kernels with model to get absolute traveltime of "phase"
  if (parameters%int_over_background) then
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') 'Kernels times background model'
    write(lu_out,'(A)') '***************************************************************'
    do ikernel = 1, parameters%nkernel
      total_time = sum(K_x(:, ikernel) * Bg_Model(:, parameters%kernel(ikernel)%model_parameter_index))
      print '(A,": ",E15.5," s")', parameters%kernel(ikernel)%name, total_time
    end do
  end if

  ! Multiply relative kernels with heterogeneities to get traveltime shift
  if (parameters%int_over_hetero) then     
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') 'Kernels times 3D model'
    write(lu_out,'(A)') '***************************************************************'
     do ikernel = 1, parameters%nkernel
        if (parameters%kernel(ikernel)%hetero_parameter_index>0) then
          delay_time = sum(K_x(:, ikernel) * Het_Model(:, parameters%kernel(ikernel)%hetero_parameter_index))/100.d0
          print '(A,": ",E15.5," s")', parameters%kernel(ikernel)%name, delay_time
        else
          print '("Kernel ", A, " has model parameter not in het file")', parameters%kernel(ikernel)%name
        end if
     end do
  end if

     

end subroutine finalize
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Create intermediate NetCDF file
subroutine create_intermediate(filename)
  character(len=*), intent(in)  :: filename

  call nc_create_file(filename = filename, ncid = ncid_intermediate, overwrite=.true.)

  call nc_putvar_by_name(ncid    = ncid_intermediate, &
                         varname = 'K_x',             & 
                         values  = real(K_x, kind=sp) )
  call nc_putvar_by_name(ncid    = ncid_intermediate, &
                         varname = 'Var',             & 
                         values  = real(Var, kind=sp) )
  call nc_putvar_by_name(ncid    = ncid_intermediate, &
                         varname = 'Bg_Model',        & 
                         values  = real(Bg_Model, kind=sp) )
  call nc_putvar_by_name(ncid    = ncid_intermediate, &
                         varname = 'Het_Model',        & 
                         values  = real(Het_Model, kind=sp) )
  call nc_putvar_by_name(ncid    = ncid_intermediate, &
                         varname = 'niterations',     & 
                         values  = niterations )
  call nc_putvar_by_name(ncid    = ncid_intermediate, &
                         varname = 'element_proc',    & 
                         values  = element_proc )

end subroutine create_intermediate                       
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Read intermediate results and set to global variables
subroutine read_intermediate(filename, completed)
  character(len=*), intent(in)             :: filename
  logical, allocatable, intent(out)        :: completed(:)

  real(kind=sp), allocatable               :: real_temp(:,:)
  integer                                  :: ncid_old_intermediate

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Intermediate results found, loading'
  write(lu_out,'(A)') '***************************************************************'
  call nc_open_for_read(filename = filename, ncid = ncid_old_intermediate)

  call nc_getvar_by_name(ncid    = ncid_old_intermediate, &
                         varname = 'K_x',             & 
                         values  = real_temp )
  K_x = real(real_temp, kind=dp)
  deallocate(real_temp)
  call nc_getvar_by_name(ncid    = ncid_old_intermediate, &
                         varname = 'Var',             & 
                         values  = real_temp)
  Var = real(real_temp, kind=dp)
  deallocate(real_temp)
  call nc_getvar_by_name(ncid    = ncid_old_intermediate, &
                         varname = 'Bg_Model',        & 
                         values  = real_temp)
  Bg_Model = real(real_temp, kind=dp)
  deallocate(real_temp)
  call nc_getvar_by_name(ncid    = ncid_old_intermediate, &
                         varname = 'Het_Model',        & 
                         values  = real_temp)
  Het_Model = real(real_temp, kind=dp)
  deallocate(real_temp)
  call nc_getvar_by_name(ncid    = ncid_old_intermediate, &
                         varname = 'niterations',     & 
                         values  = niterations )
  call nc_getvar_by_name(ncid    = ncid_old_intermediate, &
                         varname = 'element_proc',    & 
                         values  = element_proc )

  ! Since master has proc 0, all slaves have proc ids >= 1.
  completed = (element_proc .ge. 1)

  call nc_close_file(ncid_old_intermediate)

end subroutine read_intermediate
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_intermediate(itask)
  use netcdf,            only : nf90_sync
  use clocks_mod,        only : tick
  use global_parameters, only : id_dump, long, int4
  integer, intent(in)        :: itask
  integer(kind=int4)         :: iel, ielement, ibasisfunc, ipoint, status
  integer(kind=long)         :: iclockold

! The check, whether to dump intermediate results at all needs to be done here, since the
! calling routine in master_mod.f90 has no access to the parameter object
if (parameters%create_intermediate) then
    iclockold = tick()

    do iel = 1, parameters%nelems_per_task
      ielement = elems_in_task(itask, iel)
      if (ielement.eq.-1) cycle

      ! are we in volumetric or vertex mode?
      select case(trim(parameters%int_type))
      case('onvertices')         
        do ibasisfunc = 1, inv_mesh%nbasisfuncs_per_elem
          ipoint = connectivity(ibasisfunc, ielement)
          call nc_putvar_by_name(ncid    = ncid_intermediate, &
                                 varname = 'K_x',             & 
                                 values  = real(K_x(ipoint:ipoint, :), kind=sp),    &
                                 start   = [ipoint, 1],       &
                                 count   = [1, parameters%nkernel] )
          call nc_putvar_by_name(ncid    = ncid_intermediate, &
                                 varname = 'Var',             & 
                                 values  = real(Var(ipoint, :), kind=sp),    &
                                 start   = [ipoint, 1],       &
                                 count   = [1, parameters%nkernel] )
          call nc_putvar_by_name(ncid    = ncid_intermediate, &
                                 varname = 'Bg_Model',        & 
                                 values  = real(Bg_Model(ipoint, :), kind=sp),&
                                 start   = [ipoint, 1],        &
                                 count   = [1, nmodel_parameters] )
          call nc_putvar_by_name(ncid    = ncid_intermediate, &
                                 varname = 'Het_Model',        & 
                                 values  = real(Het_Model(ipoint, :), kind=sp),&
                                 start   = [ipoint, 1],        &
                                 count   = [1, nmodel_parameters_hetero] )
        end do

      case('volumetric')
        call nc_putvar_by_name(ncid    = ncid_intermediate, &
                               varname = 'K_x',             & 
                               values  = real(K_x(ielement, :), kind=sp),    &
                               start   = [ielement, 1],       &
                               count   = [1, parameters%nkernel] )
        call nc_putvar_by_name(ncid    = ncid_intermediate, &
                               varname = 'Var',             & 
                               values  = real(Var(ielement, :), kind=sp),    &
                               start   = [ielement, 1],       &
                               count   = [1, parameters%nkernel] )
        call nc_putvar_by_name(ncid    = ncid_intermediate, &
                               varname = 'Bg_Model',        & 
                               values  = real(Bg_Model(ielement, :), kind=sp),&
                               start   = [ielement, 1],        &
                               count   = [1, nmodel_parameters] )
        call nc_putvar_by_name(ncid    = ncid_intermediate, &
                               varname = 'Het_Model',        & 
                               values  = real(Het_Model(ielement, :), kind=sp),&
                               start   = [ielement, 1],        &
                               count   = [1, nmodel_parameters_hetero] )

      end select ! int_type

      call nc_putvar_by_name(ncid    = ncid_intermediate, &
                             varname = 'niterations',        & 
                             values  = niterations(ielement, :),&
                             start   = [ielement, 1],        &
                             count   = [1, parameters%nkernel] )
    end do

    call nc_putvar_by_name(ncid    = ncid_intermediate, &
                           varname = 'element_proc',    & 
                           values  = element_proc)
    
    status = nf90_sync(ncid_intermediate)
    iclockold = tick(id=id_dump, since=iclockold)
  end if

end subroutine dump_intermediate
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine delete_intermediate
  character(len=512)                :: cmdmsg, sys_cmd
  integer                           :: exitstat

  if (parameters%create_intermediate) then
    ! Delete intermediate file
    sys_cmd = 'rm intermediate_results.nc'
    call execute_command_line(command = trim(sys_cmd), wait = .true., exitstat = exitstat, &
                              cmdmsg = cmdmsg)

    if (exitstat.ne.0) then
      write(*,*) 'WARNING: Deleting the intermediate file failed with status ', exitstat
      write(*,*) '         and message:'
      write(*,*) cmdmsg
    end if
  end if

end subroutine delete_intermediate
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
  if (wait) continue
end subroutine execute_command_line
# endif
end module
!=========================================================================================
