module kernel_meshes

    use global_parameters
    use commpi
!    use netcdf
    use mesh_types
    implicit none

    contains
!-----------------------------------------------------------------------------------------

subroutine init_read_decomp_save_rot_kernel_mesh(params, kernel_mesh, &
    backgroundmodel)

    
    use parameter_types, only            : paramtype
    type(paramtype), intent(in)         :: params
    type(kernelmeshtype), intent(inout) :: kernel_mesh
    type(modeltype), intent(inout)      :: backgroundmodel
    type(meshtype)                      :: sem_mesh

    call read_sem_meshes(sem_mesh, backgroundmodel, params%io%fwd(1))

    call read_external_mesh(kernel_mesh, backgroundmodel%maxr, params%io%ext_mesh_name)

    call decompose_kernel_mesh(kernel_mesh)

    call find_gridpoints(sem_mesh, kernel_mesh, params%k%thr, params%k%phr)
    
    call save_kernelmesh(kernel_mesh)

end subroutine
!-----------------------------------------------------------------------------------------

subroutine read_sem_meshes(mesh, model, nc)
    use netcdf
    use parameter_types, only:             ncparamtype
    use nc_routines, only:            nc_read_att_char, check, getvarid

    type(modeltype), intent(inout) :: model
    type(meshtype), intent(inout)  :: mesh
    type(ncparamtype), intent(in)  :: nc
    integer                        :: ncvarid_mesh_s, ncvarid_mesh_z, ncvarid_idom
    integer                        :: nstrain
    integer                        :: ncvarid_disc
    !include 'mesh_params.h'
    !include 'mesh_params_kernel.h'
   
    !!@TODO: This is a mess. Fix it, SCS!
    !mesh%npoints  = (iend-ibeg+1)**2 * nelem

    call check( nf90_get_att(ncid = nc%ncid, &
                             name = 'npoints', &
                             varid = NF90_GLOBAL, &
                             values = mesh%npoints) ) 

    !if (nstrain.ne.mesh%npoints) then
    !    write(6,*) 'Inconsistent number of points: ', nstrain, ' in NetCDF file',&
    !        mesh%npoints, ' in mesh_params.h'
    !    stop
    !end if

    if (mynum.eq.0) write(6,*) 'Read SEM mesh from first forward simulation'
    
    allocate(mesh%s(mesh%npoints))
    allocate(mesh%z(mesh%npoints))
    allocate(mesh%idom(mesh%npoints))
    
    call nc_read_att_char(model%background_model, 'background model', nc)
    
    call  getvarid( ncid  = nc%mesh, &
                    name  = "mesh_S", &
                    varid = ncvarid_mesh_s) 
    call  getvarid( ncid  = nc%mesh, &
                    name  = "mesh_Z", &
                    varid = ncvarid_mesh_z) 
    call  getvarid( ncid  = nc%mesh, &
                    name  = "model_domain", &
                    varid = ncvarid_idom) 
    call  getvarid( ncid  = nc%mesh, &
                    name  = "disc_depths", &
                    varid = ncvarid_disc) 

    call check( nf90_get_var(ncid  = nc%mesh, &
                             varid = ncvarid_mesh_s, &
                             values = mesh%s )) 
    call check( nf90_get_var(ncid  = nc%mesh, &
                             varid = ncvarid_mesh_z, &
                             values = mesh%z )) 
    call check( nf90_get_var(ncid  = nc%mesh, &
                             varid = ncvarid_idom, &
                             values = mesh%idom )) 
    
    call check( nf90_get_att(ncid = nc%mesh, &
                             name = 'ndisc', &
                             varid = NF90_GLOBAL, &
                             values = model%ndisc) ) 
    allocate(model%discs(model%ndisc))
    call check( nf90_get_var(ncid  = nc%mesh, &
                             varid = ncvarid_disc, &
                             values = model%discs )) 
    model%maxr = maxval(model%discs)

end subroutine
 

!------------------------------------------------------------------------------

subroutine read_external_mesh(mesh, maxr, ext_mesh_name)
    use coord_trafo, only : sphi2xy,rthetaphi2xyz
    !use data_mesh, only: ext_mesh_name, lfext, maxr
    !use data_mesh, only: nproc, mynum
    type(kernelmeshtype), intent(inout) :: mesh
    real, intent(in)                    :: maxr
    character(len=200), intent(in)      :: ext_mesh_name
    real(kind=realkind)   :: x1, x2, x3, r
    integer               :: ipt, lfext
  
    if (mynum == 0) write(6,*) 'loading external mesh...'
  
    if (mynum == 0) write(6,*) ' External mesh:', trim(ext_mesh_name)
    open(unit=91, file=trim(ext_mesh_name), status='old')
    
    read(91,*) mesh%npoints_all
  
    if (mynum == 0) then
       write(6,*)
       write(6,*) 'External mesh total grid points:', mesh%npoints_all
       write(6,*) 'Kernel procs:', nproc
    endif
  
    allocate(mesh%x_all(1:mesh%npoints_all))
    allocate(mesh%y_all(1:mesh%npoints_all))
    allocate(mesh%z_all(1:mesh%npoints_all))
    allocate(mesh%fwd_pointid(1:mesh%npoints_all))
    allocate(mesh%bwd_pointid(1:mesh%npoints_all))
    allocate(mesh%idom(1:mesh%npoints_all))
  
    lfext = index(ext_mesh_name,' ')-1
 
    select case(ext_mesh_name(lfext-2:lfext))
    case('sph')
        if(mynum == 0) write(6,*)'reading external mesh in spherical coordinates...'
        do ipt=1, mesh%npoints_all
            read(91,*) x1, x2, x3
            if (x1 > maxr .or. x1 < 0.) then 
                write(6,*) mynum, 'problem with radius in external mesh:', ipt, x1
                stop
            endif
            if (x2 > pi .or. x2 < 0.) then 
                write(6,*) mynum, 'problem with theta in external mesh:', ipt, x2
                stop
            endif
            if (x3 > 2.*pi .or. x3 < 0.) then 
                write(6,*) mynum, 'problem with phi in external mesh:', ipt, x3
                stop
            endif
  
            call rthetaphi2xyz(mesh%x_all(ipt), mesh%y_all(ipt), mesh%z_all(ipt), x1, x2, x3)
        enddo
  
    case('cyl') 
    !elseif (ext_mesh_name(lfext-2:lfext) == 'cyl') then 
        if (mynum == 0) write(6,*) 'reading external mesh in cylindrical coordinates..'
        do ipt=1, mesh%npoints_all
            read(91,*) x1, x2, x3
            if (x1 > maxr .or. x1 < 0.) then 
                write(6,*) mynum, 'problem with coordinate s in external mesh:', ipt, x1
                stop
            endif
            if (abs(x2) > maxr) then 
                write(6,*) mynum, 'problem with coordinate z in external mesh:', ipt, x2
                stop
            endif
            if (x3 > 2.*pi .or. x3 < 0.) then 
                write(6,*) mynum, 'problem with phi in external mesh:', ipt, x3
                stop
            endif
  
            call sphi2xy(mesh%x_all(ipt), mesh%y_all(ipt), x1, x3)
            mesh%z_all(ipt) = x2
        enddo
  
    case('xyz')
     !elseif (ext_mesh_name(lfext-2:lfext) == 'xyz') then 
        if(mynum == 0) write(6,*) 'reading external mesh in cartesian coordinates...'
        do ipt=1, mesh%npoints_all
            read(91,*) mesh%x_all(ipt), mesh%y_all(ipt), mesh%z_all(ipt)
            r = sqrt(mesh%x_all(ipt)**2 + mesh%y_all(ipt)**2 + mesh%z_all(ipt)**2)
            if (r > 1.0001 * maxr) then
                write(6,*) mynum, 'problem with coordinates in external mesh, radius:', ipt, r
                stop
            elseif (abs(mesh%x_all(ipt)) > maxr) then
                write(6,*) mynum, 'problem with x in external mesh:', ipt, mesh%x_all(ipt)
                stop
            elseif (abs(mesh%y_all(ipt)) > maxr) then
                write(6,*) mynum, 'problem with y in external mesh:', ipt, mesh%y_all(ipt)
                stop
            elseif (abs(mesh%z_all(ipt)) > maxr) then
                write(6,*) mynum, 'problem with z in external mesh:', ipt, mesh%z_all(ipt)
                stop
            endif
        enddo
    
    case default
        if(mynum == 0) then 
            write(6,*) 'external mesh type not recognized!'
            write(6,*) 'mesh file name must end with:'
            write(6,*) " 1) 'xyz' -- cartesian coordinates (x[m], y[m], z[m])"
            write(6,*) " 2) 'sph' -- spherical coordinates (r[m], theta[rad], phi[rad])"
            write(6,*) " 3) 'cyl' -- cylindrical coordinates (s[m], z[m], phi[rad])"
        endif
        stop
    end select
  
    close(91)
  
end subroutine read_external_mesh

!------------------------------------------------------------------------------

subroutine find_gridpoints(sem_mesh, kernel_mesh, theta, phi)
    use kdtree2_module
    use coord_mapping, only: rotate_frame_rd

    type(meshtype),       intent(in)    :: sem_mesh
    type(kernelmeshtype), intent(inout) :: kernel_mesh
    real,                 intent(in)    :: theta, phi  

    type(kdtree2_result), allocatable  :: nextpoint(:)
    type(kdtree2),pointer              :: tree
    integer                            :: ipt
    real, dimension(kernel_mesh%npoints) :: rotmesh_s, rotmesh_phi, rotmesh_z

    tree => kdtree2_create(reshape([sem_mesh%s, sem_mesh%z],[2, sem_mesh%npoints]),&
                           dim = 2, &
                           sort = .false.,&
                           rearrange = .true.)
    
    allocate(nextpoint(1))
    do ipt = 1, kernel_mesh%npoints
        call kdtree2_n_nearest(tree, &
                               [kernel_mesh%x(ipt), kernel_mesh%z(ipt)], &
                               nn = 1, &
                               results = nextpoint )
        kernel_mesh%fwd_pointid(ipt) = nextpoint(1)%idx
        kernel_mesh%idom(ipt) = sem_mesh%idom(nextpoint(1)%idx) !!@HACK
    end do

    call rotate_frame_rd(kernel_mesh%npoints, rotmesh_s, rotmesh_phi, rotmesh_z, &
                         kernel_mesh%x, kernel_mesh%y, kernel_mesh%z, phi, theta)

    do ipt = 1, kernel_mesh%npoints
        call kdtree2_n_nearest(tree, &
                               [rotmesh_s(ipt), rotmesh_z(ipt)], &
                               nn = 1, &
                               results = nextpoint )
        
        kernel_mesh%bwd_pointid(ipt) = nextpoint(1)%idx
    end do

    call kdtree2_destroy(tree)

end subroutine

!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------
!> Spreads kernel mesh over processors
subroutine decompose_kernel_mesh(mesh)
  
!    use input_output, only : write_vtk_bin_scal_topology
    use commpi
  
    type(kernelmeshtype), intent(inout) :: mesh
    integer :: ipt, nptsperproc, npadd, istart

! decompose the mesh in an embarrassingly, almost 100%-balanced fashion...
    nptsperproc = mesh%npoints_all / nproc

    npadd = 0
    if (mod(mesh%npoints_all, nproc).ne.0) then 
       if (mynum<nproc-1) then 
           mesh%npoints = nptsperproc
       else
           npadd = mod(mesh%npoints_all, nproc)
           mesh%npoints = nptsperproc + npadd
       endif
    else
       mesh%npoints = nptsperproc 
    endif

    if (nproc>1) call barrier 
    write(6,*) mynum, 'reporting: total kernel points:', mesh%npoints_all
    write(6,*) mynum, 'has', mesh%npoints, 'grid points'
    if (nproc>1) call barrier

    allocate(mesh%x(mesh%npoints)) 
    allocate(mesh%y(mesh%npoints)) 
    allocate(mesh%z(mesh%npoints)) 
    
    istart = mynum * nptsperproc
    if (mynum==nproc-1 .and. npadd.ne.0 ) then
        mesh%x = mesh%x_all(istart+1:)
        mesh%y = mesh%y_all(istart+1:)
        mesh%z = mesh%z_all(istart+1:)
    else
        mesh%x = mesh%x_all(istart+1:istart+nptsperproc)
        mesh%y = mesh%y_all(istart+1:istart+nptsperproc)
        mesh%z = mesh%z_all(istart+1:istart+nptsperproc)
    end if
    
    !write(6,*)'minmax xgd:',mynum,minval(xgd),maxval(xgd)
    !write(6,*)'minmax xgd_tot:',mynum,minval(xgd_tot),maxval(xgd_tot)
    !write(6,*)'minmax ygd:',mynum,minval(ygd),maxval(ygd)
    !write(6,*)'minmax ygd_tot:',mynum,minval(ygd_tot),maxval(ygd_tot)
    !write(6,*)'minmax zgd:',mynum,minval(zgd),maxval(zgd)
    !write(6,*)'minmax zgd_tot:',mynum,minval(zgd_tot),maxval(zgd_tot)

    deallocate(mesh%x_all)
    deallocate(mesh%y_all)
    deallocate(mesh%z_all)

    if (nproc>1) call barrier 
    if (mynum==0) write(6,*) 'finished the domain decomposition of the kernel mesh.'
    if (nproc>1) call barrier 

    !call define_io_appendix(appiproc,mynum)
    !filename1='Data/mesh_all_cell_'//appiproc    
    !call write_VTK_bin_scal_topology(real(xgd),real(ygd),real(zgd),zgd,npts/4,filename1)

end subroutine decompose_kernel_mesh
!----------------------------------------------------------------------------
subroutine save_kernelmesh(mesh)
    use coord_trafo, only : xyz2rthetaphi

    type(kernelmeshtype), intent(in) :: mesh
    integer                          :: ipt
    real(kind=realkind)              :: r,th,ph
    character(len=4)                 :: appmynum

    call define_io_appendix(appmynum, mynum)
    open(unit=19,file='Data/kernelmesh_xyz_'//appmynum//'.dat')
    open(unit=20,file='Data/kernelmesh_rthph_'//appmynum//'.dat')
    
    do ipt=1,mesh%npoints
       write(19,*) real(mesh%x(ipt)), real(mesh%y(ipt)), real(mesh%z(ipt))
       call xyz2rthetaphi(r, th, ph, mesh%x(ipt), mesh%y(ipt), mesh%z(ipt))
       write(20,*) real(r), real(th), real(ph), mesh%idom(ipt)
    enddo
    close(20)
    close(19)
    
end subroutine save_kernelmesh

!subroutine decompose_mesh(fullmesh, meshpart, nparts)
!    type(kernelmeshtype), intent(in)  :: fullmesh
!    type(kernelmeshtype), allocatable, intent(out)  :: meshpart(:)
!    integer                                         :: nparts
!
!    integer                                         :: ifilememory
!    real(8)                                         :: memory_max, memory_needed 
!
!    open(ifilememory, file='input_memory.dat' status='old')  
!    read(ifilememory, *) memory_max
!    close(ifilememory
!
!    memory_needed = real(fullmesh%npoints * nomega)  
!
!
!end subroutine

end module
