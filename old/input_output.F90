!========================
module input_output
!========================

    use data_mesh
    use global_parameters

    implicit none

    public :: save_kernel, save_snapshot, save_misfit_kernel, save_devstrain
    public :: read_velocities_src, read_velocities_rec
    public :: read_summedstraintrace_src, read_summedstraintrace_rec
    public :: dump2d, dump1d
    public :: write_VTK_bin_scal, write_avs_file_scal
    public :: write_VTK_bin_scal_topology, write_VTK_bin_scal_mesh3d
    public :: read_deviatoricstrain_rec, read_deviatoricstrain_src, read_strain_src
#ifdef unc
    public :: open_netcdf,  close_netcdf
#endif
    integer, parameter                 :: nvar = 10
    character(len=16), dimension(nvar) :: varnamelist
    integer, dimension(nvar)           :: nc_varid
    integer                            :: ncid_snapshot_bw, ncid_snapshot_fw
  
    private
    contains 

!-----------------------------------------------------------------------------------------
subroutine read_summedstraintrace_src(u, appidump, dir_sim, it)

    use data_rota, only: iproc_d_src, ipt_d_src 
    include 'mesh_params.h'
    include 'mesh_params_kernel.h'
  
    character(len=4), intent(in)              :: appidump
    character(len=200), intent(in)            :: dir_sim
    integer, intent(in)                       :: it
  
    real(kind=realkind), dimension(npts), intent(out) :: u
  
    real(kind=realkind), dimension(nsize_sol) :: u_sol
    real(kind=realkind), dimension(nsize_flu) :: u_flu
    real(kind=realkind), dimension(nsize)     :: u1
    real(kind=realkind), dimension(nsize)     :: u2
    integer                                   :: ipt,imesh
    character(len=4)                          :: appimesh
    character(len=100)                        :: filename
    
    u = -123456
    do imesh = 1, num_meshes_fwd
        if(.not.(any(iproc_d_src==meshes_fwd(imesh)))) cycle
        call define_io_appendix(appimesh, meshes_fwd(imesh))
        call read_wavefield_1d(dir_sim, .true.,'straintrace_flu', & ! Fluid domain
                               nsize_flu, u_flu, appimesh, appidump)
        call read_wavefield_1d(dir_sim, .true., 'straintrace_sol', & ! Solid domain
                               nsize_sol, u_sol, appimesh, appidump)
  
        u1(1:nsize_flu) = u_flu(:)
        u1(nsize_flu+1:nsize) = u_sol(:)
  
        if (save_snaps .and. mod(it,10) == 0) then 
           filename = 'Data/sem_Eii_src_'//appimesh//'_'//appidump
           call write_vtk_bin_scal_mesh(u1, mesh1d(:,meshes_fwd(imesh),:), nsize, filename)
        endif
  
        ! apply mapping to the kernel mesh
        do ipt = 1, npts
            if(iproc_d_src(ipt).eq.meshes_fwd(imesh)) then
                u(ipt) = u1(ipt_d_src(ipt))
            end if
        end do
    enddo

    if (any(u==-123456)) then
        write(6,*) count(u==-123456), ' element was not filled in read_summedstraintrace_src'
        stop
    end if
end subroutine read_summedstraintrace_src
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_summedstraintrace_rec(u, appidump, dir_sim, it)

    use data_rota, only: iproc_d_rec, ipt_d_rec 
    include 'mesh_params.h'
    include 'mesh_params_kernel.h'
  
    character(len=4), intent(in)              :: appidump
    character(len=200), intent(in)            :: dir_sim
    integer, intent(in)                       :: it
  
    real(kind=realkind), dimension(npts), intent(out) :: u

    real(kind=realkind), dimension(nsize_sol) :: u_sol
    real(kind=realkind), dimension(nsize_flu) :: u_flu
    real(kind=realkind), dimension(nsize)     :: u1
    integer                                   :: ipt,imesh
    character(len=4)                          :: appimesh
    character(len=100)                        :: filename
  
    u = -123456
    do imesh=1, num_meshes
        if(.not.(any(iproc_d_rec==meshes(imesh)))) cycle
        
        call define_io_appendix(appimesh, meshes(imesh))
        call read_wavefield_1d(dir_sim, .false., 'straintrace_flu', & ! Fluid domain
                               nsize_flu, u_flu, appimesh, appidump)
        call read_wavefield_1d(dir_sim, .false., 'straintrace_sol', & ! Solid domain
                               nsize_sol, u_sol, appimesh, appidump)
  
        u1(1:nsize_flu) = u_flu(:) 
        u1(nsize_flu+1:nsize) = u_sol(:)
  
        if (save_snaps .and. mod(it,10) == 0) then 
           filename='Data/sem_Eii_rec_'//appimesh//'_'//appidump
           call write_vtk_bin_scal_mesh(u1, mesh1d(:,meshes(imesh),:), nsize, filename)
        endif
  
        do ipt = 1, npts
            if(iproc_d_rec(ipt).eq.meshes(imesh)) then
                u(ipt) = u1(ipt_d_rec(ipt))
            end if
        end do
        !! apply mapping to the kernel mesh
        !do ipt = 1, npts_iproc(imesh)
        !   u(iptkern_proc_rec(ipt,imesh)) = u1(iptsem_proc_rec(ipt,imesh))
        !enddo
    end do
    if (any(u==-123456)) then
        write(6,*) count(u==-123456), ' element was not filled in read_summedstraintrace_rec'
        stop
    end if

end subroutine read_summedstraintrace_rec
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_velocities_src(u, appidump, dir_sim)

    use data_rota, only: iproc_d_src, ipt_d_src 
    include 'mesh_params.h'
    include 'mesh_params_kernel.h'
  
    character(len=200), intent(in)    :: dir_sim
    character(len=4), intent(in)      :: appidump
  
    real(kind=realkind), dimension(1:npts,1:dim_fwd), intent(out) :: u
  
    real(kind=realkind), dimension(nsize_sol,1:dim_fwd) :: g
    real(kind=realkind), dimension(nsize_flu,1:dim_fwd) :: h
    real(kind=realkind), dimension(1:nsize,1:dim_fwd)   :: u1
    integer                           :: ipt, imesh
    character(len=4)                  :: appimesh
 
    u = -123456
    do imesh=1, num_meshes_fwd
  
        if(.not.(any(iproc_d_src==meshes_fwd(imesh)))) cycle
        call define_io_appendix(appimesh,meshes_fwd(imesh))
  
        ! Fluid domain
        call read_wavefields_1d(dir_sim, .true., 'velo_flu', &
                                nsize_flu, dim_fwd, h, appimesh, appidump)
  
        ! Solid domain
        call read_wavefields_1d(dir_sim, .true., 'velo_sol', &
                                nsize_sol, dim_fwd, g, appimesh, appidump)
  
        u1(1:nsize_flu,1:dim_fwd) = h(:,1:dim_fwd)
        u1(nsize_flu+1:nsize,1:dim_fwd) = g(:,1:dim_fwd)
  
        ! apply mapping to the kernel mesh
  
        do ipt = 1, npts
            if(iproc_d_src(ipt).eq.meshes_fwd(imesh)) then
                u(ipt,1:dim_fwd) = u1(ipt_d_src(ipt),1:dim_fwd)
            end if
        end do
    enddo
  
    if (any(u==-123456)) then
        write(6,*) count(u==-123456), ' element was not filled in read_velocities_src'
        stop
    end if
end subroutine read_velocities_src
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_velocities_rec(u, appidump, dir_sim)

    use data_rota, only: iproc_d_rec, ipt_d_rec 
    include 'mesh_params.h'
    include 'mesh_params_kernel.h'
  
    character(len=200), intent(in)    :: dir_sim
    character(len=4), intent(in)      :: appidump
  
    real(kind=realkind), dimension(1:npts,1:dim_bwd), intent(out) :: u
  
    real(kind=realkind), dimension(1:nsize_sol,1:dim_bwd) :: g
    real(kind=realkind), dimension(1:nsize_flu,1:dim_bwd) :: h
    real(kind=realkind), dimension(1:nsize,1:dim_bwd)     :: u1
    integer                           :: ipt, iproc
    character(len=4)                  :: appiproc
 
    u = -123456
    do iproc = 1, num_meshes
        call define_io_appendix(appiproc,meshes(iproc))
  
        ! Fluid domain
        call read_wavefields_1d(dir_sim, .false., 'velo_flu',  nsize_flu, dim_bwd, h, &
                                appiproc, appidump)
  
        ! Solid domain
        call read_wavefields_1d(dir_sim, .false., 'velo_sol',  nsize_sol, dim_bwd, g, &
                                appiproc, appidump)
  
        u1(1:nsize_flu,1:dim_bwd) = h(:,1:dim_bwd)
        u1(nsize_flu+1:nsize,1:dim_bwd) = g(:,1:dim_bwd)
  
        ! apply mapping to the kernel mesh
        do ipt = 1, npts
            if(iproc_d_rec(ipt).eq.meshes(iproc)) then
                u(ipt,1:dim_bwd) = u1(ipt_d_rec(ipt),1:dim_bwd)
            end if
        end do
    enddo

    if (any(u==-123456)) then
        write(6,*) count(u==-123456), ' element was not filled in read_velocities_rec'
        stop
    end if
end subroutine read_velocities_rec
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_deviatoricstrain_src(u, appidump, dir_sim, it)

  use data_rota, only: iproc_d_src, ipt_d_src 
  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  character(len=200), intent(in)        :: dir_sim
  character(len=4), intent(in)          :: appidump
  integer, intent(in)                   :: it

  real(kind=realkind), dimension(1:npts,1:dim_fwd*2), intent(out) :: u
  
  real(kind=realkind), dimension(1:nsize_sol)           :: g
  real(kind=realkind), dimension(1:nsize_flu)           :: h
  real(kind=realkind), dimension(1:nsize,1:dim_fwd*2)   :: u1
  integer                               :: ipt, iproc
  character(len=4)                      :: appiproc
  character(len=100)                    :: filename

    u = -123456
  do iproc=1, num_meshes
     if(.not.(any(iproc_d_src==meshes(iproc)))) cycle
     call define_io_appendix(appiproc,meshes(iproc))

     ! Fluid domain diagonal
     call read_wavefield_1d(dir_sim, .true., 'strain_dsus_flu',  nsize_flu, h, appiproc, &
                            appidump)
     u1(1:nsize_flu,1) = h
     call read_wavefield_1d(dir_sim, .true., 'strain_dpup_flu',  nsize_flu, h, appiproc, &
                            appidump)
     u1(1:nsize_flu,2) = h
     call read_wavefield_1d(dir_sim, .true., 'straintrace_flu', nsize_flu, h, appiproc, &
                            appidump)

     u1(1:nsize_flu,1) = u1(1:nsize_flu,1) - h / 3.
     u1(1:nsize_flu,2) = u1(1:nsize_flu,2) - h / 3.
     u1(1:nsize_flu,3) = - u1(1:nsize_flu,1) - u1(1:nsize_flu,2)

     ! Fluid domain off-diagonal
     call read_wavefield_1d(dir_sim, .true., 'strain_dsuz_flu',  nsize_flu, h, appiproc, &
                            appidump)
     u1(1:nsize_flu,4) = h
     if (dim_fwd == 3) then 
        call read_wavefield_1d(dir_sim, .true., 'strain_dsup_flu',  nsize_flu, h, &
                               appiproc, appidump)
        u1(1:nsize_flu,5) = h
        call read_wavefield_1d(dir_sim, .true., 'strain_dzup_flu',  nsize_flu, h, &
                               appiproc, appidump)
        u1(1:nsize_flu,6) = h
     endif

     ! Solid domain diagonal
     call read_wavefield_1d(dir_sim, .true., 'strain_dsus_sol', nsize_sol, g, appiproc, &
                            appidump)
     u1(nsize_flu+1:nsize,1) = g
     call read_wavefield_1d(dir_sim, .true., 'strain_dpup_sol', nsize_sol, g, appiproc, &
                            appidump)
     u1(nsize_flu+1:nsize,2) = g
     call read_wavefield_1d(dir_sim, .true. , 'straintrace_sol',  nsize_sol, g, &
                            appiproc, appidump)

     u1(nsize_flu+1:nsize,1) = u1(nsize_flu+1:nsize,1) - g / 3.
     u1(nsize_flu+1:nsize,2) = u1(nsize_flu+1:nsize,2) - g / 3.
     u1(nsize_flu+1:nsize,3) = - u1(nsize_flu+1:nsize,1) - u1(nsize_flu+1:nsize,2) 

     ! Solid domain off-diagonal
     call read_wavefield_1d(dir_sim, .true., 'strain_dsuz_sol', nsize_sol, g, appiproc, &
                            appidump)
     u1( nsize_flu+1:nsize,4) = g
     if (dim_fwd == 3) then
        call read_wavefield_1d(dir_sim, .true., 'strain_dsup_sol', nsize_sol, g, &
                               appiproc, appidump)
        u1(nsize_flu+1:nsize,5) = g
        call read_wavefield_1d(dir_sim, .true., 'strain_dzup_sol', nsize_sol, g, &
                               appiproc, appidump)
        u1(nsize_flu+1:nsize,6) = g
     endif

     if (save_snaps .and. mod(it,20) == 0) then 
        filename='Data/sem_Esz_src_'//appiproc//'_'//appidump
        call write_vtk_bin_scal_mesh(u1(:,4), mesh1d(:,meshes_fwd(iproc),:), nsize, &
                                     filename)
     endif

    ! ! apply mapping to the kernel mesh
    ! do ipt=1, npts_iproc(iproc)
    !    u(iptkern_proc_src(ipt,iproc),1:dim_fwd*2) = &
    !                                        u1(iptsem_proc_src(ipt,iproc),1:dim_fwd*2)
    ! enddo
      do ipt = 1, npts
          if(iproc_d_src(ipt).eq.meshes_fwd(iproc)) then
              u(ipt,1:dim_fwd*2) = u1(ipt_d_src(ipt),1:dim_fwd*2)
          end if
      end do
  enddo


    if (any(u==-123456)) then
        write(6,*) count(u==-123456), ' element was not filled in read_deviatoricstrain_src'
        stop
    end if
end subroutine read_deviatoricstrain_src
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_deviatoricstrain_rec(u, appidump, dir_sim, it)

  use data_rota, only: iproc_d_rec, ipt_d_rec 
  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  character(len=200), intent(in)    :: dir_sim
  character(len=4), intent(in)      :: appidump
  integer, intent(in)               :: it

  real(kind=realkind), dimension(1:npts,1:dim_bwd*2), intent(out) :: u
  
  real(kind=realkind), dimension(1:nsize_sol)           :: g
  real(kind=realkind), dimension(1:nsize_flu)           :: h
  real(kind=realkind), dimension(1:nsize,1:dim_bwd*2)   :: u1
  integer                           :: ipt, iproc
  character(len=4)                  :: appiproc
  character(len=100)                :: filename

    u = -123456
  do iproc=1, num_meshes
     call define_io_appendix(appiproc, meshes(iproc))

     ! Fluid domain diagonal
     call read_wavefield_1d(dir_sim, .false., 'strain_dsus_flu', nsize_flu, h, &
                            appiproc, appidump)
     u1(1:nsize_flu,1) = h

     call read_wavefield_1d(dir_sim, .false., 'strain_dpup_flu', nsize_flu, h, &
                            appiproc, appidump)
     u1(1:nsize_flu,2) = h
     ! NEED TO CHECK ALL THIS ONCE MORE!!!
     u1(1:nsize_flu,3) = - u1(1:nsize_flu,1) - u1(1:nsize_flu,2)

     call read_wavefield_1d(dir_sim, .false., 'straintrace_flu', nsize_flu, h, &
                            appiproc, appidump)
     u1(1:nsize_flu,1) = u1(1:nsize_flu,1) - h / 3.
     u1(1:nsize_flu,2) = u1(1:nsize_flu,2) - h / 3.

     ! Fluid domain off-diagonal
     call read_wavefield_1d(dir_sim, .false., 'strain_dsuz_flu', nsize_flu, h, &
                            appiproc, appidump)
     u1(1:nsize_flu,4) = h
     if (dim_bwd == 3) then 
        call read_wavefield_1d(dir_sim, .false., 'strain_dsup_flu', nsize_flu, h, &
                               appiproc,appidump)
        u1(1:nsize_flu,5) = h
        call read_wavefield_1d(dir_sim, .false., 'strain_dzup_flu', nsize_flu, h, &
                               appiproc, appidump)
        u1(1:nsize_flu,6) = h
     endif

     ! Solid domain diagonal
     call read_wavefield_1d(dir_sim, .false., 'strain_dsus_sol', nsize_sol, g, &
                            appiproc, appidump)
     u1( nsize_flu+1:nsize,1) = g

     call read_wavefield_1d(dir_sim, .false., 'strain_dpup_sol', nsize_sol, g, &
                            appiproc, appidump)
     u1( nsize_flu+1:nsize,2) = g
     ! NEED TO CHECK
     u1( nsize_flu+1:nsize,3) = - u1(nsize_flu+1:nsize,1) - u1(nsize_flu+1:nsize,2)
   
     call read_wavefield_1d(dir_sim, .false., 'straintrace_sol', nsize_sol, g, &
                            appiproc, appidump)
     u1(nsize_flu+1:nsize,1) = u1(nsize_flu+1:nsize,1) - g / 3.
     u1(nsize_flu+1:nsize,2) = u1(nsize_flu+1:nsize,2) - g / 3.

     ! Solid domain off-diagonal
     call read_wavefield_1d(dir_sim, .false., 'strain_dsuz_sol', nsize_sol, g, &
                            appiproc, appidump)
     u1(nsize_flu+1:nsize,4) = g

     if (dim_bwd == 3) then
        call read_wavefield_1d(dir_sim, .false., 'strain_dsup_sol', nsize_sol, g, &
                               appiproc, appidump)
        u1(nsize_flu+1:nsize,5) = g

        call read_wavefield_1d(dir_sim, .false., 'strain_dzup_sol', nsize_sol, g, &
                               appiproc, appidump)
        u1(nsize_flu+1:nsize,6) = g
     endif

     if (save_snaps .and. mod(it,20) == 0) then 
        filename = 'Data/sem_Esz_rec_'//appiproc//'_'//appidump
        call write_vtk_bin_scal_mesh(u1(:,4), mesh1d(:,meshes(iproc),:), nsize, filename)
     endif

    ! ! apply mapping to the kernel mesh
    ! do ipt = 1,npts_iproc(iproc)
    !    u(iptkern_proc_rec(ipt,iproc),1:dim_bwd*2) = &
    !                                        u1(iptsem_proc_rec(ipt,iproc),1:dim_bwd*2)
    ! enddo
      do ipt = 1, npts
          if(iproc_d_rec(ipt).eq.meshes(iproc)) then
              u(ipt,1:dim_bwd*2) = u1(ipt_d_rec(ipt),1:dim_bwd*2)
          end if
      end do
  enddo

    if (any(u==-123456)) then
        write(6,*) count(u==-123456), ' element was not filled in read_deviatoricstrain_rec'
        stop
    end if

end subroutine read_deviatoricstrain_rec
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_strain_src(u, appidump, dir_sim)

  use data_rota, only : iproc_srckern, ipt_srckern

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  character(len=200), intent(in)    :: dir_sim
  character(len=4), intent(in)      :: appidump

  real(kind=realkind), dimension(1:n_src,1:dim_bwd*2), intent(out) :: u

  real(kind=realkind), dimension(1:nsize_sol)           :: g
  real(kind=realkind), dimension(1:nsize_flu)           :: h
  real(kind=realkind), dimension(1:nsize,1:dim_bwd*2)   :: u1
  integer                           :: ipt, iproc
  character(len=4)                  :: appiproc
  character(len=100)                :: filename

  do ipt=1, n_src
     iproc = iproc_srckern(ipt)
     call define_io_appendix(appiproc, iproc)

     ! Fluid domain
     call read_wavefield_1d(dir_sim, .false., 'strain_dsus_flu', nsize_flu, h, &
                            appiproc, appidump)
     u1(1:nsize_flu,1) = h

     call read_wavefield_1d(dir_sim, .false., 'strain_dpup_flu', nsize_flu, h, &
                            appiproc, appidump)
     u1(1:nsize_flu,2) = h

     call read_wavefield_1d(dir_sim, .false., 'straintrace_flu', nsize_flu, h, &
                            appiproc, appidump)
     ! NEED TO CHECK ALL THIS ONCE MORE!!!
     u1(1:nsize_flu,3)= h - u1(1:nsize_flu,1) - u1(1:nsize_flu,2) 

     call read_wavefield_1d(dir_sim, .false., 'strain_dsuz_flu', nsize_flu, h, &
                            appiproc, appidump)
     u1(1:nsize_flu,4) = h

     if (dim_bwd == 3) then 
        call read_wavefield_1d(dir_sim, .false., 'strain_dsup_flu', nsize_flu, h, &
                               appiproc, appidump)
        u1(1:nsize_flu,5) = h

        call read_wavefield_1d(dir_sim, .false., 'strain_dzup_flu', nsize_flu, h, &
                               appiproc, appidump)
        u1(1:nsize_flu,6) = h
     endif

     ! Solid domain
     call read_wavefield_1d(dir_sim, .false., 'strain_dsus_sol', nsize_sol, g, &
                            appiproc, appidump)
     u1( nsize_flu+1:nsize,1) = g

     call read_wavefield_1d(dir_sim, .false., 'strain_dpup_sol', nsize_sol, g, &
                            appiproc, appidump)
     u1( nsize_flu+1:nsize,2) = g

     call read_wavefield_1d(dir_sim, .false., 'straintrace_sol', nsize_sol, g, &
                            appiproc, appidump)
     ! NEED TO CHECK ALL THIS ONCE MORE!!!
     u1( nsize_flu+1:nsize,3) = g - u1(nsize_flu+1:nsize,1) - u1(nsize_flu+1:nsize,2)

     call read_wavefield_1d(dir_sim, .false., 'strain_dsuz_sol', nsize_sol, g, &
                            appiproc, appidump)
     u1( nsize_flu+1:nsize,4) = g

     if (dim_bwd == 3) then
        call read_wavefield_1d(dir_sim, .false., 'strain_dsup_sol', nsize_sol, g, &
                               appiproc, appidump)
        u1(nsize_flu+1:nsize,5) = g

        call read_wavefield_1d(dir_sim, .false., 'strain_dzup_sol', nsize_sol, g, &
                               appiproc, appidump)
        u1(nsize_flu+1:nsize,6) = g
     endif

     u(ipt,1:dim_bwd*2) = u1(ipt_srckern(ipt),1:dim_bwd*2)
  enddo

end subroutine read_strain_src
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
#ifdef unc
subroutine close_netcdf
  use netcdf
  use data_mesh,only: ncid_fw_in, ncid_bw_in

  call check( nf90_close(ncid=ncid_fw_in) )
  call check( nf90_close(ncid=ncid_bw_in) )

end subroutine
#endif
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
#ifdef unc
subroutine open_netcdf

  use netcdf
  use data_mesh, only: ncid_fw_in, ncid_bw_in, dir_fwdmesh, dir_bwdmesh

  integer                   :: ivar, status 
  character(len=200)        :: format20, format21

  format20 = "('Trying to open NetCDF file ', A, ' on CPU ', I5)"
  format21 = "('Succeded with  NetCDF file ', A, ' has NCID ', I6, "// &
             "' Snapshots group NCID:', I6)"

  varnamelist = (/'strain_dsus_sol', 'strain_dsuz_sol', 'strain_dpup_sol', &
                  'straintrace_sol', 'velo_sol       ', 'strain_dsus_flu', &
                  'strain_dsuz_flu', 'strain_dpup_flu', 'straintrace_flu', &
                  'velo_flu       '/)

  ! Forward wavefield
  write(6,format20) trim(dir_fwdmesh)//'Data/axisem_output.nc4', mynum 
  call check( nf90_open(path=trim(dir_fwdmesh)//'Data/axisem_output.nc4', &
                        mode=nf90_nowrite, ncid=ncid_fw_in) )

  ! The following block tries to find the variable and group ids from the 
  ! netcdf output file of the SOLVER. This is done only for the forward
  ! file. However, there is absolutely no legitimate reason, why they should
  ! be different in the backward file, unless it was produced with a different
  ! revision of Axisem. In this case, you get what you deserve.
  status =  nf90_inq_ncid(ncid_fw_in, name="Snapshots", grp_ncid=ncid_snapshot_fw)
  do ivar = 1, nvar
    status =  nf90_inq_varid(ncid_snapshot_fw, trim(varnamelist(ivar)), nc_varid(ivar))
    if (mynum.eq.0) write(6,*) 'fw:', ivar, trim(varnamelist(ivar)), nc_varid(ivar)
    if (status /= nf90_NoErr) then
      print *,'Error in wavefield file, Variable ', trim(varnamelist(ivar)),' not found.'
      stop
    end if
  end do

  write(6,format21) trim(dir_fwdmesh)//'Data/axisem_output.nc4', ncid_fw_in, &
                    ncid_snapshot_fw
  
  ! Backward wavefield
  write(6,format20) trim(dir_bwdmesh)//'Data/axisem_output.nc4', mynum 
  call check( nf90_open(path=trim(dir_bwdmesh)//'Data/axisem_output.nc4', &
                        mode=nf90_nowrite, ncid=ncid_bw_in) )
  status =  nf90_inq_ncid(ncid_bw_in, name="Snapshots", grp_ncid=ncid_snapshot_bw)

  do ivar = 1, nvar
    status =  nf90_inq_varid(ncid_snapshot_bw, trim(varnamelist(ivar)), nc_varid(ivar))
    if (mynum.eq.0) write(6,*) 'bw:',ivar, trim(varnamelist(ivar)), nc_varid(ivar)
    if (status /= nf90_NoErr) then
      print *,'Error in wavefield file, Variable ', trim(varnamelist(ivar)),' not found.'
      stop
    end if
  end do

  write(6,format21) trim(dir_bwdmesh)//'Data/axisem_output.nc4', ncid_bw_in, &
                    ncid_snapshot_bw
  call flush(6) 

end subroutine
#endif
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
#ifdef unc
subroutine check(status)
! Translates netcdf error codes into error messages

  use netcdf

  implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop 2
  end if
end subroutine check  
#endif
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_wavefield_1d(datapath, lforward, filename, n, field, appiproc, appidump)

#ifdef unc
  use netcdf
#endif
  use data_mesh, only: ncid_fw_in, ncid_bw_in, use_netcdf

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  character(len=200), intent(in)        :: datapath
  integer, intent(in)                   :: n
  logical, intent(in)                   :: lforward   ! Are we loading the forward wavefield (true)?
  character(len=15), intent(in)         :: filename
  character(len=4), intent(in)          :: appiproc,appidump

  real(kind=realkind), dimension(1:n), intent(out) :: field

  integer                               :: idump_loc, iproc_loc 
  integer                               :: ncid_in, ivar

#ifdef unc
  if (use_netcdf) then
     if (lforward) then
       ncid_in = ncid_snapshot_fw
     else
       ncid_in = ncid_snapshot_bw
     end if
     
     do ivar=1, nvar
         ! Check whether this Variable actually exists in file ncid_out
         if (trim(varnamelist(ivar)) == trim(filename)) exit
     end do
     if (ivar > nvar) then
         write(6,*) 'read_wavefield_1d: Trying to access variable: ', trim(filename), &
                    ' which is not in varnamelist. Contact a developer and shout at him!'
         stop 1
     end if

     read(appidump,*) idump_loc
     read(appiproc,*) iproc_loc

     call check( nf90_get_var(ncid_in, nc_varid(ivar), &
                              start=(/1, iproc_loc+1, idump_loc/), &
                              count=(/n, 1, 1/), values=field) )
  end if
#endif

  if (.not. use_netcdf) then

     open(unit=10000, file=trim(datapath)//'/Data/'//filename//'_' &
                               //appiproc//'_'//appidump//'.bindat', &
          FORM="UNFORMATTED", ACTION='read', POSITION="REWIND")
     read(10000) field
     close(10000)
  end if

end subroutine read_wavefield_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_wavefields_1d(datapath, lforward, filename, n, dim2, field, appiproc, &
                              appidump)

#ifdef unc
  use netcdf
#endif
  use data_mesh, only: ncid_fw_in, ncid_bw_in, use_netcdf

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  character(len=200), intent(in)            :: datapath
  logical, intent(in)                       :: lforward
  integer, intent(in)                       :: n, dim2
  character(len=8), intent(in)              :: filename
  character(len=4), intent(in)              :: appiproc, appidump

  real(kind=realkind), dimension(1:n,1:dim2), intent(out) :: field

  integer                                   :: ierr, idump_loc, iproc_loc
  integer                                   :: status, ncid_in, ivar

#ifdef unc
  if (use_netcdf) then
     if (lforward) then
       ncid_in = ncid_snapshot_fw
     else
       ncid_in = ncid_snapshot_bw
     end if

     do ivar = 1, nvar
         ! Check whether this Variable actually exists in file ncid_out
         if (trim(varnamelist(ivar)) == trim(filename)) exit
     end do
     if (ivar>nvar) then
         write(6,*) 'read_wavefields_1d: Trying to access variable: ', trim(filename), &
                    ' which is not in varnamelist. Contact a developer and shout at him!'
         stop 1
     end if
     
     read(appidump,*) idump_loc
     read(appiproc,*) iproc_loc

     if (status /= nf90_NoErr) then
       print *,'Error in wavefield file, Variable ', filename,' not found.'
       stop
     end if

     call check( nf90_get_var(ncid_in, nc_varid(ivar), &
                              start=(/1, iproc_loc+1, idump_loc/), &
                              count=(/n*dim2, 1, 1/), values=field) )
  end if
#endif

  if (.not. use_netcdf) then
     open(unit=100000, file=trim(datapath)//'/Data/'//filename//'_' &
                            //appiproc//'_'//appidump//'.bindat', &
          FORM="UNFORMATTED", ACTION='read', status='old', POSITION="REWIND")

     read(100000) field
     close(100000)
  end if

end subroutine read_wavefields_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine save_kernel(n, kern, kerntype, appidump1, globmax)

  use data_arrays, only : kern_cell
  !use mpi

  include 'mpif.h'

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  integer, intent(in)               :: n
  real,dimension(1:n),intent(in)    :: kern
  character(len=20),intent(in)      :: kerntype
  character(len=4),intent(in)       :: appidump1
  real, intent(in), optional        :: globmax

  real,dimension(1:n)               :: kernloc
  integer                           :: ipt
  real                              :: globmaxloc
  character(len=200)                :: filename

  if (.not. present(globmax)) then 
     globmaxloc = 1.
  else
     if (globmax < epsi) then 
        globmaxloc = 1.
     else
        globmaxloc = globmax
     endif
  endif
  kernloc = kern / globmaxloc

  filename = 'Data/'//trim(kerntype)//'_'//appmynum//'_'//appidump1

  if (dump_avs) call write_avs_file_scal(kernloc, n, filename)
  if (dump_bin) call write_bin_scal(kernloc, n, filename)
  if (dump_ascii) call write_ascii_scal(kernloc, n, filename)
  if (dump_vtk) then 
     if (cell_topology) then 
        call field2cell(npts, kernloc, npts_cell, kern_cell)
        call write_VTK_bin_scal_topology(xgd_cell, ygd_cell, zgd_cell, kern_cell, &
                                         npts_cell/4, filename)
     else 
        call write_VTK_bin_scal(kernloc, n, filename)
     endif
  endif

end subroutine save_kernel
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine save_devstrain(npts, devstrain, fname, ndim1, appidump)

  use data_arrays, only : kern_cell
  
  implicit none
  
  real(kind=realkind), dimension(npts,ndim1),intent(in) :: devstrain
  integer,intent(in)                   :: ndim1, npts
  character(len=4),intent(in)          :: appidump
  character(len=200),intent(in)        :: fname

  integer                              :: i
  character(len=200)                   :: fname1,filename
  character(len=3)                     :: fname2(ndim1)
  real                                 :: shiftx(ndim1), shiftz(ndim1), dshiftx, dshiftz
  
  fname2 = ['Exx', 'Eyy', 'Ezz', 'Exz']
  if (ndim1 == 6) then
     fname2(5) = 'Exy'
     fname2(6) = 'Eyz'
  endif
  
  !MvD: what is going on here? XXX
  if (cell_topology) then 
     shiftx = 0.
     shiftz = 0.
     dshiftx = 1.2 * (maxval(xgd_cell) - minval(xgd_cell))
     dshiftz = 1.2 * (maxval(zgd_cell) - minval(zgd_cell))
     if (ndim1 == 4) then 
        shiftx(2) = dshiftx
        shiftz(3:4) = -dshiftz
        shiftx(4) = dshiftx
     else
        shiftx(2) = dshiftx
        shiftx(3) = 2. * dshiftx
        shiftz(4:6) = -dshiftz
        shiftx(5) = dshiftx
        shiftx(6) = 2. * dshiftx
     endif
  endif
  
  do i=1, ndim1
     fname1 = trim(fname2(i))//trim(fname)
     filename = 'Data/'//trim(fname1)//'_'//appmynum//'_'//appidump
    if (dump_vtk) then 
       if (cell_topology) then 
          call field2cell(npts, devstrain(1:npts,i), npts_cell, kern_cell)
          call write_VTK_bin_scal_topology(xgd_cell+shiftx(i), ygd_cell, &
                                           zgd_cell+shiftz(i), kern_cell, npts_cell/4, &
                                           filename)
       else 
          call write_VTK_bin_scal(devstrain(1:npts,i), npts, filename)
       endif
     endif
  enddo

end subroutine save_devstrain
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine field2cell(npts, field, npts_cell, field_cell)

  integer, intent(in)                                       :: npts,npts_cell
  real(kind=realkind), dimension(npts), intent(in)          :: field
  real(kind=realkind), dimension(npts_cell), intent(out)    :: field_cell
  real(kind=realkind), dimension(npts)                      :: buf
  integer                                                   :: i, j, iproc
  
  field_cell = -1.e30
  do i=1, npts
     do j=1, valence_tot2cell(i)
        field_cell( ind_tot2cell(i,j) ) = field(i)
     enddo
  enddo
  
  !if (nproc>1) then
  !   call MPI_SEND(buf,npts_cell,mpi_real,0,mynum,MPI_COMM_WORLD,IERROR)
  !   if (mynum==0) then 
  !     do iproc=1,nproc-2
  !         call MPI_RECV(buf,npts_cell,mpi_real,iproc,iproc,MPI_COMM_WORLD,IERROR)
  !         field_cell( iproc*npts+1:iproc*npts+npts ) = buf
  !      enddo
  !   endif
  !endif
  
  !if (nproc>1) &
  ! call MPI_REDUCE(field_cell,field_cell,npts_cell,MPI_REAL,MPI_SUM, 0,&
  !                    MPI_COMM_WORLD,IERROR)
  
  if (minval(field_cell) < -1e25) then 
     write(6,*) mynum, 'something went wrong in the mapping from global to cell!', &
                minval(field_cell) 
     stop 2
  endif

end subroutine field2cell
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine save_misfit_kernel(n, kern1, kerntype, appidump1)

  !use mpi
  include 'mpif.h'
  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  integer, intent(in)                     :: n
  real,dimension(1:n), intent(in)         :: kern1
  character(len=30), intent(in)           :: kerntype
  character(len=4), intent(in), optional  :: appidump1

  real,dimension(1:n)                     :: kern
  integer                                 :: ipt
  real                                    :: globmax
  character(len=100)                      :: filename

  kern = kern1

  if (norm_kernels) then
     globmax = maxval(abs(kern))
     if (globmax < epsi) globmax = 1.
     call MPI_ALLREDUCE(globmax, globmax, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, IERROR) 
     kern = kern / globmax
  endif

  if (present(appidump1)) then 
     filename = 'Data/'//trim(kerntype)//'_kern_'//appmynum//'_'//appidump1
  else
     filename = 'Data/'//trim(kerntype)//'_kern_'//appmynum
  endif

  if (dump_vtk) call write_VTK_bin_scal(kern, n, filename)
  if (dump_avs) call write_avs_file_scal(kern, n, filename)
  if (dump_bin) call write_bin_scal(kern, n, filename)
  if (dump_ascii) call write_ascii_scal(kern, n, filename)

end subroutine save_misfit_kernel
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_ascii_scal(kern,n,filename)

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  integer, intent(in)               :: n
  real,dimension(1:n),intent(in)    :: kern
  character(len=100), intent(in)    :: filename

  integer :: ipt

  open(unit=150+mynum, file=trim(filename)//'.ascii')

  do ipt=1, n
     write(150+mynum,*) kern(ipt)
  enddo

  close(150+mynum)

end subroutine write_ascii_scal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_bin_scal(kern, n, filename)

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  integer, intent(in)               :: n
  real,dimension(1:n),intent(in)    :: kern
  character(len=100), intent(in)    :: filename

  integer                           :: ipt

  open(unit=150+mynum, file=trim(filename)//'.bin', FORM="UNFORMATTED", &
       STATUS="UNKNOWN", POSITION="REWIND")
  write(150+mynum) kern
  close(150+mynum)

end subroutine write_bin_scal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine save_snapshot(n, usrc, urec, appidump1, snaptype)

  implicit none

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  integer, intent(in)               :: n
  real,dimension(1:n),intent(in)    :: usrc
  real,dimension(1:n),intent(in)    :: urec
  character(len=4),intent(in)       :: appidump1
  character(len=3),intent(in)       :: snaptype

  integer :: ipt

  open(unit=250,file='Data/'//'snapshot_'//snaptype//'_' &
                     //appmynum//'_'//appidump1//'.dat')
  do ipt=1, n
     write(250,*) real(usrc(ipt)), real(urec(ipt))
  enddo
  close(250)

end subroutine save_snapshot
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump2d(fname, npts, r2d, theta2d, f2d)

  implicit none

  character(len=200)                 :: fname
  integer, intent(in)                :: npts
  real, intent(in),dimension(1:npts) :: r2d
  real, intent(in),dimension(1:npts) :: theta2d
  real, intent(in),dimension(1:npts) :: f2d

  integer :: ipt

  open(10,file=trim(fname), STATUS="UNKNOWN")
  do ipt = 1, npts
     write(10,"(3(1pe12.5,2x))") r2d(ipt), theta2d(ipt), f2d(ipt)
  end do

  close(10)

end subroutine dump2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump1d(fname, npts, r2d, theta2d, f2d)

  implicit none

  character(len=200)                 :: fname
  integer, intent(in)                :: npts
  real, intent(in),dimension(1:npts) :: r2d
  real, intent(in),dimension(1:npts) :: theta2d
  real, intent(in),dimension(1:npts) :: f2d

  integer                            :: ipt

  open(10, file=trim(fname), STATUS="UNKNOWN")
  do ipt = 1, npts
     write(10,"(3(1pe12.5,2x))") r2d(ipt),f2d(ipt)
  end do

  close(10)

end subroutine dump1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_VTK_bin_scal(u2, rows, filename)

  implicit none

  integer*4, intent(in)                 :: rows
  real, dimension(1:rows), intent(in)   :: u2
  character (len=55)                    :: filename

  integer*4                             :: i, t, ioerr, dims
  real, dimension(1:rows)               :: u1
  integer*4, dimension(1:rows*2)        :: cell
  integer*4, dimension(1:rows)          :: cell_type
  character (len=50)                    :: ss
 
  !points structure
  do i=2, rows*2, 2
   cell(i-1) = 1
   cell(i) = (i/2) - 1
  enddo

  cell_type = 1
  
  u1 = real(u2)
  do i=1, rows
     if (abs(u1(i)) < 1.e-25) u1(i) = 0.
  enddo

  write(6789,*) size(u1), maxval(u1), minval(u1)
  
  open(100+mynum, file=trim(filename)//'.vtk', access='stream', status='replace', &
                convert='big_endian')
  
  write(100+mynum) '# vtk DataFile Version 4.0'//char(10)
  write(100+mynum) 'mittico'//char(10)
  write(100+mynum) 'BINARY'//char(10)
  write(100+mynum) 'DATASET UNSTRUCTURED_GRID'//char(10)
  write(ss,fmt='(A6,I10,A5)') 'POINTS',rows,'float'
  write(100+mynum) ss//char(10)

  !points
  do i=1, rows
     write(100+mynum) real(W_vtk(i,1)), real(W_vtk(i,2)), real(W_vtk(i,3))
  enddo
  write(100+mynum) char(10)

  !cell topology
  write(ss,fmt='(A5,2I10)') 'CELLS', rows, rows*2
  write(100+mynum) char(10)//ss//char(10)
  write(100+mynum) cell
  write(100+mynum) char(10)
  
  !cell type
  write(ss,fmt='(A10,2I10)') 'CELL_TYPES', rows
  write(100+mynum) char(10)//ss//char(10)
  write(100+mynum) cell_type
  write(100+mynum) char(10)

  !data
  write(ss,fmt='(A10,I10)') 'CELL_DATA',rows
  write(100+mynum) char(10)//ss//char(10)
  write(100+mynum) 'SCALARS Displ_u1 float 1'//char(10)
  write(100+mynum) 'LOOKUP_TABLE default'//char(10)
  write(100+mynum) real(u1)
  close(100+mynum)

end subroutine write_VTK_bin_scal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_VTK_bin_scal_mesh(u2, mesh, rows, filename)

  implicit none

  integer*4, intent(in)                     :: rows
  real, dimension(1:rows), intent(in)       :: u2
  real, dimension(1:rows,1:2), intent(in)   :: mesh
  character (len=55)                        :: filename

  integer*4                                 :: i, t, ioerr, dims
  real, dimension(1:rows)                   :: u1
  integer*4, dimension(1:rows*2)            :: cell
  integer*4, dimension(1:rows)              :: cell_type
  character (len=50)                        :: ss


  do i=2, rows*2, 2
     cell(i-1) = 1
     cell(i) = (i / 2) - 1;
  enddo

  do i=1, rows
     cell_type(i) = 1
  enddo

  u1 = real(u2)
  do i=1, rows
     if (abs(u1(i)) < 1.e-25) u1(i) = 0.001 * maxval(abs(u1))
  enddo

  open(1000+mynum, file=trim(filename)//'.vtk', access='stream', status='replace', &
        convert='big_endian')
  write(1000+mynum) '# vtk DataFile Version 4.0'//char(10)
  write(1000+mynum) 'mittico'//char(10)
  write(1000+mynum) 'BINARY'//char(10)
  write(1000+mynum) 'DATASET UNSTRUCTURED_GRID'//char(10)
  write(ss,fmt='(A6,I10,A5)') 'POINTS',rows,'float'
  write(1000+mynum) ss//char(10)

  !points
  do i=1, rows
     write(1000+mynum) real(mesh(i,1)), 0., real(mesh(i,2))
  enddo
  write(1000+mynum) char(10)

  !cell topology
  write(ss,fmt='(A5,2I10)') 'CELLS', rows, rows*2
  write(1000+mynum) char(10)//ss//char(10)
  write(1000+mynum) cell
  write(1000+mynum) char(10)

  !cell type
  write(ss,fmt='(A10,2I10)') 'CELL_TYPES', rows
  write(1000+mynum) char(10)//ss//char(10)
  write(1000+mynum) cell_type
  write(1000+mynum) char(10)

  !data
  write(ss,fmt='(A10,I10)') 'CELL_DATA', rows
  write(1000+mynum) char(10)//ss//char(10)
  write(1000+mynum) 'SCALARS Displ_u1 float 1'//char(10)
  write(1000+mynum) 'LOOKUP_TABLE default'//char(10)
  write(1000+mynum) real(u1)
  close(1000+mynum)

end subroutine write_VTK_bin_scal_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_VTK_bin_scal_mesh3d(x, y, z, u1, rows, nelem_disk, filename1)

  implicit none

  integer, intent(in)                   :: rows, nelem_disk
  real, dimension(1:rows), intent(in)   :: x,y,z,u1
  character (len=100), intent(in)       :: filename1

  integer                               :: t, ioerr, dims, i
  integer, dimension(1:rows,2)          :: cell
  integer, dimension(1:rows)            :: cell_type
  real, dimension(1:rows,3)             :: W
  
  character (len=30)                    :: celltype
  character (len=50)                    :: ss
  
  W(1:rows,1) = x
  W(1:rows,2) = y
  W(1:rows,3) = z

  !points structure
  do i=1, rows
     cell(i,1) = 1
     cell(i,2) = i
  enddo

  cell_type = 1
  
  write(6,*) 'computing VTK bin file ', trim(filename1)//'.vtk  ...'
  
  ! 1 IS WRONG FOR OUTDIR !!!!! 
  open(100, file=trim(filename1)//'.vtk', access='stream', &
                           status='replace', convert='big_endian')
  
  write(100) '# vtk DataFile Version 3.0'//char(10)
  write(100) 'Cell Fractions'//char(10)
  write(100) 'BINARY'//char(10)
  write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
  write(ss,fmt='(A8,I12,A10)') 'POINTS', rows, ' float'
  write(100) ss//char(10)

  !points
  do i=1, rows
     write(100) W(i,1:3)
  enddo
  write(100) char(10)

  !cell topology
  write(ss,fmt='(A5,2I12)') 'CELLS', rows, rows*2
  write(100) char(10)//ss//char(10)
  do i=1, rows
     write(100) cell(i,1:2)
  enddo
  write(100) char(10)
  
  !cell type
  write(ss,fmt='(A10,2I12)') 'CELL_TYPES', rows
  write(100) char(10)//ss//char(10)
  do i=1,rows
     write(100) cell_type(i)
  enddo
  write(100) char(10)
  
  !data
  write(ss,fmt='(A10,I12)') 'CELL_DATA',rows
  write(100) char(10)//ss//char(10)
  write(100) 'SCALARS Displ_u1 float 1'//char(10)
  write(100) 'LOOKUP_TABLE default'//char(10)
  do i=1, rows
     write(100) u1(i)
  enddo
  close(100)
  
  write(6,*)'...saved ',trim(filename1)//'.vtk'

end subroutine write_VTK_bin_scal_mesh3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_vtk_bin_scal_topology(x, y, z, u1, elems, filename)

  implicit none

  integer*4, intent(in)                     :: elems
  real*4, dimension(1:elems*4), intent(in)  :: x, y, z, u1
  character (len=200)                       :: filename

  integer*4                                 :: i, t
  integer*4, dimension(1:elems*5)           :: cell
  integer*4, dimension(1:elems)             :: cell_type
  character (len=50)                        :: ss 

  !points structure
  t = 0
  do i=5, elems*5, 5
    t = t + 4
    cell(i-4) = 4
    cell(i-3) = t - 4
    cell(i-2) = t - 3
    cell(i-1) = t - 2
    cell(i) = t - 1
  enddo
  
  cell_type = 9
  
  open(100,file=trim(filename)//'.vtk', access='stream', status='replace', &
       convert='big_endian')
  write(100) '# vtk DataFile Version 4.0'//char(10)
  write(100) 'mittico'//char(10)
  write(100) 'BINARY'//char(10)
  write(100) 'DATASET UNSTRUCTURED_GRID'//char(10)
  write(ss,fmt='(A6,I10,A5)') 'POINTS', elems*4, 'float'
  write(100) ss//char(10)

  !points
  write(100) (x(i),y(i),z(i), i=1, elems*4)
  write(100) char(10)

  !cell topology
  write(ss,fmt='(A5,2I10)') 'CELLS', elems, elems*5
  write(100) char(10)//ss//char(10)
  write(100) cell
  write(100) char(10)

  !cell type
  write(ss,fmt='(A10,2I10)') 'CELL_TYPES', elems
  write(100) char(10)//ss//char(10)
  write(100) cell_type
  write(100) char(10)

  !data
  write(ss,fmt='(A10,I10)') 'POINT_DATA', elems*4
  write(100) char(10)//ss//char(10)
  write(100) 'SCALARS data float 1'//char(10)
  write(100) 'LOOKUP_TABLE default'//char(10)
  write(100) u1
  close(100)
end subroutine write_vtk_bin_scal_topology
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_avs_file_scal(u1, rows, filename1)

  implicit none
  
  integer, intent(in)                   :: rows
  real, dimension(1:rows), intent(in)   :: u1
  character (len=40),intent(in)         :: filename1

  integer                               :: i, t, ielem, ipol, jpol, ioerr
  real                                  :: phi
  character (len=30)                    :: celltype
  real, dimension(1:rows)               :: u2
 
  u2 = u1
 
  open(unit=200+mynum,file=trim(filename1)//'.inp',status='replace',form='formatted')
 
  celltype = 'pt'
  !avs file headers
  write(200+mynum,*) rows, rows, 1, 0, 0

  ! write coordinates and numbering
  do i=1, rows    
     write(200+mynum,'(i10,3f14.2)') i, W_vtk(i,1), W_vtk(i,2), W_vtk(i,3)
  enddo   
   
  ! connectivity structure
  do i=1,rows
     write(200+mynum,'(i12,i9,a4,i12)') i,1,trim(celltype),i
  enddo                   

  !property headers
  write(200+mynum,*) 1, 1
  write(200+mynum,*) 'u1,', 'm';

  !displacement fields
  do i=1, rows
     if (abs(u2(i)) < 1.e-25) u2(i) = 0.0
     write(200+mynum,'(i12,e14.3)') i, u2(i)
  enddo               
 
  close(200+mynum)

end subroutine write_avs_file_scal
!-----------------------------------------------------------------------------------------


!============================
end module input_output
!============================
