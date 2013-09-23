!========================
module kernel_meshes
!========================
 
    use data_mesh
    use global_parameters
    use input_output,  only: write_VTK_bin_scal, write_avs_file_scal
    use coord_mapping
  
    implicit none
  
    public :: init_read_decomp_save_rot_kernel_mesh, init_src_kernel_index_arrays
    private

contains 


!-----------------------------------------------------------------------------------------
subroutine init_read_decomp_save_rot_kernel_mesh(thr, phr)

  use data_arrays, only : kernmesh_type

  real(kind=realkind)           :: thr,phr
  character(len=10)             :: fraction
  integer                       :: ipol,jpol,iel,iproc,i

  cell_topology = .false.
  cell_mesh = 'coarse'

  if (have_fluid) then 
     call read_sem_meshes(thr)
  else
     if (mynum == 0) write(6,*) '------------ NO FLUID -------------'
     call read_sem_meshes_solid(thr)
  end if

  if (kernmesh_type == 'sem_mesh') then 
     if (mynum == 0) write(6,*)'computing kernels on SEM mesh.'
     call get_semkernel_mesh(globmesh)
  
  elseif (kernmesh_type == 'sem2mesh') then 
     if (mynum == 0) write(6,*)'computing kernels on SEM mesh with phi=pi mirror.'
     call get_semkernel_mesh2(globmesh)
  
  elseif (kernmesh_type == 'cmb_mesh') then 
     if (mynum == 0) write(6,*)'computing kernels on CMB fraction of SEM mesh .'
     fraction = 'cmb'
     call get_semkernel_fraction(globmesh, fraction, thr)
  
  elseif (kernmesh_type == 'mtz_mesh') then 
     if (mynum == 0) &
        write(6,*) 'computing kernels on transition zone fraction of SEM mesh .'
     fraction = 'mtz'
     call get_semkernel_fraction(globmesh, fraction, thr)
  
  elseif (kernmesh_type == 'mantlmsh') then 
     if (mynum == 0) write(6,*)'computing kernels on mantle fraction of SEM mesh .'
     fraction = 'mantle'
     call get_semkernel_fraction(globmesh, fraction, thr)
  
  elseif (kernmesh_type == '3dslices') then 
     if (ibeg == 0 .and. iend == npol) then 
        if (mynum == 0) write(6,*)'computing kernels on 3D cell geometry mesh on the fly.'
        if (cell_mesh == 'coarse') then
           call compute_3dslices_cell4(globmesh, phr) 
                                ! using ibeg and iend for the cell (4/el)
        else
           call compute_3dslices_cell16(globmesh, phr) 
                                ! using ibeg and iend for the cell (4/el)
        endif
     else
        write(6,*) 'Please run SOLVER with ibeg=iend=0 in inparam, if you'
        write(6,*) 'want to use 3dslices.'
        stop 2 
        ! Is not working reliably and not really important.
        !if (mynum==0) write(6,*)'computing kernels on 3D grid on the fly.'
        !if (cell_mesh=='coarse') then
        !   call compute_3dslices_4(globmesh,phr) ! using ibeg and iend (4 times as large!)
        !else
        !   call compute_3dslices(globmesh)  ! using gll point index 2 only
        !endif
     endif

  elseif (kernmesh_type == 'ext_mesh') then 
     if (mynum == 0) write(6,*) 'computing kernels on external mesh.'
     call read_external_mesh
  
  elseif (kernmesh_type == 'cub_mesh') then 
     if (mynum == 0) write(6,*) 'computing kernels on cubed sphere mesh.'
     call get_cs_mesh_part
     call prepare_dump_cs_part

  elseif (kernmesh_type == 'analytic') then 
     if (mynum == 0) write(6,*) 'computing kernels on analytical cross-section mesh.'
     call define_mymesh(thr)

  elseif (kernmesh_type == 'dbg_mesh') then 
     if (mynum == 0) write(6,*) 'computing kernels on debugging mesh.'
     call define_mymesh_dbg(thr)
  
  else
     write(6,*) 'unknown kernel mesh type:', kernmesh_type
     stop 2
  endif ! kernmesh_type

  if (.not. cell_topology) call decompose_kernel_mesh
 
  if (kernmesh_type == 'sem_mesh') then 
     call find_gridpoints_semmesh
  else
     call find_gridpoints_1d(globmesh)
  endif
  call save_kernelmesh
  call find_rotated_gridpoints_1d(maxr, thr, phr, globmesh)

  if (compute_src_kernels) &
       call init_src_kernel_index_arrays(maxr, thr, phr, globmesh)

  if (save_snaps) then 
     allocate(mesh1d(nsize,0:nproc_mesh-1,2))
     mesh1d = 0.
     do iproc=0, nproc_mesh-1
         i = 0
         do iel=1, nelem
             do jpol=ibeg, iend
                 do ipol=ibeg, iend
                     i= i + 1
                     mesh1d(i,iproc,1:2) = globmesh(ipol,jpol,iel,iproc,1:2)
                 enddo
             enddo
         enddo
     enddo
     write(6,*)mynum,'mesh1d:', i, nsize, maxval(mesh1d), maxval(globmesh)
     write(6,*)mynum,'mesh1d:', minval(mesh1d), minval(globmesh)
  endif
   
  deallocate(globmesh)
    
  end subroutine init_read_decomp_save_rot_kernel_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_src_kernel_index_arrays(maxr,thr,phr,mesh)

  use data_arrays, only : srcmij_kern,srcloc_kern,srcstf_kern
  use data_arrays, only : srcmij_kern_spec,srcloc_kern_spec
  use data_rota, only : ipt_srckern,iproc_srckern,iel_srckern,s_srckern,z_srckern
  use data_fft, only : nomega,ntimes
  use coord_mapping, only : rotate_frame_rd,define_kern_sem_mesh_mapping_arrays

  integer :: i
  real(kind=realkind), intent(in) :: thr,phr,maxr
  real(kind=realkind),intent(in) ::mesh(ibeg:iend,ibeg:iend,1:nelem,0:nproc_mesh-1,1:2)
  real(kind=realkind) :: rsrc,thsrc
  integer :: iproc,iel,iel_min,iproc_min,ipt
  real(kind=realkind), allocatable, dimension(:) :: ssrc,phisrc,zsrc,phird_src,srd_src,zrd_src

  integer :: ipol,jpol,ii,nelpts
  real(kind=realkind) :: dr,r,th
  integer :: ipol_min,jpol_min,ipt_count,ipt_count_min
  integer, dimension(:), allocatable :: iipt_count
  real(kind=realkind), dimension(:), allocatable :: mindist
  integer :: thetadisttmp(2),iel3(3)

  write(6,*)'Preparing source kernel computation...'

!!!!!!!!!!!!!!!!!!!!!
  n_src=1
!!!!!!!!!!!!!!!!!!!!!

! allocate arrays
  allocate(srcmij_kern(1:ntimes,6),srcloc_kern(1:ntimes,3),srcstf_kern(1:ntimes))
  allocate(srcmij_kern_spec(0:nomega,6),srcloc_kern_spec(0:nomega,3))

  allocate(ssrc(n_src),phisrc(n_src),zsrc(n_src),phird_src(n_src),srd_src(n_src),zrd_src(n_src))
  allocate(azim1_src(n_src),azim2_src(n_src))

! find index that corresponds to source location
   open(unit=20000,file='sourceparams_fwd.dat')
   read(20000,*)realjunk
   read(20000,*)src_type1fwd
   read(20000,*)src_type2fwd
   read(20000,*)
   read(20000,*)
   read(20000,*)src_depth
   close(20000)
  write(6,*)'maxr, src_depth:',maxr,src_depth
  zsrc=maxr-src_depth*1000.
  write(6,*)'source depth:',zsrc
  ssrc=0.
  phisrc=0.

! rotate source coordinates according to receiver location
  call rotate_frame_rd(n_src,srd_src,phird_src,zrd_src,ssrc,phisrc,zsrc,phr,thr)
  do i=1,n_src
     rsrc=sqrt(srd_src(i)**2 + zrd_src(i)**2)
     thsrc=atan2(srd_src(i),zrd_src(i))
!     write(6,*)'source coordinates:',i,'s,z,phi:',srd_src(i)/1000.,zrd_src(i)/1000.,phird_src(i)*180./pi
     write(6,*)'rotated source coordinates:',i,'r,th:',rsrc/1000.,thsrc*180./pi
  enddo

! find rotated source location index in unrotated receiver frame
  allocate(mindist(n_src))
  allocate(ipt_srckern(n_src),iproc_srckern(n_src),iel_srckern(n_src))
  nelpts=(iend-ibeg+1)**2
  allocate(s_srckern(nelpts,n_src),z_srckern(nelpts,n_src)  )
  do ipt=1,n_src
    mindist(ipt) = maxr
     do iproc=0, nproc_mesh-1
        do iel = 1, nelem
           dr = sqrt(( mesh(npol/2,npol/2,iel,iproc,1) - srd_src(ipt)  )**2+ &
                      ( mesh(npol/2,npol/2,iel,iproc,2) - zrd_src(ipt)  )**2)
           if (dr<mindist(ipt)) then
              iel_min=iel ; iproc_min = iproc; mindist(ipt) = dr
           end if 
          if (mynum==0 .and. ipt==1) write(4000+iproc,*) &
               mesh(npol/2,npol/2,iel,iproc,1),mesh(npol/2,npol/2,iel,iproc,2) 
        end do !iel
     end do !iproc
     
     ipt_count = (iel_min-1)*(iend-ibeg+1)**2
     ii=0
     do jpol=ibeg,iend
        do ipol=ibeg,iend
           ii=ii+1
           ipt_count = ipt_count + 1
              dr = sqrt(( mesh(ipol,jpol,iel_min,iproc_min,1) - srd_src(ipt) )**2 + &
                         ( mesh(ipol,jpol,iel_min,iproc_min,2) - zrd_src(ipt) )**2)
           if (dr<=mindist(ipt)) then
              mindist(ipt) = dr
              ipol_min = ipol; jpol_min = jpol ; ipt_count_min = ipt_count
           end if
           s_srckern(ii,ipt)=mesh(ipol,jpol,iel_min,iproc_min,1) 
           z_srckern(ii,ipt)=mesh(ipol,jpol,iel_min,iproc_min,2) 

        end do
     end do

     write(9000+mynum,*)iel_min,iproc_min,ipt_count_min

     write(1000+mynum,*)srd_src(ipt),zrd_src(ipt),mindist(ipt)
     ipt_srckern(ipt) = ipt_count_min
     iproc_srckern(ipt) = iproc_min
     iel_srckern(ipt)=iel_min

     rsrc=sqrt( mesh(ipol_min,jpol_min,iel_srckern(ipt),iproc_srckern(ipt),1)**2 + &
                       mesh(ipol_min,jpol_min,iel_srckern(ipt),iproc_srckern(ipt),2)**2)
     thsrc=atan2(mesh(ipol_min,jpol_min,iel_srckern(ipt),iproc_srckern(ipt),1), &
                           mesh(ipol_min,jpol_min,iel_srckern(ipt),iproc_srckern(ipt),2))
     write(6,*)ipt,'src loc mesh:',rsrc/1000.,thsrc*180./pi

     do ii=1,nelpts
        rsrc=sqrt(s_srckern(ii,ipt)**2 + z_srckern(ii,ipt)**2)
        thsrc=atan2(s_srckern(ii,ipt),z_srckern(ii,ipt))
        write(77777,*)rsrc/1000.,thsrc*180./pi
     enddo

  enddo !n_src

   deallocate(mindist,srd_src,zrd_src)

! azimuthal prefactors

     ! backward fields
     if (src_type2bwd=='mxz' .or. src_type2bwd=='xforce') then
        if (mynum==0) write(6,*) &
             '  Calculating prefactors for correct azimuth for',src_type2bwd
        azim1_src =  cos(phird_src)
        if (allocated(azim2_bwd)) azim2_src =  -sin(phird_src)

     elseif (src_type2bwd=='myz' .or. src_type2bwd=='yforce') then 
        if (mynum==0) write(6,*) &
             '  Calculating prefactors for correct azimuth for',src_type2bwd
        azim1_src =  sin(phird_src)
        if (allocated(azim2_src)) azim2_src=  cos(phird_src)
        
     elseif (src_type2bwd=='mxx_m_myy') then 
        if (mynum==0) write(6,*) &
             '  Calculating prefactors for correct azimuth for',src_type2bwd
        azim1_src =  cos(2.d0*phird_src)
        if (allocated(azim2_bwd )) azim2_src  =  -sin(2.d0*phird_src)
        
     elseif (src_type2bwd=='mxy') then 
        if (mynum==0) write(6,*) &
             '  Calculating prefactors for correct azimuth for',src_type2bwd
        azim1_src =  sin(2.d0*phird_src)
        if (allocated(azim2_bwd)) azim2_src =  cos(2.d0*phird_src)
     endif

write(6,*)'.... done preparing source kernel indices etc.'

end subroutine init_src_kernel_index_arrays
!----------------------------------------------------------------------------------------------------------

!dk read_sem_mesh-------------------------------------
 subroutine read_sem_meshes(thr)

   include 'mesh_params.h'
   include 'mesh_params_kernel.h'

   real(kind=realkind), intent(in)          :: thr
   real(8), allocatable, dimension(:,:,:,:) :: mesh_sol, mesh_flu, mesh
   integer                                  :: it, ndumpstmp
   integer                                  :: ipol, jpol, iel, iproc
   character(len=4)                         :: appiproc
   real(kind=realkind)                      :: r, th

   real(kind=realkind) :: discont1,discont2,r1
   integer :: iel2,iidom

   nsize     = (iend-ibeg+1)**2 * nelem
   nsize_sol = (iend-ibeg+1)**2 * nel_solid
   nsize_flu = (iend-ibeg+1)**2 * nel_fluid

   write(6,*)'number of global points in SEM mesh:',nsize
   write(6,*)'number of solid,fluid points in SEM mesh:',nsize_sol,nsize_flu

   if (nel_solid+nel_fluid /= nelem) then 
      write(6,*)'solid & fluid elements dont sum up to total!'
      write(6,*)'nel_solid,nel_fluid:',nel_solid,nel_fluid
      write(6,*)'nelem',nelem
      stop
   endif

   if (mynum==0) write(6,*)'loading sem mesh from',trim(dir_fwdmesh1(1))//'/Data'//'/strain_mesh_sol_*'

   do iproc=0,nproc_mesh-1
      call define_io_appendix(appiproc,iproc)

 !     allocate(mesh_sol(ibeg:iend,ibeg:iend,1:nel_solid,2))
      allocate(mesh_sol(0:iend,0:iend,1:nel_solid,2))
      if (have_fluid) allocate(mesh_flu(0:iend,0:iend,1:nel_fluid,2))
      
      open(unit=88,file=trim(dir_fwdmesh1(1))//'/Data'// &
                             & '/strain_mesh_sol_'//appiproc//'.dat', &
           FORM="UNFORMATTED",STATUS="OLD")
      read(88) mesh_sol(ibeg:iend,ibeg:iend,1:nel_solid,1), &
               mesh_sol(ibeg:iend,ibeg:iend,1:nel_solid,2)

      close(88)

      if (have_fluid) then
         open(unit=78,file=trim(dir_fwdmesh1(1))//'/Data'// &
                                & '/strain_mesh_flu_'//appiproc//'.dat', &
              FORM="UNFORMATTED",STATUS="OLD")
         read(78) mesh_flu(ibeg:iend,ibeg:iend,1:nel_fluid,1), &
                  mesh_flu(ibeg:iend,ibeg:iend,1:nel_fluid,2)

         close(78)
      endif

! Saving mesh
      allocate(mesh(ibeg:iend,ibeg:iend,1:nelem,2))
      if (have_fluid) then
          mesh(ibeg:iend,ibeg:iend,1:nel_fluid,1:2) = &
               mesh_flu(ibeg:iend,ibeg:iend,1:nel_fluid,1:2)
      end if
      mesh(ibeg:iend,ibeg:iend,nel_fluid+1:nelem,1:2) = &
          mesh_sol(ibeg:iend,ibeg:iend,1:nel_solid,1:2)
      if (mynum==0) write(6,*)'saving mesh..........'
      open(unit=50,file='Data/'//'mesh'//appiproc//'.dat')
      do iel=1,nelem
         do jpol=ibeg,iend
            do ipol=ibeg,iend
               write(50,*) mesh(ipol,jpol,iel,1), mesh(ipol,jpol,iel,2)
            enddo
         enddo
      enddo
      close(50)

! Run some mesh tests
      if (mynum==0) write(6,*)'some mesh tests...'
      ! call sem_mesh_rotation_tests(mesh)
      deallocate(mesh,mesh_sol)
      if (have_fluid) deallocate(mesh_flu)
   enddo ! iproc

   allocate(mesh_sol(ibeg:iend,ibeg:iend,1:nel_solid,2))
   if (have_fluid) allocate(mesh_flu(ibeg:iend,ibeg:iend,1:nel_fluid,2))

! load wavefied information
   if (mynum==0) write(6,*)'loading wavefield info..........'
   open(unit=97,file=trim(dir_fwdmesh1(1))//'/Data'//'/strain_info.dat0000')
   read(97,*)ndumpstmp
   if (ndumpstmp/=ndumps) then
      write(6,*)'  Problem with number of dumps!'
      write(6,*)'  mesh_params_kernel.h:',ndumps
      write(6,*)'  strain_info.dat0000 :',ndumpstmp
      stop
   endif
 
   if (mynum==0) then
      write(6,*)
      write(6,*)'number of wavefield snapshots:',ndumps
   endif
   call flush(6)
   
   do it=1,ndumps
      read(97,*)time(it),timestep(it)
   enddo
   close(97)
 
   if (mynum==0) then
      write(6,*)
      write(6,*)'loading global mesh to each processor........'
   endif
   call flush(6)
   
   allocate(globmesh(ibeg:iend,ibeg:iend,1:nelem,0:nproc_mesh-1,2))
   
   allocate(iidom_glob(nelem,0:nproc_mesh-1))
 
   do iproc=0,nproc_mesh-1
 
      call define_io_appendix(appiproc,iproc)
 
      open(unit=88,file=trim(dir_fwdmesh1(1))//'/Data'// &
                        '/strain_mesh_sol_'//appiproc//'.dat', &
                        FORM="UNFORMATTED",STATUS="OLD")
      read(88) mesh_sol(ibeg:iend,ibeg:iend,:,1), mesh_sol(ibeg:iend,ibeg:iend,:,2)
      close(88)
 
      if (have_fluid) then
          open(unit=78,file=trim(dir_fwdmesh1(1))//'/Data'// &
                            '/strain_mesh_flu_'//appiproc//'.dat', &
               FORM="UNFORMATTED",STATUS="OLD")
          read(78) mesh_flu(ibeg:iend,ibeg:iend,:,1), mesh_flu(ibeg:iend,ibeg:iend,:,2)
          close(78)
      endif
 
      if (have_fluid) then
          globmesh(ibeg:iend,ibeg:iend, 1:nel_fluid      , iproc, : ) = mesh_flu
          globmesh(ibeg:iend,ibeg:iend, nel_fluid+1:nelem, iproc, : ) = mesh_sol
      else
          globmesh(ibeg:iend,ibeg:iend, 1:nel_solid,       iproc, : ) = mesh_sol
      endif
 
 !     need domain for medium parameters
      open(unit=65,file=trim(dir_fwdmesh1(1))//'/Info' &
           //'/elems_bkgrdmodel_domain.dat'//appiproc,status='old')
      do iel=1,nelem
         read(65,10) iel2, r1, iidom, discont1, discont2
         iidom_glob(iel2,iproc) = iidom
      enddo
      close(65)
 
    enddo ! iproc = 0, nproc_mesh

10  format(i9,1pe11.3,i3,2(1pe11.3))

    deallocate(mesh_sol)
    if (have_fluid) deallocate (mesh_flu)
 
    if (mynum==0) then
       write(6,*)
       write(6,*)'computing rotated receiver matrix & mesh........??????' 
    endif
    call flush(6)
 
!   read background model domains and discontinuities
    open(unit=65,file=trim(dir_fwdmesh1(1))//'/Info'// &
                      '/discontinuities_solver.dat0000',status='old')
    do iidom=1,ndisc
       read(65,*)discont(iidom)
       if (mynum==0)write(6,*)iidom,'discontinuity:',discont(iidom)
    enddo
    close(65)
  
    maxr = maxval(discont)
    ar = 1/maxr
 
    write(6,*)'maximum radius maxr:',maxr

 end subroutine read_sem_meshes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_external_mesh
  use coord_trafo, only : sphi2xy,rthetaphi2xyz

  real(kind=realkind)   :: x1, x2, x3
  integer               :: ipt

  if (mynum == 0) write(6,*) 'loading external mesh...'

  if (mynum == 0) write(6,*) ' External mesh:', ext_mesh_name(1:lfext)
  open(unit=91, file=ext_mesh_name(1:lfext), status='old')
  read(91,*) npts_tot

  if (mynum == 0) then
     write(6,*)
     write(6,*) 'External mesh total grid points:', npts_tot
     write(6,*) 'Kernel procs:', nproc
     write(6,*) 'SEM procs:', nproc_mesh
  endif

  allocate(xgd_tot(1:npts_tot), ygd_tot(1:npts_tot), zgd_tot(npts_tot))

  if (ext_mesh_name(lfext-2:lfext) == 'sph') then 
     if(mynum == 0) write(6,*)'reading external mesh in spherical coordinates...'
     do ipt=1, npts_tot
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

        call rthetaphi2xyz(xgd_tot(ipt), ygd_tot(ipt), zgd_tot(ipt), x1, x2, x3)
     enddo

  elseif (ext_mesh_name(lfext-2:lfext) == 'cyl') then 
     if (mynum == 0) write(6,*) 'reading external mesh in cylindrical coordinates..'
     do ipt=1, npts_tot
        read(91,*) x1, zgd(ipt), x3
        if (x1 > maxr .or. x1 < 0.) then 
           write(6,*) mynum, 'problem with coordinate s in external mesh:', ipt, x1
           stop
        endif
        if (abs(zgd(ipt)) > maxr) then 
           write(6,*) mynum, 'problem with coordinate z in external mesh:', ipt, x2
           stop
        endif
        if (x3 > 2.*pi .or. x3 < 0.) then 
           write(6,*) mynum, 'problem with phi in external mesh:', ipt, x3
           stop
        endif

        call sphi2xy(xgd_tot(ipt), ygd_tot(ipt), x1, x3)
     enddo

  elseif (ext_mesh_name(lfext-2:lfext) == 'xyz') then 
    if(mynum == 0) write(6,*) 'reading external mesh in cartesian coordinates...'
    do ipt=1, npts_tot
       read(91,*) xgd_tot(ipt), ygd_tot(ipt), zgd_tot(ipt)
       if (sqrt(xgd_tot(ipt)**2 + ygd_tot(ipt)**2 + zgd_tot(ipt)**2) > 1.0001 * maxr) then
          write(6,*) mynum, 'problem with coordinates in external mesh, radius:', ipt, &
                    sqrt(xgd_tot(ipt)**2 + ygd_tot(ipt)**2 + zgd_tot(ipt)**2) 
          stop
       elseif (abs(xgd_tot(ipt)) > maxr) then
          write(6,*) mynum, 'problem with x in external mesh:', ipt, xgd_tot(ipt)
          stop
       elseif (abs(ygd_tot(ipt)) > maxr) then
          write(6,*) mynum, 'problem with y in external mesh:', ipt, ygd_tot(ipt)
          stop
       elseif (abs(zgd_tot(ipt)) > maxr) then
          write(6,*) mynum, 'problem with z in external mesh:', ipt, zgd_tot(ipt)
          stop
       endif
     enddo
  else 
     if(mynum == 0) then 
        write(6,*) 'external mesh type not recognized!'
        write(6,*) 'mesh file name must end with:'
        write(6,*) " 1) 'xyz' -- cartesian coordinates (x[m], y[m], z[m])"
        write(6,*) " 2) 'sph' -- spherical coordinates (r[m], theta[rad], phi[rad])"
        write(6,*) " 3) 'cyl' -- cylindrical coordinates (s[m], z[m], phi[rad])"
     endif
     stop
  endif

  close(91)

end subroutine read_external_mesh
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_semkernel_mesh(mesh)
  
  use data_rota
  use mpi

  real(kind=realkind), intent(in)    :: mesh(ibeg:iend,ibeg:iend,nelem,0:nproc_mesh-1,1:2)
  real(kind=realkind), dimension(:), allocatable :: sgd_tot
  integer                            :: ipt, ipol, jpol, iel
  integer                            :: iproc, ipt_in_proc

  if (mynum==0) write(6,*)'define kernel mesh upon initial SEM mesh....'

  npts_tot = nproc_mesh * nelem * (iend-ibeg+1)**2
  allocate(sgd_tot(npts_tot))
  allocate(xgd_tot(npts_tot))
  allocate(ygd_tot(npts_tot))
  allocate(zgd_tot(npts_tot))
  
  allocate(iproc_proc_tot(npts_tot))
  allocate(iel_proc_tot(npts_tot))
  allocate(ipt_proc_tot(npts_tot))

  if (mynum==0) then
     write(6,*)
     write(6,*)'Total grid points:', npts_tot
     write(6,*)'Kernel procs:',      nproc
     write(6,*)'SEM procs:',         nproc_mesh
  endif

  sgd_tot = 0. 
  zgd_tot = 0.
  iproc_proc_tot = 0
  ipt=0 
  do iproc=0, nproc_mesh-1
     ipt_in_proc = 0
     do iel = 1, nelem
        do jpol=ibeg,iend
           do ipol=ibeg,iend
              ipt=ipt+1
              ipt_in_proc = ipt_in_proc + 1
              sgd_tot(ipt) = mesh(ipol,jpol,iel,iproc,1)
              zgd_tot(ipt) = mesh(ipol,jpol,iel,iproc,2)

            ! define fwd mapping at the global level
              iproc_proc_tot(ipt) = iproc
              iel_proc_tot(ipt) = iel
              ipt_proc_tot(ipt) = ipt_in_proc
           enddo
        enddo
     enddo
  enddo

  xgd_tot = sgd_tot
  ygd_tot = 0. ! this is only to avoid rotated domains when searching s,z
               ! later, we will use phi0 to compute things at the right y.

  deallocate(sgd_tot)

end subroutine get_semkernel_mesh
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine get_semkernel_fraction(mesh,fraction,thr)
 
  use data_rota
  use coord_trafo, only : xyz2rthetaphi
  use input_output, only : write_vtk_bin_scal_topology,write_vtk_bin_scal_mesh3d
  use data_arrays, only : kern_cell
  use mpi

  real, intent(in) ::  thr
  character(len=10), intent(in) :: fraction
  real(kind=realkind),intent(in) ::mesh(ibeg:iend,ibeg:iend,1:nelem,0:nproc_mesh-1,1:2)
  real(kind=realkind), dimension(:), allocatable :: sgd_tot,xgd_tottmp,zgd_tottmp
  integer :: ipt,ipol,jpol,iel, iproc, ipt_in_proc,ii,i,j,imid
  real :: theta,r,phi,y,s,z,r1,r2,th1,th2,rtmp(16)
  character(len=200) :: filename1
  character(len=4) :: appiproc
  real(kind=realkind), dimension(:), allocatable :: vptot

  if (trim(fraction)=='cmb') then 
     if (mynum==0)write(6,*)'cmb fractional mesh'
     r1=3482000. ;  r2=3631000.
     th1=thr/2.*0.8; th2=thr/2.*1.2
     lfext=17
     ext_mesh_name(1:lfext)='Data/CMB_mesh.xyz'
  elseif (trim(fraction)=='mtz') then 
     if (mynum==0)write(6,*)'transition zone fractional mesh'
     r1=5711000. ;  r2=5961000. 
     th1=thr/2.*0.97;  th2= thr/2.*1.03
     lfext=17
     ext_mesh_name(1:lfext)='Data/MTZ_mesh.xyz'
  elseif (trim(fraction)=='mantle') then 
     if (mynum==0)write(6,*)'whole-mantle fractional mesh'
     if (thr<=60./180.*pi) then  
        r1=4000000. ;  r2=6372000.
     elseif (thr<=97./180.*pi) then  
        r1=3480000. ;  r2=6372000.
     elseif (thr<=130./180.*pi) then  
        r1=1200000. ;  r2=6372000.
     else
        r1=0. ;  r2=6372000.
     endif
     write(6,*)'kernel mesh radii:',r1,r2
     th1=0.;  th2= thr*1.1
     lfext=20
     ext_mesh_name(1:lfext)='Data/MANTLE_mesh.xyz'
  endif

!====================================================================
if (ibeg==0 .and. iend==npol) then 
!====================================================================
  cell_topology=.true.

  if (trim(cell_mesh)=='coarse') then 
     if (mynum==0) write(6,*)'define fractional kernel mesh upon initial SEM mesh: 4-point cell topology....'
     npts_cell=nproc_mesh*nelem*4
     allocate(xgd_tottmp(npts_cell))
     allocate(zgd_tottmp(npts_cell))
     allocate(vptot(npts_cell))
     xgd_tottmp = 0.
     y=0.
     zgd_tottmp = 0.
     ii=0
     
     if (mynum==0)write(6,*)'big loop to find fractional mesh... max cell points:',npts_cell
     do iproc=0, nproc_mesh-1
        if (mynum==0)write(6,*)'processor',iproc
        ipt_in_proc = 0
        do iel = 1, nelem
           s = (mesh(ibeg,ibeg,iel,iproc,1) + mesh(iend,iend,iel,iproc,1))/2.
           z = (mesh(ibeg,ibeg,iel,iproc,2) + mesh(iend,iend,iel,iproc,2))/2.
           call xyz2rthetaphi(r,theta,phi,s,y,z)
           if (r>=r1 .and. r <=r2  .and. theta>=th1 .and. theta <= th2) then
              xgd_tottmp(ii+1) = mesh(ibeg,ibeg,iel,iproc,1)
              zgd_tottmp(ii+1) = mesh(ibeg,ibeg,iel,iproc,2)
              xgd_tottmp(ii+2) = mesh(iend,ibeg,iel,iproc,1)
              zgd_tottmp(ii+2) = mesh(iend,ibeg,iel,iproc,2)
              xgd_tottmp(ii+3) = mesh(iend,iend,iel,iproc,1)
              zgd_tottmp(ii+3) = mesh(iend,iend,iel,iproc,2)
              xgd_tottmp(ii+4) = mesh(ibeg,iend,iel,iproc,1)
              zgd_tottmp(ii+4) = mesh(ibeg,iend,iel,iproc,2)
              rtmp(1:4) = sqrt(xgd_tottmp(ii+1:ii+4)**2 + zgd_tottmp(ii+1:ii+4)**2 )
              rtmp(1:4) = rtmp(1:4) + 1.e-4*(r-rtmp(1:4))
              do i=1,4
                 vptot(ii+i)=prem(rtmp(i),'v_p')
              enddo
              ii=ii+4
              if (theta<=th2-thr .and. th1<=abs(th2-thr)) then !add some parts from "behind" the source
                 xgd_tottmp(ii+1) =-mesh(ibeg,ibeg,iel,iproc,1)
                 zgd_tottmp(ii+1) = mesh(ibeg,ibeg,iel,iproc,2)
                 xgd_tottmp(ii+2) =-mesh(iend,ibeg,iel,iproc,1)
                 zgd_tottmp(ii+2) = mesh(iend,ibeg,iel,iproc,2)
                 xgd_tottmp(ii+3) =-mesh(iend,iend,iel,iproc,1)
                 zgd_tottmp(ii+3) = mesh(iend,iend,iel,iproc,2)
                 xgd_tottmp(ii+4) =-mesh(ibeg,iend,iel,iproc,1)
                 zgd_tottmp(ii+4) = mesh(ibeg,iend,iel,iproc,2) 
                 do i=1,4
                    vptot(ii+i) = prem(rtmp(i),'v_p')
                 enddo
                 ii=ii+4
              endif
           endif

        enddo
     enddo
  else
     if (mynum==0) write(6,*)'define fractional kernel mesh upon initial SEM mesh: 16-point cell topology....'
     npts_cell=nproc_mesh*nelem*16
     allocate(xgd_tottmp(npts_cell),zgd_tottmp(1:npts_cell),vptot(1:npts_cell))
     xgd_tottmp = 0.; y=0.; zgd_tottmp = 0.; ii=0
     if (mynum==0)write(6,*)'big loop to find fractional mesh... max cell points:',npts_cell
     imid=iend/2
     write(6,*)'mid point index (16-point topology):',imid
     do iproc=0, nproc_mesh-1
        if (mynum==0)write(6,*)'processor',iproc
        ipt_in_proc = 0
        do iel = 1, nelem
           s = (mesh(ibeg,ibeg,iel,iproc,1)+mesh(iend,iend,iel,iproc,1))/2.
           z = (mesh(ibeg,ibeg,iel,iproc,2)+mesh(iend,iend,iel,iproc,2))/2.
           call xyz2rthetaphi(r,theta,phi,s,y,z)
           if (r>=r1 .and. r <=r2  .and. theta>=th1 .and. theta <= th2) then
              xgd_tottmp(ii+1)=mesh(ibeg,ibeg,iel,iproc,1);zgd_tottmp(ii+1)=mesh(ibeg,ibeg,iel,iproc,2)
              xgd_tottmp(ii+2)=mesh(imid,ibeg,iel,iproc,1);zgd_tottmp(ii+2)=mesh(imid,ibeg,iel,iproc,2)
              xgd_tottmp(ii+3)=mesh(imid,imid,iel,iproc,1);zgd_tottmp(ii+3)=mesh(imid,imid,iel,iproc,2)
              xgd_tottmp(ii+4)=mesh(ibeg,imid,iel,iproc,1);zgd_tottmp(ii+4)=mesh(ibeg,imid,iel,iproc,2)

              xgd_tottmp(ii+5)=mesh(imid,ibeg,iel,iproc,1);zgd_tottmp(ii+5)=mesh(imid,ibeg,iel,iproc,2)
              xgd_tottmp(ii+6)=mesh(iend,ibeg,iel,iproc,1);zgd_tottmp(ii+6)=mesh(iend,ibeg,iel,iproc,2)
              xgd_tottmp(ii+7)=mesh(iend,imid,iel,iproc,1);zgd_tottmp(ii+7)=mesh(iend,imid,iel,iproc,2)
              xgd_tottmp(ii+8)=mesh(imid,imid,iel,iproc,1);zgd_tottmp(ii+8)=mesh(imid,imid,iel,iproc,2)

              xgd_tottmp(ii+9)=mesh(imid,imid,iel,iproc,1);zgd_tottmp(ii+9)=mesh(imid,imid,iel,iproc,2)
              xgd_tottmp(ii+10)=mesh(iend,imid,iel,iproc,1);zgd_tottmp(ii+10)=mesh(iend,imid,iel,iproc,2)
              xgd_tottmp(ii+11)=mesh(iend,iend,iel,iproc,1);zgd_tottmp(ii+11)=mesh(iend,iend,iel,iproc,2)
              xgd_tottmp(ii+12)=mesh(imid,iend,iel,iproc,1);zgd_tottmp(ii+12)=mesh(imid,iend,iel,iproc,2)

              xgd_tottmp(ii+13)=mesh(ibeg,imid,iel,iproc,1);zgd_tottmp(ii+13)=mesh(ibeg,imid,iel,iproc,2)
              xgd_tottmp(ii+14)=mesh(imid,imid,iel,iproc,1);zgd_tottmp(ii+14)=mesh(imid,imid,iel,iproc,2)
              xgd_tottmp(ii+15)=mesh(imid,iend,iel,iproc,1);zgd_tottmp(ii+15)=mesh(imid,iend,iel,iproc,2)
              xgd_tottmp(ii+16)=mesh(ibeg,iend,iel,iproc,1);zgd_tottmp(ii+16)=mesh(ibeg,iend,iel,iproc,2)

              rtmp(1:16)=sqrt(xgd_tottmp(ii+1:ii+16)**2 + zgd_tottmp(ii+1:ii+16)**2 )
              rtmp(1:16)=rtmp(1:16) + 1.e-4*(r-rtmp(1:16))
              do i=1,16
                 vptot(ii+i)=prem(rtmp(i),'v_p')
              enddo
              ii=ii+16
              if (theta<=th2-thr .and. th1<=abs(th2-thr)) then !add some parts from "behind" the source
                 xgd_tottmp(ii+1)=-mesh(ibeg,ibeg,iel,iproc,1);zgd_tottmp(ii+1)=mesh(ibeg,ibeg,iel,iproc,2)
                 xgd_tottmp(ii+2)=-mesh(iend,ibeg,iel,iproc,1);zgd_tottmp(ii+2)=mesh(iend,ibeg,iel,iproc,2)
                 xgd_tottmp(ii+3)=-mesh(iend,iend,iel,iproc,1);zgd_tottmp(ii+3)=mesh(iend,iend,iel,iproc,2)
                 xgd_tottmp(ii+4)=-mesh(ibeg,iend,iel,iproc,1);zgd_tottmp(ii+4)=mesh(ibeg,iend,iel,iproc,2) 

                 xgd_tottmp(ii+5)=-mesh(imid,ibeg,iel,iproc,1);zgd_tottmp(ii+5)=mesh(imid,ibeg,iel,iproc,2)
                 xgd_tottmp(ii+6)=-mesh(iend,ibeg,iel,iproc,1);zgd_tottmp(ii+6)=mesh(iend,ibeg,iel,iproc,2)
                 xgd_tottmp(ii+7)=-mesh(iend,imid,iel,iproc,1);zgd_tottmp(ii+7)=mesh(iend,imid,iel,iproc,2)
                 xgd_tottmp(ii+8)=-mesh(imid,imid,iel,iproc,1);zgd_tottmp(ii+8)=mesh(imid,imid,iel,iproc,2)

                 xgd_tottmp(ii+9)=-mesh(imid,imid,iel,iproc,1);zgd_tottmp(ii+9)=mesh(imid,imid,iel,iproc,2)
                 xgd_tottmp(ii+10)=-mesh(iend,imid,iel,iproc,1);zgd_tottmp(ii+10)=mesh(iend,imid,iel,iproc,2)
                 xgd_tottmp(ii+11)=-mesh(iend,iend,iel,iproc,1);zgd_tottmp(ii+11)=mesh(iend,iend,iel,iproc,2)
                 xgd_tottmp(ii+12)=-mesh(imid,iend,iel,iproc,1);zgd_tottmp(ii+12)=mesh(imid,iend,iel,iproc,2)
                 
                 xgd_tottmp(ii+13)=-mesh(ibeg,imid,iel,iproc,1);zgd_tottmp(ii+13)=mesh(ibeg,imid,iel,iproc,2)
                 xgd_tottmp(ii+14)=-mesh(imid,imid,iel,iproc,1);zgd_tottmp(ii+14)=mesh(imid,imid,iel,iproc,2)
                 xgd_tottmp(ii+15)=-mesh(imid,iend,iel,iproc,1);zgd_tottmp(ii+15)=mesh(imid,iend,iel,iproc,2)
                 xgd_tottmp(ii+16)=-mesh(ibeg,iend,iel,iproc,1);zgd_tottmp(ii+16)=mesh(ibeg,iend,iel,iproc,2)

                 do i=1,16
                    vptot(ii+i)=prem(rtmp(i),'v_p')
                 enddo
                 ii=ii+16
              endif
           endif

        enddo
     enddo
  endif

  write(6,*)'number of cell points:',npts_cell
  npts_cell = ii
  filename1='Data/mesh_all_cell'
  if (mynum==nproc-1) &
  call write_VTK_bin_scal_topology(real(xgd_tottmp(1:npts_cell)),real(0.*xgd_tottmp(1:npts_cell)),&
                  real(zgd_tottmp(1:npts_cell)),real(vptot(1:npts_cell)),npts_cell/4,filename1)

! reducing overlapping coordinates: double loop
     allocate(ind_tot2cell(1:npts_cell,10))
     allocate(ind_cell2tot(1:npts_cell))
     allocate(valence_tot2cell(1:npts_cell))
     allocate(xgd_cell(1:npts_cell))
     allocate(ygd_cell(1:npts_cell))
     allocate(zgd_cell(1:npts_cell))
     xgd_cell=xgd_tottmp(1:npts_cell)
     ygd_cell=0.
     zgd_cell=zgd_tottmp(1:npts_cell)
 
     ! THIS INCLUDES DOMAIN DECOMPOSITION!!
     call reduce_coord_overlap(npts_cell,xgd_cell,ygd_cell,zgd_cell,npts,&
                               ind_tot2cell,ind_cell2tot,valence_tot2cell,nloc)

     write(6,*)'allocating xgd, size:',npts
     allocate(xgd(1:npts),ygd(1:npts),zgd(1:npts))
     do i=1,nloc
        xgd(ind_cell2tot(i)) = xgd_cell(i)
        ygd(ind_cell2tot(i)) = ygd_cell(i)
        zgd(ind_cell2tot(i)) = zgd_cell(i)
     enddo
     deallocate(xgd_cell,ygd_cell,zgd_cell)
     call define_io_appendix(appiproc,mynum)
     filename1='Data/mesh_all_totpts_'//appiproc
     call write_VTK_bin_scal_mesh3d(xgd,ygd,zgd,zgd,npts,0,filename1)

! test the mapping cell to total and back ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     allocate(xgd_cell(1:npts_cell),ygd_cell(1:npts_cell),zgd_cell(1:npts_cell))
     xgd_cell=0; ygd_cell=0;zgd_cell=0.
     write(6,*)'max ind:',maxval(ind_tot2cell),maxval(ind_cell2tot)
     do i=1,npts
        do j=1,valence_tot2cell(i)
           xgd_cell(ind_tot2cell(i,j)) = xgd(i)
           ygd_cell(ind_tot2cell(i,j)) = ygd(i)
           zgd_cell(ind_tot2cell(i,j)) = zgd(i)
        enddo
     enddo
     filename1='Data/mesh_all_cell_retry_'//appiproc
     call write_VTK_bin_scal_topology(real(xgd_cell(1:nloc)),real(ygd_cell(1:nloc)),real(zgd_cell(1:nloc)),&
          vptot(1:nloc),nloc/4,filename1)
     ! end test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     write(6,*)mynum,'Done with 3D slices in 4-corner cell topology.'
     allocate(kern_cell(nloc))

!====================================================================
else  
!====================================================================
  if (mynum==0) write(6,*)'define fractional kernel mesh upon initial SEM mesh....'

  npts_tot=nproc_mesh*nelem*(iend-ibeg+1)**2
  allocate(xgd_tottmp(npts_tot),zgd_tottmp(1:npts_tot))

  y=0.
  xgd_tottmp = 0.; zgd_tottmp = 0.; ii=0
  if (mynum==0)write(6,*)'big loop to find fractional mesh'
  do iproc=0, nproc_mesh-1
     if (mynum==0)write(6,*)'processor',iproc
     ipt_in_proc = 0
     do iel = 1, nelem
        do jpol=ibeg,iend
           do ipol=ibeg,iend
              ipt_in_proc = ipt_in_proc + 1
              s = mesh(ipol,jpol,iel,iproc,1);  z = mesh(ipol,jpol,iel,iproc,2)
              call xyz2rthetaphi(r,theta,phi,s,y,z)
              if (r>=r1 .and. r <=r2  .and. theta>=th1 .and. theta <= th2) then
                 ii=ii+1
                 xgd_tottmp(ii)=s;  zgd_tottmp(ii)=z

                 !add some parts from "behind" the source
                 if (theta<=th2-thr .and. th1<=abs(th2-thr)) then 
                    ii=ii+1
                    xgd_tottmp(ii)=-xgd_tottmp(ii-1)
                    zgd_tottmp(ii)=zgd_tottmp(ii-1)
                 endif 

              endif
           enddo
        enddo
     enddo
  enddo
  if (mynum==0) write(6,*)'Total grid points in SEM mesh:',npts_tot
  npts_tot=ii
  allocate(xgd_tot(1:npts_tot),zgd_tot(npts_tot))
  if (mynum==0) then
     write(6,*)
     write(6,*)'Total kernel grid points:',npts_tot
     write(6,*)'Kernel procs:',nproc
     write(6,*)'SEM procs:',nproc_mesh
  endif

  xgd_tot = xgd_tottmp(1:npts_tot);  zgd_tot = zgd_tottmp(1:npts_tot)
  deallocate(xgd_tottmp,zgd_tottmp)
  if (mynum==0) then
     open(unit=9999,file=ext_mesh_name(1:lfext))
     write(9999,*)npts_tot
     do ii=1,npts_tot
        write(9999,*)xgd_tot(ii),y,zgd_tot(ii)
     enddo
     close(9999)
  endif 
  if (nproc>1) call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
  deallocate(xgd_tot,zgd_tot)
  call read_external_mesh
 
endif

end subroutine get_semkernel_fraction
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine get_semkernel_mesh2(mesh)
  
  use data_rota
  use mpi

  real(kind=realkind),intent(in) ::mesh(ibeg:iend,ibeg:iend,1:nelem,0:nproc_mesh-1,1:2)
  real(kind=realkind), dimension(:), allocatable :: sgd_tot
  integer :: ipt,ipol,jpol,iel, iproc, ipt_in_proc

  if (mynum==0) write(6,*)'define kernel mesh upon initial SEM mesh....'

  npts_tot=2*nproc_mesh*nelem*(iend-ibeg+1)**2
  allocate(sgd_tot(npts_tot),xgd_tot(1:npts_tot))
  allocate(ygd_tot(1:npts_tot),zgd_tot(npts_tot))
  allocate(iproc_proc_tot(npts_tot),iel_proc_tot(npts_tot))
  allocate(ipt_proc_tot(npts_tot))

  if (mynum==0) then
     write(6,*)
     write(6,*)'Total grid points:',npts_tot
     write(6,*)'Kernel procs:',nproc
     write(6,*)'SEM procs:',nproc_mesh
  endif

  sgd_tot = 0.; zgd_tot = 0.; iproc_proc_tot = 0; ipt=0; 
  do iproc=0, nproc_mesh-1
     ipt_in_proc = 0
     do iel = 1, nelem
        do jpol=ibeg,iend
           do ipol=ibeg,iend
              ipt=ipt+1
              ipt_in_proc = ipt_in_proc + 1
              sgd_tot(ipt) = mesh(ipol,jpol,iel,iproc,1)
              zgd_tot(ipt) = mesh(ipol,jpol,iel,iproc,2)

            ! define fwd mapping at the global level
              iproc_proc_tot(ipt) = iproc
              iel_proc_tot(ipt) = iel
              ipt_proc_tot(ipt) = ipt_in_proc
           enddo
        enddo
     enddo
  enddo

  xgd_tot(1:npts_tot/2) = sgd_tot(1:npts_tot/2)
  ygd_tot = 0. ! this is only to avoid rotated domains when searching s,z
               ! later, we will use phi0 to compute things at the right y.

! mirror the D-shaped domain to phi=pi
  xgd_tot(npts_tot/2+1:npts_tot) = -xgd_tot(1:npts_tot/2)
  zgd_tot(npts_tot/2+1:npts_tot) = zgd_tot(1:npts_tot/2)

  deallocate(sgd_tot)

end subroutine get_semkernel_mesh2
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine compute_3dslices(mesh)
  
  use data_rota
  use coord_trafo, only : sphi2xy_arr,sphi2xy,rthetaphi2xyz
  use mpi

  real(kind=realkind),intent(in) ::mesh(0:iend,1:iend,1:nelem,0:nproc_mesh-1,1:2)
  real(kind=realkind), dimension(:), allocatable :: sgd_tot2,xgd_tot2,ygd_tot2,zgd_tot2,vp,vptop,vpbot
  real(kind=realkind), dimension(:), allocatable :: vptot
  integer :: ipt,ipol,jpol,iel, iproc, ipt_in_proc,nphi,ntheta,nsim,k1,k2,npts_slice
  integer :: npts_top,npts_bot,j,isim,snap1,snap2,skipfact,i,startind,snapskip
  integer, dimension(:), allocatable :: ind_top,ind_bot
  real(kind=realkind)  :: rbot,rtop,dtheta,phi00,phi_incr,dphi,r,smallval_bot,smallval_top,theta_meri,r0
  character(len=50) :: filename1
  character(len=34), allocatable :: simdir(:)
  logical :: use_top,use_bot,use_meri

  if (mynum==0) write(6,*)'define kernel mesh upon initial SEM mesh....'

open(unit=99,file='param_snaps')
read(99,*)phi00
write(6,*)'starting azimuth/phi for cross section on the right [deg]:',phi00
read(99,*)dphi
write(6,*)'ending azimuth/phi for cross section on the left [deg]:',dphi
read(99,*)rtop
write(6,*)'top surface [km]:',rtop
read(99,*)rbot
write(6,*)'bottom surface [km]:',rbot
read(99,*)theta_meri
write(6,*)'colatitude of meridional cross section:',theta_meri
read(99,*)snap1,snap2,snapskip
write(6,*)'1st,last snap, skipfactor:',snap1,snap2,snapskip
read(99,*)use_meri
write(6,*)'consider meridional cross section?',use_meri
read(99,*)use_top
write(6,*)'consider top surface?',use_top
read(99,*)use_bot
write(6,*)'consider bottom surface?',use_bot
read(99,*)skipfact
write(6,*)'skipping factor in mesh:',skipfact
read(99,*)startind
write(6,*)'starting index:',startind
read(99,*)nsim
if (nsim==1) write(6,*)'no need to sum, one simulation only!'
allocate(simdir(nsim))
do isim=1,nsim
   read(99,*)simdir(isim)
   write(6,*)isim,'simulation directory: ',trim(simdir(isim))
enddo
close(99)

phi00=phi00/180.*pi; dphi=dphi/180.*pi; rtop=rtop*1000.; rbot=rbot*1000.
theta_meri=theta_meri*pi/180.

npts_tot=nproc_mesh*nelem*10
npts_slice=nproc_mesh*nelem

write(6,*)'preliminary npts_tot:',npts_tot

write(6,*)mynum,'allocating preliminary global arrays...'
  allocate(sgd_tot2(npts_tot),xgd_tot2(1:npts_tot))
  allocate(ygd_tot2(1:npts_tot),zgd_tot2(npts_tot))

  if (mynum==0) then
     write(6,*)
     write(6,*)'Total grid points:',npts_tot
     write(6,*)'Kernel procs:',nproc
     write(6,*)'SEM procs:',nproc_mesh
  endif

  sgd_tot2 = 0.; zgd_tot2 = 0.; ipt=0; 
  do iproc=0, nproc_mesh-1
     ipt_in_proc = 0
     do iel = 1, nelem
              ipt=ipt+1
              ipt_in_proc = ipt_in_proc + 1
              sgd_tot2(ipt) = mesh(2,2,iel,iproc,1)
              zgd_tot2(ipt) = mesh(2,2,iel,iproc,2)
     enddo
  enddo

! first slice
write(6,*)mynum,'computing first slice...'
call sphi2xy_arr(xgd_tot2(1:npts_slice),ygd_tot2(1:npts_slice),sgd_tot2(1:npts_slice),phi00,npts_slice)

! second slide
write(6,*)mynum,'computing second slice...'
call sphi2xy_arr(xgd_tot2(npts_slice+1:2*npts_slice), &
                      ygd_tot2(npts_slice+1:2*npts_slice),sgd_tot2,dphi+phi00,npts_slice)
zgd_tot2(npts_slice+1:2*npts_slice)=zgd_tot2(1:npts_slice)

! load p wave velocity
write(6,*)mynum,'loading p wave velocity...'
allocate(vp(2*npts_slice))
do i=1,2*npts_slice  
   r= sqrt(xgd_tot2(i)**2+ygd_tot2(i)**2+zgd_tot2(i)**2)
   vp(i)=prem(r,'v_p')
enddo

! rescale vp
vp=vp/maxval(vp)*0.2-0.1

filename1='mesh_xsect'
!call write_VTK_bin_scal(xgd_tot2(1:2*npts_slice),ygd_tot2(1:2*npts_slice),&
!     zgd_tot2(1:2*npts_slice), vp(1:2*npts_slice),2*npts_slice,0,filename1)
!call write_avs_file_scal(xgd_tot2(1:2*npts_slice),ygd_tot2(1:2*npts_slice),&
!     zgd_tot2(1:2*npts_slice), vp(1:2*npts_slice),2*npts_slice,0,filename1)

! rtop & rbot
write(6,*)mynum,'computing top and bottom surfaces...'
allocate(ind_top(floor(real(npts_slice)/10.)),ind_bot(floor(real(npts_slice)/10.)))
smallval_top = minval(abs(sqrt(sgd_tot2(1:npts_slice)**2+zgd_tot2(1:npts_slice)**2)-rtop))*1.001
smallval_bot = minval(abs(sqrt(sgd_tot2(1:npts_slice)**2+zgd_tot2(1:npts_slice)**2)-rbot))*1.001
write(6,*)'smallval top:',smallval_top
write(6,*)'smallval bot:',smallval_bot

ipt=0; k1=0;k2=0
do iproc=0, nproc_mesh-1
   ipt_in_proc = 0
   do iel = 1, nelem
      ipt=ipt+1
      r0=sqrt(sgd_tot2(ipt)**2+zgd_tot2(ipt)**2)

      if ( abs(r0-rtop) < smallval_top ) then ! rtop
         k1=k1+1 
         ind_top(k1)=ipt
         write(61,*)iel,r0,k1
      endif
      
      if ( abs(r0-rbot) < smallval_bot ) then ! rbot
         k2=k2+1 
         ind_bot(k2)=ipt
      endif
   enddo
enddo

npts_top=k1
npts_bot=k2
write(6,*)'points on rtop (input mesh):',npts_top
write(6,*)'points on rbot (input mesh):',npts_bot

k1=0
! top surface
write(6,*)mynum,'top!'
dtheta=rtop*pi/npts_top
write(6,*)mynum,'dtheta [km] dphi [deg] rtop:',dtheta/1000.,dphi*180./pi
do i=1,npts_top,4
   nphi = max(floor(sgd_tot2(ind_top(i))*(2*pi-(dphi-phi00))/dtheta/4.),1)
   phi_incr = (2.*pi-(dphi-phi00))/nphi
!   write(6,*)mynum,'nphi, phi_incr:',nphi,phi_incr*180./pi
   do j=1,nphi
      k1=k1+1
      write(67,*)k1,2*npts_slice+k1,npts_tot
      xgd_tot2(2*npts_slice+k1) = sgd_tot2(ind_top(i))*cos(phi00-(j-1)*phi_incr)
      ygd_tot2(2*npts_slice+k1) = sgd_tot2(ind_top(i))*sin(phi00-(j-1)*phi_incr)
      zgd_tot2(2*npts_slice+k1) = zgd_tot2(ind_top(i))
   enddo
enddo

! extract vp
allocate(vptop(k1))
vptop(1:k1) = vp(ind_top(1))

! save into vtk
filename1='mesh_top'
!call write_VTK_bin_scal(xgd_tot2(2*npts_slice+1:2*npts_slice+k1) ,&
!                                           ygd_tot2(2*npts_slice+1:2*npts_slice+k1), &
!                                           zgd_tot2(2*npts_slice+1:2*npts_slice+k1),vptop,k1,0,filename1)
!
!call write_avs_file_scal(xgd_tot2(2*npts_slice+1:2*npts_slice+k1) ,&
!                                           ygd_tot2(2*npts_slice+1:2*npts_slice+k1), &
!                                           zgd_tot2(2*npts_slice+1:2*npts_slice+k1),vptop,k1,0,filename1)
! bottom surface
k2=0
write(6,*)mynum,'bottom!'
dtheta=rbot*pi/npts_bot
write(6,*)mynum,'dtheta [km] dphi [deg] rbot:',dtheta/1000.,dphi*180./pi
do i=1,npts_bot
   nphi = max(floor(sgd_tot2(ind_bot(i))*(dphi-phi00)/dtheta/1.),1)
   phi_incr = (dphi-phi00)/nphi
   do j=1,nphi
      k2=k2+1
      write(68,*)k2,2*npts_slice+k1+k2,npts_tot
      xgd_tot2(2*npts_slice+k1+k2) = sgd_tot2(ind_bot(i))*cos(phi00+(j-1)*phi_incr)
      ygd_tot2(2*npts_slice+k1+k2) = sgd_tot2(ind_bot(i))*sin(phi00+(j-1)*phi_incr)
      zgd_tot2(2*npts_slice+k1+k2) = zgd_tot2(ind_bot(i))
   enddo
enddo

! extract vp
allocate(vpbot(k2))
!vpbot(1:k2) = vp(ind_bot(1))
vpbot(1:k2) = minval(vp)*0.95

! save into vtk
filename1='mesh_bot'
!call write_VTK_bin_scal(xgd_tot2(2*npts_slice+k1+1:2*npts_slice+k1+k2) , &
!                                           ygd_tot2(2*npts_slice+k1+1:2*npts_slice+k1+k2), &
!                                           zgd_tot2(2*npts_slice+k1+1:2*npts_slice+k1+k2),vpbot,k2,0,filename1)

!call write_avs_file_scal(xgd_tot2(2*npts_slice+k1+1:2*npts_slice+k1+k2) , &
!                                           ygd_tot2(2*npts_slice+k1+1:2*npts_slice+k1+k2), &
!                                           zgd_tot2(2*npts_slice+k1+1:2*npts_slice+k1+k2),vpbot,k2,0,filename1)

deallocate(sgd_tot2)

write(6,*)mynum,'points in top surface:',k1
write(6,*)mynum,'points in bottom surface:',k2

write(6,*)mynum,'total points (preliminary):',npts_tot
npts_tot=2*npts_slice+k1+k2
write(6,*)mynum,'total points (actual):',npts_tot

allocate(xgd_tot(1:npts_tot),ygd_tot(1:npts_tot),zgd_tot(npts_tot))
xgd_tot(1:npts_tot)=xgd_tot2(1:npts_tot)
ygd_tot(1:npts_tot)=ygd_tot2(1:npts_tot)
zgd_tot(1:npts_tot)=zgd_tot2(1:npts_tot)

deallocate(xgd_tot2,ygd_tot2,zgd_tot2)

allocate(vptot(1:npts_tot))
vptot(1:2*nelem*nproc_mesh) = vp
vptot(2*nelem*nproc_mesh+1:2*nelem*nproc_mesh+k1)=vptop
vptot(2*nelem*nproc_mesh+k1+1:2*nelem*nproc_mesh+k1+k2)=vpbot

! save into AVS
filename1='Data/mesh_all'
!call write_VTK_bin_scal(xgd_tot,ygd_tot,zgd_tot,vptot,npts_tot,0,filename1)
!call write_avs_file_scal(xgd_tot,ygd_tot,zgd_tot,vptot,npts_tot,0,filename1)

deallocate(vptot,vp,vpbot,vptop)

end subroutine compute_3dslices
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine compute_3dslices_4(mesh,phr)
  
  use data_rota
  use coord_trafo, only : sphi2xy_arr,sphi2xy,rthetaphi2xyz
  use mpi
  real(kind=realkind), intent(in) :: phr
  real(kind=realkind), intent(in) ::mesh(0:iend,1:iend,1:nelem,0:nproc_mesh-1,1:2)
  real(kind=realkind), dimension(:), allocatable :: sgd_tot2,xgd_tot2,ygd_tot2,zgd_tot2
  real(kind=realkind), dimension(:), allocatable :: vp,vptop,vpbot
  real(kind=realkind), dimension(:), allocatable :: vptot
  integer :: ipt,ipol,jpol,iel, iproc, ipt_in_proc,nphi,ntheta,nsim,k1,k2,npts_slice
  integer :: npts_top,npts_bot,j,isim,snap1,snap2,skipfact,i,startind,snapskip
  integer, dimension(:), allocatable :: ind_top,ind_bot
  real(kind=realkind)  :: rbot,rtop,dtheta,phi00,phi_incr,dphi,r,smallval_bot,smallval_top,theta_meri,r0
  character(len=50) :: filename1
  character(len=34), allocatable :: simdir(:)
  logical :: use_top,use_bot,use_meri

  if (mynum==0) write(6,*)'define kernel mesh upon initial SEM mesh....'

open(unit=99,file='param_snaps')
read(99,*)phi00
write(6,*)'starting azimuth/phi for cross section on the right [deg]:',phi00
read(99,*)dphi
write(6,*)'ending azimuth/phi for cross section on the left [deg]:',dphi
read(99,*)rtop
write(6,*)'top surface [km]:',rtop
read(99,*)rbot
write(6,*)'bottom surface [km]:',rbot
read(99,*)theta_meri
write(6,*)'colatitude of meridional cross section:',theta_meri
read(99,*)snap1,snap2,snapskip
write(6,*)'1st,last snap, skipfactor:',snap1,snap2,snapskip
read(99,*)use_meri
write(6,*)'consider meridional cross section?',use_meri
read(99,*)use_top
write(6,*)'consider top surface?',use_top
read(99,*)use_bot
write(6,*)'consider bottom surface?',use_bot
read(99,*)skipfact
write(6,*)'skipping factor in mesh:',skipfact
read(99,*)startind
write(6,*)'starting index:',startind
read(99,*)nsim
if (nsim==1) write(6,*)'no need to sum, one simulation only!'
allocate(simdir(nsim))
do isim=1,nsim
   read(99,*)simdir(isim)
   write(6,*)isim,'simulation directory: ',trim(simdir(isim))
enddo
close(99)

phi00=phr

phi00=phi00/180.*pi; dphi=dphi/180.*pi; rtop=rtop*1000.; rbot=rbot*1000.
theta_meri=theta_meri*pi/180.

npts_tot=nproc_mesh*nelem*10*(iend-ibeg+1)**2
npts_slice=nproc_mesh*nelem*(iend-ibeg+1)**2

write(6,*)'preliminary npts_tot:',npts_tot

write(6,*)mynum,'allocating preliminary global arrays...'
  allocate(sgd_tot2(npts_tot),xgd_tot2(1:npts_tot))
  allocate(ygd_tot2(1:npts_tot),zgd_tot2(npts_tot))

  if (mynum==0) then
     write(6,*)
     write(6,*)'Total grid points:',npts_tot
     write(6,*)'Kernel procs:',nproc
     write(6,*)'SEM procs:',nproc_mesh
  endif

  sgd_tot2 = 0.; zgd_tot2 = 0.; ipt=0; 
  do iproc=0, nproc_mesh-1
     ipt_in_proc = 0
     do iel = 1, nelem
        do i=ibeg,iend!,iend-ibeg
           do j=ibeg,iend!,iend-ibeg
!        i=2; j=2
              ipt=ipt+1
              ipt_in_proc = ipt_in_proc + 1
              sgd_tot2(ipt) = mesh(i,j,iel,iproc,1)
              zgd_tot2(ipt) = mesh(i,j,iel,iproc,2)
              write(888,*)i,j,iel,ipt
          enddo
       enddo
     enddo
  enddo

! first slice
write(6,*)mynum,'computing first slice...'
call sphi2xy_arr(xgd_tot2(1:npts_slice),ygd_tot2(1:npts_slice),sgd_tot2(1:npts_slice),phi00,npts_slice)

! second slide
write(6,*)mynum,'computing second slice...'
call sphi2xy_arr(xgd_tot2(npts_slice+1:2*npts_slice), &
                      ygd_tot2(npts_slice+1:2*npts_slice),sgd_tot2,dphi+phi00,npts_slice)
zgd_tot2(npts_slice+1:2*npts_slice)=zgd_tot2(1:npts_slice)

! load p wave velocity
write(6,*)mynum,'loading p wave velocity...'
allocate(vp(2*npts_slice))
do i=1,2*npts_slice  
   r= sqrt(xgd_tot2(i)**2+ygd_tot2(i)**2+zgd_tot2(i)**2)
   write(6,*) i, r, xgd_tot2(i), ygd_tot2(i), zgd_tot2(i), sgd_tot2(i)
   vp(i)=prem(r,'v_p')
enddo

! rescale vp
vp=vp/maxval(vp)*0.2-0.1

! rtop & rbot
write(6,*)mynum,'computing top and bottom surfaces...'
allocate(ind_top(floor(real(npts_slice)/1.)),ind_bot(floor(real(npts_slice)/1.)))
smallval_top = minval(abs(sqrt(sgd_tot2(1:npts_slice)**2+zgd_tot2(1:npts_slice)**2)-rtop))*1.05
smallval_bot = minval(abs(sqrt(sgd_tot2(1:npts_slice)**2+zgd_tot2(1:npts_slice)**2)-rbot))*1.05
write(6,*)'smallval top:',smallval_top
write(6,*)'smallval bot:',smallval_bot

ipt=0; k1=0;k2=0
do iproc=0, nproc_mesh-1
   ipt_in_proc = 0
   do iel = 1, nelem
      do i=ibeg,iend
           do j=ibeg,iend!,iend-ibeg
!      i=2;j=2
              ipt=ipt+1
              r0=sqrt(sgd_tot2(ipt)**2+zgd_tot2(ipt)**2)
              if ( abs(r0-rtop) < smallval_top ) then ! rtop
                 k1=k1+1 
                 ind_top(k1)=ipt
                 write(61,*)iel,r0,k1
              endif
      
              if ( abs(r0-rbot) < smallval_bot ) then ! rbot
                 k2=k2+1 
                 ind_bot(k2)=ipt
              endif
           enddo
        enddo
     enddo
enddo

npts_top=k1
npts_bot=k2
write(6,*)'points on rtop (input mesh):',npts_top
write(6,*)'points on rbot (input mesh):',npts_bot

k1=0
! top surface
write(6,*)mynum,'top!'
dtheta=rtop*pi/npts_top
write(6,*)mynum,'dtheta [km] dphi [deg] rtop:',dtheta/1000.,dphi*180./pi
do i=1,npts_top,1
   nphi = max(floor(sgd_tot2(ind_top(i))*(2.*pi-(dphi-phi00))/dtheta/1.),1)
   phi_incr = (2.*pi-(dphi-phi00))/nphi
   do j=1,nphi
      if (abs(phi00+(j-1)*phi_incr-phi00)< 40.*pi/180. )then
         k1=k1+1
         write(67,*)k1,2*npts_slice+k1,npts_tot
         xgd_tot2(2*npts_slice+k1) = sgd_tot2(ind_top(i))*cos(phi00-(j-1)*phi_incr)
         ygd_tot2(2*npts_slice+k1) = sgd_tot2(ind_top(i))*sin(phi00-(j-1)*phi_incr)
         zgd_tot2(2*npts_slice+k1) = zgd_tot2(ind_top(i))
      endif
   enddo
enddo

! extract vp
allocate(vptop(k1))
vptop(1:k1) = vp(ind_top(1))

! bottom surface
k2=0
write(6,*)mynum,'bottom!'
dtheta=rbot*pi/npts_bot
write(6,*)mynum,'dtheta [km] dphi [deg] rbot:',dtheta/1000.,dphi*180./pi
do i=1,npts_bot
   nphi = max(floor(sgd_tot2(ind_bot(i))*(dphi-phi00)/dtheta/1.),1)
   phi_incr = (dphi-phi00)/nphi
   do j=1,nphi
      k2=k2+1
      write(68,*)k2,2*npts_slice+k1+k2,npts_tot
      xgd_tot2(2*npts_slice+k1+k2) = sgd_tot2(ind_bot(i))*cos(phi00+(j-1)*phi_incr)
      ygd_tot2(2*npts_slice+k1+k2) = sgd_tot2(ind_bot(i))*sin(phi00+(j-1)*phi_incr)
      zgd_tot2(2*npts_slice+k1+k2) = zgd_tot2(ind_bot(i))
   enddo
enddo

! extract vp
allocate(vpbot(k2))
vpbot(1:k2) = minval(vp)*0.95

! save into vtk
filename1='mesh_bot'
deallocate(sgd_tot2)

write(6,*)mynum,'points in top surface:',k1
write(6,*)mynum,'points in bottom surface:',k2

write(6,*)mynum,'total points (preliminary):',npts_tot
npts_tot=2*npts_slice+k1+k2
write(6,*)mynum,'total points (actual):',npts_tot

allocate(xgd_tot(1:npts_tot),ygd_tot(1:npts_tot),zgd_tot(npts_tot))
xgd_tot(1:npts_tot)=xgd_tot2(1:npts_tot)
ygd_tot(1:npts_tot)=ygd_tot2(1:npts_tot)
zgd_tot(1:npts_tot)=zgd_tot2(1:npts_tot)

deallocate(xgd_tot2,ygd_tot2,zgd_tot2)

allocate(vptot(1:npts_tot))
vptot(1:2*npts_slice) = vp
vptot(2*npts_slice+1:2*npts_slice+k1)=vptop
vptot(2*npts_slice+k1+1:2*npts_slice+k1+k2)=vpbot

deallocate(vptot,vp,vpbot,vptop)


end subroutine compute_3dslices_4
!-------------------------------------------------------------------------


!=========================================================

subroutine compute_3dslices_cell4(mesh,phr)
  
  use data_rota
  use coord_trafo, only : sphi2xy_arr, sphi2xy, rthetaphi2xyz
  use input_output, only: write_VTK_bin_scal_mesh3d, write_VTK_bin_scal_topology
  use data_mesh, only : ind_cell2tot, valence_tot2cell, ind_tot2cell,npts
  use data_arrays, only : kern_cell
  use mpi

  real(kind=realkind),intent(in) :: phr
  real(kind=realkind), intent(in) :: mesh(ibeg:iend,ibeg:iend,1:nelem,0:nproc_mesh-1,1:2)
  integer iproc, npts_top, npts_bot, npts_meri, nphi, snapskip, snap1, snap2
  real, dimension(:,:), allocatable :: coord, disp, disptot, disptot_sum, azi_xsect2, azi_xsect1
  real, dimension(:,:,:), allocatable :: azi, fluid_prefact, azi_meri, mij_snap
  real(kind=realkind), dimension(:), allocatable :: vp, x, y, z, vptop, vpbot, vpmeri, xtot, ytot, ztot, azi_phi_meri, vptot
  double precision, dimension(:), allocatable :: xtop, ytop, ztop, xbot, ybot, zbot, azi_phi_top, azi_phi_bot
  real, dimension(:), allocatable :: xmeri, ymeri, zmeri, longit_snap
  real :: dphi, phi0, r, theta_meri, smallval_meri, dr, theta0, theta_bot, theta_top, phi
  integer, dimension(:), allocatable :: ind_proc_top_tmp, ind_pts_top_tmp, ind_proc_bot_tmp, ind_pts_bot_tmp
  integer, dimension(:), allocatable :: ind_proc_top, ind_pts_top, ind_proc_bot
  integer, dimension(:), allocatable :: ind_pts_bot, azi_ind_top, azi_ind_bot, azi_ind_meri
  integer, dimension(:), allocatable :: ind_proc_meri_tmp, ind_pts_meri_tmp, ind_proc_meri, ind_pts_meri
  real :: smallval_north_top, smallval_south_top, smallval_north_bot, smallval_south_bot, r0, coord9(9,2), disp9(9,3)
  character(len=200) :: filename1
  logical :: use_meri, use_top, use_bot
  character(len=4) :: appmynum2, appiproc
  integer :: nptstot, k1, k2, npts_fluid, skipfact, j, naang, i, iii, iel
  real :: rbot, rtop, dtheta, phi_incr, smallval
  double precision :: r8, theta8, phi8

  skipfact = 2
  smallval = 10000.

  ! read snap plot parameters
  open(unit=99,file='param_snaps')
  read(99,*) phi0
  write(6,*) 'starting azimuth/phi for cross section on the left [deg]:', phi0
  read(99,*) dphi
  write(6,*) 'ending azimuth/phi for cross section on the right [deg]:', dphi
  read(99,*) rtop
  write(6,*) 'top surface [km]:',rtop
  read(99,*) rbot
  write(6,*) 'bottom surface [km]:', rbot
  read(99,*) theta_meri
  write(6,*) 'colatitude of meridional cross section:', theta_meri
  read(99,*) snap1,snap2,snapskip
  write(6,*) '1st,last snap, skipfactor:', snap1, snap2, snapskip
  read(99,*) use_meri
  write(6,*) 'consider meridional cross section?', use_meri
  read(99,*) use_top
  write(6,*) 'consider top surface?', use_top
  read(99,*) use_bot
  write(6,*) 'consider bottom surface?', use_bot
  close(99)

  phi0 = phi0 / 180. * pi
  dphi = dphi / 180. * pi
  rtop = rtop * 1000.
  rbot = rbot * 1000.
  theta_meri = theta_meri * pi / 180.

  npts = nelem * 4
  nptstot = npts * nproc_mesh
  write(6,*) 'number of points per proc, total points:', npts, nptstot

     ! load and construct global mesh (one semi-disk)
     write(6,*)'reading partitioned mesh...'
     allocate(coord(nptstot,2))
     smallval_north_top = rtop
     smallval_south_top = rtop
     smallval_north_bot = rbot
     smallval_south_bot = rbot
     write(6,*)'Smallest distance to rtop (North,South) BEFORE [km]:', &
          real(smallval_north_top/1000.),real(smallval_south_top/1000.)
     write(6,*)'Smallest distance to rbot (North,South) BEFORE [km]:', &
          real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)
     do iproc=0, nproc_mesh-1
        i=0
        do iel=1, nelem
           iii = 0
           coord(iproc*npts+i+1,:) = mesh(ibeg,ibeg,iel,iproc,1:2)
           coord(iproc*npts+i+2,:) = mesh(ibeg,iend,iel,iproc,1:2)
           coord(iproc*npts+i+3,:) = mesh(iend,iend,iel,iproc,1:2)
           coord(iproc*npts+i+4,:) = mesh(iend,ibeg,iel,iproc,1:2)      

           ! determine minimal distance from rtop and rbot
            do iii=1, 4
               r0 = sqrt(coord(iproc*npts+i+iii,1)**2 + coord(iproc*npts+i+iii,2)**2) 
               if (coord(iproc*npts+i+iii,2) > 0.0 .and. abs(r0-rtop) < smallval_north_top) then ! north
                  smallval_north_top = abs(r0-rtop)
               elseif (coord(iproc*npts+i+iii,2) < 0.0 .and. abs(r0-rtop) < smallval_south_top) then ! south
                  smallval_south_top = abs(r0-rtop)
               endif
               
               if (coord(iproc*npts+i+iii,2)>0.0 .and. abs(r0-rbot)< smallval_north_bot) then ! north            
                  smallval_north_bot = abs(r0-rbot)
               elseif (coord(iproc*npts+i+iii,2)<0.0 .and. abs(r0 -rbot)< smallval_south_bot) then ! south
                  smallval_south_bot = abs(r0-rbot)
               endif
           enddo
            i = i + 4
        enddo
     enddo

     smallval_north_top = smallval_north_top + smallval
     smallval_south_top = smallval_south_top + smallval
     smallval_north_bot = smallval_north_bot + smallval
     smallval_south_bot = smallval_south_bot + smallval 
     
     write(6,*)'Smallest distance to rtop (North,South) AFTER [km]:', &
          real(smallval_north_top/1000.),real(smallval_south_top/1000.)
     write(6,*)'Smallest distance to rbot (North,South) AFTER [km]:', &
          real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)
     
     if (use_top ) allocate(ind_proc_top_tmp(floor(real(nptstot)/10.)), &
                            ind_pts_top_tmp(floor(real(nptstot)/10.)))
     if (use_bot ) allocate(ind_proc_bot_tmp(floor(real(nptstot)/10.)), &
                            ind_pts_bot_tmp(floor(real(nptstot)/10.)))
     
     k1 = 0
     k2 = 0 
     do iproc=0, nproc_mesh-1
        do i=1, npts, skipfact
           ! check for top and bottom radii
           r0 = sqrt(coord(iproc*npts+i,1)**2 + coord(iproc*npts+i,2)**2) 

           if (use_top ) then
              if ( (coord(iproc*npts+i,2) >= 0.0 .and.  &
                   abs(r0-rtop) <= smallval_north_top) .or. &
                   (coord(iproc*npts+i,2) < 0.0 .and.  &
                   abs(r0-rtop) <= smallval_south_top)) then 
                 k1 = k1 + 1         
                 ind_proc_top_tmp(k1) = iproc
                 ind_pts_top_tmp(k1) = i
              endif
           endif

           if (use_bot) then 
              if ( (coord(iproc*npts+i,2) >= 0.0 .and.  &
                   abs(r0-rbot) <= smallval_north_bot) .or. &
                   (coord(iproc*npts+i,2) < 0.0 .and.  &
                   abs(r0-rbot) <= smallval_south_bot)) then 
                 k2 = k2 + 1
                 ind_proc_bot_tmp(k2) = iproc
                 ind_pts_bot_tmp(k2) = i
              endif
           endif
        enddo
     enddo

     npts_top = k1
     npts_bot=k2
     write(6,*) '# points on top,bottom:', npts_top, npts_bot
     write(6,*) 'allocating index arrays for surfaces....'

     if (use_top ) then
        allocate(ind_proc_top(npts_top), ind_pts_top(npts_top))
        ind_proc_top = ind_proc_top_tmp(1:npts_top)
        ind_pts_top = ind_pts_top_tmp(1:npts_top)
        deallocate(ind_proc_top_tmp,ind_pts_top_tmp)
     endif

     if (use_bot) then
        allocate(ind_proc_bot(npts_bot), ind_pts_bot(npts_bot))
        ind_proc_bot = ind_proc_bot_tmp(1:npts_bot)
        ind_pts_bot = ind_pts_bot_tmp(1:npts_bot)
        deallocate(ind_proc_bot_tmp, ind_pts_bot_tmp)
     endif

     ! xyz coordinates---------------------------------------------------------------------------------------------
     write(6,*)'defining xyz...'
     allocate(x(1:2*nptstot), y(1:2*nptstot), z(1:2*nptstot))
     x = 0.
     y = 0.
     z = 0.

     ! left cross section---------------------------------------------------------------------------------------------
     call sphi2xy_arr(x(1:nptstot), y(1:nptstot), coord(:,1), phi0, nptstot)
     z(1:nptstot) = coord(1:nptstot,2)
     write(6,*) 'max s,x,y:', maxval(coord(:,1)), maxval(x), maxval(y)

     ! right cross section---------------------------------------------------------------------------------------------
     call sphi2xy_arr(x(nptstot+1:2*nptstot), y(nptstot+1:2*nptstot), coord(:,1), phi0+dphi, nptstot)
     z(nptstot+1:2*nptstot) = coord(1:nptstot,2)
     write(6,*) 'max s,x,y:', maxval(coord(:,1)), maxval(x), maxval(y)

     allocate(vp(2*nptstot))
     do i=1, 2*nptstot  
        r = sqrt(x(i)**2 + y(i)**2 + z(i)**2)
        vp(i) = prem(r,'v_p')
     enddo
    
     !XXX Why this rescaling??
     vp = vp / maxval(vp) * 0.002 - 0.001 ! rescale vp
     filename1 = 'Data/mesh_xsect'
     call write_VTK_bin_scal_mesh3d(x, y, z, vp, 2*nptstot, 0, filename1)
     filename1= 'Data/mesh_xsect_cell'
     call write_VTK_bin_scal_topology(x, y, z, vp, 2*nptstot/4, filename1)

     ! top surface---------------------------------------------------------------------------------------------
     if (use_top) then
        k1=0
        write(6,*)'defining top surface...'
        dtheta=rtop*pi/npts_top
        naang = floor(real(npts_top)/real(2.))**2*6*4
        write(6,*)'points on top surface:',naang
        allocate(xtop(1:naang),ytop(1:naang),ztop(1:naang))
        allocate(azi_phi_top(1:naang),azi_ind_top(1:naang))
        call construct_surface_cubed_sphere(npts_top,npts,dble(rtop),ind_proc_top,ind_pts_top,nptstot,dble(coord),k1, &
                                  dble(dphi),dble(phi0),'outside',naang,xtop,ytop,ztop,azi_phi_top,azi_ind_top)
        write(6,*)'allocating vptop...',k1
        allocate(vptop(k1))
        vptop(1:k1) = vp(ind_proc_top(1)*npts+ind_pts_top(1))
        write(6,*)'save into cell vtk...',size(xtop),k1
        filename1='Data/mesh_top_cell'
        call write_VTK_bin_scal_topology(real(xtop(1:k1)),real(ytop(1:k1)),real(ztop(1:k1)),vptop(1:k1),k1/4,filename1)
     endif

     ! bottom surface ---------------------------------------------------------------------------------------------
     if (use_bot) then
        k2=0
        write(6,*)'defining bot surface...'
        dtheta=rbot*pi/npts_bot
        naang = floor(real(npts_bot)/real(2.))**2*6*4
        write(6,*)'points on bottom surface:',naang
        allocate(xbot(1:naang),ybot(1:naang),zbot(1:naang))
        allocate(azi_phi_bot(1:naang),azi_ind_bot(1:naang))
        call construct_surface_cubed_sphere(npts_bot,npts,dble(rbot),ind_proc_bot,ind_pts_bot,nptstot,dble(coord),k2, &
                                  dble(dphi),dble(phi0),'innside',naang,xbot,ybot,zbot,azi_phi_bot,azi_ind_bot)
        allocate(vpbot(k2))
        vpbot(1:k2) = vp(ind_proc_bot(1)*npts+ind_pts_bot(1))
        filename1='Data/mesh_bot_cell'
        call write_VTK_bin_scal_topology(real(xbot(1:k2)),real(ybot(1:k2)),real(zbot(1:k2)),vpbot,k2/4,filename1)
     endif
     deallocate(coord)

!=========================================================================
     write(6,*)'assembling different surfaces to one coordinate array...'
     npts_cell = 2 * nptstot + k1 + k2
     write(6,*)'allocating global velocity array...'
     allocate(vptot(npts_cell))
     vptot(1:2*nptstot) = vp
     if (use_top) vptot(2*nptstot+1:2*nptstot+k1) = vptop
     if (use_bot) vptot(2*nptstot+k1+1:2*nptstot+k1+k2) = vpbot
     deallocate(vp)

     write(6,*)'allocating global coordinates xgd...'
     allocate(xgd_cell(1:npts_cell), ygd_cell(1:npts_cell), zgd_cell(1:npts_cell))

     write(6,*)'Filling xgd_cell...'
     xgd_cell(1:2*nptstot) = x(1:2*nptstot)
     ygd_cell(1:2*nptstot) = y(1:2*nptstot)
     zgd_cell(1:2*nptstot) = z(1:2*nptstot)
     
     write(6,*)'  ... done with cross sections ...'

     if (use_top) then 
        xgd_cell(2*nptstot+1:2*nptstot+k1) = real(xtop(1:k1))
        ygd_cell(2*nptstot+1:2*nptstot+k1) = real(ytop(1:k1))
        zgd_cell(2*nptstot+1:2*nptstot+k1) = real(ztop(1:k1))
     endif

     if (use_bot) then 
        xgd_cell(2*nptstot+k1+1:2*nptstot+k1+k2) = real(xbot(1:k2))
        ygd_cell(2*nptstot+k1+1:2*nptstot+k1+k2) = real(ybot(1:k2))
        zgd_cell(2*nptstot+k1+1:2*nptstot+k1+k2) = real(zbot(1:k2))
     endif

     cell_topology = .true.

     filename1 = 'Data/mesh_all_cell'
     call write_VTK_bin_scal_topology(real(xgd_cell), real(ygd_cell), real(zgd_cell), &
                                      vptot, npts_cell/4, filename1)

! reducing overlapping coordinates: double loop
     allocate(ind_tot2cell(1:npts_cell,10), ind_cell2tot(1:npts_cell), &
              valence_tot2cell(1:npts_cell))
     call reduce_coord_overlap(npts_cell, xgd_cell, ygd_cell, zgd_cell,npts, & ! THIS NOW INCLUDES DOMAIN DECOMPOSITION!!
                               ind_tot2cell, ind_cell2tot, valence_tot2cell, nloc)
     call define_io_appendix(appiproc, mynum)

     write(6,*) 'allocating xgd, size:', npts
     allocate(xgd(1:npts), ygd(1:npts), zgd(1:npts))
     do i=1, nloc
        xgd(ind_cell2tot(i)) = xgd_cell(i)
        ygd(ind_cell2tot(i)) = ygd_cell(i)
        zgd(ind_cell2tot(i)) = zgd_cell(i)
     enddo

     deallocate(xgd_cell, ygd_cell, zgd_cell)
     allocate(xgd_cell(1:npts_cell), ygd_cell(1:npts_cell), zgd_cell(1:npts_cell))

     filename1 = 'Data/mesh_all_totpts_'//appiproc
     call write_VTK_bin_scal_mesh3d(xgd, ygd, zgd, zgd, npts, 0, filename1)

! test the mapping cell to total and back ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     xgd_cell = 0.
     ygd_cell = 0.
     zgd_cell = 0.

     write(6,*) 'max ind:', maxval(ind_tot2cell), maxval(ind_cell2tot)
     do i=1, npts
        do j=1, valence_tot2cell(i)
           xgd_cell(ind_tot2cell(i,j)) = xgd(i)
           ygd_cell(ind_tot2cell(i,j)) = ygd(i)
           zgd_cell(ind_tot2cell(i,j)) = zgd(i)
        enddo
     enddo
     filename1 = 'Data/mesh_all_cell_retry_'//appiproc
     call write_VTK_bin_scal_topology(real(xgd_cell(1:nloc)), real(ygd_cell(1:nloc)), &
                                      real(zgd_cell(1:nloc)), vptot(1:nloc), &
                                      nloc/4, filename1)
! end test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
     write(6,*) mynum, 'Done with 3D slices in 4-corner cell topology.'
     allocate(kern_cell(nloc))

end subroutine compute_3dslices_cell4
!---------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------
subroutine reduce_coord_overlap(n, x, y, z, n_tot, ind_tot2cell, ind_cell2tot, &
                                valence_tot2cell, nloc)

    use input_output, only : write_VTK_bin_scal_topology
    use mpi
    
    integer, intent(in)                                 :: n
    real(kind=realkind), dimension(n), intent(inout)    :: x, y, z

    integer, intent(out)                                :: n_tot, nloc
    integer, dimension(n), intent(out)                  :: valence_tot2cell, ind_cell2tot
    integer, dimension(n,10), intent(out)               :: ind_tot2cell

    real(kind=realkind), dimension(:), allocatable      :: xloc, yloc, zloc
    real(kind=realkind)                                 :: tiny
    integer                                             :: ipt, icount, ival, &
                                                           ii, jj, i, npadd
    logical                                             :: already_found(n)
    character(len=200)                                  :: filename1
    character(len=4)                                    :: appiproc
    
    if (mynum == 0) write(6,*) 'Decomposing the cell mesh....'

    ! Needs to be multiple of 4, because of implicit connectivity in cell output
    if (mod(n, nproc * 4) /= 0) then 
        if (mynum < nproc-1) then 
           nloc = n / (nproc * 4) * 4
        else
           npadd = mod(n, nproc * 4)
           nloc = n / (nproc * 4) * 4 + npadd
        endif
    else
       nloc = n / nproc
    endif
    
    if (nproc > 1) CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    write(6,*) mynum, 'reporting: total cell points:', n
    write(6,*) mynum, 'has', nloc, ' cell points'
    write(6,*) mynum, 'mynum*n/nproc: ', mynum*n/(nproc*4)*4
    
    allocate(xloc(nloc), yloc(nloc), zloc(nloc))
    xloc = 0.
    yloc = 0.
    zloc = 0.

    do ipt=1, n / (nproc * 4) * 4
       xloc(ipt) = x(mynum*n/(nproc*4)*4+ipt)
       yloc(ipt) = y(mynum*n/(nproc*4)*4+ipt)
       zloc(ipt) = z(mynum*n/(nproc*4)*4+ipt)
    enddo
    
    if (mynum == nproc-1 .and. mod(n,nproc * 4) /= 0) then
       xloc(n/(nproc*4)*4+1:nloc) = x(n-npadd+1:n)
       yloc(n/(nproc*4)*4+1:nloc) = y(n-npadd+1:n)
       zloc(n/(nproc*4)*4+1:nloc) = z(n-npadd+1:n)
    endif
    
    write(6,*) 'minmax xgd:', mynum, minval(xloc), maxval(xloc)
    write(6,*) 'minmax xgd_tot:', mynum, minval(x), maxval(x)
    write(6,*) 'minmax ygd:', mynum, minval(yloc), maxval(yloc)
    write(6,*) 'minmax ygd_tot:', mynum, minval(y), maxval(y)
    write(6,*) 'minmax zgd:', mynum, minval(zloc), maxval(zloc)
    write(6,*) 'minmax zgd_tot:', mynum, minval(z), maxval(z)
    
    if (nproc > 1) CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    if (mynum == 0) write(6,*) 'finished the domain decomposition of the kernel mesh.'
    
    call define_io_appendix(appiproc,mynum)
    filename1 = 'Data/mesh_all_cell_'//appiproc    
    call write_VTK_bin_scal_topology(real(xloc), real(yloc), real(zloc), zloc, &
                                     nloc/4, filename1)
    
    !===============================================================================
    write(6,*) mynum, 'eliminating overlapping points in kernel grid....'
    !XXX: MvD - is it smart to do this after loadbalancing?
    !               -> each processor just checks its own points
    !               -> leeds to slightly uneven distibution of points
    !               -> same point could still be in two different processes

    tiny = maxval(sqrt(xloc**2 + yloc**2 + zloc**2)) * smallval
    ii = 0
    jj = 0
    icount = 0
    already_found = .false.
    valence_tot2cell = -5
    ind_tot2cell = -5
    ind_cell2tot = -5
    
    do ipt=1, nloc
       ival = 0
       if (mod(ipt,nloc/10) == 0 .and. mynum == 0) &
          write(6,"(i2, '0% done in finding overlapping points')") 10 * ipt / nloc
       if (.not. already_found(ipt) ) then
          icount = icount + 1
          ind_tot2cell(icount,1) = ipt
          ind_cell2tot(ipt) = icount
          valence_tot2cell(icount) = 1
          if (ipt < nloc) then
             do i=ipt, nloc
                if ( ((xloc(ipt) - xloc(i))**2 + (yloc(ipt) - yloc(i))**2 + &
                        (zloc(ipt) - zloc(i))**2)  < tiny**2 ) then 
                   ival = ival + 1
                   already_found(i) = .true.
                   ind_tot2cell(icount,ival) = i
                   valence_tot2cell(icount) = ival
                   ind_cell2tot(i) = icount
                endif
             enddo
          endif
       endif
    enddo
    
    n_tot = icount
    write(6,*) mynum, ' REDUCE: original/reduced # points:', nloc, n_tot
    write(6,*) mynum, ' REDUCE: min/max valence:', minval(valence_tot2cell(1:n_tot)), &
               maxval(valence_tot2cell(1:n_tot))
    
    x(nloc+1:n) = 0.
    y(nloc+1:n) = 0.
    z(nloc+1:n) = 0.
    x(1:nloc) = xloc
    y(1:nloc) = yloc
    z(1:nloc) = zloc
    
    npts_cell = nloc

end subroutine reduce_coord_overlap
!-----------------------------------------------------------------------------------------


!=========================================================

subroutine compute_3dslices_cell16(mesh,phr)
  
  use data_rota
  use coord_trafo, only : sphi2xy_arr,sphi2xy,rthetaphi2xyz
  use input_output, only: write_VTK_bin_scal_mesh3d,write_VTK_bin_scal_topology
  use data_mesh, only : ind_cell2tot,valence_tot2cell,ind_tot2cell,npts
  use data_arrays, only : kern_cell
  use mpi

  real(kind=realkind),intent(in) :: phr
  real(kind=realkind), intent(in) ::mesh(ibeg:iend,ibeg:iend,1:nelem,0:nproc_mesh-1,1:2)
  integer iproc,npts_top,npts_bot,npts_meri,nphi,snapskip,snap1,snap2
  real, dimension(:,:), allocatable :: coord,disp,disptot,disptot_sum,azi_xsect2,azi_xsect1
  real, dimension(:,:,:), allocatable :: azi,fluid_prefact,azi_meri,mij_snap
  real(kind=realkind), dimension(:), allocatable :: vp,x,y,z,vptop,vpbot,vpmeri,xtot,ytot,ztot,azi_phi_meri,vptot
  double precision, dimension(:), allocatable :: xtop,ytop,ztop,xbot,ybot,zbot,azi_phi_top,azi_phi_bot
  real, dimension(:), allocatable :: xmeri,ymeri,zmeri,longit_snap
  real :: dphi,phi0,r,theta_meri,smallval_meri,dr,theta0,theta_bot,theta_top,phi
  integer, dimension(:), allocatable :: ind_proc_top_tmp,ind_pts_top_tmp,ind_proc_bot_tmp,ind_pts_bot_tmp
  integer, dimension(:), allocatable :: ind_proc_top,ind_pts_top,ind_proc_bot
  integer, dimension(:), allocatable :: ind_pts_bot,azi_ind_top,azi_ind_bot,azi_ind_meri
  integer, dimension(:), allocatable :: ind_proc_meri_tmp,ind_pts_meri_tmp,ind_proc_meri,ind_pts_meri
  real :: smallval_north_top, smallval_south_top,smallval_north_bot, smallval_south_bot,r0,coord9(9,2),disp9(9,3)
  character(len=200) :: filename1
  logical :: use_meri,use_top,use_bot
  character(len=4) :: appmynum2,appiproc
  integer :: nptstot,k1,k2,npts_fluid,skipfact,j,naang,i,iii,iel
  real :: rbot,rtop,dtheta,phi_incr,smallval
  double precision :: r8,theta8,phi8

  skipfact=2
  smallval=10000.

  ! read snap plot parameters
  open(unit=99,file='param_snaps')
  read(99,*)phi0
  write(6,*)'starting azimuth/phi for cross section on the left [deg]:',phi0
  read(99,*)dphi
  write(6,*)'ending azimuth/phi for cross section on the right [deg]:',dphi
  read(99,*)rtop
  write(6,*)'top surface [km]:',rtop
  read(99,*)rbot
  write(6,*)'bottom surface [km]:',rbot
  read(99,*)theta_meri
  write(6,*)'colatitude of meridional cross section:',theta_meri
  read(99,*)snap1,snap2,snapskip
  write(6,*)'1st,last snap, skipfactor:',snap1,snap2,snapskip
  read(99,*)use_meri
  write(6,*)'consider meridional cross section?',use_meri
  read(99,*)use_top
  write(6,*)'consider top surface?',use_top
  read(99,*)use_bot
  write(6,*)'consider bottom surface?',use_bot
  close(99)

  phi0=phi0/180.*pi; dphi=dphi/180.*pi; rtop=rtop*1000.; rbot=rbot*1000.
  theta_meri=theta_meri*pi/180.

  npts=nelem*16
  nptstot=npts*nproc_mesh
  write(6,*)'number of points per proc, total points:',npts,nptstot

     ! load and construct global mesh (one semi-disk)
     write(6,*)'reading partitioned mesh...'
     allocate(coord(nptstot,2))
     smallval_north_top = rtop; smallval_south_top = rtop
     smallval_north_bot = rbot; smallval_south_bot = rbot
     write(6,*)'Smallest distance to rtop (North,South) BEFORE [km]:', &
          real(smallval_north_top/1000.),real(smallval_south_top/1000.)
     write(6,*)'Smallest distance to rbot (North,South) BEFORE [km]:', &
          real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)
     do iproc=0,nproc_mesh-1
        i=0
        do iel=1,nelem
           iii=0
           coord(iproc*npts+i+1,:) = mesh(ibeg,ibeg,iel,iproc,1:2)
           coord(iproc*npts+i+2,:) = mesh(ibeg,iend/2,iel,iproc,1:2)
           coord(iproc*npts+i+3,:) = mesh(iend/2,iend/2,iel,iproc,1:2)
           coord(iproc*npts+i+4,:) = mesh(iend/2,ibeg,iel,iproc,1:2)

           coord(iproc*npts+i+5,:) = mesh(ibeg,iend/2,iel,iproc,1:2)
           coord(iproc*npts+i+6,:) = mesh(ibeg,iend,iel,iproc,1:2)
           coord(iproc*npts+i+7,:) = mesh(iend/2,iend,iel,iproc,1:2)
           coord(iproc*npts+i+8,:) = mesh(iend/2,iend/2,iel,iproc,1:2)

           coord(iproc*npts+i+9,:) =   mesh(iend/2,ibeg,iel,iproc,1:2)
           coord(iproc*npts+i+10,:) = mesh(iend/2,iend/2,iel,iproc,1:2)
           coord(iproc*npts+i+11,:) = mesh(iend,iend/2,iel,iproc,1:2)
           coord(iproc*npts+i+12,:) = mesh(iend,ibeg,iel,iproc,1:2)

           coord(iproc*npts+i+13,:) = mesh(iend/2,iend/2,iel,iproc,1:2)
           coord(iproc*npts+i+14,:) = mesh(iend/2,iend,iel,iproc,1:2)
           coord(iproc*npts+i+15,:) = mesh(iend,iend,iel,iproc,1:2)
           coord(iproc*npts+i+16,:) = mesh(iend,iend/2,iel,iproc,1:2)

           ! determine minimal distance from rtop and rbot
            do iii=1,16
               r0 = sqrt(coord(iproc*npts+i+iii,1)**2+coord(iproc*npts+i+iii,2)**2) 
               if (coord(iproc*npts+i+iii,2)>0.0 .and. abs(r0-rtop)< smallval_north_top) then ! north
                  smallval_north_top = abs(r0-rtop)
               elseif (coord(iproc*npts+i+iii,2)<0.0 .and. abs(r0-rtop)< smallval_south_top) then ! south
                  smallval_south_top = abs(r0-rtop)
               endif
               
               if (coord(iproc*npts+i+iii,2)>0.0 .and. abs(r0-rbot)< smallval_north_bot) then ! north            
                  smallval_north_bot = abs(r0-rbot)
               elseif (coord(iproc*npts+i+iii,2)<0.0 .and. abs(r0 -rbot)< smallval_south_bot) then ! south
                  smallval_south_bot = abs(r0-rbot)
               endif
           enddo
            i=i+16
        enddo
     enddo
     smallval_north_top = smallval_north_top + smallval
     smallval_south_top = smallval_south_top + smallval
     smallval_north_bot = smallval_north_bot + smallval
     smallval_south_bot = smallval_south_bot + smallval 
     write(6,*)'Smallest distance to rtop (North,South) AFTER [km]:', &
          real(smallval_north_top/1000.),real(smallval_south_top/1000.)
     write(6,*)'Smallest distance to rbot (North,South) AFTER [km]:', &
          real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)
     if (use_top ) allocate(ind_proc_top_tmp(floor(real(nptstot)/10.)),ind_pts_top_tmp(floor(real(nptstot)/10.)))
     if (use_bot ) allocate(ind_proc_bot_tmp(floor(real(nptstot)/10.)),ind_pts_bot_tmp(floor(real(nptstot)/10.)))
     k1=0; k2=0; 
     do iproc=0,nproc_mesh-1
        do i=1,npts,skipfact
           ! check for top and bottom radii
           r0 = sqrt(coord(iproc*npts+i,1)**2+coord(iproc*npts+i,2)**2) 

           if (use_top ) then
              if ( (coord(iproc*npts+i,2)>=0.0 .and.  &
                   abs(r0-rtop)<= smallval_north_top) .or. &
                   (coord(iproc*npts+i,2)<0.0 .and.  &
                   abs(r0-rtop) <= smallval_south_top)) then 
                 k1=k1+1         
                 ind_proc_top_tmp(k1) = iproc
                 ind_pts_top_tmp(k1) = i
              endif
           endif

           if (use_bot) then 
              if ( (coord(iproc*npts+i,2)>=0.0 .and.  &
                   abs(r0-rbot) <= smallval_north_bot) .or. &
                   (coord(iproc*npts+i,2)<0.0 .and.  &
                   abs(r0-rbot) <= smallval_south_bot)) then 
                 k2=k2+1
                 ind_proc_bot_tmp(k2) = iproc
                 ind_pts_bot_tmp(k2) = i
              endif
           endif

        enddo
     enddo
     npts_top=k1;     npts_bot=k2
     write(6,*)'# points on top,bottom:',npts_top,npts_bot
     write(6,*)'allocating index arrays for surfaces....'
     if (use_top ) then
        allocate(ind_proc_top(npts_top),ind_pts_top(npts_top))
        ind_proc_top=ind_proc_top_tmp(1:npts_top); ind_pts_top=ind_pts_top_tmp(1:npts_top)
        deallocate(ind_proc_top_tmp,ind_pts_top_tmp)
     endif
     if (use_bot) then
        allocate(ind_proc_bot(npts_bot),ind_pts_bot(npts_bot))
        ind_proc_bot=ind_proc_bot_tmp(1:npts_bot); ind_pts_bot=ind_pts_bot_tmp(1:npts_bot)
        deallocate(ind_proc_bot_tmp,ind_pts_bot_tmp)
     endif

     ! xyz coordinates---------------------------------------------------------------------------------------------
     write(6,*)'defining xyz...'
     allocate(x(1:2*nptstot),y(1:2*nptstot),z(1:2*nptstot));x=0.;y=0.;z=0.

     ! left cross section---------------------------------------------------------------------------------------------
     call sphi2xy_arr(x(1:nptstot),y(1:nptstot),coord(:,1),phi0,nptstot)
     z(1:nptstot)=coord(1:nptstot,2)
     write(6,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)

     ! right cross section---------------------------------------------------------------------------------------------
     call sphi2xy_arr(x(nptstot+1:2*nptstot),y(nptstot+1:2*nptstot),coord(:,1),phi0+dphi,nptstot)
     z(nptstot+1:2*nptstot)=coord(1:nptstot,2)
     write(6,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)
     allocate(vp(2*nptstot))
     do i=1,2*nptstot  
        r= sqrt(x(i)**2+y(i)**2+z(i)**2)
        vp(i)=prem(r,'v_p')
     enddo
     vp=vp/maxval(vp)*0.002-0.001 ! rescale vp
     filename1='Data/mesh_xsect'
     call write_VTK_bin_scal_mesh3d(x,y,z,vp,2*nptstot,0,filename1)
     filename1='Data/mesh_xsect_cell'
     call write_VTK_bin_scal_topology(x,y,z,vp,2*nptstot/4,filename1)

     ! top surface---------------------------------------------------------------------------------------------
     if (use_top) then
        k1=0
        write(6,*)'defining top surface...'
        dtheta=rtop*pi/npts_top
        naang = floor(real(npts_top)/real(2.))**2*6*4
        write(6,*)'points on top surface:',naang
        allocate(xtop(1:naang),ytop(1:naang),ztop(1:naang))
        allocate(azi_phi_top(1:naang),azi_ind_top(1:naang))
        call construct_surface_cubed_sphere(npts_top,npts,dble(rtop),ind_proc_top,ind_pts_top,nptstot,dble(coord),k1, &
                                  dble(dphi),dble(phi0),'outside',naang,xtop,ytop,ztop,azi_phi_top,azi_ind_top)
        write(6,*)'allocating vptop...',k1
        allocate(vptop(k1))
        vptop(1:k1) = vp(ind_proc_top(1)*npts+ind_pts_top(1))
        write(6,*)'save into cell vtk...',size(xtop),k1
        filename1='Data/mesh_top_cell'
        call write_VTK_bin_scal_topology(real(xtop(1:k1)),real(ytop(1:k1)),real(ztop(1:k1)),vptop(1:k1),k1/4,filename1)
     endif

     ! bottom surface ---------------------------------------------------------------------------------------------
     if (use_bot) then
        k2=0
        write(6,*)'defining bot surface...'
        dtheta=rbot*pi/npts_bot
        naang = floor(real(npts_bot)/real(2.))**2*6*4
        write(6,*)'points on bottom surface:',naang
        allocate(xbot(1:naang),ybot(1:naang),zbot(1:naang))
        allocate(azi_phi_bot(1:naang),azi_ind_bot(1:naang))
        call construct_surface_cubed_sphere(npts_bot,npts,dble(rbot),ind_proc_bot,ind_pts_bot,nptstot,dble(coord),k2, &
                                  dble(dphi),dble(phi0),'innside',naang,xbot,ybot,zbot,azi_phi_bot,azi_ind_bot)
        allocate(vpbot(k2))
        vpbot(1:k2) = vp(ind_proc_bot(1)*npts+ind_pts_bot(1))
        filename1='Data/mesh_bot_cell'
        call write_VTK_bin_scal_topology(real(xbot(1:k2)),real(ybot(1:k2)),real(zbot(1:k2)),vpbot,k2/4,filename1)
     endif
     deallocate(coord)

!=========================================================================
     write(6,*)'assembling different surfaces to one coordinate array...'
     npts_cell=2*nptstot+k1+k2
     write(6,*)'allocating global velocity array...'
     allocate(vptot(npts_cell))
     vptot(1:2*nptstot)=vp
     if (use_top)     vptot(2*nptstot+1:2*nptstot+k1)=vptop
     if (use_bot)     vptot(2*nptstot+k1+1:2*nptstot+k1+k2)=vpbot
     deallocate(vp)

     write(6,*)'allocating global coordinates xgd...'
     allocate(xgd_cell(1:npts_cell),ygd_cell(1:npts_cell),zgd_cell(1:npts_cell))
     write(6,*)'Filling xgd_cell...'
     xgd_cell(1:2*nptstot)=x(1:2*nptstot); ygd_cell(1:2*nptstot)=y(1:2*nptstot); zgd_cell(1:2*nptstot)=z(1:2*nptstot)
     write(6,*)'  ... done with cross sections ...'
     if (use_top) then 
        xgd_cell(2*nptstot+1:2*nptstot+k1)=real(xtop(1:k1))
        ygd_cell(2*nptstot+1:2*nptstot+k1)=real(ytop(1:k1))
        zgd_cell(2*nptstot+1:2*nptstot+k1)=real(ztop(1:k1))
     endif
     if (use_bot) then 
        xgd_cell(2*nptstot+k1+1:2*nptstot+k1+k2)=real(xbot(1:k2))
        ygd_cell(2*nptstot+k1+1:2*nptstot+k1+k2)=real(ybot(1:k2))
        zgd_cell(2*nptstot+k1+1:2*nptstot+k1+k2)=real(zbot(1:k2))
     endif

     cell_topology=.true.

     filename1='Data/mesh_all_cell'
     call write_VTK_bin_scal_topology(real(xgd_cell),real(ygd_cell),real(zgd_cell),vptot,npts_cell/4,filename1)

! reducing overlapping coordinates: double loop
     allocate(ind_tot2cell(1:npts_cell,10),ind_cell2tot(1:npts_cell),valence_tot2cell(1:npts_cell))
     call reduce_coord_overlap(npts_cell,xgd_cell,ygd_cell,zgd_cell,npts,& ! THIS NOW INCLUDES DOMAIN DECOMPOSITION!!
                               ind_tot2cell,ind_cell2tot,valence_tot2cell,nloc)
      call define_io_appendix(appiproc,mynum)

     write(6,*)'allocating xgd, size:',npts
     allocate(xgd(1:npts),ygd(1:npts),zgd(1:npts))
     do i=1,nloc
        xgd(ind_cell2tot(i)) = xgd_cell(i)
        ygd(ind_cell2tot(i)) = ygd_cell(i)
        zgd(ind_cell2tot(i)) = zgd_cell(i)
     enddo
     deallocate(xgd_cell,ygd_cell,zgd_cell)
     allocate(xgd_cell(1:npts_cell),ygd_cell(1:npts_cell),zgd_cell(1:npts_cell))

     filename1='Data/mesh_all_totpts_'//appiproc
     call write_VTK_bin_scal_mesh3d(xgd,ygd,zgd,zgd,npts,0,filename1)

! test the mapping cell to total and back ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xgd_cell=0; ygd_cell=0;zgd_cell=0.
write(6,*)'max ind:',maxval(ind_tot2cell),maxval(ind_cell2tot)
do i=1,npts
   do j=1,valence_tot2cell(i)
      xgd_cell(ind_tot2cell(i,j)) = xgd(i)
      ygd_cell(ind_tot2cell(i,j)) = ygd(i)
      zgd_cell(ind_tot2cell(i,j)) = zgd(i)
   enddo
enddo
filename1='Data/mesh_all_cell_retry_'//appiproc
call write_VTK_bin_scal_topology(real(xgd_cell(1:nloc)),real(ygd_cell(1:nloc)),real(zgd_cell(1:nloc)),&
                                 vptot(1:nloc),nloc/4,filename1)
! end test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write(6,*)mynum,'Done with 3D slices in 4-plus-center-corner cell topology.'
allocate(kern_cell(nloc))

end subroutine compute_3dslices_cell16
!---------------------------------------------------------------------------------------



!=======================================================================================
subroutine compute_3dslices_cell16_old(mesh,phr)
  
  use data_rota
  use coord_trafo, only : sphi2xy_arr,sphi2xy,rthetaphi2xyz
  use input_output, only: write_VTK_bin_scal_mesh3d,write_VTK_bin_scal_topology
  use mpi

  real(kind=realkind),intent(in) :: phr
  real(kind=realkind), intent(in) ::mesh(0:iend,1:iend,1:nelem,0:nproc_mesh-1,1:2)
  integer iproc,npts,npts_top,npts_bot,npts_meri,nphi,snapskip,snap1,snap2
  real, dimension(:,:), allocatable :: coord,disp,disptot,disptot_sum,azi_xsect2,azi_xsect1
  real, dimension(:,:,:), allocatable :: azi,fluid_prefact,azi_meri,mij_snap
  real(kind=realkind), dimension(:), allocatable :: vp,x,y,z,vptop,vpbot,vpmeri,azi_phi_meri
  double precision, dimension(:), allocatable :: xtop,ytop,ztop,xbot,ybot,zbot,azi_phi_top,azi_phi_bot
  real, dimension(:), allocatable :: xmeri,ymeri,zmeri,longit_snap
  real :: dphi,phi0,r,theta_meri,smallval_meri,dr,theta0,theta_bot,theta_top,phi
  integer, dimension(:), allocatable :: ind_proc_top_tmp,ind_pts_top_tmp,ind_proc_bot_tmp,ind_pts_bot_tmp
  integer, dimension(:), allocatable :: ind_proc_top,ind_pts_top,ind_proc_bot
  integer, dimension(:), allocatable :: ind_pts_bot,azi_ind_top,azi_ind_bot,azi_ind_meri
  integer, dimension(:), allocatable :: ind_proc_meri_tmp,ind_pts_meri_tmp,ind_proc_meri,ind_pts_meri
  real :: smallval_north_top, smallval_south_top,smallval_north_bot, smallval_south_bot,r0,coord9(9,2),disp9(9,3)
  character(len=200) :: filename1
  logical :: use_meri,use_top,use_bot
  character(len=4) :: appmynum2
  integer :: nptstot,k1,k2,npts_fluid,skipfact,j,naang,npts_fluid_read,i,iii,iel,ipol,jpol
  real :: rbot,rtop,dtheta,phi_incr,smallval
  double precision :: r8,theta8,phi8

  skipfact=1
  smallval=10000.

  ! read snap plot parameters
  open(unit=99,file='param_snaps')
  read(99,*)phi0
  write(6,*)'starting azimuth/phi for cross section on the right [deg]:',phi0
  read(99,*)dphi
  write(6,*)'ending azimuth/phi for cross section on the left [deg]:',dphi
  read(99,*)rtop
  write(6,*)'top surface [km]:',rtop
  read(99,*)rbot
  write(6,*)'bottom surface [km]:',rbot
  read(99,*)theta_meri
  write(6,*)'colatitude of meridional cross section:',theta_meri
  read(99,*)snap1,snap2,snapskip
  write(6,*)'1st,last snap, skipfactor:',snap1,snap2,snapskip
  read(99,*)use_meri
  write(6,*)'consider meridional cross section?',use_meri
  read(99,*)use_top
  write(6,*)'consider top surface?',use_top
  read(99,*)use_bot
  write(6,*)'consider bottom surface?',use_bot
  close(99)

  phi0=phi0/180.*pi; dphi=dphi/180.*pi; rtop=rtop*1000.; rbot=rbot*1000.
  theta_meri=theta_meri*pi/180.

  npts=nelem*16
  nptstot=npts*nproc_mesh
  write(6,*)'number of points per proc, total points:',npts,nptstot

     ! load and construct global mesh (one semi-disk)
     write(6,*)'reading partitioned mesh...'
     allocate(coord(nptstot,2))
     smallval_north_top = rtop; smallval_south_top = rtop
     smallval_north_bot = rbot; smallval_south_bot = rbot
     write(6,*)'Smallest distance to rtop (North,South) BEFORE [km]:', &
          real(smallval_north_top/1000.),real(smallval_south_top/1000.)
     write(6,*)'Smallest distance to rbot (North,South) BEFORE [km]:', &
          real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)
     do iproc=0,nproc_mesh-1
        i=0
        do iel=1,nelem
           iii=0
           do ipol=ibeg,iend
              do jpol=ibeg,iend
                 iii=iii+1
                 coord9(iii,1:2) = mesh(ipol,jpol,iel,iproc,1:2)
              enddo
           enddo
           coord(iproc*npts+i+1:iproc*npts+i+2,:) = coord9(1:2,:) ! (0,0), (0,npol/2)
           coord(iproc*npts+i+3,:) = coord9(5,:)    ! (npol/2,npol/2)
           coord(iproc*npts+i+4,:) = coord9(4,:)    ! (npol/2,0)
           coord(iproc*npts+i+5:iproc*npts+i+6,:) = coord9(2:3,:)    ! (0,npol/2),(0,npol)
           coord(iproc*npts+i+7,:) = coord9(6,:)    ! (npol/2,npol)
           coord(iproc*npts+i+8,:) = coord9(5,:)    ! (npol/2,npol/2)
           coord(iproc*npts+i+9:iproc*npts+i+10,:) = coord9(4:5,:)    ! (npol/2,0),(npol/2,npol/2)
           coord(iproc*npts+i+11,:) = coord9(8,:)    ! (npol,npol/2)
           coord(iproc*npts+i+12,:) = coord9(7,:)    ! (npol,0)            
           coord(iproc*npts+i+13:iproc*npts+i+14,:) = coord9(5:6,:)    ! (npol/2,npol/2),(npol/2,npol)
           coord(iproc*npts+i+15,:) = coord9(9,:)    ! (npol,npol)
           coord(iproc*npts+i+16,:) = coord9(8,:)    ! (npol,npol/2)            

           ! determine minimal distance from rtop and rbot
            do iii=1,16
               r0 = sqrt(coord(iproc*npts+i+iii,1)**2+coord(iproc*npts+i+iii,2)**2) 
               if (coord(iproc*npts+i+iii,2)>0.0 .and. abs(r0-rtop)< smallval_north_top) then ! north
                  smallval_north_top = abs(r0-rtop)
               elseif (coord(iproc*npts+i+iii,2)<0.0 .and. abs(r0-rtop)< smallval_south_top) then ! south
                  smallval_south_top = abs(r0-rtop)
               endif
               
               if (coord(iproc*npts+i+iii,2)>0.0 .and. abs(r0-rbot)< smallval_north_bot) then ! north            
                  smallval_north_bot = abs(r0-rbot)
               elseif (coord(iproc*npts+i+iii,2)<0.0 .and. abs(r0 -rbot)< smallval_south_bot) then ! south
                  smallval_south_bot = abs(r0-rbot)
               endif
           enddo
            i=i+16
        enddo
     enddo
     smallval_north_top = smallval_north_top + smallval
     smallval_south_top = smallval_south_top + smallval
     smallval_north_bot = smallval_north_bot + smallval
     smallval_south_bot = smallval_south_bot + smallval 
     write(6,*)'Smallest distance to rtop (North,South) AFTER [km]:', &
          real(smallval_north_top/1000.),real(smallval_south_top/1000.)
     write(6,*)'Smallest distance to rbot (North,South) AFTER [km]:', &
          real(smallval_north_bot/1000.),real(smallval_south_bot/1000.)
     if (use_top ) allocate(ind_proc_top_tmp(floor(real(nptstot)/10.)),ind_pts_top_tmp(floor(real(nptstot)/10.)))
     if (use_bot ) allocate(ind_proc_bot_tmp(floor(real(nptstot)/10.)),ind_pts_bot_tmp(floor(real(nptstot)/10.)))
     k1=0; k2=0; 
     do iproc=0,nproc_mesh-1
        do i=1,npts,skipfact
           ! check for top and bottom radii
           r0 = sqrt(coord(iproc*npts+i,1)**2+coord(iproc*npts+i,2)**2) 

           if (use_top ) then
              if ( (coord(iproc*npts+i,2)>=0.0 .and.  &
                   abs(r0-rtop)<= smallval_north_top) .or. &
                   (coord(iproc*npts+i,2)<0.0 .and.  &
                   abs(r0-rtop) <= smallval_south_top)) then 
                 k1=k1+1         
                 ind_proc_top_tmp(k1) = iproc
                 ind_pts_top_tmp(k1) = i
              endif
           endif

           if (use_bot) then 
              if ( (coord(iproc*npts+i,2)>=0.0 .and.  &
                   abs(r0-rbot) <= smallval_north_bot) .or. &
                   (coord(iproc*npts+i,2)<0.0 .and.  &
                   abs(r0-rbot) <= smallval_south_bot)) then 
                 k2=k2+1
                 ind_proc_bot_tmp(k2) = iproc
                 ind_pts_bot_tmp(k2) = i
              endif
           endif

        enddo
     enddo
     npts_top=k1;     npts_bot=k2
     write(6,*)'# points on top,bottom:',npts_top,npts_bot
     write(6,*)'allocating index arrays for surfaces....'
     if (use_top ) then
        allocate(ind_proc_top(npts_top),ind_pts_top(npts_top))
        ind_proc_top=ind_proc_top_tmp(1:npts_top); ind_pts_top=ind_pts_top_tmp(1:npts_top)
        deallocate(ind_proc_top_tmp,ind_pts_top_tmp)
     endif
     if (use_bot) then
        allocate(ind_proc_bot(npts_bot),ind_pts_bot(npts_bot))
        ind_proc_bot=ind_proc_bot_tmp(1:npts_bot); ind_pts_bot=ind_pts_bot_tmp(1:npts_bot)
        deallocate(ind_proc_bot_tmp,ind_pts_bot_tmp)
     endif

     ! xyz coordinates---------------------------------------------------------------------------------------------
     write(6,*)'defining xyz...'
     allocate(x(1:2*nptstot),y(1:2*nptstot),z(1:2*nptstot));x=0.;y=0.;z=0.

     ! left cross section---------------------------------------------------------------------------------------------
     call sphi2xy_arr(x(1:nptstot),y(1:nptstot),coord(:,1),phi0,nptstot)
     z(1:nptstot)=coord(1:nptstot,2)
     write(6,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)

     ! right cross section---------------------------------------------------------------------------------------------
     call sphi2xy_arr(x(nptstot+1:2*nptstot),y(nptstot+1:2*nptstot),coord(:,1),phi0+dphi,nptstot)
     z(nptstot+1:2*nptstot)=coord(1:nptstot,2)
     write(6,*)'max s,x,y:',maxval(coord(:,1)),maxval(x),maxval(y)
     ! save xyz
     allocate(vp(2*nptstot))
     do i=1,2*nptstot  
        r= sqrt(x(i)**2+y(i)**2+z(i)**2)
        vp(i)=prem(r,'v_p')
     enddo

     ! rescale vp
     vp=vp/maxval(vp)*0.002-0.001
     filename1='Data/mesh_xsect'
     call write_VTK_bin_scal_mesh3d(x,y,z,vp,2*nptstot,0,filename1)
     filename1='Data/mesh_xsect_cell'
     call write_VTK_bin_scal_topology(x,y,z,vp,2*nptstot/4,filename1)

     ! top surface---------------------------------------------------------------------------------------------
     if (use_top) then
        k1=0
        write(6,*)'defining top surface...'
        dtheta=rtop*pi/npts_top
        naang = floor(real(npts_top)/real(2.))**2*6*4
        write(6,*)'points on top surface:',naang
        allocate(xtop(1:naang),ytop(1:naang),ztop(1:naang))
        allocate(azi_phi_top(1:naang),azi_ind_top(1:naang))
        call construct_surface_cubed_sphere(npts_top,npts,dble(rtop),ind_proc_top,ind_pts_top,nptstot,dble(coord),k1, &
                                  dble(dphi),dble(phi0),'outside',naang,xtop,ytop,ztop,azi_phi_top,azi_ind_top)
        ! extract vp
        write(6,*)'allocating vptop...',k1
        allocate(vptop(k1))
        vptop(1:k1) = vp(ind_proc_top(1)*npts+ind_pts_top(1))
        write(6,*)'save into cell vtk...',size(xtop),k1
        filename1='Data/mesh_top_cell'
        call write_VTK_bin_scal_topology(real(xtop(1:k1)),real(ytop(1:k1)),real(ztop(1:k1)),vptop(1:k1),k1/4,filename1)
        open(unit=99,file='Data/xyz_top.dat')
        do i=1,k1
           write(99,*)xtop(i),ytop(i),ztop(i)
        enddo
        close(99)
     endif

     ! bottom surface ---------------------------------------------------------------------------------------------
     if (use_bot) then
        k2=0
        write(6,*)'defining bot surface...'
        dtheta=rbot*pi/npts_bot
        naang = floor(real(npts_bot)/real(2.))**2*6*4
        write(6,*)'points on bottom surface:',naang
        allocate(xbot(1:naang),ybot(1:naang),zbot(1:naang))
        allocate(azi_phi_bot(1:naang),azi_ind_bot(1:naang))
        call construct_surface_cubed_sphere(npts_bot,npts,dble(rbot),ind_proc_bot,ind_pts_bot,nptstot,dble(coord),k2, &
                                  dble(dphi),dble(phi0),'innside',naang,xbot,ybot,zbot,azi_phi_bot,azi_ind_bot)

        ! extract vp
        allocate(vpbot(k2))
        vpbot(1:k2) = vp(ind_proc_bot(1)*npts+ind_pts_bot(1))

        ! save into cell vtk
        filename1='Data/mesh_bot_cell'
        call write_VTK_bin_scal_topology(real(xbot(1:k2)),real(ybot(1:k2)),real(zbot(1:k2)),vpbot,k2/4,filename1)
     endif
     deallocate(coord)
     deallocate(vp)

     ! assembling everything to one coordinate array-----------------------------------------------------------------
     allocate(xgd_tot(2*nptstot+k1+k2),ygd_tot(2*nptstot+k1+k2),zgd_tot(2*nptstot+k1+k2))

     xgd_tot(1:2*nptstot)=x(1:2*nptstot); ygd_tot(1:2*nptstot)=y(1:2*nptstot); zgd_tot(1:2*nptstot)=z(1:2*nptstot)
     if (use_top) then 
        xgd_tot(2*nptstot+1:2*nptstot+k1)=real(xtop(1:k1))
        ygd_tot(2*nptstot+1:2*nptstot+k1)=real(ytop(1:k1))
        zgd_tot(2*nptstot+1:2*nptstot+k1)=real(ztop(1:k1))
     endif
     if (use_bot) then 
        xgd_tot(2*nptstot+k1+1:2*nptstot+k1+k2)=real(xbot(1:k2))
        ygd_tot(2*nptstot+k1+1:2*nptstot+k1+k2)=real(ybot(1:k2))
        zgd_tot(2*nptstot+k1+1:2*nptstot+k1+k2)=real(zbot(1:k2))
     endif

     ! save mesh for kerner
     open(unit=99,file='Data/mesh_tot.xyz')
     do i=1,2*nptstot+k1+k2
        write(99,*)xgd_tot(i),ygd_tot(i),zgd_tot(i)
     enddo
     close(99)

     cell_topology=.true.

end subroutine compute_3dslices_cell16_old
!---------------------------------------------------------------------------------------

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine construct_surface_cubed_sphere(npts_surf,npts,rsurf,ind_proc_surf,ind_pts_surf,&
                            nptstot,coord1,kpts,dphi,phi0,in_or_out,n,xsurf,ysurf,zsurf,azi_phi_surf,azi_ind_surf)

use coord_trafo, only : xyz2rthetaphi_dble,get_r_theta_dble
use splib, only : get_welegl,zelegl

!!!!! BEG CUBED SPHERE
! af2tnm: along a great circle, there's the equivalent of two chunks of the cubed
! sphere. According to the convention we use in defining the mapping in the cubed sphere
! module, that amounts to 2*nang spectral elements.
! Assuming we will be using npol_cs=1 in the following, nang has to be even.
! We therefore want nang=.5*(npts_surf-1) if npts_surf is odd
! We therefore want nang=.5*(npts_surf) if npts_surf is even
!
!use data_all
implicit none

integer, intent(in) :: npts_surf,npts,nptstot,n
double precision, intent(in) :: rsurf,dphi,phi0
character(len=7), intent(in) :: in_or_out
integer, dimension(1:npts_surf), intent(in) :: ind_proc_surf,ind_pts_surf
double precision, dimension(nptstot,2), intent(in) :: coord1
double precision, dimension(1:n), intent(out) :: xsurf,ysurf,zsurf,azi_phi_surf
integer, dimension(1:n), intent(out) :: azi_ind_surf
integer, intent(out) :: kpts

double precision, allocatable,dimension(:) :: xi_el,eta_el,r_el
integer, allocatable, dimension(:,:,:,:) :: number_el
double precision, allocatable,dimension(:) :: wt_i,xi_i,xi_j,wt_j,wt_k,dxi,dxk,xi_k
double precision, allocatable,dimension(:,:,:,:) :: xcol,ycol,zcol
double precision, allocatable,dimension(:,:,:,:) :: x_el,y_el,z_el
double precision :: dang,C,D,re,ri,teta,tksi,tr,Xe,Ye,delta
integer :: npol_cs,ii,iii,izone,nang,nelt,nr,nrpol,iel,nel_surf,jj,ipol,jpol,kpol,j,i
double precision ::  dist,r_ref,th_ref,th,dr,r,xc,yc,zc,ph
double precision, parameter :: pi = 3.1415926535898

    write(6,*)'computing cubed sphere for surface at r=',rsurf
    write(6,*)'pts surf,total:',npts_surf,npts
    nang=floor(sqrt(real(n))/6./2.)
    write(6,*)'naang,nang,npts_surf:',n,nang,npts_surf
    write(6,*)'phi,dphi:',phi0*180./pi,dphi*180./pi
    write(6,*)'in or outside:',in_or_out

     nr=1 ! one radial layer only
     write(6,*)'NANG:',nang
     allocate(xi_el(0:nang),eta_el(0:nang),r_el(0:nr))
     xi_el = 0.d0;  eta_el = 0.d0 ; r_el = 0.d0
     re=rsurf
     ri=rsurf-1.d0 ! (in km)
     !
     dang=pi/(2.d0*dble(nang))
     dr=(re-ri)/dble(nr)
     do i=0,nang
      xi_el(i)=-pi*0.25d0+dang*dble(i)
      eta_el(i)=-pi*0.25d0+dang*dble(i)
     enddo
     do i=0,nr
      r_el(i)=ri+dr*dble(i)
     enddo
!     write(6,*)'allocating element coords...'
     allocate(x_el(0:nang,0:nang,0:nr,1:6))
     allocate(y_el(0:nang,0:nang,0:nr,1:6))
     allocate(z_el(0:nang,0:nang,0:nr,1:6))
     x_el = 0.d0; y_el=0.d0; z_el = 0.d0
     allocate(number_el(0:nang,0:nang,0:nr,1:6))
     number_el = 0
!     The mapping is defined following Emmanuel Chaljub (PhD thesis, 2000).
     do izone = 1, 6
     do iii=0,nr-1     ! loop over  r
      do ii=0,nang-1   ! loop over eta
       do i=0,nang-1   ! loop over ksi
        number_el(i,ii,iii,izone) = (izone-1)*(nang*nang*nr)+&
                              iii*(nang*nang)+((ii*nang)+i+1)
       end do
      end do
     end do
     end do
     nelt = maxval(number_el)
     nrpol=1
     npol_cs=1
     allocate(xi_i(0:npol_cs),xi_j(0:npol_cs),xi_k(0:nrpol))
     allocate(wt_i(0:npol_cs),wt_j(0:npol_cs),wt_k(0:nrpol))
!
!     dummies needed by spectral library
     allocate(dxi(0:npol_cs),dxk(0:nrpol))
!
     call ZELEGL(npol_cs,xi_i,dxi)
     call get_welegl(npol_cs,xi_i,wt_i)
!     xi_i -> npol_cs+1 pts de quadrature pour l'ordre ploynomial npol_cs
     call ZELEGL(npol_cs,xi_j,dxi)
     call get_welegl(npol_cs,xi_j,wt_j)
     call ZELEGL(nrpol,xi_k,dxk)
     call get_welegl(nrpol,xi_k,wt_k)
!     get rid of dummies
     deallocate(dxk,dxi)
     write(6,*)'Now define collocation points for the SEM-cubed sphere grid'
      allocate(xcol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
      allocate(ycol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
      allocate(zcol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
      kpts=0
      nel_surf=0

      write(6,*)'ZONE NANG:',nr,nang,npol_cs,nelt
      write(6,*)'ZONE NANG 2 nel_surf:',6*nr*nang**2
      
      do izone=1,6         ! loop over the chunks
!         write(6,*)'working on chunk',izone
       do iii=0,nr-1       ! loop over r   (the global r)
        do ii=0,nang-1     ! loop over eta (the global eta)
         do i=0,nang-1     ! loop over xi  (the global xi)
          iel = number_el(i,ii,iii,izone)
          nel_surf=nel_surf+1
          do kpol = 0, 0 ! TNM nrpol  ! loop over the elemental k index (r direction)
           do jpol = 0, npol_cs  ! loop over the elemental j index (eta direction)
            do ipol = 0, npol_cs ! loop over the elemental i index (xi direction)
             tksi= xi_el(  i) +(1.d0+ xi_i(ipol))*.5d0*dang
             teta=eta_el( ii) +(1.d0 + xi_j(jpol))*.5d0*dang
             tr=    r_el(iii) +(1.d0 + xi_k(kpol))*.5d0*dr
             Xe=tan(tksi)
             Ye=tan(teta)
             C=(1.d0+Xe**2)**(0.5d0)
             D=(1.d0+Ye**2)**(0.5d0)
             delta=1.d0+Xe**2+Ye**2
             kpts=kpts+1
             if(izone==1) then
              Xcol(ipol,jpol,kpol,iel)=tr*delta**(-0.5d0)
              Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
              Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
             endif
             if(izone==2) then
              Xcol(ipol,jpol,kpol,iel)=-tr*Xe*delta**(-0.5d0)
              Ycol(ipol,jpol,kpol,iel)=tr*delta**(-0.5d0)
              Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
             endif
             if(izone==3) then
              Xcol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5d0)
              Ycol(ipol,jpol,kpol,iel)=-tr*Xe*delta**(-0.5d0)
              Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
             endif
             if(izone==4) then
              Xcol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
              Ycol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5d0)
              Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
             endif
             if(izone==5) then
              Xcol(ipol,jpol,kpol,iel)=-tr*Ye*delta**(-0.5d0)
              Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
              Zcol(ipol,jpol,kpol,iel)=tr*delta**(-0.5d0)
             endif
             if(izone==6) then
              Xcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5d0)
              Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5d0)
              Zcol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5d0)
             endif
            end do
           end do
          end do
         end do
        end do
       end do
      end do

!     At this stage we know Xcol, Ycol, Zcol for the cubed sphere
!     These arrays are four dimensional (ipol,jpol,kpol,iel)
!     Their knowledge suffice to define the vtk output that will properly take
!     into account the cubed sphere topology (see the paraview.f90 module)
!!!!! END CUBED SPHERE

      write(6,*)'number of surface pts,elems,tot pts:',npts_surf,nel_surf,kpts
      write(6,*)'max ind_proc, ind_pts:',maxval(ind_proc_surf),maxval(ind_pts_surf)
!      allocate(xsurf(1:kpts),ysurf(1:kpts),zsurf(1:kpts))
      write(6,*)size(xsurf)
      xsurf=0.d0; ysurf=0.;zsurf=0.d0
!      write(6,*)'allocating azi ind phi...'
!      allocate(azi_ind_surf(1:kpts),azi_phi_surf(1:kpts))
      iii=0
      write(6,*)'constructing 1d array for surface coordinates...'
      do iel=1,nel_surf
         if ( mod(iel,floor(nel_surf/10.))==0  )write(6,*)'percentage done:',ceiling(real(iel)/real(nel_surf)*100.)
         xc=sum(xcol(:,:,0,iel))/4.d0
         yc=sum(ycol(:,:,0,iel))/4.d0
         zc=sum(zcol(:,:,0,iel))/4.d0
         
         call xyz2rthetaphi_dble(r,th,ph,xc,yc,zc)

         if ( (in_or_out=='innside' .and. (ph>=phi0 .and. ph<=phi0+dphi) ) .or. &
                  (in_or_out=='outside' .and. (ph<=phi0 .or.  ph>=phi0+dphi)   .and. &
               ( (abs(ph-phi0)<=pi/6. .or. ph>=2.*pi-pi/6.-phi0) .or. abs(ph-(phi0+dphi))<=pi/10. ) ) ) then

            xsurf(iii+1)=xcol(0,0,0,iel)
            ysurf(iii+1)=ycol(0,0,0,iel)
            zsurf(iii+1)=zcol(0,0,0,iel)
            
            xsurf(iii+2)=xcol(npol_cs,0,0,iel)
            ysurf(iii+2)=ycol(npol_cs,0,0,iel)
            zsurf(iii+2)=zcol(npol_cs,0,0,iel)
            
            xsurf(iii+3)=xcol(npol_cs,npol_cs,0,iel)
            ysurf(iii+3)=ycol(npol_cs,npol_cs,0,iel)
            zsurf(iii+3)=zcol(npol_cs,npol_cs,0,iel)
            
            xsurf(iii+4)=xcol(0,npol_cs,0,iel)
            ysurf(iii+4)=ycol(0,npol_cs,0,iel)
            zsurf(iii+4)=zcol(0,npol_cs,0,iel)
            
            ! determine the corresponding point in the D-shape domain
            do jj=1,4
               call xyz2rthetaphi_dble(r,th, azi_phi_surf(iii+jj),xsurf(iii+jj),ysurf(iii+jj),zsurf(iii+jj))
               dist=2.d0*pi
               do i=1,npts_surf
                  call get_r_theta_dble(coord1(ind_proc_surf(i)*npts+ind_pts_surf(i),1), & 
                       coord1(ind_proc_surf(i)*npts+ind_pts_surf(i),2),r_ref,th_ref)
                  if (abs(th-th_ref)< dist) then 
                     dist=abs(th-th_ref)
                     azi_ind_surf(iii+jj)=ind_proc_surf(i)*npts+ind_pts_surf(i)
                  endif
               enddo
            enddo
            iii=iii+4
            endif ! in_or_out
         enddo

         kpts=iii
         write(6,*)'total points in surface:',kpts
         
      write(6,*)'done with cubed sphere for r=',rsurf
      deallocate(xcol,ycol,zcol,x_el,y_el,z_el)

    end subroutine construct_surface_cubed_sphere
!---------------------------------------------------------------------------------------



!!$!---------------------------------------------------------------------------------------
!!$subroutine compute_3dslices_snaps
!!$
!!$     ! load snaps ===============================================================
!!$     write(6,*)'loading snaps...'
!!$     allocate(disp(2*nptstot+k1+k2+k3,3))
!!$!     allocate(disptot(2*nptstot+k1+k2+k3,3))
!!$
!!$     do j=snap1,snap2,snapskip
!!$
!!$        disp=0.
!!$!        disptot=0.
!!$        if (any_sum_seis_true) disptot_sum=0.
!!$
!!$        do isim=1,nsim
!!$           do iproc=0,nproc_mesh-1
!!$              call define_io_appendix(appmynum,iproc)
!!$              call define_io_appendix(appmynum2,j)      
!!$              open(unit=99,file=trim(simdir(isim))//'/Data/snap_'//appmynum//'_'//appmynum2//'.dat', &
!!$                        FORM="UNFORMATTED",STATUS="OLD",POSITION="REWIND")
!!$              i=0
!!$              do ii=1,npts_read/9,skipfact
!!$                 do iii=1,9
!!$                    read(99)disp9(iii,1),disp9(iii,2),disp9(iii,3)
!!$                 enddo
!!$                 disp(iproc*npts+i+1:iproc*npts+i+2,:) = disp9(1:2,:) ! (0,0), (0,npol/2)
!!$                 disp(iproc*npts+i+3,:) = disp9(5,:)    ! (npol/2,npol/2)
!!$                 disp(iproc*npts+i+4,:) = disp9(4,:)    ! (npol/2,0)
!!$                 disp(iproc*npts+i+5:iproc*npts+i+6,:) = disp9(2:3,:)    ! (0,npol/2),(0,npol)
!!$                 disp(iproc*npts+i+7,:) = disp9(6,:)    ! (npol/2,npol)
!!$                 disp(iproc*npts+i+8,:) = disp9(5,:)    ! (npol/2,npol/2)
!!$                 disp(iproc*npts+i+9:iproc*npts+i+10,:) = disp9(4:5,:)    ! (npol/2,0),(npol/2,npol/2)
!!$                 disp(iproc*npts+i+11,:) = disp9(8,:)    ! (npol,npol/2)
!!$                 disp(iproc*npts+i+12,:) = disp9(7,:)    ! (npol,0)            
!!$                 disp(iproc*npts+i+13:iproc*npts+i+14,:) = disp9(5:6,:)    ! (npol/2,npol/2),(npol/2,npol)
!!$                 disp(iproc*npts+i+15,:) = disp9(9,:)    ! (npol,npol)
!!$                 disp(iproc*npts+i+16,:) = disp9(8,:)    ! (npol,npol/2)            
!!$                 i=i+16
!!$              enddo
!!$              close(99)
!!$
!!$              disp(iproc*npts+1:iproc*npts+npts_fluid,1)= &
!!$                   disp(iproc*npts+1:iproc*npts+npts_fluid,1) * fluid_prefact(:,1,isim)
!!$              disp(iproc*npts+1:iproc*npts+npts_fluid,2)= &
!!$                   disp(iproc*npts+1:iproc*npts+npts_fluid,2) * fluid_prefact(:,1,isim) * fluid_prefact(:,2,isim)
!!$              disp(iproc*npts+1:iproc*npts+npts_fluid,3)= &
!!$                   disp(iproc*npts+1:iproc*npts+npts_fluid,3) * fluid_prefact(:,1,isim)
!!$           enddo
!!$
!!$           disp(nptstot+1:2*nptstot,1)=disp(1:nptstot,1)
!!$           disp(nptstot+1:2*nptstot,2)=disp(1:nptstot,2)
!!$           disp(nptstot+1:2*nptstot,3)=disp(1:nptstot,3)
!!$
!!$           do i=1,k1
!!$              disp(2*nptstot+i,1) = disp(azi_ind_top(i),1)
!!$              disp(2*nptstot+i,2) = disp(azi_ind_top(i),2)
!!$              disp(2*nptstot+i,3) = disp(azi_ind_top(i),3)
!!$           enddo
!!$           do i=1,k2
!!$              disp(2*nptstot+k1+i,1) = disp(azi_ind_bot(i),1)
!!$              disp(2*nptstot+k1+i,2) = disp(azi_ind_bot(i),2)
!!$              disp(2*nptstot+k1+i,3) = disp(azi_ind_bot(i),3)
!!$           enddo
!!$           if (use_meri) then 
!!$              do i=1,k3
!!$                 disp(2*nptstot+k1+k2+i,1) = disp(azi_ind_meri(i),1)   
!!$                 disp(2*nptstot+k1+k2+i,2) = disp(azi_ind_meri(i),2)
!!$                 disp(2*nptstot+k1+k2+i,3) = disp(azi_ind_meri(i),3)
!!$              enddo
!!$           endif
!!$
!!$           disp(1:nptstot,1)=disp(1:nptstot,1)
!!$           disp(1:nptstot,2)=disp(1:nptstot,2)
!!$           disp(1:nptstot,3)=disp(1:nptstot,3)
!!$
!!$        enddo
!!$           filename1='Data/snap_mij_cell_'//appmynum2//'_z'
!!$           write(6,*)'filename out vtk :',filename1
!!$           call write_VTK_bin_scal_topology(xtot(1:2*nptstot+k1+k2),ytot(1:2*nptstot+k1+k2),ztot(1:2*nptstot+k1+k2), & 
!!$                disptot_sum(1:2*nptstot+k1+k2,3),(2*nptstot+k1+k2)/4,filename1)           
!!$     enddo
!!$
!!$end subroutine compute_3dslices_snaps
!!$!---------------------------------------------------------------------------------------


real(kind=realkind) function prem(r0,param)
!
! prem model in terms of domains separated by discontinuities
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

real(kind=realkind) , intent(in) :: r0
real(kind=realkind)            :: r,x_prem
real(kind=realkind)             :: ro_prem,vp_prem,vs_prem
character(len=3), intent(in) :: param !rho, vs,vp

  r=r0/1000.
  
  x_prem=r/6371.     ! Radius (normalized to x(surface)=1 )

  IF(r>6356. .and. r .le. 6371.01)THEN        ! upper crustal layer
     ro_prem=2.6
     vp_prem=5.8
     vs_prem=3.2
  ELSEIF(r>6346.6 .and. r .le. 6356.)THEN
     ro_prem=2.9                       ! lower crustal layer
     vp_prem=6.8
     vs_prem=3.9
  ELSEIF(r>6151. .and. r .le. 6346.6)THEN
     ro_prem=2.691+.6924*x_prem             ! upper mantle
     vp_prem=4.1875+3.9382*x_prem
     vs_prem=2.1519+2.3481*x_prem
  ELSEIF(r>5971. .and. r .le. 6151. )THEN
     ro_prem=7.1089-3.8045*x_prem
     vp_prem=20.3926-12.2569*x_prem
     vs_prem=8.9496-4.4597*x_prem
  ELSEIF(r>5771. .and. r .le. 5971.)THEN
     ro_prem=11.2494-8.0298*x_prem
     vp_prem=39.7027-32.6166*x_prem
     vs_prem=22.3512-18.5856*x_prem
  ELSEIF(r>5701. .and. r .le. 5771. )THEN
     ro_prem=5.3197-1.4836*x_prem
     vp_prem=19.0957-9.8672*x_prem
     vs_prem=9.9839-4.9324*x_prem
  ELSEIF(r>5600. .and. r .le. 5701. )THEN   !lower mantle
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=29.2766-23.6027*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=22.3459-17.2473*x_prem-2.0834*x_prem**2+0.9783*x_prem**3
  ELSEIF(r>3630. .and. r .le. 5600. )THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=24.9520-40.4673*x_prem+51.4832*x_prem**2-26.6419*x_prem**3
     vs_prem=11.1671-13.7818*x_prem+17.4575*x_prem**2-9.2777*x_prem**3
  ELSEIF(r>3480. .and. r .le. 3630.)THEN
     ro_prem=7.9565-6.4761*x_prem+5.5283*x_prem**2-3.0807*x_prem**3
     vp_prem=15.3891-5.3181*x_prem+5.5242*x_prem**2-2.5514*x_prem**3
     vs_prem=6.9254+1.4672*x_prem-2.0834*x_prem**2+.9783*x_prem**3
  ELSEIF(r>1221.5 .and. r .le. 3480. )THEN  ! outer core
     ro_prem=12.5815-1.2638*x_prem-3.6426*x_prem**2-5.5281*x_prem**3
     vp_prem=11.0487-4.0362*x_prem+4.8023*x_prem**2-13.5732*x_prem**3
     vs_prem=0.00
  ELSEIF(r .le. 1221.5)THEN                        ! inner core
     ro_prem=13.0885-8.8381*x_prem**2
     vp_prem=11.2622-6.3640*x_prem**2
     vs_prem=3.6678-4.4475*x_prem**2
  ELSE 
     write(6,*)'wrong radius!',r; stop 2
  ENDIF

  if (param=='rho') then
     prem=ro_prem*1000.
  elseif (param=='v_p') then
     prem=vp_prem*1000.
  elseif (param=='v_s') then
     prem=vs_prem*1000.
  else
     write(6,*)'ERROR IN PREM FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function prem
!=============================================================================


!----------------------------------------------------------------------------
subroutine sem_mesh_rotation_tests(mesh)

   use coord_trafo, ONLY : get_r_theta
   include 'mesh_params.h'
!   include 'mesh_params_kernel.h'

  real(kind=realkind), intent(in)  :: mesh(ibeg:iend,ibeg:iend,1:nelem,1:2)
  real(kind=realkind), allocatable, dimension(:,:,:)  :: testu
  real(kind=realkind), allocatable, dimension(:,:,:)  :: testu2
  integer :: iel,ipol,jpol
  real(kind=realkind) :: r,th

 if (have_fluid) then
! testing mesh via simple sine field and solid/fluid distinction
  allocate(testu(ibeg:iend,ibeg:iend,nelem))
  open(unit=52,file='Data/'//'mesh_test'//appmynum//'.dat')
  do iel=1,nelem
     do jpol=ibeg,iend
        do ipol=ibeg,iend
           call get_r_theta(mesh(ipol,jpol,iel,1),mesh(ipol,jpol,iel,2),r,th)
           if (iel<=nel_fluid) then 
              testu(ipol,jpol,iel) = 2.*sin(8.*pi*r/maxr)
           else
              testu(ipol,jpol,iel) = sin(8.*pi*r/maxr)  
           endif
           write(52,*)testu(ipol,jpol,iel)
        enddo
     enddo
  enddo
  close(52)
  open(unit=52,file='Data/'//'mesh_test_flu'//appmynum//'.bindat', &
                                         FORM="UNFORMATTED",STATUS="NEW")
  write(52)testu(:,:,1:nel_fluid)
  close(52)

  open(unit=52,file='Data/'//'mesh_test_sol'//appmynum//'.bindat', &
                                         FORM="UNFORMATTED",STATUS="NEW")
  write(52)testu(:,:,nel_fluid+1:nelem)
  close(52)

  deallocate(testu)


  allocate(testu(ibeg:iend,ibeg:iend,nelem))
  open(unit=52,file='Data/'//'mesh_test_flu'//appmynum//'.bindat', &
                                         FORM="UNFORMATTED",STATUS="OLD")
  read(52)testu(:,:,1:nel_fluid)
  close(52)

  open(unit=52,file='Data/'//'mesh_test_sol'//appmynum//'.bindat', &
                                         FORM="UNFORMATTED",STATUS="OLD")
  read(52)testu(:,:,nel_fluid+1:nelem)
  close(52)

  open(unit=52,file='Data/'//'mesh_test_reread_'//appmynum//'.dat')
  do iel=1,nelem
     do jpol=ibeg,iend
        do ipol=ibeg,iend
           write(52,*)testu(ipol,jpol,iel)
        enddo
     enddo
  enddo
  close(52)

  deallocate(testu)

  allocate(testu2(ibeg:iend,ibeg:iend,nelem))
 
  open(unit=52,file=trim(dir_fwdmesh1(1))//'/Data/'//'velo_sol'//'_'&
                            //appmynum//'_0037.bindat',&
                            FORM="UNFORMATTED",STATUS="OLD",POSITION="REWIND")
  read(52)testu2(:,:,nel_fluid+1:nelem)

  close(52)

  if (have_fluid) then
  open(unit=52,file=trim(dir_fwdmesh1(1))//'/Data/'//'velo_flu'//'_'&
                            //appmynum//'_0037.bindat',&
                            FORM="UNFORMATTED",STATUS="OLD",POSITION="REWIND")
  read(52)testu2(:,:,1:nel_fluid)
  close(52)
  endif 

  open(unit=52,file='Data/'//'mesh_test_fromsolver_'//appmynum//'.dat')
  do iel=1,nelem
     do jpol=ibeg,iend
        do ipol=ibeg,iend
           write(52,*)testu2(ipol,jpol,iel)
        enddo
     enddo
  enddo
  close(52)
  deallocate(testu2)
  endif

end subroutine sem_mesh_rotation_tests

!----------------------------------------------------------------------------

!dk read_sem_mesh_solid-------------------------------------
 subroutine read_sem_meshes_solid(thr)

   include 'mesh_params.h'

   real(kind=realkind), intent(in)                  :: thr
   real(kind=realkind), allocatable, dimension(:,:,:,:) :: mesh_sol
   integer :: it,ndumpstmp
   integer :: ipol,jpol,iel,iproc
   character(len=4) :: appiproc
   real(kind=realkind):: r,th
   real(kind=realkind), allocatable, dimension(:,:,:)  :: testu
   real(kind=realkind), allocatable, dimension(:,:,:)  :: testu2

   nsize=(iend-ibeg+1)**2*nelem

   if (nel_solid+nel_fluid /= nelem) then 
      write(6,*)'solid & fluid elements dont sum up to total!'
      write(6,*)'nel_solid,nel_fluid:',nel_solid,nel_fluid
      write(6,*)'nelem',nelem
      stop
   endif

   open(unit=98,file="input_sem.dat")
   read(98,*)dir_fwdmesh
   read(98,*)dir_bwdmesh
   close(98)

   lffwd = index(dir_fwdmesh,' ')-1
   lfbwd = index(dir_bwdmesh,' ')-1

   if (mynum==0) then 
      write(6,*)'forward calc :',dir_fwdmesh(1:lffwd)
      write(6,*)'backward calc:',dir_bwdmesh(1:lfbwd)
      write(6,*)
      write(6,*)'loading mesh..........'
   endif
   call flush(6)

   allocate(mesh_sol(ibeg:iend,ibeg:iend,1:nel_solid,2))
   open(unit=88,file=trim(dir_fwdmesh1(1))//'/Data'//'/strain_mesh_sol_'&
                             //appmynum//'.dat', &
                             FORM="UNFORMATTED",STATUS="OLD",POSITION="REWIND")
   read(88)mesh_sol
   close(88)

! saving mesh
   if (mynum==0) then 
      write(6,*)'saving mesh..........'
   endif
   call flush(6)

   open(unit=50,file='Data/'//'mesh'//appmynum//'.dat')
   do iel=1,nelem
       do jpol=ibeg,iend
         do ipol=ibeg,iend
            write(50,*)mesh_sol(ipol,jpol,iel,1),mesh_sol(ipol,jpol,iel,2)
         enddo
      enddo
   enddo
   close(50)

   if (mynum==0) write(6,*)'loading dumped wavefield info..........'

   open(unit=97,file=trim(dir_fwdmesh1(1))//'/Data'//'/strain_info.dat0000')
   read(97,*)ndumpstmp
   if (ndumpstmp/=ndumps) then
      write(6,*)'  Problem with number of dumps!'
      write(6,*)'  mesh_params_kernel.h:',ndumps
      write(6,*)'  strain_info.dat0000 :',ndumpstmp
      stop
   endif

   if (mynum==0) then
      write(6,*)
      write(6,*)'number of wavefield snapshots:',ndumps
   endif
   call flush(6)
   
   do it=1,ndumps
      read(97,*)time(it),timestep(it)
   enddo
   close(97)

   if (mynum==0) then
      write(6,*)
      write(6,*)'loading global mesh to each processor........'
   endif
   call flush(6)
   
   allocate(globmesh(ibeg:iend,ibeg:iend,1:nelem,0:nproc_mesh-1,1:2))
   
   do iproc=0,nproc-1

     call define_io_appendix(appiproc,iproc)

     open(unit=88,file=trim(dir_fwdmesh1(1))//'/Data'//'/strain_mesh_sol_'&
                       //appiproc//'.dat', &
                       FORM="UNFORMATTED",STATUS="OLD")
      read(88)mesh_sol(ibeg:iend,ibeg:iend,:,1),mesh_sol(ibeg:iend,ibeg:iend,:,2)
      close(88)

      globmesh(ibeg:iend,ibeg:iend,nel_fluid+1:nelem,iproc,:)=mesh_sol
      
   enddo

   deallocate(mesh_sol)

   if (mynum==0) then
      write(6,*)
      write(6,*)'computing rotated receiver matrix & mesh........' 
   endif
   call flush(6)

   write(6,*)"NEEDS TO BE REDONE!!!!"; stop

 end subroutine read_sem_meshes_solid
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
subroutine save_kernelmesh
    use data_rota, only : elem_d_src,iproc_d_src,ipt_d_src
    use coord_trafo, only : xyz2rthetaphi

    integer :: ipt
    real(kind=realkind) :: r,th,ph

    allocate(idom_ipt(npts))
    do ipt=1,npts
       idom_ipt(ipt) = iidom_glob(elem_d_src(ipt),iproc_d_src(ipt))
    enddo
    deallocate(iidom_glob, elem_d_src)
    !deallocate(ipt_d_src, iproc_d_src)

    open(unit=19,file='Data/kernelmesh_xyz_'//appmynum//'.dat')
    open(unit=20,file='Data/kernelmesh_rthph_'//appmynum//'.dat')
    
    do ipt=1,npts
       write(19,*) real(xgd(ipt)), real(ygd(ipt)), real(zgd(ipt))
       call xyz2rthetaphi(r, th, ph, xgd(ipt), ygd(ipt), zgd(ipt))
       write(20,*) real(r), real(th), real(ph), idom_ipt(ipt)
    enddo
    close(20);close(19)
    deallocate(idom_ipt)
    
end subroutine save_kernelmesh
!----------------------------------------------------------------------------

!-------------------------------------------------------------------------
!> Spreads kernel mesh over processors
subroutine decompose_kernel_mesh
  
  use input_output, only : write_vtk_bin_scal_topology
  use mpi

  integer :: ipt
  character(len=200) :: filename1
  character(len=4) :: appiproc

! decompose the mesh in an embarrassingly, almost 100%-balanced fashion...
  if (mod(npts_tot,nproc)/=0) then 
    if (mynum<nproc-1) then 
       npts=npts_tot/nproc
    else
       npadd=mod(npts_tot,nproc)
       npts=npts_tot/nproc + npadd
    endif
  else
     npts=npts_tot/nproc
  endif

  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  write(6,*)mynum,'reporting: total kernel points:',npts_tot
  write(6,*)mynum,'has',npts,'grid points'
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

  allocate(xgd(npts))
  allocate(ygd(npts))
  allocate(zgd(npts))
  xgd = 0.
  ygd = 0.
  zgd = 0.
  do ipt=1,npts_tot/nproc
     xgd(ipt) = xgd_tot(mynum*(npts_tot/nproc) + ipt)
     ygd(ipt) = ygd_tot(mynum*(npts_tot/nproc) + ipt)
     zgd(ipt) = zgd_tot(mynum*(npts_tot/nproc) + ipt)
  enddo

  if (mynum==nproc-1 .and. mod(npts_tot,nproc)/=0 ) then
     xgd(npts_tot/nproc+1:npts) = xgd_tot(npts_tot-npadd+1:npts_tot)
     ygd(npts_tot/nproc+1:npts) = ygd_tot(npts_tot-npadd+1:npts_tot)
     zgd(npts_tot/nproc+1:npts) = zgd_tot(npts_tot-npadd+1:npts_tot)
  endif
  
  write(6,*)'minmax xgd:',mynum,minval(xgd),maxval(xgd)
  write(6,*)'minmax xgd_tot:',mynum,minval(xgd_tot),maxval(xgd_tot)
  write(6,*)'minmax ygd:',mynum,minval(ygd),maxval(ygd)
  write(6,*)'minmax ygd_tot:',mynum,minval(ygd_tot),maxval(ygd_tot)
  write(6,*)'minmax zgd:',mynum,minval(zgd),maxval(zgd)
  write(6,*)'minmax zgd_tot:',mynum,minval(zgd_tot),maxval(zgd_tot)

  deallocate(xgd_tot,ygd_tot,zgd_tot)

  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  if (mynum==0) write(6,*) 'finished the domain decomposition of the kernel mesh.'
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

  call define_io_appendix(appiproc,mynum)
  filename1='Data/mesh_all_cell_'//appiproc    
  call write_VTK_bin_scal_topology(real(xgd),real(ygd),real(zgd),zgd,npts/4,filename1)

end subroutine decompose_kernel_mesh
!----------------------------------------------------------


!============ CUBED SPHERE ================================

!dk get_cs_mesh_part-------------------------------------------------------
!af here I should add a bookkeeping array to dump any field in the full cubed
!sphere done

!--------------------------------------------------------------------------------------
 subroutine get_cs_mesh_part
  use data_mesh!, only : npts_tot,xgd_tot,ygd_tot,zgd_tot
  use coord_trafo, only : xyz2rthetaphi

  integer :: ipt
  integer :: iel,ipol,jpol,kpol
  real(kind=realkind) :: x, y, z 
  real(kind=realkind) :: r,theta,phi
  real(kind=realkind) :: pi

  pi = 2.*asin(1.)
  call init_params
  call init_cubed_sphere
  call init_sem_grid

  ipt = 0 
  do iel=1,nelt
   do kpol=0,nrpol
    do jpol=0, npol_cs
     do ipol=0, npol_cs
      x = xcol(ipol,jpol,kpol,iel)
      y = ycol(ipol,jpol,kpol,iel) 
      z = zcol(ipol,jpol,kpol,iel)
      call xyz2rthetaphi(r,theta,phi,x,y,z)
! here one defines the actual conditions on the coordinates that must be
! extracted 
      ipt=ipt+1
!      if ( abs(x) < 1. .or. abs(y)<1. .or. abs(z)<1.)  ipt=ipt+1
!     if ( abs(phi) < 1.e-3 .AND. (abs(z) < 1.e-3)) ipt = ipt + 1
     end do
    end do
   end do
  end do
  npts_tot = ipt
  write(6,*) npts_tot
  allocate(xgd_tot(1:npts_tot))
  allocate(ygd_tot(1:npts_tot))
  allocate(zgd_tot(1:npts_tot))
!af Zu
  allocate(iel_cs(1:npts_tot))
  allocate(ipol_cs(1:npts_tot))
  allocate(jpol_cs(1:npts_tot))
  allocate(kpol_cs(1:npts_tot))
!end af Zu 
!
  ipt=0
  do iel=1,nelt
   do kpol=0,nrpol
    do jpol=0, npol_cs
     do ipol=0, npol_cs
      x = xcol(ipol,jpol,kpol,iel)
      y = ycol(ipol,jpol,kpol,iel) 
      z = zcol(ipol,jpol,kpol,iel)
      call xyz2rthetaphi(r,theta,phi,x,y,z)
! af same here 
!     if ( abs(phi) < 1.e-3  ) then
!     if ( abs(phi) < 1.e-3 .AND. (abs(z) < 1.e-3)) then 
!     if ( abs(x) < 1. .or. abs(y)<1. .or. abs(z)<1.)  then
       ipt = ipt+1
       xgd_tot(ipt) = x
       ygd_tot(ipt) = y
       zgd_tot(ipt) = z 
!
       iel_cs(ipt) = iel
       ipol_cs(ipt)=ipol
       jpol_cs(ipt)=jpol
       kpol_cs(ipt)=kpol 
!     end if
     end do
    end do
   end do
  end do
! deallocate(xcol)
! deallocate(ycol)
! deallocate(zcol)
! Not deallocated anymore here, since subsequently used for 3d vtk dump 
  end subroutine get_cs_mesh_part
!--------------------------------------------------------------------------
! 
!dk prepare_dump_cs_part---------------------------------------------------
  subroutine prepare_dump_cs_part
  use data_mesh
  use coord_trafo, only : xyz2rthetaphi

  implicit none
  integer :: ipt
  real(kind=realkind) :: x,y,z
  real(kind=realkind) :: r,theta,phi
!  allocate(thetagd_tot(1:npts_tot))
!  allocate(phigd_tot(1:npts_tot))
  do ipt = 1, npts_tot
   x = xgd_tot(ipt)
   y = ygd_tot(ipt)
   z = zgd_tot(ipt)
   call xyz2rthetaphi(r,theta,phi,x,y,z) 
!   thetagd_tot(ipt) = theta
!   phigd_tot(ipt) = phi
  end do
  end subroutine prepare_dump_cs_part
!--------------------------------------------------------------------------
! 
!dk init_cubed_sphere------------------------------------------------------
  subroutine init_cubed_sphere
! Here we define the skeleton of the grid, as I usually put it.  
  use global_parameters
    use data_mesh

  integer :: i,ii,iii,izone
  real(kind=realkind) :: Xe,Ye,delta
  allocate(xi_el(0:nang),eta_el(0:nang),r_el(0:nr))
  xi_el = 0.
  eta_el = 0.
  r_el = 0.
!
  dang=pi/(2.*real(nang))
  dr=(re-ri)/real(nr)
  do i=0,nang
   xi_el(i)=-pi*0.25+dang*real(i)
   eta_el(i)=-pi*0.25+dang*real(i)
  enddo
  do i=0,nr
   r_el(i)=ri+dr*i
  enddo
  allocate(x_el(0:nang,0:nang,0:nr,1:6))
  allocate(y_el(0:nang,0:nang,0:nr,1:6))
  allocate(z_el(0:nang,0:nang,0:nr,1:6))
  x_el = 0; y_el=0.; z_el = 0.
  allocate(number_el(0:nang,0:nang,0:nr,1:6))
  number_el = 0 
! The mapping is defined following Chaljub (2000). 
  do iii=0,nr    ! loop over r
   do ii=0,nang  ! loop over eta
    do i=0,nang  ! loop over xi
     Xe=tan(xi_el(i))
     Ye=tan(eta_el(ii))
     delta=1+Xe**2+Ye**2
! zone 1
     x_el(i,ii,iii,1)=r_el(iii)*delta**(-0.5)
     y_el(i,ii,iii,1)=r_el(iii)*Xe*delta**(-0.5)
     z_el(i,ii,iii,1)=r_el(iii)*Ye*delta**(-0.5)
! zone 2
     x_el(i,ii,iii,2)=-r_el(iii)*Xe*delta**(-0.5)
     y_el(i,ii,iii,2)=r_el(iii)*delta**(-0.5)
     z_el(i,ii,iii,2)=r_el(iii)*Ye*delta**(-0.5)
! zone 3
     x_el(i,ii,iii,3)=-r_el(iii)*delta**(-0.5)
     y_el(i,ii,iii,3)=-r_el(iii)*Xe*delta**(-0.5)
     z_el(i,ii,iii,3)=r_el(iii)*Ye*delta**(-0.5)
! zone 4
     x_el(i,ii,iii,4)=r_el(iii)*Xe*delta**(-0.5)
     y_el(i,ii,iii,4)=-r_el(iii)*delta**(-0.5)
     z_el(i,ii,iii,4)=r_el(iii)*Ye*delta**(-0.5)
! zone 5
     x_el(i,ii,iii,5)=-r_el(iii)*Ye*delta**(-0.5)
     y_el(i,ii,iii,5)=r_el(iii)*Xe*delta**(-0.5)
     z_el(i,ii,iii,5)=r_el(iii)*delta**(-0.5)
! zone 6
     x_el(i,ii,iii,6)=r_el(iii)*Ye*delta**(-0.5)
     y_el(i,ii,iii,6)=r_el(iii)*Xe*delta**(-0.5)
     z_el(i,ii,iii,6)=-r_el(iii)*delta**(-0.5)
    enddo
   enddo
  enddo
! element global numbering (lexicographical ordering) 
  do izone = 1, 6
   do iii=0,nr-1     ! loop over  r
    do ii=0,nang-1   ! loop over eta
     do i=0,nang-1   ! loop over ksi
      number_el(i,ii,iii,izone) = (izone-1)*(nang*nang*nr)+&
                            iii*(nang*nang)+((ii*nang)+i+1)
     end do
    end do
   end do
  end do
  nelt = maxval(number_el)
  write(6,*)
  write(6,"(10x,'mesh info: ')") 
  write(6,"(10x,'number of spectral elements: ',i8)") nelt
  write(6,"(10x,'number of grid points per spectral element: ',i8)") (nrpol+1)*(npol_cs+1)**2
  write(6,"(10x,'total number of grid points on the order of: ',i8)") nelt*(nrpol+1)*(npol_cs+1)**2
  write(6,*)
  end subroutine init_cubed_sphere
!--------------------------------------------------------------------------------------
!
!dk exit_cubed_sphere
!--------------------------------------------------------------------------------------
  subroutine exit_cubed_sphere
    use data_mesh, only : number_el,Z_el,Y_el,X_el,r_el,eta_el,xi_el

  deallocate(number_el) 
  deallocate(Z_el,Y_el,X_el) 
  deallocate(r_el,eta_el,xi_el)
  end subroutine exit_cubed_sphere
!--------------------------------------------------------------------------------------

!
!dk init_sem_grid
!--------------------------------------------------------------------------------------
  subroutine init_sem_grid

! Here the full grid is defined. 
  use splib 
  use data_mesh

  real(kind=realkind), dimension(:), allocatable :: dxi,dxk
  integer :: izone,i,ii,iii
  integer :: ipol, jpol, kpol, iel
  real(kind=realkind) :: tksi,teta,tr,C,D,delta,Xe,Ye,J
  allocate(xi_i(0:npol_cs),xi_j(0:npol_cs),xi_k(0:nrpol))
  allocate(wt_i(0:npol_cs),wt_j(0:npol_cs),wt_k(0:nrpol))
!
! dummies needed by spectral library
  allocate(dxi(0:npol_cs),dxk(0:nrpol))
!
  call ZELEGL(npol_cs,dble(xi_i),dble(dxi))
  call get_welegl(npol_cs,dble(xi_i),dble(wt_i))
! xi_i -> npol_cs+1 pts de quadrature pour l'ordre ploynomial npol_cs
  call ZELEGL(npol_cs,dble(xi_j),dble(dxi))
  call get_welegl(npol_cs,dble(xi_j),dble(wt_j))
  call ZELEGL(nrpol,dble(xi_k),dble(dxk))
  call get_welegl(nrpol,dble(xi_k),dble(wt_k))
! get rid of dummies
  deallocate(dxk,dxi)
!
! Now define collocation points for the SEM-cubed sphere grid
  allocate(xcol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
  allocate(ycol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
  allocate(zcol(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
  allocate(jacob(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
  allocate(massmat(0:npol_cs,0:npol_cs,0:nrpol,1:nelt))
  do izone=1,6         ! loop over the chunks
   do iii=0,nr-1       ! loop over r   (the global r)
    do ii=0,nang-1     ! loop over eta (the global eta)
     do i=0,nang-1     ! loop over xi  (the global xi)
      iel = number_el(i,ii,iii,izone)
      do kpol = 0, nrpol  ! loop over the elemental k index (r direction)
       do jpol = 0, npol_cs  ! loop over the elemental j index (eta direction)
        do ipol = 0, npol_cs ! loop over the elemental i index (xi direction)
         tksi= xi_el(  i) +(1 + xi_i(ipol))*.5*dang
         teta=eta_el( ii) +(1 + xi_j(jpol))*.5*dang
         tr=    r_el(iii) +(1 + xi_k(kpol))*.5*dr
         Xe=tan(tksi)
         Ye=tan(teta)
         C=(1+Xe**2)**(0.5)
         D=(1+Ye**2)**(0.5)
         delta=1+Xe**2+Ye**2
         J=(tr**2)*(C**2)*(D**2)*(delta**(-1.5))
         J=J*((.5*dang)**2)*.5*dr
         jacob(ipol,jpol,kpol,iel) = J
         massmat(ipol,jpol,kpol,iel) = wt_i(ipol)*wt_j(jpol)*wt_k(kpol)*J
         if(izone==1) then
          Xcol(ipol,jpol,kpol,iel)=tr*delta**(-0.5)
          Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5)
          Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5)
         endif
         if(izone==2) then
          Xcol(ipol,jpol,kpol,iel)=-tr*Xe*delta**(-0.5)
          Ycol(ipol,jpol,kpol,iel)=tr*delta**(-0.5)
          Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5)
         endif
         if(izone==3) then
          Xcol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5)
          Ycol(ipol,jpol,kpol,iel)=-tr*Xe*delta**(-0.5)
          Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5)
         endif
         if(izone==4) then
          Xcol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5)
          Ycol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5)
          Zcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5)
         endif
         if(izone==5) then
          Xcol(ipol,jpol,kpol,iel)=-tr*Ye*delta**(-0.5)
          Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5)
          Zcol(ipol,jpol,kpol,iel)=tr*delta**(-0.5)
         endif
         if(izone==6) then
          Xcol(ipol,jpol,kpol,iel)=tr*Ye*delta**(-0.5)
          Ycol(ipol,jpol,kpol,iel)=tr*Xe*delta**(-0.5)
          Zcol(ipol,jpol,kpol,iel)=-tr*delta**(-0.5)
         endif
        end do
       end do
      end do
     end do
    end do
   end do
  end do

  deallocate(massmat)
  deallocate(jacob)
  deallocate(zcol,ycol,xcol)
  deallocate(wt_k,wt_j,wt_i)
  deallocate(xi_k,xi_j,xi_i)

  end subroutine init_sem_grid
!--------------------------------------------------------------------------------------
!
!dk init_params
 !--------------------------------------------------------------------------------------
 subroutine init_params
   use data_mesh

  implicit none
  open(5,file='inparam_cs')
  read(5,*) ri
  read(5,*) re
  read(5,*) nang
  read(5,*) nr
  read(5,*) npol_cs
  read(5,*) nrpol
  close(5)
  write(6,*)
  write(6,"(10x,'cubed sphere mesh parameters')")
  write(6,*)
  write(6,"(10x,'inner radius of the shell  : ri = ',1pe12.5)") ri
  write(6,"(10x,'outer radius of the shell  : re = ',1pe12.5)") re
  write(6,*)
  write(6,"(10x,'number of horizontal segments along a cube edge: nang = ',i4)") nang
  write(6,"(10x,'number of segments in the  radial direction    :   nr = ',i4)") nr
  write(6,"(10x,'horizontal polynomial order : npol_cs = ',i4)") npol_cs
  write(6,"(10x,'radial polynomial order : nrpol = ',i4)") nrpol
  write(6,*)
  end subroutine init_params
!--------------------------------------------------------------------------------------
!
!============  END CUBED SPHERE ================================


 
!  ~~~~~~~~~~~~~~~ ANALYTICAL MESH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!---------------------------------------------------------------------
  subroutine define_mymesh(thr)

  use global_parameters
  use data_mesh, only : npts_tot,xgd_tot,ygd_tot,zgd_tot
  use coord_trafo, only : rthetaphi2xyz

  real(kind=realkind), intent(in) :: thr
  real(kind=realkind) :: r,theta,phi,x,y,z
  integer :: ir, itheta, iphi
  real(kind=realkind) :: ri,ro
  integer :: thetar
  real(kind=realkind) :: phibeg,phiend
  real(kind=realkind) :: thetabeg, thetaend
  integer :: mynr,mynphi,myntheta
  real(kind=realkind) :: dr, dphi, dtheta
  integer :: ipt 
  character(len=200) :: meshtype

  meshtype='halfplane'

  select case(meshtype)

  case('halfplane')
     mynr = 80 
     mynphi = 180
     phibeg = -.5*asin(1.)
     phiend =  .5*asin(1.)
     ro = 6371.e3
     ri = 2371.e3
     phi=0
     dr = (ro -ri)/real(mynr)
     dphi=(phiend-phibeg)/real(mynphi)
     npts_tot = (mynr+1)*(mynphi+1)        
     write(6,*) ' number of points ', npts_tot 
     allocate(xgd_tot(1:npts_tot)) ; xgd_tot = 0. 
     allocate(ygd_tot(1:npts_tot)) ; ygd_tot = 0. 
     allocate(zgd_tot(1:npts_tot)) ; zgd_tot = 0. 
     ipt = 0
     theta = thr/2. 
     do ir = 0, mynr
        r = ri + real(ir) * dr
        do iphi = 0, mynphi
           phi = phibeg + real(iphi)*dphi
           ipt = ipt + 1 
           call rthetaphi2xyz(x,y,z,r,theta,phi)
           xgd_tot(ipt) = x
           ygd_tot(ipt) = y
           zgd_tot(ipt) = z
        end do
     end do

  case('meridional')
     mynr = 80 
     myntheta = 360
     thetabeg = 0.
     thetaend = 2.*asin(1.) 
     ro = 6371.e3
     ri = 2371.e3
     phi = 0. 
     dr = (ro -ri)/real(mynr)
     dtheta = (thetaend-thetabeg)/real(myntheta)
     npts_tot = (mynr+1)*(myntheta+1) 
     write(6,*) ' number of points ', npts_tot 
     allocate(xgd_tot(1:npts_tot)) ; xgd_tot = 0. 
     allocate(ygd_tot(1:npts_tot)) ; ygd_tot = 0. 
     allocate(zgd_tot(1:npts_tot)) ; zgd_tot = 0. 
     do ir = 0, mynr
        r = ri + real(ir) * dr
        do itheta = 0, myntheta
           theta = thetabeg + real(itheta)*dtheta
           ipt=ipt+1
           call rthetaphi2xyz(x,y,z,r,theta,phi)
           xgd_tot(ipt) = x
           ygd_tot(ipt) = y
           zgd_tot(ipt) = z
        end do
     end do
     
  end select
  
  end subroutine define_mymesh
!---------------------------------------------------------------------

!---------------------------------------------------------------------
  subroutine define_mymesh_dbg(thr)
  use global_parameters
  use data_mesh, only : npts_tot,xgd_tot,ygd_tot,zgd_tot
  use coord_trafo, only : rthetaphi2xyz

  real(kind=realkind), intent(in) :: thr
  real(kind=realkind) :: r,theta,phi,x,y,z
  integer :: ir, itheta, iphi
  real(kind=realkind) :: ri,ro
  integer :: thetar
  real(kind=realkind) :: phibeg,phiend
  real(kind=realkind) :: thetabeg, thetaend
  integer :: mynr,mynphi,myntheta
  real(kind=realkind) :: dr, dphi, dtheta
  integer :: ipt 
  character(len=200) :: meshtype

  meshtype='onesinglepoint'

  select case(meshtype)

  case('onesinglepoint')
     ro = 6371.e3
     ri = 2371.e3
     npts_tot = 1        
     write(6,*) ' number of points ', npts_tot 
     allocate(xgd_tot(1:npts_tot)) ; xgd_tot = 0. 
     allocate(ygd_tot(1:npts_tot)) ; ygd_tot = 0. 
     allocate(zgd_tot(1:npts_tot)) ; zgd_tot = 0. 
     r = 5500.e3
     phi = 0. 
     theta = .5*thr
     call rthetaphi2xyz(x,y,z,r,theta,phi)
     xgd_tot(1) = x
     ygd_tot(1) = y
     zgd_tot(1) = z
  end select
  
  end subroutine define_mymesh_dbg
!---------------------------------------------------------------------
!  ~~~~~~~~~~~~~~~  END ANALYTICAL MESH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!========================
end module kernel_meshes
!========================

 
