!====================
module coord_mapping
!====================

  use data_mesh
  use global_parameters
  use data_rota
  use coord_trafo

  implicit none


  public :: find_rotated_gridpoints_1d, rthetaphi2xyz, rotate_frame_rd
  public :: find_gridpoints_semmesh, find_gridpoints_1d, compute_srdphirdzrd
  public :: define_kern_sem_mesh_mapping_arrays

  private

  contains 
!/////////////////////////////////////////////////

!-------------------------------------------------------------------------
subroutine find_gridpoints_semmesh

  use mpi
  integer               :: ipt
  real(kind=realkind)   :: sgd

  allocate(iproc_d_src(npts))
  allocate(elem_d_src(npts))
  allocate(ipt_d_src(npts))

  ! completing parallelization here in this specific SEM mesh case,
  ! for which it is easier to have had some index mapping predefined globally 

  do ipt=1, npts_tot / nproc
     iproc_d_src(ipt) = iproc_proc_tot(mynum * npts_tot / nproc + ipt)
     elem_d_src(ipt) = iel_proc_tot(mynum * npts_tot / nproc + ipt)
     ipt_d_src(ipt) = ipt_proc_tot(mynum * npts_tot / nproc + ipt)
  enddo

  if (mynum == nproc-1 .and. mod(npts_tot,nproc) /= 0 ) then
     iproc_d_src(npts_tot/nproc+1:npts) = iproc_proc_tot(npts_tot-npadd+1:npts_tot)
     elem_d_src(npts_tot/nproc+1:npts) = iel_proc_tot(npts_tot-npadd+1:npts_tot)
     ipt_d_src(npts_tot/nproc+1:npts) = ipt_proc_tot(npts_tot-npadd+1:npts_tot)
  endif

  ! denounce remnants of the global world
  deallocate(iproc_proc_tot, ipt_proc_tot, iel_proc_tot)
  allocate(phigd(npts))
  phigd(:) = phi0

  ! define index arrays for mapping from kernel to sem mesh
  call define_kern_sem_mesh_mapping_arrays('fwd', ipt_d_src, iproc_d_src, &
                                           npts_iproc_fwd, num_meshes_fwd, &
                                           meshes_fwd, iptkern_proc_src, &
                                           iptsem_proc_src)
  !deallocate(ipt_d_src,iproc_d_src)
  
  write(6,*)mynum,'finished search for fwd/src grid points.'
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

end subroutine find_gridpoints_semmesh
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!> Makes the mapping arrays iptkern_proc_src and iptsem_proc_src, used to identify
!! from which SOLVER processor to read. 
!! Probably this could be scrapped.
subroutine define_kern_sem_mesh_mapping_arrays(casestrg, ipt_proc, iproc_map, &
                                               npts_iproc1, num_meshes1, &
                                               meshes1, iptkern, iptsem)

  use mpi

  character(len=3),intent(in)                       :: casestrg    !< src or rec, just for Output
  integer, intent(in), dimension(1:npts)            :: ipt_proc    !< Element number of ipt on its SEM processor
  integer, intent(in), dimension(1:npts)            :: iproc_map   !< Processor number of point ipt
  integer, intent(out)                              :: num_meshes1 !< Number of SEM meshes from which this KERNER processor needs elements
  integer, intent(out), allocatable, dimension(:)   :: npts_iproc1 !< How many elements has this KERNER proc on each SEM mesh? dimension(num_meshes1)
  integer, intent(out), dimension(nproc_mesh)       :: meshes1     !
  integer, intent(out), allocatable, dimension(:,:) :: iptkern     !< Element number of point on SEM mesh, known as iptkern_proc_rec or _src in the rest of the program
  integer, intent(out), allocatable, dimension(:,:) :: iptsem      !< Element number of point on SEM mesh, known as iptsem_proc_rec or _src in the rest of the program
  integer, dimension(:), allocatable                :: iipt_count
  integer                                           :: iproc, ipt, i, imesh
  integer, dimension(0:nproc_mesh-1)                :: invmeshes
  
  ! defining mapping arrays to access the different decomposed SEM meshes
  num_meshes1 = 0
  meshes1 = -1
  invmeshes = -1
  do iproc=0, nproc_mesh-1
      if ( minval(abs(iproc_map-iproc)) == 0 ) then 
          ! Read: If at least one element of
          ! iproc_map and iproc are identical
          ! Has SEM mesh iproc elements of the
          ! kerner mesh?
          num_meshes1 = num_meshes1 + 1
          meshes1(num_meshes1) = iproc 
          invmeshes(iproc) = num_meshes1 
      endif
  enddo

  allocate(npts_iproc1(num_meshes1))
  npts_iproc1 = 0
  do ipt=1, npts
      npts_iproc1(invmeshes(iproc_map(ipt))) = npts_iproc1(invmeshes(iproc_map(ipt))) + 1
  enddo

  ! just some info to stdio
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
  if (mynum==0) write(6,*)

  do iproc=0,nproc-1
     if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
     if (iproc == mynum) then
        write(6,87) mynum, npts, num_meshes1, casestrg
        do i=1, num_meshes1
           write(6,88) mynum, npts_iproc1(i), casestrg, meshes1(i)
        enddo
        write(6,*)
     endif
   enddo
   
87 format('Processor',i4,' has', i9,' points in',i4,' SEM ',a4,' mesh(es)')
88 format('Processor',i4,' has', i9,' points in SEM ',a4,' mesh #',i4)

  allocate(iptkern(maxval(npts_iproc1), num_meshes1))
  allocate(iptsem( maxval(npts_iproc1), num_meshes1))
  allocate(iipt_count(num_meshes1))

  iptkern = -1
  iptsem = -1
  iipt_count = 0
  do imesh=1, num_meshes1
      do ipt=1, npts
          if ( iproc_map(ipt) == meshes1(imesh) ) then
              iipt_count(imesh) = iipt_count(imesh) + 1

              if (iipt_count(imesh) > maxval(npts_iproc1)) then
                 write(6,*) ' a bug'
                 write(6,*) mynum, iproc, ipt, iipt_count(imesh)
                 stop
              endif
              
              ! iptkern is defined as follows:
              ! Arguments are 
              ! 1) the count through the amount of points 
              !    for a given SEM mesh (as decomposed by the forward simulation)
              ! 2) the counting number of that mesh (up until the total number of 
              !    SEM meshes from which points are read)
              ! Value is: index of the point in the kernel mesh
              iptkern(iipt_count(imesh),imesh) = ipt

              ! iptsem is defined as follows:
              ! Argument is: the count through the amount of points for a given 
              !               SEM mesh (as decomposed by the forward simulation)
              ! Value is: index of the point in the SEM mesh
              iptsem(iipt_count(imesh),imesh) = ipt_proc(ipt)

          endif
      enddo
  enddo

  do iproc=1, num_meshes1
     if (minval(iptkern(1:npts_iproc1(iproc),iproc)) < 0) then
        write(6,*) mynum, 'Problem with definition of kernel ipt mapping!', casestrg
        write(6,*) mynum, '...undefined point:', &
             minloc(iptkern(1:npts_iproc1(iproc),iproc))
        write(6,*) mynum, ipt, npts_iproc1(iproc), iproc
        stop
     endif
     if (minval(iptsem(1:npts_iproc1(iproc),iproc)) < 0) then
        write(6,*) mynum, 'Problem with definition of sem ipt mapping!', casestrg
        write(6,*) mynum, '...undefined point:', &
             minloc(iptsem(1:npts_iproc1(iproc),iproc))
        write(6,*) mynum, ipt, npts_iproc1(iproc), iproc
        stop
     endif
  enddo

  deallocate(iipt_count)

  if (maxval(iptkern) > npts) then 
      write(6,*) mynum, 'Problem with definition of kernel ipt mapping:', &
                 maxval(iptkern), maxloc(iptkern), npts, casestrg
      stop
  endif
  if (maxval(iptsem) > nsize) then 
      write(6,*) mynum, 'Problem with definition of sem ipt mapping:', &
                 maxval(iptsem), maxloc(iptsem), nsize, casestrg
      stop
  endif

  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

end subroutine define_kern_sem_mesh_mapping_arrays
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine find_gridpoints_1d(mesh)

  use data_rota
  use mpi
  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  real(kind=realkind),intent(in)    :: mesh(ibeg:iend,ibeg:iend,1:nelem,0:nproc_mesh-1,1:2)
  real(kind=realkind)               :: xgr,ygr,zgr ! global coordinates of the receiver
  integer                           :: ipt,ipol,jpol,iel,iproc,i
  real(kind=realkind)               :: sgd,dr,r,th
  real(kind=realkind)               :: dr2, deltas, deltaz
  integer                           :: ipol_min,jpol_min,iel_min,iproc_min,ipt_count,ipt_count_min
  integer                           :: thetadisttmp(2),iel3(3)
  logical                           :: foundit
  integer, dimension(:), allocatable                :: iipt_count
  real(kind=realkind), dimension(:), allocatable    :: mindist, mindist2

!  write(6,*)mynum,' finding fwd grid points....'
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

  allocate(mindist(npts))
  allocate(ipt_d_src(npts))
  allocate(iproc_d_src(npts))
  allocate(elem_d_src(npts))
  allocate(mindist2(npts))
  allocate(phigd(npts))

  open(unit=1000+mynum,file='Info/min_dist_fwd_'//appmynum//'.dat')

  write(6,*)'nproc_mesh:',nproc_mesh

! find corresponding points in the unrotated source mesh
!=============================
  do ipt=1,npts
!=============================
     if ((npts>10) .and. mod(ipt,npts/10) == 0 .and. mynum==0) &
        write(6,*) 100 * ipt / npts, ' % done in finding fwd grid points'
     call xyz2sphiz(sgd,phigd(ipt),xgd(ipt),ygd(ipt),zgd(ipt))

     mindist2(ipt) = maxr**2
     mindist(ipt) = maxr
     foundit=.false.
     ! Loop over all (forward) processors
     do iproc=0, nproc_mesh-1
        ! Loop over all elements 
        do iel = 1, nelem
            ! If ds or dz is larger than size of largest element hmax
            ! this is not the closest element.
            deltas = abs(mesh(npol/2,npol/2,iel,iproc,1) - sgd)
            if (deltas>hmax*5) cycle
            deltaz = abs(mesh(npol/2,npol/2,iel,iproc,2) - zgd(ipt))
            if (deltaz>hmax*5) cycle
           
            dr2 = deltas**2 + deltaz**2

!            dr = sqrt(( mesh(npol/2,npol/2,iel,iproc,1) - sgd  )**2+ &
!                       ( mesh(npol/2,npol/2,iel,iproc,2) - zgd(ipt)  )**2)
            if (dr2<mindist2(ipt)) then
               iel_min=iel  
               iproc_min = iproc 
               mindist2(ipt) = dr2
               !mindist(ipt) = sqrt(dr2)
               !write(666,*)mynum,ipt,mindist(ipt),iel_min,iproc_min
               foundit=.true.
            end if 
            if (mynum==0 .and. ipt==1) then
                write(4000+iproc,*) &
                mesh(npol/2,npol/2,iel,iproc,1),0.,mesh(npol/2,npol/2,iel,iproc,2) 
            end if
        end do !iel
     end do !iproc
     if (.not. foundit) &
        write(6,*)mynum,'hasnt found a point for',sgd,zgd(ipt),mindist(ipt),ipt

     ! Loop over all GLL points of closest element
     mindist(ipt) = maxr
     ipt_count = (iel_min-1)*(iend-ibeg+1)**2
     do jpol=ibeg,iend
        do ipol=ibeg,iend
           ipt_count = ipt_count + 1
           dr = sqrt(( mesh(ipol,jpol,iel_min,iproc_min,1) - sgd      )**2 + &
                     ( mesh(ipol,jpol,iel_min,iproc_min,2) - zgd(ipt) )**2)
           if (dr<=mindist(ipt)) then
              mindist(ipt) = dr
              ipol_min = ipol
              jpol_min = jpol  
              ipt_count_min = ipt_count
           end if
        end do
     end do

     write(1000+mynum,*)sgd,zgd(ipt),mindist(ipt), &
     mesh(ipol_min,jpol_min,iel_min,iproc_min,1), &
     mesh(ipol_min,jpol_min,iel_min,iproc_min,2)

     ipt_d_src(ipt) = ipt_count_min
     iproc_d_src(ipt) = iproc_min
     elem_d_src(ipt) = iel_min
  end do ! loop over the externally defined grid points 
  
  deallocate(mindist)

! defining mapping arrays to access the different decomposed SEM meshes
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
  write(6,*)
  write(6,*)mynum,'defining mapping between kernel and fwd/src mesh'
  write(6,*)mynum,'number of kernel grid points:',npts
  write(6,*)mynum,'number of cell points:',nloc

  call define_kern_sem_mesh_mapping_arrays('fwd', ipt_d_src, iproc_d_src,  &
                                           npts_iproc_fwd, num_meshes_fwd, &
                                           meshes_fwd, iptkern_proc_src,   &
                                           iptsem_proc_src)

  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
  write(6,*)mynum,'finished search for fwd grid points.'
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

end subroutine find_gridpoints_1d
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
subroutine find_rotated_gridpoints_1d(rgr,thetagr,phigr,mesh)

  use data_rota
  use mpi
  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  real(kind=realkind),intent(in) ::rgr,thetagr,phigr
  real(kind=realkind),intent(in) ::mesh(ibeg:iend,ibeg:iend,1:nelem,0:nproc_mesh-1,1:2)

  real(kind=realkind)   :: xgr, ygr, zgr ! global coordinates of the receiver
  integer               :: ipt, ipol, jpol, iel, iproc, i
  real(kind=realkind)   :: dr, r, th, dr2, deltas, deltaz
  real(kind=realkind), allocatable :: srd(:),zrd(:)
  integer               :: ipol_min, jpol_min, iel_min, iproc_min
  integer               :: ipt_count, ipt_count_min
  integer, dimension(:), allocatable :: iipt_count
  real(kind=realkind), dimension(:), allocatable :: mindist, mindist2
  integer               :: thetadisttmp(2),iel3(3)
  logical               :: foundit

!  write(6,*)mynum,' finding rotated bwd grid points....'
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  allocate(mindist(npts))
  allocate(mindist2(npts))
  allocate(ipt_d_rec(npts))
  allocate(iproc_d_rec(npts))
  allocate(phird(1:npts), srd(1:npts), zrd(1:npts))

  call rthetaphi2xyz(xgr,ygr,zgr,rgr,thetagr,phigr)
  open(unit=1000+mynum,file='Info/min_dist_bwd_'//appmynum//'.dat')

! W contains the coordinates for the vtk binary outputs
  allocate(W_vtk(1:npts,3))
  W_vtk(1:npts,1)=xgd 
  W_vtk(1:npts,2)=ygd
  W_vtk(1:npts,3)=zgd
  call rotate_frame_rd(npts,srd,phird,zrd,xgd,ygd,zgd,phigr,thetagr)
  if (.not. cell_topology)  deallocate(xgd,ygd,zgd)

! find corresponding points in the rotated receiver mesh
!=============================
  do ipt=1,npts ! Points in the kerner mesh
!=============================
!af
     if ( (npts>10) .AND. mod(ipt,npts/10)==0 .and. mynum==0) &
     write(6,*) 100*ipt/npts, ' % done in finding rotated bwd grid points'
     mindist(ipt) = maxr
     mindist2(ipt) = maxr**2
     foundit=.false.
     do iproc=0, nproc_mesh-1
         do iel = 1, nelem ! Points in the solver mesh
            ! If ds or dz is larger than size of largest element hmax
            ! this is not the closest element.
            deltas = abs(mesh(npol/2, npol/2, iel, iproc, 1) - srd(ipt))
            if (deltas.gt.hmax*5) cycle
            deltaz = abs(mesh(npol/2, npol/2, iel, iproc, 2) - zrd(ipt))
            if (deltaz.gt.hmax*5) cycle
           
            dr2 = deltas**2 + deltaz**2
           ! dr = sqrt(( mesh(npol/2,npol/2,iel,iproc,1) - srd(ipt)  )**2+ &
           !            ( mesh(npol/2,npol/2,iel,iproc,2) - zrd(ipt)  )**2)
            if (dr2<mindist2(ipt)) then
                iel_min = iel 
                iproc_min = iproc
                mindist2(ipt) = dr2
                mindist(ipt) = sqrt(dr2)
                foundit=.true.
            end if 
            if (mynum==0 .and. ipt==1) then
               write(4000+iproc,*) &
               mesh(npol/2,npol/2,iel,iproc,1),mesh(npol/2,npol/2,iel,iproc,2) 
            end if
         end do !iel
     end do !iproc

     if (.not. foundit) &
        write(6,*)mynum,'hasnt found a point for',srd(ipt),zrd(ipt),mindist(ipt),ipt
     
     ! Loop over all GLL points of closest element
     ipt_count = (iel_min-1)*(iend-ibeg+1)**2
     ipt_count_min = ipt_count
     mindist(ipt) = maxr

     do jpol=ibeg,iend
        do ipol=ibeg,iend
           ipt_count = ipt_count + 1
              dr = sqrt(( mesh(ipol,jpol,iel_min,iproc_min,1) - srd(ipt) )**2 + &
                        ( mesh(ipol,jpol,iel_min,iproc_min,2) - zrd(ipt) )**2)
           if (dr<=mindist(ipt)) then
              mindist(ipt) = dr
              ipol_min = ipol
              jpol_min = jpol 
              ipt_count_min = ipt_count
           end if
        end do
     end do
    ! if ((1+(ipol_min-ibeg)+ (jpol_min-ibeg)*(iend-ibeg+1) + &
    ! (iel_min-1)*(iend-ibeg+1)**2).ne.ipt_count_min) then
    !     write(6,*) ipt, iproc_min, ((ipol_min-ibeg)+ (jpol_min-ibeg)*(iend-ibeg+1) + &
    !                (iel_min-1)*(iend-ibeg+1)**2), ipt_count_min
    ! end if
     write(9000+mynum,*) iel_min, iproc_min, ipt_count_min

     write(1000+mynum,*) srd(ipt), zrd(ipt), mindist(ipt), &
     mesh(ipol_min,jpol_min,iel_min,iproc_min,1), &
     mesh(ipol_min,jpol_min,iel_min,iproc_min,2)

     ipt_d_rec(ipt)   = ipt_count_min
     iproc_d_rec(ipt) = iproc_min
     write(8000+mynum,*)iproc_d_rec(ipt)

  end do ! loop over the externally defined grid points 
  

! defining mapping arrays to access the different decomposed SEM meshes
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
  write(6,*)
  write(6,*)mynum,'defining mapping between kernel and bwd/rec mesh'
  call define_kern_sem_mesh_mapping_arrays('bwd', ipt_d_rec, iproc_d_rec, &
                                           npts_iproc, num_meshes, meshes,&
                                           iptkern_proc_rec, iptsem_proc_rec)
!  Some output to test the SEM mesh/Kerner mesh assignments. May be deleted in the next
!  Revision, once we are sure that it is working. SCS 14.10.12
!  do ipt = 1, npts
!      print *, ipt, ipt_d_rec(ipt), iptsem_proc_rec(ipt), iptkern_proc_rec(ipt,1),iptkern_proc_rec(ipt,2)
!  end do

!     allocate(mesh1d(nsize,0:nproc_mesh-1,2))
!     mesh1d=0.
!     do iproc=0, nproc_mesh - 1
!         i=0
!         do iel=1, nelem
!             do jpol=ibeg, iend
!                 do ipol=ibeg, iend
!                     i = i + 1
!                     mesh1d(i,iproc,1:2) = globmesh(ipol,jpol,iel,iproc,1:2)
!                 enddo
!             enddo
!         enddo
!     enddo
!
!  do ipt=1,npts_iproc(1)
!      if (iptkern_proc_rec(ipt,iproc_d_rec(ipt)+1)<1) cycle
!      write(12345,*) ipt, & 
!                     mesh1d(iptsem_proc_rec(ipt,iproc_d_rec(ipt)+1), iproc_d_rec(ipt) , 1:2), &
!                     srd(iptkern_proc_rec(ipt,iproc_d_rec(ipt) + 1)), &
!                     zrd(iptkern_proc_rec(ipt,iproc_d_rec(ipt) + 1))
!  end do
!  
!  deallocate(ipt_d_rec,iproc_d_rec)
  deallocate(mindist,srd,zrd)

  write(6,*)mynum,'finished search for bwd grid points.'
  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

end subroutine find_rotated_gridpoints_1d
!-------------------------------------------------------------------------

!!!!-------------------------------------------------------!!!
subroutine rotate_frame_rd(npts,srd,phird,zrd,xgd,ygd,zgd,phigr,thetagr)

use input_output, only : write_VTK_bin_scal

implicit none
integer, intent(in) :: npts
real(kind=realkind), dimension(npts), intent(in) :: xgd,ygd,zgd
real(kind=realkind), intent(in) :: phigr,thetagr
real(kind=realkind), dimension(npts) :: xp,yp,zp,xp_cp,yp_cp,zp_cp
real(kind=realkind), dimension(npts), intent(out) :: srd,zrd,phird
real(kind=realkind) :: phi_cp,rgd,thetagd
integer :: i
character(len=55) :: filename
double precision :: thgr_dble,phgr_dble

thgr_dble=dble(thetagr)
phgr_dble=dble(phigr)

!!first rotation (longitude)
xp_cp=xgd*dcos(phgr_dble)+ygd*dsin(phgr_dble)
yp_cp=-xgd*dsin(phgr_dble)+ygd*dcos(phgr_dble)
zp_cp=zgd

!second rotation (colat)
xp=xp_cp*dcos(thgr_dble)-zp_cp*dsin(thgr_dble)
yp=yp_cp
zp=xp_cp*dsin(thgr_dble)+zp_cp*dcos(thgr_dble)

srd=dsqrt(dble(xp)**2+dble(yp)**2)
zrd=zp
do i=1,npts
   phi_cp=atan2(yp(i),xp(i))
   if (phi_cp<0.d0) then
      phird(i)=2.d0*pi+phi_cp
   else
      phird(i)=phi_cp
   endif
   if (phgr_dble==0.0 .and. ygd(i)==0.0)  phird(i)=0.
enddo

write(6,*)'Done with rotating frame rd.'

end subroutine rotate_frame_rd

!-----------------------------------------------------------------------
subroutine compute_srdphirdzrd(srd,phird,zrd,xgd,ygd,zgd, &
                               xorig,yorig,zorig,phigr,thetagr)
  real(kind=realkind), intent(out)::  srd,phird,zrd
  real(kind=realkind), intent(in) :: xgd,ygd,zgd
  real(kind=realkind), intent(in) :: xorig,yorig,zorig
  real(kind=realkind), intent(in) :: phigr,thetagr
  real(kind=realkind) :: rgd,phigd,thetagd,alpha
  real(kind=realkind) :: arr

!  thetagd and phigd as computed from
!  xgd(ipt),ygd(ipt),zgd(ipt), using the xyz2rthetaphi routine. 
  call xyz2rthetaphi(rgd,thetagd,phigd,xgd,ygd,zgd)

!  In the case of a source at the NP .AND. the z-component of the receiver, 
!  things are fine.

!  find closest gll point for that specific diffracting point in the D domain
!  for that specific diffracting point

  arr = (xorig**2+yorig**2+zorig**2)**(-.5)

  srd  =arr*( (ygd*zorig-zgd*yorig)**2 + &
             (zgd*xorig-xgd*zorig)**2 + &
             (xgd*yorig-ygd*xorig)**2     )**(.5)
  zrd = arr*(xgd*xorig+ygd*yorig+zgd*zorig)

! phird
  alpha = abs(phigr-phigd)
  phird = atan(sin(alpha)/(sin(thetagr)/tan(thetagd) -cos(thetagr)*cos(alpha)))
  if (abs(phird)<smallval_dble) phird = 0.d0

end subroutine compute_srdphirdzrd

end module coord_mapping
