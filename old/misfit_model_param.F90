!========================
module misfit_model_param
!========================

  use data_mesh
  use global_parameters

  implicit none

  public :: init_time_window_misfit_medium,allocate_arrays
  public :: init_wavefields,misfit_kernels,waveform_kernels

  private 

  contains

!-----------------------------------------------------------------------------------------------------
subroutine init_time_window_misfit_medium(thr,phr)

  use data_misfits
  use data_posteriori, only : post_per,filter_type_posteriori,do_filt
! use posteriori_misfit_operations, only : read_posteriori
  use data_fft, only : uspec1d,nomega,omega,filter_what,filter_period_hi,filter_period_low,filter_type
  use filtering, only : define_filter,time_shift_vector,filter_vector
  use fft, only : fftf1d,fftb1d
  
  real(kind=realkind), intent(in)          :: thr,phr  
  complex(kind=8) :: uspec(0:nomega),vspec(0:nomega),dspec(0:nomega)
  character(len=200) :: appifilt1,appifilt2
  integer :: ifilt,i,iwin,it

if (mynum==0) write(6,*)'  Loading/computing synthetics, data, time window, medium factors................'

  if (times_medium) then 
! calculate medium parameters to later multiply with
     if (mynum==0) write(6,*)' calculate medium prefactor...'
     call calc_medium_factor
  endif

! Causality: Set beginning of time to just before straight P ray
  if (mynum==0) write(6,*)' determine stuff regarding time series/windows....'
  call calc_beg_end_time(thr)

  if (mynum==0) write(6,*)'loading synthetic seismograms........' 
  call read_disp_velocity_seis(thr,phr)

  allocate(timewin(ndumps),misfit_fct(ndumps,num_misfits),chi(nwin,num_misfits))

! for now just taking the component of the singular backward field source....
!  if (mynum==0) write(6,*)'vector-to-scalar reduction for synthetics........' 
 ! call vec2scal_seismogram(u0,v0,uscal,uscal)

   if (mynum==0) write(6,*)"loading data seismograms......"
  call load_data

! for now just taking the component of the singular backward field source....
!  if (mynum==0) write(6,*)'vector-to-scalar reduction for data........' 
!  call vec2scal_seismogram

! needed to add this to allocate post_per... NOT OPTIMAL! 
! ultimately the posteriori input should be an option for the full workflow too...
  call read_posteriori

! shift synthetic time series to zero
  call fftf1d(u0(:,reccomp),uspec)
  call time_shift_vector(uspec,shift_fact_fwd)
  call fftb1d(uspec,u0(:,reccomp))
  call fftf1d(v0(:,reccomp),vspec)
  call time_shift_vector(vspec,shift_fact_fwd)
  call fftb1d(vspec,v0(:,reccomp))

  if (mynum==nproc-1) then
     write(6,*)'Receiver component::',reccomp
     write(6,*)"save unfiltered time series..."
     open(unit=100,file='Data/usyn_unfilt.dat')
     open(unit=101,file='Data/vsyn_unfilt.dat')
     open(unit=102,file='Data/udat_unfilt.dat')
     do i=1,ndumps
        write(100,*)time(i),u0(i,reccomp)
        write(101,*)time(i),v0(i,reccomp)
        write(102,*)time(i),displ_data(i)
     enddo
     close(100);close(101);close(102)
  endif

  if (do_filt .and. filter_what==0) then ! this is not optimal!! incorporate multiple filters into general workflow
     allocate(usyn_filt(ndumps,num_filt),vsyn_filt(ndumps,num_filt),udat_filt(ndumps,num_filt))
     call fftf1d(displ_data,dspec)
     if (mynum==0) write(6,*)"filter synthetics & data multiple times",num_filt
     do ifilt=1,num_filt
        call define_filter(filter_type_posteriori(ifilt),post_per(ifilt,1:2))
        uspec1d = uspec
        call filter_vector(uspec1d)
        call fftb1d(uspec1d,usyn_filt(:,ifilt))
        uspec1d = vspec
        call filter_vector(uspec1d)
        call fftb1d(uspec1d,vsyn_filt(:,ifilt))
        uspec1d = dspec
        call filter_vector(uspec1d)
        call fftb1d(uspec1d,udat_filt(:,ifilt))
        if (mynum==nproc-1) then
           call define_io_appendix(appifilt1,post_per(ifilt,1))
           call define_io_appendix(appifilt2,post_per(ifilt,2))
           open(unit=100,file='Data/usyn_filt_'//appifilt1//'_'//appifilt2//'.dat')
           open(unit=101,file='Data/vsyn_filt_'//appifilt1//'_'//appifilt2//'.dat')
           open(unit=102,file='Data/udat_filt_'//appifilt1//'_'//appifilt2//'.dat')
           do i=1,ndumps
              write(100,*)time(i),usyn_filt(i,ifilt)
              write(101,*)time(i),vsyn_filt(i,ifilt)
              write(102,*)time(i),udat_filt(i,ifilt)
           enddo
           close(100);close(101);close(102)
        endif
     enddo
     u0(:,reccomp) = usyn_filt(:,1)
     v0(:,reccomp) = vsyn_filt(:,1)
     displ_data = udat_filt(:,1)
  elseif (filter_what>0) then
     allocate(usyn_filt(ndumps,1),vsyn_filt(ndumps,1),udat_filt(ndumps,1))
     if (mynum==0) write(6,*)"filter synthetics & data once "
     call fftf1d(displ_data,dspec)
     call define_filter(filter_type)
     uspec1d = uspec
     call filter_vector(uspec1d)
     call fftb1d(uspec1d,usyn_filt(:,1))
     uspec1d = vspec
     call filter_vector(uspec1d)
     call fftb1d(uspec1d,vsyn_filt(:,1))
     uspec1d = dspec
     call filter_vector(uspec1d)
     call fftb1d(uspec1d,udat_filt(:,1))
     if (mynum==nproc-1) then
        write(65,*)'period filter hi low:',filter_period_hi,filter_period_low
        call define_io_appendix(appifilt1,filter_period_low)
        call define_io_appendix(appifilt2,filter_period_hi)
        write(65,*)'period filter hi low:',appifilt1,appifilt2

        open(unit=100,file='Data/usyn_filt_'//appifilt1//'_'//appifilt2//'.dat')
        open(unit=101,file='Data/vsyn_filt_'//appifilt1//'_'//appifilt2//'.dat')
        open(unit=102,file='Data/udat_filt_'//appifilt1//'_'//appifilt2//'.dat')
        do i=1,ndumps
           write(100,*)time(i),usyn_filt(i,1)
           write(101,*)time(i),vsyn_filt(i,1)
           write(102,*)time(i),udat_filt(i,1)
        enddo
        close(100);close(101);close(102)
     endif
     u0(:,reccomp) = usyn_filt(:,1)
     v0(:,reccomp) = vsyn_filt(:,1)
     displ_data = udat_filt(:,1)
  endif

! save time series used in kernels, along with their windowed sistas
  if (mynum==nproc-1) then
     call define_io_appendix(appifilt1,filter_period_low)
     call define_io_appendix(appifilt2,filter_period_hi)
     open(unit=100,file='Data/usyn_final.dat')
     open(unit=101,file='Data/vsyn_final.dat')
     open(unit=102,file='Data/udat_final.dat')
     write(6,*)'NWIN, NWT:',nwin,ntw
     do iwin=1,ntw
        open(unit=1000+iwin,file='Data/usyn_window_'//trim(phase_name(iwin))//'.dat')
        open(unit=2000+iwin,file='Data/vsyn_window_'//trim(phase_name(iwin))//'.dat')
        open(unit=3000+iwin,file='Data/data_window_'//trim(phase_name(iwin))//'.dat')
     enddo
     do iwin = 1, nwin ! number of time windows. 1 if waveform kernels are saved.
        do it=begdumps(iwin),enddumps(iwin)
           write(100,*)time(it),u0(it,reccomp)
           write(101,*)time(it),v0(it,reccomp)
           write(102,*)time(it),displ_data(it)
           if (inside_win(it,iwin)) then
              write(1000+ind_win(it,iwin),*)time(it),u0(it,reccomp)
              write(2000+ind_win(it,iwin),*)time(it),v0(it,reccomp)
              write(3000+ind_win(it,iwin),*)time(it),displ_data(it)
           endif
        enddo
     enddo
     close(100);close(101);close(102)
     do i=1,ntw
        close(1000+i); close(2000+i); close(3000+i)
     enddo
     write(6,*)'min max displacement time series:',minval(u0(:,reccomp)),maxval(u0(:,reccomp))
     write(6,*)'min max velocity time series:',minval(v0(:,reccomp)),maxval(v0(:,reccomp))
  endif

  if (mynum==nproc-1) call time_dependent_seismogram_vtk(u0(:,reccomp),ndumps,thr,phr)

! end time series analysis ----------------------------------------------------

  if (mynum==0) write(6,*)'define shape of time window........' 
  call def_time_window_shape

  allocate(denomu(ntw),denomv(ntw))   
  do i = 1,ntw
     call compute_denominator(t1(i),t2(i),time,timestep,denomu(i),denomv(i))
  enddo

end subroutine init_time_window_misfit_medium
!----------------------------------------------------------------------------------------------------------

subroutine time_dependent_seismogram_vtk(u0,n,thr,phr)

use input_output, only : write_VTK_bin_scal_mesh3d,write_VTK_bin_scal_topology
use coord_trafo, only : rthetaphi2xyz
use data_arrays, only : kernmesh_type

integer, intent(in) :: n
real,dimension(n), intent(in) :: u0
real(kind=realkind), intent(in)  :: thr,phr  
real :: x(n),y(n),z(n),u1(n),maxu,dx,dz,maxx,x0,y0,z0,maxz,r
real :: xc(n*4),yc(n*4),zc(4*n),uc(4*n) 
character(len=100) :: filename1
character(len=200) :: filename
integer :: iwin,i
character(len=4) :: appidump 

r=6371000.
if (trim(kernmesh_type)=='3dslides') then 
   call rthetaphi2xyz(x0,y0,z0,r,thr,phr)
else
   x0=0.
   y0=0.
endif
z0=6371000.*1.1

maxu=maxval(abs(u0))
maxx=.75*6371000.
dx=(maxx-x0)/(real(n))
maxz=maxx/5.
do i=1,n
   if (trim(kernmesh_type)=='3dslides') then 
      x(i) = x0-real(i)*dx
   else
      x(i) = x0+real(i)*dx
   endif
   y(i) = 0.
   z(i) = z0 + u0(i)/maxu*maxz
enddo
u1=-1.e20
do iwin = 1, nwin ! number of time windows. 1 if waveform kernels are saved.
   do i=begdumps(iwin),enddumps(iwin)
      u1(i) = 1.e20
      call define_io_appendix(appidump,i)
      filename1='Data/seis_'//appidump
      call write_VTK_bin_scal_mesh3d(x,y,z,u1,n,1,filename1)
   enddo
enddo

!!! CELL TOPOLOGY
do i=1,n-1
   xc(4*(i-1) + 1)= x(i)-2.*dx
   xc(4*(i-1) + 2)= x(i+1)-2.*dx
   xc(4*(i-1) + 3)= x(i+1)+2.*dx
   xc(4*(i-1) + 4)= x(i)+2.*dx
   yc(4*(i-1) + 1)= y(i)
   yc(4*(i-1) + 2)= y(i+1)
   yc(4*(i-1) + 3)= y(i+1)
   yc(4*(i-1) + 4)= y(i)
   zc(4*(i-1) + 1)= z(i)+2.*dx
   zc(4*(i-1) + 2)= z(i+1)+2.*dx
   zc(4*(i-1) + 3)= z(i+1)-2.*dx
   zc(4*(i-1) + 4)= z(i)-2.*dx
enddo
xc(4*(n-1)+1:4*n)=x(n)
yc(4*(n-1)+1:4*n)=y(n)
zc(4*(n-1)+1:4*n)=z(n)
uc=-1.e20
do iwin = 1, nwin ! number of time windows. 1 if waveform kernels are saved.
   do i=begdumps(iwin),enddumps(iwin)
      if (i>1) then 
         uc(begdumps(1):4*(i-1) + 1) = 1.e20
         uc(begdumps(1):4*(i-1) + 4) = 1.e20
      endif
      call define_io_appendix(appidump,i)
      filename='Data/seis_cell_'//appidump
      call write_VTK_bin_scal_topology(xc,yc,zc,uc,n,filename)
   enddo
enddo

write(6,*)'Saved seismograms as evolving series in vtk to plot alongside waveform kernels.'
write(6,*)'load tem into paraview: Data/seis_*.vtk'

end subroutine time_dependent_seismogram_vtk

!----------------------------------------------------------------------------------------------------------
subroutine read_posteriori
use data_posteriori

use data_fft, only : filter_type,ntimes,nomega,uphys1d,uspec1d
use parameters, only: define_sem_kernel_meshes
use data_misfits

include 'mesh_params.h'

integer :: i,ii,j,iwin,iidom,npad
character(len=100) :: junk
real(kind=realkind) :: src_depthbwd
real(kind=realkind) :: thr,phr
character(len=4) :: model_reconsttmp(2),model_runtmp(20)
character(len=5) :: misfittmp(20)
logical :: travelpicktmp(20)

open(unit=99,file='input_posteriori.dat')
read(99,*)junk
write(6,*)'Misfits to compute:'
read(99,*)xc_tt ; if (xc_tt) write(6,*)'  ... cross-correlation traveltime'
read(99,*)xc_am; if (xc_am) write(6,*)'  ... cross-correlation amplitude'
read(99,*)pk_tt ; if (pk_tt) write(6,*)'  ... peak-arrival traveltime'
read(99,*)on_tt ; if (on_tt) write(6,*)'  ... onset traveltime'
read(99,*)insta ; if (insta) write(6,*)'  ... instantaneous waveform pick'
read(99,*)wf_dd; if (wf_dd) write(6,*)'  ... waveform misfit (data needed)'
read(99,*)instph; if (instph) write(6,*)'  ... instantaneous phase misfit (data needed)'
read(99,*)insten ; if (insten) write(6,*)'  ... instanteneous envelope misfit (data needed)'
read(99,*)phase; if (phase) write(6,*)'  ... phase misfit (data needed)'
read(99,*)envel ; if (envel) write(6,*)'  ... envelope misfit (data needed)'

read(99,*)junk

read(99,*)num_filt
write(6,*)'Number of filters:',num_filt
need_filt=0
if (num_filt/=0) then
   allocate(post_per(num_filt,2),filter_type_posteriori(num_filt),app_period1(num_filt),app_period2(num_filt),filter_char(num_filt))
   do i=1,num_filt
      read(99,*)filter_type_posteriori(i),post_per(i,1),post_per(i,2)
      write(6,*)'filter type and min/max period: ',filter_type_posteriori(i),post_per(i,1),post_per(i,2)
      if (filter_type_posteriori(i)/='non') need_filt=need_filt+1
      write(app_period1,'(I3)') nint(post_per(i,1))
      write(app_period2,'(I3)') nint(post_per(i,2))
      !call define_io_appendix3_real(app_period1(i),post_per(i,1))
      !call define_io_appendix3_real(app_period2(i),post_per(i,2))
      filter_char(i) = trim(filter_type_posteriori(i))//'_'//app_period1(i)//'_'//app_period2(i)
      write(6,*)'filter char:',filter_char(i)
   enddo
   do_filt=.true.
else
   num_filt=1  ! just for convenience in the loop in main.f90
   do_filt=.false.
endif

read(99,*)junk
write(6,*)'medium parameterization:'
read(99,*)do_lam_post ; if (do_lam_post) write(6,*)' ... lambda'
read(99,*)do_mu_post; if (do_mu_post) write(6,*)' ... mu'
read(99,*)do_rho_post; if (do_rho_post) write(6,*)' ...rho'
read(99,*)do_vp_post; if (do_vp_post) write(6,*)' ... vp'
read(99,*)do_vs_post; if (do_vs_post) write(6,*)' ... vs'
read(99,*)do_imped_post; if (do_imped_post) write(6,*)' ... imped'
read(99,*)model_norm; if (model_norm) write(6,*)' ... model normalization'

read(99,*)junk

read(99,*)num_dir
write(6,*)'Number of Mij-Pn couples (input waveform kernels):',num_dir
allocate(run_dir(num_dir))
do i=1,num_dir
   read(99,*)run_dir(i)
   write(6,*)run_dir(i)
enddo
read(99,*)(Mij_post(i),i=1,6)
write(6,*)'Moment tensor:'
write(6,*)(Mij_post(i),i=1,6)

read(99,*)junk
read(99,*)consider_data; write(6,*)'Consider data?',consider_data
read(99,*)data_dir; write(6,*)'Data directory:',trim(data_dir)
read(99,*)num_tot_picks; write(6,*)'number of time picks:',num_tot_picks
allocate(time_pick_read(1:num_tot_picks),pick_misfit(1:num_tot_picks))
time_pick_read=0.
do j=1,num_tot_picks
   read(99,*)time_pick_read(j),pick_misfit(j)
   write(6,*)'time pick and type:',time_pick_read(j),pick_misfit(j)
enddo
close(99)

end subroutine read_posteriori
!---------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------
subroutine init_wavefields
  use data_arrays

 if (do_lam) allocate(urec(1:npts),usrc(1:npts))
 if (do_rho) allocate(vrec(1:npts,1:dim_bwd),vsrc(1:npts,1:dim_fwd))

 end subroutine init_wavefields
!---------------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------------
subroutine calc_medium_factor

  use data_arrays, only : do_rho,do_lam,do_mu,do_vp,do_vs

  use coord_trafo, ONLY: get_r_theta
  use background_models, ONLY : velocity

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  real(kind=realkind) :: r,theta,phi,vp,vs,rhotmp
  integer :: ipt,idom

  if (do_rho) allocate(rho(1:npts))
  if (do_lam) allocate(lam(1:npts))
  if (do_mu .or. do_vp) allocate(mu(1:npts))
  if (do_vp) allocate(vp_fac(1:npts))
  if (do_vs) allocate(vs_fac(1:npts))

!  deallocate(idom_ipt)

! read mesh as saved from input_output.f90
  if (calc_or_load_kernels=='load_wfkernels') then 
       open(unit=50,file=trim(run_dir(1))//'/Data/'//'kernelmesh_rthph_'//appmynum//'.dat',status='old')
   else
      open(unit=50,file='Data/'//'kernelmesh_rthph_'//appmynum//'.dat',status='old')
  endif 
  do ipt=1,npts
    ! check for each element which domain it belongs to
     read(50,*)r,theta,phi,idom
     vp=velocity(dble(r),'v_p',idom,bkgrdmodel,lfbkgrdmodel)
     vs=velocity(dble(r),'v_s',idom,bkgrdmodel,lfbkgrdmodel)
     rhotmp=velocity(dble(r),'rho',idom,bkgrdmodel,lfbkgrdmodel)

     if (do_rho) rho(ipt)=rhotmp ! mass density
     if (do_lam) lam(ipt) = & ! bulk mod, NOT Lame!
          rhotmp*(vp*vp-two*vs*vs)+2./3.*rhotmp*vs*vs
     if (do_mu) mu(ipt) = rhotmp*vs*vs ! Shear modulus
  enddo 
  close(50)
  if (do_lam) write(6,*)'max,min lam:',maxval(lam),minval(lam)

  if (do_vp) vp_fac = 2.d0+8.d0/3.d0*mu/lam 
  if (do_vs) vs_fac = -8.d0/3.d0*mu/lam     
  if (do_vp .and. .not. do_mu) deallocate(mu)

10 format(i9,1pe11.3,i3,2(1pe11.3))
  
end subroutine calc_medium_factor
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
subroutine calc_beg_end_time(thr)
                             
  use coord_trafo, ONLY: get_r_theta
  use background_models, ONLY : velocity
  use mpi

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  real(kind=realkind), intent(in) :: thr
  integer                :: begdumps1, i, it, iel2
  real(kind=realkind)    :: r1, th, ph
  real(kind=realkind)    :: straight_dist, straight_time
  integer                :: ipt, idom, nrad
  real(kind=realkind)    :: src_shift
  character(len=20)      :: stf_type

  ! get vpmax from forward run
  ! Is now all in 'mesh_params_kernel.h'

  !vpmax = 0.0
  !rmax = 0.0

  !open(unit=50,file='backgroundmodel_kmscale.dat0000',status='old')

  !read(50,*) nrad
  !do ipt=1,nrad
  !  read(50,*) r, vp, vs, rho
  !  if (r > rmax) rmax = r
  !  if (vp > vpmax) vpmax = vp
  !   enddo
  !close(50)
  
  !use source radius to have a stict lower bound:
  straight_dist = dsqrt(2.d0) * (rmax - src_depth * 1e3) * dsqrt(1.d0-dcos(dble(thr)))
  straight_time = straight_dist/vpmax
  begdumps1 = floor(straight_time/((dble(nt)*deltat)/dble(ndumps)))-1

  if (mynum==0) then 
     write(6,*)''
     write(6,*)'Max. radius [km]:',real(rmax/1000.)
     write(6,*)'Max. radius - src_depth [km]:',real(rmax/1000. - src_depth)
     write(6,*)'Max. velocity [km/s]:',real(vpmax/1000.)
     write(6,*)'epicentral distance [deg]:',real(thr*180./pi)
     write(6,*)'min. causal distance [km]:  ',real(straight_dist/1000.)
     write(6,*)'min. causal time [s]:  ',real(straight_time)
  endif

  save_wfkern=.true.

! define arrays for identifying which time window (if any) a given time is in
  select case(char_wfkern)

! saving the waveform kernel at all times (specified by time window in input)
  case('allt')
     nwin = 1
     if (calc_or_load_kernels=='calc_wfkernels') allocate(begdumps(nwin),enddumps(nwin))
     if (twf1<0.0) then  ! from the first causal frame
        begdumps(1) = begdumps1
     else ! from the specified time
        begdumps(1) = minloc(abs(time-twf1),1)
     endif

     if (twf2<0.0) then  ! until the last frame
        enddumps(1) = ndumps
     else ! until the specified last time 
        enddumps(1) = min(minloc(abs(time-twf2),1),ndumps)
     endif

!    define arrays for identifying which time window (if any) a given time is in
     allocate(ind_win(ndumps,1),inside_win(ndumps,1),save_n_erase(ndumps,1))
     ind_win = 1
     inside_win = .false.
     save_n_erase = .false.

     do it=begdumps(1),enddumps(1)
        do ipt = 1,ntw
           if ( time(it)>=t1(ipt) .and. time(it)<=t2(ipt) ) then 
              ind_win(it,1)=ipt
              inside_win(it,1)=.true.
              if (mynum==0) then 
              endif
              write(6,*)'t1,time,t2:',t1(ipt),time(it),t2(ipt)
              if  ( t2(ipt)-time(it)<time(2)-time(1) .and. time(it)<=t2(ipt) ) save_n_erase(it,1) = .true.
           endif
           
           if ( ipt>1 ) then 
              if ( t1(ipt)<t2(ipt-1) ) then 
              if (mynum==0) then
                 write(6,*)' WE DO NOT ALLOW FOR OVERLAPPING TIME WINDOWS IF',&
                        '      ALL WAVEFORM KERNELS ARE TO BE SAVED!'
                 write(6,*)' Either change waveform kernel settings to only'
                 write(6,*)'   saving them inside the time window, or define'
                 write(6,*)'   non-overlapping time windows'
                 write(6,*)' .... changing it such that they are exactly non-overlapping:'
                 write(6,*)'old time window:',ipt,t1(ipt),t2(ipt)
              endif
              t1(ipt)=t2(ipt-1)
              if (mynum==0) write(6,*)'new time window,phase:',ipt,t1(ipt),t2(ipt),trim(phase_name(ipt))
!                 stop
              endif
           endif

        enddo
     enddo

! only saving waveform and static kernels in windows [t1,t2]
  case('wint')
     nwin = ntw
 if (calc_or_load_kernels=='calc_wfkernels')    allocate(begdumps(nwin),enddumps(nwin))
     do i=1,nwin
        begdumps(i) = minloc(abs(time-t1(i)),1)
        enddumps(i) = min(minloc(abs(time-t2(i)),1),ndumps)
        if (mynum==0) then 
           write(6,*)'wint time window, t1,t2:',i,real(t1(i)),real(t2(i))
           write(6,*)'  wint beg/end dumps:',begdumps(i),enddumps(i)
           write(6,*)'  wint beg/end dumps t:',real(time(begdumps(i))), &
                                               real(time(enddumps(i)))
        endif
     enddo

     allocate(ind_win(ndumps,nwin),inside_win(ndumps,nwin))
     allocate(save_n_erase(ndumps,nwin))
     ind_win = 0; inside_win = .false.; save_n_erase = .false.
     do i=1,nwin
        do it=begdumps(i),enddumps(i)
!           if ( time(it)>=t1(i) .and. time(it)<=t2(i) ) then 
            if ( time(it)>=time(begdumps(i)) .and. time(it)<=time(enddumps(i)) ) then
              ind_win(it,i)=i
              inside_win(it,i)=.true.
              if (it==enddumps(i)) save_n_erase(it,i) = .true.
           else
              write(6,*)' Something is so wrong here.'
              write(6,*)'  time window #,t1,t2:',i,t1(i),t2(i)
              write(6,*)'  time:',time(it)
              stop
           endif
        enddo
     enddo

! save no waveform kernels at all
  case('none')
     save_wfkern=.false.
     nwin = ntw
 if (calc_or_load_kernels=='calc_wfkernels')    allocate(begdumps(nwin),enddumps(nwin))
     do i=1,nwin
        begdumps(i) = minloc(abs(time-t1(i)),1)
        enddumps(i) = min(minloc(abs(time-t2(i)),1),ndumps)
     enddo

     allocate(ind_win(ndumps,nwin),inside_win(ndumps,nwin))
     allocate(save_n_erase(ndumps,nwin))
     ind_win = 0; inside_win = .false.; save_n_erase = .false.
     do i=1,nwin
        do it=begdumps(i),enddumps(i)
           if ( time(it)>=t1(i) .and. time(it)<=t2(i) ) then 
              ind_win(it,i)=i
              inside_win(it,i)=.true.
              if (it==enddumps(i)) save_n_erase(it,i) = .true.
           else
              write(6,*)' Something is so wrong here.'; stop
           endif
           
        enddo
     enddo
 
  case default
     write(6,*)'WRONG definition of what to do with waveform kernels!', &
               char_wfkern
     stop

  end select     

  dt=time(2)-time(1)
  if (mynum==0) then 
     write(6,*)' number of separate windows:',nwin
     do i=1,nwin
        write(6,*)i,'time window requested:',real(t1(i)),real(t2(i))
        write(6,*)i,'time window used     :',real(time(begdumps(i))), &
                                             real(time(enddumps(i)))
        write(6,*)i,'window index         :',begdumps(i),enddumps(i)
     enddo
     write(6,*)''
     write(6,*)'time step:', dt
  endif

10 format(i9,1pe11.3,i3,2(1pe11.3))

! source shift for the backward field
! SHIFTING THE TIME OF THE BACKWARD WAVEFIELD TO MAKE SURE SOURCE 
! RUPTURES AT TIME ZERO
  if (mynum==0) &
       write(6,*)'  SHIFTING BACKWARD TIME to put "source" at time zero!'

  open(unit=20000,file='sourceparams_bwd.dat')
     read(20000,*)realjunk
     read(20000,*)charjunk
     read(20000,*)charjunk
     read(20000,*)
     read(20000,*)
     read(20000,*)realjunk
     read(20000,*)realjunk
     read(20000,*)realjunk
     read(20000,*)stf_type
  close(20000)
  if (stf_type=='dirac_0' .or. stf_type=='quheavi') then 
     src_shift=0.
  elseif (stf_type=='heavis') then 
     src_shift=300.
  else
     src_shift = 1.5*strain_samp*dt ! 1.5 times the dominant source period
  endif
  write(6,*)mynum,'source shift:',src_shift,dt

  ibwd_shift = floor(src_shift/dt)
  if (mynum==0) then
     write(6,*) 'dt,source shift,real(shift),int(shift)'
     write(6,*)real(time(2)-time(1)),real(src_shift), &
               real(src_shift/(time(2)-time(1))),ibwd_shift
     write(6,*)"1.5*strain_samp (should be equal to src shift):", &
               floor(1.5*strain_samp)
  endif
  if (stf_type .ne. 'dirac_0' .and. stf_type .ne. 'heavis' .and. floor(1.5*strain_samp)/=ibwd_shift) then 
     write(6,*)'ERROR: time shift not consistent between strain sampling'
     write(6,*)'       and calculation based on time series:'
     write(6,*)'       strain sampling units:',floor(1.5*strain_samp)
     write(6,*)'       ibwd_shift:',ibwd_shift
     stop
  endif

  num_it=0
  do i=1,nwin
     do it=begdumps(i),enddumps(i)
        num_it = num_it+1
     enddo
  enddo

  if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
  do i=0,nproc-1
     if (nproc>1) CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
     if (mynum==i) &
          write(6,*)mynum,' Total number of kernel calculations:',num_it
  enddo

end subroutine calc_beg_end_time
!------------------------------------------------------------------------

!------------------------------------------------------------------------
subroutine def_time_window_shape
use data_misfits 
integer :: it,iwin

! define boxcar time window (for now...)
timewin=0.
do iwin = 1, nwin ! number of time windows. 1 if all waveform kernels are saved.
   do it=begdumps(iwin),enddumps(iwin) 
      if (inside_win(it,iwin)) timewin(it)=1.
   enddo
enddo
open(unit=90,file='Data/timewindow.dat')
do it=1,ndumps
   write(90,*)time(it),timewin(it)
enddo
close(90)

end subroutine def_time_window_shape
!------------------------------------------------------------------------

!------------------------------------------------------------------------
subroutine vec2scal_seismogram(uvec,vvec,uscal,vscal)

use data_misfits
use data_posteriori, only : reccomp_post,src_type2bwd_post  ! should not be in here....

! define scalar seismogram to be used in the misfits
real,dimension(1:ndumps,3), intent(in) :: uvec,vvec
real,dimension(1:ndumps), intent(out) :: uscal,vscal
integer :: i

uscal=uvec(:,reccomp_post(1))
vscal=vvec(:,reccomp_post(1))
do i=2,num_dir
   if (src_type2bwd_post(i)/=src_type2bwd_post(i-1)) then
      ! stack components
      if (i==2) then
      vscal=vscal + vvec(:,reccomp_post(i))
      uscal=uscal + uvec(:,reccomp_post(i))
      else
        if (src_type2bwd_post(i)/=src_type2bwd_post(i-2)) then 
           vscal=vscal + vvec(:,reccomp_post(i))
           uscal=uscal + uvec(:,reccomp_post(i))
     endif
  endif
  endif
enddo

end subroutine vec2scal_seismogram
!------------------------------------------------------------------------

!------------------------------------------------------------------------
subroutine load_data

use data_misfits

integer :: it

! load data: data seismogram needs to be in a file with 1st column=time, 2nd column=displacement component
if (consider_data) then 
write(6,*)'loading data from',trim(data_dir)
open(unit=98,file=trim(data_dir))
open(unit=99,file='Data/data_displ.dat')
do it=1,ndumps
   read(98,*)time_data(it),displ_data(it)
   write(99,*)time_data(it),displ_data(it)
enddo
close(99);close(98)

if (abs(time_data(2)-time_data(1) - dt) >smallval) then 
   write(6,*)'Problem with data time sampling: not equal to synthetics:'
   write(6,*)'data:',time_data(2)-time_data(1)
   write(6,*)'synthetics:',dt
   stop
endif
else
   time_data=time
   displ_data=0.0
endif

end subroutine load_data
!------------------------------------------------------------------------

!-------------------------------------------------------------------------
  subroutine compute_denominator(t1,t2,time,timestep,denomu1,denomv1)
    
  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

    real(kind=realkind), intent(in)  :: t1,t2
    real(kind=realkind), intent(out) :: denomu1,denomv1
    real(kind=realkind), dimension(ndumps) :: time,timestep
    integer :: it

    denomu1=zero; denomv1=zero

    do it=1,ndumps
       if (time(it)>=t1 .and. time(it)<=t2) then
          if (src_type=='monopole') then 
             denomv1=denomv1+v0(it,1)*v0(it,1)+v0(it,2)*v0(it,2)
             denomu1=denomu1+u0(it,1)*u0(it,1)+u0(it,2)*u0(it,2)
          else
            denomv1=denomv1+v0(it,1)*v0(it,1)+v0(it,2)*v0(it,2)+v0(it,3)*v0(it,3)
            denomu1=denomu1+u0(it,1)*u0(it,1)+u0(it,2)*u0(it,2)+u0(it,3)*u0(it,3)
          endif
       endif
    enddo

! If zero-trace time window was picked, at least don't normalize with infty
    if (denomv1==0.d0) denomv1 = 1.d0
    if (denomu1==0.d0) denomu1 = 1.d0

    denomv1=1.d0/denomv1
    denomu1=1.d0/denomu1

  end subroutine compute_denominator
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
subroutine read_disp_velocity_seis(thr,phr)

#ifdef unc
   use netcdf
#endif

   use data_mesh, only: use_netcdf, ncid_fw_in, ncid_fw_snapshot
   use mpi
   implicit none
   include 'mesh_params.h'
   include 'mesh_params_kernel.h'

   real(kind=realkind), intent(in)             :: thr,phr
   real(kind=realkind), allocatable, dimension(:) :: surftheta
   integer :: it,recel,nsurf,iel
   real(kind=realkind) :: minth,buff,buff2,prefac1,prefac2
   character(len=4) :: appielem
   logical :: have_rec
   real(kind=realkind), allocatable :: utmp(:,:)
   integer    :: ncvarid_surfdisp, ncvarid_surfvelo, status
! find closest receiver location for thr
   minth=361.
   
   open(unit=666,file=trim(dir_fwdmesh1(1))//'/Data'//'/surfelem_coords.dat')
   read(666,*)nsurf
   allocate(surftheta(nsurf))
   do iel=1,nsurf
      read(666,*)surftheta(iel)
      if (abs(surftheta(iel)-thr*180./pi)<minth) then 
         recel=iel
         minth=abs(surftheta(iel)-thr*180./pi)
      endif
   enddo
   close(666)  
   if (mynum==0) then
     write(6,*)
     write(6,*)'Desired receiver location:',thr*180./pi
     write(6,*)'Closest location found   :',surftheta(recel)
     write(6,*)'Receiver number:',recel
     write(6,*)
  endif
    if (src_type2 == 'mxx_p_myy' .or. src_type2=='mzz' &
        .or. src_type2=='explosion' .or. src_type2=='vertforce')  then
       prefac1 = 1.d0
    elseif (src_type2=='mxz' .or. src_type2=='xforce') then 
       prefac1 =  dcos(dble(phr))
       prefac2 =  -dsin(dble(phr))
    elseif (src_type2=='myz' .or. src_type2=='yforce') then 
       prefac1 = dsin(dble(phr))
       prefac2 = dcos(dble(phr))
    elseif (src_type2=='mxx_m_myy') then 
       prefac1 = dcos(2.*dble(phr))
       prefac2 = -dsin(2.*dble(phr))
    elseif (src_type2=='mxy') then 
       prefac1 = dsin(2.*dble(phr))
       prefac2 = dcos(2.*dble(phr))
    else
       write(6,*)'radiation type',src_type2,'unknown!'
       stop
    endif

     call define_io_appendix(appielem,recel)
     
#ifdef unc
     if (use_netcdf) then
         call check( nf90_inq_ncid(ncid=ncid_fw_in, name="Surface", &
                                   grp_ncid=ncid_fw_snapshot) )
         write(6,*) 'Netcdf group Surface has NCID ', ncid_fw_snapshot
         call check( nf90_inq_varid(ncid=ncid_fw_snapshot, name="displacement", &
                                    varid=ncvarid_surfdisp) )
         call check( nf90_inq_varid(ncid=ncid_fw_snapshot, name="velocity", &
                                    varid=ncvarid_surfvelo) )
         if (src_type=='monopole') then
             call check( nf90_get_var(ncid=ncid_fw_snapshot, varid=ncvarid_surfdisp, &
                                      start = (/1,1,recel/), count=(/ndumps,2,1/), &
                                      values = u0 )) 
             call check( nf90_get_var(ncid=ncid_fw_snapshot, varid=ncvarid_surfvelo, &
                                      start = (/1,1,recel/), count=(/ndumps,2,1/), &
                                      values = v0 )) 
         else
             call check( nf90_get_var(ncid=ncid_fw_in, varid=ncvarid_surfdisp, &
                                      start = (/1,1,recel/), count=(/ndumps,3,1/), &
                                      values = u0 )) 
             call check( nf90_get_var(ncid=ncid_fw_in, varid=ncvarid_surfvelo, &
                                      start = (/1,1,recel/), count=(/ndumps,3,1/), &
                                      values = v0 )) 
         end if
     end if
#endif 

     if (.not. use_netcdf) then
         open(unit=50,file=trim(dir_fwdmesh1(1))//&
              '/Data'//'/surfelem_disp.dat'//appielem)
         open(unit=51,file=trim(dir_fwdmesh1(1))//&
              '/Data'//'/surfelem_velo.dat'//appielem)

         do it=1,ndumps
            if (src_type=='monopole') then
               read(50,*)u0(it,1:2)
               read(51,*)v0(it,1:2)
            else 
               read(50,*)u0(it,1:3)
               read(51,*)v0(it,1:3)
            endif
         end do
         close(50)
         close(51)
     end if
     
     if (mynum==nproc-1) then
        open(unit=61,file='Data/rec_disp_s.dat',POSITION="REWIND")
        open(unit=62,file='Data/rec_disp_z.dat',POSITION="REWIND")     
        if (src_type/='monopole') &
             open(unit=63,file='Data/rec_disp_ph.dat',status='new')
        open(unit=64,file='Data/rec_dispr.dat',POSITION="REWIND")
        open(unit=65,file='Data/rec_dispth.dat',POSITION="REWIND")     
        
        open(unit=71,file='Data/rec_velo_s.dat',POSITION="REWIND")
        open(unit=72,file='Data/rec_velo_z.dat',POSITION="REWIND") 
        if (src_type/='monopole') &
             open(unit=73,file='Data/rec_velo_ph.dat',POSITION="REWIND")
        open(unit=74,file='Data/rec_velo_r.dat',POSITION="REWIND")
        open(unit=75,file='Data/rec_velo_th.dat',POSITION="REWIND") 
        endif

! Note that this dumps at the desired, not actual theta, phi!
    do it=1,ndumps
        !if (src_type=='monopole') then
        !   read(50,*)u0(it,1:2)
        !   read(51,*)v0(it,1:2)
        !else 
        !   read(50,*)u0(it,1:3)
        !   read(51,*)v0(it,1:3)
        !endif

        if (mynum==nproc-1) then       
           write(61,*)time(it),u0(it,1)*prefac1
           write(62,*)time(it),u0(it,2)*prefac1
           if (src_type/='monopole') &
                write(63,*)time(it),u0(it,3)*prefac2
           
           write(64,*)time(it),prefac1*(sin(thr)*u0(it,1) + cos(thr)*u0(it,2))
           write(65,*)time(it),prefac1*(cos(thr)*u0(it,1) - sin(thr)*u0(it,2))
           
           write(71,*)time(it),v0(it,1)*prefac1
           write(72,*)time(it),v0(it,2)*prefac1
           if (src_type/='monopole') &
                write(73,*)time(it),v0(it,3)*prefac2
           
           write(74,*)time(it),prefac1*(sin(thr)*v0(it,1) + cos(thr)*v0(it,2))
           write(75,*)time(it),prefac1*(cos(thr)*v0(it,1) - sin(thr)*v0(it,2))
        endif
     enddo

     close(61); close(62); close(63); close(64); close(65)
     close(71); close(72); close(73); close(74); close(75)

! rotate seismograms to spherical system
allocate(utmp(ndumps,3))
utmp=u0
u0(:,1)=prefac1*(cos(thr)*utmp(:,1) - sin(thr)*utmp(:,2))
u0(:,3)=prefac1*(sin(thr)*utmp(:,1) + cos(thr)*utmp(:,2))
deallocate(utmp)

end subroutine read_disp_velocity_seis
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
#ifdef unc
subroutine check(status)
! Translates netcdf error codes into error messages
  use netcdf
  implicit none
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop 1
  end if
end subroutine check  
#endif
!-------------------------------------------------------------------------

!--------------------------------------------------------------
subroutine allocate_arrays

use data_arrays
use data_mesh, only : npts,dim_fwd,dim_bwd
use global_parameters

implicit none
   if (mynum==0)write(6,*)mynum,' allocating arrays........'; call flush(6)

  allocate(kernel(1:npts))

  if (do_lam) allocate(lamttkern(1:npts),lamampkern(1:npts))
  if (do_mu)  allocate(muttkern(1:npts),muampkern(1:npts)) 
  if (do_rho) allocate(rhottkern(1:npts),rhoampkern(1:npts),&
                       vsrc(1:npts,1:dim_fwd),vrec(1:npts,1:dim_bwd))

  if (do_vp) allocate(vpkernel(1:npts),vpttkern(1:npts),vpampkern(1:npts))
  if (do_vs)  allocate(vskernel(1:npts),vsttkern(1:npts),vsampkern(1:npts)) 
  if (do_imped)allocate(impkernel(1:npts),impttkern(1:npts),impampkern(1:npts))

  if (do_lam) lamttkern = zero; if (do_mu) muttkern = zero
  if (do_rho) rhottkern = zero; if (do_vp) vpttkern = zero
  if (do_vs) vsttkern = zero; if (do_imped) impttkern = zero

  if (do_lam) lamampkern = zero; if (do_mu) muampkern = zero
  if (do_rho) rhoampkern = zero; if (do_vp) vpampkern = zero
  if (do_vs) vsampkern = zero; if (do_imped) impampkern = zero

end subroutine allocate_arrays
!--------------------------------------------------------------

!-------------------------------------------------------------------------
subroutine waveform_kernels(appidump1,t_now,it)

  use data_mesh, only : mynum,npts,times_medium
  use data_mesh, only : rho,lam,mu,vp_fac,vs_fac
  use data_mesh, only : dt,save_wfkern,do_azim
  use data_mesh, only : globmax_lam,globmax_rho,globmax_mu
  use data_mesh, only : globmax_imp,globmax_vp,globmax_vs
  use data_arrays
  use input_output, only : save_kernel

  implicit none
  integer :: it
  character(len=4),intent(in) :: appidump1
  real(kind=realkind), intent(in) :: t_now
  character(len=20) :: fname
  
  if (do_vp) vpkernel = vp_fac * lamkern_kxt(it,1:npts)
  if (do_vs) vskernel = 2.d0*mukern_kxt(it,1:npts) + vs_fac*lamkern_kxt(it,1:npts)
  if (do_imped) impkernel = rhokern_kxt(it,1:npts) + mukern_kxt(it,1:npts) + &
                            lamkern_kxt(it,1:npts)

  if (save_wfkern) then
     if (do_vp) then 
        fname='Kxt_vp'
        call save_kernel(npts,vpkernel,fname,appidump1,globmax_vp)
        write(6,124)t_now,mynum,'vp',maxval(abs(vpkernel))
     endif

     if (do_vs) then
        fname='Kxt_vs'
        call save_kernel(npts,vskernel,fname,appidump1,globmax_vs)
        write(6,124)t_now,mynum,'vs',maxval(abs(vskernel))
     endif

     if (do_imped) then 
        fname='Kxt_imp'
        call save_kernel(npts,impkernel,fname,appidump1,globmax_imp)
        write(6,124)t_now,mynum,'imp',maxval(abs(impkernel))
     endif
  endif

124 format('<<Runtime info>> t=',f6.1,'s, proc=',i3,', full ',a4,' kernel:',1pe10.2)

end subroutine waveform_kernels
!-------------------------------------------------------------------------

! -----------------------------------------------------------------------
subroutine misfit_kernels(it,iwin,appidump1)

  use data_mesh, only : time,t1,t2,ind_win,inside_win,mynum,npts,times_medium
  use data_mesh, only : save_n_erase,rho,lam,mu,vp_fac,vs_fac,reccomp
  use data_mesh, only : dt,save_wfkern,denomu,denomv
  use data_arrays
  use global_parameters, only : zero,two,epsi
  use input_output, only : save_kernel
  use mpi

  implicit none
  integer, intent(in) :: it,iwin
  character(len=4),intent(in) :: appidump1
  character(len=4) :: appmywin
  character(len=20) :: fname

124 format('<<Runtime info>> t=',f6.1,'s, proc=',i3,', full ',a4,' kernel:',1pe10.2)
  appmywin=trim(phase_name(ind_win(it,iwin)))

!---------------------------------------------------------------------------
     if (do_lam) then
           lamttkern  = lamttkern  + v0(it,reccomp) * lamkern_kxt(it,1:npts)
           lamampkern = lamampkern + u0(it,reccomp) * lamkern_kxt(it,1:npts)
!          last point window: scale kernels, save
           if (save_n_erase(it,iwin)) then
              lamttkern = dt*lamttkern*denomv(ind_win(it,iwin))
              lamampkern = dt*lamampkern*denomu(ind_win(it,iwin))
              fname='Kx_lam_tt' ! should include filtering etc.....
              call save_kernel(npts,lamttkern,fname,appmywin)
              fname='Kx_lam_am' ! should include filtering etc.....
              call save_kernel(npts,lamampkern,fname,appmywin)

              write(6,*)mynum,'========================================'
              write(6,*)mynum,'Saved lam tt/amp static kernels at t=',time(it)
              write(6,*)mynum,'Saved lam for time window #', ind_win(it,iwin)
              write(6,*)mynum,'Saved lam between:',&
                               t1(ind_win(it,iwin)),t2(ind_win(it,iwin)) 
              write(6,*)mynum,'========================================'
           endif
        endif

!---------------------------------------------------------------------------
     if (do_mu) then
           muttkern  = muttkern  + v0(it,reccomp) * mukern_kxt(it,1:npts)
           muampkern = muampkern + u0(it,reccomp) * mukern_kxt(it,1:npts)
!          last point window: scale kernels, save
           if (save_n_erase(it,iwin)) then
              muttkern = dt*muttkern*denomv(ind_win(it,iwin))
              muampkern = dt*muampkern*denomu(ind_win(it,iwin))
              fname='Kx_mu_tt' ! should include filtering etc.....
              call save_kernel(npts,muttkern,fname,appmywin)
              fname='Kx_mu_am' ! should include filtering etc.....
              call save_kernel(npts,muampkern,fname,appmywin)
              write(6,*)mynum,'========================================'
              write(6,*)mynum,'Saved mu tt/amp static kernels at t=',time(it)
              write(6,*)mynum,'Saved mu for time window #', ind_win(it,iwin)
              write(6,*)mynum,'Saved mu between:',&
                               t1(ind_win(it,iwin)),t2(ind_win(it,iwin)) 
              write(6,*)mynum,'========================================'
           endif
        endif
     
!---------------------------------------------------------------------------
     if (do_rho) then 
           rhottkern =  rhottkern  + v0(it,reccomp) * rhokern_kxt(it,1:npts)
           rhoampkern = rhoampkern + u0(it,reccomp) * rhokern_kxt(it,1:npts)
!          last point window: scale kernels, save
           if (save_n_erase(it,iwin)) then
              rhottkern = dt*rhottkern*denomv(ind_win(it,iwin))
              rhoampkern = dt*rhoampkern*denomu(ind_win(it,iwin))
              fname='Kx_rho_tt' ! should include filtering etc.....
              call save_kernel(npts,rhottkern,fname,appmywin)
              fname='Kx_rho_am' ! should include filtering etc.....
              call save_kernel(npts,rhoampkern,fname,appmywin)
              write(6,*)mynum,'========================================'
              write(6,*)mynum,'Saved rho tt/amp static kernels at t=',time(it)
              write(6,*)mynum,'Saved rho for time window #', ind_win(it,iwin)
              write(6,*)mynum,'Saved rho between:',&
                               t1(ind_win(it,iwin)),t2(ind_win(it,iwin))
              write(6,*)mynum,'========================================'
           endif
        endif

!---------------------------------------------------------------------------
     if (do_imped) then 
           impttkern =  impttkern  + v0(it,reccomp) * impkernel
           impampkern = impampkern + u0(it,reccomp) * impkernel
!          last point window: scale kernels, save, erase
           if (save_n_erase(it,iwin)) then
              impttkern = dt*impttkern*denomv(ind_win(it,iwin))
              impampkern = dt*impampkern*denomu(ind_win(it,iwin))
              fname='Kx_imp_tt' ! should include filtering etc.....
              call save_kernel(npts,impttkern,fname,appmywin)
              fname='Kx_imp_am' ! should include filtering etc.....
              call save_kernel(npts,impampkern,fname,appmywin)
              write(6,*)mynum,'========================================'
              write(6,*)mynum,'Saved imp tt/amp static kernels at t=',time(it)
              write(6,*)mynum,'Saved imp for time window #', ind_win(it,iwin)
              write(6,*)mynum,'Saved imp between:',&
                               t1(ind_win(it,iwin)),t2(ind_win(it,iwin)) 
              write(6,*)mynum,'========================================'         
           endif
        endif

!---------------------------------------------------------------------------
     if (do_vp) then 
           vpttkern  = vpttkern  + v0(it,reccomp) * vpkernel
           vpampkern = vpampkern + u0(it,reccomp) * vpkernel
!          last point window: scale kernels, save, erase
           if (save_n_erase(it,iwin)) then
              vpttkern = dt*vpttkern*denomv(ind_win(it,iwin))
              vpampkern = dt*vpampkern*denomu(ind_win(it,iwin))
              fname='Kx_vp_tt' ! should include filtering etc.....
              call save_kernel(npts,vpttkern,fname,appmywin)
              fname='Kx_vp_am' ! should include filtering etc.....
              call save_kernel(npts,vpampkern,fname,appmywin)
              write(6,*)mynum,'========================================'
              write(6,*)mynum,'Saved vp tt/amp static kernels at t=',time(it)
              write(6,*)mynum,'Saved vp for time window #', ind_win(it,iwin)
              write(6,*)mynum,'Saved vp between:',&
                               t1(ind_win(it,iwin)),t2(ind_win(it,iwin)) 
              write(6,*)mynum,'========================================'    
           endif
        endif
     
!---------------------------------------------------------------------------
     if (do_vs) then 
           vsttkern  = vsttkern  + v0(it,reccomp) * vskernel
           vsampkern = vsampkern + u0(it,reccomp) * vskernel
!          last point window: scale kernels, save, erase
           if (save_n_erase(it,iwin)) then
              vsttkern = dt*vsttkern*denomv(ind_win(it,iwin))
              vsampkern = dt*vsampkern*denomu(ind_win(it,iwin))
              fname='Kx_vs_tt' ! should include filtering etc.....
              call save_kernel(npts,vsttkern,fname,appmywin)
              fname='Kx_vs_am' ! should include filtering etc.....
              call save_kernel(npts,vsampkern,fname,appmywin)
              write(6,*)mynum,'========================================'
              write(6,*)mynum,'Saved vs tt/amp static kernels at t=',time(it)
              write(6,*)mynum,'Saved vs for time window #', ind_win(it,iwin)
              write(6,*)mynum,'Saved vs between:',&
                               t1(ind_win(it,iwin)),t2(ind_win(it,iwin)) 
              write(6,*)mynum,'========================================'    
           endif
        endif
!---------------------------------------------------------------------------

! setting static kernels to zero at the last time point of a time window
    if (save_n_erase(it,iwin)) then
       if (do_lam) lamttkern = zero; if (do_mu) muttkern = zero
       if (do_rho) rhottkern = zero; if (do_vp) vpttkern = zero
       if (do_vs) vsttkern = zero; if (do_imped) impttkern = zero
       
       if (do_lam) lamampkern = zero; if (do_mu) muampkern = zero
       if (do_rho) rhoampkern = zero; if (do_vp) vpampkern = zero
       if (do_vs) vsampkern = zero; if (do_imped) impampkern = zero
    endif
    if (mynum==0)  write(6,*)''

end subroutine misfit_kernels
!--------------------------------------------------------------

!========================
end module misfit_model_param
!========================
