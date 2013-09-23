!========================
  module posteriori_misfit_operations
!========================
 
  use global_parameters
  use data_mesh
  use data_arrays
  use data_posteriori

  implicit none

  public :: read_posteriori,read_simulation_info
! public :: load_sum_waveform_kernels,frequency_kernels,change_wfkern,prepare_misfits
  public :: load_sum_waveform_kernels,frequency_kernels,prepare_misfits
  public :: reconst_model_param_kernel,compute_misfit_kernel, filter_kernels_post
  private

  contains 

!----------------------------------------------------------------------------------------------------------
subroutine read_posteriori

use data_fft, only : filter_type,ntimes,nomega,uphys1d,uspec1d
use parameters, only: define_sem_kernel_meshes
use data_misfits

!include 'mesh_params_kernel.h'
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
   allocate(post_per(num_filt,2))
   allocate(filter_type_posteriori(num_filt))
   allocate(app_period1(num_filt))
   allocate(app_period2(num_filt))
   allocate(filter_char(num_filt))   
   do i=1,num_filt
      read(99,*)filter_type_posteriori(i),post_per(i,1),post_per(i,2)
      write(6,*)'filter type and min/max period: ',filter_type_posteriori(i),post_per(i,1),post_per(i,2)
      if (filter_type_posteriori(i)/='non') need_filt=need_filt+1
      call define_io_appendix3_real(app_period1(i),post_per(i,1))
      call define_io_appendix3_real(app_period2(i),post_per(i,2))
      filter_char(i) = trim(filter_type_posteriori(i))//'_'//app_period1(i)//'_'//app_period2(i)
      write(6,*)'filter char:',filter_char(i)
   enddo
endif
if (num_filt==0) num_filt=1  ! just for convenience in the loop in main.f90

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
subroutine read_simulation_info
use data_fft, only : filter_type,ntimes,nomega,uphys1d,uspec1d
use parameters, only: define_sem_kernel_meshes
use data_misfits

!include 'mesh_params_kernel.h'
include 'mesh_params.h'

integer :: i,ii,j,iwin,iidom,npad
character(len=100) :: junk
real(kind=realkind) :: src_depthbwd
real(kind=realkind) :: thr,phr
character(len=4) :: model_reconsttmp(2),model_runtmp(20)
character(len=5) :: misfittmp(20)
logical :: travelpicktmp(20)

write(6,*)'load general simulation info...'

allocate(src_type1fwd_post(num_dir),src_type2fwd_post(num_dir))
allocate(src_type1bwd_post(num_dir),src_type2bwd_post(num_dir))
allocate(fwdstf_type_post(num_dir),bwdstf_type_post(num_dir))
allocate(dim_fwd_post(num_dir),dim_bwd_post(num_dir))
do i=1,num_dir
   open(unit=100,file=trim(run_dir(i))//'/Data/simulation.info')
   read(100,*)junk,nproc; write(6,*)'nprocs=',nproc
   read(100,*)junk,npts; write(6,*)'npts=',npts
   read(100,*)junk,src_type1fwd_post(i),src_type2fwd_post(i)
   read(100,*)junk,src_depth
   read(100,*)junk,fwdstf_type_post(i)
   read(100,*)junk,dim_fwd_post(i)
   read(100,*)junk,src_type1bwd_post(i),src_type2bwd_post(i)
   read(100,*)junk,src_depthbwd
   read(100,*)junk,bwdstf_type_post(i)
   read(100,*)junk,dim_bwd_post(i)
   read(100,*)junk,thr
   read(100,*)junk,phr
   write(6,*)'theta,phi:',thr,phr
   read(100,*)junk,save_wfkern; write(6,*)'save wf kern:',save_wfkern
   read(100,*)junk,do_rho,do_lam,do_mu; write(6,*)'rho,lam,mu:',do_rho,do_lam,do_mu
   read(100,*)junk,do_vp,do_vs,do_imped; write(6,*)'rho,lam,mu:',do_vp,do_vs,do_imped
   read(100,*)junk,dt; write(6,*)'dt=',dt
   read(100,*)junk,twf1,twf2; write(6,*)'twf1,twf2:',twf1,twf2
   read(100,*)junk,nwin; write(6,*)'number of time windows:',nwin
   allocate(wft1(nwin),wft2(nwin),begdumps(nwin),enddumps(nwin))
   write(6,*)' Waveform kernel time window:'
   do iwin=1,nwin
      read(100,*)junk,ii,begdumps(iwin),enddumps(iwin),wft1(iwin),wft2(iwin)
      write(6,*)'wfkern window times:',wft1(iwin),wft2(iwin)
      write(6,*)'wfkern window indices:',begdumps(iwin),enddumps(iwin)
   enddo
   read(100,*)junk,ntw
   write(6,*)'Time windows:',ntw
   do j=1,ntw
      read(100,*)junk,ii,t1(j),t2(j)
      write(6,*)t1(j),t2(j)
   enddo
   read(100,*)junk,num_it
   write(6,*)'number of snaps:',num_it
   ii=0
   if (i==1) allocate(dump_ind(ndumps,num_dir),dump_time(ndumps,num_dir))
   dump_ind=0; dump_time=0.
   do iwin=1,nwin
      do j=begdumps(iwin),enddumps(iwin)
!         write(6,*)'time window index:',iwin,j
         read(100,*)junk,dump_ind(j,i),dump_time(j,i)
      enddo
   enddo
close(100)
write(6,*)i,'source type fwd:',src_type1fwd_post(i),src_type2fwd_post(i)
write(6,*)i,'source type bwd:',src_type1bwd_post(i),src_type2bwd_post(i)

write(6,*)' determine receiver component for run:',run_dir(i)
allocate(reccomp_post(num_dir))
reccomp_post=3 ! radial by default
if (src_type1fwd_post(i)=='monopole') reccomp_post(i)=2

if (src_type2bwd_post(i)=='vertforce' .or. src_type2bwd_post(i)=='explosion') then
   if (src_type1fwd_post(i)=='monopole') then
      reccomp_post(i)=2
   else
      reccomp_post(i)=3
    endif
elseif (src_type2bwd_post(i)=='xforce') then
   reccomp_post(i)=1 ! theta
elseif (src_type2bwd_post(i)=='yforce') then
      if (src_type1fwd_post(i)=='monopole') then
         write(6,*) 'This combination does not make sense:'
         write(6,*)'monopole fwd source will not have any signal on the transverse bwd field...'
         stop
   else
      reccomp_post(i)=2 !phi
    endif
 endif
write(6,*)'receiver component:',reccomp_post(i)
enddo

! model parameterization
write(6,*)'model parameterization.....'
num_model_param=0
num_model_reconst=0
if (do_lam_post .and. .not. do_lam) then 
   write(6,*)'cannot compute lambda kernels - not existent!'
   do_lam_post=.false.
endif

if (do_mu_post .and. .not. do_mu) then 
   write(6,*)'cannot compute mu kernels - not existent!'
   do_mu_post=.false.
endif

if (do_rho_post .and. .not. do_rho) then 
   write(6,*)'cannot compute rho kernels - not existent!'
   do_rho_post=.false.
endif

if (do_vp_post .and. .not. do_vp) then 
   if (.not. do_lam) then
      write(6,*)'cannot compute vp kernels - lam nor vp existent!'
      do_vp_post=.false.
   else
      write(6,*)'need to construct vp kernels from lambda kernels'
      num_model_reconst=num_model_reconst+1
      model_reconsttmp(num_model_reconst)='v_p'
   endif
endif

if (do_vs_post .and. .not. do_vs ) then 
   if (.not. do_mu .or. .not. do_lam) then
      write(6,*)'cannot compute vs kernels - lam,mu,vs not existent!'
      do_vs_post=.false.
   elseif (do_mu .and. do_lam) then
      write(6,*)'need to construct vs kernels from lambda and mu kernels'
      num_model_reconst=num_model_reconst+1
      model_reconsttmp(num_model_reconst)='v_s'
   endif
endif

if (do_imped_post .and. .not. do_imped) then 
   if (.not. do_lam .or. .not. do_mu .or. .not. do_rho) then
      write(6,*)'cannot compute imped kernels - imped,rho,lam,mu not existent!'
      do_imped_post=.false.
   else
      write(6,*)'need to construct imped kernels from lam,mu,rho kernels'
      num_model_reconst=num_model_reconst+1
      model_reconsttmp(num_model_reconst)='imp'
   endif
endif

if (do_lam_post) then 
   num_model_param=num_model_param+1
   model_runtmp(num_model_param)='lam'
endif
if (do_mu_post) then 
   num_model_param=num_model_param+1
   model_runtmp(num_model_param)='mu_'
endif
if (do_rho_post) then 
   num_model_param=num_model_param+1
   model_runtmp(num_model_param)='rho'
endif
if (do_vp_post) then 
   num_model_param=num_model_param+1
   model_runtmp(num_model_param)='v_p'
endif
if (do_vs_post) then 
   num_model_param=num_model_param+1
   model_runtmp(num_model_param)='v_p'
endif
if (do_imped_post) then 
   num_model_param=num_model_param+1
   model_runtmp(num_model_param)='imp'
endif

allocate(model_run(1:num_model_param),model_reconst(1:num_model_reconst))
model_run(1:num_model_param)=model_runtmp(1:num_model_param)
model_reconst(1:num_model_reconst)=model_reconsttmp(1:num_model_reconst)

! time
do i=1,ndumps
   time(i)=dt*real(i)
enddo

write(6,*)' relate moment tensor to the different simulations....'
allocate(Mijsum(num_dir))

! this includes the case of having e.g. the same Mij, but different receiver wavefields, 
! i.e. stacking components. For now, they are simply 
do i=1,num_dir
   select case(src_type2fwd_post(i))
   case('explosion') 
      Mijsum(i) = Mij_post(1)+Mij_post(2)+Mij_post(3)
   case('mzz') 
      Mijsum(i) = Mij_post(1)
   case('mxx') 
      Mijsum(i) = Mij_post(2)
   case('myy') 
      Mijsum(i) = Mij_post(3)
   case('mxz') 
      Mijsum(i) = Mij_post(4)
   case('myz') 
      Mijsum(i) = Mij_post(5)
   case('mxy') 
      Mijsum(i) = Mij_post(6)
   case default
      write(6,*)trim(run_dir(i)),'has an undefined source type:',src_type2fwd_post(i)
      stop
   end select
enddo

! number of misfits to be computed
num_misfits=0
travelpicktmp(1:20)=.false.
hilbert_needed=.false.
gabor_needed=.false.
if (xc_tt) then 
   num_misfits=num_misfits+1
   misfittmp(num_misfits)='xc_tt'
endif
if (xc_am) then 
   num_misfits=num_misfits+1
   misfittmp(num_misfits)='xc_am'
endif
if (on_tt) then 
   num_misfits=num_misfits+1
   misfittmp(num_misfits)='on_tt'
   travelpicktmp(num_misfits)=.true.
endif
if (pk_tt) then 
   num_misfits=num_misfits+1
   misfittmp(num_misfits)='pk_tt'
   travelpicktmp(num_misfits)=.true.
endif
if (insta) then 
   num_misfits=num_misfits+1
   misfittmp(num_misfits)='insta'
   travelpicktmp(num_misfits)=.true.
endif
if (wf_dd) then 
   if (consider_data) then 
      num_misfits=num_misfits+1
      misfittmp(num_misfits)='wf_dd'
   else
      write(6,*)'need data for waveform misfit!'
      write(6,*)'... therefore IGNORING this misfit'
      wf_dd=.false.
   endif
endif
if (instph) then 
   if (consider_data) then 
      num_misfits=num_misfits+1
      misfittmp(num_misfits)='hilph'
      hilbert_needed=.true.
   else
      write(6,*)'need data for instantaneous phase misfit!'
      write(6,*)'... therefore IGNORING this misfit'
      instph=.false.
   endif
endif
if (insten) then 
   if (consider_data) then 
      num_misfits=num_misfits+1
      misfittmp(num_misfits)='hilen'
      hilbert_needed=.true.
   else
      write(6,*)'need data for instantaneous envelope misfit!'
      write(6,*)'... therefore IGNORING this misfit'
      insten=.false.
   endif
endif
if (phase) then 
   if (consider_data) then
      num_misfits=num_misfits+1
      misfittmp(num_misfits)='phase'
      gabor_needed=.true.
   else
      write(6,*)'need data for phase misfit!'
      write(6,*)'... therefore IGNORING this misfit'
      phase=.false.
   endif
endif
if (envel) then 
   if (consider_data) then 
      num_misfits=num_misfits+1
      misfittmp(num_misfits)='envel'
      gabor_needed=.true.
   else
      write(6,*)'need data for envelope misfit!'
      write(6,*)'... therefore IGNORING this misfit'
      envel=.false.
   endif
endif

allocate(misfit(num_misfits),travelpick(1:num_misfits))
misfit(1:num_misfits)=misfittmp(1:num_misfits)
travelpick(1:num_misfits)=travelpicktmp(1:num_misfits)

write(6,*)'determine number of picks for each misfit...'
allocate(num_picks(num_misfits),time_pick(1:num_tot_picks,num_misfits))
num_picks=0.
time_pick=0.
do j=1,num_misfits
   do ii=1,num_tot_picks
      if (misfit(j)==pick_misfit(ii))then 
         num_picks(j)=num_picks(j)+1
         time_pick(num_picks(j),j)=time_pick_read(ii)
         write(6,*)'misfits & picks:',misfit(j),' ',pick_misfit(ii),num_picks(j)
      endif
   enddo
enddo

! decide whether frequency-domain kernels are necessary
if (pk_tt .or. on_tt .or. gabor_needed .or. hilbert_needed .or. need_filt>0) then 
!   if (twf1>1. .and. twf2 < time(ndumps)) then
!      write(6,*)'.... waveform kernels have not been saved for the entire time history!'
!      write(6,*)'.... therefore NOT doing any FFTW and frequency-domain operations!'
!      write(6,*)'.... Now turning filtering off, and not computing misfits pk_tt, on_tt'
!      need_omkern=.false.
!      filter_type_posteriori(1:num_filt)='non'
!      need_filt=0
!      pk_tt=.false.
!      on_tt=.false.
!   else
      need_omkern=.true.
!   endif
   write(6,*)'We need to go into frequency domain...'
else
   need_omkern=.false.
endif

! get kernel mesh directories
call define_sem_kernel_meshes

write(6,*)'read background model domains and discontinuities....'
  write(6,*)'dir fwd mesh:',dir_fwdmesh(1:lffwd)
  open(unit=65,file=trim(dir_fwdmesh1(1))//'/Info' &
              //'/discontinuities_solver.dat0000',status='old')
 do iidom=1,ndisc
        read(65,*)discont(iidom)
      if (mynum==0)write(6,*)iidom,'discontinuity:',discont(iidom)
  enddo
  close(65)

write(6,*)'read kernel mesh for output in vtk...'
allocate(W_vtk(1:npts,3))
open(unit=50,file=trim(run_dir(1))//'/Data/'//'kernelmesh_xyz_'//appmynum//'.dat',status='old')
do i=1,npts
   read(50,*)W_vtk(i,1),W_vtk(i,2),W_vtk(i,3)
enddo
close(50)

! allocate various arrays
  npad = ndumps !/2 ! TNM added /2 to make it cheaper...
  ntimes=ndumps+1 + npad
  nomega=(ntimes-1)/2 ! AF from the definition provided in the fftw manual

  if (mynum==0) write(6,*)'number of frequencies:',nomega
  if (mynum==0) write(6,*)'number of points:',npts
!
end subroutine read_simulation_info
!----------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------
subroutine load_sum_waveform_kernels(modelname,ipar)

use data_arrays
use data_fft, only : ntimes
use data_mesh, only : ndumps
use data_misfits 

character(len=3), intent(in) :: modelname
integer, intent(in) :: ipar
integer :: i,ii,iwin,it
character(len=20) :: junk
character(len=4) :: appidump1
character(len=100) :: filename

! Load waveform kernels
if (ipar==1) allocate(kern_all_phys(1:ntimes,1:npts))

write(6,*)'loading waveform kernels......',npts

do i=1,num_dir
   do iwin=1,nwin
      do it=begdumps(iwin),enddumps(iwin)
         call define_io_appendix(appidump1,it)
!         write(6,*)'INDICES:',i,iwin,it,begdumps(1)
!         write(6,*)'time window:',it,appidump1; call flush(6)
         if (do_lam) then
!            open(unit=100,file=trim(run_dir(i))//'/Data/'//modelname//'_wf_kern_'//appmynum//'_'//appidump1//'.bin',&
!                 status='old',position='rewind',form='unformatted')
!            read(100)kern_all_phys(it,1:npts)
            filename=trim(run_dir(i))//'/Data/'//modelname//'_wf_kern_'//appmynum//'_'//appidump1//'.bin'
!               write(6,*)'waveform kernel ',trim(filename)
               write(6666,*)i,iwin,it,ii
            open(unit=100,file=trim(filename),form='unformatted',POSITION="REWIND",status='old')
               read(100)kern_all_phys(it,1:npts)
            close(100)
         endif
      enddo
   enddo ! num_it: waveform kernel snapshots
   
! sum 
  kern_all_phys = Mijsum(i)*kern_all_phys + kern_all_phys

enddo ! num_dir: different simulations

end subroutine load_sum_waveform_kernels
!----------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------------------
subroutine frequency_kernels(thr,phr)
use data_fft, only : nomega
use fft
include 'mesh_params_kernel.h'

real(kind=realkind), intent(in) :: thr,phr

write(6,*)'computing kernels in frequency domain....'

kerphys= kern_all_phys
write(6,*)'fftf...'
call fftf_ker(kern_all_spec(0:nomega,1:npts))
write(6,*)'...done.'

end subroutine frequency_kernels
!----------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------------------
subroutine filter_kernels_post(ifilt)

use data_fft
use filtering

integer, intent(in) :: ifilt
write(6,*)'filtering with ',filter_type_posteriori(ifilt),' and periods',post_per(ifilt,1),post_per(ifilt,2)
kerspec=kern_all_spec
write(6,*)'calling the filter...'
!call af_filt_dble(filter_type_posteriori(ifilt),1./post_per(ifilt,2),1./post_per(ifilt,1),6,npts,nomega,omega,kerspec)

call filter_field(kerspec)
write(6,*)'done with the filter...'

end subroutine filter_kernels_post
!----------------------------------------------------------------------------------------------------------
!
!subroutine change_wfkern(wfkerntype)
!use frequency_domain, only : compute_deriv_omkern_pk_tt,compute_3deriv_omkern_on_tt
!!this is needed if the misfit kernel does not depend on the raw K(x,t),
!! but e.g. a time derivative or Gabor transform of it

!character(len=3), intent(in) :: wfkerntype

!! choose correct kernel type
!if (wfkerntype=='1st') then 
!   call compute_deriv_omkern_pk_tt(kerspec,kerspec)
!elseif(wfkerntype=='3rd') then 
!   call compute_3deriv_omkern_on_tt(kerspec,kerspec)
!elseif(wfkerntype=='gab') then 
!   write(6,*)'not done yet....'
!   stop
!endif

!end subroutine change_wfkern

!----------------------------------------------------------------------------------------------------------
subroutine prepare_misfits(thr,phr)
use data_fft, only : v0spec,nomega,omega,ntimes
use data_misfits
use fft, only : fftf1d,fftb1d
use misfit_model_param, only : init_time_window_misfit_medium
use filtering

!include 'mesh_params_kernel.h'

real(kind=realkind) :: thr,phr
integer :: j,it,i,iw,iwin,ifilt
real(kind=realkind), allocatable :: inst_phase_syn(:),inst_phase_dat(:),utmp(:)
real(kind=realkind), allocatable :: E_syn(:),E_dat(:),hilb_syn(:),hilb_dat(:)
complex(kind=8), allocatable :: gab_syn(:,:),gab_dat(:,:)
real(kind=realkind), allocatable :: tf_env_syn(:,:),tf_env_dat(:,:),tf_phase_syn(:,:),tf_phase_dat(:,:)

character(len=200) :: filename
character(len=200) :: appifilt1,appifilt2

allocate(wfkernel_type(1:num_misfits))
wfkernel_type(:) = 'reg'

write(6,*)'determine picks in the discrete time series...........'
write(6,*)'# misfits:',num_misfits
write(6,*)'# max picks',maxval(num_picks)
write(6,*)'picks:',(num_picks(j),j=1,num_misfits)
write(6,*)'travelpick:',(travelpick(j),j=1,num_misfits)

allocate(pick_ind(1:maxval(num_picks),num_misfits))
pick_ind=0
do j=1,num_misfits
   if (travelpick(j)) then 
      do it=1,num_picks(j)
         write(6,*)'picks',j,it,num_picks(j)
         pick_ind(it,j)=min(max(minloc(abs(time(1:ndumps)-time_pick(it,j)),1),1),ndumps)
         write(6,*)'Time pick expected:',time_pick(it,j)
         write(6,*)'Time pick offered:',time(pick_ind(it,j))
      enddo
   endif
enddo

! specific to each misfit
allocate(misfit_prefactor(1:ndumps,1:num_misfits,1:num_filt))
num_windows=max(ntw,maxval(num_picks))
allocate(norm_fac(1:num_windows,1:num_misfits,1:num_filt),ttdelay(num_windows,num_misfits))
allocate(ampdiff(num_windows,num_misfits))
ttdelay=0.; ampdiff=0.
misfit_prefactor=0.
norm_fac=1.

if (need_omkern .or. hilbert_needed) then
   allocate(a0(1:ntimes), v0spec(0:nomega),utmp(1:ndumps))
   a0=0.0; v0spec=0.; utmp=0.
endif

if (gabor_needed) then
   allocate(gab_syn(ndumps,0:nomega),gab_dat(ndumps,0:nomega))
   allocate(tf_env_syn(ndumps,0:nomega),tf_env_dat(ndumps,0:nomega))
   allocate(tf_phase_syn(ndumps,0:nomega),tf_phase_dat(ndumps,0:nomega))
endif

if (hilbert_needed) then
   allocate(E_syn(ndumps),E_dat(ndumps),hilb_syn(ndumps),hilb_dat(ndumps))
   allocate(inst_phase_syn(ndumps),inst_phase_dat(ndumps))
endif

write(6,*)'calculating misfit prefactors......'
write(6,*)size(misfit_prefactor),size(v0post)
misfit_fct=0.
open(unit=50,file='Data/cost_functions.dat')

!=========================================
do ifilt=1,num_filt   ! START BIG LOOP
!=========================================

   call define_io_appendix(appifilt1,post_per(ifilt,1))
   call define_io_appendix(appifilt2,post_per(ifilt,2))

   ! Hilbert transform, instantaneous phase and envelope
   if (hilbert_needed) then 
      write(6,*)'Calculating Hilbert transform & instantaneous phase/envelope for filter',ifilt
      call hilbert_insta_phase_env(usyn_filt(:,ifilt),ndumps,nomega,hilb_syn,inst_phase_syn,E_syn)
      call hilbert_insta_phase_env(udat_filt(:,ifilt),ndumps,nomega,hilb_dat,inst_phase_dat,E_dat)
      write(6,*)'...done with Hilbert transform.'
      if (mynum==0) then 
         open(unit=100,file='Data/inst_phase_syn_'//appifilt1//'_'//appifilt2//'.dat')
         open(unit=101,file='Data/inst_phase_dat_'//appifilt1//'_'//appifilt2//'.dat')
         open(unit=102,file='Data/inst_envel_syn_'//appifilt1//'_'//appifilt2//'.dat')
         open(unit=103,file='Data/inst_envel_dat_'//appifilt1//'_'//appifilt2//'.dat')
         do i=1,ndumps
            write(100,*)time(i),inst_phase_syn(i)
            write(101,*)time(i),inst_phase_dat(i)
            write(102,*)time(i),E_syn(i)
            write(103,*)time(i),E_dat(i)
         enddo
         close(100);close(101);close(102);close(103)
      endif
   endif

   ! Gabor transform for time-frequency misfits
   if (gabor_needed) then 
      write(6,*)'Calculating Gabor transform for time-frequency phase/envelope...'
      call gabor_trafo(udat_filt(:,ifilt),ndumps,nomega,post_per(ifilt,1),gab_dat,tf_phase_syn,tf_env_syn)
      call gabor_trafo(usyn_filt(:,ifilt),ndumps,nomega,post_per(ifilt,1),gab_syn,tf_phase_dat,tf_env_dat)
      if (mynum==0 ) then
         open(unit=90,file='Data/gabor_trafo_displ_syn_'//appifilt1//'_'//appifilt2//'.dat')
         open(unit=93,file='Data/gabor_trafo_displ_dat_'//appifilt1//'_'//appifilt2//'.dat')
         if (ifilt==1) then
            open(unit=91,file='Data/time.dat')
            open(unit=92,file='Data/omega.dat')
         endif
         open(unit=100,file='Data/tf_phase_syn_'//appifilt1//'_'//appifilt2//'.dat')
         open(unit=101,file='Data/tf_phase_dat_'//appifilt1//'_'//appifilt2//'.dat')
         open(unit=102,file='Data/tf_envel_syn_'//appifilt1//'_'//appifilt2//'.dat')
         open(unit=103,file='Data/tf_envel_dat_'//appifilt1//'_'//appifilt2//'.dat')
         do i=1,ndumps
            if (ifilt==1) write(91,*)time(i)
            do iw=0,nomega
               write(90,*)real(gab_syn(i,iw))
               write(93,*)real(gab_dat(i,iw))
               write(100,*)time(i),tf_phase_syn(i,iw)
               write(101,*)time(i),tf_phase_dat(i,iw)
               write(102,*)time(i),tf_env_syn(i,iw)
               write(103,*)time(i),tf_env_dat(i,iw)
               if (i==1 .and. ifilt==1)write(92,*)omega(iw)
            enddo
         enddo
         close(90)
         if (ifilt==1) then 
            close(91);close(92)
         endif
         close(93)
         close(100);close(101);close(102);close(103)
         write(6,*)'Done with Gabor transform.'
      endif
   endif

!------------------------------------
   do j=1,num_misfits
!------------------------------------
      write(6,*)'MISFIT & FILTERING: ',misfit(j),' filter periods:',post_per(ifilt,1),post_per(ifilt,2)

      if (misfit(j)=='xc_tt') then 
         misfit_prefactor(1:ndumps,j,ifilt)=vsyn_filt(:,ifilt)
         do it=1,ntw               
            ttdelay(it,j)=1.
            if (consider_data) call traveltime_delay(t1(it),t2(it),time,usyn_filt(:,ifilt),vsyn_filt(i:,ifilt),&
                                                                              udat_filt(:,ifilt),ndumps,ttdelay(it,j))
            norm_fac(it,j,ifilt)=denomv(it)*ttdelay(it,j)
            misfit_fct(it,j)=ttdelay(it,j)
         enddo
      endif

      if (misfit(j)=='xc_am') then 
         misfit_prefactor(1:ndumps,j,ifilt)=usyn_filt(:,ifilt)
         do it=1,ntw
            ampdiff(it,j)=1.
            if (consider_data)  call amplitude_difference(t1(it),t2(it),time,usyn_filt(:,ifilt),&
                                                                                        udat_filt(:,ifilt),ndumps,ampdiff(it,j))
            norm_fac(it,j,ifilt)=denomu(it)*ampdiff(it,j)
            misfit_fct(it,j)=ampdiff(it,j)
         enddo
      endif

      if (misfit(j)=='wf_dd') then 
         misfit_prefactor(1:ndumps,j,ifilt) = usyn_filt(:,ifilt)-udat_filt(:,ifilt)
         norm_fac(it,j,ifilt) = denomu(it)
         misfit_fct(:,j) = usyn_filt(:,ifilt)-udat_filt(:,ifilt)
      endif

      if (misfit(j)=='hilen') then
         misfit_prefactor(1:ndumps,j,ifilt)= log((E_dat+epsi)/(E_syn+epsi)) * usyn_filt(:,ifilt)/(E_syn+smallval_dble)**2
         utmp(1:ndumps) = log((E_dat+epsi)/(E_syn+epsi)) *hilb_syn/(E_syn+smallval_dble)**2
         call hilbert_insta_phase_env(utmp,ndumps,nomega,hilb_dat,utmp,utmp)
         misfit_prefactor(1:ndumps,j,ifilt) = misfit_prefactor(1:ndumps,j,ifilt) - hilb_dat
         misfit_fct(:,j) = log((E_dat+epsi)/(E_syn+epsi))
      endif
      
      if (misfit(j)=='hilph' ) then
         misfit_prefactor(1:ndumps,j,ifilt)= (inst_phase_dat-inst_phase_syn)*hilb_syn/(E_syn+smallval_dble)**2 
         utmp(1:ndumps) = (inst_phase_dat-inst_phase_syn)*usyn_filt(:,ifilt)/(E_syn+smallval_dble)**2
         call hilbert_insta_phase_env(utmp,ndumps,nomega,hilb_dat,utmp,utmp)
         misfit_prefactor(1:ndumps,j,ifilt) = misfit_prefactor(1:ndumps,j,ifilt) + hilb_dat
         misfit_fct(:,j) = inst_phase_dat-inst_phase_syn
      endif

      if (misfit(j)=='phase') then
         wfkernel_type(j)='gab'
         write(6,*)misfit(j),' not yet implemented';stop
      endif
      
      if (misfit(j)=='envel') then
         wfkernel_type(j)='gab'
         write(6,*)misfit(j),' not yet implemented';stop
      endif

      if (misfit(j)=='pk_tt' .or. misfit(j)=='on_tt' .or. misfit(j)=='insta') then
         write(6,*)'calculating traveltime pick misfits and prefactors......'

         do it=1,num_picks(j) ! NEED TO CHECK WHETHER THIS NUMBERING IS ALL CORRECT.....

         ! TRAVELTIME DELAY SHOULD BE CALCULATED DIFFERENTLY FOR EACH MISFIT!!!!!

            if (misfit(j)=='pk_tt') then
               wfkernel_type(j)='1st'
               misfit_prefactor(pick_ind(it,j),j,ifilt) = 1.
               if (consider_data) call traveltime_delay(t1(it),t2(it),time,usyn_filt(:,ifilt),vsyn_filt(:,ifilt),&
                                                               udat_filt(:,ifilt),ndumps,misfit_prefactor(pick_ind(it,j),j,ifilt) )
               write(6,*)'misfit,pick index, ttdelay: ',misfit(j),' ',pick_ind(it,j),misfit_prefactor(pick_ind(it,j),j,ifilt) 
               call fftf1d(vsyn_filt(:,ifilt),v0spec)
               do i=1,nomega
                  v0spec(i)=cmplx(0.,1.)*omega(i)*v0spec(i) ! acceleration
               enddo
               call fftb1d(v0spec,a0)
               if (a0(pick_ind(it,j))/=0.) norm_fac(it,j,ifilt)= 1./(a0(pick_ind(it,j))**2)
            endif
            
            if (misfit(j)=='on_tt') then
               wfkernel_type(j)='3rd'
               misfit_prefactor(pick_ind(it,j),j,ifilt) = 1.
               if (consider_data) call traveltime_delay(t1(it),t2(it),time,usyn_filt(:,ifilt),vsyn_filt(:,ifilt),&
                                                              udat_filt(:,ifilt),ndumps,misfit_prefactor(pick_ind(it,j),j,ifilt) )
               write(6,*)'misfit,pick index, ttdelay: ',misfit(j),' ',pick_ind(it,j),misfit_prefactor(pick_ind(it,j),j,ifilt) 
               call fftf1d(vsyn_filt(:,ifilt),v0spec)
               do i=0,nomega
                  v0spec(i)=-cmplx(0.,1.)*omega(i)**3*v0spec(i) ! 4th derivative !!
               enddo
               call fftb1d(v0spec,a0)
               if (a0(pick_ind(it,j))/=0.) norm_fac(it,j,ifilt)= 1./(a0(pick_ind(it,j))**2)
            endif
            
            if (misfit(j)=='insta') then
               misfit_prefactor(pick_ind(it,j),j,ifilt)=usyn_filt(pick_ind(it,j),ifilt)-udat_filt(pick_ind(it,j),ifilt)
               write(6,*)'INSTA:',pick_ind(it,j),time(pick_ind(it,j)),maxval(misfit_prefactor(:,j,ifilt))
               if (usyn_filt(pick_ind(it,j),ifilt)/=0.) norm_fac(it,j,ifilt)=1./(usyn_filt(pick_ind(it,j),ifilt)**2)
            endif
            misfit_fct(pick_ind(it,j),j)=misfit_prefactor(pick_ind(it,j),j,ifilt) 

         enddo !num_picks
      endif ! picking misfits
      
      ! write misfit function and misfit prefactor to file (both time-dependent)
      open(unit=60,file='Data/misfit_fct_'//trim(misfit(j))//'_'//appifilt1//'_'//appifilt2//'.dat')
      filename='Data/misfit_prefact_'//trim(misfit(j))//'_'//&
                      trim(filter_type_posteriori(ifilt))//'_'//appifilt1//'_'//appifilt2//'.dat'
      open(unit=61,file=trim(filename))
      do it=1,ndumps
         write(60,*)time(it),misfit_fct(it,j)
         write(61,*)time(it),timewin(it)*misfit_prefactor(it,j,ifilt)/maxval(abs(timewin*misfit_prefactor(:,j,ifilt)))
      enddo
      close(60);close(61)

      ! calculate misfit scalar (i.e. the actual value of Chi) for each time window
      chi(:,j) = 0.      
      do iwin = 1, nwin ! number of time windows. 1 if all waveform kernels are saved.
         do it=begdumps(iwin),enddumps(iwin) 
            chi(iwin,j) = chi(iwin,j) + timewin(it)*misfit_fct(it,j)**2
         enddo
         chi(iwin,j)=0.5*norm_fac(iwin,j,ifilt)*chi(iwin,j)*dt
      enddo
!------------------------------------
      enddo !misfits 
!------------------------------------
      if (mynum==0 .and. ifilt==1) write(50,*)(misfit(j),j=1,num_misfits)
      if (mynum==0) write(50,*)(chi(iwin,j),j=1,num_misfits)

!=========================================
enddo !filt
!=========================================

close(50)
if (hilbert_needed) deallocate(E_syn,E_dat,hilb_syn,hilb_dat,inst_phase_syn,inst_phase_dat)
if (gabor_needed)  deallocate(tf_env_syn,tf_env_dat,gab_syn,gab_dat,tf_phase_syn,tf_phase_dat)
if (need_omkern) deallocate(a0,v0spec,utmp)

norm_fac=norm_fac*dt
write(6,*)'..... done prepping misfits.'; call flush(6)

allocate(kernel(1:npts))
if (need_omkern) allocate(kern_all_spec(0:nomega,1:npts))
end subroutine prepare_misfits
!----------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------
subroutine hilbert_insta_phase_env(uin,n,nom,hilbert,inst_phase,inst_env)

use fft, only : fftf1d,fftb1d

implicit none

integer, intent(in) :: n,nom
real(kind=realkind), intent(in) :: uin(n)
real(kind=realkind) :: h(n)
real(kind=realkind), intent(out) :: hilbert(n)
real(kind=realkind), intent(out) :: inst_phase(n),inst_env(n)
complex(kind=8) :: analyt_sign(n),hspec(0:nom),uspec(0:nom)
integer :: i

  do i=1,n
     h(i)=1./(pi*i*dt)
  enddo
  call fftf1d(h,hspec)
  
  call fftf1d(uin,uspec)
  do i=0,nom
     uspec(i)=hspec(i)*uspec(i) 
  enddo
  call fftb1d(uspec,hilbert)
  
  analyt_sign = dble(uin) + cmplx(0.,1.)*dble(hilbert)
  inst_env = sqrt( real(analyt_sign)**2 + aimag(analyt_sign)**2)
  inst_phase = atan2(aimag(analyt_sign),real(analyt_sign))

end subroutine hilbert_insta_phase_env
!----------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------
subroutine gabor_trafo(uin,nt,nfreq,T0,uout,phase,env)

use data_fft, only : omega 

implicit none

integer, intent(in) :: nt,nfreq
real(kind=realkind), intent(in) :: uin(nt),T0
real(kind=realkind), intent(out) :: phase(nt,0:nfreq),env(nt,0:nfreq)
complex(kind=8), intent(out) :: uout(nt,0:nfreq)
integer :: it,iw,itt
real(kind=realkind) :: prefac,sigma,tau,t

sigma=T0
prefac=(4.*pi**3*sigma)**(-0.25)
do iw=0,nfreq
   do it=1,nt
      t=(it-1)*dt
      do itt=1,nt
         tau=(itt-1)*dt
         uout(it,iw)= uout(it,iw)+uin(it)*exp(-cmplx(0.,1.)*omega(iw)*tau-((tau-t)/(2.*sigma))**2)
      enddo
      uout(it,iw)=uout(it,iw)*dt
   enddo
enddo

env=abs(uout)
phase=log(uout/env)

end subroutine gabor_trafo

!----------------------------------------------------------------------------------------------------------
subroutine compute_misfit_kernel(ipar,ifilt,j)

use input_output, only : save_kernel,save_misfit_kernel
use data_misfits 

integer, intent(in) :: ipar,ifilt,j
integer :: it,iwin,i
character(len=4) :: appidump1,appmywin
character(len=6) :: filename
character(len=30) :: filename2

!!!! this still doesn't honor the case of keeping primitive kernels from the earlier 
!!!!! loop to combine into, e.g. vs-kernel...

! for misfits based on waveforms (i.e. time windows)
   write(6,*)'computing misfit kernel for ',model_run(ipar),' ',misfit(j),' ',filter_char(ifilt)

   ! traveltime picks
   if (travelpick(j)) then 

      do iwin=1,num_picks(j)
         it=pick_ind(iwin,j)
         write(6,*)'Misfit & prefactor:',misfit(j),norm_fac(iwin,j,ifilt),& 
                          misfit_prefactor(it,j,ifilt),maxval(kerphys(it,:)),maxval(kernel)
         kernel = norm_fac(iwin,j,ifilt)*misfit_prefactor(it,j,ifilt)*real(kerphys(it,1:npts))
         call define_io_appendix(appmywin,it)
         if (num_filt>0) then 
            filename2=trim(model_run(ipar))//'_'//trim(misfit(j))//'_'//trim(filter_char(ifilt))
         else
            filename2=trim(model_run(ipar))//'_'//trim(misfit(j))
         endif
         call save_misfit_kernel(npts,real(kernel)/maxval(real(abs(kernel))),filename2,appmywin)
      enddo

   elseif (misfit(j)=='xc_tt' .or. misfit(j)=='xc_am' .or. misfit(j)=='wf_dd'  .or. misfit(j)=='hilph' .or. &
             misfit(j)=='hilen' .or. misfit(j)=='phase' .or. misfit(j)=='envel') then 
      do iwin = 1, nwin ! number of time windows. 1 if all waveform kernels are saved.
         do it=begdumps(iwin),enddumps(iwin) 
            call define_io_appendix(appidump1,it) 
            if (num_filt>0) then 
               filename2=trim(model_run(ipar))//'_wf'//'_'//trim(filter_char(ifilt))
            else
               filename2=trim(model_run(ipar))//'_wf'//'_'//trim(filter_char(ifilt))
            endif
            write(6,*)misfit(j),iwin,it,'saving waveform kernel...',inside_win(it,iwin)
            if (save_wfkern) call save_misfit_kernel(npts,real(kerphys(it,1:npts)),filename2,appidump1) ! waveform kernel
            if (inside_win(it,iwin)) then
               kernel = kernel + misfit_prefactor(it,j,ifilt) * real(kerphys(it,1:npts))
               !call sum_1d(misfit_prefactor(it,j,ifilt),real(kerphys(it,1:npts)),kernel,npts)
               write(6,*)'Misfit & prefactor:', &
                                misfit(j),norm_fac(iwin,j,ifilt),misfit_prefactor(it,j,ifilt),maxval(kerphys(it,:)),maxval(kernel)

               if (save_n_erase(it,iwin)) then
                  kernel = norm_fac(iwin,j,ifilt)*kernel
                  write(6,*)'save n erase ',misfit(j),norm_fac(iwin,j,ifilt),maxval(kernel)

                  call define_io_appendix(appmywin,ind_win(it,iwin))
                  if (num_filt>0) then 
                     filename2=trim(model_run(ipar))//'_'//trim(misfit(j))//'_'//trim(filter_char(ifilt))
                  else
                     filename2=trim(model_run(ipar))//'_'//trim(misfit(j))             
                  endif
                  write(6,*)'saving time window kernel ',trim(filename2)
                  call save_misfit_kernel(npts,real(kernel)/maxval(real(abs(kernel))),filename2,appmywin)
                  kernel = 0.
               endif ! save n erase

            endif ! inside win
         enddo ! time insidecwindow
      enddo ! loop over time windows

   else
      write(6,*)misfit(j),' is an unknown misfit!!'
      stop

   endif

end subroutine compute_misfit_kernel
!----------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------
subroutine traveltime_delay(tstart,tend,time,u0,v0,udata,n,delay)

integer, intent(in) :: n
real(kind=realkind), intent(in) :: tstart,tend,time(:)
real(kind=realkind), intent(in) :: u0(:),v0(:),udata(:)
real(kind=realkind), intent(out) :: delay
real(kind=realkind) :: denom
integer :: it

delay=0.; denom=0.
do it=1,n
   if (time(it)>=tstart .and. time(it)<=tend) then
      delay = v0(it)*(u0(it)-udata(it)) + delay
      denom= v0(it)**2 + denom
   endif
enddo
delay =delay/denom ! dt cancels out from numerator and denominator

end subroutine traveltime_delay
!----------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------------------
subroutine amplitude_difference(tstart,tend,time,u0,udata,n,ampdiff)

integer, intent(in) :: n
real(kind=realkind), intent(in) :: tstart,tend,time(:)
real(kind=realkind), intent(in) :: u0(:),udata(:)
real(kind=realkind), intent(out) :: ampdiff
real(kind=realkind) :: denom
integer :: it

ampdiff=0.; denom=0.
do it=1,n
   if (time(it)>=tstart .and. time(it)<=tend) then
      ampdiff = u0(it)*(u0(it)-udata(it)) + ampdiff
      denom= u0(it)**2 + denom
   endif
enddo
ampdiff  = ampdiff/denom ! dt cancels out from numerator and denomi

end subroutine amplitude_difference
!----------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------------------
subroutine reconst_model_param_kernel(ipar,ifilt,j)

integer, intent(in) :: ipar,ifilt,j

write(6,*)'NOT DONE YET!!!' 
write(6,*)'.... need to modify the organization of the nested loops by keeping some kernel fields unchanged...'
stop

end subroutine reconst_model_param_kernel
!----------------------------------------------------------------------------------------------------------


!========================
end module posteriori_misfit_operations
!========================


!--------------------------------------------------------------
!dk define_io_appendix
  subroutine define_io_appendix3(app,iproc)
!
! Defines the 4 digit character string appended to any 
! data or io file related to process myid. 
!
  integer :: iproc
  character(len=3) :: app
  character(len=1) :: milp,cenp,dizp,unip

  cenp = char(48+iproc/100)
  dizp = char(48+mod(iproc/10,10))
  unip = char(48+mod(iproc,10))
  
  app = cenp//dizp//unip

  end subroutine define_io_appendix3
!---------------------------------------------------------------

!--------------------------------------------------------------
!dk define_io_appendix
  subroutine define_io_appendix3_real(app,iproc_real)

    use global_parameters

! Defines the 4 digit character string appended to any 
! data or io file related to process myid. 
!
  real(kind=realkind), intent(in) :: iproc_real
  integer :: iproc
  character(len=3) :: app
  character(len=1) :: milp,cenp,dizp,unip

  iproc=floor(iproc_real)

  cenp = char(48+iproc/100)
  dizp = char(48+mod(iproc/10,10))
  unip = char(48+mod(iproc,10))
  
  app = cenp//dizp//unip

  write(6,*)'app real:',iproc_real,iproc,cenp,dizp,unip

  end subroutine define_io_appendix3_real
!---------------------------------------------------------------
