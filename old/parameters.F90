!========================
module parameters
!========================
 
  use global_parameters

  implicit none

  public :: read_kernel_input, initialize_stuff,define_sem_kernel_meshes
  public :: save_simulation_info
  private

  contains 

!-------------------------------------------------------------------------------------------------
subroutine initialize_stuff
  use data_mesh, only : nproc,nproc_mesh,mynum,appmynum,appnproc,ierror

  use mpi
  integer                   :: it

  call MPI_INIT( ierror)
  call MPI_COMM_RANK( MPI_COMM_WORLD, mynum, ierror )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierror )

  if (mynum == 0) then 
     write(6,*)
     write(6,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
     write(6,*)'<><><><><>< COMPUTING KERNELS FROM SEM-WAVEFORMS ><><><><><>'
     write(6,*)'    (c) Tarje Nissen-Meyer & Alexandre Fournier 2007-2010'
     write(6,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
     write(6,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
     write(6,*)
  endif

  call define_io_appendix(appnproc,nproc)
  call define_io_appendix(appmynum,mynum)

  do it=0, nproc-1
     if (nproc>1) call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
     if (mynum == it) write(6,*)'Hello in kernery from proc', it; call flush(6)
  enddo

  if (mynum == 0) write(6,*) 'SEM mesh number of procs:', nproc_mesh

end subroutine initialize_stuff
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
 subroutine read_kernel_input(thr,phr)
   use data_mesh
   use data_arrays, only : kernmesh_type
   use data_arrays, only : do_rho,do_lam,do_mu,do_vp,do_vs,do_imped
   use data_fft, only : filter_what, filter_period,filter_type,filter_period_hi,filter_period_low
   use data_fft,only : omegamin,omegamax,fmin,fmax

   real(kind=realkind), intent(inout) ::  thr,phr
   integer :: i,ierr,idim
   character(len=20) :: stf_type
   logical :: taup_exists
   character(len=4) :: appmywin

   open(unit=12,file='input_kernel.dat',status='old',iostat=ierr)
   read(12,*)calc_or_load_kernels
   read(12,'(a8)')kernmesh_type
   read(12,*)thr
   read(12,*)phr
   read(12,*)phi0
   read(12,*)char_wfkern
   read(12,*)twf1
   read(12,*)twf2
   read(12,*)td_fd
   read(12,*)filter_what
   read(12,*)filter_type
   read(12,*)filter_period_low,filter_period_hi
   fmin=1./filter_period_hi
   fmax=1./filter_period_low
   omegamin=2.*pi/filter_period_hi
   omegamax=2.*pi/filter_period_low
   !   endif
   if (filter_type=='gau' .or. filter_type=='low') filter_period=filter_period_low
   if (filter_type=='hig') filter_period=filter_period_hi
   
   read(12,*)compute_src_kernels
   read(12,*)do_rho
   read(12,*)do_lam
   read(12,*)do_mu
   read(12,*)do_vp
   read(12,*)do_vs
   read(12,*)do_imped
   read(12,*)save_snaps
   read(12,*)times_medium
   read(12,*)norm_kernels
   read(12,*)dump_vtk
   read(12,*)dump_avs
   read(12,*)dump_bin
   read(12,*)dump_ascii
   read(12,*)nsim
   read(12,*)(Mij(i),i=1,6)
   close(12)
   if (mynum==0) write(6,*)'Number of simulations:',nsim

   if ( .not. dump_vtk .and. .not. dump_avs .and. .not. dump_bin .and. .not. dump_ascii) then 
     if (mynum==0) then 
        write(6,*)'Surely want to output *some* kind of kernel?'
        write(6,*)'Check input_kernel.dat....'
        write(6,*)'... but for now changing settings to dump binary...'
     endif
     dump_bin=.true.
  endif

   if (do_vp .or. do_vs .or. do_imped .and. .not. times_medium) then 
      if (mynum==0) write(6,*)'Need medium parameters for wavespeed/impedance!'
      times_medium=.true.
   endif
   
   if (do_vp .and. .not. do_lam) then
      if (mynum==0) write(6,*)'Need bulk kernel for compressional vel. kernel!'
      do_lam=.true.
   endif
 
   if (do_vs .and. .not. ( do_lam .or. do_mu ) ) then 
      if (mynum==0)write(6,*)'need bulk & shear kernels for shear vel. kernel!'
      do_lam=.true.; do_mu=.true.
   endif
   
   if (do_imped .and. .not. (do_lam .or. do_mu .or. do_rho) ) then 
      if (mynum==0) write(6,*)'need all primitive kernels for impedance kernel'
      do_lam=.true.; do_mu=.true.; do_rho=.true.
   endif
   
   thr=thr*pi/180.d0
   phr=phr*pi/180.d0
 
   inquire(file='input_timewindows_taup.dat',exist=taup_exists)
   if (taup_exists) then
      if (mynum==0) write(6,*) 'reading input_timewindows_taup.dat'
      open(unit=13,file='input_timewindows_taup.dat') 
   else
      if (mynum==0) write(6,*) 'reading input_timewindows.dat'
      open(unit=13,file='input_timewindows.dat')
   endif
   read(13,*) ntw ! number of time windows
   if (mynum==0) write(6,*)'Number of time windows:',ntw
   if (ntw==0) then
      if (mynum==0) write(6,*)'No time window? doing waveform kernels only.'
      ntw=1
      allocate(t1(ntw),t2(ntw))
      t1(ntw)=twf1 
      t2(ntw)=twf2
      char_wfkern='allt'
   else
      allocate(t1(ntw),t2(ntw))
      allocate(phase_name(ntw))
      do i=1,ntw
         if (taup_exists) then 
            read(13,*)t1(i),t2(i),phase_name(i)
         else
            read(13,*)t1(i),t2(i)
            call define_io_appendix(appmywin,i)
            phase_name(i)=appmywin
         endif
         if (mynum==0) write(6,*)'time window #,from/to/phase',i,t1(i),t2(i),trim(phase_name(i))

         if (t1(i)<twf1) then 
            if (mynum==0) write(6,*)'time window starts before wf kernel!'
            twf1=t1(i)
         elseif (t2(i)>twf2) then 
            if (mynum==0) write(6,*)'time window extends beyond wf kernel!'
            twf2=t2(i)
         endif
      enddo
   endif
   close(13)

   if (mynum==0) then 
      write(6,*)
      write(6,*)' receiver theta,phi:',thr/pi*180.,phr/pi*180.
      write(6,*)' waveform time window t1,t2:',twf1,twf2
   endif
   call flush(6)

   if (calc_or_load_kernels=='load_wfkernels') then 
      posteriori=.true.
   elseif (calc_or_load_kernels=='calc_wfkernels') then 
      posteriori=.false.
  else 
     write(6,*)'unknown task',calc_or_load_kernels
     stop
  endif

! define component of the receiver-seismogram

!! TNM Nov 2010......... kinda stupid: xforce and yforce, being the 
! same simulations... should be handled jointly, NOT like below.... plus rotations are necessary!
if (src_type1bwd=='monopole') then
   if (src_type1fwd=='monopole') then 
      reccomp=2
   else
      reccomp=3
    endif
elseif (src_type2bwd=='xforce') then 
   reccomp=1
elseif (src_type2bwd=='yforce') then 
      if (src_type1fwd=='monopole') then 
         write(6,*) 'This combination does not make sense:'
         write(6,*)'monopole fwd source will not have any signal on the transverse bwd field...'
         stop
   else
      reccomp=2
    endif
  else 
    reccomp=2
 endif

open(unit=55,file='simulation.info_fwd',status='old')
     read(55,*)junk !trim(bkgrdmodel) ! background model'
     read(55,*)junk !deltat ! time step [s]'
     read(55,*)junk !niter ! number of time steps'
     read(55,*)junk !trim(src_type(1)) ! source type'
     read(55,*)junk !trim(src_type(2)) ! source type'
     read(55,*)junk !trim(stf_type) ! source time function'
     read(55,*)junk !trim(src_file_type) ! source file type'
     read(55,*)junk !period ! dominant source period'
     read(55,*)junk !src_depth/1000. ! source depth [km]'
     read(55,*)magnitude ! scalar source magnitude'
     read(55,*)junk !num_rec_tot ! number of receivers'
     read(55,*)junk !floor(real(niter)/real(seis_it)) ! length of seismogram [time samples]'
     read(55,*)junk !deltat*seis_it ! seismogram sampling [s]'
     read(55,*)junk !correct_azi ! compute seismograms at correct azimuth?'
     read(55,*)junk !floor(real(niter)/real(strain_it)) ! number of strain dumps'
     read(55,*)junk !period/real(strain_samp) ! strain dump sampling rate [s]'
     read(55,*)junk !floor(real(niter)/real(snap_it)) ! number of snapshot dumps'
     read(55,*)junk !deltat*real(snap_it) ! snapshot dump sampling rate [s]'      
     read(55,*)junk !rot_rec ! receiver components '
     read(55,*)ibeg_fwd !   ibeg: beginning gll index for wavefield dumps'
     read(55,*)iend_fwd ! iend: end gll index for wavefield dumps'
     read(55,*)shift_fact_fwd ! source shift factor [s]'
     read(55,*)shift_fact_deltat_fwd ! source shift factor for deltat'
     read(55,*)shift_fact_seis_dt_fwd ! source shift factor for seis_dt'
     read(55,*)shift_fact_deltat_coarse_fwd ! source shift factor for deltat_coarse'
     read(55,*)junk  ! receiver file type
     read(55,*)junk  ! receiver spacing (0 if not even)
     read(55,*)use_netcdf  ! use netcdf for wavefield output?
close(55)
open(unit=55,file='simulation.info_bwd')
     read(55,*)junk !trim(bkgrdmodel) ! background model'
     read(55,*)junk !deltat ! time step [s]'
     read(55,*)junk !niter ! number of time steps'
     read(55,*)junk !trim(src_type(1)) ! source type'
     read(55,*)junk !trim(src_type(2)) ! source type'
     read(55,*)junk !trim(stf_type) ! source time function'
     read(55,*)junk !trim(src_file_type) ! source file type'
     read(55,*)junk !period ! dominant source period'
     read(55,*)junk !src_depth/1000. ! source depth [km]'
     read(55,*)junk !magnitude ! scalar source magnitude'
     read(55,*)junk !num_rec_tot ! number of receivers'
     read(55,*)junk !floor(real(niter)/real(seis_it)) ! length of seismogram [time samples]'
     read(55,*)junk !deltat*seis_it ! seismogram sampling [s]'
     read(55,*)junk !correct_azi ! compute seismograms at correct azimuth?'
     read(55,*)junk !floor(real(niter)/real(strain_it)) ! number of strain dumps'
     read(55,*)junk !period/real(strain_samp) ! strain dump sampling rate [s]'
     read(55,*)junk !floor(real(niter)/real(snap_it)) ! number of snapshot dumps'
     read(55,*)junk !deltat*real(snap_it) ! snapshot dump sampling rate [s]'      
     read(55,*)junk !rot_rec ! receiver components '
     read(55,*)ibeg_bwd !   ibeg: beginning gll index for wavefield dumps'
     read(55,*)iend_bwd ! iend: end gll index for wavefield dumps'
     read(55,*)shift_fact_bwd ! source shift factor [s]'
     read(55,*)shift_fact_deltat_bwd ! source shift factor for deltat'
     read(55,*)shift_fact_seis_dt_bwd ! source shift factor for seis_dt'
     read(55,*)shift_fact_deltat_coarse_bwd ! source shift factor for deltat_coarse'
     read(55,*)junk  ! receiver file type
     read(55,*)junk  ! receiver spacing (0 if not even)
     read(55,*)use_netcdf  ! use netcdf for wavefield output?
close(55)

   if (mynum.eq.0) then
       if ((use_netcdf)) then
#ifdef unc
          write(6,*) 'Using NetCDF input'
#else
          write(6,*) 'trying to use NetCDF input, '
          write(6,*) 'but KERNER is compiled without netcdf support'
          stop 2
#endif
       else
          write(6,*) 'Using binary input'
       end if
   end if

   open(unit=20000,file='sourceparams_fwd.dat')
   read(20000,*)realjunk
   read(20000,*)src_type1fwd
   read(20000,*)src_type2fwd
   read(20000,*)
   read(20000,*)
   read(20000,*)src_depth
   read(20000,*)realjunk
   read(20000,*)realjunk
   read(20000,*)stf_type
   close(20000)

   dim_fwd = 3; if (src_type1fwd=='monopole') dim_fwd = 2; 

   call define_sem_kernel_meshes

   if (nsim>1) then 
      allocate(dir_fwdmesh1(4),srctype(4))
      dir_fwdmesh1(1) = trim(dir_fwdmesh)//'MZZ/'
      dir_fwdmesh1(2) = trim(dir_fwdmesh)//'MXX_P_MYY/'
      dir_fwdmesh1(3) = trim(dir_fwdmesh)//'MXZ_MYZ/'
      dir_fwdmesh1(4) = trim(dir_fwdmesh)//'MXY_MXX_M_MYY/'
      srctype = ['mzz','mxx','mxz','mxy']
   else
      allocate(dir_fwdmesh1(1),srctype(1))
      dir_fwdmesh1(1) = trim(dir_fwdmesh)
      srctype = trim(src_type2fwd)
      Mij = magnitude
   endif

   if (mynum.eq.0) then
       do idim=1,nsim
         write(6,*)'DIRECTORY OF FWD MESH:',trim(dir_fwdmesh1(idim))
       enddo
       write(6,*)'DIRECTORY OF BWD MESH:',trim(dir_bwdmesh)
   end if

   Mij=Mij/magnitude

 end subroutine read_kernel_input
!-------------------------------------------------

!-------------------------------------------------
subroutine define_sem_kernel_meshes
  use data_mesh
  use data_arrays, only : kernmesh_type
  use data_misfits

   open(unit=98,file="input_sem.dat")
   read(98,*)dir_fwdmesh
   read(98,*)dir_bwdmesh
   read(98,*)ext_mesh_name
   read(98,*)data_dir
   close(98)
   lffwd = index(dir_fwdmesh,' ')-1
   lfbwd = index(dir_bwdmesh,' ')-1
   lfext = index(ext_mesh_name,' ')-1

  if (mynum==0) then 
     write(6,*)'forward calc :',dir_fwdmesh(1:lffwd)
     write(6,*)'backward calc:',dir_bwdmesh(1:lfbwd)
     if (kernmesh_type=='ext_mesh') &
          write(6,*)'external mesh:',ext_mesh_name(1:lfext)
     write(6,*)
  endif

end subroutine define_sem_kernel_meshes
!------------------------------------------------

!-------------------------------------------------
subroutine save_simulation_info(thr,phr)

use data_arrays
use data_mesh
implicit none

real(kind=realkind), intent(in) :: thr,phr
character(len=6) :: nnproc
integer :: i,iwin,nnnproc
character(len=20) :: stf_type
real(kind=realkind) :: src_depthbwd

  if (mynum==0) write(6,*)'saving Data/simulation.info ........' 
!af changed all read into write for simulation info
!tnm changed all af-changed write into read for simulation info.

   open(unit=45,file='Data/simulation.info',status='unknown')
!   read(45,*)nnproc,nnnproc
   write(45,15)'nproc=',nproc
   write(45,15)'npts=',npts

12 format(a20,a15,a15)
13 format(a20,1pe11.3)
14 format(a20,a15)
15 format(a20,i10)
16 format(a20,3(a5))
17 format(a20,2(1pe11.3))
18 format(a20,i10,2(1pe11.3))
19 format(a20,i10,1pe11.3)
20 format(a20,l10)
21 format(a20,3(l10))
28 format(a20,i10,2(i10),2(1pe11.3))

   write(45,12)'fwdsrctype= ',src_type1fwd,src_type2fwd
   write(45,13)'fwdsrcdepth=',src_depth
   write(45,14)'fwdstf_type= ',stf_type
   write(45,15)'dim_fwd= ',dim_fwd

   open(unit=20000,file='sourceparams_bwd.dat')
   read(20000,*)realjunk
   read(20000,*)src_type1bwd
   read(20000,*)src_type2bwd
   read(20000,*)
   read(20000,*)
   read(20000,*)src_depthbwd
   read(20000,*)realjunk
   read(20000,*)realjunk
   read(20000,*)stf_type
   close(20000)

   dim_bwd = 3; if (src_type1bwd=='monopole') dim_bwd = 2; 

   write(45,12)'bwdsrctype= ',src_type1bwd,src_type2bwd
   write(45,13)'bwdsrcdepth=',src_depthbwd
   write(45,14)'bwdstf_type= ',stf_type
   write(45,15)'dim_bwd= ',dim_bwd

   write(45,13)'thetar=',thr*180./pi
   write(45,13)'phir=',phr*180./pi
   write(45,20)'save_wfkern=',save_wfkern
   write(45,21)'rho_lam_mu=',do_rho,do_lam,do_mu
   write(45,21)'vp_vs_imped=',do_vp,do_vs,do_imped
   write(45,13)'dt=',dt
   write(45,17)'twf1_twf2=',twf1,twf2
   write(45,15)'nwindows=',nwin
   do i=1,nwin
      write(45,28)'iwin_t1_t2=',i,begdumps(i),enddumps(i),time(begdumps(i)),time(enddumps(i))
   enddo
   write(45,15)'ntw=',ntw
   do i=1,ntw
      write(45,18)'ttwin_t1_t2=',i,t1(i),t2(i)
   enddo
   write(45,15)'nwfkernels=',num_it
   do iwin=1,nwin
      do i=begdumps(iwin),enddumps(iwin)
         write(45,19)'kernid_time=',i,time(i)
      enddo
   enddo
   close(45)
   write(6,*)'saved simulation info into Data/simulation.info'



end subroutine save_simulation_info
!-------------------------------------------------

!========================
end module parameters
!========================
