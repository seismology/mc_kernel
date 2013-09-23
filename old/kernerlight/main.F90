!===================
program kerner
!===================


  !use data_mesh
  !use data_arrays
  use commpi
  use global_parameters
  !use data_posteriori
  !use data_fft, only : nomega, omega
  !use data_misfits 

  use parameters
  use mesh_types
  use parameter_types

  !use kernel_meshes
  !use misfit_model_param, only : init_time_window_misfit_medium, misfit_kernels
  !use misfit_model_param, only : allocate_arrays, waveform_kernels
  !use coord_trafo, only: init_rotations
  !use frequency_domain
!  use posteriori_misfit_operations
  !use fft, only : fftb_ker_dble,init_fft
  !use input_output 
  use nc_routines, ONLY: open_netcdf, close_netcdf
  use kernel_meshes, ONLY: init_read_decomp_save_rot_kernel_mesh


  implicit none

  type(kernelmeshtype) :: mesh
  type(paramtype)      :: params
  type(modeltype)      :: bgmodel
  integer                           :: it, ifilt, iwin, icount, j, ipar, nom_step, iom
  integer                           :: num_thr, i
  real(kind=realkind)               :: thr, phr
  real(kind=realkind),allocatable   :: thr_arr(:)
  character(len=4)                  :: appidump1

  ! initialize MPI
  call initialize_mpi
  
  ! Read directories of forward and backward wavefield
  call define_sem_kernel_meshes(params%io, params%src_fwd%nsim, params%src_bwd%nsim)

  ! Open NetCDF files for reading
  call open_netcdf(params%io)

  ! Read parameters from input_kernel and parameters of forward runs
  call read_parameters(params)
  
  write(6,*)'Computing waveform kernels.....'

  ! initialize, read, and decompose kernel mesh
  call init_read_decomp_save_rot_kernel_mesh(params, mesh, bgmodel)

!!!! ALL GLOBAL VARIABLES REPLACED UNTIL HERE
  
!  ! initialize everything related to frequency-domain operations
!  if (td_fd=='fd' ) call init_fft(time, thr, phr)
!
!  ! initialize everything related to measurement, misfit, model param. (X-Factor)
!  call init_time_window_misfit_medium(thr, phr)
!
!  ! information to store for later use (e.g. when loading waveform kernels only)
!  call save_simulation_info(thr, phr)
!
!  ! compute all a priori rotations
!  call init_rotations(thr, phr)
!
!  ! frequency-domain computation of convolutions
!  if (td_fd=='fd' .or. td_fd=='tf') then 
!      !! @todo: Add domain separation for very large meshes, which exceed system memory, here
!      call define_io_appendix(appmynum, mynum )
!      call frequency_domain_waveform_kernels
!      call allocate_arrays
!      if (mynum == 0) write(6,*) 'Starting time loop for # kernels:', num_it
!      icount = 0
!      nom_step = nomega + 1
!      if (td_fd=='tf') call time_frequency_kernel_norm_init(nom_step)
!   
!      do iom=0, nomega, nom_step ! This is just a dummy loop for regular
!                                 ! frequ-dom. kernels
!         if (td_fd=='tf') call time_frequency_kernel_norm_mult(iom)
!   
!         do iwin = 1, nwin ! number of time windows. 1 if waveform kernels are saved.
!            do it=begdumps(iwin), enddumps(iwin) 
!               call define_io_appendix(appidump1, it) 
!               icount = icount + 1
!               call waveform_kernels(appidump1, time(it), it)              
!               if (td_fd=='tf') call time_frequency_kernel_norm_write_kern(iom, it)
!               if (inside_win(it,iwin)) call misfit_kernels(it, iwin, appidump1)
!               call runtime_info(it, iwin, icount)
!            enddo
!         enddo
!      enddo
!
!      if (td_fd=='tf') then
!          close(999)
!          close(998)
!          close(997)
!      endif
!
!  endif
!  
!#ifdef unc
!  if (use_netcdf) call close_netcdf
!#endif

  call pend 

  write(6,*) mynum, '=========PROGRAM kerner FINISHED============='; call flush(6)

!=======================
end program kerner
!=======================


!-------------------------------------------------------------------------------------------------
subroutine initialize_mpi

  use commpi
  integer                   :: iproc

  call pinit

  if (mynum == 0) then 
     write(6,*)
     write(6,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
     write(6,*)'<><><><><>< COMPUTING KERNELS FROM SEM-WAVEFORMS ><><><><><>'
     write(6,*)'    (c) Tarje Nissen-Meyer & Alexandre Fournier 2007-2010'
     write(6,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
     write(6,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
     write(6,*)
  endif

  call define_io_appendix(appnproc, nproc)
  call define_io_appendix(appmynum, mynum)

  do iproc=0, nproc-1
     if (nproc>1) call barrier 
     if (mynum == iproc) then
         write(6,*)'Hello in kernery from proc', iproc
         call flush(6)
     end if
  enddo

end subroutine initialize_mpi
!-------------------------------------------------------------------------------------------------

!------------------------------------------------------------
!subroutine runtime_info(it, iwin, icount)
!
!  use data_mesh
!  use data_arrays
!
!  integer, intent(in) :: it,iwin
!
!       if (mynum==0) write(6,88) iwin, it, time(it), real(icount)/real(num_it)*100.
!       !if (save_wfkern .and. mynum==0)  write(6,*)' ...saving waveform kernel at t=',time(it)
!
!       if (inside_win(it,iwin) .and. mynum==nproc-1) &
!            write(6,13) t1(ind_win(it,iwin)), time(it), t2(ind_win(it,iwin))
!       !if (mynum==0) write(6,12) time(it),maxval(abs(usrc)),maxval(abs(urec))
!
!       if (mynum==nproc-1 .and. do_lam) write(6,121) time(it), maxval(abs(lamkern_kxt(it,:)))
!       if (mynum==nproc-1 .and. do_mu) write(6,122) time(it), maxval(abs(mukern_kxt(it,:)))
!       if (mynum==0 .and. do_rho) write(6,123) time(it), maxval(abs(rhokern_kxt(it,:)))
!
!88 format('<<Runtime info>> time window:',i2,' , it=',i4, &
!          ', t=',f6.1,'s, i.e.',f5.1,'% done.')
!13 format(' ...updating traveltime kernel:',f8.2,' < ',f8.2,' < ',f8.2)
!12 format('<<Runtime info>> t=',f6.1,'s, ufwd:',1pe9.2,' ,ubwd=',1pe10.2)
!121 format('<<Runtime info>> t=',f6.1,'s, lam kernel:',1pe10.2)
!122 format('<<Runtime info>> t=',f6.1,'s, mu kernel:',1pe10.2)
!123 format('<<Runtime info>> t=',f6.1,'s, ro kernel:',1pe10.2)
!
!end subroutine runtime_info
!--------------------------------------------------------------

!--------------------------------------------------------------
!dk define_io_appendix
subroutine define_io_appendix(app,iproc)
!
! Defines the 4 digit character string appended to any 
! data or io file related to process myid. 
!
  integer :: iproc
  character(len=4) :: app
  character(len=1) :: milp,cenp,dizp,unip

  milp = char(48+    iproc/1000)
  cenp = char(48+mod(iproc/100,10))
  dizp = char(48+mod(iproc/10,10))
  unip = char(48+mod(iproc,10))
  
  app = milp//cenp//dizp//unip

end subroutine define_io_appendix
!---------------------------------------------------------------,
!
!dk flush-------------------------------------------------------
! subroutine flush(iunit)
!integer :: iunit
!	Pseudo flush routine, in case flush is not supported
! by compiler. Comment out/Remove if flush supported by compiler.

!end subroutine flush
!---------------------------------------------------------------
