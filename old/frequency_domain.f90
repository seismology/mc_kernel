!========================
module frequency_domain
!========================

  use global_parameters
  use data_mesh
  use data_arrays
  use fft

  implicit none

  public :: frequency_domain_waveform_kernels
  public :: time_frequency_kernel_norm_init, time_frequency_kernel_norm_mult
  public :: time_frequency_kernel_norm_write_kern
  public :: compute_source_kernels
  public :: compute_deriv_omkern_pk_tt, compute_3deriv_omkern_on_tt
  public :: change_wfkern
  private

  contains 

!----------------------------------------------------------------------------------------
subroutine frequency_domain_waveform_kernels

  use data_rota, only : sinph,cosph


  write(6,*)mynum,'starting frequency domain kerner....'

  allocate(time_mean(1:npts))
  time_mean(1:npts) = 0.0

  if (compute_src_kernels) call compute_source_kernels

  if (do_lam) call wavefields_time2frequency_lam
  if (do_rho) call wavefields_time2frequency_rho
  if (do_mu)  call wavefields_time2frequency_mu

  deallocate(kerspec, kerphys, time_mean)
  deallocate(azim1_fwd, azim2_fwd, azim1_bwd, azim2_bwd, sinph,cosph)
 
end subroutine frequency_domain_waveform_kernels
!----------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
subroutine time_frequency_kernel_norm_init(nom_step)

  use data_arrays, only : ufwd1, ubwd1
  use data_fft, only    : nomega, ntimes

  integer, intent(inout) :: nom_step

  kerphys(1:ntimes,1:npts) = 0.
  kerspec(0:nomega,1:npts) = 0.

  allocate(ufwd1(0:nomega,1:npts,1))
  allocate(ubwd1(0:nomega,1:npts,1))

  open(unit=999, file='Data/lamkern_time_frequency.dat')
  open(unit=998, file='Data/omega.dat')
  open(unit=997, file='Data/time.dat')

  nom_step = 1

end subroutine time_frequency_kernel_norm_init
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine time_frequency_kernel_norm_mult(iomega)

  use data_arrays, only : ufwd1, ubwd1
  use data_fft, only    : omega,nomega, ntimes

  integer, intent(in)   :: iomega

  write(998,*) omega(iomega)
  ufwd1 = 0.
  ufwd1(iomega,1:npts,:) = ugd_src_spec(iomega,1:npts,:)
  ubwd1 = 0.
  ubwd1(iomega,1:npts,:) = ugd_rec_spec(iomega,1:npts,:)
  lamkern_kxt = 0.
  kerphys(1:ntimes,1:npts) = 0.
  kerspec(0:nomega,1:npts) = 0.
  
  call fftb_ker_multiply(lamkern_kxt(1:ntimes,1:npts), ufwd1(0:nomega,1:npts,:), &
                         ubwd1(0:nomega,1:npts,:), 1)

end subroutine time_frequency_kernel_norm_mult
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine time_frequency_kernel_norm_write_kern(iomega, it)

  integer, intent(in) :: iomega, it

  if (iomega==0) write(997,*) time(it)
  write(999,*) sum(lamkernel) / npts, maxval(abs(lamkernel))

end subroutine time_frequency_kernel_norm_write_kern
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine wavefields_time2frequency_lam

  use input_output
  use data_arrays
  use data_fft
  use fft
  use filtering
  use coord_trafo, only : azimuthal_prefactor_fwd_isim
  use mpi
  integer               :: it, i, iomega, isim, iwin
  character(len=4)      :: appiproc, appidump1, appidump2
  real(kind=realkind)   :: mytimewidth
  character(len=200)    :: fname
  real(kind=realkind)   :: globmax
  integer               :: rand_id
  real                  :: real_rand

  write(6,*)mynum,'starting lambda kernel read in frequency domain..'
  ndim = 1 

  if (.not. allocated(kerphys)) allocate(kerphys(1:ntimes,1:npts))
  if (.not. allocated(kerspec)) allocate(kerspec(0:nomega,1:npts))

  ! Read forward field
  kerphys(1:ntimes,1:npts) = 0. 
  do isim=1, nsim
      fname = 'Ekk_fw_'//trim(srctype(isim))
      write(6,*) mynum, 'reading forward field of type: ', srctype(isim)
      call azimuthal_prefactor_fwd_isim(isim, azim1_fwd, azim2_fwd)
      write(6,*) mynum, 'azimuthal fwd prefactor 1 min/max:', minval(azim1_fwd), &
                 maxval(azim1_fwd)
      do it=1, ndumps
          usrc = 0.
          if (mod(it,ndumps/10) == 0 .and. mynum == 0) &
              write(6,'(i2,"0% done in reading lambda fwd field")') 10 * it / ndumps

          call define_io_appendix(appidump1, it)    ! fwd wavefield time index
          call read_summedstraintrace_src(usrc, appidump1, dir_fwdmesh1(isim), it)
          usrc(1:npts) =  azim1_fwd * usrc(1:npts)

          ! write wavefield isim to file Ekk_fw_<srctype>
          if (save_snaps .and. mod(it,20)==0) &
              call save_kernel(npts, usrc(1:npts), fname, appidump1)

          kerphys(it+1,1:npts) = kerphys(it+1,1:npts) + usrc(1:npts)
      enddo
  enddo
  deallocate(usrc)
 
  if (save_snaps .and. nsim>1) then 
     ! write summed wavefields to file Ekkmij
     fname = 'Ekkmij'
     do it=2, ndumps+1
        call define_io_appendix(appidump1, it)   
        call save_kernel(npts, real(kerphys(it,1:npts)), fname, appidump1)   
     enddo
  endif
 
  if (filter_what /= 0) then
     ! demean
     write(6,*) mynum, 'src demean is deactivated...'
     call prepare_kernel(kerphys(1:ntimes,1:npts), time_mean)
  endif
 
  write(6,*) mynum, ' npts = ', npts
   
  ! write some random time series for plotting before fft for debugging
  !do it = 1, 10 
  !  call random_number(real_rand) 
  !  rand_id = int(real_rand * npts) + 1 
  !  write(fname,'(A,I8.8,A)') 'Data/seismo_', rand_id, '.dat' 
  !  write(6,*) fname 
  !  open(unit=4321,file=fname) 
  !  write(4321,*) kerphys(:,rand_id) 
  !  close(4321) 
  !enddo 
  
  ! allocate memory for fourier domain fwd field
  ! XXX This is a bottleneck atm
  write(6,*) mynum, 'src forward fourier transform...'
  write(6,*) mynum, 'allocating memory for fourier domain fwd field'
  
  allocate(ugd_src_spec(0:nomega,1:npts,1))
  ugd_src_spec(0:nomega,1:npts,1) = cmplx(0.,0.)
  call fftf_ker(ugd_src_spec(0:nomega,1:npts,1))
 
  write(6,*) mynum, 'src time wavefield maxmin:', maxval(kerphys), minval(kerphys)
  write(6,*) mynum, 'src frequency wavefield maxmin:', maxval(real(ugd_src_spec)), &
             minval(real(ugd_src_spec))

  ! read receiver wavefield
  write(6,*) 'starting to read rec wavefields from ', trim(dir_bwdmesh)
  kerphys(1:ntimes,1:npts) = 0.
  fname = 'Ekk_bw_'//trim(src_type2bwd)

  do it=1, ndumps
     urec = 0.
     if (mod(it,ndumps/10) == 0 .and. mynum == 0) &
         write(6,"(i2,'0% done in reading lambda bwd field')") 10 * it / ndumps

     call define_io_appendix(appidump1, it) 
     call read_summedstraintrace_rec(urec, appidump1, dir_bwdmesh, it)

     ! write wavefield to file Ekk_bw_<srctype>
     if (save_snaps .and. mod(it,20) == 0) &
         call save_kernel(npts,urec(1:npts), fname, appidump1) 
     
     kerphys(it+1,1:npts) = azim1_bwd * urec(1:npts)
  enddo
  deallocate(urec)
  ! done reading receiver wavefields

  if (filter_what/=0) then
     write(6,*) mynum, 'rec demean is deactivated...'
     !call prepare_kernel(kerphys, time_mean)
     time_mean = 0.
  endif

  write(6,*) mynum, 'rec allocating spectral wavefields...', nomega, npts; call flush(6)
  
  allocate(ugd_rec_spec(0:nomega,1:npts,1))
  ugd_rec_spec(0:nomega,1:npts,1) = cmplx(0.,0.)
  
  write(6,*) mynum, 'rec forward fourier transform...'; call flush(6)
  call fftf_ker(ugd_rec_spec(0:nomega,1:npts,1))
  
  write(6,*) mynum, 'rec time wavefield maxmin:', maxval(kerphys), minval(kerphys)
  write(6,*) mynum, 'rec frequency wavefield maxmin:', maxval(real(ugd_rec_spec)), &
             minval(real(ugd_rec_spec))

  ! starting with wavefield kernel
  allocate(lamkern_kxt(1:ndumps,1:npts))
  lamkern_kxt(:,:) = 0. 
  
  if (td_fd == 'fd') then ! MvD: what if not fd?? TNM: only time will tell...
     call filter_and_mult_kernels(lamkern_kxt, ugd_src_spec, ugd_rec_spec)
     deallocate(ugd_src_spec, ugd_rec_spec)
     fname = 'Kxw_lam_filt'
     if (save_wfkern) call save_omega_kernel(fname, kerspec)
  else
     write(6,*) 'Only ''fd'' is implemented so far'
     stop  
  endif

  write(6,*) mynum, 'max lam time_mean:', maxval(time_mean)

 if (norm_kernels) then 
     globmax_lam = maxval(abs(lamkern_kxt))
     call MPI_ALLREDUCE(globmax_lam, globmax_lam, 1, MPI_REAL, MPI_MAX, &
          MPI_COMM_WORLD, IERROR)
  else 
     globmax_lam = 1.
  endif

  write(6,*) mynum, 'global maximum of lambda Kxt:', globmax_lam

  if (filter_what == 0) time_mean = 0.
  
  call define_io_appendix(appiproc,mynum) 
  open(unit=101,file='Data/maxabs_lam_wfkern_'//appiproc//'.dat')
  do it=1, ndumps
     if (times_medium) then 
        lamkern_kxt(it,1:npts) = -lam * ( lamkern_kxt(it,1:npts) + time_mean )
     else
        lamkern_kxt(it,1:npts) = -1. * ( lamkern_kxt(it,1:npts) + time_mean )
     endif
     write(101,*) time(it), maxval(abs(lamkern_kxt(it,1:npts)))
  enddo
  close(101)

  if (times_medium) deallocate(lam)

  write(6,*) mynum, 'global maximum of final lambda Kxt:', maxval(abs(lamkern_kxt))

  if (save_wfkern) then 
     fname = 'Kxt_lam'
     do iwin=1, nwin ! number of time windows. 1 if waveform kernels are saved.
        do it = begdumps(iwin), enddumps(iwin) 
           call define_io_appendix(appidump1, it) 
           call save_kernel(npts, lamkern_kxt(it,1:npts), fname, appidump1, globmax_lam)
        enddo
     enddo
  endif

  write(6,*) mynum, 'done with frequency-domain lambda kernel computation'

end subroutine wavefields_time2frequency_lam
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine wavefields_time2frequency_rho

  use coord_trafo, only: rotate_velocity_components, rotate_vec_cyl2cart
  use input_output, only : read_velocities_rec, read_velocities_src, save_kernel
  use data_arrays
  use data_fft
  use fft
  use mpi

  integer               :: it,i,idim,isim,iwin
  character(len=4)      :: appiproc,appidump1,appidump2
  real(kind=realkind)   :: mytimewidth
  character(len=20)     :: fname
  real(kind=realkind)   :: globmax
  write(6,*)'starting rho kernel read in frequency domain..'

  ! since we rotate to xyz, we always need to consider 3 components...
  ndim = 3

  ! Read forward field
  allocate(field_phys_dim(1:ntimes,1:npts,1:ndim)) 
  field_phys_dim(1:ntimes,1:npts,1:ndim) = 0. 
  kerphys(1:ntimes,1:npts) = 0. 

  allocate(vrot(1:npts,1:ndim))
  do isim=1, nsim !!!!!!!! NEEEEEEEEED TOOOOOOOOO FIIIIIIIIIIXXXXXXXXXXXXX!!!!!!!!!!!
     do it=1, ndumps
        vsrc = 0.
        vrot = 0.
        if (mod(it,ndumps/10)==0 .and. mynum==0) &
           write(6,"(i2,'0% done in reading rho fwd field')") 10 * it / ndumps
        call define_io_appendix(appidump1,it)    ! fwd wavefield time index
        call read_velocities_src(vsrc, appidump1, dir_fwdmesh1(isim))
        call rotate_vec_cyl2cart(vsrc, vrot, dim_fwd, azim1_fwd, azim2_fwd)
        field_phys_dim(it+1,1:npts,1:ndim) = field_phys_dim(it+1,1:npts,1:ndim) + &
                                             vrot(1:npts,1:ndim)
     enddo
  enddo
  deallocate(vsrc)

  allocate(ugd_src_spec(0:nomega,1:npts,1:ndim))
  ugd_src_spec(0:nomega,1:npts,1:ndim) = cmplx(0.,0.)
  
  do idim=1, ndim
     kerphys = field_phys_dim(:,:,idim)
     call fftf_ker(ugd_src_spec(0:nomega,1:npts,idim))
  enddo

  write(6,*) 'done with fftf_ker'
  write(6,*) 'src kerphys maxmin:', maxval(field_phys_dim), minval(field_phys_dim)
  write(6,*) 'src ugd_src_spec maxmin:', maxval(real(ugd_src_spec)), &
             minval(real(ugd_src_spec))

  ! receiver wavefield
  write(6,*)'starting to read rec wavefields...';call flush(6)
  globmax = 0.
  field_phys_dim(1:ntimes,1:npts,1:ndim) = 0. 
  do it=1, ndumps
     vrec = 0.
     vrot = 0.
     if (mod(it,ndumps/10)==0 .and. mynum==0) &
        write(6,"(i2,'0% done in reading rho bwd field')") 10 * it / ndumps
     if (it+ibwd_shift > ndumps) then
        call define_io_appendix(appidump1, ndumps) 
     else
        call define_io_appendix(appidump1, it+ibwd_shift) 
     endif
     call read_velocities_rec(vrec, appidump1, dir_bwdmesh)
     call rotate_vec_cyl2cart(vrec, vrot, dim_bwd, azim1_bwd, azim2_bwd)
     call rotate_velocity_components(vrot)
     field_phys_dim(it+1,1:npts,1:ndim) = vrot(1:npts,1:ndim) 

     if (maxval(abs(field_phys_dim(it+1,1:npts,1:ndim))) > globmax) &
          globmax = maxval(abs(field_phys_dim(it+1,1:npts,1:ndim)))
     !write(6,*)'dump:',it,maxval(vrec),maxval(field_phys_dim(it+1,1:npts,1:3)),globmax
  end do
  deallocate(vrec,vrot)

  write(6,*)'Forward Fourier transform for each component'
  allocate(ugd_rec_spec(0:nomega,1:npts,1:ndim))
  ugd_rec_spec(0:nomega,1:npts,1:ndim) = cmplx(0.,0.)
  kerphys(1:ntimes,1:npts) = 0.

  do idim=1, ndim
     kerphys = field_phys_dim(:,:,idim)
     call fftf_ker(ugd_rec_spec(0:nomega,1:npts,idim))
  enddo

  write(6,*) 'rec kerphys maxmin:', maxval(field_phys_dim), minval(field_phys_dim)
  write(6,*) 'rec ugd_rec_spec maxmin:', maxval(real(ugd_rec_spec)), &
             minval(real(ugd_rec_spec))
  deallocate(field_phys_dim)

  allocate(rhokern_kxt(1:ntimes,1:npts))
  rhokern_kxt(:,:) = 0. 

  write(6,*) 'filtering rho kerns:', maxval(real(ugd_src_spec)), &
             maxval(real(ugd_rec_spec))
  call filter_and_mult_kernels(rhokern_kxt, ugd_src_spec, ugd_rec_spec)
  deallocate(ugd_rec_spec, ugd_src_spec)
  write(6,*) 'done filtering', maxval(rhokern_kxt)

  fname = 'Kxw_rho_filt'
  call save_omega_kernel(fname, kerspec)
  write(6,*) mynum, 'max rho time_mean:', maxval(time_mean)

  if (norm_kernels) then 
     globmax_rho = maxval(abs(rhokern_kxt))
     call MPI_ALLREDUCE(globmax_rho,globmax_rho,1,MPI_REAL,MPI_MAX, &
          MPI_COMM_WORLD,IERROR)
  else 
     globmax_rho = 1.
  endif
  write(6,*) mynum, 'global maximum of rho Kxt:', globmax_rho

  if (filter_what == 0) time_mean = 0.

  call define_io_appendix(appiproc, mynum) 
  open(unit=101, file='Data/maxabs_lam_wfkern_'//appiproc//'.dat')

  do it=1, ndumps
     if (times_medium) then 
        rhokern_kxt(it,1:npts) = -rho * (rhokern_kxt(it,1:npts) + time_mean)
     else
        rhokern_kxt(it,1:npts) = -1. * (rhokern_kxt(it,1:npts) + time_mean)
     endif
     write(101,*) time(it), maxval(abs(rhokern_kxt(it,1:npts)))
  enddo
  close(101)
  if (times_medium) deallocate(rho)
  write(6,*) mynum, 'global maximum of final rho Kxt:', maxval(abs(rhokern_kxt))

  if (save_wfkern) then 
     fname = 'Kxt_rho'
     do iwin=1, nwin ! number of time windows. 1 if waveform kernels are saved.
        do it = begdumps(iwin), enddumps(iwin) 
           call define_io_appendix(appidump1, it) 
           call save_kernel(npts, rhokern_kxt(it,1:npts), fname, appidump1, globmax_rho)
        enddo
     enddo
  endif
  write(6,*) mynum, 'done with frequency-domain rho kernel computation'

end subroutine wavefields_time2frequency_rho
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine wavefields_time2frequency_mu

  use coord_trafo,  only: rotate_strain_components, rotate_tens_cyl2cart, &
                          azimuthal_prefactor_fwd_isim
  use input_output, only: read_deviatoricstrain_src, read_deviatoricstrain_rec, &
                          save_kernel, save_devstrain
  use data_arrays
  use data_fft
  use fft
  use filtering
  use mpi

  integer               :: it, i, idim, isim, iwin
  character(len=4)      :: appiproc, appidump1, appidump2
  real(kind=realkind)   :: mytimewidth, globmax
  character(len=200)    :: fname1, fname2, fname3, fname4, fname5
  character(len=20)     :: fname
  real(kind=realkind), allocatable :: devstrain_rot(:,:)

  write(6,*)'starting mu kernel read in frequency domain..'
  ndim = 6
  
  ! Read forward fields
  allocate(field_phys_dim(1:ntimes,1:npts,1:ndim)) 
  kerphys(1:ntimes,1:npts) = 0.
  field_phys_dim(1:ntimes,1:npts,1:ndim) = 0. 
  allocate(devstrain(1:npts,dim_fwd*2))
  allocate(devstrain_rot(1:npts,ndim))

  do isim=1, nsim
     fname1 = 'Exz'//trim(srctype(isim))
     call azimuthal_prefactor_fwd_isim(isim, azim1_fwd, azim2_fwd)
     do it=1, ndumps
        devstrain = 0.
        devstrain_rot = 0.
        if (mod(it,ndumps/10)==0 .and. mynum==0) &
            write(6,"(i2,'0% done in reading mu fwd field')") 10 * it / ndumps

        call define_io_appendix(appidump1, it)    ! fwd wavefield time index
        call read_deviatoricstrain_src(devstrain, appidump1, dir_fwdmesh1(isim), it) 
        call rotate_tens_cyl2cart(npts,devstrain, devstrain_rot, azim1_fwd, azim2_fwd, &
                                  dim_fwd*2 )

        if (save_snaps .and. mod(it,20)==0) &
             call save_kernel(npts, devstrain_rot(1:npts,4), fname1, appidump1) ! Exz/Ert

        field_phys_dim(it+1,1:npts,1:ndim) = field_phys_dim(it+1,1:npts,1:ndim) + &
                                             devstrain_rot(1:npts,1:ndim)
     enddo
  enddo
  ! done reading forward fields

  deallocate(devstrain)
  if (save_snaps .and. nsim > 1) then
     fname1 = 'Exzmij'
     do it=2, ndumps+1
        call define_io_appendix(appidump1, it)   
        call save_kernel(npts, field_phys_dim(it,1:npts,4), fname1, appidump1)
     enddo
  endif

  allocate(ugd_src_spec(0:nomega,1:npts,1:ndim))
  ugd_src_spec(0:nomega,1:npts,1:ndim) = cmplx(0.,0.)

  do idim=1, ndim
     kerphys = field_phys_dim(:,:,idim)
     call fftf_ker(ugd_src_spec(0:nomega,1:npts,idim))
  enddo
  fname = 'Kxw_fwd_Exz'
  if (save_wfkern .and. save_snaps) &
     call save_omega_kernel(fname,ugd_src_spec(0:nomega,1:npts,4))

  ! receiver wavefields
  globmax = 0.
  write(6,*)'starting to read rec wavefields...'; call flush(6)

  allocate(devstrain(1:npts,dim_bwd*2))
  kerphys(1:ntimes,1:npts) = 0.0
  field_phys_dim(1:ntimes,1:npts,1:ndim) = 0. 

  fname1 = 'ExzPz0'
  fname2 = 'ExzPzr'
  fname3 = 'EszPz'
  fname4 = 'P0'
  fname5 = 'Pr'

  do it=1, ndumps
     devstrain = 0.
     devstrain_rot = 0.
     if (mod(it,ndumps/10) == 0 .and. mynum == 0) &
        write(6,"(i2,'0% done in reading mu bwd field')") 10 * it / ndumps
     call define_io_appendix(appidump1, it)
     call read_deviatoricstrain_rec(devstrain,appidump1,dir_bwdmesh,it)

     if (save_snaps .and. mod(it,20) == 0) &
        call save_kernel(npts, devstrain(1:npts,4), fname3, appidump1) ! Esz/Ert
     call rotate_tens_cyl2cart(npts, devstrain, devstrain_rot, azim1_bwd, azim2_bwd, &
                               dim_bwd*2 )
     !if (save_snaps .and. mod(it,20) == 0) &
     !     call save_kernel(npts, devstrain_rot(1:npts,4), fname1, appidump1) ! Exz/Ert
     if (save_snaps .and. mod(it,20) == 0) &
        call save_devstrain(npts, devstrain_rot, fname4, 6, appidump1)

     call rotate_strain_components(npts, devstrain_rot)
     if (save_snaps .and. mod(it,20) == 0) &
        call save_devstrain(npts, devstrain_rot, fname5, 6, appidump1)
     !if (save_snaps .and. mod(it,20) == 0) &
     !     call save_kernel(npts, devstrain_rot(1:npts,4), fname2, appidump1) 
                                                                    ! Exz/Ert rotated
     field_phys_dim(it+1,1:npts,1:ndim) = magnitude * devstrain_rot(1:npts,1:ndim) 
                                                             !MAGNITUDE SHOULDNT BE HERE

     if (maxval(abs(devstrain_rot(1:npts,1:ndim))) > globmax) &
          globmax = maxval(abs(devstrain_rot(1:npts,1:ndim)))
  end do
  ! end receiver wavefields

  deallocate(devstrain_rot, devstrain)
  if (filter_what /= 0) then
     write(6,*) mynum,'rec compute & subtract mean...'
     do idim=1, ndim
        call prepare_kernel_real(field_phys_dim(:,:,idim), time_mean)
     enddo
  endif

  write(6,*)mynum,'Forward Fourier transform for receiver-side deviatoric wavefields...'
  allocate(ugd_rec_spec(0:nomega,1:npts,1:ndim))
  ugd_rec_spec(0:nomega,1:npts,1:ndim) = cmplx(0.,0.)
  
  do idim=1, ndim
     kerphys = field_phys_dim(:,:,idim)
     call fftf_ker(ugd_rec_spec(0:nomega,1:npts,idim))
  enddo
  deallocate(field_phys_dim)
 
  fname='Kxw_bwd_Exz'
  if (save_wfkern .and. save_snaps) &
     call save_omega_kernel(fname, ugd_rec_spec(0:nomega,1:npts,4))

  allocate(mukern_kxt(1:ndumps,1:npts))
  mukern_kxt(:,:) = 0. 

  if (td_fd == 'fd') then
     call filter_and_mult_kernels(mukern_kxt, ugd_src_spec, ugd_rec_spec)
     deallocate(ugd_rec_spec, ugd_src_spec)
     fname = 'Kxw_mu_filt'
     if (save_wfkern .and. save_snaps) call save_omega_kernel(fname, kerspec)
  endif
  write(6,*) mynum, 'max mu time_mean:', maxval(time_mean)

  if (norm_kernels) then 
     globmax_mu = maxval(abs(mukern_kxt))
     call MPI_ALLREDUCE(globmax_mu, globmax_mu, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, &
                        IERROR)
  else 
     globmax_mu = 1.
  endif
  write(6,*) mynum, 'global maximum of mu Kxt:', globmax_mu

  if (filter_what==0) time_mean = 0.

  call define_io_appendix(appiproc, mynum) 
  open(unit=101, file='Data/maxabs_mu_wfkern_'//appiproc//'.dat')

  do it=1, ndumps
     if (times_medium) then 
        mukern_kxt(it,1:npts) = -2. * mu * ( mukern_kxt(it,1:npts) + time_mean)
     else
        mukern_kxt(it,1:npts) = -2. * ( mukern_kxt(it,1:npts) + time_mean)
     endif
     write(101,*) time(it), maxval(abs(mukern_kxt(it,1:npts)))
  enddo
  if (times_medium) deallocate(mu)
  close(101)

  write(6,*) mynum, 'global maximum of final mu Kxt:', maxval(abs(mukern_kxt))
  if (save_wfkern) then 
     fname = 'Kxt_mu'
     do iwin=1, nwin ! number of time windows. 1 if waveform kernels are saved.
        do it=begdumps(iwin), enddumps(iwin) 
           call define_io_appendix(appidump1, it) 
           call save_kernel(npts, mukern_kxt(it,1:npts), fname, appidump1, globmax_mu)
        enddo
     enddo
  endif
  write(6,*) mynum, 'done with frequency-domain mu kernel computation'

end subroutine wavefields_time2frequency_mu
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine filter_and_mult_kernels(kern_all, ufwd_spec, ubwd_spec)
  use data_arrays
  use data_fft
  use fft
  use filtering, only : filter_field,time_shift_field
  use input_output, only : save_kernel
  
  implicit none

  real(kind=4), dimension(1:ndumps,1:npts), intent(out) :: kern_all
  complex(kind=8), dimension(0:nomega,1:npts,1:ndim), intent(inout) :: ufwd_spec, ubwd_spec

  integer                       :: idim, npts_om, iomega, i, it
  complex(kind=8), allocatable  :: timeshift_fourier_fwd(:), timeshift_fourier_bwd(:)
  character(len=20)             :: fname
  character(len=4)              :: appidump1

  npts_om = npts * (nomega + 1)
  kerspec(0:nomega,1:npts) = 0.

  ! exact time shift
  write(6,*) mynum, 'max ufwd,ubwd spec before time shift:', &
             maxval(abs(real(ufwd_spec))), maxval(abs(real(ubwd_spec)))
  call time_shift_field(ufwd_spec, ndim,shift_fact_fwd)
  call time_shift_field(ubwd_spec, ndim,shift_fact_bwd)
  write(6,*) mynum, 'max ufwd,ubwd spec after time shift:', &
             maxval(abs(real(ufwd_spec))), maxval(abs(real(ubwd_spec)))

  write(6,*) 'Filtering and backward Fourier transform'

  if (filter_what == 0) then
     write(6,*) mynum, 'NOT filtering'
     call fftb_ker_multiply(kern_all(1:ndumps,1:npts), ufwd_spec(0:nomega,1:npts,:), &
                            ubwd_spec(0:nomega,1:npts,:), ndim)

  elseif (filter_what == 1) then
     write(6,*) mynum, 'filtering fwd field'
     do idim=1, ndim
        call filter_field(ufwd_spec(0:nomega,1:npts,idim))
        kerspec(0:nomega,1:npts) = kerspec(0:nomega,1:npts) + &
                                   ufwd_spec(0:nomega,1:npts,idim) * &
                                   ubwd_spec(0:nomega,1:npts,idim)
     enddo
     call fftb_ker(kern_all(1:ndumps,1:npts), kerspec(0:nomega,1:npts))
     call prepare_kernel_real(kern_all(1:ndumps,1:npts), time_mean(1:npts))

  elseif (filter_what == 2) then 
     write(6,*) mynum, 'filtering bwd field'
     do idim=1, ndim
        call filter_field(ubwd_spec(0:nomega,1:npts,idim))
        kerspec(0:nomega,1:npts) = kerspec(0:nomega,1:npts) + &
                                   ufwd_spec(0:nomega,1:npts,idim) * &
                                   ubwd_spec(0:nomega,1:npts,idim)
     enddo
     call fftb_ker(kern_all(1:ndumps,1:npts), kerspec(0:nomega,1:npts))
     call prepare_kernel_real(kern_all(1:ndumps,1:npts), time_mean(1:npts))

  elseif (filter_what == 3) then
      write(6,*)mynum,'filtering fwd & bwd field'
      do idim=1, ndim
         call filter_field(ufwd_spec(0:nomega,1:npts,idim))
         call filter_field(ubwd_spec(0:nomega,1:npts,idim))
         kerspec(0:nomega,1:npts) = kerspec(0:nomega,1:npts) + &
                                    ufwd_spec(0:nomega,1:npts,idim) * &
                                    ubwd_spec(0:nomega,1:npts,idim)

      enddo
      call fftb_ker(kern_all(1:ndumps,1:npts), kerspec(0:nomega,1:npts))
      call prepare_kernel_real(kern_all(1:ndumps,1:npts), time_mean(1:npts))

!!$      if (save_wfkern) then 
!!$         if (ndim == 1) then 
!!$            fname = 'Kxw_fwd_Exz_fil'
!!$            call save_omega_kernel(fname, ufwd_spec(0:nomega,1:npts,1))
!!$            fname = 'Kxw_bwd_Exz_fil'
!!$            call save_omega_kernel(fname, ubwd_spec(0:nomega,1:npts,1))
!!$         elseif (ndim > 3) then
!!$            fname = 'Kxw_fwd_Exz_fil'
!!$            call save_omega_kernel(fname, ufwd_spec(0:nomega,1:npts,4))
!!$            fname = 'Kxw_bwd_Exz_fil'
!!$            call save_omega_kernel(fname, ubwd_spec(0:nomega,1:npts,4))
!!$         endif
!!$
!!$         call fftb_ker(kern_all(1:ndumps,1:npts), ufwd_spec(0:nomega,1:npts,ndim))
!!$         fname = 'fwdfil'
!!$         do it=1, ndumps
!!$            call define_io_appendix(appidump1,it)   
!!$            if (mod(it,20)==0 ) call save_kernel(npts,kern_all(it,1:npts),fname,appidump1)
!!$         enddo
!!$         
!!$         call fftb_ker(kern_all(1:ndumps,1:npts),ubwd_spec(0:nomega,1:npts,ndim))
!!$         fname = 'bwdfil'
!!$         do it=1, ndumps
!!$            call define_io_appendix(appidump1,it)   
!!$            if (mod(it,20)==0 ) call save_kernel(npts,kern_all(it,1:npts),fname,appidump1)
!!$         enddo
!!$      endif
         
  elseif (filter_what == 4) then
     write(6,*) mynum, 'filtering waveform kernel'
     do idim=1, ndim
        kerspec(0:nomega,1:npts) = kerspec(0:nomega,1:npts) + &
                                   ufwd_spec(0:nomega,1:npts,idim) * &
                                   ubwd_spec(0:nomega,1:npts,idim)
     enddo

     fname = 'Kxw_unfil'
     if (save_wfkern) call save_omega_kernel(fname, kerspec)
     call filter_field(kerspec)
     call fftb_ker(kern_all(1:ndumps,1:npts), kerspec(0:nomega,1:npts))
     call prepare_kernel_real(kern_all(1:ndumps,1:npts), time_mean(1:npts))
     
  else
     write(6,*) 'unknown filter option!', filter_what
     stop 2
  endif

  write(6,*) mynum, 'Filtering done. Max values fwd, bwd:', &
             maxval(abs(real(ufwd_spec))), maxval(abs(real(ubwd_spec)))
  write(6,*) mynum, 'Filtering done. Max values kern om time:', &
             maxval(abs(real(kerspec))), maxval(abs(kern_all))

end subroutine filter_and_mult_kernels
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_source_kernels

  use coord_trafo,  only: rotate_strain_components, rotate_tens_cyl2cart
  use input_output, only: read_strain_src
  use data_arrays
  use data_fft
  use fft

  integer               :: it, i, idim, iomega
  character(len=4)      :: appiproc, appidump1, appidump2
  real(kind=realkind)   :: mytimewidth, t
  character(len=200)    :: fname
  real(kind=realkind)   :: globmax
  real(kind=realkind), allocatable  :: strain_src(:,:), strain_src_rot(:,:), stf(:), &
                                       uphys(:)
  complex(kind=8), allocatable      :: stf_spec(:), field_spec(:,:,:), uspec(:)

  write(6,*)mynum,'starting src kernel read in frequency domain..'

  ndim = 6

  ! read strain wavefield and retrieve only time series at source element
  globmax = 0.
  write(6,*) 'starting to read strain tensor wavefields for src kernels...', n_src
  call flush(6)

  allocate(field_phys_dim(1:ntimes,1:n_src,1:ndim)) 
  allocate(strain_src(1:n_src,dim_bwd*2))
  allocate(strain_src_rot(1:n_src,ndim))

  strain_src_rot = 0.
  kerphys(1:ntimes,1:n_src) = 0.0
  field_phys_dim(1:ntimes,1:n_src,1:ndim) = 0. 
  do it=1, ndumps
     strain_src = 0.
     if (mod(it,ndumps/10) == 0 .and. mynum == 0) &
        write(6,"(i2,'0% done in reading bwd strain tensor for src kernels')") &
            10 * it / ndumps

     if (it+ibwd_shift > ndumps) then
        call define_io_appendix(appidump1, ndumps) 
     else
        call define_io_appendix(appidump1, it+ibwd_shift) 
     endif
     !write(6,*)'reading deviatoric strain...';call flush(6)
     call read_strain_src(strain_src, appidump1, dir_bwdmesh)
     !write(6,*)'rotating tens cyl...';call flush(6)
     call rotate_tens_cyl2cart(n_src, strain_src, strain_src_rot, azim1_src, azim2_src, &
                               dim_bwd*2 )
     !write(6,*)'rotating strain ...';call flush(6)
     call rotate_strain_components(n_src, strain_src_rot)
     field_phys_dim(it+1,1:n_src,1:ndim) = strain_src_rot(1:n_src,1:ndim)
     if (maxval(abs(strain_src_rot(1:n_src,1:ndim))) > globmax) &
          globmax = maxval(abs(strain_src_rot(1:n_src,1:ndim)))
  end do
  deallocate(strain_src,strain_src_rot)

  ! write strain elements to file
  open(unit=101, file='Data/src_wfkern_strain_el.dat')
  do i=1, n_src
     do it=1, ndumps
        write(101,20) time(it), (field_phys_dim(it,i,idim), idim=1, ndim)
     enddo
  enddo
  close(101)

  write(6,*) 'construct source time function'
  open(unit=101, file='Data/stf.dat')
  allocate(stf(ntimes), stf_spec(0:nomega))
  stf = 0.
  stf_spec = cmplx(0.,0.)

  do i=1, ntimes
     t = dble(i) * dt
     stf(i) = dexp(-( (3.5d0 / filter_period * (t - 1.5d0 * filter_period))**2 ))
     write(101,*) t, stf(i)
  enddo
  close(101)
  call fftf1d(stf(1:ndumps), stf_spec(0:nomega))

  open(unit=101, file='Data/stf_spec.dat')
  do iomega=0, nomega
     write(101,*) omega(iomega), real(stf_spec(iomega))
  enddo
  close(101)

  write(6,*) 'compute omega kernel'
  allocate(field_spec(0:nomega,n_src,ndim))
  field_spec = cmplx(0.,0.)
  do idim=1, ndim
     do i=1, n_src
        call fftf1d(field_phys_dim(1:ndumps,i,idim), field_spec(0:nomega,i,idim))
     enddo
  enddo

  write(6,*) 'write omega kernel to file'
  open(unit=101, file='Data/src_wfkern_strain_el_spec.dat')
  do i=1, n_src
     do iomega=0, nomega
        write(101,20) omega(iomega), (real(field_spec(iomega,i,idim)), idim=1, ndim)
     enddo
  enddo
  close(101)

  write(6,*) 'convolve with source time function'
  do idim=1, ndim
     do i=1, n_src
        field_spec(0:nomega,i,idim) = field_spec(0:nomega,i,idim) * stf_spec
        call fftb1d(field_spec(0:nomega,i,idim), field_phys_dim(:,i,idim))
     enddo
  enddo

  write(6,*) 'write omega kernel convolved with stf to file'
  open(unit=101, file='Data/src_wfkern_strain_el_convstf_spec.dat')
  do i=1, n_src
     do iomega=0, nomega
        write(101,20) omega(iomega), (real(field_spec(iomega,i,idim)), idim=1, ndim)
     enddo
  enddo
  close(101)

  write(6,*) 'write strain convolved with stf to file'
  open(unit=101, file='Data/src_wfkern_strain_el_convstf.dat')
  do i=1, n_src
     do it=1, ndumps
        write(101,20) time(it), (field_phys_dim(it,i,idim), idim=1, ndim)
     enddo
  enddo
  close(101)

20 format(7(1pe11.3))

  !write(6,*)'rec kerphys maxmin:',maxval(field_phys_dim),minval(field_phys_dim)
  !write(6,*)'rec ugd_rec_spec maxmin:',maxval(real(ugd_rec_spec)),minval(real(ugd_rec_spec))
  !deallocate(field_phys_dim)
  
  !allocate(srckern_all_phys(1:ntimes,1:n_src)); srckern_all_phys(:,:)= 0. 
  !ugd_src_spec=1.
  !call filter_and_mult_kernels(srckern_all_phys,ugd_src_spec,ugd_rec_spec)
  deallocate(field_spec,stf,stf_spec)

  !if (filter_what/=0) then
  !   do it=1, ntimes
  !      srckern_all_phys(it,1:n_src) = srckern_all_phys(it,1:n_src) + time_mean(1:n_src)
  !   end do
  !   deallocate(time_mean)
  !endif
  ! azimuth factor should be here somewhere.... at least once F_x,y are considered!!!!!!

  write(6,*) mynum, 'done with frequency-domain source kernel computation'  

  stop

end subroutine compute_source_kernels
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine save_omega_kernel(kerntype, wspec)

  use input_output
  use data_arrays
  use data_fft
  use mpi

  complex(kind=8), dimension(:,:), intent(in)   :: wspec
  character(len=20), intent(in)                 :: kerntype
  
  integer               :: iomega
  character(len=4)      :: appidump1
  character(len=200)    :: fname
  real                  :: globmax

  if (save_wfkern) then
     write(6,*) mynum, 'saving omega wavefield/kernel....'
     globmax = 0.
     do iomega=5, nomega
        if (maxval(abs(real(wspec(iomega,1:npts)))) > globmax) &
           globmax = maxval(abs(real(wspec(iomega,1:npts))))
     enddo
     !write(6,*)mynum,'globmax kerspec local',globmax

     call MPI_ALLREDUCE(globmax, globmax, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, IERROR)  

     if (globmax == 0.) &
        globmax = 1.

     !write(6,*)mynum,'globmax ',globmax
     !open(unit=101,file='Data/omega_in_om_period.dat')
     fname = trim(kerntype)
     do iomega=5, nomega
        !write(101,*)iomega,omega(iomega),2.*pi/omega(iomega)
        call define_io_appendix(appidump1, iomega)
        call save_kernel(npts, sngl(real(wspec(iomega,:))), fname, appidump1)
     enddo
  endif
  !close(101)

end subroutine save_omega_kernel
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_deriv_omkern_pk_tt(omkern_in,omkern_deriv)
  use data_fft,only : omega, nomega

  complex(kind=8), dimension(0:nomega,1:npts), intent(in)   :: omkern_in
  complex(kind=8), dimension(0:nomega,1:npts), intent(out)  :: omkern_deriv
  
  integer :: i

  write(6,*)'computing first derivative of spectral kernel...'

  do i=1, nomega
     omkern_deriv(i,:) = cmplx(0.,1.) * omega(i) * omkern_in(i,:)
  enddo 

end subroutine compute_deriv_omkern_pk_tt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_3deriv_omkern_on_tt(omkern_in,omkern_deriv)
  use data_fft,only : omega,nomega

  complex(kind=8), dimension(0:nomega,1:npts), intent(in)   :: omkern_in
  complex(kind=8), dimension(0:nomega,1:npts), intent(out)  :: omkern_deriv

  integer :: i
  
  write(6,*) 'computing third derivative of spectral kernel...'
  
  do i=1, nomega
     omkern_deriv(i,:) = -cmplx(0.,1.) * omega(i)**3 * omkern_in(i,:)
  enddo 

end subroutine compute_3deriv_omkern_on_tt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!ipgp october 2011
subroutine change_wfkern(wfkerntype)
  !this is needed if the misfit kernel does not depend on the raw K(x,t),
  ! but e.g. a time derivative or Gabor transform of it
  
  character(len=3), intent(in) :: wfkerntype
  
  ! choose correct kernel type
  if (wfkerntype == '1st') then
     call compute_deriv_omkern_pk_tt(kerspec, kerspec)
  elseif(wfkerntype == '3rd') then
     call compute_3deriv_omkern_on_tt(kerspec, kerspec)
  elseif(wfkerntype == 'gab') then
     write(6,*) 'not done yet....'
     stop 2
  endif

end subroutine change_wfkern

!========================
end module frequency_domain
!========================
