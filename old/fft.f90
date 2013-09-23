!============
module fft
!============
  use global_parameters, only : realkind
  use data_fft
  use data_mesh, only : mynum, ierror, nproc, npts, u0, v0
  use data_arrays, only : kerphys, kerspec
  implicit none

  public :: init_fft, exit_fft
  public :: fftf1d, fftb1d, fftf1d_dble, fftb1d_dble
  public :: fftf_ker, fftb_ker, fftb_ker_multiply, fftb_ker_dble
  public :: prepare_kernel, prepare_kernel_realkind, prepare_kernel_real
  public :: prepare_fwd_field
  public :: prepare_bwd_field
  private

  contains

!-----------------------------------------------------------------------------------------
subroutine init_fft(time,thr,phr)
  use filtering,  only : define_filter
  use data_arrays,  only: uzrec, vzrec
  use global_parameters
  use,  intrinsic :: iso_c_binding
  
  include 'fftw3.f03'
  include 'mesh_params_kernel.h'
  
  real(kind=realkind),  intent(in) :: time(ndumps), thr, phr
  integer               :: itimes, iomega, it, nextpow2
  integer               :: rank, n, howmany, inembed, istride, idist
  integer               :: onembed, ostride, odist
  real(kind=realkind)   :: dt
  integer, parameter    :: fftw_mode = FFTW_ESTIMATE !  quicker for testing (see
                                            !manual page 25), else use FFTW_PATIENT

  ! af make sure to include t=0 as well in the timeseries 

  if (mynum==0) write(6,*) 'Initializing FFT related things...'

  ! calculate next power of 2 for most efficient fft zero padding
  nextpow2 = 2
  do while (nextpow2 < ndumps) 
      nextpow2 = nextpow2 * 2
  end do
  if (mynum==0) write(6,*) 'nextpow2: ', nextpow2
  npad = nextpow2 * 2 - ndumps - 1

  !npad = ndumps! MvD removed /2 to make it correct ! TNM added /2 to make it cheaper...
  ntimes = ndumps + 1 + npad
  write(6,*) 'number of time samples without padding: ', ndumps
  write(6,*) 'number of time samples after padding: ', ntimes
  nomega = ntimes / 2  ! MvD from the definition provided in the fftw manual
                       ! (n.b. all arrays over omega are defined 0:nomega, so
                       ! they have dimension ntimes / 2 + 1)
  write(6,*) 'number of freq samples: ', nomega

  if (mynum==0) write(6,*) 'number of frequencies:', nomega
  if (mynum==0) write(6,*) 'number of points:', npts

  allocate(uphys1d(1:ntimes), uspec1d(0:nomega))
  call dfftw_plan_dft_r2c_1d(plan_fftf1d, ntimes, uphys1d, uspec1d, fftw_mode) 
  call dfftw_plan_dft_c2r_1d(plan_fftb1d, ntimes, uspec1d, uphys1d, fftw_mode)

  allocate(kerphys(1:ntimes,1:npts)) 
  kerphys(1:ntimes,1:npts) = 0. 
  
  allocate(kerspec(0:nomega,1:npts))
  kerspec(0:nomega,1:npts) = cmplx(0.,0.) 

  ! Forward plan for full kernel
  rank = 1
  n = ntimes
  howmany = npts
  inembed = npts
  idist = ntimes
  istride = 1
  odist = nomega + 1
  ostride = 1
  onembed = npts
      
  call dfftw_plan_many_dft_r2c(plan_fftf_ker, rank,n,howmany,kerphys, inembed, istride, &
                               idist, kerspec, onembed, ostride, odist, fftw_mode)

  ! Backward plan for full kernel
  rank = 1
  n = ntimes
  howmany = npts
  idist = nomega + 1
  istride = 1
  inembed = npts
  onembed = npts
  odist = ntimes
  ostride = 1
  call dfftw_plan_many_dft_c2r(plan_fftb_ker, rank,n,howmany,kerspec, inembed, istride, &
                               idist, kerphys, onembed, ostride, odist, fftw_mode)
       
  if (mynum == 0) write(6,*) 'done with preparing fftw plans'

  allocate(times(1:ntimes))
  times(:) = 0. 
  allocate(omega(0:nomega))
  omega(0:nomega) = 0. 
  allocate(period(0:nomega))
  period(0:nomega) = 0. 
  allocate(frequency(0:nomega))
  frequency(0:nomega) = 0. 

  times(2:ntimes-npad) = time(1:ndumps)
  dt = time(ndumps) / real(ndumps)
  if (mynum==0) write(6,*) 'the numerical timestep is ', dt, ' s'

  do itimes = ntimes-npad+1, ntimes
     times(itimes) = (itimes - 1) * dt
  end do
  timewidth = times(ntimes)-times(1)
  if (mynum==0) write(6,*) 'the width of the (potentially padded) time window:', &
                           timewidth, 's.' 

  do iomega=0, nomega
    omega(iomega) = 2. * pi * real(iomega) / timewidth
    frequency(iomega) = real(iomega) / timewidth ! iomega cycles per sample
    if (iomega > 0) &
       period(iomega) = timewidth / real(iomega) 
       !period(0) has no meaning (corresponds to the mean) 
  end do

  allocate(uzrec(1:ntimes))
  uzrec(:) = 0. 
  uzrec(2:ndumps+1) = dsin(dble(thr)) * u0(1:ndumps,1) + &
                      dcos(dble(thr)) * u0(1:ndumps,2) ! for amp

  allocate(vzrec(1:ntimes))
  vzrec(:) = 0. 
  vzrec(2:ndumps+1) = dsin(dble(thr)) * v0(1:ndumps,1) + &
                      dcos(dble(thr)) * v0(1:ndumps,2) ! for TT

  write(6,*) 'Defining the filter.......'
  allocate(specfilt(0:nomega))
  call define_filter(filter_type)
  write(6,*) 'FINISHED filter definition'

  if (mynum==0) write(6,*) 'Done with initializing stufft' 

end subroutine init_fft
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine exit_fft

  call dfftw_destroy_plan(plan_fftb_ker)
  call dfftw_destroy_plan(plan_fftf_ker)
  call dfftw_destroy_plan(plan_fftb1d)
  call dfftw_destroy_plan(plan_fftf1d)

end subroutine exit_fft
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fftf1d_dble(uphys, uspec)
  use data_mesh, only : ndumps

  real(kind=8), dimension(1:ndumps), intent(in)     :: uphys
  complex(kind=8), dimension(0:nomega), intent(out) :: uspec
  
  uphys1d(:) = 0.d0
  uphys1d(1:ndumps) = dble(uphys(1:ndumps))
  write(6,*) 'inside fftf1d_dble', maxval(uphys), maxval(uphys1d)

  call dfftw_execute(plan_fftf1d)
  uspec(:) = uspec1d(:)
  write(6,*) 'executed with fftf1d_dble', maxval(real(uspec))
  write(6,*) 'done with fftf1d_dble'

end subroutine fftf1d_dble
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fftf1d(uphys,uspec)
  use data_mesh, only : ndumps
  
  real(kind=realkind), dimension(1:ndumps), intent(in)  :: uphys
  complex(kind=8), dimension(0:nomega), intent(out)     :: uspec
  
  uphys1d(:) = 0.d0
  uphys1d(1:ndumps) = dble(uphys(1:ndumps))
  write(6,*) 'inside fftf1d', maxval(uphys), maxval(uphys1d)

  call dfftw_execute(plan_fftf1d)
  uspec(:) = uspec1d(:)
  write(6,*) 'executed with fftf1d', maxval(real(uspec))
  write(6,*) 'done with fftf1d'
 
end subroutine fftf1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fftb1d_dble(uspec,uphys)
  use data_mesh, only : ndumps

  complex(kind=8), dimension(0:nomega), intent(in)  :: uspec
  real(kind=8) , dimension(1:ndumps), intent(out)   :: uphys

  real(kind=8)  :: arn

  arn = 1.d0 / real(ntimes)
  uspec1d(0:nomega) = uspec(0:nomega)
  write(6,*)'executing fftb1d',size(uspec1d)

  call dfftw_execute(plan_fftb1d)
  uphys(:) = arn * uphys1d(1:ndumps) 
  write(6,*) 'done executing fftb1d'

end subroutine fftb1d_dble
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fftb1d(uspec,uphys)
  use data_mesh, only : ndumps

  complex(kind=8), dimension(0:nomega), intent(in)        :: uspec
  real(kind=realkind), dimension(1:ndumps), intent(out)   :: uphys

  real(kind=realkind)  :: arn

  arn = 1.d0 / real(ntimes)
  uspec1d(0:nomega) = uspec(0:nomega)
  write(6,*) 'executing fftb1d', size(uspec1d)

  call dfftw_execute(plan_fftb1d)
  uphys(:) = arn * uphys1d(1:ndumps) 
  write(6,*) 'done executing fftb1d'

end subroutine fftb1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fftf_ker(fspec)
  complex(kind=8), dimension(0:nomega,1:npts,1),intent(out) :: fspec
  integer :: i

  !write(6,*) npts,nomega, ntimes, 'size fspec:', size(fspec)
  !write(6,*) size(kerphys), maxval(kerphys)
!  if (filter_what > 0) then 
!     write(6,*) mynum, 'tapering input signal...'
!     do i=1, npts
!        call taper(ntimes,kerphys(:,i))
!     enddo
!  endif

  write(6,*) mynum, 'executing forward plan...'
  call dfftw_execute(plan_fftf_ker)

  write(6,*) size(kerspec), maxval(real(kerspec))

  fspec(0:nomega,1:npts,1) = kerspec(0:nomega,1:npts)

end subroutine fftf_ker
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine taper(npts, signal)

  implicit none
  integer, intent(in)           :: npts
  real(kind=8), intent(inout)   :: signal(npts)

  integer               :: nabs, i 
  real, allocatable     :: win(:)
  real, parameter       :: D = 3.3
  
  nabs = floor(real(npts) / real(4.))
  allocate(win(1:nabs))
  do i=1, nabs
     win(i) = win(i) * exp( -(D * real(i) / real(nabs)) * (D * real(i) / real(nabs)) )
     !write(66,*) win(i) !MvD: file 66 is never opened
  enddo
  
  signal(npts-nabs+1:npts) = signal(npts-nabs+1:npts) * win(1:nabs)

end subroutine taper
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fftb_ker(fphys,fspec)
  implicit none
  include 'mesh_params_kernel.h'
  
  real(kind=4) , dimension(1:ndumps,1:npts),intent(out)  :: fphys
  complex(kind=8), dimension(0:nomega,1:npts),intent(in) :: fspec
  real(kind=4) :: arn

  arn = 1. / real(ntimes)

  !kerspec(0:nomega,1:npts)=conjg(transpose(fspec(0:nomega,1:npts)))
  kerspec = fspec
  !write(6,*)'dfftw';call flush(6)
  
  call dfftw_execute(plan_fftb_ker)
  fphys(1:ndumps,1:npts) = arn * kerphys(1:ndumps,1:npts)

end subroutine fftb_ker
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fftb_ker_dble(fphys,fspec)
  implicit none
  include 'mesh_params_kernel.h'

  real(kind=8) , dimension(1:ndumps,1:npts),intent(out)  :: fphys
  complex(kind=8), dimension(0:nomega,1:npts),intent(in) :: fspec
  double precision :: arn

  arn = 1.d0 / dble(ntimes)

  !kerspec(0:nomega,1:npts)=conjg(transpose(fspec(0:nomega,1:npts)))
  kerspec = fspec
  !write(6,*)'dfftw';call flush(6)

  call dfftw_execute(plan_fftb_ker)
  fphys(1:ndumps,1:npts) = arn * kerphys(1:ndumps,1:npts)

end subroutine fftb_ker_dble
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fftb_ker_multiply(fphys, ffwd, fbwd, n)
  implicit none
  include 'mesh_params_kernel.h'

  real(kind=4), dimension(1:ndumps,1:npts), intent(out)       :: fphys
  complex(kind=8), dimension(0:nomega,1:npts,1:n), intent(in) :: ffwd, fbwd
  integer, intent(in)                                         :: n

  real(kind=4) :: arn
  integer      :: idim

  write(6,*) 'fftb including multiply ndim:' ,n
  arn = 1.d0 / dble(ntimes)
  kerspec(0:nomega,1:npts) = 0.
  do idim=1, n
     kerspec(0:nomega,1:npts) = kerspec(0:nomega,1:npts) + &
                                ffwd(0:nomega,1:npts,idim) * fbwd(0:nomega,1:npts,idim)
  enddo

  call dfftw_execute(plan_fftb_ker)
  fphys(1:ndumps,1:npts) = arn * kerphys(1:ndumps,1:npts)

end subroutine fftb_ker_multiply
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prepare_kernel(kernel, tmean)
  implicit none
  include 'mesh_params_kernel.h'

  double precision, intent(inout) :: kernel(1:ndumps,1:npts)
  double precision, intent(out)   :: tmean(npts)

  integer :: it, isize

  ! remove mean 
  do isize=1, npts
     tmean(isize) = sum(kernel(1:ndumps,isize), dim=1) / real(ndumps)
     kernel(1:ndumps,isize) = kernel(1:ndumps,isize) - tmean(isize)
  end do

end subroutine prepare_kernel
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prepare_kernel_realkind(kernel,tmean)
  implicit none
  include 'mesh_params_kernel.h'

  real(kind=realkind), intent(inout) :: kernel(1:ndumps,1:npts)
  double precision, intent(inout)    :: tmean(npts)

  integer :: it, isize

  ! remove mean 
  do isize=1, npts
     tmean(isize) = sum(kernel(1:ndumps,isize), dim=1) / real(ndumps)
     kernel(1:ndumps,isize) = kernel(1:ndumps,isize) - tmean(isize)
  end do

end subroutine prepare_kernel_realkind
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prepare_kernel_real(kernel,tmean)
  implicit none
  include 'mesh_params_kernel.h'
  real(kind=4), intent(inout)       :: kernel(1:ndumps,1:npts)
  double precision, intent(inout)   :: tmean(npts)

  integer :: it, isize

  ! remove mean 
  do isize=1, npts
     tmean(isize) = SUM(kernel(1:ndumps,isize),dim=1) / real(ndumps)
     kernel(1:ndumps,isize) = kernel(1:ndumps,isize) - tmean(isize)
  end do

end subroutine prepare_kernel_real
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prepare_fwd_field(ufwd,tmean)
  implicit none
  include 'mesh_params_kernel.h'

  real(kind=realkind), intent(out)  :: ufwd(npts,ndumps)
  double precision, intent(inout)   :: tmean(npts)

  integer :: it, isize

  ! remove mean 
  do isize=1, npts
     tmean(isize) = SUM(ufwd(isize,:),dim=1) / real(ndumps)
     ufwd(isize,:) = ufwd(isize,:) - tmean(isize)
  end do

end subroutine prepare_fwd_field
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prepare_bwd_field(ubwd,tmean)
  implicit none
  include 'mesh_params_kernel.h'

  real(kind=realkind) , intent(inout) :: ubwd(npts,ndumps)
  double precision, intent(inout)     :: tmean(npts)

  integer :: it, isize

  ! remove mean 
  do isize=1, npts
     tmean(isize) = SUM(ubwd(isize,:),dim=1) / real(ndumps)
     ubwd(isize,:) = ubwd(isize,:) - tmean(isize)
  end do

end subroutine prepare_bwd_field
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=realkind) function bh4(t,timew)
  !blackman harris window
  implicit none

  real(kind=realkind) , intent(in) :: t, timew
  real(kind=realkind), parameter   :: a1 =  .35875
  real(kind=realkind), parameter   :: a2 = -.48829
  real(kind=realkind), parameter   :: a3 =  .14128
  real(kind=realkind), parameter   :: a4 =-0.01168
  real(kind=realkind)  :: pi, arg
  
  pi = 2. * asin(1.)
  arg = t / timew
  bh4 = a1 + &
        a2 * cos(2. * pi * arg) + &
        a3 * cos(4. * pi * arg) + &
        a4 * cos(6. * pi * arg)
end function bh4
!-----------------------------------------------------------------------------------------

!=================
end module fft
!=================
