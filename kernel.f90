!=========================================================================================
module kernel
use global_parameters,                    only: sp, dp, verbose, lu_out, firstslave
use filtering,                            only: filter_type
use commpi,                               only: pabort
implicit none

    type kernelspec_type
        character(len=32), public            :: name
        real(kind=dp), dimension(2)          :: time_window
        real(kind=dp), allocatable           :: seis(:)           ! Seismogram (Velocity or displacement)
                                                                  ! in the time window of 
                                                                  ! this kernel
        real(kind=dp), allocatable           :: t(:)
        real(kind=dp)                        :: dt
        real(kind=dp)                        :: normalization
        integer                              :: filter_type
        character(len=4)                     :: misfit_type
        character(len=16)                    :: model_parameter
        type(filter_type), pointer           :: filter
        logical                              :: initialized = .false.
        contains 
           procedure, pass                   :: init 
           procedure, pass                   :: calc_misfit_kernel
           procedure, pass                   :: isinitialized
           procedure, pass                   :: freeme
           procedure, pass                   :: apply_filter_2d
           procedure, pass                   :: apply_filter_3d
           generic                           :: apply_filter => apply_filter_2d, apply_filter_3d
    end type

contains 

!-------------------------------------------------------------------------------
subroutine init(this, name, time_window, filter, misfit_type, model_parameter, &
                seis, dt, timeshift_fwd, deconv_stf, write_smgr)
   use fft,                           only  : rfft_type, taperandzeropad
   use filtering,                     only  : timeshift_type

   class(kernelspec_type)                  :: this
   character(len=*), intent(in)            :: name
   real(kind=dp), intent(in)               :: time_window(2)
   type(filter_type), target, intent(in)   :: filter
   character(len=4), intent(in)            :: misfit_type
   character(len=*), intent(in)            :: model_parameter
   real(kind=dp), intent(in)               :: seis(:)
   real(kind=dp), intent(in)               :: dt
   real(kind=dp), intent(in), optional     :: timeshift_fwd
   logical, intent(in), optional           :: deconv_stf
   logical, intent(in), optional           :: write_smgr
   
   type(rfft_type)                         :: fft_data
   type(timeshift_type)                    :: timeshift
   complex(kind=dp), allocatable           :: seis_fd(:,:)
   real(kind=dp),    allocatable           :: seis_td(:,:), t_cut(:)
   real(kind=dp),    allocatable           :: seis_filtered(:,:)
   character(len=32)                       :: fmtstring
   logical                                 :: deconv_stf_loc, write_smgr_loc

   integer                                 :: ntimes, ntimes_ft, nomega, isample

   if(this%initialized) then
      write(*,*) 'This kernel is already initialized'
      call pabort 
   end if
   this%name              = trim(name)
   this%time_window       = time_window 
   this%filter            => filter
   this%misfit_type       = misfit_type
   this%model_parameter   = model_parameter
   this%initialized       = .true.
   this%dt                = dt

   deconv_stf_loc = .false.
   if (present(deconv_stf)) then
     deconv_stf_loc = deconv_stf
   end if

   write_smgr_loc = .false.
   if (present(write_smgr)) then
     write_smgr_loc = write_smgr
   end if

   if ( (.not.deconv_stf_loc).and.(.not.present(timeshift_fwd)) ) then
     write(*,*) 'Initialize Kernel needs a value for timeshift_fwd, if '
     write(*,*) 'deconv_stf is set to false'
     call pabort()
   end if
     

   if (verbose>0) then
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,'(2(A))')         '   Initialize kernel ', this%name
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,'(A,2(F8.1))')    '   Time window:  ', this%time_window
      write(lu_out,'(2(A))')         '   Misfit type:  ', this%misfit_type
      write(lu_out,'(2(A))')         '   Model param:  ', this%model_parameter
      write(lu_out,'(A,L)')          '   Initialized:  ', this%initialized
      write(lu_out,'(2(A))')         '   Filter type:  ', this%filter%filterclass
      write(lu_out,'(A,4(F8.3))')    '   Filter freq:  ', this%filter%frequencies
      if (present(timeshift_fwd)) then
        write(lu_out,'(A,F8.3)')       '   Time shift :  ', timeshift_fwd
      end if
      write(lu_out,'(A,L)')          '   Deconvolve :  ', deconv_stf_loc
      write(lu_out,'(A,L)')          '   Write smgrs:  ', write_smgr_loc
   end if


   ! Save seismogram in the timewindow of this kernel

   ! Allocate temporary variable
   ntimes = size(seis, 1)
   allocate(seis_td(ntimes, 1))
   seis_td(:,1) = seis

   ! Initialize FFT for the seismogram
   call fft_data%init(ntimes, 1, 1, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   if (verbose>0) then
      fmtstring = '(A, I8, A, I8)'
      write(lu_out,fmtstring) '   ntimes: ',  ntimes,     '  , nfreq: ', nomega
   end if

   allocate(seis_fd(nomega, 1))
   allocate(seis_filtered(ntimes_ft, 1))
   allocate(this%t(ntimes_ft))
   this%t = fft_data%get_t()


   ! FFT, timeshift and filter the seismogram
   call fft_data%rfft(taperandzeropad(seis_td, ntimes_ft), seis_fd)
 
   if (.not.deconv_stf_loc) then
     ! It's slightly ineffective to init that every time, but it is only 
     ! called once per receiver at initialization.
     call timeshift%init_ts(fft_data%get_f(), dtshift = timeshift_fwd)
     call timeshift%apply(seis_fd)
     call timeshift%freeme()
   end if

   seis_fd = this%filter%apply_2d(seis_fd)
   call fft_data%irfft(seis_fd, seis_filtered)
   call fft_data%freeme()

   ! Cut timewindow of the kernel from the seismogram
   call cut_timewindow(this%t,             &
                       seis_filtered(:,1), &
                       this%time_window,   &
                       this%seis )
   call cut_timewindow(this%t,             &
                       this%t,             &
                       this%time_window,   &
                       t_cut )

   if (sum(this%seis**2).lt.1.d-100) then
       this%normalization = 0
   else
       this%normalization = 1.d0/sum(this%seis**2)
   end if
   
   if (verbose>0) then
      write(lu_out,*) '  Normalization coefficient: ', this%normalization
      write(lu_out,*) '  Length of seismogram: ', size(this%seis,1), ' samples'
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,*) ''
   end if

   if (firstslave.and.write_smgr_loc) then
      open(unit=100,file='seismogram_raw_'//trim(this%name), action='write')
      do isample = 1, size(seis,1)
         write(100,*) this%t(isample), seis(isample)
      end do
      close(100)

      open(unit=100,file='seismogram_'//trim(this%name), action='write')
      do isample = 1, size(seis_filtered,1)
         write(100,*) this%t(isample), seis_filtered(isample,1)
      end do
      close(100)

      open(unit=100,file='seismogram_cut_'//trim(this%name), action='write')
      do isample = 1, size(this%seis,1)
         write(100,*) t_cut(isample), this%seis(isample)
      end do
      close(100)
   end if

end subroutine init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
logical function isinitialized(this)
   class(kernelspec_type)                   :: this
   isinitialized = this%initialized
end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine freeme(this)
   class(kernelspec_type)                   :: this
   this%filter => NULL()
   this%initialized = .false.
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function calc_misfit_kernel(this, timeseries)
   !< This routine can take multiple time series and calculate the kernel at each
   class(kernelspec_type)                  :: this
   real(kind=dp), intent(in)               :: timeseries(:,:)
   real(kind=dp)                           :: calc_misfit_kernel(size(timeseries,2))
   real(kind=dp), allocatable              :: cut_timeseries(:)
   integer                                 :: itrace, ntrace, lu_errorlog
   character(len=64)                       :: fmtstring

   ntrace = size(timeseries,2)

   select case(trim(this%misfit_type))
   case('CC')

      do itrace = 1, ntrace
         
         ! print*,"TEST99",timeseries(:,itrace)
         ! print*,"TEST88",this%seis

         call cut_timewindow(this%t,                 &
                             timeseries(:, itrace),  &
                             this%time_window,       &
                             cut_timeseries)

         calc_misfit_kernel(itrace) = integrate( cut_timeseries * this%seis, this%dt ) &
                                      * this%normalization
      end do

   case('AM')

      ! @TODO: Amplitude kernels not yet implemented

      ! do itrace = 1, ntrace

      !    call cut_timewindow(this%t,                 &
      !                        timeseries(:, itrace),  &
      !                        this%time_window,       &
      !                        cut_timeseries)
         
      !    calc_misfit_kernel(itrace) = integrate( cut_timeseries * this%seis, this%dt ) &
      !                                 * this%normalization
      ! end do

   end select

   if (any(calc_misfit_kernel.ne.calc_misfit_kernel)) then ! Kernel is NaN
       print *, 'ERROR Kernel has value NaN!'
       print *, 'Values of Kernel:'
       write(fmtstring, "('(',I6,'(ES11.3))')") ntrace
       print fmtstring, calc_misfit_kernel
       print *, 'Value of normalization:'
       print *, this%normalization
       print *, 'Convolved time series is written into "errorlog_timeseries",'
       print *, 'seismogram is written into "errorlog_seismogram"'
       open(newunit=lu_errorlog, file='errorlog_timeseries', status='replace')
       write(lu_errorlog, *) cut_timeseries
       close(lu_errorlog)
       open(newunit=lu_errorlog, file='errorlog_seismogram', status='replace')
       write(lu_errorlog, *) this%seis
       close(lu_errorlog)
       call pabort
   end if

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function assemble_basekernel(basekernel_id, strain_type, &
                             fw_field_fd, bw_field_fd) result(conv_field_fd)

!< Assembles the basic (waveform) kernels from the forward and the backward field
!! in frequency domain.
!! readfields.f90 returns strain tensor in voigt notation
!! elements 1,...,6 contain tt,pp,rr,pr,tr,tp. See Fichtner
!! book, p. 168-169 how the fundamental kernels K_a, K_b and
!! K_c are defined

  character(len=*), intent(in)      :: basekernel_id, strain_type
  complex(kind=dp), intent(in)      :: fw_field_fd(:,:,:), bw_field_fd(:,:,:)
  complex(kind=dp)                  :: conv_field_fd(size(fw_field_fd,1),size(fw_field_fd,3))
  
  select case(trim(basekernel_id))
  case('k_lambda')
     select case(trim(strain_type))
     case('straintensor_trace')
        conv_field_fd = sum(fw_field_fd * bw_field_fd, 2)
     case('straintensor_full')
        ! Only convolve trace of fw and bw fields
        conv_field_fd = fw_field_fd(:,1,:) * &
                        fw_field_fd(:,2,:) * &
                        fw_field_fd(:,3,:) * &
                        bw_field_fd(:,1,:) * &
                        bw_field_fd(:,2,:) * &
                        bw_field_fd(:,3,:)
     end select
  case('k_mu')
     conv_field_fd = sum(fw_field_fd * bw_field_fd, 2) * 2.d0
  case('k_a')
     conv_field_fd = ( bw_field_fd(:,2,:) + bw_field_fd(:,1,:) ) * &
                     ( fw_field_fd(:,2,:) * bw_field_fd(:,1,:) )
  case('k_b')
     conv_field_fd = ( ( bw_field_fd(:,5,:) * fw_field_fd(:,5,:) ) + &
                     ( bw_field_fd(:,4,:) * fw_field_fd(:,4,:) ) ) * 4.d0
  case('k_c')
     conv_field_fd = ( fw_field_fd(:,1,:) + fw_field_fd(:,2,:) ) * &
                       bw_field_fd(:,3,:) + &
                     ( bw_field_fd(:,1,:) * bw_field_fd(:,2,:) ) * &
                       fw_field_fd(:,3,:)
  case('k_rho')
     print*,"ERROR: Density kernels not yet implemented"
     stop
  end select

end function assemble_basekernel
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Here the we tabulate the base kernels required to assemble 
!! the physical kernels for a desired model parameter
!! TODO: Move this to own module, together with assemble_basekernel, but wait for 
!!       merge with Ludwigs branch - SCS 240914
subroutine tabulate_kernels(model_param, nbasekernels, basekernel_id, strain_type)
  character(len=*), intent(in)   :: model_param
  integer, intent(out)           :: nbasekernels
  character(len=16), allocatable :: basekernel_id(:)
  character(len=32), intent(out) :: strain_type

  select case(trim(model_param))     
  case('lambda')
     nbasekernels     = 1
     allocate(basekernel_id(nbasekernels))

     basekernel_id(1) = 'k_lambda'
     strain_type      = 'straintensor_trace'

  case('vp')
     nbasekernels     = 1
     allocate(basekernel_id(nbasekernels))
     
     basekernel_id(1) = 'k_lambda'     
     strain_type      = 'straintensor_trace'

  case('vs')
     nbasekernels     = 2
     allocate(basekernel_id(nbasekernels))

     basekernel_id(1) = 'k_lambda'
     basekernel_id(2) = 'k_mu'      
     strain_type      = 'straintensor_full'
  
  case('mu')
     nbasekernels     = 1
     allocate(basekernel_id(nbasekernels))

     basekernel_id(1) = 'k_mu'
     strain_type      = 'straintensor_full'
  
  case('vsh')
     nbasekernels     = 4
     allocate(basekernel_id(nbasekernels))

     basekernel_id(1) = 'k_lambda'
     basekernel_id(2) = 'k_mu'
     basekernel_id(3) = 'k_b'
     basekernel_id(4) = 'k_c'
     strain_type      = 'straintensor_full'
  
  case('vsv')
     nbasekernels     = 1
     allocate(basekernel_id(nbasekernels))

     basekernel_id(1) = 'k_b'
     strain_type      = 'straintensor_full'
  
  case('vph')
     nbasekernels     = 2
     allocate(basekernel_id(nbasekernels))

     basekernel_id(1) = 'k_a'
     basekernel_id(2) = 'k_b'
     strain_type      = 'straintensor_full'

  case('vpv')
     nbasekernels     = 3
     allocate(basekernel_id(nbasekernels))

     basekernel_id(1) = 'k_lambda'
     basekernel_id(2) = 'k_a'
     basekernel_id(3) = 'k_c'
     strain_type      = 'straintensor_full'

  case('kappa')
     write(*,*) "Error: Kappa kernels not yet implemented"
     call pabort
  case('eta')
     write(*,*) "Error: Eta kernels not yet implemented"
     call pabort
  case('xi')
     write(*,*) "Error: Xi kernels not yet implemented"
     call pabort
  case('rho')
     write(*,*) "Error: Density kernels not yet implemented"
     call pabort
  case default
     write(*,*) "Error in tabulate_kernels: Unknown model parameter", trim(model_param)
     call pabort      
  end select
end subroutine tabulate_kernels
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function calc_physical_kernels(model_param, base_kernel, bg_model)  &
                               result(physical_kernel)

  use backgroundmodel, only               :  backgroundmodel_type
  character(len=*), intent(in)           :: model_param
  real(kind=dp),    intent(in)           :: base_kernel(:,:)
  type(backgroundmodel_type), intent(in) :: bg_model
  real(kind=dp)                          :: physical_kernel(size(base_kernel,1))

  ! See Fichtner p. 169 how kernels are defined and type_parameter.f90
  ! how basekernels inside kernelvalue_weighted are ordered
  select case(trim(model_param))
  case('lambda')
     physical_kernel(:) =                                 &
            base_kernel(:, 1)
  case('mu')
     physical_kernel(:) =                                 & 
            base_kernel(:, 1)
  case('vp')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vp *           &
            base_kernel(:, 1)
  case('vs')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vs *           &
            base_kernel(:, 2) - 4.d0 * bg_model%c_rho * bg_model%c_vs * &
            base_kernel(:, 1)
  case('vsh')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vsh *          &
          ( base_kernel(:, 1) * 2 +                   &
            base_kernel(:, 2) +                       &
            base_kernel(:, 3) +                       &
            base_kernel(:, 4) * 2*(1-bg_model%c_eta))
  case('vsv')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vsv *          &
            base_kernel(:, 1)
  case('vph')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vph *          &
          ( base_kernel(:, 1) +                       &
            base_kernel(:, 2) * bg_model%c_eta )
  case('vpv')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vpv *          &
          ( base_kernel(:, 1) +                       &
            base_kernel(:, 2) +                       &
            base_kernel(:, 3) )
  case('kappa')
     write(*,*) "Error: Kappa kernels not yet implemented"
     stop
  case('eta')
     write(*,*) "Error: Eta kernels not yet implemented"
     stop
  case('xi')
     write(*,*) "Error: Xi kernels not yet implemented"
     stop
  case('rho')
     write(*,*) "Error: Density kernels not yet implemented"
     stop
  end select

end function
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine cut_timewindow(t, x, timewindow, cut_tw, ierror)
   real(kind=dp), intent(in)        :: t(:), x(:)
   real(kind=dp), intent(in)        :: timewindow(2)
   integer, intent(out), optional   :: ierror
   real(kind=dp), allocatable       :: cut_tw(:)
   real(kind=dp), allocatable       :: cut_timewindow_temp(:)
   integer                          :: ntimes, i, iintimewindow
   real(kind=dp)                    :: dt, tlast
  

   tlast = -1d308
   do i=1, size(t)
      if (t(i)<=tlast) then
         write(*,*) 'Error in cut_timewindow:'
         write(*,*) 'times t(:) must be monotonously increasing'
         if (present(ierror)) then
            ierror = -1
            return
         else
            call pabort
         end if
      else
         tlast = t(i)
      end if
   end do

   if (present(ierror)) ierror = 0

   ntimes = size(t)
   dt = t(2) - t(1)
   if (timewindow(1).lt.t(1)) then
      write(*,*) 'Time window starts before beginning of time series'
      write(*,*) 'time window:', timewindow, '; t(1):', t(1)
      if (present(ierror)) then
         ierror = -1
         return
      else
         call pabort
      end if
   end if
   if (timewindow(2).gt.t(ntimes)) then
      write(*,*) 'Time window ends after end of time series'
      write(*,*) 'time window:', timewindow, '; t(ntimes):', t(ntimes)
      if (present(ierror)) then
         ierror = -1
         return
      else
         call pabort
      end if
   end if

   if ((timewindow(2) - timewindow(1)).le.0) then
      write(*,*) 'length of time window is negative'
      write(*,*) 'Beginning: ', timewindow(1), '; end: ', timewindow(2)
      if (present(ierror)) then
         ierror = -1
         return 
      else
         call pabort
      end if
   end if

   allocate(cut_timewindow_temp(ntimes))
   iintimewindow = 0
   do i = 1, ntimes
       if (t(i).le.timewindow(2).and.t(i).ge.timewindow(1)) then
          iintimewindow = iintimewindow + 1
          cut_timewindow_temp(iintimewindow) = x(i)
       end if
       if (t(i).gt.timewindow(2)) exit
   end do

   if (iintimewindow.eq.0) then
       print *, 'Time window: ', timewindow, ', t(1):', t(1), ', t(ntimes):', t(ntimes)
       write(*,*) 'Time window length was zero'
       if (present(ierror)) then
          ierror = -1
          return
       else
          call pabort
       end if
   end if

   if (allocated(cut_tw)) then
       deallocate(cut_tw)
   end if
   allocate(cut_tw(iintimewindow))

   cut_tw(:) = cut_timewindow_temp(1:iintimewindow)
   

end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function apply_filter_3d(this, freq_series)
   class(kernelspec_type)           :: this
   complex(kind=dp), intent(in)     :: freq_series(:,:,:)
   complex(kind=dp)                 :: apply_filter_3d(size(freq_series,1), &
                                                       size(freq_series,2), &
                                                       size(freq_series,3))


   apply_filter_3d = this%filter%apply_3d(freq_series)
   apply_filter_3d = this%filter%apply_3d(apply_filter_3d)

end function apply_filter_3d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function apply_filter_2d(this, freq_series)
   class(kernelspec_type)           :: this
   complex(kind=dp), intent(in)     :: freq_series(:,:)
   complex(kind=dp)                 :: apply_filter_2d(size(freq_series,1), &
                                                       size(freq_series,2))


   apply_filter_2d = this%filter%apply_2d(freq_series)
   apply_filter_2d = this%filter%apply_2d(apply_filter_2d)

end function apply_filter_2d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
pure function integrate(timeseries, dt)
    real(kind=dp), intent(in)   :: timeseries(:), dt
    real(kind=dp)               :: integrate
    integer                     :: npoints

    npoints = size(timeseries, 1)

    ! Trapezoidal rule: I = dt/2 * (f(x1) + 2f(x2) + ... + 2f(xN-1) + f(xN))
    integrate = timeseries(1) + timeseries(npoints) + 2*sum(timeseries(2:npoints-1), 1)
    integrate = integrate * dt * 0.5

end function integrate
!-------------------------------------------------------------------------------

end module
!=========================================================================================
