!=========================================================================================
module kernel
use global_parameters,                    only: sp, dp, verbose, lu_out, firstslave, myrank
use filtering,                            only: filter_type
use fft,                                  only: rfft_type, taperandzeropad
use commpi,                               only: pabort
implicit none

  type kernelspec_type
    private
    character(len=32), public            :: name
    real(kind=dp), dimension(2)          :: time_window
    real(kind=dp), allocatable           :: seis_cut(:)       ! Seismogram (Velocity or displacement)
                                                              ! in the time window of 
                                                              ! this kernel
    complex(kind=dp), allocatable        :: seis_cut_fd(:,:)  ! Seismogram (Velocity or displacement)
                                                              ! in the time window of 
                                                              ! this kernel in frequency domain
    real(kind=dp), allocatable           :: t(:)
    real(kind=dp), allocatable, public   :: t_cut(:)
    real(kind=dp)                        :: dt
    real(kind=dp)                        :: normalization
    integer                              :: filter_type
    character(len=4)                     :: misfit_type
    character(len=16), public            :: model_parameter
    integer,           public            :: model_parameter_index
    integer,           public            :: hetero_parameter_index
    character(len=32)                    :: strain_type 
    type(filter_type), pointer           :: filter
    integer, public                      :: ntimes            !< Length of time window in samples
    type(rfft_type)                      :: fft_data          !< FFT type for Parseval integration
    real(kind=dp), allocatable           :: datat(:,:)        !< input vector for FFT
    complex(kind=dp), allocatable        :: dataf(:,:)        !< result vector for FFT
    integer                              :: ntimes_ft         !< Length of FFT input vector
    integer                              :: nomega            !< Length of FFT result vector
    logical                              :: initialized = .false.

    integer, public                      :: nbasekernels
    logical, public                      :: needs_basekernel(6) = .false.
                                                                !< Which of the base
                                                                !! kernels does this
                                                                !! phys. kernel need
                                                                !! lam, mu, rho, a, b, c 
    contains 
       procedure, pass                   :: init 
       procedure, pass                   :: cut_and_add_seismogram
       procedure, pass                   :: calc_misfit_kernel
       procedure, pass                   :: isinitialized
       procedure, pass                   :: freeme
       procedure, pass                   :: apply_filter_2d
       procedure, pass                   :: apply_filter_3d
       generic                           :: apply_filter => apply_filter_2d, apply_filter_3d
       procedure, pass                   :: integrate_parseval_both_real
       procedure, pass                   :: integrate_parseval_b_complex
       generic                           :: integrate_parseval => integrate_parseval_both_real, &
                                                                  integrate_parseval_b_complex
  end type

contains 

!-------------------------------------------------------------------------------
subroutine init(this, name, time_window, filter, misfit_type, model_parameter, &
                seis, dt, timeshift_fwd, deconv_stf, write_smgr)
   use backgroundmodel, only                : get_parameter_index
   use heterogeneities, only                : get_hetparam_index
              
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
   
   character(len=32)                       :: fmtstring
   logical                                 :: deconv_stf_loc, write_smgr_loc
   real(kind=dp)                           :: timeshift_fwd_loc

   integer                                 :: ntimes, isample, iparam

   if(this%initialized) then
      write(*,*) 'This kernel is already initialized'
      call pabort 
   end if
   this%name              = trim(name)
   this%time_window       = time_window 
   this%filter            => filter
   this%misfit_type       = misfit_type
   this%model_parameter   = model_parameter
   this%dt                = dt

   ! Test argument consistency
   deconv_stf_loc = .false.
   if (present(deconv_stf)) then
     deconv_stf_loc = deconv_stf
     timeshift_fwd_loc = 0
   else
     if (.not.present(timeshift_fwd)) then
       write(*,*) 'Initialize Kernel needs a value for timeshift_fwd, if '
       write(*,*) 'deconv_stf is set to false'
       call pabort()
     else
       timeshift_fwd_loc = timeshift_fwd
     end if
   end if

   write_smgr_loc = .false.
   if (present(write_smgr)) then
     write_smgr_loc = write_smgr
   end if

   ! Cut seismogram to kernel time window and add to kernel type
   call this%cut_and_add_seismogram(seis, deconv_stf_loc, write_smgr_loc, timeshift_fwd_loc)

   ! Check and tabulate, which base kernels the model parameter for this 
   ! specific kernel needs
   call tabulate_kernels(this%model_parameter, this%needs_basekernel, this%strain_type)

   this%model_parameter_index = get_parameter_index(this%model_parameter)
   this%hetero_parameter_index = get_hetparam_index(this%model_parameter)

   if (verbose>0) then
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,'(2(A))')         '   Initialize kernel ', this%name
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,'(A,2(F8.1))')    '   Time window:  ', this%time_window
      write(lu_out,'(A,I5)')         '   nsamples   :  ', this%ntimes
      write(lu_out,'(2(A))')         '   Misfit type:  ', this%misfit_type
      write(lu_out,'(2(A))')         '   Model param:  ', this%model_parameter
      write(lu_out,'(A, I3)')        '   Param index:  ', this%model_parameter_index
      write(lu_out,'(A, I3)')        '   Het3d index:  ', this%hetero_parameter_index
      write(lu_out,'(2(A))')         '   Filter type:  ', this%filter%filterclass
      write(lu_out,'(A,4(F8.3))')    '   Filter freq:  ', this%filter%frequencies
      if (present(timeshift_fwd)) then
        write(lu_out,'(A,F8.3)')     '   Time shift :  ', timeshift_fwd
      end if
      write(lu_out,'(A,L1)')         '   Deconvolve :  ', deconv_stf_loc
      write(lu_out,'(A,L1)')         '   Write smgrs:  ', write_smgr_loc
   end if

   this%initialized       = .true.

end subroutine init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Cuts seismogram to kernel time window and adds it to the kernel type
!! Is called by init subroutine
subroutine cut_and_add_seismogram(this, seis, deconv_stf, write_smgr, timeshift_fwd)
   use filtering,                     only  : timeshift_type
   class(kernelspec_type)                  :: this
   real(kind=dp), intent(in)               :: seis(:)
   logical, intent(in)                     :: deconv_stf
   logical, intent(in)                     :: write_smgr
   real(kind=dp), intent(in)               :: timeshift_fwd

   type(rfft_type)                         :: fft_data
   type(timeshift_type)                    :: timeshift
   complex(kind=dp), allocatable           :: seis_fd(:,:)
   real(kind=dp),    allocatable           :: seis_td(:,:), t_cut(:)
   real(kind=dp),    allocatable           :: seis_filtered(:,:)
   real(kind=dp)                           :: normalization_term
   integer                                 :: ntimes, ntimes_ft, nomega, isample

   ! Save seismogram in the timewindow of this kernel

   ! Allocate temporary variable
   ntimes = size(seis, 1)
   allocate(seis_td(ntimes, 1))
   seis_td(:,1) = seis

   ! Initialize FFT type for the seismogram 
   call fft_data%init(ntimes, 1, 1, this%dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   !if (verbose>0) then
   !   fmtstring = '(A, I8, A, I8)'
   !   write(lu_out,fmtstring) '   ntimes: ',  ntimes,     '  , nfreq: ', this%nomega
   !end if

   allocate(seis_fd(nomega, 1))
   allocate(seis_filtered(ntimes_ft, 1))
   allocate(this%t(ntimes_ft))
   this%t = fft_data%get_t()


   ! FFT, timeshift and filter the seismogram
   call fft_data%rfft(taperandzeropad(seis_td, ntimes_ft), seis_fd)
 
   if (.not.deconv_stf) then
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
                       this%seis_cut )
   call cut_timewindow(this%t,             &
                       this%t,             &
                       this%time_window,   &
                       this%t_cut )

   ! Set length of kernel time window
   this%ntimes = size(this%seis_cut,1)

   ! Init FFT type for kernel time window, here for Parseval integration
   call this%fft_data%init(this%ntimes, 1, 1, this%dt)
   this%ntimes_ft = this%fft_data%get_ntimes()
   this%nomega    = this%fft_data%get_nomega()

   allocate(this%datat(this%ntimes,1))
   allocate(this%dataf(this%nomega,1))
   allocate(this%seis_cut_fd(this%nomega,1))

   ! FFT the cut and filtered seismogram
   call this%fft_data%rfft(taperandzeropad(array = reshape(this%seis_cut, [this%ntimes, 1]),  &
                                           ntimes = this%ntimes_ft, &
                                           ntaper = 0 ),            &
                            this%seis_cut_fd)

   ! Integrate over squared seismogram for normalization term in 5.13, TNM thesis
   normalization_term = this%integrate_parseval(this%seis_cut, this%seis_cut)                  
   if (normalization_term.lt.1.d-100) then
       this%normalization = 0
   else
       this%normalization = 1.d0 / normalization_term !/sum(this%seis_cut**2)
   end if
   
   if (verbose>0) then
      write(lu_out,*) '  Normalization coefficient: ', this%normalization
      write(lu_out,*) '  Length of seismogram: ', size(this%seis_cut,1), ' samples'
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,*) ''
   end if

   ! Write seismogram to disk (raw, filtered and cut to time window)
   if (firstslave.and.write_smgr) then
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
      do isample = 1, size(this%seis_cut,1)
         write(100,*) this%t_cut(isample), this%seis_cut(isample)
      end do
      close(100)
   end if
end subroutine cut_and_add_seismogram
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
logical function isinitialized(this)
   class(kernelspec_type)                   :: this

   isinitialized = this%initialized

end function isinitialized
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine freeme(this)
   class(kernelspec_type)                   :: this

   this%filter => NULL()
   this%initialized = .false.
   call this%fft_data%freeme()

end subroutine freeme
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function calc_misfit_kernel(this, timeseries, int_scheme)
   use simple_routines, only                : lowtrim
   !< This routine can take multiple time series and calculate the kernel at each
   class(kernelspec_type)                  :: this
   real(kind=dp), intent(in)               :: timeseries(:,:)
   character(len=*), intent(in)            :: int_scheme
   real(kind=dp)                           :: calc_misfit_kernel(size(timeseries,2))
   real(kind=dp), allocatable              :: timeseries_cut(:)
   integer                                 :: itrace, ntrace, lu_errorlog, it
   character(len=64)                       :: fmtstring, fnam_errorlog

   ntrace = size(timeseries,2)

   select case(trim(this%misfit_type))
   case('CC')

      do itrace = 1, ntrace
         
         call cut_timewindow(this%t,                 &
                             timeseries(:, itrace),  &
                             this%time_window,       &
                             timeseries_cut)

         select case(lowtrim(int_scheme))
         case('trapezoidal')
           calc_misfit_kernel(itrace) = - integrate_trapezoidal( timeseries_cut * this%seis_cut, &
                                                                 this%dt ) &
                                        * this%normalization
         case('parseval')
           calc_misfit_kernel(itrace) = - this%integrate_parseval( timeseries_cut, this%seis_cut_fd) &
                                        * this%normalization
         end select
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

   ! Handle the case that the Kernel is NaN
   if (any(calc_misfit_kernel.ne.calc_misfit_kernel)) then ! Kernel is NaN
       print *, 'ERROR Kernel has value NaN!'
       write(lu_out, *) 'ERROR Kernel has value NaN! Aborting calculation...'
       print *, 'Values of Kernel:'
       write(fmtstring, "('(',I6,'(ES11.3))')") ntrace
       print fmtstring, calc_misfit_kernel
       print *, 'Value of normalization:'
       print *, this%normalization
       print *, 'Convolved time series is written into "errorlog_timeseries",'
       print *, 'seismogram is written into "errorlog_seismogram"'
       do itrace = 1, ntrace
         if (calc_misfit_kernel(itrace).ne.calc_misfit_kernel(itrace)) then !is NaN
           call cut_timewindow(this%t,                 &
                               timeseries(:, itrace),  &
                               this%time_window,       &
                               timeseries_cut)

           calc_misfit_kernel(itrace) = integrate_trapezoidal( timeseries_cut * this%seis_cut, this%dt ) &
                                        * this%normalization

           write(fnam_errorlog,'(I06, "_errorlog_timeseries_full")') myrank                             
           open(newunit=lu_errorlog, file=fnam_errorlog, status='replace')
           do it = 1, size(this%t)
             write(lu_errorlog, '(ES11.3, ES15.8)') this%t(it), timeseries(it,itrace)
           end do
           close(lu_errorlog)

           write(fnam_errorlog,'(I06, "_errorlog_timeseries_cut")') myrank                             
           open(newunit=lu_errorlog, file=fnam_errorlog, status='replace')
           write(lu_errorlog, '(ES15.8)') timeseries_cut
           close(lu_errorlog)

           exit ! Only the first time series is written out, exit the loop
         end if
       end do
       write(fnam_errorlog,'(I06, "_errorlog_seismogram")') myrank                             
       open(newunit=lu_errorlog, file=fnam_errorlog, status='replace')
       write(lu_errorlog, *) this%seis_cut
       close(lu_errorlog)
       call pabort
   end if

end function calc_misfit_kernel
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function calc_basekernel(ibasekernel, strain_type_fwd, strain_type_bwd, &
                         fw_field_fd, bw_field_fd) result(conv_field_fd)

!< Calculates the basic (waveform) kernels from the forward and the backward field
!! in frequency domain.
!! readfields.f90 returns strain tensor in voigt notation
!! elements 1,...,6 contain tt,pp,rr,pr,tr,tp. See Fichtner
!! book, p. 168-169 how the fundamental kernels K_a, K_b and
!! K_c are defined

  integer,          intent(in)      :: ibasekernel
  character(len=*), intent(in)      :: strain_type_fwd, strain_type_bwd
  complex(kind=dp), intent(in)      :: fw_field_fd(:,:,:), bw_field_fd(:,:,:)
  complex(kind=dp)                  :: conv_field_fd(size(fw_field_fd,1),size(fw_field_fd,3))
  
  select case(ibasekernel)
  case(1) ! k_lambda
     select case(trim(strain_type_fwd))
     case('straintensor_trace')
        conv_field_fd = sum(fw_field_fd * bw_field_fd, 2)

     case('straintensor_full')
        select case(trim(strain_type_bwd))
        case('straintensor_trace') !This receiver has only trace base vectors
           conv_field_fd = (fw_field_fd(:,1,:) +    &
                            fw_field_fd(:,2,:) +    &
                            fw_field_fd(:,3,:)  ) * &
                            bw_field_fd(:,1,:) ! This row contains the trace
        case('straintensor_full')
           ! Only convolve trace of fw and bw fields
           conv_field_fd = (fw_field_fd(:,1,:) +    &
                            fw_field_fd(:,2,:) +    &
                            fw_field_fd(:,3,:)  ) * &
                           (bw_field_fd(:,1,:) +    &
                            bw_field_fd(:,2,:) +    &
                            bw_field_fd(:,3,:)  )
        end select
     end select
  case(2) ! k_mu
     conv_field_fd = sum(fw_field_fd * bw_field_fd, 2) * 2.d0
  case(3) ! k_rho
     print*,"ERROR: Density kernels not yet implemented"
     stop
  case(4) ! k_a
     conv_field_fd = ( bw_field_fd(:,2,:) + bw_field_fd(:,1,:) ) * &
                     ( fw_field_fd(:,2,:) * bw_field_fd(:,1,:) )
  case(5) ! k_b
     conv_field_fd = ( ( bw_field_fd(:,5,:) * fw_field_fd(:,5,:) ) + &
                     ( bw_field_fd(:,4,:) * fw_field_fd(:,4,:) ) ) * 4.d0
  case(6) ! k_c
     conv_field_fd = ( fw_field_fd(:,1,:) + fw_field_fd(:,2,:) ) * &
                       bw_field_fd(:,3,:) + &
                     ( bw_field_fd(:,1,:) * bw_field_fd(:,2,:) ) * &
                       fw_field_fd(:,3,:)
  end select

end function calc_basekernel
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Here the we tabulate the base kernels required to assemble 
!! the physical kernels for a desired model parameter
!! TODO: Move this to own module, together with assemble_basekernel, but wait for 
!!       merge with Ludwigs branch - SCS 240914
subroutine tabulate_kernels(model_param, needs_basekernel, strain_type)
  use simple_routines, only       : to_lower
  character(len=*),  intent(in)  :: model_param
  logical,           intent(out) :: needs_basekernel(6)
  character(len=32), intent(out) :: strain_type

  needs_basekernel = .false.
  select case(trim(to_lower(model_param)))
  case('lam')
     needs_basekernel(1)   = .true.   ! Lambda
     strain_type      = 'straintensor_trace'

  case('vp')
     needs_basekernel(1)   = .true.   ! Lambda
     strain_type      = 'straintensor_trace'

  case('vs')
     needs_basekernel(1:2) = .true.   ! Lambda, mu
     strain_type      = 'straintensor_full'
  
  case('mu')
     needs_basekernel(2)   = .true.   ! mu
     strain_type      = 'straintensor_full'
  
  case('vsh')
     needs_basekernel([1,2,5,6]) = .true. ! Lambda, mu, b, c
     strain_type      = 'straintensor_full'
  
  case('vsv')
     needs_basekernel(5)   = .true.   ! b
     strain_type      = 'straintensor_full'
  
  case('vph')
     needs_basekernel([4,6]) = .true.   ! a, c
     strain_type      = 'straintensor_full'

  case('vpv')
     needs_basekernel([1,4,6]) = .true.! Lambda, a, c
     strain_type      = 'straintensor_full'

  case('eta')
     needs_basekernel([1,4,6]) = .true. !Lambda, a, c
     strain_type      = 'straintensor_full'

  case('rho')
     needs_basekernel(1:6) = .true. !Lambda, mu, rho, a, b, c
     strain_type      = 'straintensor_full'

  case default
     write(*,*) "Error in tabulate_kernels: Unknown model parameter: ", trim(model_param)
     call pabort      
  end select
end subroutine tabulate_kernels
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function calc_physical_kernels(model_param, base_kernel, bg_model, relative_kernel)  &
                               result(physical_kernel)

  use backgroundmodel, only               :  backgroundmodel_type
  character(len=*), intent(in)           :: model_param
  real(kind=dp),    intent(in)           :: base_kernel(:,:)
  type(backgroundmodel_type), intent(in) :: bg_model
  logical, intent(in)                    :: relative_kernel
  real(kind=dp)                          :: physical_kernel(size(base_kernel,1))

  ! See Fichtner p. 169 how kernels are defined and type_parameter.f90
  ! how basekernels inside kernelvalue_weighted are ordered
  ! base_kernel contains the following order: 
  ! lambda, mu, rho, a, b, c
  select case(trim(model_param))
  case('lam')
     physical_kernel(:) =                                         &
            base_kernel(:, 1)                                ! lambda
     if (relative_kernel) physical_kernel = physical_kernel * bg_model%c_lam

  case('mu')
     physical_kernel(:) =                                         & 
            base_kernel(:, 2)                                ! Mu
     if (relative_kernel) physical_kernel = physical_kernel * bg_model%c_mu 

  case('vp')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vp * &
            base_kernel(:, 1)                                ! Lambda
     if (relative_kernel) physical_kernel = physical_kernel * bg_model%c_vp 

  case('vs')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vs * &
          ( base_kernel(:, 2)                             &  ! Mu
          - base_kernel(:, 1) * 2.d0)                        ! Lambda
     if (relative_kernel) physical_kernel = physical_kernel * bg_model%c_vs 

  case('vsh')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vsh * &
          (  2 * base_kernel(:, 1)                        &  ! Lambda
           +     base_kernel(:, 2)                        &  ! Mu
           +     base_kernel(:, 5)                        &  ! B
           + 2 * base_kernel(:, 6) * (1-bg_model%c_eta))     ! C 
     if (relative_kernel) physical_kernel = physical_kernel * bg_model%c_vsh

  case('vsv')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vsv * &
            base_kernel(:, 5)                                ! B
     if (relative_kernel) physical_kernel = physical_kernel * bg_model%c_vsv

  case('vph')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vph * &
          ( base_kernel(:, 4) +                           &  ! A
            base_kernel(:, 6) * bg_model%c_eta )             ! C
     if (relative_kernel) physical_kernel = physical_kernel * bg_model%c_vph

  case('vpv')
     physical_kernel(:) = 2.d0 * bg_model%c_rho * bg_model%c_vpv * &
          ( base_kernel(:, 1) +                           &  ! Lambda
            base_kernel(:, 4) +                           &  ! A
            base_kernel(:, 6) )                              ! C 
     if (relative_kernel) physical_kernel = physical_kernel * bg_model%c_vpv

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
   

end subroutine cut_timewindow
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
pure function integrate_trapezoidal(timeseries, dt) result(integrate)
    real(kind=dp), intent(in)   :: timeseries(:), dt
    real(kind=dp)               :: integrate
    integer                     :: npoints

    npoints = size(timeseries, 1)

    ! Trapezoidal rule: I = dt/2 * (f(x1) + 2f(x2) + ... + 2f(xN-1) + f(xN))
    integrate = timeseries(1) + timeseries(npoints) + 2*sum(timeseries(2:npoints-1), 1)
    integrate = integrate * dt * 0.5

end function integrate_trapezoidal
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function integrate_parseval_both_real(this, a, b) result(integrate)
    class(kernelspec_type)                        :: this
    real(kind=dp), intent(in)                     :: a(:), b(:)
    real(kind=dp)                                 :: integrate
    complex(kind=dp), dimension(this%nomega,1)    :: af, bf, atimesb

    this%datat(:,1) = a
    call this%fft_data%rfft(taperandzeropad(array = this%datat,      &
                                            ntimes = this%ntimes_ft, &
                                            ntaper = 0 ),            &
                            af)
    
    this%datat(:,1) = b
    call this%fft_data%rfft(taperandzeropad(array = this%datat,      &
                                            ntimes = this%ntimes_ft, &
                                            ntaper = 0 ),            &
                            bf)

    atimesb = af * conjg(bf)                      

    ! Integrate by Trapezoidal rule in frequency domain
    integrate = real((atimesb(1, 1) +                          &
                      atimesb(this%nomega, 1) +                &
                      2 * sum(atimesb(2:this%nomega-1, 1)))  * &
                     this%fft_data%get_df() , kind=dp)

end function integrate_parseval_both_real
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function integrate_parseval_b_complex(this, a, bf) result(integrate)
    class(kernelspec_type)                        :: this
    real(kind=dp), intent(in)                     :: a(:)
    complex(kind=dp), intent(in)                  :: bf(:,:)
    real(kind=dp)                                 :: integrate
    complex(kind=dp), dimension(this%nomega,1)    :: af, atimesb

    this%datat(:,1) = a
    call this%fft_data%rfft(taperandzeropad(array = this%datat,      &
                                            ntimes = this%ntimes_ft, &
                                            ntaper = 0 ),            &
                            af)
    
    atimesb = af * conjg(bf)                      

    ! Integrate by Trapezoidal rule in frequency domain
    integrate = real((atimesb(1, 1) +                          &
                      atimesb(this%nomega, 1) +                &
                      2 * sum(atimesb(2:this%nomega-1, 1)))  * &
                     this%fft_data%get_df() , kind=dp)

end function integrate_parseval_b_complex
!-------------------------------------------------------------------------------

end module kernel
!=========================================================================================
