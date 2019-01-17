!******************************************************************************
!
!    This file is part of:
!    MC Kernel: Calculating seismic sensitivity kernels on unstructured meshes
!    Copyright (C) 2016 Simon Staehler, Martin van Driel, Ludwig Auer
!
!    You can find the latest version of the software at:
!    <https://www.github.com/tomography/mckernel>
!
!    MC Kernel is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MC Kernel is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MC Kernel. If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

!=========================================================================================
module kernel
use global_parameters,                    only: sp, dp, verbose, lu_out, firstslave, myrank, testing
use filtering,                            only: filter_type
use fft,                                  only: rfft_type, taperandzeropad
use commpi,                               only: pabort
implicit none

  type kernelspec_type
    private
    character(len=32), public            :: name
    real(kind=dp), dimension(2)          :: time_window
    real(kind=dp), allocatable           :: seis_disp_cut(:), &      ! Seismogram (Velocity or displacement)
                                            seis_velo_cut(:)         ! in the time window of 
                                                                     ! this kernel
    complex(kind=dp), allocatable        :: seis_disp_cut_fd(:,:), & ! Seismogram (Velocity or displacement)
                                            seis_velo_cut_fd(:,:)    ! in the time window of 
                                                                     ! this kernel in frequency domain
    real(kind=dp), allocatable           :: t(:)
    real(kind=dp), allocatable, public   :: t_cut(:)
    real(kind=sp), allocatable, public   :: seis(:)
    real(kind=dp)                        :: dt
    real(kind=dp), public                :: normalization
    character(len=4)                     :: misfit_type
    character(len=16), public            :: model_parameter
    integer,           public            :: model_parameter_index
    integer,           public            :: hetero_parameter_index
    character(len=32)                    :: strain_type 
    type(filter_type), public, pointer   :: filter
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

  interface calc_physical_kernels
    module procedure                     :: calc_physical_kernels_scalar
    module procedure                     :: calc_physical_kernels_time_series
  end interface calc_physical_kernels



contains 

!-------------------------------------------------------------------------------
subroutine init(this, name, time_window, filter, misfit_type, model_parameter, &
                seis, dt, timeshift_fwd, deconv_stf, write_smgr)
   use background_model, only               : get_parameter_index
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
   
   logical                                 :: deconv_stf_loc, write_smgr_loc
   real(kind=dp)                           :: timeshift_fwd_loc

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
   use simple_routines,               only  : check_nan, cumsum_trapezoidal
   class(kernelspec_type)                  :: this
   real(kind=dp), intent(in)               :: seis(:)
   logical, intent(in)                     :: deconv_stf
   logical, intent(in)                     :: write_smgr
   real(kind=dp), intent(in)               :: timeshift_fwd

   type(rfft_type)                         :: fft_data
   type(timeshift_type)                    :: timeshift
   real(kind=dp),    allocatable           :: seis_disp_td(:,:)
   real(kind=dp),    allocatable           :: seis_velo_td(:,:)
   complex(kind=dp), allocatable           :: seis_disp_fd(:,:)
   complex(kind=dp), allocatable           :: seis_velo_fd(:,:)
   complex(kind=dp), allocatable           :: seis_disp_filtered_fd(:,:)
   complex(kind=dp), allocatable           :: seis_velo_filtered_fd(:,:)
   real(kind=dp),    allocatable           :: seis_disp_filtered(:,:)
   real(kind=dp),    allocatable           :: seis_velo_filtered(:,:), f(:)
   real(kind=dp)                           :: normalization_term
   integer                                 :: ntimes, ntimes_ft, nomega, isample, nan_loc(2)
   logical                                 :: isnan 

   ! Save seismogram in the timewindow of this kernel

   ! Allocate temporary variable
   ntimes = size(seis, 1)
   allocate(seis_disp_td(ntimes, 1))
   allocate(seis_velo_td(ntimes, 1))
   seis_disp_td(:,1) = cumsum_trapezoidal(seis, this%dt)
   seis_velo_td(:,1) = seis


   ! Initialize FFT type for the seismogram 
   call fft_data%init(ntimes, 1, 1, this%dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   allocate(seis_disp_fd(nomega, 1))
   allocate(seis_velo_fd(nomega, 1))
   allocate(seis_disp_filtered_fd(nomega, 1))
   allocate(seis_velo_filtered_fd(nomega, 1))
   allocate(seis_disp_filtered(ntimes_ft, 1))
   allocate(seis_velo_filtered(ntimes_ft, 1))
   allocate(this%seis(ntimes_ft))
   allocate(this%t(ntimes_ft))
   this%t = fft_data%get_t()

   ! FFT the seismograms
   call fft_data%rfft(taperandzeropad(seis_disp_td, ntimes_ft, ntaper=0), seis_disp_fd)
   call fft_data%rfft(taperandzeropad(seis_velo_td, ntimes_ft, ntaper=0), seis_velo_fd)

   if (.not.deconv_stf) then
     ! It's slightly ineffective to init that every time, but it is only 
     ! called once per receiver at initialization.
     call timeshift%init_ts(fft_data%get_f(), dtshift = timeshift_fwd)
     call timeshift%apply(seis_disp_fd)
     call timeshift%apply(seis_velo_fd)
     call timeshift%freeme()
   end if

   seis_disp_filtered_fd = this%filter%apply_2d(seis_disp_fd, kind='fwd')
   seis_velo_filtered_fd = this%filter%apply_2d(seis_velo_fd, kind='fwd')

   call fft_data%irfft(seis_disp_filtered_fd, seis_disp_filtered)
   call fft_data%irfft(seis_velo_filtered_fd, seis_velo_filtered)

   call check_NaN(seis_velo_filtered, isnan, nan_loc)
   if (isnan) then
     print *, myrank, ': ERROR: NaN found in velocity seismogram:'
     print *, 'NaN index: ', nan_loc
     print *, 'Wrote seismogram to seis_velo_nan.txt'
     open(unit=100, file='seis_velo_nan.txt', action='write')
     do isample = 1, size(seis_velo_filtered,1)
         write(100,*) this%t(isample), seis_velo_filtered(isample,1)
     end do
     close(100)
     open(unit=100, file='seis_velo_nan_fd.txt', action='write')
     f = fft_data%get_f()
     do isample = 1, size(seis_velo_fd,1)
         write(100,*) f(isample), seis_velo_fd(isample,1)
     end do
     close(100)
     stop
   end if

   call check_NaN(seis_disp_filtered, isnan, nan_loc)
   if (isnan) then
     print *, myrank, ': ERROR: NaN found in displacement seismogram:'
     print *, 'NaN index: ', nan_loc
     print *, 'Wrote seismogram to seis_disp_nan.txt'
     open(unit=100, file='seis_disp_nan.txt', action='write')
     do isample = 1, size(seis_disp_filtered,1)
         write(100,*) this%t(isample), seis_disp_filtered(isample,1)
     end do
     close(100)
     open(unit=100, file='seis_disp_nan_fd.txt', action='write')
     f = fft_data%get_f()
     do isample = 1, size(seis_disp_fd,1)
         write(100,*) f(isample), seis_disp_fd(isample,1)
     end do
     close(100)
     stop
   end if

   call fft_data%freeme()

   ! Cut timewindow of the kernel from the seismogram
   call cut_timewindow(this%t,             &
                       seis_disp_filtered(:,1), &
                       this%time_window,   &
                       this%seis_disp_cut )
   call cut_timewindow(this%t,             &
                       seis_velo_filtered(:,1), &
                       this%time_window,   &
                       this%seis_velo_cut )
   call cut_timewindow(this%t,             &
                       this%t,             &
                       this%time_window,   &
                       this%t_cut )

   ! Set length of kernel time window
   this%ntimes = size(this%t_cut,1)

   ! Init FFT type for kernel time window, here for Parseval integration, but 
   ! will be used later for integration over waveform kernels
   call this%fft_data%init(this%ntimes, 1, 1, this%dt)
   this%ntimes_ft = this%fft_data%get_ntimes()
   this%nomega    = this%fft_data%get_nomega()

   allocate(this%datat(this%ntimes_ft,1))
   allocate(this%dataf(this%nomega,1))
   allocate(this%seis_disp_cut_fd(this%nomega,1))
   allocate(this%seis_velo_cut_fd(this%nomega,1))

   ! FFT the cut and filtered seismogram
   call this%fft_data%rfft(taperandzeropad(array = reshape(this%seis_disp_cut, [this%ntimes, 1]),  &
                                           ntimes = this%ntimes_ft, &
                                           ntaper = 0 ),            &
                            this%seis_disp_cut_fd)
   call this%fft_data%rfft(taperandzeropad(array = reshape(this%seis_velo_cut, [this%ntimes, 1]),  &
                                           ntimes = this%ntimes_ft, &
                                           ntaper = 0 ),            &
                            this%seis_velo_cut_fd)


   ! Integrate over squared seismogram for normalization term in 5.13, TNM thesis
   select case(this%misfit_type)
   case('CC')
     ! For travel-time kernels, the velocity seismogram has to be used
     normalization_term = this%integrate_parseval(this%seis_velo_cut, this%seis_velo_cut)                  
     ! Save the complete seismogram once. This is only used for waveform kernel 
     ! plotting, not for kernel calculation.
     this%seis = real(seis_velo_filtered(:,1), kind=sp)
   
   case('AM')
     ! For travel-time kernels, the velocity seismogram has to be used
     normalization_term = this%integrate_parseval(this%seis_disp_cut, this%seis_disp_cut)                  
     ! Save the complete seismogram once. This is only used for waveform kernel 
     ! plotting, not for kernel calculation.
     this%seis = real(seis_disp_filtered(:,1), kind=sp)

   case default
     print *, 'Unkown misfit type: ', trim(this%misfit_type)
     call pabort()

   end select

   if (normalization_term.gt.1.d-100) then
       this%normalization = 1.d0 / normalization_term 
   else
       this%normalization = 0
   end if
   
   if (verbose>0) then
      write(lu_out,*) '  Normalization coefficient: ', this%normalization
      write(lu_out,*) '  Length of seismogram: ', size(this%t_cut,1), ' samples'
      write(lu_out,*) '  ---------------------------------------------------------'
      write(lu_out,*) ''
   end if

   ! Write seismogram to disk (raw, filtered and cut to time window)
   ! First column: time
   ! 2nd column: displacement seismogram in meters
   ! 3rd column: velocity seismogram in meters/sec
   if ((firstslave.and.write_smgr).or.testing) then
      open(unit=100,file='./Seismograms/seism_raw_'//trim(this%name), action='write')
      do isample = 1, size(seis,1)
         write(100,*) this%t(isample), seis(isample)
      end do
      close(100)

      open(unit=100,file='./Seismograms/seism_'//trim(this%name), action='write')
      do isample = 1, size(seis,1)
         write(100,*) this%t(isample), seis_disp_filtered(isample,1), seis_velo_filtered(isample,1)
      end do
      close(100)

      open(unit=100,file='./Seismograms/seism_cut_'//trim(this%name), action='write')
      do isample = 1, size(this%t_cut,1)
         write(100,*) this%t_cut(isample), this%seis_disp_cut(isample), this%seis_velo_cut(isample)
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
   integer                                 :: itrace, ntrace !, lu_errorlog, it
   !character(len=64)                       :: fmtstring, fnam_errorlog

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
           calc_misfit_kernel(itrace) = integrate_trapezoidal( timeseries_cut * this%seis_velo_cut, &
                                                                 this%dt ) &
                                        * this%normalization * -1.d0
         case('parseval')
           calc_misfit_kernel(itrace) = this%integrate_parseval( timeseries_cut, this%seis_velo_cut_fd) &
                                        * this%normalization * -1.d0
         end select
      end do

   case('AM')

      do itrace = 1, ntrace
         
         call cut_timewindow(this%t,                 &
                             timeseries(:, itrace),  &
                             this%time_window,       &
                             timeseries_cut)

         select case(lowtrim(int_scheme))
         case('trapezoidal')
           calc_misfit_kernel(itrace) = integrate_trapezoidal( timeseries_cut * this%seis_disp_cut, &
                                                                 this%dt ) &
                                        * this%normalization * -1.d0
         case('parseval')
           calc_misfit_kernel(itrace) = this%integrate_parseval( timeseries_cut, this%seis_disp_cut_fd) &
                                        * this%normalization * -1.d0
         end select
      end do

   end select

   ! Handle the case that the Kernel is NaN
   !if (any(calc_misfit_kernel.ne.calc_misfit_kernel)) then ! Kernel is NaN
   !    print *, 'ERROR Kernel has value NaN!'
   !    write(lu_out, *) 'ERROR Kernel has value NaN! Aborting calculation...'
   !    print *, 'Values of Kernel:'
   !    write(fmtstring, "('(',I6,'(ES11.3))')") ntrace
   !    print fmtstring, calc_misfit_kernel
   !    print *, 'Value of normalization:'
   !    print *, this%normalization
   !    print *, 'Convolved time series is written into "errorlog_timeseries",'
   !    print *, 'seismogram is written into "errorlog_seismogram"'
   !    do itrace = 1, ntrace
   !      if (calc_misfit_kernel(itrace).ne.calc_misfit_kernel(itrace)) then !is NaN
   !        call cut_timewindow(this%t,                 &
   !                            timeseries(:, itrace),  &
   !                            this%time_window,       &
   !                            timeseries_cut)

   !        calc_misfit_kernel(itrace) = integrate_trapezoidal( timeseries_cut * this%seis_cut, this%dt ) &
   !                                     * this%normalization

   !        write(fnam_errorlog,'(I06, "_errorlog_timeseries_full")') myrank                             
   !        open(newunit=lu_errorlog, file=fnam_errorlog, status='replace')
   !        do it = 1, size(this%t)
   !          write(lu_errorlog, '(ES11.3, ES15.8)') this%t(it), timeseries(it,itrace)
   !        end do
   !        close(lu_errorlog)

   !        write(fnam_errorlog,'(I06, "_errorlog_timeseries_cut")') myrank                             
   !        open(newunit=lu_errorlog, file=fnam_errorlog, status='replace')
   !        write(lu_errorlog, '(ES15.8)') timeseries_cut
   !        close(lu_errorlog)

   !        exit ! Only the first time series is written out, exit the loop
   !      end if
   !    end do
   !    write(fnam_errorlog,'(I06, "_errorlog_seismogram")') myrank                             
   !    open(newunit=lu_errorlog, file=fnam_errorlog, status='replace')
   !    write(lu_errorlog, *) this%seis_cut
   !    close(lu_errorlog)
   !    call pabort
   !end if

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
        conv_field_fd = - sum(fw_field_fd * bw_field_fd, 2)

     case('straintensor_full')
        select case(trim(strain_type_bwd))
        case('straintensor_trace') !This receiver has only trace base vectors
           conv_field_fd = - (fw_field_fd(:,1,:) +    &
                            fw_field_fd(:,2,:) +    &
                            fw_field_fd(:,3,:)  ) * &
                            bw_field_fd(:,1,:) ! This row contains the trace
        case('straintensor_full')
           ! Only convolve trace of fw and bw fields
           conv_field_fd =  - (fw_field_fd(:,1,:) +    &
                            fw_field_fd(:,2,:) +    &
                            fw_field_fd(:,3,:)  ) * &
                           (bw_field_fd(:,1,:) +    &
                            bw_field_fd(:,2,:) +    &
                            bw_field_fd(:,3,:)  )
        end select
     end select
  case(2) ! k_mu
     conv_field_fd = - sum(fw_field_fd * bw_field_fd, 2) * 2.d0
  case(3) ! k_rho
     print*,"ERROR: Density kernels not yet implemented"
     stop
  case(4) ! k_a
     conv_field_fd = ( fw_field_fd(:,2,:) + fw_field_fd(:,1,:) ) * &
                     ( bw_field_fd(:,2,:) + bw_field_fd(:,1,:) ) 
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
!> Just a wrapper for calc_physical_kernels_time_series
function calc_physical_kernels_scalar(model_param, base_kernel, bg_model, &
                                      relative_kernel)  &
                                      result(physical_kernel)

  use background_model, only              :  backgroundmodel_type
  character(len=*), intent(in)           :: model_param
  real(kind=dp),    intent(in)           :: base_kernel(:,:)
  real(kind=dp), allocatable             :: temp_in(:,:,:)
  real(kind=dp)                          :: temp_out(1,size(base_kernel,1))
  type(backgroundmodel_type), intent(in) :: bg_model
  logical, intent(in)                    :: relative_kernel
  real(kind=dp)                          :: physical_kernel(size(base_kernel,1))

  if (size(base_kernel,2).ne.6) then
    print *, 'Error in calc_physical_kernels: base_kernel has wrong size:'
    print *, 'size(base_kernel): ', size(base_kernel, 1), size(base_kernel, 2)
  end if

  allocate(temp_in(1, size(base_kernel, 1), size(base_kernel, 2)))
  temp_in(1,:,:) = base_kernel(:,:)

  temp_out = calc_physical_kernels_time_series(model_param, temp_in, bg_model, &
                                               relative_kernel)
  physical_kernel = temp_out(1,:)

end function calc_physical_kernels_scalar
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Calculate physical kernels from the base kernels. An adaptation of eq. 9.21 
!! in the Fichtner book.
!! Works on whole seismograms, along the first dimension of base_kernel
function calc_physical_kernels_time_series(model_param, base_kernel, bg_model, &
                                           relative_kernel)  &
                                           result(physical_kernel)

  use background_model, only              :  backgroundmodel_type
  character(len=*), intent(in)           :: model_param
  real(kind=dp),    intent(in)           :: base_kernel(:,:,:)
  type(backgroundmodel_type), intent(in) :: bg_model
  logical, intent(in)                    :: relative_kernel
  real(kind=dp)                          :: physical_kernel(size(base_kernel, 1), &
                                                            size(base_kernel, 2))
  integer                                :: it, nt

  nt = size(base_kernel, 1)
  if (size(base_kernel,3).ne.6) then
    print *, 'Error in calc_physical_kernels: base_kernel has wrong size:'
    print *, 'size(base_kernel): ', size(base_kernel, 1), size(base_kernel, 2)
  end if

  ! See Fichtner p. 169 how kernels are defined and type_parameter.f90
  ! how basekernels inside kernelvalue_weighted are ordered
  ! base_kernel contains the following order: 
  ! lambda, mu, rho, a, b, c
  select case(trim(model_param))
  case('lam')
    if (relative_kernel) then
      do it = 1, nt
        physical_kernel(it, :) =                                         &
               base_kernel(it, :, 1)                            &    ! lambda
             * bg_model%c_lam(:)  
      end do
    else
      do it = 1, nt
        physical_kernel(it, :) =                                &         
               base_kernel(it, :, 1)                                ! lambda
      end do
    end if

  case('mu')
    if (relative_kernel) then
      do it = 1, nt
        physical_kernel(it, :) =                                         &
               base_kernel(it, :, 2)                            &    ! Mu
             * bg_model%c_mu(:)  
      end do
    else
      do it = 1, nt
        physical_kernel(it, :) =                                &         
               base_kernel(it, :, 2)                                ! Mu
      end do
    end if

  case('vp')
    if (relative_kernel) then
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vp**2 * &
               base_kernel(it, :, 1)                                ! Lambda
      end do
    else
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vp * &
               base_kernel(it, :, 1)                                ! Lambda
      end do
    end if

  case('vs')
    if (relative_kernel) then
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vs**2 * &
             (      base_kernel(it, :, 2)                        &  ! Mu
              - 2 * base_kernel(it, :, 1))                          ! Lambda
      end do
    else
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vs * &
             (      base_kernel(it, :, 2)                        &  ! Mu
              - 2 * base_kernel(it, :, 1))                          ! Lambda
      end do
    end if

  case('vsh')
    if (relative_kernel) then
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vsh**2 * &
             (      base_kernel(it, :, 2)                        &  ! Mu
              - 2 * base_kernel(it, :, 1)                        &  ! Lambda
              -     base_kernel(it, :, 5)                        &  ! B
              + 2 * base_kernel(it, :, 6) * (1-bg_model%c_eta))     ! C 
      end do
    else
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vsh * &
             (      base_kernel(it, :, 2)                        &  ! Mu
              - 2 * base_kernel(it, :, 1)                        &  ! Lambda
              -     base_kernel(it, :, 5)                        &  ! B
              + 2 * base_kernel(it, :, 6) * (1-bg_model%c_eta))     ! C 
      end do
    end if

  case('vsv')
    if (relative_kernel) then
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vsv**2 * &
               base_kernel(it, :, 5)                                ! B
      end do
    else
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vsv * &
               base_kernel(it, :, 5)                                ! B
      end do
    end if

  case('vph')
    if (relative_kernel) then
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vph**2 * &
             (  base_kernel(it, :, 4)                            &  ! A
              + base_kernel(it, :, 6) * bg_model%c_eta )            ! C
      end do
    else
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vph * &
             (  base_kernel(it, :, 4)                            &  ! A
              + base_kernel(it, :, 6) * bg_model%c_eta )            ! C
      end do
    end if

  case('vpv')
    if (relative_kernel) then
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vpv**2 * &
             (  base_kernel(it, :, 1)                            &  ! Lambda
              - base_kernel(it, :, 4)                            &  ! A
              - base_kernel(it, :, 6) )                             ! C 
      end do
    else
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho * bg_model%c_vpv * &
             (  base_kernel(it, :, 1)                            &  ! Lambda
              - base_kernel(it, :, 4)                            &  ! A
              - base_kernel(it, :, 6) )                             ! C 
      end do
    end if

  case('eta')
    if (relative_kernel) then
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho *         &
               (bg_model%c_vph**2 - bg_model%c_vsh**2) *         &
               base_kernel(it, :, 6) * bg_model%c_eta              ! C 
      end do
    else
      do it = 1, nt
        physical_kernel(it, :) = 2.d0 * bg_model%c_rho *         &
               (bg_model%c_vph**2 - bg_model%c_vsh**2) *         &
               base_kernel(it, :, 6)                               ! C 
      end do
    end if

  case('rho')
    write(*,*) "Error: Density kernels not yet implemented"
    stop

  case default
    write(*,*) "Error: Unknown model parameter: ", trim(model_param)
  end select

end function calc_physical_kernels_time_series
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


   apply_filter_3d = this%filter%apply_3d(freq_series, kind='fwd')
   apply_filter_3d = this%filter%apply_3d(apply_filter_3d, kind='bwd')

end function apply_filter_3d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function apply_filter_2d(this, freq_series)
   class(kernelspec_type)           :: this
   complex(kind=dp), intent(in)     :: freq_series(:,:)
   complex(kind=dp)                 :: apply_filter_2d(size(freq_series,1), &
                                                       size(freq_series,2))


   apply_filter_2d = this%filter%apply_2d(freq_series, kind='fwd')
   apply_filter_2d = this%filter%apply_2d(apply_filter_2d, kind='bwd')

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

    this%datat(:,1) = 0
    this%datat(1:size(a),1) = a
    call this%fft_data%rfft(this%datat, af)
    
    this%datat(:,1) = 0
    this%datat(1:size(b),1) = b
    call this%fft_data%rfft(this%datat, bf)

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

    this%datat(:,1) = 0
    this%datat(1:size(a),1) = a
    call this%fft_data%rfft(this%datat, af)
    
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
