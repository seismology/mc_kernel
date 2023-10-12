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
module filtering
   
   use global_parameters
   use simple_routines,  only           : mult2d_1d, mult3d_1d
   use commpi,           only           : pabort
   
   implicit none
   private

   public   filter_type, timeshift_type

   type filter_type
       private
       character(len=32), public       :: name
       complex(kind=dp), allocatable   :: transferfunction(:)
       complex(kind=dp), allocatable   :: transferfunction_fwd(:)
       complex(kind=dp), allocatable   :: transferfunction_bwd(:)
       integer                         :: nfreq                !< Number of frequencies
       real(kind=dp), allocatable      :: f(:)
       logical                         :: initialized = .false.
       logical                         :: stf_added = .false.
       character(len=32), public       :: filterclass
       real(kind=dp),     public       :: frequencies(4)

       contains
       procedure, pass   :: create
       procedure, pass   :: add_stfs
       procedure, pass   :: deleteme
       procedure, pass   :: apply_1d
       procedure, pass   :: apply_2d
       procedure, pass   :: apply_3d
       procedure, pass   :: isinitialized
       procedure, pass   :: get_transferfunction
       generic           :: apply => apply_1d, apply_2d, apply_3d
   end type

   type timeshift_type
       complex(kind=dp), allocatable  :: shift_fd(:) !< Array with complex multiplicators 
                                                     !! to apply the timeshift
       logical                        :: isinitialized = .false.
       contains
       procedure, pass   :: init_ts
       procedure, pass   :: timeshift_1d
       procedure, pass   :: timeshift_md
       procedure, pass   :: freeme
       generic           :: apply => timeshift_1d, timeshift_md
   end type

contains

!-----------------------------------------------------------------------------------------
!> Create a filter object with the specified parameters
subroutine create(this, name, dfreq, nfreq, filterclass, frequencies, norder)

    class(filter_type)              :: this
    integer, intent(in)             :: nfreq
    real(kind=dp), intent(in)       :: dfreq, frequencies(4)
    integer                         :: ifreq
    character(len=32), intent(in)   :: name
    character(len=32), intent(in)   :: filterclass
    integer                         :: ifreq

    character(len=64)               :: fmtstring
    character(len=64)               :: fnam

    fmtstring = '("Class: ", A, ", freq: ", 4(F8.4))'

    if (this%initialized) then
       write(*,*) 'ERROR: This filter is already initialized as'
       write(*, fmtstring) trim(this%filterclass), this%frequencies
       write(*,*) 'delete it first before using'
       call pabort
    end if
    allocate(this%transferfunction(nfreq))
    allocate(this%transferfunction_fwd(nfreq))
    allocate(this%transferfunction_bwd(nfreq))
    allocate(this%f(nfreq))

    this%name              = name
    this%transferfunction  = 0.0
    this%nfreq             = nfreq
    this%filterclass       = filterclass
    this%frequencies       = frequencies

    do ifreq = 1, nfreq
       this%f(ifreq) = dfreq * (ifreq-1)
    end do

    select case(trim(filterclass))
    case('ident')
       this%transferfunction = ident(this%f)

    case('Gabor')
       ! For the gabor filter, the three parameters are 
       ! 1. center period
       ! 2. sigma divided by the center period
       ! 3. time shift in seconds
       this%transferfunction = loggabor(this%f, frequencies(1), frequencies(2), frequencies(3))

    case('Butterw_LP')
       this%transferfunction = butterworth_lowpass(this%f, frequencies(1), int(frequencies(3)))
    case('Butterw_HP')
       this%transferfunction = butterworth_highpass(this%f, frequencies(1), int(frequencies(3)))
    case('Butterw_BP')
       this%transferfunction = butterworth_highpass(this%f, frequencies(1), int(frequencies(3))) &
                             * butterworth_lowpass(this%f, frequencies(2), int(frequencies(3)))

    case('Butterw_LP_ZP')
       this%transferfunction = butterworth_lowpass(this%f, frequencies(1), int(frequencies(3)))
       this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    case('Butterw_HP_ZP')
       this%transferfunction = butterworth_highpass(this%f, frequencies(1), int(frequencies(3)))
       this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    case('Butterw_BP_ZP')
       this%transferfunction = butterworth_highpass(this%f, frequencies(1), int(frequencies(3))) &
                             * butterworth_lowpass(this%f, frequencies(2), int(frequencies(3)))
       this%transferfunction = this%transferfunction * conjg(this%transferfunction)

    !  case('Butterw_LP_O4')
    !     this%transferfunction = butterworth_lowpass(this%f, frequencies(1), 4)
    !  case('Butterw_HP_O4')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 4)
    !  case('Butterw_BP_O4')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 4) &
    !                           * butterworth_lowpass(this%f, frequencies(2), 4)
    !  case('Butterw_LP_O6')
    !     this%transferfunction = butterworth_lowpass(this%f, frequencies(1), 6)
    !  case('Butterw_HP_O6')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 6)
    !  case('Butterw_BP_O6')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 6) &
    !                           * butterworth_lowpass(this%f, frequencies(2), 6)
    !  case('Butterw_LP_O8')
    !     this%transferfunction = butterworth_lowpass(this%f, frequencies(1), 8)
    !  case('Butterw_HP_O8')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 8)
    !  case('Butterw_BP_O8')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 8) &
    !                           * butterworth_lowpass(this%f, frequencies(2), 8)

    !  case('Butterw_LP_O2_ZP')
    !     this%transferfunction = butterworth_lowpass(this%f, frequencies(1), 2)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_HP_O2_ZP')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 2)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_BP_O2_ZP')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 2) &
    !                           * butterworth_lowpass(this%f, frequencies(2), 2)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_LP_O4_ZP')
    !     this%transferfunction = butterworth_lowpass(this%f, frequencies(1), 4)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_HP_O4_ZP')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 4)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_BP_O4_ZP')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 4) &
    !                           * butterworth_lowpass(this%f, frequencies(2), 4)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_LP_O6_ZP')
    !     this%transferfunction = butterworth_lowpass(this%f, frequencies(1), 6)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_HP_O6_ZP')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 6)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_BP_O6_ZP')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 6) &
    !                           * butterworth_lowpass(this%f, frequencies(2), 6)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_LP_O8_ZP')
    !     this%transferfunction = butterworth_lowpass(this%f, frequencies(1), 8)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_HP_O8_ZP')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 8)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    !  case('Butterw_BP_O8_ZP')
    !     this%transferfunction = butterworth_highpass(this%f, frequencies(1), 8) &
    !                           * butterworth_lowpass(this%f, frequencies(2), 8)
    !     this%transferfunction = this%transferfunction * conjg(this%transferfunction)
    case default
       print *, 'ERROR: Unknown filter type: ', trim(filterclass)
       call pabort
    end select

    ! Set forward and backward transfer functions.
    ! STF of fwd and bwd earthquake may be added later
    this%transferfunction_fwd       = this%transferfunction
    this%transferfunction_bwd       = this%transferfunction

    if (firstslave) then
20     format('./Filters/filterresponse_', A, 2('_', F0.6))
       write(fnam,20) trim(filterclass), frequencies(1:2)

       open(10, file=trim(fnam), action='write')
       do ifreq = 1, nfreq
          write(10,*) this%f(ifreq), real(this%transferfunction(ifreq)), &
                                     imag(this%transferfunction(ifreq))
       end do
       close(10)
    end if   

    this%initialized = .true.
    this%stf_added = .false.

end subroutine create
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Devides the transferfunction of the filter with the complex spectra of 
!! the Source time function of the forward SEM simulation, to cancel its effect.
!! The filter uses only the forward STF, since fwd and bwd simulation are checked to
!! have the same STF in readfields.f90
!! We could get around this by having separate filters for the forward and the backward 
!! field. 
!! If we had nothing else to do.
subroutine add_stfs(this, stf_sem_fwd, sem_dt, amplitude_fwd,   &
                    istf, stf_source, stf_dt, stf_shift)
    use fft,                     only: rfft_type, taperandzeropad
    use simple_routines,         only: absreldiff
    use lanczos,                 only: lanczos_resample
    class(filter_type)              :: this
    real(kind=dp)   , intent(in)    :: stf_sem_fwd(:) ! STF of the AxiSEM simulation
    real(kind=dp)   , intent(in)    :: sem_dt         ! time step of STF
    real(kind=dp)   , intent(in)    :: amplitude_fwd
    integer,          intent(in)    :: istf           ! index of STF
    real(kind=dp)   , intent(in)    :: stf_source(:)  ! STF of the actual earthquake
    real(kind=dp)   , intent(in)    :: stf_dt         ! time step of STF
    real(kind=dp)   , intent(in)    :: stf_shift      ! time shift of STF

    ! The FFT routines need 2D arrays. Second dimension will be size 1
    real(kind=dp)   , allocatable   :: stf_sem(:,:), stf_src(:,:), t(:)
    real(kind=dp)   , allocatable   :: stf_src_td(:,:)
    real(kind=dp)   , allocatable   :: stf_sem_td(:,:)
    complex(kind=dp), allocatable   :: stf_src_fd(:,:), stf_sem_fd(:,:)
    complex(kind=dp), allocatable   :: lowpass(:), shift_fd(:)

    real(kind=dp)   , allocatable   :: stf_resampled(:)

    real(kind=dp)                   :: ratio_fwd, ratio_bwd

    type(rfft_type)                 :: fft_stf
    character(len=64)               :: fnam
    integer                         :: ifreq, i

    if (.not.this%initialized) then
       write(*,*) 'ERROR: Filter is not initialized yet'
       call pabort
    end if

    if (this%stf_added) then
       write(*,*) 'ERROR: STF has already been added to filter ', trim(this%name)
       call pabort
    end if

    call fft_stf%init(ntimes_in = size(stf_sem_fwd), &
                      ndim      = 1,                 &
                      ntraces   = 1,                 &           
                      dt        = sem_dt)


    allocate(stf_sem(size(stf_sem_fwd), 1))
    allocate(stf_sem_fd(this%nfreq, 1))
    allocate(stf_src(size(stf_sem_fwd), 1))
    allocate(stf_src_fd(this%nfreq, 1))
    allocate(shift_fd(this%nfreq))

    ! Allocate variables to store FTd and iFTd STFs
    allocate(stf_src_td(fft_stf%get_ntimes(), 1))
    allocate(stf_sem_td(fft_stf%get_ntimes(), 1))

    ! Allocate time variable
    allocate(t(fft_stf%get_ntimes()))

    stf_sem(:,1) = stf_sem_fwd 
    t = fft_stf%get_t()

    ! Resample earthquake STF to sampling rate of the SEM simulation
    stf_src(:,:) = 0.0d0

    if (absreldiff(stf_dt, sem_dt)>1d-8) then
      allocate(stf_resampled(int(size(stf_source) * stf_dt / sem_dt)))
      stf_resampled = lanczos_resample(stf_source, stf_dt, sem_dt, a=8)
      stf_src(1:size(stf_resampled,1),1) = stf_resampled
    else
      stf_src(1:size(stf_source,1),1) = stf_source
    end if

    ! Write out original Source STF
    if (firstslave) then
16     format('./Filters/stf_source_', A, 2('_', F0.3))
17     format(3(E16.8))
       write(fnam,16) trim(this%filterclass), this%frequencies(1:2)

       open(10, file=trim(fnam), action='write')
       do ifreq = 1, size(stf_source)
         write(10,17) stf_dt*(ifreq-1), stf_source(ifreq)
       end do
       close(10)
    end if

    ! FT STF of AxiSEM
    call fft_stf%rfft(taperandzeropad(stf_sem,              &
                                      fft_stf%get_ntimes(), &
                                      ntaper = 2),          &
                      stf_sem_fd)

    ! FT STF of Earthquake
    call fft_stf%rfft(taperandzeropad(stf_src,              &
                                      fft_stf%get_ntimes(), &
                                      ntaper = 2,           &
                                      end_only=.true.),     &
                      stf_src_fd)

    if (firstslave.or.testing) then
18     format('./Filters/stf_spectrum_', A, 2('_', F0.3))
19     format(5(E16.8))
       write(fnam,18) trim(this%filterclass), this%frequencies(1:2)

       open(10, file=trim(fnam), action='write')
       do ifreq = 1, this%nfreq
           write(10,19) this%f(ifreq), real(stf_sem_fd(ifreq,1)),    &
                                       imag(stf_sem_fd(ifreq,1)),    &
                                       real(stf_src_fd(ifreq,1)), &
                                       imag(stf_src_fd(ifreq,1))
       end do
       close(10)
    end if

    ! Divide filter by spectrum of forward STF. It is established at initialization
    ! that forward and backward STF are identical.
    this%transferfunction_fwd = this%transferfunction / stf_sem_fd(:,1)
    this%transferfunction_fwd = this%transferfunction_fwd * amplitude_fwd 

    this%transferfunction_bwd = this%transferfunction / stf_sem_fd(:,1)
    this%transferfunction_bwd = this%transferfunction_bwd * amplitude_fwd 

    ! Apply Earthquake STF, but only to forward transfer function
    this%transferfunction_fwd = this%transferfunction_fwd * stf_src_fd(:,1)

    ! Time shift STF
    shift_fd = exp( 2*pi * cmplx(0., 1.) * this%f(:) * stf_shift)
    this%transferfunction_fwd = this%transferfunction_fwd * shift_fd

    ! Apply high order butterworth filter to delete frequencies above mesh frequency
    allocate(lowpass(this%nfreq))
    ! print *, 'Butterworth at ', this%f(this%nfreq)*0.70d0, ' Hz'
    ! print *, 'Butterworth at ', 1./(this%f(this%nfreq)*0.70d0), ' sec'
    lowpass = butterworth_lowpass(this%f, 1.d0/(this%f(this%nfreq)*0.75d0), 12)
    lowpass = lowpass * conjg(lowpass)
  
    this%transferfunction_fwd       = this%transferfunction_fwd       * lowpass 
    this%transferfunction_bwd       = this%transferfunction_bwd       * lowpass 

    ! Replace NaNs with zero
    where(abs(this%transferfunction_fwd).ne.abs(this%transferfunction_fwd)) 
      this%transferfunction_fwd = 0
    end where

    where(abs(this%transferfunction_bwd).ne.abs(this%transferfunction_bwd)) 
      this%transferfunction_bwd = 0
    end where

    call fft_stf%irfft(stf_sem_fd, stf_sem_td)
    call fft_stf%irfft(stf_src_fd, stf_src_td)

    call fft_stf%freeme()

    if (firstslave.or.testing) then
20     format('./Filters/filterresponse_stf_', I0.2, '_', A, 2('_', F0.3))
       write(fnam,20) istf, trim(this%filterclass), this%frequencies(1:2)
       open(10, file=trim(fnam), action='write')
       do ifreq = 1, this%nfreq
          write(10,*) this%f(ifreq),  real(this%transferfunction(ifreq)), &
                                      imag(this%transferfunction(ifreq)), &
                                      real(this%transferfunction_fwd(ifreq)), &
                                      imag(this%transferfunction_fwd(ifreq)), &
                                      real(this%transferfunction_bwd(ifreq)), &
                                      imag(this%transferfunction_bwd(ifreq)), &
                                      real(lowpass(ifreq)), &
                                      imag(lowpass(ifreq))

       end do
       close(10)
       
21     format('./Filters/stf_spectrum_deriv_', I0.2, '_', A, 2('_', F0.3))
22     format(5(E16.8))
       write(fnam,21) istf, trim(this%filterclass), this%frequencies(1:2)

       open(10, file=trim(fnam), action='write')
       do ifreq = 1, this%nfreq
           write(10,22) this%f(ifreq), real(stf_sem_fd(ifreq,1)), &
                                        imag(stf_sem_fd(ifreq,1)), &
                                        real(stf_src_fd(ifreq,1)), &
                                        imag(stf_src_fd(ifreq,1))
       end do
       close(10)
       
23     format('./Filters/stf_in_', I0.2, '_', A, 2('_', F0.3))
24     format(3(E16.8))
       write(fnam,23) istf, trim(this%filterclass), this%frequencies(1:2)

       open(10, file=trim(fnam), action='write')
       do i = 1, size(stf_sem_fwd)
         write(10,24) t(i), stf_sem_fwd(i), stf_src(i,1)
       end do
       close(10)

       ! Contains the STFs after FFTs
25     format('./Filters/stf_out_', I0.2, '_', A, 2('_', F0.3))
26     format(3(E16.8))
       write(fnam,25) istf, trim(this%filterclass), this%frequencies(1:2)

       open(10, file=trim(fnam), action='write')
       do i = 1, size(stf_sem_fwd)
         write(10,26) t(i), stf_sem_td(i,1), stf_src_td(i,1)
       end do
       close(10)
       
    end if   

    ! Check whether transfer functions at highest frequency are close to zero
    ! or at least less than 1% of maximum value.
    ratio_fwd = abs(this%transferfunction_fwd(this%nfreq)) / &
                maxval(abs(this%transferfunction_fwd(:)))
    ratio_bwd = abs(this%transferfunction_bwd(this%nfreq)) / &
                maxval(abs(this%transferfunction_bwd(:)))

    if ((ratio_fwd > 1e-2) .or. (ratio_bwd > 1e-2)) then
      if (firstslave.or.testing) then
         print *, 'ERROR: Filter ', trim(this%name), ' is not vanishing fast enough for '
         print *, 'high frequencies.'
         print *, 'Numerical noise from frequencies above mesh limit will be propagated'
         print *, 'into the kernels. Check the files'
         print *, 'filterresponse*'
         print *, 'filterresponse_stf_*'
         print *, 'Relative value of forward TF at half mesh period : ', ratio_fwd
         print *, 'Relative value of backward TF at half mesh period :', ratio_bwd
         print *, 'Cutoff period: ', 1./this%f(this%nfreq)
       end if
       stop
    end if

    this%stf_added = .true.

end subroutine add_stfs
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Delete this filter and free the memory
subroutine deleteme(this)

    class(filter_type)              :: this

    this%filterclass = ''
    deallocate(this%transferfunction)
    deallocate(this%transferfunction_fwd)
    deallocate(this%transferfunction_bwd)
    deallocate(this%f)
    this%initialized = .false.

end subroutine deleteme
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Apply this filter to one trace (in the frequency domain)
function apply_1d(this, freq_series, kind)

   class(filter_type)               :: this
   complex(kind=dp), intent(in)     :: freq_series(:)
   character(len=*), intent(in)     :: kind
   complex(kind=dp)                 :: apply_1d(size(freq_series))

   if (.not.this%initialized) then
      write(*,*) 'ERROR: Filter is not initialized yet'
      call pabort(do_traceback=.false.)
   end if

   if (size(freq_series, 1).ne.this%nfreq) then
      write(*,*) 'ERROR: Filter length: ', this%nfreq, ', data length: ', &
                    size(freq_series, 1)
      call pabort(do_traceback=.false.)
   end if

   select case(kind)
   case('fwd')
     apply_1d = freq_series * this%transferfunction_fwd
   case('bwd')
     apply_1d = freq_series * this%transferfunction_bwd
   case default
     apply_1d = freq_series * this%transferfunction
   end select
end function apply_1d 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Apply this filter to multiple traces (in the frequency domain)
function apply_2d(this, freq_series, kind)

   class(filter_type)               :: this
   complex(kind=dp), intent(in)     :: freq_series(:,:)
   character(len=*), intent(in)     :: kind
   complex(kind=dp)                 :: apply_2d(size(freq_series,1), size(freq_series,2))
   integer                          :: itrace

   if (.not.this%initialized) then
      write(*,*) 'ERROR: Filter is not initialized yet'
      call pabort(do_traceback=.false.)
   end if

   if (size(freq_series, 1).ne.this%nfreq) then
      write(*,*) 'ERROR: Filter length: ', this%nfreq, ', data length: ', &
                    size(freq_series, 1)
      call pabort(do_traceback=.false.)
   end if
   
   select case(trim(kind))
   case('fwd')
     do itrace = 1, (size(freq_series, 2))
        apply_2d(:,itrace) = freq_series(:,itrace) * this%transferfunction_fwd(:)
     end do
   case('bwd')
     do itrace = 1, (size(freq_series, 2))
        apply_2d(:,itrace) = freq_series(:,itrace) * this%transferfunction_bwd(:)
     end do
   case default
     do itrace = 1, (size(freq_series, 2))
        apply_2d(:,itrace) = freq_series(:,itrace) * this%transferfunction(:)
     end do
   end select

end function apply_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Apply this filter to multiple dimensions and traces (in the frequency domain)
function apply_3d(this, freq_series, kind)

   use simple_routines, only       :  mult3d_1d
   
   class(filter_type)              :: this
   complex(kind=dp), intent(in)    :: freq_series(:,:,:)
   character(len=*), intent(in)    :: kind
   complex(kind=dp)                :: apply_3d(size(freq_series,1), size(freq_series,2), &
                                               size(freq_series,3))

   if (.not.this%initialized) then
      write(*,*) 'ERROR: Filter is not initialized yet'
      call pabort(do_traceback=.false.)
   end if

   if (size(freq_series, 1).ne.this%nfreq) then
      write(*,*) 'ERROR: Filter length: ', this%nfreq, ', data length: ', &
                    size(freq_series, 1)
      call pabort
   end if

   select case(kind)
   case('fwd')
     apply_3d = mult3d_1d(freq_series, this%transferfunction_fwd)
   case('bwd')
     apply_3d = mult3d_1d(freq_series, this%transferfunction_bwd)
   case default
     apply_3d = mult3d_1d(freq_series, this%transferfunction)
   end select

end function apply_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Returns the transfer function of this filter
function get_transferfunction(this)
    class(filter_type)           :: this
    real(kind=dp)                :: get_transferfunction(this%nfreq, 3)

    if (.not.this%initialized) then
       write(*,*) 'ERROR: Filter is not initialized yet'
       call pabort
    end if
    get_transferfunction(:,1) = this%f(:)
    get_transferfunction(:,2) = real(this%transferfunction(:))
    get_transferfunction(:,3) = imag(this%transferfunction(:))

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function isinitialized(this)
   class(filter_type)            :: this
   logical                       :: isinitialized

   isinitialized = this%initialized
end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Init a timeshift object
subroutine init_ts(this, freq, dtshift)
   class(timeshift_type)            :: this
   real(kind=dp),    intent(in)     :: freq(:)     !< 1D array with the frequencies of 
                                                   !! the entries in field
                                                   !! Dimension: nfreq
   real(kind=dp),    intent(in)     :: dtshift     !< Time shift to apply (in seconds)

   if (verbose>1) &
      write(lu_out, *) 'Initializing timeshift object with ', dtshift, ' sec'
   allocate(this%shift_fd(size(freq,1)) )
   this%shift_fd = exp( 2*pi * cmplx(0., 1.) * freq(:) * dtshift)
   this%isinitialized = .true.

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Apply a timeshift of dtshift on the frequency domain traces in field
subroutine timeshift_md(this, field)
   class(timeshift_type)            :: this
   complex(kind=dp), intent(inout)  :: field(:,:,:) !< Frequency domain traces to apply 
                                                   !! time shift on. 
                                                   !! Dimension: nfreq x ndim x ntraces
   if (this%isinitialized) then
     field(:,:,:) = mult3d_1d(field, this%shift_fd)
   else
     write(*,*) 'Timeshift is not initialized yet!'
     call pabort()
   end if

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Apply a timeshift of dtshift on the frequency domain traces in field
subroutine timeshift_1d(this, field)
   class(timeshift_type)            :: this
   complex(kind=dp), intent(inout)  :: field(:,:)  !< Frequency domain traces to apply 
                                                   !! time shift on. 
                                                   !! Dimension: nfreq x ntraces
   if (this%isinitialized) then
     field(:,:) = mult2d_1d(field, this%shift_fd)
   else
     write(*,*) 'Timeshift is not initialized yet!'
     call pabort()
   end if

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Delete this timeshift object
subroutine freeme(this)
   class(timeshift_type)           :: this

   if (allocated(this%shift_fd)) deallocate(this%shift_fd)

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------!
!                 BLOCK WITH SEPARATE FILTER ROUTINES                         !
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------
function ident(f)                        !< identical filter, returns a vector of ones the 
                                         !! length of f. Exists just for testing reasons.

    real(kind=dp), intent(in)   :: f(:)  !< Frequency array
    complex(kind=dp)            :: ident(size(f))

    ident(:) = cmplx(1.0, 0.0, kind=dp)

end function ident
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function loggabor(f, p_c, sigma, tshift)             !< Log-Gabor filter as employed by K. Sigloch

    real(kind=dp), intent(in)   :: f(:)  !< Frequency array
    real(kind=dp), intent(in)   :: p_c   !< Center period of the filter
    real(kind=dp), intent(in)   :: sigma !< (fixed) ratio of sigma divided by p_c, where
                                         !! sigma is the standard deviation of the
                                         !! Gaussian that describes the log-Gabor filter
                                         !! in the *time* domain, and fc is the filter's
                                         !! center frequency.  Hence larger sigmaIfc means
                                         !! narrower bandwidth sigmaIfc = .5 means the
                                         !! one-sigma interval i equals one octave
    real(kind=dp), intent(in)   :: tshift !< Adds a time shift to the impulse response 
                                         !! Useful to avoid wrap around on tsin that start at
                                         !! or near zero. Choose tshift on the order of
                                         !! p_c, where p_c is the center period.
                                         !! Defined in seconds

    complex(kind=dp)            :: loggabor(size(f))
    real(kind=dp)               :: f_c
    complex(kind=dp)            :: phase_shift(size(f))

    f_c = 1 / p_c
    loggabor = 0
    
    ! Calculate phase shift operator
    phase_shift = exp(cmplx(0d0, -2d0*pi, kind=dp) * f * tshift)

    ! Calculate filter itself
    where (f>0.0d0) loggabor = exp( -(log(f/f_c))**2 / ( 2 * log(sigma)**2) )

    ! Combine filter and phase shift  
    loggabor = loggabor * phase_shift  
    ! G(f,k) = exp( -(ln(f/f_k))^2 / (2*ln(sigmaIfc)^2)  )

end function loggabor
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function butterworth_lowpass(f, p_c, n)
    !http://en.wikipedia.org/wiki/Butterworth_filter

    real(kind=dp), intent(in)   :: f(:)  !< Frequency array
    real(kind=dp), intent(in)   :: p_c   !< cutoff period of the filter
    integer, intent(in)         :: n     !< order

    complex(kind=dp)            :: butterworth_lowpass(size(f))
    real(kind=dp)               :: omega_c
    complex(kind=dp)            :: s(1:n), j
    integer                     :: k

    omega_c = 2 * pi / p_c
    j = cmplx(0, 1)

    do k=1, n
        s(k) = omega_c * exp(j * (2 * k + n - 1) * pi / (2 * n))
    enddo

    butterworth_lowpass = 1

    do k=1, n
        butterworth_lowpass = butterworth_lowpass / ((2 * pi * j * f - s(k)) / omega_c)
    enddo

end function butterworth_lowpass
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function butterworth_highpass(f, p_c, n)
    !http://en.wikipedia.org/wiki/Butterworth_filter
    !http://en.wikipedia.org/wiki/Prototype_filter#Lowpass_to_highpass

    real(kind=dp), intent(in)   :: f(:)  !< Frequency array
    real(kind=dp), intent(in)   :: p_c   !< cutoff period of the filter
    integer, intent(in)         :: n     !< order

    complex(kind=dp)            :: butterworth_highpass(size(f))
    real(kind=dp)               :: omega_c
    complex(kind=dp)            :: s(1:n), j
    integer                     :: k

    omega_c = 2 * pi / p_c
    j = cmplx(0, 1)

    do k=1, n
        s(k) = omega_c * exp(j * (2 * k + n - 1) * pi / (2 * n))
    enddo

    butterworth_highpass = 1

    do k=1, n
        butterworth_highpass(2:) &
            = butterworth_highpass(2:) / (omega_c / (2 * pi * j * f(2:)) - s(k) / omega_c)
    enddo
    butterworth_highpass(1) = 0

end function butterworth_highpass
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
