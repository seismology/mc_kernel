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
module test_filter

   use global_parameters, only: sp, dp, pi
   use fft, only: rfft_type, taperandzeropad
   use filtering
   use ftnunit
   implicit none
   public
contains

!-----------------------------------------------------------------------------------------
subroutine test_filter_ident
   ! Test identical filter. Is later used to test STF deconvolution, so we should be 
   ! certain that it works
   type(filter_type)             :: ident
   type(rfft_type)               :: fft_data
   integer, parameter            :: ntimes = 1000, npoints = 1
   integer                       :: ntimes_ft, nomega
   character(len=32)             :: filtername, filterclass
   real(kind=dp)                 :: data_in(ntimes,1)
   real(kind=dp), allocatable    :: data_filtered_ref(:,:), data_filtered(:,:), tf(:,:), tf_ref(:,:)
   complex(kind=dp), allocatable :: data_fd(:,:)
   real(kind=dp)                 :: df, thalf
   real(kind=dp), parameter      :: dt = 0.1d0
   integer                       :: x

   filterclass = 'ident'
   filtername  = 'identical'

   call fft_data%init(ntimes, 1, npoints, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()
   df = fft_data%get_df()

   allocate(data_filtered(ntimes_ft,1))
   allocate(data_filtered_ref(ntimes_ft, 1))
   allocate(data_fd(nomega,1))
   allocate(tf(nomega, 3))
   allocate(tf_ref(nomega,2))

   ! Reference transfer function: Real part: 1, imaginary part: 0
   tf_ref(:,1) = 1.0d0
   tf_ref(:,2) = 0.0d0

   ! Data to filter: Fumiko-style STF:
   thalf = 20.d0
   data_in(1:ntimes/4,1) = 0
   data_in(ntimes/4+1:ntimes,1) = [(2.d0 * exp(-(x/thalf)**2)*x/thalf**2, x=1, ntimes-ntimes/4 )]
   data_filtered_ref(:, 1) = 0d0
   data_filtered_ref(1:ntimes, 1) = data_in(:,1)

   ! Test identical filter
   call ident%create(filtername, df, nomega, filterclass, [0.0d0, 0.0d0, 0.d0, 0.d0])

   call assert_true(ident%isinitialized(), 'filter is initialized after creation')

   tf = ident%get_transferfunction() 
  
   call assert_comparable_real(real(tf(2,1)-tf(1,1)), real(df), 1e-10, &
                               'df returned by filter equal to df')

   call assert_comparable_real1d(real(tf_ref(1:1024,1)), real(tf(1:1024,2)), 1e-5, &
                                 ' real part of transfer function is one')

   call assert_comparable_real(real(sum(tf(:,3))), 0.0, 1e-10, &
                               'complex part of transfer function is zero')

   call fft_data%rfft(taperandzeropad(data_in, ntimes_ft), data_fd)
   data_fd = ident%apply_2d(data_fd, kind='def')
   call fft_data%irfft(data_fd, data_filtered)

   call assert_comparable(data_filtered(:,1)+1.d0, &
                          data_filtered_ref(:,1)+1.d0, 1d-15, &
                          'Response of ident filter ')

   call ident%deleteme()

   call assert_true(.not.ident%isinitialized(), 'filter is not initialized after deletion')

end subroutine test_filter_ident
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_filter_gabor_response
   use fft, only: rfft_type, taperandzeropad
   ! Compares transfer and response functions with solutions from Karin's 
   ! Matlab code. Tries log-Gabor filters with 5s and 2.5s central period
   ! and a width of one octave

   type(filter_type)             :: gabor
   type(rfft_type)               :: fft_data
   integer, parameter            :: ntimes = 1000, npoints = 1, ntests = 12
   integer                       :: ntimes_ft, nomega, itest, i
   real(kind=dp)                 :: data_in(ntimes,1)
   real(kind=dp), allocatable    :: data_filtered_ref(:,:), data_filtered(:,:), tf(:,:)
   complex(kind=dp), allocatable :: data_fd(:,:)
   real(kind=dp)                 :: df
   real(kind=dp), parameter      :: dt = 0.1d0
   character(len=32)             :: filtername, filterclass, filename_ref
   character(len=64)             :: assert_text
   ! Filter settings
   real(kind=dp)                 :: pmax(ntests), tshift(ntests), sigmaIfc(ntests)

   ! Center periods, time shifts and sigmas for all reference files
   pmax =     [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
   tshift =   [0.0, 2.5, 0.0, 2.5, 0.0, 2.5, 0.0, 5.0, 0.0, 5.0, 0.0, 5.0]
   sigmaIfc = [0.3, 0.3, 0.5, 0.5, 0.7, 0.7, 0.3, 0.3, 0.5, 0.5, 0.7, 0.7]

   call fft_data%init(ntimes, 1, npoints, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()
   df = fft_data%get_df()

   call assert_equal(ntimes_ft, 2048, 'ntimes==expected value')

   open(100, file='./gaborinput.txt')
   read(100,*) data_in
   close(100)

   allocate(data_filtered(ntimes_ft,1))
   allocate(data_filtered_ref(ntimes_ft, 1))
   allocate(data_fd(nomega,1))
   allocate(tf(nomega, 3))

   do itest = 1, ntests
     write(filename_ref,'("./gaborresponse_", I2.2, ".txt")') itest
     open(100, file=trim(filename_ref))
     read(100,*) data_filtered_ref
     close(100)
   
     filterclass = 'Gabor'
     filtername  = 'Gabor test'
     call gabor%create(filtername, df, nomega, filterclass, &
                       [pmax(itest), sigmaIfc(itest), tshift(itest), 0d0])

     call assert_true(gabor%isinitialized(), 'filter is initialized after creation')

     tf = gabor%get_transferfunction() 
     
     call assert_comparable_real(real(tf(2,1)-tf(1,1)), real(df), 1e-10, &
                                 'df returned by filter equal to df')


     call fft_data%rfft(taperandzeropad(data_in, ntimes_ft), data_fd)
     data_fd = gabor%apply_2d(data_fd, kind='def')
     call fft_data%irfft(data_fd, data_filtered)

     write(filename_ref,'("./gaboroutput_", I2.2, ".txt")') itest
     open(100, file=trim(filename_ref), action='write')
     do i=1, size(data_filtered_ref)
       write(100,*) data_filtered_ref(i,1)
     end do
     close(100)

     write(assert_text, '("Gabor filter with (", 3(F8.3),")")') pmax(itest), &
             sigmaIfc(itest), tshift(itest)
     call assert_comparable(1e0+real(data_filtered(:,1)), &
                            1e0+real(data_filtered_ref(:,1)), 1e-06, &
                            assert_text)

     call gabor%deleteme()
   end do

   call fft_data%freeme()
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_filter_butterworth_lp_response
   use fft, only: rfft_type, taperandzeropad
   
   ! reference data was produced with the same code, because a comparison to obspy.lowpass
   ! revealed some Gibbs effect for the input dirac but otherwise good agreement

   type(filter_type)             :: butter
   type(rfft_type)               :: fft_data
   integer, parameter            :: ntimes = 1000, npoints = 1
   integer                       :: ntimes_ft, nomega
   real(kind=dp)                 :: data_in(ntimes,1)
   real(kind=dp), allocatable    :: data_filtered_ref(:,:), data_filtered(:,:)
   complex(kind=dp), allocatable :: data_fd(:,:)
   real(kind=dp)                 :: df
   real(kind=dp), parameter      :: dt = 0.1d0
   character(len=32)             :: filtername, filterclass

   call fft_data%init(ntimes, 1, npoints, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   call assert_equal(ntimes_ft, 2048, 'ntimes==expected value')

   open(100, file='./gaborinput.txt')
   read(100,*) data_in
   close(100)

   allocate(data_filtered(ntimes_ft,1))
   allocate(data_filtered_ref(ntimes_ft,1))
   allocate(data_fd(nomega,1))

   open(100, file='./butterresponse.txt')
   read(100,*) data_filtered_ref
   close(100)
   
   df = fft_data%get_df()

   ! Test filter with 5s period
   filterclass = 'Butterw_LP_O2'
   filtername  = 'Butterw test'
   call butter%create(filtername, df, nomega, filterclass, [5.0d0, 0.d0, 0.d0, 0.d0])

   call assert_true(butter%isinitialized(), 'filter is initialized after creation')

   call fft_data%rfft(taperandzeropad(data_in, ntimes_ft), data_fd)
   data_fd = butter%apply_2d(data_fd, kind='def')
   call fft_data%irfft(data_fd, data_filtered)

   call assert_comparable_real1d(real(data_filtered(:,1) + 0.1d0), &
                                 real(data_filtered_ref(:,1) + 0.1d0), 1e-7, &
                                 'Butterworth LP filter with 5s cutoff period')

   call butter%deleteme()

   call assert_true(.not.butter%isinitialized(), &
                    'filter is not initialized after deletion')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_filter_butterworth_hp_response
   use fft, only: rfft_type, taperandzeropad
   
   ! reference data was produced with the same code, because a comparison to obspy.highpass
   ! revealed some Gibbs effect for the input dirac but otherwise good agreement

   type(filter_type)             :: butter
   type(rfft_type)               :: fft_data
   integer, parameter            :: ntimes = 1000, npoints = 1
   integer                       :: ntimes_ft, nomega
   real(kind=dp)                 :: data_in(ntimes,1)
   real(kind=dp), allocatable    :: data_filtered_ref(:,:), data_filtered(:,:)
   complex(kind=dp), allocatable :: data_fd(:,:)
   real(kind=dp)                 :: df
   real(kind=dp), parameter      :: dt = 0.1d0
   character(len=32)             :: filtername, filterclass

   call fft_data%init(ntimes, 1, npoints, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   call assert_equal(ntimes_ft, 2048, 'ntimes==expected value')

   open(100, file='./gaborinput.txt')
   read(100,*) data_in
   close(100)

   allocate(data_filtered(ntimes_ft,1))
   allocate(data_filtered_ref(ntimes_ft,1))
   allocate(data_fd(nomega,1))

   open(100, file='./butterresponse_hp.txt')
   read(100,*) data_filtered_ref
   close(100)
   
   df = fft_data%get_df()

   ! Test filter with 5s period
   filterclass = 'Butterw_HP_O2'
   filtername  = 'Butterw test'
   call butter%create(filtername, df, nomega, filterclass, [5.0d0, 0.d0, 0.d0, 0.d0])

   call assert_true(butter%isinitialized(), 'filter is initialized after creation')

   call fft_data%rfft(taperandzeropad(data_in, ntimes_ft, ntaper=0), data_fd)
   data_fd = butter%apply_2d(data_fd, kind='def')
   call fft_data%irfft(data_fd, data_filtered)

   call assert_comparable_real1d(real(data_filtered(:,1) + 0.1d0), &
                                 real(data_filtered_ref(:,1) + 0.1d0), 1e-7, &
                                 'Butterworth HP filter with 5s cutoff period')

   call butter%deleteme()

   call assert_true(.not.butter%isinitialized(), &
                    'filter is not initialized after deletion')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_filter_butterworth_bp_response
   use fft, only: rfft_type, taperandzeropad
   
   ! reference data was produced with the same code, though using 4th order here reduces 
   ! the Gibbs effect from the previous two tests by a lot

   type(filter_type)             :: butter
   type(rfft_type)               :: fft_data
   integer, parameter            :: ntimes = 1000, npoints = 1
   integer                       :: ntimes_ft, nomega
   real(kind=dp)                 :: data_in(ntimes,1)
   real(kind=dp), allocatable    :: data_filtered_ref(:,:), data_filtered(:,:)
   complex(kind=dp), allocatable :: data_fd(:,:)
   real(kind=dp)                 :: df
   real(kind=dp), parameter      :: dt = 0.1d0
   character(len=32)             :: filtername, filterclass

   call fft_data%init(ntimes, 1, npoints, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   call assert_equal(ntimes_ft, 2048, 'ntimes==expected value')

   open(100, file='./gaborinput.txt')
   read(100,*) data_in
   close(100)

   allocate(data_filtered(ntimes_ft,1))
   allocate(data_filtered_ref(ntimes_ft,1))
   allocate(data_fd(nomega,1))

   open(100, file='./butterresponse_bp.txt')
   read(100,*) data_filtered_ref
   close(100)
   
   df = fft_data%get_df()

   ! Test filter with 10s - 5s passband
   filterclass = 'Butterw_BP_O4'
   filtername  = 'Butterw test'
   call butter%create(filtername, df, nomega, filterclass, [10.0d0, 5.d0, 0.d0, 0.d0])

   call assert_true(butter%isinitialized(), 'filter is initialized after creation')

   call fft_data%rfft(taperandzeropad(data_in, ntimes_ft, ntaper=0), data_fd)
   data_fd = butter%apply_2d(data_fd, kind='def')
   call fft_data%irfft(data_fd, data_filtered)

   call assert_comparable_real1d(real(data_filtered(:,1) + 0.1d0), &
                                 real(data_filtered_ref(:,1) + 0.1d0), 1e-7, &
                                 'Butterworth BP filter with 10s - 5s period passband')

   call butter%deleteme()

   call assert_true(.not.butter%isinitialized(), &
                    'filter is not initialized after deletion')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_filter_timeshift()
   use fft, only: rfft_type, taperandzeropad
   use filtering, only: timeshift_type

   type(timeshift_type)             :: timeshift
   type(rfft_type)                  :: fft_data
   integer, parameter               :: ntimes=32, npoints=1
   integer                          :: ntimes_ft, nomega
   real(kind=dp)                    :: dt
   real(kind=dp)                    :: data_in(ntimes, npoints)
   real(kind=dp)                    :: data_out(ntimes*2, npoints), data_ref(ntimes*2)
   complex(kind=dp), allocatable    :: data_fd(:,:)

   dt = 0.1

   call fft_data%init(ntimes, 1, npoints, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()
   allocate(data_fd(nomega, npoints))

   call assert_equal(ntimes_ft, ntimes*2, 'ntimes==expected value')

   ! Time shift forward 1s (10 samples)
   data_in       = 0.0d0
   data_in(15,1) = 1.0

   data_ref      = 0.0d0
   data_ref(5)  = 1.0

   call fft_data%rfft(taperandzeropad(data_in, ntimes_ft), data_fd)

   call timeshift%init_ts(fft_data%get_f(), 1.0d0)
   call timeshift%apply(data_fd)
   call timeshift%freeme()
   call fft_data%irfft(data_fd, data_out)

   call assert_comparable_real1d(1+real(data_out(:,1), kind=sp), &
                                 1+real(data_ref, kind=sp),      &
                                 1e-6,                         &
                                 'Time shift forward 1 sec')


   ! Time shift backward 1s (10 samples)

   data_in       = 0.0d0
   data_in(15,1) = 1.0

   data_ref      = 0.0d0
   data_ref(25)  = 1.0

   call fft_data%rfft(taperandzeropad(data_in, ntimes_ft), data_fd)
   
   call timeshift%init_ts(fft_data%get_f(), -1.0d0)
   call timeshift%apply(data_fd)
   call timeshift%freeme()

   call fft_data%irfft(data_fd, data_out)

   call assert_comparable_real1d(1+real(data_out(:,1), kind=sp), &
                                 1+real(data_ref, kind=sp),      &
                                 1e-6,                         &
                                 'Time shift backward 1 sec')


   ! Time shift forward 6.4s (64 samples, should remain the same)

   data_in       = 0.0d0
   data_in(15,1) = 1.0

   data_ref      = 0.0d0 
   data_ref(15)  = 1.0

   call fft_data%rfft(taperandzeropad(data_in, ntimes_ft), data_fd)
   
   call timeshift%init_ts(fft_data%get_f(), 6.4d0)
   call timeshift%apply(data_fd)
   call timeshift%freeme()

   call fft_data%irfft(data_fd, data_out)

   call assert_comparable_real1d(1+real(data_out(:,1), kind=sp), &
                                 1+real(data_ref, kind=sp),      &
                                 1e-6,                         &
                                 'Time shift 2 pi')


   call fft_data%freeme()

end subroutine test_filter_timeshift
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
