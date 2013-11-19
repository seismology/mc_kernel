!=========================================================================================
module test_filter

   use global_parameters, only: sp, pi
   use filtering
   use ftnunit
   implicit none
   public
contains

!------------------------------------------------------------------------------
subroutine test_filter_gabor_response
   use fft, only: rfft_type, taperandzeropad
   ! Compares transfer and response functions with solutions from Karin's 
   ! Matlab code. Tries log-Gabor filters with 5s and 2.5s central period
   ! and a width of one octave

   type(filter_type)             :: gabor
   type(rfft_type)               :: fft_data
   integer, parameter            :: ntimes = 1000, npoints = 1
   integer                       :: ntimes_ft, nomega, i 
   real(kind=sp)                 :: data_in(ntimes,1)
   real(kind=dp), allocatable    :: data_filtered_ref(:,:), data_filtered(:,:), tf(:,:), tf_ref(:,:)
   complex(kind=dp), allocatable :: data_fd(:,:)
   real(kind=sp)                 :: windowlength
   real(kind=dp)                 :: df
   real(kind=dp), parameter      :: dt = 0.1d0
   character(len=32)             :: filtername, filterclass

   call fft_data%init(ntimes, npoints, dt)
   ntimes_ft = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()

   call assert_equal(ntimes_ft, 2048, 'ntimes==expected value')

   open(100, file='unit_tests/gaborinput.txt')
   read(100,*) data_in
   close(100)

   allocate(data_filtered(ntimes_ft,1))
   allocate(data_filtered_ref(ntimes_ft, 2))
   allocate(data_fd(nomega,1))
   allocate(tf(nomega, 3))
   allocate(tf_ref(nomega,2))

   open(100, file='unit_tests/gaborresponse.txt')
   read(100,*) data_filtered_ref
   close(100)
   
   open(100, file='unit_tests/gabortransferfunction.txt')
   read(100,*) tf_ref(1:1024,1)
   read(100,*) tf_ref(1:1024,2)
   close(100)

   df = fft_data%get_df()

   ! Test filter with 5s period
   filterclass = 'Gabor'
   filtername  = 'Gabor test'
   call gabor%create(filtername, df, nomega, filtername, [5.0d0, 0.5d0, 0.d0, 0.d0])

   call assert_true(gabor%isinitialized(), 'filter is initialized after creation')

   tf = gabor%get_transferfunction() 
   
  
   call assert_comparable_real(real(tf(2,1)-tf(1,1)), real(df), 1e-10, &
                               'df returned by filter equal to df')

   call assert_comparable_real1d(real(tf_ref(1:1024,1)), real(tf(1:1024,2)), 1e-5, &
                                 ' real part of transfer function is correct')

   call assert_comparable_real(real(sum(tf(:,3))), 0.0, 1e-10, &
                               'complex part of transfer function is zero')

   call fft_data%rfft(taperandzeropad(data_in, ntimes_ft), data_fd)
   data_fd = gabor%apply_2d(data_fd)
   call fft_data%irfft(data_fd, data_filtered)

   call assert_comparable_real1d(real(data_filtered(:,1) + 0.1d0), &
                                 real(data_filtered_ref(:,1) + 0.1d0), 1e-7, &
                                 'Gabor filter with 5s center period')

   call gabor%deleteme()

   call assert_true(.not.gabor%isinitialized(), 'filter is not initialized after deletion')

   ! Test filter with 2.5s period
   filtername = 'Gabor'
   filtername  = 'Gabor test'
   call gabor%create(filtername, df, nomega, filterclass, [2.5d0, 0.5d0, 0.d0, 0.d0])

   tf = gabor%get_transferfunction() 
   
   call assert_comparable_real1d(real(tf_ref(1:1024,2)), real(tf(1:1024,2)), 1e-5, &
                                 ' real part of transfer function is correct')

   call assert_comparable_real(real(sum(tf(:,3))), 0.0, 1e-10, &
                               'complex part of transfer function is zero')

   call fft_data%rfft(taperandzeropad(data_in, ntimes_ft), data_fd)
   data_fd = gabor%apply_2d(data_fd)
   call fft_data%irfft(data_fd, data_filtered)

   call assert_comparable_real1d(real(data_filtered(:,1) + 0.1d0), &
                                 real(data_filtered_ref(:,2) + 0.1d0), 1e-7, &
                                 'Gabor filter with 5s center period')

   call gabor%deleteme()
   call fft_data%freeme()
end subroutine
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine test_filter_timeshift()
   use fft, only: rfft_type, taperandzeropad

   type(rfft_type)                  :: fft_data
   integer, parameter               :: ntimes=32, npoints=1
   integer                          :: ntimes_ft, nomega, i
   real(kind=dp)                    :: dt
   real(kind=sp)                    :: data_in(ntimes, npoints)
   real(kind=dp)                    :: data_out(ntimes*2, npoints), data_ref(ntimes*2)
   complex(kind=dp), allocatable    :: data_fd(:,:)

   dt = 0.1

   call fft_data%init(ntimes, npoints, dt)
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
   call timeshift(data_fd, fft_data%get_f(), 1.0)
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
   call timeshift(data_fd, fft_data%get_f(), -1.0)
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
   call timeshift(data_fd, fft_data%get_f(), 6.4)
   call fft_data%irfft(data_fd, data_out)

   call assert_comparable_real1d(1+real(data_out(:,1), kind=sp), &
                                 1+real(data_ref, kind=sp),      &
                                 1e-6,                         &
                                 'Time shift 2 pi')


   call fft_data%freeme()

end subroutine test_filter_timeshift
!------------------------------------------------------------------------------


end module

