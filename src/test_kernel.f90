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
module test_kernel

  use global_parameters
  use kernel
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------

subroutine test_kernel_cut_timewindow
   integer, parameter         :: ntimes = 100
   real(kind=dp)              :: timeseries(ntimes), dataseries(ntimes)
   real(kind=dp)              :: dt = 0.1d0
   real(kind=dp), allocatable :: cut_dataseries(:), cut_dataseries_ref(:)
   integer                    :: i

   ! Time window with normal length
   do i=1, ntimes
      timeseries(i) = i * dt
      dataseries(i) = i 
      ! print *, timeseries(i), dataseries(i)
   end do

   
   ! Test 1: Time window bounds are between samples
   allocate(cut_dataseries_ref(6))
   cut_dataseries_ref = [15, 16, 17, 18, 19, 20]

   call cut_timewindow(timeseries, dataseries, [1.45d0, 2.01d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 1)')

   deallocate(cut_dataseries_ref)
  
   ! Test 2: Upper time window bound is exactly on sample
   !         Correct behaviour would be that this sample is included in the resulting time window.
   !         This test may fail due to floating point comparison, however, 2 should have an exact 
   !         floating point representation 
   allocate(cut_dataseries_ref(6))
   cut_dataseries_ref = [15, 16, 17, 18, 19, 20]

   call cut_timewindow(timeseries, dataseries, [1.45d0, 2.0d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 2)')

   deallocate(cut_dataseries_ref)

   ! Test 3: Upper time window bound is narrowly (but within real(8) resolution) below sample
   allocate(cut_dataseries_ref(5))
   cut_dataseries_ref = [15, 16, 17, 18, 19]

   call cut_timewindow(timeseries, dataseries, [1.45d0, 1.9999999d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 3)')


   deallocate(cut_dataseries_ref)

   ! Test 4: Time window with length 1
   allocate(cut_dataseries_ref(1))
   cut_dataseries_ref = [80]

   call cut_timewindow(timeseries, dataseries, [7.95d0, 8.05d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 4)' )



   ! Test 5: Time window with negative starting time
   do i=1, ntimes
      timeseries(i) = i * dt - 5
      dataseries(i) = i - 50
   end do

   deallocate(cut_dataseries_ref)
   allocate(cut_dataseries_ref(6))
   
   cut_dataseries_ref = [-3, -2, -1, 0, 1, 2]

   call cut_timewindow(timeseries, dataseries, [-0.39d0, 0.29d0], cut_dataseries)

   call assert_comparable_real1d(real(cut_dataseries, kind=sp), real(cut_dataseries_ref, kind=sp), &
                                 1e-9, 'cut_dataseries == cut_dataseries_ref (Test 5)')


end subroutine test_kernel_cut_timewindow
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine test_kernel_init

   type(kernelspec_type)       :: kernel
   type(filter_type), target   :: gabor

   character(len=32)           :: filtername, filterclass
   real(kind=dp)               :: veloseis(256)

   call assert_false(kernel%isinitialized(), 'Kernel not initialized')
   ! Just some arbitrary filter to create
   filtername  = 'Gabor'
   filterclass = 'Gabor'
   call gabor%create(filtername, 1.0d0, 257, filterclass, [5.0d0, 0.5d0, 0.d0, 0.d0])

   veloseis = 0.0

   call kernel%init(name            = 'Testkernel      ',&
                    time_window     = [1.0d0, 2.0d0],    &
                    filter          = gabor,             &
                    misfit_type     = 'CC  ',            &  
                    model_parameter = 'vp  ',            &
                    seis            = veloseis,          &
                    dt              = 0.1d0,             &
                    timeshift_fwd   = 1.0d0,             &
                    deconv_stf      = .false.)

   call assert_true(kernel%isinitialized(), 'Kernel initialized')
   call kernel%freeme()
   
   call assert_false(kernel%isinitialized(), 'Kernel not initialized')

end subroutine test_kernel_init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine test_integrate_parseval_real

   type(kernelspec_type)       :: kernel
   type(filter_type), target   :: gabor

   character(len=32)           :: filtername, filterclass
   real(kind=dp)               :: veloseis(256), integral_parseval, integral_trapezoidal
   real(kind=dp), allocatable  :: func_to_int1(:), func_to_int2(:)

   ! Just create some kernel with a length of 256 samples
   ! Same as in test_kernel_init
   filtername  = 'Gabor'
   filterclass = 'Gabor'
   call gabor%create(filtername, 1.0d0, 257, filterclass, [5.0d0, 0.5d0, 0.d0, 0.d0])

   veloseis = 0.0

   call kernel%init(name            = 'Testkernel      ',&
                    time_window     = [1.0d0, 2.0d0],    &
                    filter          = gabor,             &
                    misfit_type     = 'CC  ',            &  
                    model_parameter = 'vp  ',            &
                    seis            = veloseis,          &
                    dt              = 0.1d0,             &
                    timeshift_fwd   = 1.0d0,             &
                    deconv_stf      = .false.)

   allocate(func_to_int1(kernel%ntimes))
   allocate(func_to_int2(kernel%ntimes))

   func_to_int1(:) = 1 - 4*(kernel%t_cut - 1.5)**2
   func_to_int2(:) = 1

   integral_parseval    = kernel%integrate_parseval(func_to_int1, func_to_int2)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-2,  &
                          'Integral of 1 * (1-4(x-1.5)**2) from 1 to 2')


   ! Now create a function to integrate (sum of two sines)

   func_to_int1(:) = sin(kernel%t_cut * 2 * pi)
   func_to_int2(:) = sin(kernel%t_cut * 2 * pi)

   integral_parseval    = kernel%integrate_parseval(func_to_int1, func_to_int2)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-7,  &
                          'Integral of sin(2pix)**2 from 1 to 2')


   func_to_int1(:) = sin(kernel%t_cut * 2 * pi)
   func_to_int2(:) = cos(kernel%t_cut * 3 * pi)

   integral_parseval = kernel%integrate_parseval(func_to_int1, func_to_int2)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-7,  &
                          'Integral of sin(2pix)*cos(3pix) from 1 to 2')
   !call assert_comparable(integral_parseval, 0.26568757573375174d0, 1d-7,  &


   func_to_int1(:) = sin(kernel%t_cut * 2.2 * pi)
   func_to_int2(:) = cos(kernel%t_cut * 1.7 * pi)

   integral_parseval = kernel%integrate_parseval(func_to_int1, func_to_int2)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   !call assert_comparable(integral_parseval, 0.32528761700943976d0, 1d-7,  &
   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-2,  &
                          'Integral of sin(2.2pix)*cos(1.7pix) from 1 to 2')

   call kernel%freeme()

end subroutine test_integrate_parseval_real
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine test_integrate_parseval_complex
   use fft,                    only: rfft_type, taperandzeropad

   type(kernelspec_type)          :: kernel
   type(filter_type), target      :: gabor
   type(rfft_type)                :: fft_data
  

   character(len=32)              :: filtername, filterclass
   real(kind=dp)                  :: veloseis(256), integral_parseval, integral_trapezoidal
   real(kind=dp), allocatable     :: func_to_int1(:), func_to_int2(:)
   complex(kind=dp), allocatable  :: func_to_int1_fd(:,:), func_to_int2_fd(:,:)
   real(kind=dp), parameter       :: dt=0.1d0
   integer                        :: nomega, ntimes_ft, ntimes

   ! Just create some kernel with a length of 256 samples
   ! Same as in test_kernel_init
   filtername  = 'Gabor'
   filterclass = 'Gabor'
   call gabor%create(filtername, 1.0d0, 257, filterclass, [5.0d0, 0.5d0, 0.d0, 0.d0])

   veloseis = 0.0

   call kernel%init(name            = 'Testkernel      ',&
                    time_window     = [1.0d0, 2.0d0],    &
                    filter          = gabor,             &
                    misfit_type     = 'CC  ',            &  
                    model_parameter = 'vp  ',            &
                    seis            = veloseis,          &
                    dt              = dt,                &
                    timeshift_fwd   = 1.0d0,             &
                    deconv_stf      = .false.)

   ! Create a FFT object
   call fft_data%init(kernel%ntimes, 1, 1, dt)
   nomega    = fft_data%get_nomega()
   ntimes_ft = fft_data%get_ntimes()
   ntimes    = kernel%ntimes

   ! Now create a function to integrate (sum of two sines)

   allocate(func_to_int1(ntimes))
   allocate(func_to_int2(ntimes))
   allocate(func_to_int1_fd(nomega,1))
   allocate(func_to_int2_fd(nomega,1))

   func_to_int1(:) = 1 - 4*(kernel%t_cut - 1.5)**2
   func_to_int2(:) = 1

   call fft_data%rfft(taperandzeropad(array  = reshape(func_to_int2, [ntimes, 1]), &
                                      ntimes = ntimes_ft,   &
                                      ntaper = 0),          &
                      func_to_int2_fd)

   integral_parseval    = kernel%integrate_parseval(func_to_int1, func_to_int2_fd)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-2,  &
                          'Integral of 1 * (1-4(x-1.5)**2) from 1 to 2')

   func_to_int1(:) = sin(kernel%t_cut * 2 * pi)
   func_to_int2(:) = sin(kernel%t_cut * 2 * pi)

   call fft_data%rfft(taperandzeropad(array  = reshape(func_to_int2, [ntimes, 1]), &
                                      ntimes = ntimes_ft,   &
                                      ntaper = 0),          &
                      func_to_int2_fd)

   integral_parseval    = kernel%integrate_parseval(func_to_int1, func_to_int2_fd)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-7,  &
                          'Integral of sin(2pix)**2 from 1 to 2')


   func_to_int1(:) = sin(kernel%t_cut * 2 * pi)
   func_to_int2(:) = cos(kernel%t_cut * 3 * pi)

   call fft_data%rfft(taperandzeropad(array  = reshape(func_to_int2, [ntimes, 1]), &
                                      ntimes = ntimes_ft,   &
                                      ntaper = 0),          &
                      func_to_int2_fd)

   integral_parseval    = kernel%integrate_parseval(func_to_int1, func_to_int2_fd)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-7,  &
                          'Integral of sin(2pix)*cos(3pix) from 1 to 2')


   func_to_int1(:) = sin(kernel%t_cut * 2.2 * pi)
   func_to_int2(:) = cos(kernel%t_cut * 1.7 * pi)

   call fft_data%rfft(taperandzeropad(array  = reshape(func_to_int2, [ntimes, 1]), &
                                      ntimes = ntimes_ft,   &
                                      ntaper = 0),          &
                      func_to_int2_fd)

   integral_parseval    = kernel%integrate_parseval(func_to_int1, func_to_int2_fd)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   ! In fact, the trapezoidal integration has a problem here
   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-2,  &
                          'Integral of sin(2.2pix)*cos(1.7pix) from 1 to 2')


   func_to_int1(:) = sin(kernel%t_cut * 2 * pi)
   func_to_int2(:) = cos(kernel%t_cut * 1.5 * pi) + 1

   call fft_data%rfft(taperandzeropad(array  = reshape(func_to_int2, [ntimes, 1]), &
                                      ntimes = ntimes_ft,   &
                                      ntaper = 0),          &
                      func_to_int2_fd)

   integral_parseval    = kernel%integrate_parseval(func_to_int1, func_to_int2_fd)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   ! In fact, the trapezoidal integration has a problem here
   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-5,  &
                          'Integral of sin(2pix)*(cos(1.5pix)+1) from 1 to 2')

   func_to_int1(:) = sin(kernel%t_cut * 2 * pi) 
   func_to_int2(:) = cos(kernel%t_cut * 1.5 * pi) - 1

   call fft_data%rfft(taperandzeropad(array  = reshape(func_to_int2, [ntimes, 1]), &
                                      ntimes = ntimes_ft,   &
                                      ntaper = 0),          &
                      func_to_int2_fd)

   integral_parseval    = kernel%integrate_parseval(func_to_int1, func_to_int2_fd)
   integral_trapezoidal = integrate_trapezoidal(func_to_int1 * func_to_int2, 0.1d0)

   ! In fact, the trapezoidal integration has a problem here
   call assert_comparable(integral_parseval, integral_trapezoidal, 1d-5,  &
                          'Integral of (sin(2pix)-1)*cos(1.5pix) from 1 to 2')

   call kernel%freeme()

end subroutine test_integrate_parseval_complex
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine test_tabulate_kernel()
   character(len=6)  :: model_param
   logical           :: needs_basekernel(6)
   character(len=32) :: strain_type

   ! Taking the kernel definition from the Fichtner book and checking them 
   ! one by one. N.B. They might not be applicable to our problem, but this
   ! is a separate question

   ! Lambda
   model_param = 'lam'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_true(needs_basekernel(1), model_param//' kernel needs lambda base kernel')
   call assert_false(needs_basekernel(2), model_param//' kernel needs no mu base kernel')
   call assert_false(needs_basekernel(3), model_param//' kernel needs no rho base kernel')
   call assert_false(needs_basekernel(4), model_param//' kernel needs no a base kernel')
   call assert_false(needs_basekernel(5), model_param//' kernel needs no b base kernel')
   call assert_false(needs_basekernel(6), model_param//' kernel needs no c base kernel')

   ! Mu 
   model_param = 'mu'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_false(needs_basekernel(1), model_param//' kernel needs no lambda base kernel')
   call assert_true(needs_basekernel(2), model_param//' kernel needs mu base kernel')
   call assert_false(needs_basekernel(3), model_param//' kernel needs no rho base kernel')
   call assert_false(needs_basekernel(4), model_param//' kernel needs no a base kernel')
   call assert_false(needs_basekernel(5), model_param//' kernel needs no b base kernel')
   call assert_false(needs_basekernel(6), model_param//' kernel needs no c base kernel')

   ! vp 
   model_param = 'vp'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_true(needs_basekernel(1), model_param//' kernel needs lambda base kernel')
   call assert_false(needs_basekernel(2), model_param//' kernel needs no mu base kernel')
   call assert_false(needs_basekernel(3), model_param//' kernel needs no rho base kernel')
   call assert_false(needs_basekernel(4), model_param//' kernel needs no a base kernel')
   call assert_false(needs_basekernel(5), model_param//' kernel needs no b base kernel')
   call assert_false(needs_basekernel(6), model_param//' kernel needs no c base kernel')

   ! vs 
   model_param = 'vs'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_true(needs_basekernel(1), model_param//' kernel needs lambda base kernel')
   call assert_true(needs_basekernel(2), model_param//' kernel needs mu base kernel')
   call assert_false(needs_basekernel(3), model_param//' kernel needs no rho base kernel')
   call assert_false(needs_basekernel(4), model_param//' kernel needs no a base kernel')
   call assert_false(needs_basekernel(5), model_param//' kernel needs no b base kernel')
   call assert_false(needs_basekernel(6), model_param//' kernel needs no c base kernel')

   ! vpv 
   model_param = 'vpv'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_true(needs_basekernel(1), model_param//' kernel needs lambda base kernel')
   call assert_false(needs_basekernel(2), model_param//' kernel needs no mu base kernel')
   call assert_false(needs_basekernel(3), model_param//' kernel needs no rho base kernel')
   call assert_true(needs_basekernel(4), model_param//' kernel needs a base kernel')
   call assert_false(needs_basekernel(5), model_param//' kernel needs no b base kernel')
   call assert_true(needs_basekernel(6), model_param//' kernel needs c base kernel')

   ! vsv 
   model_param = 'vsv'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_false(needs_basekernel(1), model_param//' kernel needs no lambda base kernel')
   call assert_false(needs_basekernel(2), model_param//' kernel needs no mu base kernel')
   call assert_false(needs_basekernel(3), model_param//' kernel needs no rho base kernel')
   call assert_false(needs_basekernel(4), model_param//' kernel needs no a base kernel')
   call assert_true(needs_basekernel(5), model_param//' kernel needs b base kernel')
   call assert_false(needs_basekernel(6), model_param//' kernel needs no c base kernel')

   ! vph 
   model_param = 'vph'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_false(needs_basekernel(1), model_param//' kernel needs no lambda base kernel')
   call assert_false(needs_basekernel(2), model_param//' kernel needs no mu base kernel')
   call assert_false(needs_basekernel(3), model_param//' kernel needs no rho base kernel')
   call assert_true(needs_basekernel(4), model_param//' kernel needs a base kernel')
   call assert_false(needs_basekernel(5), model_param//' kernel needs no b base kernel')
   call assert_true(needs_basekernel(6), model_param//' kernel needs c base kernel')

   ! vsh 
   model_param = 'vsh'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_true(needs_basekernel(1), model_param//' kernel needs lambda base kernel')
   call assert_true(needs_basekernel(2), model_param//' kernel needs mu base kernel')
   call assert_false(needs_basekernel(3), model_param//' kernel needs no rho base kernel')
   call assert_false(needs_basekernel(4), model_param//' kernel needs no a base kernel')
   call assert_true(needs_basekernel(5), model_param//' kernel needs b base kernel')
   call assert_true(needs_basekernel(6), model_param//' kernel needs c base kernel')

   ! eta 
   model_param = 'eta'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_true(needs_basekernel(1), model_param//' kernel needs lambda base kernel')
   call assert_false(needs_basekernel(2), model_param//' kernel needs no mu base kernel')
   call assert_false(needs_basekernel(3), model_param//' kernel needs no rho base kernel')
   call assert_true(needs_basekernel(4), model_param//' kernel needs a base kernel')
   call assert_false(needs_basekernel(5), model_param//' kernel needs no b base kernel')
   call assert_true(needs_basekernel(6), model_param//' kernel needs c base kernel')

   ! rho 
   model_param = 'rho'
   call tabulate_kernels(model_param, needs_basekernel, strain_type)
   call assert_true(needs_basekernel(1), model_param//' kernel needs lambda base kernel')
   call assert_true(needs_basekernel(2), model_param//' kernel needs mu base kernel')
   call assert_true(needs_basekernel(3), model_param//' kernel needs rho base kernel')
   call assert_true(needs_basekernel(4), model_param//' kernel needs a base kernel')
   call assert_true(needs_basekernel(5), model_param//' kernel needs b base kernel')
   call assert_true(needs_basekernel(6), model_param//' kernel needs c base kernel')


end subroutine test_tabulate_kernel
!-------------------------------------------------------------------------------

end module test_kernel
!=========================================================================================
