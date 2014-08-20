!=========================================================================================
module test_resampling

  use global_parameters
  use resampling
  use ftnunit
  implicit none
  public
contains

!-----------------------------------------------------------------------------------------
subroutine test_resampling_const
  
  type(resampling_type)                         :: resamp

  integer                                       :: ntimes_in, ntimes_out, ntraces
  integer                                       :: i
  integer                                       :: myunit
  real(kind=dp), dimension(:,:), allocatable    :: data_in, data_out, data_ref
  real(kind=dp), dimension(:,:), allocatable    :: T_in, T_out
  real(kind=dp)                                 :: dt_in, dt_out
  character(len=32)                             :: fnam

  ntimes_in = 100
  ntimes_out = 200
  ntraces = 1

  allocate(data_in(1:ntimes_in, 1:ntraces))
  allocate(data_out(1:ntimes_out, 1:ntraces))
  allocate(data_ref(1:ntimes_out, 1:ntraces))
  allocate(T_in(1:ntimes_in, 1:ntraces))
  allocate(T_out(1:ntimes_out, 1:ntraces))

  dt_in = 0.1
  dt_out = dt_in * ntimes_in / ntimes_out

  do i = 1, ntimes_in
     T_in(i,:) = dt_in * (i - 1)
  end do

  do i = 1, ntimes_out
     T_out(i,:) = dt_out * (i - 1)
  end do

  data_in = 1
  data_ref = 1

  fnam = 'unit_tests/test_resampling_const_in'
  open(newunit=myunit, file=fnam, status='replace')

  do i = 1, ntimes_in
     write(myunit,*) T_in(i,1), data_in(i,:)
  enddo

  close(myunit)

  call resamp%init(ntimes_in, ntimes_out, ntraces)

  call resamp%resample(data_in, data_out)

  do i=1, ntraces
     call assert_comparable_real1d(real(data_out(:,i), sp), &
                                   real(data_ref(:,i), sp), &
                                   1e-5, 'resampling constant function')
  enddo

  fnam = 'unit_tests/test_resampling_const_out'
  open(newunit=myunit, file=fnam, status='replace')

  do i = 1, ntimes_out
     write(myunit,*) T_out(i,1), data_out(i,:)
  enddo

  close(myunit)

  call resamp%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_resampling_const_ntraces
  
  type(resampling_type)                         :: resamp

  integer                                       :: ntimes_in, ntimes_out, ntraces
  integer                                       :: i
  integer                                       :: myunit
  real(kind=dp), dimension(:,:), allocatable    :: data_in, data_out, data_ref
  real(kind=dp), dimension(:,:), allocatable    :: T_in, T_out
  real(kind=dp)                                 :: dt_in, dt_out
  character(len=64)                             :: fnam

  ntimes_in = 100
  ntimes_out = 200
  ntraces = 5

  allocate(data_in(1:ntimes_in, 1:ntraces))
  allocate(data_out(1:ntimes_out, 1:ntraces))
  allocate(data_ref(1:ntimes_out, 1:ntraces))
  allocate(T_in(1:ntimes_in, 1:ntraces))
  allocate(T_out(1:ntimes_out, 1:ntraces))

  dt_in = 0.1
  dt_out = dt_in * ntimes_in / ntimes_out

  do i = 1, ntimes_in
     T_in(i,:) = dt_in * (i - 1)
  end do

  do i = 1, ntimes_out
     T_out(i,:) = dt_out * (i - 1)
  end do

  do i=1, ntraces
     data_in(:,i) = i
     data_ref(:,i) = i
  enddo

  fnam = 'unit_tests/test_resampling_const_ntraces_in'
  open(newunit=myunit, file=fnam, status='replace')

  do i = 1, ntimes_in
     write(myunit,*) T_in(i,1), data_in(i,:)
  enddo

  close(myunit)

  call resamp%init(ntimes_in, ntimes_out, ntraces)

  call resamp%resample(data_in, data_out)

  do i=1, ntraces
     call assert_comparable_real1d(real(data_out(:,i), sp), &
                                   real(data_ref(:,i), sp), 1e-5, &
                                   'resampling constant function with multiple traces')
  enddo

  fnam = 'unit_tests/test_resampling_const_ntraces_out'
  open(newunit=myunit, file=fnam, status='replace')

  do i = 1, ntimes_out
     write(myunit,*) T_out(i,1), data_out(i,:)
  enddo

  close(myunit)

  call resamp%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_resampling_triangle
  
  type(resampling_type)                         :: resamp

  integer                                       :: ntimes_in, ntimes_out, ntraces
  integer                                       :: i
  integer                                       :: myunit
  real(kind=dp), dimension(:,:), allocatable    :: data_in, data_out, data_ref
  real(kind=dp), dimension(:,:), allocatable    :: T_in, T_out
  real(kind=dp)                                 :: dt_in, dt_out
  character(len=32)                             :: fnam

  ntimes_in = 100
  ntimes_out = 200
  ntraces = 1

  allocate(data_in(1:ntimes_in, 1:ntraces))
  allocate(data_out(1:ntimes_out, 1:ntraces))
  allocate(data_ref(1:ntimes_out, 1:ntraces))
  allocate(T_in(1:ntimes_in, 1:ntraces))
  allocate(T_out(1:ntimes_out, 1:ntraces))

  dt_in = 0.1
  dt_out = dt_in * ntimes_in / ntimes_out

  do i = 1, ntimes_in
     T_in(i,:) = dt_in * (i - 1)
  end do

  do i = 1, ntimes_out
     T_out(i,:) = dt_out * (i - 1)
  end do

  data_in(1:ntimes_in/2+1,1) = T_in(1:ntimes_in/2+1,1)
  data_in(ntimes_in/2+1:ntimes_in,1) = T_in(ntimes_in/2:1:-1,1)

  data_ref(1:ntimes_out/2+1,1) = T_out(1:ntimes_out/2+1,1)
  data_ref(ntimes_out/2+1:,1) = T_out(ntimes_out/2:1:-1,1)


  fnam = 'unit_tests/test_resampling_triangle_in'
  open(newunit=myunit, file=fnam, status='replace')

  do i = 1, ntimes_in
     write(myunit,*) T_in(i,1), data_in(i,:)
  enddo

  close(myunit)

  call resamp%init(ntimes_in, ntimes_out, ntraces)

  call resamp%resample(data_in, data_out)

  do i=1, ntraces
     call assert_comparable_real1d(real(data_out(:,i), sp) + 10, &
                                   real(data_ref(:,i), sp) + 10, &
                                   1e-2, 'resampling triangle function')
  enddo

  fnam = 'unit_tests/test_resampling_triangle_out'
  open(newunit=myunit, file=fnam, status='replace')

  do i = 1, ntimes_out
     write(myunit,*) T_out(i,1), data_out(i,:), data_ref(i,:)
  enddo

  close(myunit)

  call resamp%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_resampling_timeshift_triangle
  
  use source_class
  type(resampling_type)                         :: resamp

  type(src_param_type)                          :: source(1)

  integer                                       :: ntimes_in, ntimes_out, ntraces
  integer                                       :: i
  integer                                       :: myunit
  real(kind=dp), dimension(:,:), allocatable    :: data_in, data_out, data_ref
  real(kind=dp), dimension(:,:), allocatable    :: T_in, T_out
  real(kind=dp)                                 :: dt_in, dt_out
  character(len=64)                             :: fnam

  ntimes_in = 100
  ntimes_out = 100
  ntraces = 1

  allocate(data_in(1:ntimes_in, 1:ntraces))
  allocate(data_out(1:ntimes_out, 1:ntraces))
  allocate(data_ref(1:ntimes_out, 1:ntraces))
  allocate(T_in(1:ntimes_in, 1:ntraces))
  allocate(T_out(1:ntimes_out, 1:ntraces))

  dt_in = 0.1
  dt_out = dt_in * ntimes_in / ntimes_out

  do i = 1, ntimes_in
     T_in(i,:) = dt_in * (i - 1)
  end do

  do i = 1, ntimes_out
     T_out(i,:) = dt_out * (i - 1)
  end do

  data_in(1:ntimes_in/2+1,1) = T_in(1:ntimes_in/2+1,1)
  data_in(ntimes_in/2+1:ntimes_in,1) = T_in(ntimes_in/2:1:-1,1)

  data_ref = 4.9d0 - data_in

  fnam = 'unit_tests/test_resampling_timeshift_triangle_in'
  open(newunit=myunit, file=fnam, status='replace')

  do i = 1, ntimes_in
     write(myunit,*) T_in(i,1), data_in(i,:)
  enddo

  close(myunit)

  call source(1)%init(0d0, 0d0, [0d0, 0d0, 0d0, 0d0, 0d0, 0d0], 0d0)
  source(1)%shift_time_sample = 5d0 / dt_in

  call resamp%init(ntimes_in, ntimes_out, ntraces)

  call resamp%resample_timeshift(data_in, data_out, source)

  do i=1, ntraces
     call assert_comparable_real1d(real(data_out(:,i), sp) + 10, &
                                   real(data_ref(:,i), sp) + 10, &
                                   1e-2, 'timeshifting triangle function')
  enddo

  fnam = 'unit_tests/test_resampling_timeshift_triangle_out'
  open(newunit=myunit, file=fnam, status='replace')

  do i = 1, ntimes_out
     write(myunit,*) T_out(i,1), data_out(i,:), data_ref(i,:)
  enddo

  close(myunit)

  call resamp%freeme()

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
