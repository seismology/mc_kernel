!=========================================================================================
program test_resampling

  use commpi
  use global_parameters
  use resampling

  implicit none

  type(resampling_type)                :: resamp

  integer                                       :: ntimes_in, ntimes_out, ntraces
  integer                                       :: i, j
  real(kind=dp), dimension(:,:), allocatable    :: data_in, data_out
  real(kind=dp), dimension(:,:), allocatable    :: T_in, T_out
  real(kind=dp)                                 :: dt_in, dt_out, f_ref

  ntimes_in = 100
  ntimes_out = 150
  ntraces = 1

  allocate(data_in(1:ntimes_in, 1:ntraces))
  allocate(data_out(1:ntimes_out, 1:ntraces))
  allocate(T_in(1:ntimes_in, 1:ntraces))
  allocate(T_out(1:ntimes_out, 1:ntraces))

  dt_in = 0.2
  dt_out = dt_in * ntimes_in / ntimes_out

  write(6,*) dt_out

  do i = 1, ntimes_in
     T_in(i,:) = dt_in * (i - 1)
  end do

  do i = 1, ntimes_out
     T_out(i,:) = dt_out * (i - 1)
  end do

  f_ref = .5
  data_in = sin(2 * pi * T_in * f_ref)

  do i = 1, ntimes_in
     write(222,*) T_in(i,1), data_in(i,:)
  enddo

  call resamp%init(ntimes_in, ntimes_out, ntraces)

  call resamp%resample(data_in, data_out)

  do i = 1, ntimes_out
     write(333,*) T_out(i,1), data_out(i,:)
  enddo

end program
!=========================================================================================
