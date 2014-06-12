!=========================================================================================
module global_parameters

  implicit none
  public
  integer, parameter         :: sp = selected_real_kind(6, 37)
  integer, parameter         :: dp = selected_real_kind(15, 307)
  integer, parameter         :: qp = selected_real_kind(33, 4931)

  real(kind=dp), parameter   :: pi = 3.1415926535898D0
  real(kind=dp), parameter   :: deg2rad = pi / 180.d0
  real(kind=dp), parameter   :: rad2deg = 180.d0 / pi
  integer                    :: verbose = 0

  integer, parameter         :: WORKTAG = 1
  integer, parameter         :: DIETAG  = 2
  
  logical                    :: master, firstslave
  integer                    :: myrank, nproc
  integer                    :: lu_out !< LOgical unit for output. 
                                       !! 6 (Screen) for master
                                       !! FIle 'OUTPUT_#rank' for slaves

  integer                    :: id_read, id_fft, id_fwd, id_bwd, id_mc, id_mpi,&
                                id_filter_conv, id_inv_mesh, id_kernel, id_init, &
                                id_buffer, id_netcdf, id_rotate
  contains


!-----------------------------------------------------------------------------------------
subroutine init_random_seed()
   integer :: i, n, clock
   integer, dimension(:), allocatable :: seed
                                                  
   call random_seed(size = n)
   allocate(seed(n))

   call system_clock(count=clock)

   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call random_seed(put = seed)

   deallocate(seed)
end subroutine
!-----------------------------------------------------------------------------------------


end module
!=========================================================================================
