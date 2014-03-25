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
  integer                    :: verbose 

  integer, parameter         :: WORKTAG = 1
  integer, parameter         :: DIETAG  = 2
  
  logical                    :: master, firstslave
  integer                    :: myrank, nproc
  integer                    :: lu_out !< Logical unit for output. 
                                       !! 6 (Screen) for master
                                       !! file 'OUTPUT_#rank' for slaves

  integer                    :: id_read, id_fft, id_fwd, id_bwd, id_mc, id_mpi,&
                                id_filter_conv, id_inv_mesh, id_kernel, id_init, &
                                id_buffer, id_netcdf, id_rotate
  contains


SUBROUTINE init_random_seed()
   INTEGER :: i, n, clock
   INTEGER, DIMENSION(:), ALLOCATABLE :: seed
                                                  
   CALL RANDOM_SEED(size = n)
   ALLOCATE(seed(n))

   CALL SYSTEM_CLOCK(COUNT=clock)

   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)

   DEALLOCATE(seed)
END SUBROUTINE


end module
!=========================================================================================
