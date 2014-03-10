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
  integer                    :: myrank, verbose 

  integer, parameter         :: WORKTAG = 1
  integer, parameter         :: DIETAG  = 2
  
  integer                    :: lu_out !< Logical unit for output. 
                                       !! 6 (Screen) for master
                                       !! file 'OUTPUT_#rank' for slaves
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
