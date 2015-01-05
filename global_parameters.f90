!=========================================================================================
module global_parameters

  implicit none
  public
  integer, parameter         :: sp = selected_real_kind(6, 37)
  integer, parameter         :: dp = selected_real_kind(15, 307)
  integer, parameter         :: qp = selected_real_kind(33, 4931)

  integer,parameter          :: int4 = selected_int_kind(5)
  integer,parameter          :: long = selected_int_kind(15)

  real(kind=dp), parameter   :: pi = 3.1415926535898D0
  real(kind=dp), parameter   :: deg2rad = pi / 180.d0
  real(kind=dp), parameter   :: rad2deg = 180.d0 / pi
  integer                    :: verbose = 0

  integer, parameter         :: WORKTAG = 1
  integer, parameter         :: DIETAG  = 2
  
  logical                    :: master, firstslave
  logical                    :: testing = .false. !< Set to true only for unit test, 
                                                  !! because some routines require action
                                                  !! from master or slave, which would not 
                                                  !! be tested otherwise
  integer, protected         :: myrank, nproc
  integer, protected         :: lu_out !< Logical unit for output. 
                                       !! 6 (Screen) for master
                                       !! File 'OUTPUT_#rank' for slaves

  integer                    :: id_read, id_fft, id_fwd, id_bwd, id_mc, id_mpi,&
                                id_filter_conv, id_inv_mesh, id_kernel, id_init,   &
                                id_buffer, id_netcdf, id_rotate, id_load_strain,   &
                                id_kdtree, id_calc_strain, id_find_point_fwd,      &
                                id_find_point_bwd, id_lagrange, id_int_model,      &
                                id_read_params, id_create_tasks, id_get_next_task, &
                                id_extract, id_write_kernel, id_mult_kernel,       &
                                id_init_fft
  integer                    :: id_load, id_resamp, id_out

  contains


!----------------------------------------------------------------------------------------
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

!----------------------------------------------------------------------------------------
subroutine set_myrank(myrank_value)
  integer, intent(in)   :: myrank_value

  myrank = myrank_value
end subroutine
!----------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
subroutine set_nproc(nproc_value)
  integer, intent(in)   :: nproc_value

  nproc = nproc_value

end subroutine
!----------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
subroutine set_lu_out(lu_out_value)
  integer, intent(in)   :: lu_out_value

  lu_out = lu_out_value

end subroutine
!----------------------------------------------------------------------------------------


end module
!=========================================================================================
