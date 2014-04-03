program rdbm

  use readfields

  write(6,*) 'kuckuck'

contains


!-----------------------------------------------------------------------------
subroutine start_clock
  !
  ! Driver routine to start the timing, using the clocks_mod module.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  use global_parameters
  use clocks_mod, only : clock_id, clocks_init
  
  implicit none
  
  character(len=8)  :: mydate
  character(len=10) :: mytime

  call date_and_time(mydate,mytime) 
  write(lu_out,11) mydate(5:6), mydate(7:8), mydate(1:4), mytime(1:2), mytime(3:4)

11 format('     Kerner started on ', A2,'/',A2,'/',A4,' at ', A2,'h ',A2,'min',/)


  call clocks_init(0)

  id_fft         = clock_id('FFT routines')
  id_fwd         = clock_id('Reading fwd field')
  id_netcdf      = clock_id(' - NetCDF routines')
  id_rotate      = clock_id(' - Rotate fields')
  id_buffer      = clock_id(' - Buffer routines')
  id_bwd         = clock_id('Reading bwd field')
  id_mc          = clock_id('Monte Carlo routines')
  id_filter_conv = clock_id('Filtering and convolution')
  id_inv_mesh    = clock_id('Inversion mesh routines')
  id_kernel      = clock_id('Kernel routines')
  id_init        = clock_id('Initialization per task')
  id_mpi         = clock_id('MPI communication with Master')

end subroutine start_clock
!=============================================================================

!-----------------------------------------------------------------------------
subroutine end_clock 
  !
  ! Wapper routine to end timing and display clock informations.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  use global_parameters, only : lu_out
  use clocks_mod,        only : clocks_exit

  implicit none

  write(lu_out,*)
  write(lu_out,"(10x,'Summary of timing measurements:')")
  write(lu_out,*)

  call clocks_exit(0)

  write(lu_out,*)

end subroutine end_clock
!=============================================================================
end program
