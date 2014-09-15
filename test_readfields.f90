!=========================================================================================
module test_readfields

   use global_parameters
   use readfields
   use ftnunit

   implicit none
   public
contains

!-----------------------------------------------------------------------------------------
subroutine test_readfields_load_seismogram()
   use type_parameter, only : parameter_type
   type(semdata_type)      :: sem_data
   type(parameter_type)    :: parameters

   call parameters%read_parameters('unit_tests/inparam_test')
   call parameters%read_source()
   call parameters%read_receiver()
   call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir, 100, 100)
   call sem_data%open_files()
   call sem_data%read_meshes()
   call sem_data%load_seismogram(parameters%receiver, parameters%source)
   
   call sem_data%close_files()

end subroutine test_readfields_load_seismogram
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_set_params()
   type(semdata_type)   :: sem_data
   character(len=512)   :: fwd_dir, bwd_dir

   fwd_dir = './test_wavefield/fwd'
   bwd_dir = './test_wavefield/bwd'
   
   call sem_data%set_params(fwd_dir, bwd_dir, 100, 100)

   call assert_equal_int(sem_data%nsim_fwd, 4, 'nsim_fwd == 4')
   call assert_equal_int(sem_data%nsim_bwd, 2, 'nsim_bwd == 2')

end subroutine test_readfields_set_params
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_open_files()
   type(semdata_type)   :: sem_data
   character(len=512)   :: fwd_dir, bwd_dir

   fwd_dir = './test_wavefield/fwd'
   bwd_dir = './test_wavefield/bwd'
   
   call sem_data%set_params(fwd_dir, bwd_dir, 100, 100) 

   call sem_data%open_files()

   call sem_data%close_files()

end subroutine  test_readfields_open_files
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
