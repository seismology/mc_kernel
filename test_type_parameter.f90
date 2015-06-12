!=========================================================================================
module test_type_parameter

  use global_parameters
  use type_parameter, only : parameter_type
  use ftnunit
  implicit none
  public
contains


!-----------------------------------------------------------------------------------------
subroutine test_parameter_reading
   use readfields,     only : semdata_type
   type(parameter_type)    :: parameters
   type(semdata_type)      :: sem_data
   logical                 :: ref_var(6)
   
   call parameters%read_parameters('unit_tests/inparam_test')
   call parameters%read_source()
   call parameters%read_receiver()

   ! Do checks on receivers
   ! Check, whether number of receivers is correct
   call assert_equal(parameters%nrec, 2, '2 receivers read in')

   call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir, 100, 100,  &
                            parameters%strain_type_fwd, parameters%source%depth)
   call sem_data%open_files()
   call sem_data%read_meshes()
   !call sem_data%build_kdtree()
   call sem_data%load_seismogram_rdbm(parameters%receiver, parameters%source)

   call parameters%read_filter(nomega=129, df=0.1d0)
     
   ! Do checks on filters
   call assert_equal(parameters%nfilter, 2, '2 Filters in input file')
   call assert_true(parameters%filter(1)%name=='Gabor_40', 'Filter 1, name: Gabor_40')
   call assert_true(parameters%filter(2)%name=='Gabor_30', 'Filter 2, name: Gabor_30')
   call assert_comparable(parameters%filter(1)%frequencies, [40.d0, 0.5d0, 0.d0, 0.d0], &
                          1d-10, 'Frequencies of filter 1 correct')
   call assert_comparable(parameters%filter(2)%frequencies, [30.d0, 0.5d0, 0.d0, 0.d0], &
                          1d-10, 'Frequencies of filter 2 correct')
   call assert_true(parameters%filter(1)%filterclass=='Gabor', 'Filter 1, type: Gabor')
   call assert_true(parameters%filter(2)%filterclass=='Gabor', 'Filter 2, type: Gabor')



   call parameters%read_kernel(sem_data, parameters%filter)

   ! Check, whether receiver and kernel file has been read in correctly
   call assert_equal(parameters%receiver(1)%nkernel, 3, 'Receiver 1 has 3 kernels')
   call assert_equal(parameters%receiver(2)%nkernel, 4, 'Receiver 2 has 4 kernel')

   ! Check, whether firstkernel and lastkernel are set correctly
   call assert_equal(parameters%receiver(1)%firstkernel, 1, 'Rec 1: first kernel: 1')
   call assert_equal(parameters%receiver(2)%firstkernel, 4, 'Rec 2: first kernel: 4')
   call assert_equal(parameters%receiver(1)%lastkernel, 3, 'Rec 1: last kernel: 1')
   call assert_equal(parameters%receiver(2)%lastkernel, 7, 'Rec 1: last kernel: 1')

   ! Check, whether needs_basekernel has been set correctly
   ! First receiver should need lambda, mu and A and C base kernels
   ref_var = [.true., .true., .false., .true., .false., .true.] 
   call assert_true(parameters%receiver(1)%needs_basekernel.eqv.ref_var, &
                    'First receiver should need lambda, mu, A and C base kernels')

   ref_var = [.true., .true., .true., .true., .true., .true.] 
   call assert_true(parameters%receiver(2)%needs_basekernel.eqv.ref_var, &
                    'Second receiver should need all base kernels')

   ! Check whether needs_basekernel has been set correctly on the kernels themselves

   ! Kernels of receiver 1
   ref_var = [.true., .false., .false., .false., .false., .false.] 
   call assert_true(parameters%kernel(1)%needs_basekernel.eqv.ref_var, &
                    'Kernel 1 should need only lambda base kernels (lambda)')

   ref_var = [.true., .false., .false., .true., .false., .true.] 
   call assert_true(parameters%kernel(2)%needs_basekernel.eqv.ref_var, &
                    'Kernel 2 should need lambda, a and c base kernels (eta)')

   ref_var = [.false., .true., .false., .false., .false., .false.] 
   call assert_true(parameters%kernel(3)%needs_basekernel.eqv.ref_var, &
                    'Kernel 3 should need only mu base kernels (mu)')

   ! Kernels of receiver 2
   ref_var = [.true., .true., .true., .true., .true., .true.] 
   call assert_true(parameters%kernel(4)%needs_basekernel.eqv.ref_var, &
                    'Kernel 4 should need all base kernels (rho)')

   ref_var = [.true., .false., .false., .false., .false., .false.] 
   call assert_true(parameters%kernel(5)%needs_basekernel.eqv.ref_var, &
                    'Kernel 5 should need only lambda base kernels (vp)')

   ref_var = [.true., .true., .false., .false., .false., .false.] 
   call assert_true(parameters%kernel(6)%needs_basekernel.eqv.ref_var, &
                    'Kernel 6 should need lambda and mu base kernels (vs)')

   ref_var = [.true., .false., .false., .true., .false., .true.] 
   call assert_true(parameters%kernel(7)%needs_basekernel.eqv.ref_var, &
                    'Kernel 7 should need lambda, a and c base kernels (eta)')


end subroutine test_parameter_reading
!-----------------------------------------------------------------------------------------

end module test_type_parameter                 
!=========================================================================================
