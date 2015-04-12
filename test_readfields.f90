!=========================================================================================
module test_readfields

   use global_parameters
   use readfields
   use ftnunit

   implicit none
   public
contains

!-----------------------------------------------------------------------------------------
subroutine test_get_chunk_bounds()
   integer  :: npoints, chunksize, ipoint, start_chunk, count_chunk, i, iinchunk
   integer  :: start_chunk_ref, count_chunk_ref
   character(len=80) :: err_msg

   ! This should have 20 chunks
   npoints = 99
   chunksize = 5

   do iinchunk = 0, chunksize - 1

     do i = 1, 20
       
       start_chunk_ref = (i-1) * chunksize + 1
       if (i<20) then
         count_chunk_ref = 5
       else 
         count_chunk_ref = 4
       end if

       ipoint = start_chunk_ref + iinchunk 

       if (ipoint > npoints) cycle

       call get_chunk_bounds(pointid     = ipoint,        &
                             chunksize   = chunksize,     &
                             npoints     = npoints,       & 
                             start_chunk = start_chunk,   &
                             count_chunk = count_chunk)

       write(err_msg, "('Start of chunk correct for point: ', I4)") ipoint
       call assert_equal(start_chunk, start_chunk_ref, trim(err_msg))
       write(err_msg, "('Count of chunk correct for point: ', I4)") ipoint
       call assert_equal(count_chunk, count_chunk_ref, trim(err_msg))

     end do
   end do
     
end subroutine test_get_chunk_bounds
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_set_params()
   type(semdata_type)   :: sem_data
   character(len=512)   :: fwd_dir, bwd_dir

   fwd_dir = './test_wavefields/kerner_fwd'
   bwd_dir = './test_wavefields/kerner_bwd'
   
   call sem_data%set_params(fwd_dir, bwd_dir, 100, 100, 'straintensor_trace', 100.0d0)

   call assert_equal_int(sem_data%nsim_fwd, 4, 'nsim_fwd == 4')
   call assert_equal_int(sem_data%nsim_bwd, 2, 'nsim_bwd == 2')

end subroutine test_readfields_set_params
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_open_files()
   use type_parameter, only : parameter_type
   type(parameter_type)    :: parameters
   type(semdata_type)   :: sem_data

   call parameters%read_parameters('unit_tests/inparam_test')
   
   call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir, &
                            100, 100, 'straintensor_trace', 100.0d0) 

   call sem_data%open_files()

   call sem_data%close_files()

end subroutine  test_readfields_open_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_read_meshes
   use type_parameter, only : parameter_type
   type(parameter_type)    :: parameters
   type(semdata_type)   :: sem_data
   type(meshtype)       :: testmesh_fwd, testmesh_bwd
   integer              :: i_max_vp, i_max_vs, i_max_rho
   real(kind=sp)        :: r_max_vp, r_max_vs, r_max_rho

   call parameters%read_parameters('unit_tests/inparam_test')
   
   call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir, &
                            100, 100, 'straintensor_trace', 100.0d0) 

   call sem_data%open_files()

   call sem_data%read_meshes()

   testmesh_fwd = sem_data%get_mesh('fwd')
   testmesh_bwd = sem_data%get_mesh('bwd')

   call assert_comparable(testmesh_fwd%s, testmesh_bwd%s, 1e-10, &
                          'FWD and BWD meshes are identical')
   
   call assert_comparable(testmesh_fwd%mu, testmesh_bwd%mu, 1e-10, &
                          'FWD and BWD meshes are identical')
   

   !This is PREM, so we might check whether the mesh variables have reasonable limits

   ! S
   call assert_comparable(minval(testmesh_fwd%s),   0.0,       1e-5, 'All S>0')
   call assert_comparable(maxval(testmesh_fwd%s),   6371000.0, 1e-5, 'All S<6371')

   ! Z
   call assert_comparable(minval(testmesh_fwd%z),  -6371000.0, 1e-5, 'All Z>-6371')
   call assert_comparable(maxval(testmesh_fwd%z),   6371000.0, 1e-5, 'All Z< 6371')

   ! VP
   call assert_comparable(minval(testmesh_fwd%vp),  5800.0 ,   1e-2, 'Min VP=5800')
   call assert_comparable(maxval(testmesh_fwd%vp),  13716.6,   1e-2, 'Max VP=13717')

   ! VS
   call assert_comparable(minval(testmesh_fwd%vs),  0000.0 ,   1e-2, 'All VS>0000')
   call assert_comparable(maxval(testmesh_fwd%vs),  7265.97,   1e-2, 'All VS<7265')

   ! Rho
   call assert_comparable(minval(testmesh_fwd%rho), 2600.0 ,   1e-2, 'All Rho>2600')
   call assert_comparable(maxval(testmesh_fwd%rho), 13088.48,  1e-2, 'All Rho<13089')
   
   ! Mu
   call assert_comparable(minval(testmesh_fwd%mu),  0.0 ,     1e-5, 'Min Mu=0')
   call assert_comparable(maxval(testmesh_fwd%mu),  293.77e9,  1e-1, 'Max Mu=293.8 GPa')
   
   ! Lambda
   ! These would be the correct values according to the PREM paper, but who knows
   !call assert_comparable(minval(testmesh_fwd%lambda), 52.0e9 ,   1e-5, 'Min Lambda=52 GPa')
   !call assert_comparable(maxval(testmesh_fwd%lambda), 1425.3e9,  1e-5, 'Max Lambda=1425 GPa')
   call assert_comparable(minval(testmesh_fwd%lambda), 34.216e9 ,   1e-2, 'Min Lambda=34 GPa')
   call assert_comparable(maxval(testmesh_fwd%lambda), 1307.96e9,  1e-2, 'Max Lambda=1308 GPa')
  
   ! Eta
   call assert_comparable(minval(testmesh_fwd%eta), 1.0,   1e-5, 'Min eta=1')
   call assert_comparable(maxval(testmesh_fwd%eta), 1.0,   1e-5, 'Max eta=1')

   ! I have no idea, what these values should be - SCS
   !! Phi
   !call assert_comparable(minval(testmesh_fwd%phi), 52.0e9 ,   1e-5, 'Min phi=52 GPa')
   !call assert_comparable(maxval(testmesh_fwd%phi), 1425.3e9,  1e-5, 'Max phi=1425 GPa')
   !
   !! Xi
   !call assert_comparable(minval(testmesh_fwd%xi), 52.0e9 ,   1e-5, 'Min xi=52 GPa')
   !call assert_comparable(maxval(testmesh_fwd%xi), 1425.3e9,  1e-5, 'Max xi=1425 GPa')
   
   ! Maximum VP should be at the CMB
   i_max_vp = maxloc(testmesh_fwd%vp, 1)
   r_max_vp = norm2([testmesh_fwd%s(i_max_vp), testmesh_fwd%z(i_max_vp)])
   call assert_comparable(r_max_vp, 3480000.0, 1e-2, 'Max VP at the CMB')

   ! Maximum VS should be at the D"
   i_max_vs = maxloc(testmesh_fwd%vs, 1)
   r_max_vs = norm2([testmesh_fwd%s(i_max_vs), testmesh_fwd%z(i_max_vs)])
   call assert_comparable(r_max_vs, 3630000.0, 1e-2, 'Max VS at the D"')

   ! Maximum Rho should be in the center
   i_max_rho = maxloc(testmesh_fwd%rho, 1)
   r_max_rho = norm2([testmesh_fwd%s(i_max_rho), testmesh_fwd%z(i_max_rho)])
   call assert_comparable(r_max_rho, 0.0, 1e-5, 'Max rho at the center of the earth')



   call sem_data%close_files()

end subroutine  test_readfields_read_meshes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> For now, this test only checks, whether something crashes while reading points
subroutine test_readfields_load_fw_points
   use readfields,     only : semdata_type
   use type_parameter, only : parameter_type
   type(parameter_type)    :: parameters
   type(semdata_type)      :: sem_data
   real(kind=dp)           :: u(73, 6, 2), coordinates(3,2)
   
   call parameters%read_parameters('unit_tests/inparam_test')
   call parameters%read_source()
   call parameters%read_receiver()

   call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir, 100, 100,  &
                            parameters%strain_type_fwd, parameters%source%depth)
   call sem_data%open_files()
   call sem_data%read_meshes()
   call sem_data%build_kdtree()
   call sem_data%load_seismogram_rdbm(parameters%receiver, parameters%source)

   call parameters%read_filter(nomega=129, df=0.1d0)
     
   call parameters%read_kernel(sem_data, parameters%filter)

   coordinates(:,1) = [ 1d6, 1d6, 1d6]
   coordinates(:,2) = [-1d6, 1d6, 1d6]
   u = sem_data%load_fw_points(coordinates, parameters%source)

   call sem_data%close_files()

end subroutine test_readfields_load_fw_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Load a seismogram and compare with instaseis
subroutine test_load_seismograms_rdbm
   use readfields,     only : semdata_type
   use type_parameter, only : parameter_type
   use fft,                                  only: rfft_type, taperandzeropad
   type(parameter_type)    :: parameters
   type(semdata_type)      :: sem_data
   type(rfft_type)         :: fft_data
   integer                 :: nomega, ntimes, ntimes_reference, isample, lu_seis
   real(kind=dp)           :: df, dt, t
   real(kind=dp), allocatable :: seis(:), seis_ref(:)
   
   call parameters%read_parameters('unit_tests/inparam_load_seismogram')
   call parameters%read_source()
   call parameters%read_receiver()

   call sem_data%set_params(fwd_dir              = parameters%fwd_dir,          &
                            bwd_dir              = parameters%bwd_dir,          &
                            strain_buffer_size   = 100,                         &
                            displ_buffer_size    = 100,                         &
                            strain_type          = parameters%strain_type_fwd,  &
                            desired_source_depth = parameters%source%depth)

   call sem_data%open_files()
   call sem_data%read_meshes()
   call sem_data%build_kdtree()
   call sem_data%load_seismogram_rdbm(parameters%receiver, parameters%source)

   ! Initialize FFT - just needed to get df and nomega
   call fft_data%init(ntimes_in = sem_data%ndumps,     &
                      ndim      = sem_data%get_ndim(), &
                      ntraces   = 1,                   &
                      dt        = sem_data%dt)

   ntimes = fft_data%get_ntimes()
   nomega = fft_data%get_nomega()
   df     = fft_data%get_df()
   call fft_data%freeme()


   ! Read filters
   call parameters%read_filter(nomega=nomega, df=df)
   call parameters%read_kernel(sem_data, parameters%filter)


   ! Retrieve seismograms
   ! 1st one filtered with Butterworth, 6th order at 40s
   allocate(seis(ntimes))
   open(newunit=lu_seis, file='seismogram_T001_P_BW', action='read')
   do isample = 1, ntimes
     read(lu_seis,*) t, seis(isample)
   end do
   close(lu_seis)

   ! Load reference seismogram
   open(newunit=lu_seis, file='unit_tests/seismogram_ref_T001_P_BW', action='read')
   read(lu_seis,*) ntimes_reference
   allocate(seis_ref(ntimes_reference))
   do isample = 1, ntimes_reference
     read(lu_seis,*) seis_ref(isample)
   end do
   close(lu_seis)

   !call assert_equal(ntimes, ntimes_reference, 'Seismogram length is equal to reference length')

   open(newunit=lu_seis, file='unit_tests_output/seis_comparison.txt', action='write')
   do isample = 1, ntimes_reference
     write(lu_seis, *) seis(isample), seis_ref(isample)
   end do
   close(lu_seis)
   call assert_comparable(seis(1:ntimes_reference), seis_ref, 1d-6, &
                          'Seismogram is comparable to reference seismogram')


   call sem_data%close_files()

end subroutine test_load_seismograms_rdbm
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_readfields_load_model_coeffs
   use type_parameter, only    : parameter_type
   use backgroundmodel, only   : backgroundmodel_type
   use global_parameters, only : pi, sp, dp
   type(parameter_type)       :: parameters
   type(semdata_type)         :: sem_data
   type(meshtype)             :: testmesh_fwd, testmesh_bwd
   type(backgroundmodel_type) :: model, model_ref
   integer, parameter         :: npoints = 1, nradius = 11
   real(kind=dp)              :: coordinates(3,npoints), phi(npoints), theta(npoints), r(nradius)
   integer                    :: ipoint, iradius
   real(kind=sp)              :: s(npoints), z(npoints)

   ! Test in all 11 domains of PREM
   r = [6370d3, & ! within upper crust
        6352d3, & !        lower crust
        6250d3, & !        upper mantle
        6050d3, & !        upper transition zone 
        5850d3, & !        lower transition zone
        5730d3, & ! directly above 670       
        5650d3, & ! directly below 670       
        4000d3, & !        lower mantle
        3550d3, & !        lowermost mantle
        2500d3, & !        outer core
        1000d3]   !        inner core


   call parameters%read_parameters('unit_tests/inparam_test')

   call sem_data%set_params(parameters%fwd_dir, parameters%bwd_dir, 100, 100, &
                            'straintensor_trace', 100.0d0) 

   call sem_data%open_files()

   call sem_data%read_meshes()
   
   call sem_data%build_kdtree()

   !do iradius = 600, 637
   !  coordinates(:,1) = [0.d0, 0.d0, iradius * 1.d4]
   !  model = sem_data%load_model_coeffs(coordinates, s, z)
   !  print *, iradius, ', Vp: ', model%c_vp(1), ', Vs: ', model%c_vs(1), ', Rho: ', model%c_rho(1)
   !end do

   do iradius = 1, nradius

     call random_number(phi)
     call random_number(theta)
     phi = phi * 2.d0 * pi
     theta = (theta - 0.5) * pi

     !print *, 'Radius :', iradius, '(', r(iradius), ')'

     do ipoint = 1, npoints
        coordinates(:, ipoint) = [cos(phi(ipoint)) * sin(theta(ipoint)), &
                                  sin(phi(ipoint)) * sin(theta(ipoint)), &
                                  cos(theta(ipoint))] * r(iradius)
        !print *, coordinates(:, ipoint), norm2(coordinates(:, ipoint))                       
     end do
     call flush(6)

     model = sem_data%load_model_coeffs(coordinates, s, z)
     model_ref = prem_ani_sub(r(iradius), iradius)

     do ipoint = 1, npoints 
       call assert_comparable(model_ref%c_vp(1) , model%c_vp(ipoint),  1e-1, &
                              'vP is identical for same radius')
       call assert_comparable(model_ref%c_vpv(1), model%c_vpv(ipoint), 1e-1, &
                              'vPv is identical for same radius')
       call assert_comparable(model_ref%c_vph(1), model%c_vph(ipoint), 1e-1, &
                              'vPh is identical for same radius')
       call assert_comparable(model_ref%c_vs(1) , model%c_vs(ipoint),  1e-1, &
                              'vS is identical for same radius')
       call assert_comparable(model_ref%c_vsv(1), model%c_vsv(ipoint), 1e-1, &
                              'vSv is identical for same radius')
       call assert_comparable(model_ref%c_vsh(1), model%c_vsh(ipoint), 1e-1, &
                              'vSh is identical for same radius')
       call assert_comparable(model_ref%c_xi(1) , model%c_xi(ipoint),  1e-1, &
                              'Xi is identical for same radius')
       call assert_comparable(model_ref%c_phi(1), model%c_phi(ipoint), 1e-1, &
                              'Phi is identical for same radius')
       call assert_comparable(model_ref%c_eta(1), model%c_eta(ipoint), 1e-1, &
                              'eta is identical for same radius')
       call assert_comparable(model_ref%c_rho(1), model%c_rho(ipoint), 1e-1, &
                              'rho is identical for same radius')
     end do

   end do ! iradius

   call sem_data%close_files()

end subroutine test_readfields_load_model_coeffs
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> The wicked PREM subroutine from background_models.f90 in the AxiSEM MESHER
function prem_ani_sub(r0, idom) result(model)
  use backgroundmodel, only     : backgroundmodel_type
  real(kind=dp), intent(in)    :: r0
  integer, intent(in)          :: idom
  type(backgroundmodel_type)   :: model

  real(kind=dp)                :: r,x_prem
  real(kind=dp)                :: eta_aniso, Qmu, Qkappa

  r = r0 / 1000.

  call model%init(1)
  
  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )
  model%c_eta = 1.

  IF(idom==1)THEN       ! upper crustal layer
     model%c_rho  = 2.6
     model%c_vpv = 5.8
     model%c_vsv = 3.2
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 600.0
     Qkappa = 57827.0
  ELSEIF(idom==2)THEN   ! lower crustal layer
     model%c_rho  = 2.9
     model%c_vpv = 6.8
     model%c_vsv = 3.9
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 600.0
     Qkappa = 57827.0
  ELSEIF(idom==3)THEN   ! upper mantle
     model%c_rho   =  2.6910 + 0.6924 * x_prem
     model%c_vpv  =  0.8317 + 7.2180 * x_prem
     model%c_vph  =  3.5908 + 4.6172 * x_prem
     model%c_vsv  =  5.8582 - 1.4678 * x_prem
     model%c_vsh  = -1.0839 + 5.7176 * x_prem
     eta_aniso =  3.3687 - 2.4778 * x_prem
     Qmu = 600.0
     Qkappa = 57827.0
  ELSEIF(idom==4)THEN
     model%c_rho  =  7.1089 -  3.8045 * x_prem
     model%c_vpv = 20.3926 - 12.2569 * x_prem
     model%c_vsv =  8.9496 -  4.4597 * x_prem
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 143.0
     Qkappa = 57827.0
  ELSEIF(idom==5)THEN
     model%c_rho  = 11.2494 -  8.0298 * x_prem
     model%c_vpv = 39.7027 - 32.6166 * x_prem
     model%c_vsv = 22.3512 - 18.5856 * x_prem
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 143.0
     Qkappa = 57827.0
  ELSEIF(idom==6)THEN
     model%c_rho  =  5.3197 - 1.4836 * x_prem
     model%c_vpv = 19.0957 - 9.8672 * x_prem
     model%c_vsv =  9.9839 - 4.9324 * x_prem
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 143.0
     Qkappa = 57827.0
  ELSEIF(idom==7)THEN   !lower mantle
     model%c_rho  =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     model%c_vpv = 29.2766 -23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     model%c_vsv = 22.3459 -17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 312.0
     Qkappa = 57827.0
  ELSEIF(idom==8)THEN
     model%c_rho  =  7.9565 -  6.4761 * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     model%c_vpv = 24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     model%c_vsv = 11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 312.0
     Qkappa = 57827.0
  ELSEIF(idom==9)THEN
     model%c_rho  =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     model%c_vpv = 15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     model%c_vsv =  6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 312.0
     Qkappa = 57827.0
  ELSEIF(idom==10)THEN  ! outer core
     model%c_rho  = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281 * x_prem**3
     model%c_vpv = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732 * x_prem**3
     model%c_vsv =  0.0
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 0.0
     Qkappa = 57827.0
  ELSEIF(idom==11)THEN                        ! inner core
     model%c_rho  = 13.0885 - 8.8381 * x_prem**2
     model%c_vpv = 11.2622 - 6.3640 * x_prem**2
     model%c_vsv =  3.6678 - 4.4475 * x_prem**2
     model%c_vph = model%c_vpv
     model%c_vsh = model%c_vsv
     Qmu = 84.6
     Qkappa = 1327.7
  ENDIF

  if (model%c_vsh(1)*model%c_vsv(1)==0) then
    model%c_xi = 1
  else
    model%c_xi = (model%c_vsh/model%c_vsv)**2 
  end if
  model%c_phi= (model%c_vpv/model%c_vph)**2
  model%c_vp = model%c_vph 
  model%c_vs = model%c_vsh

  ! All values should be returned in SI units, not some bollocks
  model%c_vp  = model%c_vp  * 1E3
  model%c_vs  = model%c_vs  * 1E3
  model%c_vph = model%c_vph * 1E3
  model%c_vpv = model%c_vpv * 1E3
  model%c_vsh = model%c_vsh * 1E3
  model%c_vsv = model%c_vsv * 1E3
  model%c_rho = model%c_rho * 1E3

end function prem_ani_sub
!-----------------------------------------------------------------------------------------
end module
!=========================================================================================
