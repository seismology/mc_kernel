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
   use iso_c_binding, only: c_new_line
   integer  :: npoints, chunksize, ipoint, start_chunk, count_chunk, i, iinchunk
   integer  :: start_chunk_ref, count_chunk_ref
   character(len=512) :: err_msg
   character(len=1) :: cnl = c_new_line

   ! This should have 20 chunks
   ! The first 19 are 5 points long, the last one has 4.
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

   ! Second test. start_chunk should never be smaller than one and 
   ! start_chunk+count_chunk-1 should be smaller than npoints

   do chunksize = 3, 11
     do npoints = 100, 110
       do ipoint = 1, npoints
         call get_chunk_bounds(pointid     = ipoint,        &
                               chunksize   = chunksize,     &
                               npoints     = npoints,       & 
                               start_chunk = start_chunk,   &
                               count_chunk = count_chunk)

         write(err_msg, 100) cnl, start_chunk, cnl, count_chunk, cnl, ipoint, cnl, npoints, cnl, chunksize
         call assert_true(start_chunk >= 0, trim(err_msg))

         write(err_msg, 101) cnl, start_chunk, cnl, count_chunk, cnl, ipoint, cnl, npoints, cnl, chunksize 
         call assert_true(start_chunk + count_chunk - 1 <= npoints, trim(err_msg))
       end do
     end do
   end do

100      format('Start of chunk larger zero:', A, &        
                'start_chunk: ', I4, A, &
                'count_chunk: ', I4, A, &
                'ipoint:      ', I4, A, &
                'npoints:     ', I4, A, &
                'chunksize:   ', I4)
101      format('Start + count of chunk smaller npoints:', A, &        
                'start_chunk: ', I4, A, &
                'count_chunk: ', I4, A, &
                'ipoint:      ', I4, A, &
                'npoints:     ', I4, A, &
                'chunksize:   ', I4)
end subroutine test_get_chunk_bounds
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_set_params()
   use type_parameter, only : parameter_type
   type(semdata_type)      :: sem_data
   type(parameter_type)    :: parameters
   character(len=512)      :: fwd_dir, bwd_dir

   call parameters%read_parameters('./inparam_test')
   call parameters%read_source()

   call sem_data%set_params(fwd_dir              = parameters%fwd_dir,          &
                            bwd_dir              = parameters%bwd_dir,          &
                            strain_buffer_size   = 100,                         &
                            displ_buffer_size    = 100,                         &
                            strain_type          = 'straintensor_trace',        &
                            desired_source_depth = parameters%source%depth)

   call assert_equal_int(sem_data%nsim_fwd, 4, 'nsim_fwd == 4')
   call assert_equal_int(sem_data%nsim_bwd, 2, 'nsim_bwd == 2')

end subroutine test_readfields_set_params
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_open_files()
   use type_parameter, only : parameter_type
   type(parameter_type)    :: parameters
   type(semdata_type)      :: sem_data

   call parameters%read_parameters('./inparam_test')
   call parameters%read_source()
   
   call sem_data%set_params(fwd_dir              = parameters%fwd_dir,          &
                            bwd_dir              = parameters%bwd_dir,          &
                            strain_buffer_size   = 100,                         &
                            displ_buffer_size    = 100,                         &
                            strain_type          = 'straintensor_trace',        &
                            desired_source_depth = parameters%source%depth)

   call sem_data%open_files()

   call sem_data%close_files()

end subroutine  test_readfields_open_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Just show that at least it does not crash
subroutine test_readfields_reopen_files()
   use type_parameter, only : parameter_type
   type(parameter_type)    :: parameters
   type(semdata_type)      :: sem_data

   call parameters%read_parameters('./inparam_test')
   call parameters%read_source()
   
   call sem_data%set_params(fwd_dir              = parameters%fwd_dir,          &
                            bwd_dir              = parameters%bwd_dir,          &
                            strain_buffer_size   = 100,                         &
                            displ_buffer_size    = 100,                         &
                            strain_type          = 'straintensor_trace',        &
                            desired_source_depth = parameters%source%depth)

   call sem_data%open_files()

   call sem_data%reopen_files()

   call sem_data%close_files()

end subroutine  test_readfields_reopen_files
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Reads straintrace from a coordinate and compares with a reference version
!! Reference version was manually extracted from a AxiSEM pure MRR run.
subroutine test_readfields_load_fw_points
   use readfields,     only    : semdata_type
   use type_parameter, only    : parameter_type
   use fft,            only: rfft_type, taperandzeropad
   type(parameter_type)       :: parameters
   type(semdata_type)         :: sem_data
   type(rfft_type)            :: fft_data
   real(kind=dp), allocatable :: u(:,:,:), straintrace_ref(:), straintrace_res(:)
   real(kind=dp), allocatable :: straintrace_res_filt(:,:)
   real(kind=dp)              :: coordinates(3,2), df
   integer                    :: ntimes, nomega, i, lu_refstrain
   
   call parameters%read_parameters('./inparam_load_wavefield')
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

   call parameters%read_filter(nomega=nomega, df=df)
   call parameters%read_kernel(sem_data, parameters%filter)

   ! Read reference strain trace and point coordinates
   open(newunit=lu_refstrain, file='./read_strain_point1.txt', action='read')
   read(lu_refstrain, *) coordinates(:,1)
   allocate(straintrace_ref(sem_data%ndumps))
   do i=1, ntimes
     read(lu_refstrain, *) straintrace_ref(i)
   end do

   coordinates(:,2) = [-1d6, 1d6, 1d6]

   u = sem_data%load_fw_points(coordinates, parameters%source)

   straintrace_res(:) = u(:, 1, 1)
   !straintrace_res_filt = parameters%kernel(1)%filter%apply(u(:,1,:), kind='fwd')

   call assert_comparable(straintrace_ref, straintrace_res, 1d-7, &
                          'strain trace is equal to reference, point one')

   call sem_data%close_files()

end subroutine test_readfields_load_fw_points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Load a seismogram and compare with instaseis, Z-component
subroutine test_load_seismograms_rdbm_Z
   use readfields,     only : semdata_type
   use type_parameter, only : parameter_type
   use fft,            only: rfft_type, taperandzeropad
   type(parameter_type)    :: parameters
   type(semdata_type)      :: sem_data
   type(rfft_type)         :: fft_data
   integer                 :: nomega, ntimes, ntimes_reference, ikernel
   real(kind=dp)           :: df, t
   real(kind=dp), allocatable :: seis(:), seis_ref(:)
   character(len=128)       :: kernel_name
   
   call parameters%read_parameters('./inparam_load_seismogram_Z')
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

   ! Read Kernels. (Filtered) seismograms are written out here.
   call parameters%read_kernel(sem_data, parameters%filter)

   call sem_data%close_files()

   do ikernel = 1, parameters%nkernel
     ! Compare seismograms 
     kernel_name = parameters%kernel(ikernel)%name
     call compare_seismograms(stat_name = trim(kernel_name), &
                              message   = 'Seismograms '//trim(kernel_name)//', misfit: ') 
   end do

end subroutine test_load_seismograms_rdbm_Z
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Load a seismogram and compare with instaseis, R-component
subroutine test_load_seismograms_rdbm_R
   use readfields,     only : semdata_type
   use type_parameter, only : parameter_type
   use fft,            only: rfft_type, taperandzeropad
   type(parameter_type)    :: parameters
   type(semdata_type)      :: sem_data
   type(rfft_type)         :: fft_data
   integer                 :: nomega, ntimes, ntimes_reference, ikernel
   real(kind=dp)           :: df, t
   real(kind=dp), allocatable :: seis(:), seis_ref(:)
   character(len=128)       :: kernel_name
   
   call parameters%read_parameters('./inparam_load_seismogram_R')
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

   ! Read Kernels. (Filtered) seismograms are written out here.
   call parameters%read_kernel(sem_data, parameters%filter)

   call sem_data%close_files()

   do ikernel = 1, parameters%nkernel
     ! Compare seismograms 
     kernel_name = parameters%kernel(ikernel)%name
     call compare_seismograms(stat_name = trim(kernel_name), &
                              message   = 'Seismograms '//trim(kernel_name)//', misfit: ') 
   end do

end subroutine test_load_seismograms_rdbm_R
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Load a seismogram and compare with instaseis, R-component
subroutine test_load_seismograms_rdbm_T
   use readfields,     only : semdata_type
   use type_parameter, only : parameter_type
   use fft,            only: rfft_type, taperandzeropad
   type(parameter_type)    :: parameters
   type(semdata_type)      :: sem_data
   type(rfft_type)         :: fft_data
   integer                 :: nomega, ntimes, ntimes_reference, ikernel
   real(kind=dp)           :: df, t
   real(kind=dp), allocatable :: seis(:), seis_ref(:)
   character(len=128)       :: kernel_name
   
   call parameters%read_parameters('./inparam_load_seismogram_T')
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

   ! Read Kernels. (Filtered) seismograms are written out here.
   call parameters%read_kernel(sem_data, parameters%filter)

   call sem_data%close_files()

   do ikernel = 1, parameters%nkernel
     ! Compare seismograms 
     kernel_name = parameters%kernel(ikernel)%name
     call compare_seismograms(stat_name = trim(kernel_name), &
                              message   = 'Seismograms '//trim(kernel_name)//', misfit: ') 
   end do

end subroutine test_load_seismograms_rdbm_T
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compare_seismograms(stat_name, message)
  character(len=*), intent(in)  :: stat_name
  character(len=*), intent(in)  :: message
  character(len=256)            :: filename, message_full
  real(kind=dp), allocatable    :: seis_raw(:), seis_ref_raw(:)
  real(kind=dp), allocatable    :: seis_disp(:), seis_ref_disp(:)
  real(kind=dp), allocatable    :: seis_velo(:), seis_ref_velo(:)
  integer                       :: ntimes_reference, isample, lu_seis
  real(kind=dp)                 :: temp, misfit_raw, misfit_velo, misfit_disp

  ! Compare RAW seismogram                     
  ! Load reference seismogram
  filename = './seism_ref_raw_'//trim(stat_name)
  open(newunit=lu_seis, file=filename, action='read', status='old')
  read(lu_seis,*) ntimes_reference
  allocate(seis_ref_raw(ntimes_reference))
  do isample = 1, ntimes_reference
    read(lu_seis,*) seis_ref_raw(isample)
  end do
  close(lu_seis)

  ! Retrieve seismograms
  ! 1st one filtered with Butterworth, 6th order at 40s
  filename = './Seismograms/seism_raw_'//trim(stat_name)
  open(newunit=lu_seis, file=filename, action='read', status='old')
  allocate(seis_raw(ntimes_reference))
  do isample = 1, ntimes_reference
    read(lu_seis,*) temp, seis_raw(isample)
  end do
  close(lu_seis)

  ! Compare FILTERED seismograms                     
  ! Load reference seismograms
  filename = './seism_ref_'//trim(stat_name)
  open(newunit=lu_seis, file=filename, action='read', status='old')
  read(lu_seis,*) ntimes_reference
  allocate(seis_ref_disp(ntimes_reference))
  allocate(seis_ref_velo(ntimes_reference))
  do isample = 1, ntimes_reference
    read(lu_seis,*) seis_ref_disp(isample), seis_ref_velo(isample)
  end do
  close(lu_seis)

  ! Retrieve seismograms
  ! 1st one filtered with Butterworth, 6th order at 40s
  filename = './Seismograms/seism_'//trim(stat_name)
  open(newunit=lu_seis, file=filename, action='read', status='old')
  allocate(seis_disp(ntimes_reference))
  allocate(seis_velo(ntimes_reference))
  do isample = 1, ntimes_reference
    read(lu_seis,*) temp, seis_disp(isample), seis_velo(isample)
  end do
  close(lu_seis)

  misfit_raw = norm2(seis_raw - seis_ref_raw) / norm2(seis_ref_raw)
  misfit_disp = norm2(seis_disp - seis_ref_disp) / norm2(seis_ref_disp)
  misfit_velo = norm2(seis_velo - seis_ref_velo) / norm2(seis_ref_velo)

  write(message_full, '(A, " (raw):  ", E15.8)') message, misfit_raw
  call assert_true(misfit_raw < 1d-4, message_full)
  ! The limits are very high here, I know. Seems to be some inconsistency
  ! in definition of filters compared to instaseis.
  write(message_full, '(A, " (disp): ",  E15.8)') message, misfit_disp
  call assert_true(misfit_disp < 1d-0, message_full)
  write(message_full, '(A, " (velo): ",  E15.8)') message, misfit_velo
  call assert_true(misfit_velo < 1d-0, message_full)

end subroutine compare_seismograms
!-----------------------------------------------------------------------------------------
!
!!   ! Compare SEM Source Time Function
!!   ! Load reference STF
!!   open(newunit=lu_seis, file='./SEM_stf_ref.dat', action='read')
!!   read(lu_seis,*) ntimes_reference
!!   allocate(seis_ref(ntimes_reference))
!!   do isample = 1, ntimes_reference
!!     read(lu_seis,*) seis_ref(isample)
!!   end do
!!   close(lu_seis)
!!
!!   ! Retrieve seismograms
!!   ! 1st one filtered with Butterworth, 6th order at 40s
!!   open(newunit=lu_seis, file='Filters/stf_out_Butterw_LP_O2_60.000_.000')
!!   allocate(seis(ntimes_reference))
!!   do isample = 1, ntimes_reference
!!     read(lu_seis,*) t, seis(isample)
!!   end do
!!   close(lu_seis)
!!
!!   open(newunit=lu_seis, file='./output/sem_stf_comparison.txt', action='write')
!!   do isample = 1, ntimes_reference
!!     write(lu_seis, *) seis(isample), seis_ref(isample)
!!   end do
!!   close(lu_seis)
!!   ! The threshold is very high.
!!   call assert_comparable(seis, seis_ref, 5d-4, &  
!!                          'SEM STF is comparable to reference')
!!   deallocate(seis_ref, seis) 
!!
!
!end subroutine test_load_seismograms_rdbm_T
!!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_readfields_load_model_coeffs
   use type_parameter, only    : parameter_type
   use backgroundmodel, only   : backgroundmodel_type
   use global_parameters, only : pi, sp, dp
   type(parameter_type)       :: parameters
   type(semdata_type)         :: sem_data
   type(backgroundmodel_type) :: model, model_ref
   integer, parameter         :: npoints = 1, nradius = 11
   real(kind=dp)              :: coordinates(3,npoints), phi(npoints), theta(npoints), r(nradius)
   integer                    :: ipoint, iradius

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


   call parameters%read_parameters('./inparam_test')
   call parameters%read_source()
   call parameters%read_receiver()

   call sem_data%set_params(fwd_dir              = parameters%fwd_dir,          &
                            bwd_dir              = parameters%bwd_dir,          &
                            strain_buffer_size   = 100,                         &
                            displ_buffer_size    = 100,                         &
                            strain_type          = 'straintensor_trace',        &
                            desired_source_depth = parameters%source%depth)

   call sem_data%open_files()

   call sem_data%read_meshes()
   
   do iradius = 1, nradius

     call random_number(phi)
     call random_number(theta)
     phi = phi * 2.d0 * pi
     theta = (theta - 0.5) * pi

     do ipoint = 1, npoints
        coordinates(:, ipoint) = [cos(phi(ipoint)) * sin(theta(ipoint)), &
                                  sin(phi(ipoint)) * sin(theta(ipoint)), &
                                  cos(theta(ipoint))] * r(iradius)
     end do
     call flush(6)

     model = sem_data%load_model_coeffs(coordinates)
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
subroutine test_dampen_field
  integer, parameter :: npoints = 1000
  real(kind=dp)      :: field(npoints,2,6), field_ref(npoints,2,6)
  real(kind=dp)      :: r_points(3,npoints), r_src(3), r_max, pre_fac
  integer            :: ipoint

  call random_number(field)
  field_ref = field

  call random_number(r_points)

  r_src = [0.5d0, 0.5d0, 0.5d0]

  r_max = 0.2d0

  call dampen_field(field, r_points, r_src, r_max)

  do ipoint = 1, npoints
    if (norm2(r_src-r_points(:, ipoint))>r_max) then
      call assert_comparable(field(ipoint,:,:), field_ref(ipoint,:,:), &
                             1d-5, 'Value outside r_max equal')
    else
      pre_fac = norm2(r_src-r_points(:, ipoint))/r_max
      call assert_comparable(field(ipoint,:,:), field_ref(ipoint,:,:)*pre_fac, &
                             1d-5, 'Value inside r_max corrected')
    end if
  end do

end subroutine test_dampen_field
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
