!=========================================================================================
module source_class

    use global_parameters,               only : sp, dp, qp, pi, deg2rad, rad2deg, verbose, &
                                                lu_out
    use commpi,                          only : pabort
    implicit none

    private
    public   :: src_param_type

    type src_param_type
        real(kind=dp)               :: mij(6)              ! Mrr Mtt Mpp Mrt Mrp Mtp
        real(kind=dp)               :: mij_voigt(6)        ! Mtt Mpp Mrr Mrp Mrt Mtp
        real(kind=dp)               :: colat, lat, lon     ! in radians
        real(kind=dp)               :: colatd, latd, lond  ! in degrees
        real(kind=dp), dimension(3) :: r                   ! cartesian coordinates in m
        real(kind=dp)               :: depth, radius       ! in km
        real(kind=dp)               :: shift_time          ! in seconds
        real(kind=dp), allocatable  :: shift_time_sample   ! in samples (based on sampling
                                                           ! rate of the netcdf file)
        logical                              :: have_stf = .false.
        real(kind=dp), allocatable           :: stf(:), stf_resampled(:)
        real(kind=dp)                        :: stf_dt, stf_dt_resampled
        complex(kind=dp), allocatable        :: stf_fd(:), stf_reconv_fd(:)

        real(kind=dp), dimension(3,3)        :: rot_mat, trans_rot_mat
        contains
           procedure, pass                   :: init
           procedure, pass                   :: read_cmtsolution
           procedure, pass                   :: def_rot_matrix
           procedure, pass                   :: set_shift_time_sample
           procedure, pass                   :: read_stf
    end type
contains

!-----------------------------------------------------------------------------------------
!> This routine initializes the source object
subroutine init(this, lat, lon, mij, depth)
   class(src_param_type)      :: this

   real(kind=dp), intent(in)  :: lat, lon, mij(6), depth ! MIJ in Nm here

   this%latd   = lat
   this%lond   = lon
   this%colatd = 90 - this%latd

   this%colat  = this%colatd * deg2rad
   this%lon    = this%lond   * deg2rad
   this%lat    = this%latd   * deg2rad
   this%depth  = depth

   this%shift_time = 0

   !TODO hardcoded earth radius for now until I know where to get earth's radius from (MvD)
   this%radius = 6371 - depth

   this%r(1) = cos(this%lat) * cos(this%lon) * this%radius * 1d3
   this%r(2) = cos(this%lat) * sin(this%lon) * this%radius * 1d3
   this%r(3) = sin(this%lat)                 * this%radius * 1d3

   this%mij    = mij

   ! CMTSOLUTION : Mrr Mtt Mpp Mrt Mrp Mtp
   ! voigt in tpr: Mtt Mpp Mrr Mrp Mrt Mtp
   this%mij_voigt(1) = this%mij(2)
   this%mij_voigt(2) = this%mij(3)
   this%mij_voigt(3) = this%mij(1)
   this%mij_voigt(4) = this%mij(5)
   this%mij_voigt(5) = this%mij(4)
   this%mij_voigt(6) = this%mij(6)

   call this%def_rot_matrix()
  
   write(lu_out, '(" Source:    ")') 
   write(lu_out, '("  Depth:     ", F8.3)')   depth
   write(lu_out, '("  Latitude:  ", F8.3)')   this%latd
   write(lu_out, '("  Longitude: ", F8.3)')   this%lond
   write(lu_out, '("  M_rr:      ", E15.8)')  this%mij(1)
   write(lu_out, '("  M_tt:      ", E15.8)')  this%mij(2)
   write(lu_out, '("  M_pp:      ", E15.8)')  this%mij(3)
   write(lu_out, '("  M_rt:      ", E15.8)')  this%mij(4)
   write(lu_out, '("  M_rp:      ", E15.8)')  this%mij(5)
   write(lu_out, '("  M_tp:      ", E15.8)')  this%mij(6)

   this%have_stf = .false.

end subroutine init
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> This routine initializes the source object
subroutine read_cmtsolution(this, fname)
   class(src_param_type)      :: this

   character(len=*), intent(in), optional :: fname
   integer                    :: lu_cmtsolution, ioerr
   character(len=256)         :: cmtsolution_file, junk
   real(kind=dp)              :: mij_dyncm(6)

   if (present(fname)) then
      cmtsolution_file = trim(fname)
   else
      cmtsolution_file = 'CMTSOLUTION'
   endif

   open(newunit=lu_cmtsolution, file=trim(cmtsolution_file), status='old', &
        action='read', iostat=ioerr)

   if (ioerr /= 0) then
      print *, 'ERROR: Check input file ''', trim(cmtsolution_file), '''! Is it still there?' 
      call pabort
   end if

   read(lu_cmtsolution,*) junk ! first crap line
   read(lu_cmtsolution,*) junk ! event name
   read(lu_cmtsolution,*) junk, junk, this%shift_time
   read(lu_cmtsolution,*) junk ! half duration
   read(lu_cmtsolution,*) junk, this%latd
   read(lu_cmtsolution,*) junk, this%lond
   read(lu_cmtsolution,*) junk, this%depth

   this%colatd = 90 - this%latd

   this%colat  = this%colatd * deg2rad
   this%lon    = this%lond   * deg2rad
   this%lat    = this%latd   * deg2rad

   !TODO hardcoded earth radius for now until I know where to get earth's radius from (MvD)
   this%radius = 6371 - this%depth

   this%r(1) = cos(this%lat) * cos(this%lon) * this%radius * 1d3
   this%r(2) = cos(this%lat) * sin(this%lon) * this%radius * 1d3
   this%r(3) = sin(this%lat)                 * this%radius * 1d3


   read(lu_cmtsolution,*) junk, mij_dyncm(1)
   read(lu_cmtsolution,*) junk, mij_dyncm(2)
   read(lu_cmtsolution,*) junk, mij_dyncm(3)
   read(lu_cmtsolution,*) junk, mij_dyncm(4)
   read(lu_cmtsolution,*) junk, mij_dyncm(5)
   read(lu_cmtsolution,*) junk, mij_dyncm(6)

   this%mij = mij_dyncm / 1e7 ! dyn cm -> Nm

   ! CMTSOLUTION : Mrr Mtt Mpp Mrt Mrp Mtp
   ! voigt in tpr: Mtt Mpp Mrr Mrp Mrt Mtp
   this%mij_voigt(1) = this%mij(2)
   this%mij_voigt(2) = this%mij(3)
   this%mij_voigt(3) = this%mij(1)
   this%mij_voigt(4) = this%mij(5)
   this%mij_voigt(5) = this%mij(4)
   this%mij_voigt(6) = this%mij(6)

   call this%def_rot_matrix()

   this%have_stf = .false.

   close(lu_cmtsolution)

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!>  This function defines the rotation matrix to rotate coordinates to the 
!!  forward source system. Taken from AxiSEM solver.
subroutine def_rot_matrix(this)
   class(src_param_type)          :: this
   real(kind=dp)                  :: srccolat, srclon
   character(len=256)             :: fmtstring

   srccolat = this%colat
   srclon   = this%lon

   fmtstring = '("  Source colatitude: ", F8.3, "; source longitude: ", F8.3)'
  
   if (verbose > 1) write(lu_out,fmtstring) srccolat/deg2rad, srclon/deg2rad
   

   ! This is the rotation matrix of Nissen-Meyer, Dahlen, Fournier, GJI 2007.
   this%rot_mat(1,1) =  dcos(srccolat) * dcos(srclon)
   this%rot_mat(2,1) =  dcos(srccolat) * dsin(srclon)
   this%rot_mat(3,1) = -dsin(srccolat)
   this%rot_mat(1,2) = -dsin(srclon)
   this%rot_mat(2,2) =  dcos(srclon)
   this%rot_mat(3,2) =  0.d0
   this%rot_mat(1,3) =  dsin(srccolat) * dcos(srclon)
   this%rot_mat(2,3) =  dsin(srccolat) * dsin(srclon)
   this%rot_mat(3,3) =  dcos(srccolat)

   this%trans_rot_mat = transpose(this%rot_mat)

   if (verbose > 1) then
      write(lu_out,*)             '  Rotation matrix:'
      write(lu_out,'(3("  ", 3(ES11.3)/))') this%rot_mat
      write(lu_out,*)             '  Rotation matrix, transposed:'
      write(lu_out,'(3("  ", 3(ES11.3)/))') this%trans_rot_mat
   end if

end subroutine def_rot_matrix
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine set_shift_time_sample(this, dt)
   class(src_param_type)          :: this
   real(kind=dp), intent(in)      :: dt

   if (.not. allocated(this%shift_time_sample)) &
      allocate(this%shift_time_sample)
   this%shift_time_sample = this%shift_time / dt

end subroutine 
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_stf(this, filename) !, dt)
   use lanczos,             only   : lanczos_resample
   class(src_param_type)          :: this
   character(len=*)               :: filename  ! file from which to read the STF
   !real(kind=dp), intent(in)      :: dt        ! dt to which to resample the STF

   !real(kind=dp), allocatable     :: time_orig(:), time_new(:)
   !real(kind=dp)                  :: dt_orig
   integer                        :: nsamp_orig
   integer                        :: lu_stf, isamp, ioerr

   ! Read STF from input file 'filename'
   open(newunit=lu_stf, file=trim(filename), status='old', &
        action='read', iostat=ioerr)

   if (ioerr /= 0) then
      print *, 'ERROR: Check STF input file ''', trim(filename), '''! Is it still there?' 
      call pabort
   end if

   read(lu_stf,*) nsamp_orig, this%stf_dt
   allocate(this%stf(nsamp_orig))
   do isamp = 1, nsamp_orig
     read(lu_stf,*) this%stf(isamp)
   end do
   close(lu_stf)

   !this%stf_dt_resampled = dt

   !this%stf_resampled = lanczos_resample(this%stf, this%stf_dt, dt, a=8)

   this%have_stf = .true.

end subroutine read_stf
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
