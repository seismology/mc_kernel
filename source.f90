!=========================================================================================
module source_class

    use global_parameters,               only : sp, dp, pi, deg2rad, rad2deg, verbose, &
                                                lu_out
    use commpi,                          only : pabort
    implicit none

    private
    public   :: src_param_type
    public   :: fft_stf

    type src_param_type
        real(kind=dp)               :: mij(6)              ! Mrr Mtt Mpp Mrt Mrp Mtp
        real(kind=dp)               :: mij_voigt(6)        ! Mtt Mpp Mrr Mrp Mrt Mtp
        real(kind=dp)               :: colat, lat, lon     ! in radians
        real(kind=dp)               :: colatd, latd, lond  ! in degrees
        real(kind=dp)               :: x, y, z             ! cartesian coordinates in km
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
           procedure, pass                   :: resample_stf
    end type
contains

!-----------------------------------------------------------------------------------------
!> This routine initializes the source object
subroutine init(this, lat, lon, mij, depth)
   class(src_param_type)      :: this

   real(kind=dp), intent(in)  :: lat, lon, mij(6), depth

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

   this%x = dcos(this%lat) * dcos(this%lon) * this%radius
   this%y = dcos(this%lat) * dsin(this%lon) * this%radius
   this%z = dsin(this%lat) * this%radius

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

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> This routine initializes the source object
subroutine read_cmtsolution(this, fname)
   class(src_param_type)      :: this

   character(len=*), intent(in), optional :: fname
   integer                    :: lu_cmtsolution, ioerr
   character(len=256)         :: cmtsolution_file, junk

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

   this%x = dcos(this%lat) * dcos(this%lon) * this%radius
   this%y = dcos(this%lat) * dsin(this%lon) * this%radius
   this%z = dsin(this%lat) * this%radius

   read(lu_cmtsolution,*) junk, this%mij(1)
   read(lu_cmtsolution,*) junk, this%mij(2)
   read(lu_cmtsolution,*) junk, this%mij(3)
   read(lu_cmtsolution,*) junk, this%mij(4)
   read(lu_cmtsolution,*) junk, this%mij(5)
   read(lu_cmtsolution,*) junk, this%mij(6)

   this%mij = this%mij / 1e7 ! dyn cm -> Nm

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
subroutine resample_stf(this, dt, nsamp)
   class(src_param_type)          :: this
   real(kind=dp), intent(in)      :: dt
   integer, intent(in)            :: nsamp

   real(kind=dp), allocatable     :: time_orig(:), time_new(:)
   real(kind=dp)                  :: dt_orig
   integer                        :: nsamp_orig
   integer                        :: i, j

   if (this%have_stf) then

      nsamp_orig = size(this%stf)
      dt_orig = this%stf_dt

      allocate(time_orig(nsamp_orig))
      allocate(time_new(nsamp))

      if (allocated(this%stf_resampled)) deallocate(this%stf_resampled)
      allocate(this%stf_resampled(nsamp))

      this%stf_dt_resampled = dt

      do i = 1, nsamp_orig
         time_orig(i) = (i-1) * dt_orig
      enddo

      do i = 1, nsamp
         time_new(i) = (i-1) * dt
      enddo

      this%stf_resampled(:) = 0

      outer: do i = 1, nsamp
         inner: do j = 1, nsamp_orig
            if (time_new(i) <= time_orig(j)) then
               if (j < 2) then
                  this%stf_resampled(i) = 0
               else
                  this%stf_resampled(i) = &
                          this%stf(j-1) * (time_orig(j) - time_new(i))      / dt_orig &
                        + this%stf(j)   * (time_new(i)  - time_orig(j - 1)) / dt_orig 
               endif

               cycle outer
            endif
         enddo inner
      enddo outer

   else
      ! no stf provided, so we assume a dirac
      if (allocated(this%stf_resampled)) deallocate(this%stf_resampled)
      allocate(this%stf_resampled(nsamp))

      this%stf_resampled = 0
      this%stf_resampled(1) = 1d0 / dt
      this%stf_dt_resampled = dt
   endif

end subroutine resample_stf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fft_stf(sources, stf_fwd)

   use fft, only : rfft_type
   type(src_param_type), intent(inout)  :: sources(:)
   real(kind=dp), optional              :: stf_fwd(:)

   type(rfft_type)                  :: fftt
   integer                          :: nsources, nsamp
   real(kind=dp)                    :: dt
   integer                          :: isource
   real(kind=dp), allocatable       :: stf_buf(:,:)
   complex(kind=dp), allocatable    :: stf_buf_fd(:,:)
   complex(kind=dp), allocatable    :: stf_fwd_fd(:,:)

   nsources = size(sources)

   do isource=1, nsources
      if (.not. allocated(sources(isource)%stf_resampled)) then
         write(6,*) 'ERROR: stf needs to be resampled before fft'
         call pabort
      endif
   enddo

   nsamp = size(sources(1)%stf_resampled)
   dt = sources(1)%stf_dt_resampled

   do isource=2, nsources
      if (size(sources(isource)%stf_resampled) /= nsamp) then
         write(6,*) 'ERROR: stfs have incompatible number of samples'
         call pabort
      endif

      if (sources(isource)%stf_dt_resampled /= dt) then
         write(6,*) 'ERROR: stfs have incompatible dt'
         call pabort
      endif
   enddo

   if (present(stf_fwd)) then
      if (size(stf_fwd) /= nsamp) then
         write(6,*) 'ERROR: stf_fwd has incompatible number of samples'
         call pabort
      endif
   endif
    
   call fftt%init(ntimes_in=nsamp*2, ndim=1, ntraces=1, dt=dt, nfft=nsamp*2)

   allocate(stf_buf(nsamp*2,1))
   allocate(stf_buf_fd(fftt%get_nomega(),1))

   if (present(stf_fwd)) then
       allocate(stf_fwd_fd(fftt%get_nomega(),1))
       stf_buf = 0
       stf_buf(1:nsamp,1) = stf_fwd

       call fftt%rfft(stf_buf, stf_fwd_fd)
   endif

   do isource=1, nsources
      if (allocated(sources(isource)%stf_fd)) deallocate(sources(isource)%stf_fd)
      allocate(sources(isource)%stf_fd(fftt%get_nomega()))

      stf_buf = 0
      stf_buf(1:nsamp,1) = sources(isource)%stf_resampled

      call fftt%rfft(stf_buf, stf_buf_fd)

      sources(isource)%stf_fd = stf_buf_fd(:,1)

      if (present(stf_fwd)) then
         if (allocated(sources(isource)%stf_reconv_fd)) &
             deallocate(sources(isource)%stf_reconv_fd)
         allocate(sources(isource)%stf_reconv_fd(fftt%get_nomega()))

         sources(isource)%stf_reconv_fd = sources(isource)%stf_fd / stf_fwd_fd(:,1)
      endif
   enddo

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
