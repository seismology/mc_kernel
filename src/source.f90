!******************************************************************************
!
!    This file is part of:
!    MC Kernel: Calculating seismic sensitivity kernels on unstructured meshes
!    Copyright (C) 2016 Simon Staehler, Martin van Driel, Ludwig Auer
!
!    You can find the latest version of the software at:
!    <https://www.github.com/tomography/mckernel>
!
!    MC Kernel is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MC Kernel is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MC Kernel. If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

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
        real(kind=dp)               :: depth               ! in km
        real(kind=dp), allocatable  :: stf_shift(:)        ! in seconds
        real(kind=dp)               :: shift_time          ! in seconds
        real(kind=dp), allocatable  :: shift_time_sample   ! in samples (based on sampling
                                                           ! rate of the netcdf file)
        integer                     :: nstf                ! number of stfs
        logical                         :: have_stf = .false.
        real(kind=dp), allocatable      :: stf(:,:), stf_resampled(:,:)
        real(kind=dp), allocatable      :: stf_dt(:)
        real(kind=dp)                   :: stf_dt_resampled
        complex(kind=dp), allocatable   :: stf_fd(:,:), stf_reconv_fd(:,:)
        character(len=32), allocatable  :: src_name(:)

        real(kind=dp), dimension(3,3)   :: rot_mat, trans_rot_mat
        contains
           procedure, pass              :: init
           procedure, pass              :: get_r
           procedure, pass              :: def_rot_matrix
           procedure, pass              :: set_shift_time_sample
           procedure, pass              :: read_stf
           procedure, pass              :: get_src_index
    end type
contains

!-----------------------------------------------------------------------------------------
!> This routine initializes the source object
subroutine init(this, lat, lon, mij, depth, planet_radius)
   class(src_param_type)      :: this

   real(kind=dp), intent(in)  :: lat, lon, mij(6), depth ! MIJ in Nm here
   real(kind=dp), intent(in), optional :: planet_radius  ! Default value is 6371 (Earth) 

   this%latd   = lat
   this%lond   = lon
   this%colatd = 90 - this%latd

   this%colat  = this%colatd * deg2rad
   this%lon    = this%lond   * deg2rad
   this%lat    = this%latd   * deg2rad
   this%depth  = depth

   this%shift_time = 0

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
function get_r(this, planet_radius) result(r)
   class(src_param_type)                :: this
   real(kind=dp), intent(in), optional  :: planet_radius
   real(kind=dp)                        :: radius
   real(kind=dp)                        :: r(3)

   if (present(planet_radius)) then
     radius = (planet_radius - this%depth * 1d3) 
   else 
     radius = (6371d3 - this%depth * 1d3)
   end if

   r(1) = cos(this%lat) * cos(this%lon) * radius
   r(2) = cos(this%lat) * sin(this%lon) * radius
   r(3) = sin(this%lat)                 * radius

end function get_r
!----------------------------------------------------------------------------------------

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
!> Read STF from input file 'filename'
subroutine read_stf(this, filename)
   use lanczos,             only   : lanczos_resample
   class(src_param_type)          :: this
   character(len=*)               :: filename  ! file from which to read the STF
   integer, allocatable           :: nsamp_orig(:)
   integer                        :: lu_stf, isamp, isrc, ioerr

   open(newunit=lu_stf, file=trim(filename), status='old', &
        action='read', iostat=ioerr)

   if (ioerr /= 0) then
      print *, 'ERROR: Check STF input file ''', trim(filename), '''! Is it still there?' 
      call pabort
   end if
   
   read(lu_stf, *) this%nstf
   allocate(this%src_name(this%nstf))
   allocate(nsamp_orig(this%nstf))
   allocate(this%stf_dt(this%nstf))
   allocate(this%stf_shift(this%nstf))

   do isrc = 1, this%nstf
     read(lu_stf, *) this%src_name(isrc), nsamp_orig(isrc), &
                     this%stf_dt(isrc), this%stf_shift(isrc)
     do isamp = 1, nsamp_orig(isrc)
       read(lu_stf,*) 
     end do
   end do

   allocate(this%stf(maxval(nsamp_orig), this%nstf))

   rewind(lu_stf)
   read(lu_stf, *) ! this%nstf

   do isrc = 1, this%nstf
     read(lu_stf, *) ! src_name, nsamp, dt, shift
     do isamp = 1, nsamp_orig(isrc)
       read(lu_stf,*) this%stf(isamp, isrc)
     end do
   end do
   close(lu_stf)

   this%have_stf = .true.

end subroutine read_stf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function get_src_index(this, src_name_in)
   class(src_param_type)          :: this
   character(len=*)               :: src_name_in
   integer                        :: isrc

   do isrc = 1, this%nstf
     if (trim(src_name_in) == trim(this%src_name(isrc))) exit
   end do

   if (isrc>this%nstf) then
     print *, 'ERROR: Source ', trim(src_name_in), 'from receiver input file '
     print *, '       not defined in STF input file'
     print *, 'Available STFs: ', [(this%src_name(isrc), isrc= 1, this%nstf)]
     call pabort(do_traceback=.false.)
   end if

   get_src_index = isrc

end function get_src_index
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
