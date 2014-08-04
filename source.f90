!=========================================================================================
module source_class

    use global_parameters,               only : sp, dp, pi, deg2rad, verbose, lu_out
    use commpi,                          only : pabort
    implicit none

    type src_param_type
        real(kind=dp)                        :: mij(6)              ! Mrr Mtt Mpp Mrt Mrp Mtp
        real(kind=dp)                        :: mij_voigt(6)        ! Mtt Mpp Mrr Mrp Mrt Mtp
        real(kind=dp)                        :: colat, lat, lon     ! in radians
        real(kind=dp)                        :: colatd, latd, lond  ! in degrees
        real(kind=dp)                        :: x, y, z             ! cartesian 
                                                                    ! coordinates in km
        real(kind=dp)                        :: depth, radius       ! in km
        real(kind=dp), dimension(3,3)        :: rot_mat, trans_rot_mat
        contains
           procedure, pass                   :: init
           procedure, pass                   :: read_cmtsolution
           procedure, pass                   :: def_rot_matrix
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

   !TODO hardcoded radius for now until I know where to get earth's radius from (MvD)
   this%radius = 6371 - depth

   this%x = dcos(this%lat) * dcos(this%lon) * this%radius
   this%y = dcos(this%lat) * dsin(this%lon) * this%radius
   this%z = dsin(this%lat) * this%radius

   this%mij    = mij
   call this%def_rot_matrix()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> This routine initializes the source object
subroutine read_cmtsolution(this)
   class(src_param_type)      :: this

   integer                    :: lu_cmtsolution, ioerr
   character(len=256)         :: cmtsolution_file, junk

   !@TODO: hardcoded for now, should be intent(in) at some point
   cmtsolution_file = 'CMTSOLUTION'

   open(newunit=lu_cmtsolution, file=trim(cmtsolution_file), status='old', &
        action='read', iostat=ioerr)

   if (ioerr /= 0) then
      print *, 'ERROR: Check input file ''', trim(cmtsolution_file), '''! Is it still there?' 
      call pabort
   end if

   read(lu_cmtsolution,*) junk ! first crap line
   read(lu_cmtsolution,*) junk ! event name
   read(lu_cmtsolution,*) junk ! time shift
   read(lu_cmtsolution,*) junk ! half duration
   read(lu_cmtsolution,*) junk, this%latd
   read(lu_cmtsolution,*) junk, this%lond
   read(lu_cmtsolution,*) junk, this%depth

   this%colatd = 90 - this%latd

   this%colat  = this%colatd * deg2rad
   this%lon    = this%lond   * deg2rad
   this%lat    = this%latd   * deg2rad

   !TODO hardcoded radius for now until I know where to get earth's radius from (MvD)
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

end module
!=========================================================================================
