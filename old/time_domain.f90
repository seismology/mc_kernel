!========================
  module time_domain
!========================
 
  use global_parameters
  use data_mesh
  use data_arrays

  implicit none

  public :: kxt_convolution
  private

  contains 

!----------------------------------------------------------------------------------------------------------
subroutine kxt_convolution(it,iwin)

  use input_output
  use operations
  use coord_trafo
  use misfit_model_param

  integer, intent(in) :: it,iwin
  integer :: jt,i
  character(len=4) ::  appidump2,appidump1

    if (do_lam) lamkernel = zero; if (do_mu) mukernel = zero
     if (do_rho) rhokernel = zero; if (do_vp) vpkernel = zero
     if (do_vs) vskernel = zero; if (do_imped) impkernel = zero
 allocate(vrot2(1:npts,1:3))
!  ===============
     do jt=1,it ! The time in the convolution integral
!  ===============

        call define_io_appendix(appidump1,jt)       ! fwd wavefield time index
        if (it-jt+1+ibwd_shift> ndumps) then
          call define_io_appendix(appidump2,ndumps) ! bwdwavefield time index
        else
          call define_io_appendix(appidump2,it-jt+1+ibwd_shift) 
        endif

! lambda kernel ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if (do_lam) then
       call read_summedstraintrace_src(usrc,appidump1,dir_fwdmesh1(1),jt) !!!!!!!!! NEED TO FIX FOR MIJ !!!!!!!!!!!!
       call read_summedstraintrace_rec(urec,appidump2,dir_bwdmesh,it-jt+1+ibwd_shift)

       ! multiply the two wavefields and add to kernel
       call collocate_sum_1d(usrc,urec,lamkernel,npts)
!af
       ! saving snaps of forward/backward wavefields
       if( save_snaps .and. mod(jt,10)==0 ) then
          call define_io_appendix(appidump2,jt)
          if (mynum==0) write(6,*)'saving lambda snaps...'
          call save_snapshot(npts,real(usrc),real(urec),appidump2,'lam')
       endif
endif
! end lambda kernel ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! rho kernel ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if (do_rho) then
       call read_velocities_src(vsrc,appidump1,dir_fwdmesh1(1))
       call rotate_vec_cyl2cart(vsrc,vrot2,dim_fwd,azim1_fwd,azim2_fwd)

       call read_velocities_rec(vrec,appidump2,dir_bwdmesh)

      ! rotate bwd velocity components 
       call rotate_vec_cyl2cart(vrec,vrot,dim_bwd,azim1_bwd,azim2_bwd)
       call rotate_velocity_components(vrot)

      ! multiply the two wavefields and add to kernel
      do i=1,ndim
       call collocate_sum_1d(vrot2(:,i),vrot(:,i),rhokernel,npts)
      enddo
          ! saving snaps of forward/backward wavefields
          ! need to think about what kind of wavefield snapshots make sense
endif
! end rhokernel ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


! mu kernel ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!af
if (do_mu) then 

!!!!!! NOT DONE YET !!!!!!!!
stop
!       call read_deviatoricstrain_src(usrc,appidump1)
!       call read_deviatoricstrain_rec(urec,appidump2)
!
!       do i = 1, dim_fwd ! needs to be a different integer.....

          ! rotate bwd strain deviator elements 
!     call rotate_tens_cyl2cart(devstrain,devstrain_rot,azim1_fwd,azim2_fwd,dim_fwd*2)
!          call rotate_strain_components(devstrain)

          ! multiply the two wavefields and add to kernel
 !         call collocate_sum_1d(usrc,urec,mukernel,npts)

       ! saving snaps of forward/backward wavefields
       ! need to think about what kind of wavefield snapshots make sense

!       enddo
endif
! end mu kernel ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!  ===============
     end do 
!  ===============

! usual data output----------------------

!  if (1+ibwd_shift> ndumps) then
   if (ndumps-it+1+ibwd_shift> ndumps) then
    call define_io_appendix(appidump2,ndumps) ! bwdwavefield time index
   else
!   call define_io_appendix(appidump2,it-jt+1+ibwd_shift) 
    call define_io_appendix(appidump2,ndumps-it+1+ibwd_shift) 
   endif
   call read_summedstraintrace_src(usrc,appidump1,dir_fwdmesh1(1),ndumps)
   call read_summedstraintrace_rec(urec,appidump2,dir_bwdmesh,ndumps-it+1+ibwd_shift)

end subroutine kxt_convolution
!----------------------------------------------------------------------------------------------------------

!========================
end module time_domain
!========================
